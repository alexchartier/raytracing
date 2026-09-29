"""Accepted-return ionogram loading and scoring for the current retrievals."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class Ionogram:
    frequencies: np.ndarray
    records: np.ndarray
    density_scale: float | None
    hmf2_shift_km: float | None
    settings: tuple
    f2_width_scale: float = 1.0
    topside_width_ratio: float = 1.0

    @classmethod
    def read(cls, path: Path) -> "Ionogram":
        with np.load(path, allow_pickle=False) as data:
            frequencies = np.asarray(data["frequencies_mhz"], dtype=float)
            records = np.asarray(data["records"], dtype=float)
            density = float(data["density_scale"]) if "density_scale" in data else None
            shift = float(data["hmf2_shift_km"]) if "hmf2_shift_km" in data else 0.0
            width = float(data["f2_width_scale"]) if "f2_width_scale" in data else 1.0
            top_ratio = float(data["topside_width_ratio"]) if "topside_width_ratio" in data else 1.0
            settings = (str(data["method"]), str(data["vertical_fan_layout"]),
                        float(data["vertical_outer_ray_fraction"]),
                        int(data["vertical_guard_seed_limit"]),
                        float(data["homing_tolerance_m"]))
        if records.ndim != 2 or records.shape[1] < 3:
            raise ValueError(f"Bad records in {path}")
        if len(frequencies) != 81 or not np.allclose(frequencies, np.arange(2, 10.0001, .1)):
            raise ValueError(f"Expected the complete 2–10 MHz, 100 kHz sweep in {path}")
        if len(records) and (np.min(records[:, 0]) < 0 or np.max(records[:, 0]) >= len(frequencies)):
            raise ValueError(f"Bad frequency indices in {path}")
        if settings not in (("adaptive", "equal_area_guarded", .5, 4, 1000.0),
                            ("adaptive_with_dense_gap_recovery", "equal_area_guarded",
                             .5, 4, 1000.0),
                            ("adaptive_with_dense_gap_recovery_and_20khz_nose_continuation",
                             "equal_area_guarded", .5, 4, 1000.0),
                            ("adaptive_with_dense_gap_recovery_and_20khz_nose_continuation_and_above_nose_probes",
                             "equal_area_guarded", .5, 4, 1000.0),
                            ("oblique_adaptive", "along_track", 1.0, -1, 1000.0)):
            raise ValueError(f"Unexpected generator settings or homing gate in {path}: {settings}")
        return cls(frequencies, records, density, shift, settings, width, top_ratio)

    def nose(self, mode: int) -> float | None:
        indices = self.records[self.records[:, 1] == mode, 0].astype(int)
        return float(self.frequencies[np.max(indices)]) if len(indices) else None

    def ridge(self, mode: int) -> dict[int, float]:
        records = self.records[self.records[:, 1] == mode]
        return {int(i): float(np.median(records[records[:, 0] == i, 2]))
                for i in np.unique(records[:, 0]).astype(int)}


def _return_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Symmetric clipped distance; all accepted returns enter the score."""
    if not len(a) and not len(b):
        return 0.0
    if not len(a) or not len(b):
        return 1.0
    # Frequency uses 0.2 MHz and group range uses 25 km. A clipped cost
    # keeps an occasional missed homing frequency from dominating a ridge.
    a_xy = np.column_stack((a[:, 0] / 2.0, a[:, 2] / 25.0))
    b_xy = np.column_stack((b[:, 0] / 2.0, b[:, 2] / 25.0))
    d_ab = cKDTree(b_xy).query(a_xy)[0]
    d_ba = cKDTree(a_xy).query(b_xy)[0]
    return float((np.mean(np.minimum(d_ab / 3.0, 1.0))
                  + np.mean(np.minimum(d_ba / 3.0, 1.0))) / 2.0)


def score(observed: Ionogram, predicted: Ionogram) -> dict[str, float]:
    if not np.allclose(observed.frequencies, predicted.frequencies, atol=1e-6):
        raise ValueError("Frequency axes differ")
    if observed.settings != predicted.settings:
        raise ValueError("Ionogram generator settings differ")
    distances = []
    noses = []
    ridges = []
    for mode in (1, -1):
        obs = observed.records[observed.records[:, 1] == mode]
        pred = predicted.records[predicted.records[:, 1] == mode]
        distances.append(_return_distance(obs, pred))
        on, pn = observed.nose(mode), predicted.nose(mode)
        if on is None or pn is None:
            noses.append(1.0)
        elif on >= observed.frequencies[-1] or pn >= observed.frequencies[-1]:
            # A return at the sweep edge gives only a lower bound on the nose.
            noses.append(min(abs(min(on, 10.0) - min(pn, 10.0)) / .8, 1.0)
                         if (on < 10.0 or pn < 10.0) else 0.0)
        else:
            noses.append(min(abs(on - pn) / .8, 1.0))
        oridge, pridge = observed.ridge(mode), predicted.ridge(mode)
        common = sorted(set(oridge) & set(pridge))
        if common:
            residual = np.array([oridge[i] - pridge[i] for i in common])
            ridges.append(float(np.mean(np.minimum(np.abs(residual) / 60.0, 1.0))))
        else:
            ridges.append(0.0 if not oridge and not pridge else 1.0)
    return {
        "return_distance": float(np.mean(distances)),
        "nose": float(np.mean(noses)),
        "ridge": float(np.mean(ridges)),
        "total": float(.55 * np.mean(distances) + .25 * np.mean(noses)
                       + .20 * np.mean(ridges)),
    }
