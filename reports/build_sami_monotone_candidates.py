"""Build pass-varying monotone topside candidates from ionograms and in-situ data."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from ionogram_metrics import Ionogram  # noqa: E402
from python_raytrace.general_field_inverse import LocalDensity  # noqa: E402
from python_raytrace.grid import (  # noqa: E402
    load_ionosphere_grid_netcdf, save_ionosphere_grid_netcdf,
)
from python_raytrace.monotone_topside import monotone_topside_grid  # noqa: E402

CASE = ROOT / "reports/data/sami3_wave_20170111_0600"
OLD = ROOT / "reports/data/general_field_sami3_wave"
RUN = ROOT / "reports/data/sami3_monotone_topside"
INDICES = (1, 5, 9, 12, 14, 17, 20)


def ridge_need(index: int) -> dict:
    observed = Ionogram.read(CASE / "vertical_truth" / f"ionogram_{index:02d}.npz")
    modeled = Ionogram.read(OLD / "forward/nose_init/ionograms" /
                           f"ionogram_{index:02d}.npz")
    modes = {}
    all_residuals = []
    for mode, label in ((1, "O"), (-1, "X")):
        a, b = observed.ridge(mode), modeled.ridge(mode)
        cutoff = min(observed.nose(mode), modeled.nose(mode)) - .2
        common = [i for i in a.keys() & b.keys()
                  if 3.0 <= observed.frequencies[i] <= min(5.0, cutoff)]
        residual = [a[i] - b[i] for i in common]
        if len(residual) < 8:
            raise ValueError(f"Insufficient {label}-mode ridge at profile {index}")
        modes[label] = {"points": len(residual),
                        "median_observed_minus_modeled_km": float(np.median(residual))}
        all_residuals.append(float(np.median(residual)))
    return {"index": index, "modes": modes,
            "mean_mode_ridge_need_km": float(np.mean(all_residuals))}


def main() -> None:
    os.umask(0o077)
    RUN.mkdir(mode=0o700, parents=True, exist_ok=True)
    insitu = json.loads((RUN / "insitu_1km.json").read_text())
    local = LocalDensity(np.asarray(insitu["lat_lon_deg"]),
                         insitu["altitude_km"],
                         np.asarray(insitu["density_cm3"]),
                         insitu["relative_uncertainty"])
    profiles = json.loads((CASE / "manifest.json").read_text())["profiles"]
    sounding = [profiles[index - 1] for index in INDICES]
    latitudes = np.array([row["latitude_deg"] for row in sounding])
    base = load_ionosphere_grid_netcdf(OLD / "nose_init/grid.nc")
    needs = [ridge_need(index) for index in INDICES]
    residual = np.array([row["mean_mode_ridge_need_km"] for row in needs])
    height = np.clip(-18.0 - .08 * (residual - 170.0), -30.0, -10.0)
    candidates = []
    for name, gain, height_gain in (
            ("varying_g04", .4, 1.0),
            ("varying_g07", .7, 1.0),
            ("varying_g10", 1.0, 1.0),
            ("varying_g13", 1.3, 1.0),
            ("varying_g10_h0", 1.0, 0.0)):
        q1 = np.clip(.65 + gain * .0007 * residual, .68, .87)
        q2 = np.clip(.15 + gain * .0010 * residual, .20, .49)
        if np.any(q1 <= q2 + .15):
            raise ValueError("Topside shape fractions are too close")
        trial = monotone_topside_grid(
            base, local, latitudes, q1, q2, height_gain * height)
        output = RUN / name
        output.mkdir(mode=0o700, parents=True, exist_ok=True)
        path = output / "grid.nc"
        save_ionosphere_grid_netcdf(path, trial)
        path.chmod(0o600)
        candidates.append({"name": name, "grid": str(path.relative_to(ROOT)),
                           "fraction_at_quarter_by_profile": q1.tolist(),
                           "fraction_at_three_fifths_by_profile": q2.tolist(),
                           "peak_height_shift_km_by_profile":
                               (height_gain * height).tolist()})
    plan = {
        "selection_uses_truth_density": False,
        "source_prior_grid": str((OLD / "nose_init/grid.nc").relative_to(ROOT)),
        "source_ionograms": str((OLD / "forward/nose_init/ionograms").relative_to(ROOT)),
        "observed_ionograms": str((CASE / "vertical_truth").relative_to(ROOT)),
        "in_situ_observations": str((RUN / "insitu_1km.json").relative_to(ROOT)),
        "manifest": str((CASE / "manifest.json").relative_to(ROOT)),
        "observed_profile_indices": list(INDICES),
        "pilot_profile_indices": [9, 14, 17],
        "ridge_need_by_profile": needs,
        "candidates": candidates,
        "ionogram_method": "raw adaptive O/X, option D, 2–10 MHz at 100 kHz",
    }
    path = RUN / "plan.json"
    path.write_text(json.dumps(plan, indent=2) + "\n")
    path.chmod(0o600)
    print(json.dumps({"ridge_need_km": residual.tolist(),
                      "candidate_count": len(candidates)}, indent=2))


if __name__ == "__main__":
    main()
