"""Check the model-free spline against analytic shapes excluded from its basis.

These are idealized scalar ranges, not a substitute for a full O/X ray check.
No NeQuick or IRI profile is loaded by this script.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from python_raytrace.general_topside_inverse import (
    PLASMA_MHZ_PER_SQRT_CM3, fit_spline_ionogram, vertical_group_range_km,
)


@dataclass(frozen=True)
class AnalyticF2:
    fof2_mhz: float
    hmf2_km: float
    shape: str

    def density_cm3(self, altitudes_km: np.ndarray | float) -> np.ndarray:
        u = np.asarray(altitudes_km, dtype=float) - self.hmf2_km
        bottom = 0.5 * (1 - u / 55 - np.exp(np.clip(-u / 55, -80, 80)))
        if self.shape == "epstein":
            top = -2 * np.log(np.cosh(np.maximum(u, 0) / 150))
        elif self.shape == "rational_gaussian":
            v = np.maximum(u, 0)
            top = -0.5 * (v / 90) ** 2 / (1 + v / 160)
        else:
            raise ValueError(self.shape)
        return (self.fof2_mhz / PLASMA_MHZ_PER_SQRT_CM3) ** 2 * np.exp(
            np.where(u >= 0, top, bottom))


def main() -> None:
    frequency_axis = np.arange(2, 10.0001, 0.1)
    result = {}
    for shape, fo, hm in (("epstein", 5.5, 265.0),
                          ("rational_gaussian", 6.2, 310.0)):
        truth = AnalyticF2(fo, hm, shape)
        rows = []
        for index, frequency in enumerate(frequency_axis):
            if frequency >= fo - 0.04:
                continue
            group_range = vertical_group_range_km(truth, frequency, 800.0)
            if np.isfinite(group_range):
                rows.append([index, 1, group_range])
        z = np.arange(hm, 600.1)
        source = truth.density_cm3(z)
        result[shape] = {"truth_fof2_mhz": fo, "truth_hmf2_km": hm}
        for name, regularize_tail in (("legacy", False), ("regularized", True)):
            fitted = fit_spline_ionogram(np.asarray(rows), frequency_axis, 800.0,
                                         regularize_tail=regularize_tail)
            retrieved = fitted.layer.density_cm3(z)
            result[shape][name] = {
                "fof2_error_mhz": fitted.layer.fof2_mhz - fo,
                "hmf2_error_km": fitted.layer.hmf2_km - hm,
                "topside_peak_normalized_rms_percent": float(
                    100 * np.sqrt(np.mean((retrieved - source) ** 2)) / source.max()),
                "range_mae_km": fitted.range_mae_km[1],
                "fitted_log_slopes_per_km": fitted.layer.log_slope_per_km.tolist(),
            }
    path = Path("reports/data/regularized_spline_scalar_validation.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
