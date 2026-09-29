"""Generate a synthetic generalized-Chapman truth grid independently of retrieval code."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports/data/chapman_model_validation"
LAT = np.arange(-80.0, -28.0 + 0.1, 2.0)
LON = np.arange(-24.0, 40.0 + 0.1, 4.0)
ALT = np.arange(60.0, 900.0 + 0.1, 20.0)
PLASMA_MHZ_PER_SQRT_CM3 = 0.00898
PARAMETERS = {
    "fof2_mhz": 5.7,
    "hmf2_km": 300.0,
    "bottomside_scale_km": 70.0,
    "topside_scale_km": 95.0,
    "tail_curvature": 0.45,
}


def analytic_density(altitude_km: np.ndarray) -> np.ndarray:
    p = PARAMETERS
    distance = np.asarray(altitude_km, dtype=float) - p["hmf2_km"]
    scale = np.where(distance < 0, p["bottomside_scale_km"],
                     p["topside_scale_km"])
    normalized_height = distance / scale
    positive = np.maximum(normalized_height, 0.0)
    modified_height = np.where(
        distance < 0, normalized_height,
        normalized_height / (1 + p["tail_curvature"] * positive / (positive + 3)))
    chapman_exponent = (1 - modified_height - np.exp(-modified_height)) / 2
    return (p["fof2_mhz"] / PLASMA_MHZ_PER_SQRT_CM3)**2 * np.exp(chapman_exponent)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    profile = analytic_density(ALT)
    grid = np.broadcast_to(profile, (len(LAT), len(LON), len(ALT))).copy()
    destination = OUT / "truth_density.npz"
    np.savez_compressed(destination, latitudes_deg=LAT, longitudes_deg=LON,
                        altitudes_km=ALT, electron_density_cm3=grid,
                        model=np.array("analytic generalized Chapman test truth"),
                        time_utc=np.array("2010-01-01T12:00:00"))
    (OUT / "truth_parameters.json").write_text(json.dumps({
        "source": "analytic generalized Chapman",
        **PARAMETERS,
        "horizontal_field": "uniform",
    }, indent=2) + "\n")
    print(destination)


if __name__ == "__main__":
    main()
