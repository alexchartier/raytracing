"""Score 600 km oblique O/X ionograms with uncapped range and 1 Hz Doppler."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from oblique_wave_pass import doppler_score  # noqa: E402
from sami3_monotone_score import score_pair  # noqa: E402

CASE = ROOT / "reports/data/sami3_wave_20170111_0600/oblique_wave_600km"
RUN = ROOT / "reports/data/sami3_monotone_oblique"
OLD_FIT = "gain1_sigma60_height0.15_spline"


def verify_geometry(truth: Path, modeled: Path, index: int) -> None:
    keys = ("tx_lat_deg", "tx_lon_deg", "tx_alt_km", "rx_lat_deg",
            "rx_lon_deg", "rx_alt_km", "satellite_separation_km",
            "spacecraft_speed_mps", "fan_launch_directions")
    with np.load(truth, allow_pickle=False) as observed, np.load(modeled, allow_pickle=False) as fitted:
        if int(observed["profile_index"]) != index or int(fitted["profile_index"]) != index:
            raise ValueError(f"Wrong link index at {truth} or {modeled}")
        for key in keys:
            if not np.isclose(float(observed[key]), float(fitted[key]), atol=1e-6):
                raise ValueError(f"Different oblique geometry: {key} at link {index}")
        if not np.isclose(float(observed["satellite_separation_km"]), 600.0, atol=.01):
            raise ValueError(f"Expected 600 km link at {index}")
        if int(observed["fan_launch_directions"]) != 75:
            raise ValueError(f"Expected the same 75-ray oblique fan at {index}")


def ionograms(name: str) -> Path:
    if name == "old_baseline":
        return CASE / "baseline"
    if name == "old_fit":
        return CASE / OLD_FIT / "ionograms"
    return RUN / "forward" / name / "ionograms"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", nargs="+", type=int, required=True)
    parser.add_argument("--candidates", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    os.umask(0o077)
    rows = []
    for name in args.candidates:
        directory = ionograms(name)
        parts = []
        for index in args.indices:
            truth = CASE / "truth" / f"ionogram_{index:02d}.npz"
            modeled = directory / f"ionogram_{index:02d}.npz"
            verify_geometry(truth, modeled, index)
            parts.append({"index": index, **score_pair(truth, modeled),
                          "doppler_1hz_score": doppler_score(truth, modeled)})
        rows.append({"name": name,
                     "mean_total_equivalent_km": float(np.mean([
                         item["total_equivalent_km"] for item in parts])),
                     "mean_doppler_1hz_score": float(np.mean([
                         item["doppler_1hz_score"] for item in parts])),
                     "by_profile": parts})
    result = {
        "profile_indices": args.indices,
        "satellite_separation_km": 600.0,
        "selection_uses_truth_density": False,
        "ionogram_score": "uncapped symmetric nearest O/X group-range distance across all accepted returns, 150 km missing-frequency penalty, 100 km/MHz mean O/X nose penalty",
        "doppler_score": "1 Hz quantization, nearest O/X return within 40 km, 20 Hz saturated error; reported separately",
        "candidates": rows,
        "lowest_ionogram_score": min(rows, key=lambda row: row["mean_total_equivalent_km"])["name"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.output.chmod(0o600)
    print(json.dumps({"scores_km": {row["name"]: row["mean_total_equivalent_km"]
                                    for row in rows},
                      "doppler_scores": {row["name"]: row["mean_doppler_1hz_score"]
                                         for row in rows}}, indent=2))


if __name__ == "__main__":
    main()
