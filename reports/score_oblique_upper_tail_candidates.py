"""Select upper-topside proposals using only oblique observables and in-situ data.

The 2.0--3.0 MHz range receives an explicit extra weight because its rays
probe the upper topside. Do not open the withheld density grid here. The
proposal design was informed by earlier post hoc IRI density diagnostics;
this remains an iterative development case, not a blind validation case.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from oblique_model_family_cases import ROOT, case_root, save_json
from oblique_wave_pass import doppler_score
from sami3_monotone_oblique_score import verify_geometry
from sami3_monotone_score import score_pair
from score_oblique_model_families import insitu_error, low_frequency_range_error

CASE = "iri2016"
SOURCE = "q2m06_gq2p08"
LOW_BAND_MHZ = (2.0, 3.0)
LOW_BAND_WEIGHT = 1.0
DOPPLER_WEIGHT_KM = 20.0
INSITU_TOLERANCE_PERCENT = 0.5


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def score() -> None:
    os.umask(0o077)
    root = case_root(CASE)
    design_path = root / "upper_tail_design.json"
    design = json.loads(design_path.read_text())
    if design["case"] != CASE or design["candidate_selection_uses_truth_density"]:
        raise ValueError("Unexpected candidate design")
    names = [SOURCE, *(row["name"] for row in design["candidates"])]
    if len(names) != len(set(names)):
        raise ValueError("Repeated candidate name")
    truth = root / "forward/truth/ionogram_01.npz"
    rows = []
    for name in names:
        modeled = root / "forward" / name / "ionogram_01.npz"
        if not modeled.is_file():
            raise FileNotFoundError(f"Trace all candidates before selection: {modeled}")
        verify_geometry(truth, modeled, 1)
        grid = root / name / "grid.nc"
        grid_hash = sha256(grid)
        if name != SOURCE:
            with np.load(modeled, allow_pickle=False) as ray:
                if "grid_sha256" not in ray or str(ray["grid_sha256"]) != grid_hash:
                    raise ValueError(f"Ray output was not traced from {grid}")
        whole = score_pair(truth, modeled)
        low = low_frequency_range_error(truth, modeled, *LOW_BAND_MHZ)
        doppler = doppler_score(truth, modeled)
        insitu = insitu_error(CASE, grid)
        rows.append({
            "name": name,
            "grid_sha256": grid_hash,
            "ray_sha256": sha256(modeled),
            "whole_ionogram_equivalent_km": whole["total_equivalent_km"],
            "upper_topside_2_to_3_mhz_range_error_km": low,
            "doppler_1hz_mismatch": doppler,
            "insitu_max_absolute_error_percent": insitu,
            "combined_equivalent_km": whole["total_equivalent_km"]
            + LOW_BAND_WEIGHT * low + DOPPLER_WEIGHT_KM * doppler,
        })
    eligible = [row for row in rows
                if row["insitu_max_absolute_error_percent"] <= INSITU_TOLERANCE_PERCENT]
    if not eligible:
        raise ValueError("No candidate satisfies the 800 km track constraint")
    chosen = min(eligible, key=lambda row: row["combined_equivalent_km"])
    output = root / "upper_tail_scores.json"
    save_json(output, {
        "case": CASE,
        "uses_full_truth_density": False,
        "design_informed_by_prior_truth_diagnostics": True,
        "design_sha256": sha256(design_path),
        "truth_ionogram": str(truth.relative_to(ROOT)),
        "score": "whole O/X accepted-return and nose score + 1x symmetric 2.0--3.0 MHz range score + 20 km times Doppler 1 Hz mismatch",
        "low_band_weight": LOW_BAND_WEIGHT,
        "doppler_weight_km": DOPPLER_WEIGHT_KM,
        "insitu_tolerance_percent": INSITU_TOLERANCE_PERCENT,
        "candidates": rows,
        "selected": chosen["name"],
    })
    save_json(root / "upper_tail_selection.json", {
        "case": CASE,
        "selected": chosen["name"],
        "uses_full_truth_density": False,
        "design_informed_by_prior_truth_diagnostics": True,
        "scores_sha256": sha256(output),
        "selected_grid_sha256": chosen["grid_sha256"],
        "selected_ray_sha256": chosen["ray_sha256"],
    })
    print(json.dumps({"selected": chosen["name"],
                      "combined_equivalent_km": chosen["combined_equivalent_km"]}))


def evaluate() -> None:
    """Read full IRI density only after checking the saved observable selection."""
    from evaluate_general_field_sami import metrics, profiles

    root = case_root(CASE)
    scores_path = root / "upper_tail_scores.json"
    scores = json.loads(scores_path.read_text())
    selection = json.loads((root / "upper_tail_selection.json").read_text())
    design_path = root / "upper_tail_design.json"
    chosen = selection["selected"]
    selected_grid = root / chosen / "grid.nc"
    selected_ray = root / "forward" / chosen / "ionogram_01.npz"
    if (scores["uses_full_truth_density"] or selection["uses_full_truth_density"]
            or scores["design_sha256"] != sha256(design_path)
            or selection["scores_sha256"] != sha256(scores_path)
            or selection["selected_grid_sha256"] != sha256(selected_grid)
            or selection["selected_ray_sha256"] != sha256(selected_ray)
            or chosen != scores["selected"]):
        raise ValueError("The observable selection is not frozen")
    with np.load(root / "forward/truth/ionogram_01.npz", allow_pickle=False) as ray:
        transmitter = float(ray["tx_lat_deg"])
        receiver = float(ray["rx_lat_deg"])
        longitude = float(ray["tx_lon_deg"])
    locations = np.column_stack((np.linspace(transmitter, receiver, 21),
                                 np.full(21, longitude)))
    altitude, truth = profiles(root / "truth_grid.nc", locations)
    rows = []
    for name in (SOURCE, chosen):
        _, field = profiles(root / name / "grid.nc", locations)
        quality = metrics(altitude, truth, field)
        upper = (altitude >= 500) & (altitude <= 700)
        relative_rms = 100 * np.sqrt(np.mean((field[:, upper] / truth[:, upper] - 1) ** 2))
        observable = next(row for row in scores["candidates"] if row["name"] == name)
        rows.append({"name": name,
                     "topside_nrmse_350_to_800_percent":
                     quality["topside_nrmse_350_to_800_percent"],
                     "relative_rms_500_to_700_percent": float(relative_rms),
                     "fof2_mae_mhz": quality["fof2_mae_mhz"],
                     "hmf2_mae_km": quality["hmf2_mae_km"],
                     "ionogram_combined_equivalent_km": observable["combined_equivalent_km"],
                     "ionogram_low_2_to_3_mhz_error_km":
                     observable["upper_topside_2_to_3_mhz_range_error_km"]})
    save_json(root / "upper_tail_evaluation.json", {
        "case": CASE,
        "selected": chosen,
        "positions_along_600km_link": 21,
        "selection_uses_full_truth_density": False,
        "design_informed_by_prior_truth_diagnostics": True,
        "full_truth_density_read_for_this_evaluation": True,
        "selection_sha256": sha256(root / "upper_tail_selection.json"),
        "rows": rows,
    })
    print(json.dumps({"selected": chosen, "rows": rows}, indent=2))


def main() -> None:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", nargs="?", choices=("score", "evaluate"),
                        default="score")
    args = parser.parse_args()
    if args.action == "evaluate":
        evaluate()
    else:
        score()


if __name__ == "__main__":
    main()
