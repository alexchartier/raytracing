"""Build IRI oblique upper-tail candidates from ionograms and 800 km samples.

This is an observable-only candidate proposal. The optional ``diagnose`` action
opens the withheld IRI density *after* all proposals have been written; it is a
post hoc representation check, not candidate selection or an ionogram result.
Full 2–10 MHz O/X rays must be traced before claiming retrieval improvement.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from oblique_model_family_cases import (  # noqa: E402
    case_root, load_ionosphere_grid_netcdf, local_observations,
    save_private_grid, save_json,
)
from python_raytrace.monotone_topside import monotone_topside_grid  # noqa: E402

CASE = "iri2016"
SOURCE = "q2m06_gq2p08"
TRACK_FRACTION = .82
PROPOSALS = {
    "ut_h20": (20.0, 0.0, 0.0, None),
    "ut_h20_g075": (20.0, .75, 0.0, None),
    "ut_h20_g075_q3p03": (20.0, .75, 0.0, .03),
    "ut_h20_g075_q3p06": (20.0, .75, 0.0, .06),
    "ut_h20_g075_q3p09": (20.0, .75, 0.0, .09),
    "ut_h20_g075_q2p04_q3p06": (20.0, .75, .04, .06),
    "ut_h20_g075_q2p08_q3p06": (20.0, .75, .08, .06),
    "ut_h20_g075_q2p12_q3p06": (20.0, .75, .12, .06),
    "ut_h20_g075_q2p16_q3p06": (20.0, .75, .16, .06),
    "ut_h20_g075_q2flat08_q3p06": (20.0, .75, .08, .06, 0.0),
    "ut_h20_g075_q2flat12_q3p06": (20.0, .75, .12, .06, 0.0),
    "ut_h20_g075_q2flat16_q3p06": (20.0, .75, .16, .06, 0.0),
    "ut_h20_g075_q2flat12_q3p00": (20.0, .75, .12, .00, 0.0),
    "ut_h20_g075_q2flat12_q3p03": (20.0, .75, .12, .03, 0.0),
    "ut_h20_g075_q2flat16_q3p00": (20.0, .75, .16, .00, 0.0),
    "ut_h20_g075_q2flat16_q3p03": (20.0, .75, .16, .03, 0.0),
    "ut_h10_g075_q2p04_q3p06": (10.0, .75, .04, .06),
    "ut_h30_g075_q2p04_q3p06": (30.0, .75, .04, .06),
    "ut_h20_g100_q2p04_q3p06": (20.0, 1.0, .04, .06),
}


def implied_upper_fraction(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    return np.array([float(PchipInterpolator(
        [0.0, .25, .60, 1.0], [1.0, a, b, 0.0])(TRACK_FRACTION))
        for a, b in zip(q1, q2)])


def build() -> None:
    root = case_root(CASE)
    previous = json.loads((root / "selection.json").read_text())
    if previous["chosen_candidate"] != SOURCE or previous["uses_full_truth_density"]:
        raise ValueError("Expected the frozen observable-only IRI start")
    source = json.loads((root / SOURCE / "parameters.json").read_text())
    start = json.loads((root / "start_parameters.json").read_text())
    latitudes = np.asarray(start["latitude_controls_deg"], dtype=float)
    q1 = np.asarray(source["q1_controls"], dtype=float)
    q2 = np.asarray(source["q2_controls"], dtype=float)
    local = local_observations(CASE)
    track_lat = np.asarray(local.lat_lon_deg)[:, 0]
    observed_800 = np.asarray(local.electron_density_cm3)
    at_controls = np.interp(latitudes, track_lat, observed_800)
    ratio = at_controls / at_controls[0]
    base = load_ionosphere_grid_netcdf(root / "unanchored_grid.nc")
    records = []
    for name, proposal in PROPOSALS.items():
        height, gradient, q2_change, q3_change = proposal[:4]
        q2_gradient_scale = proposal[4] if len(proposal) > 4 else 1.0
        trial_q2 = (q2 + q2_change if q2_gradient_scale == 1.0 else
                    q2.mean() + q2_change + q2_gradient_scale * (q2 - q2.mean()))
        trial_q3 = (None if q3_change is None else
                    implied_upper_fraction(q1, trial_q2) + q3_change)
        peak_scale = ratio ** gradient
        destination = root / name
        parameters = {
            "source": "frozen IRI oblique selection plus observable-based upper-tail proposal",
            "uses_full_truth_density": False,
            "source_candidate": SOURCE,
            "source_selection_sha256": previous["scores_sha256"],
            "height_shift_km": height,
            "peak_gradient_exponent_from_800km_density": gradient,
            "peak_density_scale_controls": peak_scale.tolist(),
            "q1_controls": q1.tolist(),
            "q2_controls": trial_q2.tolist(),
            "q3_fractional_height": TRACK_FRACTION if trial_q3 is not None else None,
            "q3_controls": trial_q3.tolist() if trial_q3 is not None else None,
            "800km_track_density_source": str((root / "insitu_1km.json").relative_to(ROOT)),
        }
        if q2_gradient_scale != 1.0:
            parameters["q2_gradient_scale_from_source"] = q2_gradient_scale
        grid_path = destination / "grid.nc"
        parameters_path = destination / "parameters.json"
        if grid_path.exists() or parameters_path.exists():
            if (not grid_path.is_file() or not parameters_path.is_file() or
                    json.loads(parameters_path.read_text()) != parameters):
                raise ValueError(f"Existing candidate differs from proposal: {name}")
        else:
            field = monotone_topside_grid(
                base, local, latitudes, q1, trial_q2,
                np.full(len(latitudes), height),
                fraction_at_eighty_two_percent=trial_q3,
                peak_density_scale=peak_scale)
            save_private_grid(grid_path, field)
            save_json(parameters_path, parameters)
        records.append({"name": name, "grid": str((destination / "grid.nc").relative_to(ROOT)),
                        "parameters": parameters})
    save_json(root / "upper_tail_design.json", {
        "case": CASE,
        "candidate_selection_uses_truth_density": False,
        "design_informed_by_prior_truth_diagnostics": True,
        "candidate_source": SOURCE,
        "rationale": "extra shape knot near 700 km; independent along-link upper-shape gradient; peak scale inferred from observed 800 km density; bounded peak-height and middle-topside checks",
        "frequency_sweep_mhz": [2.0, 10.0, .1],
        "candidates": records,
    })
    save_json(root / "upper_tail_trace_tasks.json", [
        {"case": CASE, "candidate": row["name"]} for row in records])
    print(json.dumps({"candidates_built": len(records),
                      "peak_scale_from_800km_track": ratio.tolist()}))


def diagnose() -> None:
    """Measure representation after building; never use this to choose rays."""
    from evaluate_general_field_sami import metrics, profiles

    root = case_root(CASE)
    design = json.loads((root / "upper_tail_design.json").read_text())
    with np.load(root / "forward/truth/ionogram_01.npz", allow_pickle=False) as ray:
        tx = float(ray["tx_lat_deg"])
        rx = float(ray["rx_lat_deg"])
        longitude = float(ray["tx_lon_deg"])
    locations = np.column_stack((np.array([tx, (tx + rx) / 2, rx]),
                                 np.full(3, longitude)))
    link_locations = np.column_stack((np.linspace(tx, rx, 21),
                                      np.full(21, longitude)))
    altitude, truth = profiles(root / "truth_grid.nc", locations)
    _, link_truth = profiles(root / "truth_grid.nc", link_locations)
    rows = []
    for name in ("start", SOURCE, *(row["name"] for row in design["candidates"])):
        _, field = profiles(root / name / "grid.nc", locations)
        _, link_field = profiles(root / name / "grid.nc", link_locations)
        evaluation = metrics(altitude, truth, field)
        link_evaluation = metrics(altitude, link_truth, link_field)
        mask = (altitude >= 500) & (altitude <= 700)
        upper_mid = 100 * np.sqrt(np.mean((field[:, mask] - truth[:, mask]) ** 2,
                                            axis=1)) / truth.max(axis=1)
        rows.append({"name": name, "density_truth_read_after_design": True,
                     "evaluation_positions": ["transmitter", "midpoint", "receiver"],
                     "topside_nrmse_350_to_800_percent": evaluation["topside_nrmse_350_to_800_percent"],
                     "link_21_position_topside_nrmse_350_to_800_percent":
                     link_evaluation["topside_nrmse_350_to_800_percent"],
                     "upper_500_to_700_nrmse_percent_by_position": upper_mid.tolist(),
                     "hmf2_mae_km": evaluation["hmf2_mae_km"],
                     "fof2_mae_mhz": evaluation["fof2_mae_mhz"]})
    save_json(root / "upper_tail_density_diagnostic.json", {
        "case": CASE,
        "note": "Post hoc density representation only; candidates lack new oblique ray traces and no retrieval improvement is established",
        "rows": rows,
    })
    for row in rows:
        print(row["name"], round(row["topside_nrmse_350_to_800_percent"], 3),
              round(float(np.mean(row["upper_500_to_700_nrmse_percent_by_position"])), 3),
              round(row["hmf2_mae_km"], 1))


def main() -> None:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("build", "diagnose"))
    args = parser.parse_args()
    if args.action == "build":
        build()
    else:
        diagnose()


if __name__ == "__main__":
    main()
