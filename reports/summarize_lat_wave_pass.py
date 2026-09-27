"""Validate and summarize all 20 latitude-wave truth ionograms."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from fit_d_ionogram import Ionogram


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "reports/data/lat_wave_pass_manifest.json"
BASELINE = ROOT / "reports/data/d_inverse_iri_truth_density.npz"
WAVE = ROOT / "reports/data/lat_wave_pass_truth_density.npz"
RAW_RESULTS = ROOT / "reports/data/lat_wave_pass_ionograms"
SUMMARY = ROOT / "reports/data/lat_wave_pass_summary.json"
FIGURE = ROOT / "reports/figures/lat_wave_pass_observables.png"


def summarize(results: Path) -> None:
    manifest = json.loads(MANIFEST.read_text())
    positions = manifest["profiles"]
    with np.load(BASELINE, allow_pickle=False) as data:
        axes = (data["latitudes_deg"], data["longitudes_deg"])
        baseline = data["electron_density_cm3"]
    with np.load(WAVE, allow_pickle=False) as data:
        wave = data["electron_density_cm3"]
    points = np.array([[row["latitude_deg"], row["longitude_deg"]] for row in positions])
    baseline_profiles = RegularGridInterpolator(axes, baseline)(points)
    wave_profiles = RegularGridInterpolator(axes, wave)(points)
    rows = []
    for row, base_profile, wave_profile in zip(positions, baseline_profiles, wave_profiles):
        path = results / f"ionogram_{row['index']:02d}.npz"
        ionogram = Ionogram.read(path)
        with np.load(path, allow_pickle=False) as data:
            actual = np.array([data["tx_lat_deg"], data["tx_lon_deg"], data["tx_alt_km"]], dtype=float)
            if int(data["profile_index"]) != row["index"] or not np.allclose(
                    actual, [row["latitude_deg"], row["longitude_deg"], row["altitude_km"]],
                    rtol=0.0, atol=1e-8):
                raise ValueError(f"Profile metadata does not match the manifest: {path}")
            runtime = float(data["runtime_seconds"])
            recovered = int(data["dense_recovered_return_count"]) if "dense_recovered_return_count" in data else 0
            checks = int(data["dense_recovery_frequency_checks"]) if "dense_recovery_frequency_checks" in data else 0
            recovery_runtime = float(data["dense_recovery_runtime_seconds"]) if "dense_recovery_runtime_seconds" in data else 0.0
        raw = Ionogram.read(RAW_RESULTS / f"ionogram_{row['index']:02d}.npz")
        rows.append({
            **row,
            "baseline_peak_cm3": float(np.max(base_profile)),
            "truth_peak_cm3": float(np.max(wave_profile)),
            "truth_peak_plasma_frequency_mhz": float(0.00898 * np.sqrt(np.max(wave_profile))),
            "o_nose_mhz": ionogram.nose(1),
            "x_nose_mhz": ionogram.nose(-1),
            "o_returns": int(np.count_nonzero(ionogram.records[:, 1] == 1)),
            "x_returns": int(np.count_nonzero(ionogram.records[:, 1] == -1)),
            "runtime_seconds": runtime,
            "dense_recovery_runtime_seconds": recovery_runtime,
            "dense_recovery_checked_bins": checks,
            "dense_recovered_returns": recovered,
            "raw_o_nose_mhz": raw.nose(1),
            "raw_x_nose_mhz": raw.nose(-1),
            "file": str(path.resolve().relative_to(ROOT)),
        })
    summary = {
        "source": "Independent IRI-2016 density plus imposed wave; shared forward ray tracer",
        "density_frozen_at_utc": manifest["model_epoch_utc"],
        "wave": manifest["wave"],
        "profile_count": len(rows),
        "total_accepted_returns": sum(row["o_returns"] + row["x_returns"] for row in rows),
        "total_dense_recovered_returns": sum(row["dense_recovered_returns"] for row in rows),
        "total_dense_checked_bins": sum(row["dense_recovery_checked_bins"] for row in rows),
        "profiles": rows,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    latitudes = np.array([row["latitude_deg"] for row in rows])
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    axes[0].plot(latitudes, [row["truth_peak_plasma_frequency_mhz"] for row in rows],
                 color="#555555", linewidth=2, label="Plasma frequency at density peak")
    if summary["total_dense_recovered_returns"]:
        axes[0].plot(latitudes, [row["raw_o_nose_mhz"] for row in rows], ":",
                     color="#94afd8", linewidth=1.5, label="O limit before gap recovery")
        axes[0].plot(latitudes, [row["raw_x_nose_mhz"] for row in rows], ":",
                     color="#e49b95", linewidth=1.5, label="X limit before gap recovery")
    axes[0].plot(latitudes, [row["o_nose_mhz"] for row in rows], "o-", color="#2455a4",
                 markersize=4, label="O accepted-return limit")
    axes[0].plot(latitudes, [row["x_nose_mhz"] for row in rows], "o-", color="#c43b31",
                 markersize=4, label="X accepted-return limit")
    axes[0].set(xlabel="Latitude (degrees)", ylabel="Frequency (MHz)",
                title="Mode-resolved limits after gap recovery")
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    axes[1].plot(latitudes, [row["o_returns"] for row in rows], "o-", color="#2455a4",
                 markersize=4, label="O returns")
    axes[1].plot(latitudes, [row["x_returns"] for row in rows], "o-", color="#c43b31",
                 markersize=4, label="X returns")
    axes[1].set(xlabel="Latitude (degrees)", ylabel="Accepted returns",
                title="Return coverage along the pass")
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=220)
    plt.close(fig)
    print(json.dumps({"profile_count": summary["profile_count"],
                      "total_accepted_returns": summary["total_accepted_returns"],
                      "total_dense_recovered_returns": summary["total_dense_recovered_returns"],
                      "summary": str(SUMMARY), "figure": str(FIGURE)}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    summarize(parser.parse_args().results)


if __name__ == "__main__":
    main()
