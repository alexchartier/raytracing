"""Compare the original and coverage-qualified IRI retrievals with truth.

This figure is a post-fit diagnostic. Truth density is never used for candidate
selection; the two candidate parameter sets come from the summary JSON files.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_d_inverse_fit import _fitted_profile, _profile


def plot_progress(background: Path, truth_grid: Path, original_summary: Path,
                  improved_summary: Path, output: Path) -> None:
    original = json.loads(original_summary.read_text())
    improved = json.loads(improved_summary.read_text())
    location = improved["density_at_sounder"]
    latitude, longitude = location["latitude_deg"], location["longitude_deg"]
    altitudes, truth = _profile(truth_grid, "electron_density_cm3", latitude, longitude)
    old_altitudes, old = _fitted_profile(background, latitude, longitude,
                                         original["retrieved"])
    new_altitudes, new = _fitted_profile(background, latitude, longitude,
                                         improved["retrieved"])
    if not np.array_equal(altitudes, old_altitudes) or not np.array_equal(altitudes, new_altitudes):
        raise ValueError("Density altitude grids differ")

    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    ax.plot(truth / 1e5, altitudes, color="#a34c23", linewidth=2.7,
            label="IRI-2016 truth")
    ax.plot(old / 1e5, altitudes, color="#727b86", linewidth=2.2,
            linestyle="--", label="Original ionogram retrieval")
    ax.plot(new / 1e5, altitudes, color="#263f67", linewidth=2.7,
            label="Coverage-qualified retrieval")
    ax.axhline(location["truth_peak_altitude_km"], color="#888888", linewidth=1,
               linestyle=":")
    ax.set(xlabel=r"Electron density ($10^5$ cm$^{-3}$)", ylabel="Altitude (km)",
           xlim=(0, max(float(truth.max()), float(old.max()), float(new.max())) / 1e5 * 1.06),
           ylim=(150, 600), title="IRI density profiles at the sounder")
    ax.grid(alpha=.23)
    ax.legend(loc="upper right")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("background", "truth_grid", "original_summary", "improved_summary", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    plot_progress(args.background, args.truth_grid, args.original_summary,
                  args.improved_summary, args.output)


if __name__ == "__main__":
    main()
