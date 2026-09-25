"""Plot accepted returns and convergence for a completed D-ionogram fit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, BoundaryNorm

from fit_d_ionogram import Ionogram, score


def _image(ionogram: Ionogram, mode: int) -> np.ndarray:
    image = np.zeros((1450, 81), dtype=np.uint8)
    records = ionogram.records[ionogram.records[:, 1] == mode]
    for row in records:
        index = int(row[0])
        bin_km = int(np.floor(row[2])) - 150
        if 0 <= index < 81 and 0 <= bin_km < 1450:
            image[bin_km, index] = min(int(image[bin_km, index]) + 1, 3)
    return image


def plot(state_path: Path, output_prefix: Path) -> None:
    state = json.loads(state_path.read_text())
    evaluations = state["evaluations"]
    best = min(evaluations, key=lambda row: row["score"]["total"])
    observed = Ionogram.read(Path(state["observed"]))
    retrieved = Ionogram.read(Path(best["path"]))
    colors = ListedColormap(["white", "#253a5e", "#a25422", "#5e2319"])
    norm = BoundaryNorm([-.5, .5, 1.5, 2.5, 3.5], colors.N)
    fig, axes = plt.subplots(2, 2, figsize=(16, 15), sharex=True, sharey=True,
                             constrained_layout=True)
    for column, (ionogram, title) in enumerate(((observed, "Synthetic truth"), (retrieved, "Retrieved"))):
        for row, (mode, label) in enumerate(((1, "O mode"), (-1, "X mode"))):
            ax = axes[row, column]
            ax.imshow(_image(ionogram, mode), origin="lower", aspect="auto",
                      interpolation="nearest", cmap=colors, norm=norm,
                      extent=(1.95, 10.05, 150, 1600))
            ax.set_title(f"{title}: {label} ({len(ionogram.records[ionogram.records[:, 1] == mode])} returns)")
            ax.set_xlim(2, 10)
            ax.set_ylim(150, 1600)
            ax.set_xlabel("Frequency (MHz)")
            ax.set_ylabel("Group range (km)")
    fig.suptitle("D ionogram: 0.1 MHz × 1 km bins; blue = 1, ochre = 2, red = 3+ accepted returns", fontsize=15)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    comparison = output_prefix.with_name(output_prefix.name + "_ionograms.png")
    fig.savefig(comparison, dpi=300)
    plt.close(fig)

    costs = np.array([item["score"]["total"] for item in evaluations])
    best_so_far = np.minimum.accumulate(costs)
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    ax.plot(np.arange(1, len(costs)+1), best_so_far, color="#253a5e", linewidth=2)
    for boundary in (40, 80, 120):
        if boundary < len(costs):
            ax.axvline(boundary + .5, color="#888888", linestyle="--", linewidth=1)
    ax.set_yscale("log")
    ax.set(xlabel="Completed forward ionograms", ylabel="Best return score (log scale)",
           title=(f"Fit convergence; final density {best['density_scale']:.4f}, "
                  f"height shift {best['hmf2_shift_km']:+.1f} km"))
    ax.grid(alpha=.2)
    convergence = output_prefix.with_name(output_prefix.name + "_convergence.png")
    fig.savefig(convergence, dpi=200)
    plt.close(fig)
    print(json.dumps({"ionograms": str(comparison), "convergence": str(convergence),
                      "best": best, "score_check": score(observed, retrieved)}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("state", type=Path)
    parser.add_argument("output_prefix", type=Path)
    args = parser.parse_args()
    plot(args.state, args.output_prefix)


if __name__ == "__main__":
    main()
