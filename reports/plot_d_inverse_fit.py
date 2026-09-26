"""Plot accepted returns and convergence for a completed D-ionogram fit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fit_d_ionogram import Ionogram, score
from plot_ionogram_pair import plot_pair


def plot(state_path: Path, output_prefix: Path, *, truth_label: str = "Synthetic truth") -> None:
    state = json.loads(state_path.read_text())
    evaluations = state["evaluations"]
    best = min(evaluations, key=lambda row: row["score"]["total"])
    observed = Ionogram.read(Path(state["observed"]))
    retrieved = Ionogram.read(Path(best["path"]))
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    comparison = output_prefix.with_name(output_prefix.name + "_ionograms.png")
    plot_pair(Path(state["observed"]), Path(best["path"]), comparison,
              truth_label=truth_label, retrieved_label="Retrieved PyIRI")

    costs = np.array([item["score"]["total"] for item in evaluations])
    best_so_far = np.minimum.accumulate(costs)
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    ax.plot(np.arange(1, len(costs)+1), best_so_far, color="#253a5e", linewidth=2)
    for boundary in (40, 80, 120):
        if boundary < len(costs):
            ax.axvline(boundary + .5, color="#888888", linestyle="--", linewidth=1)
    ax.set_yscale("log")
    width_label = (f", width {best['f2_width_scale']:.3f}"
                   if "f2_width_scale" in best else "")
    ax.set(xlabel="Completed forward ionograms", ylabel="Best return score (log scale)",
           title=(f"Fit convergence; final density {best['density_scale']:.4f}, "
                  f"height shift {best['hmf2_shift_km']:+.1f} km{width_label}"))
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
    parser.add_argument("--truth-label", default="Synthetic truth")
    args = parser.parse_args()
    plot(args.state, args.output_prefix, truth_label=args.truth_label)


if __name__ == "__main__":
    main()
