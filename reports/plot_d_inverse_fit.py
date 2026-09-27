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
    best = (next(row for row in evaluations if row["path"] == state["selected_path"])
            if "selected_path" in state else
            min(evaluations, key=lambda row: row["score"]["total"]))
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
    for boundary in range(40, len(costs), 40):
        if boundary < len(costs):
            ax.axvline(boundary + .5, color="#888888", linestyle="--", linewidth=1)
    ax.set_yscale("log")
    if "selected_path" in state:
        ax.axhline(best["score"]["total"], color="#b14a2d", linestyle=":",
                   linewidth=1.5, label="Selected coverage-qualified fit")
        ax.legend()
    width_label = (f", width {best['f2_width_scale']:.3f}"
                   if "f2_width_scale" in best else "")
    top_label = (f", top ratio {best['topside_width_ratio']:.3f}"
                 if "topside_width_ratio" in best else "")
    title = (f"Fit convergence; selected density {best['density_scale']:.4f}, "
             f"height shift {best['hmf2_shift_km']:+.1f} km{width_label}{top_label}")
    if top_label:
        title = (f"Fit convergence; selected scale {best['density_scale']:.4f}, "
                 f"shift {best['hmf2_shift_km']:+.1f} km\n"
                 f"F2 width {best['f2_width_scale']:.3f}, "
                 f"topside ratio {best['topside_width_ratio']:.3f}")
    ax.set(xlabel="Completed forward ionograms", ylabel="Best return score (log scale)",
           title=title)
    ax.grid(alpha=.2)
    convergence = output_prefix.with_name(output_prefix.name + "_convergence.png")
    fig.savefig(convergence, dpi=200)
    plt.close(fig)
    if state.get("search_model") == "scale_height_width_topside":
        from fit_d_ionogram_topside import count_aware_score
        score_check = count_aware_score(observed, retrieved)
    else:
        score_check = score(observed, retrieved)
    print(json.dumps({"ionograms": str(comparison), "convergence": str(convergence),
                      "best": best, "score_check": score_check}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("state", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--truth-label", default="Synthetic truth")
    args = parser.parse_args()
    plot(args.state, args.output_prefix, truth_label=args.truth_label)


if __name__ == "__main__":
    main()
