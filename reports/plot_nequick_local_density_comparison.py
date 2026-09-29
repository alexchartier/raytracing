"""Compare the April NeQuick retrieval with and without a local density datum."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ionogram_metrics import Ionogram
from nequick_heldout_fit import layer_from_row

ROOT = Path(__file__).resolve().parents[1] / "reports/data"
BASELINE = ROOT / "nequick_prospective_case"
LOCAL = ROOT / "nequick_local_density_case"


def main() -> None:
    base_evaluation = json.loads((BASELINE / "evaluation.json").read_text())
    local_evaluation = json.loads((LOCAL / "evaluation.json").read_text())
    base_name = base_evaluation["selected"]
    local_name = local_evaluation["selected"]
    base_fit = json.loads((BASELINE / "fits.json").read_text())["candidates"][base_name]
    local_fit = json.loads((LOCAL / "fits.json").read_text())["candidates"][local_name]
    with np.load(LOCAL / "truth_profile_1km.npz", allow_pickle=False) as data:
        altitude = np.asarray(data["altitudes_km"], dtype=float)
        truth = np.asarray(data["electron_density_cm3"], dtype=float)
    baseline = layer_from_row(base_fit).density_cm3(altitude)
    retrieved = layer_from_row(local_fit).density_cm3(altitude)
    observation = Ionogram.read(LOCAL / "truth_ionogram_recovered.npz")
    old_ionogram = Ionogram.read(BASELINE / f"full_ray_{base_name}_recovered.npz")
    new_ionogram = Ionogram.read(LOCAL / f"full_ray_{local_name}.npz")
    ionograms = (observation, old_ionogram, new_ionogram)
    titles = ("Independent truth", "Ionogram only", "Ionogram + local density")
    fig = plt.figure(figsize=(13.8, 8.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 6, height_ratios=(1.1, 0.9))
    maximum_range = max(2150.0, *(float(np.max(item.records[:, 2])) + 60.0
                                  for item in ionograms))
    for column, (ionogram, title) in enumerate(zip(ionograms, titles)):
        ax = fig.add_subplot(gs[0, 2 * column:2 * column + 2])
        for mode, color, label in ((1, "#1261ac", "O"), (-1, "#c52525", "X")):
            rows = ionogram.records[ionogram.records[:, 1] == mode]
            ax.scatter(ionogram.frequencies[rows[:, 0].astype(int)], rows[:, 2],
                       color=color, s=14, alpha=0.85, label=label)
        ax.set(xlim=(2, 10), ylim=(150, maximum_range),
               xlabel="Frequency (MHz)", title=title)
        if column == 0:
            ax.set_ylabel("Two-way group range (km)")
        ax.grid(alpha=0.2)
        ax.legend(frameon=False, ncol=2, loc="upper left")
    ax = fig.add_subplot(gs[1, :3])
    ax.plot(truth / 1e5, altitude, color="black", linewidth=2.5, label="NeQuick-G truth")
    ax.plot(baseline / 1e5, altitude, color="#d68513", linewidth=2.0,
            label="Ionogram only")
    ax.plot(retrieved / 1e5, altitude, color="#168450", linewidth=2.0,
            label="Ionogram + local density")
    ax.set(xlabel="Electron density (100,000 cm$^{-3}$)",
           ylabel="Altitude (km)", ylim=(150, 800), title="Vertical density cut")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    ax = fig.add_subplot(gs[1, 3:])
    valid = (altitude >= 250) & (truth > 0)
    ax.plot(100 * (baseline[valid] / truth[valid] - 1), altitude[valid],
            color="#d68513", linewidth=2.0, label="Ionogram only")
    ax.plot(100 * (retrieved[valid] / truth[valid] - 1), altitude[valid],
            color="#168450", linewidth=2.0, label="Ionogram + local density")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set(xlabel="Density error relative to truth (%)", ylim=(250, 800),
           title="Topside density error")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    fig.suptitle("April NeQuick-G: assumed exact density at the 800 km spacecraft",
                 fontsize=14)
    output = LOCAL / "nequick_local_density_comparison.png"
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
