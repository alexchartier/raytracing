"""Plot paired O/X ionograms and midpoint density after selection is frozen."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from evaluate_general_field_sami import profiles  # noqa: E402

RUN = ROOT / "reports/data/oblique_model_families"
LABELS = {"nequick": "NeQuick-G", "chapman": "Generalized Chapman",
          "iri2016": "IRI-2016"}


def draw(ax, path: Path, title: str, limits: tuple[int, int]) -> None:
    with np.load(path, allow_pickle=False) as source:
        frequency = source["frequencies_mhz"]
        records = source["records"]
    low, high = limits
    extent = (float(frequency[0] - .05), float(frequency[-1] + .05),
              low - .5, high + .5)
    for mode, color in ((1, "#2166ac"), (-1, "#b2182b")):
        chosen = records[records[:, 1] == mode]
        image = np.zeros((high - low + 1, len(frequency), 4), dtype=np.float32)
        ranges = np.rint(chosen[:, 2]).astype(int)
        bins = chosen[:, 0].astype(int)
        valid = (ranges >= low) & (ranges <= high)
        image[ranges[valid] - low, bins[valid]] = to_rgba(color, alpha=1.0)
        ax.imshow(image, origin="lower", extent=extent, interpolation="nearest",
                  aspect="auto")
    ax.set(xlim=(2, 10), ylim=limits, xlabel="Frequency (MHz)",
           ylabel="Group range (km)", title=title)
    ax.grid(alpha=.18)


def main() -> None:
    os.umask(0o077)
    for case, label in LABELS.items():
        root = RUN / case
        choice = json.loads((root / "selection.json").read_text())["chosen_candidate"]
        paths = [root / "forward" / name / "ionogram_01.npz"
                 for name in ("truth", choice)]
        ranges = []
        with np.load(paths[0], allow_pickle=False) as ray:
            midpoint = (float(ray["tx_lat_deg"]) + float(ray["rx_lat_deg"])) / 2
            longitude = float(ray["tx_lon_deg"])
        for path in paths:
            with np.load(path, allow_pickle=False) as source:
                ranges.extend(source["records"][:, 2].tolist())
        low = max(700, int(np.floor((min(ranges) - 30) / 50) * 50))
        high = int(np.ceil((max(ranges) + 30) / 50) * 50)
        figure = plt.figure(figsize=(11.6, 9.0), constrained_layout=True)
        panel = figure.add_gridspec(2, 2, height_ratios=(1.45, .8))
        axes = [figure.add_subplot(panel[0, 0]), figure.add_subplot(panel[0, 1])]
        draw(axes[0], paths[0], f"{label} truth", (low, high))
        draw(axes[1], paths[1], "Retrieved", (low, high))
        axes[0].legend(handles=[Patch(color="#2166ac", label="O"),
                                Patch(color="#b2182b", label="X")])
        locations = np.array([[midpoint, longitude]])
        altitude, truth = profiles(root / "truth_grid.nc", locations)
        _, start = profiles(root / "start/grid.nc", locations)
        _, selected = profiles(root / choice / "grid.nc", locations)
        ax = figure.add_subplot(panel[1, :])
        ax.plot(truth[0] / 1e5, altitude, color="black", linewidth=2.2,
                label=f"{label} truth")
        ax.plot(start[0] / 1e5, altitude, color="gray", linestyle="--",
                linewidth=1.8, label="Vertical/in-situ start")
        ax.plot(selected[0] / 1e5, altitude, color="#2166ac", linewidth=2,
                label="Oblique-selected")
        ax.set(xlabel="Electron density ($10^5$ cm$^{-3}$)",
               ylabel="Altitude (km)", ylim=(150, 800),
               title="Density at the oblique link midpoint")
        ax.grid(alpha=.2)
        ax.legend()
        target = root / "comparison.png"
        figure.savefig(target, dpi=300)
        target.chmod(0o600)
        plt.close(figure)
        print(target)


if __name__ == "__main__":
    main()
