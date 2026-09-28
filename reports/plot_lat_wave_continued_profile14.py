"""High-resolution O/X truth and retrieved ionograms for wave profile 14."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TRUTH = ROOT / "reports/data/lat_wave_continued_truth_ionograms/ionogram_14.npz"
RETRIEVED = ROOT / "reports/data/lat_wave_continued_retrieved_ionograms/ionogram_14.npz"
FIGURES = ROOT / "reports/figures"


def draw(ax, path: Path, zoom: bool) -> tuple[float, float]:
    with np.load(path, allow_pickle=False) as source:
        records = np.asarray(source["records"], dtype=float)
        frequencies = np.asarray(source["frequencies_mhz"], dtype=float)
        if int(source["profile_index"]) != 14:
            raise ValueError(f"Expected wave profile 14 in {path}")
    noses = []
    for mode, color, label in ((1, "#135b9a", "O"), (-1, "#b63836", "X")):
        selected = records[records[:, 1] == mode]
        noses.append(float(frequencies[int(np.max(selected[:, 0]))]))
        ax.scatter(frequencies[selected[:, 0].astype(int)],
                   np.round(selected[:, 2]),
                   s=22 if zoom else 11, marker="s", linewidths=0,
                   color=color, alpha=0.85, label=label)
    ax.set(xlim=(4.2, 6.4) if zoom else (2, 10),
           ylim=(1100, 2100) if zoom else (150, 2100),
           ylabel="Group range (km)")
    ax.grid(alpha=0.16)
    if zoom:
        ax.set_xticks(np.arange(4.2, 6.41, 0.2))
    return tuple(noses)


def main(retrieved_path: Path = RETRIEVED,
         output_stem: str = "lat_wave_continued_profile14_ionograms") -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 9),
                             constrained_layout=True, sharey="row")
    for col, (name, path) in enumerate((("Truth: IRI-2016 + wave", TRUTH),
                                       ("Retrieved: PyIRI + fitted wave", retrieved_path))):
        noses = draw(axes[0, col], path, False)
        draw(axes[1, col], path, True)
        axes[0, col].set_title(f"{name}\nO nose {noses[0]:.1f} MHz · X nose {noses[1]:.1f} MHz")
        axes[1, col].set_xlabel("Frequency (MHz)")
    axes[0, 1].legend(loc="upper right")
    fig.suptitle("Latitude-wave profile 14 · all accepted O/X returns · 0.1 MHz × 1 km")
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{output_stem}.png", dpi=300)
    fig.savefig(FIGURES / f"{output_stem}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retrieved", type=Path, default=RETRIEVED)
    parser.add_argument("--output-stem", default="lat_wave_continued_profile14_ionograms")
    args = parser.parse_args()
    main(args.retrieved, args.output_stem)
