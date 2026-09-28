"""Show profile-14 option-A truth homing before and after nose recovery."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reports/data"
FIGURES = ROOT / "reports/figures"
SOURCES = (
    ("A adaptive", DATA / "lat_wave_profile14_option_a_truth.npz"),
    ("A + full-fan gap checks", DATA / "lat_wave_profile14_option_a_truth_recovered.npz"),
    ("A + 20 kHz continuation", DATA / "lat_wave_profile14_option_a_truth_nose_refined.npz"),
)


def draw(axis, records: np.ndarray, frequencies: np.ndarray, size: float) -> None:
    for mode, color, label in ((1, "#135b9a", "O"), (-1, "#b63836", "X")):
        selected = records[records[:, 1] == mode]
        axis.scatter(frequencies[selected[:, 0].astype(int)],
                     np.round(selected[:, 2]), s=size, marker="s",
                     linewidths=0, color=color, alpha=0.84, label=label)
    axis.grid(alpha=0.15)


def main() -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharey="row",
                             constrained_layout=True)
    for column, (name, path) in enumerate(SOURCES):
        with np.load(path, allow_pickle=False) as source:
            records = np.asarray(source["records"], dtype=float)
            frequencies = np.asarray(source["frequencies_mhz"], dtype=float)
            if int(source["profile_index"]) != 14:
                raise ValueError(f"Expected profile 14 in {path}")
        o_records = records[records[:, 1] == 1]
        o_nose = frequencies[int(np.max(o_records[:, 0]))]
        draw(axes[0, column], records, frequencies, 11)
        axes[0, column].set(title=f"{name}\nO nose {o_nose:.1f} MHz",
                            xlim=(2, 10), ylim=(150, 2100))
        draw(axes[1, column], records, frequencies, 31)
        axes[1, column].set(xlabel="Frequency (MHz)", xlim=(4.4, 5.5),
                            ylim=(1350, 2050))
        axes[1, column].set_xticks(np.arange(4.4, 5.51, 0.1))
    axes[0, 0].set_ylabel("Group range (km)")
    axes[1, 0].set_ylabel("Group range (km)")
    axes[0, -1].legend(loc="upper right", frameon=True)
    fig.suptitle("Truth profile 14 · option A · O blue / X red · 0.1 MHz × 1 km bins")
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / "lat_wave_profile14_option_a_truth.png", dpi=300)
    fig.savefig(FIGURES / "lat_wave_profile14_option_a_truth.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
