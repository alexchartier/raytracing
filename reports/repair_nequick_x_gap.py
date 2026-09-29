"""Add independently homed X returns to the first NeQuick truth ionogram.

The 2.5–2.9 MHz X rays were traced at their exact plotted frequencies after
20 kHz branch continuation. They satisfy the original 1 km homing gate.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

RUN = ROOT / "reports/data/nequick_independent_case"
ORIGINAL = RUN / "truth_ionogram.npz"
CONTINUED = RUN / "continued_x_gap_2p5_to_2p9.npz"
CORRECTED = RUN / "truth_ionogram_x_gap_recovered.npz"


def main() -> None:
    with np.load(ORIGINAL, allow_pickle=False) as data:
        payload = {key: data[key] for key in data.files}
    with np.load(CONTINUED, allow_pickle=False) as data:
        recovered = {key: data[key] for key in data.files}
    frequencies = np.asarray(payload["frequencies_mhz"], dtype=float)
    subset = np.asarray(recovered["frequencies_mhz"], dtype=float)
    if not np.allclose(subset, frequencies[5:10], atol=1e-10, rtol=0):
        raise ValueError("Recovered frequencies do not match the 100 kHz bins")
    original = np.asarray(payload["records"], dtype=float)
    if np.any((original[:, 0] >= 5) & (original[:, 0] < 10) &
              (original[:, 1] == -1)):
        raise ValueError("The source X gap is not empty; review before merging")
    added = np.asarray(recovered["records"], dtype=float).copy()
    if (len(added) != 10 or np.any(added[:, 1] != -1)
            or np.any(added[:, 0] < 0) or np.any(added[:, 0] >= 5)
            or np.any(~np.isfinite(added))
            or np.any(added[:, 2] < 150.0)
            or np.any(added[:, 3] > 1000.0)):
        raise ValueError("Recovered rays fail frequency, mode, range, or homing checks")
    per_bin = np.bincount(added[:, 0].astype(int), minlength=5)
    if not np.array_equal(per_bin, [2, 2, 2, 2, 2]):
        raise ValueError("Expected both accepted X returns at every gap frequency")
    added[:, 0] += 5
    payload["records"] = np.asarray(sorted(np.concatenate((original, added)),
                                          key=lambda row: (row[0], row[1], row[2])))
    payload["count_array"] = np.asarray(payload["count_array"], dtype=int).copy()
    payload["count_array"][5:10, 1] += per_bin
    payload["method"] = np.array("adaptive_with_dense_gap_recovery")
    provenance = {
        "source": ORIGINAL.name,
        "continued_return_file": CONTINUED.name,
        "added_x_returns": int(len(added)),
        "frequency_mhz": subset.tolist(),
        "mode": -1,
        "maximum_accepted_miss_m": float(np.max(added[:, 3])),
        "homing_tolerance_m": 1000.0,
        "continuation_step_mhz": 0.02,
        "seed_frequencies_mhz": [2.4, 2.6, 2.8, 3.0],
        "truth_density_used_only_by_forward_ray_tracer": True,
    }
    payload["x_gap_recovery_json"] = np.array(json.dumps(provenance))
    np.savez_compressed(CORRECTED, **payload)
    (RUN / "x_gap_recovery.json").write_text(json.dumps(provenance, indent=2) + "\n")
    plot(original, payload["records"], frequencies)
    print(json.dumps(provenance, indent=2))


def plot(original: np.ndarray, corrected: np.ndarray, frequencies: np.ndarray) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True,
                             height_ratios=(1.1, 0.9))
    for column, (records, title) in enumerate((
            (original, "Original adaptive truth"),
            (corrected, "Truth with recovered X returns"))):
        for row in (0, 1):
            ax = axes[row, column]
            for mode, color, label in ((1, "#1261ac", "O"),
                                       (-1, "#c52525", "X")):
                rays = records[records[:, 1] == mode]
                ax.scatter(frequencies[rays[:, 0].astype(int)], rays[:, 2],
                           color=color, s=16 if row == 0 else 25, alpha=0.8,
                           label=f"{label} ({len(rays)})")
            if row == 0:
                ax.set(xlim=(2.0, 6.0), ylim=(950, 2200), title=title,
                       ylabel="Two-way group range (km)")
                ax.legend(frameon=False)
            else:
                ax.set(xlim=(2.2, 3.2), ylim=(1050, 1190),
                       xlabel="Frequency (MHz)",
                       ylabel="Two-way group range (km)",
                       title="2.5–2.9 MHz detail")
            ax.grid(alpha=0.2)
    fig.savefig(RUN / "nequick_truth_x_gap_recovered.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
