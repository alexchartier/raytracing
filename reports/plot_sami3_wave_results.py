"""Plot paired O/X ionograms and latitude–altitude density sections for a wave pilot."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))
from python_raytrace.grid import load_ionosphere_grid_netcdf
from retrieve_sami3_wave import _peak_density_and_height

CASE = ROOT / "reports/data/sami3_wave_20170111_0600"
FIGURES = ROOT / "reports/figures"


def _points(indices: list[int], oblique: bool) -> np.ndarray:
    manifest = json.loads((CASE / "manifest.json").read_text())
    lat = np.array([manifest["profiles"][i - 1]["latitude_deg"] for i in indices], dtype=float)
    if oblique:
        os.environ["SOUNDER_CASE_ROOT"] = str(CASE)
        from oblique_wave_pass import endpoints
        lat = np.array([(endpoints(i)[0].lat_deg + endpoints(i)[1].lat_deg) / 2
                        for i in indices])
    return np.column_stack((lat, np.full(len(lat), 190.0)))


def plot(truth_ions: Path, model_ions: Path, truth_grid: Path,
         model_grid: Path, prefix: str, indices: list[int], oblique: bool) -> dict:
    FIGURES.mkdir(parents=True, exist_ok=True)
    colors = {1: "#1769aa", -1: "#c73632"}
    selected = [i for i in (9, 14, 20) if i in indices]
    fig, axes = plt.subplots(len(selected), 2, figsize=(11, 3.1 * len(selected)),
                             sharex=True, constrained_layout=True, squeeze=False)
    for row, index in enumerate(selected):
        for col, (directory, label) in enumerate(((truth_ions, "SAMI3 truth"),
                                                   (model_ions, "Retrieved"))):
            with np.load(directory / f"ionogram_{index:02d}.npz", allow_pickle=False) as source:
                records = np.asarray(source["records"])
                frequencies = np.asarray(source["frequencies_mhz"])
            ax = axes[row, col]
            for mode, name in ((1, "O"), (-1, "X")):
                subset = records[records[:, 1] == mode]
                ax.scatter(frequencies[subset[:, 0].astype(int)], np.rint(subset[:, 2]),
                           s=11, color=colors[mode], alpha=.75, linewidths=0, label=name)
            ax.set(xlim=(2, 10), ylim=(700 if oblique else 150, 2200 if oblique else 1850),
                   title=f"{'Link' if oblique else 'Profile'} {index:02d} · {label}")
            ax.grid(alpha=.12)
            if col == 0:
                ax.set_ylabel("Group range (km)")
            if row == len(selected) - 1:
                ax.set_xlabel("Frequency (MHz)")
            if row == 0 and col == 0:
                ax.legend(loc="upper right", title="Mode")
    ionogram_path = FIGURES / f"{prefix}_ionograms.png"
    fig.savefig(ionogram_path, dpi=220)
    plt.close(fig)

    a = load_ionosphere_grid_netcdf(truth_grid)
    b = load_ionosphere_grid_netcdf(model_grid)
    if not (np.array_equal(a.latitudes_deg, b.latitudes_deg)
            and np.array_equal(a.longitudes_deg, b.longitudes_deg)
            and np.array_equal(a.altitudes_km, b.altitudes_km)):
        raise ValueError("Truth and retrieved density grids must share axes")
    points = _points(indices, oblique)
    truth = RegularGridInterpolator((a.latitudes_deg, a.longitudes_deg),
                                    a.iono_en_grid)(points)
    model = RegularGridInterpolator((b.latitudes_deg, b.longitudes_deg),
                                    b.iono_en_grid)(points)
    positions = points[:, 0]
    alt = a.altitudes_km
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    vmax = float(np.quantile(truth, .98) / 1e5)
    for ax, density, title in zip(axes[:2], (truth, model),
                                  ("SAMI3 truth", "Retrieved")):
        image = ax.contourf(positions, alt, density.T / 1e5, levels=np.linspace(0, vmax, 22),
                            extend="max", cmap="viridis")
        ax.set(title=title, xlabel="Latitude (°)", ylabel="Altitude (km)", ylim=(150, 600))
    fig.colorbar(image, ax=axes[:2], label="Electron density (10⁵ cm⁻³)")
    residual = (model - truth).T / 1e5
    bound = float(np.quantile(abs(residual), .98))
    image2 = axes[2].contourf(positions, alt, residual, levels=np.linspace(-bound, bound, 23),
                               extend="both", cmap="RdBu_r")
    axes[2].set(title="Retrieved − truth", xlabel="Latitude (°)",
                ylabel="Altitude (km)", ylim=(150, 600))
    fig.colorbar(image2, ax=axes[2], label="Density difference (10⁵ cm⁻³)")
    density_path = FIGURES / f"{prefix}_density.png"
    fig.savefig(density_path, dpi=220)
    plt.close(fig)
    peak_truth, height_truth = _peak_density_and_height(alt, truth)
    peak_model, height_model = _peak_density_and_height(alt, model)
    fig = plt.figure(figsize=(12, 5), constrained_layout=True)
    layout = fig.add_gridspec(2, 2, width_ratios=(.9, 1.1))
    profile_ax = fig.add_subplot(layout[:, 0])
    frequency_ax = fig.add_subplot(layout[0, 1])
    height_ax = fig.add_subplot(layout[1, 1], sharex=frequency_ax)
    cut_index = 14 if 14 in indices else indices[len(indices) // 2]
    cut_row = indices.index(cut_index)
    profile_ax.plot(truth[cut_row] / 1e5, alt, color="#222222", label="SAMI3 truth")
    profile_ax.plot(model[cut_row] / 1e5, alt, color="#1769aa", linestyle="--",
                    label="Retrieved")
    profile_ax.set(xlabel="Electron density (10⁵ cm⁻³)", ylabel="Altitude (km)",
                   ylim=(150, 600),
                   title=f"{'Link midpoint' if oblique else 'Sounder'} {cut_index:02d} altitude cut")
    profile_ax.legend()
    frequency_ax.plot(positions, .00898 * np.sqrt(peak_truth), "o-", color="#222222",
                      label="SAMI3 truth")
    frequency_ax.plot(positions, .00898 * np.sqrt(peak_model), "s--", color="#1769aa",
                      label="Retrieved")
    frequency_ax.set(ylabel="foF2 (MHz)", title="Local F2 peak along the pass")
    height_ax.plot(positions, height_truth, "o-", color="#222222")
    height_ax.plot(positions, height_model, "s--", color="#1769aa")
    height_ax.set(xlabel="Latitude (°)", ylabel="hmF2 (km)")
    for ax in (profile_ax, frequency_ax, height_ax):
        ax.grid(alpha=.2)
    peak_path = FIGURES / f"{prefix}_peaks.png"
    fig.savefig(peak_path, dpi=220)
    plt.close(fig)
    return {"ionograms": str(ionogram_path), "density": str(density_path),
            "peaks": str(peak_path)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-ions", type=Path, required=True)
    parser.add_argument("--model-ions", type=Path, required=True)
    parser.add_argument("--truth-grid", type=Path, required=True)
    parser.add_argument("--model-grid", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--indices", nargs="+", type=int, default=[1, 5, 9, 12, 14, 17, 20])
    parser.add_argument("--oblique", action="store_true")
    args = parser.parse_args()
    print(json.dumps(plot(args.truth_ions, args.model_ions, args.truth_grid,
                          args.model_grid, args.prefix, args.indices, args.oblique), indent=2))


if __name__ == "__main__":
    main()
