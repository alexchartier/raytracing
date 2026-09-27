"""Plot the latitude-altitude density section through the wave pass."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "reports/data/lat_wave_pass_manifest.json"
PNG = ROOT / "reports/figures/lat_wave_pass_density_contour.png"
PDF = ROOT / "reports/figures/lat_wave_pass_density_contour.pdf"


def plot() -> None:
    manifest = json.loads(MANIFEST.read_text())
    with np.load(ROOT / manifest["density_grid"], allow_pickle=False) as data:
        latitudes = np.asarray(data["latitudes_deg"], dtype=float)
        longitudes = np.asarray(data["longitudes_deg"], dtype=float)
        altitudes = np.asarray(data["altitudes_km"], dtype=float)
        density = np.asarray(data["electron_density_cm3"], dtype=float)
    longitude = float(manifest["longitude_deg"])
    j = int(np.searchsorted(longitudes, longitude) - 1)
    fraction = (longitude - longitudes[j]) / (longitudes[j + 1] - longitudes[j])
    section = (1 - fraction) * density[:, j, :] + fraction * density[:, j + 1, :]

    pass_latitudes = np.array([row["latitude_deg"] for row in manifest["profiles"]])
    display_latitudes = np.linspace(pass_latitudes[0], pass_latitudes[-1], 401)
    display_altitudes = np.arange(150.0, 600.1, 5.0)
    lat_mesh, alt_mesh = np.meshgrid(display_latitudes, display_altitudes)
    display = RegularGridInterpolator((latitudes, altitudes), section)(
        np.column_stack((lat_mesh.ravel(), alt_mesh.ravel()))
    ).reshape(lat_mesh.shape) / 1e5

    fig, ax = plt.subplots(figsize=(11, 6))
    levels = np.arange(0.0, 5.01, 0.2)
    filled = ax.contourf(lat_mesh, alt_mesh, display, levels=levels,
                         cmap="viridis", extend="max")
    lines = ax.contour(lat_mesh, alt_mesh, display, levels=(1.0, 2.0, 3.0, 4.0),
                       colors="white", linewidths=0.7, alpha=0.75)
    ax.clabel(lines, inline=True, fmt="%.0f", fontsize=8)
    ax.vlines(pass_latitudes, 150, 155, color="white", linewidth=1.0, alpha=0.9)
    ax.set(xlim=(pass_latitudes[0], pass_latitudes[-1]), ylim=(150, 600),
           xlabel="Latitude (degrees)", ylabel="Altitude (km)",
           title="Electron density across the 20° latitude wave pass")
    colorbar = fig.colorbar(filled, ax=ax, pad=0.02)
    colorbar.set_label("Electron density (100,000 cm$^{-3}$)")
    fig.text(0.12, 0.02,
             f"IRI-2016 + imposed wave at {longitude:.2f}° longitude, 12:00 UTC. "
             "White ticks mark sounders; contours interpolate a 2° × 20 km density grid.",
             fontsize=9, color="#444444")
    fig.subplots_adjust(bottom=0.16, right=0.92)
    PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG, dpi=250)
    fig.savefig(PDF)
    plt.close(fig)
    print(PNG)
    print(PDF)


if __name__ == "__main__":
    plot()
