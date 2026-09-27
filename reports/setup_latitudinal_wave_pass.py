"""Create a 20-position vertical sounder pass through an IRI density wave.

The IRI-2016 density is held at the 2010-01-01 12:00 UTC model epoch so
changes along the pass are spatial. The imposed wave is north-south, with a
900 km wavelength that the existing 2-degree density grid resolves.
"""

from __future__ import annotations

import datetime as dt
import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from python_raytrace.absorption import effective_collision_frequency  # noqa: E402
from python_raytrace.grid import build_pyiri_grid_from_axes, save_ionosphere_grid_netcdf  # noqa: E402


TRUTH = ROOT / "reports/data/d_inverse_iri_truth_density.npz"
REFERENCE = ROOT / "reports/data/d_inverse_iri_pyiri_background.npz"
GRID_OUTPUT = ROOT / "reports/data/lat_wave_pass_truth_density.npz"
FORWARD_GRID_OUTPUT = ROOT / "reports/data/lat_wave_pass_forward_grid.nc"
MANIFEST_OUTPUT = ROOT / "reports/data/lat_wave_pass_manifest.json"
FIGURE_OUTPUT = ROOT / "reports/figures/lat_wave_pass_density.png"
PROFILE_COUNT = 20
LATITUDE_SPAN_DEG = 20.0
SOUNDER_ALT_KM = 800.0
WAVE_AMPLITUDE = 0.18
WAVE_WAVELENGTH_KM = 900.0
WAVE_CENTER_ALT_KM = 280.0
WAVE_VERTICAL_SIGMA_KM = 60.0


def build() -> None:
    with np.load(TRUTH, allow_pickle=False) as source:
        latitudes = np.asarray(source["latitudes_deg"], dtype=float)
        longitudes = np.asarray(source["longitudes_deg"], dtype=float)
        altitudes = np.asarray(source["altitudes_km"], dtype=float)
        baseline = np.asarray(source["electron_density_cm3"], dtype=float)
        model_epoch = str(source["time_utc"])
    with np.load(REFERENCE, allow_pickle=False) as reference:
        center_lat = float(reference["tx_lat_deg"])
        longitude = float(reference["tx_lon_deg"])
    if not (latitudes[0] < center_lat - 10 and center_lat + 10 < latitudes[-1]
            and longitudes[0] < longitude < longitudes[-1]):
        raise ValueError("The 20-degree pass does not fit within the IRI density grid")
    if np.any(baseline <= 0) or not np.all(np.isfinite(baseline)):
        raise ValueError("IRI baseline contains invalid densities")

    north_km = 6371.0088 * np.deg2rad(latitudes - center_lat)
    horizontal = np.cos(2 * np.pi * north_km / WAVE_WAVELENGTH_KM)
    vertical = np.exp(-0.5 * ((altitudes - WAVE_CENTER_ALT_KM)
                              / WAVE_VERTICAL_SIGMA_KM) ** 2)
    factor = 1.0 + WAVE_AMPLITUDE * horizontal[:, None, None] * vertical[None, None, :]
    density = baseline * factor
    positions = np.linspace(center_lat - LATITUDE_SPAN_DEG / 2,
                            center_lat + LATITUDE_SPAN_DEG / 2, PROFILE_COUNT)
    center_time = dt.datetime.fromisoformat(model_epoch)
    times = [center_time + dt.timedelta(seconds=16 * (index - (PROFILE_COUNT - 1) / 2))
             for index in range(PROFILE_COUNT)]
    manifest = {
        "description": "Frozen-epoch vertical sounder latitude pass through independent IRI-2016 density with an imposed wave",
        "density_grid": str(GRID_OUTPUT.relative_to(ROOT)),
        "forward_grid": str(FORWARD_GRID_OUTPUT.relative_to(ROOT)),
        "model_epoch_utc": model_epoch,
        "latitude_center_deg": center_lat,
        "longitude_deg": longitude,
        "latitude_span_deg": LATITUDE_SPAN_DEG,
        "profile_count": PROFILE_COUNT,
        "sounder_altitude_km": SOUNDER_ALT_KM,
        "wave": {
            "amplitude_fraction_at_280_km": WAVE_AMPLITUDE,
            "wavelength_km": WAVE_WAVELENGTH_KM,
            "bearing_deg": 0.0,
            "phase_rad_at_center": 0.0,
            "vertical_center_km": WAVE_CENTER_ALT_KM,
            "vertical_sigma_km": WAVE_VERTICAL_SIGMA_KM,
        },
        "frequency_sweep_mhz": [2.0, 10.0, 0.1],
        "homing_tolerance_m": 1000.0,
        "vertical_fan_layout": "equal_area_guarded_D",
        "profiles": [
            {"index": index + 1, "latitude_deg": float(latitude),
             "longitude_deg": longitude, "altitude_km": SOUNDER_ALT_KM,
             "time_utc": when.isoformat()}
            for index, (latitude, when) in enumerate(zip(positions, times))
        ],
    }
    GRID_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    FIGURE_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(GRID_OUTPUT, latitudes_deg=latitudes, longitudes_deg=longitudes,
                        altitudes_km=altitudes, electron_density_cm3=density,
                        model=np.array("IRI-2016 with imposed north-south wave"),
                        time_utc=np.array(model_epoch))
    forward_base = build_pyiri_grid_from_axes(
        center_time, latitudes, longitudes, altitudes,
        f107=120.0, ap_daily=8.0, d_region_model="none",
    )
    collision = forward_base.collision_freq
    if (forward_base.electron_temp_k is not None
            and forward_base.ion_temp_k is not None
            and forward_base.neutral_species_cm3 is not None):
        collision = effective_collision_frequency(
            forward_base.electron_temp_k, forward_base.ion_temp_k,
            density * 1e6, forward_base.neutral_species_cm3,
        )
    forward_grid = replace(forward_base, iono_en_grid=density,
                           iono_en_grid_5=density, collision_freq=collision)
    save_ionosphere_grid_netcdf(FORWARD_GRID_OUTPUT, forward_grid)
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n")

    points = np.column_stack((positions, np.full(PROFILE_COUNT, longitude)))
    base_profiles = RegularGridInterpolator((latitudes, longitudes), baseline)(points)
    wave_profiles = RegularGridInterpolator((latitudes, longitudes), density)(points)
    peak_base = np.max(base_profiles, axis=1)
    peak_wave = np.max(wave_profiles, axis=1)
    difference = 100 * (wave_profiles / base_profiles - 1.0)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    axes[0].plot(positions, peak_base / 1e5, color="#888888", linewidth=2,
                 label="IRI-2016 baseline")
    axes[0].plot(positions, peak_wave / 1e5, color="#254f83", linewidth=2,
                 marker="o", markersize=3, label="IRI-2016 with wave")
    axes[0].set(xlabel="Latitude (degrees)", ylabel="Peak density (100,000 cm$^{-3}$)",
                title="20 sounder positions over 20° latitude")
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    image = axes[1].pcolormesh(positions, altitudes, difference.T, cmap="RdBu_r",
                               vmin=-WAVE_AMPLITUDE * 100, vmax=WAVE_AMPLITUDE * 100,
                               shading="nearest")
    axes[1].set(xlabel="Latitude (degrees)", ylabel="Altitude (km)",
                ylim=(150, 500), title="Imposed density change")
    fig.colorbar(image, ax=axes[1], label="Change from IRI baseline (%)")
    fig.savefig(FIGURE_OUTPUT, dpi=220)
    plt.close(fig)
    print(json.dumps({"manifest": str(MANIFEST_OUTPUT), "density_grid": str(GRID_OUTPUT),
                      "forward_grid": str(FORWARD_GRID_OUTPUT),
                      "figure": str(FIGURE_OUTPUT), "profile_count": PROFILE_COUNT,
                      "peak_density_change_pct": [float(v) for v in
                                                  (100 * (peak_wave / peak_base - 1))]}, indent=2))


if __name__ == "__main__":
    build()
