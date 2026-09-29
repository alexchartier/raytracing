"""Prepare a frozen SAMI3 wave cut and an independent 2017 PyIRI prior."""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from python_raytrace.grid import build_pyiri_grid_from_axes, save_ionosphere_grid_netcdf

SOURCE = Path("/Users/chartat1/data/sami3/2017_tid/sami_mat/2017-01-11_0600.mat")
RUN = ROOT / "reports/data/sami3_wave_20170111_0600"
FIGURE = ROOT / "reports/figures/sami3_wave_20170111_0600_truth.png"
EPOCH = dt.datetime(2017, 1, 11, 6, tzinfo=dt.timezone.utc)


def _vertical_resample(source_alt: np.ndarray, values: np.ndarray,
                       target_alt: np.ndarray) -> np.ndarray:
    """Interpolate within the snapshot and exponentially extend its end slopes."""
    flat = values.reshape(-1, len(source_alt))
    out = np.empty((len(flat), len(target_alt)))
    interior = (target_alt >= source_alt[0]) & (target_alt <= source_alt[-1])
    below = target_alt < source_alt[0]
    above = target_alt > source_alt[-1]
    for index, profile in enumerate(flat):
        out[index, interior] = np.interp(target_alt[interior], source_alt, profile)
        low_rate = np.clip((np.log(profile[5]) - np.log(profile[0])) /
                           (source_alt[5] - source_alt[0]), -0.2, 0.2)
        high_rate = np.clip((np.log(profile[-1]) - np.log(profile[-6])) /
                            (source_alt[-1] - source_alt[-6]), -0.2, 0.0)
        out[index, below] = profile[0] * np.exp(low_rate * (target_alt[below] - source_alt[0]))
        out[index, above] = profile[-1] * np.exp(high_rate * (target_alt[above] - source_alt[-1]))
    return out.reshape(*values.shape[:-1], len(target_alt))


def build() -> None:
    RUN.mkdir(parents=True, exist_ok=True)
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    source = loadmat(SOURCE, variable_names=["lat", "lon", "alt", "dene", "iono_en_grid", "time"])
    matlab_time = float(source["time"].item())
    source_epoch = (dt.datetime.fromordinal(int(matlab_time))
                    + dt.timedelta(days=matlab_time % 1)
                    - dt.timedelta(days=366)).replace(tzinfo=dt.timezone.utc)
    if abs((source_epoch - EPOCH).total_seconds()) > 1:
        raise ValueError(f"Source time {source_epoch} differs from requested epoch {EPOCH}")
    digest = hashlib.sha256()
    with SOURCE.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    src_lat = source["lat"].ravel().astype(float)
    src_lon = source["lon"].ravel().astype(float)
    src_alt = source["alt"].ravel().astype(float)
    src_density = np.transpose(source["dene"], (1, 2, 0))
    # The stored raytracing density is dene reordered from 0..360 to -180..180.
    assert np.array_equal(src_density[:, np.flatnonzero(src_lon == 190)[0], :],
                          source["iono_en_grid"][:, np.flatnonzero(np.arange(-180, 181, 5) == -170)[0], :])
    lat = np.arange(-84.0, -31.9, 2.0)
    lon = np.arange(150.0, 230.1, 5.0)
    alt = np.arange(80.0, 820.1, 5.0)
    subset = src_density[np.ix_(np.isin(src_lat, lat), np.isin(src_lon, lon),
                                np.arange(len(src_alt), dtype=int))]
    if subset.shape != (len(lat), len(lon), len(src_alt)) or not np.all(np.isfinite(subset)) or np.any(subset <= 0):
        raise ValueError("SAMI3 regional density has missing or nonpositive values")
    truth = _vertical_resample(src_alt, subset, alt)
    pyiri = build_pyiri_grid_from_axes(EPOCH, lat, lon, alt, f107=75.0,
                                      ap_daily=5.0, d_region_model="none")
    if not np.all(np.isfinite(pyiri.iono_en_grid)):
        raise ValueError("PyIRI prior contains invalid densities")
    truth_grid = replace(pyiri, iono_en_grid=truth, iono_en_grid_5=truth.copy(),
                         metadata={**pyiri.metadata, "density_source": "SAMI3 2017-01-11 06:00 snapshot"})
    save_ionosphere_grid_netcdf(RUN / "truth_grid.nc", truth_grid)
    save_ionosphere_grid_netcdf(RUN / "pyiri_prior_grid.nc", pyiri)
    np.savez_compressed(RUN / "truth_density.npz", latitudes_deg=lat, longitudes_deg=lon,
                        altitudes_km=alt, electron_density_cm3=truth,
                        model="SAMI3 saved dene snapshot", time_utc=EPOCH.isoformat())
    positions = np.linspace(-70.0, -50.0, 20)
    manifest = {
        "description": "Frozen remapped SAMI3/HIAMCM wave-rich latitude cut; PyIRI is an independent density prior",
        "provenance": "Repository remap_sami.m identifies 2017_tid/sami_mat as its HIAMCM-version output; original source NetCDF is not present locally",
        "source_mat": str(SOURCE), "source_sha256": digest.hexdigest(),
        "epoch_utc": EPOCH.isoformat(),
        "latitude_center_deg": -60.0, "longitude_deg": 190.0,
        "latitude_span_deg": 20.0, "profile_count": 20,
        "sounder_altitude_km": 800.0, "satellite_separation_km": 600.0,
        "density_units": "cm^-3, as in source iono_en_grid raytracing field",
        "source_altitude_km": [92, 800], "working_altitude_km": [80, 820, 5],
        "pyiri_f107": 75.0, "pyiri_ap_daily": 5.0,
        "frequency_sweep_mhz": [2.0, 10.0, 0.1],
        "homing_tolerance_m": 1000.0, "vertical_fan_option": "D",
        "profiles": [{"index": i + 1, "latitude_deg": float(p),
                       "longitude_deg": 190.0, "altitude_km": 800.0,
                       "time_utc": EPOCH.isoformat()}
                      for i, p in enumerate(positions)],
    }
    (RUN / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    section = RegularGridInterpolator((lat, lon), truth)(
        np.column_stack((positions, np.full(20, 190.0))))
    peak = section.max(1)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    image = axes[0].contourf(positions, alt, section.T / 1e5, levels=24, cmap="viridis")
    axes[0].set(xlabel="Latitude (°)", ylabel="Altitude (km)", ylim=(150, 500),
                title="SAMI3 truth density at 190°E")
    fig.colorbar(image, ax=axes[0], label="Electron density (10⁵ cm⁻³)")
    axes[1].plot(positions, .00898 * np.sqrt(peak), "o-", color="#1f5a91")
    axes[1].set(xlabel="Latitude (°)", ylabel="Peak plasma frequency (MHz)",
                title="Trough and crest on the selected cut")
    axes[1].grid(alpha=.25)
    fig.savefig(FIGURE, dpi=220)
    plt.close(fig)
    print(json.dumps({"run": str(RUN), "figure": str(FIGURE),
                      "fof2_range_mhz": [.00898 * float(np.sqrt(peak.min())),
                                         .00898 * float(np.sqrt(peak.max()))],
                      "pyiri_prior_peak_range_mhz":
                      [.00898 * float(np.sqrt(pyiri.iono_en_grid.max(axis=2).min())),
                       .00898 * float(np.sqrt(pyiri.iono_en_grid.max(axis=2).max()))]}, indent=2))


if __name__ == "__main__":
    build()
