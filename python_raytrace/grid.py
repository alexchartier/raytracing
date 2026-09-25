from __future__ import annotations

import datetime as dt
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import netCDF4
import numpy as np
import PyIRI
from PyIRI import igrf_library, main_library

from .absorption import build_msis_atmosphere, effective_collision_frequency
from .geometry import GeoPoint, coerce_longitude_for_grid, enu_to_ecef, make_regional_lat_lon_grids, wrap_longitudes
from .indices import resolve_space_weather_indices
from .iri2020_model import Iri2020Bridge, blend_low_altitude_density, fill_temperature_profile, r12_from_f107


@dataclass
class IonosphereGrid:
    latitudes_deg: np.ndarray
    longitudes_deg: np.ndarray
    altitudes_km: np.ndarray
    iono_en_grid: np.ndarray
    iono_en_grid_5: np.ndarray
    collision_freq: np.ndarray
    iono_grid_parms: list[float]
    Bx: np.ndarray
    By: np.ndarray
    Bz: np.ndarray
    geomag_grid_parms: list[float]
    electron_temp_k: np.ndarray | None = None
    ion_temp_k: np.ndarray | None = None
    neutral_temp_k: np.ndarray | None = None
    neutral_species_cm3: np.ndarray | None = None
    metadata: dict[str, float | int | str] = field(default_factory=dict)


def save_ionosphere_grid(path: str | Path, grid: IonosphereGrid) -> Path:
    target = Path(path).expanduser()
    target.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, object] = {
        "format_version": np.array([1], dtype=np.int64),
        "latitudes_deg": np.asarray(grid.latitudes_deg, dtype=float),
        "longitudes_deg": np.asarray(grid.longitudes_deg, dtype=float),
        "altitudes_km": np.asarray(grid.altitudes_km, dtype=float),
        "iono_en_grid": np.asarray(grid.iono_en_grid, dtype=float),
        "iono_en_grid_5": np.asarray(grid.iono_en_grid_5, dtype=float),
        "collision_freq": np.asarray(grid.collision_freq, dtype=float),
        "iono_grid_parms": np.asarray(grid.iono_grid_parms, dtype=float),
        "Bx": np.asarray(grid.Bx, dtype=float),
        "By": np.asarray(grid.By, dtype=float),
        "Bz": np.asarray(grid.Bz, dtype=float),
        "geomag_grid_parms": np.asarray(grid.geomag_grid_parms, dtype=float),
        "metadata_json": np.array(json.dumps(grid.metadata, sort_keys=True)),
    }

    def _store_optional(name: str, value: np.ndarray | None) -> None:
        payload[f"{name}_present"] = np.array([value is not None], dtype=bool)
        if value is not None:
            payload[name] = np.asarray(value, dtype=float)

    _store_optional("electron_temp_k", grid.electron_temp_k)
    _store_optional("ion_temp_k", grid.ion_temp_k)
    _store_optional("neutral_temp_k", grid.neutral_temp_k)
    _store_optional("neutral_species_cm3", grid.neutral_species_cm3)
    np.savez_compressed(target, **payload)
    return target


def load_ionosphere_grid(path: str | Path) -> IonosphereGrid:
    source = Path(path).expanduser()
    with np.load(source, allow_pickle=False) as data:
        version = int(np.asarray(data["format_version"]).ravel()[0])
        if version != 1:
            raise ValueError(f"unsupported ionosphere-grid cache version: {version}")

        def _load_optional(name: str) -> np.ndarray | None:
            present = bool(np.asarray(data[f"{name}_present"]).ravel()[0])
            if not present:
                return None
            return np.asarray(data[name], dtype=float)

        metadata_text = str(np.asarray(data["metadata_json"]).item())
        metadata = json.loads(metadata_text) if metadata_text else {}
        return IonosphereGrid(
            latitudes_deg=np.asarray(data["latitudes_deg"], dtype=float),
            longitudes_deg=np.asarray(data["longitudes_deg"], dtype=float),
            altitudes_km=np.asarray(data["altitudes_km"], dtype=float),
            iono_en_grid=np.asarray(data["iono_en_grid"], dtype=float),
            iono_en_grid_5=np.asarray(data["iono_en_grid_5"], dtype=float),
            collision_freq=np.asarray(data["collision_freq"], dtype=float),
            iono_grid_parms=np.asarray(data["iono_grid_parms"], dtype=float).tolist(),
            Bx=np.asarray(data["Bx"], dtype=float),
            By=np.asarray(data["By"], dtype=float),
            Bz=np.asarray(data["Bz"], dtype=float),
            geomag_grid_parms=np.asarray(data["geomag_grid_parms"], dtype=float).tolist(),
            electron_temp_k=_load_optional("electron_temp_k"),
            ion_temp_k=_load_optional("ion_temp_k"),
            neutral_temp_k=_load_optional("neutral_temp_k"),
            neutral_species_cm3=_load_optional("neutral_species_cm3"),
            metadata=metadata,
        )


def save_ionosphere_grid_netcdf(path: str | Path, grid: IonosphereGrid) -> Path:
    target = Path(path).expanduser()
    target.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(target, "w") as dataset:
        dataset.setncattr("format_version", 1)
        dataset.setncattr("metadata_json", json.dumps(grid.metadata, sort_keys=True))
        dataset.createDimension("lat", len(grid.latitudes_deg))
        dataset.createDimension("lon", len(grid.longitudes_deg))
        dataset.createDimension("alt", len(grid.altitudes_km))
        dataset.createDimension("parm", len(grid.iono_grid_parms))
        dataset.createDimension("geomag_parm", len(grid.geomag_grid_parms))
        if grid.neutral_species_cm3 is not None and np.asarray(grid.neutral_species_cm3).ndim == 4:
            neutral_species = np.asarray(grid.neutral_species_cm3, dtype=float)
            if neutral_species.shape[:3] == (len(grid.latitudes_deg), len(grid.longitudes_deg), len(grid.altitudes_km)):
                dataset.createDimension("species", int(neutral_species.shape[3]))
            elif neutral_species.shape[1:] == (len(grid.latitudes_deg), len(grid.longitudes_deg), len(grid.altitudes_km)):
                dataset.createDimension("species", int(neutral_species.shape[0]))
            else:
                raise ValueError(
                    "neutral_species_cm3 must have shape (lat, lon, alt, species) "
                    "or (species, lat, lon, alt)"
                )

        dataset.createVariable("latitudes_deg", "f8", ("lat",))[:] = np.asarray(grid.latitudes_deg, dtype=float)
        dataset.createVariable("longitudes_deg", "f8", ("lon",))[:] = np.asarray(grid.longitudes_deg, dtype=float)
        dataset.createVariable("altitudes_km", "f8", ("alt",))[:] = np.asarray(grid.altitudes_km, dtype=float)
        dataset.createVariable("iono_grid_parms", "f8", ("parm",))[:] = np.asarray(grid.iono_grid_parms, dtype=float)
        dataset.createVariable("geomag_grid_parms", "f8", ("geomag_parm",))[:] = np.asarray(grid.geomag_grid_parms, dtype=float)

        for name, value in (
            ("iono_en_grid", grid.iono_en_grid),
            ("iono_en_grid_5", grid.iono_en_grid_5),
            ("collision_freq", grid.collision_freq),
            ("Bx", grid.Bx),
            ("By", grid.By),
            ("Bz", grid.Bz),
        ):
            dataset.createVariable(name, "f8", ("lat", "lon", "alt"), zlib=True, complevel=3)[:] = np.asarray(value, dtype=float)

        for name, value in (
            ("electron_temp_k", grid.electron_temp_k),
            ("ion_temp_k", grid.ion_temp_k),
            ("neutral_temp_k", grid.neutral_temp_k),
            ("neutral_species_cm3", grid.neutral_species_cm3),
        ):
            if value is None:
                continue
            array = np.asarray(value, dtype=float)
            if name == "neutral_species_cm3" and array.ndim == 4:
                if array.shape[:3] == (len(grid.latitudes_deg), len(grid.longitudes_deg), len(grid.altitudes_km)):
                    dims = ("lat", "lon", "alt", "species")
                elif array.shape[1:] == (len(grid.latitudes_deg), len(grid.longitudes_deg), len(grid.altitudes_km)):
                    dims = ("species", "lat", "lon", "alt")
                else:
                    raise ValueError(
                        "neutral_species_cm3 must have shape (lat, lon, alt, species) "
                        "or (species, lat, lon, alt)"
                    )
            else:
                dims = ("lat", "lon", "alt")
            dataset.createVariable(name, "f8", dims, zlib=True, complevel=3)[:] = array
    return target


def load_ionosphere_grid_netcdf(path: str | Path) -> IonosphereGrid:
    source = Path(path).expanduser()
    with netCDF4.Dataset(source) as dataset:
        version = int(dataset.getncattr("format_version"))
        if version != 1:
            raise ValueError(f"unsupported ionosphere-grid netCDF cache version: {version}")

        def _load_optional(name: str) -> np.ndarray | None:
            if name not in dataset.variables:
                return None
            array = np.asarray(dataset.variables[name][:], dtype=float)
            if name == "neutral_species_cm3":
                if array.ndim == 4:
                    return array
                # Backward compatibility: some older caches wrote a malformed 3-D
                # neutral-species field. Treat it as unavailable and reuse the
                # cached collision-frequency field instead of forcing a rebuild.
                return None
            return array

        metadata_json = str(dataset.getncattr("metadata_json")) if "metadata_json" in dataset.ncattrs() else ""
        metadata = json.loads(metadata_json) if metadata_json else {}
        return IonosphereGrid(
            latitudes_deg=np.asarray(dataset.variables["latitudes_deg"][:], dtype=float),
            longitudes_deg=np.asarray(dataset.variables["longitudes_deg"][:], dtype=float),
            altitudes_km=np.asarray(dataset.variables["altitudes_km"][:], dtype=float),
            iono_en_grid=np.asarray(dataset.variables["iono_en_grid"][:], dtype=float),
            iono_en_grid_5=np.asarray(dataset.variables["iono_en_grid_5"][:], dtype=float),
            collision_freq=np.asarray(dataset.variables["collision_freq"][:], dtype=float),
            iono_grid_parms=np.asarray(dataset.variables["iono_grid_parms"][:], dtype=float).tolist(),
            Bx=np.asarray(dataset.variables["Bx"][:], dtype=float),
            By=np.asarray(dataset.variables["By"][:], dtype=float),
            Bz=np.asarray(dataset.variables["Bz"][:], dtype=float),
            geomag_grid_parms=np.asarray(dataset.variables["geomag_grid_parms"][:], dtype=float).tolist(),
            electron_temp_k=_load_optional("electron_temp_k"),
            ion_temp_k=_load_optional("ion_temp_k"),
            neutral_temp_k=_load_optional("neutral_temp_k"),
            neutral_species_cm3=_load_optional("neutral_species_cm3"),
            metadata=metadata,
        )


def _normalize_utc(when: dt.datetime) -> dt.datetime:
    if when.tzinfo is None:
        return when
    return when.astimezone(dt.timezone.utc).replace(tzinfo=None)


def _grid_parms(latitudes_deg: np.ndarray, longitudes_deg: np.ndarray, altitudes_km: np.ndarray) -> list[float]:
    lat_step = float(latitudes_deg[1] - latitudes_deg[0]) if len(latitudes_deg) > 1 else 0.0
    lon_step = float(longitudes_deg[1] - longitudes_deg[0]) if len(longitudes_deg) > 1 else 0.0
    alt_step = float(altitudes_km[1] - altitudes_km[0]) if len(altitudes_km) > 1 else 0.0
    return [
        float(latitudes_deg[0]), lat_step, float(len(latitudes_deg)),
        float(longitudes_deg[0]), lon_step, float(len(longitudes_deg)),
        float(altitudes_km[0]), alt_step, float(len(altitudes_km)),
    ]


def _build_global_sample_axes(
    lat_step_deg: float,
    lon_step_deg: float,
    regional_longitudes_deg: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    latitudes_deg = np.arange(-90.0, 90.0 + 0.5 * lat_step_deg, lat_step_deg, dtype=float)
    latitudes_deg = np.clip(latitudes_deg, -90.0, 90.0)
    latitudes_deg = np.unique(latitudes_deg)

    lon_start = math.floor(float(np.min(regional_longitudes_deg)) / lon_step_deg) * lon_step_deg
    lon_stop = lon_start + 360.0
    longitudes_deg = np.arange(lon_start, lon_stop + 0.5 * lon_step_deg, lon_step_deg, dtype=float)
    return latitudes_deg, longitudes_deg


def _find_axis_indices(source_axis: np.ndarray, requested_axis: np.ndarray, *, atol: float = 1e-6) -> np.ndarray:
    indices = []
    for value in np.asarray(requested_axis, dtype=float):
        matches = np.flatnonzero(np.isclose(source_axis, value, atol=atol, rtol=0.0))
        if matches.size == 0:
            raise ValueError(f"requested axis value {value} is not present in the source grid")
        indices.append(int(matches[0]))
    return np.asarray(indices, dtype=int)


def _extract_regional_subgrid(
    global_grid: np.ndarray,
    global_latitudes_deg: np.ndarray,
    global_longitudes_deg: np.ndarray,
    regional_latitudes_deg: np.ndarray,
    regional_longitudes_deg: np.ndarray,
) -> np.ndarray:
    lat_indices = _find_axis_indices(global_latitudes_deg, regional_latitudes_deg)
    lon_indices = _find_axis_indices(global_longitudes_deg, regional_longitudes_deg)
    return np.take(np.take(global_grid, lat_indices, axis=0), lon_indices, axis=1)


def extract_ionosphere_subgrid(
    grid: IonosphereGrid,
    latitudes_deg: np.ndarray,
    longitudes_deg: np.ndarray,
    altitudes_km: np.ndarray | None = None,
) -> IonosphereGrid:
    latitudes_deg = np.asarray(latitudes_deg, dtype=float)
    longitudes_deg = np.asarray(longitudes_deg, dtype=float)
    selected_altitudes_km = np.asarray(grid.altitudes_km if altitudes_km is None else altitudes_km, dtype=float)

    lat_indices = _find_axis_indices(np.asarray(grid.latitudes_deg, dtype=float), latitudes_deg)
    mapped_lons = np.asarray(
        [coerce_longitude_for_grid(float(value), np.asarray(grid.longitudes_deg, dtype=float)) for value in longitudes_deg],
        dtype=float,
    )
    lon_indices = _find_axis_indices(np.asarray(grid.longitudes_deg, dtype=float), mapped_lons)
    alt_indices = _find_axis_indices(np.asarray(grid.altitudes_km, dtype=float), selected_altitudes_km)

    def _subset(array: np.ndarray | None) -> np.ndarray | None:
        if array is None:
            return None
        values = np.asarray(array, dtype=float)
        if values.ndim == 4 and values.shape[1:] == (len(grid.latitudes_deg), len(grid.longitudes_deg), len(grid.altitudes_km)):
            return np.take(np.take(np.take(values, lat_indices, axis=1), lon_indices, axis=2), alt_indices, axis=3)
        subset = np.take(np.take(np.take(values, lat_indices, axis=0), lon_indices, axis=1), alt_indices, axis=2)
        return subset

    return IonosphereGrid(
        latitudes_deg=latitudes_deg,
        longitudes_deg=longitudes_deg,
        altitudes_km=selected_altitudes_km,
        iono_en_grid=_subset(grid.iono_en_grid),
        iono_en_grid_5=_subset(grid.iono_en_grid_5),
        collision_freq=_subset(grid.collision_freq),
        iono_grid_parms=_grid_parms(latitudes_deg, longitudes_deg, selected_altitudes_km),
        Bx=_subset(grid.Bx),
        By=_subset(grid.By),
        Bz=_subset(grid.Bz),
        geomag_grid_parms=_grid_parms(latitudes_deg, longitudes_deg, selected_altitudes_km),
        electron_temp_k=_subset(grid.electron_temp_k),
        ion_temp_k=_subset(grid.ion_temp_k),
        neutral_temp_k=_subset(grid.neutral_temp_k),
        neutral_species_cm3=_subset(grid.neutral_species_cm3),
        metadata=dict(grid.metadata),
    )


def build_pyiri_grid_from_axes(
    when: dt.datetime,
    latitudes_deg: np.ndarray,
    longitudes_deg: np.ndarray,
    altitudes_km: np.ndarray,
    *,
    f107: float | None = None,
    d_region_model: str = "fpt2018",
    ap_daily: float | None = None,
    blend_bottom_km: float = 120.0,
    blend_top_km: float = 140.0,
    msis_version: float | str = 2.1,
    refresh_indices: bool = False,
) -> IonosphereGrid:
    when = _normalize_utc(when)
    indices = resolve_space_weather_indices(
        when,
        f107=f107,
        ap_daily=ap_daily,
        refresh=refresh_indices,
    )
    latitudes_deg = np.asarray(latitudes_deg, dtype=float)
    longitudes_deg = np.asarray(longitudes_deg, dtype=float)
    altitudes_km = np.asarray(altitudes_km, dtype=float)
    if latitudes_deg.ndim != 1 or longitudes_deg.ndim != 1 or altitudes_km.ndim != 1:
        raise ValueError("latitudes_deg, longitudes_deg, and altitudes_km must each be one-dimensional arrays")
    if latitudes_deg.size < 2 or longitudes_deg.size < 2 or altitudes_km.size < 2:
        raise ValueError("latitudes_deg, longitudes_deg, and altitudes_km must each contain at least two samples")

    if len(latitudes_deg) > 701 or len(longitudes_deg) > 701 or len(altitudes_km) > 401:
        raise ValueError("grid exceeds pylap raytrace_3d limits (701 lat, 701 lon, 401 alt)")

    lat_step_deg = float(latitudes_deg[1] - latitudes_deg[0]) if len(latitudes_deg) > 1 else 1.0
    lon_step_deg = float(longitudes_deg[1] - longitudes_deg[0]) if len(longitudes_deg) > 1 else 1.0
    global_latitudes_deg, global_longitudes_deg = _build_global_sample_axes(
        lat_step_deg,
        lon_step_deg,
        longitudes_deg,
    )
    global_wrapped_longitudes_deg = wrap_longitudes(global_longitudes_deg)
    lon_mesh, lat_mesh = np.meshgrid(global_wrapped_longitudes_deg, global_latitudes_deg, indexing="xy")
    alon = lon_mesh.ravel()
    alat = lat_mesh.ravel()
    ut_hours = (
        when.hour
        + when.minute / 60.0
        + when.second / 3600.0
        + when.microsecond / 3.6e9
    )
    _, _, _, _, _, _, edens_m3 = main_library.IRI_density_1day(
        when.year,
        when.month,
        when.day,
        np.array([ut_hours], dtype=float),
        alon,
        alat,
        altitudes_km,
        indices.f107,
        PyIRI.coeff_dir,
    )

    density_alt_global_lat_lon_m3 = edens_m3[0].reshape(len(altitudes_km), len(global_latitudes_deg), len(global_longitudes_deg))
    density_global_lat_lon_alt_m3 = np.transpose(density_alt_global_lat_lon_m3, (1, 2, 0))
    density_lat_lon_alt_m3 = _extract_regional_subgrid(
        density_global_lat_lon_alt_m3,
        global_latitudes_deg,
        global_longitudes_deg,
        latitudes_deg,
        longitudes_deg,
    )

    neutral_temp_k = None
    neutral_species_cm3 = None
    electron_temp_k = None
    ion_temp_k = None
    wrapped_lons = wrap_longitudes(longitudes_deg)

    if d_region_model.lower() != "none":
        msis_atmosphere = build_msis_atmosphere(
            when,
            latitudes_deg,
            longitudes_deg,
            altitudes_km,
            f107=indices.f107,
            f107a=indices.f107a,
            ap_vector=indices.ap_vector,
            ap_daily=indices.ap_daily,
            version=msis_version,
        )
        neutral_temp_k = msis_atmosphere.neutral_temperature_k
        neutral_species_cm3 = msis_atmosphere.neutral_species_cm3
        electron_temp_k = np.empty_like(density_lat_lon_alt_m3)
        ion_temp_k = np.empty_like(density_lat_lon_alt_m3)
        bridge = Iri2020Bridge()
        r12 = r12_from_f107(indices.f107)

        merged_density_m3 = density_lat_lon_alt_m3.copy()
        for lat_idx, lat_deg in enumerate(latitudes_deg):
            for lon_idx, lon_deg in enumerate(wrapped_lons):
                iri_profile = bridge.profile(
                    when,
                    float(lat_deg),
                    float(lon_deg),
                    r12=r12,
                    altitudes_km=altitudes_km,
                    d_region_model=d_region_model,
                )
                fallback_temp = neutral_temp_k[lat_idx, lon_idx]
                electron_temp_k[lat_idx, lon_idx] = fill_temperature_profile(
                    altitudes_km,
                    iri_profile.electron_temperature_k,
                    fallback_temp,
                )
                ion_temp_k[lat_idx, lon_idx] = fill_temperature_profile(
                    altitudes_km,
                    iri_profile.ion_temperature_k,
                    fallback_temp,
                )
                merged_density_m3[lat_idx, lon_idx] = blend_low_altitude_density(
                    altitudes_km,
                    density_lat_lon_alt_m3[lat_idx, lon_idx],
                    iri_profile.electron_density_m3,
                    blend_bottom_km=blend_bottom_km,
                    blend_top_km=blend_top_km,
                )

        density_lat_lon_alt_m3 = merged_density_m3
        collision_freq = effective_collision_frequency(
            electron_temp_k,
            ion_temp_k,
            density_lat_lon_alt_m3,
            msis_atmosphere.neutral_species_cm3,
        )
    else:
        density_lat_lon_alt_m3 = density_lat_lon_alt_m3.copy()
        collision_freq = np.zeros_like(density_lat_lon_alt_m3)

    density_lat_lon_alt_cm3 = density_lat_lon_alt_m3 / 1e6

    lat3d, lon3d, alt3d = np.meshgrid(
        latitudes_deg,
        wrapped_lons,
        altitudes_km,
        indexing="ij",
    )
    date_decimal = main_library.decimal_year(when)
    _, north_nt, east_nt, down_nt, _, _, _ = igrf_library.inclination(
        PyIRI.coeff_dir,
        date_decimal,
        lon3d.ravel(),
        lat3d.ravel(),
        alt=alt3d.ravel(),
        only_inc=False,
    )

    b_ecef = enu_to_ecef(
        east=np.asarray(east_nt, dtype=float) * 1e-9,
        north=np.asarray(north_nt, dtype=float) * 1e-9,
        up=-np.asarray(down_nt, dtype=float) * 1e-9,
        lat_deg=lat3d.ravel(),
        lon_deg=lon3d.ravel(),
    )
    bx = b_ecef[:, 0].reshape(lat3d.shape)
    by = b_ecef[:, 1].reshape(lat3d.shape)
    bz = b_ecef[:, 2].reshape(lat3d.shape)

    grid_parms = _grid_parms(latitudes_deg, longitudes_deg, altitudes_km)

    return IonosphereGrid(
        latitudes_deg=latitudes_deg,
        longitudes_deg=longitudes_deg,
        altitudes_km=altitudes_km,
        iono_en_grid=density_lat_lon_alt_cm3,
        iono_en_grid_5=density_lat_lon_alt_cm3.copy(),
        collision_freq=collision_freq,
        iono_grid_parms=grid_parms,
        Bx=bx,
        By=by,
        Bz=bz,
        geomag_grid_parms=grid_parms.copy(),
        electron_temp_k=electron_temp_k,
        ion_temp_k=ion_temp_k,
        neutral_temp_k=neutral_temp_k,
        neutral_species_cm3=neutral_species_cm3,
        metadata={
            "source_grid": "global_pyiri_then_regional_subset",
            "source_lat_count": len(global_latitudes_deg),
            "source_lon_count": len(global_longitudes_deg),
            "d_region_model": d_region_model,
            "f107": indices.f107,
            "f107a": indices.f107a,
            "ap_daily": indices.ap_daily,
            "indices_source": indices.source,
            "refresh_indices": int(refresh_indices),
            "blend_bottom_km": float(blend_bottom_km),
            "blend_top_km": float(blend_top_km),
            "msis_version": str(msis_version),
        },
    )


def build_pyiri_grid(when: dt.datetime, tx: GeoPoint, rx: GeoPoint, *,
                     f107: float | None = None,
                     alt_min_km: float = 60.0,
                     alt_max_km: float = 500.0,
                     alt_step_km: float = 2.0,
                     lat_step_deg: float = 1.0,
                     lon_step_deg: float = 1.0,
                     lat_margin_deg: float = 5.0,
                     lon_margin_deg: float = 5.0,
                     d_region_model: str = "fpt2018",
                     ap_daily: float | None = None,
                     blend_bottom_km: float = 120.0,
                     blend_top_km: float = 140.0,
                     msis_version: float | str = 2.1,
                     refresh_indices: bool = False) -> IonosphereGrid:
    latitudes_deg, longitudes_deg = make_regional_lat_lon_grids(
        tx,
        rx,
        lat_step_deg=lat_step_deg,
        lon_step_deg=lon_step_deg,
        lat_margin_deg=lat_margin_deg,
        lon_margin_deg=lon_margin_deg,
    )
    altitudes_km = np.arange(alt_min_km, alt_max_km + 0.5 * alt_step_km, alt_step_km, dtype=float)
    return build_pyiri_grid_from_axes(
        when,
        latitudes_deg,
        longitudes_deg,
        altitudes_km,
        f107=f107,
        d_region_model=d_region_model,
        ap_daily=ap_daily,
        blend_bottom_km=blend_bottom_km,
        blend_top_km=blend_top_km,
        msis_version=msis_version,
        refresh_indices=refresh_indices,
    )
