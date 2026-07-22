from __future__ import annotations

import ctypes
import datetime as dt
import os
import platform
from dataclasses import dataclass
from pathlib import Path

import numpy as np


_OUTF_ROWS = 15
_OUTF_DTYPE = np.float32
_OARR_DTYPE = np.float32
_ERROR_BUFFER_LEN = 512


class Iri2020BridgeError(RuntimeError):
    pass


@dataclass
class Iri2020Profile:
    electron_density_m3: np.ndarray
    neutral_temperature_k: np.ndarray
    ion_temperature_k: np.ndarray
    electron_temperature_k: np.ndarray
    extras: np.ndarray


def r12_from_f107(f107: float) -> float:
    if f107 <= 63.75:
        return 1.0
    a = 0.00089
    b = 0.728
    c = 63.75 - float(f107)
    root = (-b + np.sqrt(b * b - 4.0 * a * c)) / (2.0 * a)
    return float(np.clip(root, 1.0, 200.0))


def fill_temperature_profile(
    altitudes_km: np.ndarray,
    profile_k: np.ndarray,
    fallback_k: np.ndarray,
) -> np.ndarray:
    alts = np.asarray(altitudes_km, dtype=float)
    profile = np.asarray(profile_k, dtype=float)
    fallback = np.asarray(fallback_k, dtype=float)
    valid = np.isfinite(profile) & (profile > 0.0)
    if np.any(valid):
        filled = np.interp(alts, alts[valid], profile[valid], left=profile[valid][0], right=profile[valid][-1])
    else:
        filled = fallback.copy()
    filled[~np.isfinite(filled) | (filled <= 0.0)] = fallback[~np.isfinite(filled) | (filled <= 0.0)]
    return filled


def repair_d_region_density_profile(altitudes_km: np.ndarray, density_m3: np.ndarray) -> np.ndarray:
    alts = np.asarray(altitudes_km, dtype=float)
    density = np.asarray(density_m3, dtype=float).copy()
    valid = np.isfinite(density) & (density > 0.0)
    if not np.any(valid):
        return np.zeros_like(density)

    first_valid = np.flatnonzero(valid)[:5]
    if first_valid.size > 0:
        anchor_alts = np.concatenate(([40.0], alts[first_valid]))
        anchor_log_density = np.log(np.concatenate(([1e-10], density[first_valid])))
        invalid = ~valid
        if np.any(invalid):
            density[invalid] = np.exp(np.interp(alts[invalid], anchor_alts, anchor_log_density))

    density[alts < 40.0] = 0.0
    density[~np.isfinite(density)] = 0.0
    return density


def blend_low_altitude_density(
    altitudes_km: np.ndarray,
    upper_density_m3: np.ndarray,
    lower_density_m3: np.ndarray,
    *,
    blend_bottom_km: float = 120.0,
    blend_top_km: float = 140.0,
) -> np.ndarray:
    alts = np.asarray(altitudes_km, dtype=float)
    upper = np.maximum(np.asarray(upper_density_m3, dtype=float), 1e-10)
    lower = np.maximum(repair_d_region_density_profile(alts, lower_density_m3), 1e-10)
    merged = upper.copy()

    if blend_top_km <= blend_bottom_km:
        merged[alts <= blend_bottom_km] = lower[alts <= blend_bottom_km]
        return merged

    low_mask = alts <= blend_bottom_km
    mid_mask = (alts > blend_bottom_km) & (alts < blend_top_km)
    merged[low_mask] = lower[low_mask]

    if np.any(mid_mask):
        low_anchor = np.exp(np.interp(blend_bottom_km, alts, np.log(lower)))
        high_anchor = np.exp(np.interp(blend_top_km, alts, np.log(upper)))
        blend_axis = np.array([blend_bottom_km, blend_top_km], dtype=float)
        blend_values = np.log(np.array([low_anchor, high_anchor], dtype=float))
        merged[mid_mask] = np.exp(np.interp(alts[mid_mask], blend_axis, blend_values))

    return merged


def _default_library_path() -> Path:
    suffix = ".dylib" if platform.system() == "Darwin" else ".so"
    return Path(__file__).resolve().parent / "_lib" / f"libiri2020_bridge{suffix}"


def _default_reference_data_dir() -> Path | None:
    pharlap_home = Path(os.environ.get("PHARLAP_HOME", "/Users/chartat1/pharlap")).expanduser()
    candidate = pharlap_home / "dat"
    return candidate if candidate.is_dir() else None


def _ensure_reference_data_env() -> None:
    if os.environ.get("DIR_MODELS_REF_DAT"):
        return
    default_dir = _default_reference_data_dir()
    if default_dir is not None:
        os.environ["DIR_MODELS_REF_DAT"] = str(default_dir)


class Iri2020Bridge:
    def __init__(self, library_path: str | os.PathLike[str] | None = None):
        _ensure_reference_data_env()
        lib_path = Path(library_path) if library_path is not None else _default_library_path()
        if not lib_path.is_file():
            raise Iri2020BridgeError(
                f"IRI2020 bridge library not found at {lib_path}. "
                "Build it with `python3 python_raytrace/iri2020_build/build.py`."
            )
        self._lib = ctypes.CDLL(str(lib_path))
        self._func = self._lib.python_raytrace_iri2020_profile
        self._func.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_char),
            ctypes.c_int,
        ]
        self._func.restype = ctypes.c_int

    def profile(
        self,
        when: dt.datetime,
        lat_deg: float,
        lon_deg: float,
        *,
        r12: float,
        altitudes_km: np.ndarray,
        d_region_model: str = "fpt2018",
    ) -> Iri2020Profile:
        alts = np.asarray(altitudes_km, dtype=float)
        if alts.ndim != 1 or alts.size < 2:
            raise ValueError("altitudes_km must be a one-dimensional array with at least two elements")
        alt_step = np.diff(alts)
        if not np.allclose(alt_step, alt_step[0], atol=1e-6, rtol=0.0):
            raise ValueError("altitudes_km must be uniformly spaced for the IRI2020 bridge")

        d_model = d_region_model.lower()
        use_fpt = 1 if d_model in {"fpt2018", "firi", "firi2020"} else 0

        outf = np.empty((_OUTF_ROWS * alts.size,), dtype=_OUTF_DTYPE)
        oarr = np.empty((100,), dtype=_OARR_DTYPE)
        error_message = ctypes.create_string_buffer(_ERROR_BUFFER_LEN)
        status = self._func(
            float(lat_deg),
            float(lon_deg),
            float(r12),
            int(when.year),
            int(when.month),
            int(when.day),
            int(when.hour),
            int(when.minute),
            float(alts[0]),
            float(alt_step[0]),
            int(alts.size),
            use_fpt,
            outf.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            oarr.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            error_message,
            _ERROR_BUFFER_LEN,
        )
        if status != 0:
            raise Iri2020BridgeError(error_message.value.decode("utf-8", errors="replace") or "unknown IRI2020 bridge error")

        profiles = outf.reshape((alts.size, _OUTF_ROWS)).astype(float, copy=False)
        density_m3 = profiles[:, 0].copy()
        if use_fpt:
            density_m3 = repair_d_region_density_profile(alts, density_m3)

        return Iri2020Profile(
            electron_density_m3=density_m3,
            neutral_temperature_k=profiles[:, 1].copy(),
            ion_temperature_k=profiles[:, 2].copy(),
            electron_temperature_k=profiles[:, 3].copy(),
            extras=oarr.astype(float, copy=True),
        )
