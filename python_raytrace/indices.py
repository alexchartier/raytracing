from __future__ import annotations

import datetime as dt
import ssl
import subprocess
import urllib.request
from dataclasses import dataclass

import certifi
import numpy as np
import pymsis


@dataclass(frozen=True)
class SpaceWeatherIndices:
    f107: float
    f107a: float
    ap_daily: float
    ap_vector: np.ndarray
    source: str


def _coerce_datetime64_minute(when: dt.datetime) -> np.datetime64:
    return np.datetime64(when.replace(tzinfo=None), "m")


def refresh_pymsis_indices() -> None:
    try:
        pymsis.utils.download_f107_ap()
        return
    except Exception:
        pass

    data_url = getattr(pymsis.utils, "_F107_AP_URL", None)
    data_path = getattr(pymsis.utils, "_F107_AP_PATH", None)
    if data_url is None or data_path is None:
        raise RuntimeError("pymsis did not expose its F10.7/Ap cache path and URL")

    context = ssl.create_default_context(cafile=certifi.where())
    try:
        with urllib.request.urlopen(data_url, context=context) as response:
            data = response.read()
    except Exception:
        result = subprocess.run(
            ["curl", "-fsSL", str(data_url)],
            check=True,
            capture_output=True,
        )
        data = result.stdout
    data_path.parent.mkdir(parents=True, exist_ok=True)
    data_path.write_bytes(data)


def resolve_space_weather_indices(
    when: dt.datetime,
    *,
    f107: float | None = None,
    f107a: float | None = None,
    ap_daily: float | None = None,
    ap_vector: np.ndarray | None = None,
    refresh: bool = False,
) -> SpaceWeatherIndices:
    if f107 is not None and f107a is None:
        f107a = float(f107)
    if ap_daily is not None and ap_vector is None:
        ap_vector = np.full((7,), float(ap_daily), dtype=float)

    fetched = None
    need_fetched_data = f107 is None or f107a is None or ap_daily is None or ap_vector is None

    if refresh:
        refresh_pymsis_indices()

    if f107 is None or f107a is None or ap_daily is None or ap_vector is None:
        try:
            fetched_f107, fetched_f107a, fetched_ap = pymsis.utils.get_f107_ap(
                np.array([_coerce_datetime64_minute(when)], dtype="datetime64[m]")
            )
        except Exception:
            refresh_pymsis_indices()
            fetched_f107, fetched_f107a, fetched_ap = pymsis.utils.get_f107_ap(
                np.array([_coerce_datetime64_minute(when)], dtype="datetime64[m]")
            )
        fetched = (
            float(fetched_f107[0]),
            float(fetched_f107a[0]),
            np.asarray(fetched_ap[0], dtype=float),
        )

    if fetched is None:
        resolved_f107 = float(f107)
        resolved_f107a = float(f107a)
        resolved_ap_vector = np.asarray(ap_vector, dtype=float)
        source = "manual"
    else:
        fetched_f107_value, fetched_f107a_value, fetched_ap_vector = fetched
        resolved_f107 = float(f107) if f107 is not None else fetched_f107_value
        resolved_f107a = float(f107a) if f107a is not None else fetched_f107a_value
        resolved_ap_vector = np.asarray(ap_vector, dtype=float) if ap_vector is not None else fetched_ap_vector
        source = "manual+pymsis" if any(value is not None for value in (f107, f107a, ap_daily, ap_vector)) else "pymsis"

    if resolved_ap_vector.shape != (7,):
        raise ValueError("ap_vector must have shape (7,)")

    resolved_ap_daily = float(ap_daily) if ap_daily is not None else float(resolved_ap_vector[0])
    resolved_ap_vector = resolved_ap_vector.copy()
    resolved_ap_vector[0] = resolved_ap_daily

    return SpaceWeatherIndices(
        f107=resolved_f107,
        f107a=resolved_f107a,
        ap_daily=resolved_ap_daily,
        ap_vector=resolved_ap_vector,
        source=source,
    )
