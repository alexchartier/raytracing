"""Serialize memory-intensive local PyLap jobs across separate commands."""

from __future__ import annotations

import fcntl
import os
import sys
from contextlib import contextmanager
from pathlib import Path

import netCDF4

LOCK = Path(os.environ.get(
    "RAYTRACING_RAY_LOCK",
    "/private/tmp/raytracing_chartat1_local_ray.lock" if sys.platform == "darwin"
    else "/tmp/raytracing_chartat1_local_ray.lock"))


def require_remote_for_large_grid(grid_path: Path) -> None:
    """Reject full high-resolution regional grids on the 16 GB Mac."""
    if sys.platform != "darwin":
        return
    with netCDF4.Dataset(grid_path) as dataset:
        shape = tuple(len(dataset.dimensions[name]) for name in ("lat", "lon", "alt"))
    if all(actual >= minimum for actual, minimum in zip(shape, (20, 10, 100))):
        raise RuntimeError(
            f"Grid {shape} is too large for local PyLap tracing; the SAMI3 "
            "case measured 26.7 GiB peak RSS. Run it on a private Cartman compute node.")


@contextmanager
def local_ray_lock():
    if sys.platform == "darwin" and os.environ.get("RAYTRACING_ALLOW_LOCAL_RAYS") != "1":
        raise RuntimeError(
            "Local PyLap tracing is disabled on this 16 GB workstation. "
            "Use private Cartman compute jobs, or set RAYTRACING_ALLOW_LOCAL_RAYS=1 "
            "only after the AGENTS.md memory and swap checks.")
    fd = os.open(LOCK, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        os.fchmod(fd, 0o600)
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)
