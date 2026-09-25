from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from setuptools import Extension, setup


PHARLAP_HOME = Path(os.environ.get("PHARLAP_HOME", "/Users/chartat1/pharlap")).expanduser().resolve()
GFORTRAN_LIB = Path(os.environ.get("GFORTRAN_LIB", Path(os.popen("gfortran -print-file-name=libgfortran.dylib").read().strip()).parent)).resolve()
PYLAP_SOURCE = Path(__file__).resolve().parent / "vendor" / "pylap"

if not PHARLAP_HOME.is_dir():
    raise OSError(f"PHARLAP_HOME not found: {PHARLAP_HOME}")
if not (PYLAP_SOURCE / "modules" / "source" / "raytrace_3d.c").is_file():
    raise OSError(f"vendored PyLap source not found: {PYLAP_SOURCE}")

pharlap_lib_dir = PHARLAP_HOME / "lib" / "maca"
pharlap_include_dir = PHARLAP_HOME / "src" / "C"
common_sources = sorted(str(path) for path in (PYLAP_SOURCE / "modules" / "source" / "common").glob("*.c"))

extra_link_args = [
    f"-Wl,-rpath,{GFORTRAN_LIB}",
]

raytrace_3d = Extension(
    "pylap.raytrace_3d",
    sources=[str(PYLAP_SOURCE / "modules" / "source" / "raytrace_3d.c"), *common_sources],
    include_dirs=[np.get_include(), str(pharlap_include_dir), str(PYLAP_SOURCE / "modules" / "include")],
    library_dirs=[str(pharlap_lib_dir), str(GFORTRAN_LIB)],
    libraries=["propagation", "maths", "iri2020", "gfortran", "gomp", "quadmath"],
    extra_link_args=extra_link_args,
)


setup(
    name="pylap",
    version="0.1.0+macos.local",
    description="Minimal local macOS build of pylap raytrace_3d against PHaRLAP",
    license_files=["vendor/pylap/LICENSE"],
    packages=["pylap"],
    ext_modules=[raytrace_3d],
)
