from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


PHARLAP_HOME = Path(os.environ.get("PHARLAP_HOME", "/Users/chartat1/pharlap")).expanduser().resolve()
PYLAP_SOURCE = Path(os.environ.get("PYLAP_SOURCE", "/tmp/PyLap")).expanduser().resolve()
GFORTRAN_LIB = Path(os.environ.get("GFORTRAN_LIB", Path(os.popen("gfortran -print-file-name=libgfortran.dylib").read().strip()).parent)).resolve()
PATCH_FILE = Path(__file__).resolve().parent / "patches" / "cached_state_vector.patch"

if not PHARLAP_HOME.is_dir():
    raise OSError(f"PHARLAP_HOME not found: {PHARLAP_HOME}")
if not PYLAP_SOURCE.is_dir():
    raise OSError(f"PYLAP_SOURCE not found: {PYLAP_SOURCE}")

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


class PatchedBuildExt(build_ext):
    def build_extensions(self):
        build_root = Path(self.build_temp).resolve() / "pylap_patched"
        source_dir = build_root / "modules" / "source"
        include_dir = build_root / "modules" / "include"
        source_dir.mkdir(parents=True, exist_ok=True)
        include_dir.mkdir(parents=True, exist_ok=True)
        patched_source = source_dir / "raytrace_3d.c"
        shutil.copy2(PYLAP_SOURCE / "modules" / "source" / "raytrace_3d.c", patched_source)
        shutil.copy2(PYLAP_SOURCE / "modules" / "include" / "pharlap.h", include_dir / "pharlap.h")
        subprocess.run(
            ["patch", "-t", "-p1", "-d", str(build_root), "-i", str(PATCH_FILE)],
            check=True,
        )
        raytrace_3d.sources[0] = str(patched_source)
        super().build_extensions()

setup(
    name="pylap",
    version="0.1.0+macos.local",
    description="Minimal local macOS build of pylap raytrace_3d against PHaRLAP",
    packages=["pylap"],
    ext_modules=[raytrace_3d],
    cmdclass={"build_ext": PatchedBuildExt},
)
