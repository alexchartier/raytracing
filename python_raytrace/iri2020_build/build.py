from __future__ import annotations

import os
import platform
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path(__file__).resolve().parent / "iri2020_bridge.c"
OUTPUT_DIR = ROOT / "python_raytrace" / "_lib"


def _gfortran_library_dir() -> Path:
    candidates = ["libgfortran.dylib", "libgfortran.so"]
    for name in candidates:
        result = subprocess.run(
            ["gfortran", f"-print-file-name={name}"],
            check=True,
            capture_output=True,
            text=True,
        )
        resolved = Path(result.stdout.strip())
        if resolved.is_file():
            return resolved.parent
    raise RuntimeError("could not locate libgfortran via gfortran -print-file-name")


def main() -> None:
    pharlap_home = Path(os.environ.get("PHARLAP_HOME", "/Users/chartat1/pharlap")).expanduser().resolve()
    if not pharlap_home.is_dir():
        raise RuntimeError(f"PHARLAP_HOME not found: {pharlap_home}")

    system = platform.system()
    if system == "Darwin":
        lib_subdir = "maca"
        output_name = "libiri2020_bridge.dylib"
        compile_mode = ["-dynamiclib"]
    elif system == "Linux":
        lib_subdir = "linux"
        output_name = "libiri2020_bridge.so"
        compile_mode = ["-shared", "-fPIC"]
    else:
        raise RuntimeError(f"unsupported platform for local IRI2020 bridge build: {system}")

    pharlap_lib_dir = pharlap_home / "lib" / lib_subdir
    if not pharlap_lib_dir.is_dir():
        raise RuntimeError(f"PHaRLAP library directory not found: {pharlap_lib_dir}")

    gfortran_lib_dir = _gfortran_library_dir()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / output_name

    command = [
        "cc",
        "-O2",
        "-std=c99",
        *compile_mode,
        str(SOURCE),
        "-o",
        str(output_path),
        "-L",
        str(pharlap_lib_dir),
        "-L",
        str(gfortran_lib_dir),
        "-Wl,-rpath," + str(gfortran_lib_dir),
        "-liri2020",
        "-lmaths",
        "-lgfortran",
        "-lgomp",
        "-lquadmath",
        "-lm",
    ]

    subprocess.run(command, check=True)
    print(output_path)


if __name__ == "__main__":
    main()
