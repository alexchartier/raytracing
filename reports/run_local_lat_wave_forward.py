"""Generate a complete option-D latitude-pass ionogram set on local workers.

Completed, validated ionograms are reused so an interrupted pass can resume.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "reports/data/lat_wave_pass_manifest.json"


def generate(grid: Path, output_dir: Path, workers: int, indices: list[int]) -> dict:
    if not 1 <= workers <= 12:
        raise ValueError("workers must be between 1 and 12")
    profiles = {int(row["index"]): row
                for row in json.loads(MANIFEST.read_text())["profiles"]}
    if not indices or any(index not in profiles for index in indices):
        raise ValueError("Expected profile indices in 1..20")
    if not grid.is_file():
        raise FileNotFoundError(grid)
    output_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    started = time.perf_counter()

    def one(index: int) -> dict:
        destination = output_dir / f"ionogram_{index:02d}.npz"
        if destination.is_file():
            try:
                with np.load(destination, allow_pickle=False) as data:
                    if (int(data["profile_index"]) == index
                            and str(data["method"]) == "adaptive"
                            and str(data["vertical_fan_option"]) == "D"
                            and len(data["frequencies_mhz"]) == 81
                            and len(data["spacecraft_doppler_hz"]) == len(data["records"])):
                        return {"index": index, "returns": len(data["records"]),
                                "reused": True}
            except (OSError, ValueError, KeyError):
                pass
        row = profiles[index]
        command = [sys.executable, str(ROOT / "reports/generate_lat_wave_ionogram.py"),
                   "--grid", str(grid), "--latitude-deg", str(row["latitude_deg"]),
                   "--longitude-deg", str(row["longitude_deg"]),
                   "--altitude-km", str(row["altitude_km"]),
                   "--profile-index", str(index), "--output", str(destination),
                   "--density-source", "second-round joint wave candidate"]
        result = subprocess.run(command, cwd=ROOT, env=environment,
                                capture_output=True, text=True)
        if result.returncode or not destination.is_file():
            raise RuntimeError(f"Profile {index} failed ({result.returncode}):\n"
                               f"{result.stdout}\n{result.stderr}")
        destination.chmod(0o600)
        with np.load(destination, allow_pickle=False) as data:
            if (int(data["profile_index"]) != index
                    or str(data["method"]) != "adaptive"
                    or len(data["frequencies_mhz"]) != 81):
                raise ValueError(f"Invalid output: {destination}")
            return {"index": index, "returns": len(data["records"]),
                    "runtime_seconds": float(data["runtime_seconds"]), "reused": False}

    completed = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one, index): index for index in indices}
        for future in as_completed(futures):
            item = future.result()
            completed.append(item)
            print(f"forward: {item['index']:02d} complete, {item['returns']} returns",
                  flush=True)
    result = {"grid": str(grid), "output_dir": str(output_dir),
              "workers": workers, "elapsed_seconds": time.perf_counter() - started,
              "profiles": sorted(completed, key=lambda row: row["index"])}
    (output_dir / "local_forward.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--indices", nargs="+", type=int, default=list(range(1, 21)))
    arguments = parser.parse_args()
    print(json.dumps(generate(arguments.grid, arguments.output_dir,
                              arguments.workers, arguments.indices), indent=2))
