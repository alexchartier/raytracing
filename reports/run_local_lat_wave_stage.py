"""Run 20 independent ionogram recovery or continuation tasks locally.

    python3 reports/run_local_lat_wave_stage.py recover --source-dir DIR --grid GRID --output-dir DIR
    python3 reports/run_local_lat_wave_stage.py continue --source-dir DIR --grid GRID --output-dir DIR
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


def run(stage: str, source_dir: Path, grid: Path, output_dir: Path,
        workers: int, above_only: bool, indices: list[int] | None = None) -> dict:
    if workers != 1:
        raise ValueError("local memory safety requires exactly one ray worker")
    if stage == "recover" and above_only:
        raise ValueError("above-only applies to continuation")
    indices = list(range(1, 21)) if indices is None else sorted(set(indices))
    if not indices or any(index < 1 or index > 20 for index in indices):
        raise ValueError("profile indices must be between 1 and 20")
    source = [source_dir / f"ionogram_{index:02d}.npz" for index in indices]
    if any(not path.is_file() for path in source) or not grid.is_file():
        raise FileNotFoundError("Expected 20 source ionograms and the forward grid")
    output_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    output_dir.chmod(0o700)
    script = ROOT / "reports" / ("recover_lat_wave_ionogram.py" if stage == "recover"
                                  else "continue_lat_wave_ionogram.py")
    if stage == "recover":
        expected_method = "adaptive_with_dense_gap_recovery"
    else:
        expected_method = ("adaptive_with_dense_gap_recovery_and_20khz_nose_continuation"
                           + ("_and_above_nose_probes" if above_only else ""))
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    started = time.perf_counter()

    def one(index: int) -> dict:
        name = f"ionogram_{index:02d}.npz"
        destination = output_dir / name
        if destination.is_file():
            try:
                with np.load(destination, allow_pickle=False) as data:
                    if (int(data["profile_index"]) == index
                            and len(data["frequencies_mhz"]) == 81
                            and str(data["method"]) == expected_method
                            and len(data["spacecraft_doppler_hz"]) == len(data["records"])):
                        return {"profile": index,
                                "accepted_returns": len(data["records"]),
                                "message": "reused completed output"}
            except (OSError, ValueError, KeyError):
                pass
        command = [sys.executable, str(script), "--source", str(source_dir / name),
                   "--grid", str(grid), "--output", str(destination)]
        if above_only:
            command.append("--above-only")
        result = subprocess.run(command, cwd=ROOT, env=environment,
                                text=True, capture_output=True)
        if result.returncode or not destination.is_file():
            raise RuntimeError(f"Profile {index} failed ({result.returncode}):\n"
                               f"{result.stdout}\n{result.stderr}")
        destination.chmod(0o600)
        with np.load(destination, allow_pickle=False) as data:
            if (int(data["profile_index"]) != index or len(data["frequencies_mhz"]) != 81
                    or str(data["method"]) != expected_method
                    or len(data["spacecraft_doppler_hz"]) != len(data["records"])):
                raise ValueError(f"Invalid output for profile {index}")
            accepted = len(data["records"])
        message = result.stdout.strip()
        if message.startswith("{"):
            detail = json.loads(message)
            message = f"added {detail.get('added_returns', 0)} returns"
        return {"profile": index, "accepted_returns": accepted,
                "message": message}

    completed = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one, index): index for index in indices}
        for future in as_completed(futures):
            item = future.result()
            completed.append(item)
            print(f"{stage}: {item['profile']:02d} complete, {item['accepted_returns']} returns",
                  flush=True)
    result = {
        "stage": stage,
        "above_only": above_only,
        "source_dir": str(source_dir),
        "grid": str(grid),
        "output_dir": str(output_dir),
        "workers": workers,
        "elapsed_seconds": time.perf_counter() - started,
        "profiles": sorted(completed, key=lambda item: item["profile"]),
    }
    (output_dir / "local_stage.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("recover", "continue"))
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--above-only", action="store_true")
    parser.add_argument("--indices", nargs="+", type=int)
    arguments = parser.parse_args()
    print(json.dumps(run(arguments.stage, arguments.source_dir, arguments.grid,
                         arguments.output_dir, arguments.workers,
                         arguments.above_only, arguments.indices), indent=2))
