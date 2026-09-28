"""Run private 20-profile ionogram nose continuation on Cartman.

    python3 reports/run_cartman_lat_wave_continuation.py submit RUN --source-dir DIR --grid GRID
    python3 reports/run_cartman_lat_wave_continuation.py status RUN
    python3 reports/run_cartman_lat_wave_continuation.py fetch RUN --output-dir DIR
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from tools.cartman_mcp.server import _call_tool  # noqa: E402

REMOTE_ROOT = "/homes/chartat1/private_raytracing"
REMOTE_PYTHON = "/homes/chartat1/rt_superdarn_cartman/.venv/bin/python"

JOB = """#!/bin/bash
#$ -S /bin/bash
set -euo pipefail
umask 077
run=__RUN__
repo=__ROOT__/repo
task_id=${SGE_TASK_ID:?}
printf -v task '%02d' "$task_id"
mkdir -p -m 0700 "$run/logs" "$run/results" "$run/status" "$run/tmp/$task" "$run/cache/$task"
exec > "$run/logs/$task.out" 2> "$run/logs/$task.err"
trap 'code=$?; printf "%s\\n" "$code" > "$run/status/$task.exit"' EXIT
export TMPDIR="$run/tmp/$task" XDG_CACHE_HOME="$run/cache/$task"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$repo"
date -u +%s.%N > "$run/status/$task.start"
__PYTHON__ reports/continue_lat_wave_ionogram.py \\
  --source "$run/input/ionogram_$task.npz" \\
  --grid "$run/grid.nc" \\
  --output "$run/results/ionogram_$task.npz" __ABOVE_ONLY__
date -u +%s.%N > "$run/status/$task.finish"
"""


def run_path(name: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_]+", name):
        raise ValueError("Run name must contain only letters, digits, and underscores")
    return f"{REMOTE_ROOT}/runs/{name}"


def remote(command: str, name: str) -> str:
    result = _call_tool("cartman_exec_remote", {
        "zone": "sandbox", "cwd": name, "command": command,
        "timeout_seconds": 120,
    })["structuredContent"]
    if result["returncode"]:
        raise RuntimeError(result["stderr"] or result["stdout"])
    return result["stdout"]


def submit(name: str, source_dir: Path, grid: Path, max_concurrent: int,
           above_only: bool = False) -> None:
    if not 1 <= max_concurrent <= 20:
        raise ValueError("max-concurrent must be between 1 and 20")
    expected = [source_dir / f"ionogram_{index:02d}.npz" for index in range(1, 21)]
    if any(not path.is_file() for path in expected):
        raise ValueError("Expected exactly 20 existing source ionograms")
    if not grid.is_file():
        raise FileNotFoundError(grid)
    run = run_path(name)
    _call_tool("cartman_mkdir", {"zone": "sandbox", "path": name})
    for subdir in ("input", "results", "logs", "status", "tmp", "cache"):
        _call_tool("cartman_mkdir", {"zone": "sandbox", "path": f"{name}/{subdir}"})
    with tempfile.TemporaryDirectory(prefix="lat_wave_continue_stage_") as temporary:
        staged = Path(temporary)
        os.chmod(staged, 0o700)
        for path in expected + [grid]:
            target = staged / ("grid.nc" if path == grid else path.name)
            shutil.copyfile(path, target)
            os.chmod(target, 0o600)
        for sources, target in (([staged / path.name for path in expected],
                                  f"cartman:{run}/input/"),
                                 ([staged / "grid.nc"], f"cartman:{run}/grid.nc")):
            result = subprocess.run(["rsync", "-a", *map(str, sources), target],
                                    text=True, capture_output=True)
            if result.returncode:
                raise RuntimeError(f"Private input staging failed: {result.stderr}")
    for path in ("reports/continue_lat_wave_ionogram.py",
                 "reports/generate_lat_wave_ionogram.py",
                 "python_raytrace/spacecraft_doppler.py"):
        _call_tool("cartman_write_file", {"zone": "repo", "path": path,
                                           "text": (ROOT / path).read_text()})
    check = remote(
        f"stat -c '%U %a %n' {run} {run}/grid.nc "
        f"{REMOTE_ROOT}/repo/reports/continue_lat_wave_ionogram.py "
        f"{REMOTE_ROOT}/repo/reports/generate_lat_wave_ionogram.py "
        f"{REMOTE_ROOT}/repo/python_raytrace/spacecraft_doppler.py "
        f"{run}/input/ionogram_*.npz", name)
    lines = check.splitlines()
    if (len(lines) != 25 or lines[0] != f"chartat1 700 {run}"
            or any(line.split(" ", 2)[:2] != ["chartat1", "600"]
                   for line in lines[1:])):
        raise RuntimeError(f"Private staging verification failed:\n{check}")
    script = (JOB.replace("__RUN__", run).replace("__ROOT__", REMOTE_ROOT)
              .replace("__PYTHON__", REMOTE_PYTHON)
              .replace("__ABOVE_ONLY__", "--above-only" if above_only else ""))
    result = _call_tool("cartman_qsub_submit", {
        "zone": "sandbox", "cwd": name, "script_path": f"{name}/job.sh",
        "script_text": script,
        "qsub_args": ["-terse", "-cwd", "-S", "/bin/bash", "-N", "lat_wave_continue",
                       "-m", "n", "-t", "1-20", "-tc", str(max_concurrent),
                       "-o", "/dev/null", "-e", "/dev/null"],
        "timeout_seconds": 120,
    })["structuredContent"]
    print(json.dumps({"run": name, "job_id": result["job_id"],
                      "above_only": above_only,
                      "remote_path": run, "private_staging_entries": len(lines)}, indent=2))


def status(name: str) -> dict:
    run = run_path(name)
    command = f"""{REMOTE_PYTHON} - <<'PYREMOTE'
import json, os, stat
from pathlib import Path
import numpy as np
run = Path({run!r})
starts = {{int(p.stem): float(p.read_text()) for p in (run/'status').glob('*.start')}}
finishes = {{int(p.stem): float(p.read_text()) for p in (run/'status').glob('*.finish')}}
exits = {{int(p.stem): int(p.read_text()) for p in (run/'status').glob('*.exit')}}
valid = []
additions = {{}}
for path in (run/'results').glob('ionogram_*.npz'):
    task = int(path.stem.split('_')[-1])
    with np.load(path, allow_pickle=False) as z:
        records = z['records']
        if (1 <= task <= 20 and int(z['profile_index']) == task
                and len(z['frequencies_mhz']) == 81
                and len(z['spacecraft_doppler_hz']) == len(records)
                and int(z['above_nose_added_return_count'] if 'above_nose_added_return_count' in z
                        else z['nose_refinement_added_return_count']) >= 0
                and np.all(records[:,3] <= 1000.0 + 1e-6)
                and np.all(records[:,2] >= 150.0)):
            valid.append(task)
            additions[task] = int(z['above_nose_added_return_count'] if 'above_nose_added_return_count' in z
                                  else z['nose_refinement_added_return_count'])
violations = []
for parent, dirs, names in os.walk(run):
    for path in [Path(parent)] + [Path(parent)/item for item in dirs+names]:
        info = path.lstat()
        if info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) & 0o077:
            violations.append(str(path))
print(json.dumps({{'expected':20,'started':len(starts),'finished':len(finishes),
                  'exit_codes':exits,'valid_outputs':sorted(valid),
                  'added_returns_by_profile':additions,
                  'privacy_violations':violations,
                  'first_start_to_last_finish_seconds':
                    max(finishes.values())-min(starts.values()) if len(finishes)==20 else None}}))
PYREMOTE"""
    payload = json.loads(remote(command, name))
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def fetch(name: str, output_dir: Path) -> None:
    state = status(name)
    if (state["valid_outputs"] != list(range(1, 21)) or state["privacy_violations"]
            or any(code != 0 for code in state["exit_codes"].values())):
        raise RuntimeError("Cannot fetch incomplete or nonprivate run")
    output_dir.mkdir(parents=True, exist_ok=True)
    os.chmod(output_dir, 0o700)
    result = subprocess.run(["rsync", "-a",
                             f"cartman:{run_path(name)}/results/", str(output_dir) + "/"],
                            text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr)
    files = list(output_dir.glob("ionogram_*.npz"))
    if len(files) != 20 or any(path.stat().st_mode & 0o077 for path in files):
        raise RuntimeError("Fetched results failed file count or permission check")
    print(json.dumps({"output_dir": str(output_dir), "files": len(files)}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    submit_parser = sub.add_parser("submit")
    submit_parser.add_argument("run_name")
    submit_parser.add_argument("--source-dir", type=Path, required=True)
    submit_parser.add_argument("--grid", type=Path, required=True)
    submit_parser.add_argument("--max-concurrent", type=int, default=20)
    submit_parser.add_argument("--above-only", action="store_true")
    sub.add_parser("status").add_argument("run_name")
    fetch_parser = sub.add_parser("fetch")
    fetch_parser.add_argument("run_name")
    fetch_parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.action == "submit":
        submit(args.run_name, args.source_dir, args.grid, args.max_concurrent,
               args.above_only)
    elif args.action == "status":
        status(args.run_name)
    else:
        fetch(args.run_name, args.output_dir)


if __name__ == "__main__":
    main()
