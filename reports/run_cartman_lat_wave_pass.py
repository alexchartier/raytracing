"""Run the 20-profile latitude-wave pass privately on Cartman.

    python3 reports/run_cartman_lat_wave_pass.py submit RUN_NAME
    python3 reports/run_cartman_lat_wave_pass.py status RUN_NAME
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
MANIFEST = ROOT / "reports/data/lat_wave_pass_manifest.json"

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
line=$(sed -n "${task_id}p" "$run/positions.tsv")
read -r latitude longitude altitude <<< "$line"
test -n "$latitude" && test -n "$longitude" && test -n "$altitude"
cd "$repo"
date -u +%s.%N > "$run/status/$task.start"
__PYTHON__ reports/generate_lat_wave_ionogram.py \\
  --grid "$run/wave_forward_grid.nc" \\
  --latitude-deg "$latitude" --longitude-deg "$longitude" --altitude-km "$altitude" \\
  --profile-index "$task_id" \\
  --density-source "__DENSITY_SOURCE__" \\
  --output "$run/results/ionogram_$task.npz"
date -u +%s.%N > "$run/status/$task.finish"
"""

RECOVERY_JOB = """#!/bin/bash
#$ -S /bin/bash
set -euo pipefail
umask 077
run=__RUN__
repo=__ROOT__/repo
task_id=${SGE_TASK_ID:?}
printf -v task '%02d' "$task_id"
mkdir -p -m 0700 "$run/recovery_logs" "$run/recovered" "$run/recovery_status" "$run/recovery_tmp/$task" "$run/recovery_cache/$task"
exec > "$run/recovery_logs/$task.out" 2> "$run/recovery_logs/$task.err"
trap 'code=$?; printf "%s\\n" "$code" > "$run/recovery_status/$task.exit"' EXIT
export TMPDIR="$run/recovery_tmp/$task" XDG_CACHE_HOME="$run/recovery_cache/$task"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$repo"
date -u +%s.%N > "$run/recovery_status/$task.start"
__PYTHON__ reports/recover_lat_wave_ionogram.py \\
  --source "$run/results/ionogram_$task.npz" \\
  --grid "$run/wave_forward_grid.nc" \\
  --output "$run/recovered/ionogram_$task.npz"
date -u +%s.%N > "$run/recovery_status/$task.finish"
"""


def _run_path(name: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_]+", name):
        raise ValueError("Run name must contain only letters, digits, and underscores")
    return f"{REMOTE_ROOT}/runs/{name}"


def submit(name: str, grid_override: Path | None = None,
           max_concurrent: int = 20) -> None:
    if not 1 <= max_concurrent <= 20:
        raise ValueError("max_concurrent must be between 1 and 20")
    manifest = json.loads(MANIFEST.read_text())
    profiles = manifest["profiles"]
    if len(profiles) != 20 or [row["index"] for row in profiles] != list(range(1, 21)):
        raise ValueError("Expected 20 sequential pass positions")
    run = _run_path(name)
    _call_tool("cartman_mkdir", {"zone": "sandbox", "path": name})
    table = "".join(f"{row['latitude_deg']:.12f} {row['longitude_deg']:.12f} "
                    f"{row['altitude_km']:.12f}\n" for row in profiles)
    _call_tool("cartman_write_file", {"zone": "sandbox", "path": f"{name}/positions.tsv",
                                       "text": table})
    source = (ROOT / "reports/generate_lat_wave_ionogram.py").read_text()
    _call_tool("cartman_write_file", {"zone": "repo",
                                       "path": "reports/generate_lat_wave_ionogram.py",
                                       "text": source})
    _call_tool("cartman_write_file", {"zone": "repo",
                                       "path": "python_raytrace/spacecraft_doppler.py",
                                       "text": (ROOT / "python_raytrace/spacecraft_doppler.py").read_text()})
    grid = grid_override if grid_override is not None else ROOT / manifest["forward_grid"]
    with tempfile.TemporaryDirectory(prefix="lat_wave_stage_") as temporary:
        staged = Path(temporary) / "wave_forward_grid.nc"
        shutil.copyfile(grid, staged)
        os.chmod(staged, 0o600)
        transfer = subprocess.run(
            ["rsync", "-a", str(staged), f"cartman:{run}/wave_forward_grid.nc"],
            text=True, capture_output=True,
        )
    if transfer.returncode:
        raise RuntimeError(f"Grid transfer failed: {transfer.stderr}")
    check = _call_tool("cartman_exec_remote", {
        "zone": "sandbox", "cwd": name,
        "command": (f"stat -c '%U %a %n' {run} {run}/wave_forward_grid.nc "
                    f"{REMOTE_ROOT}/repo/reports/generate_lat_wave_ionogram.py "
                    f"{REMOTE_ROOT}/repo/python_raytrace/spacecraft_doppler.py"),
        "timeout_seconds": 120,
    })["structuredContent"]["stdout"]
    if check.splitlines() != [f"chartat1 700 {run}",
                              f"chartat1 600 {run}/wave_forward_grid.nc",
                              f"chartat1 600 {REMOTE_ROOT}/repo/reports/generate_lat_wave_ionogram.py",
                              f"chartat1 600 {REMOTE_ROOT}/repo/python_raytrace/spacecraft_doppler.py"]:
        raise RuntimeError(f"Private staging verification failed:\n{check}")
    source_label = ("Ionogram-selected PyIRI latitude-wave fit" if grid_override is not None
                    else "IRI-2016 with imposed wave")
    script = (JOB.replace("__RUN__", run).replace("__ROOT__", REMOTE_ROOT)
              .replace("__PYTHON__", REMOTE_PYTHON)
              .replace("__DENSITY_SOURCE__", source_label))
    result = _call_tool("cartman_qsub_submit", {
        "zone": "sandbox", "cwd": name, "script_path": f"{name}/job.sh",
        "script_text": script,
        "qsub_args": ["-terse", "-cwd", "-S", "/bin/bash", "-N", "lat_wave_20",
                       "-m", "n", "-t", "1-20", "-tc", str(max_concurrent),
                       "-o", "/dev/null", "-e", "/dev/null"],
        "timeout_seconds": 120,
    })["structuredContent"]
    print(json.dumps({"name": name, "job_id": result["job_id"],
                      "remote_path": run, "staging": check.splitlines()}, indent=2))


def submit_recovery(name: str, task_ids: list[int] | None = None,
                    max_concurrent: int = 20) -> None:
    if not 1 <= max_concurrent <= 20:
        raise ValueError("max_concurrent must be between 1 and 20")
    run = _run_path(name)
    initial = status(name)
    requested = sorted(set(task_ids if task_ids is not None else range(1, 21)))
    if task_ids is not None and (not task_ids or requested != task_ids):
        raise ValueError("Recovery task IDs must be unique and sorted")
    if (not requested or any(not 1 <= task <= 20 for task in requested)
            or initial["expected"] != 20 or initial["privacy_violations"]
            or any(task not in initial["valid_outputs"] for task in requested)):
        raise RuntimeError("Requested forward outputs are incomplete or not private")
    for filename in ("generate_lat_wave_ionogram.py", "recover_lat_wave_ionogram.py"):
        _call_tool("cartman_write_file", {
            "zone": "repo", "path": f"reports/{filename}",
            "text": (ROOT / "reports" / filename).read_text(),
        })
    _call_tool("cartman_write_file", {
        "zone": "repo", "path": "python_raytrace/spacecraft_doppler.py",
        "text": (ROOT / "python_raytrace/spacecraft_doppler.py").read_text(),
    })
    check = _call_tool("cartman_exec_remote", {
        "zone": "sandbox", "cwd": name,
        "command": (f"stat -c '%U %a %n' {REMOTE_ROOT}/repo/reports/generate_lat_wave_ionogram.py "
                    f"{REMOTE_ROOT}/repo/reports/recover_lat_wave_ionogram.py "
                    f"{REMOTE_ROOT}/repo/python_raytrace/spacecraft_doppler.py"),
        "timeout_seconds": 120,
    })["structuredContent"]["stdout"]
    if check.splitlines() != [f"chartat1 600 {REMOTE_ROOT}/repo/reports/generate_lat_wave_ionogram.py",
                              f"chartat1 600 {REMOTE_ROOT}/repo/reports/recover_lat_wave_ionogram.py",
                              f"chartat1 600 {REMOTE_ROOT}/repo/python_raytrace/spacecraft_doppler.py"]:
        raise RuntimeError(f"Private recovery staging verification failed:\n{check}")
    script = RECOVERY_JOB.replace("__RUN__", run).replace("__ROOT__", REMOTE_ROOT).replace("__PYTHON__", REMOTE_PYTHON)
    # This SGE installation accepts only one contiguous range per -t flag.
    ranges = []
    first = previous = requested[0]
    for task in requested[1:]:
        if task != previous + 1:
            ranges.append((first, previous))
            first = task
        previous = task
    ranges.append((first, previous))
    submissions = []
    for first, last in ranges:
        task_range = str(first) if first == last else f"{first}-{last}"
        result = _call_tool("cartman_qsub_submit", {
            "zone": "sandbox", "cwd": name, "script_path": f"{name}/recovery_job.sh",
            "script_text": script,
            "qsub_args": ["-terse", "-cwd", "-S", "/bin/bash", "-N", "lat_wave_recover",
                           "-m", "n", "-t", task_range, "-tc", str(max_concurrent),
                           "-o", "/dev/null", "-e", "/dev/null"],
            "timeout_seconds": 120,
        })["structuredContent"]
        submissions.append({"range": task_range, "job_id": result["job_id"]})
    print(json.dumps({"name": name, "recovery_jobs": submissions,
                      "task_ids": requested}, indent=2))


def status(name: str) -> dict:
    run = _run_path(name)
    command = f"""{REMOTE_PYTHON} - <<'PYREMOTE'
import json, os, stat
from pathlib import Path
import numpy as np
run = Path({run!r})
starts = {{int(p.stem): float(p.read_text()) for p in (run/'status').glob('*.start')}}
finishes = {{int(p.stem): float(p.read_text()) for p in (run/'status').glob('*.finish')}}
exits = {{int(p.stem): int(p.read_text()) for p in (run/'status').glob('*.exit')}}
positions = [tuple(map(float, row.split())) for row in (run/'positions.tsv').read_text().splitlines()
             if row.strip()]
valid = []
return_counts = {{}}
for path in (run/'results').glob('ionogram_*.npz'):
    task = int(path.stem.split('_')[-1])
    if not 1 <= task <= len(positions):
        continue
    with np.load(path, allow_pickle=False) as z:
        coordinates = (float(z['tx_lat_deg']), float(z['tx_lon_deg']), float(z['tx_alt_km']))
        doppler_valid = ('spacecraft_doppler_hz' not in z or
                         (len(z['spacecraft_doppler_hz']) == len(z['records'])
                          and z['launch_angles_deg'].shape == (len(z['records']), 2)
                          and z['arrival_angles_deg'].shape == (len(z['records']), 2)
                          and np.all(np.isfinite(z['spacecraft_doppler_hz']))))
        if (len(z['frequencies_mhz']) == 81 and int(z['profile_index']) == task
                and np.allclose(coordinates, positions[task-1], atol=1e-8)
                and str(z['vertical_fan_layout']) == 'equal_area_guarded'
                and doppler_valid):
            valid.append(task)
            return_counts[task] = len(z['records'])
violations = []
for parent, dirs, names in os.walk(run):
    for path in [Path(parent)] + [Path(parent)/item for item in dirs+names]:
        info = path.lstat()
        if info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) & 0o077:
            violations.append(str(path))
print(json.dumps({{'expected':len(positions),'started':len(starts),'finished':len(finishes),
                  'exit_codes':exits,'valid_outputs':sorted(valid),
                  'return_counts':return_counts,'privacy_violations':violations,
                  'first_start_to_last_finish_seconds':
                    max(finishes.values())-min(starts.values()) if len(finishes)==len(positions) else None}}))
PYREMOTE"""
    result = _call_tool("cartman_exec_remote", {"zone": "sandbox", "cwd": name,
                                                "command": command, "timeout_seconds": 120})["structuredContent"]
    payload = json.loads(result["stdout"])
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def recovery_status(name: str) -> dict:
    run = _run_path(name)
    command = f"""{REMOTE_PYTHON} - <<'PYREMOTE'
import json, os, stat
from pathlib import Path
import numpy as np
run = Path({run!r})
starts = {{int(p.stem): float(p.read_text()) for p in (run/'recovery_status').glob('*.start')}}
finishes = {{int(p.stem): float(p.read_text()) for p in (run/'recovery_status').glob('*.finish')}}
exits = {{int(p.stem): int(p.read_text()) for p in (run/'recovery_status').glob('*.exit')}}
valid = []
recovered_counts = {{}}
for path in (run/'recovered').glob('ionogram_*.npz'):
    task = int(path.stem.split('_')[-1])
    with np.load(path, allow_pickle=False) as z:
        doppler_valid = ('spacecraft_doppler_hz' not in z or
                         (len(z['spacecraft_doppler_hz']) == len(z['records'])
                          and z['launch_angles_deg'].shape == (len(z['records']), 2)
                          and z['arrival_angles_deg'].shape == (len(z['records']), 2)
                          and np.all(np.isfinite(z['spacecraft_doppler_hz']))))
        if (1 <= task <= 20 and int(z['profile_index']) == task
                and str(z['method']) == 'adaptive_with_dense_gap_recovery'
                and len(z['frequencies_mhz']) == 81 and doppler_valid):
            valid.append(task)
            recovered_counts[task] = int(z['dense_recovered_return_count'])
violations = []
for parent, dirs, names in os.walk(run):
    for path in [Path(parent)] + [Path(parent)/item for item in dirs+names]:
        info = path.lstat()
        if info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) & 0o077:
            violations.append(str(path))
print(json.dumps({{'expected':20,'started':len(starts),'finished':len(finishes),
                  'exit_codes':exits,'valid_outputs':sorted(valid),
                  'recovered_return_counts':recovered_counts,
                  'privacy_violations':violations,
                  'first_start_to_last_finish_seconds':
                    max(finishes.values())-min(starts.values()) if len(finishes)==20 else None}}))
PYREMOTE"""
    result = _call_tool("cartman_exec_remote", {"zone": "sandbox", "cwd": name,
                                                "command": command, "timeout_seconds": 120})["structuredContent"]
    payload = json.loads(result["stdout"])
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    submit_parser = sub.add_parser("submit")
    submit_parser.add_argument("run_name")
    submit_parser.add_argument("--grid", type=Path,
                               help="Independent candidate forward grid")
    submit_parser.add_argument("--max-concurrent", type=int, default=20)
    sub.add_parser("status").add_argument("run_name")
    recover_parser = sub.add_parser("recover")
    recover_parser.add_argument("run_name")
    recover_parser.add_argument("--task-ids", type=int, nargs="+",
                                help="Recover completed forward tasks while the rest finish")
    recover_parser.add_argument("--max-concurrent", type=int, default=20)
    sub.add_parser("recovery-status").add_argument("run_name")
    args = parser.parse_args()
    if args.action == "submit":
        submit(args.run_name, args.grid, args.max_concurrent)
    elif args.action == "recover":
        submit_recovery(args.run_name, args.task_ids, args.max_concurrent)
    elif args.action == "recovery-status":
        recovery_status(args.run_name)
    else:
        status(args.run_name)


if __name__ == "__main__":
    main()
