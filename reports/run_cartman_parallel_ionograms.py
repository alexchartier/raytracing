"""Submit and summarize 40 private D-setting Cartman ionograms.

Usage:
  python3 reports/run_cartman_parallel_ionograms.py submit
  python3 reports/run_cartman_parallel_ionograms.py status RUN_NAME

The SGE array requests up to 40 simultaneous one-slot tasks. Each task uses a
different synthetic density scale (0.90–1.29) and writes only under the
chartat1-owned, mode-0700 Cartman run directory.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.cartman_mcp.server import _call_tool  # noqa: E402

REMOTE_ROOT = "/homes/chartat1/private_raytracing"
REMOTE_PYTHON = "/homes/chartat1/rt_superdarn_cartman/.venv/bin/python"

JOB_TEMPLATE = """#!/bin/bash
#$ -S /bin/bash
set -euo pipefail
umask 077
run=__RUN__
repo=__ROOT__/repo
task_id=${SGE_TASK_ID:?}
printf -v task '%02d' "$task_id"
mkdir -p -m 0700 "$run/logs" "$run/results" "$run/started" "$run/finished" "$run/status"
mkdir -p -m 0700 "$run/tmp/$task" "$run/cache/$task"
exec > "$run/logs/$task.out" 2> "$run/logs/$task.err"
trap 'code=$?; printf "%s\\n" "$code" > "$run/status/$task.exit"' EXIT
export TMPDIR="$run/tmp/$task"
export XDG_CACHE_HOME="$run/cache/$task"
export MPLCONFIGDIR="$run/cache/$task/matplotlib"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
scale=$(awk -v id="$task_id" 'BEGIN {printf "%.2f", 0.89 + 0.01 * id}')
cd "$repo"
date -u +%s.%N > "$run/started/$task.epoch"
__PYTHON__ reports/generate_synthetic_truth_returns.py --start 0 --stop 81 \\
  --density-scale "$scale" --output "$run/results/ionogram_$task.npz"
date -u +%s.%N > "$run/finished/$task.epoch"
"""


def run_path(name: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_]+", name):
        raise ValueError("Run name must contain only letters, digits, and underscores")
    return f"{REMOTE_ROOT}/runs/{name}"


def remote(command: str, *, run_name: str, timeout_seconds: int = 120) -> str:
    result = _call_tool("cartman_exec_remote", {
        "zone": "sandbox", "cwd": run_name, "command": command,
        "timeout_seconds": timeout_seconds,
    })["structuredContent"]
    return str(result["stdout"])


def submit(name: str) -> None:
    run = run_path(name)
    _call_tool("cartman_mkdir", {"zone": "sandbox", "path": name})
    remote(f"test ! -e {run}/submitted.epoch", run_name=name)
    remote(f"date -u +%s.%N > {run}/submitted.epoch", run_name=name)
    script = (JOB_TEMPLATE.replace("__RUN__", run)
              .replace("__ROOT__", REMOTE_ROOT)
              .replace("__PYTHON__", REMOTE_PYTHON))
    result = _call_tool("cartman_qsub_submit", {
        "zone": "sandbox", "cwd": name, "script_path": f"{name}/job.sh",
        "script_text": script,
        "qsub_args": [
            "-terse", "-cwd", "-S", "/bin/bash", "-N", "ionogram_D_40",
            "-m", "n", "-t", "1-40", "-tc", "40",
            "-o", "/dev/null", "-e", "/dev/null",
        ],
        "timeout_seconds": 120,
    })["structuredContent"]
    print(json.dumps({"run_name": name, "remote_path": run,
                      "job_id": result["job_id"], "task_count": 40,
                      "max_simultaneous_tasks": 40}, indent=2))


def status(name: str) -> dict:
    run = run_path(name)
    inspector = f"""{REMOTE_PYTHON} - <<'PYREMOTE'
import json, os, stat
from pathlib import Path
import numpy as np
run = Path({run!r})
def times(subdir):
    return {{int(path.stem): float(path.read_text()) for path in (run / subdir).glob('*.epoch')}}
started = times('started')
finished = times('finished')
exit_codes = {{int(path.stem): int(path.read_text()) for path in (run / 'status').glob('*.exit')}}
submitted = float((run / 'submitted.epoch').read_text())
files = sorted((run / 'results').glob('ionogram_*.npz'))
checks = []
for path in files:
    with np.load(path, allow_pickle=False) as data:
        task = int(path.stem.split('_')[-1])
        checks.append({{'task': task, 'returns': len(data['records']),
                       'density_scale': float(data['density_scale']),
                       'fan': int(data['fan_launch_directions']),
                       'outer_fraction': float(data['vertical_outer_ray_fraction']),
                       'guard_seed_limit': int(data['vertical_guard_seed_limit']),
                       'frequencies': len(data['frequencies_mhz']),
                       'runtime_seconds': float(data['runtime_seconds'])}})
events = sorted([(value, 1) for value in started.values()] +
                [(value, -1) for value in finished.values()], key=lambda item: (item[0], -item[1]))
active = peak = 0
for _, change in events:
    active += change
    peak = max(peak, active)
privacy_violations = []
for parent, dirs, names in os.walk(run):
    for path in [Path(parent)] + [Path(parent) / item for item in dirs + names]:
        info = path.lstat()
        if info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) & 0o077:
            privacy_violations.append(str(path))
payload = {{'run_name': {name!r}, 'remote_path': str(run),
           'started_tasks': len(started), 'finished_tasks': len(finished),
           'exit_codes': exit_codes, 'output_count': len(files),
           'peak_overlap': peak, 'outputs': checks,
           'privacy_violations': privacy_violations,
           'submission_to_last_finish_seconds':
               max(finished.values()) - submitted if len(finished) == 40 else None,
           'first_start_to_last_finish_seconds':
               max(finished.values()) - min(started.values()) if len(finished) == 40 else None}}
payload['validated_outputs'] = sorted(item['task'] for item in checks
    if item['frequencies'] == 81 and item['fan'] == 162
    and item['outer_fraction'] == 0.5 and item['guard_seed_limit'] == 4
    and abs(item['density_scale'] - (0.89 + 0.01 * item['task'])) < 1e-8)
print(json.dumps(payload, sort_keys=True))
PYREMOTE"""
    result = json.loads(remote(inspector, run_name=name, timeout_seconds=120))
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("submit", "status"))
    parser.add_argument("run_name", nargs="?")
    args = parser.parse_args()
    if args.action == "submit":
        name = args.run_name or "parallel_D_40_" + datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        submit(name)
    elif args.run_name is None:
        parser.error("status requires RUN_NAME")
    else:
        status(args.run_name)


if __name__ == "__main__":
    main()
