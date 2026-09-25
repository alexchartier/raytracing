"""Run one 40-member D-ionogram inverse population privately on Cartman.

    python3 reports/run_cartman_inverse_population.py submit FIT_DIR ROUND RUN_NAME
    python3 reports/run_cartman_inverse_population.py status RUN_NAME

After all 40 outputs finish, copy RUN_NAME/results into a local private folder
and call ``fit_d_ionogram.py advance`` with that folder.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
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
export TMPDIR="$run/tmp/$task"
export XDG_CACHE_HOME="$run/cache/$task"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
line=$(sed -n "${task_id}p" "$run/candidates.tsv")
read -r density shift <<< "$line"
test -n "$density" && test -n "$shift"
cd "$repo"
date -u +%s.%N > "$run/status/$task.start"
__PYTHON__ reports/generate_synthetic_truth_returns.py --start 0 --stop 81 \
  --density-scale "$density" --hmf2-shift-km "$shift" \
  --output "$run/results/ionogram_$task.npz"
date -u +%s.%N > "$run/status/$task.finish"
"""


def _run_path(name: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_]+", name):
        raise ValueError("Run name must contain only letters, digits, and underscores")
    return f"{REMOTE_ROOT}/runs/{name}"


def submit(workdir: Path, round_number: int, name: str) -> None:
    candidates = json.loads((workdir / f"candidates_round_{round_number}.json").read_text())
    if len(candidates) != 40 or [c["task"] for c in candidates] != list(range(1, 41)):
        raise ValueError("Expected tasks 1–40")
    path = _run_path(name)
    _call_tool("cartman_mkdir", {"zone": "sandbox", "path": name})
    table = "".join(f"{c['density_scale']:.12f} {c['hmf2_shift_km']:.12f}\n" for c in candidates)
    _call_tool("cartman_write_file", {"zone": "sandbox", "path": f"{name}/candidates.tsv", "text": table})
    # Keep the exact local D generator in the private remote checkout.
    source = (ROOT / "reports" / "generate_synthetic_truth_returns.py").read_text()
    _call_tool("cartman_write_file", {"zone": "repo", "path": "reports/generate_synthetic_truth_returns.py", "text": source})
    script = JOB.replace("__RUN__", path).replace("__ROOT__", REMOTE_ROOT).replace("__PYTHON__", REMOTE_PYTHON)
    result = _call_tool("cartman_qsub_submit", {
        "zone": "sandbox", "cwd": name, "script_path": f"{name}/job.sh",
        "script_text": script,
        "qsub_args": ["-terse", "-cwd", "-S", "/bin/bash", "-N", "inverse_D_40",
                      "-m", "n", "-t", "1-40", "-tc", "40",
                      "-o", "/dev/null", "-e", "/dev/null"],
        "timeout_seconds": 120,
    })["structuredContent"]
    print(json.dumps({"name": name, "job_id": result["job_id"], "remote_path": path}, indent=2))


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
outputs = sorted((run/'results').glob('ionogram_*.npz'))
valid = []
for path in outputs:
    with np.load(path, allow_pickle=False) as z:
        if len(z['frequencies_mhz']) == 81 and 'hmf2_shift_km' in z:
            valid.append(int(path.stem.split('_')[-1]))
violations = []
for parent, dirs, names in os.walk(run):
    for path in [Path(parent)] + [Path(parent)/name for name in dirs+names]:
        info = path.lstat()
        if info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) & 0o077:
            violations.append(str(path))
print(json.dumps({{'started':len(starts),'finished':len(finishes), 'exit_codes':exits,
                  'valid_outputs':valid, 'privacy_violations':violations,
                  'first_start_to_last_finish_seconds':
                    max(finishes.values())-min(starts.values()) if len(finishes)==40 else None}}))
PYREMOTE"""
    result = _call_tool("cartman_exec_remote", {"zone": "sandbox", "cwd": name,
                                                "command": command, "timeout_seconds": 120})["structuredContent"]
    payload = json.loads(result["stdout"])
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    sub_submit = sub.add_parser("submit")
    sub_submit.add_argument("workdir", type=Path)
    sub_submit.add_argument("round_number", type=int, choices=(1, 2, 3))
    sub_submit.add_argument("run_name")
    sub_status = sub.add_parser("status")
    sub_status.add_argument("run_name")
    args = parser.parse_args()
    if args.action == "submit":
        submit(args.workdir, args.round_number, args.run_name)
    else:
        status(args.run_name)


if __name__ == "__main__":
    main()
