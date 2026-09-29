#!/bin/bash
#$ -S /bin/bash
set -euo pipefail
umask 077
ulimit -c 0
run=/homes/chartat1/private_raytracing/runs/sami3_wave_20170111_0600
code="$run/code"
task=${SGE_TASK_ID:?}
indices=(1 5 9 12 14 17 20)
index=${indices[$((task - 1))]}
stage=${SOUNDER_STAGE:-truth}
printf -v label 'oblique_%s_%02d' "$stage" "$index"
mkdir -p -m 0700 "$run/logs" "$run/status" "$run/tmp/$label" "$run/cache/$label"
exec > "$run/logs/$label.out" 2> "$run/logs/$label.err"
trap 'result=$?; printf "%s\n" "$result" > "$run/status/$label.exit"; chmod 0600 "$run/status/$label.exit"' EXIT
export TMPDIR="$run/tmp/$label"
export XDG_CACHE_HOME="$run/cache/$label"
export MPLCONFIGDIR="$run/cache/$label/matplotlib"
export RAYTRACING_RAY_LOCK="$run/tmp/$label/ray.lock"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SOUNDER_CASE_ROOT="$code/reports/data/sami3_wave_20170111_0600"
export SOUNDER_STAGE="$stage" SOUNDER_INDEX="$index"
cd "$code"
date -u +%s.%N > "$run/status/$label.start"
echo "starting $label on $(hostname) at $(date -u)" >&2
/usr/bin/time -v /homes/chartat1/rt_superdarn_cartman/.venv/bin/python -u - <<'PY'
import json
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd() / 'reports'))
from oblique_wave_pass import trace

case = Path(os.environ['SOUNDER_CASE_ROOT'])
stage = os.environ['SOUNDER_STAGE']
index = int(os.environ['SOUNDER_INDEX'])
base = case / 'oblique_wave_600km'
if stage == 'truth':
    grid = case / 'truth_grid.nc'
    output_dir = base / 'truth'
elif stage == 'baseline':
    grid = Path(os.environ.get('SOUNDER_START_GRID', str(case / 'vertical_selected_grid.nc')))
    output_dir = base / 'baseline'
else:
    grid = base / stage / 'grid.nc'
    output_dir = base / stage / 'ionograms'
output_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
output = output_dir / f'ionogram_{index:02d}.npz'
if not output.exists():
    print(json.dumps(trace(grid, index, output, f'SAMI3 pilot {stage}')))
output.chmod(0o600)
print(output)
PY
date -u +%s.%N > "$run/status/$label.finish"
chmod 0600 "$run/status/$label.start" "$run/status/$label.finish"
