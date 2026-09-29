#!/bin/bash
#$ -S /bin/bash
set -euo pipefail
umask 077
ulimit -c 0
run=/homes/chartat1/private_raytracing/runs/sami3_wave_20170111_0600
code="$run/code"
task=${SGE_TASK_ID:?}
indices=(1 5 9 12 14 17 20)
if (( task <= 7 )); then
    kind=truth
    offset=$((task - 1))
else
    kind=prior
    offset=$((task - 8))
fi
index=${indices[$offset]}
printf -v label '%s_%02d' "$kind" "$index"
mkdir -p -m 0700 "$run/logs" "$run/status" "$run/tmp/$label" "$run/cache/$label"
exec > "$run/logs/$label.out" 2> "$run/logs/$label.err"
trap 'result=$?; printf "%s\n" "$result" > "$run/status/$label.exit"; chmod 0600 "$run/status/$label.exit"' EXIT
export TMPDIR="$run/tmp/$label"
export XDG_CACHE_HOME="$run/cache/$label"
export MPLCONFIGDIR="$run/cache/$label/matplotlib"
export RAYTRACING_RAY_LOCK="$run/tmp/$label/ray.lock"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SOUNDER_KIND="$kind" SOUNDER_INDEX="$index"
cd "$code"
date -u +%s.%N > "$run/status/$label.start"
echo "starting $label on $(hostname) at $(date -u)" >&2
/usr/bin/time -v /homes/chartat1/rt_superdarn_cartman/.venv/bin/python -u - <<'PY'
import json
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd() / 'reports'))
from generate_lat_wave_ionogram import generate

base = Path('reports/data/sami3_wave_20170111_0600')
index = int(os.environ['SOUNDER_INDEX'])
kind = os.environ['SOUNDER_KIND']
row = json.loads((base / 'manifest.json').read_text())['profiles'][index - 1]
output_dir = base / ('vertical_truth' if kind == 'truth' else 'vertical_prior')
output_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
output = output_dir / f'ionogram_{index:02d}.npz'
if not output.exists():
    generate(base / ('truth_grid.nc' if kind == 'truth' else 'pyiri_prior_grid.nc'),
             float(row['latitude_deg']), float(row['longitude_deg']),
             float(row['altitude_km']), index, output,
             'SAMI3 2017-01-11 snapshot' if kind == 'truth' else 'independent 2017 PyIRI prior')
output.chmod(0o600)
print(output)
PY
date -u +%s.%N > "$run/status/$label.finish"
chmod 0600 "$run/status/$label.start" "$run/status/$label.finish"
