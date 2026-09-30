#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -tc 2
set -euo pipefail
umask 077
ulimit -c 0

run=/homes/chartat1/private_raytracing/runs/general_field_iri_wave
code="$run/code"
name=${SOUNDER_CANDIDATE:?Set a candidate from the private plan}
index=${SOUNDER_INDEX:-${SGE_TASK_ID:-}}
[[ -n "$index" ]] || { echo "Set SOUNDER_INDEX or submit an array" >&2; exit 2; }
[[ "$name" =~ ^[A-Za-z0-9_]+$ ]] || exit 2
[[ "$index" =~ ^[0-9]+$ ]] || exit 2
label=$(printf '%s_%02d' "$name" "$index")
filename=$(printf 'ionogram_%02d.npz' "$index")
if [[ -s "$run/status/$label.exit" && -s "$run/$name/ionograms/$filename" ]]; then
    [[ $(cat "$run/status/$label.exit") == 0 ]] && exit 0
fi

mkdir -p -m 0700 "$run/logs" "$run/status" "$run/tmp/$label" "$run/cache/$label"
for stage in raw recovered continued ionograms; do
    mkdir -p -m 0700 "$run/$name/$stage"
done
exec > "$run/logs/$label.out" 2> "$run/logs/$label.err"
trap 'result=$?; printf "%s\n" "$result" > "$run/status/$label.exit"; chmod 0600 "$run/status/$label.exit"' EXIT
export TMPDIR="$run/tmp/$label"
export XDG_CACHE_HOME="$run/cache/$label"
export MPLCONFIGDIR="$run/cache/$label/matplotlib"
export RAYTRACING_RAY_LOCK="$run/tmp/$label/ray.lock"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SOUNDER_CANDIDATE="$name" SOUNDER_INDEX="$index" SOUNDER_RUN="$run"

cd "$code"
date -u +%s.%N > "$run/status/$label.start"
/usr/bin/time -v /homes/chartat1/rt_superdarn_cartman/.venv/bin/python -u - <<'PY'
import json
import os
import subprocess
import sys
from pathlib import Path

code = Path.cwd()
plan = json.loads((code / 'reports/data/general_field_iri_wave_structured/plan.json').read_text())
name = os.environ['SOUNDER_CANDIDATE']
index = int(os.environ['SOUNDER_INDEX'])
candidate = next((row for row in plan['candidates'] if row['name'] == name), None)
if candidate is None and name == 'combined':
    proposal_path = code / 'reports/data/general_field_iri_wave_structured/combined/proposal.json'
    if proposal_path.is_file():
        candidate = json.loads(proposal_path.read_text())
if candidate is None:
    raise ValueError(f'{name} is not in the private candidate plan')
profiles = json.loads((code / plan['manifest']).read_text())['profiles']
profile = next((row for row in profiles if int(row['index']) == index), None)
if profile is None:
    raise ValueError(f'Profile {index} is not in the manifest')
grid = code / candidate['grid']
if not grid.is_file():
    raise FileNotFoundError(grid)
run = Path(os.environ['SOUNDER_RUN'])
paths = {stage: run / name / stage / f'ionogram_{index:02d}.npz'
         for stage in ('raw', 'recovered', 'continued', 'ionograms')}
python = sys.executable
jobs = [
    [python, 'reports/generate_lat_wave_ionogram.py', '--grid', str(grid),
     '--latitude-deg', str(profile['latitude_deg']),
     '--longitude-deg', str(profile['longitude_deg']),
     '--altitude-km', str(profile['altitude_km']),
     '--profile-index', str(index), '--output', str(paths['raw']),
     '--density-source', 'general log-density field candidate'],
    [python, 'reports/recover_lat_wave_ionogram.py',
     '--source', str(paths['raw']), '--grid', str(grid),
     '--output', str(paths['recovered'])],
    [python, 'reports/continue_lat_wave_ionogram.py',
     '--source', str(paths['recovered']), '--grid', str(grid),
     '--output', str(paths['continued'])],
    [python, 'reports/continue_lat_wave_ionogram.py',
     '--source', str(paths['continued']), '--grid', str(grid),
     '--output', str(paths['ionograms']), '--above-only'],
]
for command, path in zip(jobs, paths.values()):
    subprocess.run(command, cwd=code, check=True)
    path.chmod(0o600)
import numpy as np
with np.load(paths['ionograms'], allow_pickle=False) as result:
    if (len(result['frequencies_mhz']) != 81
            or len(result['records']) != len(result['spacecraft_doppler_hz'])
            or str(result['method']) !=
            'adaptive_with_dense_gap_recovery_and_20khz_nose_continuation_and_above_nose_probes'):
        raise ValueError('Incomplete final ionogram')
PY
date -u +%s.%N > "$run/status/$label.finish"
chmod 0600 "$run/status/$label.start" "$run/status/$label.finish" "$run/logs/$label.out" "$run/logs/$label.err"
bad_owner=$(find "$run" ! -user chartat1 -print -quit)
bad_file=$(find "$run" -type f -perm /077 -print -quit)
bad_dir=$(find "$run" -type d -perm /077 -print -quit)
[[ -z "$bad_owner" && -z "$bad_file" && -z "$bad_dir" ]] || {
    echo "Private-run ownership or permissions failed" >&2
    exit 1
}
