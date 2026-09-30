#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -tc 2
set -euo pipefail
umask 077
ulimit -c 0

run=/homes/chartat1/private_raytracing/runs/sami3_monotone_topside
code="$run/code"
task=${SGE_TASK_ID:?Submit as a task array}
[[ "$task" =~ ^[0-9]+$ ]] || exit 2
read -r name index < <(/homes/chartat1/rt_superdarn_cartman/.venv/bin/python - "$run/tasks.json" "$task" <<'PY'
import json, sys
tasks = json.load(open(sys.argv[1]))
task = tasks[int(sys.argv[2]) - 1]
print(task['candidate'], task['profile_index'])
PY
)
[[ "$name" =~ ^[A-Za-z0-9_]+$ && "$index" =~ ^[0-9]+$ ]] || exit 2
label=$(printf '%s_%02d' "$name" "$index")
filename=$(printf 'ionogram_%02d.npz' "$index")
if [[ -s "$run/status/$label.exit" && -s "$run/$name/ionograms/$filename" ]]; then
    [[ $(cat "$run/status/$label.exit") == 0 ]] && exit 0
fi
mkdir -p -m 0700 "$run/logs" "$run/status" "$run/tmp/$label" "$run/cache/$label" "$run/$name/ionograms"
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
import json, os, subprocess, sys
from pathlib import Path
import numpy as np

code = Path.cwd()
plan = json.loads((code / 'reports/data/sami3_monotone_topside/plan.json').read_text())
name = os.environ['SOUNDER_CANDIDATE']
index = int(os.environ['SOUNDER_INDEX'])
candidate = next((row for row in plan['candidates'] if row['name'] == name), None)
if candidate is None:
    raise ValueError(f'Unknown candidate {name}')
profile = next((row for row in json.loads((code / plan['manifest']).read_text())['profiles']
                if int(row['index']) == index), None)
if profile is None:
    raise ValueError(f'Unknown profile {index}')
grid = code / candidate['grid']
output = Path(os.environ['SOUNDER_RUN']) / name / 'ionograms' / f'ionogram_{index:02d}.npz'
subprocess.run([sys.executable, 'reports/generate_lat_wave_ionogram.py',
                '--grid', str(grid), '--latitude-deg', str(profile['latitude_deg']),
                '--longitude-deg', str(profile['longitude_deg']),
                '--altitude-km', str(profile['altitude_km']),
                '--profile-index', str(index), '--output', str(output),
                '--density-source', 'spatially varying monotone topside retrieval'], check=True)
output.chmod(0o600)
with np.load(output, allow_pickle=False) as result:
    if (len(result['frequencies_mhz']) != 81
            or len(result['records']) != len(result['spacecraft_doppler_hz'])
            or str(result['method']) != 'adaptive'):
        raise ValueError('Incomplete raw adaptive ionogram')
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
