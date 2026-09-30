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
read -r candidate index < <(/homes/chartat1/rt_superdarn_cartman/.venv/bin/python - "$run/tasks_oblique.json" "$task" <<'PY'
import json, sys
row = json.load(open(sys.argv[1]))[int(sys.argv[2]) - 1]
print(row['candidate'], int(row['profile_index']))
PY
)
[[ "$candidate" =~ ^[A-Za-z0-9_]+$ && "$index" =~ ^[0-9]+$ ]] || exit 2
label=$(printf 'oblique_%s_%02d' "$candidate" "$index")
filename=$(printf 'ionogram_%02d.npz' "$index")
output="$run/oblique/$candidate/ionograms/$filename"
if [[ -s "$run/status/$label.exit" && -s "$output" ]]; then
    [[ $(cat "$run/status/$label.exit") == 0 ]] && exit 0
fi
mkdir -p -m 0700 "$run/logs" "$run/status" "$run/tmp/$label" "$run/cache/$label" "$(dirname "$output")"
exec > "$run/logs/$label.out" 2> "$run/logs/$label.err"
trap 'result=$?; printf "%s\n" "$result" > "$run/status/$label.exit"; chmod 0600 "$run/status/$label.exit"' EXIT
export TMPDIR="$run/tmp/$label"
export XDG_CACHE_HOME="$run/cache/$label"
export MPLCONFIGDIR="$run/cache/$label/matplotlib"
export RAYTRACING_RAY_LOCK="$run/tmp/$label/ray.lock"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SOUNDER_CASE_ROOT="$code/reports/data/sami3_wave_20170111_0600"
export SOUNDER_CANDIDATE="$candidate" SOUNDER_INDEX="$index" SOUNDER_OUTPUT="$output"
cd "$code"
date -u +%s.%N > "$run/status/$label.start"
/usr/bin/time -v /homes/chartat1/rt_superdarn_cartman/.venv/bin/python -u - <<'PY'
import json, os, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path.cwd() / 'reports'))
from oblique_wave_pass import trace

index = int(os.environ['SOUNDER_INDEX'])
output = Path(os.environ['SOUNDER_OUTPUT'])
name = os.environ['SOUNDER_CANDIDATE']
grid = Path.cwd() / 'reports/data/sami3_monotone_oblique' / name / 'grid.nc'
if not grid.is_file():
    grid = Path.cwd() / 'reports/data/sami3_monotone_topside' / name / 'grid.nc'
print(json.dumps(trace(grid, index, output,
                       'SAMI3 spatially varying monotone topside retrieval')))
output.chmod(0o600)
with np.load(output, allow_pickle=False) as result:
    if (len(result['frequencies_mhz']) != 81
            or str(result['method']) != 'oblique_adaptive'
            or abs(float(result['satellite_separation_km']) - 600.0) > 0.01
            or len(result['records']) != len(result['spacecraft_doppler_hz'])):
        raise ValueError('Incomplete 600 km oblique ionogram')
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
