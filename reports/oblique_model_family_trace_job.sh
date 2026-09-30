#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -tc 2
set -euo pipefail
umask 077
ulimit -c 0

run=/homes/chartat1/private_raytracing/runs/oblique_model_families
code="$run/code"
task=${SGE_TASK_ID:?Submit as a task array}
[[ "$task" =~ ^[0-9]+$ ]] || exit 2
read -r case_name candidate < <(/homes/chartat1/rt_superdarn_cartman/.venv/bin/python - "$run/tasks_trace.json" "$task" <<'PY'
import json, sys
row = json.load(open(sys.argv[1]))[int(sys.argv[2]) - 1]
print(row['case'], row['candidate'])
PY
)
[[ "$case_name" =~ ^[a-z0-9_]+$ && "$candidate" =~ ^[a-z0-9_]+$ ]] || exit 2
label="${case_name}_${candidate}"
case_root="$code/reports/data/oblique_model_families/$case_name"
grid="$case_root/$candidate/grid.nc"
[[ "$candidate" == truth ]] && grid="$case_root/truth_grid.nc"
output="$case_root/forward/$candidate/ionogram_01.npz"
if [[ -s "$run/status/$label.exit" && -s "$output" ]]; then
    [[ $(cat "$run/status/$label.exit") == 0 ]] && exit 0
fi
mkdir -p -m 0700 "$run/logs" "$run/status" "$run/tmp/$label" \
    "$run/cache/$label" "$(dirname "$output")"
exec > "$run/logs/$label.out" 2> "$run/logs/$label.err"
trap 'result=$?; printf "%s\n" "$result" > "$run/status/$label.exit"; chmod 0600 "$run/status/$label.exit"' EXIT
export TMPDIR="$run/tmp/$label"
export XDG_CACHE_HOME="$run/cache/$label"
export MPLCONFIGDIR="$run/cache/$label/matplotlib"
export RAYTRACING_RAY_LOCK="$run/tmp/$label/ray.lock"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SOUNDER_CASE_ROOT="$case_root"
export SOUNDER_GRID="$grid" SOUNDER_OUTPUT="$output"
cd "$code"
date -u +%s.%N > "$run/status/$label.start"
/usr/bin/time -v /homes/chartat1/rt_superdarn_cartman/.venv/bin/python -u - <<'PY'
import json, os, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path.cwd() / 'reports'))
from oblique_wave_pass import trace

output = Path(os.environ['SOUNDER_OUTPUT'])
print(json.dumps(trace(Path(os.environ['SOUNDER_GRID']), 1, output,
                       'independent model-profile oblique test')))
output.chmod(0o600)
with np.load(output, allow_pickle=False) as result:
    if (len(result['frequencies_mhz']) != 81
            or str(result['method']) != 'oblique_adaptive'
            or abs(float(result['satellite_separation_km']) - 600.0) > .01
            or len(result['records']) != len(result['spacecraft_doppler_hz'])):
        raise ValueError('Incomplete oblique ionogram')
PY
date -u +%s.%N > "$run/status/$label.finish"
chmod 0600 "$run/status/$label.start" "$run/status/$label.finish" \
    "$run/logs/$label.out" "$run/logs/$label.err"
[[ -z $(find "$run" ! -user chartat1 -print -quit) &&
   -z $(find "$run" -type f -perm /077 -print -quit) &&
   -z $(find "$run" -type d -perm /077 -print -quit) ]] || exit 1
