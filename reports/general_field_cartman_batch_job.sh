#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -tc 4
set -euo pipefail
umask 077

run=/homes/chartat1/private_raytracing/runs/general_field_iri_wave
code="$run/code"
task=${SGE_TASK_ID:?Submit 1-36 after the single-profile memory check}
[[ "$task" =~ ^[0-9]+$ ]] || exit 2
(( task >= 1 && task <= 36 )) || exit 2
names=(baseline plus_peak minus_peak plus_height minus_height plus_width minus_width plus_wave minus_wave)
indices=(3 8 14 18)
export SOUNDER_CANDIDATE=${names[$(((task - 1) / 4))]}
export SOUNDER_INDEX=${indices[$(((task - 1) % 4))]}
label=$(printf '%s_%02d' "$SOUNDER_CANDIDATE" "$SOUNDER_INDEX")
filename=$(printf 'ionogram_%02d.npz' "$SOUNDER_INDEX")
if [[ -s "$run/status/$label.exit" && -s "$run/$SOUNDER_CANDIDATE/ionograms/$filename" ]]; then
    [[ $(cat "$run/status/$label.exit") == 0 ]] && exit 0
fi
exec /bin/bash "$code/reports/general_field_cartman_job.sh"
