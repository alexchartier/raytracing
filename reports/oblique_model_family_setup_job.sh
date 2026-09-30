#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -tc 2
set -euo pipefail
umask 077
ulimit -c 0

run=/homes/chartat1/private_raytracing/runs/oblique_model_families
case "${SGE_TASK_ID:?}" in
    1) label=support; args=(support) ;;
    2) label=truth_nequick; args=(truth --case nequick) ;;
    3) label=truth_chapman; args=(truth --case chapman) ;;
    4) label=truth_iri2016; args=(truth --case iri2016) ;;
    5) label=start_nequick; args=(start --case nequick) ;;
    6) label=start_chapman; args=(start --case chapman) ;;
    7) label=start_iri2016; args=(start --case iri2016) ;;
    *) exit 2 ;;
esac
mkdir -p -m 0700 "$run/logs" "$run/status" "$run/tmp/$label" "$run/cache/$label"
exec > "$run/logs/$label.out" 2> "$run/logs/$label.err"
trap 'result=$?; printf "%s\n" "$result" > "$run/status/$label.exit"; chmod 0600 "$run/status/$label.exit"' EXIT
export TMPDIR="$run/tmp/$label"
export XDG_CACHE_HOME="$run/cache/$label"
export MPLCONFIGDIR="$run/cache/$label/matplotlib"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$run/code"
date -u +%s.%N > "$run/status/$label.start"
/usr/bin/time -v /homes/chartat1/rt_superdarn_cartman/.venv/bin/python -u \
    reports/oblique_model_family_cases.py "${args[@]}"
date -u +%s.%N > "$run/status/$label.finish"
chmod 0600 "$run/status/$label.start" "$run/status/$label.finish" \
    "$run/logs/$label.out" "$run/logs/$label.err"
[[ -z $(find "$run" ! -user chartat1 -print -quit) &&
   -z $(find "$run" -type f -perm /077 -print -quit) &&
   -z $(find "$run" -type d -perm /077 -print -quit) ]] || exit 1
