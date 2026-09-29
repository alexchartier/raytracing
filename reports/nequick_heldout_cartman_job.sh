#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
set -euo pipefail
umask 077
ulimit -c 0
run=/homes/chartat1/private_raytracing/runs/nequick_heldout_2026_case
case "${SGE_TASK_ID:?}" in
    1) name=truth ;;
    2) name=legacy_spline ;;
    3) name=tail_spline ;;
    4) name=x_anchor_low ;;
    5) name=x_anchor_mid ;;
    6) name=x_anchor_high ;;
    *) exit 2 ;;
esac
mkdir -p -m 0700 "$run/tmp/$name" "$run/cache/$name" "$run/logs" "$run/status"
exec > "$run/logs/$name.out" 2> "$run/logs/$name.err"
trap 'result=$?; printf "%s\n" "$result" > "$run/status/$name.exit"; chmod 0600 "$run/status/$name.exit"' EXIT
export TMPDIR="$run/tmp/$name"
export XDG_CACHE_HOME="$run/cache/$name"
export MPLCONFIGDIR="$run/cache/$name/matplotlib"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /homes/chartat1/private_raytracing/repo
date -u +%s.%N > "$run/status/$name.start"
/usr/bin/time -v /homes/chartat1/rt_superdarn_cartman/.venv/bin/python -u \
    reports/generate_synthetic_truth_returns.py \
    --start 0 --stop 81 \
    --output "$run/full_ray_$name.npz" \
    --density-grid-npz "$run/${name}_density.npz" \
    --density-scale 1.0
chmod 0600 "$run/full_ray_$name.npz"
date -u +%s.%N > "$run/status/$name.finish"
chmod 0600 "$run/status/$name.start" "$run/status/$name.finish"
