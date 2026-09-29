#!/bin/bash
#$ -S /bin/bash
# One candidate per SGE task; submit with -tc 2 and -o/-e /dev/null.
set -euo pipefail
umask 077
ulimit -c 0
run=/homes/chartat1/private_raytracing/runs/general_topside_iri_case
case "${SGE_TASK_ID:?}" in
    1) name=basis_fo0995 ;;
    2) name=basis_fo0990 ;;
    3) name=basis_fo0985 ;;
    4) name=basis_fo0990_hmplus5 ;;
    5) name=chapman_o_ridge ;;
    6) name=basis_gradient_050 ;;
    7) name=basis_gradient_000 ;;
    8) name=basis_gradient_015 ;;
    9) name=basis_gradient_030 ;;
    10) name=basis_flat_hm075 ;;
    11) name=basis_flat_hm100 ;;
    12) name=basis_g015_hm075 ;;
    13) name=basis_g015_hm100 ;;
    *) exit 2 ;;
esac
mkdir -p -m 0700 "$run/tmp/$name" "$run/cache/$name"
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
    --density-grid-npz "$run/candidate_${name}_density.npz" \
    --density-scale 1.0
chmod 0600 "$run/full_ray_$name.npz"
date -u +%s.%N > "$run/status/$name.finish"
chmod 0600 "$run/status/$name.start" "$run/status/$name.finish"
