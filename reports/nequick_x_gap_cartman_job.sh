#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
set -euo pipefail
umask 077
ulimit -c 0
run=/homes/chartat1/private_raytracing/runs/nequick_x_gap_recovery
mkdir -p -m 0700 "$run/tmp/continuation" "$run/cache/continuation" "$run/logs" "$run/status"
exec > "$run/logs/continuation.out" 2> "$run/logs/continuation.err"
trap 'result=$?; printf "%s\n" "$result" > "$run/status/continuation.exit"; chmod 0600 "$run/status/continuation.exit"' EXIT
export TMPDIR="$run/tmp/continuation"
export XDG_CACHE_HOME="$run/cache/continuation"
export MPLCONFIGDIR="$run/cache/continuation/matplotlib"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /homes/chartat1/private_raytracing/repo
date -u +%s.%N > "$run/status/continuation.start"
/usr/bin/time -v /homes/chartat1/rt_superdarn_cartman/.venv/bin/python -u - <<'PY'
import sys
import numpy as np
from reports import generate_synthetic_truth_returns as generator
from python_raytrace.multisat_topside_inverse_demo import (
    _home_frequency_returns, _continue_homed_returns, _deduplicate_homed_returns,
)
from python_raytrace.tracer import PointToPointRayTracer

generator.FREQUENCIES = 2.5 + 0.1 * np.arange(5)

def continue_x_branch(**kwargs):
    targets = np.asarray(kwargs["frequencies_mhz"], dtype=float)
    if kwargs["ox_mode"] == 1:
        return tuple(() for _ in targets)
    tracer = PointToPointRayTracer(cache_native_grid=True)
    recovered = [[] for _ in targets]
    for source, direction in ((2.4, 1), (2.6, -1), (2.8, 1), (3.0, -1)):
        paths = tuple(ray for ray in _home_frequency_returns(
            tracer, tx=kwargs["tx"], rx=kwargs["rx"], grid=kwargs["grid"],
            fan_elevations_deg=kwargs["fan_elevations_deg"],
            fan_bearings_deg=kwargs["fan_bearings_deg"],
            frequency_mhz=source, ox_mode=-1, config=kwargs["config"],
            optimizer_method="Powell") if ray.group_range_km >= 150)
        print("seed",source,len(paths),flush=True)
        for step in range(1, 31):
            frequency = round(source + direction * 0.02 * step, 2)
            if frequency < 2.5 - 1e-9 or frequency > 2.9 + 1e-9:
                if (direction < 0 and frequency < 2.5) or (direction > 0 and frequency > 2.9):
                    break
            paths = _continue_homed_returns(
                tracer, paths, tx=kwargs["tx"], rx=kwargs["rx"], grid=kwargs["grid"],
                frequency_mhz=frequency, ox_mode=-1, config=kwargs["config"],
                range_min_km=150.0)
            if not paths:
                print("branch stopped",source,frequency,flush=True)
                break
            match = np.flatnonzero(np.isclose(targets, frequency, atol=1e-8))
            for index in match:
                recovered[int(index)].extend(paths)
    result = tuple(_deduplicate_homed_returns(
        rays, kwargs["config"].homed_max_returns_per_frequency)
                   for rays in recovered)
    print("recovered per bin",[(float(f),len(rays)) for f,rays in zip(targets,result)],flush=True)
    return result

generator.home_frequency_sweep_adaptive = continue_x_branch
sys.argv = ["generate_synthetic_truth_returns.py",
            "--start", "0", "--stop", "5", "--method", "adaptive",
            "--output", "/homes/chartat1/private_raytracing/runs/nequick_x_gap_recovery/continued_x_2p5_to_2p9.npz",
            "--density-grid-npz",
            "/homes/chartat1/private_raytracing/runs/nequick_independent_case/truth_density.npz",
            "--density-scale", "1.0"]
generator.main()
PY
chmod 0600 "$run/continued_x_2p5_to_2p9.npz"
date -u +%s.%N > "$run/status/continuation.finish"
chmod 0600 "$run/status/continuation.start" "$run/status/continuation.finish"
