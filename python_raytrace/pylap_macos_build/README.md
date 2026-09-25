This is a minimal local build helper for `pylap` on macOS.

It only builds the `pylap.raytrace_3d` extension because the upstream `PyLap` package is Linux-only and hardcodes PHaRLAP library names that do not match the local macOS install in this environment.

The build helper applies [`patches/cached_state_vector.patch`](patches/cached_state_vector.patch)
to a temporary copy of the upstream source. The patch makes the 10-argument
call accept a ray state vector after the grids have been loaded, matching
MATLAB's `raytrace_3d(..., tol, ray_state_vec_in)` calling sequence. It also
corrects two upstream return-type errors required by current macOS Clang.
The upstream checkout is left unchanged.

The adaptive sweep reuses the native grid through `PyLapRaytraceBackend`.
PyLap stores one grid per process, so the backend reloads it when the grid
object changes. If code calls `pylap.raytrace_3d` directly in the same process
or edits a grid object in place, call
`PyLapRaytraceBackend.invalidate_native_grid_cache()` afterward before resuming
a cached sweep.

Expected environment:

- `PHARLAP_HOME`: defaults to `/Users/chartat1/pharlap`
- `PYLAP_SOURCE`: defaults to `/tmp/PyLap`
- `GFORTRAN_LIB`: optional override for the GCC runtime library directory

Build/install example:

```bash
export PHARLAP_HOME=/Users/chartat1/pharlap
export PYLAP_SOURCE=/tmp/PyLap
git clone https://github.com/HamSCI/PyLap.git "$PYLAP_SOURCE"
python3 -m pip install --user --no-build-isolation --no-deps ./python_raytrace/pylap_macos_build
```
