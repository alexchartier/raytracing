This is a local build of `pylap.raytrace_3d` on macOS.

The required PyLap C source is included in [`vendor/pylap`](vendor/pylap).
It is pinned to HamSCI/PyLap commit `4b32d0a` and distributed under its MIT
license. This build only includes the `raytrace_3d` extension because the
upstream build assumes Linux PHaRLAP library names.

The vendored `raytrace_3d.c` accepts a 10-argument call with a ray state vector
after the grids have been loaded, matching MATLAB's
`raytrace_3d(..., tol, ray_state_vec_in)` sequence. It also has two return-type
corrections required by current macOS Clang. See [`vendor/pylap/UPSTREAM.md`](vendor/pylap/UPSTREAM.md)
for the source list and changes.

The adaptive sweep reuses the native grid through `PyLapRaytraceBackend`.
PyLap stores one grid per process, so the backend reloads it when the grid
object changes. If code calls `pylap.raytrace_3d` directly in the same process
or edits a grid object in place, call
`PyLapRaytraceBackend.invalidate_native_grid_cache()` afterward before resuming
a cached sweep.

Expected environment:

- `PHARLAP_HOME`: defaults to `/Users/chartat1/pharlap`
- `GFORTRAN_LIB`: optional override for the GCC runtime library directory

Build/install example:

```bash
export PHARLAP_HOME=/Users/chartat1/pharlap
python3 -m pip install --user --no-build-isolation --no-deps ./python_raytrace/pylap_macos_build
```
