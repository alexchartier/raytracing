This is a minimal local build helper for `pylap` on macOS.

It only builds the `pylap.raytrace_3d` extension because the upstream `PyLap` package is Linux-only and hardcodes PHaRLAP library names that do not match the local macOS install in this environment.

Expected environment:

- `PHARLAP_HOME`: defaults to `/Users/chartat1/pharlap`
- `PYLAP_SOURCE`: defaults to `/tmp/PyLap`
- `GFORTRAN_LIB`: optional override for the GCC runtime library directory

Build/install example:

```bash
export PHARLAP_HOME=/Users/chartat1/pharlap
export PYLAP_SOURCE=/tmp/PyLap
python -m pip install ./python_raytrace/pylap_macos_build
```
