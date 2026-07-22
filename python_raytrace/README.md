This directory contains a Python port of the point-to-point HF raytracing flow used by the MATLAB code in this repository.

It mirrors the `raytrace_3dhome` approach:

- build a regional ionosphere and geomagnetic grid
- launch a coarse fan of rays
- pick the closest candidate
- numerically home the launch elevation and bearing onto the receiver

The port swaps the original dependencies as requested:

- `pylap` in place of PHaRLaP
- `PyIRI` in place of SAMI3 / PHaRLaP IRI grid generation

The current grid builder also adds the missing absorption inputs:

- `pymsis` for neutral atmosphere densities
- local `iri2020`/`FPT-2018` profiles for D-region electron density plus `Ti`/`Te`
- a PHaRLAP-compatible collision-frequency calculation instead of the previous zeroed grid

Files:

- `grid.py`: build a regional PyIRI electron density grid, merge in low-altitude IRI2020 D-region structure, and compute collision frequencies
- `geometry.py`: WGS84 helpers, great-circle sampling, and point-to-ray distance
- `tracer.py`: the point-to-point solver and a thin `pylap.raytrace_3d` backend
- `cli.py`: command line entrypoint for single-link and frequency-sweep runs
- `webapp.py`: prototype local web visualizer and JSON API
- `iri2020_build/`: local bridge build for the PHaRLAP `iri2020` Fortran core used for D-region/temperature profiles

Minimal usage:

```bash
python3 -m python_raytrace.cli \
  --time 2020-01-15T12:00:00 \
  --tx -77.8 166.4 0.0 \
  --rx -89.9 166.4 1.0 \
  --freq 4.1 --freq 5.1 --freq 6.0 \
  --f107 100 \
  --alt-min-km 60 \
  --alt-max-km 500 \
  --d-region-model fpt2018
```

Prototype web UI:

```bash
python3 -m python_raytrace.webapp --host 127.0.0.1 --port 8000
```

Then open `http://127.0.0.1:8000/`.

Programmatic usage:

```python
from datetime import datetime

from python_raytrace import GeoPoint, PointToPointRayTracer

tracer = PointToPointRayTracer()
results = tracer.trace_frequencies(
    when=datetime(2020, 1, 15, 12, 0, 0),
    tx=GeoPoint(-77.8, 166.4, 0.0),
    rx=GeoPoint(-89.9, 166.4, 1.0),
    frequencies_mhz=[4.1, 5.1, 6.0],
    f107=100.0,
)
```

Ground-to-space homing example:

```bash
python3 -m python_raytrace.ground_to_space_homing_example --mode synthetic
```

This example launches from a ground station, homes onto a 550 km target in full 3D, and prints the solved launch angles plus miss distance as JSON. Use `--mode real` to swap in the existing `PyIRI` + `pylap` pipeline when those local dependencies are available.

Collm Skiymet time-series predictor:

```bash
python3 -m python_raytrace.collm_skiymet_timeseries \
  --tle-file ~/superdarn/digital_rf_tools/artifacts/iss_stations.tle \
  --measurement-csv ~/superdarn/digital_rf_tools/tmp/pulsed_meteor_radar_detection_METnwDEU_Collm_Skiymet_36.2000MHz.csv \
  --fit-geometry \
  --fit-snr \
  --fit-start 2026-06-01T15:58:50.5+00:00 \
  --start 2026-06-01T15:58:35.5+00:00 \
  --seconds 88 \
  --step-seconds 1 \
  --f107 100 \
  --f107a 100
```

This generates a Collm-specific three-panel prediction matching the detector conventions:

- top: one-way range folded by the 625 Hz PRI
- middle: Doppler aliased into the PRF band
- bottom: PHaRLAP-based peak SNR using homed-ray path loss, a local ray-tube expansion metric, and D-region/F-region absorption

It also writes a CSV with the underlying per-time geometry and homing-ray solution, a separate foF2 `pcolormesh` map at closest approach with the ISS ground track and transmitter overlaid, and a CPA along-track versus altitude swath plot with the homed ray path.

The Collm predictor now always uses the real `PyIRI` + `pylap`/PHaRLAP-style raytracing path. There is no synthetic straight-line mode in this workflow.

Single-wave TID fitting:

The ionospheric perturbation is a single oriented Hooke-style wave applied multiplicatively to electron density:

```text
deltaNe/Ne = A * exp(-0.5 * ((z - z0)/sigma_z)^2) * cos(2*pi*s/lambda - 2*pi*t/T + phi)
```

- `A`: fractional electron-density amplitude
- `lambda`: horizontal wavelength in km
- `s`: horizontal coordinate measured along the chosen wave-bearing direction
- `T`: temporal period in seconds
- `phi`: phase offset in radians
- `z0`, `sigma_z`: Gaussian vertical center and thickness in km

The wave orientation can be given either as an absolute azimuth with `--hooke-bearing-deg` or as an offset from the ISS ground-track bearing with `--hooke-bearing-offset-deg`. The offset form is useful when you want the same relative geometry to remain valid across nearby passes.

The current best fit-window-only TID recipe used for the Collm/ISS example after the dropout at `2026-06-01T15:58:50.5+00:00` is:

```bash
python3 -m python_raytrace.collm_skiymet_timeseries \
  --tle-file ~/superdarn/digital_rf_tools/artifacts/iss_stations.tle \
  --measurement-csv ~/superdarn/digital_rf_tools/tmp/pulsed_meteor_radar_detection_METnwDEU_Collm_Skiymet_36.2000MHz.csv \
  --fit-geometry \
  --fit-snr \
  --fit-start 2026-06-01T15:58:50.5+00:00 \
  --start 2026-06-01T15:58:50.5+00:00 \
  --seconds 53 \
  --step-seconds 1 \
  --f107 100 \
  --f107a 100 \
  --hooke-amplitude-fraction 0.50 \
  --hooke-horizontal-wavelength-km 82 \
  --hooke-bearing-offset-deg -45 \
  --hooke-period-seconds 12.6 \
  --hooke-phase-rad 0 \
  --hooke-vertical-center-km 285 \
  --hooke-vertical-sigma-km 40
```

That fit uses:

- `time_shift = 18.439 s`
- `range_offset = -149.933 km`
- `doppler_offset = 50.228 Hz`
- `SNR fit RMS = 0.758 dB` on the post-dropout fit window

The TID fit is intentionally restricted to the post-dropout interval, but the script can still simulate and plot the full observed time span by keeping `--start` earlier than `--fit-start`.

Notes:

- `PyIRI` densities are converted from `m^-3` to `cm^-3` before being passed to `pylap`, because `pylap.raytrace_3d` expects electron density in electrons per cubic centimeter.
- The `PyIRI` source field is evaluated on a full-world lattice aligned to the requested horizontal step size, then the regional raytracing subgrid is extracted in the original unwrapped longitude frame.
- If `f107` and `Ap` are not supplied, the code now resolves them from `pymsis`’s cached space-weather dataset for the requested UTC timestamp; `--refresh-indices` forces a fresh download of that cache before the run.
- For absorption work, `pymsis` neutral densities are converted from `m^-3` to `cm^-3` before the collision-frequency calculation, matching the PHaRLAP formulas.
- The low-altitude density blend uses local `iri2020` output with the `FPT-2018` D-region option and a 120-140 km transition back to the `PyIRI` grid.
- The default model volume used by the solver is now `60-500 km`, and the prototype web app searches up to `2` hops.
- `pylap` is imported lazily. Grid-building tests can run without a working compiled `pylap` install.
- The local `iri2020` bridge must be built once with `python3 python_raytrace/iri2020_build/build.py`. It links against the existing PHaRLAP static libraries and writes the shared library under `python_raytrace/_lib/`.
