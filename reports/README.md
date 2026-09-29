# Active retrieval scripts

The Python files in this directory support the current latitude-wave retrieval
and the 600 km two-satellite oblique experiment. Earlier one-off report,
benchmark, and Cartman scripts were removed after commit `29ca74a`; that commit
contains the historical source needed to reproduce those older investigations.
The older Markdown reports and saved results remain as records.

- `setup_latitudinal_wave_pass.py` creates the independent IRI-2016 wave truth
  and forward grid.
- `generate_lat_wave_ionogram.py`, `recover_lat_wave_ionogram.py`, and
  `continue_lat_wave_ionogram.py` produce the vertical O/X observations.
- `run_local_lat_wave_forward.py` and `run_local_lat_wave_stage.py` run those
  independent jobs on local workers.
- `build_lat_wave_doppler_peak_round3.py`,
  `select_lat_wave_doppler_peak_round3.py`, and
  `evaluate_lat_wave_doppler_peak_round3.py` reproduce the current vertical
  retrieval from saved prior inputs.
- `ionogram_metrics.py` supplies the accepted-return score used by both
  selectors.
- `oblique_wave_pass.py` simulates, fits, selects, and evaluates the 600 km
  two-satellite links using the saved wave and vertical retrieval as inputs.

The current vertical report is
[`lat_wave_doppler_peak_round3_2026-09-28.md`](lat_wave_doppler_peak_round3_2026-09-28.md).
The oblique experiment and its commands are documented in
[`oblique_wave_600km.md`](oblique_wave_600km.md).
