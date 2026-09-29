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

The combined [vertical and oblique report](vertical_oblique_retrieval_report.pdf)
documents the density prior, fitting stages, scoring rules, and provenance,
then shows paired truth/retrieved ionograms, latitude–altitude density sections,
representative altitude cuts, and foF2/hmF2 evaluations. Its source is
[`vertical_oblique_retrieval_report.md`](vertical_oblique_retrieval_report.md);
`build_vertical_oblique_report_figures.py` regenerates the altitude cuts and
peak evaluation from the frozen selections and saved density grids.

## SAMI3 wave pilot

`setup_sami3_wave_case.py` prepares the frozen 2017-01-11 06:00 SAMI3 trough and
crest near 190°E as density truth, alongside an independent, date-matched
PyIRI prior. `retrieve_sami3_wave.py` infers peak and topside corrections from
O/X vertical ionograms and selects them by forward-traced ionogram score.
`oblique_wave_pass.py` supports the same case via `SOUNDER_CASE_ROOT` and fits
600 km two-satellite links after the vertical stage. The Cartman array scripts
`sami3_wave_cartman_vertical_job.sh`, `sami3_wave_cartman_candidate_job.sh`,
and `sami3_wave_cartman_oblique_job.sh` run this large grid in `chartat1`'s
private server area. `plot_sami3_wave_results.py` makes paired O/X ionograms
and density sections after selection. The full SAMI3 grid
must not be ray-traced on the 16 GB laptop; `local_ray_lock.py` enforces this.
The [case report](sami3_wave_retrieval_report.pdf) gives the truth and retrieved
ionograms, density cuts, peak metrics, method, and validation limits; its
[Markdown source](sami3_wave_retrieval_report.md) is kept beside the PDF.
Rebuild it from the repository root with
`pandoc reports/sami3_wave_retrieval_report.md --pdf-engine=tectonic --resource-path=reports -o reports/sami3_wave_retrieval_report.pdf`.
