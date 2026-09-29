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

## General topside inverse

`python_raytrace/general_topside_inverse.py` fits foF2, hmF2, and four
IRI-derived topside shape coefficients, with a generalized Chapman branch as
an alternative. `general_topside_fit.py` builds the 100-profile global
IRI-2020 basis, fits the saved O/X ionogram, generates candidate density
grids, and evaluates candidates only after full 3-D ray tracing. The basis
builder does not read the IRI-2016 validation density. The existing private
Cartman generator runs the candidates; `general_topside_cartman_job.sh`
documents the restricted two-at-a-time refinement job. The reproducible
scalar shape and family check is `validate_general_topside_basis.py`.

From the repository root, run `python3 reports/general_topside_fit.py
build-basis`, then `python3 reports/general_topside_fit.py fit` and
`python3 reports/validate_general_topside_basis.py`. After private full-ray
candidate jobs finish, `refine` proposes peak and horizontal-gradient checks;
`height-refine` uses saved O/X range residuals for a bounded height step.
Run `evaluate` only after the final full-ray candidate jobs finish.
The [method and results report](general_topside_inverse_report.md) separates
idealized vertical-range validation from full 3-D O/X validation.

## Analytic Chapman cross-family check

`prepare_chapman_truth.py` generates a uniform generalized-Chapman truth grid
without importing the retrieval or using IRI density. After its private full-ray
ionogram is copied back, `chapman_model_validation.py fit` builds independent
Chapman and IRI-basis candidates from O/X returns; `refine` derives a height
step from the first Chapman full-ray range residual; `evaluate` selects using
saved full-ray ionograms before opening truth density. The jobs are defined in
`chapman_cartman_job.sh`. The [Chapman validation report](chapman_model_validation_report.pdf)
shows scalar and full-ray errors, paired O/X ionograms, and density cuts;
[Markdown source](chapman_model_validation_report.md) is beside it.

## Independent NeQuick-G profile

`prepare_nequick_truth.py` freezes a public NeQuick-G model profile from
`tpl2go/NequickG` commit `1d1783412ed9811540b56d409a1dcf27d2413740`.
The upstream Python 2 source is converted mechanically with `lib2to3` in an
ignored cache; no NeQuick code or profile enters the retrieval's IRI basis,
Chapman branch, or monotone spline. The generated density and provenance are
saved under `data/nequick_independent_case/`. `nequick_independent_fit.py fit`
reads the **ionogram only**, `refine` reads saved ray returns, and `evaluate`
opens the NeQuick density after selecting by O/X score. Full rays run through
`nequick_cartman_job.sh` in `chartat1`'s private server area.

The first out-of-family 3-D check exposed a limitation in the topside fit.
Its initial O/X scores and 10.92% topside RMS are **provisional**: the adaptive
truth sweep missed every X return at 2.5–2.9 MHz. Dense-fan seeds followed at
20 kHz steps recovered ten accepted X rays at the exact plotted 100 kHz
frequencies. The [before/after truth figure](data/nequick_independent_case/nequick_truth_x_gap_recovered.png)
and [gap recovery note](nequick_x_gap_recovery.md) show the correction. The
original `evaluation.json` used the incomplete truth and must not be cited as
a score against the corrected truth. `generate_synthetic_truth_returns.py` is
restored here because it is the active full-ray forward generator.

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
