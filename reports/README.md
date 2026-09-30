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
`iri_wave_spline_check.py` is a later, profile-level regression screen using
the saved IRI-wave ionograms and a synthetic 800 km in-situ density scalar at
each of the 20 positions. Run its `prepare`, `fit`, `tail`, and `evaluate`
stages in that order. The eight-slope O-mode spline is **rejected** on this
three-dimensional wave: peak-density MAE rises from 1.73% to 3.17%, and four
fits miss the local density by over 40%. A conservative 600–800 km tail anchor
preserves the previous peak result and lowers 300–800 km peak-normalized RMS
from 1.291% to 1.261%, but has not received full 3-D ray selection.
Replaying the frozen 20-position ionogram selector still chooses `wide` with
the original 0.211682 combined score. The saved final ionograms retain 29
truth and 34 retrieved interior empty mode/frequency bins after the earlier
recovery passes; the newer inline gap repair has not been ray-checked on this
wave case.
The [screen data and figure](data/iri_wave_spline_local_density/) are separate
from the frozen report result.
The oblique experiment and its commands are documented in
[`oblique_wave_600km.md`](oblique_wave_600km.md).

## General density-field retrieval in development

`python_raytrace/general_field_inverse.py` applies a smooth, bounded correction
to **log electron density** on any supplied forward grid. It uses Gaussian
functions across the pass and at heights relative to each prior F2 peak, so
the same parameterization can change peak density, peak height, and topside
shape. `reports/general_field_retrieval.py` prepares paired pilot grids,
compares every accepted O/X return after complete ray tracing, and proposes a
regularized trust-region step. At each candidate, an observed in-situ density
at the spacecraft is imposed with a 600–800 km taper (for the current 800 km
pass). The fit reads no truth density beyond those explicit local scalars.

The first IRI-wave development plan has **9 grids** (one anchored baseline
and paired 5% log-density changes along peak, height, width, and an
ionogram-inferred wave direction) and four pilot positions: **36 ionograms**.
The measured low-Doppler O/X cutoff pattern inferred a **902 km** wavelength
without reading the wave truth parameters. The plan is saved at
`data/general_field_iri_wave_structured/plan.json`. All 36 pilot ionograms
and both 20-position finalist passes were traced on Cartman with the same
complete O/X, gap-recovery, and 20 kHz nose-continuation pipeline as the
observations. The private Cartman scripts are `general_field_cartman_job.sh`
and `general_field_cartman_batch_job.sh`. The run tree passed ownership and
permission checks; peak resident memory was 8.9 GiB per ray job.

The [pilot scores](data/general_field_iri_wave_structured/pilot_scores.json)
chose `minus_wave`; the paired scores proposed a combined step. The
[20-position scores](data/general_field_iri_wave_structured/full_scores.json)
and [frozen ionogram-only selection](data/general_field_iri_wave_structured/selection.json)
chose `minus_wave` before opening full truth density. Lower scores are better.
The density results are in the two evaluation JSON files in that directory.

| Grid | Pilot score | Full score | Peak-density MAE | 220–600 km NRMSE | 220–800 km NRMSE | hmF2 MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Previous `wide` | 0.22908 | 0.21168 | 1.73% | 4.19% | 3.44% | 12 km |
| Selected `minus_wave` | 0.19326 | 0.21074 | 1.58% | 3.83% | 3.14% | 12 km |
| Rejected combined step | 0.20671 | 0.22379 | 2.88% | 3.33% | 2.73% | 12 km |

The selected score gain is **0.00094**, with a paired-position bootstrap 95%
interval of **−0.0245 to +0.0227**; its peak-density MAE gain is also small.
The mean hmF2 bias remains **+12 km**, reaching **+40 km** at one position.
The combined step improves broad profile RMS but worsens the peak and the
ionogram score, so it was rejected. The 800 km density is an assimilated
synthetic in-situ observation; its 0.056% error in the new grids is constraint
consistency, not independent validation. Full density truth is separately
generated Fortran IRI-2016 with an imposed wave, while candidates use a
PyIRI-derived prior and shared PyLap ray physics. Earlier development examined
this same truth case, so these are exploratory regression results, not blind
evidence that the method generalizes to another ionosphere model.

The selected [density comparison](data/general_field_iri_wave_structured/density_comparison_minus_wave.png)
and [paired ionograms](data/general_field_iri_wave_structured/ionogram_comparison_minus_wave.png)
show representative positions. Additional rounds can use a denser basis with
a selected subset of active spatial centers per round.

`check_general_field_capacity.py` is an **oracle capacity diagnostic**, not a
retrieval: it fits coefficients directly to full density truth. With a denser
12 × 8 basis, peak-normalized 220–800 km RMS falls from 3.44% to 0.32% for
the IRI wave and from 13.24% to 1.28% for an independent analytic Chapman
wave. The corresponding [capacity data](data/general_field_iri_wave_structured/capacity.json)
show that the representation can fit both shapes. These oracle fits use full
truth density, so their numbers do not establish observable retrieval accuracy.

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

The [independent-model report](independent_nequick_retrieval_report.pdf) and
[Markdown source](independent_nequick_retrieval_report.md) document a separate,
prospectively selected April NeQuick-G profile. The eight-slope monotone
spline was fit from saved O/X returns without using NeQuick density, then eight
candidate profiles were checked with full O/X rays. The minimum-score candidate
has foF2 error −0.032 MHz, hmF2 error −13.9 km, peak density error −1.1%, and
peak-normalized topside RMS 7.9% through 600 km; its 800 km density is 42.5%
low. The [paired ionograms and density cut](data/nequick_prospective_case/nequick_prospective_case_validation.png)
show both the near-peak fit and upper-tail limitation. This controlled test
replicates one independent vertical profile over the ray grid; it does not
test horizontal structure. `nequick_heldout_fit.py` with `--run
reports/data/nequick_prospective_case` fits or evaluates the saved case;
`recover_missing_ionogram_returns.py` performs the exact-frequency gap check
before full-ray scoring. The July NeQuick-G case in
`data/nequick_heldout_2026_case/` was used during method development, so its
scores are not a prospective validation result.

The complete adaptive forward generator now checks for empty interior
mode-frequency bins before writing an ionogram. It continues already accepted
rays through those gaps at 20 kHz steps and uses a full fan only where
continuation fails. The 1 km homing gate and exact 100 kHz output frequencies
are unchanged. `--no-gap-recovery` disables the check for a baseline timing
run; short frequency chunks retain their existing behavior.
In a [controlled five-bin X-gap test](data/nequick_independent_case/inline_gap_recovery_validation.json),
the inline pass recovered 15 exact-frequency accepted returns in 2.8 seconds
of a 132-second full sweep, with no dense-fan fallback. With no interior gaps,
the pass took less than 1 ms on the April candidate sweeps.

`prepare_nequick_local_density_case.py` freezes one synthetic 800 km in situ
density observation from the April NeQuick-G truth. The local-density branch
of `nequick_heldout_fit.py fit --run reports/data/nequick_local_density_case`
reads that scalar and the saved O/X ionogram, without opening the profile or
peak parameters. `height-refine` uses the first round's full-ray O/X range
residual to propose bounded peak-height steps while preserving the measured
800 km density. `evaluate` selects by full-ray score before opening the
remaining truth. `nequick_local_density_cartman_job.sh` traces the candidates
privately on Cartman. The [paired result](data/nequick_local_density_case/nequick_local_density_case_validation.png)
and [direct A/B comparison](data/nequick_local_density_case/nequick_local_density_comparison.png)
show the selected profile. Score improved from 0.181 to 0.151; peak-to-600 km
normalized topside RMS fell from 7.9% to 1.8%, and hmF2 error moved from
−13.9 to +3.4 km. The 800 km agreement is imposed by the synthetic scalar
measurement; relative density RMS from 600 to 800 km is still 6.5%. This
April case was already inspected during the prior study, so these updated
results are a development-case check rather than a new held-out validation.

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
