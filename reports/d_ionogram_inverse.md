# D ionogram inverse fit

`fit_d_ionogram.py` fits two parameters of the D forward model: electron
density scale and a uniform F2-height shift. A single vertical ionogram has
little independent information about the wave's horizontal bearing and phase,
so the fit fixes wave amplitude at zero. The forward generator accepts those
wave parameters for later multi-view experiments.

## Search

The initial population is 40 precomputed D ionograms with density scales
0.90–1.29 and zero height shift. Their O/X maximum reflected frequencies
calibrate the observed frequency nose against the *same* ray tracer. The
observed nose is a proxy for critical frequency: a return at the 10 MHz sweep
edge is only a lower bound, and homing gaps can move the accepted-return nose
by a frequency bin. A median low-frequency range residual supplies a first
height guess using the approximate two-way 2 km range change per 1 km of
reflector shift.

Three more populations of 40 candidates refine those guesses. The first is
stratified around the nose and height estimates, with eight additional local
draws. The next two use the best eight candidates and shrink their spread;
each retains eight wider draws. All 160 forward ionograms are evaluated by
the same fixed score and the best evaluated candidate is kept.

The score compares **every accepted return** for each O/X mode. It combines
a symmetric nearest-return distance in frequency/range (55%), nose agreement
(25%), and the median range curve (20%). Distances use 0.2 MHz and 25 km
scales and are clipped to keep an occasional missed homing return from
dominating. The range curve uses every frequency with a return; display bins
are 0.1 MHz by 1 km. The frequency sweep is 2–10 MHz and the homing gate is
1,000 m.

## Run

Start with 40 private baseline NPZ files from
`run_cartman_parallel_ionograms.py`. Use an observed NPZ containing returns
and D settings; its generating parameters are not needed by the fit. For a
synthetic validation, remove the parameter fields from the observed copy
before initializing.

```sh
python3 reports/fit_d_ionogram.py init OBSERVED.npz LIBRARY_DIR FIT_DIR
python3 reports/run_cartman_inverse_population.py submit FIT_DIR 1 RUN_NAME_1
python3 reports/run_cartman_inverse_population.py status RUN_NAME_1
rsync -a cartman:/homes/chartat1/private_raytracing/runs/RUN_NAME_1/results/ FIT_DIR/results_round_1/
python3 reports/fit_d_ionogram.py advance FIT_DIR FIT_DIR/results_round_1
```

Repeat submit, status, copy, and advance for rounds 2 and 3. The Cartman
runner creates owner-only job, log, scratch, cache, and result paths in
`/homes/chartat1/private_raytracing/runs`. Its status command checks
completion and privacy. `plot_d_inverse_fit.py FIT_DIR/state.json PREFIX`
then writes ionogram and convergence figures.

The score is a fit to the chosen forward model. It does not establish that
the model captures profile shapes outside its density/height family; a
separately perturbed profile is needed to measure that model error.

## Held-out synthetic validation, 2026-09-25

The generating parameters were removed from the observed NPZ before the
search. This tests parameter recovery from an unseen ionogram, while still
using the same D model family and Cartman PyLap build on both sides.

| Quantity | Synthetic truth | Retrieved | Error |
| --- | ---: | ---: | ---: |
| Density scale | 1.137000 | 1.136385 | −0.000615 |
| F2-height shift (km) | +12.000 | +12.0466 | +0.0466 |

The best score fell from 0.24857 in the 40-member density library to 0.05773,
0.01558, and 0.00265 after the three refinement rounds. All 120 refinement
jobs finished, with no server ownership or permission violations. The three
40-job elapsed times were 248, 626, and 315 seconds; shifted profiles had a
wide runtime spread.

The [truth/retrieved ionograms](figures/d_inverse_fit_ionograms.png) show every
accepted return in 0.1 MHz by 1 km bins. The
[convergence figure](figures/d_inverse_fit_convergence.png) and
[score table](data/d_inverse_fit_scores.csv) show all 160 evaluations. The
[validation summary](data/d_inverse_fit_validation.json) and two NPZ files
hold the exact inputs and best prediction. This is an in-family synthetic
recovery result, not an independent physical truth test.
