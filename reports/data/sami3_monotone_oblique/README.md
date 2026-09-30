# SAMI3/HIAMCM two-satellite oblique test

This test uses seven archived O/X reflected-return ionograms from the SAMI3/HIAMCM
wave case. The transmitter and receiver are both at 800 km, 600 km apart along
track. The frequency sweep is 2–10 MHz in 100 kHz steps, with a 1 km homing
gate and the same 75-direction oblique fan for truth and candidates. The
quasi-direct path is excluded by the generator, so it is absent from the plots
and scores.

The starting `varying_g04` field is the frozen result of the preceding
**vertical** retrieval. Its F2 height and two monotone-topside shape parameters
vary independently along latitude. `extended_anchor` rebuilds this same shape
from the PyIRI-derived `nose_init` base, replacing the 800 km anchor with
synthetic in-situ samples every 1 km over both spacecraft tracks. There are
2,759 samples from −70° to −45.204°. They were sampled from SAMI density as
the measurement assumed in this simulation; interpolation of the SAMI grid
does not create new spatial resolution. No oblique ionogram residual was used
to adjust the shape parameters in this test. The two new fields therefore
measure transfer to oblique sounding and the effect of extending the in-situ
constraint, rather than an oblique-specific shape optimization.

The full 3-D SAMI density field was not used to choose between candidates.
`selection.json` freezes the choice from the seven oblique ionogram scores and
the 800 km track samples, with hashes of the scores, plan, and selected grid.
The selection requires every in-situ sample within 0.5%, then minimizes the
uncapped O/X ionogram score. The 1 Hz Doppler score is reported separately and
was not used in selection. The truth electron density comes from SAMI3/HIAMCM,
independent of the PyIRI density prior. The magnetic and collision fields are
shared PyIRI inputs. Previous inspection of this SAMI case influenced method
development, so this is not blind validation.

| Field | Ionogram score, equivalent km ↓ | 1 Hz Doppler mismatch ↓ | Maximum full-track in-situ error |
| --- | ---: | ---: | ---: |
| PyIRI baseline | 138.9 | 0.966 | 78.80% |
| Earlier oblique fit | 113.3 | 0.962 | 54.73% |
| Frozen vertical field, transferred directly | **72.3** | 0.754 | 1.62% |
| Full-track in-situ anchor, selected | 74.2 | **0.750** | **0.283%** |

The score compares all accepted O/X group ranges at each occupied frequency,
charges 150 km for a frequency occupied on only one side, and adds 100 km per
MHz of mean O/X nose error. The Doppler number averages 1 Hz-quantized errors
against the nearest matched return, capped at 20 Hz; zero is ideal. The
selected field improves the ionogram score by 35% relative to the earlier
oblique fit, while satisfying the assumed in-situ tolerance. It does not
match every return: link 14 retains a clear low-frequency branch-shape error
in `ionogram_14.png`, and link 12 scores slightly worse than the earlier fit.

After selection was frozen, the full SAMI electron-density field was opened
for evaluation at the seven observed link midpoints and all twenty pass
midpoints. Values below are for the seven observed links.

| Density metric at link midpoints | Earlier oblique fit | Selected field |
| --- | ---: | ---: |
| Peak-density MAE | **1.37%** | 2.20% |
| foF2 MAE | **0.035 MHz** | 0.055 MHz |
| hmF2 MAE | 10.0 km | **7.86 km** |
| Density NRMSE, 150–800 km | 11.47% | **8.27%** |
| Topside NRMSE, 350–800 km | 5.82% | **2.99%** |

Across all twenty link midpoints, the selected field has 2.41% peak-density
MAE, 0.058 MHz foF2 MAE, 6.0 km hmF2 MAE, and 2.83% topside NRMSE. The
earlier oblique fit has better peak density, but the selected field has
better height and topside shape. The selected forward grid and all twenty
sampled midpoint columns have zero upward steps above their F2 peaks.

`ionogram_14.png` puts truth and retrieved O/X returns in side-by-side panels,
with a 0.1 MHz × 1 km plotting grid. `density_cut.png` shows the truth and
retrieved latitude–altitude cuts at link midpoints and their difference.
`density_profiles.png` shows midpoint profiles for links 9 and 14.
`evaluation.json` contains per-link peak and height errors, the full-pass
comparison, and the in-situ and monotonicity checks.

All fourteen new oblique traces ran privately on Cartman in jobs 3396446–
3396449 with at most two large jobs active. Every trace exited zero. The
largest measured resident set was 19,920,496 kB; no SAMI rays ran on the
16 GB laptop. The private run had owner `chartat1`, mode 0700 directories,
and mode 0600 files after completion.
