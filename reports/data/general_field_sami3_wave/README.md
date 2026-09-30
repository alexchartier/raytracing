# General-field retrieval on the SAMI3/HIAMCM pass

This run tests a PyIRI-started density-field retrieval on seven archived
SAMI3/HIAMCM vertical ionograms (indices 1, 5, 9, 12, 14, 17, 20). The
observations and candidates use the same raw adaptive O/X generator, 2–10 MHz
at 100 kHz spacing, option D fan, and 1 km homing tolerance. Scoring uses all
accepted O/X returns, ridge ranges, and O/X noses. Doppler gating is disabled
because its near-vertical subset truncates the observed nose at several
profiles. The saved Doppler values remain available for future fitting.

The density prior is `../sami3_wave_20170111_0600/pyiri_prior_grid.nc`, built
independently of SAMI3. `prepare_general_field_sami.py` opens the SAMI truth
once to make **seven scalar 800 km local-density observations**. The candidate
builder, nose initialization, score, and selection do not open the SAMI truth
density. The full truth grid is opened only by `evaluate_general_field_sami.py`
after selection is frozen.

The SAMI truth grid carries SAMI electron density remapped onto the same
latitude/longitude/altitude axes as PyIRI. The setup code retains the PyIRI
magnetic and collision background in that ray grid, so the independent part
of this test is the **electron-density field**, not every forward-model input.
The original source NetCDF is unavailable locally; the saved MATLAB remap has
its SHA-256 recorded in the case manifest.

The four paired pilot directions alter peak density, peak height, topside
width, and an along-track peak wave inferred from observed O/X cutoffs. Nine
spatial centers cover the 20° latitude track. The nose-initialized candidate
fits relative O/X cutoff shifts against the independent PyIRI forward
ionograms with the same log-density basis. Every selected candidate must pass
the full-ray forward generator.

The location was previously chosen for visible SAMI wave structure, and the
retrieval method was developed after examining related truth diagnostics.
This is an independent **model** test, not a blind validation. Only seven
positions have archived truth ionograms; the other 13 positions in the density
cut assess interpolation and have no direct sounding constraints.

## Result

All 45 new Cartman traces completed with exit status zero. The private run was
checked for `chartat1` ownership and modes 0700 (directories) and 0600
(files). Peak resident memory among these jobs was 24,374,900 kB on a Cartman
compute node; no PyLap rays ran on the laptop. Ionogram-only selection chose
`nose_init` before loading the withheld SAMI density. The seven-position mean
score was 0.694, compared with 0.735 for
untouched PyIRI, 0.724 for `minus_peak`, and 0.707 for `combined`. The selected
candidate improved the score at six of seven positions. The paired-position
bootstrap 95% interval for its score change versus PyIRI was -0.079 to -0.005.
These are same-data fitting scores, not held-out ionogram scores.

The private Cartman run is
`/homes/chartat1/private_raytracing/runs/general_field_sami3_wave` (job IDs
3396418, 3396419, 3396437). The 45 traces consist of 27 three-position
pilots and 18 additional traces for complete seven-position candidates.

| Density metric | PyIRI prior | Selected fit |
| --- | ---: | ---: |
| Peak-density MAE, seven observed positions | 7.84% | 2.23% |
| foF2 MAE, seven observed positions | 0.190 MHz | 0.054 MHz |
| hmF2 MAE, seven observed positions | 18.57 km | 18.57 km |
| Density NRMSE, 150–800 km, seven positions | 17.97% | 17.36% |
| Topside NRMSE, 350–800 km, seven positions | 12.57% | 11.77% |

The 20-position density cut gives a 2.30% peak-density MAE and a 19.25 km
hmF2 MAE for the selected fit. After removal of a linear latitude background,
the along-track foF2 residual correlates with truth at 0.988, versus 0.570
for PyIRI. Its RMS amplitude is 74% of the truth residual amplitude, so the
wave is still too weak. The selected grid leaves every peak 5–35 km too high
and underestimates much of the upper topside. At profile 14, its 3–5 MHz
O-mode ridge lies about 186 km below the truth ridge on average. The peak
improvement therefore does **not** establish an accurate full-profile
retrieval. The figures are `density_cut.png`, `density_profiles.png`, and
`ionogram_14.png`; exact per-position metrics are in `evaluation.json`.
