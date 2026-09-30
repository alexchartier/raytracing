# Spatially varying monotone-topside SAMI3 retrieval

This run fits seven archived vertical O/X ionograms from a SAMI3/HIAMCM
electron-density snapshot. The candidate generator starts from the earlier
`nose_init` grid. At 800 km, it uses 2,225 synthetic in-situ density samples
spaced by at most 1 km along the 20° latitude pass. The samples were extracted
once from the saved SAMI grid by `prepare_sami_dense_insitu.py`; candidate
construction, ionogram scoring, and selection use these one-dimensional
samples and the archived truth **ionograms**, not the full SAMI density field.
The 1 km sampling here comes from interpolation of the saved model grid, so it
does not add spatial detail absent from that grid.

At each of the seven sounding latitudes, a candidate has its own F2 peak-height
shift and two topside log-density fractions. The fractions specify log density at
25% and 60% of the distance between the shifted peak and the 800 km
spacecraft. Shape parameters are interpolated along latitude using PCHIP. A
second PCHIP builds each vertical log-density profile through the F2 peak,
the two shape points, and the measured spacecraft density. Its strictly
descending knot values enforce a monotone topside; every grid column is checked
for upward steps above its peak. Every 1 km in-situ sample is checked against
the resulting forward-grid interpolation within a 0.5% tolerance.

The initial shape at each latitude is inferred from the previous model's
3–5 MHz O/X ridge residual. Four amplitudes and a no-height-shift control
were traced on pilot soundings 9, 14, and 17. Each trace uses the same option-D
fan, 2–10 MHz at 100 kHz, and 1 km homing gate as truth. The score uses every
accepted O/X return: symmetric nearest group-range distance at each occupied
frequency, 150 km for a frequency occupied on only one side, plus 100 km per
MHz of mean O/X nose error. Group-range errors are not clipped.

| Pilot candidate | Three-sounding score (equivalent km) |
| --- | ---: |
| Earlier `nose_init` | 196.4 |
| Varying shape, gain 0.4 | **89.0** |
| Varying shape, gain 0.7 | 99.3 |
| Varying shape, gain 1.0 | 135.2 |
| Varying shape, gain 1.3 | 190.3 |
| Gain 1.0, no peak-height shift | 143.2 |

At profile 14, the previous 3–5 MHz O-mode ridge averaged 186 km below the
observed ridge. The gain-0.4 candidate averages 28 km above it. Its pilot
ionogram score is 84.8 versus 160.3 for the previous field. This is a fit to
the same ionograms used to tune the candidate; it is not a held-out prediction.

The truth's electron-density field comes from SAMI3/HIAMCM rather than PyIRI.
The saved ray grid keeps the PyIRI magnetic and collision backgrounds, so only
the electron density is independent. The SAMI location was previously chosen
for a visible wave; this is an independent-model test, not blind validation.
The full density field was opened only after the seven-ionogram selection was
frozen in `selection.json`, which records SHA-256 hashes of the scores, plan,
and selected grid. Earlier examination of truth-density diagnostics motivated
the method and the peak-height adjustment, so this is not a blind test.

## Seven-sounding result

Both pilot finalists were traced at all seven observed positions. The mean
uncapped ionogram score was 178.2 equivalent km for the previous `nose_init`
field, 90.2 for gain 0.7, and **78.7 for gain 0.4**. Gain 0.4 was selected
before inspecting the full truth density. All seven gain-0.4 ionograms scored
better than the earlier field. At profile 14, its ionogram score improved from
160.3 to 84.8 equivalent km; its 3–5 MHz O-mode ridge went from a mean 186 km
below truth to 28 km above truth. The retrieved O/X branches still disagree
with the observed low-frequency shape and upper ridge near the nose.

| Density metric at seven observed positions | Earlier field | Selected varying topside |
| --- | ---: | ---: |
| Peak-density MAE | 2.23% | 2.22% |
| foF2 MAE | 0.054 MHz | 0.054 MHz |
| hmF2 MAE | 18.57 km | 4.29 km |
| Density NRMSE, 150–800 km | 17.36% | 8.42% |
| Topside NRMSE, 350–800 km | 11.77% | 2.96% |

Across all twenty latitude positions in the cut, the selected field has 2.25%
peak-density MAE, 4.75 km hmF2 MAE, and 2.77% topside NRMSE. Its maximum
error against the 2,225 synthetic in-situ samples is 0.0157%, and it has no
upward density steps above the peak in any grid column or sampled track
profile. At profile 14 specifically, topside NRMSE decreased from 16.32% to
5.57%. The selected profile remains low through part of the middle topside,
and its bottomside below about 230 km is too thin; the topside sounder cannot
directly observe that lower region. Peak density changed little because this
iteration mainly changed height and topside curvature, with the earlier peak
amplitude retained.

The figures are `density_cut.png` (truth, retrieved, and difference),
`density_profiles.png` (profiles 9 and 14), `ionogram_14.png` (accepted O/X
returns in truth and retrieval), and `topside_parameters.png` (the seven
latitude-varying shape and peak-height parameters). `evaluation.json` gives
exact errors at the seven soundings and all twenty positions.

The full SAMI ray grid exceeded this laptop's memory. All new rays ran in
private Cartman jobs as `chartat1`, with at most two active at a time. The
private run is `/homes/chartat1/private_raytracing/runs/sami3_monotone_topside`
(job arrays 3396442, 3396443, and 3396444). All 23 new traces finished with
zero exit status and private owner-only permissions. The largest measured
resident set was 32,969,948 kB on Cartman; no rays ran on the 16 GB laptop.
