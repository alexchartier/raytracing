---
title: "Retrieving a wave-rich SAMI3/HIAMCM ionosphere"
subtitle: "Vertical sounder and 600 km two-satellite oblique pilot"
date: "29 September 2026"
geometry: margin=0.72in
fontsize: 10pt
colorlinks: true
---

## Summary

We tested the retrieval on a frozen **11 January 2017, 06:00 UTC** SAMI3 density snapshot at **190°E**, from **70°S to 50°S**. Seven positions sample a natural F2 trough and crest. The truth density has local foF2 from roughly **4.44 to 5.40 MHz** along the full 20-position cut. An independently generated, date-matched PyIRI density field starts the inversion. Both simulations use the same ray tracer, 800 km spacecraft altitude, separate O and X modes, and frequencies **2–10 MHz in 0.1 MHz steps**.

The vertical inversion improves the **F2 peak density** at the seven sounder positions, but the ionogram group ranges and lower density profile still differ markedly from truth. The oblique test starts from that vertical solution, so its result measures the additional value of the 600 km links rather than a separate inversion from PyIRI.

| Seven-position evaluation | PyIRI vertical prior | Vertical selected | Oblique starting field at link midpoints | Oblique selected at link midpoints |
|:--|--:|--:|--:|--:|
| Mean O/X ionogram score, lower is better | 0.7348 | **0.5837** | 0.6374 | **0.5886** |
| Mean absolute foF2 error | 0.190 MHz | **0.051 MHz** | 0.041 MHz | **0.035 MHz** |
| Mean absolute F2 peak-density error | 7.86% | **2.07%** | 1.60% | **1.36%** |
| Mean absolute hmF2 error | 19.9 km | 20.2 km | 23.1 km | 10.8 km |
| Density RMS, 150–600 km, normalized by truth-section peak | 17.9% | **16.2%** | 17.34% | **12.11%** |

The vertical sounders and oblique link midpoints occupy different latitudes, so their density-error columns are not a direct ranking of the two geometries. FoF2 and hmF2 above are calculated from the **density profiles**, not from the final ionogram return frequency. The ionogram score is a dimensionless discrepancy, not a percentage accuracy.

\clearpage

## 1. Truth, prior, and what is independent

The saved MATLAB snapshot `2017-01-11_0600.mat` contains the `dene` electron-density field and a prepared ray-tracing density grid. Its time value converts to **11 January 2017, 06:00 UTC**. We checked that those densities agree over the regional grid after longitude reordering; the source SHA-256 is recorded in the case manifest. The repository's `remap_sami.m` identifies this snapshot path as the HIAMCM version of the SAMI3 output. Its original source NetCDF is unavailable here, so that upstream label cannot be independently rechecked from the original file. We interpolated density onto a **5 km altitude grid from 80 to 820 km** and traced rays through a regional latitude–longitude grid. The latitudinal wave and trough are present in the snapshot; they were not imposed for this experiment.

PyIRI supplied the independent **starting electron density** at the same epoch, with F10.7 = 75 and daily Ap = 5. The truth grid retains PyIRI's magnetic, neutral, and thermal supporting fields while replacing its electron density with SAMI3. Truth and retrieval ionograms therefore share the PyLap propagation and homing physics; this is an independent **density-model** comparison, not independent ray-physics or instrument validation. The wave location was selected after inspecting the SAMI3 density, so this is an exploratory case rather than a blind random sample.

The seven tested transmitter latitudes correspond to positions **1, 5, 9, 12, 14, 17, and 20** of a 20-position pass. Every plotted ionogram displays **all accepted reflected returns** at the exact traced 100 kHz frequencies, with group ranges rounded to **1 km** for display. Blue is O mode and red is X mode. A receiver miss of at most **1 km** is accepted; homing gaps can therefore affect the last-return frequency, especially near the nose.

## 2. Vertical inversion

At each of the seven positions we measured the last accepted O and X return frequencies in the SAMI3 and PyIRI ionograms. Their squared frequency ratios proposed a local peak-density multiplier. A shape-preserving cubic interpolation made one smooth latitude field from those seven multipliers. The median low-frequency O-mode group-range residual proposed a topside stretch about each prior F2 peak. We forward-traced the prior and a peak-plus-topside candidate with the full three-dimensional O/X ray generator. A second peak correction, confined by a **60 km Gaussian** around the current local peak, used the remaining O/X nose residuals and received another full forward check. Candidate selection used only these ionograms.

The ray generator used the option-D equal-area guarded vertical fan. It recovered empty frequency/mode bins with extra seeds and used 20 kHz intermediate frequencies near the nose for continuation; each **saved** return was traced again at its displayed 100 kHz grid frequency. The score averaged a symmetric nearest-return distance, O/X last-return-frequency error, and common-frequency median group-range error with weights **0.55, 0.25, and 0.20**. The full-set minimum selected `peak_round2`: score **0.5837**, versus **0.7348** for PyIRI. The fitted first-stage peak multipliers span **0.857–1.060**; the selected full topside proposal spans **1.315–2.299**. The second local peak factors span **0.923–1.031**.

\clearpage

### Vertical O/X ionograms

![SAMI3 truth (left) and selected vertical retrieval (right) for positions 9, 14, and 20. Every accepted O return is blue and every accepted X return is red. Frequencies are the traced 0.1 MHz grid, and displayed group ranges are rounded to 1 km. The large low-frequency range mismatch at positions 14 and 20 is visible despite similar noses.](figures/sami3_wave_vertical_ionograms.png){width=96%}

The retrieval follows the **frequency cutoff** better than the **range ridge**. At sounder 14, truth has accepted low-frequency O and X returns at substantially longer group ranges than the selected model. At sounder 20 the modes overlap well around 4–5 MHz, then separate in range toward lower frequencies. These differences are not a plotting artifact: both panels display every saved accepted return under the same homing gate.

\clearpage

### Vertical density and F2 peak

![SAMI3 truth, vertically retrieved electron density, and retrieved-minus-truth residual at the seven sounder latitudes. Truth and retrieved panels share one scale. The retrieval is deficient below about 250 km and puts its peak too high over much of the cut.](figures/sami3_wave_vertical_density.png){width=100%}

![Left: truth and retrieved altitude profiles at vertical sounder 14. Right: density-derived local foF2 and hmF2 at all seven vertical sounders. These are profile properties evaluated against independent SAMI3 density, not fitted ionogram noses.](figures/sami3_wave_vertical_peaks.png){width=96%}

The selected vertical field has **0.051 MHz foF2 MAE** and **2.07% peak-density MAE**, but its **20.2 km hmF2 MAE** is essentially no better than PyIRI's 19.9 km. Its 150–600 km normalized density RMS remains **16.2%**. The paired ionograms reveal why the peak result alone is insufficient: the truth group range at lower frequencies can be over 100 km longer than the modeled return.

\clearpage

## 3. Two-satellite oblique inversion

For each tested position, the transmitter and receiver are both at **800 km**, with a **600 km straight-line separation** along the latitude track. We traced O and X modes from the same 2–10 MHz grid. The starting oblique density is exactly the selected vertical grid above. The saved oblique observables keep reflected rays with group range at least **700 km** whose paths descend at least **100 km** below spacecraft altitude. This filter removes the roughly 600 km quasi-direct path; it cannot appear in either plotted oblique ionogram.

The seven pairs of truth and starting-field oblique ionograms supplied two proposals. O/X last-return-frequency residuals suggested a local peak correction through $\Delta N/N\simeq2\Delta f/f$. Median group-range residuals over **2.5–4.5 MHz** suggested a reflector-height shift through $\Delta h\simeq-\Delta R/1.7$. A regularized constant, latitude trend, sine, and cosine basis at a **fixed 1300 km** wavelength smoothed the peak proposal; a cubic spline smoothed the height proposal. We tested a height-only candidate and one with the same height shift plus peak correction. Their allowed peak change was at most about **3%**, and their downward altitude shift ranged **8–20 km**. These are conservative steps from the vertical solution: the raw range residuals imply a much larger displacement, outside the pre-set **25 km** height trust region.

Each candidate was forward-traced on all seven links. Selection minimized **0.8 times the O/X ionogram score plus 0.2 times a Doppler score**. We assumed spacecraft speed **8 km/s** and rounded simulated Doppler to **1 Hz**. An observed return is paired to the nearest modeled return of the same mode and frequency when their ranges differ by at most **40 km**; its rounded Doppler difference is scaled by **20 Hz** and capped at one. Unmatched rays cost one. This score can saturate when paths or group ranges differ substantially, so its value is reported separately. The selected candidate is **`gain1_sigma60_height0.15_spline`**.

\clearpage

### Oblique O/X ionograms

![SAMI3 truth (left) and selected oblique retrieval (right) for links 9, 14, and 20. Blue and red mark O and X accepted reflected returns at the traced 0.1 MHz frequencies, with 1 km display bins. The absent quasi-direct branch was excluded before saving.](figures/sami3_wave_oblique_ionograms.png){width=96%}

The downward height correction moves the modeled oblique ridge closer to truth, but links 9 and 20 still have low-frequency range deficits and incomplete near-nose agreement. These are the reflected branches used in the inversion. A quasi-direct trace would begin near 600 km and would require a separate forward simulation and score; adding a line to this plot would not recover an excluded ray.

\clearpage

### Oblique density and F2 peak

![SAMI3 truth, oblique-selected electron density, and retrieved-minus-truth residual at the seven link midpoints. The first two panels share one density scale.](figures/sami3_wave_oblique_density.png){width=100%}

![Left: truth and retrieved altitude profiles at the midpoint of oblique link 14. Right: density-derived local foF2 and hmF2 at all seven link midpoints. Truth and retrieved values use the same midpoint coordinates.](figures/sami3_wave_oblique_peaks.png){width=96%}

The selected oblique candidate changes the seven-link mean O/X score from **0.6374** to **0.5886** and its 1 Hz Doppler score from **0.9657** to **0.9620**. Doppler pairing remains almost fully penalized, so that number is weak evidence about angle recovery. The combined score is **0.6633**, only **0.0018** below the height-only candidate's score. In the independent density evaluation, midpoint RMS falls from **17.34%** to **12.11%**, foF2 MAE from **0.041** to **0.035 MHz**, peak-density MAE from **1.60%** to **1.36%**, and hmF2 MAE from **23.1** to **10.8 km**. The remaining hmF2 bias is **+10.7 km**, and the lower profile and parts of the ionogram ridge remain visibly wrong.

\clearpage

## 4. Interpretation and limits

The height-only candidate was almost tied with the selected height-plus-peak candidate on the observable score (**0.6650** versus **0.6633**). Its independent midpoint density RMS was **12.15%**, foF2 MAE **0.042 MHz**, and hmF2 MAE **10.8 km**. Adding the small peak correction improved foF2 and peak-density accuracy, while the height result barely changed. The apparent advantage is modest compared with the large remaining ionogram ridge mismatch.

The SAMI3 bottomside structure and latitude-dependent peak height challenge this low-dimensional correction. Vertical sounding recovers F2 peak frequency much better than the full altitude profile. Oblique links add a second viewing geometry, but inherit the vertical field and sample only seven midpoints. Large range residuals, receiver-homing gaps, and saturated Doppler pairings limit the current score's ability to constrain the lower profile.

Candidate generation and the stated selection scores use saved truth **ionograms**, modeled ionograms, and the candidate grids. The SAMI3 truth **density** is used for the displayed errors and contour cuts. A baseline diagnostic render loaded that truth grid before the final oblique score selection, but it was neither inspected nor used to change the two predeclared oblique candidates or the selection rule. The earlier choice of this wave-rich location from the truth density and the common ray physics remain the principal limits on claims of independent validation.

### Reproducibility record

The [case manifest](data/sami3_wave_20170111_0600/manifest.json) records the epoch, positions, sweep, homing gate, and source checksum. The [vertical selection](data/sami3_wave_20170111_0600/vertical_selection.json) and [density evaluation](data/sami3_wave_20170111_0600/vertical_evaluation.json), plus the [oblique selection](data/sami3_wave_20170111_0600/oblique_wave_600km/pilot_selection.json) and [selected density evaluation](data/sami3_wave_20170111_0600/oblique_evaluation_gain1_sigma60_height0.15_spline.json), hold the exact reported values. `plot_sami3_wave_results.py` regenerates the paired ionograms, density sections, altitude cuts, and peak curves from saved outputs. Full-grid rays ran in `chartat1`'s private Cartman area; the 16 GB laptop was used only for lightweight scoring, plotting, and PDF generation.
