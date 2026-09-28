# Latitude-wave truth and retrieval rerun with 20 kHz continuation

## Result

The 20-profile pass was rerun from the saved, 100 kHz O/X ionograms. Near the
measured nose, full-fan solutions were seeded at the nose and 0.5 and 1.0 MHz
below it, then homed in 20 kHz increments up and down by as much as 0.54 MHz.
The final pass also checked full-fan seeds 0.2, 0.4, and 0.6 MHz above each
mode's nose. Every saved point was raytraced at its plotted 100 kHz frequency;
the 20 kHz intermediate solutions were used only to continue paths. All
accepted returns meet the 1 km homing gate and 150 km minimum group range.

The truth pass gained **365 accepted solutions** relative to its original
2,053 returns. These filled **122 previously empty frequency/mode bins** and
raised the last accepted frequency for 19 mode/profile combinations. The
above-nose probe added no further truth returns. Multiple accepted solutions
can still describe nearly the same physical path; the counts are not counts
of unique propagation branches.

The new joint wave fit inferred a **915.2 km** wavelength and **14.43%**
fractional modulation at the local F2 peak. The original fit inferred
916.7 km and 14.84%; the imposed truth wave has a 900 km latitude wavelength
and 18% modulation at 280 km, which is a different amplitude definition.
Both candidates were forward traced with the same option D fan, O/X separation,
empty-bin checks, continuation, and above-nose probes. The ionogram-only
comparison selected the **previous fit**: mean score **0.231** versus **0.252**
for the new refit. [The saved selection](data/lat_wave_final_candidate_selection.json)
contains all 20 profile scores.

For the selected retrieval, peak-density mean absolute error is **3.06%** and
foF2 mean absolute error is **0.084 MHz** against the withheld density. The
150–600 km mean profile normalized RMS error is **6.98%**. The new refit was
worse on this same check: **3.56%** peak-density error and **0.098 MHz** foF2
error. Both density errors were calculated after ionogram-based candidate
selection. The selected retrieval has 2,310 accepted returns versus 2,418 in
truth. Its mean absolute O/X nose errors are **0.19/0.19 MHz**. Common-frequency
median group-range errors are **25.0/21.8 km** for O/X.

## Figures

[Profile 14 truth and selected retrieval, with a nose zoom](figures/lat_wave_final_profile14_ionograms.png)
shows every accepted return on a 0.1 MHz × 1 km display, O in blue and X in
red. The full panel starts at 150 km group range. The O nose now reaches
**5.2 MHz in both truth and retrieval**; the above-nose probe found a candidate
O branch that the initial continuation missed. The X nose is 5.6 MHz in truth
and 5.4 MHz in retrieval. The retrieved group ranges remain too short near
both noses, so equal O cutoffs do not imply matching profile shape. A vector
[PDF of the pair](figures/lat_wave_final_profile14_ionograms.pdf) is also saved.

The [four-position O/X ionogram comparison](figures/lat_wave_final_ionograms.png),
[all-pass cutoff and peak plot](figures/lat_wave_final_peaks.png), and
[latitude-altitude density comparison](figures/lat_wave_final_density.png)
show the pass-wide result. The figures place truth and retrieval side by side;
they do not overlay a third curve on the ionograms.

## Provenance and limits

Truth electron density came from Fortran IRI-2016 with an imposed wave. The
retrieved density is a transformed PyIRI background. Both ionogram sets use
the same PyLap ray physics, and the PyIRI prior was calibrated on an earlier
IRI ionogram at pass center. The truth pass density and imposed wave parameters
were not read by the fitting or candidate-selection code. The previous fit,
however, had already been evaluated against this same truth density in earlier
work. This rerun is therefore a numerical coverage and retrieval check, not a
fully blind independent validation case.

All 20 jobs completed in each Cartman batch, and every remote batch passed
owner-only permission checks. The truth continuation took **259 s**, and its
above-nose check took **124 s**. The selected candidate's continuation took
**307 s**, and its above-nose check took **118 s**. The saved
[remote-run record](data/lat_wave_continuation_remote_runs.json),
[final evaluation](data/lat_wave_final_summary.json),
[truth ionograms](data/lat_wave_final_truth_ionograms), and
[selected retrieval ionograms](data/lat_wave_final_previous_fit_ionograms)
contain the underlying values.
