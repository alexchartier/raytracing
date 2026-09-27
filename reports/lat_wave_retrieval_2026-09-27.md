# Joint retrieval of the 20-profile latitude-wave pass

## Result

The 20 O/X ionograms between −63.68° and −43.68° latitude were fitted jointly.
The fitted latitude wavelength is **917 km**; the withheld wave has a **900 km**
wavelength. The fitted fractional density modulation is **14.8% at the local F2
peak**, with a phase of **0.305 rad** at pass center. The truth wave is specified
as **18% at 280 km**, so the two amplitude numbers have different altitude
definitions. The retrieval assumed a 65 km vertical Gaussian width, and a
north-south wave bearing; one fixed-longitude pass cannot determine bearing.

The [density comparison](figures/lat_wave_retrieval_density.png) shows the
independently generated IRI-2016 density, the retrieved PyIRI-based density,
and their difference. Across the 20 sounder positions, peak-density mean
absolute error is **3.06%**, and the associated foF2 mean absolute error is
**0.084 MHz**. The unperturbed PyIRI prior has **9.12%** and **0.246 MHz**
errors on the same pass. Mean density-profile normalized RMS error from
150–600 km falls from **8.14%** for the unperturbed prior to **6.98%** for the
retrieval. The topside (260–600 km) error is **2.18%**, while the 150–240 km
bottomside error is **14.31%**. The bottomside bias is visible in the difference
panel and remains a limitation of this topside sounder inversion.

The [side-by-side ionograms](figures/lat_wave_retrieval_ionograms.png) show four
representative positions. O returns are blue and X returns red, with every
accepted ray shown at 0.1 MHz frequency resolution and its group range rounded
to the nearest kilometer. The display starts at 150 km group range. The
[cutoff and peak plot](figures/lat_wave_retrieval_peaks.png) covers all 20
positions. Using the common-frequency median range, the mean absolute
ionogram difference is **19.95 km for O** and **19.52 km for X**. Cutoff mean
absolute errors are **0.265 MHz for O** and **0.220 MHz for X**. The truth has
2,053 accepted returns and the retrieved forward run has 1,911. These counts
and cutoff errors retain sensitivity to missed homing solutions; neither is
an exact measurement of density error.

## Method and independence

The search used only the 20 recovered ionograms and the earlier PyIRI profile
fit as its background. It fitted O and X noses separately, plus the latitude
variation in their median group ranges. A fast plane-stratified group-delay
approximation supplied the search objective; frequency-specific offsets
absorbed its systematic difference from the full magnetoionic ray tracer.
After selection, the candidate density was built from the saved PyIRI
background and new PyIRI physical fields. The 20 candidate ionograms were
traced using the same option D equal-area fan, separate O/X modes, 1 km homing
gate, 2–10 MHz sweep at 0.1 MHz steps, and full-fan empty-bin checks used for
the truth ionograms.

The withheld density field was opened **only after** the candidate was fixed
and its full forward run finished. Its electron density came from Fortran
IRI-2016; candidate density came from PyIRI. Both sets share PyLap forward
physics and similarly prepared magnetic/collision fields. The prior PyIRI
profile was calibrated on an earlier IRI ionogram at the pass center, so this
pass is independent in density generation and withheld wave parameters, but
not a completely untouched cross-validation case. The fixed-longitude pass
also leaves wave bearing unidentifiable.

The full 20-ionogram candidate run completed in **222.7 s** on Cartman; the
empty-bin recovery completed in **89.1 s**. All 40 jobs succeeded and the
remote run passed owner-only permission checks. See the
[initial fit](data/lat_wave_pass_initial_fit.json),
[evaluation summary](data/lat_wave_pass_retrieval_summary.json),
[retrieved density](data/lat_wave_pass_retrieved_density.npz), and
[retrieved ionograms](data/lat_wave_pass_retrieved_ionograms).
