# Spacecraft Doppler on the latitude-wave topside pass

## Setup

The 20-profile pass uses an assumed 8 km/s northbound spacecraft and a
temporally frozen ionosphere. For each accepted O and X return, the carrier
Doppler from spacecraft motion is calculated from PyLap's launch and arrival
phase momenta:

\[
f_D = \frac{f}{c}\,\mathbf{v}\cdot
      (\mathbf{p}_{\rm launch}-\mathbf{p}_{\rm arrival}).
\]

The arrival momentum is interpolated at the receiver. The forward runs use the
same option-D equal-area fan, 1 km homing gate, 2–10 MHz sweep at 0.1 MHz, and
empty-bin recovery as the earlier pass. The example [truth and retrieved
Doppler curves](figures/lat_wave_doppler_examples.png) show the modeled
frequency dependence at two positions; O is blue and X is red.

## Does it add a measurable constraint?

The [Doppler diagnostics](figures/lat_wave_doppler_diagnostics.png) show a
3.41 Hz median absolute Doppler in the synthetic truth, with 10.49 Hz at the
90th percentile. Of 2,053 accepted truth rays, 1,449 remain after grouping
near-duplicate returns within 1 km and 0.1 Hz in the same frequency and mode
bin. All accepted rays remain in the saved ionograms. Among 527 distinct
truth/retrieved paths matched within 10 km in group range, 412 (78.2%) differ
by more than the assumed **1 Hz** Doppler precision. The residual correlation
between group range and Doppler is 0.12 for matches within 30 km. These
numbers show useful angular information even where group range is similar;
they do not establish that the present density retrieval improves when it is
used. For scale, a 1 Hz error at 5 MHz corresponds to about 0.21° in the
effective two-way along-track angle near nadir at 8 km/s.

## Candidate-selection experiment

Four nearby density candidates were fixed before scoring: global density scale
±0.02 and wave-amplitude fraction ±0.03 around the existing ionogram-only
fit. Each was traced at all 20 positions and recovered in empty frequency
bins. The conventional O/X ionogram score was combined with a Doppler cost.
The Doppler cost pairs distinct returns at the same frequency and mode by
group range (30 km gate), clips residuals at three assumed 1 Hz standard
deviations, and assigns the maximum penalty to unmatched returns. The primary
joint score is ionogram score + 0.2 × Doppler cost.

| Candidate | Ionogram score | Doppler cost | Matched distinct returns |
| --- | ---: | ---: | ---: |
| Existing fit | 0.2344 | 0.7054 | 1,032 |
| Density scale −0.02 | 0.2425 | 0.7015 | 984 |
| Density scale +0.02 | 0.2431 | 0.7146 | 1,050 |
| Wave amplitude −0.03 | 0.2435 | 0.7215 | 974 |
| Wave amplitude +0.03 | **0.2236** | **0.6940** | **1,086** |

The higher-amplitude candidate wins with Doppler weights 0, 0.1, 0.2, and 0.4.
**Doppler did not change the chosen candidate in this experiment.** It agrees
with the ionogram preference. The [selected candidate ionograms](figures/lat_wave_doppler_retrieval_ionograms.png)
show O (blue) and X (red) together in paired truth and retrieval panels, with
every accepted return displayed. Only after the choice was fixed was the
withheld density opened. The [peak-density comparison](figures/lat_wave_doppler_retrieval_peaks.png)
shows mean absolute peak-density error of 3.06% for the existing fit and 2.99%
for the chosen candidate; mean absolute foF2 error changes from 0.084 to
0.081 MHz. The small improvement cannot be attributed to Doppler because the
ionogram score alone chose the same candidate.

## Scope and provenance

The withheld pass density comes from Fortran IRI-2016 with an imposed wave;
the candidate density is a fitted PyIRI-based grid. Candidate selection reads
only O/X group-range ionograms and synthetic spacecraft Doppler. Both use the
same PyLap forward physics. The PyIRI background was calibrated against one
earlier IRI ionogram at the pass center, so this is a withheld pass/wave test,
not fully independent cross-validation. The synthetic ionosphere is frozen;
plasma motion, instrument noise, Doppler extraction error, and uncertainty
from the 1 km homing gate were not simulated. A single velocity projection
cannot uniquely determine both endpoint angles or wave bearing.

The [selection record](data/lat_wave_doppler_selection.json),
[Doppler assessment](data/lat_wave_doppler_assessment.json), and
[withheld-density evaluation](data/lat_wave_doppler_retrieval_evaluation.json)
contain the underlying numbers. All Cartman jobs completed and the remote
owner-only permission checks passed.
