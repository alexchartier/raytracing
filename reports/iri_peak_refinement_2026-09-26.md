# IRI ionogram peak refinement

## Method

The previous 240-candidate retrieval estimated an F2 peak of 350,853 cm⁻³.
Its O-mode accepted-return limit was 5.2 MHz, while the observed limit was
5.4 MHz. The X-mode limit already matched at 5.4 MHz. Both modes are measured
and scored separately.

A local Gaussian multiplier was added to the candidate electron-density grid
at each F2 peak. Its amplitude was centered on the observed O-mode limit using
the plasma-frequency relation, then bracketed with eight amplitudes at 40 km
width. Eight more candidates varied the width between 20 and 30 km. All 16
new candidates were ray traced over 2–10 MHz at 0.1 MHz spacing with the same
option D fan and 1 km homing gate as the earlier retrieval.

Selection used the ionogram only: each O/X frequency limit had to match the
observed 0.1 MHz bin, each mode's accepted-return count had to be within 10%
of observed, and the count-aware return/range score broke ties. Six of the 256
total candidates passed these conditions. The selected correction is 3.0646%
at a 20 km Gaussian width. The IRI-2016 truth density was opened only after
selection. The 40 km trial alone matched both limits but worsened the overall
score, motivating the narrower second sweep.

## Results

| Quantity | Previous fit | Peak refinement | Independent truth |
| --- | ---: | ---: | ---: |
| F2 peak density (cm⁻³) | 350,853 | 360,675 | 362,539 |
| F2 peak error | −3.22% | −0.51% | — |
| O / X accepted-return limits (MHz) | 5.2 / 5.4 | 5.4 / 5.4 | 5.4 / 5.4 |
| O / X accepted returns | 54 / 58 | 54 / 59 | 56 / 60 |
| Count-aware ionogram score | 0.1353 | 0.1118 | — |
| Topside density RMS / truth peak | 1.76% | 1.68% | — |
| Bottomside density RMS / truth peak | 15.67% | 15.47% | — |

The [O/X ionograms](figures/d_inverse_iri_peak_ionograms.png) place every
accepted return in a 0.1 MHz by 1 km cell, with O blue and X red overlaid in
each of the two panels. The [density profiles](figures/d_inverse_iri_peak_density_profile.png)
show the selected peak and the remaining bottomside discrepancy. The full
[summary](data/d_inverse_iri_peak_summary.json) and
[candidate scores](data/d_inverse_iri_peak_evaluations.json) retain the exact
numbers.

## Independence and limits

The withheld truth density comes from Fortran IRI-2016; candidates use PyIRI
densities transformed by the retrieval parameters. The ionograms share the
same geometry, magnetic-field treatment, ray tracer, and homing algorithm.
This therefore checks density-model independence, not an independent forward
code or a real sounder observation. The near-peak error fell by about a factor
of six, while bottomside shape error remains about 15% of the truth peak.
The O-mode median range error at common frequencies rose from 10.54 to
12.90 km, although the overall count-aware ionogram score improved. Further
bottomside refinement should be evaluated against multiple independent
profiles before treating this as a general accuracy estimate.
