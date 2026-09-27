# Latitude-wave sounder pass

## Setup

Twenty vertical sounder positions run from −63.6786° to −43.6786° latitude
in 20° total, spaced by 20/19° at a fixed 7.7401° longitude and 800 km
altitude. Position times span 2010-01-01 11:57:28 to 12:02:32 UTC at 16 s
intervals. The density field is frozen at the 12:00 UTC model epoch, so
variation along the pass is spatial.

The baseline electron density is the separately generated Fortran IRI-2016
grid. A north-south cosine wave multiplies it, with 18% maximum amplitude
at 280 km altitude, 900 km horizontal wavelength, and a 60 km vertical
Gaussian width. The 900 km wavelength is resolved by the grid's 2° latitude
spacing. The 20 interpolated sounder profiles have peak-density changes from
−15.7% to +15.0% relative to the IRI baseline. The
[density figure](figures/lat_wave_pass_density.png) shows the wave and the
sample positions; the [manifest](data/lat_wave_pass_manifest.json) records
their coordinates, times, and wave parameters.

Each position has a 2–10 MHz vertical ionogram in 0.1 MHz steps. The option D
equal-area fan, separate O/X modes, 1 km homing gate, and all accepted returns
are retained. To check adaptive homing gaps, each empty mode/frequency bin
at or below 0.6 MHz above the local peak plasma frequency received a full-fan
search; only reflected returns at 150 km group range or greater were added.
The full-fan check tested 336 empty bins and recovered 43 returns, raising
the total from 2,010 to 2,053. It changed the O-mode cutoff in 10 profiles
and the X-mode cutoff in nine. The mean absolute O cutoff difference from
the local peak plasma frequency fell from 0.46 to 0.18 MHz. This cutoff is
an accepted-return diagnostic, not an exact measurement of peak density in
the spatially varying plasma. Full-fan recovery checks empty bins; it does
not certify that every additional multipath return was found.

The [pass observables](figures/lat_wave_pass_observables.png) show mode-resolved
cutoffs and counts. A [profile 14 comparison](figures/lat_wave_pass_profile14_gap_recovery.png)
shows the nine added accepted returns using the same O-blue/X-red ionogram
display. Raw ionograms are in [the adaptive set](data/lat_wave_pass_ionograms);
the working truth ionograms are in [the recovered set](data/lat_wave_pass_ionograms_recovered).
The [summary](data/lat_wave_pass_summary.json) gives per-profile results.

## Runtime and provenance

The 20 initial ionograms completed in 222.5 s from first start to last
finish on Cartman. The 20 full-fan gap checks completed in another 134.8 s.
All 40 jobs exited successfully; the remote run, inputs, outputs, and logs
were owned by `chartat1` with owner-only permissions. The staged forward
grid's SHA-256 matched the local file.

Only electron density comes from IRI-2016 and the imposed wave. The magnetic
and collision fields are prepared with PyIRI and the same PyLap ray tracer
used by candidate ionograms. This pass is a synthetic truth set ready for a
multi-profile retrieval; no inversion has been run on it yet.
