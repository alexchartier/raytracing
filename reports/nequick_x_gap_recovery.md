# NeQuick-G truth ionogram: X-mode gap recovery

The first NeQuick-G truth ionogram (−53.6786°, 7.7401°, January 12 UT,
Galileo broadcast Az = 64) had no accepted X-mode return at 2.5, 2.6, 2.7,
2.8, or 2.9 MHz. It had X returns at 2.4 and 3.0 MHz. This was a homing-search
gap in the adaptive equal-area fan, not a modeled physical gap.

![Original and corrected NeQuick-G truth ionograms, with a detailed view of the gap](data/nequick_independent_case/nequick_truth_x_gap_recovered.png)

A dense-fan diagnostic recovered X rays at 2.6–2.8 MHz. Starting from full-fan
solutions at 2.4, 2.6, 2.8, and 3.0 MHz, `nequick_x_gap_cartman_job.sh` then
followed the branches in 20 kHz steps. It traced and homed two X rays at
**each exact 100 kHz display frequency** from 2.5 through 2.9 MHz. All ten
new returns met the original 1 km homing gate; the largest miss was 780.3 m.
Their group ranges continue smoothly from 1134.1 km at 2.5 MHz to 1145.0 km
at 2.9 MHz. `repair_nequick_x_gap.py` checks the frequencies, mode, range,
homing miss, and count before writing `truth_ionogram_x_gap_recovered.npz`.
It retains every accepted return, including both new X rays per bin.

The independent NeQuick-G density was used only by the forward ray tracer.
The recovery did not read retrieval candidates or modify the density truth.
The first retrieval's `evaluation.json` and paired truth/retrieval plot were
made against the incomplete ionogram; their score ranking and density errors
are preliminary until the candidate ionograms receive the same gap-recovery
pass and are scored again. The original file is retained as a record of the
failure.

All Cartman files were placed under `chartat1`'s private run directory with
directories at `0700`, files at `0600`, and scheduler output sent to
`/dev/null`; the job's own logs remained in the private directory. Both jobs
exited successfully, and a final ownership/permission check found no
violations.
