# Server privacy

- Anything this project creates or changes on a server must be private to `chartat1`.
- Use `chartat1`'s own account and private home or scratch directories. Do not place project code, jobs, logs, results, or credentials in public or group-accessible trees.
- Set `umask 077`; create directories with mode `0700` and files with mode `0600`. Verify ownership and permissions after staging and after jobs finish.
- Keep scheduler output, temporary files, and copied dependencies under the same private directory. If a server or scheduler cannot enforce this, stop before submitting work.

# Local memory safety

- This workstation has 16 GB RAM. Run at most **one** PyLap ionogram ray-tracing process at a time on this computer, across all terminals, scripts, and agent sessions. Do not launch overlapping truth, prior, vertical, or oblique batches. Use one worker by default and reject requests for more workers locally. Prefer the private Cartman compute nodes for full ionogram passes.
- Local PyLap CLI runs are disabled by default. Set `RAYTRACING_ALLOW_LOCAL_RAYS=1` only after checking memory pressure and swap use for a justified single test. Do not start if available memory is below 6 GB or swap use is growing; stop a run immediately if swapping begins or rises. Before increasing concurrency for any future workload, measure one complete representative ionogram's peak resident memory, account for all other running jobs and the OS, and leave at least 6 GB of RAM headroom. CPU count alone is not grounds for parallel ray jobs.
- **Never ray trace the full SAMI3 wave grid locally.** A single representative ionogram measured on Cartman peaked at 27,955,576 kB (about 26.7 GiB) resident memory, exceeding this laptop's 16 GB RAM even with one worker. The earlier crashes were on this laptop, not Cartman. Use private compute jobs with measured memory capacity for this grid; keep at most two large ray jobs active until a lower peak is demonstrated. Loading the small saved regional grids for scoring or figures is allowed.
- Use a shared local lock around each ray-tracing process so separate commands cannot accidentally overlap. After interrupted work, check for surviving ray processes before restarting.
