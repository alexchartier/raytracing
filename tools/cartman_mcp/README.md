# Private Cartman MCP

This is the Cartman MCP copied from `~/ampere-analysis/alex/pipeline/tools/cartman_mcp`
and narrowed for this raytracing project. `.codex/config.toml` launches it locally.
It connects with `ssh -T cartman` as `chartat1` and exposes only two remote zones:

- `repo`: `/homes/chartat1/private_raytracing/repo`
- `sandbox`: `/homes/chartat1/private_raytracing/runs`

The remote shell starts with `umask 077`. The MCP's directory and file writers
and qsub script staging retain owner-only permissions. Public, website, archive,
rsync, and admin tools from the source server are not exported. Remote commands
can run arbitrary shell code, so callers must still honor `AGENTS.md` and keep
all created files in the private root.

The local MCP log is under `.cache/cartman_mcp/`. Run
`tools/cartman_mcp/launch.sh --self-check` to inspect the exposed surface.

For SGE jobs, pass `-S /bin/bash` explicitly. Route scheduler stdout and stderr
to `/dev/null` and redirect inside the script after `umask 077`; SGE otherwise
creates log files with mode `0644` even in a private directory.
