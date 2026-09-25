from __future__ import annotations

from dataclasses import dataclass
from pathlib import PurePosixPath
import shlex
import subprocess

from tools.cartman_mcp.policy import REMOTE_HOST
from tools.cartman_mcp.policy import REMOTE_PYTHON
from tools.cartman_mcp.policy import REMOTE_PYTHON_COMPAT_ROOT
from tools.cartman_mcp.policy import REMOTE_SGE_CELL
from tools.cartman_mcp.policy import REMOTE_SGE_ROOT

REMOTE_PATH_PREFIXES = (
    str(PurePosixPath(REMOTE_PYTHON).parent),
    "/opt/ge-GE2011.11-11p1/bin/linux-x64",
)


@dataclass
class RemoteResult:
    command: list[str]
    returncode: int
    stdout: str
    stderr: str


class RemoteCommandError(RuntimeError):
    def __init__(self, result: RemoteResult):
        message = result.stderr.strip() or result.stdout.strip() or f"Remote command failed with {result.returncode}"
        super().__init__(message)
        self.result = result


def _bootstrap_remote_shell() -> str:
    path_prefix = ":".join(shlex.quote(path) for path in dict.fromkeys(REMOTE_PATH_PREFIXES))
    lines = ["umask 077", f"export PATH={path_prefix}:$PATH"]
    lines.extend(
        [
            "if [ -z \"${HOSTNAME:-}\" ]; then",
            "    export HOSTNAME=\"$(hostname -s 2>/dev/null || hostname 2>/dev/null || uname -n)\"",
            "fi",
            f"if [ -z \"${{SGE_ROOT:-}}\" ] && [ -d {shlex.quote(REMOTE_SGE_ROOT)} ]; then",
            f"    export SGE_ROOT={shlex.quote(REMOTE_SGE_ROOT)}",
            "fi",
            "if [ -n \"${SGE_ROOT:-}\" ] && [ -z \"${SGE_CELL:-}\" ]; then",
            f"    export SGE_CELL={shlex.quote(REMOTE_SGE_CELL)}",
            "fi",
            f"if [ -d {shlex.quote(REMOTE_PYTHON_COMPAT_ROOT)} ]; then",
            f"    export PYTHONPATH={shlex.quote(REMOTE_PYTHON_COMPAT_ROOT)}${{PYTHONPATH:+:$PYTHONPATH}}",
            "fi",
        ]
    )
    return "\n".join(lines)


def run_remote_bash(script: str, timeout_seconds: int = 60) -> RemoteResult:
    script = f"{_bootstrap_remote_shell()}\n{script}"
    command = ["ssh", "-T", REMOTE_HOST, "/bin/bash", "--noprofile", "--norc", "-s", "--"]
    completed = subprocess.run(
        command,
        input=script,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
        check=False,
    )
    result = RemoteResult(
        command=command,
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
    )
    if completed.returncode != 0:
        raise RemoteCommandError(result)
    return result
