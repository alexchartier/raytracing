from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.cartman_mcp.policy import REMOTE_REPO_ROOT
from tools.cartman_mcp.policy import REMOTE_HOST
from tools.cartman_mcp.policy import ZONES
from tools.cartman_mcp.policy import dump_redacted_json
from tools.cartman_mcp.policy import policy_summary_markdown
from tools.cartman_mcp.policy import resolve_bind_mount_destination
from tools.cartman_mcp.policy import resolve_bind_mount_source
from tools.cartman_mcp.policy import resolve_read_path
from tools.cartman_mcp.policy import resolve_system_file_path
from tools.cartman_mcp.policy import resolve_write_path
from tools.cartman_mcp.policy import validate_systemctl_action
from tools.cartman_mcp.policy import validate_systemctl_service
from tools.cartman_mcp.remote import RemoteCommandError
from tools.cartman_mcp.remote import run_remote_bash


SERVER_NAME = "raytracing-cartman-mcp"
SERVER_VERSION = "0.3.0"
SUPPORTED_PROTOCOL_VERSIONS = (
    "2025-11-25",
    "2025-06-18",
    "2025-03-26",
)

LOCAL_REPO_ROOT = Path(__file__).resolve().parents[2]
LOCAL_AGENTS_PATH = LOCAL_REPO_ROOT / "AGENTS.md"
PRIVATE_TOOLS = frozenset({
    "cartman_read_file", "cartman_tail_log", "cartman_list_dir",
    "cartman_stat_path", "cartman_df_path", "cartman_du_path",
    "cartman_mkdir", "cartman_write_file", "cartman_exec_remote",
    "cartman_qsub_submit", "cartman_git_status_repo",
    "cartman_git_pull_repo", "cartman_find", "cartman_search_text",
})
TRANSPORT_MODE = "framed"
WRITE_ENABLED_ZONES = tuple(sorted(name for name, zone in ZONES.items() if zone.write_roots))
EXEC_ENABLED_ZONES = tuple(zone for zone in ("repo", "sandbox") if zone in ZONES)
QSUB_ENABLED_ZONES = tuple(zone for zone in ("repo", "sandbox") if zone in ZONES)
LOCAL_RSYNC_WRITE_ROOTS = tuple(
    root.resolve()
    for root in (
        Path.home(),
        Path("/tmp"),
        Path("/private/tmp"),
    )
    if root.exists()
)


def _select_shared_runtime_root() -> Path:
    env_dir = os.environ.get("RAYTRACING_CARTMAN_MCP_LOG_DIR")
    if env_dir:
        return Path(env_dir).expanduser()
    return LOCAL_REPO_ROOT / ".cache" / "cartman_mcp"


LOG_PATH = _select_shared_runtime_root() / "cartman_mcp.log"


def _log_event(kind: str, **fields: Any) -> None:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "ts": datetime.now(timezone.utc).isoformat(),
        "kind": kind,
        **fields,
    }
    with LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def _read_local_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _read_local_json(path: Path) -> dict[str, Any]:
    return json.loads(_read_local_text(path))


def _tool_definitions() -> list[dict[str, Any]]:
    write_zones = list(WRITE_ENABLED_ZONES)
    exec_zones = list(EXEC_ENABLED_ZONES)
    qsub_zones = list(QSUB_ENABLED_ZONES)
    definitions = [
        {
            "name": "cartman_read_file",
            "description": "Read a text file from an allowed cartman zone.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": sorted(ZONES)},
                    "path": {"type": "string"},
                    "start_line": {"type": "integer", "minimum": 1, "default": 1},
                    "end_line": {"type": "integer", "minimum": 1, "default": 200},
                },
                "required": ["zone", "path"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_tail_log",
            "description": "Tail a text file from an allowed cartman zone.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": sorted(ZONES)},
                    "path": {"type": "string"},
                    "lines": {"type": "integer", "minimum": 1, "maximum": 2000, "default": 200},
                },
                "required": ["zone", "path"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_list_dir",
            "description": "List one directory level from an allowed cartman zone.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": sorted(ZONES)},
                    "path": {"type": "string"},
                    "max_entries": {"type": "integer", "minimum": 1, "maximum": 500, "default": 200},
                },
                "required": ["zone", "path"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_stat_path",
            "description": "Return basic stat metadata for a path in an allowed cartman zone.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": sorted(ZONES)},
                    "path": {"type": "string"},
                },
                "required": ["zone", "path"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_df_path",
            "description": "Report filesystem capacity for the filesystem containing a path in an allowed cartman zone.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": sorted(ZONES)},
                    "path": {"type": "string"},
                },
                "required": ["zone", "path"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_du_path",
            "description": "Report total disk usage for a path in an allowed cartman zone.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": sorted(ZONES)},
                    "path": {"type": "string"},
                },
                "required": ["zone", "path"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_mkdir",
            "description": "Create a directory in an allowed writable zone.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": sorted(ZONES)},
                    "path": {"type": "string"},
                },
                "required": ["zone", "path"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_write_file",
            "description": "Write or append text to a file in an allowed writable zone.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": sorted(ZONES)},
                    "path": {"type": "string"},
                    "text": {"type": "string"},
                    "append": {"type": "boolean", "default": False},
                },
                "required": ["zone", "path", "text"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_bind_mount_ro",
            "description": "Create or replace a read-only bind mount for an allowlisted AMPERE public mirror path.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "source_path": {"type": "string"},
                    "target_path": {"type": "string"},
                    "replace_existing": {"type": "boolean", "default": False},
                },
                "required": ["source_path", "target_path"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_install_system_file",
            "description": "Install text into an allowlisted system file path using sudo install.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "path": {"type": "string"},
                    "text": {"type": "string"},
                    "mode": {"type": "string", "default": "0644"},
                },
                "required": ["path", "text"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_systemctl",
            "description": "Run an allowlisted systemctl action on an allowlisted service.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "service": {"type": "string"},
                    "action": {"type": "string"},
                    "timeout_seconds": {"type": "integer", "minimum": 1, "maximum": 600, "default": 120},
                },
                "required": ["service", "action"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_exec_remote",
            "description": "Run a remote shell command from an allowed writable working directory on cartman. Detached runs require explicit log paths.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": exec_zones},
                    "cwd": {"type": "string", "default": "."},
                    "command": {"type": "string"},
                    "timeout_seconds": {"type": "integer", "minimum": 1, "maximum": 86400, "default": 300},
                    "detach": {"type": "boolean", "default": False},
                    "stdout_path": {"type": "string"},
                    "stderr_path": {"type": "string"},
                },
                "required": ["zone", "command"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_qsub_submit",
            "description": "Submit a qsub job from repo or sandbox. Optionally write the job script first under an allowed writable path.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": qsub_zones},
                    "cwd": {"type": "string", "default": "."},
                    "script_path": {"type": "string"},
                    "script_text": {"type": "string"},
                    "qsub_args": {
                        "type": "array",
                        "items": {"type": "string"},
                        "default": [],
                    },
                    "timeout_seconds": {"type": "integer", "minimum": 1, "maximum": 600, "default": 120},
                },
                "required": ["zone", "script_path"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_rsync_from_remote",
            "description": "Copy files or directories from an allowed remote Cartman read path to a guarded local destination using rsync.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": sorted(ZONES)},
                    "remote_path": {"type": "string"},
                    "local_path": {"type": "string"},
                    "dry_run": {"type": "boolean", "default": False},
                    "delete": {"type": "boolean", "default": False},
                    "timeout_seconds": {"type": "integer", "minimum": 1, "maximum": 86400, "default": 3600},
                },
                "required": ["zone", "remote_path", "local_path"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_git_status_repo",
            "description": "Report git status for the cartman pipeline checkout or a subdirectory within it.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "path": {"type": "string", "default": "."},
                },
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_git_pull_repo",
            "description": "Run `git pull --ff-only` in the cartman pipeline checkout or a subdirectory within it.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "path": {"type": "string", "default": "."},
                    "remote": {"type": "string", "default": "origin"},
                    "branch": {"type": "string", "default": "main"},
                },
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_find",
            "description": "Find files or directories by case-insensitive basename pattern.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": sorted(ZONES)},
                    "root": {"type": "string"},
                    "pattern": {"type": "string"},
                    "max_results": {"type": "integer", "minimum": 1, "maximum": 500, "default": 200},
                },
                "required": ["zone", "root", "pattern"],
                "additionalProperties": False,
            },
        },
        {
            "name": "cartman_search_text",
            "description": "Search text files recursively for a literal string under an allowed read path.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "zone": {"type": "string", "enum": sorted(ZONES)},
                    "root": {"type": "string"},
                    "query": {"type": "string", "minLength": 1},
                    "max_results": {"type": "integer", "minimum": 1, "maximum": 500, "default": 200},
                },
                "required": ["zone", "root", "query"],
                "additionalProperties": False,
            },
        },
    ]
    return [definition for definition in definitions if definition["name"] in PRIVATE_TOOLS]


def _resources() -> list[dict[str, Any]]:
    return [
        {
            "uri": "policy://summary",
            "name": "Cartman policy summary",
            "description": "Effective path boundaries and write rules for the Cartman MCP server.",
            "mimeType": "text/markdown",
        },
        {
            "uri": "notes://agents",
            "name": "Project AGENTS.md",
            "description": "Local server privacy rule for the raytracing project.",
            "mimeType": "text/markdown",
        },
    ]


def _prompts() -> list[dict[str, Any]]:
    return []


def _resource_contents(uri: str) -> dict[str, Any]:
    if uri == "policy://summary":
        return {"contents": [{"uri": uri, "mimeType": "text/markdown", "text": policy_summary_markdown()}]}
    if uri == "notes://agents":
        return {"contents": [{"uri": uri, "mimeType": "text/markdown", "text": _read_local_text(LOCAL_AGENTS_PATH)}]}
    raise ValueError(f"Unknown resource: {uri}")


def _prompt_contents(name: str, arguments: dict[str, Any]) -> dict[str, Any]:
    if name == "triage_failed_job":
        job_id = arguments["job_id"]
        text = (
            f"Investigate cartman job `{job_id}`.\n"
            f"Start with `cartman_find` or `cartman_list_dir` under `repo` for `qsub_logs`, then use "
            f"`cartman_tail_log` and `cartman_read_file` on matching stdout/stderr files. "
            f"Summarize the failure mode, impacted dates, and the next safe read-only checks."
        )
    elif name == "prepare_backfill_submission":
        start_date = arguments["start_date"]
        end_date = arguments["end_date"]
        text = (
            f"Plan a backfill for `{start_date}` through `{end_date}` without submitting anything.\n"
            f"Review `config://remote`, `notes://agents`, and [qsub_run_pipeline_multi.sh]({REMOTE_REPO_ROOT}/qsub_run_pipeline_multi.sh) "
            f"before proposing chunking, host exclusions, and any pipeline overrides."
        )
    else:
        raise ValueError(f"Unknown prompt: {name}")
    return {
        "description": text,
        "messages": [
            {
                "role": "user",
                "content": {"type": "text", "text": text},
            }
        ],
    }


def _normalize_line_range(start_line: Any, end_line: Any) -> tuple[int, int]:
    start = int(start_line or 1)
    end = int(end_line or max(start, 200))
    if start < 1 or end < start:
        raise ValueError("Line bounds must satisfy 1 <= start_line <= end_line")
    return start, min(end, start + 1999)


def _make_text_result(payload: dict[str, Any]) -> dict[str, Any]:
    text = json.dumps(payload, indent=2, sort_keys=True)
    return {
        "content": [{"type": "text", "text": text}],
        "structuredContent": payload,
        "isError": False,
    }


def _parse_df_output(text: str) -> dict[str, Any]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError(f"Unexpected df output: {text!r}")
    fields = lines[-1].split(None, 5)
    if len(fields) != 6:
        raise ValueError(f"Unexpected df fields: {lines[-1]!r}")
    total_kib = int(fields[1])
    used_kib = int(fields[2])
    available_kib = int(fields[3])
    return {
        "filesystem": fields[0],
        "total_kib": total_kib,
        "used_kib": used_kib,
        "available_kib": available_kib,
        "total_bytes": total_kib * 1024,
        "used_bytes": used_kib * 1024,
        "available_bytes": available_kib * 1024,
        "capacity_percent": int(fields[4].rstrip("%")),
        "mount_point": fields[5],
    }


def _parse_du_output(text: str) -> dict[str, Any]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"Unexpected du output: {text!r}")
    fields = lines[-1].split(None, 1)
    if len(fields) != 2:
        raise ValueError(f"Unexpected du fields: {lines[-1]!r}")
    size_kib = int(fields[0])
    return {
        "size_kib": size_kib,
        "size_bytes": size_kib * 1024,
        "reported_path": fields[1],
    }


def _basename_pattern(pattern: str) -> str:
    if any(char in pattern for char in "*?[]"):
        return pattern
    return f"*{pattern}*"


def _shell_single_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def _join_shell_args(values: list[str]) -> str:
    return " ".join(_shell_single_quote(str(value)) for value in values)


def _normalize_timeout(value: Any, *, default: int, min_value: int, max_value: int) -> int:
    timeout = default if value is None else int(value)
    if timeout < min_value or timeout > max_value:
        raise ValueError(f"timeout_seconds must be between {min_value} and {max_value}")
    return timeout


def _parse_job_id(text: str) -> str | None:
    stripped = text.strip()
    if not stripped:
        return None
    match = re.match(r"^(\S+)", stripped)
    return match.group(1) if match else None


def _is_local_path_under(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def _resolve_local_rsync_path(path_text: str) -> Path:
    if not path_text:
        raise ValueError("local_path must be non-empty")
    path = Path(path_text).expanduser()
    if not path.is_absolute():
        raise ValueError("local_path must be an absolute local path")
    resolved = path.resolve(strict=False)
    if not any(_is_local_path_under(resolved, root) for root in LOCAL_RSYNC_WRITE_ROOTS):
        roots = ", ".join(str(root) for root in LOCAL_RSYNC_WRITE_ROOTS)
        raise ValueError(f"local_path {resolved} is outside allowed local write roots: {roots}")
    return resolved


def _run_rsync_from_remote(
    remote_path: str,
    local_path: Path,
    *,
    dry_run: bool,
    delete: bool,
    timeout_seconds: int,
) -> subprocess.CompletedProcess[str]:
    if local_path.exists() and local_path.is_dir():
        local_path.mkdir(parents=True, exist_ok=True)
    elif str(local_path).endswith(os.sep):
        local_path.mkdir(parents=True, exist_ok=True)
    else:
        local_path.parent.mkdir(parents=True, exist_ok=True)

    command = ["rsync", "-a", "--human-readable", "--stats"]
    if dry_run:
        command.extend(["--dry-run", "--itemize-changes"])
    if delete:
        command.append("--delete")
    command.extend([f"{REMOTE_HOST}:{remote_path}", str(local_path)])
    return subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
        check=False,
    )


def _call_tool(name: str, arguments: dict[str, Any]) -> dict[str, Any]:
    if name not in PRIVATE_TOOLS:
        raise ValueError(f"Tool {name} is unavailable in the private raytracing MCP")
    zone = arguments.get("zone", "repo")
    if name == "cartman_read_file":
        path = resolve_read_path(zone, arguments["path"])
        start, end = _normalize_line_range(arguments.get("start_line"), arguments.get("end_line"))
        command = f"sed -n {shlex.quote(f'{start},{end}p')} -- {shlex.quote(path)}"
        result = run_remote_bash(command)
        return _make_text_result(
            {
                "zone": zone,
                "path": path,
                "start_line": start,
                "end_line": end,
                "text": result.stdout,
            }
        )
    if name == "cartman_tail_log":
        path = resolve_read_path(zone, arguments["path"])
        lines = min(max(int(arguments.get("lines", 200)), 1), 2000)
        command = f"tail -n {lines} -- {shlex.quote(path)}"
        result = run_remote_bash(command)
        return _make_text_result({"zone": zone, "path": path, "lines": lines, "text": result.stdout})
    if name == "cartman_list_dir":
        path = resolve_read_path(zone, arguments["path"])
        max_entries = min(max(int(arguments.get("max_entries", 200)), 1), 500)
        command = (
            f"find {shlex.quote(path)} -mindepth 1 -maxdepth 1 "
            "-printf '%f\t%y\t%s\t%TY-%Tm-%TdT%TH:%TM:%TS\n' | sort | head -n "
            f"{max_entries}"
        )
        result = run_remote_bash(command)
        entries = []
        for line in result.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) != 4:
                continue
            entries.append(
                {
                    "name": parts[0],
                    "kind": parts[1],
                    "size_bytes": int(parts[2]),
                    "mtime": parts[3],
                }
            )
        return _make_text_result({"zone": zone, "path": path, "entries": entries, "count": len(entries)})
    if name == "cartman_stat_path":
        path = resolve_read_path(zone, arguments["path"])
        command = f"stat -c '%F\t%s\t%Y\t%n' -- {shlex.quote(path)}"
        result = run_remote_bash(command)
        fields = result.stdout.rstrip("\n").split("\t", 3)
        payload = {
            "zone": zone,
            "path": path,
            "file_type": fields[0] if len(fields) > 0 else "",
            "size_bytes": int(fields[1]) if len(fields) > 1 and fields[1] else 0,
            "mtime_epoch": int(fields[2]) if len(fields) > 2 and fields[2] else 0,
            "reported_path": fields[3] if len(fields) > 3 else path,
        }
        return _make_text_result(payload)
    if name == "cartman_df_path":
        path = resolve_read_path(zone, arguments["path"])
        command = f"LC_ALL=C df -Pk -- {shlex.quote(path)}"
        result = run_remote_bash(command)
        payload = {
            "zone": zone,
            "path": path,
            **_parse_df_output(result.stdout),
        }
        return _make_text_result(payload)
    if name == "cartman_du_path":
        path = resolve_read_path(zone, arguments["path"])
        command = f"LC_ALL=C du -sk -- {shlex.quote(path)}"
        result = run_remote_bash(command)
        payload = {
            "zone": zone,
            "path": path,
            **_parse_du_output(result.stdout),
        }
        return _make_text_result(payload)
    if name == "cartman_mkdir":
        path = resolve_write_path(zone, arguments["path"])
        command = f"mkdir -p -- {shlex.quote(path)}"
        run_remote_bash(command)
        return _make_text_result({"zone": zone, "path": path, "created": True})
    if name == "cartman_write_file":
        path = resolve_write_path(zone, arguments["path"])
        text = str(arguments["text"])
        append = bool(arguments.get("append", False))
        path_parent = str(Path(path).parent).replace("\\", "/")
        redirect = ">>" if append else ">"
        command = (
            f"mkdir -p -- {shlex.quote(path_parent)}\n"
            f"cat {redirect} {shlex.quote(path)} <<'__CARTMAN_EOF__'\n"
            f"{text}\n"
            "__CARTMAN_EOF__"
        )
        run_remote_bash(command)
        return _make_text_result({"zone": zone, "path": path, "bytes_written": len(text.encode('utf-8')), "append": append})
    if name == "cartman_bind_mount_ro":
        source_path = resolve_bind_mount_source(arguments["source_path"])
        target_path = resolve_bind_mount_destination(arguments["target_path"])
        replace_existing = bool(arguments.get("replace_existing", False))
        target_parent = str(Path(target_path).parent).replace("\\", "/")
        commands = [f"sudo mkdir -p -- {shlex.quote(target_parent)} {shlex.quote(target_path)}"]
        commands.append(f"if mountpoint -q -- {shlex.quote(target_path)}; then")
        if replace_existing:
            commands.extend(
                [
                    f"  sudo umount -- {shlex.quote(target_path)}",
                    "else",
                    "  true",
                    "fi",
                ]
            )
        else:
            commands.extend(
                [
                    "  echo 'mountpoint already active' >&2",
                    "  exit 2",
                    "else",
                    "  true",
                    "fi",
                ]
            )
        commands.extend(
            [
                f"sudo mount --bind {shlex.quote(source_path)} {shlex.quote(target_path)}",
                f"sudo mount -o remount,bind,ro {shlex.quote(target_path)}",
                f"mount | grep -F -- {_shell_single_quote(' on ' + target_path + ' ')}",
            ]
        )
        result = run_remote_bash("\n".join(commands), timeout_seconds=120)
        return _make_text_result(
            {
                "source_path": source_path,
                "target_path": target_path,
                "replace_existing": replace_existing,
                "text": result.stdout,
            }
        )
    if name == "cartman_install_system_file":
        path = resolve_system_file_path(arguments["path"])
        text = str(arguments["text"])
        mode = str(arguments.get("mode", "0644"))
        if not re.fullmatch(r"0?[0-7]{3,4}", mode):
            raise ValueError("mode must be an octal string such as 0644")
        command = (
            "tmp=$(mktemp)\n"
            "trap 'rm -f \"$tmp\"' EXIT\n"
            "cat > \"$tmp\" <<'__CARTMAN_EOF__'\n"
            f"{text}\n"
            "__CARTMAN_EOF__\n"
            f"sudo install -D -m {shlex.quote(mode)} \"$tmp\" {shlex.quote(path)}\n"
            f"sudo stat -c '%a\t%n' -- {shlex.quote(path)}"
        )
        result = run_remote_bash(command, timeout_seconds=120)
        return _make_text_result(
            {
                "path": path,
                "mode": mode,
                "bytes_written": len(text.encode("utf-8")),
                "text": result.stdout,
            }
        )
    if name == "cartman_systemctl":
        service = validate_systemctl_service(arguments["service"])
        action = validate_systemctl_action(arguments["action"])
        timeout = _normalize_timeout(arguments.get("timeout_seconds"), default=120, min_value=1, max_value=600)
        if action == "enable-now":
            command = f"sudo systemctl enable --now {shlex.quote(service)}"
        elif action == "daemon-reload":
            command = "sudo systemctl daemon-reload"
        else:
            command = f"sudo systemctl {shlex.quote(action)} {shlex.quote(service)}"
        result = run_remote_bash(command, timeout_seconds=timeout)
        return _make_text_result(
            {
                "service": service,
                "action": action,
                "text": result.stdout,
            }
        )
    if name == "cartman_exec_remote":
        cwd = resolve_write_path(zone, arguments.get("cwd", "."))
        command_text = str(arguments["command"])
        detach = bool(arguments.get("detach", False))
        timeout = _normalize_timeout(arguments.get("timeout_seconds"), default=300, min_value=1, max_value=86400)
        if detach:
            stdout_arg = arguments.get("stdout_path")
            stderr_arg = arguments.get("stderr_path")
            if not stdout_arg or not stderr_arg:
                raise ValueError("Detached runs require both stdout_path and stderr_path")
            stdout_path = resolve_write_path(zone, stdout_arg)
            stderr_path = resolve_write_path(zone, stderr_arg)
            parent_dirs = sorted({str(Path(stdout_path).parent).replace("\\", "/"), str(Path(stderr_path).parent).replace("\\", "/")})
            prelude = [f"cd -- {shlex.quote(cwd)}"] + [f"mkdir -p -- {shlex.quote(parent)}" for parent in parent_dirs]
            prelude.extend(
                [
                    f"nohup /bin/bash --noprofile --norc -c {_shell_single_quote(command_text)} > {shlex.quote(stdout_path)} 2> {shlex.quote(stderr_path)} < /dev/null &",
                    "pid=$!",
                    "printf '%s\\n' \"$pid\"",
                ]
            )
            result = run_remote_bash("\n".join(prelude), timeout_seconds=min(timeout, 60))
            pid_text = result.stdout.strip()
            pid = int(pid_text) if pid_text.isdigit() else None
            return _make_text_result(
                {
                    "zone": zone,
                    "cwd": cwd,
                    "command": command_text,
                    "detached": True,
                    "pid": pid,
                    "stdout_path": stdout_path,
                    "stderr_path": stderr_path,
                }
            )
        result = run_remote_bash(f"cd -- {shlex.quote(cwd)}\n{command_text}", timeout_seconds=timeout)
        return _make_text_result(
            {
                "zone": zone,
                "cwd": cwd,
                "command": command_text,
                "detached": False,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "returncode": result.returncode,
            }
        )
    if name == "cartman_qsub_submit":
        cwd = resolve_write_path(zone, arguments.get("cwd", "."))
        script_path = resolve_write_path(zone, arguments["script_path"])
        script_text = arguments.get("script_text")
        qsub_args = [str(value) for value in arguments.get("qsub_args", [])]
        timeout = _normalize_timeout(arguments.get("timeout_seconds"), default=120, min_value=1, max_value=600)
        commands = [f"cd -- {shlex.quote(cwd)}"]
        wrote_script = "script_text" in arguments
        if wrote_script:
            script_body = "" if script_text is None else str(script_text)
            script_parent = str(Path(script_path).parent).replace("\\", "/")
            commands.extend(
                [
                    f"mkdir -p -- {shlex.quote(script_parent)}",
                    f"cat > {shlex.quote(script_path)} <<'__CARTMAN_EOF__'\n{script_body}\n__CARTMAN_EOF__",
                    f"chmod 0700 -- {shlex.quote(script_path)}",
                ]
            )
        qsub_command = "qsub"
        if qsub_args:
            qsub_command += f" {_join_shell_args(qsub_args)}"
        qsub_command += f" {shlex.quote(script_path)}"
        commands.append(qsub_command)
        result = run_remote_bash("\n".join(commands), timeout_seconds=timeout)
        return _make_text_result(
            {
                "zone": zone,
                "cwd": cwd,
                "script_path": script_path,
                "wrote_script": wrote_script,
                "bytes_written": len(str(script_text).encode("utf-8")) if wrote_script and script_text is not None else 0,
                "qsub_args": qsub_args,
                "job_id": _parse_job_id(result.stdout),
                "text": result.stdout,
            }
        )
    if name == "cartman_rsync_from_remote":
        remote_path = resolve_read_path(zone, arguments["remote_path"])
        local_path = _resolve_local_rsync_path(str(arguments["local_path"]))
        dry_run = bool(arguments.get("dry_run", False))
        delete = bool(arguments.get("delete", False))
        timeout = _normalize_timeout(arguments.get("timeout_seconds"), default=3600, min_value=1, max_value=86400)
        result = _run_rsync_from_remote(
            remote_path,
            local_path,
            dry_run=dry_run,
            delete=delete,
            timeout_seconds=timeout,
        )
        payload = {
            "zone": zone,
            "remote_path": remote_path,
            "local_path": str(local_path),
            "dry_run": dry_run,
            "delete": delete,
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
        if result.returncode != 0:
            return {
                "content": [{"type": "text", "text": json.dumps(payload, indent=2, sort_keys=True)}],
                "structuredContent": payload,
                "isError": True,
            }
        return _make_text_result(payload)
    if name == "cartman_git_status_repo":
        path = resolve_read_path("repo", arguments.get("path", "."))
        command = (
            f"cd -- {shlex.quote(path)}\n"
            "git status --short --branch"
        )
        result = run_remote_bash(command)
        return _make_text_result({"zone": "repo", "path": path, "text": result.stdout})
    if name == "cartman_git_pull_repo":
        path = resolve_write_path("repo", arguments.get("path", "."))
        remote = str(arguments.get("remote", "origin"))
        branch = str(arguments.get("branch", "main"))
        command = (
            f"cd -- {shlex.quote(path)}\n"
            f"git pull --ff-only {_shell_single_quote(remote)} {_shell_single_quote(branch)}"
        )
        result = run_remote_bash(command, timeout_seconds=300)
        return _make_text_result(
            {
                "zone": "repo",
                "path": path,
                "remote": remote,
                "branch": branch,
                "text": result.stdout,
            }
        )
    if name == "cartman_find":
        root = resolve_read_path(zone, arguments["root"])
        pattern = _basename_pattern(arguments["pattern"])
        max_results = min(max(int(arguments.get("max_results", 200)), 1), 500)
        command = (
            f"find {shlex.quote(root)} \\( -type f -o -type d \\) "
            f"-iname {shlex.quote(pattern)} | sort | head -n {max_results}"
        )
        result = run_remote_bash(command)
        matches = result.stdout.splitlines()
        return _make_text_result(
            {
                "zone": zone,
                "root": root,
                "pattern": pattern,
                "matches": matches,
                "count": len(matches),
            }
        )
    if name == "cartman_search_text":
        root = resolve_read_path(zone, arguments["root"])
        query = str(arguments["query"])
        if not query:
            raise ValueError("query must be non-empty")
        max_results = min(max(int(arguments.get("max_results", 200)), 1), 500)
        command = (
            "LC_ALL=C grep -rInI -F --exclude-dir=.git --exclude='*.o' --exclude='*.a' --exclude='*.so' "
            f"-- {_shell_single_quote(query)} {_shell_single_quote(root)} | head -n {max_results}"
        )
        result = run_remote_bash(command)
        matches = result.stdout.splitlines()
        return _make_text_result(
            {
                "zone": zone,
                "root": root,
                "query": query,
                "matches": matches,
                "count": len(matches),
            }
        )
    raise ValueError(f"Unknown tool: {name}")


def _success_response(request_id: Any, result: dict[str, Any]) -> dict[str, Any]:
    return {"jsonrpc": "2.0", "id": request_id, "result": result}


def _error_response(request_id: Any, code: int, message: str) -> dict[str, Any]:
    return {"jsonrpc": "2.0", "id": request_id, "error": {"code": code, "message": message}}


def _read_message() -> dict[str, Any] | None:
    global TRANSPORT_MODE
    first_line = sys.stdin.buffer.readline()
    if not first_line:
        _log_event("eof")
        return None

    _log_event("readline", sample=first_line[:200].decode("utf-8", "replace"))
    stripped = first_line.lstrip()
    if stripped.startswith(b"{"):
        TRANSPORT_MODE = "line"
        _log_event("transport", mode=TRANSPORT_MODE)
        chunks = [first_line]
        decoder = json.JSONDecoder()
        while True:
            text = b"".join(chunks).decode("utf-8", "replace")
            try:
                obj, index = decoder.raw_decode(text)
            except json.JSONDecodeError:
                next_line = sys.stdin.buffer.readline()
                if not next_line:
                    raise ValueError("Unexpected EOF while reading newline-delimited JSON-RPC message")
                _log_event("readline", sample=next_line[:200].decode("utf-8", "replace"))
                chunks.append(next_line)
                continue
            remainder = text[index:].strip()
            if remainder:
                _log_event("json_remainder", sample=remainder[:200])
            return obj

    headers = {}
    TRANSPORT_MODE = "framed"
    _log_event("transport", mode=TRANSPORT_MODE)
    line = first_line
    while True:
        if line in (b"\r\n", b"\n"):
            break
        decoded = line.decode("utf-8")
        if ":" not in decoded:
            raise ValueError(f"Malformed MCP header line: {decoded!r}")
        key, value = decoded.split(":", 1)
        headers[key.strip().lower()] = value.strip()
        line = sys.stdin.buffer.readline()
        if not line:
            raise ValueError("Unexpected EOF while reading MCP headers")
        _log_event("readline", sample=line[:200].decode("utf-8", "replace"))
    length = int(headers["content-length"])
    payload = sys.stdin.buffer.read(length)
    _log_event("payload", sample=payload[:200].decode("utf-8", "replace"), size=len(payload))
    return json.loads(payload.decode("utf-8"))


def _write_message(message: dict[str, Any]) -> None:
    data = json.dumps(message, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    if TRANSPORT_MODE == "line":
        sys.stdout.buffer.write(data)
        sys.stdout.buffer.write(b"\n")
    else:
        header = f"Content-Length: {len(data)}\r\n\r\n".encode("ascii")
        sys.stdout.buffer.write(header)
        sys.stdout.buffer.write(data)
    sys.stdout.buffer.flush()
    _log_event("response", id=message.get("id"), has_error="error" in message, mode=TRANSPORT_MODE)


def _handle_request(message: dict[str, Any]) -> dict[str, Any] | None:
    request_id = message.get("id")
    method = message.get("method")
    params = message.get("params", {})
    if method == "initialize":
        _log_event("request", method=method, request_id=request_id, params=params)
    else:
        _log_event("request", method=method, request_id=request_id)
    try:
        if method == "initialize":
            requested_version = str(params.get("protocolVersion", "") or "")
            if requested_version in SUPPORTED_PROTOCOL_VERSIONS:
                negotiated_version = requested_version
            else:
                negotiated_version = SUPPORTED_PROTOCOL_VERSIONS[0]
            return _success_response(
                request_id,
                {
                    "protocolVersion": negotiated_version,
                    "serverInfo": {"name": SERVER_NAME, "title": "Cartman MCP", "version": SERVER_VERSION},
                    "capabilities": {
                        "tools": {"listChanged": False},
                        "resources": {"subscribe": False, "listChanged": False},
                        "prompts": {"listChanged": False},
                    },
                    "instructions": "Use the cartman tools for structured remote inspection, scoped writes, remote command execution from writable directories, and qsub submission through the Cartman MCP policy.",
                },
            )
        if method == "ping":
            return _success_response(request_id, {})
        if method == "tools/list":
            return _success_response(request_id, {"tools": _tool_definitions()})
        if method == "tools/call":
            return _success_response(request_id, _call_tool(params["name"], params.get("arguments", {})))
        if method == "resources/list":
            return _success_response(request_id, {"resources": _resources()})
        if method == "resources/read":
            return _success_response(request_id, _resource_contents(params["uri"]))
        if method == "prompts/list":
            return _success_response(request_id, {"prompts": _prompts()})
        if method == "prompts/get":
            return _success_response(request_id, _prompt_contents(params["name"], params.get("arguments", {})))
        if method in {"notifications/initialized", "initialized"}:
            _log_event("initialized")
            return None
        if method == "notifications/cancelled":
            _log_event("cancelled", params=params)
            return None
        return _error_response(request_id, -32601, f"Method not found: {method}")
    except RemoteCommandError as exc:
        return _success_response(
            request_id,
            {
                "content": [{"type": "text", "text": str(exc)}],
                "structuredContent": {
                    "error": str(exc),
                    "returncode": exc.result.returncode,
                    "stdout": exc.result.stdout,
                    "stderr": exc.result.stderr,
                },
                "isError": True,
            },
        )
    except Exception as exc:
        return _error_response(request_id, -32000, str(exc))


def _serve_stdio() -> int:
    while True:
        message = _read_message()
        if message is None:
            return 0
        response = _handle_request(message)
        if response is not None and "id" in response:
            _write_message(response)


def _self_check() -> int:
    resource_uris = [resource["uri"] for resource in _resources()]
    payload = {
        "server": SERVER_NAME,
        "version": SERVER_VERSION,
        "protocol_versions": list(SUPPORTED_PROTOCOL_VERSIONS),
        "repo_root": str(LOCAL_REPO_ROOT),
        "tool_count": len(_tool_definitions()),
        "resource_uris": resource_uris,
        "prompt_names": [prompt["name"] for prompt in _prompts()],
        "log_path": str(LOG_PATH),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


def main(argv: list[str] | None = None) -> int:
    _log_event("startup", argv=list(argv or sys.argv[1:]), cwd=str(Path.cwd()), pid=os.getpid())
    parser = argparse.ArgumentParser(description="Cartman MCP server")
    parser.add_argument("--self-check", action="store_true", help="Validate local server metadata and exit.")
    args = parser.parse_args(argv)
    if args.self_check:
        return _self_check()
    return _serve_stdio()


if __name__ == "__main__":
    raise SystemExit(main())
