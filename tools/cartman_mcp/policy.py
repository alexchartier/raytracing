from __future__ import annotations

from dataclasses import dataclass
from pathlib import PurePosixPath
from typing import Iterable
import json


REMOTE_HOST = "cartman"
PRIVATE_ROOT = "/homes/chartat1/private_raytracing"
REMOTE_REPO_ROOT = f"{PRIVATE_ROOT}/repo"
REMOTE_HOME_ROOT = "/homes/chartat1"
REMOTE_PYTHON = "/project/ampere/software/python/bin/python3"
REMOTE_SGE_ROOT = "/opt/ge-GE2011.11-11p1"
REMOTE_SGE_CELL = "default"
REMOTE_PYTHON_COMPAT_ROOT = f"{REMOTE_REPO_ROOT}/tools/cartman_mcp/python_compat"
REMOTE_SPICEYPY_FALLBACK = (
    "/project/ampere/software/python/venv_3.11/lib/python3.11/site-packages/spiceypy"
)

ARCHIVE_ROOT = "/project/ampere/data"
TEST_ROOT = "/project/ampere/test"
PUBLIC_ROOT = "/project/ampere/public"
SOFTWARE_ROOT = "/project/ampere/software"
SANDBOX_ROOTS = (f"{PRIVATE_ROOT}/runs",)
WWW_ROOT = "/project/ampere/www"
SUPERMAG_ROOT = "/disks/d0510/project/supermag/deployed/data"
ALLOWED_BIND_MOUNT_SOURCES = (
    "/project/ampere/sandbox/fits",
)
ALLOWED_BIND_MOUNT_DESTINATIONS = (
    f"{PUBLIC_ROOT}/hapi",
)
ALLOWED_SYSTEM_FILE_PATHS = (
    "/etc/systemd/system/ampere-hapi.service",
    "/etc/httpd/conf.d/ampere-hapi.conf",
    "/etc/apache2/conf-enabled/ampere-hapi.conf",
    "/etc/apache2/sites-available/ampere-hapi.conf",
)
ALLOWED_SYSTEMCTL_SERVICES = (
    "ampere-hapi.service",
    "httpd",
    "apache2",
)
ALLOWED_SYSTEMCTL_ACTIONS = (
    "daemon-reload",
    "enable",
    "enable-now",
    "start",
    "stop",
    "restart",
    "reload",
    "status",
)


@dataclass(frozen=True)
class Zone:
    name: str
    description: str
    read_roots: tuple[str, ...]
    write_roots: tuple[str, ...]


_SOURCE_ZONES: dict[str, Zone] = {
    "repo": Zone(
        name="repo",
        description="Pipeline repo under the cartman checkout.",
        read_roots=(REMOTE_REPO_ROOT,),
        write_roots=(REMOTE_REPO_ROOT,),
    ),
    "homes": Zone(
        name="homes",
        description="Read-only /homes/ampere2 tree.",
        read_roots=(REMOTE_HOME_ROOT,),
        write_roots=(),
    ),
    "archive": Zone(
        name="archive",
        description="Read-only archive inputs under /project/ampere/data.",
        read_roots=(ARCHIVE_ROOT,),
        write_roots=(),
    ),
    "test": Zone(
        name="test",
        description="Read-only test products under /project/ampere/test.",
        read_roots=(TEST_ROOT,),
        write_roots=(),
    ),
    "software": Zone(
        name="software",
        description="Read-only legacy AMPERE software under /project/ampere/software.",
        read_roots=(SOFTWARE_ROOT,),
        write_roots=(),
    ),
    "sandbox": Zone(
        name="sandbox",
        description="Generated products and writable sandbox data roots.",
        read_roots=SANDBOX_ROOTS,
        write_roots=SANDBOX_ROOTS,
    ),
    "public": Zone(
        name="public",
        description="Public-facing non-website data roots, including read-only mirror targets.",
        read_roots=(PUBLIC_ROOT,),
        write_roots=(PUBLIC_ROOT,),
    ),
    "www": Zone(
        name="www",
        description="Website tree; reads are broad and writes are allowed.",
        read_roots=(WWW_ROOT,),
        write_roots=(WWW_ROOT,),
    ),
    "supermag": Zone(
        name="supermag",
        description="Read-only SuperMAG deployed data tree.",
        read_roots=(SUPERMAG_ROOT,),
        write_roots=(),
    ),
}

# This project's Cartman access is confined to chartat1's private tree.
ZONES: dict[str, Zone] = {
    name: _SOURCE_ZONES[name] for name in ("repo", "sandbox")
}

def get_zone(name: str) -> Zone:
    zone = ZONES.get(name)
    if zone is None:
        raise ValueError(f"Unknown zone: {name}")
    return zone


def _normalize(path: str) -> str:
    if not path:
        raise ValueError("Path must be non-empty")
    if "\x00" in path:
        raise ValueError("NUL bytes are not allowed in paths")
    posix = PurePosixPath(path)
    if ".." in posix.parts:
        raise ValueError("Parent-directory traversal is not allowed in paths")
    if posix.is_absolute():
        return str(posix)
    return str(posix)


def _is_under(path: str, root: str) -> bool:
    root_path = PurePosixPath(root)
    path_obj = PurePosixPath(path)
    try:
        path_obj.relative_to(root_path)
        return True
    except ValueError:
        return False


def _resolve_relative(zone: Zone, path: str) -> str:
    if PurePosixPath(path).is_absolute():
        return path
    return str(PurePosixPath(zone.read_roots[0]) / path)


def resolve_read_path(zone_name: str, path: str) -> str:
    zone = get_zone(zone_name)
    candidate = _resolve_relative(zone, _normalize(path))
    if any(_is_under(candidate, root) for root in zone.read_roots):
        return candidate
    roots = ", ".join(zone.read_roots)
    raise ValueError(f"Read path {candidate} is outside allowed roots for zone {zone_name}: {roots}")


def resolve_write_path(zone_name: str, path: str) -> str:
    zone = get_zone(zone_name)
    if not zone.write_roots:
        raise ValueError(f"Zone {zone_name} is read-only")
    candidate = _resolve_relative(zone, _normalize(path))
    if not any(_is_under(candidate, root) for root in zone.write_roots):
        roots = ", ".join(zone.write_roots)
        raise ValueError(f"Write path {candidate} is outside allowed roots for zone {zone_name}: {roots}")
    if not can_write_path(candidate):
        raise ValueError(f"Write path {candidate} is not permitted by policy")
    return candidate


def can_write_path(path: str) -> bool:
    normalized = _normalize(path)
    return any(_is_under(normalized, root) for root in
               (REMOTE_REPO_ROOT, *SANDBOX_ROOTS))


def resolve_bind_mount_source(path: str) -> str:
    candidate = _normalize(path)
    if any(_is_under(candidate, root) for root in ALLOWED_BIND_MOUNT_SOURCES):
        return candidate
    roots = ", ".join(ALLOWED_BIND_MOUNT_SOURCES)
    raise ValueError(f"Bind-mount source {candidate} is outside allowed roots: {roots}")


def resolve_bind_mount_destination(path: str) -> str:
    candidate = _normalize(path)
    if any(_is_under(candidate, root) for root in ALLOWED_BIND_MOUNT_DESTINATIONS):
        return candidate
    roots = ", ".join(ALLOWED_BIND_MOUNT_DESTINATIONS)
    raise ValueError(f"Bind-mount destination {candidate} is outside allowed roots: {roots}")


def resolve_system_file_path(path: str) -> str:
    candidate = _normalize(path)
    if candidate in ALLOWED_SYSTEM_FILE_PATHS:
        return candidate
    allowed = ", ".join(ALLOWED_SYSTEM_FILE_PATHS)
    raise ValueError(f"System file path {candidate} is not in the allowlist: {allowed}")


def validate_systemctl_service(name: str) -> str:
    service = str(name)
    if service in ALLOWED_SYSTEMCTL_SERVICES:
        return service
    allowed = ", ".join(ALLOWED_SYSTEMCTL_SERVICES)
    raise ValueError(f"Service {service} is not in the allowlist: {allowed}")


def validate_systemctl_action(name: str) -> str:
    action = str(name)
    if action in ALLOWED_SYSTEMCTL_ACTIONS:
        return action
    allowed = ", ".join(ALLOWED_SYSTEMCTL_ACTIONS)
    raise ValueError(f"systemctl action {action} is not in the allowlist: {allowed}")


def policy_summary_markdown() -> str:
    lines = [
        "# Cartman MCP Policy",
        "",
        f"- Remote host: `{REMOTE_HOST}`",
        f"- Private root: `{PRIVATE_ROOT}` (0700, owned by chartat1)",
        f"- Remote repo root: `{REMOTE_REPO_ROOT}`",
        f"- Remote Python: `{REMOTE_PYTHON}`",
        f"- Remote scheduler root: `{REMOTE_SGE_ROOT}`",
        "- Remote shells use `umask 077`.",
        "",
        "## Zones",
    ]
    for zone in ZONES.values():
        lines.append(f"- `{zone.name}`: {zone.description}")
        for root in zone.read_roots:
            lines.append(f"  read: `{root}`")
        if zone.write_roots:
            for root in zone.write_roots:
                lines.append(f"  write: `{root}`")
        else:
            lines.append("  write: denied")
    lines.extend(("", "Only the repo and run zones are exposed. Public, website, archive, and admin tools are unavailable."))
    return "\n".join(lines) + "\n"


def redact_mapping(value):
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            lower = key.lower()
            if any(token in lower for token in ("token", "secret", "password", "passwd", "key")):
                result[key] = "***REDACTED***"
            else:
                result[key] = redact_mapping(item)
        return result
    if isinstance(value, list):
        return [redact_mapping(item) for item in value]
    return value


def dump_redacted_json(data: dict) -> str:
    return json.dumps(redact_mapping(data), indent=2, sort_keys=True) + "\n"
