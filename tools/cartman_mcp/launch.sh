#!/bin/sh
set -eu
umask 077

ROOT="$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)"
SERVER="${ROOT}/tools/cartman_mcp/server.py"
PYTHON_BIN="${CARTMAN_MCP_PYTHON:-python3}"

SHARED_RUNTIME_ROOT="${ROOT}/.cache/cartman_mcp"
export RAYTRACING_CARTMAN_MCP_LOG_DIR="${SHARED_RUNTIME_ROOT}"
LOG="${SHARED_RUNTIME_ROOT}/wrapper.log"

mkdir -p -m 0700 "${SHARED_RUNTIME_ROOT}"
{
  printf '%s kind=wrapper-start pid=%s ppid=%s cwd=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$$" "${PPID:-}" "$(pwd)"
} >>"${LOG}"

export PYTHONUNBUFFERED=1
exec "${PYTHON_BIN}" "${SERVER}" "$@" 2>>"${LOG}"
