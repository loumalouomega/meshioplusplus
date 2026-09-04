"""MCP server: expose every meshio++ operation to AI agents.

Serves the whole library surface — reading/writing 40+ mesh formats, format
conversion, and every mesh and data operation — over the Model Context
Protocol so agent clients (Claude Code, Claude Desktop, any MCP client) can
drive it. Tools are **stateless and file-path based**: input path(s) in,
output path(s) out, a strict-JSON report back; no mesh lives between calls.

Nothing here is part of the C++ core. Split exactly like
:mod:`meshioplusplus._interop` and for the same reason: the entire behaviour
lives in the pure payload layer :mod:`._tools` (imports only meshioplusplus +
numpy + stdlib, works on every supported Python, tested in the default CI
matrix with no ``mcp`` installed), while the thin registration layer
:mod:`._server` is the only thing importing the ``mcp`` SDK — which needs the
``[mcp]`` extra and Python >= 3.10 (the SDK's floor, not meshio++'s).

Run it with the ``meshioplusplus-mcp`` console script or
``python -m meshioplusplus.mcp`` (both take ``--root DIR`` /
``MESHIOPLUSPLUS_MCP_ROOT`` to sandbox paths). Docs: ``doc/mcp.md``.

``--http`` turns the same process into the browser dataset manager's
companion (a JSON API over the tool registry plus MCP over HTTP); that
front-end lives in :mod:`._http`, the only module importing Starlette and
uvicorn (the ``[dashboard]`` extra). Docs: ``doc/dashboard.md``.
"""

from ._tools import TOOL_REGISTRY, formats_payload, set_root  # noqa: F401

__all__ = [
    "TOOL_REGISTRY",
    "create_server",
    "formats_payload",
    "has_dashboard",
    "has_mcp",
    "main",
    "set_root",
]


def has_mcp() -> bool:
    """Whether the ``mcp`` SDK is importable (the ``[mcp]`` extra)."""
    try:
        import mcp  # noqa: F401

        return True
    except ImportError:
        return False


def has_dashboard() -> bool:
    """Whether the HTTP front-end can run (the ``[dashboard]`` extra:
    the ``mcp`` SDK plus Starlette and uvicorn)."""
    if not has_mcp():
        return False
    try:
        import starlette  # noqa: F401
        import uvicorn  # noqa: F401

        return True
    except ImportError:
        return False


def _require_http():
    if not has_dashboard():
        raise ImportError(
            "meshio++: mcp server: the HTTP front-end needs starlette and "
            "uvicorn; install them with `pip install meshioplusplus[dashboard]` "
            "(requires Python >= 3.10)"
        )
    from . import _http

    return _http


def _require_server():
    if not has_mcp():
        raise ImportError(
            "meshio++: mcp server: mcp is not installed; install it with "
            "`pip install meshioplusplus[mcp]` (requires Python >= 3.10)"
        )
    from . import _server

    return _server


def create_server(root=None):
    """Build the FastMCP server (raises the named error without the extra)."""
    return _require_server().create_server(root=root)


def main(argv=None) -> int:
    """Console entry point (``meshioplusplus-mcp``); runs on stdio."""
    try:
        server_mod = _require_server()
    except ImportError as e:
        import sys

        print(str(e), file=sys.stderr)
        return 1
    return server_mod.main(argv)
