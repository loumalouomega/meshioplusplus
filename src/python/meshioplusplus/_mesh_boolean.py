"""Watertight surface mesh boolean operations (roadmap 9.1).

Dispatches ``"union"``, ``"difference"`` and ``"intersection"`` on two closed
triangle meshes through the C++ core's registered :cpp:class:`BooleanProvider`
(for example the optional Apache-2.0 Manifold provider).
"""

from __future__ import annotations

__all__ = ["mesh_boolean", "has_boolean_provider"]


def has_boolean_provider() -> bool:
    """Return ``True`` when a surface boolean provider is registered in the C++ core."""
    from . import _core

    return bool(_core.has_boolean_provider())


def mesh_boolean(
    a,
    b,
    op: str = "union",
    *,
    backend: str = "",
    return_report: bool = False,
):
    """Compute a watertight surface boolean (``"union"``, ``"difference"`` or ``"intersection"``).

    :returns: the resulting triangle :class:`~meshioplusplus.Mesh`, or
        ``(mesh, {"backend": str, "num_triangles": int})`` when ``return_report=True``.
    """
    from . import _core

    res = _core.mesh_boolean(a, b, str(op), str(backend))
    out = res.pop("mesh")
    return (out, res) if return_report else out
