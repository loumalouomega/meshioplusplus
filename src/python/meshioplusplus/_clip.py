"""Cut a surface by a plane, keep a side and cap it (roadmap 7.8.3). C++-core
only (see ``doc/clip.md``)."""

from __future__ import annotations

__all__ = ["clip"]

_PREFIX = "meshio++: clip: "


def clip(
    mesh, origin, normal, *, keep: str = "below", cap: bool = True, return_report=False
):
    """Clip ``mesh`` by the plane through ``origin`` normal to ``normal``.

    :param keep: ``"below"`` keeps ``(x - origin) . normal < 0``, ``"above"``
        the rest.
    :param cap: triangulate the cut's closed loops so a closed part stays closed.
    :returns: one ``triangle`` block with ``clip:cap`` cell data; with
        ``return_report`` ``(mesh, report)``: ``num_kept``, ``num_split``,
        ``num_cap_triangles``, ``num_loops``, ``num_open_chains`` and ``warnings``.
    :raises ValueError: on a zero normal, a volume mesh or a bad ``keep``.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    if keep not in ("below", "above"):
        raise ValueError(f"{_PREFIX}keep is 'below' or 'above', not {keep!r}")
    o = tuple(float(v) for v in origin)
    n = tuple(float(v) for v in normal)
    if len(o) != 3 or len(n) != 3:
        raise ValueError(f"{_PREFIX}origin and normal have three components")
    res = _core_op("clip")(mesh, o, n, keep == "below", bool(cap))
    out = res.pop("mesh")
    if not return_report:
        return out
    report = {
        k: int(res[k])
        for k in (
            "num_kept",
            "num_split",
            "num_cap_triangles",
            "num_loops",
            "num_open_chains",
        )
    }
    report["warnings"] = [str(w) for w in res["warnings"]]
    return out, report
