"""The outline of a part seen along a direction (roadmap 7.6.3).

The surface projected onto a plane and unioned exactly, as closed chains with
its holes, optionally offset by a clearance or made a prism (a nest, a jig, a
foam insert). C++-core only (see ``doc/silhouette.md``).
"""

from __future__ import annotations

__all__ = ["silhouette"]

_PREFIX = "meshio++: silhouette: "


def silhouette(
    mesh,
    normal=(0.0, 0.0, 1.0),
    *,
    clearance: float = 0.0,
    round_angle: float = 180.0,
    arc_segments: int = 64,
    extrude: bool = False,
    region: str | None = None,
    merge_tolerance: float = 1e-9,
    return_report: bool = False,
):
    """The outline of ``mesh`` seen along ``normal``.

    :param mesh: a surface (triangles, quads, polygons) or a volume mesh,
        outlined through its skin.
    :param normal: the viewing direction; the outline lies in the plane normal
        to it at the part's lowest point along it.
    :param clearance: > 0 grows the outline by this distance (a loose fit),
        < 0 shrinks it (a press fit).
    :param round_angle: a corner the offset opens is rounded with an arc when
        it turns by less than this many degrees, else kept sharp (mitred).
    :param arc_segments: polygon sides a full turn of arc takes.
    :param extrude: return the closed triangle surface of the prism over the
        outline, from the part's lowest to its highest point, instead.
    :param region: outline only this cell region of a surface.
    :param merge_tolerance: projected points closer than this times the
        projected box's diagonal are one.
    :param return_report: also return the report.
    :returns: a chained line mesh (``section:chain``, ``section:closed``,
        ``section:offset`` and ``silhouette:hole`` cell data, ready for
        ``fit_curves``), or the prism; both with ``silhouette:normal`` and
        ``silhouette:extent`` field data. With ``return_report``
        ``(mesh, report)``: ``num_loops``, ``num_outer``, ``num_holes``,
        ``num_triangles``, ``area``, ``perimeter``, ``min_height``,
        ``max_height`` and ``warnings``.
    :raises ValueError: on a bad option or a mesh with no surface cells.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    n = tuple(float(v) for v in normal)
    if len(n) != 3:
        raise ValueError(f"{_PREFIX}normal must have three components")
    res = _core_op("silhouette")(
        mesh,
        n,
        float(clearance),
        float(round_angle),
        int(arc_segments),
        bool(extrude),
        region or "",
        float(merge_tolerance),
    )
    out = res.pop("mesh")
    if not return_report:
        return out
    report = {
        k: int(res[k]) for k in ("num_loops", "num_outer", "num_holes", "num_triangles")
    }
    for k in ("area", "perimeter", "min_height", "max_height"):
        report[k] = float(res[k])
    report["warnings"] = [str(w) for w in res["warnings"]]
    return out, report
