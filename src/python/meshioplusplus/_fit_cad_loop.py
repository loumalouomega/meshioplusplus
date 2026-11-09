"""Closing the reverse-engineering loop (roadmap 7.7.1-7.7.3): the scan's
support of each face of a fitted solid, a driver that refits towards a target,
and primitive constraints under control. C++-core only (see
``doc/fit_cad_loop.md`` and ``doc/constrain_primitives.md``)."""

from __future__ import annotations

__all__ = ["face_support", "fit_cad_iterate", "constrain_primitives"]


def face_support(
    solid,
    scan,
    *,
    tolerance: float = 0.0,
    measured_share: float = 0.9,
    inferred_share: float = 0.1,
    hole_data: str = "repair:hole",
    return_report: bool = False,
):
    """Class each face of a ``fit_cad`` re-mesh by how much of it the scan
    supports: ``measured`` (at least ``measured_share`` of its area within
    ``tolerance`` of the scan), ``inferred`` (below ``inferred_share``: no
    scan near) or ``partial``. Scan cells whose ``hole_data`` is at least 0
    (``repair``'s fills) count as unmeasured.

    :returns: the solid with ``fit_cad:face_support`` (0 measured, 1 partial,
        2 inferred), ``fit_cad:face_measured`` (field data, per face) and
        ``fit_cad:support`` (cell data); with ``return_report`` ``(solid,
        report)``: ``face_share``, ``face_class``, ``num_measured``,
        ``num_partial``, ``num_inferred``, ``share_within`` (the scan's points
        within tolerance of the solid) and ``tolerance``.
    :raises ValueError: on a solid without ``fit_cad:face``.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    res = _core_op("face_support")(
        solid,
        scan,
        float(tolerance),
        float(measured_share),
        float(inferred_share),
        str(hole_data),
    )
    out = res.pop("mesh")
    if not return_report:
        return out
    res["face_share"] = [float(v) for v in res["face_share"]]
    res["face_class"] = [int(v) for v in res["face_class"]]
    return out, res


def fit_cad_iterate(
    scan,
    step_path: str | None = None,
    *,
    segment_tolerance: float = 0.0,
    fit_tolerance: float = 0.0,
    target_within: float = 0.0,
    target_explained: float = 0.0,
    part_class: str = "any",
    max_rounds: int = 4,
    worst_faces: int = 3,
    refine_factor: float = 0.25,
    min_patch_cells: int = 1,
    backend: str = "auto",
    return_report: bool = False,
):
    """Segment (``segment``'s primitive method), fit (``fit_primitives``) and
    build (``fit_cad``) a closed scan in rounds: each round measures the
    share of the scan within tolerance of the solid (and, with
    ``target_explained``, the share of the skin's area analytic faces
    explain), stops at the target, at no improvement or after ``max_rounds``,
    and otherwise re-segments the patches of the worst faces finer by
    ``refine_factor``. ``part_class`` (``"analytic"``, ``"cast"``,
    ``"freeform"``) sets the target (0.90, 0.85, 0.75) and the tolerances
    (``1e-3``, ``2e-3``, ``4e-3`` of the diagonal).

    :returns: the best round's solid with ``face_support``'s arrays; with
        ``return_report`` ``(solid, report)``: ``skin`` (its fitted skin),
        ``rounds`` (each ``share_within``, ``explained``, ``hausdorff_max``,
        ``segment_tolerance``, ``num_patches``, ``num_faces``,
        ``resegmented``), ``best_round``, ``target_within``, ``target_met``,
        ``stop_reason``, the support counts, ``valid`` and ``num_faces``.
    :raises ValueError: on an option out of range.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    res = _core_op("fit_cad_iterate")(
        scan,
        str(step_path or ""),
        float(segment_tolerance),
        float(fit_tolerance),
        float(target_within),
        float(target_explained),
        str(part_class),
        int(max_rounds),
        int(worst_faces),
        float(refine_factor),
        int(min_patch_cells),
        str(backend),
    )
    out = res.pop("mesh")
    return (out, res) if return_report else out


def _relations(items, what):
    out = []
    for item in items or ():
        a, b, kind = item
        out.append((int(a), int(b), str(kind)))
    return out


def constrain_primitives(
    mesh,
    *,
    face_tolerance: float = 0.0,
    angle_tolerance: float = -1.0,
    distance_tolerance: float = -1.0,
    min_fraction: float = -1.0,
    force=None,
    forbid=None,
    return_report: bool = False,
):
    """Fit related primitives of a ``fit_primitives`` skin together: coaxial
    groups, parallel and perpendicular planes and equal radii, as
    ``fit_cad(constrain_primitives=True)`` does, with ``force`` and ``forbid``
    lists of ``(patch_a, patch_b, relation)`` (``"coaxial"``, ``"parallel"``,
    ``"perpendicular"``, ``"equal_radius"``; patches by ``segment:patch``
    label) and angle (degrees) and distance tolerances in place of the pass's
    own (negative keeps them).

    :returns: the skin with its ``fit:patch_*`` table rewritten; with
        ``return_report`` ``(mesh, report)``: ``num_constrained``,
        ``num_groups``, ``face_tolerance``, ``relations`` (each ``patch_a``,
        ``patch_b``, ``relation``, ``angle_deg``, ``distance``, ``forced``)
        and ``unmet`` (forced relations not kept).
    :raises ValueError: on a mesh without the table, an unknown patch or
        relation.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    res = _core_op("constrain_primitives")(
        mesh,
        float(face_tolerance),
        float(angle_tolerance),
        float(distance_tolerance),
        float(min_fraction),
        _relations(force, "force"),
        _relations(forbid, "forbid"),
    )
    out = res.pop("mesh")
    if not return_report:
        return out
    res["unmet"] = [tuple(u) for u in res["unmet"]]
    return out, res
