"""Volume decimation (``meshioplusplus.decimate_volume``).

A dependency-free mesh *operation* (not a file format): reduce a tetrahedral
mesh's cell count by greedy quadric-error-metric tet-edge collapse. The
volume-mesh sibling of surface :func:`decimate` -- a **separate operation**,
not a mode on it: ``decimate`` keeps raising by name on any 3D volume block,
pointing here, and this module never touches ``decimate``'s own machinery.

Boundary vertices **participate** in decimation with a real quadric-error
objective (``preserve_boundary`` defaults to ``False`` here, the opposite of
``decimate``'s own default) -- pass ``preserve_boundary=True`` to reproduce
the simpler pinned-boundary behaviour.

**This operation is C++-core only, with no numpy fallback at all** -- like
:func:`subdivide`/:func:`agglomerate` and unlike surface :func:`decimate`
(whose own guards are safely twinned only because every floating-point
expression feeding them is transcribed token for token). The volume-specific
validity guards here (the vertex-link set-equality check, the duplicate-tet
guard, and the tet-inversion check) sit on top of the mesh's own outer skin,
whose incremental maintenance through hundreds of collapses is exactly the
kind of intricate bookkeeping a second, independently written implementation
would risk silently diverging from -- the same reasoning ``_subdivide.py``
and ``_agglomerate.py`` already document for their own winding repairs.
Rather than ship a second, subtly different implementation of that
bookkeeping, this module raises when ``_core`` cannot be used, for any input.

Named regions (and so ``point_sets``/``cell_sets``) are carried by the C++
core itself, exactly as surface ``decimate`` carries them on its own C++
path -- no shim-side remapping is needed here at all.

Public API:

* :func:`decimate_volume` -- decimate a tetrahedral mesh.
"""

from __future__ import annotations

import numpy as np

__all__ = ["decimate_volume"]


def decimate_volume(
    mesh,
    ratio=None,
    target_cells=None,
    max_error=None,
    placement="optimal",
    preserve_boundary=False,
    preserve_features=True,
    feature_angle=30.0,
    frozen=None,
    return_report=False,
):
    """Decimate a tetrahedral mesh by greedy quadric-error-metric tet-edge
    collapse.

    Exactly one of ``ratio``, ``target_cells``, ``max_error`` must be given.
    Tet-only: any hexahedron/wedge/pyramid/polyhedron/ragged/higher-order 3D
    block, or any non-3D block mixed in, raises -- run
    :func:`convert_cells` with ``mode="simplexify"`` first. The block
    structure stays 1:1 with the input.

    Every vertex accumulates a Garland-Heckbert plane quadric from its
    incident **boundary triangles only** (the mesh's own outer skin) -- a
    purely interior vertex's quadric is exactly zero, so interior-only edges
    are scored by squared length instead and always rank behind
    boundary-touching ones in the greedy queue. ``max_error`` is only really
    meaningful for boundary-touching (real quadric-error) collapses;
    ``ratio``/``target_cells`` remain exactly well-defined either way.

    :param mesh: the tet mesh to decimate (never modified).
    :param ratio: fraction of the tets to KEEP, in ``(0, 1]``.
    :param target_cells: absolute tet count to stop at. A collapse removes
        one or more tets (every tet sharing the collapsed edge), so the
        result lands within one collapse of the target.
    :param max_error: collapse only while the cheapest boundary-touching
        candidate's quadric error is at most this (squared mesh units).
    :param placement: where the surviving vertex goes -- ``"optimal"`` (the
        quadric minimizer, midpoint when ill-conditioned or purely
        interior), ``"midpoint"``, or ``"endpoint"`` (the endpoint with the
        lower error).
    :param preserve_boundary: pin every boundary vertex outright (the
        once-used-face test on the mesh's own outer skin), reproducing
        :func:`decimate`'s own default instead of letting boundary vertices
        participate.
    :param preserve_features: pin boundary vertices whose incident
        boundary-triangle normals differ by more than ``feature_angle``, so
        corners and creases of the outer surface survive.
    :param feature_angle: that angle, in degrees, between boundary-triangle
        normals.
    :param frozen: optional extra vertex ids to pin, as an index array or the
        name of one of ``mesh.point_sets``.
    :param return_report: also return the run summary dict
        (``tets_removed``, ``points_removed``, ``collapses_rejected``,
        ``max_error_applied``).
    :returns: the decimated mesh, or ``(mesh, report)`` if ``return_report``.
    :raises ValueError: when not exactly one stopping criterion is given, on
        a mis-sized frozen mask, on a non-manifold boundary face, and on any
        block outside the tet-only scope.
    :raises NotImplementedError: when the compiled C++ core
        (``meshioplusplus._core``) is unavailable -- this operation has no
        pure-Python fallback.

    ``point_data`` at a surviving vertex is blended between the collapse's
    two endpoints (clamped edge parameter) -- except **integer** arrays,
    which keep the survivor's own value. Each surviving tet keeps its own
    ``cell_data`` row; ``field_data`` passes through verbatim.
    """
    try:
        from . import _core
    except ImportError:
        _core = None

    if _core is None:
        raise NotImplementedError(
            "meshio++: decimate_volume: this operation is C++-core only and "
            "has no pure-numpy fallback (the validity guards depend on "
            "incrementally maintaining the mesh's own outer skin, intricate "
            "bookkeeping a second implementation could disagree with) -- "
            "install a build with the compiled meshioplusplus._core extension"
        )

    if isinstance(frozen, str):
        try:
            frozen = np.asarray(mesh.point_sets[frozen], dtype=np.int64)
        except KeyError:
            raise ValueError(
                f"meshio++: decimate_volume: no point_set named {frozen!r} "
                f"(available: {sorted(mesh.point_sets)})"
            ) from None

    target_ratio = -1.0 if ratio is None else float(ratio)
    cells = -1 if target_cells is None else int(target_cells)
    error = -1.0 if max_error is None else float(max_error)

    res = _core.decimate_volume(
        mesh,
        target_ratio,
        cells,
        error,
        placement,
        bool(preserve_boundary),
        bool(preserve_features),
        float(feature_angle),
        None if frozen is None else np.asarray(frozen, dtype=np.int64),
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {
        "tets_removed": res["tets_removed"],
        "points_removed": res["points_removed"],
        "collapses_rejected": res["collapses_rejected"],
        "max_error_applied": res["max_error_applied"],
    }
    return out, report
