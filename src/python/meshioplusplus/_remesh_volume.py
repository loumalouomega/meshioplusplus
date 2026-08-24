"""Volumetric retetrahedralization by isosurface stuffing (``meshioplusplus.remesh_volume``).

A dependency-free mesh *operation* (not a file format): retetrahedralize a
**volume** mesh (or a closed surface) at a caller-chosen resolution, by
isosurface stuffing over a body-centered cubic (BCC) lattice. It is the
volumetric sibling of :func:`remesh` -- the same relationship
:func:`decimate_volume` has to :func:`decimate`. Nothing else in this repo
can *raise* the quality of a tetrahedral mesh at a chosen target resolution:
:func:`refine` subdivides the input's own (possibly bad) cells,
:func:`decimate_volume` can only remove elements, :func:`smooth` (see
``method="odt"``) only moves points. This operation instead discards the
input's own tetrahedra entirely and generates a fresh mesh, so output
element quality is a property of the lattice construction, not of the input.

**Attribution.** Implemented from the PUBLISHED DESCRIPTION of Labelle &
Shewchuk, "Isosurface Stuffing: Fast Tetrahedral Meshes with Good Dihedral
Angles" (SIGGRAPH 2007) only. Neither the paper's own reference
implementation (Stellar, non-commercial) nor TetGen (AGPL-3.0, also by
Shewchuk) is read or vendored here -- the same clean-room posture
``doc/remesh.md`` documents for its ACVD attribution. See
``doc/remesh_volume.md`` for the algorithm and its measured warp/quality
tradeoff.

**This operation is C++-core only, with no numpy fallback at all** -- like
:func:`remesh`/:func:`subdivide`/:func:`agglomerate`/:func:`decimate_volume`.
Warp and cut are discrete branches on a sign near a threshold, so a second,
independently-written implementation could land on the other side of a
near-tie and diverge into a different mesh, not a last-ulp difference.
Rather than ship a second implementation that could silently disagree, this
module raises when ``_core`` cannot be used, for any input.

The output is a brand-new mesh with new points and new connectivity, so --
like :func:`remesh` -- there is no meaningful point or cell map:
``point_data``, ``cell_data`` and named regions (``point_sets``/
``cell_sets``) are dropped. Transfer a field onto the result with
:func:`interpolate` or :func:`conservative_interpolate`. ``field_data``
passes through verbatim.

Public API:

* :func:`remesh_volume` -- retetrahedralize a volume mesh or closed surface.
"""

from __future__ import annotations

__all__ = ["remesh_volume"]


def remesh_volume(
    mesh,
    resolution=None,
    cell_size=None,
    bounds=None,
    padding=0.0,
    padding_relative=0.1,
    max_cells=20000000,
    max_tets=20000000,
    warp_fraction=0.35,
    sign="pseudonormal",
    weight="angle",
    watertight_check="warn",
    surface_region="",
    grid_cell_size=0.0,
    max_winding_work=2.0e9,
    return_report=False,
):
    """Retetrahedralize ``mesh`` at a caller-chosen resolution.

    :param mesh: a volume mesh (its boundary is extracted internally via
        :func:`extract_surface`, which is closed by construction) or a
        closed surface directly (triangle, quad, rectangular polygon --
        exactly :func:`compute_sdf`'s own scope). Never modified.
    :param resolution: cell counts per axis of the root (pre-cut) BCC
        lattice. Exactly one of ``resolution``/``cell_size`` must be set.
    :param cell_size: cubic cell size of the root lattice.
    :param bounds: explicit ``(lo, hi)`` or a flat 6-tuple; ``None`` uses the
        surface's own bounding box.
    :param padding: padding added to every side, in world units.
    :param padding_relative: padding added to every side, as a fraction of
        the bounding-box diagonal. Defaults to 0.1: a lattice that stops
        exactly at the surface leaves nowhere for a fully-outside cell.
    :param max_cells: refuse to generate a root lattice above this many
        cells, checked before any lattice is built. ``<= 0`` lifts the
        limit.
    :param max_tets: refuse an output with more tets than this, checked
        after cutting (the true output count depends on how much of the
        lattice the surface excludes). Unlike ``max_cells``, this has no
        "lifts the limit" value -- the check is unconditional.
    :param warp_fraction: fraction of a lattice vertex's own shortest
        incident edge length within which it may be warped onto the
        surface. ``0`` disables warping (a plain, unwarped cut -- watertight
        and correctly oriented, but with arbitrarily bad boundary dihedral
        angles); see ``doc/remesh_volume.md`` for the measured sweep this
        default was chosen from.
    :param sign, weight, watertight_check, surface_region, grid_cell_size,
        max_winding_work: passed through to the surface-distance query
        every lattice vertex is classified with (:func:`sample_distance`'s
        own vocabulary; ``location``/``band``/``record_closest_cell``/
        ``record_inside`` do not apply here, since this operation queries
        raw points directly rather than through a mesh location).
    :param return_report: also return a run summary dict (``input_quality``
        -- the verdict for the INPUT surface, not the output --
        ``num_tets``, ``num_vertices_warped``, ``num_tets_rejected``,
        ``num_non_manifold_edges``).
    :returns: the output tet mesh, or ``(mesh, report)`` if ``return_report``.
    :raises ValueError: when neither or both of ``resolution``/``cell_size``
        are set, on a non-positive resolution/cell size, on a negative
        ``warp_fraction``, when the root lattice would exceed ``max_cells``,
        when the cut output would exceed ``max_tets``, and on any input the
        surface extraction itself refuses (a higher-order or ragged
        surface block, an empty surface).
    :raises NotImplementedError: when the compiled C++ core
        (``meshioplusplus._core``) is unavailable -- this operation has no
        pure-Python fallback.

    A caller wanting a genuinely high-quality boundary composes with
    :func:`extract_surface` + :func:`remesh` on the result -- this operation
    only ever approximates the boundary at lattice scale. See
    ``mNumNonManifoldEdges``'s own doc comment in ``remesh_volume.hpp`` (and
    ``doc/remesh_volume.md``) for the measured, bounded manifoldness effect
    warping can have on the OUTPUT boundary specifically.
    """
    try:
        from . import _core
    except ImportError:
        _core = None

    if _core is None:
        raise NotImplementedError(
            "meshio++: remesh_volume: this operation is C++-core only and "
            "has no pure-numpy fallback (warp and cut are discrete branches "
            "on a sign near a threshold, so a second implementation could "
            "silently diverge into a different mesh) -- install a build "
            "with the compiled meshioplusplus._core extension"
        )

    res = _core.remesh_volume(
        mesh,
        None if resolution is None else [int(v) for v in resolution],
        None if cell_size is None else float(cell_size),
        None if bounds is None else [float(v) for v in bounds],
        float(padding),
        float(padding_relative),
        int(max_cells),
        int(max_tets),
        float(warp_fraction),
        sign,
        weight,
        watertight_check,
        surface_region,
        float(grid_cell_size),
        float(max_winding_work),
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {
        "input_quality": res["input_quality"],
        "num_tets": res["num_tets"],
        "num_vertices_warped": res["num_vertices_warped"],
        "num_tets_rejected": res["num_tets_rejected"],
        "num_non_manifold_edges": res["num_non_manifold_edges"],
    }
    return out, report
