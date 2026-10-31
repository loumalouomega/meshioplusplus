"""Tetrahedral <-> polyhedral conversion.

A dependency-free mesh *operation* (not a file format):

* :func:`tet_to_poly` builds the **median dual** of a tetrahedral mesh -- the
  tet-to-poly conversion CFD meshers perform. Every vertex becomes one
  polyhedral cell, bounded by the ring of tetrahedron and triangle centroids
  around each edge leaving it and, on the boundary, by the planar pieces of
  the boundary triangles around it (fused while their normals agree within
  ``feature_angle``, so feature edges and corners survive).
* :func:`poly_to_tet` decomposes a polyhedral (or mixed) mesh into
  tetrahedra: ``convert_cells(mode="simplexify")``, except that a polyhedron
  that already is a tetrahedron becomes exactly that one ``tetra``, with no
  new points.

**Both are C++-core only, with no numpy fallback** -- the same reasoning
``_agglomerate.py`` documents: the construction depends on winding and
orientation decisions (discrete branches) that a second, independently
written implementation could resolve differently near degenerate cells.

Public API:

* :func:`tet_to_poly` -- tetrahedral mesh to its polyhedral dual.
* :func:`poly_to_tet` -- polyhedral mesh to tetrahedra.
"""

from __future__ import annotations

__all__ = ["tet_to_poly", "poly_to_tet"]


def _core_or_raise(op):
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None:
        raise NotImplementedError(
            f"meshio++: {op}: this operation is C++-core only and has no "
            "pure-numpy fallback (it depends on winding decisions a second "
            "implementation could resolve differently) -- install a build "
            "with the compiled meshioplusplus._core extension"
        )
    return _core


def tet_to_poly(
    mesh,
    feature_angle: float = 1.0,
    point_data_to_cell_data: bool = False,
    record_vertex_ids: bool = False,
    return_report: bool = False,
):
    """Convert a tetrahedral mesh into its polyhedral (median) dual: one
    polyhedron per vertex.

    Each vertex's cell is bounded by one face per edge leaving it (the ring of
    tetrahedron and triangle centroids around the edge, closed through the
    edge midpoint where the edge reaches the boundary) and, if the vertex is
    on the boundary, by its share of the boundary triangles around it. Those
    pieces are fused into one polygon per run of consecutive triangles whose
    normals differ by at most ``feature_angle`` degrees (and that lie on the
    same side and in the same input triangle block), so feature edges and
    corners are kept. With the default of 1 degree every boundary face is
    planar and the total volume equals the tetrahedral mesh's.

    A vertex on the interface between two ``tetra`` blocks yields one cell per
    block, so material interfaces stay sharp. Each ``tetra`` block becomes one
    ``polyhedron`` block; each ``triangle`` block becomes a ``polygon`` block of
    the dual boundary faces over it; any other block is dropped with a warning.

    The input points keep their indices (interior vertices are left
    unreferenced -- :func:`clean` with ``remove_orphans=True`` compacts them),
    followed by the new edge midpoints, triangle centroids and tetrahedron
    centroids; ``point_data`` is extended to them by averaging. ``cell_data``
    is the volume-weighted mean over a cell's tetrahedra for floating-point
    arrays and the largest tetrahedron's value for integer ones (polygon
    blocks weight by triangle area). Point regions survive; cell and side
    regions come back empty.

    :param mesh: the tetrahedral mesh (never modified).
    :param feature_angle: degrees, in ``[0, 180)``; boundary triangles whose
        normals differ by more start a new dual boundary face.
    :param point_data_to_cell_data: also copy each vertex's ``point_data`` onto
        its dual cells (and boundary polygons) as ``cell_data``; a name already
        used by a ``cell_data`` array is skipped with a warning.
    :param record_vertex_ids: attach an Int64 ``dual:vertex`` ``cell_data``
        array holding each output cell's source vertex.
    :param return_report: also return ``{"vertex_map", "num_new_points",
        "num_boundary_faces", "num_split_vertices"}``.
    :returns: the dual mesh; with ``return_report``, ``(mesh, report)``.
    :raises ValueError: when a volume block is not linear ``tetra``, the mesh
        has no tetrahedra or 3D points, a face is shared by three or more
        tetrahedra, or ``feature_angle`` is out of range.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    _core = _core_or_raise("tet_to_poly")
    res = _core.tet_to_poly(
        mesh,
        float(feature_angle),
        bool(point_data_to_cell_data),
        bool(record_vertex_ids),
    )
    if return_report:
        return res["mesh"], {
            "vertex_map": res["vertex_map"],
            "num_new_points": int(res["num_new_points"]),
            "num_boundary_faces": int(res["num_boundary_faces"]),
            "num_split_vertices": int(res["num_split_vertices"]),
        }
    return res["mesh"]


def poly_to_tet(mesh, record_parent_ids: bool = False, return_report: bool = False):
    """Decompose a polyhedral (or mixed) mesh into tetrahedra.

    ``convert_cells(mode="simplexify")`` with one addition: a polyhedron that
    already is a tetrahedron (four triangular faces over four distinct nodes)
    becomes that one positively oriented ``tetra``, with no new points. Every
    other polyhedron is decomposed by the volume-conserving corner-average fan
    (one tetrahedron per face edge, adding the cell's and each face's corner
    average as points); hexahedra, wedges and pyramids use simplexify's fixed
    templates and polygons a triangle fan.

    ``poly_to_tet(tet_to_poly(mesh))`` is a valid tetrahedral mesh of the same
    volume, not ``mesh`` itself; a tetrahedral mesh stored as ``polyhedron4``
    cells does come back exactly.

    :param mesh: the mesh to decompose (never modified).
    :param record_parent_ids: attach an Int64 ``convert:parent_cell``
        ``cell_data`` array recording each output cell's input cell.
    :param return_report: also return ``{"point_map", "cell_maps",
        "num_exact_tets", "num_fanned_cells"}``.
    :returns: the tetrahedral mesh; with ``return_report``, ``(mesh, report)``.
    :raises ValueError: when a polyhedron's faces are not a closed orientable
        surface.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    _core = _core_or_raise("poly_to_tet")
    res = _core.poly_to_tet(mesh, bool(record_parent_ids))
    if return_report:
        return res["mesh"], {
            "point_map": res["point_map"],
            "cell_maps": list(res["cell_maps"]),
            "num_exact_tets": int(res["num_exact_tets"]),
            "num_fanned_cells": int(res["num_fanned_cells"]),
        }
    return res["mesh"]
