"""Polyhedral coarsening: merge groups of cells into single larger polyhedral
cells.

A dependency-free mesh *operation* (not a file format). ``decimate`` raises by
name on a polyhedron, pointing at ``convert_cells(mode="simplexify")`` -- its
fixed-template QEM edge collapse has no analogue for merging arbitrary
polyhedral cells. ``agglomerate`` is a genuinely different algorithm: greedy
seed-and-grow over the mesh's shared-face dual, absorbing face-adjacent
neighbours into a group until it reaches a target size, then emitting one
polyhedron per group whose faces are exactly that group's *external*
boundary -- every face shared by two members of the same group cancels out of
the result, by construction, not by a tolerance.

**This operation is C++-core only, with no numpy fallback at all** -- the same
reasoning ``_subdivide.py`` (its one-to-many sibling) already documents: the
emit step depends transitively on a winding repair (a discrete branch on the
sign of an enclosed volume) that a second, independently written
implementation could disagree with near-degenerate cells.

Public API:

* :func:`agglomerate` -- polyhedrally coarsen a mesh.
"""

from __future__ import annotations

# Kept for parity with every other operation shim's import shape, even though
# nothing here actually needs numpy -- there is no pure-Python arithmetic path.
import numpy as np  # noqa: F401

__all__ = ["agglomerate"]


def agglomerate(
    mesh,
    target_group_size: int = 8,
    merge_coplanar_faces: bool = False,
    coplanar_angle: float = 1.0,
    min_sphericity: float = 0.0,
    return_report: bool = False,
):
    """Polyhedrally coarsen a mesh: merge groups of cells into single larger
    polyhedral cells.

    Greedy seed-and-grow over the mesh's shared-face dual: cells are seeded
    in ascending order and a group absorbs its unclaimed face-neighbour with
    the largest accumulated shared-face area until it reaches
    ``target_group_size``, or no unclaimed neighbour remains (a short group
    at a mesh boundary or pocket is expected, not an error).
    ``target_group_size=1`` groups every cell by itself -- an identity
    transform in everything but representation.

    Non-volume blocks (2D/1D boundary markers, and any 3D block with no face
    table) pass through unchanged. Points are never pruned or renumbered --
    a group can leave an interior node unreferenced; :func:`clean` with
    ``remove_orphans=True`` is the documented follow-up for a caller who
    wants a minimal point set. Point and Cell regions (and so
    ``point_sets``/``cell_sets``) survive; named **Side** regions do not, a
    many-to-one collapse having no facet correspondence to preserve.

    :param mesh: the mesh to coarsen (never modified).
    :param target_group_size: approximate member cells per output group;
        must be at least 1.
    :param merge_coplanar_faces: fuse the coplanar faces two groups (or a group
        and the boundary) share into single polygons, identically on both
        sides, when their outline is one simple loop and every touched
        polyhedron stays closed.
    :param coplanar_angle: the largest angle, in degrees, between face normals
        still treated as coplanar; in ``[0, 90)``.
    :param min_sphericity: refuse to absorb a cell when the union would be
        less round than this (``pi^(1/3) (6V)^(2/3) / A``: 1 for a ball, about
        0.81 for a cube); 0 disables the gate.
    :param return_report: also return ``{"cell_map", "num_faces_merged",
        "num_rejected"}``.
    :returns: the coarsened mesh; with ``return_report``, ``(mesh, report)``.
    :raises ValueError: when ``target_group_size`` is 0, or when the mesh
        contains a face shared by three or more cells (non-manifold) -- the
        owner/neighbour classification the merge relies on is only
        well-defined on a 2-manifold face.
    :raises NotImplementedError: when the compiled C++ core
        (``meshioplusplus._core``) is unavailable -- this operation has no
        pure-Python fallback.
    """
    try:
        from . import _core
    except ImportError:
        _core = None

    if _core is None:
        raise NotImplementedError(
            "meshio++: agglomerate: this operation is C++-core only and has "
            "no pure-numpy fallback (the emit step depends on a winding "
            "repair, a discrete branch a second implementation could "
            "disagree with near-degenerate cells) -- install a build with "
            "the compiled meshioplusplus._core extension"
        )

    res = _core.agglomerate(
        mesh,
        int(target_group_size),
        bool(merge_coplanar_faces),
        float(coplanar_angle),
        float(min_sphericity),
    )
    if return_report:
        return res["mesh"], {
            "cell_map": res["cell_map"],
            "num_faces_merged": int(res["num_faces_merged"]),
            "num_rejected": int(res["num_rejected"]),
        }
    return res["mesh"]
