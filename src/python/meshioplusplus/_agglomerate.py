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


def agglomerate(mesh, target_group_size: int = 8):
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
    :returns: the coarsened mesh.
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

    res = _core.agglomerate(mesh, int(target_group_size))
    return res["mesh"]
