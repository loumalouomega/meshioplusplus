"""Polyhedral refinement: one polyhedral child per face, connected to a new
interior point.

A dependency-free mesh *operation* (not a file format). ``refine`` and
``decimate`` both raise by name on a polyhedron, pointing at
``convert_cells(mode="simplexify")`` -- both are built on fixed same-type
subdivision templates, and an arbitrary polyhedron has none. ``subdivide`` is
the answer for polyhedral refinement specifically: it goes through the same
uniform face-ring machinery ``gradient``/``compute_quality`` already use, so it
needs no per-type table at all -- the same code handles every 3D cell type the
mesh already supports (tabulated, reduced to corners for a quadratic variant,
or an existing polyhedron block).

**This operation is C++-core only, with no numpy fallback at all** -- unlike
most operations in this repo, which fall back to a pure-numpy reference when
``_core`` is unavailable. The winding repair (``orient_rings``) is a discrete
branch on the sign of an enclosed volume, and a second, independently written
implementation of that branch could land on the opposite side for a
near-degenerate cell and then diverge macroscopically rather than by
round-off -- the same reasoning ``_convert_cells.py``'s polyhedron branch and
``_smooth.py``'s inversion guard already document. Rather than ship a second,
subtly different winding repair, this module raises when ``_core`` cannot be
used, for any input, tabulated or not.

Public API:

* :func:`subdivide` -- polyhedrally refine a mesh.
"""

from __future__ import annotations

# Kept for parity with every other operation shim's import shape, even though
# nothing here actually needs numpy -- there is no pure-Python arithmetic path.
import numpy as np  # noqa: F401

__all__ = ["subdivide"]


def subdivide(mesh, record_parent_ids: bool = False):
    """Polyhedrally refine a mesh: one polyhedral child per face of every
    eligible 3D cell, connected to a new interior point.

    Automatically conforming -- each child uses exactly one *original* face
    untouched, so a neighbouring cell across it still sees the same face on the
    far side. No closure, no 2:1 balance, no hanging-node bookkeeping is
    needed, unlike :func:`refine`.

    Non-3D blocks (2D/1D boundary markers) and 3D blocks with no face table
    (the full-Lagrange family) pass through unchanged. Point and Cell regions
    (and so ``point_sets``/``cell_sets``) survive; named **Side** regions do
    not, since a child's facets have no correspondence with the parent's --
    the same limitation ``convert_cells(mode="simplexify")`` already has.

    :param mesh: the mesh to subdivide (never modified).
    :param record_parent_ids: when true, attach an Int64
        ``subdivide:parent_cell`` cell_data array recording, per output cell,
        the index of the input cell it came from within its own block.
    :returns: the subdivided mesh.
    :raises ValueError: when a cell's faces are not a closed orientable
        surface (naming the cell and block), so it bounds no volume.
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
            "meshio++: subdivide: this operation is C++-core only and has no "
            "pure-numpy fallback (the winding repair is a discrete branch a "
            "second implementation could disagree with near-degenerate "
            "cells) -- install a build with the compiled meshioplusplus._core "
            "extension"
        )

    res = _core.subdivide(mesh, bool(record_parent_ids))
    return res["mesh"]
