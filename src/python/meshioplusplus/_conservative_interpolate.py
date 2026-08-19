"""Mass-preserving cross-mesh field transfer.

``conservative_interpolate(source, target, ...)`` returns a new mesh that is a
copy of the target — geometry, connectivity, and the target's own data
preserved exactly — with the requested source arrays remapped onto it so that,
over the region the two meshes share, ``sum(target value * target measure)``
equals ``sum(source value * source measure)``. This is the property
:func:`interpolate`'s ``"barycentric"`` mode does not have, since a pointwise
sampler has no notion of "how much of the source region a target sample
stands for" — CFD/FEM solver coupling and successive-remeshing workflows need
exactly this guarantee.

Both meshes are reduced to triangles (2D) or tetrahedra (3D) via
``convert_cells(mode="simplexify")`` — which accepts ragged and polyhedron
blocks for free, since Simplexify already fans them into simplices via a
shipped, tested path — and every overlapping simplex pair is measured with an
exact geometric clip; a target cell's value is the overlap-measure-weighted
mean of every source cell it intersects. ``cell_data`` is transferred by this
algorithm directly; ``point_data`` is transferred by **composition**
(``point_data_to_cell_data`` lumps the source array onto a cell proxy, the
same overlap algorithm remaps it, ``cell_data_to_point_data`` distributes the
result back onto the target's points) — a layered approximation, not exact
nodal/FEM conservation.

This operation reports no built-in integral/conservation diagnostic — use the
``data`` module's field-integration reductions to check how well conservation
held on a given mesh.

**This operation is C++-core only, with no pure-numpy fallback at all** — like
``subdivide``/``agglomerate``/``decimate_volume``. The 3D clip kernel is a
discrete-branch geometric algorithm (half-space in/out classification, cutting
-plane chord deduplication, angle-sorted cap triangulation), and a second,
independently written implementation of those branches could disagree with
the first near a degenerate or tangent overlap.

Public API:

* :func:`conservative_interpolate` -- mass-preservingly remap data between
  two meshes.
"""

from __future__ import annotations

# Kept for parity with every other operation shim's import shape, even though
# nothing here actually needs numpy -- there is no pure-Python arithmetic path.
import numpy as np  # noqa: F401

__all__ = ["conservative_interpolate"]


def conservative_interpolate(
    source, target, arrays=None, default_value=0.0, on_conflict="error"
):
    """Mass-preservingly (overlap-measure weighted) transfer data arrays from
    ``source`` onto ``target``.

    Parameters
    ----------
    source :
        the mesh whose data is sampled (unmodified).
    target :
        the mesh receiving the samples (unmodified); its geometry,
        connectivity, own data and sets are preserved exactly.
    arrays :
        source array names to transfer, or ``None`` (default) for every
        source ``point_data`` **and** ``cell_data`` array — there is one
        algorithm regardless of location, unlike :func:`interpolate`, so there
        is no "cell_data only when named" special case. An unknown name
        raises.
    default_value :
        the fill value (every component) for a target cell whose covered
        fraction of source overlap is below the coverage tolerance.
    on_conflict :
        ``"error"`` (default), ``"overwrite"``, or ``"suffix"`` (writes to
        ``name + "_interp"``) when a transferred name already exists on the
        target.

    Returns
    -------
    Mesh
        a copy of the target with the requested source arrays attached
        (always Float64 — a measure-weighted mean is not integral). The
        target's ``point_sets``/``cell_sets`` ride through untouched; source
        sets are **not** transferred; ``mesh.info``/``gmsh_periodic`` are not
        carried.

    Raises
    ------
    ValueError
        on an empty source, mismatched maximum topological dimensions between
        the two meshes, an unknown array name, or a name conflict.
    NotImplementedError
        when the compiled C++ core (``meshioplusplus._core``) is unavailable
        -- this operation has no pure-Python fallback.
    """
    try:
        from . import _core
    except ImportError:
        _core = None

    if _core is None:
        raise NotImplementedError(
            "meshio++: conservative_interpolate: this operation is C++-core "
            "only and has no pure-numpy fallback (the 3D clip kernel is a "
            "discrete-branch geometric algorithm a second implementation "
            "could disagree with near-degenerate or tangent overlaps) -- "
            "install a build with the compiled meshioplusplus._core extension"
        )

    out = _core.conservative_interpolate(
        source,
        target,
        None if arrays is None else [str(a) for a in arrays],
        float(default_value),
        str(on_conflict),
    )

    # Nothing is renumbered, so the target's sets pass through verbatim (they
    # never reach the C++ core); source sets are deliberately not transferred.
    out.point_sets = {k: v for k, v in target.point_sets.items()}
    out.cell_sets = {k: list(v) for k, v in target.cell_sets.items()}
    return out
