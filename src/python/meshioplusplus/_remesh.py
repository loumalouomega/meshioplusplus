"""Isotropic / feature-preserving surface remeshing (``meshioplusplus.remesh``).

A dependency-free mesh *operation* (not a file format): replace a surface
mesh's own triangulation with a new, near-uniformly-sized, well-shaped one at
a caller-chosen vertex count, by approximated centroidal Voronoi diagram
(ACVD) clustering. Every other resolution-changing operation --
:func:`refine`, :func:`decimate`, :func:`decimate_volume`,
:func:`subdivide`/:func:`agglomerate`, :func:`smooth` -- works on the input's
own triangulation and so cannot *raise* a badly-shaped mesh's element quality
at every target count; this one can, because output quality is a property of
the clustering rather than of the input.

**Attribution.** The isotropic clustering engine is derived from
`pyacvd <https://github.com/pyvista/pyacvd>`_ (``src/clustering.cpp``,
nanobind), MIT License, Copyright (c) 2017-2024 The PyVista Developers --
itself an independent implementation of the published research of
S. Valette and J.-M. Chassery, not of the CeCILL-B licensed ACVD project,
which this codebase never reads code from and never vendors (see
``doc/remesh.md``). The optional ``metric="quadric"`` (feature-preserving)
mode has no pyacvd counterpart and is an original synthesis built on this
project's own pre-existing Garland-Heckbert quadric machinery
(``detail/decimate_common.hpp``, shared with :func:`decimate`/
:func:`decimate_volume`).

**This operation is C++-core only, with no numpy fallback at all** -- like
:func:`subdivide`/:func:`agglomerate`/:func:`decimate_volume`/
:func:`conservative_interpolate`. The energy-minimisation loop branches on a
comparison between two candidate cluster configurations, and an accepted
move immediately mutates the accumulators the next comparison reads, so two
independently written implementations disagreeing on a single near-tie
diverge into a genuinely different clustering rather than a last-ulp
difference. Rather than ship a second implementation that could silently
disagree, this module raises when ``_core`` cannot be used, for any input.

The output is a brand-new mesh with new points and new connectivity, so --
unlike every other operation in this layer -- there is no meaningful point or
cell map: ``point_data``, ``cell_data`` and named regions (``point_sets``/
``cell_sets``) are dropped. Transfer a field onto the result with
:func:`interpolate` or :func:`conservative_interpolate`. ``field_data``
passes through verbatim.

Public API:

* :func:`remesh` -- remesh a surface mesh.
"""

from __future__ import annotations

__all__ = ["remesh"]


def remesh(
    mesh,
    num_clusters,
    subdivide=None,
    subsample_ratio=10.0,
    max_subdivide=4,
    max_iterations=100,
    max_repair_passes=10,
    metric="isotropic",
    gradation=0.0,
    preserve_boundary=True,
    return_report=False,
):
    """Remesh a surface uniformly by ACVD clustering.

    :param mesh: the surface mesh to remesh (never modified). ``quad``/
        rectangular ``polygon`` blocks are triangulated first; 3D volume
        cells, higher-order cells and ragged blocks raise, naming the
        operation to use instead (``extract_surface``,
        ``convert_cells(mode="linearize")``).
    :param num_clusters: number of clusters, i.e. the number of vertices the
        output aims for. Must be at least 4 and no more than the subdivided
        input's vertex count.
    :param subdivide: uniform :func:`refine` passes applied before
        clustering (each multiplies the triangle count by four). ``None``
        (the default) picks the smallest count reaching ``subsample_ratio``
        items per cluster, capped at ``max_subdivide``; ``0`` disables
        subdivision.
    :param subsample_ratio: items per cluster targeted by automatic
        subdivision.
    :param max_subdivide: ceiling on automatically chosen subdivision.
        Ignored when ``subdivide`` is given explicitly.
    :param max_iterations: maximum energy-minimisation sweeps per pass (the
        loop also stops early once a sweep changes nothing).
    :param max_repair_passes: how many times to alternate "split
        disconnected clusters, minimise again" -- the fix for clusters the
        energy sweep splits into disjoint pieces, which is where
        non-manifold output comes from. ``0`` skips repair entirely.
    :param metric: ``"isotropic"`` (area-weighted centroidal distance, fast
        and robust, rounds sharp features) or ``"quadric"`` (Garland-Heckbert
        quadric error, preserves sharp edges/corners at extra cost per
        candidate move).
    :param gradation: curvature-gradation exponent ``gamma`` in the item
        weight ``area * kappa**gamma``. ``0.0`` (the default) disables
        gradation entirely -- curvature is not even computed -- and
        reproduces plain area weighting byte-for-byte. Positive values
        concentrate clusters where the surface bends more sharply; ``kappa``
        comes from a local osculating-paraboloid fit over each vertex's
        1-ring.
    :param preserve_boundary: detect the input's open boundary (if any) and
        seed it before the interior, then emit a ``line`` dual cell along
        every boundary edge whose endpoints land in different clusters. A
        no-op, and so free, on a closed mesh.
    :param return_report: also return the run summary dict (``num_clusters``,
        ``num_iterations``, ``subdivide_applied``, ``num_isolated_clusters``,
        ``num_non_manifold_vertices`` -- both of the last two are non-zero
        when repair could not fully fix the output, for two distinct causes:
        disconnected clusters vs. non-manifold "bowtie" vertices; check
        both rather than assuming).
    :returns: the remeshed mesh, or ``(mesh, report)`` if ``return_report``.
    :raises ValueError: on a cluster count below 4 or above the subdivided
        input's vertex count, on a non-positive limit, and on any block
        outside the surface scope.
    :raises NotImplementedError: when the compiled C++ core
        (``meshioplusplus._core``) is unavailable -- this operation has no
        pure-Python fallback.

    Topology is not preserved -- two surface sheets closer together than a
    cluster can merge, and the genus can change; this is inherent to
    dualising a discrete Voronoi partition, not an implementation limit.
    """
    try:
        from . import _core
    except ImportError:
        _core = None

    if _core is None:
        raise NotImplementedError(
            "meshio++: remesh: this operation is C++-core only and has no "
            "pure-numpy fallback (an accepted move in the energy loop "
            "mutates the accumulators the next comparison reads, so a "
            "second implementation could silently diverge into a different "
            "clustering) -- install a build with the compiled "
            "meshioplusplus._core extension"
        )

    res = _core.remesh(
        mesh,
        int(num_clusters),
        -1 if subdivide is None else int(subdivide),
        float(subsample_ratio),
        int(max_subdivide),
        int(max_iterations),
        int(max_repair_passes),
        metric,
        float(gradation),
        bool(preserve_boundary),
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {
        "num_clusters": res["num_clusters"],
        "num_iterations": res["num_iterations"],
        "subdivide_applied": res["subdivide_applied"],
        "num_isolated_clusters": res["num_isolated_clusters"],
        "num_non_manifold_vertices": res["num_non_manifold_vertices"],
    }
    return out, report
