"""ODT remeshing of a tetrahedral mesh (``meshioplusplus.optimize_volume``).

A dependency-free mesh *operation* (not a file format): raise a **tetrahedral**
mesh's worst element quality by *ODT remeshing* -- relocating vertices AND
changing connectivity. This is the genuine "ODT remeshing" the roadmap named,
and the missing third member of a trio whose other two each do half the job:

* :func:`smooth` with ``method="odt"`` is ODT *smoothing* -- it moves each free
  interior tet vertex to the volume-weighted circumcenter average of its
  incident tets, on the mesh's **fixed** connectivity. It cannot fix a badly
  *connected* tetrahedralization (a sliver no vertex motion removes).
* :func:`remesh_volume` *generates* a fresh tet mesh from a signed-distance
  lattice (surface-in, volume-out); it discards the input's tets entirely.

:func:`optimize_volume` alternates the ODT vertex relocation above with
quality-improving topological **flips** (2-3 and 3-2), so both the vertex
positions and the connectivity change -- that is what makes it *remeshing*
rather than *smoothing*.

**Predicate-free (in-posture).** A flip is applied only when, using the
*signed volume* of the candidate tets alone (no in-sphere / Delaunay
predicate), the local configuration is convex and the **minimum** scaled
Jacobian over the new tets strictly beats the minimum over the tets it
replaces. That improvement rule (Freitag & Ollivier-Gooch, "Tetrahedral mesh
improvement using swapping and smoothing", 1997) makes the process monotone in
worst quality, hence terminating -- no Delaunay optimality argument is needed.

**Boundary is invariant by construction.** The flips act only on interior
faces/edges, so with ``preserve_boundary=True`` the output's boundary surface
is byte-identical to the input's: watertight in => watertight out.

**What survives.** The point set is invariant (relocation moves points, flips
only reconnect them), so ``point_data``, ``field_data`` and named **Point**
regions (``point_sets``) carry through unchanged; ``cell_data`` and named
**Cell**/**Side** regions are dropped (a flip has no cell correspondence). The
output is a single ``tetra`` block.

**This operation is C++-core only, with no numpy fallback at all** -- like
:func:`remesh_volume`/:func:`subdivide`/:func:`agglomerate`/
:func:`decimate_volume`. The flip acceptance is a discrete branch on the sign
of a volume and on a near-tie quality comparison; a second implementation could
land on the other side of such a tie and diverge into a different connectivity.
Rather than ship one that could silently disagree, this module raises when
``_core`` cannot be used, for any input.

Public API:

* :func:`optimize_volume` -- ODT-remesh a tetrahedral mesh.
"""

from __future__ import annotations

__all__ = ["optimize_volume"]


def optimize_volume(
    mesh,
    max_iterations=10,
    relocate=True,
    flip=True,
    preserve_boundary=True,
    min_improvement=1e-6,
    frozen=None,
    return_report=False,
):
    """ODT-remesh ``mesh``: relocate vertices and flip connectivity.

    :param mesh: a tetra-only volume mesh (one or more ``tetra`` blocks). Never
        modified.
    :param max_iterations: optimisation sweeps; each sweep is one ODT
        relocation pass followed by one flip pass. The loop stops early once a
        sweep moves no vertex and accepts no flip (a fixed point).
    :param relocate: run the ODT vertex-relocation half of each sweep.
    :param flip: run the topological-flip half of each sweep. With ``flip``
        off and ``relocate`` on, the operation reduces to ODT *smoothing*.
    :param preserve_boundary: pin boundary vertices during relocation (the
        flips never touch the boundary, so with this on the boundary surface is
        exactly preserved). ``False`` lets boundary vertices drift off the
        surface -- there is no surface re-projection here.
    :param min_improvement: the strict quality gain, in scaled-Jacobian units,
        a flip must deliver to be accepted.
    :param frozen: optional pin mask for the relocation half -- either ``None``
        or a length-``len(mesh.points)`` array-like whose non-zero entries pin
        those vertices.
    :param return_report: also return a run-summary dict (``num_flips``,
        ``num_23_flips``, ``num_32_flips``, ``num_vertices_moved``,
        ``num_tets``, ``min_quality_before``, ``min_quality_after``).
    :returns: the optimised mesh, or ``(mesh, report)`` if ``return_report``.
    :raises ValueError: when the mesh contains a non-``tetra`` 3D block (run
        ``convert_cells(mode="simplexify")`` first), a ragged/polyhedron block,
        a non-3D block alongside the tets (drop it via :func:`split`), or no
        tetra block at all.
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
            "meshio++: optimize_volume: this operation is C++-core only and "
            "has no pure-numpy fallback (a flip is accepted on a discrete "
            "sign-of-volume and near-tie quality branch, so a second "
            "implementation could silently diverge into a different "
            "connectivity) -- install a build with the compiled "
            "meshioplusplus._core extension"
        )

    res = _core.optimize_volume(
        mesh,
        int(max_iterations),
        bool(relocate),
        bool(flip),
        bool(preserve_boundary),
        float(min_improvement),
        None if frozen is None else [int(v) for v in frozen],
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {
        "num_flips": res["num_flips"],
        "num_23_flips": res["num_23_flips"],
        "num_32_flips": res["num_32_flips"],
        "num_vertices_moved": res["num_vertices_moved"],
        "num_tets": res["num_tets"],
        "min_quality_before": res["min_quality_before"],
        "min_quality_after": res["min_quality_after"],
    }
    return out, report
