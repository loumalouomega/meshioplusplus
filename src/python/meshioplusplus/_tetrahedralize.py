"""Tetrahedra filling a closed triangle surface (roadmap 6.7).

The input is a closed, manifold, consistently oriented triangle surface, such
as a :func:`meshioplusplus.mesh_cad_surface` result or a boundary layer's cap:
an open, non-manifold, misoriented or degenerate one is refused by name. Its
connected shells are oriented outward, and each shell bounds a region (the
inside of the shell minus the shells directly inside it): a box inside a box
gives two regions.

Two backends. Netgen's nglib (J. Schöberl, *NETGEN: an advancing front
2D/3D-mesh generator based on abstract rules*, 1997), linked when the core is
built with ``MESHIOPLUSPLUS_WITH_NETGEN``, keeps the surface's triangles and
points exactly, which is checked, so neighbouring regions conform. The native
octree backend, always built, is :func:`meshioplusplus.octree_mesh`'s
isosurface stuffing at the surface's own size: its boundary is fitted to the
surface within the reported ``hausdorff_to_surface``. ``backend="auto"`` takes
Netgen when ``meshioplusplus._core.__has_netgen__`` and the octree backend
otherwise, with a warning. See ``doc/tetrahedralize.md``.

The C++ core does the work; there is no numpy twin.
"""

from __future__ import annotations

import math
import warnings

from ._mesh_cad import _core_or_raise

__all__ = ["tetrahedralize"]

_PREFIX = "meshio++: tetrahedralize: "

_BACKENDS = ("auto", "netgen", "octree")

_REPORT_INTS = (
    "num_tetrahedra",
    "num_points",
    "num_new_points",
    "num_regions",
    "num_boundary_triangles",
)
_REPORT_FLOATS = (
    "hausdorff_to_surface",
    "min_dihedral_deg",
    "max_dihedral_deg",
    "volume",
)


def _int_option(name, value, lo, hi):
    if isinstance(value, bool) or int(value) != value:
        raise ValueError(f"{_PREFIX}{name} must be an integer")
    if not lo <= int(value) <= hi:
        raise ValueError(f"{_PREFIX}{name} must be between {lo} and {hi}")
    return int(value)


def tetrahedralize(
    mesh,
    *,
    backend: str = "auto",
    max_size: float = 0.0,
    max_volume: float = 0.0,
    grading: float = 0.3,
    optimize_steps: int = 3,
    max_cells: int = 50_000_000,
    face_array: str = "cad:face",
    return_report: bool = False,
):
    """Fill the closed triangle surface ``mesh`` with tetrahedra.

    :param mesh: the surface: ``triangle`` blocks (lines and vertices are
        skipped), closed, manifold, consistently oriented and without
        degenerate triangles; quads, polygons, higher-order and volume cells
        are refused.
    :param backend: ``"auto"`` (Netgen when built in, else ``"octree"``),
        ``"netgen"`` (the surface's triangles kept) or ``"octree"`` (the
        native backend, the boundary fitted).
    :param max_size: the largest edge asked of the tetrahedra (0, the
        default, is none: the surface's triangles size the inside).
    :param max_volume: the largest volume asked of a tetrahedron (0 is none);
        a regular tetrahedron's edge for it joins ``max_size``.
    :param grading: how fast the size grows away from the surface, in
        ``[0, 1]`` (Netgen).
    :param optimize_steps: passes of the backend's volume optimisation (0 to
        100).
    :param max_cells: refuse a result of more tetrahedra.
    :param face_array: the cell data naming each triangle's CAD face, carried
        to the triangles and made ``face:<v>`` Side regions when present.
    :param return_report: also return the report.
    :returns: the mesh (a ``tetra`` block and a ``triangle`` block: with
        Netgen the input's points first and unchanged and the surface's
        triangles, with the octree backend new points and the tetrahedra's
        boundary facets; the cell data ``tetrahedralize:region``;
        ``region:<k>`` Cell regions), or ``(mesh, report)``: ``num_tetrahedra``, ``num_points``,
        ``num_new_points``, ``num_regions``, ``num_boundary_triangles``,
        ``boundary_preserved``, ``hausdorff_to_surface``,
        ``min_dihedral_deg``, ``max_dihedral_deg``, ``volume``, ``backend``
        and ``warnings`` (each also emitted as a :class:`RuntimeWarning`).
    :raises ValueError: for an option out of range, a surface that has no
        triangle, another cell type, or is open, non-manifold, misoriented or
        degenerate, or a result past ``max_cells``.
    :raises NotImplementedError: for ``backend="netgen"`` in a core built
        without ``MESHIOPLUSPLUS_WITH_NETGEN``.
    :raises RuntimeError: when the backend fails (a self-intersecting
        surface, for one).
    """
    if not isinstance(backend, str) or backend not in _BACKENDS:
        raise ValueError(
            f"{_PREFIX}backend must be 'auto', 'netgen' or 'octree', not {backend!r}"
        )
    max_size = float(max_size)
    if not (math.isfinite(max_size) and max_size >= 0.0):
        raise ValueError(f"{_PREFIX}max_size must be finite and not negative")
    max_volume = float(max_volume)
    if not (math.isfinite(max_volume) and max_volume >= 0.0):
        raise ValueError(f"{_PREFIX}max_volume must be finite and not negative")
    grading = float(grading)
    if not (math.isfinite(grading) and 0.0 <= grading <= 1.0):
        raise ValueError(f"{_PREFIX}grading must lie in [0, 1]")
    optimize_steps = _int_option("optimize_steps", optimize_steps, 0, 100)
    max_cells = _int_option("max_cells", max_cells, 1, 2**63 - 1)
    _core = _core_or_raise("tetrahedralize")
    res = _core.tetrahedralize(
        mesh,
        backend=backend,
        max_size=max_size,
        max_volume=max_volume,
        grading=grading,
        optimize_steps=optimize_steps,
        max_cells=max_cells,
        face_array=str(face_array),
    )
    out = res["mesh"]
    found = [str(w) for w in res["warnings"]]
    for w in found:
        warnings.warn(f"{_PREFIX}{w}", RuntimeWarning, stacklevel=2)
    if not return_report:
        return out
    report = {k: int(res[k]) for k in _REPORT_INTS}
    report["boundary_preserved"] = bool(res["boundary_preserved"])
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report["backend"] = str(res["backend"])
    report["warnings"] = found
    return out, report
