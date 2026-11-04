"""Field-aligned quad remeshing of a triangle surface (roadmap 7.7.3).

:func:`remesh_quads` hands a triangle surface to a quad remesher and returns
its quads: QuadriFlow (Huang et al. 2018, MIT) or Instant Meshes (Jakob et
al. 2015, BSD-3), installed programs the C++ core runs as separate processes
(found on ``PATH``, or through ``MESHIOPLUSPLUS_QUADRIFLOW`` and
``MESHIOPLUSPLUS_INSTANT_MESHES``; never in WebAssembly). They align the quads
with a smooth cross field, so a sphere comes out all quads without the poor
quads at its poles that recombining a mesher's own triangles leaves. The
result is a new surface near the input; the report says how near. See
``doc/remesh_quads.md``.

The C++ core does the work; there is no numpy twin.
"""

from __future__ import annotations

import warnings

from ._mesh_cad import _core_or_raise

__all__ = ["remesh_quads", "remesh_quads_backends"]

_PREFIX = "meshio++: remesh_quads: "

_REPORT_INTS = ("num_quads", "num_triangles", "num_points")
_REPORT_FLOATS = (
    "min_scaled_jacobian",
    "mean_scaled_jacobian",
    "hausdorff_to_surface",
)


def remesh_quads_backends() -> list:
    """The ``backend`` names :func:`remesh_quads` takes besides ``"auto"``,
    sorted: every registered quad remesher, the installed ``"quadriflow"`` and
    ``"instant_meshes"`` programs among them."""
    try:
        from . import _core
    except ImportError:
        return []
    from ._plugins import load_plugins

    load_plugins()
    names = getattr(_core, "remesh_quads_backends", None)
    return sorted(set(names())) if names is not None else []


def _int_option(name, value, lo, hi):
    if isinstance(value, bool) or int(value) != value:
        raise ValueError(f"{_PREFIX}{name} must be an integer")
    if not lo <= int(value) <= hi:
        raise ValueError(f"{_PREFIX}{name} must be between {lo} and {hi}")
    return int(value)


def remesh_quads(
    mesh,
    *,
    target_faces: int = 0,
    backend: str = "auto",
    sharp: bool = False,
    boundary: bool = True,
    adaptive: bool = False,
    dominant: bool = False,
    seed: int = 0,
    max_cells: int = 50_000_000,
    return_report: bool = False,
):
    """Remesh the triangle surface ``mesh`` into field-aligned quads.

    :param mesh: the surface: ``triangle`` blocks (lines and vertices are
        skipped); quads, polygons, higher-order and volume cells are refused.
    :param target_faces: the number of faces asked for; 0 (the default) takes
        half the input's triangles. The programs treat it as a target.
    :param backend: ``"auto"`` (QuadriFlow, else Instant Meshes) or a name of
        :func:`remesh_quads_backends`.
    :param sharp: align with and keep the sharp edges.
    :param boundary: keep the boundary of an open surface.
    :param adaptive: let the face size follow the curvature (QuadriFlow).
    :param dominant: allow triangles where quads fit badly (Instant Meshes'
        quad-dominant mode; QuadriFlow is always all quads).
    :param seed: the program's random seed, so a run repeats.
    :param max_cells: refuse a target or result of more faces.
    :param return_report: also return the report.
    :returns: the mesh (a ``quad`` block, and a ``triangle`` block when the
        remesher left any), or ``(mesh, report)``: ``num_quads``,
        ``num_triangles``, ``num_points``, ``min_scaled_jacobian``,
        ``mean_scaled_jacobian``, ``hausdorff_to_surface`` (the two-sided
        distance to the input), ``backend`` and ``warnings`` (each also
        emitted as a :class:`RuntimeWarning`).
    :raises ValueError: for an option out of range, or a surface that is not
        linear triangles.
    :raises NotImplementedError: when no quad remesher is installed, or none of
        that name (the message lists the available ones).
    :raises RuntimeError: when the remesher fails or writes a bad mesh.
    """
    if not isinstance(backend, str) or not backend.strip():
        raise ValueError(
            f"{_PREFIX}backend must be 'auto' or a quad remesher's name "
            f"(remesh_quads_backends()), not {backend!r}"
        )
    target_faces = _int_option("target_faces", target_faces, 0, 2**63 - 1)
    max_cells = _int_option("max_cells", max_cells, 1, 2**63 - 1)
    seed = _int_option("seed", seed, -(2**31), 2**31 - 1)
    _core = _core_or_raise("remesh_quads")
    res = _core.remesh_quads(
        mesh,
        target_faces=target_faces,
        backend=backend,
        sharp=bool(sharp),
        boundary=bool(boundary),
        adaptive=bool(adaptive),
        dominant=bool(dominant),
        seed=seed,
        max_cells=max_cells,
    )
    out = res["mesh"]
    found = [str(w) for w in res["warnings"]]
    for w in found:
        warnings.warn(f"{_PREFIX}{w}", RuntimeWarning, stacklevel=2)
    if not return_report:
        return out
    report = {k: int(res[k]) for k in _REPORT_INTS}
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report["backend"] = str(res["backend"])
    report["warnings"] = found
    return out, report
