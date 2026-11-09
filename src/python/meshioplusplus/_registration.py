"""Rigid registration by iterative closest points (roadmap 7.5.3).

``register_rigid(source, target)`` returns the 4x4 matrix that lays the
source onto the target -- apply it with ``transform(source, matrix=m)`` --
and a report. Point-to-point (Horn's closed form) or point-to-plane
(linearised, applied as an exact rotation); closest points on the target's
surface, or its nearest points for a cloud. See ``doc/registration.md``.

**C++-core only, with no numpy fallback** (the ``geodesic_distance``
precedent): the iteration stops on a rounded test, so a second
implementation could stop a step early or late and return a different
matrix in the last digits. This module raises when ``_core`` cannot be used.
"""

from __future__ import annotations

import numpy as np

__all__ = ["register_rigid", "register_landmarks", "find_symmetry_plane", "auto_align"]

_METHODS = ("point_to_point", "point_to_plane")


def register_rigid(
    source,
    target,
    *,
    method: str = "point_to_point",
    max_iterations: int = 50,
    tolerance: float = 1e-10,
    max_distance=0.0,
    initial=None,
    symmetric_check: bool = False,
    symmetric_tolerance: float = 0.05,
    source_region: str | None = None,
    target_region: str | None = None,
    adaptive_trim: float = 0.0,
):
    """The rigid motion laying ``source`` onto ``target``, by iterative
    closest points (Besl and McKay 1992).

    :param source: the mesh or ``(n, 3)`` points to move (only its points are
        read).
    :param target: the mesh to lay it onto: closest points on its surface when
        it has surface cells (a volume cell is refused naming
        ``extract_surface``), else its nearest points (a cloud).
    :param method: ``"point_to_point"`` or ``"point_to_plane"`` (the distance
        to the target's tangent planes; usually far fewer iterations; a cloud
        target needs ``normals`` point data, see ``estimate_normals``).
    :param max_iterations: the iteration cap.
    :param tolerance: stop when a step rotates by less than this (radians)
        and moves the source's centre by less than this times its
        bounding-box diagonal.
    :param max_distance: drop pairs farther apart than this; 0 keeps all.
        ``"auto"`` trims adaptively instead: ``adaptive_trim=3``.
    :param initial: a starting 4x4 rigid matrix (row-major), ``None`` (the
        identity), or ``"auto"``: the best of the 24 principal-frame candidates,
        each refined by a coarse point-to-point ICP on a subsample (no start
        needed; see doc/registration.md).
    :param symmetric_check: with ``initial="auto"``, warn when a candidate of
        a pose more than 5 degrees away fits within ``symmetric_tolerance``
        (relative trimmed RMS) of the chosen one: a symmetric part.
    :param source_region: pair only the points of this cell region's surface
        cells of ``source`` (the part, without its fixture).
    :param target_region: match only this cell region of ``target``.
    :param adaptive_trim: drop, at each iteration, the pairs farther apart than
        this multiple of the median pair distance; 0 keeps them.
    :returns: ``(matrix, report)``: the 4x4 ``numpy`` matrix mapping source
        coordinates onto the target's (``transform(source, matrix=matrix)``
        applies it), and ``{"initial_rms", "rms", "iterations", "num_pairs",
        "converged", "num_rejected", "candidate_rms", "chosen_candidate",
        "warnings"}``. Without ``initial="auto"`` ICP finds the nearest local
        optimum: start roughly aligned.
    :raises ValueError: on an unknown method, an empty or non-finite input,
        a target with volume cells, ``point_to_plane`` against a cloud with no
        ``normals``, or no pair within ``max_distance``.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _as_mesh, _core_op

    if method not in _METHODS:
        raise ValueError(
            f"meshio++: register_rigid: unknown method '{method}' (expected "
            "'point_to_point' or 'point_to_plane')"
        )
    initial_mode = "matrix"
    if isinstance(initial, str):
        if initial != "auto":
            raise ValueError(
                f"meshio++: register_rigid: unknown initial '{initial}' (expected a "
                "4x4 matrix, None or 'auto')"
            )
        initial_mode, start = "auto", np.eye(4)
    elif initial is None:
        start = np.eye(4)
    else:
        start = np.asarray(initial, dtype=np.float64).reshape(4, 4)
    if isinstance(max_distance, str):
        if max_distance != "auto":
            raise ValueError(
                f"meshio++: register_rigid: unknown max_distance '{max_distance}' "
                "(expected a number or 'auto')"
            )
        max_distance, adaptive_trim = 0.0, adaptive_trim or 3.0
    res = _core_op("register_rigid")(
        _as_mesh(source),
        _as_mesh(target),
        method,
        int(max_iterations),
        float(tolerance),
        float(max_distance),
        [float(v) for v in start.ravel()],
        initial_mode,
        bool(symmetric_check),
        float(symmetric_tolerance),
        source_region or "",
        target_region or "",
        float(adaptive_trim),
    )
    matrix = np.asarray(res.pop("matrix"), dtype=np.float64).reshape(4, 4)
    res["candidate_rms"] = [float(v) for v in res["candidate_rms"]]
    res["warnings"] = [str(w) for w in res["warnings"]]
    if res["warnings"]:
        from ._common import warn

        for w in res["warnings"]:
            warn(f"register_rigid: {w}")
    if not res["converged"]:
        from ._common import warn

        warn(
            f"register_rigid: no convergence in {res['iterations']} iterations "
            f"(rms {res['rms']:.3g}); raise max_iterations or start closer"
        )
    return matrix, res


def _points3(values, name):
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"meshio++: register_landmarks: {name} must be (n, 3) points")
    return [tuple(float(v) for v in row) for row in arr]


def register_landmarks(source_points, target_points, weights=None):
    """The rigid motion laying picked ``source_points`` onto ``target_points``:
    the weighted Horn solve in closed form (roadmap 7.3.2), the start a viewer
    or a script gives ``register_rigid``.

    :param source_points: ``(n, 3)`` picked points, ``n >= 3``, not collinear.
    :param target_points: their ``(n, 3)`` targets.
    :param weights: ``n`` non-negative weights, or ``None`` for all ones.
    :returns: ``(matrix, report)``: the 4x4 matrix and ``{"residuals", "rms"}``.
    :raises ValueError: on fewer than three pairs, mismatched sizes, a bad
        weight or collinear source points.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    src = _points3(source_points, "source_points")
    dst = _points3(target_points, "target_points")
    w = [] if weights is None else [float(v) for v in np.asarray(weights).ravel()]
    res = _core_op("register_landmarks")(src, dst, w)
    matrix = np.asarray(res["matrix"], dtype=np.float64).reshape(4, 4)
    return matrix, {
        "residuals": [float(v) for v in res["residuals"]],
        "rms": float(res["rms"]),
    }


def find_symmetry_plane(
    mesh,
    *,
    seed=None,
    tolerance: float = 0.0,
    min_share: float = 0.9,
    max_samples: int = 4000,
    max_iterations: int = 50,
):
    """The plane ``mesh`` is most nearly symmetric about (roadmap 7.3.6).

    Each principal plane through the centroid (and ``seed``) is tried: the
    mesh is mirrored across it, the image registered to the original by ICP,
    and the plane read off the product of the reflection and the motion.

    :param mesh: a surface or a point cloud.
    :param seed: an extra candidate ``(origin, normal)``.
    :param tolerance: a sample whose mirror image lies within this of the
        surface counts as symmetric; 0 uses 0.005 of the diagonal.
    :param min_share: refuse a plane with a smaller share of symmetric samples.
    :param max_samples: the samples of the fine fit (a stride over the points).
    :param max_iterations: the iteration cap of each ICP run.
    :returns: ``{"origin", "normal", "rms", "share", "tolerance",
        "candidate_rms", "chosen_candidate", "warnings"}``; ``origin`` is the
        centroid projected onto the plane.
    :raises ValueError: on a bad option, or when no plane reaches ``min_share``.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _as_mesh, _core_op

    so, sn = (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)
    if seed is not None:
        so = tuple(float(v) for v in np.asarray(seed[0], dtype=np.float64).ravel())
        sn = tuple(float(v) for v in np.asarray(seed[1], dtype=np.float64).ravel())
        if len(so) != 3 or len(sn) != 3 or not any(sn):
            raise ValueError(
                "meshio++: find_symmetry_plane: seed must be (origin, normal), "
                "a non-zero normal"
            )
    res = _core_op("find_symmetry_plane")(
        _as_mesh(mesh),
        so,
        sn,
        float(tolerance),
        float(min_share),
        int(max_samples),
        int(max_iterations),
    )
    out = {k: float(res[k]) for k in ("rms", "share", "tolerance")}
    out["origin"] = [float(v) for v in res["origin"]]
    out["normal"] = [float(v) for v in res["normal"]]
    out["candidate_rms"] = [float(v) for v in res["candidate_rms"]]
    out["chosen_candidate"] = int(res["chosen_candidate"])
    out["warnings"] = [str(w) for w in res["warnings"]]
    return out


def auto_align(
    mesh,
    *,
    angle_tolerance: float = 10.0,
    subdivisions: int = 3,
    region: str | None = None,
    rotate_vector_data: bool = False,
    return_report: bool = False,
):
    """Pose ``mesh`` from the scan alone (roadmap 7.3.5): on its largest flat
    face, or a cylindrical part on the end of its axis, from an area-weighted
    histogram of its facet normals; an open scan's opening counts as a face.

    :param mesh: a surface (its triangles, quads and polygons vote).
    :param angle_tolerance: degrees within which a normal agrees with a
        direction, in (0, 45).
    :param subdivisions: the Gauss-sphere icosahedron's splits (642 bins at 3).
    :param region: only this cell region's surface cells vote.
    :param rotate_vector_data: rotate vector and tensor point data too.
    :param return_report: also return the report.
    :returns: the aligned mesh (its bounding box centred on the Z axis, lowest
        point at z = 0), or ``(mesh, report)`` with ``matrix`` (4x4), ``mode``
        (``"face"`` or ``"axis"``), ``primary`` and ``secondary`` (input
        directions sent to -Z or +Z and to +X), ``share``, ``face_share``,
        ``axis_share``, ``opening_area`` and ``warnings``.
    :raises ValueError: on a bad option or no surface cell of positive area.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    res = _core_op("auto_align")(
        mesh,
        float(angle_tolerance),
        int(subdivisions),
        region or "",
        bool(rotate_vector_data),
    )
    out = res["mesh"]
    if not return_report:
        return out
    report = {
        "matrix": np.asarray(res["matrix"], dtype=np.float64).reshape(4, 4),
        "mode": str(res["mode"]),
        "primary": [float(v) for v in res["primary"]],
        "secondary": [float(v) for v in res["secondary"]],
        "warnings": [str(w) for w in res["warnings"]],
    }
    for k in ("share", "face_share", "axis_share", "opening_area"):
        report[k] = float(res[k])
    return out, report
