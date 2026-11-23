"""Screened-Poisson smoothing of a point field on a triangle surface (roadmap 15.2.2).

``(M + l**2 K) g = M f`` over the lumped mass ``M`` and the cotangent stiffness
``K`` of :func:`~meshioplusplus.mesh_operator`: the field is smoothed on the
scale ``l`` and the points do not move. See ``doc/mesh_operator.md``.

**C++-core only, with no numpy fallback** (the ``geodesic_distance``
precedent): the conjugate-gradient solve stops on a rounded test.
"""

from __future__ import annotations

__all__ = ["smooth_field", "SUFFIX"]

#: The default output name is ``<array>`` + this. Twins ``kSmoothFieldSuffix``.
SUFFIX = ":smooth"


def smooth_field(
    mesh,
    array,
    length,
    *,
    output=None,
    tolerance: float = 1e-12,
    max_iterations: int = 0,
    return_report: bool = False,
):
    """Smooth a point array over a triangle surface on the scale ``length``.

    Returns the mesh with ``point_data["<array>:smooth"]`` (or ``output``): the
    solution of ``(M + length**2 K) g = M f``. Features shorter than ``length``
    are flattened and longer ones pass nearly unchanged; the ``M``-weighted
    mean is preserved and ``length = 0`` returns the field. Each component of a
    vector array is smoothed alone. A point in no triangle keeps its value.

    :param mesh: a surface (its ``triangle`` blocks), never modified.
    :param array: the name of a numeric ``point_data`` array.
    :param length: the smoothing length in the mesh's units, ``>= 0``.
    :param return_report: also return ``{"num_degenerate", "iterations",
        "converged"}``.
    :raises ValueError: on an unknown array, a negative length, or a mesh
        outside :func:`~meshioplusplus.mesh_operator`'s scope.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "smooth_field"):
        raise NotImplementedError(
            "meshio++: smooth_field: this operation is C++-core only and has no "
            "pure-numpy fallback -- install a build with the compiled "
            "meshioplusplus._core extension"
        )
    res = _core.smooth_field(
        mesh,
        str(array),
        float(length),
        output or "",
        float(tolerance),
        int(max_iterations),
    )
    out = res.pop("mesh")
    if not res["converged"]:
        from ._common import warn

        warn(
            f"smooth_field: the solve hit its iteration cap ({res['iterations']}); "
            "the field is approximate -- raise max_iterations"
        )
    return (out, res) if return_report else out
