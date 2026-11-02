"""A metric toolbox that needs no MMG (roadmap 7.4.3): intersect the
anisotropic metrics of several fields, normalise to a target vertex count,
clamp and grade, ready for ``adapt``. See ``doc/metric.md``.

**C++-core only, with no numpy fallback** (the ``metric_from_hessian``
precedent): per-point symmetric eigen-decompositions and a swept gradation
whose order and stopping test a second implementation would have to match
bit for bit.
"""

from __future__ import annotations

__all__ = ["prepare_metric"]


def prepare_metric(
    mesh,
    arrays,
    *,
    output: str = "metric",
    target_vertices: float = 0.0,
    norm: float = 2.0,
    hmin: float = 0.0,
    hmax: float = 0.0,
    gradation: float = 0.0,
    max_sweeps: int = 100,
    return_report: bool = False,
):
    """Intersect, normalise, clamp and grade anisotropic metric fields.

    Reads one or more ``point_data`` metrics (six components ``m11 m21 m22
    m31 m32 m33`` in 3-D, three ``m11 m21 m22`` in 2-D, as
    :func:`metric_from_hessian` writes) and, per point: intersects them
    (simultaneous reduction: the smallest metric asking for at least what
    each does in every direction), normalises to ``target_vertices`` (the
    continuous-mesh L^p normalisation, ``norm`` = p, 0 for infinity), clamps
    every size to ``[hmin, hmax]`` and grades the field so sizes grow by at
    most ``gradation`` per unit length (Alauzet's metric-space gradation).

    :param arrays: a metric field's name, or several.
    :param output: the ``point_data`` name of the result.
    :param target_vertices: the vertex count to normalise to (0: skip).
    :param norm: the L^p norm of the normalisation (at least 1; 0: infinity).
    :param hmin, hmax: the size range (0: that side open).
    :param gradation: the largest size growth per unit length, > 1 (0: skip).
    :param max_sweeps: gradation sweeps at most.
    :param return_report: also return ``{"complexity_before",
        "complexity_after", "num_graded", "sweeps"}`` (complexity: the
        integral of ``sqrt(det M)``, about the vertex count a mesher makes).
    :returns: the mesh with ``point_data[output]``, or ``(mesh, report)``.
    :raises ValueError: on a missing or misshapen field, or bad options.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    if isinstance(arrays, str):
        arrays = [arrays]
    res = _core_op("prepare_metric")(
        mesh,
        [str(a) for a in arrays],
        str(output),
        float(target_vertices),
        float(norm),
        float(hmin),
        float(hmax),
        float(gradation),
        int(max_sweeps),
    )
    out = res.pop("mesh")
    return (out, res) if return_report else out
