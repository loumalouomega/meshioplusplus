"""Write a named region chosen by geometry, data, quality, component or growth.

The Python face of ``operations/select.hpp``. ``select`` adds one region of
points, cells or sides to a copy of the mesh and leaves everything else --
points, connectivity, data, the other regions -- as it was, so a
boundary-condition or contact patch becomes a line of a script rather than a
pick in a GUI.

The operation is C++-core only: the patch growth walks the shared crease test
and the quality criterion the core's metrics, and a numpy twin of either would
be a second implementation that has to stay byte-identical (the
``remesh_volume`` precedent). Without the compiled core it raises
``NotImplementedError``.
"""

from __future__ import annotations

import numpy as np

_PREFIX = "meshio++: select: "
_COMPARISONS = ("<", "<=", ">", ">=", "==", "!=")
_REPORT_KEYS = ("num_candidates", "num_selected", "num_grown", "num_components")


def _vec3(value, what):
    try:
        out = [float(v) for v in value]
    except TypeError:
        raise ValueError(f"{_PREFIX}{what} must be three numbers") from None
    if len(out) != 3:
        raise ValueError(f"{_PREFIX}{what} must be three numbers")
    return out


def _box(box):
    flat = [float(v) for v in np.asarray(box, dtype=float).ravel()]
    if len(flat) != 6:
        raise ValueError(
            f"{_PREFIX}box= must be (xmin, ymin, zmin, xmax, ymax, zmax) or (lo3, hi3)"
        )
    return (flat[:3], flat[3:])


def _comparison(value, what):
    try:
        name, op, rhs = value
    except (TypeError, ValueError):
        raise ValueError(
            f"{_PREFIX}{what}= must be (name, comparison, value)"
        ) from None
    op = "==" if op == "=" else str(op)
    if op not in _COMPARISONS:
        raise ValueError(
            f"{_PREFIX}unknown comparison '{op}' (expected one of {', '.join(_COMPARISONS)})"
        )
    return (str(name), op, float(rhs))


def select(
    mesh,
    name: str,
    kind: str = "cell",
    *,
    box=None,
    plane=None,
    sphere=None,
    where=None,
    quality=None,
    component=None,
    seed=None,
    feature_angle: float = 30.0,
    normal_tolerance: float = 0.0,
    mode: str = "all",
    boundary_only: bool = True,
    tag: int = -1,
    return_report: bool = False,
):
    """Select points, cells or sides and write them as one named region.

    A candidate is selected when it passes **every** criterion given; at least
    one is required. See ``doc/select.md``.

    :param mesh: the input mesh (unmodified).
    :param name: the region's name. A region of the same name and kind is
        replaced, whatever its dimension or tag.
    :param kind: ``"cell"`` (global, block-major cell indices), ``"point"``, or
        ``"side"`` (``(cell, local facet)`` pairs: faces of the 3-D cells, or
        edges of the 2-D cells when there are none).
    :param box: ``(xmin, ymin, zmin, xmax, ymax, zmax)`` or ``(lo3, hi3)``.
    :param plane: ``(point3, normal3)``; inside is ``(p - point) . normal >= 0``.
    :param sphere: ``(center3, radius)``, the closed ball.
    :param where: ``(array, comparison, value)`` on a scalar ``point_data``
        array (points) or ``cell_data`` array (cells; a side takes its owner
        cell's value). A NaN never matches.
    :param quality: ``(metric, comparison, value)`` on a
        :func:`compute_quality` metric (``"scaled_jacobian"``,
        ``"min_dihedral"``, ...); cells and sides only.
    :param component: the k-th node-connected component, numbered as
        :func:`split` with ``by="component"`` numbers its pieces.
    :param seed: ``(x, y, z)``: grow a patch from the facet nearest this point
        across the edges :func:`feature_edges` would not draw.
    :param feature_angle: the largest dihedral angle, in degrees, the growth
        crosses; in ``[0, 180]``.
    :param normal_tolerance: the largest angle, in degrees, between a grown
        facet's normal and the seed's; ``0`` disables the test.
    :param mode: how a cell or side reduces its nodes for ``box``/``plane``/
        ``sphere`` (which together form one region, their intersection):
        ``"all"`` nodes inside, ``"any"`` node inside, or its ``"centroid"``.
    :param boundary_only: side regions consider only facets owned by one cell.
    :param tag: the region's integer tag, ``-1`` for none.
    :param return_report: also return ``num_candidates``, ``num_selected``,
        ``num_grown`` (``-1`` without a seed) and ``num_components`` (``-1``
        without a component criterion).
    :returns: the mesh with the region, or ``(mesh, report)``.
    :raises ValueError: on an empty name, no criterion, an unknown kind, mode,
        array, metric or comparison, a quality criterion on points, or a seed
        with nothing to grow over.
    :raises NotImplementedError: without the compiled C++ core.
    """
    name = str(name)
    if not name:
        raise ValueError(f"{_PREFIX}the region name must not be empty")
    if kind not in ("point", "cell", "side"):
        raise ValueError(
            f"{_PREFIX}unknown kind '{kind}' (expected 'point', 'cell' or 'side')"
        )
    if mode not in ("all", "any", "centroid"):
        raise ValueError(
            f"{_PREFIX}unknown mode '{mode}' (expected 'all', 'any' or 'centroid')"
        )
    args = {
        "box": None if box is None else _box(box),
        "plane": None,
        "sphere": None,
        "where": None if where is None else _comparison(where, "where"),
        "quality": None if quality is None else _comparison(quality, "quality"),
        "component": -1 if component is None else int(component),
        "seed": None if seed is None else _vec3(seed, "seed="),
    }
    if plane is not None:
        try:
            point, normal = plane
        except (TypeError, ValueError):
            raise ValueError(f"{_PREFIX}plane= must be (point3, normal3)") from None
        args["plane"] = (
            _vec3(point, "the plane's point"),
            _vec3(normal, "the plane's normal"),
        )
    if sphere is not None:
        try:
            center, radius = sphere
        except (TypeError, ValueError):
            raise ValueError(f"{_PREFIX}sphere= must be (center3, radius)") from None
        args["sphere"] = (_vec3(center, "the sphere's centre"), float(radius))
    if component is not None and args["component"] < 0:
        raise ValueError(f"{_PREFIX}component= must be >= 0")
    if all(
        args[k] is None for k in ("box", "plane", "sphere", "where", "quality", "seed")
    ) and (args["component"] < 0):
        raise ValueError(
            f"{_PREFIX}no criterion given (box, plane, sphere, where, quality, component or seed)"
        )

    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None:
        raise NotImplementedError(
            f"{_PREFIX}this operation is C++-core only and has no pure-numpy "
            "fallback (the patch growth and the quality metrics are the core's) "
            "-- install a build with the compiled meshioplusplus._core extension"
        )
    try:
        res = _core.select_region(
            mesh,
            name,
            kind,
            int(tag),
            mode,
            bool(boundary_only),
            args["box"],
            args["plane"],
            args["sphere"],
            args["where"],
            args["quality"],
            args["component"],
            args["seed"],
            float(feature_angle),
            float(normal_tolerance),
        )
    except TypeError as exc:
        # The core's std::invalid_argument arrives as ValueError already; a
        # value a cast refused is the caller's argument too.
        raise ValueError(f"{_PREFIX}{exc}") from exc
    out = res["mesh"]
    if not return_report:
        return out
    return out, {k: int(res[k]) for k in _REPORT_KEYS}


__all__ = ["select"]
