"""Cut a mesh by level sets, and quadrature on the cut cells (doc/cut.md).

:func:`cut_level_set` returns the part of a mesh that one or several level
sets select: the pieces of the cells they cross, on one side of them or on
their zero sets, and the cells wholly inside. :func:`cut_quadrature` returns
quadrature rules on that part, the cut-cell quadrature of unfitted finite
elements (CutFEM). Both run CutCells (S. Claus, MIT), an optional dependency
of the C++ core (``MESHIOPLUSPLUS_WITH_CUTCELLS``; :func:`cut_has_cutcells`).

The C++ core does the work; there is no numpy twin.
"""

from __future__ import annotations

import math
import warnings

from ._mesh_cad import _core_or_raise

__all__ = ["cut_has_cutcells", "cut_level_set", "cut_quadrature"]

_PREFIX = "meshio++: cut: "
_BACKENDS = ("auto", "lut", "quadrays")


def cut_has_cutcells() -> bool:
    """Whether the C++ core carries CutCells, which both cut operations need."""
    try:
        from . import _core
    except ImportError:
        return False
    return bool(getattr(_core, "__has_cutcells__", False))


def _names(level_sets):
    if level_sets is None:
        return []
    if isinstance(level_sets, str):
        return [level_sets]
    out = []
    for name in level_sets:
        if not isinstance(name, str) or not name:
            raise ValueError(
                f"{_PREFIX}level_sets are point_data array names, not {name!r}"
            )
        out.append(name)
    return out


def _planes(planes):
    if planes is None:
        return []
    out = []
    for i, plane in enumerate(planes):
        try:
            origin, normal = plane
            origin = [float(v) for v in origin]
            normal = [float(v) for v in normal]
        except (TypeError, ValueError):
            raise ValueError(
                f"{_PREFIX}plane {i} must be (origin, normal), each three numbers"
            ) from None
        while len(origin) < 3:
            origin.append(0.0)
        while len(normal) < 3:
            normal.append(0.0)
        if len(origin) != 3 or len(normal) != 3:
            raise ValueError(
                f"{_PREFIX}plane {i} must be (origin, normal), each three numbers"
            )
        if not all(math.isfinite(v) for v in origin + normal) or not any(normal):
            raise ValueError(
                f"{_PREFIX}plane {i} needs a finite origin and a nonzero normal"
            )
        out.append((tuple(origin), tuple(normal)))
    return out


def _int(name, value, lo, hi):
    if isinstance(value, bool) or int(value) != value:
        raise ValueError(f"{_PREFIX}{name} must be an integer")
    if not lo <= int(value) <= hi:
        raise ValueError(f"{_PREFIX}{name} must be between {lo} and {hi}")
    return int(value)


def _common(level_sets, planes, selection, backend):
    names = _names(level_sets)
    planes = _planes(planes)
    if not names and not planes:
        raise ValueError(
            f"{_PREFIX}give a level set: a point_data array (level_sets) or a plane (planes)"
        )
    if selection is None:
        selection = ""
    if not isinstance(selection, str):
        raise ValueError(f"{_PREFIX}selection must be a string such as 'phi < 0'")
    if backend not in _BACKENDS:
        raise ValueError(
            f"{_PREFIX}backend must be one of {', '.join(_BACKENDS)}, not {backend!r}"
        )
    return names, planes, selection


def _warn(found):
    for w in found:
        warnings.warn(
            w if w.startswith("meshio++") else f"{_PREFIX}{w}",
            RuntimeWarning,
            stacklevel=3,
        )


def cut_level_set(
    mesh,
    level_sets=None,
    *,
    planes=None,
    selection: str = "",
    include_uncut: bool = True,
    backend: str = "auto",
    degree: int = 2,
    triangulate: bool = False,
    record_parent_ids: bool = True,
    return_report: bool = False,
):
    """The part of ``mesh`` that the level sets select.

    :param mesh: the mesh. Its cells of the highest dimension are cut
        (``line``, ``triangle``, ``quad``, ``tetra``, ``hexahedron``,
        ``wedge``, ``pyramid`` and their higher-order forms, through their
        corners); lower-dimensional blocks are left out with a warning.
        Polygons and polyhedra are refused.
    :param level_sets: a scalar ``point_data`` array name, or several. On
        ``triangle6``/``tetra10`` cells a level set is quadratic.
    :param planes: plane level sets ``(origin, normal)``, ``(x - origin) .
        normal``, named ``plane0``, ``plane1``, ... in ``selection``.
    :param selection: CutCells' expression over the level sets' names, such as
        ``"phi < 0"``, ``"phi = 0"`` (the interface) or ``"wall < 0 and hole >
        0"``. Empty takes ``"<first> < 0"``.
    :param include_uncut: keep the cells wholly in the part, not only the
        pieces of the cut ones.
    :param backend: ``"lut"`` (straight pieces by lookup tables), ``"quadrays"``
        (curved Lagrange pieces of ``degree``) or ``"auto"`` (quadrays for a
        quadratic level set, the lookup tables otherwise).
    :param degree: quadrays only: the pieces' Lagrange degree, 1 to 8.
    :param triangulate: lut only: split the pieces into simplices.
    :param record_parent_ids: add ``cut:parent_cell``, the input cell each
        piece came from (global, block-major).
    :param return_report: also return the report.
    :returns: the mesh, with point data interpolated in each piece's parent,
        cell data and cell regions carried from it and ``cut:cut`` (1 on a
        piece of a cut cell), or ``(mesh, report)``: ``num_cells``,
        ``num_cut_cells``, ``num_uncut_cells``, ``num_pieces``, ``dimension``,
        ``measure`` (the part's volume, area or length), ``selection``,
        ``backend`` and ``warnings`` (each also a :class:`RuntimeWarning`).
    :raises ValueError: for no level set, a missing or non-scalar array, a bad
        selection or option.
    :raises NotImplementedError: without CutCells, or for polygons and
        polyhedra.
    """
    names, planes, selection = _common(level_sets, planes, selection, backend)
    degree = _int("degree", degree, 1, 8)
    _core = _core_or_raise("cut_level_set")
    res = _core.cut_level_set(
        mesh,
        level_sets=names,
        planes=planes,
        selection=selection,
        include_uncut=bool(include_uncut),
        backend=backend,
        degree=degree,
        triangulate=bool(triangulate),
        record_parent_ids=bool(record_parent_ids),
    )
    found = [str(w) for w in res["warnings"]]
    _warn(found)
    if not return_report:
        return res["mesh"]
    report = {
        k: int(res[k])
        for k in (
            "num_cells",
            "num_cut_cells",
            "num_uncut_cells",
            "num_pieces",
            "dimension",
        )
    }
    report["measure"] = float(res["measure"])
    report["selection"] = str(res["selection"])
    report["backend"] = str(res["backend"])
    report["warnings"] = found
    return res["mesh"], report


def cut_quadrature(
    mesh,
    level_sets=None,
    *,
    planes=None,
    selection: str = "",
    order: int = 2,
    include_uncut: bool = True,
    backend: str = "auto",
    triangulate: bool = False,
    record_parent_ids: bool = True,
    return_report: bool = False,
):
    """Quadrature rules of degree ``order`` on the part the level sets select.

    The options are :func:`cut_level_set`'s. ``order`` is the polynomial
    degree the rules integrate exactly (on straight pieces), 0 to 20.

    :returns: a mesh of one ``vertex`` per quadrature point at its physical
        position, with point data ``quadrature:weight`` (physical: summed, the
        part's measure), ``quadrature:reference`` (in the parent cell's
        reference element), ``quadrature:rule`` and ``quadrature:parent_cell``;
        or ``(mesh, report)``: ``num_rules``, ``num_points``, ``sum_weights``,
        ``dimension``, ``selection``, ``backend`` and ``warnings``.
    :raises: as :func:`cut_level_set`.
    """
    names, planes, selection = _common(level_sets, planes, selection, backend)
    order = _int("order", order, 0, 20)
    _core = _core_or_raise("cut_quadrature")
    res = _core.cut_quadrature(
        mesh,
        level_sets=names,
        planes=planes,
        selection=selection,
        include_uncut=bool(include_uncut),
        backend=backend,
        order=order,
        triangulate=bool(triangulate),
        record_parent_ids=bool(record_parent_ids),
    )
    found = [str(w) for w in res["warnings"]]
    _warn(found)
    if not return_report:
        return res["mesh"]
    report = {k: int(res[k]) for k in ("num_rules", "num_points", "dimension")}
    report["sum_weights"] = float(res["sum_weights"])
    report["selection"] = str(res["selection"])
    report["backend"] = str(res["backend"])
    report["warnings"] = found
    return res["mesh"], report
