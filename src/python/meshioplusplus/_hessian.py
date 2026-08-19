"""The Hessian (second derivative) of a scalar ``point_data`` field, for
curvature-based adaptive refinement -- :func:`meshioplusplus.gradient`'s
companion one order further: ``gradient`` differentiates a field once, this
differentiates it twice.

**Composition, not a new numerical kernel.** ``hessian`` is built entirely out
of two calls to ``gradient``, exactly the precedent
:func:`meshioplusplus.estimate_error` already sets: (1) ``gradient(array,
location="point")`` gives the field's gradient as a genuine ``point_data``
array (``(n, 3)``); (2) ``gradient`` again on THAT array, with the default
``"gradient"`` operator, differentiates a 3-component field and so produces
``(n, 9)`` -- the flattened row-major 3x3 Hessian, ``H[i][j] =
d^2f/dxi dxj`` at index ``i*3+j``. The first call must use ``"point"``
location: the second call's ``point_data``-only validation would otherwise
reject a ``cell_data`` intermediate. ``method`` simply forwards to both
internal ``gradient`` calls.

**Exactness, stated honestly.** A field that is at most LINEAR has an
exactly-zero Hessian everywhere (its gradient is a constant, and Green-Gauss
of a constant field is trivially exact) -- the one mesh-shape-independent
guarantee. For a genuinely QUADRATIC field the composition is exact on a
structured/symmetric mesh away from its boundary (measured, not assumed --
see ``tests/cpp/test_hessian.cpp``) and a good, standard, but genuinely
approximate curvature estimate on an irregular mesh, because the mandatory
intermediate ``cell_data_to_point_data(weight="uniform")`` averaging step is
only exact for a field that is constant over the averaged neighbourhood.

**Scope: scalar fields only.** The input must have exactly one component. A
vector field's Hessian is a separate quantity per component -- call
``hessian`` once per component.

**Curvature-driven refinement needs no new marking step.**
:func:`meshioplusplus.data_calc`'s ``norm(...)`` is a plain
sum-of-squares-then-sqrt over however many components its argument has, so
``norm(<array>:hessian)`` on the 9-component output is exactly its Frobenius
norm -- a scalar curvature indicator with zero new code, ready for
``refine(mesh, where=...)``. See ``doc/hessian.md``'s worked composition.

The C++ core (``_core.hessian``) does the work; this module is the thin shim
(try C++, fall back to a pure-Python composition calling the already-composed
public :func:`meshioplusplus.gradient` twice -- which itself already has a
numpy fallback, so no separate numpy Hessian kernel is needed here) plus that
fallback itself.

Public API:

* :func:`hessian` -- the Hessian composition.
"""

from __future__ import annotations

from ._gradient import gradient

_PREFIX = "meshio++: hessian: "

#: Private working name for the intermediate gradient -- never reaches the
#: returned mesh (dropped/renamed before returning), so a clash with a real
#: user array only affects an internal, discarded intermediate copy.
_RAW_GRAD_NAME = "__hessian_raw_gradient__"

_METHOD_ALIASES = {
    "": "green-gauss",
    "green-gauss": "green-gauss",
    "green_gauss": "green-gauss",
    "gg": "green-gauss",
    "least-squares": "least-squares",
    "least_squares": "least-squares",
    "lsq": "least-squares",
}
_LOCATION_ALIASES = {"": "cell", "cell": "cell", "point": "point"}


def _num_components(mesh, array):
    arr = mesh.point_data[array]
    if arr.ndim < 2:
        return 1
    n = 1
    for d in arr.shape[1:]:
        n *= d
    return n


def _validate(mesh, array, method, location, output, overwrite):
    """Raise on every user error, before either path runs.

    Shared deliberately: the C++ core validates the same things, but the
    numpy-composition fallback must reject the same inputs with the same
    message even when ``_core`` is unavailable.
    """
    if method not in _METHOD_ALIASES:
        raise ValueError(
            f"{_PREFIX}unknown method '{method}' "
            "(expected 'green-gauss' or 'least-squares')"
        )
    if location not in _LOCATION_ALIASES:
        raise ValueError(f"{_PREFIX}unknown location '{location}' (expected 'cell' or 'point')")

    if not array:
        raise ValueError(f"{_PREFIX}an array name is required")
    if array not in mesh.point_data:
        if array in mesh.cell_data:
            raise ValueError(
                f"{_PREFIX}'{array}' is cell_data, which is piecewise constant "
                "and has no derivative; convert it first with "
                "cell_data_to_point_data (CLI: meshioplusplus data to-point)"
            )
        have = ", ".join(sorted(mesh.point_data))
        raise ValueError(
            f"{_PREFIX}no point_data array named '{array}'"
            + (f" (available: {have})" if have else " (the mesh has no point_data)")
        )
    comp = _num_components(mesh, array)
    if comp != 1:
        raise ValueError(
            f"{_PREFIX}'{array}' has {comp} components; hessian currently "
            "supports scalar fields only -- call it once per component of a "
            "vector field"
        )

    loc = _LOCATION_ALIASES[location]
    out_name = output or f"{array}:hessian"
    if not overwrite:
        taken = out_name in (mesh.point_data if loc == "point" else mesh.cell_data)
        if taken:
            raise ValueError(
                f"{_PREFIX}'{out_name}' already exists in "
                f"{'point_data' if loc == 'point' else 'cell_data'} "
                "(pass overwrite=True to replace it)"
            )


def _hessian_py(mesh, array, method, location, output, overwrite):
    """Pure-Python composition: two calls to the public ``gradient``.

    No separate numpy kernel -- ``gradient`` already has its own numpy
    fallback, so this reuses it exactly as ``estimate_error``'s own Python
    composition reuses the public averaging functions.
    """
    out_name = output or f"{array}:hessian"

    g1 = gradient(
        mesh, array, operator="gradient", method=method, location="point",
        output=_RAW_GRAD_NAME, overwrite=True,
    )
    g2, report = gradient(
        g1, _RAW_GRAD_NAME, operator="gradient", method=method, location="cell",
        output=out_name, overwrite=True, return_report=True,
    )

    out = g2
    out.cell_data.pop(_RAW_GRAD_NAME, None)

    if location == "point":
        from ._data_average import cell_data_to_point_data
        from ._data_manage import data_drop

        with_points = cell_data_to_point_data(
            out, keys=[out_name], weighted=False, overwrite=True
        )
        out = data_drop(with_points, "cell", [out_name])

    return out, report["num_skipped"], report["num_fallback"]


def hessian(
    mesh,
    array,
    method="green-gauss",
    location="cell",
    output=None,
    overwrite=False,
    return_report=False,
):
    """The Hessian (second derivative) of a scalar ``point_data`` field.

    :param mesh: the source mesh (never modified).
    :param array: name of the scalar ``point_data`` array to differentiate
        twice.
    :param method: ``"green-gauss"`` (default) or ``"least-squares"``,
        forwarded to both internal ``gradient`` passes.
    :param location: ``"cell"`` (default) or ``"point"`` for the result.
    :param output: output array name; ``None`` selects ``"<array>:hessian"``.
    :param overwrite: allow replacing an existing array of the output name.
    :param return_report: also return ``{num_skipped, num_fallback}``.
    :returns: the mesh carrying the produced ``(n, 9)`` array, or
        ``(mesh, report)`` when ``return_report`` is set.
    :raises ValueError: on an unknown or ``cell_data`` array name, an array
        with more than one component, an unknown method/location, or an
        output-name collision without ``overwrite``.
    """
    _validate(mesh, array, method, location, output, overwrite)

    out = None
    report = None
    try:
        from . import _core

        res = _core.hessian(
            mesh,
            array,
            _METHOD_ALIASES[method],
            _LOCATION_ALIASES[location],
            "" if output is None else str(output),
            bool(overwrite),
        )
        out = res["mesh"]
        report = {"num_skipped": res["num_skipped"], "num_fallback": res["num_fallback"]}
    except (ValueError, TypeError):
        # A genuine user error must not fall through to the composition path,
        # which would either raise something less helpful or silently
        # succeed differently.
        raise
    except Exception:
        out = None

    if out is None:
        out, skipped, fallback = _hessian_py(
            mesh, array, _METHOD_ALIASES[method], _LOCATION_ALIASES[location], output, overwrite
        )
        report = {"num_skipped": skipped, "num_fallback": fallback}

    return (out, report) if return_report else out
