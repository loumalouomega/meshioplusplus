"""``tensor_invariants``: von Mises, principal, hydrostatic and deviatoric
fields of a symmetric (6-component) or general 3x3 (9-component) tensor array.

A dependency-free mesh *operation* (not a file format). Before this operation
the only route to these was the CalculiX ``.frd`` reader's format-specific
``derived=`` option (:mod:`meshioplusplus.frd`); it now calls
:func:`tensor_invariants` internally and keeps its existing output names.

Input layout, per row:

- 6 components: a symmetric tensor ``xx yy zz xy yz zx`` (the order used
  throughout this codebase; see ``doc/mesh_data_model.md``).
- 9 components: a general 3x3 tensor, row-major (``xx xy xz yx yy yz zx zy
  zz``), the layout :func:`meshioplusplus.gradient`/:func:`meshioplusplus.hessian`
  produce. ``mises`` and ``principal`` are only defined for a symmetric
  tensor, so both use the symmetric part ``0.5*(T + T^T)``.

A row with any non-finite input component produces NaN in every requested
output for that row.

Public API:
    tensor_invariants
"""

from __future__ import annotations

import copy

import numpy as np

from ._data_common import location_map, normalize_location, num_components, require_key
from ._fallback import core_op_declined

#: What an unqualified call computes: the four original invariants.
_ALL_OUTPUTS = ("mises", "principal", "hydrostatic", "deviatoric")
#: Every accepted name; the algebraic ones read the tensor as given and are opt-in.
_VALID_OUTPUTS = _ALL_OUTPUTS + ("trace", "determinant", "invariants", "inverse")


def _symmetric6(flat, ncomp):
    """Rows reduced to ``xx yy zz xy yz zx``, symmetrizing a 9-comp input."""
    if ncomp == 6:
        return flat
    # Row-major 3x3: xx xy xz yx yy yz zx zy zz.
    xx, xy, xz, yx, yy, yz, zx, zy, zz = (flat[:, k] for k in range(9))
    return np.stack(
        [xx, yy, zz, 0.5 * (xy + yx), 0.5 * (yz + zy), 0.5 * (xz + zx)], axis=1
    )


def _matrix9(flat, ncomp):
    """Rows as the nine entries of the row-major 3x3, as given (not symmetrised)."""
    if ncomp == 9:
        return [flat[:, k] for k in range(9)]
    xx, yy, zz, xy, yz, xz = (flat[:, k] for k in range(6))
    return [xx, xy, xz, xy, yy, yz, xz, yz, zz]


def _finite_rows(flat):
    return np.isfinite(flat).all(axis=1)


def _det9(m):
    """The triple product of the rows, in the arithmetic order of the core."""
    return (
        m[0] * (m[4] * m[8] - m[5] * m[7])
        + m[1] * (m[5] * m[6] - m[3] * m[8])
        + m[2] * (m[3] * m[7] - m[4] * m[6])
    )


def _trace(flat, ncomp):
    m = _matrix9(flat, ncomp)
    return np.where(_finite_rows(flat), m[0] + m[4] + m[8], np.nan)


def _determinant(flat, ncomp):
    m = _matrix9(flat, ncomp)
    return np.where(_finite_rows(flat), _det9(m), np.nan)


def _invariants(flat, ncomp):
    m = _matrix9(flat, ncomp)
    i1 = m[0] + m[4] + m[8]
    i2 = (
        (m[0] * m[4] - m[1] * m[3])
        + (m[4] * m[8] - m[5] * m[7])
        + (m[0] * m[8] - m[2] * m[6])
    )
    out = np.stack([i1, i2, _det9(m)], axis=1)
    return np.where(_finite_rows(flat)[:, None], out, np.nan)


def _inverse(flat, ncomp):
    """``adj(A) / det(A)`` in the input's width; NaN for a zero determinant."""
    m = _matrix9(flat, ncomp)
    det = np.where(_finite_rows(flat), _det9(m), np.nan)
    bad = ~np.isfinite(det) | (det == 0.0)
    safe = np.where(bad, 1.0, det)
    inv = [
        (m[4] * m[8] - m[5] * m[7]) / safe,
        (m[2] * m[7] - m[1] * m[8]) / safe,
        (m[1] * m[5] - m[2] * m[4]) / safe,
        (m[5] * m[6] - m[3] * m[8]) / safe,
        (m[0] * m[8] - m[2] * m[6]) / safe,
        (m[2] * m[3] - m[0] * m[5]) / safe,
        (m[3] * m[7] - m[4] * m[6]) / safe,
        (m[1] * m[6] - m[0] * m[7]) / safe,
        (m[0] * m[4] - m[1] * m[3]) / safe,
    ]
    if ncomp == 6:
        cols = [inv[0], inv[4], inv[8], inv[1], inv[5], inv[2]]
    else:
        cols = inv
    out = np.stack(cols, axis=1)
    return np.where(bad[:, None], np.nan, out)


_ALGEBRAIC = {
    "trace": _trace,
    "determinant": _determinant,
    "invariants": _invariants,
    "inverse": _inverse,
}


def _mises(sym6):
    xx, yy, zz, xy, yz, xz = (sym6[:, k] for k in range(6))
    return np.sqrt(
        0.5
        * (
            (xx - yy) ** 2
            + (yy - zz) ** 2
            + (zz - xx) ** 2
            + 6.0 * (xy**2 + yz**2 + xz**2)
        )
    )


def _principal(sym6):
    """Eigenvalues of the symmetric tensors, ascending (min, mid, max)."""
    out = np.full((len(sym6), 3), np.nan)
    ok = np.isfinite(sym6).all(axis=1)
    if ok.any():
        s = sym6[ok]
        m = np.empty((len(s), 3, 3))
        m[:, 0, 0], m[:, 1, 1], m[:, 2, 2] = s[:, 0], s[:, 1], s[:, 2]
        m[:, 0, 1] = m[:, 1, 0] = s[:, 3]
        m[:, 1, 2] = m[:, 2, 1] = s[:, 4]
        m[:, 0, 2] = m[:, 2, 0] = s[:, 5]
        out[ok] = np.linalg.eigvalsh(m)
    return out


def _hydrostatic(flat, ncomp):
    """Mean of the tensor's own (unsymmetrized) diagonal."""
    yy_idx = 1 if ncomp == 6 else 4
    zz_idx = 2 if ncomp == 6 else 8
    xx, yy, zz = flat[:, 0], flat[:, yy_idx], flat[:, zz_idx]
    with np.errstate(invalid="ignore"):
        return (xx + yy + zz) / 3.0


def _deviatoric(flat, ncomp, hydro):
    yy_idx = 1 if ncomp == 6 else 4
    zz_idx = 2 if ncomp == 6 else 8
    out = flat.copy()
    ok = np.isfinite(hydro)
    for idx in (0, yy_idx, zz_idx):
        out[:, idx] = np.where(
            ok, flat[:, idx] - np.where(ok, hydro, 0.0), flat[:, idx]
        )
    return out


def _tensor_invariants_py(mesh, location, names, outputs, prefix, suffix, overwrite):
    """Pure-Python reference for :func:`tensor_invariants`."""
    loc = normalize_location(location)
    if loc == "field_data":
        raise ValueError(
            "meshio++: tensor_invariants: field_data has no per-row tensor to reduce"
        )
    source = location_map(mesh, loc)
    outs = list(outputs) if outputs else list(_ALL_OUTPUTS)
    for o in outs:
        if o not in _VALID_OUTPUTS:
            raise ValueError(
                f"meshio++: unknown tensor invariant '{o}' (expected one of "
                f"{', '.join(_VALID_OUTPUTS)})"
            )

    if names:
        selected = list(names)
        for name in selected:
            require_key(mesh, loc, name)
            ncomp = num_components(
                source[name][0] if loc == "cell_data" else source[name]
            )
            if ncomp not in (6, 9):
                raise ValueError(
                    f"meshio++: tensor_invariants: '{name}' has {ncomp} "
                    "component(s); expected 6 (symmetric) or 9 (general 3x3)"
                )
    else:
        selected = [
            name
            for name in sorted(source.keys())
            if num_components(source[name][0] if loc == "cell_data" else source[name])
            in (6, 9)
        ]

    out = copy.deepcopy(mesh)

    def _store(target, loc_, blocks_or_array):
        if not overwrite and target in location_map(out, loc_):
            raise ValueError(
                f"meshio++: tensor_invariants: '{target}' already exists (overwrite=false)"
            )
        location_map(out, loc_)[target] = blocks_or_array

    for name in selected:
        base = prefix + name
        if loc == "cell_data":
            blocks = source[name]
            if len(blocks) != len(mesh.cells):
                raise ValueError(
                    f"meshio++: tensor_invariants: cell_data '{name}' has "
                    f"{len(blocks)} block(s) but the mesh has "
                    f"{len(mesh.cells)} cell block(s)"
                )
            ncomp = num_components(blocks[0])
            per_output = {o: [] for o in outs}
            for b in blocks:
                flat = np.asarray(b, dtype=float).reshape(-1, ncomp)
                sym6 = _symmetric6(flat, ncomp)
                hydro = (
                    _hydrostatic(flat, ncomp)
                    if "hydrostatic" in outs or "deviatoric" in outs
                    else None
                )
                if "mises" in outs:
                    per_output["mises"].append(_mises(sym6))
                if "principal" in outs:
                    per_output["principal"].append(_principal(sym6))
                if "hydrostatic" in outs:
                    per_output["hydrostatic"].append(hydro)
                if "deviatoric" in outs:
                    per_output["deviatoric"].append(
                        _deviatoric(flat, ncomp, hydro).reshape(np.asarray(b).shape)
                    )
                for o in outs:
                    if o in _ALGEBRAIC:
                        value = _ALGEBRAIC[o](flat, ncomp)
                        if o == "inverse":
                            value = value.reshape(np.asarray(b).shape)
                        per_output[o].append(value)
            for o in outs:
                _store(base + f"_{o}" + suffix, "cell_data", per_output[o])
            continue

        arr = np.asarray(source[name], dtype=float)
        ncomp = num_components(arr)
        flat = arr.reshape(-1, ncomp)
        sym6 = _symmetric6(flat, ncomp)
        if "mises" in outs:
            _store(base + "_mises" + suffix, loc, _mises(sym6))
        if "principal" in outs:
            _store(base + "_principal" + suffix, loc, _principal(sym6))
        hydro = None
        if "hydrostatic" in outs or "deviatoric" in outs:
            hydro = _hydrostatic(flat, ncomp)
        if "hydrostatic" in outs:
            _store(base + "_hydrostatic" + suffix, loc, hydro)
        if "deviatoric" in outs:
            _store(
                base + "_deviatoric" + suffix,
                loc,
                _deviatoric(flat, ncomp, hydro).reshape(arr.shape),
            )
        for o in outs:
            if o in _ALGEBRAIC:
                value = _ALGEBRAIC[o](flat, ncomp)
                if o == "inverse":
                    value = value.reshape(arr.shape)
                _store(base + f"_{o}" + suffix, loc, value)
    return out


def tensor_invariants(
    mesh,
    location: str = "point",
    keys=None,
    outputs=None,
    prefix: str = "",
    suffix: str = "",
    overwrite: bool = True,
):
    """Compute von Mises / principal / hydrostatic / deviatoric (and algebraic) fields.

    Args:
        mesh: the source mesh (unmodified).
        location: ``"point"`` or ``"cell"`` (``"field"`` is rejected).
        keys: names to process; ``None`` means every 6- or 9-component array
            at ``location``.
        outputs: any of ``"mises"``, ``"principal"``, ``"hydrostatic"``,
            ``"deviatoric"``, ``"trace"``, ``"determinant"``, ``"invariants"``
            (``I1, I2, I3``, three components) and ``"inverse"`` (the input's
            width; NaN where the determinant is zero); ``None`` means the first
            four. The last four read the tensor as given: a 9-component input is
            not symmetrised, unlike ``mises`` and ``principal``.
        prefix: prepended to every output array's name.
        suffix: appended after the invariant's own name segment.
        overwrite: throw instead of silently overwriting an existing array of
            the target name.

    Returns:
        A new mesh with the source data plus the requested invariant arrays;
        geometry is unchanged.
    """
    loc = normalize_location(location)
    names = list(keys) if keys else []
    outs = list(outputs) if outputs else []
    args = (loc, names, outs, prefix, suffix, overwrite)
    try:
        from . import _core
    except ImportError:
        return _tensor_invariants_py(mesh, *args)
    try:
        out = _core.tensor_invariants(
            mesh, loc, names, ",".join(outs), prefix, suffix, overwrite
        )
    except Exception as exc:
        if not core_op_declined(exc, "tensor_invariants"):
            raise
        return _tensor_invariants_py(mesh, *args)
    for attr in ("point_sets", "cell_sets"):
        value = getattr(mesh, attr, None)
        if value:
            setattr(out, attr, copy.deepcopy(value))
    return out
