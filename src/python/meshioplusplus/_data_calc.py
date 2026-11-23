"""Derive a new data array from an elementwise expression.

A dependency-free mesh *operation* (not a file format). The geometry is never
modified; a single derived array is added at the chosen location.

The evaluator is a hand-written tokenizer plus recursive-descent parser -- no
external parser library and no evaluation of arbitrary code. It mirrors the C++
implementation in ``src/cpp/src/operations/data_calc.cpp`` token for token, because
Windows CI builds without the native paths and runs entirely on this fallback.
Evaluation here is whole-array numpy, so broadcasting falls out for free.

Grammar::

    expr    := term (('+'|'-') term)*
    term    := unary (('*'|'/') unary)*
    unary   := ('-'|'+') unary | primary
    primary := number | ident | ident '(' expr (',' expr)* ')' | '(' expr ')'

Functions: ``abs``, ``sqrt``, ``min``, ``max``, ``norm``, ``exp``, ``log``,
``log10``, ``sin``, ``cos``, ``tan``, ``asin``, ``acos``, ``atan``, ``atan2``,
``sinh``, ``cosh``, ``tanh``, ``floor``, ``ceil``, ``sign``, ``pow`` and the
tuple functions ``dot``, ``cross``, ``det``, ``trace``, ``transpose``, ``mag``, ``matmul``, ``outer``, ``inv`` and the mesh-aware calls ``grad``,
``div``, ``curl`` (of a point_data array, by name) and ``area``, ``volume``,
``size`` (per cell). Identifiers may
contain ``:`` and ``.`` after the first character (so ``gmsh:physical`` works),
and a name with spaces or operator characters can be backtick-quoted.

Public API:
    data_calc
"""

from __future__ import annotations

import copy

import numpy as np

from ._data_common import available_keys, location_map, normalize_location
from ._fallback import core_op_declined

MAX_COMPONENTS = 16
MAX_DEPTH = 64

#: The coordinate names, and the column of ``(x, y, z)`` each one is.
_COORDS = {"x": 0, "y": 1, "z": 2, "coords": slice(0, 3)}

_FUNCS = {
    "abs": 1,
    "sqrt": 1,
    "min": 2,
    "max": 2,
    "norm": 1,
    "exp": 1,
    "log": 1,
    "log10": 1,
    "sin": 1,
    "cos": 1,
    "tan": 1,
    "asin": 1,
    "acos": 1,
    "atan": 1,
    "atan2": 2,
    "sinh": 1,
    "cosh": 1,
    "tanh": 1,
    "floor": 1,
    "ceil": 1,
    "sign": 1,
    "pow": 2,
    "dot": 2,
    "cross": 2,
    "det": 1,
    "trace": 1,
    "transpose": 1,
    "mag": 1,
    "matmul": 2,
    "outer": 2,
    "inv": 1,
    "grad": 1,
    "div": 1,
    "curl": 1,
    "area": 0,
    "volume": 0,
    "size": 0,
}

#: The mesh-aware calls: derivatives of a point_data array by name, and the
#: per-cell measures. They are evaluated once, before the element loop.
_DERIVATIVES = {"grad": "gradient", "div": "divergence", "curl": "curl"}
_MEASURES = {"area": 2, "volume": 3, "size": 0}

#: The unary math functions beyond ``abs``/``sqrt``/``norm``, width-preserving.
_UNARY = {
    "exp": np.exp,
    "log": np.log,
    "log10": np.log10,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "asin": np.arcsin,
    "acos": np.arccos,
    "atan": np.arctan,
    "sinh": np.sinh,
    "cosh": np.cosh,
    "tanh": np.tanh,
    "floor": np.floor,
    "ceil": np.ceil,
    "sign": np.sign,
}


def _fail(message: str, pos: int | None = None):
    where = "" if pos is None else f" at position {pos}"
    raise ValueError(f"meshio++: data_calc: {message}{where}")


def _tokenize(text: str) -> list:
    """Return ``[(kind, value, pos)]`` with kind in number/ident/op."""
    out = []
    i = 0
    n = len(text)
    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
            continue
        if c.isdigit() or (c == "." and i + 1 < n and text[i + 1].isdigit()):
            j = i
            while j < n and (text[j].isdigit() or text[j] == "."):
                j += 1
            if j < n and text[j] in "eE":
                k = j + 1
                if k < n and text[k] in "+-":
                    k += 1
                if k < n and text[k].isdigit():
                    j = k
                    while j < n and text[j].isdigit():
                        j += 1
            try:
                value = float(text[i:j])
            except ValueError:
                _fail("malformed number", i)
            out.append(("number", value, i))
            i = j
            continue
        if c == "`":
            close = text.find("`", i + 1)
            if close < 0:
                _fail("unterminated `-quoted name", i)
            name = text[i + 1 : close]
            if not name:
                _fail("empty `-quoted name", i)
            out.append(("ident", name, i))
            i = close + 1
            continue
        if c.isalpha() or c == "_":
            j = i + 1
            while j < n and (text[j].isalnum() or text[j] in "_:."):
                j += 1
            while j > i + 1 and text[j - 1] in ":.":
                j -= 1
            out.append(("ident", text[i:j], i))
            i = j
            continue
        if c in "+-*/(),":
            out.append(("op", c, i))
            i += 1
            continue
        _fail(f"unexpected character '{c}'", i)
    out.append(("end", None, n))
    return out


class _Parser:
    """Recursive-descent parser producing a nested tuple tree."""

    def __init__(self, tokens):
        self.tokens = tokens
        self.pos = 0

    def peek(self):
        return self.tokens[self.pos]

    def take(self):
        tok = self.tokens[self.pos]
        self.pos += 1
        return tok

    def parse_all(self):
        node = self.expr(0)
        kind, value, pos = self.peek()
        if kind != "end":
            _fail("trailing input after the expression", pos)
        return node

    def _depth(self, depth):
        if depth > MAX_DEPTH:
            _fail(f"expression nests deeper than {MAX_DEPTH} levels")

    def expr(self, depth):
        self._depth(depth)
        node = self.term(depth + 1)
        while True:
            kind, value, pos = self.peek()
            if kind != "op" or value not in "+-":
                return node
            self.take()
            node = (value, node, self.term(depth + 1), pos)

    def term(self, depth):
        self._depth(depth)
        node = self.unary(depth + 1)
        while True:
            kind, value, pos = self.peek()
            if kind != "op" or value not in "*/":
                return node
            self.take()
            node = (value, node, self.unary(depth + 1), pos)

    def unary(self, depth):
        self._depth(depth)
        kind, value, pos = self.peek()
        if kind == "op" and value == "+":
            self.take()
            return self.unary(depth + 1)
        if kind == "op" and value == "-":
            self.take()
            return ("neg", self.unary(depth + 1), None, pos)
        return self.primary(depth + 1)

    def primary(self, depth):
        self._depth(depth)
        kind, value, pos = self.peek()
        if kind == "end":
            _fail("unexpected end of expression")
        if kind == "number":
            self.take()
            return ("const", value, None, pos)
        if kind == "op" and value == "(":
            self.take()
            inner = self.expr(depth + 1)
            k2, v2, p2 = self.peek()
            if not (k2 == "op" and v2 == ")"):
                _fail("expected ')'", p2)
            self.take()
            return inner
        if kind == "ident":
            self.take()
            k2, v2, _ = self.peek()
            if not (k2 == "op" and v2 == "("):
                return ("array", value, None, pos)
            if value not in _FUNCS:
                _fail(
                    f"unknown function '{value}' (known: {', '.join(_FUNCS)})",
                    pos,
                )
            self.take()
            args = []
            k3, v3, _ = self.peek()
            if not (k3 == "op" and v3 == ")"):
                while True:
                    args.append(self.expr(depth + 1))
                    k4, v4, _ = self.peek()
                    if not (k4 == "op" and v4 == ","):
                        break
                    self.take()
            k5, v5, p5 = self.peek()
            if not (k5 == "op" and v5 == ")"):
                _fail("expected ')'", p5)
            self.take()
            arity = _FUNCS[value]
            if len(args) != arity:
                word = "argument" if arity == 1 else "arguments"
                _fail(
                    f"'{value}' takes exactly {arity} {word} (got {len(args)})",
                    pos,
                )
            if value in _DERIVATIVES:
                # A derivative differentiates a point_data array by name: an
                # expression (and so a nested call) has no such name.
                if args[0][0] != "array":
                    _fail(
                        f"'{value}' needs a point_data array name, not an expression "
                        "(for a second derivative use hessian)",
                        pos,
                    )
                return (value, args[0][1], None, pos)
            if value in _MEASURES:
                return (value, None, None, pos)
            return (value, args[0], args[1] if len(args) > 1 else None, pos)
        _fail("unexpected token", pos)


def _as_2d(values, rows):
    """Reshape to ``(rows, ncomp)``."""
    a = np.asarray(values, dtype=float)
    if a.ndim < 2:
        return a.reshape(rows, 1)
    return a.reshape(rows, -1)


def _combine(a, b, op, pos):
    """Broadcast two ``(rows, ncomp)`` operands, allowing a scalar column."""
    wa = a.shape[1]
    wb = b.shape[1]
    if wa != wb and wa != 1 and wb != 1:
        _fail(
            f"cannot combine a {wa}-component array with a "
            f"{wb}-component array in '{op}'",
            pos,
        )
    return a, b


def _eval(node, arrays, rows):
    """Evaluate ``node`` to a ``(rows, ncomp)`` float array."""
    kind, lhs, rhs, pos = node
    if kind == "const":
        return np.full((rows, 1), float(lhs))
    if kind == "array":
        # Names are resolved and reported before evaluation starts, so a miss
        # here would be an internal inconsistency rather than user error.
        return arrays[lhs]
    if kind == "neg":
        return -_eval(lhs, arrays, rows)
    if kind in _DERIVATIVES:
        return arrays[(kind, lhs)]
    if kind in _MEASURES:
        return arrays[(kind, None)]
    if kind == "abs":
        return np.abs(_eval(lhs, arrays, rows))
    if kind == "sqrt":
        with np.errstate(invalid="ignore"):
            return np.sqrt(_eval(lhs, arrays, rows))
    if kind in _UNARY:
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            return _UNARY[kind](_eval(lhs, arrays, rows))
    if kind in ("norm", "mag"):
        v = _eval(lhs, arrays, rows)
        return np.sqrt((v * v).sum(axis=1)).reshape(rows, 1)
    if kind == "inv":
        v = _eval(lhs, arrays, rows)
        if v.shape[1] not in (6, 9):
            _fail(
                "'inv' needs a 6-component (symmetric) or 9-component (3x3) "
                f"operand (got {v.shape[1]})",
                pos,
            )
        from ._tensor_invariants import _inverse

        return _inverse(v, v.shape[1])
    if kind in ("det", "trace", "transpose"):
        v = _eval(lhs, arrays, rows)
        if v.shape[1] not in (6, 9) or (kind == "transpose" and v.shape[1] != 9):
            need = (
                "a 9-component (3x3) operand"
                if kind == "transpose"
                else "a 6-component (symmetric) or 9-component (3x3) operand"
            )
            _fail(f"'{kind}' needs {need} (got {v.shape[1]})", pos)
        if v.shape[1] == 6:
            m = [
                v[:, 0],
                v[:, 3],
                v[:, 5],
                v[:, 3],
                v[:, 1],
                v[:, 4],
                v[:, 5],
                v[:, 4],
                v[:, 2],
            ]
        else:
            m = [v[:, k] for k in range(9)]
        if kind == "trace":
            return (m[0] + m[4] + m[8]).reshape(rows, 1)
        if kind == "det":
            det = (
                m[0] * (m[4] * m[8] - m[5] * m[7])
                + m[1] * (m[5] * m[6] - m[3] * m[8])
                + m[2] * (m[3] * m[7] - m[4] * m[6])
            )
            return det.reshape(rows, 1)
        return np.stack([m[j * 3 + i] for i in range(3) for j in range(3)], axis=1)
    a = _eval(lhs, arrays, rows)
    b = _eval(rhs, arrays, rows)
    if kind == "matmul":
        if a.shape[1] != 9 or b.shape[1] not in (3, 9):
            _fail(
                "'matmul' needs a 9-component (3x3) operand and a 9-component "
                f"(3x3) or 3-component (vector) operand (got {a.shape[1]} and {b.shape[1]})",
                pos,
            )
        am = a.reshape(rows, 3, 3)
        cols = 3 if b.shape[1] == 9 else 1
        bm = b.reshape(rows, 3, cols)
        out = np.zeros((rows, 3, cols))
        for k in range(3):  # k ascending, as the core sums
            out = out + am[:, :, k : k + 1] * bm[:, k : k + 1, :]
        return out.reshape(rows, 3 * cols)
    if kind == "outer":
        if a.shape[1] != 3 or b.shape[1] != 3:
            _fail(
                f"'outer' needs two 3-component operands (got {a.shape[1]} and {b.shape[1]})",
                pos,
            )
        return (a[:, :, None] * b[:, None, :]).reshape(rows, 9)
    if kind in ("dot", "cross"):
        if kind == "dot" and a.shape[1] != b.shape[1]:
            _fail(
                f"'dot' needs two operands of the same width (got {a.shape[1]} and {b.shape[1]})",
                pos,
            )
        if kind == "cross" and (a.shape[1] != 3 or b.shape[1] != 3):
            _fail(
                f"'cross' needs two 3-component operands (got {a.shape[1]} and {b.shape[1]})",
                pos,
            )
        if kind == "dot":
            return (a * b).sum(axis=1).reshape(rows, 1)
        return np.stack(
            [
                a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1],
                a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2],
                a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0],
            ],
            axis=1,
        )
    a, b = _combine(a, b, kind, pos)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        if kind == "atan2":
            return np.arctan2(a, b)
        if kind == "pow":
            return np.power(a, b)
        if kind == "+":
            return a + b
        if kind == "-":
            return a - b
        if kind == "*":
            return a * b
        if kind == "/":
            return a / b
        if kind == "min":
            return np.minimum(a, b)
        if kind == "max":
            return np.maximum(a, b)
    _fail(f"unsupported operator '{kind}'", pos)


def _collect_names(node, out):
    """Gather every array name the tree references."""
    kind, lhs, rhs, _ = node
    if kind == "array":
        out.add(lhs)
        return
    if kind == "const":
        return
    if isinstance(lhs, tuple):
        _collect_names(lhs, out)
    if isinstance(rhs, tuple):
        _collect_names(rhs, out)


def _collect_calls(node, out):
    """Gather every ``(kind, array name or None, position)`` mesh-aware call."""
    kind, lhs, rhs, pos = node
    if kind in _DERIVATIVES or kind in _MEASURES:
        out.add((kind, lhs, pos))
        return
    if isinstance(lhs, tuple):
        _collect_calls(lhs, out)
    if isinstance(rhs, tuple):
        _collect_calls(rhs, out)


def _precompute_calls(mesh, tree, loc):
    """Evaluate the mesh-aware calls once: ``{(kind, name): [array per block]}``
    (a single array at the point location)."""
    from ._data_average import _DIM, _cell_measures
    from ._gradient import _gradient_py, _validate

    calls = set()
    _collect_calls(tree, calls)
    out = {}
    for kind, name, pos in sorted(calls, key=lambda c: (c[0], str(c[1]), c[2])):
        if loc == "field_data":
            _fail(
                f"'{kind}' needs the point or cell location: field_data has no mesh to read",
                pos,
            )
        if kind in _DERIVATIVES:
            if name not in mesh.point_data:
                _fail(
                    f"'{kind}' differentiates a point_data array and '{name}' is not one",
                    pos,
                )
            where = "point" if loc == "point_data" else "cell"
            op = _DERIVATIVES[kind]
            try:
                _validate(mesh, name, op, "green-gauss", where, "\x01calc", None, True)
            except ValueError as exc:
                _fail(f"'{kind}': {exc}", pos)
            res, _skipped, _fallback = _gradient_py(
                mesh, name, op, "green-gauss", where, "\x01calc", None, True
            )
            if loc == "point_data":
                out[(kind, name)] = [np.asarray(res.point_data["\x01calc"])]
            else:
                out[(kind, name)] = [np.asarray(b) for b in res.cell_data["\x01calc"]]
            continue
        if loc != "cell_data":
            _fail(f"'{kind}' is a per-cell quantity: use the cell location", pos)
        want = _MEASURES[kind]
        blocks = []
        for block in mesh.cells:
            n = len(block.data)
            ctype = block.type.split("_")[0] if "_" in block.type else block.type
            if block.type.startswith("polyhedron"):
                dim = 3
            elif block.type.startswith("polygon") or isinstance(block.data, list):
                dim = 2
            else:
                dim = _DIM.get(ctype)
            meas = _cell_measures(mesh, block)
            if meas is None or (want != 0 and dim != want):
                meas = np.full(n, np.nan)
            blocks.append(np.asarray(meas, dtype=float))
        out[(kind, None)] = blocks
    return out


def _calc_py(mesh, expression, location, output, overwrite):
    """Pure-Python reference for :func:`data_calc`."""
    loc = normalize_location(location)
    if not output:
        _fail("an output array name is required")
    if not overwrite and output in location_map(mesh, loc):
        _fail(
            f"output name '{output}' already exists in {loc} "
            "(pass overwrite=True to replace it)"
        )

    tokens = _tokenize(expression)
    if len(tokens) == 1:
        _fail("the expression is empty")
    tree = _Parser(tokens).parse_all()

    referenced = set()
    _collect_names(tree, referenced)
    source = location_map(mesh, loc)
    # A name that no data array has but that is `x`, `y`, `z` or `coords` is a
    # coordinate (a point's, or a cell's centroid); an array of that name wins.
    coord_names = (
        {n for n in referenced if n in _COORDS and n not in source}
        if loc != "field_data"
        else set()
    )
    for name in sorted(referenced - coord_names):
        if name not in source:
            keys = available_keys(mesh, loc)
            detail = (
                f"available: {', '.join(keys)}" if keys else f"the mesh has no {loc}"
            )
            _fail(f"unknown {loc} array '{name}' ({detail})")

    out = copy.deepcopy(mesh)
    call_values = _precompute_calls(mesh, tree, loc)

    def run(rows, fetch, block=0):
        arrays = {n: _as_2d(fetch(n), rows) for n in referenced - coord_names}
        for key, values in call_values.items():
            arrays[key] = _as_2d(values[min(block, len(values) - 1)], rows)
        for n in coord_names:
            xyz = fetch(n)
            arrays[n] = xyz[:, _COORDS[n]] if n != "coords" else xyz
            arrays[n] = _as_2d(arrays[n], rows)
        for n, a in arrays.items():
            if a.shape[1] > MAX_COMPONENTS:
                _fail(
                    f"array '{n}' has {a.shape[1]} components, more than the "
                    f"supported maximum of {MAX_COMPONENTS}"
                )
        res = _eval(tree, arrays, rows)
        return res.reshape(rows) if res.shape[1] == 1 else res

    if loc == "point_data":
        rows = len(mesh.points)
        for name in referenced - coord_names:
            got = len(np.asarray(source[name]))
            if got != rows:
                _fail(f"array '{name}' has {got} rows but point_data needs {rows}")
        xyz = None
        if coord_names:
            pts = np.asarray(mesh.points, dtype=float)
            xyz = np.zeros((rows, 3))
            xyz[:, : min(pts.shape[1], 3)] = pts[:, :3]
        out.point_data[output] = run(
            rows, lambda n: xyz if n in coord_names else source[n]
        )
        return out

    if loc == "field_data":
        rows = None
        for name in sorted(referenced):
            got = len(np.asarray(source[name]).reshape(-1))
            rows = got if rows is None else rows
            if got != rows:
                _fail(f"field_data arrays disagree in length ({rows} vs {got})")
        out.field_data[output] = run(rows or 0, lambda n: source[n])
        return out

    nblocks = len(mesh.cells)
    for name in referenced - coord_names:
        if len(source[name]) != nblocks:
            _fail(
                f"cell_data '{name}' has {len(source[name])} block(s) but the "
                f"mesh has {nblocks} cell block(s)"
            )
    centroids = None
    if coord_names:
        from ._partition import _centroids

        total = sum(len(block.data) for block in mesh.cells)
        centroids = _centroids(mesh, total)
    blocks = []
    base = 0
    for b, block in enumerate(mesh.cells):
        rows = len(block.data)
        cent = centroids[base : base + rows] if centroids is not None else None
        blocks.append(
            run(rows, lambda n, b=b, c=cent: c if n in coord_names else source[n][b], b)
        )
        base += rows
    out.cell_data[output] = blocks
    return out


def data_calc(
    mesh,
    expression: str,
    location: str = "point",
    output: str = "",
    overwrite: bool = False,
):
    """Evaluate ``expression`` elementwise and store it as a new array.

    Args:
        mesh: the source mesh (unmodified).
        expression: the expression text (see the module docstring's grammar).
        location: ``"point"``, ``"cell"`` or ``"field"`` -- where operands are
            looked up and the result is stored.
        output: name of the new array. Required.
        overwrite: whether an existing array of that name may be replaced.

    Returns:
        A new mesh carrying the derived array.

    Raises:
        ValueError: any lexical, syntactic, name-resolution, arity,
            component-width or row-count error, prefixed
            ``meshio++: data_calc: `` and carrying the character position.
    """
    loc = normalize_location(location)
    try:
        from . import _core
    except ImportError:
        return _calc_py(mesh, expression, loc, output, overwrite)
    try:
        out = _core.data_calc(mesh, expression, loc, output, overwrite)
    except Exception as exc:
        if not core_op_declined(exc, "data_calc"):
            raise
        return _calc_py(mesh, expression, loc, output, overwrite)
    for attr in ("point_sets", "cell_sets"):
        value = getattr(mesh, attr, None)
        if value:
            setattr(out, attr, copy.deepcopy(value))
    return out
