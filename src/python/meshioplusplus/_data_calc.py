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

Functions: ``abs``, ``sqrt``, ``min``, ``max``, ``norm``. Identifiers may
contain ``:`` and ``.`` after the first character (so ``gmsh:physical`` works),
and a name with spaces or operator characters can be backtick-quoted.

Public API:
    data_calc
"""

from __future__ import annotations

import copy

import numpy as np

from ._data_common import available_keys, location_map, normalize_location

MAX_COMPONENTS = 16
MAX_DEPTH = 64

_FUNCS = {"abs": 1, "sqrt": 1, "min": 2, "max": 2, "norm": 1}


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
    if kind == "abs":
        return np.abs(_eval(lhs, arrays, rows))
    if kind == "sqrt":
        with np.errstate(invalid="ignore"):
            return np.sqrt(_eval(lhs, arrays, rows))
    if kind == "norm":
        v = _eval(lhs, arrays, rows)
        return np.sqrt((v * v).sum(axis=1)).reshape(rows, 1)
    a = _eval(lhs, arrays, rows)
    b = _eval(rhs, arrays, rows)
    a, b = _combine(a, b, kind, pos)
    with np.errstate(divide="ignore", invalid="ignore"):
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
    for name in sorted(referenced):
        if name not in source:
            keys = available_keys(mesh, loc)
            detail = (
                f"available: {', '.join(keys)}" if keys else f"the mesh has no {loc}"
            )
            _fail(f"unknown {loc} array '{name}' ({detail})")

    out = copy.deepcopy(mesh)

    def run(rows, fetch):
        arrays = {n: _as_2d(fetch(n), rows) for n in referenced}
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
        for name in referenced:
            got = len(np.asarray(source[name]))
            if got != rows:
                _fail(f"array '{name}' has {got} rows but point_data needs {rows}")
        out.point_data[output] = run(rows, lambda n: source[n])
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
    for name in referenced:
        if len(source[name]) != nblocks:
            _fail(
                f"cell_data '{name}' has {len(source[name])} block(s) but the "
                f"mesh has {nblocks} cell block(s)"
            )
    blocks = []
    for b, block in enumerate(mesh.cells):
        rows = len(block.data)
        blocks.append(run(rows, lambda n, b=b: source[n][b]))
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
    except Exception:
        return _calc_py(mesh, expression, loc, output, overwrite)
    try:
        out = _core.data_calc(mesh, expression, loc, output, overwrite)
    except ValueError:
        raise
    except Exception:
        return _calc_py(mesh, expression, loc, output, overwrite)
    for attr in ("point_sets", "cell_sets"):
        value = getattr(mesh, attr, None)
        if value:
            setattr(out, attr, copy.deepcopy(value))
    return out
