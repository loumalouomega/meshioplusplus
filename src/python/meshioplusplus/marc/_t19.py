"""Write MSC Marc formatted post files (``.t19``).

The Python twin of ``write_marc_t19`` in ``src/cpp/src/formats/marc.cpp``: both
write the same bytes, in the layout ``read_t19`` (``_marc.py``) reads: the
analysis title, the model blocks 502 to 515 (header and element post codes,
connectivities, coordinates, nodal codes, sets, flow lines) and one group of
blocks per increment between ``****`` and ``----`` (title, integer and real
verification data, element integration point values, nodal results), closed by
``++++``. Every integer is ``i13`` and every real ``e13.6``, six to a line, each
Fortran record starting a new line.

The model is the input deck writer's (``deck_model``): the cells that have a
Marc element type, ``marc:element`` numbers, point and cell regions as node and
element sets. Point arrays of one, two or three components per point are nodal
results (a vector is three components); cell arrays are element post codes, per
integration point when ``field_data["marc:layout:<name>"]`` says so.
"""

from __future__ import annotations

import numpy as np

from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from ._marc import _SCALAR_NAMES, _TENSOR_NAMES, TYPES, deck_model

__all__ = ["write", "SeriesWriter"]

_W = 13
_DEFAULT_SCALAR = {name: code for code, name in _SCALAR_NAMES.items()}
_DEFAULT_TENSOR = {name: code for code, name in _TENSOR_NAMES.items()}
_FREE_TENSOR_BASES = (541, 461, 421, 431, 441, 411, 391, 381, 371, 361, 351)
_FIRST_FREE_SCALAR = 600
_LABEL = 24
_NODAL_IDS = {"Displacement": 1, "Reaction Force": 5}


def _e13(x):
    """Fortran ``e13.6``: `` 0.123457E+01`` (NaN and infinities as text)."""
    x = float(x)
    if x != x:
        return f"{'NaN':>{_W}}"
    if x in (float("inf"), float("-inf")):
        return f"{'inf' if x > 0 else '-inf':>{_W}}"
    if x == 0.0:
        return " 0.000000E+00"
    mant, exp = f"{abs(x):.5E}".split("E")
    digits = mant.replace(".", "")
    e = int(exp) + 1
    sign = "-" if x < 0 else " "
    if abs(e) > 99:
        return f"{sign}0.{digits}{e:+04d}"
    return f"{sign}0.{digits}E{e:+03d}"


def _lines(values, per_line, fmt):
    return [
        "".join(fmt(v) for v in values[k : k + per_line])
        for k in range(0, len(values), per_line)
    ]


def _ints(values):
    return _lines([int(v) for v in values], 6, lambda v: f"{v:{_W}d}")


def _reals(values):
    return _lines([float(v) for v in values], 6, _e13)


def _block(number, name, body):
    return [f"=beg={number} ({name})".ljust(70)] + body + ["=end="]


def _text_line(text, width):
    return str(text).encode("latin-1", "replace").decode("latin-1")[:width].ljust(width)


class _Codes:
    """The element post codes of a mesh's cell arrays and where each array's
    columns sit in a record."""

    def __init__(self, mesh, model):
        self.entries = []  # (code, label)
        self.arrays = []  # (name, data per written cell: (points, comps))
        self.nstres = 1
        layouts = {}
        for name in mesh.field_data:
            if name.startswith("marc:layout:"):
                layouts[name[len("marc:layout:") :]] = np.asarray(
                    mesh.field_data[name]
                ).reshape(-1)
        stres = {int(v[0]) for v in layouts.values() if len(v) >= 2}
        if len(stres) > 1:
            warn("Marc .t19: arrays with different integration point counts; writing 1")
        elif stres:
            self.nstres = stres.pop()
        written = model.written
        scalar_code = _FIRST_FREE_SCALAR
        tensor_free = list(_FREE_TENSOR_BASES)
        used_codes = set()
        for name in sorted(mesh.cell_data):
            if name in ("marc:element", "marc:type") or name.startswith("tecplot:"):
                continue
            parts = mesh.cell_data[name]
            if len(parts) != len(mesh.cells):
                continue
            rows = []
            for c in written:
                b, r, _ = model.cells[c]
                rows.append(np.asarray(parts[b])[r])
            data = np.asarray(rows, dtype=np.float64)
            layout = layouts.get(name)
            if (
                layout is not None
                and len(layout) >= 2
                and int(layout[0]) == self.nstres
            ):
                comps = int(layout[1])
                if data.ndim != 2 or data.shape[1] != self.nstres * comps:
                    warn(
                        f"Marc .t19: cell data '{name}' does not match its layout; dropped"
                    )
                    continue
                values = data.reshape(len(rows), self.nstres, comps)
            elif self.nstres == 1 and data.ndim <= 2:
                values = data.reshape(len(rows), 1, -1)
                comps = values.shape[2]
            else:
                warn(
                    f"Marc .t19: cell data '{name}' has no integration point layout; dropped"
                )
                continue
            base_name, layer = name, 0
            if "@layer" in name:
                head, _, tail = name.rpartition("@layer")
                if tail.isdigit():
                    base_name, layer = head, int(tail)
            codes = None
            if comps == 6:
                if base_name in _DEFAULT_TENSOR:
                    base = _DEFAULT_TENSOR[base_name]
                    codes = [(layer * 1000 + base + k, "") for k in range(6)]
                elif tensor_free:
                    base = tensor_free.pop(0)
                    codes = [(layer * 1000 + base + k, "") for k in range(6)]
                    codes[0] = (codes[0][0], base_name[:_LABEL])
                else:
                    codes = None
            if codes is None and comps == 1:
                if base_name in _DEFAULT_SCALAR:
                    codes = [(layer * 1000 + _DEFAULT_SCALAR[base_name], "")]
                else:
                    codes = [(layer * 1000 + scalar_code, base_name[:_LABEL])]
                    scalar_code += 1
            if codes is None:
                # Any other width: one scalar code per component.
                codes = []
                for k in range(comps):
                    codes.append(
                        (layer * 1000 + scalar_code, f"{base_name}_{k}"[:_LABEL])
                    )
                    scalar_code += 1
            if any(code in used_codes for code, _ in codes):
                warn(f"Marc .t19: cell data '{name}' repeats a post code; dropped")
                continue
            used_codes.update(code for code, _ in codes)
            self.entries.extend(codes)
            self.arrays.append((name, values))
        self.columns = len(self.entries)

    def record(self, element, point):
        out = []
        for _, values in self.arrays:
            out.extend(values[element, point, :].tolist())
        return out


def _model(mesh, model, codes):
    """Blocks 502 to 515."""
    points, npts = model.points, model.npts
    elements = []
    for c in model.written:
        b, r, etype = model.cells[c]
        block = mesh.cells[b]
        nodes = [int(v) + 1 for v in np.asarray(block.data[r])]
        kind = TYPES[etype][0]
        if kind == "hexahedron" and block.type != "hexahedron":
            from ..lsdyna._lsdyna import _expand_solid

            nodes = _expand_solid(block.type, nodes)
        elif 168 <= etype <= 170 and len(nodes) == 3:
            nodes = [nodes[0], nodes[2], nodes[1]]
        elements.append((model.elem_id[c], etype, nodes))
    nnodmx = max([len(n) for _, _, n in elements] + [1])
    sets = [
        (name, 1 if fam == "node" else 0, members, None)
        for fam, name, members in model.sets
    ]
    for fam, name, items in model.facet_sets:
        pairs = [tuple(int(v) for v in it.strip().split(":")) for it in items]
        sets.append(
            (
                name,
                12 if fam == "edge" else 13,
                [p[0] for p in pairs],
                [p[1] for p in pairs],
            )
        )
    lm = [
        codes.columns,
        npts,
        len(elements),
        3,
        codes.nstres,
        9,
        1,
        0,
        3,
        nnodmx,
        2,
        0,
        0,
        12,
        0,
        len(sets),
        0,
        0,
    ]
    lm += [0] * 12
    out = []
    out += _block(50200, "Analysis Verification Data", _ints(lm))
    out += _block(50400, "Dummy", _ints([0]))
    out += _block(50500, "Domain Decomposition Information", _ints([0, 0]))
    out += _block(
        50600,
        "Element Variable Postcodes",
        [f"{code:{_W}d}" + label.ljust(24) for code, label in codes.entries],
    )
    body = []
    for ident, etype, conn in elements:
        body += _ints([ident, etype, len(conn)] + conn + [0] * (nnodmx - len(conn)))
    out += _block(50700, "Element Connectivities", body)
    body = [
        f"{p + 1:{_W}d}"
        + "".join(
            _e13(float(points[p, d]) if d < points.shape[1] else 0.0) for d in range(3)
        )
        for p in range(npts)
    ]
    out += _block(50800, "Nodal Coordinates", body)
    out += _block(51000, "Nodal Codes and Transformation ID", _ints([3] * npts))
    body = _ints([len(sets)])
    for name, kind, members, numbers in sets:
        body.append(_text_line(name, 32))
        body += _ints([len(members), kind])
        body += _ints(members)
        if numbers is not None:
            body += _ints(numbers)
    out += _block(51301, "Set Definitions", body)
    out += _block(51501, "Flow Line Data", _ints([0] * 6))
    return out


def _nodal_results(mesh, npts):
    """(name, id, components, values) of the point arrays a post file holds."""
    out = []
    for name in sorted(mesh.point_data):
        a = np.asarray(mesh.point_data[name], dtype=np.float64)
        if a.ndim == 1:
            a = a[:, None]
        if a.ndim != 2 or a.shape[0] != npts or a.shape[1] not in (1, 2, 3, 6, 9):
            warn(f"Marc .t19: point data '{name}' is not written")
            continue
        out.append((name, a))
    return out


def _increment(index, time, mesh, model, codes, nodal, new_model=None):
    """One increment; ``new_model`` (blocks 502 to 515 of a remeshed mesh)
    makes it a remeshing increment: ``newmo`` is set and the model follows
    block 518 under block 519."""
    newmo = 0 if new_model is None else 1
    out = ["****"]
    out += _block(51600, "Loadcase Title", ["          meshio++".ljust(70)])
    out += _block(
        51701,
        "Integer Increment Verification Data",
        _ints(
            [newmo, index + 1, 0, 102 if codes.columns else 1, 2, 0, 0, 1, 0, 0, 0, 0]
        ),
    )
    out += _block(
        51801,
        "Real Increment Verification Data",
        _ints([24]) + _reals([time] + [0.0] * 23),
    )
    if new_model is not None:
        out += _block(51900, "New Model", [])
        out += new_model
    if codes.columns:
        body = []
        for e in range(len(model.written)):
            for p in range(codes.nstres):
                body += _reals(codes.record(e, p))
        out += _block(52300, "Element Integration Point Values", body)
    body = _ints([len(nodal), 6])
    for name, a in nodal:
        ident = _NODAL_IDS.get(name, 0)
        body.append(_text_line(name, 48))
        body += _ints([ident, 0, 0, a.shape[1], 0, 0, -1, 0, 0, 0, 0, 0])
        body += _reals(a.reshape(-1))
    out += _block(52401, "Nodal Results", body)
    out.append("----")
    return out


def _time_of(mesh):
    fd = getattr(mesh, "field_data", {}) or {}
    if "meshio:time" in fd:
        return float(np.asarray(fd["meshio:time"]).reshape(-1)[0])
    return 0.0


class _Writer:
    """One post file: the model of the first mesh and an increment per mesh; a
    mesh that differs from the previous one repeats the model (a remeshing
    increment)."""

    def __init__(self, filename):
        if is_buffer(filename, "w"):
            raise WriteError("Marc .t19: cannot be written to a stream")
        self.path = filename
        self.header = None
        self.model = None
        self.nodal_signature = None
        self.increments = []

    def add(self, mesh, time):
        model = deck_model(mesh)
        if not model.written:
            raise WriteError("Marc .t19: the mesh has no cell with a Marc element type")
        codes = _Codes(mesh, model)
        nodal = _nodal_results(mesh, model.npts)
        nodal_signature = [(n, a.shape[1]) for n, a in nodal]
        blocks = _model(mesh, model, codes)
        new_model = None
        if self.header is None:
            self.header = blocks
            self.nodal_signature = nodal_signature
        elif nodal_signature != self.nodal_signature:
            raise WriteError(
                f"Marc .t19: series: step {len(self.increments)}'s point arrays "
                "differ from the first step's"
            )
        elif blocks != self.model:
            # A different mesh is a remeshing increment: it repeats the model.
            new_model = blocks
        self.model = blocks
        self.increments.append(
            _increment(len(self.increments), time, mesh, model, codes, nodal, new_model)
        )

    def close(self):
        if self.header is None:
            raise WriteError("Marc .t19: series: no step was written")
        out = _block(50100, "Analysis Title", ["          meshio++".ljust(70)])
        out += self.header
        for lines in self.increments:
            out += lines
        out.append("++++")
        with open(self.path, "w", encoding="latin-1", newline="\n") as f:
            f.write("\n".join(out) + "\n")


def write(filename, mesh):
    """Write one increment's post file: the mesh at ``field_data["meshio:time"]``."""
    writer = _Writer(filename)
    writer.add(mesh, _time_of(mesh))
    writer.close()


class SeriesWriter:
    """An increment per mesh, written when the writer closes; every mesh must
    have the first one's point arrays. A mesh whose points, cells or cell
    arrays differ from the previous one's is written as a remeshing increment."""

    def __init__(self, filename):
        self._writer = _Writer(filename)

    def write(self, time, mesh):
        self._writer.add(mesh, float(time))

    def close(self):
        if self._writer is not None:
            writer, self._writer = self._writer, None
            if writer.increments:
                writer.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
