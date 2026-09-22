"""
Reader for CalculiX result files (``.frd``), the ASCII format ``ccx`` writes and
``cgx`` reads.

A ``.frd`` is a stream of fixed-column records keyed by their first columns: ``1C``/
``1U`` header, ``2C`` node block, ``3C`` element block, then one ``100C`` block per
result per increment and ``9999`` to end. Data lines are ``-1`` (first line), ``-2``
(continuation) and ``-3`` (end of block); ``-4`` names a result and ``-5`` names each
of its components. Values are ``E12.5`` with *no separator* (``7-1.18144E-06`` is two
numbers), so every field is sliced by column. Node and element ids are ``I5`` in the
short format (flag 0) and ``I10`` in the long one (flag 1); flag 2 is binary and is
not read. The C++ twin is ``formats/frd.cpp``.

Increments become the steps of the sequence engine: the mesh is read once and only the
requested step's result blocks are parsed. CalculiX expands shells and beams into
solids before it writes, so the mesh here is not the ``.inp`` mesh.
"""

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._mesh import Mesh
from .._tensor_invariants import _mises as _ti_mises
from .._tensor_invariants import _principal as _ti_principal

TIME_KEY = "meshio:time"

# FRD element type -> (cell type, node count, permutation). ``connectivity[k] =
# frd_nodes[permutation[k]]``; None is the identity. The he20 and pe15 mid-edge groups
# and the be3 mid-node sit elsewhere than in Abaqus order (confirmed against ccx 2.23
# output for the same ``.inp``). Types 7-10 are cgx's own shells: ccx expands its
# shells into solids and never writes them.
_TYPES = {
    1: ("hexahedron", 8, None),
    2: ("wedge", 6, None),
    3: ("tetra", 4, None),
    4: (
        "hexahedron20",
        20,
        list(range(12)) + list(range(16, 20)) + list(range(12, 16)),
    ),
    5: ("wedge15", 15, list(range(9)) + [12, 13, 14, 9, 10, 11]),
    6: ("tetra10", 10, None),
    7: ("triangle", 3, None),
    8: ("triangle6", 6, None),
    9: ("quad", 4, None),
    10: ("quad8", 8, None),
    11: ("line", 2, None),
    12: ("line3", 3, [0, 2, 1]),
}

# Result names whose 6 components are a symmetric tensor xx yy zz xy yz zx.
_TENSOR_NAMES = ("STRESS", "TOSTRAIN", "MESTRAIN", "ZZSTR")

_VALUES_PER_LINE = 6
_VALUE_WIDTH = 12


def _err(message):
    return ReadError(f"CalculiX FRD: {message}")


def _int(text, where):
    text = text.strip()
    if not text:
        return 0
    try:
        return int(text)
    except ValueError:
        raise _err(f"invalid integer field {text!r} in {where}") from None


def _real(text, where):
    s = text.strip().replace("D", "E").replace("d", "e")
    try:
        return float(s)
    except ValueError:
        raise _err(f"invalid real field {text!r} in {where}") from None


def _flag(line, default):
    """Trailing format flag of a ``2C``/``3C`` header line (0 short, 1 long)."""
    tokens = line.split()
    if len(tokens) >= 3:
        try:
            return int(tokens[-1])
        except ValueError:
            pass
    return default


def _width(flag):
    if flag == 2:
        raise _err("binary .frd files are not supported")
    return 5 if flag == 0 else 10


class _Block:
    """One ``100C`` result block: where its lines are and what they hold."""

    __slots__ = ("name", "ncomps", "data_ncomps", "irtype", "first", "last", "width")

    def __init__(self, name, ncomps, irtype, width):
        self.name = name
        self.ncomps = ncomps
        self.data_ncomps = ncomps
        self.irtype = irtype
        self.width = width
        self.first = 0
        self.last = 0


class _Frame:
    """One increment: a ``100C`` header shared by all the result blocks under it."""

    def __init__(self, key, value, analysis, step):
        self.key = key
        self.value = value
        self.analysis = analysis
        self.step = step
        self.blocks = []


class _Frd:
    def __init__(self, text):
        self.lines = text.replace("\r", "").split("\n")
        self.node_ids = []
        self.coords = None
        self.elements = []  # (type, group, material, [file node ids])
        self.frames = []
        self.skipped_types = set()
        self._parse()

    def _parse(self):
        lines = self.lines
        n = len(lines)
        pos = 0
        flag = 1
        seen_nodes = seen_elements = False
        frame_of = {}
        while pos < n:
            line = lines[pos].rstrip("\r")
            pos += 1
            if line.startswith("    2C"):
                flag = _flag(line, flag)
                if seen_nodes:
                    warn("CalculiX FRD: a second node block was ignored")
                    pos = self._skip_block(pos)
                else:
                    pos = self._read_nodes(pos, _width(flag))
                    seen_nodes = True
            elif line.startswith("    3C"):
                flag = _flag(line, flag)
                if seen_elements:
                    warn("CalculiX FRD: a second element block was ignored")
                    pos = self._skip_block(pos)
                else:
                    pos = self._read_elements(pos, _width(flag))
                    seen_elements = True
            elif line.startswith("  100C"):
                pos = self._read_frame(pos, line, flag, frame_of)
            elif line.startswith("  9999"):
                break
        if not seen_nodes:
            raise _err("no node block (2C record): not a CalculiX result file")

    def _skip_block(self, pos):
        while pos < len(self.lines) and not self.lines[pos].startswith(" -3"):
            pos += 1
        return pos + 1

    def _read_nodes(self, pos, w):
        lines = self.lines
        ids = []
        coords = []
        end = 3 + w
        while pos < len(lines):
            line = lines[pos]
            pos += 1
            if line.startswith(" -3"):
                break
            if not line.startswith(" -1"):
                continue
            ids.append(_int(line[3:end], "a node line"))
            coords.append(
                [
                    _real(
                        line[end + k * _VALUE_WIDTH : end + (k + 1) * _VALUE_WIDTH],
                        "a node line",
                    )
                    for k in range(3)
                ]
            )
        self.node_ids = ids
        self.coords = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
        return pos

    def _read_elements(self, pos, w):
        lines = self.lines
        end = 3 + w
        current = None
        while pos < len(lines):
            line = lines[pos]
            pos += 1
            if line.startswith(" -3"):
                break
            if line.startswith(" -1"):
                etype = _int(line[end : end + 5], "an element line")
                group = _int(line[end + 5 : end + 10], "an element line")
                material = _int(line[end + 10 : end + 15], "an element line")
                current = (etype, group, material, [])
                self.elements.append(current)
            elif line.startswith(" -2") and current is not None:
                body = line[3:].rstrip()
                if len(body) % w:
                    raise _err(f"element line {line!r} is not a whole number of ids")
                current[3].extend(
                    _int(body[k : k + w], "an element line")
                    for k in range(0, len(body), w)
                )
        return pos

    def _read_frame(self, pos, header, flag, frame_of):
        lines = self.lines
        key = header[6:12]
        value_text = header[12:24]
        frame_flag = flag
        if len(header) >= 75:
            frame_flag = _int(header[73:75], "a 100C record")
        w = _width(frame_flag)
        analysis = _int(header[56:58], "a 100C record")
        step = _int(header[58:63], "a 100C record")
        ident = (key, value_text)
        frame = frame_of.get(ident)
        if frame is None:
            frame = _Frame(key, _real(value_text, "a 100C record"), analysis, step)
            frame_of[ident] = frame
            self.frames.append(frame)
        block = None
        while pos < len(lines):
            line = lines[pos]
            pos += 1
            if line.startswith(" -4"):
                block = _Block(
                    line[5:13].strip(),
                    _int(line[13:18], "a -4 record"),
                    _int(line[18:23], "a -4 record"),
                    w,
                )
                calculated = 0
                while pos < len(lines) and lines[pos].startswith(" -5"):
                    if _int(lines[pos][33:38], "a -5 record") == 1:
                        calculated += 1
                    pos += 1
                block.data_ncomps = block.ncomps - calculated
                block.first = pos
                while pos < len(lines) and not lines[pos].startswith(" -3"):
                    pos += 1
                block.last = pos
                pos += 1
                frame.blocks.append(block)
                break
            if line.startswith(("  100", "    1P", "  9999")):
                pos -= 1
                break
        return pos

    # -- mesh --------------------------------------------------------------------
    def build_cells(self):
        index = {nid: i for i, nid in enumerate(self.node_ids)}
        cells = []
        group = []
        material = []
        for etype, grp, mat, nodes in self.elements:
            spec = _TYPES.get(etype)
            if spec is None:
                self.skipped_types.add(etype)
                continue
            ctype, count, perm = spec
            if len(nodes) != count:
                raise _err(
                    f"element type {etype} ({ctype}) needs {count} nodes, found {len(nodes)}"
                )
            try:
                row = [index[nid] for nid in nodes]
            except KeyError as exc:
                raise _err(
                    f"an element references undefined node {exc.args[0]}"
                ) from None
            if perm is not None:
                row = [row[k] for k in perm]
            if not cells or cells[-1][0] != ctype:
                cells.append((ctype, []))
                group.append([])
                material.append([])
            cells[-1][1].append(row)
            group[-1].append(grp)
            material[-1].append(mat)
        if self.skipped_types:
            warn(
                "CalculiX FRD: skipped elements of unsupported type(s) "
                + ", ".join(str(t) for t in sorted(self.skipped_types))
            )
        return (
            [(t, np.asarray(rows, dtype=np.int64)) for t, rows in cells],
            [np.asarray(g, dtype=np.int64) for g in group],
            [np.asarray(m, dtype=np.int64) for m in material],
            index,
        )

    # -- results -----------------------------------------------------------------
    def read_block(self, block, index):
        lines = self.lines
        w = block.width
        nc = block.data_ncomps
        values = np.full((len(self.node_ids), nc), np.nan)
        pos = block.first
        end = 3 + w
        while pos < block.last:
            line = lines[pos]
            pos += 1
            if not line.startswith(" -1"):
                continue
            node = _int(line[3:end], "a result line")
            row = index.get(node)
            if row is None:
                raise _err(f"result for {block.name} refers to undefined node {node}")
            got = 0
            while True:
                take = min(_VALUES_PER_LINE, nc - got)
                for k in range(take):
                    lo = end + k * _VALUE_WIDTH
                    field = line[lo : lo + _VALUE_WIDTH]
                    if len(field.strip()) == 0:
                        raise _err(f"result line for {block.name} is missing values")
                    values[row, got + k] = _real(field, "a result line")
                got += take
                if got >= nc:
                    break
                if pos >= block.last or not lines[pos].startswith(" -2"):
                    raise _err(f"result line for {block.name} is missing values")
                line = lines[pos]
                pos += 1
        return values[:, 0] if nc == 1 else values


# The eigensolver used to live here as private `_mises`/`_principal`; it now
# backs the public `tensor_invariants` operation as well, so both callers
# share `meshioplusplus._tensor_invariants` instead of drifting apart. The
# names stay as thin aliases so the call sites below read unchanged.
_mises = _ti_mises
_principal = _ti_principal


def _decode(source):
    if hasattr(source, "read"):
        data = source.read()
        if isinstance(data, bytes):
            data = data.decode("utf-8", errors="replace")
        return data
    try:
        # binary, so no newline translation: a stray "\r\r\n" must not become a blank line
        with open(source, "rb") as f:
            return f.read().decode("utf-8", errors="replace")
    except OSError as exc:
        raise _err(f"could not read {source}: {exc}") from None


def _resolve_step(time_step, count):
    if count == 0:
        if time_step in (0, -1):
            return None
        raise ReadError(
            f"meshio++: time step {time_step} requested, but this file carries no time steps"
        )
    step = time_step + count if time_step < 0 else time_step
    if not 0 <= step < count:
        raise ReadError(
            f"meshio++: time step {time_step} is out of range: this file has {count} "
            + ("step" if count == 1 else "steps")
        )
    return step


def _load(filename):
    return _Frd(_decode(filename))


def read(filename, points_only=False, arrays=None, time_step=0, derived=False):
    frd = _load(filename)
    cells, group, material, index = frd.build_cells()
    mesh = Mesh(frd.coords.copy(), cells)
    if group:
        mesh.cell_data["frd:group"] = group
        mesh.cell_data["frd:material"] = material
    step = _resolve_step(time_step, len(frd.frames))
    if step is None:
        return mesh
    mesh.time_values = [f.value for f in frd.frames]
    frame = frd.frames[step]
    mesh.field_data[TIME_KEY] = np.array([frame.value])
    mesh.field_data["frd:step"] = np.array([frame.step], dtype=np.int64)
    mesh.field_data["frd:analysis"] = np.array([frame.analysis], dtype=np.int64)
    if points_only:
        return mesh
    wanted = None if arrays is None else set(arrays)

    def keep(name):
        return wanted is None or name in wanted

    used = set()

    def unique(name):
        candidate, k = name, 2
        while candidate in used:
            candidate = f"{name}_{k}"
            k += 1
        used.add(candidate)
        return candidate

    for block in frame.blocks:
        name = unique(block.name)
        wants_derived = (
            derived
            and block.name in _TENSOR_NAMES
            and block.data_ncomps == 6
            and (keep(f"{name}_mises") or keep(f"{name}_principal"))
        )
        if not (keep(name) or wants_derived):
            continue
        values = frd.read_block(block, index)
        if keep(name):
            mesh.point_data[name] = values
        if wants_derived:
            if keep(f"{name}_mises"):
                mesh.point_data[f"{name}_mises"] = _mises(values)
            if keep(f"{name}_principal"):
                mesh.point_data[f"{name}_principal"] = _principal(values)
    return mesh


def time_values(filename):
    """The value of every increment, in file order (the sequence engine's steps)."""
    return [frame.value for frame in _load(filename).frames]
