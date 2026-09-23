"""
Reader for CalculiX result files (``.frd``), the ASCII format ``ccx`` writes and
``cgx`` reads.

A ``.frd`` is a stream of fixed-column records keyed by their first columns: ``1C``/
``1U`` header, ``2C`` node block, ``3C`` element block, then one ``100C`` block per
result per increment and ``9999`` to end. In the ASCII layout, data lines are ``-1``
(first line), ``-2`` (continuation) and ``-3`` (end of block); ``-4`` names a result
and ``-5`` names each of its components. Values are ``E12.5`` with *no separator*
(``7-1.18144E-06`` is two numbers), so every field is sliced by column. Node and
element ids are ``I5`` in the short format (flag 0) and ``I10`` in the long one
(flag 1).

The binary layout ccx writes for ``*NODE OUTPUT``/``*ELEMENT OUTPUT`` keeps the same
ASCII header lines (banner, ``2C``, ``3C``, ``1PSTEP``, ``100CL``, ``-4``, ``-5``) but
replaces the ``-1``/``-2``/``-3`` data lines with a raw little-endian blob: one
fixed-size ``[int32 id][reals]`` record per node/element/result entry, back to back,
with the record count read from the header's own count field (there is no terminator
to scan for). Reals are ``float32`` (flag 2) or ``float64`` (flag 3) per each block's
*own* flag -- ccx writes node coordinates as flag 3 and result values as flag 2 by
default, so the flag is re-read per header, never assumed constant for the file. This
means raw bytes, not decoded text, so ``_load``/``_decode`` keep the file as ``bytes``
until the ASCII-vs-binary branch. The C++ twin is ``formats/frd.cpp``.

Increments become the steps of the sequence engine: the mesh is read once and only the
requested step's result blocks are parsed. CalculiX expands shells and beams into
solids before it writes, so the mesh here is not the ``.inp`` mesh.
"""

import pathlib
import re
import struct

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._mesh import Mesh
from .._node_order import node_order
from .._tensor_invariants import _mises as _ti_mises
from .._tensor_invariants import _principal as _ti_principal

TIME_KEY = "meshio:time"

# FRD element type -> (cell type, node count). The node permutations (he20, pe15 and
# be3) live in meshioplusplus/_node_order.py under "frd". Types 7-10 are cgx's own
# shells: ccx expands its shells into solids and never writes them.
_TYPES = {
    1: ("hexahedron", 8),
    2: ("wedge", 6),
    3: ("tetra", 4),
    4: ("hexahedron20", 20),
    5: ("wedge15", 15),
    6: ("tetra10", 10),
    7: ("triangle", 3),
    8: ("triangle6", 6),
    9: ("quad", 4),
    10: ("quad8", 8),
    11: ("line", 2),
    12: ("line3", 3),
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
    """ASCII id column width: I5 (flag 0) or I10 (flag 1). Never called for a binary
    flag (2 or 3) -- binary detection happens before this is reached."""
    return 5 if flag == 0 else 10


def _detect_binary(raw):
    """Whether ``raw`` (the file's bytes) is the binary layout: the flag on its
    ``2C`` header line is 2 or 3 rather than 0 or 1. Scans only the always-ASCII
    banner that precedes ``2C``, which is safe even before binary detection."""
    start = 0
    while start < len(raw):
        stop = raw.find(b"\n", start)
        end = len(raw) if stop == -1 else stop
        line = raw[start:end].rstrip(b"\r").decode("ascii", errors="replace")
        if line.startswith("    2C"):
            return _flag(line, 1) >= 2
        if stop == -1:
            break
        start = stop + 1
    return False


class _Block:
    """One ``100C`` result block: where its lines (ASCII) or raw bytes (binary) are,
    and what they hold."""

    __slots__ = (
        "name",
        "ncomps",
        "data_ncomps",
        "irtype",
        "first",
        "last",
        "width",
        "binary",
        "byte_offset",
        "real_bytes",
        "num_entries",
    )

    def __init__(self, name, ncomps, irtype, width):
        self.name = name
        self.ncomps = ncomps
        self.data_ncomps = ncomps
        self.irtype = irtype
        self.width = width
        self.first = 0
        self.last = 0
        self.binary = False
        self.byte_offset = 0
        self.real_bytes = 8
        self.num_entries = 0


class _Frame:
    """One increment: a ``100C`` header shared by all the result blocks under it."""

    def __init__(self, key, value, analysis, step):
        self.key = key
        self.value = value
        self.analysis = analysis
        self.step = step
        self.blocks = []


class _Frd:
    def __init__(self, raw):
        self.node_ids = []
        self.coords = None
        self.elements = []  # (type, group, material, [file node ids])
        self.frames = []
        self.skipped_types = set()
        self.binary = _detect_binary(raw)
        if self.binary:
            self.raw = raw
            self.lines = []
            self._parse_binary()
        else:
            self.raw = None
            text = raw.decode("utf-8", errors="replace")
            self.lines = text.replace("\r", "").split("\n")
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

    # -- binary layout -------------------------------------------------------------
    #
    # Header lines stay plain ASCII text terminated by b"\n"; a raw little-endian
    # record blob immediately follows a "2C"/"3C" header or a "100C" frame's last
    # "-5" line, with no "-1"/"-2"/"-3" markers and no line boundaries of its own (a
    # record's bytes may well contain 0x0A). The record count comes from the header's
    # own count field. See the module docstring and doc/formats/frd.md.

    def _binary_line(self, pos):
        """One header line starting at byte ``pos``, decoded and stripped of a
        trailing ``\\r``; returns ``(line, next_pos)``, ``next_pos`` past its
        ``\\n`` (or end of file if there is none)."""
        raw = self.raw
        stop = raw.find(b"\n", pos)
        end = len(raw) if stop == -1 else stop
        line = raw[pos:end].decode("ascii", errors="replace").rstrip("\r")
        return line, (len(raw) if stop == -1 else stop + 1)

    def _binary_require(self, end, what):
        if end > len(self.raw):
            raise _err(f"binary {what} runs past the end of the file")

    def _parse_binary(self):
        pos = 0
        flag = 1
        seen_nodes = seen_elements = False
        frame_of = {}
        n = len(self.raw)
        while pos < n:
            line, next_pos = self._binary_line(pos)
            if line.startswith("    2C"):
                flag = _flag(line, flag)
                count = _int(line[6:36], "a 2C record")
                real_bytes = 8 if flag == 3 else 4
                rec = 4 + 3 * real_bytes
                self._binary_require(next_pos + count * rec, "node block")
                if seen_nodes:
                    warn("CalculiX FRD: a second node block was ignored")
                else:
                    self._read_nodes_binary(next_pos, count, real_bytes)
                seen_nodes = True
                pos = next_pos + count * rec
            elif line.startswith("    3C"):
                flag = _flag(line, flag)
                count = _int(line[6:36], "a 3C record")
                pos = self._read_elements_binary(
                    next_pos, count, keep=not seen_elements
                )
                if seen_elements:
                    warn("CalculiX FRD: a second element block was ignored")
                seen_elements = True
            elif line.startswith("  100C"):
                pos = self._read_frame_binary(next_pos, line, flag, frame_of)
            elif line.startswith(("9999", " 9999", "  9999")):
                break
            else:
                pos = next_pos
        if not seen_nodes:
            raise _err("no node block (2C record): not a CalculiX result file")

    def _read_nodes_binary(self, pos, count, real_bytes):
        rec = 4 + 3 * real_bytes
        blob = self.raw[pos : pos + count * rec]
        ids = np.frombuffer(
            blob, dtype=np.dtype(f"<i4, ({3},)f{real_bytes}"), count=count
        )
        self.node_ids = ids["f0"].tolist()
        self.coords = ids["f1"].astype(np.float64)

    def _read_elements_binary(self, pos, count, keep):
        raw = self.raw
        for _ in range(count):
            self._binary_require(pos + 16, "element header")
            etype, group, material = struct.unpack_from("<3i", raw, pos + 4)
            spec = _TYPES.get(etype)
            if spec is None:
                raise _err(f"unknown FRD element type {etype}")
            _, node_count = spec
            pos += 16
            self._binary_require(pos + node_count * 4, "element node list")
            if keep:
                nodes = list(struct.unpack_from(f"<{node_count}i", raw, pos))
                self.elements.append((etype, group, material, nodes))
            pos += node_count * 4
        return pos

    def _read_frame_binary(self, pos, header, flag, frame_of):
        key = header[6:12]
        value_text = header[12:24]
        frame_flag = flag
        if len(header) >= 75:
            frame_flag = _int(header[73:75], "a 100C record")
        real_bytes = 8 if frame_flag == 3 else 4
        numnod = _int(header[24:36], "a 100C record")
        analysis = _int(header[56:58], "a 100C record")
        step = _int(header[58:63], "a 100C record")
        ident = (key, value_text)
        frame = frame_of.get(ident)
        if frame is None:
            frame = _Frame(key, _real(value_text, "a 100C record"), analysis, step)
            frame_of[ident] = frame
            self.frames.append(frame)
        line, next_pos = self._binary_line(pos)
        if not line.startswith(" -4"):
            return next_pos
        ncomps = _int(line[13:18], "a -4 record")
        calculated = 0
        after = next_pos
        while True:
            l2, peek = self._binary_line(after)
            if not l2.startswith(" -5"):
                break
            if _int(l2[33:38], "a -5 record") == 1:
                calculated += 1
            after = peek
        block = _Block(line[5:13].strip(), ncomps, _int(line[18:23], "a -4 record"), 0)
        block.data_ncomps = ncomps - calculated
        block.binary = True
        block.real_bytes = real_bytes
        block.num_entries = numnod
        block.byte_offset = after
        rec = 4 + block.data_ncomps * real_bytes
        self._binary_require(after + numnod * rec, "result block")
        frame.blocks.append(block)
        return after + numnod * rec

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
            ctype, count = spec
            order = node_order("frd", ctype)
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
            if order is not None:
                row = [row[k] for k in order.to_meshio]
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
        if block.binary:
            return self._read_block_binary(block, index)
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

    def _read_block_binary(self, block, index):
        """Binary twin of ``read_block``: raw ``[int32 id][data_ncomps reals]``
        records, one per entry, no line markers."""
        nc = block.data_ncomps
        rec_size = 4 + nc * block.real_bytes
        start = block.byte_offset
        blob = self.raw[start : start + block.num_entries * rec_size]
        values = np.full((len(self.node_ids), nc), np.nan)
        real_dtype = f"<f{block.real_bytes}"
        for i in range(block.num_entries):
            base = i * rec_size
            node = struct.unpack_from("<i", blob, base)[0]
            row = index.get(node)
            if row is None:
                raise _err(f"result for {block.name} refers to undefined node {node}")
            if nc:
                values[row] = np.frombuffer(
                    blob, dtype=real_dtype, count=nc, offset=base + 4
                )
        return values[:, 0] if nc == 1 else values


# The eigensolver used to live here as private `_mises`/`_principal`; it now
# backs the public `tensor_invariants` operation as well, so both callers
# share `meshioplusplus._tensor_invariants` instead of drifting apart. The
# names stay as thin aliases so the call sites below read unchanged.
_mises = _ti_mises
_principal = _ti_principal


def _decode(source):
    """Returns the file's raw bytes -- decoding as text (only for the ASCII layout,
    and only after the binary-vs-ASCII branch) happens in ``_Frd.__init__``, since a
    genuinely binary file must never go through a UTF-8 decode."""
    if hasattr(source, "read"):
        data = source.read()
        if isinstance(data, str):
            data = data.encode("utf-8")
        return data
    try:
        return pathlib.Path(source).read_bytes()
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


# -- .dat: the ccx tabular print, a companion file with no mesh in it ----------------
#
# *NODE PRINT and *EL PRINT write whitespace-separated, free-format tables (unlike
# .frd, there is no glued-column layout to slice), one section per print request per
# increment, framed by "S T E P n" / "INCREMENT n" banners. Each section's own
# description line names its columns in parentheses, e.g.
# "stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz) for set E and time  ...": an
# element/integration-point section always leads with "elem, integ.pnt." (two id
# columns instead of node print's one), and the *value* column order is the file's
# own -- sxx syy szz sxy sxz syz, NOT .frd's xx yy zz xy yz zx. There is no geometry
# here, so this returns plain tables, not a Mesh, and it is Python-only: a native
# reader would need the same free-format tokenizer detail/fast_number.hpp does not
# provide, and .dat already belongs to Tecplot in the format registry, so this is
# never registered as one.

_DAT_STEP_RE = re.compile(r"\s*S\s+T\s+E\s+P\s+(\d+)")
_DAT_INCREMENT_RE = re.compile(r"\s*INCREMENT\s+(\d+)")
_DAT_SECTION_RE = re.compile(
    r"\s*(\S.*?)\s*\(([^)]*)\)\s+for set\s+(\S+)\s+and time\s+([-+0-9.EeDd]+)\s*$"
)


def read_dat(filename):
    """Parses a ccx ``.dat`` tabular print (``*NODE PRINT``/``*EL PRINT``) into a list
    of tables -- there is no mesh in this file, so this returns plain data, not a
    :class:`~meshioplusplus.Mesh`.

    Each table is a dict:

    - ``step``, ``increment``: the ``S T E P`` / ``INCREMENT`` banner values.
    - ``time``: the section's own time value.
    - ``kind``: ``"node"`` or ``"element"``.
    - ``quantity``: the description before the parenthesised column list
      (``"displacements"``, ``"stresses"``, ...).
    - ``set``: the ``NSET``/``ELSET`` name the print request named.
    - ``ids``: point/node ids (``kind="node"``) or element ids (``kind="element"``).
    - ``int_points``: integration-point indices, only present when ``kind="element"``.
    - ``components``: the column names, in the file's own order -- for stress/strain
      this is ``sxx syy szz sxy sxz syz``, not ``.frd``'s ``xx yy zz xy yz zx``.
    - ``values``: ``(len(ids), len(components))`` float64 array.

    :raises ReadError: if the file cannot be read or is not a recognisable ``.dat``.
    """
    raw = _decode(filename)
    text = raw.decode("utf-8", errors="replace")
    lines = text.replace("\r", "").split("\n")
    tables = []
    step = increment = None
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        m = _DAT_STEP_RE.match(line)
        if m:
            step = _int(m.group(1), "a S T E P banner")
            i += 1
            continue
        m = _DAT_INCREMENT_RE.match(line)
        if m:
            increment = _int(m.group(1), "an INCREMENT banner")
            i += 1
            continue
        m = _DAT_SECTION_RE.match(line)
        if not m:
            i += 1
            continue
        quantity, cols_text, set_name, time_text = m.groups()
        cols = [c.strip() for c in cols_text.split(",")]
        kind = "element" if cols[:2] == ["elem", "integ.pnt."] else "node"
        components = cols[2:] if kind == "element" else cols
        time_value = _real(time_text, "a .dat section header")
        i += 1
        while i < n and not lines[i].strip():
            i += 1
        ids = []
        int_points = [] if kind == "element" else None
        rows = []
        while i < n and lines[i].strip():
            parts = lines[i].split()
            skip = 2 if kind == "element" else 1
            ids.append(_int(parts[0], "a .dat data row"))
            if kind == "element":
                int_points.append(_int(parts[1], "a .dat data row"))
            rows.append([_real(v, "a .dat data row") for v in parts[skip:]])
            i += 1
        table = {
            "step": step,
            "increment": increment,
            "time": time_value,
            "kind": kind,
            "quantity": quantity,
            "set": set_name,
            "ids": np.asarray(ids, dtype=np.int64),
            "components": components,
            "values": np.asarray(rows, dtype=np.float64).reshape(
                len(ids), len(components)
            ),
        }
        if kind == "element":
            table["int_points"] = np.asarray(int_points, dtype=np.int64)
        tables.append(table)
    if not tables:
        raise _err(f"no recognisable *NODE PRINT/*EL PRINT section in {filename}")
    return tables
