"""
I/O for COMSOL native mesh files: text ``.mphtxt`` and binary ``.mphbin``.

Both serialise the same sequence: a version pair (``0 1``), tag and type tables,
then objects (``0 0 1`` and a class name). ``Mesh`` objects (version 4 and
older, whose element records also carry parameter values and up/down pairs) and
``Selection`` objects are read; the first object of any other class stops the
read with a warning. Several Mesh objects are merged, each one's cells a region
named by its tag.

Every element's geometric entity index becomes ``cell_data["mphtxt:geom"]``
(domains count from 1, boundaries, edges and points from 0); a Selection becomes
a cell region: the elements of its dimension whose entity it lists. Node order
goes through the node-ordering registry (format key ``mphtxt``).

The C++ core (``formats/mphtxt.cpp``) is the twin of this module; both read the
same Mesh and write the same bytes.
"""

from __future__ import annotations

import struct

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import is_buffer, open_file
from .._mesh import CellBlock, Mesh, topological_dimension
from .._node_order import node_order
from .._regions import Region

_COMSOL_TO_MESHIO = {
    "vtx": "vertex",
    "edg": "line",
    "tri": "triangle",
    "quad": "quad",
    "tet": "tetra",
    "prism": "wedge",
    "pyr": "pyramid",
    "hex": "hexahedron",
    "edg2": "line3",
    "tri2": "triangle6",
    "quad2": "quad9",
    "tet2": "tetra10",
    "prism2": "wedge18",
    "hex2": "hexahedron27",
    "pyr2": "pyramid14",
}
_MESHIO_TO_COMSOL = {v: k for k, v in _COMSOL_TO_MESHIO.items()}
_NUM_NODES = {
    "vertex": 1,
    "line": 2,
    "triangle": 3,
    "quad": 4,
    "tetra": 4,
    "wedge": 6,
    "pyramid": 5,
    "hexahedron": 8,
    "line3": 3,
    "triangle6": 6,
    "quad9": 9,
    "tetra10": 10,
    "wedge18": 18,
    "hexahedron27": 27,
    "pyramid14": 14,
}


# ---------------------------------------------------------------------------
# Sources: the text and the binary serialisation hold the same values.
# ---------------------------------------------------------------------------


def _is_int_token(t):
    s = t[1:] if t[:1] in ("-", "+") else t
    return s.isdigit() and s.isascii()


class _Text:
    def __init__(self, text, fmt):
        self.text = text
        self.pos = 0
        self.fmt = fmt

    def _skip_blank(self):
        text, n = self.text, len(self.text)
        while self.pos < n:
            c = text[self.pos]
            if c == "#":
                eol = text.find("\n", self.pos)
                self.pos = n if eol < 0 else eol
            elif c in " \t\n\r":
                self.pos += 1
            else:
                break

    def _token(self):
        self._skip_blank()
        if self.pos >= len(self.text):
            raise ReadError(f"{self.fmt}: unexpected end of file")
        b = self.pos
        text, n = self.text, len(self.text)
        while self.pos < n and text[self.pos] not in " \t\n\r#":
            self.pos += 1
        return text[b : self.pos]

    def int(self):
        t = self._token()
        if not _is_int_token(t):
            raise ReadError(f"{self.fmt}: expected an integer, found '{t}'")
        return int(t)

    def ints(self, n):
        return [self.int() for _ in range(n)]

    def real(self):
        t = self._token()
        try:
            return float(t)
        except ValueError:
            raise ReadError(f"{self.fmt}: expected a number, found '{t}'")

    def reals(self, n):
        return [self.real() for _ in range(n)]

    def string(self):
        # a length, one blank, then exactly that many characters (which may
        # include blanks and '#')
        n = self.int()
        if n < 0:
            raise ReadError(f"{self.fmt}: negative string length")
        if n == 0:
            return ""
        if self.pos < len(self.text) and self.text[self.pos] in " \t":
            self.pos += 1
        if self.pos + n > len(self.text):
            raise ReadError(f"{self.fmt}: unexpected end of file in a string")
        out = self.text[self.pos : self.pos + n]
        self.pos += n
        return out

    def at_end(self):
        self._skip_blank()
        return self.pos >= len(self.text)

    def skip_value(self):
        self._token()

    def next_is_integer(self):
        keep = self.pos
        if self.at_end():
            return False
        t = self._token()
        self.pos = keep
        return _is_int_token(t)

    def next_is_name(self):
        keep = self.pos
        ok = self.next_is_integer()
        if ok:
            self._token()
            self._skip_blank()
            ok = self.pos < len(self.text) and (
                "a" <= self.text[self.pos] <= "z" or "A" <= self.text[self.pos] <= "Z"
            )
        self.pos = keep
        return ok


class _Binary:
    def __init__(self, data, fmt):
        self.data = data
        self.pos = 0
        self.fmt = fmt

    def _need(self, n):
        if self.pos + n > len(self.data):
            raise ReadError(f"{self.fmt}: unexpected end of file")

    def int(self):
        self._need(4)
        (v,) = struct.unpack_from("<i", self.data, self.pos)
        self.pos += 4
        return v

    def ints(self, n):
        self._need(4 * n)
        out = np.frombuffer(self.data, dtype="<i4", count=n, offset=self.pos)
        self.pos += 4 * n
        return out.astype(np.int64)

    def real(self):
        self._need(8)
        (v,) = struct.unpack_from("<d", self.data, self.pos)
        self.pos += 8
        return v

    def reals(self, n):
        self._need(8 * n)
        out = np.frombuffer(self.data, dtype="<f8", count=n, offset=self.pos)
        self.pos += 8 * n
        return out.astype(np.float64)

    def string(self):
        n = self.int()
        if n < 0 or n > (len(self.data) - self.pos) // 4:
            raise ReadError(f"{self.fmt}: invalid string length {n}")
        cps = self.ints(n)
        return "".join(chr(int(c)) for c in cps)

    def at_end(self):
        return self.pos >= len(self.data)

    def skip_value(self):
        self._need(8)
        self.pos += 8

    def next_is_integer(self):
        return self.pos + 4 <= len(self.data)

    def next_is_name(self):
        if self.pos + 8 > len(self.data):
            return False
        n, c = struct.unpack_from("<ii", self.data, self.pos)
        return 0 < n < 64 and (ord("a") <= c <= ord("z") or ord("A") <= c <= ord("Z"))


def _object_or_end(src):
    """Whether the next values end the file or open an object: 0 0 1 and a
    class name."""
    if src.at_end():
        return True
    for expected in (0, 0, 1):
        if not src.next_is_integer() or src.int() != expected:
            return False
    return src.next_is_name()


def _tail_fits(src, ne, last):
    """Whether the rest of an element record of a version < 4 Mesh fits once
    the parameter values are skipped: ``ne`` entity indices, then up/down
    pairs, then the next type name, the next object or the end."""
    keep = src.pos
    ok = False
    try:
        ngeom = src.int() if src.next_is_integer() else -1
        if ngeom in (ne, 0):
            ok = True
            k = 0
            while k < ngeom and ok:
                ok = src.next_is_integer() and src.int() >= 0
                k += 1
            if ok and src.next_is_integer():
                nud = src.int()
                ok = nud >= 0
                k = 0
                while k < 2 * nud and ok:
                    ok = src.next_is_integer()
                    if ok:
                        src.int()
                    k += 1
                if ok:
                    ok = _object_or_end(src) if last else src.next_is_name()
            else:
                ok = False
    except ReadError:
        ok = False
    src.pos = keep
    return ok


def _read_mesh(src, fmt):
    version = src.int()
    sdim = src.int()
    npts = src.int()
    lowest = src.int()
    if not 1 <= sdim <= 3 or npts < 0:
        raise ReadError(f"{fmt}: invalid Mesh header (sdim {sdim}, {npts} vertices)")
    points = np.asarray(src.reals(npts * sdim), dtype=np.float64).reshape(npts, sdim)
    blocks = []
    ntypes = src.int()
    for t in range(ntypes):
        ctype = src.string()
        if ctype not in _COMSOL_TO_MESHIO:
            raise ReadError(f"{fmt}: unknown element type '{ctype}'")
        cell_type = _COMSOL_TO_MESHIO[ctype]
        nep = src.int()
        ne = src.int()
        if nep != _NUM_NODES[cell_type] or ne < 0:
            raise ReadError(f"{fmt}: '{ctype}' elements with {nep} nodes")
        raw = np.asarray(src.ints(ne * nep), dtype=np.int64).reshape(ne, nep) - lowest
        bad = raw[(raw < 0) | (raw >= npts)]
        if len(bad):
            raise ReadError(
                f"{fmt}: '{ctype}' element names vertex {int(bad[0]) + lowest}"
            )
        order = node_order("mphtxt", cell_type)
        conn = raw[:, list(order.to_meshio)] if order is not None else raw
        if version < 4:
            # parameter records of npp values each, one value per parameter
            # dimension: find how many by where the rest of the record fits
            npp = src.int()
            npar = src.int()
            start = src.pos
            found = False
            for k in (1, 2, 3):
                src.pos = start
                try:
                    for _ in range(npp * npar * k):
                        src.skip_value()
                    found = _tail_fits(src, ne, t + 1 == ntypes)
                except ReadError:
                    found = False
                if found:
                    break
            if not found:
                raise ReadError(
                    f"{fmt}: cannot delimit the parameters of the '{ctype}' elements"
                )
        # no entity indices at all is allowed: every element then gets COMSOL's
        # default, domain 1 or entity 0 below the space dimension
        ngeom = src.int()
        if ngeom not in (ne, 0):
            raise ReadError(
                f"{fmt}: {ngeom} entity indices for {ne} '{ctype}' elements"
            )
        if ngeom:
            geom = np.asarray(src.ints(ngeom), dtype=np.int64)
        else:
            dim = topological_dimension[cell_type]
            geom = np.full(ne, 1 if dim == sdim else 0, dtype=np.int64)
        if version < 4:
            nud = src.int()
            src.ints(2 * nud)
        blocks.append((cell_type, conn, geom))
    return points, blocks


def _read(src, fmt):
    major, minor = src.int(), src.int()
    if (major, minor) != (0, 1):
        raise ReadError(f"{fmt}: unsupported file version {major}.{minor}")
    tags = [src.string() for _ in range(max(0, src.int()))]
    ntypes = src.int()
    for _ in range(ntypes):
        src.string()

    meshes = []  # (tag, points, blocks)
    selections = []  # (label, mesh tag, dim, entities)
    for obj in range(ntypes):
        if src.at_end():
            break
        src.int()
        src.int()
        src.int()
        cls = src.string()
        tag = tags[obj] if obj < len(tags) else ""
        if cls == "Mesh":
            points, blocks = _read_mesh(src, fmt)
            meshes.append((tag, points, blocks))
        elif cls == "Selection":
            src.int()  # version
            label = src.string()
            mesh_tag = src.string()
            dim = src.int()
            entities = [int(e) for e in src.ints(max(0, src.int()))]
            selections.append((label, mesh_tag, dim, entities))
        else:
            warn(
                f"{fmt}: object {obj} is a '{cls}', which meshio++ does not read; "
                "the objects after it are skipped too"
            )
            break
    if not meshes:
        raise ReadError(f"{fmt}: the file holds no Mesh object")

    sdim = meshes[0][1].shape[1]
    if any(m[1].shape[1] != sdim for m in meshes):
        raise ReadError(f"{fmt}: Mesh objects of different dimensions")
    points = np.concatenate([m[1] for m in meshes]).reshape(-1, sdim)
    cells, geom, dims, first_block = [], [], [], []
    offset = 0
    for _, pts, blocks in meshes:
        first_block.append(len(cells))
        for cell_type, conn, g in blocks:
            cells.append(CellBlock(cell_type, conn + offset))
            geom.append(g)
            dims.append(topological_dimension[cell_type])
        offset += len(pts)
    mesh = Mesh(points, cells)
    if cells:
        mesh.cell_data["mphtxt:geom"] = geom

    bases = np.concatenate([[0], np.cumsum([len(c) for c in cells])]).astype(np.int64)
    regions = []
    taken = set()

    def unique_name(name):
        out = name
        k = 2
        while out in taken:
            out = f"{name} ({k})"
            k += 1
        if out != name:
            warn(f"{fmt}: a second region named '{name}' is read as '{out}'")
        taken.add(out)
        return out

    # with several Mesh objects, each one's cells are a region named by its tag
    if len(meshes) > 1:
        for o, (tag, _, blocks) in enumerate(meshes):
            b0, b1 = first_block[o], first_block[o] + len(blocks)
            ids = np.arange(bases[b0], bases[b1], dtype=np.int64)
            regions.append(Region(unique_name(tag), "cell", ids, -1, -1))
    # a Selection: the elements of its dimension whose entity it lists
    for label, mesh_tag, dim, entities in selections:
        o = next((k for k, m in enumerate(meshes) if m[0] == mesh_tag), None)
        if o is None:
            warn(
                f"{fmt}: selection '{label}' refers to '{mesh_tag}', which is not "
                "a Mesh here; skipped"
            )
            continue
        ids = []
        for b in range(first_block[o], first_block[o] + len(meshes[o][2])):
            if dims[b] != dim:
                continue
            hit = np.flatnonzero(np.isin(geom[b], entities))
            ids.append(bases[b] + hit)
        ids = np.concatenate(ids).astype(np.int64) if ids else np.zeros(0, np.int64)
        regions.append(Region(unique_name(label), "cell", ids, dim, -1))
    if regions:
        mesh.regions = regions
    return mesh


def read(filename):
    if is_buffer(filename, "r"):
        data = filename.read()
    else:
        with open_file(filename, "rb") as f:
            data = f.read()
    if isinstance(data, bytes):
        data = data.decode("utf-8")
    return _read(_Text(data, "mphtxt"), "mphtxt")


def read_binary(filename):
    if is_buffer(filename, "rb"):
        data = filename.read()
    else:
        with open_file(filename, "rb") as f:
            data = f.read()
    return _read(_Binary(data, "mphbin"), "mphbin")


# ---------------------------------------------------------------------------
# Writing: one serialiser over two sinks.
# ---------------------------------------------------------------------------


class _TextSink:
    def __init__(self):
        self.out = []

    def _end(self, comment):
        self.out.append(f" # {comment}\n" if comment else "\n")

    def int(self, v, comment=None):
        self.out.append(str(int(v)))
        self._end(comment)

    def ints(self, row):
        self.out.append(" ".join(str(int(v)) for v in row) + "\n")

    def reals(self, row):
        self.out.append(" ".join(f"{float(v):.17g}" for v in row) + "\n")

    def string(self, s, comment=None):
        self.out.append(f"{len(s)} {s}")
        self._end(comment)

    def comment(self, text):
        self.out.append(f"# {text}\n")

    def blank(self):
        self.out.append("\n")


class _BinarySink:
    def __init__(self, fmt):
        self.out = []
        self.fmt = fmt

    def _put(self, values):
        a = np.asarray(values, dtype=np.int64)
        if a.size and (a.min() < -(2**31) or a.max() > 2**31 - 1):
            bad = a[(a < -(2**31)) | (a > 2**31 - 1)][0]
            raise WriteError(f"{self.fmt}: {bad} does not fit a 32-bit integer")
        self.out.append(a.astype("<i4").tobytes())

    def int(self, v, comment=None):
        self._put([v])

    def ints(self, row):
        self._put(row)

    def reals(self, row):
        self.out.append(np.asarray(row, dtype="<f8").tobytes())

    def string(self, s, comment=None):
        self._put([len(s)] + [ord(c) for c in s])

    def comment(self, text):
        pass

    def blank(self):
        pass


def _kind_key(r):
    return ({"point": 0, "cell": 1, "side": 2}[r.kind], r.name, r.dim, r.tag)


def _serialise(sink, mesh, fmt):
    sdim = mesh.points.shape[1]
    ctypes, dims = [], []
    for block in mesh.cells:
        if block.type not in _MESHIO_TO_COMSOL:
            raise WriteError(f"{fmt}: unsupported cell type {block.type}")
        ctypes.append(_MESHIO_TO_COMSOL[block.type])
        dims.append(topological_dimension[block.type])
    sizes = [len(b.data) for b in mesh.cells]
    bases = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
    ncells = int(bases[-1])
    cell_dim = np.repeat(np.array(dims, dtype=np.int64), sizes)

    # cell regions, validated, in the C++ mesh's (kind, name, dim, tag) order
    regions = []
    for reg in sorted(mesh.regions, key=_kind_key):
        if reg.kind != "cell":
            warn(
                f"{fmt}: {reg.kind} region '{reg.name}' dropped; COMSOL selections "
                "are written for cell regions only"
            )
            continue
        e = np.asarray(reg.entries, dtype=np.int64)
        bad = e[(e < 0) | (e >= ncells)]
        if len(bad):
            raise WriteError(
                f"{fmt}: cell region '{reg.name}' names cell {bad[0]} of {ncells}"
            )
        regions.append((reg, e))

    # geometric entity of every cell: mphtxt:geom when the mesh has it; else per
    # dimension, the pairwise-disjoint single-dimension cell regions in order,
    # then one entity for the cells in none. Domains (dimension sdim) count
    # from 1, lower dimensions from 0.
    if "mphtxt:geom" in mesh.cell_data:
        entity = np.concatenate(
            [np.asarray(g, dtype=np.int64) for g in mesh.cell_data["mphtxt:geom"]]
            or [np.zeros(0, np.int64)]
        )
    else:
        entity = np.zeros(ncells, dtype=np.int64)
        assigned = np.zeros(ncells, dtype=bool)
        nxt = {}

        def base(d):
            return 1 if d == sdim else 0

        for reg, e in regions:
            if len(e) == 0:
                continue
            d = int(cell_dim[e[0]])
            if (cell_dim[e] != d).any() or assigned[e].any():
                continue
            k = nxt.setdefault(d, base(d))
            entity[e] = k
            assigned[e] = True
            nxt[d] = k + 1
        for c in np.flatnonzero(~assigned):
            d = int(cell_dim[c])
            entity[c] = nxt.get(d, base(d))

    # a cell region becomes a Selection when it is exactly a union of whole
    # entities of one dimension
    selections = []
    for reg, e in regions:
        if len(e) == 0:
            selections.append((reg.name, reg.dim if reg.dim >= 0 else sdim, []))
            continue
        d = int(cell_dim[e[0]])
        ents = sorted({int(x) for x in entity[e]})
        covered = int(np.count_nonzero((cell_dim == d) & np.isin(entity, ents)))
        if (cell_dim[e] != d).any() or covered != len(e):
            warn(
                f"{fmt}: cell region '{reg.name}' is not a union of whole geometric "
                "entities of one dimension; dropped"
            )
            continue
        selections.append((reg.name, d, ents))

    nobjects = 1 + len(selections)
    sink.int(0)
    sink.int(1, "version")
    sink.int(nobjects, "number of tags")
    sink.string("mesh1")
    for k in range(1, nobjects):
        sink.string(f"mesh1_sel{k}")
    sink.int(nobjects, "number of types")
    for _ in range(nobjects):
        sink.string("obj")
    sink.blank()
    sink.comment("--------- Object 0 ----------")
    sink.int(0)
    sink.int(0)
    sink.int(1)
    sink.string("Mesh", "class")
    sink.int(4, "version")
    sink.int(sdim, "sdim")
    sink.int(len(mesh.points), "number of mesh vertices")
    sink.int(0, "lowest mesh vertex index")
    sink.blank()
    sink.comment("Mesh vertex coordinates")
    points = np.asarray(mesh.points, dtype=np.float64)
    if isinstance(sink, _BinarySink):
        sink.reals(points.ravel())
    else:
        for p in points:
            sink.reals(p)
    sink.blank()
    sink.int(len(mesh.cells), "number of element types")
    for b, block in enumerate(mesh.cells):
        data = np.asarray(block.data, dtype=np.int64)
        sink.blank()
        sink.comment(f"Type #{b}")
        sink.string(ctypes[b], "type name")
        nn = _NUM_NODES[block.type]
        sink.int(nn, "number of vertices per element")
        sink.int(len(data), "number of elements")
        sink.comment("Elements")
        order = node_order("mphtxt", block.type)
        if order is not None:
            data = data[:, list(order.from_meshio)]
        if isinstance(sink, _BinarySink):
            sink.ints(data.ravel())
        else:
            for row in data:
                sink.ints(row)
        sink.blank()
        sink.int(len(data), "number of geometric entity indices")
        sink.comment("Geometric entity indices")
        ent = entity[bases[b] : bases[b + 1]]
        if isinstance(sink, _BinarySink):
            sink.ints(ent)
        else:
            for x in ent:
                sink.ints([x])
    for k, (label, dim, ents) in enumerate(selections):
        sink.blank()
        sink.comment(f"--------- Object {k + 1} ----------")
        sink.int(0)
        sink.int(0)
        sink.int(1)
        sink.string("Selection", "class")
        sink.int(0, "Version")
        sink.string(label, "Label")
        sink.string("mesh1", "Geometry/mesh tag")
        sink.int(dim, "Dimension")
        sink.int(len(ents), "Number of entities")
        sink.comment("Entities")
        if isinstance(sink, _BinarySink):
            sink.ints(ents)
        else:
            for x in ents:
                sink.ints([x])


def _check_data(mesh, fmt):
    dropped = len(mesh.point_data) + len(mesh.field_data)
    dropped += sum(1 for name in mesh.cell_data if name != "mphtxt:geom")
    if dropped:
        warn(f"{fmt}: a COMSOL mesh holds no data arrays; {dropped} array(s) dropped")


def write(filename, mesh):
    sink = _TextSink()
    sink.out.append(_provenance.render_lines(_provenance.SlotTier.BLOCK, "# "))
    _serialise(sink, mesh, "mphtxt")
    _check_data(mesh, "mphtxt")
    text = "".join(sink.out)
    if is_buffer(filename, "w"):
        filename.write(text)
        return
    with open_file(filename, "wb") as f:
        f.write(text.encode("utf-8"))


def write_binary(filename, mesh):
    # no comment slot: consumes the pending record, and fails under
    # Mode.REQUIRED like the other slotless formats
    _provenance.lines(_provenance.SlotTier.NONE)
    sink = _BinarySink("mphbin")
    _serialise(sink, mesh, "mphbin")
    _check_data(mesh, "mphbin")
    if is_buffer(filename, "wb"):
        filename.write(b"".join(sink.out))
        return
    with open_file(filename, "wb") as f:
        f.write(b"".join(sink.out))
