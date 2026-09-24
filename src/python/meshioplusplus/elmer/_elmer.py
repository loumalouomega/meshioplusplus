"""Elmer mesh directory: the pure-Python reference engine.

The twin of ``src/cpp/src/formats/elmer.cpp``: both read the same meshes and
write the same bytes. An Elmer mesh is a directory holding ``mesh.header``,
``mesh.nodes``, ``mesh.elements``, ``mesh.boundary`` and optionally
``mesh.names``; a partitioned one adds ``partitioning.N/part.n.*``. See
``doc/formats/elmer.md``.
"""

from __future__ import annotations

import os
import pathlib
import struct

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._facets import FacetIndex, facet_nodes
from .._mesh import CellBlock, Mesh
from .._node_order import node_order
from .._regions import Region

# Elmer type code <-> meshio++ cell type (elmer.cpp's kElmTypes).
_TYPES = {
    101: "vertex",
    202: "line",
    203: "line3",
    204: "line4",
    303: "triangle",
    306: "triangle6",
    310: "triangle10",
    404: "quad",
    408: "quad8",
    409: "quad9",
    504: "tetra",
    510: "tetra10",
    605: "pyramid",
    613: "pyramid13",
    706: "wedge",
    715: "wedge15",
    718: "wedge18",
    808: "hexahedron",
    820: "hexahedron20",
    827: "hexahedron27",
}
_CODES = {t: c for c, t in _TYPES.items()}
_DIMS = {
    "vertex": 0,
    "line": 1,
    "triangle": 2,
    "quad": 2,
    "tetra": 3,
    "pyramid": 3,
    "wedge": 3,
    "hexahedron": 3,
}
_CORNERS = {
    "vertex": 1,
    "line": 2,
    "triangle": 3,
    "quad": 4,
    "tetra": 4,
    "pyramid": 5,
    "wedge": 6,
    "hexahedron": 8,
}


def _family(cell_type):
    return cell_type.rstrip("0123456789")


def _dim(cell_type):
    return _DIMS[_family(cell_type)]


def _fail(path, line, what):
    raise ReadError(f"Elmer mesh: {path}:{line}: {what}")


def _int(token):
    try:
        return int(token)
    except ValueError:
        return None


def _elm_id(token):
    """``id`` or ``id/part`` -> (id, 1-based owning part or -1)."""
    head, sep, tail = token.partition("/")
    ident = _int(head)
    if ident is None:
        return None
    if not sep:
        return ident, -1
    owner = _int(tail)
    return None if owner is None else (ident, owner)


def _lines(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for number, line in enumerate(fh, 1):
                tokens = line.split()
                if tokens:
                    yield number, tokens
    except OSError as exc:
        raise ReadError(f"Elmer mesh: cannot open {path}") from exc


class _Records:
    __slots__ = ("ids", "tags", "codes", "parts", "parents", "nodes")

    def __init__(self):
        self.ids = []
        self.tags = []
        self.codes = []
        self.parts = []
        self.parents = []
        self.nodes = []


class _ElmMesh:
    def __init__(self):
        self.node_ids = []
        self.coords = []
        self.node_index = {}
        self.bulk = _Records()
        self.boundary = _Records()
        self.body_names = {}
        self.boundary_names = {}
        self.has_parts = False


class _Binary:
    """A binary mesh file as ElmerGrid -bin writes it and ElmerSolver reads it
    (fem/src/MeshIO.F90, stream access): 32-bit native-endian integers and 64-bit
    (``.bin``) or 32-bit (``.sbin``, nodes only) reals. The byte order is taken
    from the first record's id, which is small and positive."""

    def __init__(self, path):
        self.path = path
        with open(path, "rb") as fh:
            self.data = fh.read()
        self.pos = 0
        self.record = 0
        order = "<"
        if len(self.data) >= 4:
            (le,) = struct.unpack("<I", self.data[:4])
            (be,) = struct.unpack(">I", self.data[:4])
            plausible = lambda v: 1 <= v < (1 << 30)  # noqa: E731
            # the byte order giving the smaller positive id wins (a big-endian 1
            # read little-endian is 2**24, still a possible id)
            if plausible(be) and (not plausible(le) or be < le):
                order = ">"
        self.order = order

    def at_end(self):
        return self.pos >= len(self.data)

    def _take(self, fmt, n):
        if self.pos + n > len(self.data):
            _fail(self.path, self.record + 1, "the binary file ends inside a record")
        (v,) = struct.unpack(self.order + fmt, self.data[self.pos : self.pos + n])
        self.pos += n
        return v

    def int(self):
        return self._take("i", 4)

    def real(self, single):
        return self._take("f", 4) if single else self._take("d", 8)


def _form(directory, stem, nodes):
    """Which form of ``stem`` a mesh directory holds: text, bin or sbin."""
    path = os.path.join(directory, stem)
    if os.path.isfile(path):
        return "text"
    if os.path.isfile(path + ".bin"):
        return "bin"
    if nodes and os.path.isfile(path + ".sbin"):
        return "sbin"
    raise ReadError(f"Elmer mesh: missing {path}")


def _read_nodes(directory, stem, mesh):
    def add(ident, xyz):
        # A shared node is listed by every part that uses it: keep the first.
        if ident in mesh.node_index:
            return
        mesh.node_index[ident] = len(mesh.node_ids)
        mesh.node_ids.append(ident)
        mesh.coords.append(xyz)

    form = _form(directory, stem, True)
    if form != "text":
        path = os.path.join(directory, stem + "." + form)
        b = _Binary(path)
        while not b.at_end():
            ident = b.int()
            if ident <= 0:
                _fail(path, b.record + 1, f"bad node id {ident}")
            add(ident, [b.real(form == "sbin") for _ in range(3)])
            b.record += 1
        return
    path = os.path.join(directory, stem)
    for number, tok in _lines(path):
        if len(tok) < 5:
            _fail(path, number, "a node line needs `id part x y z`")
        ident = _int(tok[0])
        if ident is None:
            _fail(path, number, f"bad node id '{tok[0]}'")
        try:
            xyz = [float(t) for t in tok[2:5]]
        except ValueError:
            _fail(path, number, "bad coordinate")
        add(ident, xyz)


def _read_elements(directory, stem, boundary, part, lenient, out, seen, seen_sides):
    """Text records ``id body type nodes`` / ``id boundary parent1 parent2 type
    nodes``; their binary forms add the owning part after the id."""
    state = {"skipped": 0}

    def accept(path, number, elem_id, owner, tag, parent, code, nodes):
        if code not in _TYPES:
            if not lenient:
                _fail(path, number, f"element type {code} has no meshio++ cell type")
            state["skipped"] += 1
            return
        elem_part = owner - 1 if owner > 0 else part
        if part >= 0:
            if boundary:
                key = (tag, tuple(sorted(nodes)))
                if key in seen_sides:
                    return
                seen_sides.add(key)
            else:
                if elem_id in seen:
                    if owner < 0:
                        out.parts[seen[elem_id]] = part
                    return
                seen[elem_id] = len(out.ids)
        out.ids.append(elem_id)
        out.tags.append(tag)
        out.codes.append(code)
        out.parts.append(elem_part)
        out.parents.append(parent if boundary else 0)
        out.nodes.append(nodes)

    if _form(directory, stem, False) == "bin":
        path = os.path.join(directory, stem + ".bin")
        b = _Binary(path)
        while not b.at_end():
            number = b.record + 1
            elem_id = b.int()
            owner = b.int()
            tag = b.int()
            parent = b.int() if boundary else 0
            if boundary:
                b.int()  # the second parent
            code = b.int()
            if elem_id <= 0 or code <= 100:
                _fail(path, number, "malformed record")
            nodes = [b.int() for _ in range(code % 100)]
            accept(
                path,
                number,
                elem_id,
                owner if owner > 0 else -1,
                tag,
                parent,
                code,
                nodes,
            )
            b.record += 1
    else:
        path = os.path.join(directory, stem)
        lead = 5 if boundary else 3
        for number, tok in _lines(path):
            if len(tok) < lead:
                _fail(
                    path,
                    number,
                    (
                        "a boundary line needs `id boundary parent1 parent2 type nodes`"
                        if boundary
                        else "an element line needs `id body type nodes`"
                    ),
                )
            ident = _elm_id(tok[0])
            tag = _int(tok[1])
            code = _int(tok[lead - 1])
            if ident is None or tag is None or code is None:
                _fail(path, number, "malformed line")
            count = code % 100
            if len(tok) != lead + count:
                _fail(path, number, f"type {code} needs {count} nodes")
            if code not in _TYPES and not lenient:
                _fail(path, number, f"element type {code} has no meshio++ cell type")
            nodes = []
            for t in tok[lead:]:
                n = _int(t)
                if n is None:
                    _fail(path, number, f"bad node id '{t}'")
                nodes.append(n)
            parent = (_int(tok[2]) or 0) if boundary else 0
            accept(path, number, ident[0], ident[1], tag, parent, code, nodes)
    if state["skipped"]:
        warn(
            f"Elmer mesh: {state['skipped']} element(s) of types with no meshio++ cell "
            f"type skipped in {os.path.join(directory, stem)}"
        )


def _element_ids(directory, stem):
    """(id, owner) of every bulk element of ``stem``, text or binary."""
    if _form(directory, stem, False) == "bin":
        b = _Binary(os.path.join(directory, stem + ".bin"))
        while not b.at_end():
            elem_id = b.int()
            owner = b.int()
            b.int()  # body
            code = b.int()
            for _ in range(code % 100):
                b.int()
            yield elem_id, (owner if owner > 0 else -1)
            b.record += 1
        return
    for _, tok in _lines(os.path.join(directory, stem)):
        ident = _elm_id(tok[0])
        if ident is not None:
            yield ident


def _read_names(path, mesh):
    if not os.path.isfile(path):
        return
    bodies = True
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            lower = line.lower()
            dollar = line.find("$")
            equals = line.find("=", dollar if dollar >= 0 else 0)
            if dollar < 0 or equals < 0:
                if "names for bound" in lower:
                    bodies = False
                elif "names for bod" in lower:
                    bodies = True
                continue
            name = line[dollar + 1 : equals].strip(" \t")
            tokens = line[equals + 1 :].split()
            ident = _int(tokens[0]) if tokens else None
            if ident is None or not name:
                continue
            (mesh.body_names if bodies else mesh.boundary_names)[ident] = name


def _count_parts(directory):
    n = 0
    while os.path.isfile(os.path.join(directory, f"part.{n + 1}.header")):
        n += 1
    return n


def _partition_dirs(directory):
    out = []
    try:
        entries = sorted(os.listdir(directory))
    except OSError:
        return out
    for name in entries:
        path = os.path.join(directory, name)
        if (
            name.startswith("partitioning.")
            and os.path.isdir(path)
            and _count_parts(path) > 0
        ):
            out.append(path)
    return out


def _read_serial(directory, lenient, mesh):
    seen, seen_sides = {}, set()
    _read_nodes(directory, "mesh.nodes", mesh)
    _read_elements(
        directory, "mesh.elements", False, -1, lenient, mesh.bulk, seen, seen_sides
    )
    _read_elements(
        directory, "mesh.boundary", True, -1, lenient, mesh.boundary, seen, seen_sides
    )
    _read_names(os.path.join(directory, "mesh.names"), mesh)


def _read_parts(directory, first, last, lenient, mesh):
    seen, seen_sides = {}, set()
    for p in range(first, last + 1):
        stem = f"part.{p + 1}"
        _read_nodes(directory, stem + ".nodes", mesh)
        _read_elements(
            directory,
            stem + ".elements",
            False,
            p,
            lenient,
            mesh.bulk,
            seen,
            seen_sides,
        )
        _read_elements(
            directory,
            stem + ".boundary",
            True,
            p,
            lenient,
            mesh.boundary,
            seen,
            seen_sides,
        )
    mesh.has_parts = True
    _read_names(
        os.path.join(os.path.dirname(os.path.normpath(directory)), "mesh.names"), mesh
    )


def _label_serial(part_dir, mesh):
    part_of = {}
    for p in range(_count_parts(part_dir)):
        stem = f"part.{p + 1}.elements"
        path = os.path.join(part_dir, stem)
        if not os.path.isfile(path) and not os.path.isfile(path + ".bin"):
            return
        for elem_id, owner in _element_ids(part_dir, stem):
            if owner < 0 or elem_id not in part_of:
                part_of[elem_id] = owner - 1 if owner > 0 else p
    labels = []
    for elem_id in mesh.bulk.ids:
        if elem_id not in part_of:
            warn(
                f"Elmer mesh: {part_dir} does not match the serial mesh; "
                "no partition labels"
            )
            return
        labels.append(part_of[elem_id])
    mesh.bulk.parts = labels
    mesh.boundary.parts = [part_of.get(p, -1) for p in mesh.boundary.parents]
    mesh.has_parts = True


def _build(src, where):
    points = np.array(src.coords, dtype=np.float64).reshape(-1, 3)
    cells = []
    part_blocks = []
    body_cells, boundary_cells = {}, {}
    body_dim, boundary_dim = {}, {}
    glob = 0
    for boundary in (False, True):
        rec = src.boundary if boundary else src.bulk
        codes = list(dict.fromkeys(rec.codes))
        for code in codes:
            cell_type = _TYPES[code]
            order = node_order("elmer", cell_type)
            dim = _dim(cell_type)
            rows = [e for e, c in enumerate(rec.codes) if c == code]
            conn = np.empty((len(rows), code % 100), dtype=np.int64)
            parts = np.empty(len(rows), dtype=np.int64)
            for r, e in enumerate(rows):
                src_nodes = rec.nodes[e]
                if order is not None:
                    src_nodes = [src_nodes[k] for k in order.to_meshio]
                try:
                    conn[r] = [src.node_index[n] for n in src_nodes]
                except KeyError as exc:
                    kind = "boundary element" if boundary else "element"
                    raise ReadError(
                        f"Elmer mesh: {where}: {kind} {rec.ids[e]} names undefined "
                        f"node {exc.args[0]}"
                    ) from None
                parts[r] = rec.parts[e]
                store = boundary_cells if boundary else body_cells
                dims = boundary_dim if boundary else body_dim
                store.setdefault(rec.tags[e], []).append(glob)
                dims[rec.tags[e]] = max(dims.get(rec.tags[e], dim), dim)
                glob += 1
            cells.append(CellBlock(cell_type, conn))
            part_blocks.append(parts)
    mesh = Mesh(points, cells)
    if src.has_parts:
        mesh.cell_data["partition:part"] = part_blocks

    body_names = {n for i, n in src.body_names.items() if i in body_cells}
    clash = {
        n
        for i, n in src.boundary_names.items()
        if i in boundary_cells and n in body_names
    }
    regions = []
    for boundary in (False, True):
        store = boundary_cells if boundary else body_cells
        names = src.boundary_names if boundary else src.body_names
        dims = boundary_dim if boundary else body_dim
        for ident in sorted(store):
            named = names.get(ident)
            if named is None:
                name = f"{'boundary' if boundary else 'body'}_{ident}"
            elif named in clash:
                name = f"{'boundary' if boundary else 'body'}:{named}"
            else:
                name = named
            regions.append(
                Region(
                    name,
                    "cell",
                    np.array(store[ident], dtype=np.int64),
                    dim=dims[ident],
                    tag=ident,
                )
            )
    regions.sort(key=lambda r: (r.name, r.dim, r.tag))
    mesh.regions = regions
    return mesh


def read(filename, points_only=False, arrays=None, piece=None, lenient=False):
    path = str(filename)
    directory = path
    if os.path.basename(path) == "mesh.header" and os.path.isfile(path):
        directory = os.path.dirname(path) or "."
    if not os.path.isdir(directory):
        raise ReadError(f"Elmer mesh: '{path}' is not a mesh directory")
    mesh = _ElmMesh()
    serial = os.path.isfile(os.path.join(directory, "mesh.header"))
    name = os.path.basename(os.path.normpath(directory))
    part_dir = None
    if name.startswith("partitioning.") and _count_parts(directory) > 0:
        part_dir = directory
    else:
        parts = _partition_dirs(directory)
        if len(parts) == 1:
            part_dir = parts[0]
        elif len(parts) > 1 and (not serial or piece is not None):
            raise ReadError(
                f"Elmer mesh: {directory} holds several partitioning directories; "
                "name one of them"
            )
        elif len(parts) > 1:
            warn(
                f"Elmer mesh: {directory} holds several partitioning directories; "
                "no partition labels"
            )
        if not serial and part_dir is None:
            raise ReadError(
                f"Elmer mesh: {directory} holds neither mesh.header nor a "
                "partitioning directory"
            )
    if piece is not None:
        if part_dir is None:
            raise ReadError(
                f"Elmer mesh: {directory} is not partitioned; it has no pieces"
            )
        count = _count_parts(part_dir)
        index = piece + count if piece < 0 else piece
        if not 0 <= index < count:
            raise ReadError(
                f"piece {piece} is out of range: the file has {count} piece(s)"
            )
        _read_parts(part_dir, index, index, lenient, mesh)
        out = _build(mesh, part_dir)
    elif serial:
        _read_serial(directory, lenient, mesh)
        if part_dir is not None:
            _label_serial(part_dir, mesh)
        out = _build(mesh, directory)
    else:
        _read_parts(part_dir, 0, _count_parts(part_dir) - 1, lenient, mesh)
        out = _build(mesh, part_dir)
    return out


# -- writing ---------------------------------------------------------------------


def _clean_name(name):
    for c in "$=\n\r!":
        name = name.replace(c, "_")
    return name


def _padded(value):
    return f"{value:<6d}"


class _Ids:
    def __init__(self):
        self.taken = set()

    def take(self, tag):
        if tag > 0 and tag not in self.taken:
            self.taken.add(tag)
            return tag
        ident = max(self.taken) + 1 if self.taken else 1
        self.taken.add(ident)
        return ident


def write(filename, mesh):
    directory = pathlib.Path(filename)
    blocks = mesh.cells
    codes, dims = [], []
    bulk_dim = -1
    for block in blocks:
        code = -1 if isinstance(block.data, list) else _CODES.get(block.type, -1)
        if code < 0:
            raise WriteError(
                f"Elmer mesh writer: cell type '{block.type}' has no Elmer element code"
            )
        codes.append(code)
        dims.append(_dim(block.type))
        if len(block.data):
            bulk_dim = max(bulk_dim, dims[-1])
    if bulk_dim < 1:
        raise WriteError(
            "Elmer mesh writer: the mesh has no cells of dimension 1 or more"
        )
    bases = [0]
    for block in blocks:
        bases.append(bases[-1] + len(block.data))
    n_cells = bases[-1]
    is_bulk = np.zeros(n_cells, dtype=bool)
    for b in range(len(blocks)):
        if dims[b] == bulk_dim:
            is_bulk[bases[b] : bases[b + 1]] = True

    cell_id = np.zeros(n_cells, dtype=np.int64)
    body_ids, boundary_ids = _Ids(), _Ids()
    bodies, boundaries = [], []
    overlaps = point_regions = 0
    side_facets, side_ids = [], []
    kind_order = {"point": 0, "cell": 1, "side": 2}
    for reg in sorted(
        mesh.regions, key=lambda r: (kind_order[r.kind], r.name, r.dim, r.tag)
    ):
        entries = np.asarray(reg.entries, dtype=np.int64)
        if reg.kind == "point":
            point_regions += 1
            continue
        if reg.kind == "side":
            if len(entries) == 0:
                continue
            ident = boundary_ids.take(reg.tag)
            boundaries.append((ident, reg.name))
            for cell, facet in entries.reshape(-1, 2).tolist():
                side_facets.append((cell, facet))
                side_ids.append(ident)
            continue
        valid = [int(g) for g in entries.tolist() if 0 <= g < n_cells]
        bulk = [g for g in valid if is_bulk[g]]
        lower = [g for g in valid if not is_bulk[g]]
        for boundary in (False, True):
            group = lower if boundary else bulk
            if not group:
                continue
            free_cells = []
            for g in group:
                if cell_id[g] != 0:
                    overlaps += 1
                else:
                    free_cells.append(g)
            if not free_cells:
                continue
            name = reg.name
            if bulk and lower and boundary:
                name += "_boundary"
            ident = (boundary_ids if boundary else body_ids).take(reg.tag)
            (boundaries if boundary else bodies).append((ident, name))
            cell_id[free_cells] = ident
    for b in range(len(blocks)):
        ident = 0
        for g in range(bases[b], bases[b + 1]):
            if cell_id[g] != 0:
                continue
            if ident == 0:
                is_body = dims[b] == bulk_dim
                ident = (body_ids if is_body else boundary_ids).take(-1)
                (bodies if is_body else boundaries).append((ident, ""))
            cell_id[g] = ident

    if overlaps:
        warn(
            f"Elmer mesh writer: {overlaps} cell(s) belong to more than one region; "
            "the first region (in name order) wins"
        )
    if point_regions:
        warn(
            f"Elmer mesh writer: {point_regions} point region(s) have no Elmer "
            "equivalent and were dropped"
        )
        _provenance.note(
            "regions-dropped",
            f"{point_regions} point region(s) have no Elmer equivalent",
        )
    partitioned = "partition:part" in mesh.cell_data
    other_cell_data = [k for k in mesh.cell_data if k != "partition:part"]
    if mesh.point_data or other_cell_data or mesh.field_data:
        warn(
            "Elmer mesh writer: an Elmer mesh holds no data arrays; point, cell and "
            "field data dropped"
        )
        _provenance.note("data-dropped", "an Elmer mesh directory holds no data arrays")

    element_no = np.zeros(n_cells, dtype=np.int64)
    element_no[is_bulk] = np.arange(1, int(is_bulk.sum()) + 1)
    n_bulk = int(is_bulk.sum())

    facets = None
    node_cells = None
    orphans = 0

    def bulk_parent(g):
        return int(element_no[g]) if g >= 0 and is_bulk[g] else 0

    def parents_of(corners, dim):
        nonlocal facets, node_cells
        if dim == bulk_dim - 1 and dim >= 1:
            if facets is None:
                facets = FacetIndex(
                    mesh, solid_faces=bulk_dim == 3, surface_edges=bulk_dim == 2
                )
            hit = facets.find(corners)
            if hit is None:
                return 0, 0
            return (
                bulk_parent(hit.first[0]),
                bulk_parent(hit.second[0]) if hit.count > 1 else 0,
            )
        if node_cells is None:
            node_cells = [[] for _ in range(len(mesh.points))]
            for b, block in enumerate(blocks):
                if dims[b] != bulk_dim:
                    continue
                for r, row in enumerate(np.asarray(block.data).tolist()):
                    g = bases[b] + r
                    for p in row:
                        lst = node_cells[p]
                        if not lst or lst[-1] != g:
                            lst.append(g)
        for g in node_cells[corners[0]]:
            if all(g in node_cells[c] for c in corners[1:]):
                return int(element_no[g]), 0
        return 0, 0

    directory.mkdir(parents=True, exist_ok=True)
    for entry in sorted(os.listdir(directory)):
        if entry.startswith("partitioning.") and not partitioned:
            warn(
                f"Elmer mesh writer: {directory / entry} is left over from an earlier "
                "mesh and no longer matches it"
            )

    points = np.asarray(mesh.points, dtype=np.float64)
    pdim = points.shape[1] if points.ndim == 2 else 0
    out = []
    for p, row in enumerate(points.tolist()):
        xyz = [row[d] if d < pdim else 0.0 for d in range(3)]
        out.append(f"{p + 1} -1 " + " ".join("%.17g" % v for v in xyz) + "\n")
    (directory / "mesh.nodes").write_bytes("".join(out).encode())

    type_counts = {}
    elements, boundary = [], []
    n_boundary = 0
    bulk_rows, boundary_rows = [], []  # kept for the partitioned copy
    for b, block in enumerate(blocks):
        order = node_order("elmer", block.type)
        n_corners = _CORNERS[_family(block.type)]
        is_body = dims[b] == bulk_dim
        for r, row in enumerate(np.asarray(block.data).tolist()):
            g = bases[b] + r
            if order is not None:
                written = [row[k] for k in order.from_meshio]
            else:
                written = row
            nodes = " ".join(str(v + 1) for v in written)
            if is_body:
                elements.append(f"{element_no[g]} {cell_id[g]} {codes[b]} {nodes}\n")
                bulk_rows.append((g, cell_id[g], codes[b], [v + 1 for v in written]))
            else:
                p1, p2 = parents_of(row[:n_corners], dims[b])
                orphans += p1 == 0
                n_boundary += 1
                boundary.append(
                    f"{n_boundary} {cell_id[g]} {p1} {p2} {codes[b]} {nodes}\n"
                )
                boundary_rows.append(
                    (cell_id[g], p1, p2, codes[b], [v + 1 for v in written])
                )
            type_counts[codes[b]] = type_counts.get(codes[b], 0) + 1
    bad_facets = 0
    for (cell, facet), ident in zip(side_facets, side_ids):
        found = facet_nodes(mesh, cell, facet)
        if found is None:
            bad_facets += 1
            continue
        ftype, row = found
        code = _CODES[ftype]
        order = node_order("elmer", ftype)
        n_corners = _CORNERS[_family(ftype)]
        p1, p2 = parents_of(row[:n_corners], _dim(ftype))
        orphans += p1 == 0
        n_boundary += 1
        written = [row[k] for k in order.from_meshio] if order is not None else row
        nodes = " ".join(str(v + 1) for v in written)
        boundary.append(f"{n_boundary} {ident} {p1} {p2} {code} {nodes}\n")
        boundary_rows.append((ident, p1, p2, code, [v + 1 for v in written]))
        type_counts[code] = type_counts.get(code, 0) + 1
    if bad_facets:
        warn(
            f"Elmer mesh writer: {bad_facets} side region entr(ies) name no facet and "
            "were dropped"
        )
    if orphans:
        warn(
            f"Elmer mesh writer: {orphans} boundary element(s) lie on no bulk "
            "element; written with parent 0"
        )
    (directory / "mesh.elements").write_bytes("".join(elements).encode())
    (directory / "mesh.boundary").write_bytes("".join(boundary).encode())

    header = [
        f"{_padded(len(points))} {_padded(n_bulk)} {_padded(n_boundary)}\n",
        f"{_padded(len(type_counts))}\n",
    ]
    for code in sorted(type_counts):
        header.append(f"{_padded(code)} {_padded(type_counts[code])}\n")
    (directory / "mesh.header").write_bytes("".join(header).encode())

    names = [
        _provenance.render_lines(_provenance.SlotTier.BLOCK, "! ").replace("$", "_")
    ]
    body_names = {_clean_name(n) for _, n in bodies if n}
    for is_boundary in (False, True):
        group = sorted(boundaries if is_boundary else bodies, key=lambda g: g[0])
        names.append(
            "! ----- names for boundaries -----\n"
            if is_boundary
            else "! ----- names for bodies -----\n"
        )
        for ident, name in group:
            if not name:
                continue
            name = _clean_name(name)
            if is_boundary and name in body_names:
                name += "_boundary"
            names.append(f"$ {name} = {ident}\n")
    (directory / "mesh.names").write_bytes("".join(names).encode())

    if partitioned:
        _write_partitioning(
            directory, mesh, bases, element_no, n_bulk, bulk_rows, boundary_rows, points
        )


def _write_partitioning(
    directory, mesh, bases, element_no, n_bulk, bulk_rows, boundary_rows, points
):
    """ElmerGrid's partitioning.N without halos (``write_elmer`` in elmer.cpp):
    ``partition:part`` (0-based) places each bulk element; a node belongs to
    every part using it and is owned by the lowest; a boundary element goes to
    each part holding one of its parents, the other parent set to 0."""
    labels = mesh.cell_data["partition:part"]
    part_of = [-1] * (n_bulk + 1)
    n_parts = 0
    for g, _, _, _ in bulk_rows:
        b = int(np.searchsorted(bases, g, side="right")) - 1
        label = int(labels[b][g - bases[b]])
        if label < 0:
            raise WriteError(
                f"Elmer mesh writer: partition:part has a negative part for cell {g}"
            )
        part_of[element_no[g]] = label
        n_parts = max(n_parts, label + 1)
    users = {}
    for g, _, _, nodes in bulk_rows:
        part = part_of[element_no[g]]
        for n in nodes:
            users.setdefault(n, set()).add(part)
    pdir = directory / f"partitioning.{n_parts}"
    pdir.mkdir(parents=True, exist_ok=True)
    pdim = points.shape[1] if points.ndim == 2 else 0
    empty = 0
    for part in range(n_parts):
        stem = f"part.{part + 1}"
        elems, sides, node_lines, shared = [], [], [], []
        bulk_types, side_types = {}, {}
        nodes = set()
        for g, tag, code, row in bulk_rows:
            no = int(element_no[g])
            if part_of[no] != part:
                continue
            elems.append(f"{no} {tag} {code} " + " ".join(str(n) for n in row) + "\n")
            nodes.update(row)
            bulk_types[code] = bulk_types.get(code, 0) + 1
        for side_no, (tag, p1, p2, code, row) in enumerate(boundary_rows, start=1):
            in1 = p1 > 0 and part_of[p1] == part
            in2 = p2 > 0 and part_of[p2] == part
            if not (in1 or in2):
                continue
            sides.append(
                f"{side_no} {tag} {p1 if in1 else 0} {p2 if in2 else 0} {code} "
                + " ".join(str(n) for n in row)
                + "\n"
            )
            side_types[code] = side_types.get(code, 0) + 1
        for n in sorted(nodes):
            row = points[n - 1].tolist()
            xyz = [row[d] if d < pdim else 0.0 for d in range(3)]
            node_lines.append(f"{n} -1 " + " ".join("%.17g" % v for v in xyz) + "\n")
            parts = sorted(users[n])
            if len(parts) > 1:
                # `id count owner others`: every part using it but the owner
                shared.append(
                    f"{n} {len(parts)} " + " ".join(str(q + 1) for q in parts) + "\n"
                )
        empty += not elems
        header = [
            f"{_padded(len(nodes))} {_padded(len(elems))} {_padded(len(sides))}\n",
            f"{_padded(len(bulk_types) + len(side_types))}\n",
        ]
        for types in (bulk_types, side_types):
            for code in sorted(types):
                header.append(f"{_padded(code)} {_padded(types[code])}\n")
        header.append(f"{_padded(len(shared))} {_padded(0)}\n")
        (pdir / f"{stem}.header").write_bytes("".join(header).encode())
        (pdir / f"{stem}.nodes").write_bytes("".join(node_lines).encode())
        (pdir / f"{stem}.elements").write_bytes("".join(elems).encode())
        (pdir / f"{stem}.boundary").write_bytes("".join(sides).encode())
        (pdir / f"{stem}.shared").write_bytes("".join(shared).encode())
    if empty:
        warn(
            f"Elmer mesh writer: {empty} of the {n_parts} parts in partition:part hold "
            "no element; ElmerSolver needs every part populated"
        )
