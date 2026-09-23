"""FEBio input files (``.feb``): the pure-Python reference engine.

The twin of ``src/cpp/src/formats/febio.cpp``: both read the same meshes and
write the same bytes. The mesh of febio_spec 2.5 (``<Geometry>``), 3.0 and 4.0
(``<Mesh>``) is read; 4.0 is written. See ``doc/formats/febio.md``.
"""

from __future__ import annotations

import math
import os
import re
import xml.etree.ElementTree as ET

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._facets import FacetIndex, facet_nodes
from .._mesh import CellBlock, Mesh
from .._node_order import node_order
from .._regions import Region
from .._skin import _CELL_FACES
from .._surface import _CELL_EDGES

# FEBio element type -> (meshio++ type, node count); febio.cpp's kFebTypes.
_TYPES = {
    "tet4": ("tetra", 4),
    "ut4": ("tetra", 4),
    "tet10": ("tetra10", 10),
    "penta6": ("wedge", 6),
    "penta15": ("wedge15", 15),
    "pyra5": ("pyramid", 5),
    "pyra13": ("pyramid13", 13),
    "hex8": ("hexahedron", 8),
    "hex20": ("hexahedron20", 20),
    "hex27": ("hexahedron27", 27),
    "quad4": ("quad", 4),
    "q4eas": ("quad", 4),
    "q4ans": ("quad", 4),
    "q4s": ("quad", 4),
    "quad8": ("quad8", 8),
    "quad9": ("quad9", 9),
    "tri3": ("triangle", 3),
    "tri3s": ("triangle", 3),
    "tri6": ("triangle6", 6),
    "tri7": ("triangle7", 7),
    "tri10": ("triangle10", 10),
    "line2": ("line", 2),
    "truss2": ("line", 2),
    "line3": ("line3", 3),
}
# FEBio types with no meshio++ cell: node count and the lenient downgrade.
_LOSSY = {"tet5": (5, "tetra"), "tet15": (15, "tetra10"), "tet20": (20, None)}
_WRITE = {
    "line": "line2",
    "line3": "line3",
    "triangle": "tri3",
    "triangle6": "tri6",
    "triangle7": "tri7",
    "quad": "quad4",
    "quad8": "quad8",
    "quad9": "quad9",
    "tetra": "tet4",
    "tetra10": "tet10",
    "pyramid": "pyra5",
    "pyramid13": "pyra13",
    "wedge": "penta6",
    "wedge15": "penta15",
    "hexahedron": "hex8",
    "hexahedron20": "hex20",
    "hexahedron27": "hex27",
}
_FAMILY = {
    "vertex": (1, 0),
    "line": (2, 1),
    "triangle": (3, 2),
    "quad": (4, 2),
    "tetra": (4, 3),
    "pyramid": (5, 3),
    "wedge": (6, 3),
    "hexahedron": (8, 3),
}
_NODES = {t: n for t, n in _TYPES.values()}
_NODES["vertex"] = 1
_KIND_ORDER = {"point": 0, "cell": 1, "side": 2}


def _corners(cell_type):
    return _FAMILY[cell_type.rstrip("0123456789")][0]


def _dim(cell_type):
    return _FAMILY[cell_type.rstrip("0123456789")][1]


def _fail(what):
    raise ReadError(f"FEBio .feb: {what}")


def _split(text):
    return [t for t in re.split(r"[,\s]+", text or "") if t]


def _int(token, what):
    if not re.fullmatch(r"-?[0-9]+", token):
        _fail(f"bad {what} '{token}'")
    return int(token)


def _float(token):
    try:
        return float(token)
    except ValueError:
        _fail(f"bad number '{token}'")


def _base_type(name):
    """``TET10G4`` -> ``tet10``: drop an integration-rule suffix."""
    t = name.lower()
    m = re.fullmatch(r"(.*[0-9])g[0-9]+", t)
    return m.group(1) if m else t


def _id_list(text):
    """Ids, with FEBio 4's ``a:b[:step]`` ranges."""
    out = []
    for t in _split(text):
        if ":" not in t:
            out.append(_int(t, "id"))
            continue
        parts = t.split(":")
        a, b = _int(parts[0], "range"), _int(parts[1], "range")
        step = _int(parts[2], "range step") if len(parts) > 2 else 1
        if step <= 0:
            _fail(f"bad range '{t}'")
        out.extend(range(a, b + 1, step))
    return out


class _Reader:
    def __init__(self, lenient):
        self.lenient = lenient
        self.coords = []
        self.node_index = {}
        self.block_types = []
        self.block_conn = []
        self.elem_index = {}
        self.num_cells = 0
        self.node_sets = {}
        self.elem_sets = {}
        self.groups = []
        self.group_index = {}
        self.surfaces = []
        self.extra = []
        self.material_ids = {}
        self.domain_material = {}
        self.lossy = 0

    def group(self, name, kind):
        key = (kind, name)
        if key not in self.group_index:
            self.group_index[key] = len(self.groups)
            self.groups.append(
                {"name": name, "kind": kind, "dim": -1, "tag": -1, "e": []}
            )
        return self.groups[self.group_index[key]]

    def node(self, ident):
        try:
            return self.node_index[ident]
        except KeyError:
            _fail(f"undefined node {ident}")

    def element(self, ident):
        try:
            return self.elem_index[ident]
        except KeyError:
            _fail(f"undefined element {ident}")

    def read_nodes(self, nodes):
        members = []
        for n in nodes:
            ident = _int(n.get("id", ""), "node id")
            tokens = _split(n.text)
            if len(tokens) != 3:
                _fail(f"node {ident} needs x,y,z")
            if ident in self.node_index:
                _fail(f"node {ident} is defined twice")
            self.node_index[ident] = len(self.coords)
            members.append(len(self.coords))
            self.coords.append([_float(t) for t in tokens])
        name = nodes.get("name", "")
        if name:
            self.group(name, "point")["e"].extend(members)
            self.node_sets.setdefault(name, []).extend(members)

    def read_elements(self, elems, index):
        raw = elems.get("type", "")
        base = _base_type(raw)
        if base in _TYPES:
            meshio_type, nodes = _TYPES[base]
        else:
            if base not in _LOSSY:
                _fail(f"unknown element type '{raw}'")
            nodes, downgrade = _LOSSY[base]
            if not self.lenient or downgrade is None:
                hint = " (read with lenient to downgrade it)" if downgrade else ""
                _fail(f"element type '{raw}' has no meshio++ cell type{hint}")
            meshio_type = downgrade
            self.lossy += 1
        keep = _NODES[meshio_type]
        order = node_order("febio", meshio_type)
        name = elems.get("name", "") or f"Part{index + 1}"
        tag = -1
        mat = elems.get("mat", "")
        if mat:
            if mat in self.material_ids:
                tag = self.material_ids[mat]
            elif re.fullmatch(r"-?[0-9]+", mat):
                tag = int(mat)
        elif name in self.domain_material:
            tag = self.domain_material[name]
        conn = []
        members = []
        for e in elems:
            ident = _int(e.get("id", ""), "element id")
            tokens = _split(e.text)
            if len(tokens) != nodes:
                _fail(f"element {ident} of type {raw} needs {nodes} nodes")
            row = [self.node(_int(t, "node id")) for t in tokens]
            conn.append(
                [row[order.to_meshio[k]] if order else row[k] for k in range(keep)]
            )
            if ident in self.elem_index:
                _fail(f"element {ident} is defined twice")
            self.elem_index[ident] = self.num_cells
            members.append(self.num_cells)
            self.num_cells += 1
        self.block_types.append(meshio_type)
        self.block_conn.append(conn)
        g = self.group(name, "cell")
        g["dim"] = max(g["dim"], _dim(meshio_type))
        if tag >= 0:
            g["tag"] = tag
        g["e"].extend(members)
        self.elem_sets.setdefault(name, []).extend(members)

    def read_node_set(self, nset):
        name = nset.get("name", "")
        members = [self.node(i) for i in _id_list(nset.text)]
        for c in nset:
            if c.tag == "node_list":
                members.extend(self.node(i) for i in _id_list(c.text))
            elif c.tag == "NodeSet":
                other = c.get("node_set", "")
                if other not in self.node_sets:
                    _fail(f"NodeSet '{name}' includes undefined set '{other}'")
                members.extend(self.node_sets[other])
            else:
                members.append(self.node(_int(c.get("id", ""), "node id")))
        self.group(name, "point")["e"].extend(members)
        self.node_sets.setdefault(name, []).extend(members)

    def read_element_set(self, eset):
        name = eset.get("name", "")
        members = [self.element(i) for i in _id_list(eset.text)]
        for c in eset:
            members.append(self.element(_int(c.get("id", ""), "element id")))
        self.group(name, "cell")["e"].extend(members)
        self.elem_sets.setdefault(name, []).extend(members)

    def facet(self, f):
        entry = _TYPES.get(_base_type(f.tag))
        if entry is None or _dim(entry[0]) > 2:
            _fail(f"unknown facet type '{f.tag}'")
        tokens = _split(f.text)
        if len(tokens) != entry[1]:
            _fail(f"a {f.tag} facet needs {entry[1]} nodes")
        return entry[0], [self.node(_int(t, "node id")) for t in tokens]

    def read_surface(self, surface):
        self.surfaces.append(
            (surface.get("name", ""), [self.facet(f) for f in surface])
        )

    def read_lines(self, lines, discrete):
        name = lines.get("name", "")
        by_type = {}
        for c in lines:
            if discrete:
                tokens = _split(c.text)
                if len(tokens) != 2:
                    _fail("a <delem> needs two nodes")
                cell_type, nodes = "line", [
                    self.node(_int(t, "node id")) for t in tokens
                ]
            else:
                cell_type, nodes = self.facet(c)
            if cell_type not in by_type:
                by_type[cell_type] = len(self.extra)
                self.extra.append([cell_type, [], name])
            self.extra[by_type[cell_type]][1].append(nodes)


def _components(data_type, data):
    if not data_type:
        for c in data:
            return max(len(_split(c.text)), 1)
        return 1
    widths = {"scalar": 1, "vec2": 2, "vec3": 3, "mat3s": 6, "mat3": 9}
    if data_type not in widths:
        _fail(f"unknown data type '{data_type}'")
    return widths[data_type]


def _fill(data, members, width, out):
    for c in data:
        lid = _int(c.get("lid", ""), "lid")
        if lid < 1 or lid > len(members):
            _fail(f"lid {lid} is outside its set")
        tokens = _split(c.text)
        if len(tokens) != width:
            _fail(f"MeshData '{data.get('name', '')}' needs {width} values per entry")
        out[members[lid - 1]] = [_float(t) for t in tokens]


def read(filename, points_only=False, arrays=None, lenient=False):
    path = str(filename)
    try:
        tree = ET.parse(path)
    except ET.ParseError as exc:
        raise ReadError(f"FEBio .feb: could not parse {path}: {exc}") from None
    except OSError as exc:
        raise ReadError(f"FEBio .feb: could not parse {path}: {exc}") from None
    root = tree.getroot()
    if root.tag != "febio_spec":
        raise ReadError(f"FEBio .feb: {path} has no <febio_spec> root")
    version = root.get("version", "")
    if version not in ("2.5", "3.0", "4.0"):
        raise ReadError(
            f"FEBio .feb: febio_spec version '{version}' is not supported "
            "(2.5, 3.0 and 4.0 are)"
        )
    spec25 = version == "2.5"
    reader = _Reader(lenient)
    materials = root.find("Material")
    for m in materials.findall("material") if materials is not None else []:
        ident = m.get("id", "")
        if m.get("name") and re.fullmatch(r"-?[0-9]+", ident):
            reader.material_ids[m.get("name")] = int(ident)
    domains = root.find("MeshDomains")
    for d in domains if domains is not None else []:
        if d.get("mat", "") in reader.material_ids:
            reader.domain_material[d.get("name", "")] = reader.material_ids[
                d.get("mat")
            ]

    if spec25:
        mesh_node = root.find("Geometry")
    else:
        mesh_node = root.find("Mesh")
        if mesh_node is None:
            mesh_node = root.find("Geometry")
    if mesh_node is None:
        section = "<Geometry>" if spec25 else "<Mesh>"
        raise ReadError(f"FEBio .feb: {path} has no {section} section")
    source = mesh_node.get("from", "")
    if source:
        other = source
        if not os.path.isabs(other):
            other = os.path.join(os.path.dirname(path), other)
        try:
            other_root = ET.parse(other).getroot()
        except (ET.ParseError, OSError):
            raise ReadError(
                f"FEBio .feb: could not read '{other}' named by from="
            ) from None
        found = (
            other_root.find(mesh_node.tag) if other_root.tag == "febio_spec" else None
        )
        if found is None:
            raise ReadError(f"FEBio .feb: '{other}' has no mesh section")
        mesh_node = found

    n_blocks = 0
    for c in mesh_node:
        tag = c.tag
        if tag == "Nodes":
            reader.read_nodes(c)
        elif tag == "Elements":
            reader.read_elements(c, n_blocks)
            n_blocks += 1
        elif tag == "NodeSet":
            reader.read_node_set(c)
        elif tag == "ElementSet":
            reader.read_element_set(c)
        elif tag == "Surface":
            reader.read_surface(c)
        elif tag == "Edge":
            reader.read_lines(c, False)
        elif tag == "DiscreteSet":
            reader.read_lines(c, True)
        elif tag in ("Part", "Instance"):
            raise ReadError(
                "FEBio .feb: the <Part>/<Instance> form is not supported; export the "
                "model from FEBio Studio as a plain mesh"
            )
        elif tag not in ("SurfacePair", "PartList", "NodeSetPair", "NodeSetSet"):
            warn(f"FEBio .feb: <{tag}> in the mesh section is not read")
    if reader.lossy:
        warn(
            f"FEBio .feb: {reader.lossy} element block(s) downgraded to a meshio++ "
            "cell type"
        )

    points = np.array(reader.coords, dtype=np.float64).reshape(-1, 3)
    cells = [
        CellBlock(t, np.array(conn, dtype=np.int64).reshape(-1, _NODES[t]))
        for t, conn in zip(reader.block_types, reader.block_conn)
    ]
    mesh = Mesh(points, cells)

    faces = None
    detached = 0
    base = reader.num_cells
    surface_blocks = []
    for name, facets in reader.surfaces:
        if faces is None:
            faces = FacetIndex(mesh, surface_edges=False)
        sides = []
        for cell_type, nodes in facets:
            hit = faces.find(nodes[: _corners(cell_type)])
            if hit is None:
                sides = None
                break
            sides.append(hit.first)
        if sides is not None:
            g = reader.group(name, "side")
            g["dim"] = 2
            g["e"].extend(sides)
            continue
        detached += 1
        by_type = {}
        for cell_type, nodes in facets:
            if cell_type not in by_type:
                by_type[cell_type] = len(surface_blocks)
                surface_blocks.append([cell_type, [], name])
            surface_blocks[by_type[cell_type]][1].append(nodes)
    if detached:
        warn(
            f"FEBio .feb: {detached} surface(s) do not lie on solid elements; read "
            "as cell blocks"
        )
    for cell_type, rows, name in surface_blocks + reader.extra:
        mesh.cells.append(
            CellBlock(
                cell_type,
                np.array(rows, dtype=np.int64).reshape(-1, _NODES[cell_type]),
            )
        )
        g = reader.group(name, "cell")
        g["dim"] = max(g["dim"], _dim(cell_type))
        g["e"].extend(range(base, base + len(rows)))
        base += len(rows)

    sizes = [len(b.data) for b in mesh.cells]
    n_cells = sum(sizes)
    mesh_data = root.find("MeshData")
    for d in mesh_data if mesh_data is not None else []:
        data_type = d.get("data_type", "") or d.get("datatype", "")
        name = d.get("name", "") or d.get("var", "")
        if d.tag == "NodeData":
            members = reader.node_sets.get(d.get("node_set", ""))
            if members is None:
                _fail(f"NodeData '{name}' names an undefined node set")
            width = _components(data_type, d)
            values = np.full((len(points), width), math.nan)
            _fill(d, members, width, values)
            mesh.point_data[name] = values[:, 0] if width == 1 else values
        elif d.tag == "ElementData" and len(d):
            members = reader.elem_sets.get(d.get("elem_set", ""))
            if members is None:
                _fail(f"ElementData '{name}' names an undefined element set")
            width = _components(data_type, d)
            values = np.full((n_cells, width), math.nan)
            _fill(d, members, width, values)
            blocks = []
            start = 0
            for n in sizes:
                chunk = values[start : start + n]
                blocks.append(chunk[:, 0].copy() if width == 1 else chunk.copy())
                start += n
            mesh.cell_data[name] = blocks
        else:
            warn(f'FEBio .feb: <{d.tag} name="{name}"> in MeshData is not read')

    regions = []
    for g in reader.groups:
        if g["kind"] == "side":
            entries = np.array(g["e"], dtype=np.int64).reshape(-1, 2)
        else:
            entries = np.array(g["e"], dtype=np.int64)
        regions.append(
            Region(g["name"], g["kind"], entries, dim=g["dim"], tag=g["tag"])
        )
    regions.sort(key=lambda r: (_KIND_ORDER[r.kind], r.name, r.dim, r.tag))
    mesh.regions = regions
    return mesh


# -- writing ---------------------------------------------------------------------


def _escape(text):
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _ids(values):
    return ",".join(str(int(v)) for v in values)


def _all_on_edges(mesh, line, dims):
    """Whether every line of block ``line`` joins two corners adjacent in a
    solid's face or along a surface cell's edge."""
    edges = set()
    for b, block in enumerate(mesh.cells):
        if dims[b] < 2:
            continue
        rows = np.asarray(block.data).tolist()
        if dims[b] == 3:
            cycles = [local[:n] for _, n, local in _CELL_FACES.get(block.type, [])]
        else:
            cycles = None
            pairs = [local[:2] for _, _, local in _CELL_EDGES.get(block.type, [])]
        for row in rows:
            if cycles is not None:
                for cyc in cycles:
                    for c in range(len(cyc)):
                        a, z = row[cyc[c]], row[cyc[(c + 1) % len(cyc)]]
                        edges.add((min(a, z), max(a, z)))
            else:
                for i, j in pairs:
                    a, z = row[i], row[j]
                    edges.add((min(a, z), max(a, z)))
    for a, z in np.asarray(mesh.cells[line].data)[:, :2].tolist():
        if (min(a, z), max(a, z)) not in edges:
            return False
    return True


def write(filename, mesh):
    blocks = mesh.cells
    dims = []
    max_dim = 0
    for block in blocks:
        if isinstance(block.data, list) or (
            block.type != "vertex" and block.type not in _WRITE
        ):
            raise WriteError(
                f"FEBio .feb writer: cell type '{block.type}' has no FEBio element type"
            )
        dims.append(_dim(block.type))
        if len(block.data):
            max_dim = max(max_dim, dims[-1])
    bases = [0]
    for block in blocks:
        bases.append(bases[-1] + len(block.data))

    faces = None
    kinds = ["elements"] * len(blocks)
    dropped_vertices = 0
    for b, block in enumerate(blocks):
        if dims[b] == 0:
            kinds[b] = "dropped"
            dropped_vertices += len(block.data)
        elif dims[b] == 1 and max_dim > 1:
            # Lines along the edges of other cells are an <Edge>; two-node lines
            # between nodes no cell connects are springs, a <DiscreteSet>.
            on_edges = _all_on_edges(mesh, b, dims)
            kinds[b] = "discrete" if block.type == "line" and not on_edges else "edge"
        elif dims[b] == 2 and max_dim == 3 and len(block.data):
            if faces is None:
                faces = FacetIndex(mesh, surface_edges=False)
            corners = _corners(block.type)
            if all(
                faces.find(row[:corners]) is not None
                for row in np.asarray(block.data).tolist()
            ):
                kinds[b] = "surface"

    regions = sorted(
        mesh.regions, key=lambda r: (_KIND_ORDER[r.kind], r.name, r.dim, r.tag)
    )
    names = [""] * len(blocks)
    names_block = [False] * len(regions)
    taken = set()
    for b in range(len(blocks)):
        for r, reg in enumerate(regions):
            if names[b]:
                break
            entries = np.asarray(reg.entries).ravel()
            if (
                reg.kind != "cell"
                or names_block[r]
                or len(entries) == 0
                or len(entries) != bases[b + 1] - bases[b]
            ):
                continue
            if (
                entries[0] == bases[b]
                and entries[-1] == bases[b + 1] - 1
                and reg.name not in taken
            ):
                names[b] = reg.name
                names_block[r] = True
        if not names[b]:
            stem = {
                "surface": "Surface",
                "edge": "Edge",
                "discrete": "DiscreteSet",
            }.get(kinds[b], "Part")
            names[b] = f"{stem}{b + 1}"
            while names[b] in taken:
                names[b] += "_"
        taken.add(names[b])

    element_no = np.zeros(bases[-1], dtype=np.int64)
    n_elements = 0
    for b in range(len(blocks)):
        if kinds[b] == "elements":
            n = bases[b + 1] - bases[b]
            element_no[bases[b] : bases[b + 1]] = np.arange(
                n_elements + 1, n_elements + n + 1
            )
            n_elements += n

    if dropped_vertices:
        warn(
            f"FEBio .feb writer: {dropped_vertices} vertex cell(s) have no FEBio "
            "element and were dropped"
        )
        _provenance.note(
            "cells-dropped",
            f"{dropped_vertices} vertex cell(s) have no FEBio element",
        )
    if mesh.point_data or mesh.cell_data or mesh.field_data:
        warn(
            "FEBio .feb writer: data arrays are not written (MeshData is not "
            "supported yet)"
        )
        _provenance.note("data-dropped", "the .feb writer does not write MeshData")

    out = [
        '<?xml version="1.0" encoding="ISO-8859-1"?>\n<febio_spec version="4.0">\n',
        _provenance.render_xml_comment(_provenance.SlotTier.BLOCK),
        '\n\t<Module type="solid"/>\n',
        "\t<Material>\n\t\t<!-- Placeholder materials: FEBio needs one per domain. "
        "Replace them. -->\n",
    ]
    mat_id = 0
    for b in range(len(blocks)):
        if kinds[b] != "elements":
            continue
        mat_id += 1
        out.append(
            f'\t\t<material id="{mat_id}" name="{_escape(names[b])}" '
            'type="isotropic elastic">\n\t\t\t<E>1</E>\n\t\t\t<v>0.3</v>\n'
            "\t\t</material>\n"
        )
    out.append("\t</Material>\n\t<Mesh>\n\t\t<Nodes>\n")
    points = np.asarray(mesh.points, dtype=np.float64)
    pdim = points.shape[1] if points.ndim == 2 else 0
    for p, row in enumerate(points.tolist()):
        xyz = ",".join("%.17g" % (row[d] if d < pdim else 0.0) for d in range(3))
        out.append(f'\t\t\t<node id="{p + 1}">{xyz}</node>\n')
    out.append("\t\t</Nodes>\n")

    for b, block in enumerate(blocks):
        if kinds[b] == "dropped":
            continue
        if kinds[b] == "discrete":
            out.append(f'\t\t<DiscreteSet name="{_escape(names[b])}">\n')
            for a, z in np.asarray(block.data)[:, :2].tolist():
                out.append(f"\t\t\t<delem>{a + 1},{z + 1}</delem>\n")
            out.append("\t\t</DiscreteSet>\n")
            continue
        ftype = _WRITE[block.type]
        order = node_order("febio", block.type)
        elements = kinds[b] == "elements"
        outer = {"elements": "Elements", "surface": "Surface", "edge": "Edge"}[kinds[b]]
        attrs = f' type="{ftype}"' if elements else ""
        out.append(f'\t\t<{outer}{attrs} name="{_escape(names[b])}">\n')
        inner = "elem" if elements else ftype
        for r, row in enumerate(np.asarray(block.data).tolist()):
            g = bases[b] + r
            ident = element_no[g] if elements else r + 1
            nodes = [row[k] for k in order.from_meshio] if order else row
            out.append(
                f'\t\t\t<{inner} id="{ident}">{_ids(v + 1 for v in nodes)}</{inner}>\n'
            )
        out.append(f"\t\t</{outer}>\n")

    skipped = 0
    for r, reg in enumerate(regions):
        if names_block[r]:
            continue
        entries = np.asarray(reg.entries)
        if reg.kind == "point":
            if len(entries) == 0:
                skipped += 1
                continue
            out.append(
                f'\t\t<NodeSet name="{_escape(reg.name)}">'
                f"{_ids(v + 1 for v in entries.ravel())}</NodeSet>\n"
            )
        elif reg.kind == "side":
            body = []
            edge = False
            for cell, facet in entries.reshape(-1, 2).tolist():
                found = facet_nodes(mesh, cell, facet)
                if found is None:
                    continue
                ftype, nodes = found
                tag = _WRITE[ftype]
                edge = _dim(ftype) == 1
                body.append(
                    f'\t\t\t<{tag} id="{len(body) + 1}">'
                    f"{_ids(v + 1 for v in nodes)}</{tag}>\n"
                )
            if not body:
                skipped += 1
                continue
            outer = "Edge" if edge else "Surface"
            out.append(
                f'\t\t<{outer} name="{_escape(reg.name)}">\n'
                + "".join(body)
                + f"\t\t</{outer}>\n"
            )
        else:
            ids = [
                int(element_no[g])
                for g in entries.ravel().tolist()
                if 0 <= g < bases[-1] and element_no[g]
            ]
            if not ids:
                skipped += 1
                continue
            out.append(
                f'\t\t<ElementSet name="{_escape(reg.name)}">{_ids(ids)}</ElementSet>\n'
            )
    if skipped:
        warn(
            f"FEBio .feb writer: {skipped} region(s) with nothing FEBio can hold "
            "were dropped"
        )

    out.append("\t</Mesh>\n\t<MeshDomains>\n")
    for b in range(len(blocks)):
        if kinds[b] != "elements":
            continue
        domain = {3: "SolidDomain", 2: "ShellDomain"}.get(dims[b], "BeamDomain")
        name = _escape(names[b])
        out.append(f'\t\t<{domain} name="{name}" mat="{name}"/>\n')
    out.append("\t</MeshDomains>\n</febio_spec>\n")
    with open(filename, "wb") as fh:
        fh.write("".join(out).encode())
