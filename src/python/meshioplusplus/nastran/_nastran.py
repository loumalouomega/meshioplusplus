"""
I/O for Nastran and OptiStruct bulk data.

The C++ core (``nastran.cpp``) is the twin of this reader and writer; both
give the same Mesh and write the same HyperMesh comment block.
"""

from __future__ import annotations

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import open_file
from .._mesh import CellBlock, Mesh, topological_dimension
from .._regions import Region

# Element card -> (meshio type, node fields). A fixed count reads exactly that
# many fields and ignores the rest of the card (THETA, ZOFFS, thicknesses, an
# orientation vector ...); 0 marks a solid whose linear or quadratic form is
# told apart by the number of node fields given.
_ELEMENTS = {
    "CELAS1": ("vertex", 1),
    "CBEAM": ("line", 2),
    "CBUSH": ("line", 2),
    "CBUSH1D": ("line", 2),
    "CROD": ("line", 2),
    "CGAP": ("line", 2),
    "CBAR": ("line", 2),
    "CTRIAR": ("triangle", 3),
    "CTRIA3": ("triangle", 3),
    "CTRAX6": ("triangle6", 6),
    "CTRIAX6": ("triangle6", 6),
    "CTRIA6": ("triangle6", 6),
    "CQUADR": ("quad", 4),
    "CSHEAR": ("quad", 4),
    "CQUAD4": ("quad", 4),
    "CQUAD8": ("quad8", 8),
    "CQUAD9": ("quad9", 9),
    "CTETRA": ("tetra", 0),
    "CPYRAM": ("pyramid", 0),
    "CPYRA": ("pyramid", 0),
    "CPENTA": ("wedge", 0),
    "CHEXA": ("hexahedron", 0),
}
# linear type -> (linear nodes, quadratic type, quadratic nodes)
_SOLIDS = {
    "tetra": (4, "tetra10", 10),
    "pyramid": (5, "pyramid13", 13),
    "wedge": (6, "wedge15", 15),
    "hexahedron": (8, "hexahedron20", 20),
}
nastran_to_meshio_type = {k: v[0] for k, v in _ELEMENTS.items()}
meshio_to_nastran_type = {
    "vertex": "CELAS1",
    "line": "CBAR",
    "triangle": "CTRIA3",
    "triangle6": "CTRIA6",
    "quad": "CQUAD4",
    "quad8": "CQUAD8",
    "quad9": "CQUAD9",
    "tetra": "CTETRA",
    "tetra10": "CTETRA",
    "pyramid": "CPYRA",
    "pyramid13": "CPYRA",
    "wedge": "CPENTA",
    "wedge15": "CPENTA",
    "hexahedron": "CHEXA",
    "hexahedron20": "CHEXA",
}

# meshio slot j holds Nastran slot P[j]. hexahedron20 and wedge15 put the
# vertical mid-edges before the top ring in Nastran (both are involutions);
# CTRIAX6/CTRAX6 list corner, mid, corner, mid, corner, mid.
_PERM_HEX20 = list(range(12)) + [16, 17, 18, 19, 12, 13, 14, 15]
_PERM_WEDGE15 = list(range(9)) + [12, 13, 14, 9, 10, 11]
_PERM_TRIAX6 = [0, 2, 4, 1, 3, 5]
_WRITE_PERM = {"hexahedron20": _PERM_HEX20, "wedge15": _PERM_WEDGE15}

# Cards read silently although they carry nothing meshio++ keeps: properties,
# materials, loads, constraints, coordinate systems, tables and solution
# parameters. Anything else not read is counted and named in one warning.
_QUIET_PREFIXES = (
    "P",
    "MAT",
    "SPC",
    "MPC",
    "FORCE",
    "MOMENT",
    "LOAD",
    "TEMP",
    "GRAV",
    "RFORCE",
    "ACCEL",
    "CORD",
    "TABLE",
    "EIGR",
    "EIGC",
    "NLPARM",
    "TSTEP",
    "FREQ",
    "SUPORT",
    "DAREA",
    "DLOAD",
    "RLOAD",
    "TLOAD",
    "SPOINT",
    "ASET",
    "OMIT",
    "INCLUDE",
)

_BLANK8 = " " * 8


def read(filename):
    with open_file(filename, "r") as f:
        out = read_buffer(f)
    return out


def _is_comment(line):
    return len(line) < 3 or line.startswith(("$", "//", "#"))


def _int(text, card):
    if text == "":
        return 0
    try:
        return int(text)
    except ValueError:
        raise ReadError(f"Nastran: invalid integer field '{text}' in a {card} card")


def _real(text, card):
    if text == "":
        return 0.0
    try:
        return _nastran_string_to_float(text)
    except (ValueError, IndexError):
        raise ReadError(f"Nastran: invalid real field '{text}' in a {card} card")


def _parse_id_list(tokens, card):
    """Explicit ids and ``a THRU b`` ranges of an id list."""
    t = [x for x in tokens if x != ""]
    ids, ranges = [], []
    k = 0
    while k < len(t):
        if k + 2 < len(t) and t[k + 1] == "THRU":
            ranges.append((_int(t[k], card), _int(t[k + 2], card)))
            k += 3
        else:
            ids.append(_int(t[k], card))
            k += 1
    return ids, ranges


def _quoted_after(line, pos):
    """``(id, name)`` from ``<spaces><digits>...\"name\"`` at ``pos``, or None."""
    k = pos
    while k < len(line) and line[k] == " ":
        k += 1
    digits = k
    while k < len(line) and "0" <= line[k] <= "9":
        k += 1
    if k == digits:
        return None
    opening = line.find('"', k)
    if opening < 0:
        return None
    closing = line.find('"', opening + 1)
    if closing < 0:
        return None
    return int(line[digits:k]), line[opening + 1 : closing], closing + 1


class _HyperMesh:
    """The HyperMesh comment cards: component membership ($HMMOVE plus the
    ``$`` id lines after it), component names ($HMNAME COMP) and SET names
    ($HMSET)."""

    def __init__(self):
        self.components = {}
        self.component_names = {}
        self.set_names = {}
        # component id -> the property id `$HMNAME COMP <id>"name" <pid> "type"`
        # records after the name, when it does
        self.component_properties = {}
        self.active = None

    @staticmethod
    def _id_line(line):
        if not line.startswith("$") or line[1:8].strip(" ") != "":
            return None
        tokens = []
        for k in range(8, len(line), 8):
            f = line[k : k + 8].strip(" \t")
            if f == "":
                continue
            if f != "THRU" and not f.isdigit():
                return None
            tokens.append(f)
        if not any(t != "THRU" for t in tokens):
            return None
        return tokens

    def feed(self, line):
        if line.startswith("$HMMOVE"):
            ident = line[7:16].strip(" \t")
            self.active = int(ident) if ident.isdigit() else None
            if self.active is not None:
                self.components.setdefault(self.active, ([], []))
            return
        if self.active is not None:
            tokens = self._id_line(line)
            if tokens is not None:
                ids, ranges = _parse_id_list(tokens, "$HMMOVE")
                self.components[self.active][0].extend(ids)
                self.components[self.active][1].extend(ranges)
                return
        self.active = None
        if line.startswith("$HMNAME COMP "):
            q = _quoted_after(line, 13)
            if q is not None:
                self.component_names[q[0]] = q[1]
                k = q[2]
                while k < len(line) and line[k] == " ":
                    k += 1
                e = k
                while e < len(line) and "0" <= line[e] <= "9":
                    e += 1
                if e > k:
                    self.component_properties[q[0]] = int(line[k:e])
        elif line.startswith("$HMSET "):
            k = 7
            while k < len(line) and line[k] == " ":
                k += 1
            e = k
            while e < len(line) and "0" <= line[e] <= "9":
                e += 1
            if e > k:
                q = _quoted_after(line, e)
                if q is not None:
                    self.set_names[int(line[k:e])] = q[1]


def _chunk_line(line):
    """Free-field lines split on commas; fixed-field lines into (at most ten)
    8-column fields, the tenth being the continuation marker."""
    if "," in line:
        return line.split(",")
    return [line[i : i + 8] for i in range(0, len(line), 8)][:10]


def _merge_large(c):
    """Re-merge each pair of 8-column chunks of a large-field line."""
    d = [c[0]]
    for k in (1, 3, 5, 7):
        if k >= len(c):
            break
        f = c[k]
        if k + 1 < len(c) and c[k + 1] is not None:
            f += c[k + 1]
        d.append(f)
    if len(c) > 9:
        d.append(c[9])
    return d


def _cards(lines):
    """The logical cards: each the flattened, stripped list of its fields,
    continuation lines merged in."""
    i = 0
    while i < len(lines):
        chunks = [_chunk_line(lines[i])]
        free = "," in lines[i]
        i += 1
        while i < len(lines):
            nxt = lines[i]
            if nxt[0] in "+*":
                if len(chunks[-1]) == 10:
                    chunks[-1][-1] = None
                c = _chunk_line(nxt)
                c[0] = None
                chunks.append(c)
                i += 1
            elif len(chunks[-1]) == 10 and chunks[-1][-1] == _BLANK8:
                # implicit continuation: a blank tenth field followed by a
                # line whose first field is blank too
                c = _chunk_line(nxt)
                if c and c[0] == _BLANK8:
                    chunks[-1][9] = None
                    c[0] = None
                    chunks.append(c)
                    i += 1
                else:
                    break
            else:
                break
        head = chunks[0][0].strip(" \t")
        if not free and head.endswith("*"):
            chunks = [_merge_large(c) for c in chunks]
        yield [f.strip(" \t") for c in chunks for f in c if f is not None]


def _field(fields, k):
    return fields[k] if k < len(fields) else ""


def _resolve(group, index, sorted_ids, missing):
    ids, ranges = group
    out = []
    for i in ids:
        if i in index:
            out.append(index[i])
        else:
            missing[0] += 1
    for a, b in ranges:
        lo = np.searchsorted(sorted_ids, a, side="left")
        hi = np.searchsorted(sorted_ids, b, side="right")
        out.extend(index[int(i)] for i in sorted_ids[lo:hi])
    return np.unique(np.array(out, dtype=np.int64))


def read_buffer(f):
    # Everything before BEGIN BULK (executive and case control, I/O options)
    # is skipped; comment lines feed the HyperMesh parser; ENDDATA ends it.
    hm = _HyperMesh()
    lines = []
    bulk = False
    for line in f:
        line = line.rstrip("\n").rstrip("\r")
        if not bulk:
            bulk = line.strip(" \t").startswith("BEGIN BULK")
            continue
        if line.startswith("ENDDATA"):
            break
        if line.startswith("$"):
            hm.feed(line)
        else:
            hm.active = None
        if not _is_comment(line):
            lines.append(line)
    if not bulk:
        raise ReadError('Nastran: "BEGIN BULK" statement not found')

    points = []
    points_id = []
    point_refs = []
    any_point_ref = False
    blocks = []  # [type, conn rows, refs, ids]
    any_cell_ref = False
    cell_index = {}
    cell_dims = []
    sets = {}
    skipped = {}

    for fields in _cards(lines):
        kw = fields[0] if fields else ""
        if kw.endswith("*"):
            kw = kw[:-1]
        if kw == "GRID":
            pid = _int(_field(fields, 1), kw)
            ref = _field(fields, 2)
            any_point_ref = any_point_ref or ref != ""
            point_refs.append(_int(ref, kw))
            points_id.append(pid)
            points.append([_real(_field(fields, c), kw) for c in (3, 4, 5)])
            continue
        if kw in _ELEMENTS:
            eid = _int(_field(fields, 1), kw)
            ref = _field(fields, 2)
            cell_type, count = _ELEMENTS[kw]
            if count > 0:
                nodes = []
                for j in range(count):
                    t = _field(fields, 3 + j)
                    if t == "":
                        raise ReadError(f"Nastran: {kw} {eid} is missing node {j + 1}")
                    nodes.append(_int(t, kw))
            else:
                nodes = [_int(t, kw) for t in fields[3:] if t != ""]
                lin, quad_type, quad = _SOLIDS[cell_type]
                if len(nodes) == quad:
                    cell_type = quad_type
                elif len(nodes) != lin:
                    raise ReadError(
                        f"Nastran: {kw} {eid} has {len(nodes)} nodes; "
                        f"expected {lin} or {quad}"
                    )
            if kw in ("CTRIAX6", "CTRAX6"):
                nodes = [nodes[i] for i in _PERM_TRIAX6]
            elif cell_type in _WRITE_PERM:
                nodes = [nodes[i] for i in _WRITE_PERM[cell_type]]
            if not blocks or blocks[-1][0] != cell_type:
                blocks.append([cell_type, [], [], []])
            blocks[-1][1].append(nodes)
            any_cell_ref = any_cell_ref or ref != ""
            blocks[-1][2].append(_int(ref, kw))
            blocks[-1][3].append(eid)
            cell_index[eid] = len(cell_dims)
            cell_dims.append(topological_dimension[cell_type])
            continue
        if kw == "SET":
            # OptiStruct: SET, id, GRID|ELEM, LIST, ids (with THRU ranges)
            sid = _int(_field(fields, 1), kw)
            kind, sub = _field(fields, 2), _field(fields, 3)
            if kind not in ("GRID", "ELEM") or sub != "LIST":
                warn(f"Nastran: SET {sid} of type '{kind} {sub}' is not read; skipped")
                continue
            ids, ranges = _parse_id_list(fields[4:], kw)
            entry = sets.setdefault(sid, [kind, [], []])
            entry[0] = kind
            entry[1].extend(ids)
            entry[2].extend(ranges)
            continue
        if not kw.startswith(_QUIET_PREFIXES):
            skipped[kw] = skipped.get(kw, 0) + 1

    if skipped:
        listed = ", ".join(f"{k} ({skipped[k]})" for k in sorted(skipped))
        warn(
            f"Nastran: skipped {sum(skipped.values())} card(s) "
            f"meshio++ does not read: {listed}"
        )

    point_index = {pid: k for k, pid in enumerate(points_id)}
    points = np.array(points, dtype=np.float64).reshape(-1, 3)
    cells = []
    cells_id = []
    cell_refs = []
    for cell_type, conn, refs, ids in blocks:
        try:
            data = np.array(
                [[point_index[n] for n in row] for row in conn], dtype=np.int64
            )
        except KeyError as exc:
            raise ReadError(
                f"Nastran: an element references grid {exc.args[0]}, "
                "which is not defined"
            )
        cells.append(CellBlock(cell_type, data))
        cell_refs.append(np.array(refs, dtype=np.int64))
        cells_id.append(np.array(ids, dtype=np.int64))

    mesh = Mesh(points, cells)
    mesh.points_id = np.array(points_id, dtype=np.int64)
    mesh.cells_id = cells_id
    if any_point_ref:
        mesh.point_data["nastran:ref"] = np.array(point_refs, dtype=np.int64)
    if any_cell_ref:
        mesh.cell_data["nastran:ref"] = cell_refs

    # Regions: HyperMesh components (tag = component id), then SET cards
    # (tag = set id). Names resolve last: $HMNAME may follow the elements.
    sorted_cells = np.array(sorted(cell_index), dtype=np.int64)
    sorted_points = np.array(sorted(point_index), dtype=np.int64)
    missing = [0]

    def cell_dim(ids):
        dims = {cell_dims[i] for i in ids}
        return dims.pop() if len(dims) == 1 else -1

    regions = []
    for comp in hm.component_names:
        hm.components.setdefault(comp, ([], []))
    members = {
        comp: _resolve(hm.components[comp], cell_index, sorted_cells, missing)
        for comp in sorted(hm.components)
    }
    # An element no $HMMOVE lists belongs to the component whose recorded
    # property is its PID (HyperMesh writes $HMMOVE only for the others).
    if hm.component_properties:
        moved = np.zeros(len(cell_dims), dtype=bool)
        for ids in members.values():
            moved[ids] = True
        by_property = {}
        for comp, prop in sorted(hm.component_properties.items()):
            by_property.setdefault(prop, comp)
        pids = np.concatenate(cell_refs) if cell_refs else np.zeros(0, np.int64)
        for prop, comp in by_property.items():
            extra = np.flatnonzero((pids == prop) & ~moved)
            if len(extra):
                members[comp] = np.union1d(members[comp], extra).astype(np.int64)
    for comp, ids in members.items():
        name = hm.component_names.get(comp, f"component_{comp}")
        regions.append(Region(name, "cell", ids, cell_dim(ids), comp))
    for sid in sorted(sets):
        kind, ids, ranges = sets[sid]
        name = hm.set_names.get(sid, f"set_{sid}")
        if kind == "GRID":
            got = _resolve((ids, ranges), point_index, sorted_points, missing)
            regions.append(Region(name, "point", got, -1, sid))
        else:
            got = _resolve((ids, ranges), cell_index, sorted_cells, missing)
            regions.append(Region(name, "cell", got, cell_dim(got), sid))
    if missing[0]:
        warn(
            f"Nastran: {missing[0]} region member id(s) name no grid or element; "
            "dropped"
        )
    if regions:
        mesh.regions = regions
    return mesh


# There are two basic categories of input data formats in NX Nastran:
#
# - "Free" format data, in which the data fields are simply separated by
#   commas. This type of data is known as free field data.
#
# - "Fixed" format data, in which your data must be aligned in columns of
#   specific width. There are two subcategories of fixed format data that differ
#   based on the size of the fixed column width:
#
#     - Small field format, in which a single line of data is divided into 10
#       fields that can contain eight characters each.
#
#     - Large field format, in which a single line of input is expanded into
#       two lines The first and last fields on each line are eight columns wide,
#       while the intermediate fields are sixteen columns wide. The large field
#       format is useful when you need greater numerical accuracy.
#
# See: https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/User.pdf


def _hypermesh_block(mesh):
    """Disjoint cell regions as HyperMesh components: ``$HMMOVE`` with the
    element ids (runs as ``a THRU b``), then ``$HMNAME COMP``. Comments only,
    so no solver sees them. Byte-identical to the C++ writer's block."""
    num_cells = sum(len(c) for c in mesh.cells)
    owned = np.zeros(num_cells, dtype=bool)
    kept = []
    # the C++ mesh's (kind, name, dim, tag) order, so both engines keep the
    # same region of an overlapping pair
    kind_order = {"point": 0, "cell": 1, "side": 2}
    for reg in sorted(
        mesh.regions, key=lambda r: (kind_order[r.kind], r.name, r.dim, r.tag)
    ):
        if reg.kind != "cell":
            warn(
                f"Nastran: {reg.kind} region '{reg.name}' dropped; "
                "HyperMesh components hold cells only"
            )
            continue
        e = np.asarray(reg.entries, dtype=np.int64)
        bad = e[(e < 0) | (e >= num_cells)]
        if len(bad):
            raise WriteError(
                f"Nastran writer: cell region '{reg.name}' names cell {bad[0]} "
                f"of {num_cells}"
            )
        if owned[e].any():
            warn(
                f"Nastran: cell region '{reg.name}' overlaps an earlier one and is "
                "dropped; a HyperMesh component owns each element once"
            )
            continue
        owned[e] = True
        kept.append(reg)
    if not kept:
        return ""
    tags = [int(r.tag) for r in kept]
    tags_ok = all(t > 0 for t in tags) and len(set(tags)) == len(tags)
    out = []
    for k, reg in enumerate(kept):
        comp = tags[k] if tags_ok else k + 1
        name = reg.name.replace('"', "'")
        out.append("$\n")
        e = [int(x) for x in np.asarray(reg.entries, dtype=np.int64)]
        if e:
            out.append(f"$HMMOVE {comp:8d}\n")
            singles = []

            def flush():
                if singles:
                    out.append("$       " + "".join(singles) + "\n")
                    singles.clear()

            j = 0
            while j < len(e):
                run = j
                while run + 1 < len(e) and e[run + 1] == e[run] + 1:
                    run += 1
                if run > j:
                    flush()
                    out.append(f"$       {e[j] + 1:8d}THRU    {e[run] + 1:8d}\n")
                else:
                    singles.append(f"{e[j] + 1:8d}")
                    if len(singles) == 8:
                        flush()
                j = run + 1
            flush()
        out.append(f'$HMNAME COMP{comp:20d}"{name}"\n')
    return "".join(out)


def write(filename, mesh, point_format="fixed-large", cell_format="fixed-small"):
    if point_format == "free":
        grid_fmt = "GRID,{:d},{:s},{:s},{:s},{:s}\n"
        float_fmt = _float_to_nastran_string
    elif point_format == "fixed-small":
        grid_fmt = "GRID    {:<8d}{:<8s}{:>8s}{:>8s}{:>8s}\n"
        float_fmt = _float_rstrip
    elif point_format == "fixed-large":
        grid_fmt = "GRID*   {:<16d}{:<16s}{:>16s}{:>16s}\n*       {:>16s}\n"
        float_fmt = _float_to_nastran_string
    else:
        raise WriteError(f'Nastran: unknown point format "{point_format}"')

    if cell_format == "free":
        int_fmt, cell_info_fmt = "{:d}", "{:s},{:d},{:s},"
        sjoin = ","
        nipl1, nipl2 = 6, 14
    elif cell_format == "fixed-small":
        int_fmt, cell_info_fmt = "{:<8d}", "{:<8s}{:<8d}{:<8s}"
        sjoin, cchar = "", "+"
        nipl1, nipl2 = 6, 14
    elif cell_format == "fixed-large":
        int_fmt, cell_info_fmt = "{:<16d}", "{:<8s}{:<16d}{:<16s}"
        sjoin, cchar = "", "*"
        nipl1, nipl2 = 2, 6
    else:
        raise WriteError(f'Nastran: unknown cell format "{cell_format}"')

    for cell_block in mesh.cells:
        if cell_block.type not in meshio_to_nastran_type:
            raise WriteError(f"Nastran writer: unsupported cell type {cell_block.type}")

    if mesh.points.shape[1] == 2:
        warn(
            "Nastran requires 3D points, but 2D points given. "
            "Appending 0 third component."
        )
        points = np.column_stack([mesh.points, np.zeros_like(mesh.points[:, 0])])
    else:
        points = mesh.points

    def ref_text(value):
        value = int(value)
        return str(value) if value != 0 else ""

    with open_file(filename, "w") as f:
        f.write(_provenance.render_lines(_provenance.SlotTier.BLOCK, "$ "))
        f.write("BEGIN BULK\n")

        # Points; a zero reference field is written blank
        point_refs = mesh.point_data.get("nastran:ref", None)
        for point_id, x in enumerate(points):
            fx = [float_fmt(k) for k in x]
            pref = ref_text(point_refs[point_id]) if point_refs is not None else ""
            string = grid_fmt.format(point_id + 1, pref, fx[0], fx[1], fx[2])
            f.write(string)

        # CellBlock
        cell_id = 0
        cell_refs = mesh.cell_data.get("nastran:ref", None)
        for ict, cell_block in enumerate(mesh.cells):
            cell_type = cell_block.type
            cells = cell_block.data
            nastran_type = meshio_to_nastran_type[cell_type]
            if cell_format.endswith("-large"):
                nastran_type += "*"
            perm = _WRITE_PERM.get(cell_type)
            for ic, cell in enumerate(cells):
                cell_ref = ref_text(cell_refs[ict][ic]) if cell_refs is not None else ""
                cell_id += 1
                cell_info = cell_info_fmt.format(nastran_type, cell_id, cell_ref)
                cell1 = cell + 1
                if perm is not None:
                    cell1 = cell1[perm]
                conn = sjoin.join(int_fmt.format(nid) for nid in cell1[:nipl1])

                if len(cell1) > nipl1:
                    if cell_format == "free":
                        cflag1 = cflag3 = ""
                        cflag2 = cflag4 = "+,"
                    else:
                        cflag1 = cflag2 = f"{cchar}1{cell_id:<6x}"
                        cflag3 = cflag4 = f"{cchar}2{cell_id:<6x}"
                    f.write(cell_info + conn + cflag1 + "\n")
                    conn = sjoin.join(int_fmt.format(nid) for nid in cell1[nipl1:nipl2])
                    if len(cell1) > nipl2:
                        f.write(cflag2 + conn + cflag3 + "\n")
                        conn = sjoin.join(int_fmt.format(nid) for nid in cell1[nipl2:])
                        f.write(cflag4 + conn + "\n")
                    else:
                        f.write(cflag2 + conn + "\n")
                else:
                    f.write(cell_info + conn + "\n")

        f.write(_hypermesh_block(mesh))
        f.write("ENDDATA\n")


def _float_rstrip(x, n=8):
    return f"{x:f}".rstrip("0")[:n]


def _float_to_nastran_string(value, length=16):
    """
    From
    <https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/User.pdf>:

    Real numbers, including zero, must contain a decimal point. You can enter
    real numbers in a variety of formats. For example, the following are all
    acceptable versions of the real number, seven:
    ```
    7.0   .7E1  0.7+1
    .70+1 7.E+0 70.-1
    ```

    This methods converts a float value into the corresponding string. Choose
    the variant with `E` to make the file less ambigious when edited by a
    human. (`5.-1` looks like 4.0, not 5.0e-1 = 0.5.)

    Examples:
        1234.56789 --> "1.23456789E+3"
        -0.1234 --> "-1.234E-1"
        3.1415926535897932 --> "3.14159265359E+0"
    """
    out = np.format_float_scientific(value, exp_digits=1, precision=11).replace(
        "e", "E"
    )
    assert len(out) <= 16
    return out
    # The following is the manual float conversion. Keep it around for a while in case
    # we still need it.

    # aux = length - 2
    # # sfmt = "{" + f":{length}s" + "}"
    # sfmt = "{" + ":s" + "}"
    # pv_fmt = "{" + f":{length}.{aux}e" + "}"

    # if value == 0.0:
    #     return sfmt.format("0.")

    # python_value = pv_fmt.format(value)  # -1.e-2
    # svalue, sexponent = python_value.strip().split("e")
    # exponent = int(sexponent)  # removes 0s

    # sign = "-" if abs(value) < 1.0 else "+"

    # # the exponent will be added later...
    # sexp2 = str(exponent).strip("-+")
    # value2 = float(svalue)

    # # the plus 1 is for the sign
    # len_sexp = len(sexp2) + 1
    # leftover = length - len_sexp
    # leftover = leftover - 3 if value < 0 else leftover - 2
    # fmt = "{" + f":1.{leftover:d}f" + "}"

    # svalue3 = fmt.format(value2)
    # svalue4 = svalue3.strip("0")
    # field = sfmt.format(svalue4 + sign + sexp2)
    # return field


def _nastran_string_to_float(string):
    try:
        return float(string)
    except ValueError:
        string = string.strip()
        return float(string[0] + string[1:].replace("+", "e+").replace("-", "e-"))


# NOTE: format registration now lives in meshioplusplus/nastran/__init__.py, which wraps
# the reader/writer above with the C++-backed fast paths.
