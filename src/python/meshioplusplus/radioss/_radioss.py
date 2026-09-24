"""I/O for OpenRadioss / Radioss starter decks (``*_0000.rad``).

The pure-Python twin of ``src/cpp/src/formats/radioss.cpp``: both engines read
the same meshes.

``#RADIOSS STARTER``, ``/BEGIN`` (run name, input version, units), then
``/KEYWORD/option/id`` blocks up to ``/END``; ``#``/``$`` lines are comments and
``#include file`` inlines another file. Fields are 10 (integer) and 20 (real)
columns from format 51 on, 8 and 16 before; a line with a comma is split on
commas. Nodes and elements make the mesh; parts, subsets, groups and ``/SURF/SEG``
surfaces become regions.
"""

import os

import numpy as np

from .._common import num_nodes_per_cell, warn
from .._exceptions import ReadError
from .._facets import FacetIndex
from .._mesh import Mesh, topological_dimension
from .._node_order import node_order
from .._regions import Region
from ..lsdyna._lsdyna import _collapse_solid

__all__ = ["read"]

_MAX_INCLUDE_DEPTH = 8

# keyword -> (group family, nodes, meshio++ type)
_ELEMENTS = {
    "BRICK": ("BRIC", 8, "hexahedron"),
    "PENTA6": ("BRIC", 6, "wedge"),
    "TETRA4": ("BRIC", 4, "tetra"),
    "TETRA10": ("BRIC", 10, "tetra10"),
    "BRIC20": ("BRIC", 20, "hexahedron20"),
    "SHELL": ("SHEL", 4, "quad"),
    "SH3N": ("SH3N", 3, "triangle"),
    "QUAD": ("QUAD", 4, "quad"),
    "TRIA": ("TRIA", 3, "triangle"),
    "BEAM": ("BEAM", 2, "line"),
    "TRUSS": ("TRUS", 2, "line"),
    "SPRING": ("SPRI", 2, "line"),
}

_GROUPS = {
    "GRNOD": "NODE",
    "GRBRIC": "BRIC",
    "GRSHEL": "SHEL",
    "GRSH3N": "SH3N",
    "GRQUAD": "QUAD",
    "GRTRIA": "TRIA",
    "GRBEAM": "BEAM",
    "GRTRUS": "TRUS",
    "GRSPRI": "SPRI",
}

_UNREAD = ("TSHELL", "TSH3N", "SHEL16", "SPHCEL", "RIVET", "XELEM")


def _collect(path, depth, out):
    """The deck's lines with includes inlined and comments dropped; True at /END."""
    if depth > _MAX_INCLUDE_DEPTH:
        raise ReadError(f"Radioss: #include nested deeper than {_MAX_INCLUDE_DEPTH}")
    try:
        with open(path, "rb") as f:
            text = f.read().decode("latin-1")
    except OSError:
        raise ReadError(f"Radioss: cannot open {path}")
    label = os.path.basename(path)
    for number, line in enumerate(text.split("\n"), start=1):
        if line.endswith("\r"):
            line = line[:-1]
        if line[:9].upper().startswith("#INCLUDE"):
            name = line[8:].strip(" \t\r")
            if len(name) >= 2 and name[0] in "\"'" and name[-1] == name[0]:
                name = name[1:-1]
            name = name.replace("\\", "/")
            inc = (
                name
                if os.path.isabs(name)
                else os.path.join(os.path.dirname(path), name)
            )
            if not os.path.isfile(inc):
                warn(f"Radioss: include file '{name}' not found; skipped")
                continue
            _collect(inc, depth + 1, out)
            continue
        if depth > 0 and line.strip(" \t\r").upper().startswith("#ENDDATA"):
            return False
        # Blank lines stay: a blank /PART or group title is still its title line.
        if line and line[0] in "#$":
            continue
        if line.strip(" \t\r").upper() == "/END":
            if depth == 0:
                return True
            continue
        out.append((line, f"{label}:{number}"))
    return False


def _fields(line, width, count):
    if "," in line:
        return [f.strip(" \t\r") for f in line.split(",")]
    return [line[k * width : (k + 1) * width].strip(" \t\r") for k in range(count)]


def _int(text, where):
    text = text.strip()
    if not text:
        return 0
    try:
        return int(text)
    except ValueError:
        try:
            v = float(text.replace("D", "E").replace("d", "e"))
        except ValueError:
            raise ReadError(f"Radioss: bad integer field '{text}' ({where})")
        if v != int(v):
            raise ReadError(f"Radioss: bad integer field '{text}' ({where})")
        return int(v)


def _real(text, where):
    text = text.strip()
    if not text:
        return 0.0
    t = text.replace("D", "E").replace("d", "e")
    try:
        return float(t)
    except ValueError:
        # Fortran's 1.5-3 (no exponent letter).
        k = max(t.rfind("+"), t.rfind("-"))
        if k > 0 and t[k - 1] not in "eE":
            try:
                return float(t[:k] + "E" + t[k:])
            except ValueError:
                pass
        raise ReadError(f"Radioss: bad real field '{text}' ({where})")


def read(filename):
    path = os.fspath(filename)
    try:
        with open(path, "rb") as f:
            head = f.read(512).decode("latin-1").upper()
    except OSError:
        raise ReadError(f"Radioss: cannot open {path}")
    if "#RADIOSS ENGINE" in head:
        raise ReadError(
            f"Radioss: {path} is an engine deck (_0001.rad); read the starter deck "
            "(_0000.rad)"
        )
    lines = []
    _collect(path, 0, lines)

    version = 2019
    iw, rw = 10, 20
    node_ids, coords = [], []
    elements = []  # (family, id, part, type, nodes, where)
    parts = {}
    part_order = []
    groups = []  # (keyword, subtype, title, id, ids)
    surfaces = []  # (id, title, segments)
    subsets = {}
    skipped_keywords = set()
    zero_springs = linear_bric20 = 0

    n = len(lines)
    i = 0
    while i < n:
        text, where = lines[i]
        if not text or text[0] != "/":
            i += 1
            continue
        path_ = [p.strip(" \t\r").upper() for p in text.strip(" \t\r")[1:].split("/")]
        key = path_[0]
        body = i + 1
        end = body
        while end < n and (not lines[end][0] or lines[end][0][0] != "/"):
            end += 1

        def last_id():
            return _int(path_[-1], where) if len(path_) >= 2 else 0

        if key == "BEGIN":
            if body + 1 < end:
                f = _fields(lines[body + 1][0], 10, 2)
                if f and f[0]:
                    version = _int(f[0], lines[body + 1][1])
            iw, rw = (10, 20) if version >= 51 else (8, 16)
        elif key == "NODE":
            for k in range(body, end):
                ln, w = lines[k]
                if "," in ln:
                    f = _fields(ln, iw, 4)
                else:
                    f = [ln[:iw].strip()] + [
                        ln[iw + rw * d : iw + rw * (d + 1)].strip() for d in range(3)
                    ]
                if not f or not f[0]:
                    continue
                node_ids.append(_int(f[0], w))
                coords.append(
                    [_real(f[d], w) if d < len(f) else 0.0 for d in (1, 2, 3)]
                )
        elif key in _ELEMENTS:
            family, count, ctype = _ELEMENTS[key]
            part = last_id()
            per_record = 11 if key == "TETRA10" else count + 1
            values = []
            rec_where = where
            for k in range(body, end):
                ln, w = lines[k]
                if not ln.strip(" \t\r"):
                    continue
                if key == "TETRA10":
                    take = 1 if not values else 10
                elif key == "BRIC20":
                    take = 9 if not values else (8 if len(values) == 9 else 4)
                else:
                    take = per_record
                f = _fields(ln, iw, take)
                if not values:
                    rec_where = w
                for j in range(take):
                    if len(values) < per_record:
                        values.append(_int(f[j], w) if j < len(f) and f[j] else 0)
                if len(values) < per_record and key in ("TETRA10", "BRIC20"):
                    continue
                nodes = values[1 : 1 + count]
                t = ctype
                keep = True
                if key == "BRICK":
                    t, nodes = _collapse_solid(nodes)
                elif key == "SHELL" and (nodes[3] == nodes[2] or nodes[3] == 0):
                    t, nodes = "triangle", nodes[:3]
                elif key == "BRIC20":
                    if 0 in nodes[8:]:
                        linear_bric20 += 1
                        t, nodes = "hexahedron", nodes[:8]
                    else:
                        order = node_order("radioss", "hexahedron20")
                        nodes = [nodes[j] for j in order.to_meshio]
                elif key == "SPRING" and nodes[1] == 0:
                    zero_springs += 1
                    keep = False
                if keep:
                    elements.append(
                        (
                            family,
                            values[0],
                            part,
                            t,
                            list(nodes),
                            rec_where,
                            key in ("TETRA4", "TETRA10"),
                        )
                    )
                values = []
            if values:
                raise ReadError(
                    f"Radioss: /{key} element {values[0]} is cut short ({rec_where})"
                )
        elif key == "PART":
            pid = last_id()
            title, prop, mat, subset = "", 0, 0, 0
            k = body
            if k < end:
                title = lines[k][0].strip(" \t\r")
                k += 1
            if k < end:
                f = _fields(lines[k][0], iw, 3)
                w = lines[k][1]
                prop = _int(f[0], w) if f and f[0] else 0
                mat = _int(f[1], w) if len(f) > 1 and f[1] else 0
                subset = _int(f[2], w) if len(f) > 2 and f[2] else 0
            if pid not in parts:
                parts[pid] = (title, prop, mat, subset)
                part_order.append(pid)
        elif key in _GROUPS and len(path_) >= 3:
            k = body
            title = ""
            if k < end:
                title = lines[k][0].strip(" \t\r")
                k += 1
            ids = []
            for kk in range(k, end):
                for f in _fields(lines[kk][0], iw, 10):
                    if f:
                        ids.append(_int(f, lines[kk][1]))
            groups.append((key, path_[1], title, last_id(), ids))
        elif key == "SURF" and len(path_) >= 3 and path_[1] == "SEG":
            k = body
            title = ""
            if k < end:
                title = lines[k][0].strip(" \t\r")
                k += 1
            segs = []
            for kk in range(k, end):
                f = _fields(lines[kk][0], iw, 5)
                w = lines[kk][1]
                seg = [
                    _int(f[j + 1], w) if j + 1 < len(f) and f[j + 1] else 0
                    for j in range(4)
                ]
                if seg[0] or seg[1] or seg[2]:
                    segs.append(seg)
            surfaces.append((last_id(), title, segs))
        elif key == "SURF":
            skipped_keywords.add("/SURF/" + (path_[1] if len(path_) > 1 else ""))
        elif key == "SUBSET":
            k = body
            title = ""
            if k < end:
                title = lines[k][0].strip(" \t\r")
                k += 1
            children = []
            for kk in range(k, end):
                for f in _fields(lines[kk][0], iw, 10):
                    if f:
                        children.append(_int(f, lines[kk][1]))
            subsets[last_id()] = (title, children)
        elif key.startswith("GR") or key in _UNREAD:
            skipped_keywords.add("/" + key)
        i = end
    if skipped_keywords:
        warn(f"Radioss: mesh keywords not read: {', '.join(sorted(skipped_keywords))}")
    if zero_springs:
        warn(f"Radioss: {zero_springs} /SPRING element(s) with a single node skipped")
    if linear_bric20:
        warn(
            f"Radioss: {linear_bric20} /BRIC20 element(s) with missing mid-edge nodes "
            "read as hexahedra"
        )

    node_index = {}
    for p, nid in enumerate(node_ids):
        if nid in node_index:
            raise ReadError(f"Radioss: node {nid} is defined twice")
        node_index[nid] = p
    points = np.array(coords, dtype=np.float64).reshape(-1, 3)

    # Tetrahedra come in either winding (every /TETRA4 of OpenRadioss's INT_25 QA
    # deck is inverted, gmsh writes them positive): an inverted one is mirrored.
    reoriented = 0
    for el in elements:
        if not el[6]:
            continue
        nodes = el[4]
        idx = [node_index.get(v) for v in nodes[:4]]
        if None in idx:
            continue
        p = points[idx]
        if np.linalg.det(np.array([p[1] - p[0], p[2] - p[0], p[3] - p[0]])) >= 0.0:
            continue
        nodes[1], nodes[2] = nodes[2], nodes[1]
        if len(nodes) == 10:  # edges 0-1 <-> 0-2 and 1-3 <-> 2-3
            nodes[4], nodes[6] = nodes[6], nodes[4]
            nodes[8], nodes[9] = nodes[9], nodes[8]
        reoriented += 1
    if reoriented:
        warn(f"Radioss: {reoriented} inverted tetrahedra reoriented")

    order_of_types = []
    by_type = {}
    for e, el in enumerate(elements):
        if el[3] not in by_type:
            by_type[el[3]] = []
            order_of_types.append(el[3])
        by_type[el[3]].append(e)
    owner = {}
    part_cells = {}
    cell_dim = []
    cell_family = []
    cell_conn = []
    cells = []
    part_blocks, prop_blocks, mat_blocks = [], [], []
    for t in order_of_types:
        members = by_type[t]
        k = num_nodes_per_cell[t]
        conn = np.empty((len(members), k), dtype=np.int64)
        pa = np.empty(len(members), dtype=np.int64)
        pr = np.empty(len(members), dtype=np.int64)
        ma = np.empty(len(members), dtype=np.int64)
        for r, e in enumerate(members):
            family, eid, part, _, nodes, w, _ = elements[e]
            for j in range(k):
                idx = node_index.get(nodes[j])
                if idx is None:
                    raise ReadError(
                        f"Radioss: element {eid} names undefined node {nodes[j]} ({w})"
                    )
                conn[r, j] = idx
            cell = len(cell_dim)
            fam = owner.setdefault(family, {})
            if eid in fam:
                raise ReadError(f"Radioss: element {eid} is defined twice ({w})")
            fam[eid] = cell
            p = parts.get(part)
            pa[r] = part
            pr[r] = p[1] if p else 0
            ma[r] = p[2] if p else 0
            part_cells.setdefault(part, []).append(cell)
            cell_dim.append(topological_dimension[t])
            cell_family.append(family)
            cell_conn.append(conn[r].tolist())
        cells.append((t, conn))
        part_blocks.append(pa)
        prop_blocks.append(pr)
        mat_blocks.append(ma)

    mesh = Mesh(points, cells)
    mesh.field_data["radioss:version"] = np.array(version, dtype=np.int64)
    if cells:
        mesh.cell_data["radioss:part"] = part_blocks
        mesh.cell_data["radioss:property"] = prop_blocks
        mesh.cell_data["radioss:material"] = mat_blocks

    regions = []
    seen = set()

    def add_region(name, kind, tag, keyword, entries):
        dim = -1
        if kind == "cell":
            dim = max((cell_dim[c] for c in entries), default=-1)
        elif kind == "side":
            dim = max((cell_dim[c] - 1 for c, _ in entries), default=-1)
        if (kind, name, tag) in seen:
            name += f" [{keyword}]"
        seen.add((kind, name, tag))
        arr = np.array(entries, dtype=np.int64)
        if kind == "side":
            arr = arr.reshape(-1, 2)
        regions.append(Region(name, kind, arr, dim, tag))

    for pid in part_cells:
        if pid not in parts:
            parts[pid] = ("", 0, 0, 0)
            part_order.append(pid)
    for pid in part_order:
        title = parts[pid][0]
        add_region(
            title or f"Part {pid}", "cell", pid, "PART", list(part_cells.get(pid, []))
        )

    def subset_closure(root):
        ids = set()
        stack = [root]
        while stack:
            c = stack.pop()
            if c in ids:
                continue
            ids.add(c)
            if c in subsets:
                stack.extend(subsets[c][1])
        return ids

    for sid in sorted(subsets):
        title, _ = subsets[sid]
        ids = subset_closure(sid)
        members = []
        for pid in part_order:
            if parts[pid][3] in ids:
                members.extend(part_cells.get(pid, []))
        add_region(title or f"Subset {sid}", "cell", sid, "SUBSET", members)

    group_of = {(g[0], g[3]): k for k, g in enumerate(groups)}
    skipped_subtypes = set()
    dropped = 0

    def resolve(g, visiting):
        nonlocal dropped
        out = set()
        if g in visiting:
            return out
        visiting.add(g)
        keyword, subtype, _, _, ids = groups[g]
        family = _GROUPS[keyword]
        nodes = family == "NODE"
        removed = set()
        for raw in ids:
            ident = abs(raw)
            hits = set()
            if subtype == family:
                idx = (
                    node_index.get(ident) if nodes else owner.get(family, {}).get(ident)
                )
                if idx is None:
                    dropped += 1
                else:
                    hits.add(idx)
            elif subtype in ("PART", "SUBSET"):
                # The parts' nodes for GRNOD, else their cells of the group's family.
                pids = [ident]
                if subtype == "SUBSET":
                    closure = subset_closure(ident)
                    pids = [pid for pid in part_order if parts[pid][3] in closure]
                for pid in pids:
                    for c in part_cells.get(pid, []):
                        if nodes:
                            hits.update(cell_conn[c])
                        elif cell_family[c] == family:
                            hits.add(c)
            elif nodes and subtype == "SURF":
                for sid, _, segs in surfaces:
                    if sid == ident:
                        for seg in segs:
                            for v in seg:
                                if v and v in node_index:
                                    hits.add(node_index[v])
            elif nodes and _GROUPS.get(subtype, "NODE") != "NODE":
                # A GRNOD of element groups: those elements' nodes.
                sub = group_of.get((subtype, ident))
                if sub is not None:
                    for c in resolve(sub, visiting):
                        hits.update(cell_conn[c])
            elif subtype == keyword:
                sub = group_of.get((keyword, ident))
                if sub is not None:
                    hits.update(resolve(sub, visiting))
            else:
                skipped_subtypes.add(f"/{keyword}/{subtype}")
                break
            (removed if raw < 0 else out).update(hits)
        out -= removed
        visiting.discard(g)
        return out

    for g, (keyword, _, title, gid, _) in enumerate(groups):
        entries = sorted(resolve(g, set()))
        add_region(
            title or f"{keyword}_{gid}",
            "point" if _GROUPS[keyword] == "NODE" else "cell",
            gid,
            keyword,
            entries,
        )
    for t in sorted(skipped_subtypes):
        warn(
            f"Radioss: {t} groups are not resolved (only entity ids, parts and other "
            "groups are); their regions are left empty or partial"
        )

    if surfaces:
        facets = FacetIndex(mesh, surface_self=True)
        for sid, title, segs in surfaces:
            entries = []
            for seg in segs:
                # 2 nodes in a 2-D analysis, 3 for a triangle (n4 blank or n3).
                if seg[2] == 0:
                    count = 2
                else:
                    count = 3 if seg[3] == 0 or seg[3] == seg[2] else 4
                idx = [node_index.get(v) for v in seg[:count]]
                hit = facets.find(idx) if None not in idx else None
                if hit is None:
                    dropped += 1
                    continue
                entries.append(hit.first)
            add_region(title or f"SURF_{sid}", "side", sid, "SURF", entries)
    if dropped:
        warn(
            f"Radioss: {dropped} group or surface entries name undefined ids or no cell "
            "facet and were dropped"
        )
    mesh.regions = regions
    return mesh
