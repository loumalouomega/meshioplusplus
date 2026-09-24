"""
I/O for LS-DYNA keyword input decks (``.k`` / ``.key`` / ``.dyn``).

The reader collects every ``*NODE``, ``*ELEMENT_*``, ``*PART`` and ``*SET_*`` card of
the deck and of the files it ``*INCLUDE``s into one state, and resolves ids only once
everything is read, because a set may refer to elements defined in a later file.
``*PART`` becomes a cell region (title as name, ``pid`` as tag); sets become point,
cell and side regions. The C++ twin is ``formats/lsdyna.cpp``.
"""

import pathlib

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._facets import FacetIndex
from .._files import open_file
from .._mesh import Mesh
from .._regions import Region
from .._skin import _CELL_FACES
from . import _cards
from ._cards import I10, IDS10, LONG, STD, format_real16, split_card, to_float, to_int

_MAX_INCLUDE_DEPTH = 8
_MAX_STD_ID = 99_999_999

# Element family -> topological dimension.
_FAMILY_DIM = {"solid": 3, "tshell": 3, "shell": 2, "beam": 1, "discrete": 1, "mass": 0}
_SET_FAMILY = {
    "NODE": "node",
    "SOLID": "solid",
    "SHELL": "shell",
    "TSHELL": "tshell",
    "BEAM": "beam",
    "DISCRETE": "discrete",
    "PART": "part",
    "SEGMENT": "segment",
}
_SHELL_OPTIONS = {"THICKNESS", "BETA", "MCID", "OFFSET"}
_PART_KEYWORDS = {"PART", "PART_CONTACT", "PART_COMPOSITE", "PART_INERTIA"}
_INCLUDE_KEYWORDS = {"INCLUDE", "INCLUDE_NO_TRANSFORM", "INCLUDE_TRANSFORM"}


def _collapse_solid(n):
    """(cell type, nodes in meshio++ order) of an 8-node LS-DYNA solid.

    LS-DYNA has no tetra, pyramid or wedge card: they are hexahedra with repeated
    nodes. The patterns are checked most-degenerate first. The twin of
    ``detail::collapse_brick``, which the Radioss reader shares.
    """
    if n[3] == n[4] == n[5] == n[6] == n[7]:
        return "tetra", [n[0], n[1], n[2], n[3]]
    if n[2] == n[3] and n[4] == n[5] == n[6] == n[7]:
        return "tetra", [n[0], n[1], n[2], n[4]]
    if n[4] == n[5] == n[6] == n[7]:
        return "pyramid", [n[0], n[1], n[2], n[3], n[4]]
    if n[2] == n[3] and n[6] == n[7]:
        return "wedge", [n[0], n[1], n[2], n[4], n[5], n[6]]
    if n[4] == n[5] and n[6] == n[7]:
        return "wedge", [n[0], n[4], n[1], n[3], n[6], n[2]]
    # One side edge collapsed in both the bottom and the top face (Radioss
    # writes 1 2 3 1 5 6 7 5): the wedge keeps the faces' cyclic order.
    for i in range(4):
        j, k, m = (i + 1) % 4, (i + 2) % 4, (i + 3) % 4
        if (
            n[i] == n[j]
            and n[i + 4] == n[j + 4]
            and n[j] != n[k]
            and n[k] != n[m]
            and n[m] != n[j]
        ):
            return "wedge", [n[j], n[k], n[m], n[j + 4], n[k + 4], n[m + 4]]
    return "hexahedron", list(n)


def _expand_solid(ctype, row):
    """The 8 LS-DYNA nodes of a meshio++ tetra, pyramid, wedge or hexahedron."""
    r = [int(v) for v in row]
    if ctype == "tetra":
        return [r[0], r[1], r[2], r[3], r[3], r[3], r[3], r[3]]
    if ctype == "pyramid":
        return [r[0], r[1], r[2], r[3], r[4], r[4], r[4], r[4]]
    if ctype == "wedge":
        return [r[0], r[2], r[5], r[3], r[1], r[1], r[4], r[4]]
    return r


class _Deck:
    def __init__(self):
        self.node_index = {}
        self.coords = []
        self.elements = []  # (family, eid, pid, cell type, [file node ids], group)
        self.group = 0  # one per element keyword, so blocks follow the deck's sections
        self.parts = {}  # pid -> title, in first-seen order
        self.sets = []  # (family, sid, title, [ids] | [(a, b) segments])
        self.include_dirs = []
        self.warned = set()
        self.param_skips = 0

    def warn_once(self, key, message):
        if key not in self.warned:
            self.warned.add(key)
            warn(message)

    # -- files -----------------------------------------------------------------
    def read_file(self, path, depth, mode):
        if depth > _MAX_INCLUDE_DEPTH:
            raise ReadError(
                f"LS-DYNA: *INCLUDE nested deeper than {_MAX_INCLUDE_DEPTH}"
            )
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                text = f.read()
        except OSError as exc:
            raise ReadError(f"LS-DYNA: could not read {path}: {exc}") from None
        self.read_text(
            text, pathlib.Path(path).resolve().parent, str(path), depth, mode
        )

    def read_text(self, text, base_dir, label, depth, mode):
        lines = [ln[:-1] if ln.endswith("\r") else ln for ln in text.split("\n")]
        pos = 0
        n = len(lines)
        while pos < n:
            line = lines[pos]
            pos += 1
            if not line.startswith("*"):
                if line.startswith("-----BEGIN PGP"):
                    pos = self.skip_pgp(lines, pos)
                continue
            keyword, card_mode, rest = _parse_keyword(line)
            block = []
            while pos < n and not lines[pos].startswith("*"):
                raw = lines[pos]
                pos += 1
                if raw.startswith("$"):
                    continue
                if raw.startswith("-----BEGIN PGP"):
                    pos = self.skip_pgp(lines, pos)
                    continue
                block.append((pos, raw))
            if keyword == "END":
                return
            if keyword == "KEYWORD":
                mode = _keyword_mode(rest, mode)
                continue
            m = card_mode or mode
            ctx = (label, m)
            if keyword == "NODE":
                self.read_nodes(block, ctx)
            elif keyword.startswith("ELEMENT_"):
                self.read_elements(keyword, block, ctx)
            elif keyword in _PART_KEYWORDS:
                self.read_parts(keyword, block, ctx)
            elif keyword.startswith("SET_"):
                self.read_set(keyword, block, ctx)
            elif keyword in _INCLUDE_KEYWORDS:
                self.read_includes(keyword, block, base_dir, depth, mode)
            elif keyword in ("INCLUDE_PATH", "INCLUDE_PATH_RELATIVE"):
                for _, raw in block:
                    d = raw.strip()
                    if d:
                        self.include_dirs.append((d, base_dir))

    def skip_pgp(self, lines, pos):
        self.warn_once("pgp", "LS-DYNA: skipped a PGP-encrypted block")
        while pos < len(lines) and not lines[pos].startswith("-----END PGP"):
            pos += 1
        return pos + 1

    def read_includes(self, keyword, block, base_dir, depth, mode):
        names = []
        cur = None
        for _, raw in block:
            s = raw.strip()
            if not s:
                continue
            if cur is not None:
                cur += s
            else:
                cur = s
            if cur.endswith("+"):
                cur = cur[:-1]
                continue
            names.append(cur)
            cur = None
            if keyword == "INCLUDE_TRANSFORM":
                self.warn_once(
                    "transform",
                    "LS-DYNA: *INCLUDE_TRANSFORM offsets are not applied; "
                    "the file is included untransformed",
                )
                break
        if cur:
            names.append(cur)
        for name in names:
            path = self.find_include(name, base_dir)
            if path is None:
                warn(f"LS-DYNA: include file {name!r} not found; skipped")
                continue
            self.read_file(path, depth + 1, mode)

    def find_include(self, name, base_dir):
        spellings = [name]
        if "\\" in name:
            spellings.append(name.replace("\\", "/"))
        for sp in spellings:
            p = pathlib.Path(sp)
            cands = [p] if p.is_absolute() else [base_dir / p]
            if not p.is_absolute():
                for d, parent in self.include_dirs:
                    dp = pathlib.Path(d)
                    cands.append(dp / p if dp.is_absolute() else parent / dp / p)
                    if not dp.is_absolute():
                        cands.append(dp / p)
                cands.append(pathlib.Path.cwd() / p)
            for c in cands:
                if c.is_file():
                    return c
        return None

    # -- cards -----------------------------------------------------------------
    def skip_param(self, line):
        if "&" in line:
            self.param_skips += 1
            return True
        return False

    def read_nodes(self, block, ctx):
        label, mode = ctx
        for lineno, raw in block:
            if not raw.strip() or self.skip_param(raw):
                continue
            where = f" (line {lineno} of {label})"
            f = split_card(raw, _cards.NODE, mode)
            nid = to_int(f[0], where)
            if nid in self.node_index:
                raise ReadError(f"LS-DYNA: duplicate node id {nid}{where}")
            self.node_index[nid] = len(self.coords)
            self.coords.append(
                (to_float(f[1], where), to_float(f[2], where), to_float(f[3], where))
            )

    def read_elements(self, keyword, block, ctx):
        label, mode = ctx
        tokens = keyword.split("_")[1:]  # after ELEMENT
        kind = tokens[0]
        opts = tokens[1:]
        extras = 0
        if kind == "SOLID" and not opts:
            pass
        elif kind == "SOLID" and opts == ["ORTHO"]:
            extras = 2
        elif kind == "SOLID" and opts[:1] in (["TET4TOTET10"], ["H8TOH20"]):
            self.warn_once(
                keyword,
                f"LS-DYNA: *{keyword} is a conversion directive and is not applied",
            )
            return
        elif kind == "SHELL" and set(opts) <= _SHELL_OPTIONS:
            extras = len(opts)
        elif kind in ("TSHELL", "BEAM", "DISCRETE", "MASS") and not opts:
            pass
        elif kind in ("SOLID", "SHELL", "TSHELL", "BEAM", "DISCRETE", "MASS"):
            self.warn_once(keyword, f"LS-DYNA: *{keyword} is not supported; skipped")
            return
        else:
            return
        family = kind.lower()
        self.group += 1
        group = self.group

        j = 0
        while j < len(block):
            lineno, raw = block[j]
            j += 1
            if not raw.strip():
                continue
            where = f" (line {lineno} of {label})"
            if self.skip_param(raw):
                j += extras
                continue
            if family == "mass":
                f = split_card(raw, _cards.ELEMENT_MASS, mode)
                self.elements.append(
                    (
                        "mass",
                        to_int(f[0], where),
                        to_int(f[3], where),
                        "vertex",
                        [to_int(f[1], where)],
                        group,
                    )
                )
                continue
            f = split_card(raw, _cards.ELEMENT, mode)
            eid = to_int(f[0], where)
            pid = to_int(f[1], where)
            if family in ("solid", "tshell"):
                if family == "solid" and not any(f[2:10]):
                    # Two-line layout: "eid pid" then up to ten node fields.
                    if j >= len(block):
                        raise ReadError(f"LS-DYNA: truncated element card{where}")
                    lineno2, raw2 = block[j]
                    j += 1
                    g = split_card(raw2, _cards.ELEMENT, mode)
                    if len(g) < 10:
                        g = g + [""] * (10 - len(g))
                    nodes = [to_int(v, where) for v in g[:10]]
                    while nodes and nodes[-1] == 0:
                        nodes.pop()
                    if len(nodes) == 10:
                        ctype = "tetra10"
                    elif len(nodes) == 8:
                        ctype, nodes = _collapse_solid(nodes)
                    elif len(nodes) == 4:
                        ctype = "tetra"
                    else:
                        raise ReadError(
                            f"LS-DYNA: solid element with {len(nodes)} nodes is not "
                            f"supported{where}"
                        )
                else:
                    nodes = [to_int(v, where) for v in f[2:10]]
                    ctype, nodes = _collapse_solid(nodes)
            elif family == "shell":
                n = [to_int(v, where) for v in f[2:10]]
                if any(n[4:]):
                    self.warn_once(
                        "shell-mid",
                        "LS-DYNA: higher-order shell nodes (n5..n8) are ignored",
                    )
                if n[3] == 0 or n[3] == n[2]:
                    ctype, nodes = "triangle", n[:3]
                else:
                    ctype, nodes = "quad", n[:4]
            else:  # beam, discrete
                ctype = "line"
                nodes = [to_int(f[2], where), to_int(f[3], where)]
                if family == "beam":
                    while j < len(block) and "." in block[j][1]:
                        j += 1  # an optional second card; element cards have no "."
            self.elements.append((family, eid, pid, ctype, nodes, group))
            j += extras

    def read_parts(self, keyword, block, ctx):
        label, mode = ctx
        lines = list(block)
        j = 0
        while j + 1 < len(lines):
            title = lines[j][1].strip()
            lineno, raw = lines[j + 1]
            if not raw.strip():
                break
            where = f" (line {lineno} of {label})"
            f = split_card(raw, IDS10, mode)
            pid = to_int(f[0], where)
            self.parts[pid] = title
            j += 2
            if keyword != "PART":
                break

    def read_set(self, keyword, block, ctx):
        label, mode = ctx
        tokens = keyword.split("_")[1:]
        fam_key = tokens[0]
        opts = set(tokens[1:])
        family = _SET_FAMILY.get(fam_key)
        if family is None:
            return
        if not opts <= {"LIST", "GENERATE", "TITLE"}:
            self.warn_once(keyword, f"LS-DYNA: *{keyword} is not supported; skipped")
            return
        lines = list(block)
        title = ""
        j = 0
        if "TITLE" in opts:
            if not lines:
                return
            title = lines[0][1].strip()
            j = 1
        if j >= len(lines):
            return
        lineno, raw = lines[j]
        j += 1
        where = f" (line {lineno} of {label})"
        sid = to_int(split_card(raw, IDS10, mode)[0], where)
        ids = []
        for lineno, raw in lines[j:]:
            if not raw.strip() or self.skip_param(raw):
                continue
            where = f" (line {lineno} of {label})"
            if family == "segment":
                f = split_card(raw, _cards.SEGMENT, mode)
                seg = [to_int(v, where) for v in f[:4]]
                if seg[3] == 0:
                    seg[3] = seg[2]
                ids.append(tuple(seg))
                continue
            f = [to_int(v, where) for v in split_card(raw, IDS10, mode)]
            if "GENERATE" in opts:
                for k in range(0, len(f) - 1, 2):
                    if f[k] > 0 and f[k + 1] >= f[k]:
                        ids.extend(range(f[k], f[k + 1] + 1))
            else:
                ids.extend(v for v in f if v != 0)
        self.sets.append((family, sid, title, ids))


def _parse_keyword(line):
    body = line[1:].strip()
    end = len(body)
    for k, ch in enumerate(body):
        if ch in " \t,":
            end = k
            break
    name = body[:end]
    rest = body[end:].strip(" \t,")
    up = name.upper()
    mode = None
    if up and up[-1] in "+-%":
        mode = {"+": LONG, "-": STD, "%": I10}[up[-1]]
        up = up[:-1]
    elif rest in ("+", "-", "%"):
        mode = {"+": LONG, "-": STD, "%": I10}[rest]
        rest = ""
    return up, mode, rest


def _keyword_mode(rest, mode):
    """The deck-wide default a ``*KEYWORD`` card sets (``LONG=Y|S|K``, ``I10=Y``)."""
    for token in rest.replace(",", " ").split():
        key, _, value = token.partition("=")
        key = key.upper()
        value = value.upper()
        if key == "LONG" and value in ("Y", "S", "K"):
            mode = LONG
        elif key == "I10" and value == "Y" and mode != LONG:
            mode = I10
    return mode


# -- reading ---------------------------------------------------------------------


def read(filename):
    """Read an LS-DYNA keyword deck."""
    deck = _Deck()
    if hasattr(filename, "read"):
        text = filename.read()
        if isinstance(text, bytes):
            text = text.decode("utf-8", "replace")
        name = getattr(filename, "name", None)
        base = pathlib.Path(name).resolve().parent if isinstance(name, str) else None
        deck.read_text(
            text, base or pathlib.Path.cwd(), str(name or "<buffer>"), 0, STD
        )
    else:
        deck.read_file(filename, 0, STD)
    return _build_mesh(deck)


def _build_mesh(deck):
    if deck.param_skips:
        warn(
            f"LS-DYNA: {deck.param_skips} mesh card(s) use *PARAMETER values (&name) "
            "and were skipped"
        )
    points = np.array(deck.coords, dtype=np.float64).reshape(-1, 3)

    # One block per (element keyword, cell type), so the blocks follow the deck's
    # sections; a keyword whose first type is the previous block's extends it.
    blocks = []  # [cell type, [element indices]]
    block_of = {}
    for e, (_, _, _, ctype, _, group) in enumerate(deck.elements):
        b = block_of.get((group, ctype))
        if b is None:
            if blocks and blocks[-1][0] == ctype:
                b = len(blocks) - 1
            else:
                blocks.append([ctype, []])
                b = len(blocks) - 1
            block_of[(group, ctype)] = b
        blocks[b][1].append(e)
    owner = {}  # (family, eid) -> global cell index
    cells = []
    cell_pid = []
    cell_family = []
    for ctype, members in blocks:
        conn = np.empty(
            (len(members), len(deck.elements[members[0]][4])), dtype=np.int64
        )
        for r, e in enumerate(members):
            family, eid, pid, _, nodes, _ = deck.elements[e]
            for c, nid in enumerate(nodes):
                idx = deck.node_index.get(nid)
                if idx is None:
                    raise ReadError(
                        f"LS-DYNA: element {eid} references undefined node {nid}"
                    )
                conn[r, c] = idx
            key = (family, eid)
            if key in owner:
                raise ReadError(f"LS-DYNA: duplicate {family} element id {eid}")
            owner[key] = len(cell_pid)
            cell_pid.append(pid)
            cell_family.append(family)
        cells.append((ctype, conn))

    regions = []
    part_cells = {}
    by_pid = {}
    for g, p in enumerate(cell_pid):
        by_pid.setdefault(p, []).append(g)
    for pid, title in deck.parts.items():
        members = by_pid.get(pid, [])
        dim = max((_FAMILY_DIM[cell_family[g]] for g in members), default=-1)
        part_cells[pid] = members
        regions.append(
            Region(title or f"Part {pid}", "cell", members, dim=dim, tag=pid)
        )

    mesh = Mesh(points, cells)
    dropped = 0
    face_map = None
    seen = {(r.kind, r.name, r.dim, r.tag) for r in regions}
    for family, sid, title, ids in deck.sets:
        name = title or f"{family.upper()} set {sid}"
        if family == "node":
            entries = [deck.node_index[i] for i in ids if i in deck.node_index]
            dropped += len(ids) - len(entries)
            kind = "point"
        elif family == "part":
            entries = [g for pid in ids for g in part_cells.get(pid, [])]
            kind = "cell"
        elif family == "segment":
            if face_map is None:
                face_map = _face_map(mesh)
            entries = []
            for seg in ids:
                idx = tuple(deck.node_index.get(i, -1) for i in seg)
                key = _segment_key(idx)
                hit = None if key is None else face_map.find(key)
                if hit is None:
                    dropped += 1
                else:
                    entries.append(hit.first)
            kind = "side"
        else:
            entries = [owner[(family, i)] for i in ids if (family, i) in owner]
            dropped += len(ids) - len(entries)
            kind = "cell"
        if (kind, name, -1, sid) in seen:
            name = f"{name} [{family}]"
            if (kind, name, -1, sid) in seen:
                continue
        seen.add((kind, name, -1, sid))
        if kind == "side":
            arr = np.array(entries, dtype=np.int64).reshape(-1, 2)
        else:
            arr = np.array(entries, dtype=np.int64)
        regions.append(Region(name, kind, arr, dim=-1, tag=sid))
    if dropped:
        warn(f"LS-DYNA: {dropped} set entries refer to undefined ids and were dropped")
    mesh.regions = regions
    return mesh


def _segment_key(nodes):
    """Corner nodes of a 3- or 4-node segment as a sorted tuple (``-1`` = undefined)."""
    if -1 in nodes:
        return None
    n = list(nodes)
    if n[3] == n[2]:
        n = n[:3]
    return tuple(sorted(n))


def _face_map(mesh):
    """Corner-node lookup over every face of every solid and every shell's own face
    (as facet 0). An interior face is shared by two cells; the lowest cell index wins.
    """
    return FacetIndex(mesh, surface_edges=False, surface_self=True)


# -- writing ---------------------------------------------------------------------

_SOLID_TYPES = {"tetra", "pyramid", "wedge", "hexahedron", "tetra10"}
_TYPE_FAMILY = {
    "vertex": "mass",
    "line": "beam",
    "triangle": "shell",
    "quad": "shell",
    "tetra": "solid",
    "pyramid": "solid",
    "wedge": "solid",
    "hexahedron": "solid",
    "tetra10": "solid",
}
_KEYWORD_OF_FAMILY = {
    "solid": "ELEMENT_SOLID",
    "shell": "ELEMENT_SHELL",
    "beam": "ELEMENT_BEAM",
    "mass": "ELEMENT_MASS",
}


def _title_line(name):
    t = str(name).replace("\r", " ").replace("\n", " ")
    return " " + t if t.startswith(("*", "$")) else t


def _assign_parts(regions, cells, block_of_cell):
    """(pid of each cell, [(pid, title)]) from the mesh's cell regions.

    A cell region with a topological dimension is a part when none of its cells is
    already claimed; every other cell region is written as a set. Cells that no part
    claims form one part per cell block.
    """
    ncells = len(block_of_cell)
    pid_of = [0] * ncells
    parts = []
    used = set()
    as_sets = []
    claimed = bytearray(ncells)
    candidates = [r for r in regions if r.kind == "cell" and r.dim >= 0]
    for r in candidates:
        if len(r.entries) and any(claimed[int(g)] for g in r.entries):
            as_sets.append(r)
            continue
        for g in r.entries:
            claimed[int(g)] = 1
        want = r.tag if r.tag > 0 and r.tag not in used else 0
        parts.append([want, r.name, r])
        if want:
            used.add(want)
    nxt = max(used, default=0)
    for entry in parts:
        if not entry[0]:
            nxt += 1
            entry[0] = nxt
        for g in entry[2].entries:
            pid_of[int(g)] = entry[0]
    result = [(pid, name) for pid, name, _ in parts]
    for b, cb in enumerate(cells):
        idx = [g for g in range(ncells) if block_of_cell[g] == b and not pid_of[g]]
        if not idx:
            continue
        nxt += 1
        for g in idx:
            pid_of[g] = nxt
        result.append((nxt, cb.type))
    return pid_of, result, as_sets


def write(filename, mesh: Mesh) -> None:
    points = np.asarray(mesh.points, dtype=np.float64)
    npts = len(points)
    if npts > _MAX_STD_ID:
        raise WriteError("LS-DYNA writer: too many nodes for 8-column ids")
    block_of_cell = []
    for b, cb in enumerate(mesh.cells):
        if cb.type not in _TYPE_FAMILY:
            raise WriteError(f"LS-DYNA writer: unsupported cell type {cb.type!r}")
        block_of_cell.extend([b] * len(cb.data))
    if len(block_of_cell) > _MAX_STD_ID:
        raise WriteError("LS-DYNA writer: too many elements for 8-column ids")
    regions = sorted(mesh.regions, key=lambda r: r.key)
    pid_of, parts, as_sets = _assign_parts(regions, mesh.cells, block_of_cell)

    with open_file(filename, "wt") as f:
        f.write("*KEYWORD\n")
        f.write(_provenance.render_lines(_provenance.SlotTier.BLOCK, "$ "))
        f.write("*NODE\n")
        for k in range(npts):
            x = points[k]
            xyz = [format_real16(x[c]) if c < len(x) else "0.0" for c in range(3)]
            f.write(f"{k + 1:>8}" + "".join(f"{v:>16}" for v in xyz) + "\n")
        g = 0
        for cb in mesh.cells:
            family = _TYPE_FAMILY[cb.type]
            two_line = cb.type == "tetra10"
            f.write(f"*{_KEYWORD_OF_FAMILY[family]}\n")
            for row in np.asarray(cb.data):
                eid = g + 1
                pid = pid_of[g]
                g += 1
                if family == "mass":
                    f.write(f"{eid:>8}{int(row[0]) + 1:>8}{'0.0':>16}{pid:>8}\n")
                    continue
                if family == "solid":
                    nodes = list(row) if two_line else _expand_solid(cb.type, row)
                elif cb.type == "triangle":
                    nodes = [row[0], row[1], row[2], row[2]]
                else:
                    nodes = list(row)
                ids = "".join(f"{int(v) + 1:>8}" for v in nodes)
                if two_line:
                    f.write(f"{eid:>8}{pid:>8}\n{ids}\n")
                else:
                    f.write(f"{eid:>8}{pid:>8}{ids}\n")
        for pid, name in parts:
            f.write(f"*PART\n{_title_line(name)}\n{pid:>10}{1:>10}{1:>10}\n")
        _write_sets(f, mesh, regions, as_sets)
        f.write("*END\n")


def _cell_families(mesh):
    out = []
    for cb in mesh.cells:
        out.extend([_TYPE_FAMILY[cb.type]] * len(cb.data))
    return out


def _write_ids(f, ids):
    for k in range(0, len(ids), 8):
        f.write("".join(f"{v:>10}" for v in ids[k : k + 8]) + "\n")


def _write_sets(f, mesh, regions, cell_regions):
    families = _cell_families(mesh)
    used = {}

    def sid_for(family, tag):
        taken = used.setdefault(family, set())
        if tag > 0 and tag not in taken:
            taken.add(tag)
            return tag
        new = max(taken, default=0) + 1
        while new in taken:
            new += 1
        taken.add(new)
        return new

    node_regions = [r for r in regions if r.kind == "point"]
    side_regions = [r for r in regions if r.kind == "side"]
    plain = [r for r in regions if r.kind == "cell" and r.dim < 0]
    for r in node_regions:
        sid = sid_for("node", r.tag)
        f.write(f"*SET_NODE_LIST_TITLE\n{_title_line(r.name)}\n{sid:>10}\n")
        _write_ids(f, [int(v) + 1 for v in r.entries])
    for r in plain + [c for c in cell_regions if c.dim >= 0]:
        by_family = {}
        for g in r.entries:
            by_family.setdefault(families[int(g)], []).append(int(g) + 1)
        if not by_family:
            by_family["solid"] = []
        if len(by_family) > 1:
            warn(f"LS-DYNA writer: region {r.name!r} spans several element families")
        for family in ("solid", "shell", "beam", "mass"):
            if family not in by_family or family == "mass":
                continue
            sid = sid_for(family, r.tag)
            f.write(
                f"*SET_{family.upper()}_LIST_TITLE\n{_title_line(r.name)}\n{sid:>10}\n"
            )
            _write_ids(f, by_family[family])
    for r in side_regions:
        sid = sid_for("segment", r.tag)
        f.write(f"*SET_SEGMENT_TITLE\n{_title_line(r.name)}\n{sid:>10}\n")
        for g, facet in r.entries:
            nodes = _segment_nodes(mesh, int(g), int(facet))
            if nodes is None:
                continue
            f.write("".join(f"{v + 1:>10}" for v in nodes) + "\n")


def _segment_nodes(mesh, g, facet):
    base = 0
    for block in mesh.cells:
        if g < base + len(block.data):
            row = [int(v) for v in block.data[g - base]]
            faces = _CELL_FACES.get(block.type)
            if faces is not None:
                if not 0 <= facet < len(faces):
                    return None
                _, ncorner, local = faces[facet]
                nodes = [row[i] for i in local[:ncorner]]
            elif block.type in ("triangle", "quad"):
                nodes = row
            else:
                return None
            return nodes + [nodes[-1]] if len(nodes) == 3 else nodes
        base += len(block.data)
    return None
