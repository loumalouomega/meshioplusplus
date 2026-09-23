"""I/O for Code_Aster's native ``.mail`` mesh (``LIRE_MAILLAGE(FORMAT='ASTER')``).

The pure-Python twin of ``src/cpp/src/formats/code_aster.cpp``: both engines read
the same meshes and write the same bytes.

The file is a sequence of blocks, each opened by a keyword and closed by
``FINSF``; it ends at ``FIN``. Tokens are separated by blanks or commas, ``%``
starts a comment, and Code_Aster reads only the first 80 columns of a line, so
this reader does too. Records are a token stream: an element's nodes may wrap
onto the next lines.

- ``COOR_1D``/``COOR_2D``/``COOR_3D`` hold named nodes.
- ``POI1`` … ``HEXA27`` are element blocks, one cell block each, with the node
  order of the ``"code_aster"`` tables in :mod:`meshioplusplus._node_order`
  (not MED's).
- ``GROUP_MA``/``GROUP_NO`` become cell/point regions, named by ``NOM=`` or by the
  block's first token.
"""

import re

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import open_file
from .._mesh import Mesh, topological_dimension
from .._node_order import node_order
from .._regions import Region

__all__ = ["read", "write"]

_COLUMNS = 80
_MAX_ENTITIES = 9999999
_MAX_GROUP_NAME = 24

# keyword -> (meshio++ cell type, node count)
_TYPES = {
    "POI1": ("vertex", 1),
    "SEG2": ("line", 2),
    "SEG3": ("line3", 3),
    "SEG4": ("line4", 4),
    "TRIA3": ("triangle", 3),
    "TRIA6": ("triangle6", 6),
    "TRIA7": ("triangle7", 7),
    "QUAD4": ("quad", 4),
    "QUAD8": ("quad8", 8),
    "QUAD9": ("quad9", 9),
    "TETRA4": ("tetra", 4),
    "TETRA10": ("tetra10", 10),
    "PENTA6": ("wedge", 6),
    "PENTA15": ("wedge15", 15),
    "PENTA18": ("wedge18", 18),
    "PYRAM5": ("pyramid", 5),
    "PYRAM13": ("pyramid13", 13),
    "HEXA8": ("hexahedron", 8),
    "HEXA20": ("hexahedron20", 20),
    "HEXA27": ("hexahedron27", 27),
}
_KEYWORD_OF = {cell_type: kw for kw, (cell_type, _) in _TYPES.items()}

# Keywords Code_Aster's reader (lrmast.F90) skips without a word.
_SKIPPED = {
    "TITRE",
    "DUMP",
    "DEBUG",
    "GROUP_FA",
    "SYS_UNIT",
    "SYS_COOR",
    "MACRO_AR",
    "MACRO_FA",
    "MACRO_EL",
    "MATERIAU",
}

_KEYWORD = re.compile(r"[A-Z][A-Z0-9_]*\Z")
_NUMBER = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?\Z")
_SPLIT = re.compile(r"[ \t,]+")


def _fail(what, line):
    raise ReadError(f"Code_Aster .mail: {what} (line {line})")


def _text(raw):
    """A token's bytes (read as latin-1, one character per byte) as text."""
    return raw.encode("latin-1").decode("utf-8", "surrogateescape")


def _tokenize(data):
    tokens = []
    warned_long = False
    for line_no, line in enumerate(data.split("\n"), start=1):
        if line.endswith("\r"):
            line = line[:-1]
        if len(line) > _COLUMNS:
            comment = line.find("%")
            beyond = line[_COLUMNS:].strip(" \t") != "" and (
                comment == -1 or comment >= _COLUMNS
            )
            if beyond and not warned_long:
                warn(
                    f"Code_Aster .mail: line {line_no} is longer than 80 columns; like "
                    "Code_Aster, the reader ignores everything past column 80"
                )
                warned_long = True
            line = line[:_COLUMNS]
        comment = line.find("%")
        if comment != -1:
            line = line[:comment]
        words = [w for w in _SPLIT.split(line) if w]
        i = 0
        while i < len(words):
            w = words[i]
            while i + 1 < len(words) and (
                w.endswith("=") or words[i + 1].startswith("=")
            ):
                i += 1
                w += words[i]
            tokens.append((w, line_no))
            i += 1
    return tokens


def _number(token, line):
    if not _NUMBER.match(token):
        _fail(f"expected a coordinate, found '{_text(token)}'", line)
    return float(token.replace("D", "E").replace("d", "E"))


def read(filename):
    with open_file(filename, "rb") as f:
        data = f.read().decode("latin-1")
    tokens = _tokenize(data)
    n = len(tokens)

    point_dim = 0
    coords = []
    node_names = []
    blocks = []  # (keyword, [names], [node names], [lines])
    groups = []  # (name, is_cell, [(member, line)])
    warned_keywords = set()
    saw_fin = False
    i = 0

    def skip_options():
        nonlocal i
        name = None
        while i < n and "=" in tokens[i][0]:
            key, _, value = tokens[i][0].partition("=")
            if key.upper() in ("NOM", "NAME"):
                name = value
            i += 1
        return name

    def require(what, line):
        if i >= n:
            _fail(f"the file ends inside a block, while reading {what}", line)

    def at_block_end():
        return tokens[i][0].upper() == "FINSF"

    while i < n:
        head, head_line = tokens[i]
        kw = head.upper()
        i += 1
        if kw == "FIN" or kw.startswith("FIN("):
            saw_fin = True
            break
        if kw == "FINSF":
            continue
        if kw in ("COOR_1D", "COOR_2D", "COOR_3D"):
            dim = int(kw[5])
            if point_dim not in (0, dim):
                _fail(f"{kw} after a COOR_{point_dim}D block", head_line)
            point_dim = dim
            skip_options()
            while True:
                require("coordinates", head_line)
                if at_block_end():
                    i += 1
                    break
                if "=" in tokens[i][0]:
                    i += 1
                    continue
                node_names.append(tokens[i][0])
                line = tokens[i][1]
                i += 1
                for _ in range(dim):
                    require("coordinates", line)
                    coords.append(_number(*tokens[i]))
                    i += 1
            continue
        if kw in _TYPES:
            count = _TYPES[kw][1]
            names, nodes, lines = [], [], []
            skip_options()
            while True:
                require("elements", head_line)
                if at_block_end():
                    i += 1
                    break
                if "=" in tokens[i][0]:
                    i += 1
                    continue
                names.append(tokens[i][0])
                line = tokens[i][1]
                lines.append(line)
                i += 1
                for _ in range(count):
                    require("element nodes", line)
                    if at_block_end():
                        _fail(
                            f"{kw} element '{_text(names[-1])}' has fewer than "
                            f"{count} nodes",
                            line,
                        )
                    nodes.append(tokens[i][0])
                    i += 1
            if names:
                blocks.append((kw, names, nodes, lines))
            continue
        if kw in ("GROUP_NO", "GROUP_MA"):
            name = skip_options() or ""
            members = []
            while True:
                require("a group", head_line)
                if at_block_end():
                    i += 1
                    break
                if "=" in tokens[i][0]:
                    i += 1
                    continue
                if not name:
                    name = tokens[i][0]
                else:
                    members.append(tokens[i])
                i += 1
            if not name:
                _fail(f"{kw} block without a name", head_line)
            groups.append((name, kw == "GROUP_MA", members))
            continue
        if _KEYWORD.match(kw):
            if kw not in _SKIPPED and kw not in warned_keywords:
                warn(
                    f"Code_Aster .mail: skipping the unsupported '{_text(head)}' block "
                    f"(line {head_line})"
                )
                warned_keywords.add(kw)
            while i < n and tokens[i][0].upper() != "FINSF":
                i += 1
            if i == n:
                _fail(f"block '{_text(head)}' has no FINSF", head_line)
            i += 1
            continue
        _fail(f"expected a keyword, found '{_text(head)}'", head_line)
    if not saw_fin:
        warn(f"Code_Aster .mail: '{filename}' has no FIN line; reading it to the end")

    node_index = {}
    for p, name in enumerate(node_names):
        if name in node_index:
            raise ReadError(f"Code_Aster .mail: node '{_text(name)}' is defined twice")
        node_index[name] = p
    pdim = point_dim or 3
    points = np.array(coords, dtype=np.float64).reshape(len(node_names), pdim)

    cells = []
    element_index = {}
    element_dim = []
    for kw, names, nodes, lines in blocks:
        cell_type, k = _TYPES[kw]
        order = node_order("code_aster", cell_type)
        to_meshio = order.to_meshio if order else range(k)
        dim = topological_dimension[cell_type]
        conn = np.empty((len(names), k), dtype=np.int64)
        for r, ename in enumerate(names):
            row = nodes[r * k : (r + 1) * k]
            for j, src in enumerate(to_meshio):
                idx = node_index.get(row[src])
                if idx is None:
                    _fail(
                        f"element '{_text(ename)}' names undefined node '{_text(row[src])}'",
                        lines[r],
                    )
                conn[r, j] = idx
            if ename in element_index:
                _fail(f"element '{_text(ename)}' is defined twice", lines[r])
            element_index[ename] = len(element_dim)
            element_dim.append(dim)
        cells.append((cell_type, conn))

    members = {}
    dims = {}
    for name, is_cell, group_members in groups:
        key = (is_cell, name)
        label = "GROUP_MA" if is_cell else "GROUP_NO"
        if key in members:
            warn(
                f"Code_Aster .mail: {label} '{_text(name)}' is defined twice; merging the two"
            )
        else:
            members[key] = []
            dims[key] = -1
        ids = members[key]
        seen = set(ids)
        index = element_index if is_cell else node_index
        warned_duplicate = False
        for member, line in group_members:
            found = index.get(member)
            if found is None:
                kind = "element" if is_cell else "node"
                _fail(
                    f"{label} '{_text(name)}' names undefined {kind} '{_text(member)}'",
                    line,
                )
            if found in seen:
                if not warned_duplicate:
                    warn(
                        f"Code_Aster .mail: {label} '{_text(name)}' lists "
                        f"'{_text(member)}' more than once"
                    )
                warned_duplicate = True
                continue
            seen.add(found)
            ids.append(found)
            if is_cell:
                dims[key] = max(dims[key], element_dim[found])

    mesh = Mesh(points, cells)
    mesh.regions = [
        Region(
            _text(name),
            "cell" if is_cell else "point",
            np.array(ids, dtype=np.int64),
            dims[(is_cell, name)] if is_cell else -1,
            -1,
        )
        for (is_cell, name), ids in members.items()
    ]
    return mesh


class _Record:
    """A record being written: tokens wrap to an indented continuation line
    rather than crossing column 80."""

    def __init__(self, out):
        self.out = out
        self.column = 0

    def add(self, token):
        if self.column == 0:
            self.out.append(token)
            self.column = len(token)
            return
        if self.column + 1 + len(token) > _COLUMNS:
            self.out.append("\n        ")
            self.column = 8
        self.out.append(" " + token)
        self.column += 1 + len(token)

    def end(self):
        self.out.append("\n")
        self.column = 0


def _group_name(name, taken, kind):
    """Letters, digits and ``_``, at most 24 characters, unique within its kind."""
    raw = name.encode("utf-8", "surrogateescape")
    clean = "".join(
        chr(b) if (chr(b).isascii() and (chr(b).isalnum() or b == ord("_"))) else "_"
        for b in raw
    )
    clean = (clean or "GROUP")[:_MAX_GROUP_NAME]
    out = clean
    k = 1
    while out in taken:
        suffix = f"_{k}"
        out = clean[: min(len(clean), _MAX_GROUP_NAME - len(suffix))] + suffix
        k += 1
    taken.add(out)
    if out != name:
        warn(f"Code_Aster .mail: {kind} '{name}' is written as '{out}'")
    return out


def write(filename, mesh):
    points = np.asarray(mesh.points)
    pdim = points.shape[1] if points.ndim == 2 else 0
    if pdim > 3:
        raise WriteError(
            f"Code_Aster .mail writer: points of dimension {pdim} (at most 3)"
        )
    num_cells = 0
    for block in mesh.cells:
        if block.type not in _KEYWORD_OF or not isinstance(block.data, np.ndarray):
            raise WriteError(
                f"Code_Aster .mail writer: unsupported cell type '{block.type}'"
            )
        num_cells += len(block.data)
    if len(points) > _MAX_ENTITIES or num_cells > _MAX_ENTITIES:
        raise WriteError(
            "Code_Aster .mail writer: more than 9,999,999 nodes or elements do not fit "
            "8-character names; write MED instead"
        )

    regions = sorted(getattr(mesh, "regions", []) or [], key=lambda r: r.key)
    side_regions = sum(1 for r in regions if r.kind == "side")
    if side_regions:
        warn(
            f"Code_Aster .mail has no facet groups; {side_regions} side region(s) dropped"
        )
        _provenance.note(
            "regions-dropped",
            f"{side_regions} side region(s) have no Code_Aster .mail group",
        )
    if mesh.point_data or mesh.cell_data or mesh.field_data:
        warn(
            "Code_Aster .mail holds no data arrays; point, cell and field data dropped"
        )
        _provenance.note("data-dropped", "a Code_Aster .mail mesh holds no data arrays")

    out = [_provenance.render_lines(_provenance.SlotTier.BLOCK, "% ")]
    rec = _Record(out)
    out.append(f"COOR_{pdim or 3}D\n")
    for p, xyz in enumerate(points.astype(np.float64, copy=False)):
        rec.add(f"N{p + 1}")
        for v in xyz:
            rec.add("%.16E" % v)
        rec.end()
    out.append("FINSF\n")

    label = 0
    for block in mesh.cells:
        keyword = _KEYWORD_OF[block.type]
        order = node_order("code_aster", block.type)
        data = np.asarray(block.data)
        if order is not None:
            data = data[:, list(order.from_meshio)]
        out.append(keyword + "\n")
        for row in data:
            label += 1
            rec.add(f"M{label}")
            for node in row:
                rec.add(f"N{int(node) + 1}")
            rec.end()
        out.append("FINSF\n")

    taken_no, taken_ma = set(), set()
    for region in regions:
        if region.kind == "side":
            continue
        cells = region.kind == "cell"
        kind = "GROUP_MA" if cells else "GROUP_NO"
        name = _group_name(region.name, taken_ma if cells else taken_no, kind)
        out.append(f"{kind} NOM={name}\n")
        prefix = "M" if cells else "N"
        for e in np.asarray(region.entries, dtype=np.int64).ravel():
            rec.add(f"{prefix}{int(e) + 1}")
        if rec.column:
            rec.end()
        out.append("FINSF\n")
    out.append("FIN\n")

    text = "".join(out)
    with open_file(filename, "w", newline="\n") as f:
        f.write(text)
