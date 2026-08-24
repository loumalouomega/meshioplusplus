"""
I/O for the Triangle file formats (.node/.ele/.poly), c.f.
<https://www.cs.cmu.edu/~quake/triangle.html> and
<https://people.sc.fsu.edu/~jburkardt/data/triangle_files/triangle_files.html>

Shewchuk's 2D mesh generator — the planar analogue of TetGen. A
``.node``/``.ele`` path selects the shared-stem sibling pair (a lone
``.node`` file reads as a point cloud); a ``.poly`` path reads the PSLG file
(vertices inline or in the sibling ``.node``, segments as a ``line`` cell
block; holes/regions skipped with a warning). Data naming mirrors tetgen:
``point_data["triangle:attr<k>"]`` / ``"triangle:ref"`` columns and matching
cell_data keys.

Note that ``.node``/``.ele`` default to the tetgen format; the dispatcher
falls through to this format when tetgen rejects a 2D file.
"""

import pathlib

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._mesh import CellBlock, Mesh


class _Tokens:
    """Whitespace token stream, '#' comments stripped to end of line."""

    def __init__(self, path):
        toks = []
        with open(path) as f:
            for line in f:
                line = line.split("#", 1)[0]
                toks.extend(line.split())
        self.toks = toks
        self.pos = 0

    def at_end(self):
        return self.pos >= len(self.toks)

    def next(self, what):
        if self.at_end():
            raise ReadError(f"Triangle: unexpected end of file reading {what}")
        tok = self.toks[self.pos]
        self.pos += 1
        return tok

    def next_int(self, what):
        try:
            return int(self.next(what))
        except ValueError:
            raise ReadError(f"Triangle: expected an integer for {what}")

    def next_float(self, what):
        try:
            return float(self.next(what))
        except ValueError:
            raise ReadError(f"Triangle: expected a number for {what}")


def _read_node_section(tokens):
    num_points = tokens.next_int("vertex count")
    dim = tokens.next_int("dimension")
    num_attrs = tokens.next_int("attribute count")
    num_markers = tokens.next_int("marker count")
    if dim != 2:
        raise ReadError("Triangle: need 2D points")

    ncol = 3 + num_attrs + num_markers
    rows = np.array(
        [tokens.next_float("vertex row") for _ in range(num_points * ncol)]
    ).reshape(num_points, ncol)

    base = int(rows[0, 0]) if num_points > 0 else 0
    if num_points > 0 and not np.all(rows[:, 0] == np.arange(base, base + num_points)):
        raise ReadError("Triangle: vertices not numbered consecutively")

    points = rows[:, 1:3]
    point_data = {}
    for k in range(num_attrs):
        point_data[f"triangle:attr{k + 1}"] = rows[:, 3 + k]
    for k in range(num_markers):
        flag = "" if k == 0 else str(k + 1)
        point_data["triangle:ref" + flag] = rows[:, 3 + num_attrs + k]
    return points, point_data, base


def _read_node_ele(stem):
    node_path = stem.with_suffix(".node")
    if not node_path.is_file():
        raise ReadError(f"Triangle: could not open file: {node_path}")
    points, point_data, base = _read_node_section(_Tokens(node_path))

    cells = []
    cell_data = {}
    ele_path = stem.with_suffix(".ele")
    if ele_path.is_file():  # optional: a lone .node file is a point cloud
        tokens = _Tokens(ele_path)
        ne = tokens.next_int("triangle count")
        npc = tokens.next_int("nodes per triangle")
        num_attrs = tokens.next_int("attribute count")
        if npc not in (3, 6):
            raise ReadError("Triangle: only 3- or 6-node triangles are supported")
        rows = np.array(
            [
                tokens.next_float("triangle row")
                for _ in range(ne * (1 + npc + num_attrs))
            ]
        ).reshape(ne, 1 + npc + num_attrs)
        conn = rows[:, 1 : 1 + npc].astype(np.int64) - base
        if ne > 0:  # an empty .ele adds no block
            if conn.min() < 0 or conn.max() >= len(points):
                raise ReadError("Triangle: connectivity index out of range")
            cells.append(CellBlock("triangle" if npc == 3 else "triangle6", conn))
            for k in range(num_attrs):
                flag = "" if k == 0 else str(k + 1)
                cell_data["triangle:ref" + flag] = [rows[:, 1 + npc + k]]

    return Mesh(points, cells, point_data=point_data, cell_data=cell_data)


def _read_poly(filename):
    tokens = _Tokens(filename)
    header_pos = tokens.pos
    nv = tokens.next_int("vertex count")
    if nv == 0:
        # Vertices live in the sibling .node file.
        tokens.next_int("dimension")
        tokens.next_int("attribute count")
        tokens.next_int("marker count")
        node_path = filename.with_suffix(".node")
        if not node_path.is_file():
            raise ReadError("Triangle: .poly refers to a missing sibling .node file")
        points, point_data, base = _read_node_section(_Tokens(node_path))
    else:
        tokens.pos = header_pos
        points, point_data, base = _read_node_section(tokens)

    ns = tokens.next_int("segment count")
    nmark = tokens.next_int("segment marker count")
    if nmark not in (0, 1):
        raise ReadError("Triangle: malformed segment header")
    rows = np.array(
        [tokens.next_int("segment row") for _ in range(ns * (3 + nmark))],
        dtype=np.int64,
    ).reshape(ns, 3 + nmark)
    conn = rows[:, 1:3] - base
    if ns > 0 and (conn.min() < 0 or conn.max() >= len(points)):
        raise ReadError("Triangle: segment endpoint out of range")
    cells = [CellBlock("line", conn)]
    cell_data = {}
    if nmark == 1:
        cell_data["triangle:ref"] = [rows[:, 3]]

    if not tokens.at_end():
        nh = tokens.next_int("hole count")
        if nh > 0:
            warn(f"Triangle: skipping {nh} hole(s) in {filename}")
        for _ in range(nh * 3):
            tokens.next_float("hole entry")
    if not tokens.at_end():
        nr = tokens.next_int("region count")
        if nr > 0:
            warn(f"Triangle: skipping {nr} regional attribute(s) in {filename}")

    return Mesh(points, cells, point_data=point_data, cell_data=cell_data)


def read(filename):
    filename = pathlib.Path(filename)
    if filename.suffix in (".node", ".ele"):
        return _read_node_ele(filename)
    if filename.suffix == ".poly":
        return _read_poly(filename)
    raise ReadError("Triangle: expected a .node, .ele, or .poly file")


def _fmt_value(v):
    return str(int(v)) if float(v) == int(v) else f"{v:.16e}"


def _split_point_keys(mesh):
    # Attribute columns plus at most one marker column (the first ":ref"
    # key, or else the first key), mirroring tetgen.
    attr_keys = sorted(mesh.point_data.keys())
    ref_keys = [k for k in attr_keys if ":ref" in k][:1]
    if attr_keys and not ref_keys:
        ref_keys = attr_keys[:1]
    attr_keys = [k for k in attr_keys if k not in ref_keys]
    return attr_keys, ref_keys


def _write_node_rows(fh, mesh, attr_keys, ref_keys):
    for i, pt in enumerate(mesh.points):
        row = [str(i), f"{pt[0]:.16e}", f"{pt[1]:.16e}"]
        row += [f"{mesh.point_data[k][i]:.16e}" for k in attr_keys]
        row += [_fmt_value(mesh.point_data[k][i]) for k in ref_keys]
        fh.write(" ".join(row) + "\n")


def _write_node_ele(stem, mesh):
    attr_keys, ref_keys = _split_point_keys(mesh)

    with open(stem.with_suffix(".node"), "w") as fh:
        fh.write(_provenance.render_lines(_provenance.SlotTier.BLOCK, "# "))
        fh.write(f"{len(mesh.points)} 2 {len(attr_keys)} {len(ref_keys)}\n")
        _write_node_rows(fh, mesh, attr_keys, ref_keys)

    tri_blocks = [c for c in mesh.cells if c.type in ("triangle", "triangle6")]
    tri_types = {c.type for c in tri_blocks}
    if len(tri_types) > 1:
        raise WriteError("Triangle: cannot mix triangle and triangle6 blocks")
    skipped = {c.type for c in mesh.cells} - tri_types
    if skipped:
        warn(
            f"Triangle only supports triangles; skipping {', '.join(sorted(skipped))}."
        )

    cell_attr_keys = sorted(mesh.cell_data.keys())
    ref = next((k for k in cell_attr_keys if ":ref" in k), None)
    if ref is not None:
        cell_attr_keys.remove(ref)
        cell_attr_keys.insert(0, ref)

    npc = 6 if "triangle6" in tri_types else 3
    with open(stem.with_suffix(".ele"), "w") as fh:
        fh.write(_provenance.render_lines(_provenance.SlotTier.BLOCK, "# "))
        fh.write(f"{sum(len(c) for c in tri_blocks)} {npc} {len(cell_attr_keys)}\n")
        eid = 0
        for ci, cell_block in enumerate(mesh.cells):
            if cell_block.type not in tri_types:
                continue
            for i, tri in enumerate(cell_block.data):
                row = [str(eid)] + [str(v) for v in tri]
                for k in cell_attr_keys:
                    blocks = mesh.cell_data[k]
                    row.append(_fmt_value(blocks[ci][i]) if ci < len(blocks) else "0")
                fh.write(" ".join(row) + "\n")
                eid += 1


def _write_poly(filename, mesh):
    attr_keys, ref_keys = _split_point_keys(mesh)
    seg_ref = next((k for k in sorted(mesh.cell_data.keys()) if ":ref" in k), None)

    with open(filename, "w") as fh:
        fh.write(_provenance.render_lines(_provenance.SlotTier.BLOCK, "# "))
        fh.write(f"{len(mesh.points)} 2 {len(attr_keys)} {len(ref_keys)}\n")
        _write_node_rows(fh, mesh, attr_keys, ref_keys)

        line_blocks = [(ci, c) for ci, c in enumerate(mesh.cells) if c.type == "line"]
        fh.write(
            f"{sum(len(c) for _, c in line_blocks)} {0 if seg_ref is None else 1}\n"
        )
        sid = 0
        for ci, cell_block in line_blocks:
            for i, seg in enumerate(cell_block.data):
                row = [str(sid), str(seg[0]), str(seg[1])]
                if seg_ref is not None:
                    blocks = mesh.cell_data[seg_ref]
                    row.append(_fmt_value(blocks[ci][i]) if ci < len(blocks) else "0")
                fh.write(" ".join(row) + "\n")
                sid += 1
        fh.write("0\n")  # holes


def write(filename, mesh):
    filename = pathlib.Path(filename)
    if mesh.points.shape[1] != 2:
        raise WriteError("Triangle: can only write 2D points")
    if filename.suffix in (".node", ".ele"):
        _write_node_ele(filename, mesh)
    elif filename.suffix == ".poly":
        _write_poly(filename, mesh)
    else:
        raise WriteError(
            f"Triangle: must specify a .node, .ele, or .poly file. Got {filename}."
        )
