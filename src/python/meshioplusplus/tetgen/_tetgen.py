"""
I/O for the TetGen file format, c.f.
<https://wias-berlin.de/software/tetgen/fformats.node.html>
"""

import pathlib

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._mesh import CellBlock, Mesh


def _header_line(f):
    """The first line that is neither blank nor a comment; ReadError at EOF."""
    # readline, not `for raw in f`: iterating a text file disables tell(),
    # which np.fromfile needs to pick up where the header ended.
    while True:
        raw = f.readline()
        if not raw:
            raise ReadError(f"TetGen: {f.name} has no header line")
        line = raw.strip()
        if line and line[0] != "#":
            return line


def _tokens(filename):
    """The non-comment tokens of a TetGen text file, and its header's."""
    header, data = None, []
    with open(filename) as f:
        for raw in f:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            if header is None:
                header = line.split()
            else:
                data.extend(line.split())
    if header is None:
        raise ReadError(f"TetGen: missing header line in {filename}")
    return header, data


def _read_smesh(filename):
    """A ``.smesh`` surface: nodes, then polygonal facets (roadmap 7.7.2)."""
    header, data = _tokens(filename)
    npoints = int(header[0])
    dim = int(header[1]) if len(header) > 1 else 3
    nattr = int(header[2]) if len(header) > 2 else 0
    nmark = int(header[3]) if len(header) > 3 else 0
    if dim != 3:
        raise ReadError("TetGen: need 3D points")
    ncol = 4 + nattr + nmark
    if npoints < 0 or npoints * ncol > len(data):
        raise ReadError("TetGen: .smesh node section is truncated")
    nodes = np.array(data[: npoints * ncol], dtype=float).reshape(npoints, ncol)
    base = int(nodes[0, 0]) if npoints else 0
    if not np.array_equal(nodes[:, 0], np.arange(base, base + npoints)):
        raise ReadError("TetGen: nodes not numbered consecutively")
    points = np.ascontiguousarray(nodes[:, 1:4])
    at = npoints * ncol

    def take():
        nonlocal at
        if at >= len(data):
            raise ReadError("TetGen: .smesh facet section is truncated")
        at += 1
        return int(data[at - 1])

    nfacets, marked = take(), take()
    if nfacets < 0 or marked not in (0, 1):
        raise ReadError("TetGen: malformed .smesh facet header")
    runs = []
    for _ in range(nfacets):
        n = take()
        if n < 3:
            raise ReadError("TetGen: a facet has fewer than three corners")
        corners = [take() - base for _ in range(n)]
        if min(corners) < 0 or max(corners) >= npoints:
            raise ReadError("TetGen: a facet names a node out of range")
        if not runs or len(runs[-1][0][0]) != n:
            runs.append(([], []))
        runs[-1][0].append(corners)
        if marked:
            runs[-1][1].append(take())
    cells = []
    for conn, _ in runs:
        n = len(conn[0])
        kind = {3: "triangle", 4: "quad"}.get(n, "polygon")
        cells.append(CellBlock(kind, np.array(conn, dtype=np.int64)))
    cell_data = {}
    if marked:
        cell_data["tetgen:ref"] = [np.array(m, dtype=np.int64) for _, m in runs]
    return Mesh(points, cells, cell_data=cell_data)


def _write_smesh(filename, mesh, float_fmt=".16e"):
    if mesh.points.shape[1] != 3:
        raise WriteError("TetGen: can only write 3D points")
    for block in mesh.cells:
        if not (
            block.type in ("triangle", "quad")
            or (
                block.type == "polygon"
                and isinstance(block.data, np.ndarray)
                and block.data.ndim == 2
            )
        ):
            raise WriteError(
                f"TetGen: a .smesh holds polygonal facets, not '{block.type}' cells"
            )
    marked = "tetgen:ref" in mesh.cell_data
    with open(filename, "w") as fh:
        fh.write(_provenance.render_lines(_provenance.SlotTier.BLOCK, "# "))
        fh.write(f"# part 1: nodes\n{mesh.points.shape[0]} 3 0 0\n")
        for k, p in enumerate(mesh.points):
            fh.write(f"{k} {p[0]:{float_fmt}} {p[1]:{float_fmt}} {p[2]:{float_fmt}}\n")
        nfacets = sum(len(block.data) for block in mesh.cells)
        fh.write(f"# part 2: facets\n{nfacets} {1 if marked else 0}\n")
        for b, block in enumerate(mesh.cells):
            data = np.asarray(block.data)
            marks = None
            if marked and b < len(mesh.cell_data["tetgen:ref"]):
                marks = np.asarray(mesh.cell_data["tetgen:ref"][b])
            for c, row in enumerate(data):
                line = f"{len(row)} " + " ".join(str(int(v)) for v in row)
                if marked:
                    line += f" {int(marks[c]) if marks is not None else 0}"
                fh.write(line + "\n")
        fh.write("# part 3: holes\n0\n# part 4: regions\n0\n")


def read(filename):
    filename = pathlib.Path(filename)
    if filename.suffix == ".smesh":
        return _read_smesh(filename)
    if filename.suffix == ".node":
        node_filename = filename
        ele_filename = filename.parent / (filename.stem + ".ele")
    elif filename.suffix == ".ele":
        node_filename = filename.parent / (filename.stem + ".node")
        ele_filename = filename
    else:
        raise ReadError()

    point_data = {}
    cell_data = {}

    # read nodes
    with open(node_filename) as f:
        line = _header_line(f)

        num_points, dim, num_attrs, num_bmarkers = (
            int(item) for item in line.split(" ") if item != ""
        )
        if dim != 3:
            raise ReadError("Need 3D points.")

        points = np.fromfile(
            f, dtype=float, count=(4 + num_attrs + num_bmarkers) * num_points, sep=" "
        ).reshape(num_points, 4 + num_attrs + num_bmarkers)

        node_index_base = int(points[0, 0])
        # make sure the nodes are numbered consecutively
        if not np.all(
            points[:, 0]
            == np.arange(node_index_base, node_index_base + points.shape[0])
        ):
            raise ReadError()
        # read point attributes
        for k in range(num_attrs):
            point_data[f"tetgen:attr{k + 1}"] = points[:, 4 + k]
        # read boundary markers, the first is "ref", the others are "ref2", "ref3", ...
        for k in range(num_bmarkers):
            flag = "" if k == 0 else str(k + 1)
            point_data["tetgen:ref" + flag] = points[:, 4 + num_attrs + k]
        # remove the leading index column, the attributes, and the boundary markers
        points = points[:, 1:4]

    # read elements
    with open(ele_filename.as_posix()) as f:
        line = _header_line(f)

        num_tets, num_points_per_tet, num_attrs = (
            int(item) for item in line.strip().split(" ") if item != ""
        )
        if num_points_per_tet != 4:
            raise ReadError()
        cells = np.fromfile(
            f, dtype=int, count=(5 + num_attrs) * num_tets, sep=" "
        ).reshape(num_tets, 5 + num_attrs)
        # read cell (region) attributes, the first is "ref", the others are "ref2",
        # "ref3", ...
        for k in range(num_attrs):
            flag = "" if k == 0 else str(k + 1)
            cell_data["tetgen:ref" + flag] = [cells[:, 5 + k]]
        # remove the leading index column and the attributes
        cells = cells[:, 1:5]
        cells -= node_index_base

    return Mesh(
        points, [CellBlock("tetra", cells)], point_data=point_data, cell_data=cell_data
    )


def write(filename, mesh, float_fmt=".16e"):
    filename = pathlib.Path(filename)
    if filename.suffix == ".smesh":
        return _write_smesh(filename, mesh, float_fmt)
    if filename.suffix == ".node":
        node_filename = filename
        ele_filename = filename.parent / (filename.stem + ".ele")
    elif filename.suffix == ".ele":
        node_filename = filename.parent / (filename.stem + ".node")
        ele_filename = filename
    else:
        raise WriteError(f"Must specify .node, .ele or .smesh file. Got {filename}.")

    if mesh.points.shape[1] != 3:
        raise WriteError("Can only write 3D points")

    # write nodes
    with open(node_filename, "w") as fh:
        # identify ":ref" key
        attr_keys = list(mesh.point_data.keys())
        ref_keys = [k for k in attr_keys if ":ref" in k]
        if len(attr_keys) > 0:
            if len(ref_keys) > 0:
                ref_keys = ref_keys[:1]
                attr_keys.remove(ref_keys[0])
            else:
                ref_keys = attr_keys[:1]
                attr_keys = attr_keys[1:]

        nattr, nref = len(attr_keys), len(ref_keys)
        fh.write(_provenance.render_lines(_provenance.SlotTier.BLOCK, "# "))
        if (nattr + nref) > 0:
            fh.write(
                "# attribute and marker names: {}\n".format(
                    ", ".join(attr_keys + ref_keys)
                )
            )
        fh.write(f"{mesh.points.shape[0]} {3} {nattr} {nref}\n")
        fmt = (
            "{} "
            + " ".join((3 + nattr) * ["{:" + float_fmt + "}"])
            + "".join((nref) * [" {}"])
            + "\n"
        )
        for k, pt in enumerate(mesh.points):
            data = (
                list(pt[:3])
                + [mesh.point_data[key][k] for key in attr_keys]
                + [mesh.point_data[key][k] for key in ref_keys]
            )
            fh.write(fmt.format(k, *data))

    if any(c.type != "tetra" for c in mesh.cells):
        string = ", ".join([c.type for c in mesh.cells if c.type != "tetra"])
        warn(f"TetGen only supports tetrahedra, but mesh has {string}. Skipping those.")

    # write cells
    with open(ele_filename, "w") as fh:
        attr_keys = list(mesh.cell_data.keys())
        ref_keys = [k for k in attr_keys if ":ref" in k]
        if len(attr_keys) > 0:
            if len(ref_keys) > 0:
                attr_keys.remove(ref_keys[0])
                attr_keys = ref_keys[:1] + attr_keys

        nattr = len(attr_keys)
        fh.write(_provenance.render_lines(_provenance.SlotTier.BLOCK, "# "))
        if nattr > 0:
            fh.write("# attribute names: {}\n".format(", ".join(attr_keys)))
        for id, c in enumerate(filter(lambda c: c.type == "tetra", mesh.cells)):
            data = c.data
            fh.write(f"{data.shape[0]} {4} {nattr}\n")
            fmt = " ".join((5 + nattr) * ["{}"]) + "\n"
            for k, tet in enumerate(data):
                data = list(tet[:4]) + [mesh.cell_data[key][id][k] for key in attr_keys]
                fh.write(fmt.format(k, *data))
