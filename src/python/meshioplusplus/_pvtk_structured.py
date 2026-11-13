"""VTK XML parallel indexes of structured data: ``.pvti``, ``.pvtr`` and ``.pvts``.

A ``<PImageData>``, ``<PRectilinearGrid>`` or ``<PStructuredGrid>`` index holds the
whole grid's ``WholeExtent`` and one ``<Piece Extent=... Source=...>`` per part, each
a standalone ``.vti``, ``.vtr`` or ``.vts`` file. Reading places every piece by its
``Extent`` in the ``WholeExtent`` lattice, so the points of a plane two pieces share
become one layer of points (nothing is welded: the placement is the weld), and the
cells are the whole lattice's. Writing writes one piece holding the whole lattice.
Python only.
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET

import numpy as np

from ._exceptions import ReadError
from ._grid import lattice_cells  # noqa: F401
from ._mesh import Mesh
from ._regions import Region

#: kind -> (index root tag, piece extension, piece package)
KINDS = {
    "pvti": ("PImageData", ".vti", "vti"),
    "pvtr": ("PRectilinearGrid", ".vtr", "vtr"),
    "pvts": ("PStructuredGrid", ".vts", "vts"),
}


def _ints(text, what, kind):
    try:
        values = [int(v) for v in text.split()]
    except (AttributeError, ValueError):
        values = []
    if len(values) != 6:
        raise ReadError(f"meshio++: {kind}: {what} needs six integers")
    return values


def _parse(filename, kind):
    root_tag, _ext, _pkg = KINDS[kind]
    try:
        tree = ET.parse(str(filename))
    except ET.ParseError as exc:
        raise ReadError(
            f"meshio++: {kind}: could not parse {filename}: {exc}"
        ) from None
    root = tree.getroot()
    if root.tag != "VTKFile" or root.get("type") != root_tag:
        raise ReadError(
            f"meshio++: {kind}: expected a VTKFile of type {root_tag}: {filename}"
        )
    body = root.find(root_tag)
    if body is None:
        raise ReadError(f"meshio++: {kind}: expected tag '{root_tag}': {filename}")
    whole = _ints(body.get("WholeExtent"), "WholeExtent", kind)
    pieces = []
    for el in body.findall("Piece"):
        source = el.get("Source")
        if not source:
            raise ReadError(f"meshio++: {kind}: a piece has no Source: {filename}")
        path = (
            source
            if os.path.isabs(source)
            else os.path.join(os.path.dirname(str(filename)), source)
        )
        if not os.path.exists(path):
            raise ReadError(
                f"meshio++: {kind}: piece file '{source}' named by {filename} does not exist"
            )
        pieces.append((_ints(el.get("Extent"), "a piece's Extent", kind), path))
    return whole, pieces


def _local_to_global(ext, whole, per_axis_cells):
    """Global lattice indices of a piece's points (or cells), ``i`` fastest."""
    n = [whole[2 * a + 1] - whole[2 * a] + 1 for a in range(3)]
    count = [ext[2 * a + 1] - ext[2 * a] + 1 for a in range(3)]
    off = [ext[2 * a] - whole[2 * a] for a in range(3)]
    if per_axis_cells:
        n = [max(v - 1, 1) for v in n]
        count = [max(v - 1, 1) for v in count]
    k, j, i = np.meshgrid(
        np.arange(count[2]) + off[2],
        np.arange(count[1]) + off[1],
        np.arange(count[0]) + off[0],
        indexing="ij",
    )
    return (i + n[0] * (j + n[1] * k)).reshape(-1)


def read(filename, kind, piece=None):
    import importlib

    whole, pieces = _parse(filename, kind)
    _root, ext_name, pkg = KINDS[kind]
    reader = importlib.import_module(f".{pkg}", __package__).read
    n = [whole[2 * a + 1] - whole[2 * a] + 1 for a in range(3)]
    if min(n) < 1:
        raise ReadError(f"meshio++: {kind}: WholeExtent holds no points")
    if piece is not None:
        index = piece + len(pieces) if piece < 0 else piece
        if not 0 <= index < len(pieces):
            raise ReadError(
                f"meshio++: {kind}: piece {piece} is out of range; the index has "
                f"{len(pieces)} piece(s)"
            )
        return reader(pieces[index][1])
    num_points = n[0] * n[1] * n[2]
    points = np.zeros((num_points, 3))
    placed = np.zeros(num_points, dtype=bool)
    point_data, cell_data, regions = {}, {}, []
    cell_type, conn = lattice_cells(*n)
    num_cells = len(conn)
    cells_global = np.zeros(num_cells, dtype=bool)
    for number, (ext, path) in enumerate(pieces):
        mesh = reader(path)
        gp = _local_to_global(ext, whole, False)
        if len(mesh.points) != len(gp):
            raise ReadError(
                f"meshio++: {kind}: piece '{path}' has {len(mesh.points)} points but "
                f"its Extent holds {len(gp)}"
            )
        points[gp] = mesh.points
        placed[gp] = True
        for name, arr in mesh.point_data.items():
            arr = np.asarray(arr)
            if name not in point_data:
                point_data[name] = np.zeros((num_points,) + arr.shape[1:], arr.dtype)
            point_data[name][gp] = arr
        cells_here = None
        if len(mesh.cells):
            gc = _local_to_global(ext, whole, True)
            if len(gc) == len(mesh.cells[0].data):
                cells_here = gc
                cells_global[gc] = True
                for name, blocks in mesh.cell_data.items():
                    arr = np.asarray(blocks[0])
                    if name not in cell_data:
                        cell_data[name] = np.zeros(
                            (num_cells,) + arr.shape[1:], arr.dtype
                        )
                    cell_data[name][gc] = arr
        if cells_here is not None:
            regions.append(
                Region(
                    f"piece_{number}",
                    "cell",
                    np.sort(cells_here).astype(np.int64),
                    tag=number,
                )
            )
    if not placed.all():
        raise ReadError(f"meshio++: {kind}: the pieces do not cover WholeExtent")
    out_cells = [(cell_type, conn.astype(np.int64))] if num_cells else []
    return Mesh(
        points,
        out_cells,
        point_data=point_data,
        cell_data={k: [v] for k, v in cell_data.items()} if num_cells else {},
        regions=regions,
    )


def write(filename, mesh, kind, **kwargs):
    """The whole lattice as one piece (``<stem>/<stem>_0000.vt?``) and an index."""
    import importlib

    root_tag, ext, pkg = KINDS[kind]
    filename = str(filename)
    stem = os.path.splitext(os.path.basename(filename))[0]
    parent = os.path.dirname(filename)
    piece_dir = os.path.join(parent, stem) if parent else stem
    os.makedirs(piece_dir, exist_ok=True)
    piece_name = f"{stem}_0000{ext}"
    piece_path = os.path.join(piece_dir, piece_name)
    importlib.import_module(f".{pkg}", __package__).write(piece_path, mesh, **kwargs)
    # The piece's own extent and grid attributes, read back from what was written.
    grid = (
        ET.parse(piece_path)
        .getroot()
        .find(
            {"pvti": "ImageData", "pvtr": "RectilinearGrid", "pvts": "StructuredGrid"}[
                kind
            ]
        )
    )
    extent = grid.get("WholeExtent")
    attrs = {k: v for k, v in grid.attrib.items() if k != "WholeExtent"}
    root = ET.Element(
        "VTKFile", type=root_tag, version="1.0", byte_order="LittleEndian"
    )
    body = ET.SubElement(root, root_tag, WholeExtent=extent, **attrs)
    first = ET.parse(piece_path).getroot().find(grid.tag).find("Piece")
    for section, decl in (("PPointData", "PointData"), ("PCellData", "CellData")):
        src = first.find(decl)
        if src is not None:
            node = ET.SubElement(body, section)
            for da in src.findall("DataArray"):
                ET.SubElement(
                    node,
                    "PDataArray",
                    {
                        k: v
                        for k, v in da.attrib.items()
                        if k in ("type", "Name", "NumberOfComponents")
                    },
                )
    pts = first.find("Points")
    if pts is not None:
        ppts = ET.SubElement(body, "PPoints")
        for da in pts.findall("DataArray"):
            ET.SubElement(
                ppts,
                "PDataArray",
                {
                    k: v
                    for k, v in da.attrib.items()
                    if k in ("type", "NumberOfComponents")
                },
            )
    coords = first.find("Coordinates")
    if coords is not None:
        pc = ET.SubElement(body, "PCoordinates")
        for da in coords.findall("DataArray"):
            ET.SubElement(
                pc, "PDataArray", {k: v for k, v in da.attrib.items() if k in ("type",)}
            )
    ET.SubElement(body, "Piece", Extent=extent, Source=f"{stem}/{piece_name}")
    ET.ElementTree(root).write(filename, xml_declaration=True, encoding="utf-8")
