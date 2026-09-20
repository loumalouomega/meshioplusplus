"""VTK XML MultiBlock, the pure-Python reference (v11.6.0, roadmap §1 tier B4,
part 3 of 3).

The index carries no geometry at all: it is a list of `<DataSet file="..."/>`
entries, each pointing at a standalone `.vtu` (or `.vtp`) piece. Writing means
carving the mesh into one piece per cell block -- pruned to only the points
that block references, via `meshioplusplus.clean` -- and delegating each
piece to the best available `.vtu` writer. Reading means parsing the index,
reading every piece with the best available reader for its own extension,
and combining them with the internal numpy `merge()` (no welding: multiblock
pieces are pre-separated by construction, not coincident-point fragments to
weld back together). Each piece becomes one `"cell"` region in the merged
mesh, named from the index's `name=` attribute, via `merge()`'s own per-input
cell index map -- correct even when same-typed pieces consolidate into one
output cell block.
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET

import numpy as np

from .. import _provenance
from .._clean import clean
from .._exceptions import ReadError
from .._merge import _merge_py
from .._mesh import Mesh
from .._pvtk_index import share_field_data
from .._regions import Region


def _collect_datasets(node, out):
    for child in node:
        if child.tag == "DataSet":
            out.append(child)
        else:
            _collect_datasets(child, out)


def _read_piece(path):
    ext = os.path.splitext(path)[1]
    if ext == ".vtu":
        from .. import vtu as pkg
    elif ext == ".vtp":
        from .. import vtp as pkg
    else:
        raise ReadError(
            f"Unsupported .vtm piece '{path}': only .vtu/.vtp pieces are read"
        )
    return pkg.read(path)


def read(filename):
    tree = ET.parse(str(filename))
    root = tree.getroot()
    if root.tag != "VTKFile":
        raise ReadError("Expected tag 'VTKFile'")
    if root.get("type") != "vtkMultiBlockDataSet":
        raise ReadError("Expected type vtkMultiBlockDataSet")
    mb = root.find("vtkMultiBlockDataSet")
    if mb is None:
        raise ReadError("Expected tag 'vtkMultiBlockDataSet'")

    datasets = []
    _collect_datasets(mb, datasets)
    if not datasets:
        return Mesh(np.zeros((0, 3), dtype=np.float64), [])

    base = os.path.dirname(str(filename))
    pieces = []
    names = []
    for i, ds in enumerate(datasets):
        file_attr = ds.get("file")
        if not file_attr:
            raise ReadError("<DataSet> is missing its 'file' attribute")
        path = os.path.join(base, file_attr) if base else file_attr
        pieces.append(_read_piece(path))
        names.append(ds.get("name") or f"block_{i}")

    out, _point_maps, cell_maps = _merge_py(
        pieces,
        weld=False,
        atol=1e-8,
        source_tag=True,
        data_policy="fill",
        drop_duplicate_cells=False,
    )
    # Dataset-global field data is the union across pieces, not namespaced.
    share_field_data(out, pieces)
    for name, cmap in zip(names, cell_maps):
        out.regions.append(Region(name, "cell", np.asarray(cmap, dtype=np.int64)))
    return out


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    filename = str(filename)
    stem = os.path.splitext(os.path.basename(filename))[0]
    parent = os.path.dirname(filename)
    piece_dir = os.path.join(parent, stem) if parent else stem
    os.makedirs(piece_dir, exist_ok=True)

    from .. import vtu as vtu_pkg

    piece_files = []
    piece_names = []
    for i, cb in enumerate(mesh.cells):
        piece = Mesh(
            mesh.points,
            [(cb.type, cb.data)],
            point_data=dict(mesh.point_data),
            cell_data={k: [v[i]] for k, v in mesh.cell_data.items() if i < len(v)},
        )
        cleaned = clean(
            piece,
            weld=False,
            remove_orphans=True,
            drop_degenerate=False,
            drop_duplicate_cells=False,
        )
        piece_file_name = f"{stem}_{i}.vtu"
        vtu_pkg.write(
            os.path.join(piece_dir, piece_file_name),
            cleaned,
            binary=binary,
            compression=compression,
            header_type=header_type,
        )
        piece_files.append(f"{stem}/{piece_file_name}")
        piece_names.append(f"block_{i}")

    lines = ['<?xml version="1.0"?>']
    lines.append(
        '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">'
    )
    lines.append(_provenance.render_xml_comment(_provenance.SlotTier.BLOCK))
    lines.append("<vtkMultiBlockDataSet>")
    lines.append('<Block index="0">')
    for i, (f, name) in enumerate(zip(piece_files, piece_names)):
        lines.append(f'<DataSet index="{i}" name="{name}" file="{f}"/>')
    lines.append("</Block>")
    lines.append("</vtkMultiBlockDataSet>")
    lines.append("</VTKFile>")
    with open(filename, "w") as f:
        f.write("\n".join(lines) + "\n")
