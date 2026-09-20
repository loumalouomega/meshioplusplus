"""Machinery shared by the ParaView index formats ``.pvtu``, ``.pvtp`` and ``.pvd``.

Each is an XML *index* that carries no geometry of its own, only the relative
paths of the piece files that do. This module holds what the three have in
common so the per-format modules stay short:

* resolving a ``file=`` / ``Source=`` attribute against the index's own
  directory, and XML-escaping one on the way out;
* the declaration check a parallel index has to satisfy (every piece declares
  identical arrays, validated *before* anything is written);
* the ``vtkGhostType`` vocabulary, and the translation of ``partition``'s halo
  layers (``partition:ghost``) into it;
* the parallel-index reader and writer themselves, parameterised on the kind
  (``pvtu`` writes ``.vtu`` pieces, ``pvtp`` writes ``.vtp`` pieces).
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from xml.sax.saxutils import quoteattr

import numpy as np

from . import _provenance
from ._clean import clean
from ._exceptions import ReadError, WriteError
from ._merge import _merge_py, _merge_sets
from ._mesh import Mesh
from ._regions import Region

PARTITION_PART = "partition:part"
PARTITION_GHOST = "partition:ghost"

# VTK's ghost vocabulary (vtkDataSetAttributes). Cells: DUPLICATECELL = 1,
# REFINEDCELL = 8, HIDDENCELL = 32. Points: DUPLICATEPOINT = 1, HIDDENPOINT = 32.
GHOST_NAME = "vtkGhostType"
DUPLICATECELL = 1
REFINEDCELL = 8
HIDDENCELL = 32
DUPLICATEPOINT = 1
HIDDENPOINT = 32

# kind -> (root element, piece extension, piece package)
KINDS = {
    "pvtu": ("PUnstructuredGrid", ".vtu", "vtu"),
    "pvtp": ("PPolyData", ".vtp", "vtp"),
}

# Which serial files an index may name. Deliberately no recursion inside a
# parallel index (a `.pvtu` names `.vtu`/`.vtp` only), so the only nesting is a
# `.pvd` over the wider set: no cycle is expressible.
PIECE_EXTENSIONS = (".vtu", ".vtp")
COLLECTION_EXTENSIONS = (".vtu", ".vtp", ".vtm", ".pvtu", ".pvtp")

GHOST_POLICIES = ("keep", "drop")


def quote(value) -> str:
    """An XML attribute value, quoted and escaped (``&``, ``<``, quotes)."""
    return quoteattr(str(value))


def check_ghosts(ghosts):
    if ghosts not in GHOST_POLICIES:
        raise ValueError(f"meshio++: ghosts must be 'keep' or 'drop', got {ghosts!r}")


def resolve_piece(piece, count):
    """Resolve ``piece`` (negative counts from the end) against ``count``."""
    k = int(piece)
    resolved = count + k if k < 0 else k
    if resolved < 0 or resolved >= count:
        raise ReadError(
            f"meshio++: piece {k} is out of range: this file has {count} "
            f"{'piece' if count == 1 else 'pieces'}"
        )
    return resolved


def resolve_path(index_path, file_attr, fmt, attr="file"):
    """Where a piece named by ``file_attr`` lives.

    Resolved against the index's own directory; an absolute path is used as
    written. There is deliberately no search for the file elsewhere: an absolute
    path from another machine that does not exist here is an error naming the
    attribute and the directory, because a fallback could read the wrong file.
    """
    index_path = str(index_path)
    if not file_attr:
        raise ReadError(
            f"meshio++: {fmt}: an entry is missing its '{attr}' attribute: {index_path}"
        )
    base = os.path.dirname(index_path)
    path = file_attr if os.path.isabs(file_attr) else os.path.join(base, file_attr)
    if not os.path.exists(path):
        raise ReadError(
            f"meshio++: {fmt}: piece file '{file_attr}' ({attr}=) named by "
            f"{index_path} does not exist (looked for {path}); a relative path is "
            "resolved against the index's own directory and an absolute path is "
            "not searched for elsewhere"
        )
    return path


def piece_layout(index_path, count, ext):
    """``(directory, [(path, index-relative posix path), ...])`` for ``count`` pieces.

    Pieces live in a sibling directory named after the index's stem, zero-padded
    to ``max(4, digits(count - 1))`` (the width ``{step}`` patterns use).
    """
    index_path = str(index_path)
    stem = os.path.splitext(os.path.basename(index_path))[0]
    parent = os.path.dirname(index_path)
    piece_dir = os.path.join(parent, stem) if parent else stem
    width = max(4, len(str(max(count - 1, 0))))
    names = [f"{stem}_{i:0{width}d}{ext}" for i in range(count)]
    return piece_dir, [(os.path.join(piece_dir, n), f"{stem}/{n}") for n in names]


def read_child(path, ghosts="keep", extensions=PIECE_EXTENSIONS):
    """Read one piece file with the best available reader for its extension."""
    ext = os.path.splitext(path)[1].lower()
    if ext not in extensions:
        raise ReadError(
            f"meshio++: unsupported piece '{path}': only "
            f"{'/'.join(extensions)} pieces are read"
        )
    from . import pvtp, pvtu, vtm, vtp, vtu

    if ext == ".vtu":
        return vtu.read(path)
    if ext == ".vtp":
        return vtp.read(path)
    if ext == ".vtm":
        return vtm.read(path)
    if ext == ".pvtu":
        return pvtu.read(path, ghosts=ghosts)
    return pvtp.read(path, ghosts=ghosts)


def empty_mesh():
    return Mesh(np.zeros((0, 3), dtype=np.float64), [])


# --------------------------------------------------------------------------- #
# declarations                                                                 #
# --------------------------------------------------------------------------- #
def _vtu_type(name, arr):
    from .vtu._vtu import numpy_to_vtu_type

    dt = np.asarray(arr).dtype.newbyteorder("=")
    if name == GHOST_NAME:
        dt = np.dtype(np.uint8)  # the piece writers always store it as UInt8
    try:
        return numpy_to_vtu_type[dt]
    except KeyError:
        raise WriteError(
            f"meshio++: '{name}' has dtype {dt}, which VTK XML cannot declare "
            f"(supported: {', '.join(sorted(numpy_to_vtu_type.values()))})"
        ) from None


def _ncomp(arr):
    arr = np.asarray(arr)
    return int(arr.shape[1]) if arr.ndim == 2 else 1


def declarations(mesh):
    """``(points, point_data, cell_data)`` declarations of what a piece will hold.

    A declaration is ``(vtk type name, number of components)``; the two dicts
    map an array name to one. This is exactly what ``<PDataArray>`` carries, so
    it is what every piece of a parallel index has to agree on.
    """
    pts = np.asarray(mesh.points)
    points = (_vtu_type("Points", pts), 3)
    point_data = {
        k: (_vtu_type(k, v), _ncomp(v)) for k, v in sorted(mesh.point_data.items())
    }
    cell_data = {}
    for k, blocks in sorted(mesh.cell_data.items()):
        decl = None
        for blk in blocks:
            d = (_vtu_type(k, blk), _ncomp(blk))
            if decl is None:
                decl = d
            elif decl != d:
                raise WriteError(
                    f"meshio++: cell_data '{k}' has differing types or component "
                    "counts across cell blocks"
                )
        if decl is not None:
            cell_data[k] = decl
    return points, point_data, cell_data


def _fmt_decl(d):
    return f"{d[0]} with {d[1]} component{'s' if d[1] != 1 else ''}"


def validate_declarations(kind, meshes):
    """Refuse a piece list whose pieces do not declare identical arrays.

    Runs before anything is emitted, so a refusal leaves no half-written tree.
    Piece 0 is the reference. A piece with no cells at all is exempt from the
    cell_data comparison only when it carries no cell arrays (an idle rank
    legitimately allocated nothing); points are never exempt.
    """
    decls = [declarations(m) for m in meshes]
    ref_points, ref_pd, ref_cd = decls[0]
    for k, (points, pd, cd) in enumerate(decls[1:], start=1):
        if points != ref_points:
            raise WriteError(
                f"meshio++: {kind}: piece {k} declares Points as {_fmt_decl(points)}, "
                f"but piece 0 declares it as {_fmt_decl(ref_points)}; every piece of "
                "a parallel index must declare identical arrays (name, type, "
                "NumberOfComponents)"
            )
        for label, ref, cur in (("point_data", ref_pd, pd), ("cell_data", ref_cd, cd)):
            if label == "cell_data" and not cur and _num_cells(meshes[k]) == 0:
                continue
            for name in sorted(set(ref) | set(cur)):
                if name not in cur:
                    raise WriteError(
                        f"meshio++: {kind}: piece {k} is missing {label} '{name}', "
                        "which piece 0 declares; every piece of a parallel index "
                        "must declare identical arrays (name, type, NumberOfComponents)"
                    )
                if name not in ref:
                    raise WriteError(
                        f"meshio++: {kind}: piece {k} declares {label} '{name}', "
                        "which piece 0 does not; every piece of a parallel index "
                        "must declare identical arrays (name, type, NumberOfComponents)"
                    )
                if ref[name] != cur[name]:
                    raise WriteError(
                        f"meshio++: {kind}: piece {k} declares {label} '{name}' as "
                        f"{_fmt_decl(cur[name])}, but piece 0 declares it as "
                        f"{_fmt_decl(ref[name])}; every piece of a parallel index "
                        "must declare identical arrays (name, type, NumberOfComponents)"
                    )
    return decls[0]


def _num_cells(mesh):
    return sum(len(cb.data) for cb in mesh.cells)


# --------------------------------------------------------------------------- #
# ghosts                                                                       #
# --------------------------------------------------------------------------- #
def _rectangular(block):
    data = block.data
    return isinstance(data, np.ndarray) and data.dtype != object and data.ndim == 2


def ghost_level(meshes):
    """``max(partition:ghost)`` over the pieces: the number of halo layers."""
    level = 0
    for m in meshes:
        for blk in m.cell_data.get(PARTITION_GHOST, []):
            blk = np.asarray(blk)
            if blk.size:
                level = max(level, int(blk.max()))
    return level


def with_ghost_arrays(mesh):
    """A copy of ``mesh`` carrying ``vtkGhostType`` derived from ``partition:ghost``.

    Cells: 0 for an owned cell (layer 0), ``DUPLICATECELL`` for any halo layer.
    Points: ``DUPLICATEPOINT`` when no owned cell of this piece references the
    point. An array the caller already supplied is passed through unchanged, and
    nothing is fabricated when ``partition:ghost`` is absent. The caller's mesh
    is never mutated.
    """
    layers = mesh.cell_data.get(PARTITION_GHOST)
    have_cell = GHOST_NAME in mesh.cell_data
    have_point = GHOST_NAME in mesh.point_data
    if layers is None or (have_cell and have_point):
        return mesh

    cell_data = dict(mesh.cell_data)
    point_data = dict(mesh.point_data)
    owned = np.zeros(len(mesh.points), dtype=bool)
    cell_gh = []
    for blk, lay in zip(mesh.cells, layers):
        lay = np.asarray(lay)
        cell_gh.append((lay > 0).astype(np.uint8) * np.uint8(DUPLICATECELL))
        if not _rectangular(blk):
            raise NotImplementedError(
                "ghost translation of polygon/polyhedron blocks needs the C++ core"
            )
        if len(blk.data):
            owned[np.asarray(blk.data)[lay == 0].ravel()] = True
    if not have_cell:
        cell_data[GHOST_NAME] = cell_gh
    if not have_point:
        point_data[GHOST_NAME] = np.where(owned, 0, DUPLICATEPOINT).astype(np.uint8)
    return Mesh(
        mesh.points,
        [(cb.type, cb.data) for cb in mesh.cells],
        point_data=point_data,
        cell_data=cell_data,
        field_data=dict(mesh.field_data),
    )


def drop_ghosts(mesh):
    """Remove ghost cells (any ``vtkGhostType`` bit set) and the points only they used.

    When a cell is dropped, every point no kept cell references is pruned (the
    ``crop`` rule the C++ core applies); empty blocks are kept, so the block
    structure is unchanged. With no flagged cell only the arrays go. The
    now-meaningless ``vtkGhostType`` arrays are removed, and ``partition:ghost``
    with them when it is all zero.
    """
    cell_gh = mesh.cell_data.get(GHOST_NAME)
    if cell_gh is None:
        return mesh
    npts = len(mesh.points)
    keeps = []
    for blk, gh in zip(mesh.cells, cell_gh):
        if not _rectangular(blk):
            raise NotImplementedError(
                "dropping ghosts from polygon/polyhedron blocks needs the C++ core"
            )
        keeps.append(np.asarray(gh) == 0)

    if all(k.all() for k in keeps):
        points = mesh.points
        cells = [(cb.type, cb.data) for cb in mesh.cells]
        point_sel = None
    else:
        referenced = np.zeros(npts, dtype=bool)
        for blk, keep in zip(mesh.cells, keeps):
            data = np.asarray(blk.data)
            if len(data):
                referenced[data[keep].ravel()] = True
        new_index = np.cumsum(referenced) - 1
        points = np.asarray(mesh.points)[referenced]
        cells = [
            (blk.type, new_index[np.asarray(blk.data)[keep]])
            for blk, keep in zip(mesh.cells, keeps)
        ]
        point_sel = referenced

    point_data = {
        k: (np.asarray(v) if point_sel is None else np.asarray(v)[point_sel])
        for k, v in mesh.point_data.items()
        if k != GHOST_NAME
    }
    layers = mesh.cell_data.get(PARTITION_GHOST)
    drop_layers = layers is not None and all(
        not np.asarray(x)[keep].any() for x, keep in zip(layers, keeps)
    )
    cell_data = {}
    for k, blocks in mesh.cell_data.items():
        if k == GHOST_NAME or (k == PARTITION_GHOST and drop_layers):
            continue
        cell_data[k] = [np.asarray(x)[keep] for x, keep in zip(blocks, keeps)]
    return Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=dict(mesh.field_data),
    )


# --------------------------------------------------------------------------- #
# carving a single mesh into parts                                             #
# --------------------------------------------------------------------------- #
def carve_by_part(mesh, key=PARTITION_PART):
    """One piece per part id ``0..max`` of the integer ``cell_data[key]``.

    Without the key the whole mesh is one piece. Parts with no cells are kept
    (they are real ranks; skipping one would renumber the rest). Each piece is
    pruned to the points its cells reference.
    """
    labels = mesh.cell_data.get(key)
    if labels is None:
        return [mesh]
    labels = [np.asarray(x) for x in labels]
    for x in labels:
        if x.dtype.kind not in "iu":
            raise WriteError(f"meshio++: cell_data '{key}' must hold integer part ids")
    sizes = [int(x.max()) for x in labels if x.size]
    lows = [int(x.min()) for x in labels if x.size]
    if lows and min(lows) < 0:
        raise WriteError(f"meshio++: cell_data '{key}' holds a negative part id")
    nparts = (max(sizes) + 1) if sizes else 1
    for blk in mesh.cells:
        if not _rectangular(blk):
            raise NotImplementedError(
                "carving polygon/polyhedron blocks by part needs the C++ core"
            )

    pieces = []
    for p in range(nparts):
        cells = []
        cell_data = {k: [] for k in mesh.cell_data}
        for b, blk in enumerate(mesh.cells):
            sel = labels[b] == p
            cells.append((blk.type, np.asarray(blk.data)[sel]))
            for k, v in mesh.cell_data.items():
                if b < len(v):
                    cell_data[k].append(np.asarray(v[b])[sel])
        piece = Mesh(
            mesh.points,
            cells,
            point_data=dict(mesh.point_data),
            cell_data=cell_data,
            field_data=dict(mesh.field_data),
        )
        pieces.append(
            clean(
                piece,
                weld=False,
                remove_orphans=True,
                drop_degenerate=False,
                drop_duplicate_cells=False,
            )
        )
    return pieces


# --------------------------------------------------------------------------- #
# the parallel index: write and read                                           #
# --------------------------------------------------------------------------- #
def write_pieces(
    filename, pieces, kind, binary=True, compression="zlib", header_type=None
):
    """Write ``pieces`` (a list of meshes) as a parallel index plus one file each.

    The primitive under ``.pvtu``/``.pvtp`` writing: ``partition``'s output goes
    straight in, halo layers included.
    """
    root, ext, pkg_name = KINDS[kind]
    pieces = list(pieces)
    if not pieces:
        raise WriteError(f"meshio++: {kind}: a parallel index needs at least one piece")
    prepared = [with_ghost_arrays(m) for m in pieces]
    points, point_decl, cell_decl = validate_declarations(kind, prepared)
    level = ghost_level(pieces)

    filename = str(filename)
    piece_dir, layout = piece_layout(filename, len(prepared), ext)
    os.makedirs(piece_dir, exist_ok=True)

    import importlib

    pkg = importlib.import_module(f".{pkg_name}", __package__)
    for (path, _rel), mesh in zip(layout, prepared):
        pkg.write(
            path, mesh, binary=binary, compression=compression, header_type=header_type
        )

    lines = ['<?xml version="1.0"?>']
    lines.append(f'<VTKFile type="{root}" version="1.0" byte_order="LittleEndian">')
    lines.append(_provenance.render_xml_comment(_provenance.SlotTier.BLOCK))
    lines.append(f'<{root} GhostLevel="{level}">')
    for tag, decl in (("PPointData", point_decl), ("PCellData", cell_decl)):
        if not decl:
            continue
        lines.append(f"<{tag}>")
        for name, (typ, ncomp) in decl.items():
            comp = f' NumberOfComponents="{ncomp}"' if ncomp != 1 else ""
            lines.append(f'<PDataArray type="{typ}" Name={quote(name)}{comp}/>')
        lines.append(f"</{tag}>")
    lines.append("<PPoints>")
    lines.append(
        f'<PDataArray type="{points[0]}" Name="Points" NumberOfComponents="3"/>'
    )
    lines.append("</PPoints>")
    for _path, rel in layout:
        lines.append(f"<Piece Source={quote(rel)}/>")
    lines.append(f"</{root}>")
    lines.append("</VTKFile>")
    with open(filename, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_mesh(
    filename,
    mesh,
    kind,
    part_key=PARTITION_PART,
    binary=True,
    compression="zlib",
    header_type=None,
):
    """Write one mesh: one piece per ``partition:part`` id, else a single piece."""
    write_pieces(
        filename,
        carve_by_part(mesh, part_key),
        kind,
        binary=binary,
        compression=compression,
        header_type=header_type,
    )


def parse_index(filename, kind):
    """``[Source, ...]`` (resolved paths) of a parallel index, in document order."""
    root_tag, _ext, _pkg = KINDS[kind]
    try:
        tree = ET.parse(str(filename))
    except ET.ParseError as e:
        raise ReadError(f"meshio++: {kind}: could not parse {filename}: {e}") from None
    root = tree.getroot()
    if root.tag != "VTKFile":
        raise ReadError(f"meshio++: {kind}: expected tag 'VTKFile': {filename}")
    if root.get("type") != root_tag:
        raise ReadError(
            f"meshio++: {kind}: expected type {root_tag}, got "
            f"{root.get('type')!r}: {filename}"
        )
    body = root.find(root_tag)
    if body is None:
        raise ReadError(f"meshio++: {kind}: expected tag '{root_tag}': {filename}")
    return [
        resolve_path(filename, el.get("Source"), kind, "Source")
        for el in body.findall("Piece")
    ]


def read_index(filename, kind, piece=None, ghosts="keep"):
    """Read a parallel index: every piece merged (one region each) or one piece."""
    check_ghosts(ghosts)
    sources = parse_index(filename, kind)
    if not sources:
        return empty_mesh()

    def one(path):
        m = read_child(path)
        return drop_ghosts(m) if ghosts == "drop" else m

    if piece is not None:
        return one(sources[resolve_piece(piece, len(sources))])
    return merge_pieces(
        [one(s) for s in sources], [f"piece_{i}" for i in range(len(sources))]
    )


def share_field_data(out, meshes):
    """Give ``out`` the union of the inputs' ``field_data``, the first input winning.

    Field data belongs to the dataset, not to a piece: a parallel writer repeats it
    in every piece (or writes it once), so the copies are the same value. ``merge``
    cannot know that and renames any key present in more than one input to
    ``0:name``, ``1:name``, which would turn a ``TimeValue`` every piece carries
    into keys nothing looks up, and make the metadata (a union of names) disagree
    with the read.
    """
    out.field_data.clear()
    for m in meshes:
        for key, value in m.field_data.items():
            out.field_data.setdefault(key, value)


def merge_pieces(meshes, names):
    """Merge without welding and add one ``cell`` region per input, named ``names[i]``.

    A single input is returned as read: wrapping it in one all-cells region would
    hide the regions it already has (a ``.pvd`` step that is one ``.pvtu``).
    """
    if len(meshes) == 1:
        return meshes[0]
    out, point_maps, cell_maps = _merge_py(
        meshes,
        weld=False,
        atol=1e-8,
        source_tag=False,
        data_policy="fill",
        drop_duplicate_cells=False,
    )
    _merge_sets(meshes, out, point_maps, cell_maps)
    share_field_data(out, meshes)
    for name, cmap in zip(names, cell_maps):
        out.regions.append(Region(name, "cell", np.asarray(cmap, dtype=np.int64)))
    return out
