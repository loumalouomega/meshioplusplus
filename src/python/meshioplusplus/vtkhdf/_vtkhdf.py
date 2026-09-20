"""VTKHDF (``.vtkhdf``), the pure-Python (h5py) reference.

Kitware's HDF5-based VTK format: one file holds geometry, fields, partitions and
time. This module reads and writes ``Type`` ``UnstructuredGrid``, ``PolyData``,
``PartitionedDataSetCollection`` and ``MultiBlockDataSet``; ``ImageData``,
``OverlappingAMR``, ``HyperTreeGrid``, ``Table``, ``RectilinearGrid`` and
``StructuredGrid`` are refused by name.

Layout facts this code depends on were measured against VTK 9.7's own
``vtkHDFWriter`` / ``vtkHDFReader`` rather than taken from the prose spec:

* ``Offsets`` holds ``NumberOfCells + 1`` *boundaries* per piece
  (``Offsets[0] == 0``), while the shared cell reconstruction wants one *end*
  offset per cell, so the reader drops each piece's leading zero.
* ``Connectivity`` and ``FaceConnectivity`` hold piece-local point ids and
  ``PolyhedronToFaces`` piece-local face ids; the reader renumbers them.
* A polyhedron's ``Connectivity`` row is its sorted unique node set, and
  ``PolyhedronOffsets`` has ``NumberOfCells + 1`` entries per piece covering
  *every* cell (a non-polyhedral cell is a zero-length span).
* A composite's ``VTKHDF`` and ``Assembly`` groups must track creation order
  (``vtkHDFReader`` aborts the process on an untracked one) and a
  ``PartitionedDataSetCollection`` block needs its ``Index`` attribute (without
  it the reader returns empty blocks).

The reader dispatches on *features* -- ``Type`` and which groups exist -- and
never on the declared minor version: files in the wild declare sloppy versions,
so the minor only produces a good error message.
"""

from __future__ import annotations

import numpy as np

from .. import _provenance
from .._common import num_nodes_per_cell, warn
from .._exceptions import ReadError, WriteError
from .._files import is_buffer
from .._mesh import CellBlock, Mesh
from .._regions import Region
from .._vtk_common import (
    meshio_to_vtk_order,
    meshio_to_vtk_type,
    vtk_cells_from_data,
    vtk_to_meshio_type,
)

ROOT = "VTKHDF"

#: The newest minor version this code knows about. A newer file is read with a
#: warning (the spec's own rule is that minors add things older readers may
#: ignore); a different *major* is refused.
KNOWN_MINOR = 8

#: The group attribute the provenance block rides in (the ``zarr`` shape).
PROVENANCE_ATTR = "meshioplusplus:provenance"

#: ``field_data`` key carrying a step's time, as every fanned-out step does.
TIME_KEY = "meshio:time"

UNSTRUCTURED_GRID = "UnstructuredGrid"
POLY_DATA = "PolyData"
PDC = "PartitionedDataSetCollection"
MULTIBLOCK = "MultiBlockDataSet"
_LEAF_TYPES = (UNSTRUCTURED_GRID, POLY_DATA)
_COMPOSITE_TYPES = (PDC, MULTIBLOCK)
_CATEGORIES = ("Vertices", "Lines", "Polygons", "Strips")

# VTK cell ids with no meshio++ type, named so the refusal says what it is.
_UNREADABLE_VTK_TYPES = {
    0: "empty cell",
    2: "poly-vertex",
    4: "poly-line",
    6: "triangle-strip",
    11: "voxel",
}
_VTK_POLYHEDRON = 42


# --------------------------------------------------------------------------- #
# small helpers                                                                #
# --------------------------------------------------------------------------- #
def _text_attr(obj, name):
    """A string attribute as ``str`` (fixed- or variable-length), or ``None``."""
    if name not in obj.attrs:
        return None
    v = obj.attrs[name]
    if isinstance(v, np.ndarray):
        v = v.ravel()[0] if v.size else ""
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace").rstrip("\x00 ")
    return str(v).rstrip("\x00 ")


def _ints(grp, name):
    """A small bookkeeping dataset (``NumberOf*``) as a flat int64 array."""
    if name not in grp:
        raise ReadError(f"meshio++: vtkhdf: '{grp.name}' has no '{name}' dataset")
    return np.asarray(grp[name][()], dtype=np.int64).ravel()


def _pick(value, fallback):
    return fallback if value is None else int(value)


def _check_version(grp):
    v = grp.attrs.get("Version")
    if v is None:
        raise ReadError(f"meshio++: vtkhdf: '{grp.name}' has no Version attribute")
    v = np.asarray(v).ravel()
    if v.size != 2:
        raise ReadError(
            f"meshio++: vtkhdf: Version must be two integers, got {v.tolist()}"
        )
    major, minor = int(v[0]), int(v[1])
    if major not in (1, 2):
        raise ReadError(
            f"meshio++: vtkhdf: unsupported VTKHDF version {major}.{minor}; "
            "this build reads 1.x and 2.x"
        )
    if major == 2 and minor > KNOWN_MINOR:
        warn(
            f"meshio++: vtkhdf: file declares VTKHDF {major}.{minor}, newer than "
            f"the {major}.{KNOWN_MINOR} this build knows; reading what it recognizes."
        )


def _friendly_type(vtk_type):
    name = _UNREADABLE_VTK_TYPES.get(vtk_type)
    return f"{name} (VTK type {vtk_type})" if name else f"VTK type {vtk_type}"


# --------------------------------------------------------------------------- #
# Steps                                                                        #
# --------------------------------------------------------------------------- #
class _Steps:
    """The ``Steps`` group: how many steps, their times, and where each lives."""

    def __init__(self, grp):
        self.group = grp["Steps"] if "Steps" in grp else None
        self.transient = self.group is not None
        self.values = None
        self.count = 1
        if self.transient:
            if "NSteps" not in self.group.attrs:
                raise ReadError(f"meshio++: vtkhdf: '{grp.name}/Steps' has no NSteps")
            self.count = int(self.group.attrs["NSteps"])
            if "Values" in self.group:
                self.values = np.asarray(
                    self.group["Values"][()], dtype=np.float64
                ).ravel()
                if self.values.size != self.count:
                    raise ReadError(
                        f"meshio++: vtkhdf: Steps/Values has {self.values.size} "
                        f"entries but NSteps is {self.count}"
                    )
            if self.count < 1:
                raise ReadError(
                    "meshio++: vtkhdf: NSteps is zero; the file has no data"
                )

    def resolve(self, time_step):
        """Step index for ``time_step`` (negative counts from the end)."""
        k = time_step + self.count if time_step < 0 else time_step
        if not 0 <= k < self.count:
            raise ReadError(
                f"meshio++: vtkhdf: time_step {time_step} is out of range; "
                f"the file has {self.count} step(s)"
            )
        return k

    def time(self, k):
        return None if self.values is None else float(self.values[k])

    def scalar(self, name, k, col=0):
        """``Steps/<name>[k]`` (or ``[k, col]``); ``None`` if absent or too short."""
        if not self.transient or name not in self.group:
            return None
        ds = self.group[name]
        if k >= ds.shape[0]:
            return None
        row = np.asarray(ds[k]).ravel()
        return int(row[min(col, row.size - 1)])

    def data_start(self, kind, name, k):
        """Start row of array ``name`` for step ``k`` (``None`` = static array)."""
        if not self.transient:
            return None
        table = f"{kind}DataOffsets"
        if table not in self.group or name not in self.group[table]:
            return None
        ds = self.group[table][name]
        return int(ds[k]) if k < ds.shape[0] else None

    def parts(self, k, total):
        """``(first piece, piece count)`` of step ``k`` among ``total`` pieces."""
        if not self.transient:
            return 0, total
        nparts = self.scalar("NumberOfParts", k)
        if nparts is None:
            nparts = (
                total // self.count
                if total >= self.count and total % self.count == 0
                else total
            )
        part0 = self.scalar("PartOffsets", k)
        if part0 is None:
            part0 = k * nparts if total >= self.count * nparts else 0
        if part0 < 0 or part0 + nparts > total:
            raise ReadError(
                f"meshio++: vtkhdf: step {k} names pieces [{part0}, {part0 + nparts}) "
                f"but the file holds {total}"
            )
        return part0, nparts


def _join_offsets(raw, counts, totals):
    """Turn per-piece ``count+1`` boundary runs into global per-cell END offsets.

    ``raw`` is the concatenation of every piece's boundaries (each starting at 0),
    ``counts`` the cells per piece and ``totals`` the ids spanned per piece.
    """
    counts = np.asarray(counts, dtype=np.int64)
    if counts.size == 0:
        return np.zeros(0, dtype=np.int64)
    raw = np.asarray(raw, dtype=np.int64)
    if raw.size != int(counts.sum()) + counts.size:
        raise ReadError(
            f"meshio++: vtkhdf: an Offsets dataset holds {raw.size} entries; the "
            f"piece counts require {int(counts.sum()) + counts.size}"
        )
    firsts = np.concatenate(([0], np.cumsum(counts + 1)[:-1]))
    keep = np.ones(raw.size, dtype=bool)
    keep[firsts] = False
    base = np.concatenate(([0], np.cumsum(np.asarray(totals, dtype=np.int64))[:-1]))
    return raw[keep] + np.repeat(base, counts)


def _starts(counts):
    """Exclusive prefix sums: where each piece starts inside a concatenation."""
    counts = np.asarray(counts, dtype=np.int64)
    return np.concatenate(([0], np.cumsum(counts)[:-1])) if counts.size else counts


# --------------------------------------------------------------------------- #
# cells                                                                        #
# --------------------------------------------------------------------------- #
def _polyhedron_run(a, b, poly, raw, blocks, cell_data, perm, at):
    """Decode polyhedra ``[a, b)`` into ``polyhedron<N>`` blocks (N = unique nodes).

    Bucketed by node count in first-seen order -- the convention the VTU, OpenFOAM,
    MED and CGNS readers all use -- so the run may be reordered; ``perm`` records
    where each file cell lands.
    """
    if poly is None:
        raise ReadError(
            "meshio++: vtkhdf: a cell has VTK type 42 (polyhedron) but the file "
            "carries no FaceConnectivity/FaceOffsets/PolyhedronToFaces/"
            "PolyhedronOffsets datasets"
        )
    face_conn, face_end, to_faces, poly_end = poly
    cells, counts = [], []
    for c in range(a, b):
        ps = int(poly_end[c - 1]) if c > 0 else 0
        pe = int(poly_end[c])
        faces = []
        for f in to_faces[ps:pe]:
            f = int(f)
            fs = int(face_end[f - 1]) if f > 0 else 0
            faces.append(np.array(face_conn[fs : int(face_end[f])], dtype=int))
        cells.append(faces)
        counts.append(
            int(np.unique(np.concatenate(faces)).size) if faces and pe > ps else 0
        )
    groups = {}
    for i, n_nodes in enumerate(counts):
        groups.setdefault(n_nodes, []).append(i)
    for n_nodes, idx in groups.items():
        blocks.append(CellBlock(f"polyhedron{n_nodes}", [cells[i] for i in idx]))
        idx = np.asarray(idx, dtype=np.int64)
        for name, arr in raw.items():
            cell_data[name].append(arr[a + idx])
        perm[a + idx] = at + np.arange(idx.size)
        at += idx.size
    return at


def _build_cells(conn, ends, types, poly, raw, lenient, piece_of_cell=None):
    """Cells from flat VTK arrays: ``(blocks, cell_data, perm)``.

    ``perm[i]`` is the block-major index file cell ``i`` ends up at, or ``-1`` if
    ``lenient`` dropped it. Polyhedra are bucketed by node count *per piece*
    (``piece_of_cell`` gives each cell's piece), so a piece's cells always occupy
    one contiguous range of the result and reading ``piece=k`` alone reproduces
    region ``k`` of the merged read exactly. Runs of any other type join across
    piece boundaries.
    """
    n = types.size
    blocks = []
    cell_data = {name: [] for name in raw}
    perm = np.empty(n, dtype=np.int64)
    if n == 0:
        return blocks, cell_data, perm
    breaks = np.flatnonzero(types[1:] != types[:-1]) + 1
    at = 0
    for a, b in zip(
        np.concatenate(([0], breaks)).tolist(), np.concatenate((breaks, [n])).tolist()
    ):
        vtk_type = int(types[a])
        if vtk_type == _VTK_POLYHEDRON:
            cuts = [a, b]
            if piece_of_cell is not None:
                cuts = [a, *(a + np.flatnonzero(np.diff(piece_of_cell[a:b])) + 1), b]
            for lo, hi in zip(cuts[:-1], cuts[1:]):
                at = _polyhedron_run(
                    int(lo), int(hi), poly, raw, blocks, cell_data, perm, at
                )
            continue
        meshio_type = vtk_to_meshio_type.get(vtk_type)
        if meshio_type is None or vtk_type == 0:
            if not lenient:
                raise ReadError(
                    f"meshio++: vtkhdf: {_friendly_type(vtk_type)} has no meshio++ "
                    "cell type; pass lenient=True to skip such cells"
                )
            warn(
                f"meshio++: vtkhdf: skipping {b - a} cell(s) of "
                f"{_friendly_type(vtk_type)}."
            )
            perm[a:b] = -1
            continue
        first = int(ends[a - 1]) if a > 0 else 0
        sub_ends = ends[a:b] - first
        fixed = num_nodes_per_cell.get(meshio_type)
        if fixed is not None and meshio_type != "polygon":
            sizes = np.diff(np.concatenate(([0], sub_ends)))
            if np.any(sizes != fixed):
                raise ReadError(
                    f"meshio++: vtkhdf: '{meshio_type}' cells must have {fixed} "
                    f"nodes, found {sorted(set(sizes.tolist()))}"
                )
        sub_blocks, sub_cd = vtk_cells_from_data(
            conn[first : int(ends[b - 1])],
            sub_ends,
            types[a:b],
            {name: arr[a:b] for name, arr in raw.items()},
        )
        blocks.extend(sub_blocks)
        for name, lst in sub_cd.items():
            cell_data[name].extend(lst)
        perm[a:b] = at + np.arange(b - a)
        at += b - a
    return blocks, cell_data, perm


# --------------------------------------------------------------------------- #
# reading one UnstructuredGrid / PolyData group                                #
# --------------------------------------------------------------------------- #
def _gather_topology(grp, steps, k, col, part0, first, cnt, pt_starts):
    """One topology's window: global ``(conn, end offsets, cells per piece, c0)``."""
    ncells_all = _ints(grp, "NumberOfCells")
    nconn_all = _ints(grp, "NumberOfConnectivityIds")
    lo, hi = part0 + first, part0 + first + cnt
    c_step = _pick(steps.scalar("CellOffsets", k, col), int(ncells_all[:part0].sum()))
    n_step = _pick(
        steps.scalar("ConnectivityIdOffsets", k, col), int(nconn_all[:part0].sum())
    )
    ncells, nconn = ncells_all[lo:hi], nconn_all[lo:hi]
    c0 = c_step + int(ncells_all[part0 : part0 + first].sum())
    n0 = n_step + int(nconn_all[part0 : part0 + first].sum())
    n_total, c_total = int(nconn.sum()), int(ncells.sum())
    conn = np.asarray(grp["Connectivity"][n0 : n0 + n_total], dtype=np.int64)
    offs = np.asarray(grp["Offsets"][c0 + lo : c0 + lo + c_total + cnt], dtype=np.int64)
    conn = conn + np.repeat(pt_starts, nconn)
    return conn, _join_offsets(offs, ncells, nconn), ncells, c0


def _gather_faces(g, steps, k, part0, first, cnt, pt_starts, ncells, c0):
    """The four polyhedron arrays for a window, renumbered to window-global ids."""
    if "FaceConnectivity" not in g:
        return None
    nfaces_all = _ints(g, "NumberOfFaces")
    nfconn_all = _ints(g, "NumberOfFaceConnectivityIds")
    npf_all = _ints(g, "NumberOfPolyhedronToFaceIds")
    lo, hi = part0 + first, part0 + first + cnt
    # Derived from the counts, not from Steps/FaceConnectivityOffsets & co: VTK's
    # writer leaves the last piece out of those tables for partitioned steps.
    f0 = int(nfaces_all[:lo].sum())
    fc0 = int(nfconn_all[:lo].sum())
    pf0 = int(npf_all[:lo].sum())
    nfaces, nfconn, npf = nfaces_all[lo:hi], nfconn_all[lo:hi], npf_all[lo:hi]
    face_conn = np.asarray(
        g["FaceConnectivity"][fc0 : fc0 + int(nfconn.sum())], dtype=np.int64
    ) + np.repeat(pt_starts, nfconn)
    face_end = _join_offsets(
        g["FaceOffsets"][f0 + lo : f0 + lo + int(nfaces.sum()) + cnt], nfaces, nfconn
    )
    to_faces = np.asarray(
        g["PolyhedronToFaces"][pf0 : pf0 + int(npf.sum())], dtype=np.int64
    ) + np.repeat(_starts(nfaces), npf)
    poly_end = _join_offsets(
        g["PolyhedronOffsets"][c0 + lo : c0 + lo + int(ncells.sum()) + cnt], ncells, npf
    )
    return face_conn, face_end, to_faces, poly_end


def _wanted(names, arrays, points_only):
    if points_only:
        return []
    if arrays is None:
        return list(names)
    keep = set(arrays)
    return [n for n in names if n in keep]


def _read_span(grp, kind, name, steps, k, in_step, count):
    """Rows of ``<Kind>Data/<name>`` for the window: ``in_step`` rows into step ``k``."""
    ds = grp[f"{kind}Data"][name]
    start = (steps.data_start(kind, name, k) or 0) + in_step
    arr = np.asarray(ds[start : start + count])
    if arr.shape[0] != count:
        raise ReadError(
            f"meshio++: vtkhdf: {kind}Data/{name} has {arr.shape[0]} rows in step "
            f"{k}; expected {count}"
        )
    return arr


def _read_field_data(g, steps, k, arrays, points_only):
    out = {}
    if "FieldData" not in g:
        return out
    fd = g["FieldData"]
    for name in _wanted(list(fd.keys()), arrays, points_only):
        ds = fd[name]
        if not hasattr(ds, "shape"):
            continue
        start = steps.data_start("Field", name, k)
        sizes = (
            steps.group["FieldDataSizes"][name]
            if steps.transient
            and "FieldDataSizes" in steps.group
            and name in steps.group["FieldDataSizes"]
            else None
        )
        if start is not None and sizes is not None and k < sizes.shape[0]:
            ncomp, ntuples = (int(x) for x in sizes[k])
            arr = np.asarray(ds[start : start + ntuples])
            if ncomp > 1 and arr.ndim == 1:
                arr = arr.reshape(ntuples, ncomp)
            out[name] = arr
        else:
            out[name] = np.asarray(ds[()])
    return out


def _read_dataset(g, kind, k, steps, only, points_only, arrays, lenient):
    """Read step ``k`` of an UnstructuredGrid/PolyData group.

    Returns ``(mesh, piece_of_cell, perm)``: the mesh, which piece each *file*
    cell came from, and the file-cell -> block-major-index permutation (``-1`` =
    dropped).
    """
    npts_all = _ints(g, "NumberOfPoints")
    part0, nparts = steps.parts(k, npts_all.size)
    first, cnt = (0, nparts) if only is None else (only, 1)
    if only is not None and not 0 <= only < nparts:  # callers resolve; belt and braces
        raise ReadError(f"meshio++: vtkhdf: piece {only} is out of range ({nparts})")
    lo, hi = part0 + first, part0 + first + cnt
    npts = npts_all[lo:hi]
    p_step = _pick(steps.scalar("PointOffsets", k), int(npts_all[:part0].sum()))
    in_step_pts = int(npts_all[part0 : part0 + first].sum())
    p0 = p_step + in_step_pts
    n_pts = int(npts.sum())
    points = np.asarray(g["Points"][p0 : p0 + n_pts])
    pt_starts = _starts(npts)

    if kind == UNSTRUCTURED_GRID:
        conn, ends, ncells, c0 = _gather_topology(
            g, steps, k, 0, part0, first, cnt, pt_starts
        )
        types = np.asarray(g["Types"][c0 : c0 + int(ncells.sum())], dtype=np.int64)
        poly = _gather_faces(g, steps, k, part0, first, cnt, pt_starts, ncells, c0)
        piece_of_cell = np.repeat(np.arange(cnt), ncells)
        c_in_step = int(_ints(g, "NumberOfCells")[part0 : part0 + first].sum())
    else:
        # Topology is stored per category (Vertices, Lines, Polygons, Strips), but
        # VTK's cell order -- and so the order of CellData -- is piece-major:
        # piece 0's cells in category order, then piece 1's (measured with
        # vtkHDFWriter). Gather each category, then interleave.
        sizes_l, starts_l, conn_l, types_l, piece_l, cat_l = [], [], [], [], [], []
        base = 0
        piece_cells = np.zeros(npts_all.size, dtype=np.int64)
        for col, cat in enumerate(_CATEGORIES):
            if cat not in g:
                continue
            sub = g[cat]
            piece_cells += _ints(sub, "NumberOfCells")
            conn_c, ends_c, ncells_c, _ = _gather_topology(
                sub, steps, k, col, part0, first, cnt, pt_starts
            )
            sizes = np.diff(np.concatenate(([0], ends_c)))
            if cat == "Vertices":
                t = np.where(sizes == 1, 1, 2)
            elif cat == "Lines":
                t = np.where(sizes == 2, 3, 4)
            elif cat == "Polygons":
                t = np.where(sizes == 3, 5, np.where(sizes == 4, 9, 7))
            else:
                t = np.full(sizes.shape, 6)
            sizes_l.append(sizes)
            starts_l.append(ends_c - sizes + base)
            conn_l.append(conn_c)
            types_l.append(t.astype(np.int64))
            piece_l.append(np.repeat(np.arange(cnt), ncells_c))
            cat_l.append(np.full(sizes.shape, col))
            base += conn_c.size
        if sizes_l:
            sizes = np.concatenate(sizes_l)
            row_starts = np.concatenate(starts_l)
            conn_all = np.concatenate(conn_l)
            order = np.lexsort((np.concatenate(cat_l), np.concatenate(piece_l)))
            sizes = sizes[order]
            ends = np.cumsum(sizes)
            gather = np.repeat(row_starts[order], sizes) + (
                np.arange(int(ends[-1]) if ends.size else 0)
                - np.repeat(ends - sizes, sizes)
            )
            conn = conn_all[gather]
            types = np.concatenate(types_l)[order]
            piece_of_cell = np.concatenate(piece_l)[order]
        else:
            conn = ends = types = np.zeros(0, dtype=np.int64)
            piece_of_cell = np.zeros(0, dtype=np.int64)
        poly = None
        c_in_step = int(piece_cells[part0 : part0 + first].sum())

    n_cells = types.size
    cell_raw = {}
    point_data = {}
    if not points_only:
        if "CellData" in g:
            for name in _wanted(list(g["CellData"].keys()), arrays, False):
                cell_raw[name] = _read_span(
                    g, "Cell", name, steps, k, c_in_step, n_cells
                )
        if "PointData" in g:
            for name in _wanted(list(g["PointData"].keys()), arrays, False):
                point_data[name] = _read_span(
                    g, "Point", name, steps, k, in_step_pts, n_pts
                )

    blocks, cell_data, perm = _build_cells(
        conn, ends, types, poly, cell_raw, lenient, piece_of_cell
    )
    field_data = _read_field_data(g, steps, k, arrays, points_only)
    if steps.transient and steps.time(k) is not None:
        field_data[TIME_KEY] = np.array([steps.time(k)])
    mesh = Mesh(
        points,
        blocks,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
    )
    return mesh, piece_of_cell, perm


# --------------------------------------------------------------------------- #
# composite datasets and the top-level reader                                  #
# --------------------------------------------------------------------------- #
def _shift(data, offset):
    """Cell block connectivity with every node id moved by ``offset``."""
    if offset == 0:
        return data
    if isinstance(data, np.ndarray):
        return data + offset
    out = []
    for cell in data:
        if len(cell) and isinstance(cell[0], (np.ndarray, list, tuple)):
            out.append([np.asarray(face) + offset for face in cell])
        else:
            out.append(np.asarray(cell) + offset)
    return out


def _can_merge(a_type, a, b_type, b):
    if a_type != b_type:
        return False
    if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        return a.ndim == 2 and b.ndim == 2 and a.shape[1] == b.shape[1]
    return isinstance(a, list) and isinstance(b, list)


def _concat_meshes(meshes):
    """Concatenate meshes without reordering, welding or regrouping by type.

    Cells keep their order (so region indices computed per input stay valid once
    offset); only *adjacent* blocks of one type are joined. Arrays present in some
    inputs but not all are dropped with a warning.
    """
    if not meshes:
        return Mesh(np.zeros((0, 3), dtype=np.float64), [])
    common_pd = set.intersection(*(set(m.point_data) for m in meshes))
    common_cd = set.intersection(*(set(m.cell_data) for m in meshes))
    dropped = (
        (
            set().union(*(set(m.point_data) for m in meshes))
            | set().union(*(set(m.cell_data) for m in meshes))
        )
        - common_pd
        - common_cd
    )
    if dropped:
        warn(
            "meshio++: vtkhdf: dropping array(s) not present in every piece: "
            + ", ".join(sorted(dropped))
        )
    points, blocks, cell_data = [], [], {name: [] for name in common_cd}
    point_data = {name: [] for name in common_pd}
    offset = 0
    for m in meshes:
        pts = np.asarray(m.points)
        points.append(pts)
        for name in common_pd:
            point_data[name].append(np.asarray(m.point_data[name]))
        for i, cb in enumerate(m.cells):
            data = _shift(cb.data, offset)
            cd = {name: np.asarray(m.cell_data[name][i]) for name in common_cd}
            if blocks and _can_merge(blocks[-1][0], blocks[-1][1], cb.type, data):
                t, prev = blocks[-1][0], blocks[-1][1]
                merged = (
                    np.concatenate([prev, data])
                    if isinstance(prev, np.ndarray)
                    else prev + data
                )
                blocks[-1] = (t, merged)
                for name in common_cd:
                    cell_data[name][-1] = np.concatenate(
                        [cell_data[name][-1], cd[name]]
                    )
            else:
                blocks.append((cb.type, data))
                for name in common_cd:
                    cell_data[name].append(cd[name])
        offset += len(pts)
    dtype = np.result_type(*(p.dtype for p in points))
    return Mesh(
        np.concatenate([p.astype(dtype, copy=False) for p in points]),
        [CellBlock(t, d) for t, d in blocks],
        point_data={n: np.concatenate(v) for n, v in point_data.items()},
        cell_data=cell_data,
    )


def _crt_keys(node):
    """Child names in creation order where the group tracks it (h5py's default)."""
    return list(node.keys())


def _composite_blocks(g, kind):
    """``[(name, block group)]`` in read order for a composite ``VTKHDF`` group."""
    import h5py

    entries = []
    if "Assembly" in g:

        def walk(node, prefix):
            for name in _crt_keys(node):
                link = node.get(name, getlink=True)
                if isinstance(link, h5py.SoftLink):
                    try:
                        target = node[name]
                    except KeyError as exc:
                        raise ReadError(
                            f"meshio++: vtkhdf: Assembly link '{prefix}{name}' points "
                            f"at '{link.path}', which does not exist"
                        ) from exc
                    entries.append((name, prefix + name, target))
                else:
                    walk(node[name], prefix + name + "/")

        walk(g["Assembly"], "")
    else:
        for name in _crt_keys(g):
            obj = g[name]
            if (
                name not in ("Assembly", "Steps", "FieldData")
                and isinstance(obj, h5py.Group)
                and "Type" in obj.attrs
            ):
                entries.append((name, name, obj))
    if len({e[0] for e in entries}) != len(entries):
        entries = [(path, path, grp) for _n, path, grp in entries]
    if kind == PDC and entries and all("Index" in e[2].attrs for e in entries):
        entries.sort(key=lambda e: int(np.asarray(e[2].attrs["Index"]).ravel()[0]))
    return [(name, grp) for name, _path, grp in entries]


def _check_no_time(g, where):
    steps = _Steps(g)
    if steps.count > 1:
        raise ReadError(
            f"meshio++: vtkhdf: {where} carries {steps.count} time steps; transient "
            "composite (PartitionedDataSetCollection / MultiBlockDataSet) datasets "
            "are not supported -- only UnstructuredGrid and PolyData files are "
            "transient"
        )
    return steps


def _read_composite(g, kind, time_step, piece, points_only, arrays, lenient):
    blocks = _composite_blocks(g, kind)
    root_steps = _check_no_time(g, ROOT)
    root_steps.resolve(time_step)  # 0 / -1 are the only valid steps of a static file
    plan = []
    for name, bg in blocks:
        btype = _text_attr(bg, "Type")
        if btype not in _LEAF_TYPES:
            raise ReadError(
                f"meshio++: vtkhdf: composite block '{name}' has Type "
                f"'{btype}'; only UnstructuredGrid and PolyData blocks are supported"
            )
        steps = _check_no_time(bg, f"block '{name}'")
        plan.append((name, bg, btype, steps, _ints(bg, "NumberOfPoints").size))

    fd = {}
    if "FieldData" in g:
        fd = _read_field_data(g, root_steps, 0, arrays, points_only)

    if piece is not None:
        total = sum(p[4] for p in plan)
        flat = piece + total if piece < 0 else piece
        if not 0 <= flat < total:
            raise ReadError(
                f"meshio++: vtkhdf: piece {piece} is out of range; the file has "
                f"{total} piece(s) across {len(plan)} block(s)"
            )
        for name, bg, btype, steps, nparts in plan:
            if flat < nparts:
                mesh, _, _ = _read_dataset(
                    bg, btype, 0, steps, flat, points_only, arrays, lenient
                )
                mesh.field_data.update(fd)
                return mesh
            flat -= nparts

    read = []
    for name, bg, btype, steps, nparts in plan:
        mesh, piece_of_cell, perm = _read_dataset(
            bg, btype, 0, steps, None, points_only, arrays, lenient
        )
        read.append((name, mesh, piece_of_cell, perm, nparts))
    out = _concat_meshes([r[1] for r in read])
    out.field_data.update(fd)
    if read:  # every block carries the same mesh-wide field_data; take the first
        out.field_data.update(read[0][1].field_data)
    cell_start = 0
    tag = 0
    for name, mesh, piece_of_cell, perm, nparts in read:
        n = sum(len(cb) for cb in mesh.cells)
        if nparts == 1:
            pieces = [(name, np.arange(n, dtype=np.int64))]
        else:
            pieces = []
            for j in range(nparts):
                idx = perm[piece_of_cell == j]
                pieces.append((f"{name}/piece_{j}", idx[idx >= 0]))
        for region_name, idx in pieces:
            # tag = the piece's position in the file: what a writer sorts by to
            # restore the order (a region list alone cannot carry it across the
            # C++ core, which stores regions sorted by name).
            out.regions.append(
                Region(
                    region_name,
                    "cell",
                    (idx + cell_start).astype(np.int64),
                    tag=tag,
                )
            )
            tag += 1
        cell_start += n
    return out


def read(
    filename,
    points_only=False,
    arrays=None,
    time_step=0,
    piece=None,
    lenient=False,
):
    """Read a VTKHDF file.

    :param points_only: geometry and connectivity, no data arrays at all.
    :param arrays: restrict point/cell/field data to these names (``None`` = all,
        an empty list = none).
    :param time_step: which step of a transient file (negative counts from the
        end); a static file has exactly one.
    :param piece: keep only this piece of a partitioned file (negative counts from
        the end). ``None`` (default) merges every piece into one mesh and attaches
        one cell region per piece, named ``piece_<i>`` (``<block>`` for a
        composite's blocks).
    :param lenient: skip cells with no meshio++ type (poly-vertex, poly-line,
        triangle strips, ...) with a warning instead of raising.
    """
    import h5py

    try:
        f = h5py.File(filename if is_buffer(filename, "r") else str(filename), "r")
    except OSError as exc:
        raise ReadError(
            f"meshio++: vtkhdf: cannot open '{filename}' as an HDF5 file ({exc})"
        ) from exc
    with f:
        if ROOT not in f:
            raise ReadError(
                f"meshio++: vtkhdf: '{filename}' is an HDF5 file with no /{ROOT} "
                "group; if it holds another HDF5 format, pass that format explicitly"
            )
        g = f[ROOT]
        kind = _text_attr(g, "Type")
        _check_version(g)
        if kind in _COMPOSITE_TYPES:
            return _read_composite(
                g, kind, time_step, piece, points_only, arrays, lenient
            )
        if kind not in _LEAF_TYPES:
            raise ReadError(
                f"meshio++: vtkhdf: Type '{kind}' is not supported; this build reads "
                f"{', '.join(_LEAF_TYPES + _COMPOSITE_TYPES)}"
            )
        steps = _Steps(g)
        k = steps.resolve(time_step)
        nparts = steps.parts(k, _ints(g, "NumberOfPoints").size)[1]
        only = None
        if piece is not None:
            only = piece + nparts if piece < 0 else piece
            if not 0 <= only < nparts:
                raise ReadError(
                    f"meshio++: vtkhdf: piece {piece} is out of range; "
                    f"the file has {nparts} piece(s)"
                )
        mesh, piece_of_cell, perm = _read_dataset(
            g, kind, k, steps, only, points_only, arrays, lenient
        )
        if steps.values is not None:
            # The side channel `read_metadata` uses when it cannot ask the C++
            # core: every step's time, not just the one this mesh holds.
            mesh.time_values = [float(v) for v in steps.values]
        if only is None and nparts > 1:
            for j in range(nparts):
                idx = perm[piece_of_cell == j]
                mesh.regions.append(
                    Region(f"piece_{j}", "cell", idx[idx >= 0].astype(np.int64), tag=j)
                )
        return mesh


# --------------------------------------------------------------------------- #
# writing                                                                      #
# --------------------------------------------------------------------------- #
def _has_polyhedra(mesh):
    return any(cb.type.startswith("polyhedron") for cb in mesh.cells)


def _needed_version(kind, has_poly):
    """The oldest VTKHDF version that can express the file being written."""
    if has_poly:
        return (2, 5)
    if kind in _COMPOSITE_TYPES:
        return (2, 1)
    if kind == POLY_DATA:
        return (2, 0)
    return (1, 0)


def _resolve_version(version, kind, has_poly):
    """``(major, minor)`` to write: the feature-derived minimum unless pinned."""
    needed = _needed_version(kind, has_poly)
    if version is None:
        # Never a fixed maximum: the widest reader compatibility is the oldest
        # version that covers the features used. UnstructuredGrid is written as
        # 2.0 (the version that added time) rather than 1.0 so a file later
        # extended with steps needs no rewrite of its Version.
        return max(needed, (2, 0))
    try:
        v = tuple(int(x) for x in version)
    except (TypeError, ValueError):
        v = ()
    if len(v) != 2 or v[0] not in (1, 2) or v[1] < 0:
        raise WriteError(
            f"meshio++: vtkhdf: version must be (major, minor) with major 1 or 2, "
            f"got {version!r}"
        )
    if v[0] == 2 and v[1] > KNOWN_MINOR:
        raise WriteError(
            f"meshio++: vtkhdf: version {v[0]}.{v[1]} is newer than the "
            f"2.{KNOWN_MINOR} this build knows"
        )
    if v < needed:
        feature = "polyhedral cells" if has_poly else f"Type '{kind}'"
        raise WriteError(
            f"meshio++: vtkhdf: {feature} needs VTKHDF {needed[0]}.{needed[1]} or "
            f"newer; version {v[0]}.{v[1]} cannot express it"
        )
    return v


def _write_type_attr(grp, value):
    """``Type`` as a fixed-length ASCII string (variable-length UTF-8 has tripped readers)."""
    import h5py

    grp.attrs.create(
        "Type", np.bytes_(value), dtype=h5py.string_dtype("ascii", len(value))
    )


#: Below this, chunking plus a deflate filter costs more in HDF5 metadata than it saves.
_GZIP_MIN_BYTES = 4096


def _ds(grp, name, data, gzip):
    """Create dataset ``name``, deflated only when big enough for that to pay."""
    data = np.asarray(data)
    if gzip is not None and data.nbytes >= _GZIP_MIN_BYTES:
        grp.create_dataset(name, data=data, compression="gzip", compression_opts=gzip)
    else:
        grp.create_dataset(name, data=data)


def _int64(values):
    return np.asarray(values, dtype=np.int64)


def _data_array(name, arr, what):
    """A ``point_data``/``cell_data``/``field_data`` array in a storable form."""
    if "/" in name or not name:
        raise WriteError(
            f"meshio++: vtkhdf: {what} name '{name}' cannot be an HDF5 dataset name "
            "(empty or containing '/')"
        )
    a = np.asarray(arr)
    if a.dtype == bool:
        a = a.astype(np.uint8)
    if a.dtype.kind not in "iuf":
        raise WriteError(
            f"meshio++: vtkhdf: {what} '{name}' has dtype {a.dtype}; only integer "
            "and floating-point arrays can be written"
        )
    if a.dtype.byteorder == ">":
        a = a.astype(a.dtype.newbyteorder("<"))
    if a.ndim == 0:
        a = a.reshape(1)
    elif a.ndim > 2:
        a = a.reshape(a.shape[0], -1)
    return a


def _raw_cell_data(mesh, order=None):
    """``cell_data`` concatenated over blocks (in ``order`` if given)."""
    out = {}
    blocks = range(len(mesh.cells)) if order is None else order
    for name, lst in mesh.cell_data.items():
        parts = []
        for bi in blocks:
            if lst[bi] is None:
                raise WriteError(
                    f"meshio++: vtkhdf: cell_data '{name}' has no data for cell "
                    f"block {bi} ('{mesh.cells[bi].type}')"
                )
            parts.append(np.asarray(lst[bi]))
        out[name] = np.concatenate(parts) if parts else np.zeros(0, dtype=np.float64)
    return out


def _write_data_groups(grp, mesh, raw_cell, gzip, field_data=True):
    pd = grp.create_group("PointData")
    for name, arr in mesh.point_data.items():
        _ds(pd, name, _data_array(name, arr, "point_data"), gzip)
    cd = grp.create_group("CellData")
    for name, arr in raw_cell.items():
        _ds(cd, name, _data_array(name, arr, "cell_data"), gzip)
    fd = grp.create_group("FieldData")
    if field_data:
        for name, arr in mesh.field_data.items():
            try:
                a = _data_array(name, arr, "field_data")
            except WriteError as exc:
                warn(f"{exc}; skipping it.")
                continue
            _ds(fd, name, a, gzip)


def _points3(mesh):
    pts = np.asarray(mesh.points)
    if pts.size == 0:
        return np.zeros((0, 3), dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] > 3:
        raise WriteError(
            f"meshio++: vtkhdf: points must be (n, 1..3), got shape {pts.shape}"
        )
    if pts.dtype.kind != "f":
        pts = pts.astype(np.float64)
    if pts.shape[1] < 3:
        pts = np.column_stack(
            [pts, np.zeros((pts.shape[0], 3 - pts.shape[1]), dtype=pts.dtype)]
        )
    return pts


def _ug_arrays(mesh):
    """Flat VTK arrays for an UnstructuredGrid: conn, sizes, types, polyhedron data."""
    conn, sizes, types = [], [], []
    face_rows, face_sizes, faces_per_cell = [], [], []
    has_poly = False
    for cb in mesh.cells:
        n = len(cb)
        t = cb.type
        if t.startswith("polyhedron"):
            has_poly = True
            for cell in cb.data:
                faces = [np.asarray(f, dtype=np.int64).ravel() for f in cell]
                nodes = (
                    np.unique(np.concatenate(faces))
                    if faces
                    else np.zeros(0, dtype=np.int64)
                )
                conn.append(nodes)
                sizes.append(np.array([nodes.size]))
                face_rows.extend(faces)
                face_sizes.extend(f.size for f in faces)
                faces_per_cell.append(len(faces))
            types.append(np.full(n, _VTK_POLYHEDRON, dtype=np.uint8))
            continue
        faces_per_cell.extend([0] * n)
        if t.startswith("polygon"):
            rows = [np.asarray(r, dtype=np.int64).ravel() for r in cb.data]
            conn.extend(rows)
            sizes.append(np.array([r.size for r in rows], dtype=np.int64))
            types.append(np.full(n, meshio_to_vtk_type["polygon"], dtype=np.uint8))
            continue
        vtk_type = meshio_to_vtk_type.get(t)
        if vtk_type is None:
            raise WriteError(
                f"meshio++: vtkhdf: cell type '{t}' has no VTK cell type id"
            )
        data = np.asarray(cb.data)
        order = meshio_to_vtk_order(t)
        if order is not None:
            data = data[:, order]
        conn.append(data.reshape(-1).astype(np.int64, copy=False))
        sizes.append(np.full(n, data.shape[1] if data.ndim == 2 else 0, dtype=np.int64))
        types.append(np.full(n, vtk_type, dtype=np.uint8))
    out = {
        "Connectivity": _int64(np.concatenate(conn)) if conn else _int64([]),
        "sizes": _int64(np.concatenate(sizes)) if sizes else _int64([]),
        "Types": np.concatenate(types) if types else np.zeros(0, dtype=np.uint8),
    }
    out["Offsets"] = _int64(np.concatenate(([0], np.cumsum(out["sizes"]))))
    if has_poly:
        out["FaceConnectivity"] = (
            _int64(np.concatenate(face_rows)) if face_rows else _int64([])
        )
        out["FaceOffsets"] = _int64(np.concatenate(([0], np.cumsum(face_sizes))))
        out["PolyhedronToFaces"] = _int64(np.arange(len(face_rows)))
        out["PolyhedronOffsets"] = _int64(
            np.concatenate(([0], np.cumsum(faces_per_cell)))
        )
    return out


def _write_unstructured(grp, mesh, gzip):
    a = _ug_arrays(mesh)
    pts = _points3(mesh)
    _ds(grp, "NumberOfPoints", _int64([len(pts)]), None)
    _ds(grp, "NumberOfCells", _int64([a["Types"].size]), None)
    _ds(grp, "NumberOfConnectivityIds", _int64([a["Connectivity"].size]), None)
    _ds(grp, "Points", pts, gzip)
    _ds(grp, "Connectivity", a["Connectivity"], gzip)
    _ds(grp, "Offsets", a["Offsets"], gzip)
    _ds(grp, "Types", a["Types"], gzip)
    if "FaceConnectivity" in a:
        _ds(grp, "NumberOfFaces", _int64([a["FaceOffsets"].size - 1]), None)
        _ds(
            grp,
            "NumberOfFaceConnectivityIds",
            _int64([a["FaceConnectivity"].size]),
            None,
        )
        _ds(
            grp,
            "NumberOfPolyhedronToFaceIds",
            _int64([a["PolyhedronToFaces"].size]),
            None,
        )
        _ds(grp, "FaceConnectivity", a["FaceConnectivity"], gzip)
        _ds(grp, "FaceOffsets", a["FaceOffsets"], gzip)
        _ds(grp, "PolyhedronToFaces", a["PolyhedronToFaces"], gzip)
        _ds(grp, "PolyhedronOffsets", a["PolyhedronOffsets"], gzip)
    return _raw_cell_data(mesh)


def _pd_category(cell_type):
    if cell_type == "vertex":
        return 0
    if cell_type == "line":
        return 1
    if cell_type in ("triangle", "quad") or cell_type.startswith("polygon"):
        return 2
    raise WriteError(
        f"meshio++: vtkhdf: PolyData cannot hold '{cell_type}' cells (only vertex, "
        "line, triangle, quad and polygon); write dataset_type='UnstructuredGrid' "
        "to keep them"
    )


def _write_polydata(grp, mesh, gzip):
    kinds = [_pd_category(cb.type) for cb in mesh.cells]
    # VTK's canonical PolyData cell order is Vertices, Lines, Polygons, Strips, so
    # the blocks are regrouped (stably) and cell_data follows the same order.
    order = [bi for want in (0, 1, 2) for bi, kind in enumerate(kinds) if kind == want]
    pts = _points3(mesh)
    _ds(grp, "NumberOfPoints", _int64([len(pts)]), None)
    _ds(grp, "Points", pts, gzip)
    rows = {0: [], 1: [], 2: []}
    for bi in order:
        cb = mesh.cells[bi]
        for r in cb.data:
            r = np.asarray(r, dtype=np.int64).ravel()
            if kinds[bi] == 0 and r.size != 1:
                raise WriteError("meshio++: vtkhdf: vertex cells must have one node")
            if kinds[bi] == 1 and r.size != 2:
                raise WriteError("meshio++: vtkhdf: line cells must have two nodes")
            rows[kinds[bi]].append(r)
    for col, cat in enumerate(_CATEGORIES):
        sub = grp.create_group(cat)
        cat_rows = rows.get(col, [])
        conn = _int64(np.concatenate(cat_rows)) if cat_rows else _int64([])
        offs = _int64(np.concatenate(([0], np.cumsum([r.size for r in cat_rows]))))
        _ds(sub, "NumberOfCells", _int64([len(cat_rows)]), None)
        _ds(sub, "NumberOfConnectivityIds", _int64([conn.size]), None)
        _ds(sub, "Connectivity", conn, gzip)
        _ds(sub, "Offsets", offs, gzip)
    return _raw_cell_data(mesh, order)


def _write_leaf(grp, mesh, kind, gzip, version):
    grp.attrs.create("Version", np.asarray(version, dtype="<i8"))
    _write_type_attr(grp, kind)
    raw = (
        _write_polydata(grp, mesh, gzip)
        if kind == POLY_DATA
        else _write_unstructured(grp, mesh, gzip)
    )
    _write_data_groups(grp, mesh, raw, gzip)


# ---- composite: carve the mesh into blocks -------------------------------- #
def _node_ids(data):
    if isinstance(data, np.ndarray):
        return data.ravel()
    parts = []
    for cell in data:
        if len(cell) and isinstance(cell[0], (np.ndarray, list, tuple)):
            parts.extend(np.asarray(f).ravel() for f in cell)
        else:
            parts.append(np.asarray(cell).ravel())
    return np.concatenate(parts) if parts else np.zeros(0, dtype=np.int64)


def _remap(data, used):
    if isinstance(data, np.ndarray):
        return np.searchsorted(used, data)
    out = []
    for cell in data:
        if len(cell) and isinstance(cell[0], (np.ndarray, list, tuple)):
            out.append([np.searchsorted(used, np.asarray(f)) for f in cell])
        else:
            out.append(np.searchsorted(used, np.asarray(cell)))
    return out


def _extract(mesh, idx):
    """The sub-mesh of cells ``idx`` (sorted, block-major), pruned to its own points."""
    idx = np.asarray(idx, dtype=np.int64)
    bounds = np.cumsum([0] + [len(cb) for cb in mesh.cells])
    blocks, cell_data = [], {name: [] for name in mesh.cell_data}
    for bi, cb in enumerate(mesh.cells):
        local = idx[(idx >= bounds[bi]) & (idx < bounds[bi + 1])] - bounds[bi]
        if local.size == 0:
            continue
        if isinstance(cb.data, np.ndarray):
            data = cb.data[local]
        else:
            data = [cb.data[i] for i in local]
        blocks.append((cb.type, data))
        for name, lst in mesh.cell_data.items():
            if lst[bi] is None:
                raise WriteError(
                    f"meshio++: vtkhdf: cell_data '{name}' has no data for cell "
                    f"block {bi} ('{cb.type}')"
                )
            cell_data[name].append(np.asarray(lst[bi])[local])
    ids = [_node_ids(d) for _t, d in blocks]
    used = np.unique(np.concatenate(ids)) if ids else np.zeros(0, dtype=np.int64)
    return Mesh(
        np.asarray(mesh.points)[used],
        [(t, _remap(d, used)) for t, d in blocks],
        point_data={n: np.asarray(v)[used] for n, v in mesh.point_data.items()},
        cell_data=cell_data,
    )


def _carve(mesh):
    """``[(block name, sub-mesh)]`` for a composite write."""
    total = sum(len(cb) for cb in mesh.cells)
    if total == 0:
        raise WriteError(
            "meshio++: vtkhdf: a composite dataset needs at least one cell to "
            "carve into blocks"
        )
    # Blocks are written in (tag, name) order: a read composite tags each region with
    # its block position, which restores the order; regions without tags fall back to
    # name order. The C++ core stores regions sorted by name, so this is the one
    # rule both engines can honour.
    regions = sorted(
        (r for r in mesh.regions if r.kind == "cell"), key=lambda r: (r.tag, r.name)
    )
    if regions:
        cover = np.concatenate([np.asarray(r.entries, dtype=np.int64) for r in regions])
        if cover.size == total and np.unique(cover).size == total:
            names = [r.name or f"block_{i}" for i, r in enumerate(regions)]
            if len(set(names)) != len(names) or any("/" in n for n in names):
                raise WriteError(
                    "meshio++: vtkhdf: cell region names must be unique and free of "
                    "'/' to name composite blocks"
                )
            return [
                (name, _extract(mesh, np.sort(np.asarray(r.entries, dtype=np.int64))))
                for name, r in zip(names, regions)
            ]
        warn(
            "meshio++: vtkhdf: cell regions overlap or do not cover every cell; "
            "writing one block per cell block instead."
        )
    bounds = np.cumsum([0] + [len(cb) for cb in mesh.cells])
    return [
        (f"block_{i}", _extract(mesh, np.arange(bounds[i], bounds[i + 1])))
        for i in range(len(mesh.cells))
        if len(mesh.cells[i])
    ]


def _write_composite(root, mesh, kind, gzip, version):
    import h5py

    pieces = _carve(mesh)
    for i, (name, sub) in enumerate(pieces):
        # A composite's root group may hold only blocks and the Assembly (vtkHDFReader
        # reads every other root group as a block and then fails), so the mesh-wide
        # field_data rides on each block's own FieldData group instead.
        sub.field_data = dict(mesh.field_data)
        bg = root.create_group(name, track_order=True)
        if kind == PDC:
            bg.attrs.create("Index", np.int64(i))
        _write_leaf(bg, sub, UNSTRUCTURED_GRID, gzip, version)
    asm = root.create_group("Assembly", track_order=True)
    for name, _sub in pieces:
        target = f"/{ROOT}/{name}"
        if kind == PDC:
            node = asm.create_group(name, track_order=True)
            node[name] = h5py.SoftLink(target)
        else:
            asm[name] = h5py.SoftLink(target)


def write(
    filename,
    mesh,
    compression="gzip",
    compression_opts=4,
    dataset_type=UNSTRUCTURED_GRID,
    version=None,
):
    """Write a static VTKHDF file.

    :param compression: ``"gzip"`` or ``None``.
    :param compression_opts: gzip level 0-9.
    :param dataset_type: ``"UnstructuredGrid"`` (default), ``"PolyData"``,
        ``"PartitionedDataSetCollection"`` or ``"MultiBlockDataSet"``. PolyData is
        never chosen automatically and refuses any cell that is not a vertex,
        line, triangle, quad or polygon. The composite types write one block per
        cell region (when the regions partition the cells) or per cell block.
    :param version: ``(major, minor)`` to declare; ``None`` writes the oldest
        version that covers the features used (2.0 for a plain mesh, 2.1 for a
        composite, 2.5 with polyhedra). A pinned version too old for the content
        raises.
    """
    import h5py

    if dataset_type not in _LEAF_TYPES + _COMPOSITE_TYPES:
        raise WriteError(
            f"meshio++: vtkhdf: unknown dataset_type '{dataset_type}'; expected one "
            f"of {', '.join(_LEAF_TYPES + _COMPOSITE_TYPES)}"
        )
    if compression not in (None, "gzip"):
        raise WriteError(
            f"meshio++: vtkhdf: compression must be 'gzip' or None, got {compression!r}"
        )
    gzip = None
    if compression == "gzip":
        gzip = 4 if compression_opts is None else int(compression_opts)
        if not 0 <= gzip <= 9:
            raise WriteError("meshio++: vtkhdf: gzip level must be 0-9")
    ver = _resolve_version(version, dataset_type, _has_polyhedra(mesh))

    with h5py.File(filename if is_buffer(filename, "w") else str(filename), "w") as f:
        root = f.create_group(ROOT, track_order=True)
        provenance = "\n".join(_provenance.lines(_provenance.SlotTier.BLOCK))
        if dataset_type in _COMPOSITE_TYPES:
            root.attrs.create("Version", np.asarray(ver, dtype="<i8"))
            _write_type_attr(root, dataset_type)
            _write_composite(root, mesh, dataset_type, gzip, ver)
        else:
            _write_leaf(root, mesh, dataset_type, gzip, ver)
        if provenance:
            root.attrs[PROVENANCE_ATTR] = provenance
