"""PLOT3D multi-block structured grids, ``.q`` solutions and ``.f`` function files.

The Python reference of ``formats/plot3d.cpp`` (roadmap 2.4.2). The file does not
state its variant, so every combination of ASCII / Fortran unformatted / raw,
single / multi-block, 2-D / 3-D, single / double precision, with / without
``iblank`` and either byte order is checked against the file's exact size (a Fortran
file's record sizes); one match is read, none or several is a ``ReadError``.
"""

from __future__ import annotations

import struct

import numpy as np

from .._exceptions import ReadError, WriteError
from .._fortran_records import fortran_record, fortran_records, sniff_fortran_records
from .._grid import lattice_cells
from .._mesh import Mesh
from .._regions import Region

DIMS_KEY = "plot3d:dims"
IBLANK_KEY = "plot3d:iblank"
_MAX_BLOCKS = 1000000
_ASCII_CHARS = set(b"0123456789 \t\r\n.-+eEdD,")


class _Stream:
    """A file as tokens (ASCII) or bytes (binary; a Fortran file's payloads joined)."""

    def __init__(self, path):
        with open(path, "rb") as f:
            data = f.read()
        self.ascii = bool(data) and set(data[:4096]) <= _ASCII_CHARS
        self.fortran = False
        self.big = False
        self.records = []
        if self.ascii:
            self.tokens = data.replace(b",", b" ").split()
            return
        self.bytes = data
        layout = sniff_fortran_records(data)
        if layout is not None:
            try:
                spans = fortran_records(data, layout, "plot3d")
            except ReadError:
                spans = None
            if spans is not None:
                marker = layout[0]
                covered = sum(size + 2 * marker for _, size in spans)
                if covered == len(data):
                    self.fortran = True
                    self.big = layout[1] == ">"
                    self.bytes = b"".join(data[o : o + s] for o, s in spans)
                    self.records = [s for _, s in spans]

    def size(self):
        return len(self.tokens) if self.ascii else len(self.bytes)

    def int_unit(self):
        return 1 if self.ascii else 4

    def real_unit(self, w):
        return 1 if self.ascii else w

    def int_at(self, at, big):
        if self.ascii:
            if at >= len(self.tokens):
                return None
            t = self.tokens[at]
            body = t[1:] if t[:1] in (b"-", b"+") else t
            if not body or not body.isdigit() or len(body) > 13:
                return None
            return int(t)
        if at + 4 > len(self.bytes):
            return None
        return struct.unpack_from(">i" if big else "<i", self.bytes, at)[0]

    def ints(self, at, count, big):
        if self.ascii:
            out = []
            for k in range(count):
                v = self.int_at(at + k, big)
                if v is None:
                    raise ReadError(
                        "meshio++: plot3d: an iblank value is not an integer"
                    )
                out.append(v)
            return np.asarray(out, dtype=np.int32)
        return np.frombuffer(
            self.bytes, dtype=">i4" if big else "<i4", count=count, offset=at
        ).astype(np.int32)

    def reals(self, at, count, w, big):
        if self.ascii:
            try:
                return np.asarray(
                    [
                        float(t.replace(b"d", b"e").replace(b"D", b"e"))
                        for t in self.tokens[at : at + count]
                    ],
                    dtype=np.float64,
                )
            except ValueError as exc:
                raise ReadError(f"meshio++: plot3d: {exc}") from None
        dtype = (">" if big else "<") + ("f4" if w == 4 else "f8")
        return np.frombuffer(self.bytes, dtype=dtype, count=count, offset=at).astype(
            np.float64
        )


def _bigs(stream, preferred=False):
    if stream.ascii or stream.fortran:
        return [stream.big]
    return [preferred, not preferred]


def _header(stream, multi, dim, big, extra):
    """``(dims rows of dim+extra ints, units consumed)`` or ``None``."""
    at = 0
    n = 1
    if multi:
        v = stream.int_at(0, big)
        if v is None or v < 1 or v > _MAX_BLOCKS:
            return None
        n = v
        at += stream.int_unit()
    rows = []
    for _ in range(n):
        row = [1, 1, 1, 0]
        for c in range(dim + extra):
            v = stream.int_at(at, big)
            if v is None or v < 1 or v > (1 << 31):
                return None
            if c < dim:
                row[c] = v
            else:
                row[3] = v
            at += stream.int_unit()
        rows.append(row)
    return rows, at


def _describe(shape):
    multi, dim, w, iblank, big = shape[:5]
    return (
        f"{'multi' if multi else 'single'}-block, {dim}-D, "
        f"{'single' if w == 4 else 'double'} precision, "
        f"{'with' if iblank else 'without'} iblank" + (", big-endian" if big else "")
    )


def _grid_shapes(stream, precision, multiblock, iblank, dimensions):
    out = []
    for big in _bigs(stream):
        for multi in (True, False):
            for dim in (3, 2):
                if multiblock is not None and multi != bool(multiblock):
                    continue
                if dimensions and dim != dimensions:
                    continue
                header = _header(stream, multi, dim, big, 0)
                if header is None:
                    continue
                rows, head = header
                dims = [r[:3] for r in rows]
                for w in (8,) if stream.ascii else (4, 8):
                    for ib in (False, True):
                        if precision and not stream.ascii and w != precision:
                            continue
                        if iblank is not None and ib != bool(iblank):
                            continue
                        total = head
                        records = ([4] if multi else []) + [4 * dim * len(dims)]
                        too_big = False
                        for d in dims:
                            np_ = d[0] * d[1] * d[2]
                            if np_ > stream.size():
                                too_big = True
                                break
                            units = np_ * (
                                dim * stream.real_unit(w)
                                + (stream.int_unit() if ib else 0)
                            )
                            total += units
                            records.append(units)
                        if too_big or total != stream.size():
                            continue
                        if stream.fortran and records != stream.records:
                            continue
                        out.append((multi, dim, w, ib, big, dims, head))
    return out


def _plausible(stream, shape):
    """Whether a reading's values look like a grid (finite, sane magnitude, no
    subnormals; small iblank values): the tie-break between readings the size alone
    cannot tell apart. The twin of ``p3d_plausible``."""
    _, dim, w, ib, big, dims, head = shape
    at = head
    for d in dims:
        np_ = d[0] * d[1] * d[2]
        try:
            v = stream.reals(at, np_ * dim, w, big)
        except ReadError:
            return False
        a = np.abs(v)
        with np.errstate(invalid="ignore"):
            bad = (a > 1e30) | ((a < 1e-30) & (a != 0))
        if not np.all(np.isfinite(v)) or bad.any():
            return False
        at += np_ * dim * stream.real_unit(w)
        if ib:
            try:
                b = stream.ints(at, np_, big)
            except ReadError:
                return False
            if (np.abs(b.astype(np.int64)) > _MAX_BLOCKS).any():
                return False
            at += np_ * stream.int_unit()
    return True


def _check_same(header_rows, dims, extra=False):
    if len(header_rows) != len(dims):
        return False
    nvar = header_rows[0][3] if extra else None
    for row, d in zip(header_rows, dims):
        if row[:3] != list(d):
            return False
        if extra and row[3] != nvar:
            return False
    return True


def _add_solution(mesh, shape, path):
    multi, dim, _, _, big0, dims, _ = shape
    stream = _Stream(path)
    nvar = dim + 2
    for big in _bigs(stream, big0):
        header = _header(stream, multi, dim, big, 0)
        if header is None or not _check_same(header[0], dims):
            continue
        rows, head = header
        for w in (8,) if stream.ascii else (4, 8):
            ru = stream.real_unit(w)
            total = head + sum((4 + d[0] * d[1] * d[2] * nvar) * ru for d in dims)
            if total != stream.size():
                continue
            rho, mom, energy, ref = [], [], [], []
            at = head
            for d in dims:
                np_ = d[0] * d[1] * d[2]
                ref.append(stream.reals(at, 4, w, big))
                at += 4 * ru
                values = stream.reals(at, np_ * nvar, w, big).reshape(nvar, np_)
                at += np_ * nvar * ru
                rho.append(values[0])
                m = np.zeros((np_, 3))
                m[:, : nvar - 2] = values[1 : nvar - 1].T
                mom.append(m)
                energy.append(values[nvar - 1])
            mesh.point_data["density"] = np.concatenate(rho)
            mesh.point_data["momentum"] = np.concatenate(mom)
            mesh.point_data["energy"] = np.concatenate(energy)
            mesh.field_data["plot3d:q_reference"] = np.asarray(ref, dtype=np.float64)
            return
    raise ReadError(
        f"meshio++: plot3d: '{path}' is not a PLOT3D solution for this grid (block "
        "count, dimensions or size disagree)"
    )


def _add_function(mesh, shape, path):
    multi, dim, _, _, big0, dims, _ = shape
    stream = _Stream(path)
    for big in _bigs(stream, big0):
        header = _header(stream, multi, dim, big, 1)
        if header is None or not _check_same(header[0], dims, extra=True):
            continue
        rows, head = header
        nv = rows[0][3]
        for w in (8,) if stream.ascii else (4, 8):
            ru = stream.real_unit(w)
            total = head + sum(d[0] * d[1] * d[2] * nv * ru for d in dims)
            if total != stream.size():
                continue
            parts = []
            at = head
            for d in dims:
                np_ = d[0] * d[1] * d[2]
                parts.append(stream.reals(at, np_ * nv, w, big).reshape(nv, np_).T)
                at += np_ * nv * ru
            f = np.concatenate(parts)
            mesh.point_data["plot3d:f"] = f[:, 0] if nv == 1 else f
            return
    raise ReadError(
        f"meshio++: plot3d: '{path}' is not a PLOT3D function file for this grid "
        "(block count, dimensions, variable count or size disagree)"
    )


def _option(value, name, table):
    if value is None:
        return None
    if value in table:
        return table[value]
    raise ReadError(
        f"meshio++: plot3d: {name} must be one of {sorted(map(str, table))}"
    )


def read(
    filename,
    solution=None,
    function=None,
    precision=None,
    multiblock=None,
    iblank=None,
    dimensions=None,
):
    """Read a PLOT3D grid; see the module docstring.

    :param solution: a ``.q`` file read as ``density``, ``momentum`` and ``energy``
        point data (its reference values to ``field_data["plot3d:q_reference"]``).
    :param function: a ``.f`` function file read as ``plot3d:f`` point data.
    :param precision: ``"single"`` or ``"double"`` (4 or 8), when the size cannot tell.
    :param multiblock: ``True`` or ``False``, likewise.
    :param iblank: ``True`` or ``False``, likewise.
    :param dimensions: 2 or 3, likewise.
    """
    w = _option(precision, "precision", {"single": 4, "double": 8, 4: 4, 8: 8})
    dims_opt = _option(dimensions, "dimensions", {2: 2, 3: 3})
    stream = _Stream(str(filename))
    shapes = _grid_shapes(stream, w, multiblock, iblank, dims_opt)
    if len(shapes) > 1:
        plausible = [s for s in shapes if _plausible(stream, s)]
        if plausible:
            shapes = plausible
    if not shapes:
        raise ReadError(
            f"meshio++: plot3d: '{filename}' is not a PLOT3D grid this reader can size "
            "(ASCII, Fortran unformatted or raw; single or multi-block; 2-D or 3-D; "
            "single or double precision; with or without iblank)"
        )
    if len(shapes) > 1:
        what = "".join(f"\n  - {_describe(s)}" for s in shapes)
        raise ReadError(
            f"meshio++: plot3d: '{filename}' reads {len(shapes)} ways; pick one with the "
            f"precision, multiblock, iblank or dimensions read option:{what}"
        )
    shape = shapes[0]
    multi, dim, w, ib, big, dims, head = shape
    points, blanks = [], []
    at = head
    for d in dims:
        np_ = d[0] * d[1] * d[2]
        coords = np.zeros((np_, 3))
        values = stream.reals(at, np_ * dim, w, big).reshape(dim, np_)
        coords[:, :dim] = values.T
        at += np_ * dim * stream.real_unit(w)
        points.append(coords)
        if ib:
            blanks.append(stream.ints(at, np_, big))
            at += np_ * stream.int_unit()
    cells, regions = [], []
    base = cell_base = 0
    for b, d in enumerate(dims):
        cell_type, conn = lattice_cells(*d)
        if len(conn):
            cells.append((cell_type, np.asarray(conn, dtype=np.int64) + base))
            regions.append(
                Region(
                    f"block_{b}", "cell", np.arange(cell_base, cell_base + len(conn))
                )
            )
            cell_base += len(conn)
        base += d[0] * d[1] * d[2]
    mesh = Mesh(
        np.concatenate(points),
        cells,
        point_data={IBLANK_KEY: np.concatenate(blanks)} if ib else {},
        field_data={DIMS_KEY: np.asarray(dims, dtype=np.int64)},
        regions=regions,
    )
    if solution is not None:
        _add_solution(mesh, shape, str(solution))
    if function is not None:
        _add_function(mesh, shape, str(function))
    return mesh


def write(filename, mesh):
    """The Fortran unformatted multi-block double-precision form, little-endian."""
    dims = mesh.field_data.get(DIMS_KEY)
    if dims is None:
        raise WriteError(
            "meshio++: plot3d: the mesh is not a set of PLOT3D blocks (it has no "
            "field_data['plot3d:dims']); read one from a PLOT3D file, or write a "
            "structured format such as vts"
        )
    dims = np.asarray(dims)
    if dims.size == 0 or dims.size % 3:
        raise WriteError(
            "meshio++: plot3d: field_data['plot3d:dims'] must be (nblocks, 3)"
        )
    dims = dims.reshape(-1, 3).astype(np.int64)
    if (dims < 1).any() or (dims > np.iinfo(np.int32).max).any():
        raise WriteError("meshio++: plot3d: a block dimension is out of range")
    npts = dims.prod(axis=1)
    total = int(npts.sum())
    if total != len(mesh.points):
        raise WriteError(
            f"meshio++: plot3d: the blocks in field_data['plot3d:dims'] hold {total} "
            f"points but the mesh has {len(mesh.points)}"
        )
    iblank = mesh.point_data.get(IBLANK_KEY)
    if iblank is not None and np.asarray(iblank).size != total:
        raise WriteError("meshio++: plot3d: plot3d:iblank needs one value per point")
    pts = np.zeros((total, 3))
    p = np.asarray(mesh.points, dtype=np.float64)
    pts[:, : min(3, p.shape[1])] = p[:, :3]
    layout = (4, "<")
    out = [
        fortran_record(struct.pack("<i", len(dims)), layout),
        fortran_record(dims.astype("<i4").tobytes(), layout),
    ]
    first = 0
    for n in npts:
        n = int(n)
        payload = pts[first : first + n].T.astype("<f8").tobytes()
        if iblank is not None:
            payload += (
                np.asarray(iblank)
                .reshape(-1)[first : first + n]
                .astype("<i4")
                .tobytes()
            )
        out.append(fortran_record(payload, layout))
        first += n
    with open(filename, "wb") as f:
        f.write(b"".join(out))
