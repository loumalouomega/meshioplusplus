"""
I/O for Medit solution files (``.sol`` ascii, ``.solb`` binary), the fields
that sit beside a Medit ``.mesh``/``.meshb`` (libMeshb's ``SolAt*`` keywords,
<https://github.com/LoicMarechal/libMeshb>) and that MMG reads and writes as
its size maps and metrics.

A solution section holds one row per vertex (``SolAtVertices``) or per element
of one type (``SolAtTriangles``, ``SolAtTetrahedra``, ...), each row the
concatenation of the section's fields. A field is a scalar (type 1), a vector
of ``dim`` components (type 2), a symmetric matrix of ``dim*(dim+1)/2``
components (type 3) or a full ``dim*dim`` matrix (type 4). A symmetric matrix
keeps the file's order, the lower triangle row by row: ``m11 m21 m22`` in 2-D,
``m11 m21 m22 m31 m32 m33`` in 3-D.

Read standalone, the file is a geometry-less :class:`Mesh` (``points`` with no
columns, one row per vertex) carrying the vertex fields as
``point_data["medit:sol:<i>"]`` and the dimension as
``field_data["medit:dimension"]``; element fields need the mesh they belong to,
and ``medit.read(path, sol=...)`` attaches both kinds onto it.
"""

import numpy as np

from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import open_file
from .._mesh import Mesh

__all__ = ["read", "write", "read_sections", "attach", "SOL_KEYWORDS"]

# keyword -> (libMeshb code, the cell type it is attached to; None: vertices)
SOL_KEYWORDS = {
    "SolAtVertices": (62, None),
    "SolAtEdges": (63, "line"),
    "SolAtTriangles": (64, "triangle"),
    "SolAtQuadrilaterals": (65, "quad"),
    "SolAtTetrahedra": (66, "tetra"),
    "SolAtPrisms": (67, "wedge"),
    "SolAtHexahedra": (68, "hexahedron"),
    "SolAtPyramids": (26, "pyramid"),
}
_SOL_CODES = {code: (kw, cell_type) for kw, (code, cell_type) in SOL_KEYWORDS.items()}
_RESERVED = (0, 2, 11, 12, 51)


def components(sol_type, dim):
    """Values per entity of a field of libMeshb solution type ``sol_type``."""
    if sol_type == 1:
        return 1
    if sol_type == 2:
        return dim
    if sol_type == 3:
        return dim * (dim + 1) // 2
    if sol_type == 4:
        return dim * dim
    raise ReadError(f"Medit: unknown solution type {sol_type}")


def sol_type_of(num_components, dim):
    """The solution type that holds ``num_components`` values in ``dim``-D, or None."""
    for sol_type in (1, 2, 3, 4):
        if components(sol_type, dim) == num_components:
            return sol_type
    return None


class Section:
    """One ``SolAt*`` section: its keyword, field types and ``(n, total)`` values."""

    def __init__(self, keyword, types, values):
        self.keyword = keyword
        self.cell_type = SOL_KEYWORDS[keyword][1]
        self.types = list(types)
        self.values = values

    def fields(self, dim):
        """The section's fields, ``(n,)`` for a scalar, else ``(n, components)``."""
        out = []
        offset = 0
        for sol_type in self.types:
            c = components(sol_type, dim)
            field = self.values[:, offset : offset + c]
            out.append(field[:, 0].copy() if c == 1 else field.copy())
            offset += c
        return out


def _check_types(count, types, dim):
    if count < 1:
        raise ReadError("Medit: a solution section declares no field")
    return sum(components(int(t), dim) for t in types)


def read_sections(filename):
    """``(dim, [Section, ...])`` of a ``.sol`` (or, ending in ``b``, ``.solb``)."""
    if str(filename)[-1] == "b":
        with open_file(filename, "rb") as f:
            return _read_binary_sections(f.read())
    with open_file(filename) as f:
        return _read_ascii_sections(f.read())


def _read_ascii_sections(text):
    tokens = []
    for line in text.splitlines():
        tokens.extend(line.split("#", 1)[0].split())
    at = 0

    def take():
        nonlocal at
        if at >= len(tokens):
            raise ReadError("Medit: the solution file is truncated")
        at += 1
        return tokens[at - 1]

    def take_int():
        token = take()
        try:
            return int(token)
        except ValueError:
            raise ReadError(f"Medit: expected an integer, not '{token}'") from None

    dim = 0
    dtype = np.float64
    sections = []
    while at < len(tokens):
        keyword = take()
        if keyword == "MeshVersionFormatted":
            dtype = np.float32 if take_int() <= 1 else np.float64
        elif keyword == "Dimension":
            dim = take_int()
            if dim not in (2, 3):
                raise ReadError(f"Medit: invalid solution dimension {dim}")
        elif keyword in SOL_KEYWORDS:
            if dim == 0:
                raise ReadError(f"Medit: {keyword} before Dimension")
            count = take_int()
            num_types = take_int()
            if count < 0 or num_types > len(tokens) - at:
                raise ReadError(f"Medit: {keyword} has an impossible size")
            types = [take_int() for _ in range(max(num_types, 0))]
            total = _check_types(num_types, types, dim)
            if count * total > len(tokens) - at:
                raise ReadError(f"Medit: {keyword} is truncated")
            try:
                values = np.array(tokens[at : at + count * total], dtype=np.float64)
            except ValueError:
                raise ReadError(f"Medit: {keyword} holds a non-number") from None
            at += count * total
            sections.append(
                Section(keyword, types, values.astype(dtype).reshape(count, total))
            )
        elif keyword in ("Time", "Iterations"):
            take()
        elif keyword == "End":
            break
        else:
            raise ReadError(f"Medit: unknown solution keyword '{keyword}'")
    if dim == 0:
        raise ReadError("Medit: expected Dimension")
    return dim, sections


class _Cursor:
    """Bounds-checked reads from the bytes of a binary file."""

    def __init__(self, data):
        self.data = data
        self.at = 0
        self.order = "<"

    def left(self):
        return len(self.data) - self.at

    def array(self, count, dtype):
        dtype = np.dtype(self.order + dtype)
        if count < 0 or count * dtype.itemsize > self.left():
            raise ReadError("Medit: the solution file is truncated")
        out = np.frombuffer(self.data, dtype, count, self.at)
        self.at += count * dtype.itemsize
        return out

    def int(self, dtype):
        return int(self.array(1, dtype)[0])


def _read_binary_sections(data):
    cursor = _Cursor(data)
    code = cursor.int("i4")
    if code == 16777216:
        cursor.order = ">"
    elif code != 1:
        raise ReadError("Medit: invalid code")
    version = cursor.int("i4")
    if version < 1 or version > 4:
        raise ReadError("Medit: invalid version")
    itype = "i8" if version == 4 else "i4"
    ftype = "f4" if version == 1 else "f8"
    postype = "i8" if version >= 3 else "i4"
    if cursor.int("i4") != 3:
        raise ReadError("Medit: expected Dimension")
    cursor.int(postype)
    dim = cursor.int("i4")
    if dim not in (2, 3):
        raise ReadError(f"Medit: invalid solution dimension {dim}")

    sections = []
    while True:
        if cursor.left() < 4:
            warn("End-of-file reached before GmfEnd keyword")
            break
        code = cursor.int("i4")
        if code < 0 or code > 198:
            raise ReadError(f"Medit: unsupported field {code}")
        if code == 54:
            break
        if code in _RESERVED:
            continue
        following = cursor.int(postype)
        section_start = cursor.at
        if code in _SOL_CODES:
            keyword = _SOL_CODES[code][0]
            count = cursor.int(itype)
            num_types = cursor.int("i4")
            if count < 0 or num_types < 0 or num_types * 4 > cursor.left():
                raise ReadError(f"Medit: {keyword} has an impossible size")
            types = [int(t) for t in cursor.array(num_types, "i4")]
            total = _check_types(num_types, types, dim)
            values = cursor.array(count * total, ftype).reshape(count, total)
            sections.append(Section(keyword, types, values.astype(ftype)))
        else:
            if following <= section_start or following > len(data):
                raise ReadError(
                    f"Medit: cannot skip keyword {code}: its next-section position "
                    "is outside the file"
                )
            warn(f"meshio++ doesn't know Medit keyword {code}. Skipping.")
            cursor.at = following
    return dim, sections


def _vertex_fields(dim, sections):
    out = {}
    count = None
    for section in sections:
        if section.cell_type is not None:
            continue
        if count is not None and len(section.values) != count:
            raise ReadError("Medit: two SolAtVertices sections of different sizes")
        count = len(section.values)
        for field in section.fields(dim):
            out[f"medit:sol:{len(out)}"] = field
    return count, out


def read(filename):
    """A geometry-less mesh carrying the file's vertex fields (see the module doc)."""
    dim, sections = read_sections(filename)
    count, point_data = _vertex_fields(dim, sections)
    for section in sections:
        if section.cell_type is not None:
            warn(
                f"Medit: {section.keyword} needs the mesh it belongs to; read it with "
                "medit.read(mesh_path, sol=...). Skipping."
            )
    return Mesh(
        np.empty((count or 0, 0)),
        [],
        point_data=point_data,
        field_data={"medit:dimension": np.array([dim])},
    )


def attach(mesh, dim, sections):
    """Put the solution's fields on ``mesh``: vertex fields as point data, element
    fields as cell data (``NaN`` on blocks of the other types). Returns ``mesh``."""
    num_points, d = mesh.points.shape
    if dim != d:
        raise ReadError(f"Medit: the solution is {dim}-D and the mesh {d}-D")
    count, point_data = _vertex_fields(dim, sections)
    if count is not None and count != num_points:
        raise ReadError(
            f"Medit: SolAtVertices has {count} rows and the mesh {num_points} points"
        )
    mesh.point_data.update(point_data)

    cell_fields = {}  # name -> {block index: values}
    for section in sections:
        if section.cell_type is None:
            continue
        blocks = [k for k, cb in enumerate(mesh.cells) if cb.type == section.cell_type]
        total = sum(len(mesh.cells[k].data) for k in blocks)
        if total != len(section.values):
            raise ReadError(
                f"Medit: {section.keyword} has {len(section.values)} rows and the "
                f"mesh {total} {section.cell_type} cells"
            )
        for j, field in enumerate(section.fields(dim)):
            name = f"medit:sol:{j}"
            known = cell_fields.get(name)
            if known is not None:
                other = next(iter(known.values()))
                if other.shape[1:] != field.shape[1:]:
                    name = f"medit:sol:{j}:{section.cell_type}"
            per_block = cell_fields.setdefault(name, {})
            start = 0
            for k in blocks:
                n = len(mesh.cells[k].data)
                per_block[k] = field[start : start + n]
                start += n
    for name, per_block in cell_fields.items():
        sample = next(iter(per_block.values()))
        mesh.cell_data[name] = [
            per_block.get(
                k, np.full((len(cb.data),) + sample.shape[1:], np.nan, sample.dtype)
            )
            for k, cb in enumerate(mesh.cells)
        ]
    return mesh


def _dimension(mesh):
    d = mesh.points.shape[1] if mesh.points.ndim == 2 else 0
    if d in (2, 3):
        return d
    dim = (getattr(mesh, "field_data", None) or {}).get("medit:dimension")
    if dim is not None and int(np.asarray(dim).reshape(-1)[0]) in (2, 3):
        return int(np.asarray(dim).reshape(-1)[0])
    raise WriteError(
        "Medit: a solution is 2-D or 3-D; give the mesh 2-D or 3-D points or "
        "field_data['medit:dimension']"
    )


def _classified(name, array, dim):
    """``(sol type, (n, components) float64 values)``, or None (with a warning)."""
    array = np.asarray(array)
    if array.dtype.kind != "f":
        if name != "medit:ref":
            warn(f"Medit: a solution holds reals; '{name}' is not written")
        return None
    values = array.reshape(len(array), -1)
    sol_type = sol_type_of(values.shape[1], dim)
    if sol_type is None:
        warn(
            f"Medit: '{name}' has {values.shape[1]} components, which no {dim}-D "
            "solution type holds; not written"
        )
        return None
    return sol_type, values.astype(np.float64)


def sections_of(mesh):
    """``(dim, [Section, ...])`` the mesh's floating-point data is written as."""
    dim = _dimension(mesh)
    sections = []
    types, columns = [], []
    for name in sorted(mesh.point_data):
        item = _classified(name, mesh.point_data[name], dim)
        if item is not None:
            types.append(item[0])
            columns.append(item[1])
    if types:
        sections.append(Section("SolAtVertices", types, np.hstack(columns)))

    by_type = {}
    for k, cb in enumerate(mesh.cells):
        by_type.setdefault(cb.type, []).append(k)
    sol_cell_types = {cell_type for _, cell_type in SOL_KEYWORDS.values()}
    names = []
    for name in sorted(mesh.cell_data):
        arrays = [np.asarray(a) for a in mesh.cell_data[name]]
        if name == "medit:ref" or not arrays:
            continue
        if any(a.dtype.kind != "f" for a in arrays):
            warn(f"Medit: a solution holds reals; '{name}' is not written")
            continue
        names.append(name)
        for k, cb in enumerate(mesh.cells):
            if cb.type not in sol_cell_types and not np.all(np.isnan(arrays[k])):
                warn(f"Medit: a solution has no {cb.type} section; '{name}' skips it")
    for keyword, (_, cell_type) in SOL_KEYWORDS.items():
        blocks = by_type.get(cell_type)
        if cell_type is None or not blocks:
            continue
        types, columns = [], []
        for name in names:
            joined = np.concatenate(
                [np.asarray(mesh.cell_data[name][k]) for k in blocks]
            )
            if len(joined) and np.all(np.isnan(joined)):
                continue  # the fill of a field other cell types carry
            item = _classified(name, joined, dim)
            if item is not None:
                types.append(item[0])
                columns.append(item[1])
        if types:
            sections.append(Section(keyword, types, np.hstack(columns)))
    if not sections:
        raise WriteError(
            "Medit: no floating-point point or cell data to write as a solution"
        )
    return dim, sections


def write(filename, mesh, float_fmt=".16e"):
    """Write the mesh's floating-point point and cell data as a ``.sol``, or
    (ending in ``b``) a little-endian ``.solb``."""
    dim, sections = sections_of(mesh)
    if str(filename)[-1] == "b":
        _write_binary(filename, dim, sections)
    else:
        _write_ascii(filename, dim, sections, float_fmt)


def _write_ascii(filename, dim, sections, float_fmt):
    with open_file(filename, "wb") as fh:
        fh.write(f"MeshVersionFormatted 2\nDimension {dim}\n".encode())
        for section in sections:
            n, total = section.values.shape
            types = " ".join(str(t) for t in section.types)
            fh.write(
                f"\n{section.keyword}\n{n}\n{len(section.types)} {types}\n".encode()
            )
            fmt = " ".join(["{:" + float_fmt + "}"] * total) + "\n"
            fh.write("".join(fmt.format(*row) for row in section.values).encode())
        fh.write(b"\nEnd\n")


def _write_binary(filename, dim, sections):
    big = any(len(s.values) > np.iinfo(np.int32).max for s in sections)
    version = 4 if big else 3
    isize = 8 if big else 4
    itype = "<i8" if big else "<i4"
    psize = 8
    pos = 4 * 4 + psize
    out = [np.array([1, version, 3], "<i4").tobytes(), np.array([pos], "<i8").tobytes()]
    out.append(np.array([dim], "<i4").tobytes())
    for section in sections:
        n, total = section.values.shape
        nbytes = 4 + 4 * len(section.types) + n * total * 8
        pos += nbytes + 4 + psize + isize
        out.append(np.array([SOL_KEYWORDS[section.keyword][0]], "<i4").tobytes())
        out.append(np.array([pos], "<i8").tobytes())
        out.append(np.array([n], itype).tobytes())
        out.append(np.array([len(section.types)] + section.types, "<i4").tobytes())
        out.append(np.ascontiguousarray(section.values, "<f8").tobytes())
    out.append(np.array([54], "<i4").tobytes())
    out.append(np.array([0], "<i8").tobytes())
    with open_file(filename, "wb") as fh:
        fh.write(b"".join(out))
