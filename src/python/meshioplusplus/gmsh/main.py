import copy
import pathlib
import struct
from warnings import warn

import numpy as np

from .._exceptions import ReadError, WriteError
from . import _gmsh22, _gmsh40, _gmsh41
from .common import _fast_forward_to_end_block, _validate_periodic

# Some mesh files out there have the version specified as version "2" when it really is
# "2.2". Same with "4" vs "4.1".
_readers = {"2": _gmsh22, "2.2": _gmsh22, "4.0": _gmsh40, "4": _gmsh41, "4.1": _gmsh41}
_writers = {"2.2": _gmsh22, "4.0": _gmsh40, "4.1": _gmsh41}


def read(filename):
    """Reads a Gmsh msh file."""
    filename = pathlib.Path(filename)
    with open(filename.as_posix(), "rb") as f:
        mesh = read_buffer(f)
    return mesh


def read_buffer(f):
    # The various versions of the format are specified at
    # <http://gmsh.info/doc/texinfo/gmsh.html#File-formats>.
    line = f.readline().decode().strip()

    # skip any $Comments/$EndComments sections
    while line == "$Comments":
        _fast_forward_to_end_block(f, "Comments")
        line = f.readline().decode().strip()

    if line != "$MeshFormat":
        raise ReadError()
    fmt_version, data_size, is_ascii = _read_header(f)

    try:
        reader = _readers[fmt_version]
    except KeyError:
        try:
            reader = _readers[fmt_version.split(".")[0]]
        except KeyError:
            raise ValueError(
                "Need mesh format in {} (got {})".format(
                    sorted(_readers.keys()), fmt_version
                )
            )
    return reader.read_buffer(f, is_ascii, data_size)


def _read_header(f):
    """Read the mesh format block

    specified as

     version(ASCII double; currently 4.1)
       file-type(ASCII int; 0 for ASCII mode, 1 for binary mode)
       data-size(ASCII int; sizeof(size_t))
     < int with value one; only in binary mode, to detect endianness >

    though here the version is left as str
    """

    # http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format

    line = f.readline().decode()
    # Split the line
    # 4.1 0 8
    # into its components.
    str_list = list(filter(None, line.split()))
    fmt_version = str_list[0]
    if str_list[1] not in ["0", "1"]:
        raise ReadError()
    is_ascii = str_list[1] == "0"
    data_size = int(str_list[2])
    if not is_ascii:
        # The next line is the integer 1 in bytes. Useful for checking endianness.
        # Just assert that we get 1 here.
        one = f.read(struct.calcsize("i"))
        if struct.unpack("i", one)[0] != 1:
            raise ReadError()
    _fast_forward_to_end_block(f, "MeshFormat")
    return fmt_version, data_size, is_ascii


# Gmsh ASCII output uses `%.16g` for floating point values,
# meshio uses same precision but exponential notation `%.16e`.
def _region_dim(mesh, region, block_ends):
    """The dimension a cell region is written under: its own, else its first cell's."""
    from .._mesh import topological_dimension

    if region.dim >= 0:
        return region.dim
    if len(region.entries) == 0:
        return -1
    block = int(np.searchsorted(block_ends, region.entries[0], side="right"))
    if block >= len(mesh.cells):
        return -1
    return topological_dimension.get(mesh.cells[block].type, -1)


def _with_region_groups(mesh):
    """A mesh whose ``Cell`` regions are gmsh physical groups.

    The C++ writer's ``gmsh_resolve_region_tags`` rule, in Python: a mesh with
    no ``gmsh:physical`` cell data of its own gets one synthesised from its
    ``Cell`` regions (a region's own tag, else a fresh one per dimension, in
    region order), and each group a ``field_data`` entry does not already
    describe is added to ``$PhysicalNames``. The mesh itself is not modified.
    """
    regions = [r for r in mesh.regions if r.kind == "cell"]
    if not regions or "gmsh:physical" in mesh.cell_data:
        return mesh
    ends = np.cumsum([len(c) for c in mesh.cells], dtype=np.int64)
    max_tag = {}
    for value in mesh.field_data.values():
        value = np.asarray(value).ravel()
        if value.size >= 2:
            dim = int(value[1])
            max_tag[dim] = max(max_tag.get(dim, 0), int(value[0]))
    resolved = {}
    for i, r in enumerate(regions):
        if r.tag >= 0:
            dim = _region_dim(mesh, r, ends)
            resolved[i] = (dim, r.tag)
            max_tag[dim] = max(max_tag.get(dim, 0), r.tag)
    for i, r in enumerate(regions):
        if r.tag >= 0:
            continue
        dim = _region_dim(mesh, r, ends)
        if dim < 0:
            continue
        max_tag[dim] = max_tag.get(dim, 0) + 1
        resolved[i] = (dim, max_tag[dim])

    field_data = dict(mesh.field_data)
    for i, r in enumerate(regions):
        if i in resolved and r.name not in field_data:
            dim, tag = resolved[i]
            field_data[r.name] = np.array([tag, dim], dtype=np.int64)

    total = int(ends[-1]) if len(ends) else 0
    flat = np.zeros(total, dtype=np.int64)
    claimed = {}
    any_tagged = False
    for i, r in enumerate(regions):
        if i not in resolved or len(r.entries) == 0:
            continue
        any_tagged = True
        for c in r.entries:
            c = int(c)
            if not 0 <= c < total:
                continue
            if c in claimed:
                warn(
                    f"Gmsh writer: cell {c} is in both region "
                    f"'{claimed[c]}' and '{r.name}'; keeping '{claimed[c]}'"
                )
                continue
            claimed[c] = r.name
            flat[c] = resolved[i][1]

    out = copy.copy(mesh)
    out.field_data = field_data
    if any_tagged:
        starts = np.concatenate([[0], ends[:-1]]) if len(ends) else ends
        out.cell_data = dict(mesh.cell_data)
        out.cell_data["gmsh:physical"] = [flat[a:b] for a, b in zip(starts, ends)]
    return out


def write(filename, mesh, fmt_version="4.1", binary=True, float_fmt=".16e"):
    """Writes a Gmsh msh file."""
    _validate_periodic(getattr(mesh, "gmsh_periodic", None), len(mesh.points))
    try:
        writer = _writers[fmt_version]
    except KeyError:
        try:
            writer = _writers[fmt_version]
        except KeyError:
            raise WriteError(
                "Need mesh format in {} (got {})".format(
                    sorted(_writers.keys()), fmt_version
                )
            )

    writer.write(
        filename, _with_region_groups(mesh), binary=binary, float_fmt=float_fmt
    )


# NOTE: format registration now lives in meshioplusplus/gmsh/__init__.py, which wraps the
# reader/writer above with the C++-backed fast paths (reads 2.2/4.0/4.1,
# writes 2.2/4.1).
