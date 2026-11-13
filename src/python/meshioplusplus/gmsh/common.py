from __future__ import annotations

import shlex

import numpy as np
from numpy.typing import ArrayLike

from .. import _node_order
from .._common import cell_data_from_raw, warn
from .._exceptions import ReadError, WriteError

c_int = np.dtype("int32")
c_double = np.dtype("float64")


def _header_int(f, tag):
    """One integer header line of a data section; anything else is ReadError."""
    line = f.readline().decode()
    try:
        return int(line)
    except ValueError:
        raise ReadError(
            f"Gmsh ${tag}: header line {line.strip()!r} is not an integer"
        ) from None


def checked_cell_data(cells, cell_data_raw):
    """``cell_data_from_raw`` for the Gmsh readers: short cell data is refused.

    The other formats share ``cell_data_from_raw``, which splits silently; a
    Gmsh ElementData array with fewer values than the elements it annotates is
    corrupt input, the same case the C++ reader refuses.
    """
    total = sum(len(c) for c in cells)
    for data in cell_data_raw.values():
        if len(data) < total:
            raise ReadError(
                "Gmsh: cell data has fewer values than the elements it annotates"
            )
    return cell_data_from_raw(cells, cell_data_raw)


def _remap_periodic(periodic, point_tags):
    """Resolve file node tags to point rows without a max-tag-sized allocation."""
    if not periodic:
        return periodic
    rows = {int(tag): i for i, tag in enumerate(point_tags)}
    if len(rows) != len(point_tags):
        raise ReadError("Gmsh $Periodic: duplicate node tag")
    for record in periodic:
        dim, (slave, master), affine, pairs = record
        if not 0 <= dim <= 3 or slave < 1 or master < 1:
            raise ReadError("Gmsh $Periodic: invalid entity dimension/tag")
        if affine is not None and (
            len(affine) not in (0, 16) or not np.isfinite(affine).all()
        ):
            raise ReadError("Gmsh $Periodic: invalid affine coefficients")
        try:
            record[3] = np.array(
                [rows[int(tag)] for tag in pairs.flat], dtype=np.int64
            ).reshape(-1, 2)
        except KeyError as exc:
            raise ReadError("Gmsh $Periodic: node tag outside $Nodes") from exc
    return periodic


def _validate_periodic(periodic, num_points):
    """Validate metadata before a reference writer opens its output file."""
    for record in periodic or []:
        if len(record) != 4 or len(record[1]) != 2:
            raise WriteError("Gmsh $Periodic: invalid record shape")
        dim, (slave, master), affine, pairs = record
        if any(not isinstance(v, (int, np.integer)) for v in [dim, slave, master]):
            raise WriteError("Gmsh $Periodic: entity tags must be integers")
        if not 0 <= dim <= 3 or not all(
            0 < v <= np.iinfo(np.int32).max for v in [slave, master]
        ):
            raise WriteError("Gmsh $Periodic: invalid entity dimension/tag")
        if affine is not None:
            affine = np.asarray(affine, dtype=float)
            if (
                affine.ndim != 1
                or affine.size not in (0, 16)
                or not np.isfinite(affine).all()
            ):
                raise WriteError(
                    "Gmsh $Periodic: affine must have 0 or 16 finite coefficients"
                )
        pairs = np.asarray(pairs)
        if pairs.ndim != 2 or pairs.shape[1] != 2 or pairs.dtype.kind not in "iu":
            raise WriteError("Gmsh $Periodic: pairs must be integer (N,2)")
        if np.any(pairs < 0) or np.any(pairs >= num_points):
            raise WriteError("Gmsh $Periodic: point index outside mesh")


def _periodic_count(f, value, stride=1):
    """Refuse negative or implausible counts before NumPy allocates."""
    value = int(value)
    position = f.tell()
    f.seek(0, 2)
    remaining = f.tell() - position
    f.seek(position)
    if value < 0 or value > remaining // stride:
        raise ReadError("Gmsh $Periodic: invalid count")
    return value


def _periodic_end(f):
    marker, _ = _fast_forward_over_blank_lines(f)
    if marker != "$EndPeriodic":
        raise ReadError("Gmsh $Periodic: expected $EndPeriodic")


def _fast_forward_to_end_block(f, block):
    """fast-forward to end of block"""
    # See also https://github.com/nschloe/pygalmesh/issues/34

    for line in f:
        try:
            line = line.decode()
        except UnicodeDecodeError:
            pass
        if line.strip() == f"$End{block}":
            break
    else:
        warn(f"${block} not closed by $End{block}.")


def _fast_forward_over_blank_lines(f):
    is_eof = False
    while True:
        line = f.readline().decode()
        if not line:
            is_eof = True
            break
        elif len(line.strip()) > 0:
            break
    # Stripped: some writers (FEconv's samples) indent every line.
    return line.strip(), is_eof


def _read_physical_names(f, field_data):
    line = f.readline().decode()
    num_phys_names = int(line)
    for _ in range(num_phys_names):
        line = shlex.split(f.readline().decode())
        key = line[2]
        value = np.array(line[1::-1], dtype=int)
        field_data[key] = value
    _fast_forward_to_end_block(f, "PhysicalNames")


#: The sections kept verbatim as ``field_data["gmsh:section:<Name>"]`` (their text body, as
#: ``uint8``) and written back: ghost elements, interpolation schemes, partitioned entities
#: and parametrizations. ASCII files only. Twins of ``kGmshKeptSections`` in ``gmsh.cpp``.
SECTION_PREFIX = "gmsh:section:"
#: Private key carrying each block's element tags from ``_read_elements`` to the reader.
ELEMENT_TAGS = "\0element_tags"
KEPT_SECTIONS = (
    "GhostElements",
    "InterpolationScheme",
    "PartitionedEntities",
    "Parametrizations",
)


def _capture_section(f, name, field_data):
    """A kept section's body, line by line to its ``$End`` line, into ``field_data``."""
    target = f"$End{name}"
    body = b""
    for line in f:
        if line.decode("latin-1").strip() == target:
            break
        body += line.rstrip(b"\r\n") + b"\n"
    field_data[SECTION_PREFIX + name] = np.frombuffer(body, dtype=np.uint8).copy()


def _ghost_cell_data(field_data, tag_blocks):
    """``{"partition:ghost": [per block]}`` from a kept ``$GhostElements`` section: 1 for
    an element the section lists (a ghost), 0 for an owned one. ``tag_blocks`` holds the
    file's element tags per cell block. Twin of ``gmsh_add_ghost_cells`` in ``gmsh.cpp``.
    """
    body = field_data.get(SECTION_PREFIX + "GhostElements")
    if body is None or tag_blocks is None:
        return {}
    words = np.asarray(body, dtype=np.uint8).tobytes().split()
    try:
        count = int(words[0])
        pos, ghosts = 1, set()
        for _ in range(count):
            tag, num = int(words[pos]), int(words[pos + 2])
            if num < 0:
                return {}
            ghosts.add(tag)
            pos += 3 + num
    except (IndexError, ValueError):
        return {}
    if count < 0:
        return {}
    return {
        "partition:ghost": [
            np.array([int(t) in ghosts for t in tags], dtype=np.int64)
            for tags in tag_blocks
        ]
    }


def _write_kept_sections(fh, field_data, binary):
    """The kept sections back as ``$Name ... $EndName``; ASCII only."""
    if binary:
        return
    for name in KEPT_SECTIONS:
        body = field_data.get(SECTION_PREFIX + name)
        if body is None:
            continue
        fh.write(f"${name}\n".encode())
        fh.write(np.asarray(body, dtype=np.uint8).tobytes())
        fh.write(f"$End{name}\n".encode())


def _read_data_header(f, tag):
    """The tags of a ``$NodeData``/``$ElementData``/``$ElementNodeData`` section:
    ``(first string tag, component count, item count)``."""
    # Read string tags
    num_string_tags = _header_int(f, tag)
    string_tags = [
        f.readline().decode().strip().replace('"', "") for _ in range(num_string_tags)
    ]
    # The real tags typically only contain one value, the time.
    # Discard it.
    num_real_tags = _header_int(f, tag)
    for _ in range(num_real_tags):
        f.readline()
    num_integer_tags = _header_int(f, tag)
    if num_integer_tags < 3:
        raise ReadError(f"Gmsh ${tag}: needs at least three integer tags")
    integer_tags = [_header_int(f, tag) for _ in range(num_integer_tags)]
    num_components = integer_tags[1]
    num_items = integer_tags[2]
    if num_components < 0 or num_items < 0:
        raise ReadError(f"Gmsh ${tag}: negative component or item count")
    return string_tags[0] if string_tags else "", num_components, num_items


#: The ``cell_data`` name prefix of a ``$ElementNodeData`` section (a value per
#: element per node), and of the ``field_data`` entry holding its component count.
#: Twins of ``kGmshElementNodePrefix`` and ``kGmshElementNodeComponents`` in ``gmsh.cpp``.
ELEMENT_NODE_PREFIX = "gmsh:element_node:"
ELEMENT_NODE_COMPONENTS = "gmsh:element_node_components:"


def _read_element_node_data(f, out, is_ascii):
    """A ``$ElementNodeData`` section into ``out[name] = (components, rows)``, one
    row per element (its ``nodes * components`` values, in the file's node order)."""
    tag = "ElementNodeData"
    name, num_components, num_items = _read_data_header(f, tag)
    rows = []
    if is_ascii:
        for _ in range(num_items):
            parts = f.readline().decode().split()
            if len(parts) < 2:
                raise ReadError(f"Gmsh ${tag}: truncated element row")
            count = int(parts[1])
            values = np.array([float(v) for v in parts[2:]])
            if len(values) != count * num_components:
                raise ReadError(f"Gmsh ${tag}: an element row has the wrong length")
            rows.append(values)
    else:
        for _ in range(num_items):
            head = np.fromfile(f, dtype=c_int, count=2)
            if len(head) != 2:
                raise ReadError(f"Gmsh ${tag}: truncated element row")
            values = np.fromfile(f, dtype=c_double, count=int(head[1]) * num_components)
            rows.append(values)
    _fast_forward_to_end_block(f, tag)
    out[name] = (num_components, rows)


def _element_node_cell_data(cells, element_node):
    """``cell_data`` and ``field_data`` entries for the ``$ElementNodeData`` sections:
    per block, the rows of its elements as ``(cells, nodes * components)``."""
    sizes = [len(c[1] if isinstance(c, tuple) else c.data) for c in cells]
    total = sum(sizes)
    cell_data, field_data = {}, {}
    for name, (num_components, rows) in element_node.items():
        if len(rows) < total:
            raise ReadError(
                "Gmsh: element-node data has fewer rows than the elements it annotates"
            )
        blocks, at = [], 0
        for size in sizes:
            chunk = rows[at : at + size]
            at += size
            widths = {len(r) for r in chunk}
            if len(widths) > 1:
                raise ReadError("Gmsh: element-node rows of one block differ in length")
            blocks.append(
                np.array(chunk, dtype=float).reshape(len(chunk), -1)
                if chunk
                else np.empty((0, 0))
            )
        cell_data[ELEMENT_NODE_PREFIX + name] = blocks
        field_data[ELEMENT_NODE_COMPONENTS + name] = np.array(
            [num_components], dtype=np.int64
        )
    return cell_data, field_data


def _write_element_node_data(fh, name, blocks, num_components, binary):
    """One ``$ElementNodeData`` section from per-block ``(cells, nodes * components)``
    arrays; elements are tagged 1 to N in block order."""
    fh.write(b"$ElementNodeData\n")
    fh.write(f'1\n"{name}"\n1\n{0.0}\n3\n0\n{num_components}\n'.encode())
    total = sum(len(b) for b in blocks)
    fh.write(f"{total}\n".encode())
    tag = 1
    for block in blocks:
        block = np.asarray(block, dtype=float).reshape(len(block), -1)
        nodes = block.shape[1] // num_components if num_components else 0
        for row in block:
            if binary:
                np.array([tag, nodes], dtype=c_int).tofile(fh)
                row.astype(c_double).tofile(fh)
            else:
                fh.write(
                    (
                        f"{tag} {nodes} " + " ".join(repr(float(v)) for v in row) + "\n"
                    ).encode()
                )
            tag += 1
    if binary:
        fh.write(b"\n")
    fh.write(b"$EndElementNodeData\n")


def _split_element_node(cell_data, field_data):
    """``(plain cell_data, {name: (blocks, components)})``: the ``gmsh:element_node:``
    arrays set apart, for the writers."""
    plain, special = {}, {}
    for key, value in cell_data.items():
        if key.startswith(ELEMENT_NODE_PREFIX):
            name = key[len(ELEMENT_NODE_PREFIX) :]
            comp = field_data.get(ELEMENT_NODE_COMPONENTS + name)
            if comp is None:
                raise WriteError(
                    f"Gmsh: '{key}' needs field_data['{ELEMENT_NODE_COMPONENTS}{name}'] "
                    "holding its component count"
                )
            special[name] = (value, int(np.asarray(comp).ravel()[0]))
        else:
            plain[key] = value
    return plain, special


def _read_data(f, tag, data_dict, data_size, is_ascii):
    string_tags, num_components, num_items = _read_data_header(f, tag)
    string_tags = [string_tags]
    if is_ascii:
        # We need to read num_items * (1 + num_components) floats.
        # np.fromfile(..., sep=" ") can be flaky if there are newlines or other
        # whitespace issues.
        # Instead, read the raw string and split it.
        data = []
        while len(data) < num_items * (1 + num_components):
            line = f.readline().decode().split()
            if not line:
                break
            try:
                data.extend([float(val) for val in line])
            except ValueError:
                raise ReadError(
                    f"Gmsh ${tag}: value line {line!r} is not numeric"
                ) from None
        data = np.array(data).reshape((num_items, 1 + num_components))
        # The first entry is the node number
        data = data[:, 1:]
    else:
        # binary
        dtype = [("index", c_int), ("values", c_double, (num_components,))]
        data = np.fromfile(f, count=num_items, dtype=dtype)
        if not (data["index"] == range(1, num_items + 1)).all():
            raise ReadError()
        data = np.ascontiguousarray(data["values"])

    _fast_forward_to_end_block(f, tag)

    # The gmsh format cannot distinguish between data of shape (n,) and (n, 1).
    # If shape[1] == 1, cut it off.
    if data.shape[1] == 1:
        data = data[:, 0]

    data_dict[string_tags[0]] = data


# Translate meshio types to gmsh codes
# http://gmsh.info//doc/texinfo/gmsh.html#MSH-file-format-version-2
_gmsh_to_meshio_type = {
    1: "line",
    2: "triangle",
    3: "quad",
    4: "tetra",
    5: "hexahedron",
    6: "wedge",
    7: "pyramid",
    8: "line3",
    9: "triangle6",
    10: "quad9",
    11: "tetra10",
    12: "hexahedron27",
    13: "wedge18",
    14: "pyramid14",
    15: "vertex",
    16: "quad8",
    17: "hexahedron20",
    18: "wedge15",
    19: "pyramid13",
    21: "triangle10",
    23: "triangle15",
    25: "triangle21",
    26: "line4",
    27: "line5",
    28: "line6",
    29: "tetra20",
    30: "tetra35",
    31: "tetra56",
    36: "quad16",
    37: "quad25",
    38: "quad36",
    42: "triangle28",
    43: "triangle36",
    44: "triangle45",
    45: "triangle55",
    46: "triangle66",
    47: "quad49",
    48: "quad64",
    49: "quad81",
    50: "quad100",
    51: "quad121",
    62: "line7",
    63: "line8",
    64: "line9",
    65: "line10",
    66: "line11",
    71: "tetra84",
    72: "tetra120",
    73: "tetra165",
    74: "tetra220",
    75: "tetra286",
    90: "wedge40",
    91: "wedge75",
    92: "hexahedron64",
    93: "hexahedron125",
    94: "hexahedron216",
    95: "hexahedron343",
    96: "hexahedron512",
    97: "hexahedron729",
    98: "hexahedron1000",
    106: "wedge126",
    107: "wedge196",
    108: "wedge288",
    109: "wedge405",
    110: "wedge550",
}
_meshio_to_gmsh_type = {v: k for k, v in _gmsh_to_meshio_type.items()}


def _gmsh_to_meshio_order(cell_type: str, idx: ArrayLike) -> np.ndarray:
    # Gmsh cells are mostly ordered like VTK; the exceptions are gmsh's
    # entries in the node-ordering registry (_node_order.py).
    return _node_order.to_meshio("gmsh", cell_type, np.asarray(idx))


def _meshio_to_gmsh_order(cell_type: str, idx: ArrayLike) -> np.ndarray:
    return _node_order.from_meshio("gmsh", cell_type, np.asarray(idx))


def _write_physical_names(fh, field_data):
    # Write physical names
    entries = []
    for phys_name in field_data:
        if phys_name.startswith((ELEMENT_NODE_COMPONENTS, SECTION_PREFIX)):
            continue
        try:
            phys_num, phys_dim = field_data[phys_name]
            phys_num, phys_dim = int(phys_num), int(phys_dim)
            entries.append((phys_dim, phys_num, phys_name))
        except (ValueError, TypeError):
            warn("Field data contains entry that cannot be processed.")
    entries.sort()
    if entries:
        fh.write(b"$PhysicalNames\n")
        fh.write(f"{len(entries)}\n".encode())
        for entry in entries:
            fh.write('{} {} "{}"\n'.format(*entry).encode())
        fh.write(b"$EndPhysicalNames\n")


# The time and step index of the data sections of a write, set by
# ``gmsh.SeriesWriter`` around it (a plain write leaves it ``None``: 0.0 and 0).
_SERIES_LABEL = None


def _write_data(fh, tag, name, data, binary):
    fh.write(f"${tag}\n".encode())
    # <http://gmsh.info/doc/texinfo/gmsh.html>:
    # > Number of string tags.
    # > gives the number of string tags that follow. By default the first
    # > string-tag is interpreted as the name of the post-processing view and
    # > the second as the name of the interpolation scheme. The interpolation
    # > scheme is provided in the $InterpolationScheme section (see below).
    fh.write(f"{1}\n".encode())
    fh.write(f'"{name}"\n'.encode())
    fh.write(f"{1}\n".encode())
    time, step = _SERIES_LABEL if _SERIES_LABEL is not None else (0.0, 0)
    fh.write((f"{0.0}\n" if _SERIES_LABEL is None else f"{time:.17g}\n").encode())
    # three integer tags:
    fh.write(f"{3}\n".encode())
    # time step
    fh.write(f"{step}\n".encode())
    # number of components
    num_components = data.shape[1] if len(data.shape) > 1 else 1
    if num_components not in [1, 3, 9]:
        raise WriteError("Gmsh only permits 1, 3, or 9 components per data field.")

    # Cut off the last dimension in case it's 1. This avoids problems with
    # writing the data.
    if len(data.shape) > 1 and data.shape[1] == 1:
        data = data[:, 0]

    fh.write(f"{num_components}\n".encode())
    # num data items
    fh.write(f"{data.shape[0]}\n".encode())
    # actually write the data
    if binary:
        if num_components == 1:
            dtype = [("index", c_int), ("data", c_double)]
        else:
            dtype = [("index", c_int), ("data", c_double, num_components)]
        tmp = np.empty(len(data), dtype=dtype)
        tmp["index"] = 1 + np.arange(len(data))
        tmp["data"] = data
        tmp.tofile(fh)
        fh.write(b"\n")
    else:
        fmt = " ".join(["{}"] + ["{!s}"] * num_components) + "\n"
        # TODO unify
        if num_components == 1:
            for k, x in enumerate(data):
                fh.write(fmt.format(k + 1, x).encode())
        else:
            for k, x in enumerate(data):
                fh.write(fmt.format(k + 1, *x).encode())

    fh.write(f"$End{tag}\n".encode())
