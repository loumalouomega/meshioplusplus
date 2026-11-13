"""
CGNS <https://cgns.github.io/> reader/writer -- a genuine CGNS/SIDS-compliant
unstructured-mesh subset stored in HDF5 (the ADF-over-HDF5 mapping every real
CGNS tool uses), readable by cgnslib/ParaView/VTK. This is the Python/h5py
twin of the C++ core (`src/cpp/src/formats/cgns.cpp`) -- kept in exact
structural sync, since the shim (`cgns/__init__.py`) catches every exception
from the C++ path and falls back here, so a gap in this module turns a C++
`WriteError` into a silently different file. See doc/formats/cgns.md for the
full on-disk layout, the supported/rejected cell-type table, and what CI can
and cannot verify.

**The `" data"` leading-space dataset name is not an ad hoc convention** --
it is cgnslib's own ADF-over-HDF5 mapping (`#define D_DATA " data"` in
`ADFH.c`); every node with a payload stores it in a child dataset by that
exact name.
"""

import warnings

import numpy as np

from .. import _node_order, _provenance
from .._common import num_nodes_per_cell
from .._exceptions import ReadError, WriteError
from .._grid import lattice_cells
from .._mesh import Mesh, topological_dimension
from .._regions import Region, canonical_regions
from .._units import (
    CGNS_LENGTH_UNITS,
    LENGTH_UNIT_KEY,
    length_unit,
    set_length_unit,
    unit_name,
)

# meshio++ <-> CGNS ElementType_t table, kept in exact sync with
# src/cpp/src/formats/cgns.cpp's cgns_type_table(). The SIDS node order of
# PENTA_15/18 and HEXA_20/27 differs from meshio++'s; those permutations are
# the "cgns" entries of the node-ordering registry (_node_order.py).
# name -> (cgns_name, ElementType_t code)
_CGNS_TYPES = {
    "vertex": ("NODE", 2),
    "line": ("BAR_2", 3),
    "line3": ("BAR_3", 4),
    "triangle": ("TRI_3", 5),
    "triangle6": ("TRI_6", 6),
    "quad": ("QUAD_4", 7),
    "quad8": ("QUAD_8", 8),
    "quad9": ("QUAD_9", 9),
    "tetra": ("TETRA_4", 10),
    "tetra10": ("TETRA_10", 11),
    "pyramid": ("PYRA_5", 12),
    "pyramid14": ("PYRA_14", 13),
    "wedge": ("PENTA_6", 14),
    "wedge15": ("PENTA_15", 15),
    "wedge18": ("PENTA_18", 16),
    "hexahedron": ("HEXA_8", 17),
    "hexahedron20": ("HEXA_20", 18),
    "hexahedron27": ("HEXA_27", 19),
    # PYRA_13's code (21) is non-monotonic -- appended to ElementType_t after
    # MIXED in the CGNS enum -- not a typo.
    "pyramid13": ("PYRA_13", 21),
}
_CODE_TO_TYPE = {code: name for name, (_cgns_name, code) in _CGNS_TYPES.items()}
_CODE_TO_CGNS_NAME = {
    2: "NODE",
    3: "BAR_2",
    4: "BAR_3",
    5: "TRI_3",
    6: "TRI_6",
    7: "QUAD_4",
    8: "QUAD_8",
    9: "QUAD_9",
    10: "TETRA_4",
    11: "TETRA_10",
    12: "PYRA_5",
    13: "PYRA_14",
    14: "PENTA_6",
    15: "PENTA_15",
    16: "PENTA_18",
    17: "HEXA_8",
    18: "HEXA_20",
    19: "HEXA_27",
    20: "MIXED",
    21: "PYRA_13",
    22: "NGON_n",
    23: "NFACE_n",
    24: "BAR_4",
    26: "TRI_10",
    28: "QUAD_16",
    30: "TETRA_20",
    36: "PENTA_40",
    39: "HEXA_64",
}


# FlowSolution_t node names. CGNS has NO component concept -- one DataArray_t
# is one scalar, with no NumberOfComponents anywhere in the SIDS -- so a k>1
# meshio++ array is split into k siblings suffixed `_0.._k-1` and re-joined on
# read. A documented meshio++ convention, not something SIDS specifies; see
# doc/formats/cgns.md. Twin of cgns.cpp's kCgnsVertexSolution/kCgnsCellSolution.
_VERTEX_SOLUTION = "FlowSolution"
_CELL_SOLUTION = "FlowSolutionCells"


def _component_name(base, k, i):
    """``<base>_<i>`` for a component of a multi-component array, else ``<base>``.

    The single owner of the suffix convention (C++ twin:
    ``cgns_component_name``).
    """
    return base if k <= 1 else f"{base}_{i}"


def _create_solution(zone, name, location):
    """A ``FlowSolution_t`` node plus its ``GridLocation_t`` child."""
    sol = _create_group(zone, name)
    _write_node_attrs(sol, name, "FlowSolution_t", "MT")
    gl = _create_group(sol, "GridLocation")
    _write_node_attrs(gl, "GridLocation", "GridLocation_t", "C1")
    # Like ZoneType's payload: the declared length carries the string, no NUL.
    gl.create_dataset(
        " data", data=np.frombuffer(_padded_bytes(location, len(location)), dtype="i1")
    )
    return sol


def _write_solution_array(sol, name, col, gzip_kwargs):
    """One ``DataArray_t`` holding a single scalar component."""
    g = _create_group(sol, name)
    _write_node_attrs(g, name, "DataArray_t", "R8")
    col = np.ascontiguousarray(col, dtype="<f8")
    kwargs = dict(gzip_kwargs)
    if kwargs and col.size:
        kwargs["chunks"] = col.shape
    g.create_dataset(" data", data=col, **kwargs)


def _strip_rind(arr, rind):
    """``arr`` without its rind planes, flat: ``Rind_t`` is ``(lo, hi)`` per index
    direction, ``i`` first, and the array's own axes are the reverse of that
    (``(nk, nj, ni)``); a 1-D array (an unstructured zone) uses the first pair.
    Twin of ``cgns_strip_rind``."""
    if not any(rind):
        return arr.reshape(-1)
    nd = max(arr.ndim, 1)
    arr = arr.reshape(arr.shape if arr.ndim else (1,))
    slices = []
    for ax in range(nd):
        pair = nd - 1 - ax
        lo = max(rind[2 * pair], 0) if 2 * pair + 1 < len(rind) else 0
        hi = max(rind[2 * pair + 1], 0) if 2 * pair + 1 < len(rind) else 0
        if lo + hi > arr.shape[ax]:
            raise ReadError("CGNS: a solution's Rind is larger than its array")
        slices.append(slice(lo, arr.shape[ax] - hi))
    return arr[tuple(slices)].reshape(-1)


def _read_solution(sol, sol_name, expected_rows):
    """Re-join a ``FlowSolution_t``'s ``DataArray_t`` children, undoing the
    ``_0.._k-1`` split. Returns a list of ``(name, (rows, k) array)`` pairs.

    Twin of ``cgns_read_solution``: a genuine multi-component array is a
    *contiguous* ``0..k-1`` run; anything else (a lone ``foo_7``, a gap) stays a
    scalar under its own literal name, since guessing would invent components.
    """
    raw_names = [
        child
        for child in sol
        if not child.startswith(" ")
        and child != "GridLocation"
        and _read_attr_str(sol[child], "label") == "DataArray_t"
    ]

    components = {}
    base_order = []
    for n in raw_names:
        base, idx = n, 0
        us = n.rfind("_")
        if us != -1 and us + 1 < len(n) and n[us + 1 :].isdigit():
            base, idx = n[:us], int(n[us + 1 :])
        if base not in components:
            base_order.append(base)
            components[base] = {}
        components[base][idx] = n

    rind = []
    if "Rind" in sol and " data" in sol["Rind"]:
        rind = [int(v) for v in np.asarray(sol["Rind"][" data"][()]).reshape(-1)]

    def _one(child):
        a = _strip_rind(np.asarray(sol[child][" data"][()], dtype=np.float64), rind)
        if a.size != expected_rows:
            raise ReadError(
                f"CGNS: '{sol_name}/{child}' has {a.size} values but the zone has "
                f"{expected_rows} entities at this GridLocation"
            )
        return a

    out = []
    for base in base_order:
        parts = components[base]
        idxs = sorted(parts)
        contiguous = idxs == list(range(len(idxs)))
        if not contiguous:
            for i in idxs:
                out.append((parts[i], _one(parts[i])))
            continue
        cols = [_one(parts[i]) for i in idxs]
        k = len(cols)
        out.append((base, cols[0] if k == 1 else np.stack(cols, axis=1)))
    return out


def _permute_conn(conn, cell_type, shift):
    """Node reorder + additive shift, the C++ `cgns_permute_conn` convention:
    shift=+1 writes (meshio++ to SIDS order), shift=-1 reads."""
    if shift > 0:
        conn = _node_order.from_meshio("cgns", cell_type, conn)
    else:
        conn = _node_order.to_meshio("cgns", cell_type, conn)
    return conn + shift


# ---- node encoding: fixed-length NULLTERM attributes + a " data" payload --


def _write_attr_str(obj, name, value, size):
    import h5py

    t = h5py.h5t.C_S1.copy()
    t.set_size(size)
    t.set_strpad(h5py.h5t.STR_NULLTERM)
    space = h5py.h5s.create(h5py.h5s.SCALAR)
    buf = value.encode("ascii")[:size].ljust(size, b"\0")
    attr_id = h5py.h5a.create(obj.id, name.encode("ascii"), t, space)
    attr_id.write(np.array(buf, dtype="S%d" % size))


def _write_attr_flags(obj):
    import h5py

    space = h5py.h5s.create_simple((1,))
    attr_id = h5py.h5a.create(obj.id, b"flags", h5py.h5t.STD_I32LE, space)
    attr_id.write(np.array([1], dtype="<i4"))


def _write_node_attrs(obj, name, label, dtype_code):
    _write_attr_str(obj, "name", name, 33)
    _write_attr_str(obj, "label", label, 33)
    _write_attr_str(obj, "type", dtype_code, 3)
    _write_attr_flags(obj)


def _create_group(parent, name):
    # track_order=True sets BOTH link and attribute creation-order tracking
    # on the new group's own GCPL -- load-bearing: cgnslib's has_child/
    # has_data (ADFH.c) iterate creation order with no name-order fallback.
    return parent.create_group(name, track_order=True)


def _type_code(dtype):
    dtype = np.dtype(dtype)
    if dtype == np.dtype("<i4"):
        return "I4"
    if dtype == np.dtype("<i8"):
        return "I8"
    if dtype == np.dtype("<f4"):
        return "R4"
    if dtype == np.dtype("<f8"):
        return "R8"
    raise WriteError(f"CGNS: unsupported dtype {dtype} for a CGNS node payload")


def _padded_bytes(s, total):
    b = s.encode("ascii")
    return b[:total].ljust(total, b"\0")


def _hdf5_version_string():
    import h5py

    maj, minr, rel = h5py.h5.get_libversion()
    return f"HDF5 Version {maj}.{minr}.{rel}"


def _read_attr_str(obj, name):
    if name not in obj.attrs:
        return None
    v = obj.attrs[name]
    if isinstance(v, bytes):
        return v.split(b"\0", 1)[0].decode("ascii", "ignore").rstrip(" ")
    return str(v).rstrip(" ")


# ---- write -----------------------------------------------------------------


_FIELD_DATA = "FieldData"
_SUBREGION_NAME = "RegionName"
_SUBREGION_TAG = "RegionTag"


def _valid_node_name(name):
    """Whether CGNS can name a node ``name`` (1 to 32 characters, no ``/``)."""
    return 0 < len(name) <= 32 and "/" not in name and not name.startswith(" ")


def _c1_node(parent, name, label, text):
    node = _create_group(parent, name)
    _write_node_attrs(node, name, label, "C1")
    node.create_dataset(
        " data", data=np.frombuffer(_padded_bytes(text, len(text)), "i1")
    )
    return node


def _read_c1(node):
    raw = np.asarray(node[" data"][()]).astype("i1")
    return bytes(raw).decode("latin-1").strip("\0 ")


_UNIT_SLOTS = 5  # Mass, Length, Time, Temperature, Angle
_UNIT_WIDTH = 32


def _write_dimensional_units(base, mesh):
    """The mesh's length unit as the base's ``DimensionalUnits_t`` (Length
    named, the other four ``Null``); False when it has none CGNS can name."""
    unit = length_unit(mesh)
    if unit is None:
        return False
    name = unit_name(unit, CGNS_LENGTH_UNITS)
    if name is None:
        warnings.warn(
            f"CGNS: a length unit of {unit!r} m has no DimensionalUnits name "
            "(Meter, Centimeter, Millimeter, Foot, Inch); kept as field data only.",
            stacklevel=3,
        )
        return False
    rows = ["Null", name, "Null", "Null", "Null"]
    data = np.frombuffer(
        b"".join(r.encode("ascii").ljust(_UNIT_WIDTH, b" ") for r in rows), "i1"
    ).reshape(_UNIT_SLOTS, _UNIT_WIDTH)
    node = _create_group(base, "DimensionalUnits")
    _write_node_attrs(node, "DimensionalUnits", "DimensionalUnits_t", "C1")
    node.create_dataset(" data", data=data)
    return True


def _read_dimensional_units(node):
    """The metres per unit a parent's ``DimensionalUnits_t`` child names, or
    ``None``."""
    for child in node:
        if child.startswith(" "):
            continue
        if _read_attr_str(node[child], "label") != "DimensionalUnits_t":
            continue
        if " data" not in node[child]:
            return None
        raw = np.asarray(node[child][" data"][()]).astype("i1").reshape(-1)
        if raw.size < 2 * _UNIT_WIDTH:
            return None
        text = bytes(raw[_UNIT_WIDTH : 2 * _UNIT_WIDTH]).decode("latin-1")
        return CGNS_LENGTH_UNITS.get(text.strip("\0 "))
    return None


def _write_field_data(zone, mesh, skip=()):
    """``field_data`` as ``DataArray_t`` children of a ``UserDefinedData_t`` node
    (the arrays keep their shape; CGNS has no per-zone global variable node)."""
    names = []
    for name in sorted(mesh.field_data):
        if name in skip:
            continue
        a = np.asarray(mesh.field_data[name])
        if not _valid_node_name(name) or a.dtype.kind not in "iuf" or a.size == 0:
            warnings.warn(
                f"CGNS: field_data '{name}' cannot be a CGNS node (name over 32 "
                "characters or with '/', or not numeric); not written.",
                stacklevel=3,
            )
            continue
        kind = "f" if a.dtype.kind == "f" else "i"
        dtype = np.dtype(
            ("<f4" if a.dtype.itemsize <= 4 else "<f8")
            if kind == "f"
            else ("<i4" if a.dtype.itemsize <= 4 and a.dtype.kind == "i" else "<i8")
        )
        names.append((name, a.astype(dtype).reshape(a.shape if a.ndim else (1,))))
    if not names:
        return
    group = _create_group(zone, _FIELD_DATA)
    _write_node_attrs(group, _FIELD_DATA, "UserDefinedData_t", "MT")
    for name, a in names:
        g = _create_group(group, name)
        _write_node_attrs(g, name, "DataArray_t", _type_code(a.dtype))
        g.create_dataset(" data", data=a)


def _read_field_data(zone):
    out = {}
    for child in zone:
        if (
            child.startswith(" ")
            or _read_attr_str(zone[child], "label") != "UserDefinedData_t"
        ):
            continue
        group = zone[child]
        if child != _FIELD_DATA:
            continue
        for name in group:
            if (
                name.startswith(" ")
                or _read_attr_str(group[name], "label") != "DataArray_t"
            ):
                continue
            if " data" in group[name]:
                out[name] = np.asarray(group[name][" data"][()])
    return out


def _write_regions(zone, mesh, starts, firsts, n_points):
    """Point and cell regions as ``ZoneSubRegion_t`` nodes (``GridLocation`` plus a
    1-based ``PointList``; a cell region lists element ids). The node value is the
    region's dimension and a ``RegionTag`` child holds its tag. Side regions and
    cells of face-based sections have no CGNS node here and are dropped."""
    used = set()
    skipped_sides = 0
    for k, region in enumerate(canonical_regions(mesh)):
        if region.kind == "side":
            skipped_sides += 1
            continue
        ids = np.asarray(region.entries, dtype=np.int64)
        if region.kind == "point":
            ids = ids + 1
            location = "Vertex"
        else:
            location = "CellCenter"
            block = np.searchsorted(starts[1:], ids, side="right")
            ok = (block < len(firsts)) & (
                firsts[np.minimum(block, len(firsts) - 1)] >= 0
            )
            if not ok.all():
                warnings.warn(
                    f"CGNS: region '{region.name}' holds cells with no element "
                    "section of their own (face-based blocks); they are not written.",
                    stacklevel=3,
                )
            ids = firsts[block[ok]] + (ids[ok] - starts[block[ok]])
        node_name = (
            region.name
            if _valid_node_name(region.name) and region.name not in used
            else f"Region{k}"
        )
        used.add(node_name)
        wide = n_points > np.iinfo(np.int32).max or (
            len(ids) and ids.max() > np.iinfo(np.int32).max
        )
        itype = np.dtype("<i8") if wide else np.dtype("<i4")
        g = _create_group(zone, node_name)
        _write_node_attrs(g, node_name, "ZoneSubRegion_t", "I4")
        g.create_dataset(" data", data=np.array([region.dim], dtype="<i4"))
        gl = _create_group(g, "GridLocation")
        _write_node_attrs(gl, "GridLocation", "GridLocation_t", "C1")
        gl.create_dataset(
            " data", data=np.frombuffer(_padded_bytes(location, len(location)), "i1")
        )
        pl = _create_group(g, "PointList")
        _write_node_attrs(pl, "PointList", "IndexArray_t", _type_code(itype))
        pl.create_dataset(" data", data=ids.astype(itype).reshape(-1, 1))
        tag = _create_group(g, _SUBREGION_TAG)
        _write_node_attrs(tag, _SUBREGION_TAG, "DataArray_t", "I4")
        tag.create_dataset(" data", data=np.array([region.tag], dtype="<i4"))
        if node_name != region.name:
            _c1_node(g, _SUBREGION_NAME, "DataArray_t", region.name)
    if skipped_sides:
        warnings.warn(
            f"CGNS: {skipped_sides} side region(s) have no ZoneSubRegion equivalent; "
            "not written.",
            stacklevel=3,
        )


def _read_regions(zone, sections, n_points, n_cells_before, point_offset=0):
    """The ``ZoneSubRegion_t`` children as regions. ``sections`` are the zone's
    ``(first element id, last, first cell index)`` triples."""
    regions = []
    dropped = 0
    for child in zone:
        if (
            child.startswith(" ")
            or _read_attr_str(zone[child], "label") != "ZoneSubRegion_t"
        ):
            continue
        node = zone[child]
        if "PointList" not in node or " data" not in node["PointList"]:
            warnings.warn(
                f"CGNS: ZoneSubRegion '{child}' has no PointList; not read.",
                stacklevel=3,
            )
            continue
        location = "Vertex"
        if "GridLocation" in node and " data" in node["GridLocation"]:
            location = _read_c1(node["GridLocation"])
        ids = np.asarray(node["PointList"][" data"][()], dtype=np.int64).reshape(-1)
        dim = -1
        if " data" in node:
            value = np.asarray(node[" data"][()]).reshape(-1)
            dim = int(value[0]) if value.size else -1
        tag = -1
        if _SUBREGION_TAG in node and " data" in node[_SUBREGION_TAG]:
            tag = int(np.asarray(node[_SUBREGION_TAG][" data"][()]).reshape(-1)[0])
        name = child
        if _SUBREGION_NAME in node and " data" in node[_SUBREGION_NAME]:
            name = _read_c1(node[_SUBREGION_NAME])
        if location == "Vertex":
            if len(ids) and (ids.min() < 1 or ids.max() > n_points):
                raise ReadError(
                    f"CGNS: ZoneSubRegion '{child}' names a vertex the zone does not have"
                )
            regions.append(
                Region(name, "point", ids - 1 + point_offset, dim=dim, tag=tag)
            )
        elif location == "CellCenter":
            cells = []
            for e in ids:
                for first, last, start in sections:
                    if first <= e <= last:
                        cells.append(n_cells_before + start + (e - first))
                        break
                else:
                    dropped += 1
            regions.append(
                Region(name, "cell", np.array(cells, dtype=np.int64), dim=dim, tag=tag)
            )
        else:
            warnings.warn(
                f"CGNS: ZoneSubRegion '{child}' has GridLocation '{location}'; only "
                "Vertex and CellCenter are supported, so it was not read.",
                stacklevel=3,
            )
    if dropped:
        warnings.warn(
            f"CGNS: {dropped} region cell(s) name elements outside the zone's "
            "fixed-size sections; dropped.",
            stacklevel=3,
        )
    return regions


def _expand_range(raw, dims):
    """The 1-based point ids of a ``PointRange``: a one-dimensional ``(first, last)``
    pair, or, in a structured zone of ``dims`` points, one pair per index direction
    (``i`` fastest). ``None`` for a range this reader cannot place. Twin of
    ``cgns_expand_range``."""
    rng = np.asarray(raw, dtype=np.int64).reshape(-1)
    if rng.size == 2:
        return np.arange(rng[0], rng[1] + 1, dtype=np.int64)
    if dims is None or rng.size % 2 or rng.size > 6:
        return None
    lo = [1, 1, 1]
    hi = [1, 1, 1]
    for a in range(rng.size // 2):
        lo[a], hi[a] = int(rng[2 * a]), int(rng[2 * a + 1])
    k, j, i = np.meshgrid(
        np.arange(lo[2], hi[2] + 1),
        np.arange(lo[1], hi[1] + 1),
        np.arange(lo[0], hi[0] + 1),
        indexing="ij",
    )
    return (i + dims[0] * ((j - 1) + dims[1] * (k - 1))).reshape(-1).astype(np.int64)


def _read_interfaces(zone, n_points, point_offset, dims):
    """The ``GridConnectivity1to1_t`` children of a zone as point regions named
    ``cgns:interface:<name>``, and the donor zone's name as
    ``field_data["cgns:interface_donor:<name>"]`` (``uint8``). Twin of
    ``cgns_read_interfaces``."""
    regions, field_data = [], {}
    for child in zone:
        if (
            child.startswith(" ")
            or _read_attr_str(zone[child], "label") != "ZoneGridConnectivity_t"
        ):
            continue
        for name in zone[child]:
            if name.startswith(" "):
                continue
            node = zone[child][name]
            if _read_attr_str(node, "label") != "GridConnectivity1to1_t":
                continue
            if "PointList" in node and " data" in node["PointList"]:
                ids = np.asarray(node["PointList"][" data"][()], dtype=np.int64)
                ids = ids.reshape(-1)
            elif "PointRange" in node and " data" in node["PointRange"]:
                ids = _expand_range(node["PointRange"][" data"][()], dims)
            else:
                ids = None
            if ids is None or not len(ids):
                warnings.warn(
                    f"CGNS: interface '{name}' lists no points this reader can "
                    "place; not read.",
                    stacklevel=3,
                )
                continue
            if ids.min() < 1 or ids.max() > n_points:
                raise ReadError(
                    f"CGNS: interface '{name}' names a vertex the zone does not have"
                )
            regions.append(
                Region(f"cgns:interface:{name}", "point", ids - 1 + point_offset)
            )
            if " data" in node:
                donor = _read_c1(node)
                if donor:
                    field_data[f"cgns:interface_donor:{name}"] = np.frombuffer(
                        donor.encode("latin-1"), np.uint8
                    )
    return regions, field_data


def _merge_regions(regions):
    """Regions of several zones as one list: those with the same name and kind
    merge, their entries joined."""
    merged = {}
    for r in regions:
        key = (r.name, r.kind)
        if key in merged:
            old = merged[key]
            merged[key] = Region(
                r.name,
                r.kind,
                np.concatenate([old.entries, r.entries]),
                dim=old.dim,
                tag=old.tag,
            )
        else:
            merged[key] = r
    return list(merged.values())


def _read_bcs(zone, sections, n_points, n_cells_before, point_offset=0, dims=None):
    """The ``ZoneBC_t/BC_t`` children as regions, and the BC types as
    ``field_data["cgns:bc_type:<region>"]`` (the type's characters as ``uint8``).

    A BC lists its points (``PointList``) or their range (``PointRange``) at its
    ``GridLocation``: ``Vertex`` names points, ``FaceCenter``, ``EdgeCenter`` and
    ``CellCenter`` name elements (the boundary elements of the file's own
    sections), so they become cell regions over those blocks. The region is named
    after the BC's ``FamilyName`` when it has one, else after the BC, and BCs
    sharing a name and kind are merged.
    """
    merged = {}  # (name, kind) -> entries
    types = {}
    dropped = 0
    for child in zone:
        if child.startswith(" ") or _read_attr_str(zone[child], "label") != "ZoneBC_t":
            continue
        for bc_name in zone[child]:
            if bc_name.startswith(" "):
                continue
            node = zone[child][bc_name]
            if _read_attr_str(node, "label") != "BC_t":
                continue
            if "PointList" in node and " data" in node["PointList"]:
                ids = np.asarray(node["PointList"][" data"][()], dtype=np.int64)
                ids = ids.reshape(-1)
            elif "PointRange" in node and " data" in node["PointRange"]:
                ids = _expand_range(node["PointRange"][" data"][()], dims)
                if ids is None:
                    warnings.warn(
                        f"CGNS: BC '{bc_name}' has a PointRange this reader cannot "
                        "place; not read.",
                        stacklevel=3,
                    )
                    continue
            else:
                warnings.warn(
                    f"CGNS: BC '{bc_name}' has no PointList or PointRange; not read.",
                    stacklevel=3,
                )
                continue
            location = "Vertex"
            if "GridLocation" in node and " data" in node["GridLocation"]:
                location = _read_c1(node["GridLocation"])
            name = bc_name
            if "FamilyName" in node and " data" in node["FamilyName"]:
                name = _read_c1(node["FamilyName"]) or bc_name
            if location == "Vertex":
                if len(ids) and (ids.min() < 1 or ids.max() > n_points):
                    raise ReadError(
                        f"CGNS: BC '{bc_name}' names a vertex the zone does not have"
                    )
                kind, entries = "point", ids - 1 + point_offset
            else:
                kind, cells = "cell", []
                for e in ids:
                    for first, last, start in sections:
                        if first <= e <= last:
                            cells.append(n_cells_before + start + (e - first))
                            break
                    else:
                        dropped += 1
                entries = np.array(cells, dtype=np.int64)
            merged.setdefault((name, kind), []).append(entries)
            if " data" in node:
                bc_type = _read_c1(node)
                if bc_type:
                    types[name] = np.frombuffer(bc_type.encode("latin-1"), np.uint8)
    if dropped:
        warnings.warn(
            f"CGNS: {dropped} BC element(s) name elements outside the zone's "
            "fixed-size sections; dropped.",
            stacklevel=3,
        )
    regions = [
        Region(name, kind, np.concatenate(parts))
        for (name, kind), parts in merged.items()
    ]
    return regions, {f"cgns:bc_type:{n}": v for n, v in types.items()}


@_provenance.slotless_writer
def write(filename, mesh, compression="gzip", compression_opts=4):
    import h5py

    gzip_kwargs = {}
    if compression not in (None, False):
        gzip_kwargs = {
            "compression": "gzip",
            "compression_opts": int(compression_opts or 4),
        }

    point_dim = mesh.points.shape[1] if mesh.points.ndim == 2 else 3
    if point_dim > 3:
        raise WriteError(f"CGNS: PhysicalDimension must be 1..3, got {point_dim}")

    # Validate every cell type up front (before any file is created) and
    # compute CellDim = the max topological dimension over all blocks.
    cell_dim = 0
    for cb in mesh.cells:
        # polygon/polyhedron* blocks are always ragged in the C++ core (see
        # CellBlock::mPolygonRows/mPolyhedronRows in the MESHIO backend)
        # regardless of whether this particular block happens to be uniform
        # -- checked by type name, not by whether `.data` is a Python list,
        # since a uniform polygon block is stored as a plain ndarray.
        if cb.type.startswith("polygon") or cb.type.startswith("polyhedron"):
            # NGON_n/NFACE_n write is C++-core-only (v9.21.0). Deliberately not
            # reimplemented here: the writer deduplicates faces and repairs each
            # cell's winding, and the repair is a discrete branch on the sign of
            # an enclosed volume -- two implementations of such a branch can
            # land on opposite sides for a near-degenerate cell and then diverge
            # macroscopically, the same reasoning that keeps `smooth`'s
            # inversion guard out of its numpy fallback. It is also unreachable
            # in practice: `_core` ships in every wheel.
            raise WriteError(
                f"CGNS: cell type '{cb.type}' needs an NGON_n/NFACE_n section, "
                "which the pure-Python writer does not produce; this requires "
                "the compiled core (meshioplusplus._core.cgns_write)"
            )
        if cb.type not in _CGNS_TYPES:
            if cb.type in topological_dimension:
                raise WriteError(
                    f"CGNS: cell type '{cb.type}' maps to a CGNS ElementType_t "
                    "but its CGNS node ordering is not yet verified in "
                    "meshio++; refusing to write a guessed ordering"
                )
            raise WriteError(
                f"CGNS: cell type '{cb.type}' has no fixed-size CGNS "
                "ElementType_t equivalent"
            )
        cell_dim = max(cell_dim, topological_dimension.get(cb.type, 0))

    phys_dim = min(3, max(1, max(point_dim, cell_dim)))
    cgns_cell_dim = cell_dim if cell_dim > 0 else phys_dim

    n_cells_at_dim = sum(
        len(cb.data)
        for cb in mesh.cells
        if topological_dimension.get(cb.type, -1) == cgns_cell_dim
    )
    n_points = mesh.points.shape[0]

    f = h5py.File(filename, "w", track_order=True)
    _write_attr_str(f, "name", "HDF5 MotherNode", 33)
    _write_attr_str(f, "label", "Root Node of HDF5 File", 33)
    _write_attr_str(f, "type", "MT", 3)
    f.create_dataset(
        " format", data=np.frombuffer(_padded_bytes("IEEE_LITTLE_32", 15), dtype="i1")
    )
    f.create_dataset(
        " hdf5version",
        data=np.frombuffer(_padded_bytes(_hdf5_version_string(), 33), dtype="i1"),
    )

    cgver = _create_group(f, "CGNSLibraryVersion")
    _write_node_attrs(cgver, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4")
    cgver.create_dataset(" data", data=np.array([3.1], dtype="<f4"))

    base = _create_group(f, "Base")
    _write_node_attrs(base, "Base", "CGNSBase_t", "I4")
    base.create_dataset(" data", data=np.array([cgns_cell_dim, phys_dim], dtype="<i4"))
    units_written = _write_dimensional_units(base, mesh)

    wide = n_points > np.iinfo(np.int32).max or n_cells_at_dim > np.iinfo(np.int32).max
    zone_np_dt = np.dtype("<i8") if wide else np.dtype("<i4")

    zone = _create_group(base, "Zone1")
    _write_node_attrs(zone, "Zone1", "Zone_t", _type_code(zone_np_dt))
    zdata = np.array([[n_points], [n_cells_at_dim], [0]], dtype=zone_np_dt)
    zone.create_dataset(" data", data=zdata)

    zt = _create_group(zone, "ZoneType")
    _write_node_attrs(zt, "ZoneType", "ZoneType_t", "C1")
    zt.create_dataset(
        " data", data=np.frombuffer(_padded_bytes("Unstructured", 12), dtype="i1")
    )

    coords = _create_group(zone, "GridCoordinates")
    _write_node_attrs(coords, "GridCoordinates", "GridCoordinates_t", "MT")
    coord_dtype = "<f4" if mesh.points.dtype == np.float32 else "<f8"
    names = ["CoordinateX", "CoordinateY", "CoordinateZ"]
    n_coords = 3 if phys_dim >= 3 else 2
    for c in range(n_coords):
        g = _create_group(coords, names[c])
        col = (
            mesh.points[:, c].astype(coord_dtype)
            if c < point_dim
            else np.zeros(n_points, dtype=coord_dtype)
        )
        _write_node_attrs(g, names[c], "DataArray_t", _type_code(coord_dtype))
        kwargs = dict(gzip_kwargs)
        if kwargs and col.size:
            kwargs["chunks"] = col.shape
        g.create_dataset(" data", data=col, **kwargs)

    # Elements: one Elements_t section per meshio++ cell block, in mesh
    # order, with contiguous non-overlapping 1-based ElementRanges. Unlike
    # MED, sections of the same type are NOT consolidated.
    next_id = 1
    type_counts = {}
    region_starts = [0]  # per block: its first cell index, and its first element id
    region_firsts = []
    for cb in mesh.cells:
        conn = np.asarray(cb.data)
        nc = conn.shape[0]
        region_starts.append(region_starts[-1] + nc)
        region_firsts.append(next_id if nc else -1)
        if nc == 0:
            continue  # a zero-length ElementRange is not representable
        cgns_name, code = _CGNS_TYPES[cb.type]

        first, last = next_id, next_id + nc - 1
        next_id = last + 1

        type_counts[cgns_name] = type_counts.get(cgns_name, 0) + 1
        section_name = f"{cgns_name}_{type_counts[cgns_name]}"

        sect = _create_group(zone, section_name)
        _write_node_attrs(sect, section_name, "Elements_t", "I4")
        sect.create_dataset(" data", data=np.array([code, 0], dtype="<i4"))

        section_wide = conn.dtype != np.int32 or n_points > np.iinfo(np.int32).max
        out_dtype = np.dtype("<i8") if section_wide else np.dtype("<i4")

        rng = _create_group(sect, "ElementRange")
        _write_node_attrs(rng, "ElementRange", "IndexRange_t", _type_code(out_dtype))
        rng.create_dataset(" data", data=np.array([first, last], dtype=out_dtype))

        permuted = _permute_conn(conn.astype(out_dtype), cb.type, 1)
        flat = permuted.reshape(-1)

        ec = _create_group(sect, "ElementConnectivity")
        _write_node_attrs(
            ec, "ElementConnectivity", "DataArray_t", _type_code(out_dtype)
        )
        kwargs = dict(gzip_kwargs)
        if kwargs and flat.size:
            kwargs["chunks"] = flat.shape
        ec.create_dataset(" data", data=flat, **kwargs)

    # FlowSolution_t: point_data at "Vertex", cell_data at "CellCenter". Twin of
    # the C++ writer -- see cgns.cpp's own comment for why a k>1 array is split
    # into `<name>_<i>` siblings (CGNS has no NumberOfComponents) and why a
    # mixed-dimension mesh's cell_data cannot be written at all. field_data has
    # no CGNS home and is not written.
    if mesh.point_data:
        sol = _create_solution(zone, _VERTEX_SOLUTION, "Vertex")
        for name in sorted(mesh.point_data):
            arr = np.asarray(mesh.point_data[name], dtype=np.float64)
            rows = arr.shape[0] if arr.ndim else 0
            k = int(np.prod(arr.shape[1:])) if arr.ndim > 1 else 1
            if rows != n_points:
                warnings.warn(
                    f"CGNS: point_data '{name}' has {rows} rows but the mesh has "
                    f"{n_points} points; not written.",
                    stacklevel=2,
                )
                continue
            flat2 = arr.reshape(rows, k)
            for i in range(k):
                _write_solution_array(
                    sol, _component_name(name, k, i), flat2[:, i], gzip_kwargs
                )

    if mesh.cell_data:
        all_at_cell_dim = len(mesh.cells) > 0 and all(
            topological_dimension.get(cb.type, -1) == cgns_cell_dim for cb in mesh.cells
        )
        if not all_at_cell_dim:
            warnings.warn(
                "CGNS: this mesh mixes cell blocks of different topological dimensions, "
                "so a zone-wide CellCenter FlowSolution cannot be distributed back "
                f"across them; {len(mesh.cell_data)} cell_data array(s) not written.",
                stacklevel=2,
            )
            _provenance.note(
                "data-dropped",
                f"{len(mesh.cell_data)} cell_data array(s) not written -- CGNS's "
                "CellCenter FlowSolution is per-zone and this mesh mixes topological "
                "dimensions",
            )
        else:
            sol = _create_solution(zone, _CELL_SOLUTION, "CellCenter")
            for name in sorted(mesh.cell_data):
                blocks = mesh.cell_data[name]
                if len(blocks) != len(mesh.cells):
                    warnings.warn(
                        f"CGNS: cell_data '{name}' covers {len(blocks)} of "
                        f"{len(mesh.cells)} blocks; not written.",
                        stacklevel=2,
                    )
                    continue
                arrs = [np.asarray(b, dtype=np.float64) for b in blocks]
                k = int(np.prod(arrs[0].shape[1:])) if arrs[0].ndim > 1 else 1
                ok = all(
                    (int(np.prod(a.shape[1:])) if a.ndim > 1 else 1) == k
                    and a.shape[0] == len(cb.data)
                    for a, cb in zip(arrs, mesh.cells)
                )
                merged = (
                    np.concatenate([a.reshape(a.shape[0], k) for a in arrs], axis=0)
                    if ok
                    else None
                )
                if not ok or merged.shape[0] != n_cells_at_dim:
                    warnings.warn(
                        f"CGNS: cell_data '{name}' does not line up with the zone's "
                        "cells; not written.",
                        stacklevel=2,
                    )
                    continue
                for i in range(k):
                    _write_solution_array(
                        sol, _component_name(name, k, i), merged[:, i], gzip_kwargs
                    )

    _write_field_data(zone, mesh, skip=(LENGTH_UNIT_KEY,) if units_written else ())
    _write_regions(
        zone,
        mesh,
        np.array(region_starts, dtype=np.int64),
        np.array(region_firsts, dtype=np.int64),
        n_points,
    )

    f.close()


# ---- read --------------------------------------------------------------


def read(filename, base=None):
    import h5py

    f = h5py.File(filename, "r")

    if _is_spec_layout(f):
        return _read_spec(f, base)
    return _read_legacy(f)


def _is_spec_layout(f):
    for name in f:
        if name.startswith(" "):
            continue
        obj = f[name]
        if _read_attr_str(obj, "label") == "CGNSBase_t":
            return True
    return False


def _join_zone_fields(zone_fields):
    """The per-zone solutions as one ``point_data`` and ``cell_data``: each array
    concatenated over the zones in block order, a zone without it filled with NaN
    (``merge``'s ``Fill`` policy)."""
    point_data, cell_data = {}, {}
    names = []
    for _, _, pd, _cd in zone_fields:
        names += [n for n in pd if n not in names]
    for name in names:
        template = next(pd[name] for _, _, pd, _ in zone_fields if name in pd)
        parts = [
            (
                pd[name]
                if name in pd
                else np.full((n_points,) + template.shape[1:], np.nan)
            )
            for n_points, _, pd, _ in zone_fields
        ]
        point_data[name] = parts[0] if len(parts) == 1 else np.concatenate(parts)
    names = []
    for _, _, _pd, cd in zone_fields:
        names += [n for n in cd if n not in names]
    for name in names:
        template = next(cd[name][0] for _, _, _, cd in zone_fields if name in cd)
        blocks = []
        for _, sizes, _, cd in zone_fields:
            if name in cd:
                blocks += cd[name]
            else:
                blocks += [np.full((n,) + template.shape[1:], np.nan) for n in sizes]
        cell_data[name] = blocks
    return point_data, cell_data


def _read_spec(f, base=None):
    base_names = [
        name
        for name in f
        if not name.startswith(" ") and _read_attr_str(f[name], "label") == "CGNSBase_t"
    ]
    n_bases = len(base_names)
    base_name = base_names[0] if base_names else None
    if base is not None:
        # A base is chosen by name, else by its position (0 is the first).
        if str(base) in base_names:
            base_name = str(base)
        elif str(base).isdigit() and int(base) < n_bases:
            base_name = base_names[int(base)]
        else:
            raise ReadError(
                f"CGNS: the file has no base '{base}' ({n_bases} base(s) found)"
            )
    elif n_bases > 1:
        warnings.warn(
            f"CGNS: file has {n_bases} CGNSBase_t nodes; only the first "
            f"('{base_name}') is read; pass a base to choose.",
            stacklevel=2,
        )
    base = f[base_name]
    unit = _read_dimensional_units(base)

    point_chunks = []
    cells = []
    zone_fields = []  # per zone: points, block sizes, point data, cell data
    field_data = {}
    regions = []
    point_offset = 0
    point_dim_out = 3

    for zname in base:
        if zname.startswith(" "):
            continue
        zone = base[zname]
        if _read_attr_str(zone, "label") != "Zone_t":
            continue

        zone_unit = _read_dimensional_units(zone)
        if zone_unit is not None:
            unit = zone_unit
        # A Structured zone is expanded to its lattice: ``Zone_t``'s value is the
        # (vertex, cell, boundary-vertex) sizes, one row each, ``i`` the fastest index.
        structured_dims = None
        if "ZoneType" in zone and " data" in zone["ZoneType"]:
            zt_raw = zone["ZoneType"][" data"][()]
            zt_name = bytes(np.asarray(zt_raw).astype("i1")).decode("ascii", "ignore")
            if zt_name == "Structured":
                zdata = (
                    np.asarray(zone[" data"][()]).reshape(-1)
                    if " data" in zone
                    else np.empty(0, dtype=np.int64)
                )
                idim = zdata.size // 3
                if not 1 <= idim <= 3:
                    raise ReadError(
                        f"CGNS: structured zone '{zname}' has a malformed size record"
                    )
                structured_dims = [int(v) for v in zdata[:idim]] + [1] * (3 - idim)
                if min(structured_dims) < 1:
                    raise ReadError(f"CGNS: structured zone '{zname}' has no points")
            elif zt_name != "Unstructured":
                raise ReadError(
                    f"CGNS: zone '{zname}' has ZoneType '{zt_name}'; only "
                    "Unstructured and Structured zones are supported."
                )

        coords = zone["GridCoordinates"]
        axes = [
            ax for ax in ("CoordinateX", "CoordinateY", "CoordinateZ") if ax in coords
        ]
        if not axes:
            raise ReadError(f"CGNS: zone '{zname}' has no GridCoordinates")
        point_dim_out = max(2, len(axes))
        cols = [np.asarray(coords[ax][" data"][()], dtype=np.float64) for ax in axes]
        if structured_dims is not None:
            cols = [c.reshape(-1) for c in cols]
        n_zone_points = cols[0].shape[0] if cols[0].ndim else 0
        zpts = np.zeros((n_zone_points, point_dim_out), dtype=np.float64)
        for d in range(point_dim_out):
            if d < len(cols):
                zpts[:, d] = cols[d]
        point_chunks.append(zpts)

        # Cell counts of the blocks this zone contributes, in emission order --
        # what a zone-wide CellCenter FlowSolution is split back across.
        zone_block_cells = []
        zone_sections = []  # (first element id, last, first cell index in the zone)
        cells_before = sum(len(c[1]) for c in cells)

        if structured_dims is not None:
            expected = structured_dims[0] * structured_dims[1] * structured_dims[2]
            if n_zone_points != expected:
                raise ReadError(
                    f"CGNS: structured zone '{zname}' lists {n_zone_points} "
                    f"coordinates for a {' x '.join(map(str, structured_dims))} lattice"
                )
            cell_type, lattice = lattice_cells(*structured_dims)
            if len(lattice):
                cells.append((cell_type, lattice.astype(np.int64) + point_offset))
                zone_block_cells.append(len(lattice))

        sects = []
        for sname in zone:
            if sname.startswith(" "):
                continue
            s = zone[sname]
            if _read_attr_str(s, "label") != "Elements_t":
                continue
            rng = np.asarray(s["ElementRange"][" data"][()])
            if rng.size < 2:
                raise ReadError(f"CGNS: section '{sname}' has a malformed ElementRange")
            sects.append((int(rng.flat[0]), sname))
        sects.sort()

        for _first_unused, sname in sects:
            s = zone[sname]
            sdata = np.asarray(s[" data"][()]).reshape(-1)
            if sdata.size < 1:
                raise ReadError(
                    f"CGNS: section '{sname}' has a malformed Elements_t descriptor"
                )
            code = int(sdata[0])

            if code == 20:
                raise ReadError(
                    f"CGNS: element section '{sname}' has ElementType "
                    f"{_CODE_TO_CGNS_NAME.get(code, '?')} ({code}); MIXED "
                    "sections are not supported."
                )
            if code in (22, 23):
                # The C++ core reads these (v9.21.0); the pure-Python twin does
                # not, so say which path is needed rather than claiming the
                # format is unsupported.
                raise ReadError(
                    f"CGNS: element section '{sname}' has ElementType "
                    f"{_CODE_TO_CGNS_NAME.get(code, '?')} ({code}); polyhedral "
                    "sections need the compiled core "
                    "(meshioplusplus._core.cgns_read)"
                )
            meshio_type = _CODE_TO_TYPE.get(code)
            if meshio_type is None:
                name = _CODE_TO_CGNS_NAME.get(code)
                if name is not None:
                    raise ReadError(
                        f"CGNS: element section '{sname}' has type {name} ({code}), "
                        "whose node ordering is not yet verified in meshio++."
                    )
                raise ReadError(
                    f"CGNS: element section '{sname}' has unknown ElementType "
                    f"code {code}."
                )
            cgns_name, _code = _CGNS_TYPES[meshio_type]

            rng = np.asarray(s["ElementRange"][" data"][()]).reshape(-1)
            first, last = int(rng[0]), int(rng[1])
            if last < first:
                raise ReadError(
                    f"CGNS: section '{sname}' has an empty/overlapping "
                    f"ElementRange [{first}, {last}]"
                )
            nc = last - first + 1
            flat = np.asarray(s["ElementConnectivity"][" data"][()]).reshape(-1)
            expected_npc = num_nodes_per_cell[meshio_type]
            if flat.size != nc * expected_npc:
                raise ReadError(
                    f"CGNS: section '{sname}' declares {nc} elements of "
                    f"{cgns_name} ({nc * expected_npc} nodes) but "
                    f"ElementConnectivity has {flat.size} entries"
                )

            conn = flat.reshape(nc, expected_npc)
            conn = _permute_conn(conn, meshio_type, -1)
            if point_offset:
                conn = conn + point_offset
            zone_sections.append((first, last, sum(zone_block_cells)))
            cells.append((meshio_type, conn.astype(np.int64)))
            zone_block_cells.append(nc)

        # FlowSolution_t (see the writer). GridLocation absent => "Vertex", the
        # SIDS default.
        sol_names = [
            child
            for child in zone
            if not child.startswith(" ")
            and _read_attr_str(zone[child], "label") == "FlowSolution_t"
        ]
        zone_pd, zone_cd = {}, {}
        if sol_names:
            zone_total_cells = sum(zone_block_cells)
            for child in sol_names:
                sol = zone[child]
                location = "Vertex"
                if "GridLocation" in sol and " data" in sol["GridLocation"]:
                    raw = np.asarray(sol["GridLocation"][" data"][()]).astype("i1")
                    location = bytes(raw).decode("ascii", "ignore").strip("\0 ")
                if location == "Vertex":
                    for name, arr in _read_solution(sol, child, n_zone_points):
                        zone_pd[name] = arr
                elif location == "CellCenter":
                    for name, arr in _read_solution(sol, child, zone_total_cells):
                        # Split the zone-wide array back across the blocks.
                        blocks, row = [], 0
                        for nb in zone_block_cells:
                            blocks.append(arr[row : row + nb])
                            row += nb
                        zone_cd[name] = blocks
                else:
                    warnings.warn(
                        f"CGNS: FlowSolution '{child}' has GridLocation '{location}'; "
                        "only Vertex and CellCenter are supported, so it was not read.",
                        stacklevel=2,
                    )

        zone_dims = structured_dims
        regions.extend(
            _read_regions(
                zone, zone_sections, n_zone_points, cells_before, point_offset
            )
        )
        bc_regions, bc_types = _read_bcs(
            zone, zone_sections, n_zone_points, cells_before, point_offset, zone_dims
        )
        regions.extend(bc_regions)
        field_data.update(bc_types)
        if_regions, if_fields = _read_interfaces(
            zone, n_zone_points, point_offset, zone_dims
        )
        regions.extend(if_regions)
        field_data.update(if_fields)
        field_data.update(_read_field_data(zone))
        zone_fields.append((n_zone_points, list(zone_block_cells), zone_pd, zone_cd))
        point_offset += n_zone_points

    point_data, cell_data = _join_zone_fields(zone_fields)
    if not point_chunks:
        raise ReadError(f"CGNS: base '{base_name}' has no Unstructured zones")
    if len(point_chunks) == 1:
        points = point_chunks[0]
    else:
        # Zones may differ in coordinate count; the widest wins and the
        # others pad with zeros, as the core reader does.
        dim = max(c.shape[1] for c in point_chunks)
        points = np.concatenate(
            [np.pad(c, ((0, 0), (0, dim - c.shape[1]))) for c in point_chunks], axis=0
        )
    out = Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
        regions=_merge_regions(regions),
    )
    if unit is not None:
        set_length_unit(out, unit)
    return out


def _read_legacy(f):
    if "Base" not in f:
        raise ReadError('Expected "Base" in file. Malformed CGNS?')
    if "Zone1" not in f["Base"]:
        raise ReadError('Expected "Zone1" in "Base". Malformed CGNS?')

    x = f["Base"]["Zone1"]["GridCoordinates"]["CoordinateX"][" data"]
    y = f["Base"]["Zone1"]["GridCoordinates"]["CoordinateY"][" data"]
    z = f["Base"]["Zone1"]["GridCoordinates"]["CoordinateZ"][" data"]
    points = np.column_stack([x, y, z])

    idx_min, idx_max = f["Base"]["Zone1"]["GridElements"]["ElementRange"][" data"]
    data = f["Base"]["Zone1"]["GridElements"]["ElementConnectivity"][" data"]
    cells = np.array(data).reshape(idx_max, -1) - 1

    if cells.shape[1] != 4:
        raise ReadError("Can only read tetrahedra.")
    cells = [("tetra", cells)]

    return Mesh(points, cells)
