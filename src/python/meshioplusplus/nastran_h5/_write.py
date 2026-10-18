"""MSC Nastran HDF5 result files, written.

The twin of the writing half of ``src/cpp/src/formats/nastran_h5.cpp``: both
engines write the same tables. The file is the reader's inverse for the core of
the schema: ``/NASTRAN/INPUT`` (the GRIDs, one table per element card, the
property and material tables the elements name) and ``/NASTRAN/RESULT`` (the
``DOMAINS`` table, the nodal tables and the element tables of the centre
values), each result table with its ``/INDEX`` twin. One mesh is one domain; a
series is one domain per step.

The file is not an MSC Nastran product: it has the ``SCHEMA`` attribute but no
``VERSION`` one, which the reader accepts and an MSC reader may not.
"""

import copy
import re

import numpy as np

from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from ..nastran_op2._write import _CARDS as _OP2_CARDS
from ..nastran_op2._write import PROPERTY_VALUES, _Model, material_ids
from ._schema import SCHEMA as _SCHEMA

TIME_KEY = "meshio:time"
SCHEMA = 20200

# The cards the writer holds a table for. Any other card of ``nastran:card`` (CBUSH, CGAP,
# CMASS1, CMASS2, CONM1: no MSC sample here) is left out, with a warning.
_CARDS = frozenset(
    "CBAR CBEAM CDAMP1 CDAMP2 CELAS1 CELAS2 CHEXA CONM2 CONROD CPENTA CPYRAM CQUAD CQUAD4 "
    "CQUAD8 CROD CSHEAR CTETRA CTRIA3 CTRIA6 CTUBE CVISC PLOTEL".split()
)
# card -> (width of its G column, the family name its result tables take). A
# solid's G holds the quadratic node count, zero-padded for a linear element.
_CARD_WIDTH = {
    "CHEXA": 20,
    "CONM2": 1,
    "CONROD": 2,
    "CPENTA": 15,
    "CPYRAM": 13,
    "CQUAD": 9,
    "CQUAD4": 4,
    "CQUAD8": 8,
    "CROD": 2,
    "CTETRA": 10,
    "CTRIA3": 3,
    "CTRIA6": 6,
    "CTUBE": 2,
    "CSHEAR": 4,
    "CVISC": 2,
    "PLOTEL": 2,
}
_FAMILY = {
    "CBAR": "BARS",
    "CBEAM": "BEAM",
    "CELAS1": "ELAS1",
    "CELAS2": "ELAS2",
    "CSHEAR": "SHEAR",
    "CHEXA": "HEXA",
    "CONROD": "CONROD",
    "CPENTA": "PENTA",
    "CPYRAM": "PYRAM",
    "CQUAD": "QUADR",
    "CQUAD4": "QUAD4",
    "CQUAD8": "QUAD8",
    "CROD": "ROD",
    "CTETRA": "TETRA",
    "CTRIA3": "TRIA3",
    "CTRIA6": "TRIA6",
    "CTUBE": "TUBE",
}
# A card without a PID column (the reader gives such an element property -1).
_NO_PID = ("CONM2", "CONROD")
_PROPERTY_MEMBER = {"PSHELL": "MID1", "PSOLID": "MID", "PROD": "MID", "PTUBE": "MID"}

_GROUP = re.compile(r"^[A-Z][A-Z0-9_]*$")
_QUIET_POINT = ("nastran:cp", "nastran:cd", "nastran:nid")
_QUIET_CELL = (
    "nastran:eid",
    "nastran:pid",
    "nastran:card",
    "nastran:c1",
    "nastran:c2",
    "nastran:x",
    "nastran:g0",
    "nastran:pin",
    "nastran:w",
    "nastran:k",
    "nastran:ge",
    "nastran:s",
    "nastran:b",
)

_KIND_DTYPE = {"i8": "<i8", "f8": "<f8", "S4": "S4", "S8": "S8"}
_schema_cache = {}


def _spec(path):
    """MSC's columns of the table at ``path`` (below ``/NASTRAN``) as ``[(name, kind,
    arity)]``, or ``None`` for a table the fixtures do not hold."""
    if path not in _schema_cache:
        text = _SCHEMA.get(path)
        columns = None
        if text is not None:
            columns = []
            for item in text.split():
                name, kind = item.split(":")
                arity = 1
                if "[" in kind:
                    kind, count = kind[:-1].split("[")
                    arity = int(count)
                columns.append((name, kind, arity))
        _schema_cache[path] = columns
    return _schema_cache[path]


def _fit(values, kind, arity, rows):
    """``values`` as a column of ``kind`` and ``arity`` for ``rows`` rows: the first entries
    of a wider array, zeros after a narrower one."""
    dtype = _KIND_DTYPE[kind]
    array = np.asarray(values)
    if rows == 0:
        return np.zeros((0,) if arity == 1 else (0, arity), dtype=dtype)
    if kind[0] == "S":
        return np.asarray(array, dtype=dtype).reshape(rows)
    if arity == 1:
        return np.asarray(array.reshape(rows, -1)[:, 0], dtype=dtype)
    wide = array.reshape(rows, -1)
    out = np.zeros((rows, arity), dtype=dtype)
    k = min(arity, wide.shape[1])
    out[:, :k] = wide[:, :k]
    return out


def _default(name, kind, arity, rows, grid_arity):
    """The value of a column the mesh gives nothing for: zeros, 'CEN/' for a shell's TERM,
    'GRID' for a solid's CTYPE and its corner count for NODEF."""
    dtype = _KIND_DTYPE[kind]
    if name == "TERM":
        return np.full(rows, b"CEN/", dtype=dtype)
    if name == "CTYPE":
        return np.full(rows, b"GRID", dtype=dtype)
    if name == "NODEF":
        return np.full(rows, max(grid_arity - 1, 0), dtype=dtype)
    if name == "DOMAIN_ID":
        return np.ones(rows, dtype=dtype)
    return np.zeros(rows if arity == 1 else (rows, arity), dtype=dtype)


def _fill(path, provided, rows):
    """The columns of the table at ``path`` as ``[(name, dtype, values)]``: MSC's columns in
    MSC's order, each from ``provided`` (``{name: values}``) or its default, then any
    provided column MSC's table lacks, before ``DOMAIN_ID``. A table the fixtures do not
    hold is its provided columns."""
    spec = _spec(path)
    if spec is None:
        return [(m, _dtype_of(v), v) for m, v in provided.items()]
    grid = next((a for n, _k, a in spec if n == "GRID"), 0)
    out = []
    names = set()
    for name, kind, arity in spec:
        names.add(name)
        if name in provided:
            values = _fit(provided[name], kind, arity, rows)
        else:
            values = _default(name, kind, arity, rows, grid)
        out.append((name, _KIND_DTYPE[kind], values))
    extras = [(m, _dtype_of(v), v) for m, v in provided.items() if m not in names]
    if extras:
        last = out.pop() if out and out[-1][0] == "DOMAIN_ID" else None
        out += extras
        if last is not None:
            out.append(last)
    return out


def _dtype_of(values):
    kind = np.asarray(values).dtype.kind
    return "<i8" if kind in "iu" else "S4" if kind == "S" else "<f8"


def _wfail(message):
    raise WriteError(f"Nastran HDF5: {message}")


def _table(group, name, columns):
    """A compound dataset from ``[(member, dtype, values)]``, values (n,) or (n, k)."""
    fields = []
    for member, dtype, values in columns:
        values = np.asarray(values)
        fields.append(
            (member, dtype, values.shape[1:]) if values.ndim > 1 else (member, dtype)
        )
    rows = len(columns[0][2])
    data = np.zeros(rows, dtype=np.dtype(fields))
    for member, _dtype, values in columns:
        data[member] = values
    group.create_dataset(name, data=data)


def _family_rows(model, card):
    """``(eid, parts, order)`` of a card's rows in the card's element order."""
    eid, _pid, parts, order = model.element_rows(card)
    return eid, parts, order


def _gather(mesh, name, parts, order):
    """A cell-data array of the cells of a card, in the card's row order; ``None``
    when no block of the card has it."""
    arrays = mesh.cell_data.get(name)
    if arrays is None:
        return None
    out = np.full(sum(len(rows) for _, rows in parts), np.nan)
    at = 0
    for blk, rows in parts:
        a = np.asarray(arrays[blk["index"]], dtype=np.float64)
        n = len(blk["eids"])
        if n and (a.size < n or a.size % n):
            _wfail(f"'{name}' must have one value per cell")
        out[at : at + len(rows)] = a.reshape(n, -1)[rows, 0]
        at += len(rows)
    return out[order]


def _nodal_members(mesh, npoints, taken):
    """``{table: {member: (npoints,)}}`` from the point data, members in order."""
    tables = {}
    names = sorted(mesh.point_data)

    def put(table, member, values):
        tables.setdefault(table, {})[member] = values

    def vector(name):
        a = np.asarray(mesh.point_data[name], dtype=np.float64)
        if a.size != npoints * 3:
            _wfail(f"'{name}' must have 3 components per point")
        taken.add(name)
        return a.reshape(npoints, 3)

    for name in names:
        if name in _QUIET_POINT or name.startswith("nastran:"):
            continue
        if name.endswith("_real") and name[:-5] + "_imag" in mesh.point_data:
            base = name[:-5]
            rot = base.endswith("_ROT")
            table = (base[:-4] if rot else base) + "_CPLX"
            prefix = "R" if rot else ""
            re_, im_ = vector(base + "_real"), vector(base + "_imag")
            for k, axis in enumerate("XYZ"):
                put(table, f"{prefix}{axis}R", re_[:, k])
            for k, axis in enumerate("XYZ"):
                put(table, f"{prefix}{axis}I", im_[:, k])
        elif name.endswith("_imag") and name[:-5] + "_real" in mesh.point_data:
            continue  # taken with its real part
    for name in names:
        if name in taken or name in _QUIET_POINT or name.startswith("nastran:"):
            continue
        a = np.asarray(mesh.point_data[name], dtype=np.float64)
        if ":" in name:
            table, member = name.split(":", 1)
            if a.size != npoints:
                continue
            put(table, member, a.reshape(npoints))
            taken.add(name)
        elif name.endswith("_ROT") and a.size == npoints * 3:
            table = name[:-4]
            for k, axis in enumerate(("RX", "RY", "RZ")):
                put(table, axis, vector(name)[:, k])
        elif a.size == npoints * 3:
            table = name
            for k, axis in enumerate("XYZ"):
                put(table, axis, vector(name)[:, k])
        elif a.size == npoints:
            put(name, "VALUE", a.reshape(npoints))
            taken.add(name)
    return tables


def _gather_wide(mesh, name, parts, order, width):
    """A ``@corner`` array of the cells of a card, in the card's row order, as
    ``width`` entries per row: entry 0 is NaN (the centre), entry j + 1 is column j;
    ``None`` when no block of the card has it."""
    arrays = mesh.cell_data.get(name)
    if arrays is None:
        return None
    out = np.full((sum(len(rows) for _, rows in parts), width), np.nan)
    at = 0
    for blk, rows in parts:
        a = np.asarray(arrays[blk["index"]], dtype=np.float64)
        n = len(blk["eids"])
        if n and (a.size < n or a.size % n):
            _wfail(f"'{name}' must have one value per cell and node")
        a = a.reshape(n, -1)
        k = min(a.shape[1], width - 1)
        out[at : at + len(rows), 1 : 1 + k] = a[rows, :k]
        at += len(rows)
    return out[order]


def _corner_members(model, mesh, taken, out):
    """Fold ``<G>:<M>@corner`` into ``out`` as array members: ``{(group, card):
    [("GRID", ints), (member, floats), ...]}`` with ``(rows, width)`` columns whose
    entry 0 is the centre value ``<G>:<M>``."""
    for name in sorted(mesh.cell_data):
        if not name.endswith("@corner") or ":" not in name:
            continue
        group, member = name[: -len("@corner")].split(":", 1)
        if not _GROUP.match(group) or not member:
            continue
        for card in sorted(model.cards):
            if card not in _FAMILY or card in _STATION_CARDS:
                continue
            eid, parts, order = _family_rows(model, card)
            width = 1 + max(blk["conn"].shape[1] for blk, _ in parts)
            wide = _gather_wide(mesh, name, parts, order, width)
            if wide is None:
                continue
            columns = out.setdefault((group, card, ""), ([], eid))[0]
            if not any(m == "GRID" for m, _ in columns):
                grid = np.zeros((wide.shape[0], width), dtype=np.int64)
                at = 0
                for blk, rows in parts:
                    conn = model.node_ids[blk["conn"][rows]]
                    grid[at : at + len(rows), 1 : 1 + conn.shape[1]] = conn
                    at += len(rows)
                columns.append(("GRID", grid[order]))
            for k, (m, column) in enumerate(columns):
                if m == member and column.ndim == 1:
                    wide[:, 0] = column
                    del columns[k]
                    break
            columns.append((member, wide))
            taken.add(name)
    for key in [k for k, v in out.items() if not v[0]]:
        del out[key]


_STATION_CARDS = ("CBAR", "CBEAM")
# the shells a composite's plies are results of, and the family their ply tables take
_PLY_FAMILY = {
    "CQUAD4": "QUAD4",
    "CTRIA3": "TRIA3",
    "CQUAD": "QUADR",
    "CQUAD8": "QUAD8",
    "CTRIA6": "TRIA6",
}


def _gather_columns(mesh, name, parts, order):
    """A ``(cells, columns)`` array of the cells of a card, in the card's row order, NaN
    where a block lacks columns; ``None`` when no block of the card has it."""
    arrays = mesh.cell_data.get(name)
    if arrays is None:
        return None
    width = 0
    shaped = []
    for blk, rows in parts:
        a = np.asarray(arrays[blk["index"]], dtype=np.float64)
        n = len(blk["eids"])
        if n and (a.size < n or a.size % n):
            _wfail(f"'{name}' must have one value per cell and column")
        a = a.reshape(n, -1)
        width = max(width, a.shape[1])
        shaped.append((a, rows))
    out = np.full((sum(len(rows) for _, rows in parts), max(width, 1)), np.nan)
    at = 0
    for a, rows in shaped:
        out[at : at + len(rows), : a.shape[1]] = a[rows]
        at += len(rows)
    return out[order]


def _split_names(mesh, suffix):
    """``{group: {member: name}}`` of the arrays ``<G>:<M><suffix>`` ('' for the plain ones,
    which have no ``@``)."""
    found = {}
    for name in sorted(mesh.cell_data):
        if name in _QUIET_CELL or ":" not in name:
            continue
        base = name[: -len(suffix)] if suffix and name.endswith(suffix) else name
        if suffix and base == name:
            continue
        group, member = base.split(":", 1)
        if not _GROUP.match(group) or not member or "@" in member:
            continue
        found.setdefault(group, {})[member] = name
    return found


def _ply_members(model, mesh, taken, out):
    """``<G>:<M>@ply`` of the shells as rows of a ``<FAMILY>_COMP`` table: one per ply with
    any value, ``PLY`` (from 1) and the members."""
    for group, members in _split_names(mesh, "@ply").items():
        for card in sorted(model.cards):
            if card not in _PLY_FAMILY:
                continue
            eid, parts, order = _family_rows(model, card)
            arrays = {}
            for member, name in members.items():
                wide = _gather_columns(mesh, name, parts, order)
                if wide is not None:
                    arrays[member] = wide
                    taken.add(name)
            if not arrays:
                continue
            width = max(a.shape[1] for a in arrays.values())
            shown = np.zeros((len(eid), width), dtype=bool)
            for a in arrays.values():
                shown[:, : a.shape[1]] |= ~np.isnan(a)
            rr, cc = np.nonzero(shown)
            columns = [("PLY", cc + 1)]
            for member, a in arrays.items():
                column = np.full(len(rr), np.nan)
                ok = cc < a.shape[1]
                column[ok] = a[rr[ok], cc[ok]]
                columns.append((member, column))
            out[(group, card, "ply")] = (columns, eid[rr])


def _station_members(model, mesh, taken, out):
    """The bars' and beams' values: a ``BEAM`` row per beam with 11-slot arrays (``GRID``
    names its two ends, ``SD`` the station along it) and a ``BARS`` row per station of a
    bar, from ``<G>:<M>@station`` (entry 0 is end A) and the plain ``<G>:<M>`` (end A).
    """
    plain = _split_names(mesh, "")
    stations = _split_names(mesh, "@station")
    for group in sorted(set(plain) | set(stations)):
        for card in _STATION_CARDS:
            if card not in model.cards:
                continue
            eid, parts, order = _family_rows(model, card)
            arrays = {}
            for member in sorted(
                set(plain.get(group, {})) | set(stations.get(group, {}))
            ):
                wide = None
                if member in stations.get(group, {}):
                    name = stations[group][member]
                    wide = _gather_columns(mesh, name, parts, order)
                    taken.add(name)
                if member in plain.get(group, {}):
                    name = plain[group][member]
                    single = _gather(mesh, name, parts, order)
                    if single is not None:
                        taken.add(name)
                        if wide is None:
                            wide = single.reshape(-1, 1)
                        else:
                            wide = wide.copy()
                            wide[:, 0] = np.where(
                                np.isnan(wide[:, 0]), single, wide[:, 0]
                            )
                if wide is not None:
                    arrays[member] = wide
            if not arrays:
                continue
            if card == "CBEAM":
                out[(group, card, "")] = (
                    _beam_columns(model, parts, order, arrays),
                    eid,
                )
            else:
                out[(group, card, "")] = _bar_columns(eid, arrays)


def _beam_columns(model, parts, order, arrays):
    """The columns of a ``BEAM`` table: ``GRID`` and the members as 11-slot arrays."""
    n = len(order)
    members = []
    shown = np.zeros((n, 11), dtype=bool)
    for member, a in arrays.items():
        full = np.full((n, 11), np.nan)
        k = min(11, a.shape[1])
        full[:, :k] = a[:, :k]
        members.append((member, full))
        shown |= ~np.isnan(full)
    grid = np.zeros((n, 11), dtype=np.int64)
    at = 0
    for blk, rows in parts:
        conn = model.node_ids[blk["conn"][rows]]
        grid[at : at + len(rows), 0] = conn[:, 0]
        grid[at : at + len(rows), 10] = conn[:, 1]
        at += len(rows)
    grid = grid[order]
    columns = [("GRID", grid)]
    if "SD" not in arrays:  # the station along the beam: 0 to 1 over the 11 slots
        columns.append(("SD", np.where(shown, np.arange(11) / 10.0, 0.0)))
    return columns + members


def _bar_columns(eid, arrays):
    """The rows of a ``BARS`` table: one per station of each bar with any value, and ``SD``
    (the station along the bar, 0 to 1) unless the mesh gives it."""
    width = max(a.shape[1] for a in arrays.values())
    shown = np.zeros((len(eid), width), dtype=bool)
    for a in arrays.values():
        shown[:, : a.shape[1]] |= ~np.isnan(a)
    rr, cc = np.nonzero(shown)
    columns = []
    if "SD" not in arrays:
        columns.append(("SD", cc / max(width - 1, 1)))
    for member, a in arrays.items():
        column = np.full(len(rr), np.nan)
        ok = cc < a.shape[1]
        column[ok] = a[rr[ok], cc[ok]]
        columns.append((member, column))
    return columns, eid[rr]


def _element_members(model, mesh, taken):
    """``{(group, card, variant): (columns, eid)}`` of the centre values, of the corner
    values as array members (variant ''), of the composite plies ('ply') and of the bars'
    and beams' stations: ``columns`` is ``[(member, values)]`` and ``eid`` the element id
    of each row."""
    out = {}
    for group, members in _split_names(mesh, "").items():
        for card in sorted(model.cards):
            if card not in _FAMILY or card in _STATION_CARDS:
                continue
            eid, parts, order = _family_rows(model, card)
            columns = []
            for member, name in members.items():
                column = _gather(mesh, name, parts, order)
                if column is not None:
                    taken.add(name)
                    columns.append((member, column))
            if columns:
                out[(group, card, "")] = (columns, eid)
    _corner_members(model, mesh, taken, out)
    _ply_members(model, mesh, taken, out)
    _station_members(model, mesh, taken, out)
    return out


def _scalar(mesh, name, default):
    array = mesh.field_data.get(name)
    if array is None:
        return default
    values = np.asarray(array).ravel()
    return values[0] if values.size else default


class _Step:
    """One domain: its row of ``DOMAINS`` and its result tables."""

    def __init__(self, model, mesh, time, index, taken):
        self.row = (
            index + 1,
            int(_scalar(mesh, "nastran:subcase", 1)),
            int(_scalar(mesh, "nastran:step", index)),
            int(_scalar(mesh, "nastran:analysis", 1 if time == 0.0 else 6)),
            float(time),
            float(_scalar(mesh, "nastran:eigi", 0.0)),
            int(_scalar(mesh, "nastran:mode", 1)),
        )
        self.nodal = _nodal_members(mesh, model.npoints, taken)
        self.elemental = _element_members(model, mesh, taken)
        self.signature = (
            sorted((t, tuple(m)) for t, m in self.nodal.items()),
            sorted((k, tuple(m for m, _ in v[0])) for k, v in self.elemental.items()),
        )


def _time_of(mesh):
    array = mesh.field_data.get(TIME_KEY)
    if array is None:
        return 0.0
    values = np.asarray(array, dtype=np.float64).ravel()
    return float(values[0]) if values.size else 0.0


def _warn_dropped(model, mesh, taken):
    names = [n for n in mesh.point_data if n not in taken and n not in _QUIET_POINT]
    names += [n for n in mesh.cell_data if n not in taken and n not in _QUIET_CELL]
    names += [
        n
        for n in mesh.field_data
        if n != TIME_KEY and not n.startswith("nastran:") and n not in taken
    ]
    if names:
        warn(
            "Nastran HDF5: arrays without a place in the file are not written: "
            + ", ".join(sorted(set(names)))
        )
    if model.skipped:
        warn(
            "Nastran HDF5: cell blocks with no element card are not written: "
            + ", ".join(sorted(set(model.skipped)))
        )
    if model.unwritten:
        warn(
            "Nastran HDF5: element cards with no table here are left out: "
            + ", ".join(sorted(model.unwritten))
        )


# Columns the MSC tables carry beyond the ones the reader uses, with the value
# this writer gives them: (member, dtype, shape or (), value). A table the schema
# holds in full is what a reader of the schema (pyNastran's, MSC's) expects.
_SHELL_T4 = ("T", "<f8", (4,), 0.0)
_SHELL_T3 = ("T", "<f8", (3,), 0.0)
_ELEMENT_EXTRA = {
    "CQUAD": [
        ("THETA", "<f8", (), 0.0),
        ("ZOFFS", "<f8", (), 0.0),
        ("TFLAG", "<i8", (), 0),
    ]
    + [_SHELL_T4, ("MCID", "<i8", (), -1)],
    "CQUAD4": [
        ("THETA", "<f8", (), 0.0),
        ("ZOFFS", "<f8", (), 0.0),
        ("TFLAG", "<i8", (), 0),
    ]
    + [_SHELL_T4, ("MCID", "<i8", (), -1)],
    "CQUAD8": [_SHELL_T4, ("THETA", "<f8", (), 0.0), ("ZOFFS", "<f8", (), 0.0)]
    + [("TFLAG", "<i8", (), 0), ("MCID", "<i8", (), -1)],
    "CTRIA3": [
        ("THETA", "<f8", (), 0.0),
        ("ZOFFS", "<f8", (), 0.0),
        ("TFLAG", "<i8", (), 0),
    ]
    + [_SHELL_T3, ("MCID", "<i8", (), -1)],
    "CTRIA6": [("THETA", "<f8", (), 0.0), ("ZOFFS", "<f8", (), 0.0), _SHELL_T3]
    + [("TFLAG", "<i8", (), 0), ("MCID", "<i8", (), -1)],
}
_PROPERTY_COLUMNS = {
    "PROD": [("MID", "<i8", (), 1), ("A", "<f8", (), 1.0), ("J", "<f8", (), 0.0)]
    + [("C", "<f8", (), 0.0), ("NSM", "<f8", (), 0.0)],
    "PSHELL": [("MID1", "<i8", (), 1), ("T", "<f8", (), 1.0), ("MID2", "<i8", (), -1)]
    + [("BK", "<f8", (), 1.0), ("MID3", "<i8", (), -1), ("TS", "<f8", (), 5.0 / 6.0)]
    + [("NSM", "<f8", (), 0.0), ("Z1", "<f8", (), -0.5), ("Z2", "<f8", (), 0.5)]
    + [("MID4", "<i8", (), -1)],
    "PSOLID": [("MID", "<i8", (), 1), ("CORDM", "<i8", (), 0), ("IN", "<i8", (), 0)]
    + [("STRESS", "<i8", (), 0), ("ISOP", "<i8", (), 0), ("FCTN", "S4", (), b"SMEC")],
    "PTUBE": [("MID", "<i8", (), 1), ("OD", "<f8", (), 1.0), ("T", "<f8", (), 0.0)]
    + [("NSM", "<f8", (), 0.0)],
}
_MATERIAL_COLUMNS = (
    [("MID", "<i8", (), 1), ("E", "<f8", (), 1.0), ("G", "<f8", (), 0.0)]
    + [("NU", "<f8", (), 0.3), ("RHO", "<f8", (), 0.0), ("A", "<f8", (), 0.0)]
    + [("TREF", "<f8", (), 0.0), ("GE", "<f8", (), 0.0), ("ST", "<f8", (), 0.0)]
    + [("SC", "<f8", (), 0.0), ("SS", "<f8", (), 0.0), ("MCSID", "<i8", (), 0)]
)
_DOMAIN_EXTRA = ("DESIGN_CYCLE", "RANDOM", "SE", "AFPM", "TRMC", "INSTANCE", "MODULE")


def _filled(rows, spec):
    """The columns of ``spec`` ((member, dtype, shape, value) entries) for ``rows`` rows."""
    out = []
    for member, dtype, shape, value in spec:
        out.append((member, dtype, np.full((rows,) + tuple(shape), value, dtype=dtype)))
    return out


def _card_rows(model, card):
    """A card's rows in element order: ``(eid, pid, nodes, comp, bar)`` with the node ids
    in Nastran's order, zero-padded to the card's width, the spring components and the
    21 columns of card fields."""
    spec = _OP2_CARDS[card]
    linear, perm = spec[2], spec[4]
    width = _CARD_WIDTH.get(card, 2)
    eid, pid, parts, order = model.element_rows(card)
    nodes = np.zeros((len(eid), width), dtype=np.int64)
    at = 0
    for blk, rows in parts:
        conn = model.node_ids[blk["conn"][rows]]
        if blk["nodes"] > linear and perm is not None:
            placed = np.zeros_like(conn)
            for i, p in enumerate(perm):
                placed[:, p] = conn[:, i]
            conn = placed
        k = min(conn.shape[1], width)
        nodes[at : at + len(rows), :k] = conn[:, :k]
        at += len(rows)
    comp = np.concatenate([blk["comps"][r] for blk, r in parts])[order]
    bar = np.concatenate([blk["bar"][r] for blk, r in parts])[order]
    return eid, pid, nodes[order], comp, bar


def _bar_card(card, eid, pid, nodes, comp, bar):
    """The provided columns of a card the writer fills from the schema (see ``_fill``)."""
    g0 = bar[:, 3].astype(np.int64)
    pins = bar[:, 4:6].astype(np.int64)
    if card == "CBAR":
        out = {"EID": eid, "PID": pid, "GA": nodes[:, 0], "GB": nodes[:, 1]}
        out["FLAG"] = np.where(g0 != 0, 2, 1)
        for k in range(3):
            out[f"X{k + 1}"] = bar[:, k]
        out.update({"GO": g0, "PA": pins[:, 0], "PB": pins[:, 1]})
        for k, name in enumerate(("W1A", "W2A", "W3A", "W1B", "W2B", "W3B")):
            out[name] = bar[:, 6 + k]
        return out
    if card == "CBEAM":
        return {
            "EID": eid,
            "PID": pid,
            "GA": nodes[:, 0],
            "GB": nodes[:, 1],
            "X": bar[:, 0:3],
            "G0": g0,
            "F": np.where(g0 != 0, 2, 1),
            "PA": pins[:, 0],
            "PB": pins[:, 1],
            "WA": bar[:, 6:9],
            "WB": bar[:, 9:12],
        }
    ends = {"G1": nodes[:, 0], "G2": nodes[:, 1], "C1": comp[:, 0], "C2": comp[:, 1]}
    if card in ("CELAS1", "CDAMP1"):
        return {"EID": eid, "PID": pid, **ends}
    if card == "CELAS2":
        return {"EID": eid, "K": bar[:, 0], **ends, "GE": bar[:, 1], "S": bar[:, 2]}
    if card == "CDAMP2":
        return {"EID": eid, "B": bar[:, 0], **ends}
    if card in ("CSHEAR", "CVISC"):
        return {"EID": eid, "PID": pid, "G": nodes}
    return {"EID": eid, "G": nodes}  # PLOTEL


_SCHEMA_CARDS = "CBAR CBEAM CDAMP1 CDAMP2 CELAS1 CELAS2 CSHEAR CVISC PLOTEL".split()
# the properties filled from the schema, with the members the field data gives
_SCHEMA_PROPERTIES = ("PBAR", "PELAS", "PDAMP", "PVISC", "PSHEAR")
# the OP2 name of a MAT2 member called otherwise in MSC's table
_MAT2_COLUMN = {"A3": "A12"}


def _write_model(f, model):
    nas = f.create_group("NASTRAN")
    n = model.npoints
    one = lambda rows: ("DOMAIN_ID", "<i8", np.ones(rows, dtype="<i8"))  # noqa: E731
    _table(
        nas.create_group("INPUT"),
        "DOMAINS",
        [("ID", "<i8", [1])]
        + [(m, "<i8", [0]) for m in ("SE", "AFPM", "TRMC")]
        + [("MODULE", "<i8", [1])],
    )
    node = nas.create_group("INPUT/NODE")
    zeros = np.zeros(n, dtype="<i8")
    _table(
        node,
        "GRID",
        [
            ("ID", "<i8", model.node_ids),
            ("CP", "<i8", zeros),
            ("X", "<f8", model.points),
            ("CD", "<i8", zeros),
            ("PS", "<i8", zeros),
            ("SEID", "<i8", zeros),
            one(n),
        ],
    )
    element = nas.create_group("INPUT/ELEMENT")
    for card in sorted(model.cards):
        eid, pid, nodes, comp, bar = _card_rows(model, card)
        count = len(eid)
        if card in _SCHEMA_CARDS:
            columns = _fill(
                f"INPUT/ELEMENT/{card}",
                _bar_card(card, eid, pid, nodes, comp, bar),
                count,
            )
            _table(element, card, columns)
            continue
        if card == "CONROD":
            columns = [("EID", "<i8", eid), ("G1", "<i8", nodes[:, 0])]
            columns += [("G2", "<i8", nodes[:, 1])]
            columns += _filled(
                count,
                [("MID", "<i8", (), 1), ("A", "<f8", (), 1.0), ("J", "<f8", (), 0.0)]
                + [("C", "<f8", (), 0.0), ("NSM", "<f8", (), 0.0)],
            )
        elif card == "CONM2":
            columns = [("EID", "<i8", eid), ("G", "<i8", nodes[:, 0])]
            columns += _filled(
                count,
                [("CID", "<i8", (), 0), ("M", "<f8", (), 1.0), ("X1", "<f8", (), 0.0)]
                + [("X2", "<f8", (), 0.0), ("X3", "<f8", (), 0.0)]
                + [("I1", "<f8", (), 0.0), ("I2", "<f8", (2,), 0.0)]
                + [("I3", "<f8", (3,), 0.0)],
            )
        else:
            columns = [("EID", "<i8", eid), ("PID", "<i8", pid), ("G", "<i8", nodes)]
            columns += _filled(count, _ELEMENT_EXTRA.get(card, []))
        columns.append(one(count))
        _table(element, card, columns)
    prop = nas.create_group("INPUT/PROPERTY")
    for ptype in ("PROD", "PSHELL", "PSOLID", "PTUBE"):
        pids = sorted(model.properties[ptype])
        if not pids:
            continue
        columns = _filled(len(pids), _PROPERTY_COLUMNS[ptype])
        for k, pid in enumerate(pids):
            for member, value in model.values.get((ptype, pid), {}).items():
                for c in columns:
                    if c[0] == member:
                        c[2][k] = value
        _table(prop, ptype, [("PID", "<i8", pids)] + columns + [one(len(pids))])
    for ptype in _SCHEMA_PROPERTIES:
        pids = sorted(model.properties[ptype])
        if not pids:
            continue
        provided = {"PID": pids}
        for member, (_word, default) in PROPERTY_VALUES[ptype].items():
            given = [
                model.values.get((ptype, pid), {}).get(member, default) for pid in pids
            ]
            provided[member] = given
        _table(prop, ptype, _fill(f"INPUT/PROPERTY/{ptype}", provided, len(pids)))
    pids = sorted(model.properties["PBEAM"])
    if pids:
        provided = {"PID": pids, **{m: [] for m in ("MID", "NSEGS", "CCF")}}
        names = "SO XXB A I1 I2 I12 J NSM C1 C2 D1 D2 E1 E2 F1 F2".split()
        for name in names:
            provided[name] = []
        for pid in pids:
            mid, nseg, ccf, slots = model.pbeam_slots(pid)
            provided["MID"].append(mid)
            provided["NSEGS"].append(nseg)
            provided["CCF"].append(ccf)
            for column, name in enumerate(names):
                provided[name].append(slots[:, column])
        _table(prop, "PBEAM", _fill("INPUT/PROPERTY/PBEAM", provided, len(pids)))
    pids = sorted(model.properties["PCOMP"])
    if pids:
        group = prop.create_group("PCOMP")
        head = {"PID": pids, "NPLIES": [], "PLY_POS": [], "PLY_LEN": []}
        for member in ("Z0", "NSM", "SB", "FT", "TREF", "GE"):
            head[member] = []
        plies, position = [], 0
        for pid in pids:
            given = model.values[("PCOMP", pid)]
            rows = np.asarray(given["PLY"], dtype=np.float64)
            head["NPLIES"].append(len(rows))
            head["PLY_POS"].append(position)
            head["PLY_LEN"].append(len(rows))
            for member in ("Z0", "NSM", "SB", "FT", "TREF", "GE"):
                head[member].append(
                    given.get(member, PROPERTY_VALUES["PCOMP"][member][1])
                )
            plies.append(rows)
            position += len(rows)
        _table(
            group, "IDENTITY", _fill("INPUT/PROPERTY/PCOMP/IDENTITY", head, len(pids))
        )
        rows = np.concatenate(plies)
        ply = {
            "MID": rows[:, 0],
            "T": rows[:, 1],
            "THETA": rows[:, 2],
            "SOUT": rows[:, 3],
        }
        _table(group, "PLY", _fill("INPUT/PROPERTY/PCOMP/PLY", ply, len(rows)))
    _write_materials(nas, model)


def _write_materials(nas, model):
    """The MAT1 table of every material id the properties name, the MAT2 and MAT4 tables of
    those whose field data gives a MAT2 or MAT4 (any other card's values are not written:
    the id is a MAT1)."""
    from ..nastran_op2._write import material_cards

    mids = material_ids(model.properties, model.values, model.beams, model.lines)
    if not mids:
        return
    cards = material_cards(model.values)
    mat = nas.create_group("INPUT/MATERIAL")
    by_card = {"MAT1": [], "MAT2": [], "MAT4": []}
    for mid in mids:
        card = cards.get(mid, "MAT1")
        by_card[card if card in by_card else "MAT1"].append(mid)
    for card, ids in by_card.items():
        if not ids:
            continue
        if card == "MAT1":
            columns = _filled(len(ids), _MATERIAL_COLUMNS)
            for k, mid in enumerate(ids):
                columns[0][2][k] = mid
                for member, value in model.values.get(("MAT1", mid), {}).items():
                    for c in columns:
                        if c[0] == member:
                            c[2][k] = value
            _table(
                mat, "MAT1", columns + [("DOMAIN_ID", "<i8", np.ones(len(ids), "<i8"))]
            )
            continue
        provided = {"MID": ids}
        for member, (_word, default) in PROPERTY_VALUES[card].items():
            column = _MAT2_COLUMN.get(member, member) if card == "MAT2" else member
            provided[column] = [
                model.values.get((card, mid), {}).get(member, default) for mid in ids
            ]
        _table(mat, card, _fill(f"INPUT/MATERIAL/{card}", provided, len(ids)))


def _write_results(f, model, steps):
    if not any(s.nodal or s.elemental for s in steps):
        return
    res = f["NASTRAN"].create_group("RESULT")
    rows = [s.row for s in steps]
    zeros = np.zeros(len(rows), dtype="<i8")
    _table(
        res,
        "DOMAINS",
        [
            ("ID", "<i8", [r[0] for r in rows]),
            ("SUBCASE", "<i8", [r[1] for r in rows]),
            ("STEP", "<i8", [r[2] for r in rows]),
            ("ANALYSIS", "<i8", [r[3] for r in rows]),
            ("TIME_FREQ_EIGR", "<f8", [r[4] for r in rows]),
            ("EIGI", "<f8", [r[5] for r in rows]),
            ("MODE", "<i8", [r[6] for r in rows]),
        ]
        + [(m, "<i8", zeros) for m in _DOMAIN_EXTRA]
        + [("SUBSTEP", "<i8", zeros), ("IMPFID", "<i8", zeros)],
    )
    index = f.create_group("INDEX/NASTRAN/RESULT")

    def emit(parent_path, name, per_step, columns_of, path=None):
        """One table over every step that has rows for it, and its INDEX twin."""
        chunks, idx, position = [], [], 0
        for step, payload in per_step:
            cols = columns_of(payload)
            count = len(cols[0][2])
            if not count:
                continue
            chunks.append((step, cols))
            idx.append((step.row[0], position, count))
            position += count
        if not chunks:
            return
        # DOMAIN_ID is the table's last column: every row's step
        names = [m for m, _d, _ in chunks[0][1] if m != "DOMAIN_ID"]
        merged = []
        for k, m in enumerate(names):
            dtype = chunks[0][1][k][1]
            merged.append(
                (m, dtype, np.concatenate([np.asarray(c[1][k][2]) for c in chunks]))
            )
        merged.append(
            (
                "DOMAIN_ID",
                "<i8",
                np.concatenate(
                    [np.full(len(c[1][0][2]), c[0].row[0], dtype="<i8") for c in chunks]
                ),
            )
        )
        group = res.require_group(parent_path) if parent_path else res
        _table(group, name, merged)
        igroup = index.require_group(parent_path) if parent_path else index
        _table(
            igroup,
            name,
            [
                ("DOMAIN_ID", "<i8", [i[0] for i in idx]),
                ("POSITION", "<i8", [i[1] for i in idx]),
                ("LENGTH", "<i8", [i[2] for i in idx]),
            ],
        )

    nodal_names = sorted({t for s in steps for t in s.nodal})
    for table in nodal_names:

        def nodal_columns(step_members, table=table):
            members = step_members[table]
            keep = np.nonzero(
                ~np.all(np.isnan(np.stack(list(members.values()))), axis=0)
            )[0]
            provided = {"ID": model.node_ids[keep]}
            provided.update({m: v[keep] for m, v in members.items()})
            return _fill(f"RESULT/NODAL/{table}", provided, len(keep))

        emit(
            "NODAL",
            table,
            [(s, s.nodal) for s in steps if table in s.nodal],
            nodal_columns,
            f"RESULT/NODAL/{table}",
        )
    elemental_keys = sorted({k for s in steps for k in s.elemental})
    for group, card, variant in elemental_keys:
        name = _table_name(group, card, variant)

        def element_columns(payload, group=group, name=name):
            columns, eid = payload
            empty = np.ones(len(eid), dtype=bool)
            for m, c in columns:
                if m not in ("GRID", "PLY"):
                    empty &= np.isnan(c).reshape(len(c), -1).all(axis=1)
            keep = np.nonzero(~empty)[0]
            provided = {"ID" if group == "ENERGY" else "EID": eid[keep]}
            provided.update({m: c[keep] for m, c in columns})
            return _fill(f"RESULT/ELEMENTAL/{group}/{name}", provided, len(keep))

        emit(
            f"ELEMENTAL/{group}",
            name,
            [
                (s, s.elemental[(group, card, variant)])
                for s in steps
                if (group, card, variant) in s.elemental
            ],
            element_columns,
            f"RESULT/ELEMENTAL/{group}/{name}",
        )


def _table_name(group, card, variant):
    """The name of a result table: the card's family (MSC calls a CQUAD4's stresses and
    strains ``QUAD_CN`` and its forces ``QUAD4_CN``), ``_COMP`` for plies."""
    if variant == "ply":
        return _PLY_FAMILY[card] + "_COMP"
    if card == "CQUAD4":
        return "QUAD4_CN" if group == "ELEMENT_FORCE" else "QUAD_CN"
    return _FAMILY[card]


def _assemble(filename, model, steps):
    import h5py

    if is_buffer(filename, "w"):
        _wfail("an HDF5 file cannot be written to a buffer")
    with h5py.File(str(filename), "w") as f:
        f.attrs["SCHEMA"] = np.array([SCHEMA])
        _write_model(f, model)
        _write_results(f, model, steps)


_LINE_PROPERTY = re.compile(r"^nastran:(PBARL|PBEAML|PCOMPG)_[0-9]+:")


def _build_model(mesh):
    # PBARL, PBEAML and PCOMPG values have no table here: their bars, beams and shells
    # are written on the plain PBAR, PBEAM and PSHELL cards, with a warning.
    lines = sorted(
        {m.group(1) for n in mesh.field_data if (m := _LINE_PROPERTY.match(n))}
    )
    if lines:
        mesh = copy.copy(mesh)
        mesh.field_data = {
            k: v for k, v in mesh.field_data.items() if not _LINE_PROPERTY.match(k)
        }
        warn(
            "Nastran HDF5: the " + ", ".join(lines) + " values have no table here; the "
            "elements are written on the plain property cards"
        )
    try:
        return _Model(mesh, per_card_ids=True, allowed=_CARDS)
    except WriteError as exc:
        raise WriteError(str(exc).replace("Nastran OP2:", "Nastran HDF5:")) from None


def write(filename, mesh):
    """Write one mesh as a model plus one result domain."""
    model = _build_model(mesh)
    taken = set()
    step = _Step(model, mesh, _time_of(mesh), 0, taken)
    _warn_dropped(model, mesh, taken)
    _assemble(filename, model, [step])


class SeriesWriter:
    """Write a time series into one results file: the model once, a domain per
    step. A step's time is the argument, not its ``meshio:time``; its subcase,
    analysis and mode are its ``nastran:*`` field data. Every step must have the
    first step's points, cell blocks and result tables. The file is assembled
    when the writer closes. Use as a context manager, or call :meth:`close`."""

    def __init__(self, filename):
        self._filename = filename
        self._model = None
        self._steps = []
        self._layout = None

    def write(self, time, mesh):
        index = len(self._steps)
        taken = set()
        if self._model is None:
            self._model = _build_model(mesh)
            self._signature = _signature(mesh)
        elif self._signature != _signature(mesh):
            _wfail(
                f"series: step {index}'s points or cells differ from the first "
                "step's; a results file holds one model -- write one file per "
                "step with '{step}'"
            )
        step = _Step(self._model, mesh, float(time), index, taken)
        if self._layout is None:
            self._layout = step.signature
            _warn_dropped(self._model, mesh, taken)
        elif step.signature != self._layout:
            _wfail(
                f"series: step {index}'s result arrays differ from the first "
                "step's; a result table holds the same members in every domain"
            )
        self._steps.append(step)

    def close(self):
        if self._filename is None:
            return
        filename, self._filename = self._filename, None
        if self._steps:
            _assemble(filename, self._model, self._steps)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if exc[0] is None:
            self.close()
        else:
            self._filename = None
        return False


def _signature(mesh):
    return (len(mesh.points), [(c.type, np.asarray(c.data).shape) for c in mesh.cells])
