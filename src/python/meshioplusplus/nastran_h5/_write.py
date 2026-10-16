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

import re

import numpy as np

from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from ..nastran_op2._write import _CARDS as _OP2_CARDS
from ..nastran_op2._write import _Model

TIME_KEY = "meshio:time"
SCHEMA = 20200

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
}
_FAMILY = {
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
_QUIET_POINT = ("nastran:cp", "nastran:cd")
_QUIET_CELL = ("nastran:eid", "nastran:pid")


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


def _element_members(model, mesh, taken):
    """``{(group, card): [(member, column)]}`` of the centre values."""
    groups = {}
    for name in sorted(mesh.cell_data):
        if name in _QUIET_CELL or ":" not in name:
            continue
        group, member = name.split(":", 1)
        if not _GROUP.match(group) or "@" in member:
            continue
        groups.setdefault(group, []).append((member, name))
    out = {}
    for group, members in groups.items():
        for card in sorted(model.cards):
            if card not in _FAMILY:
                continue
            _eid, parts, order = _family_rows(model, card)
            columns = []
            for member, name in members:
                column = _gather(mesh, name, parts, order)
                if column is not None:
                    taken.add(name)
                    columns.append((member, column))
            if columns:
                out[(group, card)] = columns
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
            sorted((k, tuple(m for m, _ in v)) for k, v in self.elemental.items()),
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
            ("ID", "<i8", np.arange(1, n + 1)),
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
        spec = _OP2_CARDS[card]
        linear, perm = spec[2], spec[4]
        width = _CARD_WIDTH[card]
        eid, pid, parts, order = model.element_rows(card)
        nodes = np.zeros((len(eid), width), dtype=np.int64)
        at = 0
        for blk, rows in parts:
            conn = blk["conn"][rows] + 1
            if blk["nodes"] > linear and perm is not None:
                placed = np.zeros_like(conn)
                for i, p in enumerate(perm):
                    placed[:, p] = conn[:, i]
                conn = placed
            nodes[at : at + len(rows), : conn.shape[1]] = conn
            at += len(rows)
        nodes = nodes[order]
        count = len(eid)
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
        _table(
            prop,
            ptype,
            [("PID", "<i8", pids)]
            + _filled(len(pids), _PROPERTY_COLUMNS[ptype])
            + [one(len(pids))],
        )
    if any(model.properties.values()):
        mat = nas.create_group("INPUT/MATERIAL")
        _table(mat, "MAT1", _filled(1, _MATERIAL_COLUMNS) + [one(1)])


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

    def emit(parent_path, name, per_step, columns_of):
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
        members = [(m, d) for m, d, _ in chunks[0][1]]
        merged = []
        for k, (m, d) in enumerate(members):
            merged.append(
                (m, d, np.concatenate([np.asarray(c[1][k][2]) for c in chunks]))
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
            cols = [("ID", "<i8", keep + 1)]
            cols += [(m, "<f8", v[keep]) for m, v in members.items()]
            return cols

        emit(
            "NODAL",
            table,
            [(s, s.nodal) for s in steps if table in s.nodal],
            nodal_columns,
        )
    elemental_keys = sorted({k for s in steps for k in s.elemental})
    for group, card in elemental_keys:

        def element_columns(payload, group=group, card=card):
            columns, eid = payload
            keep = np.nonzero(
                ~np.all(np.isnan(np.stack([c for _, c in columns])), axis=0)
            )[0]
            key = "ID" if group == "ENERGY" else "EID"
            cols = [(key, "<i8", eid[keep])]
            cols += [(m, "<f8", c[keep]) for m, c in columns]
            return cols

        eid = model.element_rows(card)[0]
        emit(
            f"ELEMENTAL/{group}",
            _FAMILY[card],
            [
                (s, (s.elemental[(group, card)], eid))
                for s in steps
                if (group, card) in s.elemental
            ],
            element_columns,
        )


def _assemble(filename, model, steps):
    import h5py

    if is_buffer(filename, "w"):
        _wfail("an HDF5 file cannot be written to a buffer")
    with h5py.File(str(filename), "w") as f:
        f.attrs["SCHEMA"] = np.array([SCHEMA])
        _write_model(f, model)
        _write_results(f, model, steps)


def _build_model(mesh):
    try:
        return _Model(mesh, per_card_ids=True)
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
