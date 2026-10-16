"""Write Ansys MAPDL binary results (``.rst``, ``.rth``).

The Python twin of ``write_ansys_rst`` in ``src/cpp/src/formats/ansys_rst.cpp``:
both write the same bytes, in the ``fdresu.inc`` layout ``_ansys_rst.py`` reads
and pymapdl-reader opens: Fortran-style records (``i32 length`` in words,
``i32 flags``, the data, a trailer word repeating the length) addressed by word
offsets, little-endian, 64-bit pointers, uncompressed.

The file holds the standard header, the result header, the node and element
equivalence tables, the geometry (element types, nodes, elements, components),
and per result set the solution header, the nodal DOF solution, the reaction
forces and the element nodal stresses and strains. The model is the ``.cdb``
writer's (``ansysInp.element_types``): a cell keeps its ``ansys:element``
routine and ``ansys:type`` slot where they fit, nodes are numbered ``1..n`` in
point order and elements ``1..m`` in cell order, point and cell regions are
components.

Not written: element nodal forces, the raw element records (``EMS``, ``ENG``...),
element rotation angles (the results are in the global axes), nodal coordinate
system angles, coordinate systems, materials, and the cyclic-symmetry records.
"""

from __future__ import annotations

import numpy as np

from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from .._regions import canonical_regions
from ..ansysInp._ansysInp import _cell_int, _layout, _num_corners, element_types
from ._ansys_rst import (
    _DOF_LABELS,
    _NO_RESULT_CELLS,
    _RAW_RECORDS,
    _REACTION_VECTORS,
    _SCALARS,
    _TENSORS,
    _TOP,
    _UNDEFINED,
    _VECTORS,
)

__all__ = ["write", "SeriesWriter"]

_VERSION = "15.0"
_NAME_WORDS = 8
# Entries of an element's pointer table (``fdresu.inc``'s ``ptrENF`` ...).
_ENS, _EEL, _EPL, _ECR, _ETH = 2, 5, 6, 7, 8
_TABLE_ENTRIES = {"S": _ENS, "EPEL": _EEL, "EPPL": _EPL, "EPCR": _ECR, "EPTH": _ETH}
_ITEMS = {"S": 6, "EPEL": 7, "EPPL": 7, "EPCR": 7, "EPTH": 8}
# AvailData bits (``resucm.inc``): item k is bit k.
_AVAIL = {"S": 3, "EPEL": 6, "EPPL": 7, "EPCR": 8, "EPTH": 9, "NSL": 27, "RF": 30}
_ENF = 1
_AVAIL["ENF"] = _ENF + 1
_AVAIL.update({name: entry + 1 for name, entry in _RAW_RECORDS})
_CODE_OF = {name: code for code, name in _SCALARS.items()}
# RF_AX ... RF_TEMP: the DOFs that are not in an RF or RMOM triple.
_REACTION_SCALAR = {
    "RF_" + label: code for code, label in _DOF_LABELS.items() if code > 6
}


def _i32(value):
    value = int(value) & 0xFFFFFFFF
    return value - (1 << 32) if value >= 1 << 31 else value


def _text_words(text, count):
    """A name as ``count`` 4-character words, each read big-endian (the reader's
    ``_name``)."""
    raw = text.encode("latin-1", "replace")[: 4 * count].ljust(4 * count, b" ")
    return [_i32(int.from_bytes(raw[4 * k : 4 * k + 4], "big")) for k in range(count)]


class _Words:
    """The file as int32 words."""

    def __init__(self):
        self.w = []

    def record(self, data, flags):
        ptr = len(self.w)
        self.w += [len(data), _i32(flags << 24), *data, len(data)]
        return ptr

    def ints(self, values, size=None, flags=0x80):
        values = [_i32(v) for v in values]
        if size is not None:
            values += [0] * (size - len(values))
        return self.record(values, flags)

    def doubles(self, values):
        a = np.ascontiguousarray(np.asarray(values, dtype="<f8").reshape(-1))
        return self.record(a.view("<i4").tolist(), 0x00)

    def patch(self, ptr, k, value):
        self.w[ptr + 2 + k] = _i32(value)

    def patch64(self, ptr, lo, hi, value):
        self.patch(ptr, lo, value & 0xFFFFFFFF)
        self.patch(ptr, hi, value >> 32)

    def tobytes(self):
        return np.asarray(self.w, dtype="<i4").tobytes()


class _Model:
    """The mesh part of the file: nodes, elements, element types, components,
    and what each element's results are placed on."""

    def __init__(self, mesh):
        points = np.asarray(mesh.points, dtype=np.float64)
        n = len(points)
        if n == 0:
            raise WriteError("Ansys .rst: the mesh has no points")
        self.points = np.zeros((n, 3))
        points = points.reshape(n, -1)
        self.points[:, : min(3, points.shape[1])] = points[:, :3]
        self.nnod = n
        for block in mesh.cells:
            data = np.asarray(block.data)
            if len(data) and (data.min() < 0 or data.max() >= n):
                raise WriteError(
                    f"Ansys .rst: a '{block.type}' cell refers to a point the mesh "
                    "does not have"
                )
        self.cell_etype, self.slot_routine = element_types(mesh, "Ansys .rst writer")
        self.has_results = [b.type not in _NO_RESULT_CELLS for b in mesh.cells]
        # Per element: (block, row, routine, slot, nodes in the element's order,
        # the element node each cell node became).
        self.elements = []
        for b, block in enumerate(mesh.cells):
            for r, row in enumerate(np.asarray(block.data).tolist()):
                routine, slot = self.cell_etype[b][r]
                layout = _layout(block.type, routine)
                first = {}
                for c, j in enumerate(layout):
                    first.setdefault(j, c)
                self.elements.append(
                    (
                        b,
                        r,
                        routine,
                        slot,
                        [row[k] + 1 for k in layout],
                        [first[j] for j in range(len(row))],
                    )
                )
        self.nelm = len(self.elements)
        # Per ET slot: nodes listed, nodes with stresses (the corner nodes).
        self.slots = {}
        for b, r, routine, slot, nodes, _ in self.elements:
            entry = self.slots.setdefault(
                slot, {"routine": routine, "nodelm": 0, "nodstr": 0}
            )
            entry["nodelm"] = max(entry["nodelm"], len(nodes))
            entry["nodstr"] = max(entry["nodstr"], _num_corners(mesh.cells[b].type))
        self.components = []
        for reg in canonical_regions(mesh):
            if reg.kind == "side":
                continue
            ids = sorted(set(int(v) + 1 for v in np.asarray(reg.entries).ravel()))
            self.components.append((reg.name, reg.kind == "point", ids))
        if any(r.kind == "side" for r in canonical_regions(mesh)):
            warn(
                "Ansys .rst writer: side regions have no component equivalent; dropped"
            )
        self.mat = [_cell_int(mesh, "ansys:mat", e[0], e[1], 1) for e in self.elements]
        self.real = [
            _cell_int(mesh, "ansys:real", e[0], e[1], 1) for e in self.elements
        ]
        self.secnum = [
            _cell_int(mesh, "ansys:secnum", e[0], e[1], 1) for e in self.elements
        ]
        # Shells with a "top" array keep both surfaces (KEYOPT(8) = 0).
        self.layered = any(name.endswith(_TOP) for name in mesh.cell_data)
        self.signature = (
            n,
            [(e[2], e[3], e[4]) for e in self.elements],
            [(c[0], c[1], c[2]) for c in self.components],
            self.layered,
        )


class _Set:
    """One result set: its nodal solution, reactions and element results."""

    def __init__(self, mesh, model, time, step):
        n = model.nnod
        self.time = float(time)
        self.step = step
        columns = {}  # DOF code -> values (n,)
        for name, codes in _VECTORS.items():
            if name in mesh.point_data:
                a = np.asarray(mesh.point_data[name], dtype=np.float64)
                if a.ndim == 2 and a.shape == (n, 3):
                    for d, code in enumerate(codes):
                        columns[code] = a[:, d]
                else:
                    warn(
                        f"Ansys .rst writer: point data '{name}' is not (points, 3); dropped"
                    )
        for code, name in _SCALARS.items():
            if name in mesh.point_data:
                a = np.asarray(mesh.point_data[name], dtype=np.float64).reshape(-1)
                if len(a) == n:
                    columns[code] = a
                else:
                    warn(
                        f"Ansys .rst writer: point data '{name}' is not (points,); dropped"
                    )
        self.reactions = []  # (node position, DOF code, value)
        for name, codes in _REACTION_VECTORS.items():
            if name in mesh.point_data:
                a = np.asarray(mesh.point_data[name], dtype=np.float64)
                if a.shape != (n, 3):
                    warn(
                        f"Ansys .rst writer: point data '{name}' is not (points, 3); dropped"
                    )
                    continue
                for d, code in enumerate(codes):
                    for p in np.nonzero(np.isfinite(a[:, d]))[0].tolist():
                        self.reactions.append((p, code, float(a[p, d])))
        for name, code in sorted(_REACTION_SCALAR.items(), key=lambda kv: kv[1]):
            if name in mesh.point_data and name not in _REACTION_VECTORS:
                a = np.asarray(mesh.point_data[name], dtype=np.float64).reshape(-1)
                if len(a) == n:
                    for p in np.nonzero(np.isfinite(a))[0].tolist():
                        self.reactions.append((p, code, float(a[p])))
        used = sorted(set(columns) | {code for _, code, _ in self.reactions})
        self.dofs = used
        self.nodal = np.full((n, len(used)), _UNDEFINED)
        for k, code in enumerate(used):
            if code in columns:
                a = columns[code]
                self.nodal[:, k] = np.where(np.isfinite(a), a, _UNDEFINED)
        # Element results.
        self.tensors = []  # (name, tensor, is the top surface, per-block cell data)
        for base, entry, items, _ in _TENSORS:
            for suffix in ("", _TOP):
                name = base + suffix
                if name not in mesh.cell_data:
                    continue
                self.tensors.append((name, base, bool(suffix), mesh.cell_data[name]))
        # raw element records, and the element nodal forces (one value per DOF)
        self.raw = [
            (name, entry, mesh.cell_data[name])
            for name, entry in _RAW_RECORDS
            if name in mesh.cell_data
        ]
        self.enf = None
        if "ENF" in mesh.cell_data:
            layout = mesh.field_data.get("ansys:layout:ENF")
            if layout is not None and int(layout[1]) == len(used) and len(used):
                self.enf = mesh.cell_data["ENF"]
            else:
                warn(
                    "Ansys .rst writer: ENF needs its layout over the set's DOFs; dropped"
                )
        self.mask = 0
        if columns or self.reactions:
            self.mask |= 1 << _AVAIL["NSL"]
        if self.reactions:
            self.mask |= 1 << _AVAIL["RF"]
        for name, base, top, _ in self.tensors:
            self.mask |= 1 << _AVAIL[base]
        for name, _, _ in self.raw:
            self.mask |= 1 << _AVAIL[name]
        if self.enf is not None:
            self.mask |= 1 << _AVAIL["ENF"]


def _file_words(model, sets):
    f = _Words()
    nnod, nelm = model.nnod, model.nelm
    nsets = len(sets)
    resmax = max(nsets, 1)
    # -- the standard header (100 words) and the result header (80) --------------
    std = [0] * 100
    std[0], std[1], std[2] = 12, -1, 120000
    std[3], std[4] = 20000101, -1
    std[9] = _text_words(_VERSION, 1)[0]
    std[10] = 20140101
    for k, w in enumerate(_text_words("LINUX x64", 3)):
        std[11 + k] = w
    for k, w in enumerate(_text_words("file", 2)):
        std[14 + k] = w
    for k, w in enumerate(_text_words("FULL", 2)):
        std[16 + k] = w
    std[22] = 0
    std[25] = 16384
    for k, w in enumerate(_text_words("file", 8)):
        std[30 + k] = w
    f.ints(std)
    header = f.ints([12, nnod, nnod, resmax, 0, nelm, nelm, 0, nsets], 80)
    # words 9.. are patched below
    avail = 0
    for s in sets:
        avail |= s.mask
    f.patch(header, 35, avail)
    f.patch(header, 39, 1)  # rstsprs: ENS holds six items per node
    f.patch(header, 20, 1)  # one sector
    f.patch(header, 48, nnod)  # the whole model's nodes: not a distributed file
    neqv = f.ints(range(1, nnod + 1))
    f.patch64(header, 14, 45, neqv)
    eeqv = f.ints(range(1, nelm + 1))
    f.patch64(header, 13, 44, eeqv)
    numdof_max = max([len(s.dofs) for s in sets] + [0])
    f.patch(header, 4, numdof_max)

    # -- the geometry --------------------------------------------------------------
    maxety = max(model.slots) if model.slots else 0
    geometry = f.ints([0, maxety, 0, nnod, nelm], 80)
    f.patch64(header, 15, 46, geometry)
    f.patch(geometry, 42, nnod)
    f.patch(geometry, 45, nnod)
    ety = f.ints([0] * max(maxety, 1))
    f.patch64(geometry, 20, 21, ety)
    for slot, info in sorted(model.slots.items()):
        record = [0] * 200
        record[0], record[1] = slot, info["routine"]
        if info["routine"] in (181, 281):
            record[9] = 0 if model.layered else 1  # KEYOPT(8)
        record[60] = info["nodelm"]
        record[62] = info["nodelm"]
        record[93] = info["nodstr"]
        ptr = f.ints(record)
        f.patch(ety, slot - 1, ptr - ety)
    first_node = None
    for k, (x, y, z) in enumerate(model.points.tolist(), start=1):
        ptr = f.doubles([float(k), x, y, z, 0.0, 0.0, 0.0])
        first_node = ptr if first_node is None else first_node
    f.patch64(geometry, 26, 27, first_node)
    eid = f.ints([0] * (2 * nelm))
    f.patch64(geometry, 28, 29, eid)
    for k, (b, r, routine, slot, nodes, _) in enumerate(model.elements):
        record = [
            model.mat[k],
            slot,
            model.real[k],
            model.secnum[k],
            0,
            1,
            0,
            0,
            k + 1,
            0,
        ]
        ptr = f.ints(record + nodes)
        f.patch(eid, 2 * k, ptr - eid)
    f.patch(geometry, 48, len(model.components))
    if model.components:
        first = None
        for name, is_nodes, ids in model.components:
            packed = []
            k = 0
            while k < len(ids):
                j = k
                while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
                    j += 1
                packed.append(ids[k])
                if j > k:
                    packed.append(-ids[j])
                k = j + 1
            ptr = f.ints(
                [1 if is_nodes else 2, *_text_words(name, _NAME_WORDS), *packed]
            )
            first = ptr if first is None else first
        f.patch64(geometry, 50, 51, first)

    # -- the result sets -----------------------------------------------------------
    bases = []
    for s in sets:
        bases.append(_set_records(f, model, s))
    dsi = f.ints(
        [b & 0xFFFFFFFF for b in bases]
        + [0] * (resmax - nsets)
        + [b >> 32 for b in bases]
        + [0] * (resmax - nsets)
    )
    tim = f.doubles([s.time for s in sets] + [0.0] * (resmax - nsets))
    lsp = f.ints([v for s in sets for v in s.step] + [0] * (3 * (resmax - nsets)))
    f.patch64(header, 10, 40, dsi)
    f.patch64(header, 11, 41, tim)
    f.patch64(header, 12, 42, lsp)
    end = len(f.w)
    f.patch(header, 9, end)
    f.patch(header, 22, end)
    return f


def _set_records(f, model, s):
    nnod = model.nnod
    numdof = len(s.dofs)
    base = f.ints([0] * 150)
    f.patch(base, 1, model.nelm)
    f.patch(base, 2, nnod)
    f.patch(base, 3, 0)
    f.patch(base, 19, numdof)
    for k, code in enumerate(s.dofs):
        f.patch(base, 20 + k, code)
    f.patch(base, 146, s.mask)
    if numdof:
        nsl = f.doubles(s.nodal.reshape(-1))
        f.patch64(base, 104, 105, nsl - base)
    if s.reactions:
        order = {code: k for k, code in enumerate(s.dofs)}
        table = []
        for p, code, _ in s.reactions:
            index = p * numdof + order[code] + 1
            table += [index & 0xFFFFFFFF, index >> 32]
        rf = f.ints(table)
        f.doubles([v for _, _, v in s.reactions])
        f.patch(base, 7, len(s.reactions))
        f.patch64(base, 106, 107, rf - base)
    if s.tensors or s.raw or s.enf is not None:
        _element_solution(f, model, s, base)
    return base


def _element_solution(f, model, s, base):
    esl = f.ints([0] * (2 * model.nelm))
    f.patch64(base, 118, 119, esl - base)
    numdof = len(s.dofs)
    for k, (b, r, routine, slot, nodes, cell_slots) in enumerate(model.elements):
        entries = {}
        pending = []
        info = model.slots[slot]
        # raw records: the values up to the last one the element has
        for name, entry, per_block in s.raw:
            arr = np.asarray(per_block[b], dtype=np.float64)
            if arr.ndim != 2 or len(arr) <= r:
                continue
            finite = np.nonzero(np.isfinite(arr[r]))[0]
            if len(finite) == 0:
                continue
            row = arr[r, : finite[-1] + 1]
            pending.append((entry, np.where(np.isfinite(row), row, _UNDEFINED)))
        if model.has_results[b]:
            for name, base_name, top, per_block in s.tensors:
                arr = np.asarray(per_block[b])
                if arr.ndim != 2 or len(arr) <= r:
                    continue
                row = arr[r]
                width = 6
                npc = row.size // width
                if npc == 0 or row.size % width:
                    continue
                entries.setdefault(base_name, {})[top] = row.reshape(npc, width)
            nodstr = info["nodstr"]
            layers = 2 if (routine in (181, 281) and model.layered) else 1
            for base_name, surfaces in entries.items():
                items = _ITEMS[base_name]
                values = np.full((nodstr * layers, items), _UNDEFINED)
                values[:, 6:] = 0.0
                any_value = False
                for top, data in surfaces.items():
                    layer = 1 if (top and layers == 2) else 0
                    if top and layers == 1:
                        continue
                    for j, slot_j in enumerate(cell_slots):
                        if j >= len(data) or slot_j >= nodstr:
                            continue
                        row6 = data[j]
                        if np.isfinite(row6).all():
                            values[layer * nodstr + slot_j, :6] = row6
                            any_value = True
                if any_value:
                    pending.append((_TABLE_ENTRIES[base_name], values.reshape(-1)))
            if s.enf is not None and info["nodelm"]:
                arr = np.asarray(s.enf[b], dtype=np.float64)
                if arr.ndim == 2 and len(arr) > r and arr.shape[1] % numdof == 0:
                    npc = arr.shape[1] // numdof
                    forces = np.zeros((info["nodelm"], numdof))
                    any_value = False
                    for j, slot_j in enumerate(cell_slots):
                        if j < npc and slot_j < info["nodelm"]:
                            row_j = arr[r, j * numdof : (j + 1) * numdof]
                            if np.isfinite(row_j).all():
                                forces[slot_j] = row_j
                                any_value = True
                    if any_value:
                        pending.append((_ENF, forces.reshape(-1)))
        if not pending:
            continue
        table = [0] * 25
        table_ptr = f.ints(table)
        for entry, values in pending:
            ptr = f.doubles(values)
            f.patch(table_ptr, entry, ptr - table_ptr)
        f.patch(esl, 2 * k, table_ptr - esl)


def _time_of(mesh):
    fd = getattr(mesh, "field_data", {}) or {}
    if "meshio:time" in fd:
        return float(np.asarray(fd["meshio:time"]).reshape(-1)[0])
    return 1.0


def _step_of(mesh, index):
    fd = getattr(mesh, "field_data", {}) or {}

    def get(name, default):
        if name in fd:
            return int(np.asarray(fd[name]).reshape(-1)[0])
        return default

    return (
        get("ansys:load_step", 1),
        get("ansys:substep", index + 1),
        get("ansys:cumulative", index + 1),
    )


def write(filename, mesh):
    """Write a results file with one result set: the mesh, its nodal DOF
    solution, reaction forces and element nodal stresses and strains."""
    if is_buffer(filename, "w"):
        raise WriteError("Ansys .rst: cannot be written to a stream")
    model = _Model(mesh)
    sets = [_Set(mesh, model, _time_of(mesh), _step_of(mesh, 0))]
    with open(filename, "wb") as fh:
        fh.write(_file_words(model, sets).tobytes())


class SeriesWriter:
    """A result set per mesh, written when the writer closes; every mesh must
    have the first one's points, cells and components."""

    def __init__(self, filename):
        if is_buffer(filename, "w"):
            raise WriteError("Ansys .rst: cannot be written to a stream")
        self._path = filename
        self._model = None
        self._sets = []

    def write(self, time, mesh):
        model = _Model(mesh)
        if self._model is None:
            self._model = model
        elif model.signature != self._model.signature:
            raise WriteError(
                f"Ansys .rst: series: step {len(self._sets)}'s points, cells or "
                "components differ from the first step's; a results file holds one "
                "model -- write one file per step with '{step}'"
            )
        self._sets.append(
            _Set(mesh, self._model, float(time), _step_of(mesh, len(self._sets)))
        )

    def close(self):
        if self._sets:
            sets, self._sets = self._sets, []
            with open(self._path, "wb") as fh:
                fh.write(_file_words(self._model, sets).tobytes())

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
