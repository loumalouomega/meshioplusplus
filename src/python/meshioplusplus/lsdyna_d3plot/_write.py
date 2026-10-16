"""Write LS-DYNA state databases (the ``d3plot`` family).

The pure-Python twin of ``write_lsdyna_d3plot`` in
``src/cpp/src/formats/lsdyna_d3plot.cpp``: both engines write the same bytes.
The layout is the one ``_d3plot.py`` reads (4-byte little-endian words): the
64-word control block, the geometry (coordinates, then solid, beam and shell
connectivity with the part index as the last word), the user ids, the end mark,
the part titles and one state per step; the states are split over ``d3plot``,
``d3plot01``... by size.

Mesh to file: hexahedra, wedges, pyramids and tetrahedra are solids (the
degenerate ones as bricks with repeated nodes), ``tetra10`` ten-node solids,
``line`` beams, ``triangle`` and ``quad`` shells. Other cell types are not
written (a warning counts them). Point data ``displacement``, ``velocity`` and
``acceleration``; cell data ``stress``, ``effective_plastic_strain`` (solids: 1
or 8 points, shells: through-thickness layers), the ``beam_*`` resultants and
the deletion mask ``lsdyna:alive``; the ``global_*`` and ``part_*`` field data.
"""

from __future__ import annotations

import os

import numpy as np

from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from .._regions import canonical_regions
from ..lsdyna._lsdyna import _expand_solid
from ._d3plot import family_files

__all__ = ["write", "SeriesWriter"]

_TITLE = b"meshio++"
_VERSION = 971.0
_END_MARK = -999999.0
_DEFAULT_FILE_BYTES = 1 << 30
_SOLID, _BEAM, _SHELL = range(3)
_FAMILY = {
    "hexahedron": _SOLID,
    "wedge": _SOLID,
    "pyramid": _SOLID,
    "tetra": _SOLID,
    "tetra10": _SOLID,
    "line": _BEAM,
    "triangle": _SHELL,
    "quad": _SHELL,
}
_BEAM_ARRAYS = (
    ("beam_axial_force", 1),
    ("beam_shear_force", 2),
    ("beam_bending_moment", 2),
    ("beam_torsion_moment", 1),
)
_GLOBALS = (
    ("global_kinetic_energy", 1),
    ("global_internal_energy", 1),
    ("global_total_energy", 1),
    ("global_velocity", 3),
)
_PART_GLOBALS = (
    ("part_internal_energy", 1),
    ("part_kinetic_energy", 1),
    ("part_velocity", 3),
    ("part_mass", 1),
    ("part_hourglass_energy", 1),
)


def _i32(values):
    return np.asarray(values, dtype="<i4").tobytes()


def _f32(values):
    return np.asarray(values, dtype="<f4").tobytes()


def _time_of(mesh):
    value = mesh.field_data.get("meshio:time")
    return float(np.asarray(value).reshape(-1)[0]) if value is not None else 0.0


def _cell_array(mesh, name, comps):
    """Per block ``(cells, points, comps)`` of a cell array, or ``None``."""
    values = mesh.cell_data.get(name)
    if values is None:
        return None
    layout = mesh.field_data.get(f"lsdyna_d3plot:layout:{name}")
    out = []
    for block, a in zip(mesh.cells, values):
        a = np.asarray(a, dtype=np.float64)
        n = len(block.data)
        if n == 0 or a.size % n:
            out.append(None)
            continue
        a = a.reshape(n, -1)
        width = a.shape[1]
        if (
            layout is not None
            and int(layout[1]) == comps
            and width == int(layout[0]) * comps
        ):
            points = int(layout[0])
        elif width % comps == 0 and width > 0:
            points = width // comps
        else:
            out.append(None)
            continue
        out.append(a.reshape(n, points, comps))
    return out


class _Model:
    """The geometry of a file and the layout of its states."""

    def __init__(self, mesh):
        points = np.asarray(mesh.points, dtype=np.float64)
        n = len(points)
        points = points.reshape(n, -1)
        if n == 0:
            raise WriteError("LS-DYNA d3plot: the mesh has no points")
        self.points = np.zeros((n, 3))
        self.points[:, : min(3, points.shape[1])] = points[:, :3]
        self.n = n
        nid = mesh.point_data.get("lsdyna:nid")
        self.nid = (
            np.asarray(nid, dtype=np.int64).reshape(-1)
            if nid is not None and len(nid) == n
            else np.arange(1, n + 1, dtype=np.int64)
        )
        # the cells of each family: (block, first row in the block, count)
        self.ranges = [[], [], []]
        self.rows = [[], [], []]  # per family: node lists in file order
        self.extra = []  # tetra10: the two extra nodes per solid
        self.cell_part = []  # per family list of part ids
        self.cell_eid = [[], [], []]
        part_of_cell = self._part_assignment(mesh)
        skipped = {}
        base = 0
        eid_arrays = mesh.cell_data.get("lsdyna:eid")
        part_arrays = mesh.cell_data.get("lsdyna:part")
        self.parts_of = [[], [], []]
        for b, block in enumerate(mesh.cells):
            data = np.asarray(block.data, dtype=np.int64)
            family = _FAMILY.get(block.type)
            if family is None:
                if len(data):
                    skipped[block.type] = skipped.get(block.type, 0) + len(data)
                base += len(data)
                continue
            if len(data) and (data.min() < 0 or data.max() >= n):
                raise WriteError(
                    f"LS-DYNA d3plot: a '{block.type}' cell refers to a point the mesh "
                    "does not have"
                )
            first = len(self.rows[family])
            for r, row in enumerate(data):
                if family == _SOLID and block.type == "tetra10":
                    nodes = [int(v) for v in row[:8]]
                    self.extra.append(
                        (len(self.rows[family]), [int(v) for v in row[8:10]])
                    )
                elif family == _SOLID:
                    nodes = _expand_solid(block.type, row)
                else:
                    nodes = [int(v) for v in row]
                    if block.type == "triangle":
                        nodes = nodes + [nodes[2]]
                self.rows[family].append(nodes)
                if (
                    part_arrays is not None
                    and b < len(part_arrays)
                    and part_arrays[b] is not None
                ):
                    pid = int(np.asarray(part_arrays[b]).reshape(-1)[r])
                else:
                    pid = part_of_cell[base + r]
                self.parts_of[family].append(pid)
                if eid_arrays is not None and b < len(eid_arrays):
                    self.cell_eid[family].append(
                        int(np.asarray(eid_arrays[b]).reshape(-1)[r])
                    )
                else:
                    self.cell_eid[family].append(len(self.rows[family]))
            self.ranges[family].append((b, first, len(data)))
            base += len(data)
        for ctype, count in skipped.items():
            warn(f"LS-DYNA d3plot: {count} '{ctype}' cell(s) are not written")
        self.counts = [len(r) for r in self.rows]
        if not any(self.counts):
            raise WriteError(
                "LS-DYNA d3plot: the mesh has no solid, beam or shell cells"
            )
        # every part, sorted; each family's number of distinct parts
        self.part_ids = sorted({p for fam in self.parts_of for p in fam})
        self.part_index = {p: i + 1 for i, p in enumerate(self.part_ids)}
        self.nummat = [len(set(fam)) for fam in self.parts_of]
        self.titles = self._titles
        self.signature = (n, [tuple(map(tuple, r)) for r in self.rows])

    def _part_assignment(self, mesh):
        """A part id for every cell, from the cell regions (the first region
        containing a cell names its part); cells in none share one part."""
        regions = [r for r in canonical_regions(mesh) if r.kind == "cell"]
        ncells = sum(len(b.data) for b in mesh.cells)
        taken = {int(r.tag) for r in regions if r.tag > 0}
        fresh = 1
        part = np.full(ncells, -1, dtype=np.int64)
        self._titles = {}
        for r in regions:
            pid = int(r.tag)
            if pid <= 0:
                while fresh in taken:
                    fresh += 1
                pid = fresh
                taken.add(pid)
            self._titles.setdefault(pid, r.name)
            idx = np.asarray(r.entries, dtype=np.int64).reshape(-1)
            idx = idx[(idx >= 0) & (idx < ncells)]
            idx = idx[part[idx] < 0]
            part[idx] = pid
        rest = (max(taken) + 1) if taken else 1
        part[part < 0] = rest
        return part


class _Plan:
    """The header: what the states of all steps need."""

    def __init__(self, model, meshes):
        self.model = model
        self.solid_points = 1
        self.shell_points = 0
        self.stress = self.pstrain = False
        self.disp = self.vel = self.acc = False
        self.alive = False
        self.beam = False
        self.globals_ = False
        self.part_globals = False
        nparts = sum(model.nummat)
        for mesh in meshes:
            self.disp = self.disp or "displacement" in mesh.point_data
            self.vel = self.vel or "velocity" in mesh.point_data
            self.acc = self.acc or "acceleration" in mesh.point_data
            self.alive = self.alive or "lsdyna:alive" in mesh.cell_data
            self.globals_ = self.globals_ or any(
                k in mesh.field_data for k, _ in _GLOBALS
            )
            if self.model.counts[_BEAM]:
                self.beam = self.beam or any(
                    k in mesh.cell_data for k, _ in _BEAM_ARRAYS
                )
            if all(
                k in mesh.field_data
                and np.asarray(mesh.field_data[k]).size == w * nparts
                for k, w in _PART_GLOBALS
            ):
                self.part_globals = True
            for name, comps in (("stress", 6), ("effective_plastic_strain", 1)):
                arrays = _cell_array(mesh, name, comps)
                if arrays is None:
                    continue
                for family in (_SOLID, _SHELL):
                    for b, _, _ in model.ranges[family]:
                        a = arrays[b]
                        if a is None:
                            continue
                        points = a.shape[1]
                        if family == _SOLID:
                            if points not in (1, 8):
                                continue
                            self.solid_points = max(self.solid_points, points)
                        else:
                            self.shell_points = max(self.shell_points, points)
                        if name == "stress":
                            self.stress = True
                        else:
                            self.pstrain = True
        self.globals_ = self.globals_ or self.part_globals
        self.nglbv = (
            (6 + (7 * nparts if self.part_globals else 0)) if self.globals_ else 0
        )
        per_point = 6 * self.stress + self.pstrain
        self.nv3d = self.solid_points * per_point if model.counts[_SOLID] else 0
        self.nv2d = self.shell_points * per_point if model.counts[_SHELL] else 0
        self.nv1d = 6 if self.beam else 0
        # a state of the time alone cannot be told from the end mark
        self.words = self._state_words()
        if self.words < 2:
            self.disp = True
            self.words = self._state_words()

    def _state_words(self):
        m = self.model
        n = 1 + self.nglbv + 3 * m.n * (self.disp + self.vel + self.acc)
        n += m.counts[_SOLID] * self.nv3d + m.counts[_BEAM] * self.nv1d
        n += m.counts[_SHELL] * self.nv2d
        if self.alive:
            n += sum(m.counts)
        return n

    def maxint(self):
        layers = max(self.shell_points, 1)
        return -(10000 + layers) if self.alive else layers


def _geometry(plan):
    m = plan.model
    words = [0] * 64
    fields = {
        11: 1,
        15: 4,
        16: m.n,
        17: 6,
        18: plan.nglbv,
        20: int(plan.disp),
        21: int(plan.vel),
        22: int(plan.acc),
        23: -m.counts[_SOLID] if m.extra else m.counts[_SOLID],
        24: m.nummat[_SOLID],
        27: plan.nv3d,
        28: m.counts[_BEAM],
        29: m.nummat[_BEAM],
        30: plan.nv1d,
        31: m.counts[_SHELL],
        32: m.nummat[_SHELL],
        33: plan.nv2d,
        36: plan.maxint(),
        43: 1000 if plan.stress else 0,
        44: 1000 if plan.pstrain else 0,
        51: len(m.part_ids),
    }
    ids = m.counts[_SOLID] + m.counts[_BEAM] + m.counts[_SHELL]
    narbs = 16 + m.n + ids + 3 * len(m.part_ids)
    fields[39] = narbs
    for i, v in fields.items():
        words[i] = v
    int_part = np.asarray(words, dtype="<i4")
    out = bytearray(int_part.tobytes())
    out[:40] = _TITLE + b" " * (40 - len(_TITLE))
    out[14 * 4 : 15 * 4] = _f32([_VERSION])
    out += _f32(m.points.reshape(-1))
    solids = np.zeros((m.counts[_SOLID], 9), dtype=np.int64)
    for e, nodes in enumerate(m.rows[_SOLID]):
        solids[e, :8] = np.asarray(nodes) + 1
        solids[e, 8] = m.part_index[m.parts_of[_SOLID][e]]
    out += _i32(solids)
    if m.extra:
        extra = np.zeros((m.counts[_SOLID], 2), dtype=np.int64)
        for e, nodes in m.extra:
            extra[e] = np.asarray(nodes) + 1
        out += _i32(extra)
    beams = np.zeros((m.counts[_BEAM], 6), dtype=np.int64)
    for e, nodes in enumerate(m.rows[_BEAM]):
        beams[e, :2] = np.asarray(nodes[:2]) + 1
        beams[e, 5] = m.part_index[m.parts_of[_BEAM][e]]
    out += _i32(beams)
    shells = np.zeros((m.counts[_SHELL], 5), dtype=np.int64)
    for e, nodes in enumerate(m.rows[_SHELL]):
        shells[e, :4] = np.asarray(nodes) + 1
        shells[e, 4] = m.part_index[m.parts_of[_SHELL][e]]
    out += _i32(shells)
    # user ids: the header (NSORT < 0), then nodes, solids, beams, shells,
    # thick shells (none) and the parts
    nparts = len(m.part_ids)
    solids_at = 17 + m.n
    beams_at = solids_at + m.counts[_SOLID]
    shells_at = beams_at + m.counts[_BEAM]
    parts_at = shells_at + m.counts[_SHELL]
    out += _i32(
        [-17, solids_at, beams_at, shells_at, parts_at]
        + [m.n, m.counts[_SOLID], m.counts[_BEAM], m.counts[_SHELL], 0]
    )
    out += _i32([parts_at, parts_at + nparts, parts_at + 2 * nparts, nparts, 0, nparts])
    out += _i32(m.nid)
    for family in (_SOLID, _BEAM, _SHELL):
        out += _i32(m.cell_eid[family])
    # the parts' ids, then (as LS-DYNA writes them) in defined order and cross-referenced
    out += _i32(m.part_ids) * 2 + _i32(range(1, nparts + 1))
    out += _f32([_END_MARK])
    titles = bytearray()
    for pid in m.part_ids:
        text = m.titles.get(pid, f"Part {pid}").encode("latin-1", "replace")[:72]
        titles += _i32([pid]) + text + b" " * (72 - len(text))
    out += _i32([90001, len(m.part_ids)]) + bytes(titles)
    out += _f32([_END_MARK])
    return bytes(out)


def _node_array(mesh, name, n):
    value = mesh.point_data.get(name)
    if value is None:
        return np.zeros((n, 3))
    a = np.asarray(value, dtype=np.float64).reshape(n, -1)
    out = np.zeros((n, 3))
    out[:, : min(3, a.shape[1])] = a[:, :3]
    return out


def _element_block(plan, mesh, family, count, nv, points):
    """The (cells, points, per point) variables of a family, NaN where absent."""
    m = plan.model
    per_point = 6 * plan.stress + plan.pstrain
    out = np.full((count, points, per_point), np.nan)
    column = 0
    for name, comps, wanted in (
        ("stress", 6, plan.stress),
        ("effective_plastic_strain", 1, plan.pstrain),
    ):
        if not wanted:
            continue
        arrays = _cell_array(mesh, name, comps)
        if arrays is not None:
            for b, first, n in m.ranges[family]:
                a = arrays[b]
                if a is None or n == 0:
                    continue
                if family == _SOLID and a.shape[1] not in (1, 8):
                    continue
                p = min(a.shape[1], points)
                out[first : first + n, :p, column : column + comps] = a[:, :p, :]
        column += comps
    return out.reshape(count, nv)


def _state(plan, mesh, time):
    m = plan.model
    n = m.n
    pieces = [_f32([time])]
    if plan.nglbv:
        g = np.zeros(plan.nglbv)
        at = 0
        for name, width in _GLOBALS:
            v = mesh.field_data.get(name)
            if v is not None:
                vals = np.asarray(v, dtype=np.float64).reshape(-1)[:width]
                g[at : at + len(vals)] = vals
            at += width
        if plan.part_globals:
            nparts = sum(m.nummat)
            for name, width in _PART_GLOBALS:
                v = mesh.field_data.get(name)
                if v is not None and np.asarray(v).size == width * nparts:
                    g[at : at + width * nparts] = np.asarray(
                        v, dtype=np.float64
                    ).reshape(-1)
                at += width * nparts
        pieces.append(_f32(g))
    if plan.disp:
        pieces.append(_f32(m.points + _node_array(mesh, "displacement", n)))
    if plan.vel:
        pieces.append(_f32(_node_array(mesh, "velocity", n)))
    if plan.acc:
        pieces.append(_f32(_node_array(mesh, "acceleration", n)))
    if plan.nv3d:
        pieces.append(
            _f32(
                _element_block(
                    plan, mesh, _SOLID, m.counts[_SOLID], plan.nv3d, plan.solid_points
                )
            )
        )
    if plan.nv1d:
        beams = np.full((m.counts[_BEAM], 6), np.nan)
        at = 0
        for name, comps in _BEAM_ARRAYS:
            arrays = _cell_array(mesh, name, comps)
            if arrays is not None:
                for b, first, count in m.ranges[_BEAM]:
                    a = arrays[b]
                    if a is not None and count:
                        beams[first : first + count, at : at + comps] = a[:, 0, :]
            at += comps
        pieces.append(_f32(beams))
    if plan.nv2d:
        pieces.append(
            _f32(
                _element_block(
                    plan, mesh, _SHELL, m.counts[_SHELL], plan.nv2d, plan.shell_points
                )
            )
        )
    if plan.alive:
        alive = mesh.cell_data.get("lsdyna:alive")
        for family in (_SOLID, _SHELL, _BEAM):
            flags = np.ones(m.counts[family])
            if alive is not None:
                for b, first, count in m.ranges[family]:
                    if b < len(alive) and count:
                        flags[first : first + count] = (
                            np.asarray(alive[b]).reshape(-1) != 0
                        )
            pieces.append(_f32(flags))
    return b"".join(pieces)


def _assemble(plan, steps, max_bytes):
    """The files of the family as ``(suffix number, bytes)``."""
    head = _geometry(plan)
    mark = _f32([_END_MARK])
    state_bytes = plan.words * 4
    files = []
    current = bytearray(head)
    used = 0
    for time, mesh in steps:
        state = _state(plan, mesh, time)
        if used and len(current) + len(state) + len(mark) > max_bytes:
            files.append(bytes(current) + mark)
            current = bytearray()
            used = 0
        current += state
        used += 1
        assert len(state) == state_bytes
    files.append(bytes(current) + mark)
    return files


def _store(path, files):
    for stale in family_files(path)[1:]:
        try:
            os.remove(stale)
        except OSError:
            pass
    for k, data in enumerate(files):
        target = path if k == 0 else f"{path}{k:02d}"
        with open(target, "wb") as fh:
            fh.write(data)


def write(filename, mesh, max_file_bytes=_DEFAULT_FILE_BYTES):
    """Write one state as a d3plot family (a single file for one state)."""
    if is_buffer(filename, "w"):
        raise WriteError("LS-DYNA d3plot: cannot be written to a stream")
    model = _Model(mesh)
    plan = _Plan(model, [mesh])
    _store(
        os.fspath(filename), _assemble(plan, [(_time_of(mesh), mesh)], max_file_bytes)
    )


class SeriesWriter:
    """A d3plot family with a state per step; the geometry is the first step's,
    which every later step must share. The header takes the union of the
    variables of all steps, so the file is written when the writer closes."""

    def __init__(self, filename, max_file_bytes=_DEFAULT_FILE_BYTES):
        if is_buffer(filename, "w"):
            raise WriteError("LS-DYNA d3plot: cannot be written to a stream")
        self._path = os.fspath(filename)
        self._max = int(max_file_bytes)
        self._model = None
        self._steps = []
        self._closed = False

    def write(self, time, mesh):
        model = _Model(mesh)
        if self._model is None:
            self._model = model
        elif model.signature != self._model.signature:
            raise WriteError(
                f"LS-DYNA d3plot: series: step {len(self._steps)}'s points or cells "
                "differ from the first step's; a family holds one mesh -- write one "
                "family per step with '{step}'"
            )
        self._steps.append((float(time), mesh))

    def close(self):
        if self._closed:
            return
        self._closed = True
        if self._model is None:
            return
        plan = _Plan(self._model, [mesh for _, mesh in self._steps])
        _store(self._path, _assemble(plan, self._steps, self._max))

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
