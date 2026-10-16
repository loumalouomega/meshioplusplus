"""Write LS-DYNA state databases (the ``d3plot`` family).

The pure-Python twin of ``write_lsdyna_d3plot`` in
``src/cpp/src/formats/lsdyna_d3plot.cpp``: both engines write the same bytes.
The layout is the one ``_d3plot.py`` reads (little-endian words of 4 bytes, or 8
in double precision): the control block, the geometry (coordinates, then solid,
beam and shell connectivity with the part index as the last word), the user ids,
the end mark, the part titles and one state per step; the states are split over
``d3plot``, ``d3plot01``... by size.

Mesh to file: hexahedra, wedges, pyramids and tetrahedra are solids (the
degenerate ones as bricks with repeated nodes), ``tetra10`` ten-node solids,
``line`` beams, ``triangle`` and ``quad`` shells. Other cell types are not
written (a warning counts them). The variables are those the reader returns, by
the same names: nodal ``displacement``, ``velocity``, ``acceleration``,
``temperature``, ``heat_flux``, ``mass_scaling``, ``temperature_gradient``,
``residual_forces`` and ``residual_moments``; the solid, shell and beam element
variables (``stress``, ``effective_plastic_strain``, ``strain``, the strain
tensors, ``history_variables``, ``thermal_variables``, the shell resultants and
extras, the beam resultants and layers), the deletion mask ``lsdyna:alive`` and
the ``global_*`` and ``part_*`` field data. A variable a block does not have is
NaN.
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
_SOLID, _BEAM, _SHELL, _SPH = range(4)
_FAMILY = {
    "hexahedron": _SOLID,
    "wedge": _SOLID,
    "pyramid": _SOLID,
    "tetra": _SOLID,
    "tetra10": _SOLID,
    "hexahedron20": _SOLID,
    "hexahedron27": _SOLID,
    "line": _BEAM,
    "triangle": _SHELL,
    "quad": _SHELL,
    "quad8": _SHELL,
    "vertex": _SPH,
}
# The default components of the element variables; 0 takes them from the layout
# array.
_COMPS = {
    "stress": 6,
    "effective_plastic_strain": 1,
    "plastic_strain_tensor": 6,
    "thermal_strain_tensor": 6,
    "strain": 6,
    "thermal_variables": 0,
    "history_variables": 0,
    "shell_bending_moment": 3,
    "shell_shear_force": 2,
    "shell_normal_force": 3,
    "thickness": 1,
    "shell_element_variables": 2,
    "strain_inner": 6,
    "strain_outer": 6,
    "internal_energy": 1,
    "beam_axial_force": 1,
    "beam_shear_force": 2,
    "beam_bending_moment": 2,
    "beam_torsion_moment": 1,
    "beam_axial_stress": 1,
    "beam_shear_stress": 2,
    "beam_axial_strain": 1,
    "sph_radius": 1,
    "sph_pressure": 1,
    "density": 1,
    "sph_neighbors": 1,
    "strain_rate": 6,
    "mass": 1,
}
# SPH particle variables in state order: (name, flag index, components); the
# strain flag is signed in the file and counts the strain rate too.
_SPH_VARS = (
    ("sph_radius", 1, 1),
    ("sph_pressure", 2, 1),
    ("stress", 3, 6),
    ("effective_plastic_strain", 4, 1),
    ("density", 5, 1),
    ("internal_energy", 6, 1),
    ("sph_neighbors", 7, 1),
)
_RIGID_FIELDS = (
    ("coordinates", 3),
    ("rotation", 9),
    ("velocity", 3),
    ("rotational_velocity", 3),
    ("acceleration", 3),
    ("rotational_acceleration", 3),
)
_BEAM_CORE = (
    ("beam_axial_force", 1),
    ("beam_shear_force", 2),
    ("beam_bending_moment", 2),
    ("beam_torsion_moment", 1),
)
# per beam layer: axial stress, shear stress (2), plastic strain, axial strain
_BEAM_LAYER = (
    ("beam_axial_stress", 1),
    ("beam_shear_stress", 2),
    ("effective_plastic_strain", 1),
    ("beam_axial_strain", 1),
)
_SHELL_FORCES = (
    ("shell_bending_moment", 3),
    ("shell_shear_force", 2),
    ("shell_normal_force", 3),
)
_SHELL_EXTRA = (("thickness", 1), ("shell_element_variables", 2))
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
_NODAL_NAMES = (
    "displacement",
    "velocity",
    "acceleration",
    "temperature",
    "heat_flux",
    "mass_scaling",
    "temperature_gradient",
    "residual_forces",
    "residual_moments",
)


class _Words:
    """Little-endian words of 4 or 8 bytes."""

    def __init__(self, double):
        self.ws = 8 if double else 4
        self.itype = "<i8" if double else "<i4"
        self.ftype = "<f8" if double else "<f4"

    def i(self, values):
        return np.asarray(values, dtype=self.itype).tobytes()

    def f(self, values):
        return np.asarray(values, dtype=self.ftype).tobytes()


def _time_of(mesh):
    value = mesh.field_data.get("meshio:time")
    return float(np.asarray(value).reshape(-1)[0]) if value is not None else 0.0


def _cell_array(mesh, name, comps=None):
    """Per block ``(cells, points, comps)`` of a cell array, or ``None``."""
    values = mesh.cell_data.get(name)
    if values is None:
        return None
    layout = mesh.field_data.get(f"lsdyna_d3plot:layout:{name}")
    default = _COMPS.get(name, 0) if comps is None else comps
    out = []
    for block, a in zip(mesh.cells, values):
        a = np.asarray(a, dtype=np.float64)
        n = len(block.data)
        if n == 0 or a.size % n:
            out.append(None)
            continue
        a = a.reshape(n, -1)
        width = a.shape[1]
        if layout is not None and width == int(layout[0]) * int(layout[1]):
            points, k = int(layout[0]), int(layout[1])
        elif default and width % default == 0 and width > 0:
            points, k = width // default, default
        elif not default and width > 0:
            points, k = 1, width
        else:
            out.append(None)
            continue
        if comps is not None and k != comps:
            out.append(None)
            continue
        out.append(a.reshape(n, points, k))
    return out


def _node_array(mesh, name, n, comps=3):
    """A point array as ``(n, comps)``: missing components and arrays are NaN."""
    value = mesh.point_data.get(name)
    out = np.full((n, comps), np.nan)
    if value is None:
        return out
    a = np.asarray(value, dtype=np.float64).reshape(n, -1)
    k = min(comps, a.shape[1])
    out[:, :k] = a[:, :k]
    return out


class _Model:
    """The geometry of a file."""

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
        self.ranges = [[], [], [], []]
        self.rows = [[], [], [], []]  # per family: node lists in file order
        self.extra = []  # tetra10: the two extra nodes per solid
        self.extra20 = []  # hexahedron20: (solid, the 12 edge nodes)
        self.extra27 = []  # hexahedron27: (solid, nodes 8 to 26)
        self.extra8 = []  # quad8: (shell, the 4 midside nodes)
        self.cell_eid = [[], [], [], []]
        part_of_cell = self._part_assignment(mesh)
        skipped = {}
        base = 0
        eid_arrays = mesh.cell_data.get("lsdyna:eid")
        part_arrays = mesh.cell_data.get("lsdyna:part")
        self.parts_of = [[], [], [], []]
        airbag = 0
        for b, block in enumerate(mesh.cells):
            data = np.asarray(block.data, dtype=np.int64)
            family = _FAMILY.get(block.type)
            if family == _SPH and len(data):
                # an airbag particle reads back on a point with no coordinates yet
                if len(data) and np.isnan(self.points[data[:, 0]]).any():
                    airbag += len(data)
                    base += len(data)
                    continue
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
                elif family == _SOLID and block.type == "hexahedron20":
                    nodes = [int(v) for v in row[:8]]
                    self.extra20.append(
                        (len(self.rows[family]), [int(v) for v in row[8:20]])
                    )
                elif family == _SOLID and block.type == "hexahedron27":
                    nodes = [int(v) for v in row[:8]]
                    self.extra27.append(
                        (len(self.rows[family]), [int(v) for v in row[8:27]])
                    )
                elif family == _SOLID:
                    nodes = _expand_solid(block.type, row)
                else:
                    nodes = [int(v) for v in row]
                    if block.type == "triangle":
                        nodes = nodes + [nodes[2]]
                    elif block.type == "quad8":
                        self.extra8.append((len(self.rows[family]), nodes[4:8]))
                        nodes = nodes[:4]
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
        if airbag:
            warn(f"LS-DYNA d3plot: {airbag} airbag particle(s) are not written")
        self.counts = [len(r) for r in self.rows]
        if not any(self.counts):
            raise WriteError(
                "LS-DYNA d3plot: the mesh has no solid, beam, shell or SPH cells"
            )
        # rigid bodies: their parts (user ids) from the reader's field data
        rigid = mesh.field_data.get("lsdyna:rigid_body_part")
        self.rigid_parts = (
            [int(v) for v in np.asarray(rigid).reshape(-1)] if rigid is not None else []
        )
        # every part, sorted; each family's number of distinct parts (SPH
        # particles count with the solids, as LS-DYNA counts them)
        self.part_ids = sorted(
            {p for fam in self.parts_of for p in fam} | set(self.rigid_parts)
        )
        self.part_index = {p: i + 1 for i, p in enumerate(self.part_ids)}
        self.nummat = [
            len(set(self.parts_of[_SOLID]) | set(self.parts_of[_SPH])),
            len(set(self.parts_of[_BEAM])),
            len(set(self.parts_of[_SHELL])),
        ]
        self.nrigid = len(self.rigid_parts)
        self.nparts = sum(self.nummat) + self.nrigid
        self.titles = self._titles
        self.signature = (
            n,
            [tuple(map(tuple, r)) for r in self.rows],
            self.extra,
            self.extra20,
            self.extra27,
            self.extra8,
            self.rigid_parts,
        )

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


class _Layout:
    """The columns of a family's element variables: the ``layer`` slots
    ``(name, components)`` repeat for each of ``layers`` layers, the ``tail``
    slots ``(name, components, count)`` follow once."""

    def __init__(self):
        self.layers = 1
        self.layer = []
        self.tail = []

    def layer_width(self):
        return sum(c for _, c in self.layer)

    def width(self):
        return self.layers * self.layer_width() + sum(c * n for _, c, n in self.tail)


class _Plan:
    """The header: what the states of all steps need."""

    def __init__(self, model, meshes, double=False):
        self.model = model
        self.words_ = _Words(double)
        # what the arrays of all steps hold: family -> name -> [points, comps]
        seen = [{}, {}, {}, {}]
        node = {}
        self.alive = self.alive_nodes = False
        self.globals_ = self.part_globals = False
        self.part_count = 0
        nparts = model.nparts
        for mesh in meshes:
            for name in _NODAL_NAMES:
                if name in mesh.point_data:
                    a = np.asarray(mesh.point_data[name])
                    node[name] = max(node.get(name, 0), a.size // max(model.n, 1))
            self.alive = self.alive or "lsdyna:alive" in mesh.cell_data
            self.alive_nodes = self.alive_nodes or "lsdyna:alive" in mesh.point_data
            self.globals_ = self.globals_ or any(
                k in mesh.field_data for k, _ in _GLOBALS
            )
            # the per-part arrays may cover more parts than the mesh uses
            count = np.asarray(mesh.field_data.get(_PART_GLOBALS[0][0], [])).size
            if (
                count >= nparts
                and not self.part_globals
                and all(
                    k in mesh.field_data
                    and np.asarray(mesh.field_data[k]).size == w * count
                    for k, w in _PART_GLOBALS
                )
            ):
                self.part_globals = True
                self.part_count = count
            for name in _COMPS:
                arrays = _cell_array(mesh, name)
                if arrays is None:
                    continue
                for family in (_SOLID, _BEAM, _SHELL, _SPH):
                    for b, _, _ in model.ranges[family]:
                        a = arrays[b]
                        if a is None:
                            continue
                        points, comps = a.shape[1], a.shape[2]
                        if (
                            family == _SOLID
                            and name != "thermal_variables"
                            and points not in (1, 8)
                        ):
                            continue
                        entry = seen[family].setdefault(name, [0, 0])
                        entry[0] = max(entry[0], points)
                        entry[1] = max(entry[1], comps)
        self.node = node
        self.seen = seen
        self.rigid = {}  # reader name -> components, over all steps
        if model.nrigid:
            for mesh in meshes:
                for name, width in _RIGID_FIELDS:
                    value = mesh.field_data.get("lsdyna:rigid_body_" + name)
                    if (
                        value is not None
                        and np.asarray(value).size == model.nrigid * width
                    ):
                        self.rigid[name] = width
        self.globals_ = self.globals_ or self.part_globals
        # the header's material counts add up to the parts the arrays cover
        self.nummat = list(model.nummat)
        if self.part_globals:
            self.nummat[0] += self.part_count - nparts
        self.nglbv = (
            (6 + (7 * self.part_count if self.part_globals else 0))
            if self.globals_
            else 0
        )
        self._element_layouts()
        self._nodal_flags()
        # a state of the time alone cannot be told from the end mark
        self.words = self._state_words()
        if self.words < 2:
            self.disp = True
            self.words = self._state_words()

    # -- element variables -------------------------------------------------

    def _has(self, family, *names):
        return any(n in self.seen[family] for n in names)

    def _element_layouts(self):
        m, seen = self.model, self.seen
        solid_seen, shell_seen, beam_seen = seen[_SOLID], seen[_SHELL], seen[_BEAM]
        # the header holds one flag of each for solids and shells together
        self.stress = self._has(_SOLID, "stress") or self._has(_SHELL, "stress")
        both = (_SOLID, _SHELL)
        self.pstrain = any(self._has(f, "effective_plastic_strain") for f in both)
        self.plastic_tensor = any(self._has(f, "plastic_strain_tensor") for f in both)
        self.thermal_tensor = any(self._has(f, "thermal_strain_tensor") for f in both)
        self.strain = self._has(_SOLID, "strain") or self._has(
            _SHELL, "strain_inner", "strain_outer"
        )
        self.forces = self._has(_SHELL, *[n for n, _ in _SHELL_FORCES])
        self.extra = self._has(
            _SHELL, "thickness", "shell_element_variables", "internal_energy"
        )

        solid = _Layout()
        if m.counts[_SOLID]:
            solid.layers = 8 if any(v[0] == 8 for v in solid_seen.values()) else 1
            hist = solid_seen.get("history_variables", [0, 0])
            if self.stress:
                solid.layer.append(("stress", 6))
            if self.pstrain:
                solid.layer.append(("effective_plastic_strain", 1))
            if self.plastic_tensor:
                solid.layer.append(("plastic_strain_tensor", 6))
            if self.thermal_tensor:
                solid.layer.append(("thermal_strain_tensor", 6))
            if hist[1]:
                solid.layer.append(("history_variables", hist[1]))
            if self.strain:
                solid.layer.append(("strain", 6))
        self.solid = solid
        self.neiph = sum(
            c for n, c in solid.layer if n not in ("stress", "effective_plastic_strain")
        )
        self.nt3d = solid_seen.get("thermal_variables", [0, 0])[1]

        shell = _Layout()
        self.neips = 0
        if m.counts[_SHELL]:
            layer_points = [
                v[0]
                for k, v in shell_seen.items()
                if k in ("stress", "effective_plastic_strain", "history_variables")
            ]
            shell.layers = max(layer_points) if layer_points else 1
            hist = shell_seen.get("history_variables", [0, 0])
            if self.stress:
                shell.layer.append(("stress", 6))
            if self.pstrain:
                shell.layer.append(("effective_plastic_strain", 1))
            if hist[1]:
                shell.layer.append(("history_variables", hist[1]))
            self.neips = hist[1]
            if self.forces:
                shell.tail += [(n, c, 1) for n, c in _SHELL_FORCES]
            if self.extra:
                shell.tail += [(n, c, 1) for n, c in _SHELL_EXTRA]
            if self.strain:
                shell.tail += [("strain_inner", 6, 1), ("strain_outer", 6, 1)]
            if self.extra:
                shell.tail.append(("internal_energy", 1, 1))
            if self.plastic_tensor:
                shell.tail.append(
                    (
                        "plastic_strain_tensor",
                        6,
                        shell_seen.get("plastic_strain_tensor", [shell.layers, 6])[0],
                    )
                )
            if self.thermal_tensor:
                shell.tail.append(("thermal_strain_tensor", 6, 1))
        self.shell = shell

        # a beam: six resultants, five values per layer, (3 + layers) * NEIPB history
        self.beam = False
        self.beam_layers = 0
        self.neipb = 0
        if m.counts[_BEAM]:
            layer_points = [beam_seen[n][0] for n, _ in _BEAM_LAYER if n in beam_seen]
            hist = beam_seen.get("history_variables", [0, 0])
            core = self._has(_BEAM, *[n for n, _ in _BEAM_CORE])
            if hist[1]:
                layer_points.append(max(hist[0] - 3, 0))
            if core or layer_points or hist[1]:
                self.beam = True
                self.beam_layers = max(layer_points) if layer_points else 0
                self.neipb = hist[1]
        # SPH particles: the flag block and one row of values per particle
        self.sph_flags = [0] * 11
        self.nsph_vars = 0
        if m.counts[_SPH]:
            sph_seen = seen[_SPH]
            flags = [11] + [0] * 10
            for name, index, comps in _SPH_VARS:
                if name in sph_seen:
                    flags[index] = comps
            if "strain_rate" in sph_seen:
                flags[8] = 12
            elif "strain" in sph_seen:
                flags[8] = 6
            flags[9] = 1 if "mass" in sph_seen else 0
            flags[10] = sph_seen.get("history_variables", [0, 0])[1]
            self.sph_flags = flags
            self.nsph_vars = 1 + sum(flags[1:8]) + abs(flags[8]) + flags[9] + flags[10]
        self.nv3d = solid.width() if m.counts[_SOLID] else 0
        self.nv2d = shell.width() if m.counts[_SHELL] else 0
        self.nv1d = (
            6 + 5 * self.beam_layers + (3 + self.beam_layers) * self.neipb
            if self.beam
            else 0
        )

    # -- nodal variables ---------------------------------------------------

    def _nodal_flags(self):
        node = self.node
        self.disp = "displacement" in node
        self.vel = "velocity" in node
        self.acc = "acceleration" in node
        # IT: temperature (1), temperature and heat flux (2), through the thickness (3)
        self.it0 = 0
        if node.get("temperature", 0) > 1:
            self.it0 = 3
        elif "heat_flux" in node:
            self.it0 = 2
        elif "temperature" in node:
            self.it0 = 1
        self.it1 = 1 if "mass_scaling" in node else 0
        self.tgrad = "temperature_gradient" in node
        self.resid = "residual_forces" in node or "residual_moments" in node

    def node_width(self):
        w = 3 * self.disp
        if self.it0:
            w += 1 if self.it0 == 1 else 4
            if self.it0 == 3:
                w += 2
        if self.it1:
            w += 1
        if self.tgrad:
            w += 1
        if self.resid:
            w += 6
        return w + 3 * self.vel + 3 * self.acc

    def _state_words(self):
        m = self.model
        n = 1 + self.nglbv + m.n * self.node_width()
        n += m.counts[_SOLID] * (self.nv3d + self.nt3d)
        n += m.counts[_BEAM] * self.nv1d + m.counts[_SHELL] * self.nv2d
        if self.alive:
            n += sum(m.counts[:3])
        elif self.alive_nodes:
            n += m.n
        n += m.counts[_SPH] * self.nsph_vars + m.nrigid * 24
        return n

    def maxint(self):
        layers = max(self.shell.layers if self.model.counts[_SHELL] else 1, 1)
        if self.alive:
            return -(10000 + layers)
        if self.alive_nodes:
            return -layers
        return layers

    def idtdt(self):
        return (
            100000
            + 10000 * self.strain
            + 1000 * self.thermal_tensor
            + 100 * self.plastic_tensor
            + 10 * self.resid
            + self.tgrad
        )

    def it(self):
        return 10 * self.it1 + self.it0

    def extra_words(self):
        m = self.model
        return 16 if (self.nt3d or self.neipb or m.extra20 or m.extra27) else 0


def _geometry(plan):
    m, w = plan.model, plan.words_
    words = [0] * (64 + plan.extra_words())
    fields = {
        11: 1,
        15: 8 if m.nrigid else 4,
        16: m.n,
        17: 6,
        18: plan.nglbv,
        19: plan.it(),
        20: int(plan.disp),
        21: int(plan.vel),
        22: int(plan.acc),
        23: -m.counts[_SOLID] if m.extra else m.counts[_SOLID],
        24: plan.nummat[_SOLID],
        27: plan.nv3d,
        28: m.counts[_BEAM],
        29: m.nummat[_BEAM],
        30: plan.nv1d,
        31: m.counts[_SHELL],
        32: m.nummat[_SHELL],
        33: plan.nv2d,
        34: plan.neiph,
        35: plan.neips,
        36: plan.maxint(),
        37: m.counts[_SPH],
        38: 1 if m.counts[_SPH] else 0,
        43: 1000 if plan.stress else 0,
        44: 1000 if plan.pstrain else 0,
        45: 1000 if plan.forces else 0,
        46: 1000 if plan.extra else 0,
        51: len(m.part_ids),
        55: len(m.extra8),
        56: plan.idtdt(),
    }
    if plan.extra_words():
        fields[57] = plan.extra_words()
        fields[64] = len(m.extra20)
        fields[65] = plan.nt3d
        fields[66] = len(m.extra27)
        fields[67] = plan.neipb
    ids = m.counts[_SOLID] + m.counts[_BEAM] + m.counts[_SHELL]
    fields[39] = 16 + m.n + ids + 3 * len(m.part_ids)
    for i, v in fields.items():
        words[i] = v
    out = bytearray(w.i(words))
    title = _TITLE + b" " * (10 * w.ws - len(_TITLE))
    out[: 10 * w.ws] = title
    out[14 * w.ws : 15 * w.ws] = w.f([_VERSION])
    if m.counts[_SPH]:
        out += w.i(plan.sph_flags)
    out += w.f(m.points.reshape(-1))
    solids = np.zeros((m.counts[_SOLID], 9), dtype=np.int64)
    for e, nodes in enumerate(m.rows[_SOLID]):
        solids[e, :8] = np.asarray(nodes) + 1
        solids[e, 8] = m.part_index[m.parts_of[_SOLID][e]]
    out += w.i(solids)
    if m.extra:
        extra = np.zeros((m.counts[_SOLID], 2), dtype=np.int64)
        for e, nodes in m.extra:
            extra[e] = np.asarray(nodes) + 1
        out += w.i(extra)
    beams = np.zeros((m.counts[_BEAM], 6), dtype=np.int64)
    for e, nodes in enumerate(m.rows[_BEAM]):
        beams[e, :2] = np.asarray(nodes[:2]) + 1
        beams[e, 5] = m.part_index[m.parts_of[_BEAM][e]]
    out += w.i(beams)
    shells = np.zeros((m.counts[_SHELL], 5), dtype=np.int64)
    for e, nodes in enumerate(m.rows[_SHELL]):
        shells[e, :4] = np.asarray(nodes) + 1
        shells[e, 4] = m.part_index[m.parts_of[_SHELL][e]]
    out += w.i(shells)
    # user ids: the header (NSORT < 0), then nodes, solids, beams, shells,
    # thick shells (none) and the parts
    nparts = len(m.part_ids)
    solids_at = 17 + m.n
    beams_at = solids_at + m.counts[_SOLID]
    shells_at = beams_at + m.counts[_BEAM]
    parts_at = shells_at + m.counts[_SHELL]
    out += w.i(
        [-17, solids_at, beams_at, shells_at, parts_at]
        + [m.n, m.counts[_SOLID], m.counts[_BEAM], m.counts[_SHELL], 0]
    )
    out += w.i(
        [parts_at, parts_at + nparts, parts_at + 2 * nparts, nparts, m.nrigid, nparts]
    )
    out += w.i(m.nid)
    for family in (_SOLID, _BEAM, _SHELL):
        out += w.i(m.cell_eid[family])
    # the parts' ids, then (as LS-DYNA writes them) in defined order and cross-referenced
    out += w.i(m.part_ids) * 2 + w.i(range(1, nparts + 1))
    # rigid bodies (part, no node lists), then the SPH particles' (node, part)
    if m.nrigid:
        out += w.i([m.nrigid])
        for part in m.rigid_parts:
            out += w.i([m.part_index[part], 0, 0])
    for e, nodes in enumerate(m.rows[_SPH]):
        out += w.i([nodes[0] + 1, m.part_index[m.parts_of[_SPH][e]]])
    # the higher-order tables: eight-node shells, then 20- and 27-node solids
    for e, nodes in m.extra8:
        out += w.i([e + 1] + [v + 1 for v in nodes])
    for e, nodes in m.extra20:
        out += w.i([e + 1] + [v + 1 for v in nodes])
    for e, nodes in m.extra27:
        out += w.i([e + 1] + [v + 1 for v in nodes])
    out += w.f([_END_MARK])
    titles = bytearray()
    for pid in m.part_ids:
        text = m.titles.get(pid, f"Part {pid}").encode("latin-1", "replace")[:72]
        titles += w.i([pid]) + text + b" " * (72 - len(text))
    out += w.i([90001, len(m.part_ids)]) + bytes(titles)
    out += w.f([_END_MARK])
    return bytes(out)


def _family_matrix(plan, mesh, family, layout):
    """A family's element variables as ``(cells, width)``, NaN where absent."""
    m = plan.model
    out = np.full((m.counts[family], layout.width()), np.nan)
    layer_width = layout.layer_width()
    offset = 0
    for name, comps in layout.layer:
        arrays = _cell_array(mesh, name)
        if arrays is not None:
            for b, first, n in m.ranges[family]:
                a = arrays[b]
                if a is None or n == 0:
                    continue
                k = min(a.shape[2], comps)
                for layer in range(min(a.shape[1], layout.layers)):
                    col = layer * layer_width + offset
                    out[first : first + n, col : col + k] = a[:, layer, :k]
        offset += comps
    base = layout.layers * layer_width
    for name, comps, count in layout.tail:
        arrays = _cell_array(mesh, name)
        if arrays is not None:
            for b, first, n in m.ranges[family]:
                a = arrays[b]
                if a is None or n == 0:
                    continue
                k = min(a.shape[2], comps)
                for point in range(min(a.shape[1], count)):
                    col = base + point * comps
                    out[first : first + n, col : col + k] = a[:, point, :k]
        base += comps * count
    return out


def _beam_matrix(plan, mesh):
    """The beams' columns: six resultants, then per layer the stress, shear
    stress (2), plastic strain and axial strain, then the history values."""
    m = plan.model
    nl, nh = plan.beam_layers, plan.neipb
    out = np.full((m.counts[_BEAM], plan.nv1d), np.nan)
    at = 0
    for name, comps in _BEAM_CORE:
        arrays = _cell_array(mesh, name)
        if arrays is not None:
            for b, first, n in m.ranges[_BEAM]:
                a = arrays[b]
                if a is not None and n:
                    k = min(a.shape[2], comps)
                    out[first : first + n, at : at + k] = a[:, 0, :k]
        at += comps
    offset = 0
    for name, comps in _BEAM_LAYER:
        arrays = _cell_array(mesh, name)
        if arrays is not None:
            for b, first, n in m.ranges[_BEAM]:
                a = arrays[b]
                if a is None or n == 0:
                    continue
                k = min(a.shape[2], comps)
                for layer in range(min(a.shape[1], nl)):
                    col = 6 + layer * 5 + offset
                    out[first : first + n, col : col + k] = a[:, layer, :k]
        offset += comps
    if nh:
        arrays = _cell_array(mesh, "history_variables")
        if arrays is not None:
            for b, first, n in m.ranges[_BEAM]:
                a = arrays[b]
                if a is None or n == 0:
                    continue
                k = min(a.shape[2], nh)
                for point in range(min(a.shape[1], 3 + nl)):
                    col = 6 + 5 * nl + point * nh
                    out[first : first + n, col : col + k] = a[:, point, :k]
    return out


def _sph_matrix(plan, mesh):
    """The SPH particles' rows: the material (negative once deleted), then the
    variables the flag block names."""
    m = plan.model
    flags = plan.sph_flags
    out = np.full((m.counts[_SPH], plan.nsph_vars), np.nan)
    out[:, 0] = [m.part_index[p] for p in m.parts_of[_SPH]]
    alive = mesh.cell_data.get("lsdyna:alive")
    if alive is not None:
        for b, first, count in m.ranges[_SPH]:
            if b < len(alive) and count:
                dead = np.asarray(alive[b]).reshape(-1) == 0
                out[first : first + count, 0] *= np.where(dead, -1.0, 1.0)
    strain = abs(flags[8])
    columns = [
        ("sph_radius", flags[1]),
        ("sph_pressure", flags[2]),
        ("stress", flags[3]),
        ("effective_plastic_strain", flags[4]),
        ("density", flags[5]),
        ("internal_energy", flags[6]),
        ("sph_neighbors", flags[7]),
        ("strain", min(strain, 6)),
        ("strain_rate", strain - min(strain, 6)),
        ("mass", flags[9]),
        ("history_variables", flags[10]),
    ]
    at = 1
    for name, width in columns:
        if width > 0:
            arrays = _cell_array(mesh, name)
            if arrays is not None:
                for b, first, count in m.ranges[_SPH]:
                    a = arrays[b]
                    if a is None or count == 0:
                        continue
                    k = min(a.shape[2], width)
                    out[first : first + count, at : at + k] = a[:, 0, :k]
        at += max(width, 0)
    return out


def _state(plan, mesh, time):
    m, w = plan.model, plan.words_
    n = m.n
    pieces = [w.f([time])]
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
            nparts = plan.part_count
            for name, width in _PART_GLOBALS:
                v = mesh.field_data.get(name)
                if v is not None and np.asarray(v).size == width * nparts:
                    g[at : at + width * nparts] = np.asarray(
                        v, dtype=np.float64
                    ).reshape(-1)
                at += width * nparts
        pieces.append(w.f(g))
    # nodal variables, in the reader's order
    if plan.disp:
        pieces.append(
            w.f(m.points + np.nan_to_num(_node_array(mesh, "displacement", n)))
        )
    if plan.it0:
        layers = plan.it0 == 3
        pieces.append(w.f(_node_array(mesh, "temperature", n, 3 if layers else 1)))
        if plan.it0 in (2, 3):
            pieces.append(w.f(_node_array(mesh, "heat_flux", n)))
    if plan.it1:
        pieces.append(w.f(_node_array(mesh, "mass_scaling", n, 1)))
    if plan.tgrad:
        pieces.append(w.f(_node_array(mesh, "temperature_gradient", n, 1)))
    if plan.resid:
        pieces.append(w.f(_node_array(mesh, "residual_forces", n)))
        pieces.append(w.f(_node_array(mesh, "residual_moments", n)))
    if plan.vel:
        pieces.append(w.f(_node_array(mesh, "velocity", n)))
    if plan.acc:
        pieces.append(w.f(_node_array(mesh, "acceleration", n)))
    if plan.nt3d:
        thermal = np.full((m.counts[_SOLID], plan.nt3d), np.nan)
        arrays = _cell_array(mesh, "thermal_variables")
        if arrays is not None:
            for b, first, count in m.ranges[_SOLID]:
                a = arrays[b]
                if a is not None and count:
                    k = min(a.shape[2], plan.nt3d)
                    thermal[first : first + count, :k] = a[:, 0, :k]
        pieces.append(w.f(thermal))
    if plan.nv3d:
        pieces.append(w.f(_family_matrix(plan, mesh, _SOLID, plan.solid)))
    if plan.nv1d:
        pieces.append(w.f(_beam_matrix(plan, mesh)))
    if plan.nv2d:
        pieces.append(w.f(_family_matrix(plan, mesh, _SHELL, plan.shell)))
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
            pieces.append(w.f(flags))
    elif plan.alive_nodes:
        value = mesh.point_data.get("lsdyna:alive")
        flags = np.ones(n) if value is None else (np.asarray(value).reshape(-1) != 0)
        pieces.append(w.f(flags))
    if m.counts[_SPH]:
        pieces.append(w.f(_sph_matrix(plan, mesh)))
    if m.nrigid:
        rigid = np.full((m.nrigid, 24), np.nan)
        at = 0
        for name, width in _RIGID_FIELDS:
            value = mesh.field_data.get("lsdyna:rigid_body_" + name)
            if value is not None and np.asarray(value).size == m.nrigid * width:
                rigid[:, at : at + width] = np.asarray(value, dtype=np.float64).reshape(
                    m.nrigid, width
                )
            at += width
        pieces.append(w.f(rigid))
    return b"".join(pieces)


def _assemble(plan, steps, max_bytes):
    """The files of the family: the geometry in the first, the end mark last."""
    w = plan.words_
    head = _geometry(plan)
    mark = w.f([_END_MARK])
    state_bytes = plan.words * w.ws
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


def write(filename, mesh, max_file_bytes=_DEFAULT_FILE_BYTES, double=False):
    """Write one state as a d3plot family (a single file for one state)."""
    if is_buffer(filename, "w"):
        raise WriteError("LS-DYNA d3plot: cannot be written to a stream")
    model = _Model(mesh)
    plan = _Plan(model, [mesh], double)
    _store(
        os.fspath(filename), _assemble(plan, [(_time_of(mesh), mesh)], max_file_bytes)
    )


class SeriesWriter:
    """A d3plot family with a state per step; the geometry is the first step's,
    which every later step must share. The header takes the union of the
    variables of all steps, so the file is written when the writer closes."""

    def __init__(self, filename, max_file_bytes=_DEFAULT_FILE_BYTES, double=False):
        if is_buffer(filename, "w"):
            raise WriteError("LS-DYNA d3plot: cannot be written to a stream")
        self._path = os.fspath(filename)
        self._max = int(max_file_bytes)
        self._double = bool(double)
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
        plan = _Plan(self._model, [mesh for _, mesh in self._steps], self._double)
        _store(self._path, _assemble(plan, self._steps, self._max))

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
