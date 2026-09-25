"""Write the generated LS-DYNA d3plot families and freeze lasso-python's reading.

The d3plot fixtures under ``tests/python/meshes/lsdyna_d3plot/`` are of two kinds:

* ``lasso/<name>/``: lasso-python's own test families (real LS-DYNA output, BSD-3),
  copied as they are;
* ``generated/<name>/``: families this script writes with lasso-python's
  ``D3plot.write_d3plot`` -- a shell + solid + beam model whose shells and solids
  yield and are then deleted, with its states in ``d3plot01`` (``shell_solid``),
  one per file (``shell_solid_split``) and in double precision
  (``shell_solid_double``); and a 20-node and a 27-node hexahedron with a
  moving rigid body (``quadratic_rigid``), for which lasso-python's writer
  needs the patches below;
* ``dyna/<name>/``: LS-DYNA families from Ansys' example data (MIT), trimmed
  to a few of their files: ``bird_strike`` (SPH particles and composite
  shells; the base file and the first state of ``d3plot01``) and
  ``projectile`` (element erosion; its ``d3plot``, ``d3plot03`` and
  ``d3plot16`` as ``d3plot``, ``d3plot01``, ``d3plot02``).

For every family the script reads each state back with lasso-python and stores
what it reads in ``lasso_reference.npz``: the tests compare meshio++ against that
file, so they need neither lasso-python nor LS-DYNA. Run it in a throwaway
environment with ``pip install lasso-python`` (tested with 2.0.4)::

    python tools/gen_d3plot_reference.py <lasso-python checkout>/test/test_data \
        [<ansys example-data checkout>/result_files]
"""

import os
import shutil
import sys

import lasso.dyna.d3plot as _lasso_d3plot
import numpy as np
from lasso.dyna import ArrayType as A
from lasso.dyna import D3plot

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "tests", "python", "meshes", "lsdyna_d3plot")

LASSO_FAMILIES = [
    "simple_d3plot",
    "d3plot_solid_int",
    "d3plot_beamip",
    "d3plot_node_temperature",
    "order_d3plot",
]


def _model(nstates=4):
    """A plate of shells over a row of solids, with two beams.

    Returns the arrays of a D3plot: nodes, one hexahedron, one wedge and one
    tetrahedron (degenerate 8-node solids), five quads and a triangle
    (degenerate quad), two beams, and ``nstates`` states in which the plate
    bends, plastic strain grows and elements are deleted one by one.
    """
    # solid block nodes 0..11 (a 2x1x1 grid), plate nodes 12..23 above it
    xs = np.array([0.0, 1.0, 2.0])
    solid_nodes = [(x, y, z) for z in (0.0, 1.0) for y in (0.0, 1.0) for x in xs]
    plate_nodes = [
        (x, y, 2.0) for y in (0.0, 1.0) for x in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5)
    ]
    coords = np.array(solid_nodes + plate_nodes, dtype=np.float32)

    def s(i, j, k):  # solid grid node index
        return k * 6 + j * 3 + i

    hexa = [
        s(0, 0, 0),
        s(1, 0, 0),
        s(1, 1, 0),
        s(0, 1, 0),
        s(0, 0, 1),
        s(1, 0, 1),
        s(1, 1, 1),
        s(0, 1, 1),
    ]
    # wedge written LS-DYNA style: n1 n2 n3 n3 n5 n6 n7 n7 -> collapsed edges
    wedge = [
        s(1, 0, 0),
        s(2, 0, 0),
        s(2, 1, 0),
        s(2, 1, 0),
        s(1, 0, 1),
        s(2, 0, 1),
        s(2, 1, 1),
        s(2, 1, 1),
    ]
    # tetrahedron: n1 n2 n3 n4 n4 n4 n4 n4
    tetra = [s(1, 0, 1), s(2, 0, 1), s(2, 1, 1), s(1, 1, 1)] + [s(1, 1, 1)] * 4
    solids = np.array([hexa, wedge, tetra], dtype=np.int32)

    def p(i, j):  # plate node index
        return 12 + j * 6 + i

    quads = [[p(i, 0), p(i + 1, 0), p(i + 1, 1), p(i, 1)] for i in range(5)]
    tria = [p(4, 0), p(5, 0), p(5, 1), p(5, 1)]  # n4 == n3 -> triangle
    shells = np.array(quads[:4] + [tria] + quads[4:], dtype=np.int32)
    beams = np.array(
        [[s(0, 0, 1), p(0, 0), 0, 0, 0], [s(2, 0, 1), p(4, 0), 0, 0, 0]], dtype=np.int32
    )

    nn, ns, nsh, nb = len(coords), len(solids), len(shells), len(beams)
    t = np.linspace(0.0, 3.0e-3, nstates).astype(np.float32)
    arrays = {
        A.node_coordinates: coords,
        A.node_ids: np.arange(101, 101 + nn, dtype=np.int32),
        A.element_solid_node_indexes: solids,
        A.element_solid_part_indexes: np.zeros(ns, dtype=np.int32),
        A.element_solid_ids: np.array([11, 12, 13], dtype=np.int32),
        A.element_shell_node_indexes: shells,
        A.element_shell_part_indexes: np.ones(nsh, dtype=np.int32),
        A.element_shell_ids: np.arange(21, 21 + nsh, dtype=np.int32),
        A.element_beam_node_indexes: beams,
        A.element_beam_part_indexes: np.full(nb, 2, dtype=np.int32),
        A.element_beam_ids: np.array([31, 32], dtype=np.int32),
        A.part_ids: np.array([100, 200, 300], dtype=np.int32),
        A.part_titles_ids: np.array([100, 200, 300], dtype=np.int32),
        A.part_titles: np.array(
            [b"block".ljust(72), b"plate".ljust(72), b"struts".ljust(72)]
        ),
        A.global_timesteps: t,
    }
    rng = np.random.default_rng(7)
    disp = np.zeros((nstates, nn, 3), dtype=np.float32)
    for k in range(nstates):
        bend = 0.05 * k * np.sin(np.pi * coords[:, 0] / 2.5)
        disp[k] = coords
        disp[k, :, 2] += bend.astype(np.float32)
    arrays[A.node_displacement] = disp
    arrays[A.node_velocity] = (rng.standard_normal((nstates, nn, 3)) * 0.1).astype(
        np.float32
    )
    nl = 3
    arrays[A.element_shell_stress] = rng.standard_normal((nstates, nsh, nl, 6)).astype(
        np.float32
    )
    eps = np.cumsum(np.abs(rng.standard_normal((nstates, nsh, nl))) * 0.01, axis=0)
    eps[0] = 0.0
    arrays[A.element_shell_effective_plastic_strain] = eps.astype(np.float32)
    arrays[A.element_shell_thickness] = np.full((nstates, nsh), 0.1, dtype=np.float32)
    arrays[A.element_shell_internal_energy] = rng.random((nstates, nsh)).astype(
        np.float32
    )
    arrays[A.element_solid_stress] = rng.standard_normal((nstates, ns, 1, 6)).astype(
        np.float32
    )
    seps = np.cumsum(np.abs(rng.standard_normal((nstates, ns, 1))) * 0.02, axis=0)
    seps[0] = 0.0
    arrays[A.element_solid_effective_plastic_strain] = seps.astype(np.float32)
    arrays[A.element_beam_axial_force] = rng.standard_normal((nstates, nb)).astype(
        np.float32
    )
    arrays[A.element_beam_shear_force] = rng.standard_normal((nstates, nb, 2)).astype(
        np.float32
    )
    arrays[A.element_beam_bending_moment] = rng.standard_normal(
        (nstates, nb, 2)
    ).astype(np.float32)
    arrays[A.element_beam_torsion_moment] = rng.standard_normal((nstates, nb)).astype(
        np.float32
    )
    # deletion: shell k dies at state k (the plate tears from its left edge), the
    # tetrahedron at the last state; beams live
    shell_alive = np.ones((nstates, nsh), dtype=np.float32)
    for k in range(1, nstates):
        shell_alive[k:, k - 1] = 0.0
    solid_alive = np.ones((nstates, ns), dtype=np.float32)
    solid_alive[-1, 2] = 0.0
    arrays[A.element_shell_is_alive] = shell_alive
    arrays[A.element_solid_is_alive] = solid_alive
    arrays[A.element_beam_is_alive] = np.ones((nstates, nb), dtype=np.float32)
    arrays[A.global_kinetic_energy] = rng.random(nstates).astype(np.float32)
    arrays[A.global_internal_energy] = rng.random(nstates).astype(np.float32)
    arrays[A.global_total_energy] = rng.random(nstates).astype(np.float32)
    arrays[A.global_velocity] = rng.random((nstates, 3)).astype(np.float32)
    return arrays


def _write(folder, arrays, double=False, single_file=True):
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    plot = D3plot()
    for key, value in arrays.items():
        if double and isinstance(value, np.ndarray):
            if value.dtype == np.float32:
                value = value.astype(np.float64)
            elif value.dtype == np.int32:
                value = value.astype(np.int64)
        plot.arrays[key] = value
    if double:
        plot.header.itype = np.int64
        plot.header.ftype = np.float64
        plot.header.wordsize = 8
    plot.write_d3plot(os.path.join(folder, "d3plot"), single_file=single_file)


# lasso-python 2.0.4's writer fails on these sections (it concatenates 1-D
# arrays along axis 1, reaches the road writer for NDIM 8, and sizes the rigid
# body states wrong): each is written as the database manual lays it out.
def _extra_nodes(self, fp, settings):
    """The 20-node solids' rows (index, 12 extra nodes) and the 27-node ones'
    (index, 19 extra nodes: QUADR = 0), 1-based."""
    n = 0
    for index, extra in (
        (
            A.element_solid_node20_element_index,
            A.element_solid_node20_extra_node_indexes,
        ),
        (
            A.element_solid_node27_element_index,
            A.element_solid_node27_extra_node_indexes,
        ),
    ):
        if index in self.arrays:
            rows = np.column_stack((self.arrays[index] + 1, self.arrays[extra] + 1))
            n += fp.write(settings.pack(rows, dtype_hint=np.integer))
    return n


_road = _lasso_d3plot.D3plot._write_geom_rigid_road_surface


def _no_road(self, fp, settings):
    if A.rigid_road_segment_node_ids not in self.arrays:
        return 0
    return _road(self, fp, settings)


def _rigid_states(self, fp, i_timestep, settings):
    """Per rigid body 24 words: coordinates, rotation matrix, velocity,
    rotational velocity, acceleration, rotational acceleration."""
    if not 8 <= settings.header["ndim"] <= 9:
        return 0
    names = [
        A.rigid_body_coordinates,
        A.rigid_body_rotation_matrix,
        A.rigid_body_velocity,
        A.rigid_body_rot_velocity,
        A.rigid_body_acceleration,
        A.rigid_body_rot_acceleration,
    ]
    rows = np.concatenate([self.arrays[n][i_timestep] for n in names], axis=1)
    return fp.write(settings.pack(rows, dtype_hint=np.floating))


def _quadratic_rigid(nstates=3):
    """A 20-node hexahedron (nodes 0-19) and a 27-node one (20-46) in VTK's node
    order (LS-DYNA's), and six more nodes a rigid body (part 20) holds."""
    rng = np.random.default_rng(3)
    corners = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        float,
    )
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
    edges += [(0, 4), (1, 5), (2, 6), (3, 7)]
    faces = [
        (0, 3, 7, 4),
        (1, 2, 6, 5),
        (0, 1, 5, 4),
        (3, 2, 6, 7),
        (0, 1, 2, 3),
        (4, 5, 6, 7),
    ]
    h20 = np.vstack([corners] + [(corners[a] + corners[b]) / 2 for a, b in edges])
    h27 = np.vstack(
        [h20]
        + [corners[list(f)].mean(0) for f in faces]
        + [corners.mean(0, keepdims=True)]
    ) + [2, 0, 0]
    extra = np.array([[5 + 0.2 * i, 0, 0] for i in range(6)], float)
    coords = np.vstack([h20, h27, extra]).astype(np.float32)
    nn = len(coords)
    arrays = {
        A.node_coordinates: coords,
        A.node_ids: np.arange(1, nn + 1, dtype=np.int32),
        A.element_solid_node_indexes: np.array([range(8), range(20, 28)], np.int32),
        A.element_solid_part_indexes: np.array([0, 0], np.int32),
        A.element_solid_ids: np.array([7, 8], np.int32),
        A.element_solid_node20_element_index: np.array([0], np.int32),
        A.element_solid_node20_extra_node_indexes: np.arange(8, 20, dtype=np.int32)[
            None, :
        ],
        A.element_solid_node27_element_index: np.array([1], np.int32),
        A.element_solid_node27_extra_node_indexes: np.arange(28, 47, dtype=np.int32)[
            None, :
        ],
        A.part_ids: np.array([10, 20], np.int32),
        A.global_timesteps: np.linspace(0, 1e-3, nstates).astype(np.float32),
        A.node_displacement: np.repeat(coords[None], nstates, 0),
        A.element_solid_stress: rng.standard_normal((nstates, 2, 1, 6)).astype(
            np.float32
        ),
        A.element_solid_effective_plastic_strain: rng.random((nstates, 2, 1)).astype(
            np.float32
        ),
        A.rigid_body_part_indexes: np.array([1], np.int32),
        A.rigid_body_n_nodes: np.array([6], np.int32),
        A.rigid_body_node_indexes_list: [np.arange(47, 53, dtype=np.int32)],
        A.rigid_body_n_active_nodes: np.array([2], np.int32),
        A.rigid_body_active_node_indexes_list: [np.array([47, 48], np.int32)],
    }
    for name, width in (
        (A.rigid_body_coordinates, 3),
        (A.rigid_body_rotation_matrix, 9),
        (A.rigid_body_velocity, 3),
        (A.rigid_body_rot_velocity, 3),
        (A.rigid_body_acceleration, 3),
        (A.rigid_body_rot_acceleration, 3),
    ):
        arrays[name] = rng.random((nstates, 1, width)).astype(np.float32)
    return arrays


def _write_patched(folder, arrays):
    saved = (
        _lasso_d3plot.D3plot._write_geom_extra_node_data,
        _lasso_d3plot.D3plot._write_geom_rigid_road_surface,
        _lasso_d3plot.D3plot._write_states_rigid_bodies,
    )
    _lasso_d3plot.D3plot._write_geom_extra_node_data = _extra_nodes
    _lasso_d3plot.D3plot._write_geom_rigid_road_surface = _no_road
    _lasso_d3plot.D3plot._write_states_rigid_bodies = _rigid_states
    try:
        _write(folder, arrays)
    finally:
        (
            _lasso_d3plot.D3plot._write_geom_extra_node_data,
            _lasso_d3plot.D3plot._write_geom_rigid_road_surface,
            _lasso_d3plot.D3plot._write_states_rigid_bodies,
        ) = saved


def copy_dyna(src):
    """The trimmed LS-DYNA families from Ansys' example data (``result_files``)."""
    dyna = os.path.join(OUT, "dyna")
    bird = os.path.join(dyna, "bird_strike")
    projectile = os.path.join(dyna, "projectile")
    for folder in (bird, projectile):
        if os.path.isdir(folder):
            shutil.rmtree(folder)
        os.makedirs(folder)
    base = os.path.join(src, "lsdyna_bird_strike")
    shutil.copyfile(os.path.join(base, "d3plot"), os.path.join(bird, "d3plot"))
    # d3plot01 holds two states: the first is kept
    one_state = D3plot(os.path.join(base, "d3plot"))._compute_n_bytes_per_state()
    data = open(os.path.join(base, "d3plot01"), "rb").read()
    with open(os.path.join(bird, "d3plot01"), "wb") as f:
        f.write(data[:one_state])
    shutil.copyfile(
        os.path.join(base, "EXAMPLE_DATA_LICENSE"), os.path.join(dyna, "LICENSE")
    )
    base = os.path.join(src, "d3plot_projectile")
    # renumbered without gaps, which lasso-python stops at
    for name, target in (
        ("d3plot", "d3plot"),
        ("d3plot03", "d3plot01"),
        ("d3plot16", "d3plot02"),
    ):
        shutil.copyfile(os.path.join(base, name), os.path.join(projectile, target))


def generate():
    gen = os.path.join(OUT, "generated")
    arrays = _model()
    _write(os.path.join(gen, "shell_solid"), arrays)
    _write(os.path.join(gen, "shell_solid_double"), arrays, double=True)
    _write(os.path.join(gen, "shell_solid_split"), arrays, single_file=False)
    _write_patched(os.path.join(gen, "quadratic_rigid"), _quadratic_rigid())


# -- the reference ---------------------------------------------------------------------

_STATE_ARRAYS = {
    "disp": A.node_displacement,
    "vel": A.node_velocity,
    "acc": A.node_acceleration,
    "temp": A.node_temperature,
    "solid_stress": A.element_solid_stress,
    "solid_eps": A.element_solid_effective_plastic_strain,
    "solid_alive": A.element_solid_is_alive,
    "shell_stress": A.element_shell_stress,
    "shell_eps": A.element_shell_effective_plastic_strain,
    "shell_alive": A.element_shell_is_alive,
    "shell_thickness": A.element_shell_thickness,
    "beam_axial_force": A.element_beam_axial_force,
    "beam_eps": A.element_beam_plastic_strain,
    "beam_alive": A.element_beam_is_alive,
    "tshell_stress": A.element_tshell_stress,
    "kinetic_energy": A.global_kinetic_energy,
    "sph_deletion": A.sph_deletion,
    "sph_radius": A.sph_radius,
    "sph_pressure": A.sph_pressure,
    "sph_stress": A.sph_stress,
    "sph_eps": A.sph_effective_plastic_strain,
    "sph_density": A.sph_density,
    "sph_internal_energy": A.sph_internal_energy,
    "sph_neighbors": A.sph_n_neighbors,
    "sph_strain": A.sph_strain,
    "sph_strainrate": A.sph_strainrate,
    "sph_mass": A.sph_mass,
    "rigid_coordinates": A.rigid_body_coordinates,
    "rigid_rotation": A.rigid_body_rotation_matrix,
    "rigid_velocity": A.rigid_body_velocity,
    "rigid_rot_velocity": A.rigid_body_rot_velocity,
    "rigid_acceleration": A.rigid_body_acceleration,
    "rigid_rot_acceleration": A.rigid_body_rot_acceleration,
}
_STATIC_ARRAYS = {
    "coords": A.node_coordinates,
    "node_ids": A.node_ids,
    "solid_ids": A.element_solid_ids,
    "shell_ids": A.element_shell_ids,
    "beam_ids": A.element_beam_ids,
    "times": A.global_timesteps,
    "sph_nodes": A.sph_node_indexes,
    "sph_material": A.sph_node_material_index,
    "rigid_part": A.rigid_body_part_indexes,
    "node20": A.element_solid_node20_extra_node_indexes,
}


def freeze():
    out = {}
    for kind in ("lasso", "generated", "dyna"):
        root = os.path.join(OUT, kind)
        for name in sorted(os.listdir(root)):
            path = os.path.join(root, name, "d3plot")
            if not os.path.isfile(path):
                continue
            plot = D3plot(path)
            key = f"{kind}/{name}"
            for short, array_type in {**_STATIC_ARRAYS, **_STATE_ARRAYS}.items():
                if array_type in plot.arrays:
                    value = plot.arrays[array_type]
                    if isinstance(value, np.ndarray) and value.dtype != object:
                        out[f"{key}:{short}"] = np.asarray(value)
            print(key, plot.n_timesteps, "states")
    np.savez_compressed(os.path.join(OUT, "lasso_reference.npz"), **out)


def main():
    if len(sys.argv) > 1:
        src = sys.argv[1]
        for name in LASSO_FAMILIES:
            dst = os.path.join(OUT, "lasso", name)
            if os.path.isdir(dst):
                shutil.rmtree(dst)
            shutil.copytree(os.path.join(src, name), dst)
    if len(sys.argv) > 2:
        copy_dyna(sys.argv[2])
    generate()
    freeze()


if __name__ == "__main__":
    main()
