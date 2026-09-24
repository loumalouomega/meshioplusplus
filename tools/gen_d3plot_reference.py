"""Write the generated LS-DYNA d3plot families and freeze lasso-python's reading.

The d3plot fixtures under ``tests/python/meshes/lsdyna_d3plot/`` are of two kinds:

* ``lasso/<name>/``: lasso-python's own test families (real LS-DYNA output, BSD-3),
  copied as they are;
* ``generated/<name>/``: families this script writes with lasso-python's
  ``D3plot.write_d3plot`` -- a shell + solid + beam model whose shells and solids
  yield and are then deleted, with its states in ``d3plot01`` (``shell_solid``),
  one per file (``shell_solid_split``) and in double precision
  (``shell_solid_double``).

For every family the script reads each state back with lasso-python and stores
what it reads in ``lasso_reference.npz``: the tests compare meshio++ against that
file, so they need neither lasso-python nor LS-DYNA. Run it in a throwaway
environment with ``pip install lasso-python`` (tested with 2.0.4)::

    python tools/gen_d3plot_reference.py <lasso-python checkout>/test/test_data
"""

import os
import shutil
import sys

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


def generate():
    gen = os.path.join(OUT, "generated")
    arrays = _model()
    _write(os.path.join(gen, "shell_solid"), arrays)
    _write(os.path.join(gen, "shell_solid_double"), arrays, double=True)
    _write(os.path.join(gen, "shell_solid_split"), arrays, single_file=False)


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
}
_STATIC_ARRAYS = {
    "coords": A.node_coordinates,
    "node_ids": A.node_ids,
    "solid_ids": A.element_solid_ids,
    "shell_ids": A.element_shell_ids,
    "beam_ids": A.element_beam_ids,
    "times": A.global_timesteps,
}


def freeze():
    out = {}
    for kind in ("lasso", "generated"):
        root = os.path.join(OUT, kind)
        for name in sorted(os.listdir(root)):
            path = os.path.join(root, name, "d3plot")
            plot = D3plot(path)
            key = f"{kind}/{name}"
            for short, array_type in {**_STATIC_ARRAYS, **_STATE_ARRAYS}.items():
                if array_type in plot.arrays:
                    out[f"{key}:{short}"] = np.asarray(plot.arrays[array_type])
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
    generate()
    freeze()


if __name__ == "__main__":
    main()
