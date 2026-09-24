"""Freeze pymapdl-reader's reading of the Ansys ``.rst``/``.rth`` fixtures.

Writes ``tests/python/meshes/ansys/rst/pymapdl_reference.npz``: for each fixture,
the result set times, and for every set the node numbers and the nodal DOF
solution pymapdl-reader returns (with the DOF labels), plus the cells of its VTK
grid (points, cell types, offsets, connectivity). Where the file has them, each
set also gets pymapdl-reader's node-averaged stresses (``nodal_stress``) and
elastic strains (``nodal_elastic_strain``, first six components) and its
reaction forces (``nodal_reaction_forces``, DOF labels). A distributed solve
(``dist_static/file0.rst``) is read as pymapdl-reader's ``DistributedResult``;
a static cyclic model also gets its full-rotor points, displacements and
stresses (``full_rotor=True``). ``test_ansys_rst.py`` compares
meshio++'s arrays against it, so the cross-check runs without pymapdl-reader
installed.

pymapdl-reader (MIT) is not a dependency. Run this in a throw-away environment:

    uv venv --python 3.12 /tmp/pymapdl && uv pip install --python /tmp/pymapdl/bin/python \
        "ansys-mapdl-reader @ git+https://github.com/ansys/pymapdl-reader@488ad8b" \
        pyvista==0.49.0
    cd /tmp && /tmp/pymapdl/bin/python /path/to/meshioplusplus/tools/gen_ansys_rst_reference.py

(The frozen file came from that commit, 0.56.dev0. Run it outside a source
checkout of pymapdl-reader, whose package directory would shadow the install.)

pymapdl-reader corrupts its heap reading the element results of SHELL181 and
SHELL281 models (it sizes their records by the doubled node count), so those are
not frozen. pymapdl-reader sizes a record's values by its length in 4-byte words, so on a
result set that holds only some nodes it reads past the record; none of these
fixtures has such a set (see doc/formats/ansys_rst.md).
"""

import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
FIXTURES = HERE.parent / "tests" / "python" / "meshes" / "ansys" / "rst"


def _element_results(rst, key, s, out):
    available = rst.available_results
    routines = set(int(v) for v in np.asarray(rst.mesh.ekey)[:, 1])
    if routines & {181, 281}:
        return
    if available["ENS"]:
        nnum, values = rst.nodal_stress(s)
        out[f"{key}/{s}/stress_nnum"] = np.asarray(nnum, dtype=np.int64)
        out[f"{key}/{s}/stress"] = np.asarray(values, dtype=np.float64)
    if available["EEL"]:
        nnum, values = rst.nodal_elastic_strain(s)
        out[f"{key}/{s}/strain_nnum"] = np.asarray(nnum, dtype=np.int64)
        out[f"{key}/{s}/strain"] = np.asarray(values, dtype=np.float64)[:, :6]
    if rst._solution_header(s)["nrf"]:
        forces, nnum, dof = rst.nodal_reaction_forces(s)
        labels = rst.result_dof(s)
        out[f"{key}/{s}/rf_nnum"] = np.asarray(nnum, dtype=np.int64)
        out[f"{key}/{s}/rf_values"] = np.asarray(forces, dtype=np.float64)
        out[f"{key}/{s}/rf_dofs"] = np.array([labels[d - 1] for d in dof], dtype=str)


def _cyclic(rst, key, s, out):
    """The full rotor of a static cyclic model, as pymapdl-reader expands it."""
    out[f"{key}/rotor_points"] = np.asarray(rst.full_rotor.points, dtype=np.float64)
    nnum, values = rst.nodal_solution(s, full_rotor=True)
    out[f"{key}/{s}/rotor_nnum"] = np.asarray(nnum, dtype=np.int64)
    out[f"{key}/{s}/rotor_values"] = np.asarray(values, dtype=np.float64)
    _, values = rst.nodal_stress(s, full_rotor=True)
    out[f"{key}/{s}/rotor_stress"] = np.asarray(values, dtype=np.float64)


def main():
    from ansys.mapdl import reader as pymapdl_reader

    out = {}
    paths = sorted(FIXTURES.glob("*.rst")) + sorted(FIXTURES.glob("*.rth"))
    paths.append(FIXTURES / "dist_static" / "file0.rst")
    paths.append(FIXTURES / "dist_static" / "file.rst")
    for path in paths:
        rst = pymapdl_reader.read_binary(str(path))
        key = str(path.relative_to(FIXTURES))
        out[f"{key}/times"] = np.asarray(rst.time_values, dtype=np.float64)
        cyclic = rst.n_sector > 1 and rst._resultheader["kan"] == 0
        for s in range(rst.n_results):
            nnum, values = rst.nodal_solution(s)
            out[f"{key}/{s}/nnum"] = np.asarray(nnum, dtype=np.int64)
            out[f"{key}/{s}/values"] = np.asarray(values, dtype=np.float64)
            out[f"{key}/{s}/dofs"] = np.array(rst.result_dof(s), dtype=str)
            _element_results(rst, key, s, out)
            if cyclic:
                _cyclic(rst, key, s, out)
        grid = rst.quadgrid
        out[f"{key}/points"] = np.asarray(grid.points, dtype=np.float64)
        out[f"{key}/celltypes"] = np.asarray(grid.celltypes, dtype=np.uint8)
        out[f"{key}/offsets"] = np.asarray(grid.cell_offsets, dtype=np.int64)
        out[f"{key}/connectivity"] = np.asarray(grid.cell_connectivity, dtype=np.int64)
    target = FIXTURES / "pymapdl_reference.npz"
    np.savez_compressed(target, **out)
    print(f"wrote {len(paths)} fixtures to {target}")


if __name__ == "__main__":
    main()
