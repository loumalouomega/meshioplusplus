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
stresses (``full_rotor=True``). Every fixture's first set also gets the element
records meshio++ reads as they are written (``element_solution_data``: ENG, EMS,
EMN, EGR, EFX, ENL, EPT, ECT, ESV). Two modal cyclic files are made from
``cyc12.rst`` by ``modal_variant`` (the analysis type set to modal, the harmonic
indices and set times rewritten: one with a mode pair, one without, whose second
half is its duplicate sector) and get pymapdl-reader's full rotor of each set.
``test_ansys_rst.py`` compares meshio++'s arrays against it, so the cross-check
runs without pymapdl-reader installed.

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
import tempfile

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


RAW_RECORDS = ("EMS", "ENG", "EGR", "EFX", "EMN", "ENL", "EPT", "ECT", "ESV")

# The modal variants of cyc12.rst: (name, harmonic index per set, set times).
MODAL_VARIANTS = (
    ("cyc12_modal_pair.rst", (4, -4, 9), (5.0, 5.0, 7.0)),
    ("cyc12_modal_dup.rst", (2, 3, 0), (1.0 / 3.0, 2.0 / 3.0, 1.0)),
)


def modal_variant(source, target, harmonic, times):
    """Write ``source`` (a static cyclic results file) to ``target`` as a modal
    one: analysis type 2, the harmonic index table and the set times rewritten.
    ``test_ansys_rst.py`` holds the same function."""
    data = bytearray(pathlib.Path(source).read_bytes())
    words = np.frombuffer(bytes(data), "<i4")
    header = int(words[0]) + 3  # the result header follows the standard one
    h = words[header + 2 : header + 2 + int(words[header])]

    def pointer(lo, hi):
        return (int(h[lo]) & 0xFFFFFFFF) | (int(h[hi]) << 32)

    def put(word, value):
        raw = value.tobytes()
        data[word * 4 : word * 4 + len(raw)] = raw

    put(header + 2 + 7, np.int32(2))
    cyc = pointer(16, 43)
    for k, v in enumerate(harmonic):
        put(cyc + 2 + k, np.int32(v))
    tim = pointer(11, 41)
    for k, v in enumerate(times):
        put(tim + 2 + 2 * k, np.float64(v))
    pathlib.Path(target).write_bytes(bytes(data))


def _raw_records(rst, key, out):
    for name in RAW_RECORDS:
        try:
            enum, data, _ = rst.element_solution_data(0, name)
        except Exception:  # noqa: BLE001 - pymapdl-reader refuses absent records
            continue
        rows = [
            np.asarray(d, dtype=np.float64) if d is not None else None for d in data
        ]
        if not any(r is not None and len(r) for r in rows):
            continue
        out[f"{key}/raw/{name}/enum"] = np.asarray(enum, dtype=np.int64)
        out[f"{key}/raw/{name}/lengths"] = np.array(
            [0 if r is None else len(r) for r in rows], dtype=np.int64
        )
        out[f"{key}/raw/{name}/values"] = np.concatenate(
            [r for r in rows if r is not None] + [np.zeros(0)]
        )


def _modal(pymapdl_reader, out):
    with tempfile.TemporaryDirectory() as tmp:
        for name, harmonic, times in MODAL_VARIANTS:
            path = pathlib.Path(tmp) / name
            modal_variant(FIXTURES / "cyc12.rst", path, harmonic, times)
            rst = pymapdl_reader.read_binary(str(path))
            for s in range(rst.n_results):
                nnum, values = rst.nodal_solution(s, full_rotor=True)
                out[f"{name}/{s}/rotor_nnum"] = np.asarray(nnum, dtype=np.int64)
                out[f"{name}/{s}/rotor_values"] = np.asarray(values, dtype=np.float64)
                _, values = rst.nodal_stress(s, full_rotor=True)
                out[f"{name}/{s}/rotor_stress"] = np.asarray(values, dtype=np.float64)


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
        _raw_records(rst, key, out)
        grid = rst.quadgrid
        out[f"{key}/points"] = np.asarray(grid.points, dtype=np.float64)
        out[f"{key}/celltypes"] = np.asarray(grid.celltypes, dtype=np.uint8)
        out[f"{key}/offsets"] = np.asarray(grid.cell_offsets, dtype=np.int64)
        out[f"{key}/connectivity"] = np.asarray(grid.cell_connectivity, dtype=np.int64)
    _modal(pymapdl_reader, out)
    target = FIXTURES / "pymapdl_reference.npz"
    np.savez_compressed(target, **out)
    print(f"wrote {len(paths)} fixtures to {target}")


if __name__ == "__main__":
    main()
