"""Z88 structure files (``z88i1.txt``) and results: both engines, the decks
written by ``tools/gen_z88_fixtures.py``, Z88R's own output for the cantilever,
basename dispatch and the writer."""

import pathlib
import shutil

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.z88 import _z88 as py_z88

from .test_node_order import _is_valid

MESHES = pathlib.Path(__file__).parent / "meshes" / "z88"
DECKS = sorted(p.parent for p in MESHES.glob("*/z88i1.txt"))


@pytest.fixture(params=["core", "python"])
def engine(request):
    if request.param == "core":
        return meshioplusplus.z88
    return py_z88


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    assert sorted(a.point_data) == sorted(b.point_data)
    for name in a.point_data:
        np.testing.assert_array_equal(a.point_data[name], b.point_data[name])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_array_equal(x, y)


def _valid(mesh):
    for block in mesh.cells:
        for row in np.asarray(block.data):
            assert _is_valid(np.asarray(mesh.points)[row], block.type), block.type


@pytest.mark.parametrize("deck", DECKS, ids=[d.name for d in DECKS])
def test_engines_agree_on_every_deck(deck):
    _same(meshioplusplus.z88.read(deck / "z88i1.txt"), py_z88.read(deck / "z88i1.txt"))


def test_cantilever_hex20_and_z88r_results(engine):
    mesh = engine.read(MESHES / "cantilever" / "z88i1.txt")
    assert [(c.type, len(c.data)) for c in mesh.cells] == [("hexahedron20", 5)]
    _valid(mesh)
    np.testing.assert_array_equal(mesh.cell_data["z88:type"][0], [10] * 5)
    u = mesh.point_data["U"]
    assert u.shape == (68, 3)
    tip = np.isclose(np.asarray(mesh.points)[:, 0], 100.0)
    # Z88R V15's tip deflection (Euler-Bernoulli gives 1.905 for this beam).
    np.testing.assert_allclose(u[tip, 2].min(), -1.8486980, rtol=1e-7)
    np.testing.assert_array_equal(u[np.asarray(mesh.points)[:, 0] == 0.0], 0.0)
    sig = mesh.cell_data["SIG"][0]
    assert sig.shape == (5, 6)
    assert mesh.cell_data["SIGV"][0].shape == (5,)
    assert (mesh.cell_data["SIGV"][0] > 0).all()


def test_results_can_be_left_out(engine):
    mesh = engine.read(MESHES / "cantilever" / "z88i1.txt", results=False)
    assert "U" not in mesh.point_data and "SIG" not in mesh.cell_data


def test_a_results_file_reads_its_structure_file(engine):
    a = engine.read(MESHES / "cantilever" / "z88o2.txt")
    b = engine.read(MESHES / "cantilever" / "z88i1.txt")
    _same(a, b)


def test_tet10_order_and_legacy_header(engine):
    tets = engine.read(MESHES / "tets" / "z88i1.txt")
    assert [c.type for c in tets.cells] == ["tetra10", "tetra"]
    _valid(tets)
    plate = engine.read(MESHES / "plate_v13" / "z88i1.txt")
    assert [c.type for c in plate.cells] == ["quad8", "triangle6"]
    _valid(plate)
    np.testing.assert_array_equal(plate.cell_data["z88:type"][0], [7])
    np.testing.assert_array_equal(plate.cell_data["z88:type"][1], [14])


def test_cylindrical_input_is_converted(engine):
    mesh = engine.read(MESHES / "polar" / "z88i1.txt")
    pts = np.asarray(mesh.points)
    r = np.hypot(pts[:, 0], pts[:, 1])
    assert np.all((r > 0.999) & (r < 2.001))
    assert np.any(np.all(np.isclose(pts, [0.0, 2.0, 0.0]), axis=1))


@pytest.mark.parametrize("deck", DECKS, ids=[d.name for d in DECKS])
def test_writers_agree_and_round_trip(deck, tmp_path):
    mesh = meshioplusplus.z88.read(deck / "z88i1.txt", results=False)
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    meshioplusplus.z88.write(tmp_path / "a" / "z88i1.txt", mesh, stubs=True)
    py_z88.write(tmp_path / "b" / "z88i1.txt", mesh, stubs=True)
    written = sorted(p.name for p in (tmp_path / "a").iterdir())
    assert written == sorted(p.name for p in (tmp_path / "b").iterdir())
    for name in written:
        assert (tmp_path / "a" / name).read_bytes() == (
            tmp_path / "b" / name
        ).read_bytes(), name
    if "z88:bc:u" not in mesh.point_data:
        assert (tmp_path / "a" / "z88i2.txt").read_text() == "0\n"
    _same(mesh, meshioplusplus.z88.read(tmp_path / "a" / "z88i1.txt"))


def test_basename_dispatch_beats_the_txt_extension(tmp_path):
    mesh = meshioplusplus.read(MESHES / "tets" / "z88i1.txt")
    assert [c.type for c in mesh.cells] == ["tetra10", "tetra"]
    out = tmp_path / "Z88I1.TXT"
    meshioplusplus.write(out, mesh)
    assert meshioplusplus.read(out).cells[0].type == "tetra10"
    assert meshioplusplus.sniff_format(out) == "z88"
    assert _core.sniff_format(str(out)) == "z88"


def test_cells_without_a_z88_type_are_dropped(engine, tmp_path):
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]),
        [("tetra", [[0, 1, 2, 3]]), ("pyramid", [[0, 1, 4, 2, 3]])],
    )
    engine.write(tmp_path / "z88i1.txt", mesh)
    back = meshioplusplus.z88.read(tmp_path / "z88i1.txt")
    assert [c.type for c in back.cells] == ["tetra"]
    only = meshioplusplus.Mesh(mesh.points, [("pyramid", [[0, 1, 4, 2, 3]])])
    with pytest.raises(meshioplusplus.WriteError):
        engine.write(tmp_path / "z88i1.txt", only)


def test_errors(engine, tmp_path):
    deck = tmp_path / "z88i1.txt"
    deck.write_text("3 1 1 3 0\n1 3 0 0 0\n1 99\n1\n")
    with pytest.raises(meshioplusplus.ReadError, match="unknown type 99"):
        engine.read(deck)
    deck.write_text(
        "3 4 1 12 0\n1 3 0 0 0\n2 3 1 0 0\n3 3 0 1 0\n4 3 0 0 1\n1 17\n1 2 3\n"
    )
    with pytest.raises(meshioplusplus.ReadError, match="cut short"):
        engine.read(deck)
    shutil.copy(MESHES / "cantilever" / "z88o2.txt", tmp_path / "z88o2.txt")
    (tmp_path / "z88i1.txt").unlink()
    with pytest.raises(meshioplusplus.ReadError, match="no z88i1.txt"):
        engine.read(tmp_path / "z88o2.txt")


def test_beam_and_truss_results(engine):
    # Z88R's own output for the frame: three 3-D beams (type 2) and a truss (4).
    mesh = engine.read(MESHES / "frame" / "z88i1.txt")
    np.testing.assert_array_equal(mesh.cell_data["z88:type"][0], [2, 2, 2, 4])
    sig = mesh.cell_data["SIGXX"][0]
    np.testing.assert_allclose(sig, [0.2086, -9.689, -9.716, 13.41])
    # Beam-only labels are NaN on the truss; no tensor for line elements.
    assert np.isnan(mesh.cell_data["SIGYY1"][0][3])
    np.testing.assert_allclose(mesh.cell_data["SIGYY1"][0][:3], [-10.30, -8.361, 6.743])
    assert "SIG" not in mesh.cell_data


@pytest.mark.parametrize(
    "deck, labels",
    [
        ("plate", ["MXX", "MYY", "MXY", "QYZ", "QZX", "SIGXX", "TAUXZ", "SIGV"]),
        ("shell", ["SIGXX", "SIGYY", "TAUXY", "SIGV"]),
        ("torus", ["SIGRR", "SIGZZ", "TAURZ", "SIGTE", "SIGV"]),
    ],
)
def test_plate_shell_and_torus_results(engine, deck, labels):
    mesh = engine.read(MESHES / deck / "z88i1.txt")
    for name in labels:
        values = mesh.cell_data[name][0]
        assert values.shape == (len(mesh.cells[0].data),)
        assert np.isfinite(values).all(), name
    assert "SIG" not in mesh.cell_data


def test_nodal_forces_and_constraints(engine):
    mesh = engine.read(MESHES / "frame" / "z88i1.txt")
    u, f, forces = (mesh.point_data[k] for k in ("z88:bc:u", "z88:bc:f", "F"))
    assert u.shape == f.shape == forces.shape == (4, 6)
    feet = np.isclose(np.asarray(mesh.points)[:, 2], 0.0)
    np.testing.assert_array_equal(u[feet], 0.0)
    assert np.isnan(u[~feet]).all()
    # The sway load is the only applied force; on free degrees of freedom the
    # nodal sums Z88R prints equal it, and the reactions balance it.
    loaded = ~np.isnan(f)
    assert loaded.sum() == 1 and f[loaded][0] == 1000.0
    # (to the precision Z88R prints them, about 1e-6 relative).
    np.testing.assert_allclose(forces[~feet], np.nan_to_num(f[~feet]), atol=1e-3)
    np.testing.assert_allclose(forces.sum(axis=0)[:3], 0.0, atol=1e-3)


def test_materials_element_parameters_and_integration(engine):
    mesh = engine.read(MESHES / "cantilever" / "z88i1.txt", results=False)
    np.testing.assert_array_equal(mesh.cell_data["z88:material"][0], [51] * 5)
    np.testing.assert_array_equal(mesh.cell_data["z88:E"][0], [210000.0] * 5)
    np.testing.assert_array_equal(mesh.cell_data["z88:nu"][0], [0.3] * 5)
    elp = mesh.cell_data["z88:elp"][0]
    assert elp.shape == (5, 12)
    # The file gives the seven section values; the rest are absent (NaN).
    np.testing.assert_array_equal(elp[:, :7], 0.0)
    assert np.isnan(elp[:, 7:]).all()
    np.testing.assert_array_equal(mesh.cell_data["z88:int"][0], [[3, 3]] * 5)
    assert "U" not in mesh.point_data and "z88:bc:u" in mesh.point_data


def test_writer_writes_the_deck(engine, tmp_path):
    mesh = engine.read(MESHES / "frame" / "z88i1.txt")
    engine.write(tmp_path / "z88i1.txt", mesh)
    names = sorted(p.name for p in tmp_path.iterdir())
    assert names == [
        "51.txt",
        "z88elp.txt",
        "z88i1.txt",
        "z88i2.txt",
        "z88int.txt",
        "z88mat.txt",
    ]
    assert (tmp_path / "z88mat.txt").read_text().split() == ["1", "1", "4", "51.txt"]
    back = meshioplusplus.z88.read(tmp_path / "z88i1.txt")
    for name in ("z88:bc:u", "z88:bc:f"):
        np.testing.assert_array_equal(back.point_data[name], mesh.point_data[name])
    for name in ("z88:E", "z88:nu", "z88:elp", "z88:int", "z88:type"):
        np.testing.assert_array_equal(back.cell_data[name][0], mesh.cell_data[name][0])


def test_writer_keeps_3d_types_in_a_3d_file(engine, tmp_path):
    # The shells lie in z = 0, but type 24 exists only in a 3-D file.
    mesh = engine.read(MESHES / "shell" / "z88i1.txt")
    engine.write(tmp_path / "z88i1.txt", mesh)
    assert (tmp_path / "z88i1.txt").read_text().split()[0] == "3"
    back = meshioplusplus.z88.read(tmp_path / "z88i1.txt")
    np.testing.assert_array_equal(back.cell_data["z88:type"][0], [24] * 4)


def test_aurora_sets_become_regions(engine, tmp_path):
    # Z88Aurora's z88sets.txt (layout of a real Aurora V2 project, checked
    # outside the repository): element and node sets become regions; surface
    # sets, which name Aurora's surface mesh rather than the elements, do not.
    shutil.copy(MESHES / "tets" / "z88i1.txt", tmp_path / "z88structure.txt")
    (tmp_path / "z88sets.txt").write_text(
        "3\n"
        '#ELEMENTS MATERIAL 1 2 "MatSet1"\n         1          2\n'
        '#SURFACE CONSTRAINT 3 1 "Set3"\n         2          0         12         7'
        "          8          9\n"
        '#NODES CONSTRAINT 4 3 "FIX"\n         1          2          3\n'
    )
    mesh = engine.read(tmp_path / "z88structure.txt")
    regions = {(r.kind, r.name): r for r in mesh.regions}
    assert sorted(regions) == [("cell", "MatSet1"), ("point", "FIX")]
    assert (regions[("cell", "MatSet1")].tag, regions[("cell", "MatSet1")].dim) == (
        1,
        3,
    )
    np.testing.assert_array_equal(regions[("cell", "MatSet1")].entries, [0, 1])
    np.testing.assert_array_equal(regions[("point", "FIX")].entries, [0, 1, 2])
