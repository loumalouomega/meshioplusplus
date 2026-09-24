"""MSC Marc input decks (``.dat``) and formatted post files (``.t19``): both
engines on the fixtures ``tools/gen_marc_fixtures.py`` writes (fixed, extended
and free format, sets, a post file with two increments), dispatch between Marc's
and Tecplot's ``.dat``, and the refusals."""

import pathlib
import shutil

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.marc import _marc as py_marc

HERE = pathlib.Path(__file__).parent
MARC = HERE / "meshes" / "marc"
DECKS = sorted(MARC.glob("*.dat"))

# The generator's models, for the expected values.
_spec = __import__("importlib.util").util.spec_from_file_location(
    "gen_marc_fixtures", HERE.parent.parent / "tools" / "gen_marc_fixtures.py"
)
gen = __import__("importlib.util").util.module_from_spec(_spec)
_spec.loader.exec_module(gen)


class _Python:
    read = staticmethod(py_marc.read)
    read_t19 = staticmethod(py_marc.read_t19)
    time_values = staticmethod(py_marc.time_values)


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read signature."""
    return meshioplusplus.marc if request.param == "core" else _Python


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    for data_a, data_b in ((a.point_data, b.point_data), (a.field_data, b.field_data)):
        assert sorted(data_a) == sorted(data_b)
        for name in data_a:
            np.testing.assert_array_equal(data_a[name], data_b[name])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_array_equal(x, y)
    key = lambda r: (r.kind, r.name, r.dim, tuple(r.entries))  # noqa: E731
    assert sorted(map(key, a.regions)) == sorted(map(key, b.regions))


@pytest.mark.parametrize("path", DECKS, ids=[p.name for p in DECKS])
def test_engines_agree_on_decks(path):
    _same(_core.marc_read(str(path)), py_marc.read(path))


def test_engines_agree_on_the_post_file():
    path = MARC / "results.t19"
    assert list(_core.marc_t19_time_values(str(path))) == py_marc.time_values(path)
    for step in range(2):
        _same(
            _core.marc_t19_read(str(path), time_step=step),
            py_marc.read_t19(path, time_step=step),
        )


def _jacobian_signs(mesh):
    out = []
    for block in mesh.cells:
        p = mesh.points[block.data]
        if block.type in ("hexahedron", "hexahedron20"):
            j = np.cross(p[:, 1] - p[:, 0], p[:, 3] - p[:, 0])
            out += list(np.einsum("ij,ij->i", j, p[:, 4] - p[:, 0]))
        elif block.type in ("tetra", "tetra10", "wedge"):
            j = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
            out += list(np.einsum("ij,ij->i", j, p[:, 3] - p[:, 0]))
        elif block.type in ("quad", "quad8", "triangle", "triangle6"):
            out += list(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])[:, 2])
    return np.array(out)


_EDGES = {
    "hexahedron20": [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
    + [(0, 4), (1, 5), (2, 6), (3, 7)],
    "tetra10": [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
    "quad8": [(0, 1), (1, 2), (2, 3), (3, 0)],
}


def _check_midpoints(mesh):
    for block in mesh.cells:
        edges = _EDGES.get(block.type)
        if not edges:
            continue
        p = mesh.points[block.data]
        corners = block.data.shape[1] - len(edges)
        for k, (a, b) in enumerate(edges):
            np.testing.assert_allclose(p[:, corners + k], 0.5 * (p[:, a] + p[:, b]))


@pytest.mark.parametrize("name", ["hex20.dat", "hex20_extended.dat"])
def test_hex20_deck(engine, name):
    """The done-when of roadmap item 1.1: a hex20 deck reads with its node order
    (positive Jacobians, mid-edge nodes at the midpoints) and element sets as
    regions."""
    mesh = engine.read(MARC / name)
    nodes, elements, x0, top = gen.hex20_model()
    assert [b.type for b in mesh.cells] == ["hexahedron20"]
    assert mesh.cell_data["marc:element"][0].tolist() == [100, 200]
    assert mesh.cell_data["marc:type"][0].tolist() == [21, 21]
    assert (_jacobian_signs(mesh) > 0).all()
    _check_midpoints(mesh)
    # Each cell's nodes are the element's, in order; node 140, defined twice,
    # keeps its later coordinates.
    for row, (_, _, ids) in zip(mesh.cells[0].data, elements):
        np.testing.assert_array_equal(mesh.points[row], [nodes[i] for i in ids])
    regions = {(r.kind, r.name): r for r in mesh.regions}
    assert regions[("cell", "left")].entries.tolist() == [0]
    assert regions[("cell", "right")].entries.tolist() == [1]
    assert regions[("cell", "both")].entries.tolist() == [0, 1]
    assert regions[("cell", "both_again")].entries.tolist() == [0, 1]
    assert regions[("cell", "both")].dim == 3
    x0_points = regions[("point", "x0_face")].entries
    np.testing.assert_array_equal(mesh.points[x0_points][:, 0], 0.0)
    assert len(x0_points) == len(x0) == 8
    assert len(regions[("point", "top_face")].entries) == len(top) == 13
    top_not_x0 = mesh.points[regions[("point", "top_not_x0")].entries]
    assert len(top_not_x0) == 10 and (top_not_x0[:, 0] > 0).all()
    assert (top_not_x0[:, 2] == 1.0).all()


def test_fixed_and_extended_decks_agree(engine):
    _same(engine.read(MARC / "hex20.dat"), engine.read(MARC / "hex20_extended.dat"))


def test_free_format_mixed_deck(engine, capfd):
    mesh = engine.read(MARC / "mixed_free.dat")
    err = capfd.readouterr().err
    assert "116" in err and "skipped" in err  # the unknown element type
    assert "FACE SET 'skin'" in err
    assert [b.type for b in mesh.cells] == [
        "hexahedron",
        "wedge",
        "tetra10",
        "quad",
        "line",
    ]
    # Element 2 repeats nodes 9 and 10: the wedge keeps its faces' order.
    assert mesh.cell_data["marc:element"][1].tolist() == [2]
    assert (_jacobian_signs(mesh)[:3] > 0).all()
    _check_midpoints(mesh)
    regions = {(r.kind, r.name): r.entries.tolist() for r in mesh.regions}
    assert regions[("cell", "solids")] == [0, 1, 2]
    assert regions[("point", "base")] == [0, 1, 2, 3]


def test_plane_deck(engine):
    mesh = engine.read(MARC / "plane_quad8.dat")
    assert [b.type for b in mesh.cells] == ["quad8"]
    assert mesh.points.shape == (13, 3) and (mesh.points[:, 2] == 0).all()
    assert (_jacobian_signs(mesh) > 0).all()
    _check_midpoints(mesh)


def test_post_file(engine):
    path = MARC / "results.t19"
    assert engine.time_values(path) == gen.T19_TIMES
    last = engine.read_t19(path, time_step=-1)
    assert last.field_data["meshio:time"].tolist() == [1.0]
    assert last.field_data["marc:increment"].tolist() == [2]
    assert [b.type for b in last.cells] == ["hexahedron"]
    assert (_jacobian_signs(last) > 0).all()
    for step in range(2):
        mesh = engine.read_t19(path, time_step=step)
        u = mesh.point_data["Displacement"]
        expected = [gen.t19_displacement(step, n) for n in range(1, 13)]
        np.testing.assert_allclose(u, expected, rtol=1e-6)
        rf = mesh.point_data["Reaction Force"]
        assert rf.shape == (12, 3) and (rf[[0, 3, 4, 7]] < 0).all()
        # Eight integration points: flattened point-major, with the layout.
        assert mesh.cell_data["Stress"][0].shape == (2, 48)
        np.testing.assert_array_equal(mesh.field_data["marc:layout:Stress"], [8, 6])
        np.testing.assert_array_equal(
            mesh.field_data["marc:layout:Equivalent Von Mises Stress"], [8, 1]
        )
        stress = mesh.cell_data["Stress"][0].reshape(2, 8, 6)
        mises = mesh.cell_data["Equivalent Von Mises Stress"][0]
        assert mises.shape == (2, 8)
        for e in range(2):
            for ip in range(8):
                t = gen.t19_stress(step, e + 1, ip)
                np.testing.assert_allclose(stress[e, ip], t, rtol=1e-6, atol=1e-6)
                assert mises[e, ip] == pytest.approx(gen.mises(stress[e, ip]), rel=1e-5)
    regions = {(r.kind, r.name): r.entries.tolist() for r in last.regions}
    assert regions == {
        ("point", "fixed"): [0, 3, 4, 7],
        ("point", "loaded"): [8, 9, 10, 11],
        ("cell", "second"): [1],
    }


def test_post_file_selective_reads(engine):
    path = MARC / "results.t19"
    only = engine.read_t19(path, points_only=True)
    assert not only.point_data
    assert sorted(k for k in only.cell_data if not k.startswith("marc:")) == []
    stress = engine.read_t19(path, arrays=["Stress"])
    assert not stress.point_data
    assert sorted(stress.cell_data) == ["Stress", "marc:element", "marc:type"]
    with pytest.raises(meshioplusplus.ReadError, match="out of range"):
        engine.read_t19(path, time_step=2)


def test_dispatch_and_sniffing(tmp_path):
    formats = meshioplusplus.formats()
    assert formats["extensions"][".t19"] == ["marc_t19"]
    assert formats["extensions"][".dat"] == ["marc", "tecplot"]
    assert {"marc", "marc_t19"} <= set(formats["readable"])
    assert not {"marc", "marc_t19"} & set(formats["writable"])
    # A Marc .dat reads as Marc; a Tecplot .dat still as Tecplot.
    mesh = meshioplusplus.read(MARC / "hex20.dat")
    assert [b.type for b in mesh.cells] == ["hexahedron20"]
    tec = tmp_path / "tecplot.dat"
    triangle = meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2]]))],
    )
    meshioplusplus.write(tec, triangle, file_format="tecplot")
    assert [b.type for b in meshioplusplus.read(tec).cells] == ["triangle"]
    with pytest.raises(meshioplusplus.ReadError, match="not a Marc input deck"):
        meshioplusplus.marc.read(tec)
    for source, expected in (
        ("hex20_extended.dat", "marc"),
        ("results.t19", "marc_t19"),
    ):
        renamed = tmp_path / f"{expected}.bin"
        shutil.copy(MARC / source, renamed)
        assert meshioplusplus.sniff_format(renamed) == expected
        assert _core.sniff_format(str(renamed)) == expected


def test_sequences_and_metadata(tmp_path):
    path = MARC / "results.t19"
    assert meshioplusplus.read_metadata(path)["time_values"] == gen.T19_TIMES
    steps = list(meshioplusplus.read_sequence(path))
    assert [time for time, _ in steps] == gen.T19_TIMES
    np.testing.assert_array_equal(
        steps[1][1].point_data["Displacement"],
        meshioplusplus.read(path, time_step=1).point_data["Displacement"],
    )


def test_refusals(engine, tmp_path):
    bad = tmp_path / "bad.dat"
    bad.write_text(
        "title   x\nend\nconnectivity\n    0    0    1\n    1    7    1    2\n"
    )
    with pytest.raises(meshioplusplus.ReadError, match="lists 2 of its 8 nodes"):
        engine.read(bad)
    bad.write_text(
        "title   x\nend\nconnectivity\n    0\n    1  134    1    2    3    9\n"
        "coordinates\n    3    3\n"
        + "".join(f"{n:5d}{0.0:10.3f}{0.0:10.3f}{float(n):10.3f}\n" for n in (1, 2, 3))
    )
    with pytest.raises(meshioplusplus.ReadError, match="undefined node"):
        engine.read(bad)
    sets = tmp_path / "sets.dat"
    sets.write_text(
        "end\ncoordinates\n    3    1\n    1        0.        0.        0.\n"
        "define              node                set                 a\n"
        "missing and 1\n"
    )
    mesh = engine.read(sets)  # a line of names that are no earlier set ends it
    assert [(r.kind, r.name, r.entries.tolist()) for r in mesh.regions] == [
        ("point", "a", [])
    ]
    sets.write_text(sets.read_text().replace("missing and 1", "    1 and 7"))
    with pytest.raises(meshioplusplus.ReadError, match="undefined node 7"):
        engine.read(sets)
    remesh = tmp_path / "remesh.t19"
    text = (MARC / "results.t19").read_text()
    head = "=beg=51701 (Integer Increment Verification Data)"
    at = text.index(head)
    line = text.index("\n", at) + 1
    remesh.write_text(text[:line] + f"{1:13d}" + text[line + 13 :])
    with pytest.raises(meshioplusplus.ReadError, match="remesh"):
        engine.read_t19(remesh)
    junk = tmp_path / "junk.t19"
    junk.write_text("hello\n")
    with pytest.raises(
        meshioplusplus.ReadError, match="not a Marc formatted post file"
    ):
        engine.read_t19(junk)
