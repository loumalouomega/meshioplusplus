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
    # a face set: (cell, Marc face number), not mapped to facets
    np.testing.assert_array_equal(mesh.field_data["marc:face_set:skin"], [[0, 1]])
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
    # The deck is written since v16.17.0, the post file since its writer joined.
    assert "marc" in formats["writable"] and "marc_t19" in formats["writable"]
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
    # remeshing flagged, but no model blocks repeated: the model stays
    assert len(engine.read_t19(remesh, time_step=0).points) == 12
    junk = tmp_path / "junk.t19"
    junk.write_text("hello\n")
    with pytest.raises(
        meshioplusplus.ReadError, match="not a Marc formatted post file"
    ):
        engine.read_t19(junk)


def test_include_files_volume_b_types_and_side_sets(engine):
    """include_main.dat INCLUDEs its coordinates, which INCLUDE its sets from a
    subdirectory; a Herrmann brick and quad keep their corners (the pressure
    nodes have no coordinates), a composite brick reads as one, and a 3-node
    rebar line lists its middle node second."""
    mesh = engine.read(MARC / "include_main.dat")
    assert [(b.type, len(b.data)) for b in mesh.cells] == [
        ("hexahedron", 2),
        ("line3", 1),
        ("quad", 1),
    ]
    pts = np.asarray(mesh.points)
    line = mesh.cells[1].data[0]
    np.testing.assert_array_equal(pts[line[2]], [0.5, 0.0, 1.0])  # the middle
    np.testing.assert_array_equal(
        mesh.field_data["marc:edge_set:rim"], [[3, 1], [3, 3]]
    )
    np.testing.assert_array_equal(
        mesh.field_data["marc:face_set:top"], [[0, 2], [1, 2]]
    )
    assert [(r.kind, r.name, r.entries.tolist()) for r in mesh.regions] == [
        ("cell", "bricks", [0, 1])
    ]


def test_remeshing_increment_brings_its_mesh(engine):
    """The second increment of remesh.t19 remeshes: its step is the new model
    (three bricks, with an edge and a face set), the first keeps the old."""
    before = engine.read_t19(MARC / "remesh.t19", time_step=0)
    after = engine.read_t19(MARC / "remesh.t19", time_step=1)
    assert [len(b.data) for b in before.cells] == [2]
    assert [len(b.data) for b in after.cells] == [3]
    assert len(after.points) == 16
    assert after.point_data["Displacement"].shape == (16, 3)
    np.testing.assert_array_equal(
        after.field_data["marc:edge_set:tip_edges"], [[2, 2], [2, 5]]
    )
    np.testing.assert_array_equal(after.field_data["marc:face_set:tip_face"], [[2, 4]])
    assert "marc:edge_set:tip_edges" not in before.field_data


def test_include_and_remesh_engines_agree():
    _same(
        meshioplusplus.marc.read(MARC / "include_main.dat"),
        _Python.read(MARC / "include_main.dat"),
    )
    for step in (0, 1):
        _same(
            meshioplusplus.marc.read_t19(MARC / "remesh.t19", time_step=step),
            _Python.read_t19(MARC / "remesh.t19", time_step=step),
        )


# --------------------------------------------------------------------------- #
# Writing (v16.17.0)                                                          #
# --------------------------------------------------------------------------- #

WRITERS = {
    "core": lambda p, m: _core.marc_write(str(p), m),
    "python": lambda p, m: py_marc.write(p, m),
}


def _set_members(mesh):
    return {
        (r.kind, r.name): sorted(np.asarray(r.entries).ravel().tolist())
        for r in mesh.regions
    }


@pytest.mark.parametrize("path", DECKS, ids=[p.name for p in DECKS])
def test_decks_round_trip_and_the_writers_agree(path, tmp_path):
    mesh = py_marc.read(path)
    written = {}
    for name, writer in WRITERS.items():
        out = tmp_path / name / "deck.dat"
        out.parent.mkdir()
        writer(out, mesh)
        written[name] = out.read_bytes()
    assert written["core"] == written["python"]
    back = py_marc.read(tmp_path / "core" / "deck.dat")
    # Types with extra nodes (Herrmann elements) come back as the default type.
    if path.name == "include_main.dat":
        mesh.cell_data.pop("marc:type")
        back.cell_data.pop("marc:type")
    _same(mesh, back)
    for a, b in zip(mesh.cell_data["marc:element"], back.cell_data["marc:element"]):
        np.testing.assert_array_equal(a, b)
    assert _set_members(mesh) == _set_members(back)
    for name, value in mesh.field_data.items():
        if "_set:" in name:
            np.testing.assert_array_equal(value, back.field_data[name])


def _cells(types, points):
    mesh = meshioplusplus.Mesh(
        points, [(t, np.arange(n)[None] + off) for t, n, off in types]
    )
    return mesh


@pytest.mark.parametrize("engine", sorted(WRITERS))
def test_default_types_follow_the_mesh(engine, tmp_path):
    # A planar mesh gets plane-strain elements, a 3-D one shells and solids;
    # a pyramid is a degenerate brick of type 7.
    flat = _cells(
        [("quad", 4, 0), ("triangle", 3, 0)],
        np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float),
    )
    WRITERS[engine](tmp_path / "flat.dat", flat)
    back = py_marc.read(tmp_path / "flat.dat")
    assert [int(t[0]) for t in back.cell_data["marc:type"]] == [11, 6]
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 1], [0, 0, 2]],
        dtype=float,
    )
    solid = _cells([("pyramid", 5, 0), ("triangle", 3, 0)], pts)
    WRITERS[engine](tmp_path / "solid.dat", solid)
    back = py_marc.read(tmp_path / "solid.dat")
    assert [c.type for c in back.cells] == ["pyramid", "triangle"]
    np.testing.assert_array_equal(back.cells[0].data, solid.cells[0].data)
    assert [int(t[0]) for t in back.cell_data["marc:type"]] == [7, 138]


@pytest.mark.parametrize("engine", sorted(WRITERS))
def test_set_names_read_back_as_sets(engine, tmp_path):
    from meshioplusplus._regions import Region

    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    mesh = meshioplusplus.Mesh(pts, [("tetra", np.array([[0, 1, 2, 3]]))])
    mesh.regions = [
        Region("end", "point", np.array([0, 1])),
        Region("two words", "cell", np.array([0])),
        Region("Two_words", "cell", np.array([0])),
    ]
    WRITERS[engine](tmp_path / "names.dat", mesh)
    back = _set_members(py_marc.read(tmp_path / "names.dat"))
    assert back == {
        ("point", "end_set"): [0, 1],
        ("cell", "Two_words"): [0],
        ("cell", "two_words_2"): [0],
    }


def test_dat_writes_tecplot_unless_marc_is_named(tmp_path):
    tet = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    meshioplusplus.write(tmp_path / "a.dat", tet)
    assert not py_marc.is_marc_deck((tmp_path / "a.dat").read_text())
    deck = py_marc.read(MARC / "hex20.dat")
    meshioplusplus.write(tmp_path / "b.dat", deck, file_format="marc")
    assert py_marc.is_marc_deck((tmp_path / "b.dat").read_text())
    from meshioplusplus._sniff import sniff_format

    assert sniff_format(tmp_path / "b.dat") == "marc"
    _same(deck, meshioplusplus.read(tmp_path / "b.dat"))


# ---------------------------------------------------------------------------
# The post-file writer: both engines write the same bytes, and either reader
# gets the increments back (the reals to the e13.6 digits the file holds).
# ---------------------------------------------------------------------------

from meshioplusplus.marc import _t19 as py_t19  # noqa: E402

T19_WRITERS = {
    "core": lambda path, mesh: _core.marc_t19_write(str(path), mesh),
    "python": py_t19.write,
}


@pytest.mark.parametrize("name", ["results.t19", "remesh.t19"])
def test_the_t19_writers_agree_and_round_trip(name, tmp_path):
    path = MARC / name
    mesh = py_marc.read_t19(path, time_step=0)
    written = {}
    for key, writer in T19_WRITERS.items():
        out = tmp_path / key / "job.t19"
        out.parent.mkdir()
        writer(out, mesh)
        written[key] = out.read_bytes()
    assert written["core"] == written["python"]
    for read in (py_marc.read_t19, meshioplusplus.marc.read_t19):
        back = read(tmp_path / "core" / "job.t19")
        np.testing.assert_allclose(back.points, mesh.points)
        assert [(c.type, c.data.tolist()) for c in back.cells] == [
            (c.type, c.data.tolist()) for c in mesh.cells
        ]
        assert set(back.point_data) == set(mesh.point_data)
        for key, value in mesh.point_data.items():
            np.testing.assert_allclose(back.point_data[key], value, rtol=1e-5)
        assert set(back.cell_data) == set(mesh.cell_data)
        for key in mesh.cell_data:
            for x, y in zip(mesh.cell_data[key], back.cell_data[key]):
                np.testing.assert_allclose(x, y, rtol=1e-5, equal_nan=True)
        assert _set_members(back) == _set_members(mesh)
        assert float(back.field_data["meshio:time"][0]) == pytest.approx(
            float(mesh.field_data["meshio:time"][0])
        )


def test_a_t19_series_is_an_increment_per_mesh(tmp_path):
    steps = []
    for k, step in enumerate((0, 1)):
        mesh = py_marc.read_t19(MARC / "results.t19", time_step=step)
        steps.append((0.5 * (k + 1), mesh))
    out = tmp_path / "series.t19"
    meshioplusplus.write_sequence(out, iter(steps), file_format="marc_t19")
    assert meshioplusplus.marc.time_values(out) == [0.5, 1.0]
    for k, (_, mesh) in enumerate(steps):
        back = meshioplusplus.marc.read_t19(out, time_step=k)
        np.testing.assert_allclose(
            back.point_data["Displacement"], mesh.point_data["Displacement"], rtol=1e-5
        )
    with py_t19.SeriesWriter(tmp_path / "py.t19") as w:
        for time, mesh in steps:
            w.write(time, mesh)
    assert (tmp_path / "py.t19").read_bytes() == out.read_bytes()
    other = py_marc.read_t19(MARC / "results.t19")
    other.point_data["extra"] = np.zeros(len(other.points))
    for make in (meshioplusplus.marc.SeriesWriter, py_t19.SeriesWriter):
        with pytest.raises(meshioplusplus.WriteError):
            with make(tmp_path / "bad.t19") as w:
                w.write(0.0, steps[0][1])
                w.write(1.0, other)


def test_a_t19_series_remeshes_when_the_mesh_changes(tmp_path):
    """The two models of remesh.t19 as a series: the second step is written as
    a remeshing increment, in both engines and byte for byte."""
    steps = [
        (0.5 * (k + 1), py_marc.read_t19(MARC / "remesh.t19", time_step=k))
        for k in (0, 1)
    ]
    written = {}
    for key, make in (
        ("core", meshioplusplus.marc.SeriesWriter),
        ("python", py_t19.SeriesWriter),
    ):
        out = tmp_path / f"{key}.t19"
        with make(out) as w:
            for time, mesh in steps:
                w.write(time, mesh)
        written[key] = out.read_bytes()
    assert written["core"] == written["python"]
    assert written["core"].count(b"=beg=51900") == 1
    for read in (py_marc.read_t19, meshioplusplus.marc.read_t19):
        out = tmp_path / "core.t19"
        for k, (_, mesh) in enumerate(steps):
            back = read(out, time_step=k)
            assert back.points.shape == mesh.points.shape
            np.testing.assert_allclose(back.points, mesh.points)
            assert [len(b.data) for b in back.cells] == [
                len(b.data) for b in mesh.cells
            ]
            np.testing.assert_allclose(
                back.point_data["Displacement"],
                mesh.point_data["Displacement"],
                rtol=1e-5,
            )
            assert _set_members(back) == _set_members(mesh)
            assert sorted(n for n in back.field_data if "_set:" in n) == sorted(
                n for n in mesh.field_data if "_set:" in n
            )
    # An unchanged model is not repeated.
    same = tmp_path / "same.t19"
    with py_t19.SeriesWriter(same) as w:
        w.write(0.0, steps[0][1])
        w.write(1.0, steps[0][1])
    assert b"=beg=51900" not in same.read_bytes()


def test_layers_custom_names_and_nan(tmp_path):
    mesh = py_marc.read_t19(MARC / "results.t19")
    stress = mesh.cell_data["Stress"]
    mesh.cell_data["Custom Scalar"] = [
        np.arange(16.0).reshape(2, 8) + 100 * k for k, _ in enumerate(stress)
    ]
    mesh.cell_data["Stress@layer2"] = [np.array(a) * 2 for a in stress]
    mesh.field_data["marc:layout:Stress@layer2"] = mesh.field_data["marc:layout:Stress"]
    mesh.field_data["marc:layout:Custom Scalar"] = np.array([8, 1])
    mesh.cell_data["Hole"] = [np.full((2, 8), np.nan) for _ in stress]
    mesh.field_data["marc:layout:Hole"] = np.array([8, 1])
    for key, writer in T19_WRITERS.items():
        out = tmp_path / f"{key}.t19"
        writer(out, mesh)
        back = py_marc.read_t19(out)
        np.testing.assert_allclose(
            back.cell_data["Custom Scalar"][0], mesh.cell_data["Custom Scalar"][0]
        )
        np.testing.assert_allclose(
            back.cell_data["Stress@layer2"][0],
            mesh.cell_data["Stress@layer2"][0],
            rtol=1e-5,
        )
        assert np.isnan(back.cell_data["Hole"][0]).all()
