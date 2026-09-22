import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._exceptions import ReadError
from meshioplusplus.unv import _unv

from . import helpers

FIXTURES = pathlib.Path(__file__).resolve().parent / "meshes" / "unv"

# wedge15 (UNV descriptor 113) and quad9 are not shared helper fixtures.
wedge15_mesh = meshioplusplus.Mesh(
    np.arange(45.0).reshape(15, 3),
    [("wedge15", np.arange(15).reshape(1, 15))],
)
quad9_mesh = meshioplusplus.Mesh(
    np.arange(27.0).reshape(9, 3),
    [("quad9", np.arange(9).reshape(1, 9))],
)


def _both_readers(path, **kwargs):
    """The C++ reader and the Python reference, on the same file."""
    return _core.unv_read(
        str(path), False, None, kwargs.get("time_step", 0)
    ), _unv.read(path, **kwargs)


# --------------------------------------------------------------------------- #
# Round trips                                                                 #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.line_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.triangle6_mesh,
        helpers.quad_mesh,
        helpers.quad8_mesh,
        quad9_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
        helpers.wedge_mesh,
        wedge15_mesh,
        helpers.pyramid_mesh,
        helpers.pyramid13_mesh,
    ],
)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.unv.write, meshioplusplus.unv.read, mesh, 1.0e-12
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.unv")
    helpers.generic_io(tmp_path / "test.0.unv")
    helpers.generic_io(tmp_path / "test.uff")


def test_groups(tmp_path):
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    mesh.point_sets = {"corners": np.array([0, 2])}
    mesh.cell_sets = {"all": [np.array([0, 1])], "first": [np.array([0])]}
    for write, read in (
        (meshioplusplus.unv.write, meshioplusplus.unv.read),
        (_unv.write, _unv.read),
    ):
        p = tmp_path / "g.unv"
        write(p, mesh)
        text = p.read_text()
        # 7 = node, 8 = finite element (SDRL 2467 entity type codes)
        assert "         7         1         0         0" in text
        assert "         8         1         0         0" in text
        out = read(p)
        assert np.array_equal(out.point_sets["corners"], [0, 2])
        assert np.array_equal(out.cell_sets["all"][0], [0, 1])
        assert np.array_equal(out.cell_sets["first"][0], [0])


def test_region_tags_and_mixed_groups_round_trip(tmp_path):
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    R = meshioplusplus.Region
    mesh.regions = [
        R("wall", "point", [0, 1], tag=7),
        R("wall", "cell", [1], dim=2, tag=7),
        R("empty", "cell", [], dim=2, tag=3),
        R("untagged", "point", [3]),
    ]
    for write in (meshioplusplus.unv.write, _unv.write):
        p = tmp_path / "t.unv"
        write(p, mesh)
        for read in (meshioplusplus.unv.read, _unv.read):
            got = {(r.name, r.kind): r for r in read(p).regions}
            assert got[("wall", "point")].tag == 7
            assert list(got[("wall", "point")].entries) == [0, 1]
            assert got[("wall", "cell")].tag == 7
            assert got[("wall", "cell")].dim == 2
            assert list(got[("wall", "cell")].entries) == [1]
            assert len(got[("empty", "cell")].entries) == 0
            assert got[("empty", "cell")].tag == 3
            # the next free number goes to the untagged group
            assert got[("untagged", "point")].tag == 1


def test_side_regions_are_dropped_with_a_warning(tmp_path, capfd):
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    mesh.regions = [meshioplusplus.Region("face", "side", [[0, 1]], dim=2)]
    p = tmp_path / "s.unv"
    for write in (_unv.write, _core.unv_write):
        write(str(p), mesh)
        assert "2467" not in p.read_text()
        assert "side region" in capfd.readouterr().err


def test_registry_metadata_carries_groups_and_steps():
    """The generic C++ path (C API, CLI, WASM) sees the groups and the steps."""
    meta = meshioplusplus.read_metadata(FIXTURES / "modes_2414.unv")
    assert meta["format"] == "unv"
    assert meta["time_values"] == [12.5, 31.0, 47.25, 1.0]
    regions = {
        (r["name"], r["kind"]): (r["tag"], r["num_entries"]) for r in meta["regions"]
    }
    assert regions[("right_face", "point")] == (11, 4)
    assert regions[("left", "cell")] == (10, 1)


def _field_mesh():
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    mesh.point_data = {
        "temp": np.array([1.0, 2, 3, 4, 5]),
        "disp": np.arange(15.0).reshape(5, 3),
    }
    mesh.cell_data = {
        "stress": [np.arange(12.0).reshape(2, 6)],
        "unv:pid": [np.array([3, 4])],
        "unv:mid": [np.array([8, 9])],
    }
    return mesh


def _assert_fields(out, rtol=1e-12):
    assert np.allclose(out.point_data["temp"], [1, 2, 3, 4, 5], rtol=rtol)
    assert np.allclose(out.point_data["disp"], np.arange(15.0).reshape(5, 3), rtol=rtol)
    assert np.allclose(
        out.cell_data["stress"][0], np.arange(12.0).reshape(2, 6), rtol=rtol
    )
    assert np.array_equal(out.cell_data["unv:pid"][0], [3, 4])
    assert np.array_equal(out.cell_data["unv:mid"][0], [8, 9])


def test_fields_roundtrip(tmp_path):
    """Dataset 2414 fields (scalar/vector/tensor) at nodes and elements."""
    mesh = _field_mesh()
    p = tmp_path / "f.unv"
    meshioplusplus.unv.write(p, mesh)
    _assert_fields(meshioplusplus.unv.read(p))


def test_fields_keep_double_precision(tmp_path):
    mesh = _field_mesh()
    mesh.point_data["temp"] = mesh.point_data["temp"] + 1.0 / 3.0
    p = tmp_path / "f.unv"
    meshioplusplus.unv.write(p, mesh)
    assert np.allclose(
        meshioplusplus.unv.read(p).point_data["temp"], mesh.point_data["temp"]
    )
    assert np.abs(meshioplusplus.unv.read(p).point_data["temp"][0] - 4.0 / 3.0) < 1e-11


def test_symmetric_tensor_is_written_in_file_order(tmp_path):
    """SDRL: a symmetric tensor is stored Sxx Sxy Syy Sxz Syz Szz."""
    mesh = meshioplusplus.Mesh(np.zeros((1, 3)), [("vertex", np.array([[0]]))])
    mesh.point_data = {
        "s": np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]])
    }  # xx yy zz xy yz zx
    p = tmp_path / "t.unv"
    meshioplusplus.unv.write(p, mesh)
    lines = p.read_text().splitlines()
    start = lines.index("  2414")
    assert lines[start + 14].strip() == "1"  # record 14: node 1
    vals = [float(v) for v in lines[start + 15].split()]
    assert vals == [1.0, 4.0, 2.0, 6.0, 5.0, 3.0]


def test_writers_are_byte_identical(tmp_path):
    mesh = _field_mesh()
    mesh.point_sets = {"corners": np.array([0, 2])}
    mesh.field_data = {
        "meshio:time": np.array([0.25]),
        "unv:analysis": np.array([4]),
        "unv:step": np.array([3]),
        "unv:units": np.array([1]),
        "unv:unit_factors": np.array([1.0, 1.0, 1.0, 273.15]),
    }
    for code_aster in (False, True):
        p_cpp = tmp_path / "cpp.unv"
        p_py = tmp_path / "py.unv"
        _core.unv_write(str(p_cpp), mesh, code_aster=code_aster)
        _unv.write(p_py, mesh, code_aster=code_aster)
        assert p_cpp.read_bytes() == p_py.read_bytes()
        for out in _both_readers(p_cpp):
            _assert_fields(out, rtol=1e-12 if not code_aster else 1e-5)
            assert out.field_data["unv:analysis"][0] == 4
            assert out.field_data["unv:step"][0] == 3
            assert np.isclose(out.field_data["meshio:time"][0], 0.25)
            assert out.field_data["unv:units"][0] == 1


def test_fields_code_aster(tmp_path):
    """Code-Aster mode emits the legacy datasets 55 (nodes) / 56 (elements)."""
    mesh = _field_mesh()
    p = tmp_path / "ca.unv"
    meshioplusplus.unv.write(p, mesh, code_aster=True)
    text = p.read_text()
    assert "\n    55\n" in text and "\n    56\n" in text and "\n    57\n" not in text
    _assert_fields(meshioplusplus.unv.read(p), rtol=1e-5)


def test_node_dataset_781(tmp_path):
    mesh = _field_mesh()
    p = tmp_path / "n.unv"
    meshioplusplus.unv.write(p, mesh, node_dataset=781)
    assert "\n   781\n" in p.read_text()
    _assert_fields(meshioplusplus.unv.read(p))


def test_nan_entities_are_not_written(tmp_path):
    mesh = _field_mesh()
    mesh.point_data = {"t": np.array([1.0, np.nan, 3.0, 4.0, 5.0])}
    mesh.cell_data = {}
    p = tmp_path / "n.unv"
    meshioplusplus.unv.write(p, mesh)
    assert "NAN" not in p.read_text().upper()
    for out in _both_readers(p):
        assert np.isnan(out.point_data["t"][1])
        assert out.point_data["t"][2] == 3.0


# --------------------------------------------------------------------------- #
# Malformed and unsupported input                                              #
# --------------------------------------------------------------------------- #

_TWO_NODES = (
    "    -1\n  2411\n"
    "         1         1         1        11\n   0.0   0.0   0.0\n"
    "         2         1         1        11\n   1.0   0.0   0.0\n"
    "    -1\n"
)


@pytest.mark.parametrize(
    "record",
    [
        # unknown descriptor
        "         1       999         1         1        11         2\n         1         2\n",
        # a linear descriptor with the wrong node count (gmsh's quad9 under 94)
        "         1        94         1         1        11         2\n         1         2\n",
    ],
)
def test_unsupported_elements_are_skipped(tmp_path, record):
    p = tmp_path / "bad.unv"
    p.write_text(_TWO_NODES + "    -1\n  2412\n" + record + "    -1\n")
    for out in _both_readers(p):
        assert len(out.points) == 2 and len(out.cells) == 0


def test_undefined_node_is_a_read_error(tmp_path):
    p = tmp_path / "bad.unv"
    p.write_text(
        _TWO_NODES
        + "    -1\n  2412\n         1        21         1         1        11         2\n"
        "         0         1         1\n         1         7\n    -1\n"
    )
    with pytest.raises(Exception, match="undefined node 7"):
        _core.unv_read(str(p))
    with pytest.raises(ReadError, match="undefined node 7"):
        _unv.read(p)


def test_time_step_out_of_range(tmp_path):
    with pytest.raises(Exception, match="out of range"):
        _core.unv_read(str(FIXTURES / "transient_2414.unv"), False, None, 3)
    with pytest.raises(ReadError, match="out of range"):
        _unv.read(FIXTURES / "transient_2414.unv", time_step=3)


# --------------------------------------------------------------------------- #
# Independent fixtures                                                         #
# --------------------------------------------------------------------------- #

# VTK edge list of each quadratic type: meshio node n_corners + k is the mid-node of
# edge k.
_EDGES = {
    "triangle6": (3, [(0, 1), (1, 2), (2, 0)]),
    "quad8": (4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    "tetra10": (4, [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]),
    "wedge15": (
        6,
        [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
    ),
    "hexahedron20": (
        8,
        [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
        + [(0, 4), (1, 5), (2, 6), (3, 7)],
    ),
}


def _signed_volume(p):
    return np.einsum(
        "ij,ij->i", p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])
    )


@pytest.mark.parametrize("name", ["gmsh_hex20", "gmsh_wedge15", "gmsh_tet10"])
def test_gmsh_files_match_their_msh_twin(name):
    """gmsh wrote the same mesh as .unv and .msh: identical cells, node for node."""
    ref = meshioplusplus.read(FIXTURES / f"{name}.msh")
    for mesh in _both_readers(FIXTURES / f"{name}.unv"):
        for cb in ref.cells:
            if cb.type not in _EDGES or cb.type == "triangle6":
                continue
            ours = [c for c in mesh.cells if c.type == cb.type][0]
            key = lambda pts: tuple(sorted(map(tuple, pts.round(9))))  # noqa: E731
            by_corners = {key(mesh.points[row]): mesh.points[row] for row in ours.data}
            for row in cb.data:
                assert np.allclose(by_corners[key(ref.points[row])], ref.points[row])
        # Straight-edged geometry: every mid-node sits at its edge's midpoint.
        for cb in mesh.cells:
            if cb.type not in _EDGES:
                continue
            nc, edges = _EDGES[cb.type]
            p = mesh.points[cb.data]
            for k, (a, b) in enumerate(edges):
                assert np.allclose(p[:, nc + k], 0.5 * (p[:, a] + p[:, b]))
            if cb.type == "tetra10":
                assert (_signed_volume(p[:, :4]) > 0).all()
        # Groups: the physical groups, each as its nodes and its elements.
        names = {r.name for r in mesh.regions}
        assert names == {n for n in ref.field_data}
        for r in mesh.regions:
            assert r.tag == ref.field_data[r.name][0]


def test_modes_2414():
    for step, freq in enumerate((12.5, 31.0, 47.25)):
        for mesh in _both_readers(FIXTURES / "modes_2414.unv", time_step=step):
            assert mesh.field_data["unv:analysis"][0] == 2
            assert mesh.field_data["unv:step"][0] == step + 1
            assert mesh.field_data["meshio:time"][0] == freq
            disp = mesh.point_data["Mode shape"]
            assert np.allclose(
                disp[:, 2], (step + 1) * 0.1 * np.arange(1, 13), rtol=1e-5
            )
            assert "stress" not in mesh.cell_data
    # The static stress is the fourth step, reordered to xx yy zz xy yz zx.
    for mesh in _both_readers(FIXTURES / "modes_2414.unv", time_step=-1):
        assert mesh.field_data["unv:analysis"][0] == 1
        assert np.allclose(
            mesh.cell_data["stress"][0], [[1, 2, 3, 4, 5, 6], [10, 20, 30, 40, 50, 60]]
        )
        assert np.isnan(mesh.cell_data["stress"][1]).all()  # the quad has no stress
    assert meshioplusplus.read_metadata(FIXTURES / "modes_2414.unv")["time_values"] == [
        12.5,
        31.0,
        47.25,
        1.0,
    ]
    for mesh in _both_readers(FIXTURES / "modes_2414.unv"):
        got = {(r.name, r.kind): (r.tag, list(r.entries)) for r in mesh.regions}
        assert got[("left", "cell")] == (10, [0])
        assert got[("right_face", "point")] == (11, [8, 9, 10, 11])
        assert got[("right_face", "cell")] == (11, [2])
        assert got[("origin", "point")] == (12, [0])  # dataset 2429: pairs
        assert np.array_equal(mesh.cell_data["unv:pid"][0], [3, 4])
        assert np.array_equal(mesh.cell_data["unv:mid"][1], [8])


def test_transient_sequence():
    seq = list(meshioplusplus.read_sequence(str(FIXTURES / "transient_2414.unv")))
    assert [t for t, _ in seq] == [0.0, 0.5, 1.0]
    last = seq[-1][1]
    assert np.allclose(last.point_data["Temperature"][:6], 100.0 + np.arange(1, 7))
    assert np.isnan(last.point_data["Temperature"][6:]).all()


def test_legacy_meshioplusplus_layout_still_reads():
    for mesh in _both_readers(FIXTURES / "legacy55.unv"):
        assert np.allclose(mesh.point_data["temp"], np.arange(1, 9))
        assert np.allclose(mesh.cell_data["pressure"][0], [7.5])


def test_units_coordinate_systems_and_legacy_datasets():
    for mesh in _both_readers(FIXTURES / "units_cs_legacy.unv"):
        assert mesh.field_data["unv:units"][0] == 5
        assert np.allclose(mesh.field_data["unv:unit_factors"], [1000, 1000, 1, 273.15])
        # node 4 is (1, 0, 0) in system 2: rotated then shifted by (10, 0, 0)
        assert np.allclose(mesh.points[3], [10.0, 1.0, 0.0])
        assert [c.type for c in mesh.cells] == ["triangle", "line"]
        assert np.array_equal(mesh.cell_data["unv:pid"][0], [6])
        assert np.array_equal(mesh.cell_data["unv:mid"][0], [9])


def test_pyuff_modes_55():
    for step, freq in enumerate((10.5, 22.0)):
        for mesh in _both_readers(FIXTURES / "pyuff_modes.uff", time_step=step):
            assert mesh.field_data["meshio:time"][0] == freq
            disp = mesh.point_data["displacement"]
            assert np.allclose(
                disp, (step + 1) * np.array([[1, 4, 7], [2, 5, 8], [3, 6, 9]])
            )


@pytest.mark.parametrize("name", ["pyuff_frf.uff", "pyuff_frf_58b.uff"])
def test_pyuff_frf_is_a_frequency_sequence(name):
    freq = np.linspace(0.0, 20.0, 9)
    assert np.allclose(
        meshioplusplus.read_metadata(FIXTURES / name)["time_values"], freq
    )
    for step in (0, 4, 8):
        for mesh in _both_readers(FIXTURES / name, time_step=step):
            assert mesh.field_data["unv:analysis"][0] == 5
            assert mesh.field_data["meshio:time"][0] == freq[step]
            h = mesh.point_data["frf_real"] + 1j * mesh.point_data["frf_imag"]
            for node in (1, 2, 3):
                for d in (1, 2, 3):
                    expected = (node + 0.1 * d) * (freq[step] + 1j * (freq[step] + 1))
                    assert np.isclose(h[node - 1, d - 1], expected)
            coh = mesh.point_data["coherence"]
            assert np.isnan(coh[[0, 2]]).all()
            assert np.isclose(coh[1, 2], np.linspace(1.0, 0.5, 9)[step], rtol=1e-6)


def test_58_without_nodes_is_an_error():
    with pytest.raises(Exception, match="no nodes"):
        _core.unv_read(str(FIXTURES / "pyuff_frf_only.uff"))
    with pytest.raises(ReadError, match="no nodes"):
        _unv.read(FIXTURES / "pyuff_frf_only.uff")


def test_pyuff_agrees_on_the_frf():
    pyuff = pytest.importorskip("pyuff")
    sets = pyuff.UFF(str(FIXTURES / "pyuff_frf_58b.uff")).read_sets()
    frfs = [s for s in sets if s["type"] == 58 and s["func_type"] == 4]
    times = meshioplusplus.read_metadata(FIXTURES / "pyuff_frf_58b.uff")["time_values"]
    assert np.allclose(frfs[0]["x"], times)
    for step in range(len(times)):
        mesh = meshioplusplus.read(FIXTURES / "pyuff_frf_58b.uff", time_step=step)
        h = mesh.point_data["frf_real"] + 1j * mesh.point_data["frf_imag"]
        for s in frfs:
            assert np.isclose(h[s["rsp_node"] - 1, s["rsp_dir"] - 1], s["data"][step])


def test_pyuff_reads_our_output(tmp_path):
    pyuff = pytest.importorskip("pyuff")
    mesh = _field_mesh()
    p = tmp_path / "ours.unv"
    meshioplusplus.unv.write(p, mesh)
    sets = pyuff.UFF(str(p)).read_sets()
    nodes = next(s for s in sets if s["type"] == 2411)
    assert np.allclose(
        np.column_stack([nodes["x"], nodes["y"], nodes["z"]]), mesh.points
    )
    disp = next(
        s for s in sets if s["type"] == 2414 and s["analysis_dataset_name"] == "disp"
    )
    assert np.allclose(disp["data_at_node"], mesh.point_data["disp"])
    stress = next(
        s for s in sets if s["type"] == 2414 and s["analysis_dataset_name"] == "stress"
    )
    assert len(stress["data_at_element"]) == 2


def test_strict_core_reads_every_fixture(monkeypatch):
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    for path in sorted(FIXTURES.glob("*.u*")):
        if path.name == "pyuff_frf_only.uff":
            continue
        meshioplusplus.unv.read(path)


def test_sniffing(tmp_path):
    from meshioplusplus._sniff import _sniff_format_py, sniff_format

    for path in sorted(FIXTURES.glob("*.u*")):
        blind = tmp_path / "noext"
        blind.write_bytes(path.read_bytes())
        assert sniff_format(blind) == "unv", path.name
        assert _sniff_format_py(blind) == "unv", path.name
