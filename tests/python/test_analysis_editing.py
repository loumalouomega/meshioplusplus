"""Tests for the roadmap §5 "Analysis and editing" operations of v16.23.0:
``feature_edges``, ``hausdorff_distance``, ``edit_regions`` and
``match_periodic_nodes`` -- the Python face of each, its numpy twin pinned
against the compiled core, both CLI verbs and the pipeline steps. The C++
behaviour is pinned in ``tests/cpp/test_feature_edges.cpp`` and
``tests/cpp/test_region_analysis.cpp``."""

import json

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _feature_edges as fe
from meshioplusplus import _hausdorff as hd
from meshioplusplus import _periodic as per
from meshioplusplus._regions import Region

from .test_curvature import icosphere, open_cylinder

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")


def cube_quads():
    pts = np.array(
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
        dtype=float,
    )
    cells = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [3, 7, 6, 2],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
    ]
    return mio.Mesh(pts, [("quad", np.array(cells))])


def quad_grid(n, z=0.0):
    xs, ys = np.meshgrid(np.arange(n + 1.0), np.arange(n + 1.0))
    pts = np.column_stack([xs.ravel(), ys.ravel(), np.full(xs.size, z)])
    quads = [
        [
            j * (n + 1) + i,
            j * (n + 1) + i + 1,
            (j + 1) * (n + 1) + i + 1,
            (j + 1) * (n + 1) + i,
        ]
        for j in range(n)
        for i in range(n)
    ]
    return mio.Mesh(pts, [("quad", np.array(quads))])


def bumpy_triangles(seed=0):
    rng = np.random.default_rng(seed)
    n = 8
    xs, ys = np.meshgrid(np.arange(n), np.arange(n))
    z = rng.normal(size=xs.shape) * 0.6
    pts = np.column_stack([xs.ravel(), ys.ravel(), z.ravel()]).astype(float)
    tris = []
    for j in range(n - 1):
        for i in range(n - 1):
            a = j * n + i
            tris += [[a, a + 1, a + n + 1], [a, a + n + 1, a + n]]
    tris = np.array(tris)
    tris[3] = tris[3][[1, 0, 2]]  # one inconsistent pair
    return mio.Mesh(pts, [("triangle", tris)])


ALL = {"feature": True, "boundary": True, "non_manifold": True, "inconsistent": True}


# --------------------------------------------------------------------------- #
# feature_edges                                                                #
# --------------------------------------------------------------------------- #
def test_cube_has_twelve_feature_edges():
    out, report = mio.feature_edges(cube_quads(), return_report=True)
    assert report == {
        "num_feature": 12,
        "num_boundary": 0,
        "num_non_manifold": 0,
        "num_inconsistent": 0,
    }
    assert out.cells[0].type == "line"
    assert len(out.cells[0].data) == 12
    np.testing.assert_allclose(out.cell_data["feature:angle"][0], 90.0)
    assert len(out.points) == 8


def test_categories_can_be_switched_off():
    out = mio.feature_edges(open_cylinder(), boundary=False)
    assert not (out.cell_data["feature:kind"][0] == fe.BOUNDARY).any()
    out = mio.feature_edges(open_cylinder(), feature=False)
    assert set(out.cell_data["feature:kind"][0].tolist()) == {fe.BOUNDARY}


def test_rejects_an_angle_out_of_range():
    with pytest.raises(ValueError, match="180"):
        mio.feature_edges(cube_quads(), feature_angle=200)


@pytest.mark.parametrize(
    "mesh",
    [cube_quads(), icosphere(1), open_cylinder(), bumpy_triangles()],
    ids=["cube", "icosphere", "cylinder", "bumpy"],
)
@pytest.mark.parametrize("angle", [0.0, 20.0, 30.0, 90.0, 180.0])
@needs_core
def test_feature_edges_twin_matches_the_core(mesh, angle):
    core, core_report = mio.feature_edges(mesh, feature_angle=angle, return_report=True)
    twin, twin_report = fe._feature_edges_py(mesh, angle, ALL, "")
    assert core_report == twin_report
    np.testing.assert_array_equal(core.cells[0].data, twin.cells[0].data)
    np.testing.assert_array_equal(
        core.cell_data["feature:kind"][0], twin.cell_data["feature:kind"][0]
    )
    np.testing.assert_allclose(
        core.cell_data["feature:angle"][0],
        twin.cell_data["feature:angle"][0],
        atol=1e-10,
    )


@needs_core
def test_feature_edges_of_a_volume_mesh_use_its_skin():
    hexa = mio.Mesh(
        cube_quads().points, [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))]
    )
    _, report = mio.feature_edges(hexa, return_report=True)
    assert report["num_feature"] == 12
    twin, twin_report = fe._feature_edges_py(hexa, 30.0, ALL, "")
    assert twin_report == report


def test_decimate_and_smooth_pin_through_the_crease_kernel():
    # A once-subdivided icosphere: adjacent faces ~21 degrees apart, faces on
    # opposite sides of a vertex ~40. The old per-vertex test pinned every
    # vertex at 30 degrees; the per-edge test pins none.
    sphere = icosphere(1)
    out = mio.decimate(sphere, ratio=0.5, feature_angle=30.0)
    assert len(out.cells[0].data) < len(sphere.cells[0].data)


# --------------------------------------------------------------------------- #
# hausdorff_distance                                                           #
# --------------------------------------------------------------------------- #
def test_hausdorff_of_a_translated_copy():
    r = mio.hausdorff_distance(quad_grid(3), quad_grid(3, 0.25))
    assert r["distance"] == pytest.approx(0.25)
    assert r["mean_a_to_b"] == pytest.approx(0.25)
    assert r["num_samples_a"] == 16


@needs_core
@pytest.mark.parametrize("samples", [0, 1, 3])
def test_hausdorff_twin_matches_the_core(samples):
    a = icosphere(1)
    b = icosphere(2, radius=1.05)
    core = mio.hausdorff_distance(a, b, face_samples=samples)
    twin = hd._hausdorff_py(a, b, samples, "", "")
    for key in ("distance", "a_to_b", "b_to_a", "mean_a_to_b", "rms_b_to_a"):
        assert core[key] == pytest.approx(twin[key], abs=1e-12), key
    assert core["num_samples_a"] == twin["num_samples_a"]


def test_hausdorff_refuses_a_mesh_without_triangles():
    lines = mio.Mesh(np.zeros((2, 3)), [("line", np.array([[0, 1]]))])
    with pytest.raises(ValueError, match="no surface triangles"):
        mio.hausdorff_distance(lines, quad_grid(1))


# --------------------------------------------------------------------------- #
# edit_regions                                                                 #
# --------------------------------------------------------------------------- #
def tagged_grid():
    m = quad_grid(3)
    m.regions = [
        Region("a", "cell", [0, 1, 2, 3], dim=2, tag=1),
        Region("b", "cell", [2, 3, 4, 5], dim=2, tag=2),
        Region("s1", "side", [[0, 0], [1, 0]]),
        Region("s2", "side", [[1, 0], [1, 1]]),
        Region("wall", "point", [0], dim=0, tag=5),
        Region("wall", "cell", [0], dim=2, tag=5),
    ]
    return m


EDITS = [
    {"op": "union", "inputs": ["a", "b"], "output": "u"},
    {"op": "intersection", "inputs": ["a", "b"], "output": "i", "tag": 9},
    {"op": "difference", "inputs": ["a", "b"], "output": "d", "keep_inputs": False},
    {"op": "union", "inputs": ["s1", "s2"], "output": "s"},
    {"op": "rename", "inputs": ["s1"], "output": "S1"},
    {"op": "retag", "inputs": [{"name": "wall", "kind": "cell"}], "tag": 7, "dim": 2},
    {"op": "delete", "inputs": [{"name": "wall", "kind": "point"}]},
]


def _described(mesh):
    return [(r.name, r.kind, r.dim, r.tag, r.entries.tolist()) for r in mesh.regions]


def test_edit_regions_results():
    m = tagged_grid()
    out = mio.edit_regions(m, EDITS)
    got = {(r.name, r.kind): r for r in out.regions}
    assert got[("u", "cell")].entries.tolist() == [0, 1, 2, 3, 4, 5]
    assert got[("i", "cell")].entries.tolist() == [2, 3]
    assert got[("i", "cell")].tag == 9
    assert got[("d", "cell")].entries.tolist() == [0, 1]
    assert ("a", "cell") not in got and ("b", "cell") not in got
    assert got[("s", "side")].entries.tolist() == [[0, 0], [1, 0], [1, 1]]
    assert ("S1", "side") in got and ("s1", "side") not in got
    assert got[("wall", "cell")].tag == 7
    assert ("wall", "point") not in got
    assert len(m.regions) == 6  # the input is untouched
    np.testing.assert_array_equal(out.points, m.points)


@needs_core
def test_edit_regions_matches_the_core():
    m = tagged_grid()
    assert _described(mio.edit_regions(m, EDITS)) == _described(
        _core.edit_regions(m, EDITS)
    )


@pytest.mark.parametrize(
    "edits, match",
    [
        ([{"op": "delete", "inputs": ["wall"]}], "matches 2 regions"),
        ([{"op": "delete", "inputs": ["nope"]}], "no region matches"),
        ([{"op": "union", "inputs": ["a", "s1"], "output": "x"}], "different kinds"),
        ([{"op": "union", "inputs": ["a"], "output": "x"}], "two or more"),
        ([{"op": "rename", "inputs": ["a"], "output": "b"}], None),
        ([{"op": "xor", "inputs": ["a"]}], "unknown operation"),
        ([{"op": "retag", "inputs": ["a"]}], "tag and/or dimension"),
    ],
)
def test_edit_regions_refuses_bad_edits(edits, match):
    m = tagged_grid()
    if match is None:
        # "b" (cell, dim 2, tag 2) differs from the renamed "a" (tag 1): no
        # collision, so this one succeeds
        mio.edit_regions(m, edits)
        return
    with pytest.raises(ValueError, match=match):
        mio.edit_regions(m, edits)
    if _core is not None:
        with pytest.raises(ValueError):
            _core.edit_regions(m, edits)


def test_edit_regions_refuses_to_overwrite_an_unrelated_region():
    m = tagged_grid()
    m.regions.append(Region("x", "cell", [8], dim=2, tag=-1))  # the union's key
    with pytest.raises(ValueError, match="would replace"):
        mio.edit_regions(m, [{"op": "union", "inputs": ["a", "b"], "output": "x"}])


def test_edit_regions_keeps_an_empty_result():
    m = quad_grid(2)
    m.regions = [Region("a", "cell", [0]), Region("b", "cell", [1])]
    out = mio.edit_regions(
        m, {"op": "intersection", "inputs": ["a", "b"], "output": "e"}
    )
    assert [r.name for r in out.regions if len(r.entries) == 0] == ["e"]


# --------------------------------------------------------------------------- #
# match_periodic_nodes                                                         #
# --------------------------------------------------------------------------- #
def periodic_grid():
    m = quad_grid(3)
    m.regions = [
        Region("left", "point", [0, 4, 8, 12]),
        Region("right", "point", [3, 7, 11, 15]),
        Region("left_side", "side", [[0, 3], [3, 3], [6, 3]]),
        Region("right_cells", "cell", [2, 5, 8]),
    ]
    return m


def test_translation_pairs():
    pairs = mio.match_periodic_nodes(
        periodic_grid(), "left", "right", translate=(3, 0, 0)
    )
    assert pairs.tolist() == [[0, 3], [4, 7], [8, 11], [12, 15]]
    pairs = mio.match_periodic_nodes(
        periodic_grid(), "left_side", "right_cells", translate=(3, 0, 0)
    )
    assert pairs.tolist() == [[0, 3], [4, 7], [8, 11], [12, 15]]


def test_rotation_about_an_origin_and_fixed_points():
    # A quarter disk around (1, 1): the +x arm maps onto the +y arm under +90
    # degrees about z through (1, 1); the centre maps onto itself.
    pts = np.array(
        [[1, 1, 0], [2, 1, 0], [3, 1, 0], [1, 2, 0], [1, 3, 0], [2, 2, 0]], dtype=float
    )
    m = mio.Mesh(
        pts, [("triangle", np.array([[0, 1, 5], [0, 5, 3], [1, 2, 5], [3, 5, 4]]))]
    )
    m.regions = [Region("x", "point", [0, 1, 2]), Region("y", "point", [0, 3, 4])]
    pairs, report = mio.match_periodic_nodes(
        m, "x", "y", rotate=("z", 90.0), origin=(1, 1, 0), atol=1e-9, return_report=True
    )
    assert pairs.tolist() == [[1, 3], [2, 4]]
    assert report["num_fixed"] == 1
    twin = per._match_py(
        m,
        {"name": "x"},
        {"name": "y"},
        per.periodic_matrix(rotate=("z", 90.0), origin=(1, 1, 0)),
        1e-9,
        True,
    )
    assert twin["slave"].tolist() == [1, 2] and twin["num_fixed"] == 1


def test_incomplete_and_double_claims():
    m = periodic_grid()
    m.regions[1] = Region("right", "point", [3, 7, 11])
    with pytest.raises(ValueError, match="no master node"):
        mio.match_periodic_nodes(m, "left", "right", translate=(3, 0, 0))
    pairs, report = mio.match_periodic_nodes(
        m,
        "left",
        "right",
        translate=(3, 0, 0),
        require_complete=False,
        return_report=True,
    )
    assert report["unmatched"].tolist() == [12]
    with pytest.raises(ValueError, match="both map onto"):
        mio.match_periodic_nodes(m, "left", "right", translate=(3, 0.5, 0), atol=1.5)
    with pytest.raises(ValueError, match="tolerance"):
        mio.match_periodic_nodes(m, "left", "right", translate=(3, 0, 0), atol=0)
    with pytest.raises(ValueError, match="give a transform"):
        mio.match_periodic_nodes(m, "left", "right")


@needs_core
def test_periodic_twin_matches_the_core():
    m = periodic_grid()
    core = mio.match_periodic_nodes(m, "left", "right", translate=(3, 0, 0))
    twin = per._match_py(
        m,
        {"name": "left"},
        {"name": "right"},
        per.periodic_matrix((3, 0, 0)),
        1e-8,
        True,
    )
    assert core[:, 0].tolist() == twin["slave"].tolist()
    assert core[:, 1].tolist() == twin["master"].tolist()


# --------------------------------------------------------------------------- #
# the CLI verbs and the pipeline steps                                         #
# --------------------------------------------------------------------------- #
def test_cli_feature_edges_and_hausdorff(tmp_path, capsys):
    src = str(tmp_path / "cube.vtu")
    mio.write(src, cube_quads())
    assert (
        mio._cli.main(["feature-edges", src, str(tmp_path / "fe.vtu"), "--json"]) == 0
    )
    report = json.loads(capsys.readouterr().out)
    assert report["edges"] == 12 and report["num_feature"] == 12

    moved = cube_quads()
    moved.points = moved.points + [0.0, 0.0, 0.25]
    other = str(tmp_path / "moved.vtu")
    mio.write(other, moved)
    assert mio._cli.main(["hausdorff", src, other, "--max", "0.5"]) == 0
    capsys.readouterr()
    assert mio._cli.main(["hausdorff", src, other, "--max", "0.1", "--json"]) == 1
    report = json.loads(capsys.readouterr().out)
    assert report["passed"] is False
    assert report["distance"] == pytest.approx(0.25)


def test_cli_regions_edits_and_periodic(tmp_path, capsys):
    src = str(tmp_path / "grid.inp")
    mio.write(src, periodic_grid())
    dst = str(tmp_path / "out.inp")
    mio._cli.main(
        [
            "regions",
            src,
            dst,
            "--union",
            "both=left,right",
            "--rename",
            "left=Left",
            "--retag",
            "point:right=4",
            "--json",
        ]
    )
    listed = {r["name"]: r for r in json.loads(capsys.readouterr().out)}
    assert listed["both"]["num_entries"] == 8
    assert "Left" in listed and "left" not in listed
    assert listed["right"]["tag"] == 4
    with pytest.raises(SystemExit):
        mio._cli.main(["regions", src, "--delete", "left"])

    csv = str(tmp_path / "pairs.csv")
    mio._cli.main(
        [
            "periodic",
            src,
            "--slave",
            "left",
            "--master",
            "right",
            "--translate",
            "3,0,0",
            "--output",
            csv,
            "--json",
        ]
    )
    report = json.loads(capsys.readouterr().out)
    assert report["num_pairs"] == 4
    assert np.loadtxt(csv, delimiter=",", dtype=int).tolist() == report["pairs"]


def _settings(tmp_path, mesh, steps, ext="vtu"):
    src = tmp_path / f"in.{ext}"
    mio.write(str(src), mesh, **({"compression": None} if ext == "vtu" else {}))
    return {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": steps,
        "Output": {"Path": str(tmp_path / f"out.{ext}")},
    }


def test_pipeline_feature_edges_and_edit_regions(tmp_path):
    s = _settings(tmp_path, cube_quads(), [{"Op": "FeatureEdges", "FeatureAngle": 45}])
    report = mio.run_pipeline(s)
    assert report["steps"][0]["NumFeature"] == 12
    assert len(mio.read(s["Output"]["Path"]).cells[0].data) == 12

    s = _settings(
        tmp_path,
        periodic_grid(),
        [
            {
                "Op": "EditRegions",
                "Edit": "union",
                "Inputs": ["left", "right"],
                "Output": "lr",
            },
            {
                "Op": "EditRegions",
                "Edit": "retag",
                "Inputs": ["lr"],
                "Kind": "point",
                "Tag": 3,
            },
        ],
        ext="inp",
    )
    report = mio.run_pipeline(s)
    assert report["steps"][1]["NumRegions"] == 5
    out = mio.read(s["Output"]["Path"])
    assert any(r.name == "lr" and len(r.entries) == 8 for r in out.regions)


@needs_core
def test_pipeline_steps_match_the_c_engine(tmp_path):
    if not hasattr(_core, "run_pipeline_json"):
        pytest.skip("this build carries no C++ pipeline engine")
    s = _settings(
        tmp_path, bumpy_triangles(), [{"Op": "FeatureEdges", "FeatureAngle": 25}]
    )
    mio.run_pipeline(s)
    py = mio.read(s["Output"]["Path"])
    s["Output"]["Path"] = str(tmp_path / "out_cpp.vtu")
    try:
        _core.run_pipeline_json(json.dumps(s))
    except RuntimeError as exc:
        if "no JSON parser" in str(exc):
            pytest.skip("this build carries no JSON parser for the C++ engine")
        raise
    cpp = mio.read(s["Output"]["Path"])
    np.testing.assert_array_equal(py.cells[0].data, cpp.cells[0].data)
    np.testing.assert_array_equal(
        py.cell_data["feature:kind"][0], cpp.cell_data["feature:kind"][0]
    )
