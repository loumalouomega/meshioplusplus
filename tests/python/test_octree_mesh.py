"""octree_mesh: a volume mesh of a skin's inside on a balanced octree.

The castellated output's regions, distances and report through the Python
wrapper, its refusals before the core, the CLI verb, and the ``OctreeMesh``
pipeline step in both engines. The geometry itself is tested in the core
(``tests/cpp/test_octree_mesh.cpp``).
"""

import json
import sys

import numpy as np
import pytest

import meshioplusplus as mio

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - the numpy-only install
    _core = None

needs_core = pytest.mark.skipif(
    not hasattr(_core, "octree_mesh"), reason="core predates octree_mesh"
)


def _sphere():
    return mio.icosa_sphere(1.0, 3)


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"output": "quad"}, "unknown output"),
        ({"fill": "half"}, "unknown fill"),
        ({"skin_grading": "fancy"}, "unknown skin_grading"),
        ({"conformity": "magic"}, "unknown conformity"),
        ({"min_size": -1.0}, "min_size"),
        ({"min_size": 1.0, "max_size": 0.5}, "min_size must not exceed"),
        ({"max_level": 20}, "max_level"),
        ({"winding_threshold": 1.0}, "winding_threshold"),
        ({"alpha_long": 0.6}, "alpha_long"),
        ({"grading_layers": 9}, "grading_layers"),
        ({"bounds": (0, 0, 0, 1)}, "six entries"),
        ({"max_cells": 0}, "max_cells"),
    ],
)
def test_bad_arguments_are_refused_before_the_core(kwargs, match):
    with pytest.raises(ValueError, match=match):
        mio.octree_mesh(_sphere(), **kwargs)


def test_without_core_the_refusal_is_named(monkeypatch):
    monkeypatch.setitem(sys.modules, "meshioplusplus._core", None)
    monkeypatch.delattr(mio, "_core", raising=False)
    with pytest.raises(NotImplementedError, match="no numpy twin"):
        mio.octree_mesh(_sphere(), output="castellated")


@needs_core
def test_castellated_regions_distances_and_report():
    out, report = mio.octree_mesh(
        _sphere(),
        output="castellated",
        max_size=0.5,
        min_size=0.0625,
        return_report=True,
    )
    assert out.cells[0].type == "hexahedron"
    names = {r.name: r for r in out.regions}
    assert {"inside", "boundary_cells", "surrogate_boundary"} <= set(names)
    assert names["surrogate_boundary"].kind == "side"
    assert report["output"] == "castellated"
    assert report["num_inside"] == report["num_cells"] == len(out.cells[0].data)
    assert report["num_surrogate_faces"] == len(names["surrogate_boundary"].entries)
    assert report["num_inverted"] == 0
    assert report["min_scaled_jacobian"] == pytest.approx(1.0)
    assert sum(report["leaves_per_level"]) == report["num_leaves"]
    assert report["quality"]["num_cells"] == report["num_cells"]
    assert "quality:scaled_jacobian" in report["quality"]["metrics"]
    assert list(report["timings"]) == [
        "skin",
        "size_field",
        "octree",
        "output",
        "quality",
    ]
    assert report["warnings"] == []
    # The distance is negative inside, the closest points lie on the skin.
    d = out.point_data["sdf:distance"]
    cp = out.point_data["sdf:closest_point"]
    assert cp.shape == (len(out.points), 3)
    r = np.linalg.norm(out.points, axis=1)
    assert np.all(d[r < 0.8] < 0)
    assert np.all(np.abs(np.linalg.norm(cp, axis=1) - 1.0) < 0.02)
    # Every inside cell's centre is inside the unit sphere's polyhedron.
    centres = out.points[out.cells[0].data].mean(axis=1)
    assert np.all(np.linalg.norm(centres, axis=1) < 1.0)


@needs_core
def test_fills_nest():
    common = dict(output="castellated", max_size=0.5, min_size=0.125)
    inside = mio.octree_mesh(_sphere(), fill="inside", **common)
    cover = mio.octree_mesh(_sphere(), fill="cover", **common)
    every = mio.octree_mesh(_sphere(), fill="all", **common)
    assert (
        len(inside.cells[0].data) < len(cover.cells[0].data) < len(every.cells[0].data)
    )


@needs_core
def test_hexahedral_output():
    out, report = mio.octree_mesh(
        _sphere(), output="hex", surface_size=0.15, return_report=True
    )
    assert [c.type for c in out.cells] == ["hexahedron"]
    assert report["output"] == "hex"
    assert report["num_hexahedra"] == len(out.cells[0].data)
    assert report["num_inverted"] == 0
    assert report["min_scaled_jacobian"] >= 0.2
    assert report["num_boundary_edges_open"] == 0
    assert report["num_non_manifold_boundary_edges"] == 0
    assert 0 < report["hausdorff_to_skin"] < 0.05
    assert report["warnings"] == []
    pillow = np.asarray(out.cell_data["octree_mesh:pillow"][0])
    names = {r.name: r for r in out.regions}
    sides = np.asarray(names["boundary"].entries)
    # The boundary is the pillow layer's outer faces.
    assert np.all(pillow[sides[:, 0]] == 1) and np.all(sides[:, 1] == 5)
    assert "sdf:distance" in out.point_data
    # Without the pillow the grid's own cells meet the skin.
    bare = mio.octree_mesh(_sphere(), output="hex", surface_size=0.15, pillow=False)
    assert not np.any(np.asarray(bare.cell_data["octree_mesh:pillow"][0]))


@needs_core
def test_tetrahedral_output():
    out, report = mio.octree_mesh(
        _sphere(), max_size=0.5, min_size=0.0625, return_report=True
    )
    assert out.cells[0].type == "tetra"
    assert report["output"] == "tet"
    assert report["num_tetrahedra"] == len(out.cells[0].data)
    assert report["num_inverted"] == 0
    assert report["min_dihedral"] > 10.0
    names = {r.name: r for r in out.regions}
    assert names["boundary"].kind == "side"
    assert len(names["inside"].entries) == len(out.cells[0].data)
    assert "sdf:distance" in out.point_data
    # The boundary nodes lie on the skin.
    assert (
        np.max(
            np.abs(
                out.point_data["sdf:distance"][
                    np.unique(out.cells[0].data[names["boundary"].entries[:, 0]])
                ]
            )
        )
        < 0.6
    )
    p = out.points
    t = out.cells[0].data
    vol = np.einsum(
        "ij,ij->i",
        p[t[:, 1]] - p[t[:, 0]],
        np.cross(p[t[:, 2]] - p[t[:, 0]], p[t[:, 3]] - p[t[:, 0]]),
    )
    assert np.all(vol > 0)
    assert vol.sum() / 6 == pytest.approx(4 / 3 * np.pi, rel=0.05)
    # The quality gate: a watertight boundary near the skin, and the
    # optimisation's counts.
    assert report["num_boundary_edges_open"] == 0
    assert report["num_non_manifold_boundary_edges"] == 0
    assert 0 < report["hausdorff_to_skin"] < 0.1
    assert report["num_optimized_moves"] > 0
    assert report["num_flips"] >= 0
    _, raw = mio.octree_mesh(
        _sphere(), max_size=0.5, min_size=0.0625, fit_iterations=0, return_report=True
    )
    assert raw["num_optimized_moves"] == 0
    assert report["min_dihedral"] >= raw["min_dihedral"] - 1e-9
    assert report["max_dihedral"] <= raw["max_dihedral"] + 1e-9


@needs_core
def test_cube_corners_are_snapped():
    lo, size = np.array([-0.61, -0.67, -0.7]), 1.37
    cube = mio.convert_cells(
        mio.extract_surface(mio.box([8, 8, 8], [size] * 3, list(lo))),
        mode="simplexify",
    )
    out, report = mio.octree_mesh(cube, max_size=0.5, min_size=0.05, return_report=True)
    assert report["num_corners"] == report["num_corners_snapped"] == 8
    assert report["num_feature_chains_missed"] == 0
    assert report["guarantee_applies"] is False
    corners = np.array(
        [
            lo + size * np.array([i, j, k])
            for i in (0, 1)
            for j in (0, 1)
            for k in (0, 1)
        ]
    )
    d = np.linalg.norm(out.points[None, :, :] - corners[:, None, :], axis=2).min(axis=1)
    assert np.all(d < 1e-12)
    names = {r.name: r for r in out.regions}
    assert len(names["feature_corners"].entries) == 8


@needs_core
def test_cli_round_trip(tmp_path, capsys):
    skin = tmp_path / "skin.stl"
    mio.write(skin, _sphere())
    out = tmp_path / "volume.vtu"
    code = mio._cli.main(
        [
            "octree-mesh",
            str(skin),
            str(out),
            "--output",
            "castellated",
            "--max-size",
            "0.5",
            "--min-size",
            "0.125",
            "--json",
        ]
    )
    assert code == 0
    report = json.loads(capsys.readouterr().out)
    assert "timings" not in report
    assert report["num_cells"] == len(mio.read(out).cells[0].data)


@needs_core
def test_pipeline_step_in_both_engines(tmp_path):
    src = tmp_path / "skin.vtu"
    mio.write(src, _sphere(), compression=None)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {
                "Op": "OctreeMesh",
                "Output": "castellated",
                "MaxSize": 0.5,
                "MinSize": 0.125,
                "Fill": "cover",
            }
        ],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    reports = [mio.run_pipeline(settings)]
    outputs = [mio.read(tmp_path / "py.vtu")]
    if hasattr(_core, "run_pipeline_json"):
        settings["Output"]["Path"] = str(tmp_path / "cpp.vtu")
        reports.append(_core.run_pipeline_json(json.dumps(settings)))
        outputs.append(mio.read(tmp_path / "cpp.vtu"))
    for report in reports:
        entry = report["steps"][-1]
        assert entry["op"] == "OctreeMesh"
        for key in (
            "NumLeaves",
            "NumCells",
            "NumInside",
            "NumSurrogateFaces",
            "NumHanging",
            "NumInverted",
            "MinDihedral",
            "MinScaledJacobian",
        ):
            assert key in entry, key
        assert entry["NumInside"] < entry["NumCells"]
    if len(reports) == 2:
        py, cpp = (r["steps"][-1] for r in reports)
        assert {k: py[k] for k in py if k[:1].isupper()} == {
            k: cpp[k] for k in cpp if k[:1].isupper()
        }
        np.testing.assert_array_equal(outputs[0].points, outputs[1].points)
        np.testing.assert_array_equal(
            outputs[0].cells[0].data, outputs[1].cells[0].data
        )
