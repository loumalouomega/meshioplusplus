"""octree_mesh: a volume mesh of a skin's inside on a balanced octree.

The castellated output's regions, distances and report through the Python
wrapper, its refusals before the core, the CLI verb, and the ``OctreeMesh``
pipeline step in both engines. The geometry itself is tested in the core
(``tests/cpp/test_octree_mesh.cpp``).
"""

import json
import os
import sys
from pathlib import Path

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
    # A closed skin: the points deep inside are flooded, not evaluated.
    assert report["num_winding_skipped"] > 0
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


# Each hexahedron corner's three edge neighbours, right-handed (VTK order).
_HEX_CORNER = np.array(
    [
        [1, 3, 4],
        [2, 0, 5],
        [3, 1, 6],
        [0, 2, 7],
        [7, 5, 0],
        [4, 6, 1],
        [5, 7, 2],
        [6, 4, 3],
    ]
)


def _corner_jacobians(mesh):
    """The scaled Jacobian of every hexahedron corner, shape (n, 8)."""
    points = np.asarray(mesh.points, dtype=float)
    hexes = np.asarray(mesh.cells[0].data, dtype=np.int64)
    out = np.empty(hexes.shape)
    for c in range(8):
        o = points[hexes[:, c]]
        a, b, d = (points[hexes[:, _HEX_CORNER[c, k]]] - o for k in range(3))
        length = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
        length *= np.linalg.norm(d, axis=1)
        triple = np.einsum("ij,ij->i", np.cross(a, b), d)
        out[:, c] = np.where(length > 0, triple / np.where(length > 0, length, 1), -1.0)
    return out


def _flat_corners(mesh):
    return int((_corner_jacobians(mesh) <= 0).sum())


def _cylinder():
    return mio.convert_cells(
        mio.extract_surface(mio.cylinder(0.5, 2.0)), mode="simplexify"
    )


@needs_core
def test_graded_hexahedral_output():
    uniform, ureport = mio.octree_mesh(
        _cylinder(), output="hex", min_size=0.05, return_report=True
    )
    out, report = mio.octree_mesh(
        _cylinder(),
        output="hex",
        min_size=0.05,
        hex_grading="graded",
        return_report=True,
    )
    assert [c.type for c in out.cells] == ["hexahedron"]
    assert report["num_hexahedra"] < ureport["num_hexahedra"]
    assert report["num_inverted"] == 0
    # The rims are chains the staircase follows by turning inside
    # quadrilaterals; the fit releases those nodes rather than flatten them.
    assert ureport["min_scaled_jacobian"] >= 0.2
    assert report["min_scaled_jacobian"] >= 0.2
    assert _flat_corners(uniform) == _flat_corners(out) == 0
    assert ureport["warnings"] == report["warnings"] == []
    assert report["num_boundary_edges_open"] == 0
    assert report["num_non_manifold_boundary_edges"] == 0
    assert ureport["hausdorff_to_skin"] < 0.01
    assert report["hausdorff_to_skin"] < 0.01
    levels = np.asarray(out.cell_data["octree:level"][0])
    assert len(np.unique(levels)) > 1
    _, conformity = mio.check_conformity(out, return_report=True)
    assert conformity["num_hanging_nodes"] == 0
    with pytest.raises(ValueError, match="hex_grading"):
        mio.octree_mesh(_cylinder(), output="hex", hex_grading="coarse")


@needs_core
@pytest.mark.parametrize("grading", ["uniform", "graded"])
@pytest.mark.parametrize("off", [{"feature_angle": 0}, {"preserve_features": False}])
def test_hexahedral_box_without_features_has_no_flat_corner(grading, off):
    # Near a box's edges and corners a graded boundary node is about as near
    # two or three faces; projected to its closest point, a face's pillow
    # copies went to the side face and flattened the pillow cells.
    box = mio.convert_cells(
        mio.extract_surface(mio.box((4, 4, 4), (1.0, 1.0, 1.0))), mode="simplexify"
    )
    out, report = mio.octree_mesh(
        box, output="hex", min_size=0.05, hex_grading=grading, return_report=True, **off
    )
    assert report["num_corners"] == 0
    assert _flat_corners(out) == 0
    assert report["min_scaled_jacobian"] >= 0.2
    assert report["warnings"] == []


@needs_core
@pytest.mark.parametrize("grading", ["uniform", "graded"])
def test_hexahedral_cleaving_gives_the_tetrahedral_materials(grading):
    # Nested shells: with conformity="cleaving" the hexahedra carry one
    # material per shell, named as the tetrahedral output names them.
    skin = mio.merge([mio.icosa_sphere(1.0, 4), mio.icosa_sphere(0.5, 4)])
    size = {"surface_size": 0.1} if grading == "uniform" else {}
    hexes, report = mio.octree_mesh(
        skin,
        output="hex",
        conformity="cleaving",
        hex_grading=grading,
        min_size=0.1,
        return_report=True,
        **size,
    )
    tets = mio.octree_mesh(skin, conformity="cleaving", min_size=0.1, max_size=0.4)
    material = np.asarray(hexes.cell_data["octree_mesh:material"][0])
    assert sorted(np.unique(material).tolist()) == [1, 2]
    names = {r.name for r in hexes.regions}
    assert {"material_1", "material_2"} <= names
    assert {n for n in names if n.startswith("material_")} == {
        r.name for r in tets.regions if r.name.startswith("material_")
    }
    assert report["min_scaled_jacobian"] >= 0.2
    assert report["num_inverted"] == 0
    assert report["warnings"] == []
    _, conformity = mio.check_conformity(hexes, return_report=True)
    assert conformity["num_hanging_nodes"] == 0
    # Without cleaving the same skin is one material, as before.
    plain = mio.octree_mesh(
        skin, output="hex", hex_grading=grading, min_size=0.1, **size
    )
    assert "octree_mesh:material" not in plain.cell_data


_BUNNY = Path(__file__).resolve().parents[2] / "example" / "data" / "Bunny.stl"


@needs_core
@pytest.mark.skipif(
    os.environ.get("MESHIOPLUSPLUS_OCTREE_MESH_BUNNY") != "1" or not _BUNNY.exists(),
    reason="the Stanford bunny's hexahedra take a minute; set MESHIOPLUSPLUS_OCTREE_MESH_BUNNY=1",
)
@pytest.mark.parametrize("grading", ["uniform", "graded"])
def test_hexahedral_bunny(grading):
    # The bunny of doc/octree_mesh.md's hexahedral table.
    out, report = mio.octree_mesh(
        mio.read(_BUNNY),
        output="hex",
        min_size=1.5,
        max_size=8.0,
        feature_angle=60,
        hex_grading=grading,
        return_report=True,
    )
    assert report["num_inverted"] == 0
    assert report["num_boundary_edges_open"] == 0
    if grading == "graded":
        assert report["min_scaled_jacobian"] >= 0.2
        assert _flat_corners(out) == 0
    else:
        # Open (doc/octree_mesh.md): a few pillow cells at the ear tips keep a
        # corner just below zero, which the warnings report.
        assert report["min_scaled_jacobian"] > -0.01
        assert any("non-positive corner" in w for w in report["warnings"])


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
    code = mio._cli.main(
        [
            "octree-mesh",
            str(skin),
            str(out),
            "--output",
            "hex",
            "--max-size",
            "0.5",
            "--hex-grading",
            "graded",
            "--json",
        ]
    )
    assert code == 0
    report = json.loads(capsys.readouterr().out)
    assert report["num_hexahedra"] == len(mio.read(out).cells[0].data) > 0


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
                "HexGrading": "graded",
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
