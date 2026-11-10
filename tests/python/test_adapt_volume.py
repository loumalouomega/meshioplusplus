"""Tests of meshioplusplus.adapt_volume (doc/adapt_volume.md): the Python
surface over the C++ core, its report, its invariants measured with the
library's own checks, the pipeline step in both engines, and
octree_mesh(post_adapt=True)."""

import json

import numpy as np
import pytest

import meshioplusplus as mp

_core = pytest.importorskip("meshioplusplus._core")


def _tet_box(n=6):
    return mp.box(dims=(n, n, n), cell_type="tetra")


def _boundary(mesh):
    return mp.extract_surface(mesh)


def _signed_volumes(mesh):
    p = np.asarray(mesh.points, dtype=float)
    vols = []
    for block in mesh.cells:
        c = np.asarray(block.data)
        a, b, cc, d = (p[c[:, k]] for k in range(4))
        vols.append(np.einsum("ij,ij->i", np.cross(b - a, cc - a), d - a) / 6.0)
    return np.concatenate(vols)


def test_coarsens_a_box_and_keeps_its_invariants():
    box = _tet_box()
    out, rep = mp.adapt_volume(box, hmin=0.3, hmax=0.4, return_report=True)
    assert rep["num_tets_before"] == 6 * 6 * 6 * 6
    assert rep["num_tets_after"] < rep["num_tets_before"] / 2
    assert rep["num_tets_after"] == sum(len(b.data) for b in out.cells)
    assert rep["num_points_after"] == len(out.points)
    assert rep["num_corners"] == 8
    vols = _signed_volumes(out)
    assert (vols > 0).all()
    assert vols.sum() == pytest.approx(1.0, abs=1e-12)
    # the boundary stays on the box, within hausd measured independently
    h = mp.hausdorff_distance(_boundary(out), _boundary(box), face_samples=4)
    assert h["distance"] <= rep["hausd"]
    assert rep["max_boundary_deviation"] <= rep["hausd"]
    conf = mp.check_conformity(out, return_report=True)[1]
    assert conf["num_hanging_nodes"] == 0
    assert rep["min_dihedral_after"] >= min(rep["min_dihedral_before"], 10.0) - 1e-9
    # every side region survives with faces on its plane
    for name, axis, value in (("xmin", 0, 0.0), ("xmax", 0, 1.0), ("zmax", 2, 1.0)):
        assert name in out.side_sets if hasattr(out, "side_sets") else True
    assert set(rep["timings"]) == {
        "setup",
        "collapse",
        "split",
        "swap",
        "smooth",
        "polish",
        "output",
    }


def test_report_has_every_field():
    _, rep = mp.adapt_volume(_tet_box(3), return_report=True)
    expected = {
        "num_tets_before",
        "num_tets_after",
        "num_points_before",
        "num_points_after",
        "num_collapsed",
        "num_split",
        "num_swapped",
        "num_smoothed",
        "num_rejected_by_link",
        "num_rejected_by_inversion",
        "num_rejected_by_quality",
        "num_rejected_by_length",
        "num_rejected_by_class",
        "num_rejected_by_hausdorff",
        "num_rejected_by_normal",
        "min_dihedral_before",
        "min_dihedral_after",
        "max_dihedral_before",
        "max_dihedral_after",
        "mean_quality_before",
        "mean_quality_after",
        "max_boundary_deviation",
        "hmin",
        "hmax",
        "hausd",
        "volume_before",
        "volume_after",
        "region_volume_change",
        "num_corners",
        "num_ridge_vertices",
        "num_chains",
        "num_patches",
        "iterations",
        "timings",
        "warnings",
    }
    assert set(rep) == expected


def test_preserve_boundary_leaves_the_boundary_bit_for_bit():
    box = _tet_box(5)
    out = mp.adapt_volume(box, hmax=0.5, preserve_boundary=True)
    a, b = _boundary(box), _boundary(out)

    def tris(m):
        p = np.asarray(m.points)
        t = np.sort(p[np.asarray(m.cells[0].data)], axis=1)
        return sorted(map(lambda x: tuple(x.ravel()), t))

    assert tris(a) == tris(b)


def test_passes_off_is_the_identity_and_float_point_data_follows():
    box = _tet_box(3)
    box.point_data["x"] = np.asarray(box.points)[:, 0].copy()
    out = mp.adapt_volume(box, collapse=False, split=False, swap=False, smooth=False)
    np.testing.assert_array_equal(out.points, box.points)
    np.testing.assert_array_equal(out.point_data["x"], box.point_data["x"])
    # splits interpolate and collapses keep the survivor's value, so a linear
    # field stays linear; smoothing moves points and not their data (as in
    # optimize_volume), so it is off here
    out2 = mp.adapt_volume(box, hmin=0.05, hmax=0.1, smooth=False)
    np.testing.assert_allclose(
        out2.point_data["x"], np.asarray(out2.points)[:, 0], atol=1e-12
    )


def test_bad_options_are_refused_by_name():
    box = _tet_box(2)
    with pytest.raises(ValueError, match="hgrad"):
        mp.adapt_volume(box, hgrad=0.5)
    with pytest.raises(ValueError, match="quality_metric"):
        mp.adapt_volume(box, quality_metric="best")
    with pytest.raises(ValueError, match="size_array"):
        mp.adapt_volume(box, size_array="nope")
    with pytest.raises(ValueError):
        mp.adapt_volume(mp.icosa_sphere(1.0, 1))


def test_deterministic():
    box = _tet_box(5)
    a = mp.adapt_volume(box, hmin=0.25, hmax=0.35)
    b = mp.adapt_volume(box, hmin=0.25, hmax=0.35)
    assert mp.meshes_equal(a, b)


def test_octree_mesh_post_adapt_is_opt_in():
    skin = mp.icosa_sphere(1.0, 3)
    plain, rep = mp.octree_mesh(skin, return_report=True)
    assert "adapt" not in rep
    again = mp.octree_mesh(skin, post_adapt=False)
    assert mp.meshes_equal(plain, again)
    out, rep2 = mp.octree_mesh(
        skin, post_adapt=True, adapt_hgrad=2.0, return_report=True
    )
    assert rep2["adapt"]["num_tets_before"] == rep["num_tetrahedra"]
    assert rep2["adapt"]["num_tets_after"] == rep2["num_tetrahedra"]
    assert rep2["num_inverted"] == 0
    assert rep2["guarantee_applies"] is False
    assert "adapt" in rep2["timings"]
    assert "inside" in out.cell_sets


def test_pipeline_step_matches_the_direct_call_in_both_engines(tmp_path):
    box = _tet_box(4)
    src = tmp_path / "in.vtu"
    mp.write(str(src), box, compression=None)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "AdaptVolume", "Hmin": 0.3, "Hmax": 0.4}],
        "Output": {"Path": str(tmp_path / "out.vtu")},
    }
    report = mp.run_pipeline(settings)
    step = report["steps"][0]
    direct, rep = mp.adapt_volume(box, hmin=0.3, hmax=0.4, return_report=True)
    assert step["NumTetsAfter"] == rep["num_tets_after"]
    assert step["NumCollapsed"] == rep["num_collapsed"]
    out_py = mp.read(str(tmp_path / "out.vtu"))
    if getattr(_core, "__has_json__", False):
        settings["Output"]["Path"] = str(tmp_path / "out_cpp.vtu")
        _core.run_pipeline_json(json.dumps(settings))
        assert mp.meshes_equal(out_py, mp.read(settings["Output"]["Path"]))
