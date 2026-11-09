"""remove_islands, flip_normals and clip (roadmap 7.8.1, 7.8.2, 7.8.3)."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "remove_islands"), reason="needs the C++ core"
)


def tri_volume(mesh):
    p = mesh.points
    t = np.concatenate([c.data for c in mesh.cells if c.type == "triangle"])
    a, b, c = p[t[:, 0]], p[t[:, 1]], p[t[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def sphere_with_crumbs(n=10):
    big = mp.icosa_sphere(1.0, 3)
    pts = [big.points]
    tris = [big.cells[0].data]
    off = len(big.points)
    for k in range(n):
        small = mp.icosa_sphere(0.01, 1, center=(3.0 + 0.1 * k, 0.0, 0.0))
        pts.append(small.points)
        tris.append(small.cells[0].data + off)
        off += len(small.points)
    return mp.Mesh(np.vstack(pts), [("triangle", np.vstack(tris))]), big


def test_remove_islands_keeps_the_sphere():
    mesh, big = sphere_with_crumbs()
    out, rep = mp.remove_islands(mesh, return_report=True)
    assert rep["num_components"] == 11
    assert rep["num_removed"] == 10
    assert rep["component_removed"][0] is False
    assert len(out.points) == len(big.points)
    np.testing.assert_array_equal(out.points, big.points)
    np.testing.assert_array_equal(out.cells[0].data, big.cells[0].data)


def test_remove_islands_dry_run_and_thresholds():
    mesh, _ = sphere_with_crumbs(3)
    lab = mp.remove_islands(mesh, dry_run=True)
    assert len(lab.points) == len(mesh.points)
    removed = lab.cell_data["islands:removed"][0]
    comp = lab.cell_data["islands:component"][0]
    assert removed.sum() == 3 * 80
    assert set(comp.tolist()) == {0, 1, 2, 3}
    out = mp.remove_islands(mesh, min_diagonal=0.0, min_cells=81)
    assert len(out.cells[0].data) == len(mesh.cells[0].data) - 3 * 80
    with pytest.raises(ValueError):
        mp.remove_islands(mesh, min_diagonal=-1.0)
    with pytest.raises(ValueError):
        mp.remove_islands(mesh, connectivity="face")


def test_flip_normals_twice_is_the_input():
    s = mp.icosa_sphere(1.0, 2)
    once, rep = mp.flip_normals(s, return_report=True)
    assert rep["num_flipped"] == len(s.cells[0].data)
    assert tri_volume(once) == pytest.approx(-tri_volume(s))
    np.testing.assert_array_equal(once.points, s.points)
    twice = mp.flip_normals(once)
    np.testing.assert_array_equal(twice.cells[0].data, s.cells[0].data)
    with pytest.raises(ValueError):
        mp.flip_normals(s, "nowhere")


def test_flip_normals_in_a_region():
    s = mp.icosa_sphere(1.0, 1)
    s.regions = [mp.Region("cap", "cell", np.array([0, 1, 2]))]
    out, rep = mp.flip_normals(s, "cap", return_report=True)
    assert rep["num_flipped"] == 3
    np.testing.assert_array_equal(out.cells[0].data[3:], s.cells[0].data[3:])
    assert not np.array_equal(out.cells[0].data[:3], s.cells[0].data[:3])


def test_clip_a_sphere_and_cap_it():
    s = mp.icosa_sphere(1.0, 4)
    below, rep = mp.clip(s, (0, 0, 0), (0, 0, 1), return_report=True)
    assert rep["num_loops"] == 1
    assert rep["num_open_chains"] == 0
    assert rep["num_cap_triangles"] > 0
    assert below.points[:, 2].max() <= 1e-12
    above = mp.clip(s, (0, 0, 0), (0, 0, 1), keep="above")
    assert above.points[:, 2].min() >= -1e-12
    v = tri_volume(s)
    assert tri_volume(below) + tri_volume(above) == pytest.approx(v, rel=1e-12)
    assert tri_volume(below) == pytest.approx(v / 2, rel=1e-9)
    assert "clip:cap" in below.cell_data
    open_part = mp.clip(s, (0, 0, 0), (0, 0, 1), cap=False)
    assert (
        len(open_part.cells[0].data)
        == len(below.cells[0].data) - rep["num_cap_triangles"]
    )
    with pytest.raises(ValueError):
        mp.clip(s, (0, 0, 0), (0, 0, 0))
    with pytest.raises(ValueError):
        mp.clip(s, (0, 0, 0), (0, 0, 1), keep="left")


def test_both_pipeline_engines(tmp_path):
    mesh, _ = sphere_with_crumbs(2)
    src = tmp_path / "scan.vtu"
    mp.write(src, mesh)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {"Op": "RemoveIslands"},
            {"Op": "FlipNormals"},
            {"Op": "Clip", "Origin": [0, 0, 0], "Normal": [0, 0, 1]},
        ],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    mp.run_pipeline(settings)
    settings["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    _core.run_pipeline_json(__import__("json").dumps(settings))
    a, b = mp.read(str(tmp_path / "py.vtu")), mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(a.points, b.points)
    np.testing.assert_array_equal(a.cells[0].data, b.cells[0].data)
    assert a.points[:, 0].max() < 2.0


def test_the_mcp_tools(tmp_path):
    from meshioplusplus.mcp import _tools

    mesh, _ = sphere_with_crumbs(2)
    src = tmp_path / "scan.vtu"
    mp.write(src, mesh)
    rep = _tools.tool_remove_islands(str(src), str(tmp_path / "a.vtu"))
    assert rep["num_removed"] == 2
    rep = _tools.tool_flip_normals(str(tmp_path / "a.vtu"), str(tmp_path / "b.vtu"))
    assert rep["num_flipped"] > 0
    rep = _tools.tool_clip(
        str(tmp_path / "b.vtu"), str(tmp_path / "c.vtu"), [0, 0, 0], [0, 0, 1]
    )
    assert rep["num_loops"] == 1


def quad_cube():
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
        + [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]],
        dtype=float,
    )
    quads = np.array(
        [
            [0, 3, 2, 1],
            [4, 5, 6, 7],
            [0, 4, 7, 3],
            [1, 2, 6, 5],
            [0, 1, 5, 4],
            [3, 7, 6, 2],
        ]
    )
    return mp.Mesh(pts, [("quad", quads)])


def shell_volume(mesh):
    p = mesh.points
    tris = []
    for c in mesh.cells:
        d = c.data
        if c.type == "triangle":
            tris.append(d)
        else:
            tris.append(d[:, [0, 1, 2]])
            tris.append(d[:, [0, 2, 3]])
    t = np.concatenate(tris)
    a, b, c = p[t[:, 0]], p[t[:, 1]], p[t[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def test_offset_a_cube_grows_it_by_the_distance():
    out, rep = mp.offset_surface(quad_cube(), 0.1, return_report=True)
    np.testing.assert_allclose(np.abs(out.points - 0.5), 0.6, atol=1e-12)
    assert rep["max_move"] == pytest.approx(0.1 * np.sqrt(3.0), abs=1e-12)
    assert rep["num_folded"] == 0 and rep["warnings"] == []
    inward = mp.offset_surface(quad_cube(), -0.1)
    np.testing.assert_allclose(np.abs(inward.points - 0.5), 0.4, atol=1e-12)
    with pytest.raises(ValueError, match="finite"):
        mp.offset_surface(quad_cube(), float("nan"))


def test_thicken_a_hemisphere_and_hollow_a_cube():
    hemi = mp.clip(
        mp.icosa_sphere(1.0, 4), (0, 0, 0), (0, 0, 1), keep="above", cap=False
    )
    for side, hi in (("inward", 1.0), ("outward", 1.05), ("both", 1.025)):
        shell, rep = mp.thicken(hemi, 0.05, side=side, return_report=True)
        assert rep["num_rim_quads"] > 0
        want = 2 * np.pi / 3 * (hi**3 - (hi - 0.05) ** 3)
        assert shell_volume(shell) == pytest.approx(want, rel=0.02)
        assert set(np.concatenate(shell.cell_data["thicken:side"]).tolist()) == {
            0,
            1,
            2,
        }
    with pytest.raises(ValueError, match="side"):
        mp.thicken(hemi, 0.05, side="sideways")
    hollowed = mp.hollow(quad_cube(), 0.1)
    assert shell_volume(hollowed) == pytest.approx(1.0 - 0.8**3, abs=1e-12)
    with pytest.raises(ValueError, match="boundary"):
        mp.hollow(hemi, 0.05)


def test_offset_pipeline_engines_agree(tmp_path):
    src = tmp_path / "cube.vtu"
    mp.write(src, quad_cube())
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {"Op": "OffsetSurface", "Distance": 0.1},
            {"Op": "Hollow", "Wall": 0.2},
        ],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep = mp.run_pipeline(settings)
    assert rep["steps"][0]["MaxMove"] == pytest.approx(0.1 * np.sqrt(3.0))
    settings["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    _core.run_pipeline_json(__import__("json").dumps(settings))
    a, b = mp.read(str(tmp_path / "py.vtu")), mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(a.points, b.points)
    settings["Operations"] = [{"Op": "Thicken"}]
    with pytest.raises(Exception, match="Thickness"):
        mp.run_pipeline(settings)


def test_offset_mcp_tools(tmp_path):
    from meshioplusplus.mcp import _tools

    src = tmp_path / "cube.vtu"
    mp.write(src, quad_cube())
    assert _tools.tool_offset_surface(str(src), str(tmp_path / "o.vtu"), 0.1)[
        "max_move"
    ] == pytest.approx(0.1 * np.sqrt(3.0))
    assert _tools.tool_hollow(str(src), str(tmp_path / "h.vtu"), 0.1)["num_folded"] == 0
    assert (
        _tools.tool_thicken(str(src), str(tmp_path / "t.vtu"), 0.1)["num_rim_quads"]
        == 0
    )


def test_decimate_max_deviation_holds_on_the_coverage():
    s = mp.icosa_sphere(1.0, 4)
    out, rep = mp.decimate(s, max_deviation=0.01, return_report=True)
    assert rep["deviation_met"] is True
    assert rep["deviation"] <= 0.01
    assert rep["rounds"] > 1
    assert len(out.cells[0].data) < len(s.cells[0].data) // 2
    h = mp.hausdorff_distance(out, s, face_samples=2)
    assert h["distance"] < 0.03
    capped, rep = mp.decimate(s, max_deviation=1e-6, max_faces=200, return_report=True)
    assert rep["deviation_met"] is False
    assert len(capped.cells[0].data) <= 202
    alone = mp.decimate(s, max_faces=300)
    assert len(alone.cells[0].data) <= 302
    with pytest.raises(ValueError, match="max_deviation"):
        mp.decimate(s, ratio=0.5, max_deviation=0.01)


def test_decimate_bounds_in_both_pipeline_engines(tmp_path):
    src = tmp_path / "s.vtu"
    mp.write(src, mp.icosa_sphere(1.0, 3))
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "Decimate", "MaxDeviation": 0.01}],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep = mp.run_pipeline(settings)
    assert rep["steps"][0]["DeviationMet"] == 1.0
    settings["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    cpp = _core.run_pipeline_json(__import__("json").dumps(settings))
    assert cpp["steps"][0]["Rounds"] == rep["steps"][0]["Rounds"]
    a, b = mp.read(str(tmp_path / "py.vtu")), mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(a.points, b.points)
