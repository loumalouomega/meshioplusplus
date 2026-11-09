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
