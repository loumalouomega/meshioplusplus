"""A UV atlas of a triangle surface (roadmap 7.3.3): the operation's probes, the
texture coordinates OBJ, glTF and USD write and read in both engines, and the
pipeline step in both engines."""

import importlib.util
import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core

pytestmark = pytest.mark.skipif(
    not hasattr(_core, "parameterize"), reason="needs the C++ core"
)


def _signed_areas(mesh):
    uv = mesh.point_data["uv"]
    t = mesh.cells[0].data
    d1 = uv[t[:, 1]] - uv[t[:, 0]]
    d2 = uv[t[:, 2]] - uv[t[:, 0]]
    return 0.5 * (d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0])


def _torus(n=24, m=12):
    u, v = np.meshgrid(
        np.linspace(0, 2 * np.pi, n, endpoint=False),
        np.linspace(0, 2 * np.pi, m, endpoint=False),
        indexing="ij",
    )
    pts = np.column_stack(
        [
            ((2 + 0.6 * np.cos(v)) * np.cos(u)).ravel(),
            ((2 + 0.6 * np.cos(v)) * np.sin(u)).ravel(),
            (0.6 * np.sin(v)).ravel(),
        ]
    )
    tris = []
    for i in range(n):
        for j in range(m):
            a, b = i * m + j, ((i + 1) % n) * m + j
            c, d = ((i + 1) % n) * m + (j + 1) % m, i * m + (j + 1) % m
            tris += [[a, b, c], [a, c, d]]
    return mp.Mesh(pts, [("triangle", np.array(tris))])


def test_a_disc_is_one_chart_with_no_flipped_triangle():
    disc = mp.disk(1.0, 8, cell_type="triangle")
    out, report = mp.parameterize(disc, return_report=True)
    assert report["num_charts"] == 1 and report["num_flipped"] == 0
    assert report["num_seam_points"] == 0 and len(out.points) == len(disc.points)
    assert report["stretch"] == pytest.approx(1.0, abs=1e-9)  # flat: isometric
    assert np.all(_signed_areas(out) > 0)
    uv = out.point_data["uv"]
    assert uv.shape == (len(out.points), 2)
    assert uv.min() >= -1e-12 and uv.max() <= 1 + 1e-12


@pytest.mark.parametrize("name", ["sphere", "torus"])
def test_closed_surfaces_get_atlases_without_flips(name):
    surface = mp.icosa_sphere(subdivisions=3) if name == "sphere" else _torus()
    out, report = mp.parameterize(surface, return_report=True)
    assert report["num_charts"] >= 2 and report["num_flipped"] == 0
    assert np.all(_signed_areas(out) > 0)
    assert 1.0 <= report["stretch"] < 2.0
    parent = out.point_data["parameterize:parent_point"]
    np.testing.assert_array_equal(
        parent[: len(surface.points)], np.arange(len(surface.points))
    )
    np.testing.assert_array_equal(out.points, surface.points[parent])
    assert len(out.points) == len(surface.points) + report["num_seam_points"]
    charts = out.cell_data["parameterize:chart"][0]
    assert len(np.unique(charts)) == report["num_charts"]


def test_data_regions_and_other_blocks_ride_along():
    sphere = mp.icosa_sphere(subdivisions=2)
    n = len(sphere.points)
    mesh = mp.Mesh(
        sphere.points,
        [("triangle", sphere.cells[0].data), ("line", np.array([[0, 1]]))],
        point_data={"h": sphere.points[:, 2].copy()},
        cell_data={"id": [np.arange(len(sphere.cells[0].data)), np.array([7])]},
    )
    out = mp.parameterize(mesh)
    parent = out.point_data["parameterize:parent_point"]
    np.testing.assert_array_equal(out.point_data["h"], mesh.point_data["h"][parent])
    np.testing.assert_array_equal(out.cells[1].data, [[0, 1]])  # first copies keep ids
    assert out.cell_data["parameterize:chart"][1].tolist() == [-1]
    assert out.cell_data["id"][1].tolist() == [7]
    assert len(out.points) >= n


def test_refusals():
    with pytest.raises(ValueError, match="volume"):
        mp.parameterize(mp.box((1, 1, 1)))
    with pytest.raises(ValueError, match="simplexify"):
        mp.parameterize(mp.disk(1.0, 4))  # quads
    tri = mp.disk(1.0, 4, cell_type="triangle")
    with pytest.raises(ValueError, match="padding"):
        mp.parameterize(tri, padding=-1.0)
    with pytest.raises(ValueError, match="max_stretch"):
        mp.parameterize(tri, max_stretch=0.5)


# --------------------------------------------------------------------------- #
# Texture coordinates in the formats                                           #
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def atlas():
    return mp.parameterize(mp.icosa_sphere(subdivisions=2))


@pytest.mark.parametrize("writer", ["core", "python"])
@pytest.mark.parametrize("reader", ["core", "python"])
def test_obj_writes_and_reads_vt(tmp_path, atlas, writer, reader):
    from meshioplusplus.obj import _obj

    path = str(tmp_path / "a.obj")
    _core.obj_write(path, atlas) if writer == "core" else _obj.write(path, atlas)
    text = open(path).read()
    assert "\nvt " in text and "/" in text.split("\nf ", 1)[1].split("\n", 1)[0]
    back = _core.obj_read(path) if reader == "core" else _obj.read(path)
    assert len(back.points) == len(atlas.points)  # one vt per point: no split
    np.testing.assert_allclose(
        back.point_data["uv"], atlas.point_data["uv"], atol=1e-15
    )
    np.testing.assert_array_equal(back.cells[0].data, atlas.cells[0].data)


@pytest.mark.parametrize("reader", ["core", "python"])
def test_obj_reader_splits_a_seam(tmp_path, reader):
    from meshioplusplus.obj import _obj

    path = tmp_path / "seam.obj"
    path.write_text(
        "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 1 1 0\n"
        "vt 0 0\nvt 1 0\nvt 0 1\nvt 1 1\nvt 0.5 0.5\n"
        "f 1/1 2/2 3/3\nf 2/5 4/4 3/3\n"
    )
    back = _core.obj_read(str(path)) if reader == "core" else _obj.read(str(path))
    assert back.cells[0].data.tolist() == [[0, 1, 2], [4, 3, 2]]
    assert back.point_data["obj:parent_point"].tolist() == [0, 1, 2, 3, 1]
    np.testing.assert_array_equal(back.points[4], back.points[1])
    assert back.point_data["uv"][4].tolist() == [0.5, 0.5]
    assert "obj:vt" not in back.point_data  # 5 rows for 4 points: uv says it
    # A face without a texture index leaves the file's old reading alone.
    path.write_text("v 0 0 0\nv 1 0 0\nv 0 1 0\nvt 0 0\nvt 1 0\nvt 0 1\nf 1/1 2 3\n")
    back = _core.obj_read(str(path)) if reader == "core" else _obj.read(str(path))
    assert "uv" not in back.point_data and len(back.points) == 3
    assert back.point_data["obj:vt"].shape == (3, 2)


def test_gltf_writes_texcoord_0_in_both_engines(tmp_path, atlas):
    from meshioplusplus.gltf import _gltf, _gltf_read

    a, b = tmp_path / "core.glb", tmp_path / "py.glb"
    mp.write(str(a), atlas)
    _gltf.write(str(b), atlas)
    assert a.read_bytes() == b.read_bytes()
    assert b"TEXCOORD_0" in a.read_bytes() and b"_UV" not in a.read_bytes()
    for back in (_core.gltf_read(str(a), False, None, False), _gltf_read.read(str(a))):
        # Every vertex comes back with the uv of the point it was written from.
        ref = np.c_[atlas.points, atlas.point_data["uv"]]
        got = np.c_[back.points, back.point_data["uv"]]
        dist = np.abs(got[:, None, :] - ref[None, :, :]).max(axis=2).min(axis=1)
        assert dist.max() < 1e-6  # float32


@pytest.mark.skipif(importlib.util.find_spec("pxr") is None, reason="needs usd-core")
def test_usd_writes_primvars_st(tmp_path, atlas):
    path = str(tmp_path / "a.usda")
    mp.write(path, atlas)
    assert "texCoord2f[] primvars:st" in open(path).read()
    back = mp.read(path)
    np.testing.assert_allclose(back.point_data["uv"], atlas.point_data["uv"], atol=1e-6)


def test_pipeline_step_in_both_engines(tmp_path):
    src = tmp_path / "s.vtu"
    mp.write(str(src), mp.icosa_sphere(subdivisions=2), compression=None)
    doc = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "Parameterize", "Padding": 0.05}],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    rep_py = mp.run_pipeline(doc)
    py = mp.read(str(tmp_path / "py.vtu"))
    assert "uv" in py.point_data
    if not getattr(_core, "__has_json__", False):
        pytest.skip("core built without the JSON pipeline")
    doc["Output"]["Path"] = str(tmp_path / "cpp.vtu")
    rep_cpp = _core.run_pipeline_json(json.dumps(doc))
    cpp = mp.read(str(tmp_path / "cpp.vtu"))
    np.testing.assert_array_equal(cpp.point_data["uv"], py.point_data["uv"])
    assert rep_py["steps"][0] == rep_cpp["steps"][0]
