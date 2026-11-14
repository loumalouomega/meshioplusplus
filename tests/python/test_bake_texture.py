"""bake_texture (roadmap 2.4.8): both engines draw the same pixels."""

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _bake_texture


def _square():
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    return mio.Mesh(
        pts,
        [("triangle", np.array([[0, 1, 2]])), ("quad", np.array([[0, 2, 3, 0]]))],
        point_data={
            "uv": pts[:, :2] * 0.9 + 0.05,
            "h": np.array([0.0, 1.0, 2.0, 3.0]),
            "vec": np.array([[3.0, 4.0], [0, 0], [1, 0], [0, 1]]),
        },
        cell_data={"c": [np.array([1.0]), np.array([np.nan])]},
    )


@pytest.mark.parametrize(
    "array, kw",
    [("h", {}), ("vec", {}), ("vec", {"component": 1}), ("c", {"vmin": 0, "vmax": 2})],
)
def test_the_engines_draw_the_same_pixels(array, kw):
    mesh = _square()
    out, image = mio.bake_texture(mesh, array, size=(37, 23), **kw)
    twin = _bake_texture._twin(
        mesh,
        array,
        37,
        23,
        "viridis",
        kw.get("vmin"),
        kw.get("vmax"),
        kw.get("component", -1),
    )
    assert image.shape == (23, 37, 4)
    assert np.array_equal(image, twin["image"])
    assert bytes(out.field_data["texture:png"])[:8] == b"\x89PNG\r\n\x1a\n"
    assert 0 < np.count_nonzero(image[..., 3]) < 37 * 23


def test_a_nan_cell_stays_transparent_and_a_bad_call_raises():
    mesh = _square()
    _, image = mio.bake_texture(mesh, "c", size=(20, 20))
    # The quad (NaN) covers the upper-left half; the triangle the lower-right.
    assert image[2, 2, 3] == 0 and image[17, 17, 3] == 255
    with pytest.raises(ValueError, match="no point or cell data"):
        mio.bake_texture(mesh, "missing")
    no_uv = mio.Mesh(mesh.points, mesh.cells, point_data={"h": mesh.point_data["h"]})
    with pytest.raises(ValueError, match="uv"):
        mio.bake_texture(no_uv, "h")


def test_the_baked_texture_rides_into_gltf(tmp_path):
    mesh, image = mio.bake_texture(_square(), "h", size=(8, 8))
    path = tmp_path / "t.glb"
    mio.write(path, mesh)
    back = mio.read(path)
    assert bytes(back.field_data["texture:png"]) == bytes(
        mesh.field_data["texture:png"]
    )
