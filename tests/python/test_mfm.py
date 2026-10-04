import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core

from . import helpers
from .text_io_helpers import mesh_arrays, native_snapshot


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.line_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.quad_mesh,
        helpers.tet_mesh,
        helpers.hex_mesh,
        helpers.wedge_mesh,
    ],
)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.mfm.write, meshioplusplus.mfm.read, mesh, 1.0e-12
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.mfm")
    helpers.generic_io(tmp_path / "test.0.mfm")


def test_reject_mixed(tmp_path):
    with pytest.raises(meshioplusplus.WriteError):
        meshioplusplus.mfm.write(tmp_path / "x.mfm", helpers.tri_quad_mesh)


@pytest.mark.parametrize("threshold", ["0", str(1 << 40)])
@pytest.mark.parametrize("ending", ["\n", "\r\n", ""])
def test_native_whitespace_and_ownership(tmp_path, monkeypatch, threshold, ending):
    monkeypatch.setenv("MESHIOPLUSPLUS_MMAP_THRESHOLD", threshold)
    path = tmp_path / "triangle.mfm"
    body = "1 2 3\t0\v0\f0 0 0 0\r\n0 0 1 0 0 1\n7"
    text = "\n\t\r\n1 3 3 2 3 3 3 1\r\n" + body + ending
    path.write_bytes(text.encode())
    reference = meshioplusplus.mfm._py_read(path)
    native = native_snapshot(path, "mfm_read", threshold, tmp_path)
    expected = mesh_arrays(reference)
    assert sorted(native) == sorted(expected)
    for name in native:
        np.testing.assert_array_equal(native[name], expected[name])


def test_native_lenient_tokens_are_preserved(tmp_path):
    path = tmp_path / "prefix.mfm"
    path.write_text(
        "description\n1 3 3 2 3 3 3 1 99 # ignored header tail\n"
        "+1suffix 2.0 3e0 ignored reference tokens are not parsed "
        "nonnumeric 0tail +1suffix nope -0suffix 1.0suffix +7suffix trailing data"
    )
    mesh = _core.mfm_read(str(path))
    np.testing.assert_array_equal(mesh.points, [[0, 0], [1, 0], [0, 1]])
    np.testing.assert_array_equal(mesh.cells[0].data, [[0, 1, 2]])
    np.testing.assert_array_equal(mesh.cell_data["mfm:ref"][0], [7])


@pytest.mark.parametrize("count", range(16))
def test_native_truncated_sections(tmp_path, count):
    path = tmp_path / "short.mfm"
    tokens = "1 2 3 0 0 0 0 0 0 0 0 1 0 0 1 7".split()
    path.write_text("1 3 3 2 3 3 3 1\n" + " ".join(tokens[:count]))
    with pytest.raises(meshioplusplus.ReadError):
        _core.mfm_read(str(path))


def test_native_empty_mesh_without_final_newline(tmp_path):
    path = tmp_path / "empty.mfm"
    path.write_text("0 0 0 3 4 4 6 4")
    mesh = _core.mfm_read(str(path))
    assert mesh.points.shape == (0, 3)
    assert mesh.cells[0].data.shape == (0, 4)
    assert mesh.cell_data["mfm:ref"][0].shape == (0,)
