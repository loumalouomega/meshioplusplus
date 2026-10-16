import pathlib

import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        # helpers.empty_mesh,
        helpers.tri_mesh
    ],
)
def test_neuroglancer(mesh, tmp_path):
    def writer(*args, **kwargs):
        return meshioplusplus.neuroglancer.write(*args, **kwargs)

    # 32bit only
    helpers.write_read(tmp_path, writer, meshioplusplus.neuroglancer.read, mesh, 1.0e-8)


@pytest.mark.parametrize("filename, ref_sum, ref_num_cells", [("simple1", 20, 4)])
def test_reference_file(filename, ref_sum, ref_num_cells):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "neuroglancer" / filename

    mesh = meshioplusplus.read(filename, "neuroglancer")
    tol = 1.0e-5
    s = np.sum(mesh.points)
    assert abs(s - ref_sum) < tol * abs(ref_sum)
    assert len(mesh.cells[0].data) == ref_num_cells


def _native(tmp_path):
    from meshioplusplus import _core
    from meshioplusplus.neuroglancer import _neuroglancer as twin

    return _core, twin


def test_engines_write_the_same_bytes(tmp_path):
    core, twin = _native(tmp_path)
    mesh = helpers.tri_mesh
    a, b = tmp_path / "core", tmp_path / "twin"
    core.neuroglancer_write(str(a), mesh)
    twin.write(str(b), mesh)
    assert a.read_bytes() == b.read_bytes()
    for path in (a, b):
        got_core, got_twin = core.neuroglancer_read(str(path)), twin.read(str(path))
        np.testing.assert_array_equal(got_core.points, got_twin.points)
        np.testing.assert_array_equal(got_core.cells[0].data, got_twin.cells[0].data)


def test_native_reads_the_reference_file():
    core, _ = _native(None)
    this_dir = pathlib.Path(__file__).resolve().parent
    mesh = core.neuroglancer_read(str(this_dir / "meshes" / "neuroglancer" / "simple1"))
    assert abs(np.sum(mesh.points) - 20) < 1e-5 * 20
    assert len(mesh.cells[0].data) == 4


@pytest.mark.parametrize(
    "payload",
    [
        b"\x01",  # shorter than the count
        b"\x05\x00\x00\x00" + bytes(12),  # count larger than the file
        b"\x01\x00\x00\x00" + bytes(12) + bytes(5),  # ragged triangle table
        b"\x01\x00\x00\x00"
        + bytes(12)
        + b"\x00\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00",
    ],
)
def test_native_rejects_malformed_files(tmp_path, payload):
    core, _ = _native(tmp_path)
    path = tmp_path / "bad"
    path.write_bytes(payload)
    with pytest.raises(meshioplusplus.ReadError):
        core.neuroglancer_read(str(path))
