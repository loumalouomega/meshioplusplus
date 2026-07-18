import numpy as np
import pytest

import meshioplusplus

from . import helpers


def test_public_attributes():
    # Just make sure this is here
    meshioplusplus.extension_to_filetypes


def test_read_unknown_extension(tmp_path):
    # An extension no format claims cannot be deduced -> ReadError.
    p = tmp_path / "mesh.nosuchext"
    p.write_text("garbage")
    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.read(p)


def test_read_unknown_format_name(tmp_path):
    p = tmp_path / "mesh.vtu"
    meshioplusplus.write(p, helpers.tri_mesh)
    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.read(p, file_format="not-a-real-format")


def test_read_missing_file(tmp_path):
    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.read(tmp_path / "does_not_exist.vtu")


def test_write_unknown_extension(tmp_path):
    # The extension-inference helper is shared with read, so an undeducible
    # extension surfaces as ReadError even on the write path.
    with pytest.raises((meshioplusplus.ReadError, meshioplusplus.WriteError)):
        meshioplusplus.write(tmp_path / "mesh.nosuchext", helpers.tri_mesh)


def test_write_unknown_format_name(tmp_path):
    with pytest.raises(meshioplusplus.WriteError):
        meshioplusplus.write(
            tmp_path / "mesh.vtu", helpers.tri_mesh, file_format="not-a-real-format"
        )


def test_read_buffer_without_format_raises():
    import io

    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.read(io.BytesIO(b"whatever"))


def test_roundtrip_via_public_read_write(tmp_path):
    # Exercise the top-level read/write dispatch (extension inference) end to end.
    p = tmp_path / "mesh.vtu"
    meshioplusplus.write(p, helpers.tri_mesh)
    mesh = meshioplusplus.read(p)
    assert np.allclose(mesh.points, helpers.tri_mesh.points)
