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


@pytest.mark.parametrize("file_format", ["tetgen", "triangle", "ensight"])
def test_read_multifile_format_from_buffer_raises(file_format):
    # tetgen/triangle/ensight are spread across several sibling files and so
    # cannot be served from a single in-memory buffer (_helpers.py:77-81).
    import io

    with pytest.raises(meshioplusplus.ReadError, match="multiple files"):
        meshioplusplus.read(io.BytesIO(b"whatever"), file_format=file_format)


@pytest.mark.parametrize("file_format", ["tetgen", "triangle", "ensight"])
def test_write_multifile_format_to_buffer_raises(file_format):
    # Symmetric guard on the write path (_helpers.py:170-174).
    import io

    with pytest.raises(meshioplusplus.WriteError, match="multiple files"):
        meshioplusplus.write(io.BytesIO(), helpers.tri_mesh, file_format=file_format)


def test_write_buffer_without_format_raises():
    import io

    with pytest.raises(meshioplusplus.WriteError):
        meshioplusplus.write(io.BytesIO(), helpers.tri_mesh)


def test_write_bad_cell_shape_raises(tmp_path):
    # A "triangle" block whose rows do not have 3 nodes must be rejected by the
    # cell-shape sanity check in write() (_helpers.py:190-198).
    bad = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        [("triangle", np.array([[0, 1, 2, 3]]))],
    )
    with pytest.raises(meshioplusplus.WriteError, match="cells array shape"):
        meshioplusplus.write(tmp_path / "mesh.vtu", bad)


def test_read_multiformat_fallback_success(tmp_path):
    # `.msh` is claimed by several formats (ansys, gmsh, freefem). A gmsh file
    # is not readable as ansys, so the dispatcher must fall through to gmsh.
    p = tmp_path / "mesh.msh"
    meshioplusplus.write(p, helpers.tri_mesh, file_format="gmsh")
    mesh = meshioplusplus.read(p)  # no explicit file_format -> tries each
    assert np.allclose(mesh.points, helpers.tri_mesh.points)


def test_read_multiformat_fallback_all_fail(tmp_path):
    # None of the `.msh` formats can read pure garbage -> the "as either of ..."
    # error lists every candidate (_helpers.py:107-113).
    p = tmp_path / "mesh.msh"
    p.write_text("this is not a mesh at all\n")
    with pytest.raises(meshioplusplus.ReadError, match="either of"):
        meshioplusplus.read(p)


def test_corruption_error_is_distinct_exception():
    # CorruptionError is raised by the VTU reader on truncated/corrupt appended
    # data (vtu/_vtu.py:567). Pin down that it exists and is its own type.
    from meshioplusplus._exceptions import CorruptionError

    assert issubclass(CorruptionError, Exception)
    assert CorruptionError is not meshioplusplus.ReadError
