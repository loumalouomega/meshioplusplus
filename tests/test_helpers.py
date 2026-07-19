from pathlib import Path

import pytest

import meshioplusplus

OBJ_PATH = Path(__file__).resolve().parent / "meshes" / "ply" / "bun_zipper_res4.ply"


def test_read_str():
    meshioplusplus.read(str(OBJ_PATH))


def test_read_pathlike():
    meshioplusplus.read(OBJ_PATH)


def test_read_buffer():
    # The PLY reader consumes bytes, so the buffer must be opened binary.
    with open(str(OBJ_PATH), "rb") as f:
        mesh = meshioplusplus.read(f, "ply")
    assert len(mesh.points) > 0


@pytest.fixture
def mesh():
    return meshioplusplus.read(OBJ_PATH)


def test_write_str(mesh, tmpdir):
    tmp_path = str(tmpdir.join("tmp.ply"))
    meshioplusplus.write(tmp_path, mesh)
    assert Path(tmp_path).is_file()


def test_write_pathlike(mesh, tmpdir):
    tmp_path = Path(tmpdir.join("tmp.ply"))
    meshioplusplus.write(tmp_path, mesh)
    assert Path(tmp_path).is_file()


def test_write_buffer(mesh):
    import io

    # Binary PLY is written to an in-memory bytes buffer.
    buf = io.BytesIO()
    meshioplusplus.write(buf, mesh, "ply")
    assert len(buf.getvalue()) > 0


# --- format-registry / auto-detection helper coverage ---
from meshioplusplus._helpers import (  # noqa: E402
    _filetypes_from_path,
    _pick_best_format,
)


def test_register_and_deregister_format_roundtrip(tmp_path):
    calls = {}

    def _reader(filename):
        calls["read"] = filename
        return meshioplusplus.Mesh([[0.0, 0.0, 0.0]], [])

    def _writer(filename, mesh, **kwargs):
        calls["write"] = filename

    p = tmp_path / "whatever.dummyext"
    p.write_text("")

    meshioplusplus.register_format(
        "dummyfmt", [".dummyext"], _reader, {"dummyfmt": _writer}
    )
    try:
        assert "dummyfmt" in meshioplusplus.extension_to_filetypes[".dummyext"]
        # dispatch by extension goes through the registered reader/writer
        meshioplusplus.read(p)
        assert calls["read"] == str(p)
    finally:
        meshioplusplus.deregister_format("dummyfmt")

    # After deregistration the extension no longer resolves.
    assert "dummyfmt" not in meshioplusplus.extension_to_filetypes.get(".dummyext", [])
    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.read(p)


def test_filetypes_from_path_multi_suffix():
    # `.msh` is claimed by several formats; a compound suffix still resolves the
    # trailing known extension.
    fmts = _filetypes_from_path(Path("mesh.foo.msh"))
    assert "gmsh" in fmts


def test_filetypes_from_path_unknown_raises():
    with pytest.raises(meshioplusplus.ReadError):
        _filetypes_from_path(Path("mesh.nosuchext"))


def test_pick_best_format_prefers_gmsh_on_gmsh_data():
    # Extension `.msh` maps to several formats; a mesh carrying gmsh-specific
    # data must resolve to gmsh (_helpers._pick_best_format).
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        [("triangle", [[0, 1, 2]])],
        cell_data={"gmsh:physical": [[1]]},
    )
    assert _pick_best_format(["ansys", "gmsh", "freefem"], mesh) == "gmsh"


def test_pick_best_format_falls_back_to_first():
    # No gmsh/med markers -> the first candidate wins.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        [("triangle", [[0, 1, 2]])],
    )
    assert _pick_best_format(["ansys", "gmsh", "freefem"], mesh) == "ansys"
