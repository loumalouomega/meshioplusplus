import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._exceptions import ReadError
from meshioplusplus._fallback import set_strict_core
from meshioplusplus.xyz import _xyz as py_xyz

from . import helpers

ALIASES = [".xyz", ".xyzn", ".xyzrgb", ".asc", ".pts", ".txt"]


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read/write signatures."""
    if request.param == "core":
        return meshioplusplus.xyz
    return py_xyz


def assert_same(a, b):
    assert a.points.dtype == b.points.dtype
    assert np.array_equal(a.points, b.points, equal_nan=True)
    assert sorted(a.point_data) == sorted(b.point_data)
    for key, value in a.point_data.items():
        other = b.point_data[key]
        assert value.dtype == other.dtype and value.shape == other.shape, key
        assert np.array_equal(value, other, equal_nan=True), key
    assert [c.type for c in a.cells] == [c.type for c in b.cells]


def write_text(path, rows, comments=()):
    path.write_text("\n".join([*comments, *rows]) + "\n")
    return path


class TestColumnSniffing:
    def test_three_columns_are_xyz(self, engine, tmp_path):
        mesh = engine.read(write_text(tmp_path / "a.xyz", ["1 2 3", "4 5 6"]))
        assert mesh.points.tolist() == [[1, 2, 3], [4, 5, 6]]
        assert mesh.point_data == {}
        assert [c.type for c in mesh.cells] == ["vertex"]

    def test_four_columns_are_a_scalar(self, engine, tmp_path):
        mesh = engine.read(write_text(tmp_path / "a.xyz", ["1 2 3 0.5"]))
        assert mesh.point_data["scalar"].tolist() == [0.5]

    def test_six_columns_by_value_range(self, engine, tmp_path):
        normals = engine.read(write_text(tmp_path / "n.xyz", ["0 0 0 0 0 1"]))
        assert normals.point_data["normals"].tolist() == [[0, 0, 1]]
        colours = engine.read(write_text(tmp_path / "c.xyz", ["0 0 0 255 0 0"]))
        assert colours.point_data["rgb"].tolist() == [[255, 0, 0]]

    def test_six_columns_by_extension(self, engine, tmp_path):
        # Unit-range values that look like neither normals nor bytes are
        # ambiguous as .xyz, but the alias decides.
        rows = ["0 0 0 0.5 0.5 0.5"]
        assert engine.read(write_text(tmp_path / "a.xyzn", rows)).point_data[
            "normals"
        ].tolist() == [[0.5, 0.5, 0.5]]
        # 0.5-valued columns are not byte-like, so .xyzrgb reads them as unit
        # floats scaled to bytes.
        assert engine.read(write_text(tmp_path / "a.xyzrgb", rows)).point_data[
            "rgb"
        ].tolist() == [[128, 128, 128]]

    def test_pts_seven_columns(self, engine, tmp_path):
        mesh = engine.read(
            write_text(
                tmp_path / "a.pts", ["2", "0 0 0 1.5 255 0 0", "1 1 1 2.5 0 0 255"]
            )
        )
        assert mesh.points.tolist() == [[0, 0, 0], [1, 1, 1]]
        assert mesh.point_data["intensity"].tolist() == [1.5, 2.5]
        assert mesh.point_data["rgb"].tolist() == [[255, 0, 0], [0, 0, 255]]

    def test_pts_count_is_validated(self, engine, tmp_path):
        with pytest.raises(ReadError, match="declares"):
            engine.read(write_text(tmp_path / "a.pts", ["3", "0 0 0"]))

    def test_ambiguous_columns_ask_for_names(self, engine, tmp_path):
        with pytest.raises(ReadError, match="columns="):
            engine.read(write_text(tmp_path / "a.xyz", ["0 0 0 1 2"]))

    def test_chemistry_xyz_is_refused_by_name(self, engine, tmp_path):
        with pytest.raises(ReadError, match="[Cc]hemistry"):
            engine.read(
                write_text(tmp_path / "mol.xyz", ["2", "comment", "C 0 0 0", "H 1 0 0"])
            )


class TestSyntax:
    @pytest.mark.parametrize("sep", [" ", ",", ";", "  ,  "])
    def test_delimiters(self, engine, tmp_path, sep):
        mesh = engine.read(write_text(tmp_path / "a.xyz", [sep.join(["1", "2", "3"])]))
        assert mesh.points.tolist() == [[1, 2, 3]]

    def test_comments_and_blank_lines(self, engine, tmp_path):
        mesh = engine.read(
            write_text(
                tmp_path / "a.xyz",
                ["# a comment", "// another", "", "1 2 3", ""],
            )
        )
        assert mesh.points.tolist() == [[1, 2, 3]]

    def test_header_comment_names_the_columns(self, engine, tmp_path):
        mesh = engine.read(
            write_text(tmp_path / "a.xyz", ["# x y z temperature", "1 2 3 300"])
        )
        assert mesh.point_data["temperature"].tolist() == [300]

    def test_explicit_columns_and_skip(self, engine, tmp_path):
        mesh = engine.read(
            write_text(tmp_path / "a.xyz", ["1 2 3 9 0.5"]),
            columns=["x", "y", "z", "_", "p"],
        )
        assert mesh.point_data.keys() == {"p"}
        assert mesh.point_data["p"].tolist() == [0.5]

    def test_explicit_delimiter(self, engine, tmp_path):
        mesh = engine.read(write_text(tmp_path / "a.xyz", ["1,2,3"]), delimiter=",")
        assert mesh.points.tolist() == [[1, 2, 3]]

    def test_empty_file(self, engine, tmp_path):
        mesh = engine.read(write_text(tmp_path / "a.xyz", []))
        assert mesh.points.shape == (0, 3)
        assert mesh.cells[0].data.shape == (0, 1)

    def test_ragged_and_nonnumeric_rows(self, engine, tmp_path):
        with pytest.raises(ReadError, match="columns"):
            engine.read(write_text(tmp_path / "a.xyz", ["1 2 3", "1 2"]))
        with pytest.raises(ReadError, match="[Nn]umeric"):
            engine.read(write_text(tmp_path / "a.xyz", ["1 2 oops"]))


class TestRoundTrip:
    def test_point_cloud(self, engine, tmp_path):
        out = tmp_path / "cloud.xyz"
        engine.write(out, helpers.point_cloud_mesh)
        back = engine.read(out)
        assert np.allclose(back.points, helpers.point_cloud_mesh.points)
        assert np.array_equal(
            back.point_data["rgb"], helpers.point_cloud_mesh.point_data["rgb"]
        )
        assert np.allclose(
            back.point_data["normals"], helpers.point_cloud_mesh.point_data["normals"]
        )
        assert np.allclose(
            back.point_data["intensity"],
            helpers.point_cloud_mesh.point_data["intensity"],
        )
        assert [c.type for c in back.cells] == ["vertex"]

    def test_scalar_name_survives(self, engine, tmp_path):
        mesh = meshioplusplus.Mesh(
            np.array([[1.0, 2.0, 3.0]]),
            [("vertex", np.array([[0]]))],
            point_data={"temperature": np.array([300.0])},
        )
        out = tmp_path / "t.xyz"
        engine.write(out, mesh)
        assert engine.read(out).point_data["temperature"].tolist() == [300.0]

    def test_engines_agree(self, tmp_path):
        out = tmp_path / "cloud.xyz"
        py_xyz.write(out, helpers.point_cloud_mesh)
        assert_same(meshioplusplus.xyz.read(out), py_xyz.read(out))

    def test_engines_write_identical_bytes(self, tmp_path):
        a, b = tmp_path / "core.xyz", tmp_path / "py.xyz"
        meshioplusplus.xyz.write(a, helpers.point_cloud_mesh)
        py_xyz.write(b, helpers.point_cloud_mesh)
        assert a.read_bytes() == b.read_bytes()

    def test_generic_io_over_every_alias(self, tmp_path):
        for suffix in ALIASES:
            path = tmp_path / f"cloud{suffix}"
            meshioplusplus.write(path, helpers.point_cloud_mesh)
            back = meshioplusplus.read(path)
            assert np.allclose(back.points, helpers.point_cloud_mesh.points)

    def test_cells_are_dropped_and_2d_is_padded(self, engine, tmp_path):
        mesh = meshioplusplus.Mesh(
            np.array([[1.0, 2.0], [3.0, 4.0]]),
            [("triangle", np.array([[0, 1, 0]]))],
        )
        out = tmp_path / "p.xyz"
        engine.write(out, mesh)
        assert engine.read(out).points.tolist() == [[1, 2, 0], [3, 4, 0]]

    def test_empty_cloud(self, engine, tmp_path):
        out = tmp_path / "e.xyz"
        engine.write(out, helpers.empty_mesh)
        back = engine.read(out)
        assert back.points.shape == (0, 3)


def test_native_path_is_used(tmp_path):
    path = tmp_path / "cloud.xyz"
    py_xyz.write(path, helpers.point_cloud_mesh)
    set_strict_core(True)
    try:
        mesh = meshioplusplus.xyz.read(path)
        assert np.allclose(mesh.points, helpers.point_cloud_mesh.points)
    finally:
        set_strict_core(None)


def test_core_has_the_functions():
    for name in ("xyz_read", "xyz_write"):
        assert hasattr(_core, name)


def test_buffers(tmp_path):
    import io

    buf = io.BytesIO()
    py_xyz.write(buf, helpers.point_cloud_mesh)
    back = meshioplusplus.xyz.read(io.BytesIO(buf.getvalue()))
    assert np.allclose(back.points, helpers.point_cloud_mesh.points)
