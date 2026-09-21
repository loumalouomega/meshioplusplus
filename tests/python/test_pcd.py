import pathlib
import struct

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core, _provenance
from meshioplusplus._exceptions import ReadError, WriteError
from meshioplusplus._fallback import set_strict_core
from meshioplusplus.pcd import _lzf
from meshioplusplus.pcd import _pcd as py_pcd

from . import helpers

DIR = pathlib.Path(__file__).resolve().parent / "meshes" / "pcd"
# PCL's own test corpus (BSD-3, see meshes/pcd/LICENSE.PCL): one file per DATA mode.
FIXTURES = ["bun0", "colored_cloud", "pcl_logo", "milk_color"]
MODES = ["ascii", "binary", "binary_compressed"]


def fixture(name):
    path = DIR / f"{name}.pcd"
    with open(path, "rb") as f:
        assert not f.read(24).startswith(
            b"version https://git-lfs"
        ), f"{path.name} is an unfetched Git-LFS pointer; run `git lfs pull`"
    return path


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read/write signatures."""
    if request.param == "core":
        return meshioplusplus.pcd
    return py_pcd


def assert_same(a, b):
    assert a.points.dtype == b.points.dtype
    assert np.array_equal(a.points, b.points, equal_nan=True)
    assert sorted(a.point_data) == sorted(b.point_data)
    for key, value in a.point_data.items():
        other = b.point_data[key]
        assert value.dtype == other.dtype and value.shape == other.shape, key
        assert np.array_equal(value, other, equal_nan=True), key
    assert sorted(a.field_data) == sorted(b.field_data)
    for key, value in a.field_data.items():
        assert np.array_equal(value, b.field_data[key]), key
    assert [c.type for c in a.cells] == [c.type for c in b.cells]


def write_pcd(
    path, fields, columns, data="ascii", width=None, height=1, viewpoint=None
):
    """Hand-write a PCD file: fields = [(name, type, size, count)], columns = rows."""
    points = len(columns)
    header = [
        "# .PCD v0.7 - Point Cloud Data file format",
        "VERSION 0.7",
        "FIELDS " + " ".join(f[0] for f in fields),
        "SIZE " + " ".join(str(f[2]) for f in fields),
        "TYPE " + " ".join(f[1] for f in fields),
        "COUNT " + " ".join(str(f[3]) for f in fields),
        f"WIDTH {width or points}",
        f"HEIGHT {height}",
        "VIEWPOINT " + " ".join(str(v) for v in (viewpoint or (0, 0, 0, 1, 0, 0, 0))),
        f"POINTS {points}",
        f"DATA {data}",
    ]
    body = b""
    if data == "ascii":
        body = (
            "\n".join(" ".join(str(v) for v in row) for row in columns) + "\n"
        ).encode()
    path.write_bytes(("\n".join(header) + "\n").encode() + body)
    return path


class TestReferenceFiles:
    """Files written by PCL itself, in each of its DATA modes."""

    def test_ascii_normals(self):
        mesh = meshioplusplus.read(fixture("bun0"))
        assert mesh.points.dtype == np.float32 and mesh.points.shape == (397, 3)
        assert [c.type for c in mesh.cells] == ["vertex"]
        assert mesh.point_data["normals"].shape == (397, 3)
        assert mesh.point_data["curvature"].shape == (397,)
        assert set(mesh.point_data) == {"normals", "curvature"}
        assert mesh.field_data == {}  # identity viewpoint, unorganised
        assert np.allclose(
            np.linalg.norm(mesh.point_data["normals"], axis=1), 1, atol=1e-4
        )

    def test_binary_organised_uint_rgb(self):
        mesh = meshioplusplus.read(fixture("colored_cloud"))
        assert mesh.points.shape == (1000, 3)
        assert mesh.field_data["pcd:width"].tolist() == [1]
        assert mesh.field_data["pcd:height"].tolist() == [1000]
        rgb = mesh.point_data["rgb"]
        assert rgb.dtype == np.uint8 and rgb.shape == (1000, 3)
        # TYPE U: read the raw uint32 straight out of the file, without the library
        raw = fixture("colored_cloud").read_bytes()
        body = raw[raw.index(b"DATA binary\n") + len(b"DATA binary\n") :]
        rows = np.frombuffer(
            body,
            dtype=[("xyz", "<f4", 3), ("rgb", "<u4"), ("nc", "<f4", 4)],
            count=1000,
        )
        packed = rows["rgb"]
        assert np.array_equal(rgb[:, 0], (packed >> 16) & 255)
        assert np.array_equal(rgb[:, 1], (packed >> 8) & 255)
        assert np.array_equal(rgb[:, 2], packed & 255)
        assert np.array_equal(mesh.points, rows["xyz"])

    def test_binary_compressed_float_rgb_and_viewpoint(self):
        path = fixture("pcl_logo")
        mesh = meshioplusplus.read(path)
        assert mesh.points.shape == (12909, 3)
        assert mesh.field_data["pcd:viewpoint"].tolist() == [0, 0, -2, 0, 1, 0, 0]
        # The stream must decode to exactly the declared size, into a smooth cloud.
        assert np.isfinite(mesh.points).all()
        assert np.abs(np.diff(mesh.points[:, 0])).max() < 0.5
        # TYPE F rgb: the uint32 bits sit in a float32 slot -- decode them by hand.
        raw = path.read_bytes()
        marker = b"DATA binary_compressed\n"
        at = raw.index(marker) + len(marker)
        comp, uncomp = struct.unpack("<II", raw[at : at + 8])
        assert uncomp == 12909 * 16
        soa = _lzf.decompress(raw[at + 8 : at + 8 + comp], uncomp)
        slot = np.frombuffer(soa, dtype="<f4", count=12909, offset=3 * 12909 * 4)
        packed = slot.view(np.uint32)
        expected = np.stack(
            [(packed >> 16) & 255, (packed >> 8) & 255, packed & 255], 1
        )
        assert np.array_equal(mesh.point_data["rgb"], expected.astype(np.uint8))
        assert mesh.point_data[
            "rgb"
        ].any()  # a value-cast (not bit-cast) would be all 0

    def test_binary_compressed_rgba(self):
        mesh = meshioplusplus.read(fixture("milk_color"))
        assert mesh.points.shape == (13704, 3)
        rgba = mesh.point_data["rgba"]
        assert rgba.dtype == np.uint8 and rgba.shape == (13704, 4)

    def test_points_are_writeable(self):
        mesh = meshioplusplus.read(fixture("bun0"))
        assert mesh.points.flags["WRITEABLE"]
        assert mesh.point_data["normals"].flags["WRITEABLE"]

    @pytest.mark.parametrize("name", FIXTURES)
    def test_engines_agree(self, name):
        assert_same(meshioplusplus.pcd.read(fixture(name)), py_pcd.read(fixture(name)))


@pytest.mark.parametrize("name", FIXTURES)
@pytest.mark.parametrize("mode", MODES)
def test_round_trip_reference(name, mode, engine, tmp_path):
    mesh = engine.read(fixture(name))
    out = tmp_path / "out.pcd"
    engine.write(out, mesh, data=mode, point_dtype="float32")
    assert_same(mesh, engine.read(out))
    assert f"DATA {mode}\n".encode() in out.read_bytes()


@pytest.mark.parametrize("mode", MODES)
def test_engines_write_identical_bytes(mode, tmp_path):
    mesh = py_pcd.read(fixture("colored_cloud"))
    a, b = tmp_path / "core.pcd", tmp_path / "py.pcd"
    meshioplusplus.pcd.write(a, mesh, data=mode)
    py_pcd.write(b, mesh, data=mode)
    assert a.read_bytes() == b.read_bytes()
    assert_same(py_pcd.read(a), meshioplusplus.pcd.read(b))


@pytest.mark.parametrize("mode", MODES)
def test_point_cloud_round_trip(mode, engine, tmp_path):
    mesh = helpers.point_cloud_mesh
    out = tmp_path / "cloud.pcd"
    engine.write(out, mesh, data=mode, point_dtype="float64")
    back = engine.read(out)
    assert back.points.dtype == np.float64
    assert np.array_equal(back.points, mesh.points)
    assert back.point_data["rgb"].dtype == np.uint8
    assert np.array_equal(back.point_data["rgb"], mesh.point_data["rgb"])
    assert back.point_data["label"].dtype == np.int32
    assert np.array_equal(back.point_data["label"], mesh.point_data["label"])
    assert np.array_equal(back.point_data["normals"], mesh.point_data["normals"])
    assert np.array_equal(back.point_data["intensity"], mesh.point_data["intensity"])
    assert [c.type for c in back.cells] == ["vertex"]


def test_generic_io(tmp_path):
    mesh = helpers.point_cloud_mesh
    for name in ("test.pcd", "test.0.pcd"):
        meshioplusplus.write(tmp_path / name, mesh)
        back = meshioplusplus.read(tmp_path / name)
        assert np.allclose(back.points, mesh.points, atol=1e-6)


def test_default_precision_is_float32(tmp_path):
    out = tmp_path / "a.pcd"
    meshioplusplus.pcd.write(out, helpers.point_cloud_mesh)
    assert b"SIZE 4 4 4" in out.read_bytes()
    assert meshioplusplus.pcd.read(out).points.dtype == np.float32
    out64 = tmp_path / "b.pcd"
    meshioplusplus.pcd.write(out64, helpers.point_cloud_mesh, point_dtype="float64")
    assert b"SIZE 8 8 8" in out64.read_bytes()


def test_binary_flag_selects_ascii(tmp_path):
    out = tmp_path / "a.pcd"
    meshioplusplus.pcd.write(out, helpers.point_cloud_mesh, binary=False)
    assert b"DATA ascii" in out.read_bytes()


class TestEncoding:
    def test_rgb_is_written_as_a_float_slot_and_rgba_as_uint(self, tmp_path):
        mesh = meshioplusplus.Mesh(
            np.zeros((2, 3)),
            [("vertex", np.arange(2).reshape(-1, 1))],
            point_data={
                "rgb": np.array([[255, 0, 0], [1, 2, 3]], dtype=np.uint8),
                "rgba": np.array([[9, 8, 7, 200], [0, 0, 0, 255]], dtype=np.uint8),
            },
        )
        out = tmp_path / "c.pcd"
        meshioplusplus.pcd.write(out, mesh, data="binary")
        text = out.read_bytes()
        assert b"FIELDS x y z rgb rgba\n" in text and b"TYPE F F F F U\n" in text
        body = text[text.index(b"DATA binary\n") + 12 :]
        rows = np.frombuffer(body, dtype="<u4").reshape(2, 5)  # xyz + rgb + rgba
        assert rows[0, 3] == 0x00FF0000 and rows[1, 3] == 0x00010203
        assert rows[0, 4] == 0xC8090807 and rows[1, 4] == 0xFF000000
        back = meshioplusplus.pcd.read(out)
        assert np.array_equal(back.point_data["rgb"], mesh.point_data["rgb"])
        assert np.array_equal(back.point_data["rgba"], mesh.point_data["rgba"])

    def test_multi_count_and_padding_and_int_types(self, engine, tmp_path):
        path = write_pcd(
            tmp_path / "h.pcd",
            [
                ("x", "F", 4, 1),
                ("y", "F", 4, 1),
                ("z", "F", 4, 1),
                ("hist", "F", 4, 3),
                ("_", "U", 1, 1),
                ("ring", "U", 2, 1),
                ("delta", "I", 4, 1),
            ],
            [[1, 2, 3, 10, 11, 12, 0, 5, 100], [4, 5, 6, 13, 14, 15, 0, 6, -200]],
        )
        mesh = engine.read(path)
        assert mesh.point_data["hist"].shape == (2, 3)
        assert mesh.point_data["hist"].tolist() == [[10, 11, 12], [13, 14, 15]]
        assert mesh.point_data["ring"].dtype == np.uint16
        assert mesh.point_data["delta"].dtype == np.int32
        assert mesh.point_data["delta"].tolist() == [100, -200]
        assert "_" not in mesh.point_data

    def test_ascii_float_rgb_is_a_bit_cast(self, engine, tmp_path):
        packed = np.array([0x00FF8040], dtype=np.uint32).view(np.float32)[0]
        path = write_pcd(
            tmp_path / "rgb.pcd",
            [("x", "F", 4, 1), ("y", "F", 4, 1), ("z", "F", 4, 1), ("rgb", "F", 4, 1)],
            [[0, 0, 0, format(float(packed), ".9g")]],
        )
        assert engine.read(path).point_data["rgb"].tolist() == [[255, 128, 64]]

    def test_ascii_rgba_uint(self, engine, tmp_path):
        path = write_pcd(
            tmp_path / "rgba.pcd",
            [("x", "F", 4, 1), ("y", "F", 4, 1), ("z", "F", 4, 1), ("rgba", "U", 4, 1)],
            [[0, 0, 0, 0x80FF8040]],
        )
        assert engine.read(path).point_data["rgba"].tolist() == [[255, 128, 64, 128]]

    def test_float64_points_keep_precision(self, engine, tmp_path):
        path = write_pcd(
            tmp_path / "d.pcd",
            [("x", "F", 8, 1), ("y", "F", 8, 1), ("z", "F", 8, 1)],
            [[0.1, 0.2, 0.3]],
        )
        mesh = engine.read(path)
        assert mesh.points.dtype == np.float64
        assert mesh.points.tolist() == [[0.1, 0.2, 0.3]]

    def test_precision_narrowing_is_recorded(self, tmp_path):
        lossy = meshioplusplus.Mesh(
            np.array([[0.1, 0.2, 0.3]]), [("vertex", np.array([[0]]))]
        )
        exact = meshioplusplus.Mesh(
            np.array([[0.5, 0.25, 1.0]]), [("vertex", np.array([[0]]))]
        )
        for write in (meshioplusplus.pcd.write, py_pcd.write):
            with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
                write(tmp_path / "a.pcd", lossy)
                assert [n.category for n in s.get().notes] == ["dtype"]
            with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
                write(tmp_path / "b.pcd", exact)  # representable in float32: no note
                write(tmp_path / "c.pcd", lossy, point_dtype="float64")
                assert s.get().notes == []

    def test_engines_record_the_same_notes(self, tmp_path):
        mesh = meshioplusplus.Mesh(
            np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]]),
            [("triangle", np.array([[0, 1, 2]]))],
            cell_data={"a": [np.array([1.0])]},
            field_data={"extra": np.array([1])},
        )
        seen = []
        for write in (_core.pcd_write, py_pcd.write):
            with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
                write(str(tmp_path / "n.pcd"), mesh)
                seen.append(s.get().notes)
        assert seen[0] == seen[1]
        assert {n.category for n in seen[0]} == {
            "point-padding",
            "dtype",
            "cells-dropped",
            "data-dropped",
        }


class TestOrganisedClouds:
    @pytest.fixture
    def organised(self, tmp_path):
        nan = float("nan")
        rows = [
            [0, 0, 1],
            [1, 0, 1],
            [nan, nan, nan],
            [1, 1, 1],
            [0, 1, 1],
            [nan, nan, nan],
        ]
        return write_pcd(
            tmp_path / "org.pcd",
            [("x", "F", 4, 1), ("y", "F", 4, 1), ("z", "F", 4, 1)],
            rows,
            width=3,
            height=2,
        )

    def test_kept_by_default(self, engine, organised):
        mesh = engine.read(organised)
        assert len(mesh.points) == 6
        assert np.isnan(mesh.points[2]).all()
        assert mesh.field_data["pcd:width"].tolist() == [3]
        assert mesh.field_data["pcd:height"].tolist() == [2]

    def test_drop_invalid(self, engine, organised):
        mesh = engine.read(organised, drop_invalid=True)
        assert len(mesh.points) == 4 and np.isfinite(mesh.points).all()
        assert (
            "pcd:width" not in mesh.field_data and "pcd:height" not in mesh.field_data
        )
        assert mesh.cells[0].data.shape == (4, 1)

    def test_drop_invalid_on_a_clean_cloud_keeps_the_organisation(
        self, engine, tmp_path
    ):
        path = write_pcd(
            tmp_path / "c.pcd",
            [("x", "F", 4, 1), ("y", "F", 4, 1), ("z", "F", 4, 1)],
            [[0, 0, 0], [1, 0, 0]],
            width=1,
            height=2,
        )
        assert engine.read(path, drop_invalid=True).field_data[
            "pcd:height"
        ].tolist() == [2]

    @pytest.mark.parametrize("mode", MODES)
    def test_round_trip_restores_the_shape(self, engine, organised, mode, tmp_path):
        mesh = engine.read(organised)
        out = tmp_path / "again.pcd"
        engine.write(out, mesh, data=mode)
        assert b"WIDTH 3\nHEIGHT 2\n" in out.read_bytes()
        assert_same(mesh, engine.read(out))


def test_viewpoint_round_trip(engine, tmp_path):
    path = write_pcd(
        tmp_path / "v.pcd",
        [("x", "F", 4, 1), ("y", "F", 4, 1), ("z", "F", 4, 1)],
        [[1, 2, 3]],
        viewpoint=(1, 2, 3, 0, 1, 0, 0),
    )
    mesh = engine.read(path)
    assert mesh.field_data["pcd:viewpoint"].tolist() == [1, 2, 3, 0, 1, 0, 0]
    assert mesh.points.tolist() == [[1, 2, 3]]  # recorded, not applied
    out = tmp_path / "w.pcd"
    engine.write(out, mesh)
    assert b"VIEWPOINT 1 2 3 0 1 0 0\n" in out.read_bytes()


def test_identity_viewpoint_is_not_stored(engine, tmp_path):
    out = tmp_path / "i.pcd"
    engine.write(out, helpers.point_cloud_mesh)
    assert b"VIEWPOINT 0 0 0 1 0 0 0\n" in out.read_bytes()
    assert "pcd:viewpoint" not in engine.read(out).field_data


class TestWriterRules:
    def test_cells_and_cell_data_are_dropped(self, engine, tmp_path):
        mesh = meshioplusplus.Mesh(
            np.random.default_rng(0).random((3, 3)),
            [("triangle", np.array([[0, 1, 2]]))],
            cell_data={"a": [np.array([1.0])]},
        )
        out = tmp_path / "t.pcd"
        engine.write(out, mesh)
        back = engine.read(out)
        assert len(back.points) == 3 and [c.type for c in back.cells] == ["vertex"]
        assert back.cell_data == {}

    def test_2d_points_are_padded(self, engine, tmp_path):
        mesh = meshioplusplus.Mesh(np.array([[1.0, 2.0], [3.0, 4.0]]), [])
        out = tmp_path / "p.pcd"
        engine.write(out, mesh)
        assert engine.read(out).points.tolist() == [[1, 2, 0], [3, 4, 0]]

    def test_empty_cloud(self, engine, tmp_path):
        for mode in MODES:
            out = tmp_path / f"e_{mode}.pcd"
            engine.write(out, helpers.empty_mesh, data=mode)
            back = engine.read(out)
            assert back.points.shape == (0, 3) and back.cells[0].data.shape == (0, 1)

    def test_awkward_point_data_names_are_made_safe(self, engine, tmp_path):
        mesh = meshioplusplus.Mesh(
            np.zeros((2, 3)),
            [("vertex", np.arange(2).reshape(-1, 1))],
            point_data={"my field": np.array([1.0, 2.0]), "x": np.array([3.0, 4.0])},
        )
        out = tmp_path / "n.pcd"
        engine.write(out, mesh)
        back = engine.read(out)
        assert back.point_data["my_field"].tolist() == [1, 2]
        assert back.point_data["x_2"].tolist() == [3, 4]

    def test_higher_rank_point_data_is_flattened(self, engine, tmp_path):
        mesh = meshioplusplus.Mesh(
            np.zeros((2, 3)),
            [("vertex", np.arange(2).reshape(-1, 1))],
            point_data={"m": np.arange(8.0).reshape(2, 2, 2)},
        )
        out = tmp_path / "r.pcd"
        engine.write(out, mesh)
        assert engine.read(out).point_data["m"].shape == (2, 4)

    def test_rejects_bad_options(self, tmp_path):
        with pytest.raises(WriteError):
            meshioplusplus.pcd.write(
                tmp_path / "a.pcd", helpers.point_cloud_mesh, data="zip"
            )
        with pytest.raises(WriteError):
            meshioplusplus.pcd.write(
                tmp_path / "a.pcd", helpers.point_cloud_mesh, point_dtype="float16"
            )

    def test_input_mesh_is_not_modified(self, engine, tmp_path):
        import copy

        mesh = copy.deepcopy(helpers.point_cloud_mesh)
        engine.write(tmp_path / "m.pcd", mesh, data="binary_compressed")
        assert np.array_equal(mesh.points, helpers.point_cloud_mesh.points)
        assert mesh.point_data.keys() == helpers.point_cloud_mesh.point_data.keys()


class TestMalformed:
    def fields(self):
        return [("x", "F", 4, 1), ("y", "F", 4, 1), ("z", "F", 4, 1)]

    @pytest.mark.parametrize(
        "mutate, message",
        [
            (lambda t: t.split("DATA")[0], "DATA"),
            (lambda t: t.replace("DATA ascii", "DATA zip"), "unknown DATA mode"),
            (lambda t: t.replace("SIZE 4 4 4", "SIZE 4 4"), "disagree"),
            (lambda t: t.replace("TYPE F F F", "TYPE F F X"), "unsupported field type"),
            (lambda t: t.replace("FIELDS x y z", "FIELDS x y w"), "no scalar 'z'"),
            (
                lambda t: t.replace("VIEWPOINT 0 0 0 1 0 0 0", "VIEWPOINT 1 2"),
                "VIEWPOINT",
            ),
            (lambda t: t.replace("1 2 3", "1 2"), "values"),
            (lambda t: t.replace("1 2 3", "1 2 abc"), "non-numeric"),
        ],
    )
    def test_bad_ascii(self, engine, tmp_path, mutate, message):
        good = write_pcd(tmp_path / "g.pcd", self.fields(), [[1, 2, 3]])
        bad = tmp_path / "bad.pcd"
        bad.write_bytes(mutate(good.read_text()).encode())
        with pytest.raises(ReadError, match=message):
            engine.read(bad)

    def test_truncated_binary(self, engine, tmp_path):
        mesh = helpers.point_cloud_mesh
        good = tmp_path / "g.pcd"
        engine.write(good, mesh, data="binary")
        bad = tmp_path / "bad.pcd"
        bad.write_bytes(good.read_bytes()[:-9])
        with pytest.raises(ReadError, match="shorter"):
            engine.read(bad)

    def test_corrupt_compressed(self, engine, tmp_path):
        good = tmp_path / "g.pcd"
        engine.write(good, helpers.point_cloud_mesh, data="binary_compressed")
        raw = good.read_bytes()
        at = raw.index(b"DATA binary_compressed\n") + len(b"DATA binary_compressed\n")
        for name, blob in {
            "no prefix": raw[: at + 3],
            "short stream": raw[: at + 8 + 5],
            "wrong size": raw[:at]
            + struct.pack("<II", *struct.unpack("<II", raw[at : at + 8])[:1], 7)
            + raw[at + 8 :],
            "garbage": raw[: at + 8] + b"\xff" * (len(raw) - at - 8),
        }.items():
            bad = tmp_path / "bad.pcd"
            bad.write_bytes(blob)
            with pytest.raises(ReadError):
                engine.read(bad)

    def test_missing_file(self):
        with pytest.raises(Exception):
            meshioplusplus.pcd.read("/nonexistent/none.pcd")


class TestLzf:
    @pytest.mark.parametrize(
        "blob",
        [
            b"",
            b"a",
            b"abc",
            b"a" * 1000,
            bytes(range(256)) * 20,
            np.random.default_rng(1).bytes(6000),
        ],
    )
    def test_round_trip(self, blob):
        assert _lzf.decompress(_lzf.compress(blob), len(blob)) == blob

    def test_compresses_repetition(self):
        assert len(_lzf.compress(b"ab" * 5000)) < 200

    def test_decodes_a_liblzf_stream(self):
        # literal "abc", then a 5-byte back reference at distance 3 ("bcabc" -> overlap)
        stream = bytes([2]) + b"abc" + bytes([(5 - 2) << 5, 2])
        assert _lzf.decompress(stream, 8) == b"abcabcab"

    def test_rejects_bad_streams(self):
        with pytest.raises(ReadError):
            _lzf.decompress(bytes([5, 1, 2]), 6)  # truncated literal run
        with pytest.raises(ReadError):
            _lzf.decompress(bytes([0x20, 5]), 4)  # reference before the start
        with pytest.raises(ReadError):
            _lzf.decompress(
                bytes([1]) + b"ab", 5
            )  # decodes to fewer bytes than declared


def test_native_path_is_used():
    set_strict_core(True)
    try:
        mesh = meshioplusplus.pcd.read(fixture("pcl_logo"))
        assert mesh.points.shape == (12909, 3)
    finally:
        set_strict_core(None)


def test_core_has_the_functions():
    for name in ("pcd_read", "pcd_write"):
        assert hasattr(_core, name)


class TestSniff:
    @pytest.mark.parametrize("name", FIXTURES)
    def test_fixtures(self, name):
        assert meshioplusplus.sniff_format(fixture(name)) == "pcd"

    def test_header_without_comment(self, tmp_path):
        path = tmp_path / "x.cloud"
        path.write_text(
            "VERSION .7\nFIELDS x y z\nSIZE 4 4 4\nTYPE F F F\nCOUNT 1 1 1\n"
            "WIDTH 1\nHEIGHT 1\nPOINTS 1\nDATA ascii\n1 2 3\n"
        )
        assert meshioplusplus.sniff_format(path) == "pcd"
        assert meshioplusplus.read(path).points.tolist() == [[1, 2, 3]]

    def test_unknown_extension_is_read_by_content(self, tmp_path):
        out = tmp_path / "cloud.bin"
        meshioplusplus.pcd.write(out, helpers.point_cloud_mesh)
        assert len(meshioplusplus.read(out).points) == 40


def test_buffers(tmp_path):
    import io

    mesh = helpers.point_cloud_mesh
    buf = io.BytesIO()
    py_pcd.write(buf, mesh, data="binary_compressed", point_dtype="float64")
    back = meshioplusplus.pcd.read(io.BytesIO(buf.getvalue()))
    assert np.array_equal(back.points, mesh.points)


class TestPipeline:
    """Roadmap 1.1: .pcd -> subsample_points -> proximity_graph -> .pcd/.vtu, from the CLI."""

    def test_cli_chain(self, tmp_path):
        sub = tmp_path / "sub.pcd"
        graph = tmp_path / "graph.vtu"
        meshioplusplus._cli.main(
            ["subsample", str(fixture("pcl_logo")), str(sub), "-n", "500"]
        )
        assert len(meshioplusplus.read(sub).points) == 500
        meshioplusplus._cli.main(
            ["proximity-graph", str(sub), str(graph), "--knn", "4", "-q"]
        )
        edges = meshioplusplus.read(graph)
        assert [c.type for c in edges.cells] == ["line"]
        assert "degree" in edges.point_data

    def test_subsample_keeps_colours_without_a_spurious_warning(
        self, tmp_path, recwarn
    ):
        mesh = meshioplusplus.read(fixture("pcl_logo"))
        out = meshioplusplus.subsample_points(mesh, 300)
        result = out[0] if isinstance(out, tuple) else out
        assert result.point_data["rgb"].shape == (300, 3)
        assert not [w for w in recwarn if "dropped" in str(w.message)]

    def test_convert_pcd_to_vtu_and_back(self, tmp_path):
        vtu = tmp_path / "c.vtu"
        back = tmp_path / "c.pcd"
        meshioplusplus._cli.main(["convert", str(fixture("bun0")), str(vtu)])
        meshioplusplus._cli.main(["convert", str(vtu), str(back)])
        assert np.allclose(
            meshioplusplus.read(back).points,
            meshioplusplus.read(fixture("bun0")).points,
        )


class TestCliVerbs:
    """ascii / binary / compress / decompress rewrite a PCD file in place."""

    def modes(self, tmp_path, mesh):
        path = tmp_path / "in.pcd"
        py_pcd.write(path, mesh, data="binary", point_dtype="float64")
        seen = []
        for verb in ("ascii", "compress", "decompress", "binary"):
            meshioplusplus._cli.main([verb, str(path)])
            head = path.read_bytes().split(b"POINTS")[1].split(b"\n")[1]
            seen.append(head.decode())
            back = meshioplusplus.read(path)
            assert back.points.dtype == np.float64  # rewritten in place, not narrowed
            assert np.array_equal(back.points, mesh.points)
            assert np.array_equal(back.point_data["rgb"], mesh.point_data["rgb"])
        return seen

    def test_verbs(self, tmp_path):
        assert self.modes(tmp_path, helpers.point_cloud_mesh) == [
            "DATA ascii",
            "DATA binary_compressed",
            "DATA binary",
            "DATA binary",
        ]

    def test_keep_follows_the_meshes_precision(self, tmp_path):
        mesh32 = meshioplusplus.Mesh(
            helpers.point_cloud_mesh.points.astype(np.float32),
            helpers.point_cloud_mesh.cells,
        )
        out = tmp_path / "k.pcd"
        meshioplusplus.pcd.write(out, mesh32, point_dtype="keep")
        assert b"SIZE 4 4 4" in out.read_bytes()
        meshioplusplus.pcd.write(out, helpers.point_cloud_mesh, point_dtype="keep")
        assert b"SIZE 8 8 8" in out.read_bytes()
