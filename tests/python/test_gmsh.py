import copy
import pathlib
import struct
from functools import partial

import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.fixture
def periodic_native_only(monkeypatch):
    from meshioplusplus import _fallback

    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    monkeypatch.setattr(_fallback, "_strict", True)


def gmsh_periodic():
    mesh = copy.deepcopy(helpers.quad_mesh)
    trns = [0] * 16  # just for io testing
    mesh.gmsh_periodic = [
        [0, (3, 1), None, [[2, 0]]],
        [0, (4, 6), None, [[3, 5]]],
        [1, (2, 1), trns, [[5, 0], [4, 1], [4, 2]]],
    ]
    return mesh


def periodic_fixture(version="4.1", binary=False, width=8):
    """Independent sparse-node file, including ordered duplicate periodic pairs."""
    if version == "4.0":
        data = gmsh40_fixture(binary, width)
    else:
        data = f"$MeshFormat\n{version} {int(binary)} {width}\n".encode()
        if binary:
            data += struct.pack("=i", 1) + b"\n"
        data += b"$EndMeshFormat\n$Nodes\n"
        tags = [30, 10, 40, 20]
        points = [[1, 0, 0], [0, 0, 0], [1, 1, 0], [0, 1, 0]]
        size = "I" if width == 4 else "Q"
        if version == "2.2":
            data += b"4\n"
            for tag, xyz in zip(tags, points):
                data += (
                    struct.pack("=i3d", tag, *xyz)
                    if binary
                    else (" ".join(map(str, [tag, *xyz])) + "\n").encode()
                )
        elif binary:
            data += struct.pack("=" + size * 4, 1, 4, 10, 40)
            data += struct.pack("=3i" + size, 2, 22, 0, 4)
            data += struct.pack("=" + size * 4, *tags)
            data += struct.pack("=12d", *np.ravel(points))
        else:
            data += b"1 4 10 40\n2 22 0 4\n30 10 40 20\n1 0 0\n0 0 0\n1 1 0\n0 1 0\n"
        data += b"\n$EndNodes\n$Elements\n"
        if version == "2.2":
            data += b"1\n"
            data += (
                struct.pack("=10i", 3, 1, 2, 1, 7, 22, 10, 30, 40, 20)
                if binary
                else b"1 3 2 7 22 10 30 40 20\n"
            )
        elif binary:
            data += struct.pack("=" + size * 4, 1, 1, 1, 1)
            data += struct.pack("=3i" + size, 2, 22, 3, 1)
            data += struct.pack("=" + size * 5, 1, 10, 30, 40, 20)
        else:
            data += b"1 1 1 1\n2 22 3 1\n1 10 30 40 20\n"
        data += b"\n$EndElements\n"
    affine = np.eye(4).ravel()
    affine[3] = 1
    pairs = [30, 10, 40, 20, 30, 10]
    data += b"$Periodic\n"
    if not binary or version == "2.2":
        data += b"2\n1 11 33\n"
        if version == "4.1":
            data += b"0\n"
        data += b"0\n1 12 33\n"
        data += b"16\n" if version == "4.1" else b"Affine "
        data += (" ".join(map(str, affine)) + "\n3\n30 10\n40 20\n30 10\n").encode()
    else:
        size = "I" if width == 4 else "Q"
        data += struct.pack("=" + (size if version == "4.1" else "i"), 2)
        data += struct.pack("=3i", 1, 11, 33)
        if version == "4.1":
            data += struct.pack("=" + size, 0)
        data += struct.pack("=" + size, 0)
        data += struct.pack("=3i", 1, 12, 33)
        data += struct.pack(
            "=" + size, 16 if version == "4.1" else 2 ** (width * 8) - 1
        )
        data += struct.pack("=16d", *affine)
        data += struct.pack("=" + size, 3)
        data += struct.pack("=" + (size if version == "4.1" else "i") * 6, *pairs)
    return data + b"\n$EndPeriodic\n"


def assert_periodic(mesh):
    assert len(mesh.gmsh_periodic) == 2
    first, link = mesh.gmsh_periodic
    assert first[0] == 1 and first[1] == (11, 33)
    assert first[2] is None or len(first[2]) == 0
    assert first[3].shape == (0, 2)
    assert link[0] == 1 and link[1] == (12, 33)
    affine = np.eye(4).ravel()
    affine[3] = 1
    np.testing.assert_array_equal(link[2], affine)
    np.testing.assert_array_equal(link[3], [[0, 1], [2, 3], [0, 1]])


@pytest.mark.parametrize("version", ["2.2", "4.0", "4.1"])
def test_native_periodic_selective_and_metadata(version, tmp_path):
    from meshioplusplus import _core

    path = tmp_path / "periodic.msh"
    path.write_bytes(periodic_fixture(version))
    mesh = _core.gmsh_read(str(path), points_only=True, arrays=[])
    assert not mesh.cells
    assert_periodic(mesh)
    metadata = meshioplusplus.read_metadata(path, file_format="gmsh")
    assert metadata["num_points"] == 4
    assert metadata["fell_back_to_full_read"]


@pytest.mark.parametrize("version", ["2.2", "4.0", "4.1"])
def test_native_periodic_time_step(version, tmp_path, periodic_native_only):
    from meshioplusplus.gmsh.main import write as reference_write

    path = tmp_path / "transient.msh"
    mesh = gmsh_periodic()
    reference_write(path, mesh, fmt_version=version, binary=False)
    with path.open("a") as f:
        for step in [0, 1]:
            f.write(
                f'$NodeData\n1\n"signal"\n1\n{step}\n3\n{step}\n1\n1\n1 {step + 10}\n$EndNodeData\n'
            )
    back = meshioplusplus.read(path, file_format="gmsh", time_step=-1)
    assert back.point_data["signal"][0] == 11
    for a, b in zip(mesh.gmsh_periodic, back.gmsh_periodic):
        np.testing.assert_array_equal(a[3], b[3])


@pytest.mark.parametrize("version", ["2.2", "4.0", "4.1"])
@pytest.mark.parametrize("binary", [False, True])
@pytest.mark.parametrize("width", [4, 8])
def test_native_periodic_sparse_layouts(
    version, binary, width, tmp_path, periodic_native_only
):
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as reference_read

    path = tmp_path / "periodic.msh"
    path.write_bytes(periodic_fixture(version, binary, width))
    assert_periodic(_core.gmsh_read(str(path)))
    assert_periodic(meshioplusplus.read(path, file_format="gmsh"))
    # The 4.0 reference uses host long; only its host-width binary is supported.
    if version != "4.0" or not binary or width == np.dtype("L").itemsize:
        assert_periodic(reference_read(path))


@pytest.mark.parametrize("version", ["2.2", "4.1"])
@pytest.mark.parametrize("binary", [False, True])
def test_native_periodic_cross_writers(version, binary, tmp_path, periodic_native_only):
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as reference_read
    from meshioplusplus.gmsh.main import write as reference_write

    path = tmp_path / "source.msh"
    path.write_bytes(periodic_fixture(version))
    mesh = _core.gmsh_read(str(path))
    mesh.cell_data["gmsh:physical"] = [
        np.full(len(c.data), 7, dtype=np.int32) for c in mesh.cells
    ]
    for writer in [reference_write, meshioplusplus.gmsh.write]:
        output = tmp_path / "output.msh"
        writer(output, mesh, fmt_version=version, binary=binary)
        assert_periodic(_core.gmsh_read(str(output)))
        assert_periodic(reference_read(output))
    # Metadata owns its buffers independently of the source mesh lifetime.
    link = mesh.gmsh_periodic[1]
    del mesh
    np.testing.assert_array_equal(link[3][-1], [0, 1])


@pytest.mark.parametrize("binary", [False, True])
def test_native_periodic41_grouped_node_order(binary, tmp_path):
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as reference_read

    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        [("line", [[0, 2]]), ("line", [[1, 3]])],
        point_data={"gmsh:dim_tags": np.array([[1, 1], [1, 2], [1, 1], [1, 2]])},
        cell_data={
            "gmsh:geometrical": [np.array([1]), np.array([2])],
            "gmsh:physical": [np.array([7]), np.array([8])],
        },
        gmsh_periodic=[[1, (2, 1), np.eye(4).ravel(), np.array([[1, 0], [3, 2]])]],
    )
    path = tmp_path / "grouped.msh"
    _core.gmsh41_write(str(path), mesh, binary, None)
    for reader in [_core.gmsh_read, reference_read]:
        back = reader(str(path))
        # Writer groups nodes by entity. Rows change, but paired coordinates do not.
        pairs = back.gmsh_periodic[0][3]
        np.testing.assert_array_equal(
            back.points[pairs], mesh.points[mesh.gmsh_periodic[0][3]]
        )
    np.testing.assert_array_equal(mesh.gmsh_periodic[0][3], [[1, 0], [3, 2]])


@pytest.mark.parametrize("version", ["2.2", "4.0", "4.1"])
@pytest.mark.parametrize("binary", [False, True])
def test_reference_periodic_writer_native_read(version, binary, tmp_path):
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import write as reference_write

    mesh = gmsh_periodic()
    path = tmp_path / "reference.msh"
    reference_write(path, mesh, fmt_version=version, binary=binary)
    back = _core.gmsh_read(str(path))
    for a, b in zip(mesh.gmsh_periodic, back.gmsh_periodic):
        np.testing.assert_array_equal(a[3], b[3])


@pytest.mark.parametrize("pairs", [[[99, 0]], [[-1, 0]], [[0.5, 1]], [0, 1]])
@pytest.mark.parametrize("version", ["2.2", "4.1"])
def test_native_periodic_writer_rejects_bad_indices(pairs, version, tmp_path):
    from meshioplusplus import _core

    mesh = gmsh_periodic()
    mesh.gmsh_periodic[0][3] = pairs
    path = str(tmp_path / "invalid.msh")
    with pytest.raises((ValueError, meshioplusplus.WriteError)):
        if version == "2.2":
            _core.gmsh22_write(path, mesh, False)
        else:
            _core.gmsh41_write(path, mesh, False, None)


@pytest.mark.parametrize("version", ["2.2", "4.0", "4.1"])
def test_periodic_missing_terminator_and_duplicates(version, tmp_path):
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as reference_read

    data = periodic_fixture(version)
    section = data[data.index(b"$Periodic") :]
    path = tmp_path / "invalid.msh"
    for invalid in [data.replace(b"$EndPeriodic", b"$Other"), data + section]:
        path.write_bytes(invalid)
        for reader in [_core.gmsh_read, reference_read]:
            with pytest.raises(meshioplusplus.ReadError):
                reader(str(path))


@pytest.mark.parametrize("version", ["2.2", "4.0", "4.1"])
def test_public_periodic_writer_never_falls_back_to_invalid_metadata(version, tmp_path):
    mesh = gmsh_periodic()
    mesh.gmsh_periodic[0][3] = [[99, 0]]
    path = tmp_path / "unchanged.msh"
    path.write_text("unchanged")
    with pytest.raises(meshioplusplus.WriteError):
        meshioplusplus.gmsh.write(path, mesh, fmt_version=version)
    assert path.read_text() == "unchanged"


@pytest.mark.parametrize(
    "mesh",
    [
        # helpers.empty_mesh,
        helpers.line_mesh,
        helpers.tri_mesh,
        helpers.triangle6_mesh,
        helpers.quad_mesh,
        helpers.quad8_mesh,
        # helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
        helpers.pyramid13_mesh,
        helpers.wedge15_mesh,
        helpers.pyramid14_mesh,
        helpers.wedge18_mesh,
        helpers.add_point_data(helpers.tri_mesh, 1),
        helpers.add_point_data(helpers.tri_mesh, 3),
        helpers.add_point_data(helpers.tri_mesh, 9),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (3,), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (9,), np.float64)]),
        helpers.add_field_data(helpers.tri_mesh, [1, 2], int),
        helpers.add_field_data(helpers.tet_mesh, [1, 3], int),
        helpers.add_field_data(helpers.hex_mesh, [1, 3], int),
        gmsh_periodic(),
    ],
)
@pytest.mark.parametrize("binary", [False, True])
def test_gmsh22(mesh, binary, tmp_path):
    writer = partial(meshioplusplus.gmsh.write, fmt_version="2.2", binary=binary)
    helpers.write_read(tmp_path, writer, meshioplusplus.gmsh.read, mesh, 1.0e-15)


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tri_mesh,
        helpers.triangle6_mesh,
        helpers.quad_mesh,
        helpers.quad8_mesh,
        # helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
        helpers.add_point_data(helpers.tri_mesh, 1),
        helpers.add_point_data(helpers.tri_mesh, 3),
        helpers.add_point_data(helpers.tri_mesh, 9),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (3,), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (9,), np.float64)]),
        helpers.add_field_data(helpers.tri_mesh, [1, 2], int),
        helpers.add_field_data(helpers.tet_mesh, [1, 3], int),
        helpers.add_field_data(helpers.hex_mesh, [1, 3], int),
    ],
)
@pytest.mark.parametrize("binary", [False, True])
def test_gmsh40(mesh, binary, tmp_path):
    writer = partial(meshioplusplus.gmsh.write, fmt_version="4.0", binary=binary)

    helpers.write_read(tmp_path, writer, meshioplusplus.gmsh.read, mesh, 1.0e-15)


def gmsh40_fixture(binary=False, count_bytes=8):
    """Independent 4.0 records with sparse tags and three entity blocks."""
    header = f"$MeshFormat\n4.0 {int(binary)} 8\n".encode()
    if binary:
        header += struct.pack("=i", 1) + b"\n"
    header += (
        b'$EndMeshFormat\n$PhysicalNames\n2\n1 8 "edge"\n'
        b'2 7 "plate"\n$EndPhysicalNames\n'
    )
    body = bytearray()

    def put(values, dtype):
        if binary:
            code = {
                "int": "i",
                "count": "I" if count_bytes == 4 else "Q",
                "double": "d",
            }[dtype]
            body.extend(struct.pack("=" + code * len(values), *values))
        else:
            body.extend((" ".join(map(str, values)) + "\n").encode())

    body.extend(b"$Entities\n")
    put([1, 2, 1, 0], "count")
    for tag in [5, 11, 33, 22]:
        put([tag], "int")
        put([0.0] * 6, "double")  # even point entities have six bbox doubles
        put([0 if tag == 5 else 1], "count")
        if tag != 5:
            put([7 if tag == 22 else 8], "int")
            put([0 if tag == 33 else 2], "count")
            if tag != 33:
                put([5, -6] if tag == 11 else [11, -33], "int")
    body.extend(b"\n$EndEntities\n$Nodes\n")
    put([2, 4], "count")
    nodes = [(30, [1, 0, 0]), (10, [0, 0, 0]), (40, [1, 1, 0]), (20, [0, 1, 0])]
    for block, (tag, dim) in enumerate([(11, 1), (22, 2)]):
        if binary:
            put([tag, dim, 0], "int")
            put([2], "count")
        else:
            put([tag, dim, 0, 2], "int")
        for node, xyz in nodes[block * 2 : block * 2 + 2]:
            if binary:
                put([node], "int")
                put(xyz, "double")
            else:
                put([node, *xyz], "int")
    body.extend(b"\n$EndNodes\n$Elements\n")
    put([3, 4], "count")
    for tag, dim, kind, rows in [
        (11, 1, 1, [[90, 10, 30]]),
        (22, 2, 2, [[80, 10, 30, 40], [60, 10, 40, 20]]),
        (33, 1, 1, [[50, 20, 10]]),
    ]:
        put([tag, dim, kind], "int")
        put([len(rows)], "count")
        for row in rows:
            put(row, "int")
    body.extend(b"\n$EndElements\n")
    return header + body


@pytest.mark.parametrize("binary,count_bytes", [(False, 8), (True, 8), (True, 4)])
def test_gmsh40_native_entities_and_flat_count_widths(
    binary, count_bytes, tmp_path, monkeypatch
):
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as py_read

    path = tmp_path / "entities40.msh"
    path.write_bytes(gmsh40_fixture(binary, count_bytes))
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    mesh = meshioplusplus.gmsh.read(path)
    assert [c.type for c in mesh.cells] == ["line", "triangle", "line"]
    np.testing.assert_array_equal(mesh.cells[0].data, [[1, 0]])
    np.testing.assert_array_equal(mesh.cells[1].data, [[1, 0, 2], [1, 2, 3]])
    np.testing.assert_array_equal(mesh.cells[2].data, [[3, 1]])
    assert "gmsh:dim_tags" not in mesh.point_data
    assert "gmsh:bounding_entities" not in mesh.cell_sets
    assert sorted((r.name, r.dim, r.tag, list(r.entries)) for r in mesh.regions) == [
        ("edge", 1, 8, [0, 3]),
        ("plate", 2, 7, [1, 2]),
    ]
    # The reference's unsigned-long width is host-dependent. The core also
    # accepts files from 32-bit producers, including when running in WASM.
    if not binary or count_bytes == np.dtype("L").itemsize:
        # The legacy 4.0 reference keeps physical tags as data but does not
        # derive cell sets/regions. The native region mapping is checked above.
        _assert_same_mesh(mesh, py_read(path), compare_sets=False)
    _assert_same_mesh(mesh, _core.gmsh_read(str(path)))
    meta = meshioplusplus.read_metadata(path, "gmsh")
    assert meta["num_points"] == 4
    assert meta["fell_back_to_full_read"] is True
    _assert_same_mesh(mesh, meshioplusplus.read(path, file_format="gmsh"))


@pytest.mark.parametrize("binary", [False, True])
@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tet10_mesh,
        helpers.hex20_mesh,
        meshioplusplus.Mesh(
            np.arange(81, dtype=float).reshape(27, 3),
            [("hexahedron27", np.arange(27).reshape(1, 27))],
        ),
        helpers.wedge15_mesh,
        helpers.wedge18_mesh,
        helpers.pyramid13_mesh,
        helpers.pyramid14_mesh,
        helpers.add_point_data(helpers.tri_mesh, 3),
        helpers.add_cell_data(helpers.tri_mesh, [("stress", (3,), np.float64)]),
    ],
)
def test_gmsh40_reference_output_matches_core(mesh, binary, tmp_path):
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as py_read
    from meshioplusplus.gmsh.main import write as py_write

    path = tmp_path / "reference40.msh"
    py_write(path, mesh, fmt_version="4.0", binary=binary)
    _assert_same_mesh(_core.gmsh_read(str(path)), py_read(path))


@pytest.mark.parametrize("binary", [False, True])
def test_gmsh40_selective_and_transient_reads(binary, tmp_path):
    from meshioplusplus import _core
    from meshioplusplus.gmsh.common import _write_data
    from meshioplusplus.gmsh.main import write as py_write

    path = tmp_path / "data40.msh"
    mesh = helpers.add_point_data(helpers.tri_mesh, 1)
    py_write(path, mesh, fmt_version="4.0", binary=binary)
    with path.open("ab") as stream:
        _write_data(stream, "NodeData", "later", np.arange(len(mesh.points)), binary)
    none = _core.gmsh_read(str(path), points_only=True)
    assert not none.point_data and not none.cell_data
    selected = _core.gmsh_read(str(path), arrays=["later"])
    assert sorted(selected.point_data) == ["later"]
    assert not selected.cell_data
    np.testing.assert_array_equal(
        selected.point_data["later"], np.arange(len(mesh.points))
    )
    with path.open("ab") as stream:
        stream.write(
            f'$NodeData\n1\n"a"\n1\n2.5\n3\n1\n1\n{len(mesh.points)}\n'.encode()
        )
        for i in range(len(mesh.points)):
            stream.write(
                struct.pack("=id", i + 1, 10.0 + i)
                if binary
                else f"{i + 1} {10 + i}\n".encode()
            )
        stream.write(b"\n$EndNodeData\n")
    later = meshioplusplus.gmsh.read(path, time_step=-1)
    np.testing.assert_array_equal(
        later.point_data["a"], 10 + np.arange(len(mesh.points))
    )
    assert "later" not in later.point_data
    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.gmsh.read(path, time_step=2)


@pytest.mark.parametrize("count_bytes", [4, 8])
def test_gmsh40_truncated_binary_sections_raise(count_bytes, tmp_path):
    from meshioplusplus import _core

    data = gmsh40_fixture(True, count_bytes)
    path = tmp_path / "bad40.msh"
    for marker in [b"$EndEntities", b"$EndNodes", b"$EndElements"]:
        path.write_bytes(data[: data.index(marker) - 3])
        with pytest.raises(meshioplusplus.ReadError):
            _core.gmsh_read(str(path))


@pytest.mark.parametrize("count_bytes", [4, 8])
def test_gmsh40_bad_binary_header_raises(count_bytes, tmp_path):
    from meshioplusplus import _core

    path = tmp_path / "header40.msh"
    data = gmsh40_fixture(True, count_bytes)
    marker = len(b"$MeshFormat\n4.0 1 8\n")
    path.write_bytes(data[:marker] + struct.pack(">i", 1) + data[marker + 4 :])
    with pytest.raises(meshioplusplus.ReadError, match="endianness"):
        _core.gmsh_read(str(path))
    path.write_bytes(data.replace(b"4.0 1 8", b"4.0 1 4", 1))
    with pytest.raises(meshioplusplus.ReadError, match="8-byte doubles"):
        _core.gmsh_read(str(path))


def test_gmsh40_partial_physical_tags_keep_block_alignment(tmp_path):
    from meshioplusplus import _core

    text = gmsh40_fixture().decode()
    # Third curve: no physical tag, even though the other entities are tagged.
    text = text.replace(
        "33\n0.0 0.0 0.0 0.0 0.0 0.0\n1\n8\n0\n", "33\n0.0 0.0 0.0 0.0 0.0 0.0\n0\n0\n"
    )
    path = tmp_path / "partial40.msh"
    path.write_text(text)
    mesh = _core.gmsh_read(str(path))
    assert [a.tolist() for a in mesh.cell_data["gmsh:physical"]] == [[8], [7, 7], [0]]
    assert sorted((r.name, list(r.entries)) for r in mesh.regions) == [
        ("edge", [0]),
        ("plate", [1, 2]),
    ]


@pytest.mark.parametrize("binary,count_bytes", [(False, 8), (True, 4), (True, 8)])
def test_gmsh40_empty_mesh(binary, count_bytes, tmp_path):
    from meshioplusplus import _core

    path = tmp_path / "empty40.msh"
    data = f"$MeshFormat\n4.0 {int(binary)} 8\n".encode()
    if binary:
        data += struct.pack("=i", 1) + b"\n"
    data += b"$EndMeshFormat\n$Nodes\n"
    counts = b"0 0\n" if not binary else bytes(count_bytes * 2)
    data += counts + b"\n$EndNodes\n$Elements\n" + counts + b"\n$EndElements\n"
    path.write_bytes(data)
    mesh = _core.gmsh_read(str(path))
    assert mesh.points.shape == (0, 3)
    assert not mesh.cells
    assert not mesh.point_data and not mesh.cell_data


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tri_mesh,
        helpers.triangle6_mesh,
        helpers.quad_mesh,
        helpers.quad8_mesh,
        # helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
        helpers.add_point_data(helpers.tri_mesh, 1),
        helpers.add_point_data(helpers.tri_mesh, 3),
        helpers.add_point_data(helpers.tri_mesh, 9),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (3,), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (9,), np.float64)]),
        helpers.add_field_data(helpers.tri_mesh, [1, 2], int),
        helpers.add_field_data(helpers.tet_mesh, [1, 3], int),
        helpers.add_field_data(helpers.hex_mesh, [1, 3], int),
        gmsh_periodic(),
    ],
)
@pytest.mark.parametrize("binary", [False, True])
def test_gmsh41(mesh, binary, tmp_path):
    writer = partial(meshioplusplus.gmsh.write, fmt_version="4.1", binary=binary)
    helpers.write_read(tmp_path, writer, meshioplusplus.gmsh.read, mesh, 1.0e-15)


def test_generic_io(tmp_path):
    """A `.msh` written without a format and without gmsh tags is Fluent's
    (the first `.msh` candidate): it reads back as the same planar triangles,
    rebuilt from their faces, so compare them as point sets."""
    mesh = helpers.tri_mesh
    # With additional, insignificant suffix too:
    for path in (tmp_path / "test.msh", tmp_path / "test.0.msh"):
        meshioplusplus.write_points_cells(path, mesh.points, mesh.cells)
        out = meshioplusplus.read(path)
        dim = out.points.shape[1]
        assert (mesh.points[:, dim:] == 0).all()

        def triangles(m, pts):
            return sorted(
                tuple(sorted(map(tuple, pts[row].tolist())))
                for block in m.cells
                if block.type == "triangle"
                for row in block.data
            )

        assert triangles(out, out.points) == triangles(mesh, mesh.points[:, :dim])


@pytest.mark.parametrize(
    "filename, ref_sum, ref_num_cells",
    [("insulated-2.2.msh", 2.001762136876221, [21, 111])],
)
@pytest.mark.parametrize("binary", [False, True])
def test_reference_file(filename, ref_sum, ref_num_cells, binary, tmp_path):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "msh" / filename
    mesh = meshioplusplus.read(filename)
    tol = 1.0e-2
    s = mesh.points.sum()
    assert abs(s - ref_sum) < tol * ref_sum
    assert [c.type for c in mesh.cells] == ["line", "triangle"]
    assert [len(c.data) for c in mesh.cells] == ref_num_cells
    assert list(map(len, mesh.cell_data["gmsh:geometrical"])) == ref_num_cells
    assert list(map(len, mesh.cell_data["gmsh:physical"])) == ref_num_cells

    writer = partial(meshioplusplus.gmsh.write, fmt_version="2.2", binary=binary)
    helpers.write_read(tmp_path, writer, meshioplusplus.gmsh.read, mesh, 1.0e-15)


@pytest.mark.parametrize(
    "filename, ref_sum, ref_num_cells, ref_num_cells_in_cell_sets",
    [
        (
            "insulated-4.1.msh",
            2.001762136876221,
            {"line": 21, "triangle": 111},
            {"line": 27, "triangle": 120},
        )
    ],
    # Note that testing on number of cells in
    # cell_sets_dict will count both cells associated with physical tags, and
    # bounding entities.
)
@pytest.mark.parametrize("binary", [False, True])
def test_reference_file_with_entities(
    filename, ref_sum, ref_num_cells, ref_num_cells_in_cell_sets, binary, tmp_path
):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "msh" / filename

    mesh = meshioplusplus.read(filename)
    tol = 1.0e-2
    s = mesh.points.sum()
    assert abs(s - ref_sum) < tol * ref_sum
    assert {k: len(v) for k, v in mesh.cells_dict.items()} == ref_num_cells
    assert {
        k: len(v) for k, v in mesh.cell_data_dict["gmsh:physical"].items()
    } == ref_num_cells

    writer = partial(meshioplusplus.gmsh.write, fmt_version="4.1", binary=binary)

    num_cells = {k: 0 for k in ref_num_cells_in_cell_sets}
    for vv in mesh.cell_sets_dict.values():
        for k, v in vv.items():
            num_cells[k] += len(v)
    assert num_cells == ref_num_cells_in_cell_sets

    # $Entities is the only place 4.1 records physical-group membership, so the
    # regions it yields carry the group's real dimension and tag -- not the -1
    # placeholders a set-derived region has.
    assert sorted((r.name, r.dim, r.tag, len(r.entries)) for r in mesh.regions) == [
        ("convection", 1, 3, 21),
        ("insulation", 2, 2, 66),
        ("wire", 2, 1, 45),
    ]

    helpers.write_read(tmp_path, writer, meshioplusplus.gmsh.read, mesh, 1.0e-15)


_ENTITY_FILES = ["tests/python/meshes/msh/insulated-4.1.msh", "example/example.msh"]


def _assert_same_mesh(a, b, compare_sets=True):
    np.testing.assert_allclose(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for ca, cb in zip(a.cells, b.cells):
        np.testing.assert_array_equal(ca.data, cb.data)
    assert sorted(a.point_data) == sorted(b.point_data)
    for k in a.point_data:
        np.testing.assert_array_equal(a.point_data[k], b.point_data[k])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for k in a.cell_data:
        for x, y in zip(a.cell_data[k], b.cell_data[k]):
            np.testing.assert_array_equal(x, y)
    assert sorted(a.field_data) == sorted(b.field_data)
    for k in a.field_data:
        np.testing.assert_array_equal(a.field_data[k], b.field_data[k])
    if not compare_sets:
        return
    assert sorted(a.cell_sets) == sorted(b.cell_sets)
    for k in a.cell_sets:
        for x, y in zip(a.cell_sets[k], b.cell_sets[k]):
            xa = np.asarray([] if x is None else x).ravel()
            ya = np.asarray([] if y is None else y).ravel()
            np.testing.assert_array_equal(xa, ya)


@pytest.mark.parametrize("filename", _ENTITY_FILES)
def test_cpp_matches_python_on_entities(filename):
    # Before $Entities landed in the C++ reader these files reached Python only
    # by way of the shim's blanket except -- so this is the gate that the two
    # readers now genuinely agree rather than one of them being unreachable.
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as py_read

    root = pathlib.Path(__file__).resolve().parents[2]
    path = str(root / filename)
    _assert_same_mesh(_core.gmsh_read(path), py_read(path))


@pytest.mark.parametrize("filename", _ENTITY_FILES)
@pytest.mark.parametrize("binary", [False, True])
def test_entities_survive_a_cpp_round_trip(filename, binary, tmp_path):
    from meshioplusplus import _core

    root = pathlib.Path(__file__).resolve().parents[2]
    mesh = _core.gmsh_read(str(root / filename))
    out = str(tmp_path / "rt.msh")
    _core.gmsh41_write(out, mesh, binary, mesh.cell_sets.get("gmsh:bounding_entities"))
    _assert_same_mesh(mesh, _core.gmsh_read(out))


def test_cpp_ascii_41_output_matches_the_python_writer(tmp_path):
    # Byte parity, not just equivalence: the entity records' spacing and line
    # discipline are easy to get subtly wrong and no read-back would notice.
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as py_read
    from meshioplusplus.gmsh.main import write as py_write

    root = pathlib.Path(__file__).resolve().parents[2]
    src = str(root / "tests/python/meshes/msh/insulated-4.1.msh")
    mesh = py_read(src)

    py_path = tmp_path / "py.msh"
    cpp_path = tmp_path / "cpp.msh"
    py_write(str(py_path), mesh, fmt_version="4.1", binary=False)
    _core.gmsh41_write(
        str(cpp_path),
        _core.gmsh_read(src),
        False,
        mesh.cell_sets.get("gmsh:bounding_entities"),
    )
    assert py_path.read_bytes() == cpp_path.read_bytes()


def test_writing_41_without_dim_tags_emits_no_entities(tmp_path):
    # No entity information to describe -> the legacy single-$Nodes-block
    # output, unchanged.
    from meshioplusplus import _core

    mesh = copy.deepcopy(helpers.tet_mesh)
    path = tmp_path / "plain.msh"
    _core.gmsh41_write(str(path), mesh, False)
    text = path.read_text()
    assert "$Entities" not in text
    # A single $Nodes block covering every point, as before.
    assert text.split("$Nodes\n")[1].splitlines()[0].split()[0] == "1"


def test_a_41_file_with_no_physical_groups_has_no_physical_cell_data():
    # example.msh tags no entity at all. Synthesizing an all-zero
    # gmsh:physical there would invent a group the file does not have.
    from meshioplusplus import _core

    root = pathlib.Path(__file__).resolve().parents[2]
    mesh = _core.gmsh_read(str(root / "example/example.msh"))
    assert "gmsh:physical" not in mesh.cell_data
    assert mesh.regions == []
    assert "gmsh:geometrical" in mesh.cell_data


def test_read_metadata_reports_both_steps_of_a_transient_file(tmp_path):
    """roadmap §1 tier B1: read_gmsh_metadata (4.1) reports the sorted union
    of every $NodeData/$ElementData section's time value, and ``time_step``
    on a real read picks exactly the one matching that resolved time --
    replacing the pre-v11.3.0 "first section per name wins" rule. A second
    step is appended by hand (ASCII text, no external tool)."""
    from meshioplusplus import _core

    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["u"] = np.array([10.0, 20.0, 30.0, 40.0])
    path = tmp_path / "transient.msh"
    _core.gmsh41_write(str(path), mesh, False)

    n = len(mesh.points)
    lines = [
        "$NodeData",
        "1",
        '"u"',
        "1",
        "2.5",
        "3",
        "1",
        "1",
        str(n),
    ]
    lines += [f"{i + 1} {11.0 + i * 10.0}" for i in range(n)]
    lines.append("$EndNodeData")
    with open(path, "a") as f:
        f.write("\n".join(lines) + "\n")

    meta = meshioplusplus.read_metadata(path, "gmsh")
    assert meta["fell_back_to_full_read"] is False
    assert meta["time_values"] == [0.0, 2.5]
    assert meta["point_data_names"].count("u") == 1

    mesh0 = meshioplusplus.gmsh.read(path, time_step=0)
    assert mesh0.point_data["u"][0] == 10.0
    mesh1 = meshioplusplus.gmsh.read(path, time_step=1)
    assert mesh1.point_data["u"][0] == 11.0
    mesh_last = meshioplusplus.gmsh.read(path, time_step=-1)
    assert mesh_last.point_data["u"][0] == 11.0

    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.gmsh.read(path, time_step=5)


def test_malformed_time_value_raises_read_error_not_value_error(tmp_path):
    """A $NodeData section whose time value isn't a valid number used to
    reach std::stod directly, propagating as a bare ValueError instead of a
    proper ReadError -- part of the roadmap's locale item (std::stod's
    exceptions are the reason A4 gave these two call sites an explicit
    check rather than a blind swap to the non-throwing parse_double)."""
    from meshioplusplus import _core

    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["u"] = np.array([10.0, 20.0, 30.0, 40.0])
    path = tmp_path / "bad_time.msh"
    _core.gmsh41_write(str(path), mesh, False)

    n = len(mesh.points)
    lines = [
        "$NodeData",
        "1",
        '"u"',
        "1",
        "not-a-number",
        "3",
        "1",
        "1",
        str(n),
    ]
    lines += [f"{i + 1} {11.0 + i * 10.0}" for i in range(n)]
    lines.append("$EndNodeData")
    with open(path, "a") as f:
        f.write("\n".join(lines) + "\n")

    # Straight through the C++ core, not the meshioplusplus.gmsh.read() shim:
    # the shim's `except Exception: if time_step: raise` only re-raises for a
    # truthy time_step, and 0 (the default/first step) is falsy, so it would
    # otherwise swallow this and silently fall back to the Python reader --
    # a separate, already-tracked gap (roadmap §1's other reader-fallback
    # item), not what this test is about.
    with pytest.raises(meshioplusplus.ReadError, match="time value"):
        _core.gmsh_read(str(path), time_step=0)


def test_untagged_region_gets_an_allocated_tag_on_write(tmp_path):
    """Roadmap §1 tier B3 (v11.5.0): a Cell region with no gmsh tag of its
    own (as Abaqus/MED/MDPA produce) gets a freshly allocated one instead of
    being dropped from $PhysicalNames/gmsh:physical. C++-core only -- the
    pure-Python fallback writer keeps its existing field_data-only behaviour.
    """
    mesh = meshioplusplus.Mesh(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        [("hexahedron", [[0, 1, 2, 3, 4, 5, 6, 7]])],
        regions=[meshioplusplus.Region("body", "cell", [0])],
    )
    path = tmp_path / "body.msh"
    meshioplusplus.write(path, mesh, file_format="gmsh22")
    back = meshioplusplus.read(path)
    cell_regions = {r.name: r for r in back.regions if r.kind == "cell"}
    assert "body" in cell_regions
    assert cell_regions["body"].tag >= 0, "an allocated tag must not round-trip as -1"
    assert list(cell_regions["body"].entries) == [0]


# --- roadmap §1 "gmsh pyramid14 (and wedge18) node ordering is not
# permuted" ---
#
# A write->read round trip cannot catch a wrong-but-consistently-inverse
# permutation pair, so this checks against gmsh's OWN edge/face tables
# independently: `_meshio_to_gmsh_order` tells us which meshio corner sits
# at each gmsh slot; gmsh's corners 0..k_corner-1 are unchanged from
# meshio's (confirmed against gmsh's own src/geo/MPyramid.h /
# src/geo/MPrism.h edge/face tables), so re-deriving each mid-edge/
# face-centre point from THOSE corners and comparing to what actually
# landed there is a check independent of meshio's own edge/face
# convention (which is what the fixture meshes in helpers.py were built
# from).


def _assert_gmsh_order_matches_gmsh_geometry(
    cell_type, points, num_corners, edges, faces
):
    from meshioplusplus.gmsh.common import _meshio_to_gmsh_order

    idx = np.arange(len(points)).reshape(1, -1)
    gmsh_conn = _meshio_to_gmsh_order(cell_type, idx)[0]
    gmsh_pts = points[gmsh_conn]

    for slot, (a, b) in enumerate(edges, start=num_corners):
        expected = (gmsh_pts[a] + gmsh_pts[b]) / 2.0
        np.testing.assert_allclose(
            gmsh_pts[slot],
            expected,
            atol=1e-12,
            err_msg=f"{cell_type}: gmsh mid-edge slot {slot} (edge {a},{b})",
        )
    edge_end = num_corners + len(edges)
    for slot, corners in enumerate(faces, start=edge_end):
        expected = gmsh_pts[list(corners)].mean(axis=0)
        np.testing.assert_allclose(
            gmsh_pts[slot],
            expected,
            atol=1e-12,
            err_msg=f"{cell_type}: gmsh face-centre slot {slot} (face {corners})",
        )


def test_pyramid14_gmsh_order_matches_gmsh_own_edge_and_face_tables():
    _assert_gmsh_order_matches_gmsh_geometry(
        "pyramid14",
        helpers.pyramid14_mesh.points,
        num_corners=5,
        edges=[(0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (2, 4), (3, 4)],
        faces=[(0, 3, 2, 1)],
    )


def test_wedge18_gmsh_order_matches_gmsh_own_edge_and_face_tables():
    _assert_gmsh_order_matches_gmsh_geometry(
        "wedge18",
        helpers.wedge18_mesh.points,
        num_corners=6,
        edges=[
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 2),
            (1, 4),
            (2, 5),
            (3, 4),
            (3, 5),
            (4, 5),
        ],
        faces=[(0, 1, 4, 3), (0, 3, 5, 2), (1, 2, 5, 4)],
    )


def test_pyramid13_gmsh_order_matches_gmsh_own_edge_table():
    # The already-shipped pyramid13 table, as an in-repo positive control
    # for the helper above.
    _assert_gmsh_order_matches_gmsh_geometry(
        "pyramid13",
        helpers.pyramid13_mesh.points,
        num_corners=5,
        edges=[(0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (2, 4), (3, 4)],
        faces=[],
    )


def test_wedge15_gmsh_order_matches_gmsh_own_edge_table():
    _assert_gmsh_order_matches_gmsh_geometry(
        "wedge15",
        helpers.wedge15_mesh.points,
        num_corners=6,
        edges=[
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 2),
            (1, 4),
            (2, 5),
            (3, 4),
            (3, 5),
            (4, 5),
        ],
        faces=[],
    )


@pytest.mark.parametrize("engine", ["core", "python"])
def test_indented_file_reads(engine):
    # FEconv's samples indent every line (tools/gen_feconv_quirk_fixtures.py);
    # the third element line is past the declared count and ignored.
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as py_read

    path = pathlib.Path(__file__).resolve().parent / "meshes" / "gmsh" / "indented.msh"
    mesh = _core.gmsh_read(str(path)) if engine == "core" else py_read(path)
    assert len(mesh.points) == 6
    assert [c.type for c in mesh.cells] == ["quad"]
    np.testing.assert_array_equal(mesh.cells[0].data, [[0, 1, 2, 3], [1, 4, 5, 2]])
