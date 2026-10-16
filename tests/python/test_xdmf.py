import pathlib

import numpy as np
import pytest

import meshioplusplus

from . import helpers

test_set_full = [
    helpers.empty_mesh,
    helpers.line_mesh,
    helpers.tri_mesh,
    helpers.line_tri_mesh,
    helpers.tri_mesh_2d,
    helpers.triangle6_mesh,
    helpers.quad_mesh,
    helpers.quad8_mesh,
    helpers.tri_quad_mesh,
    helpers.tet_mesh,
    helpers.tet10_mesh,
    helpers.hex_mesh,
    helpers.hex20_mesh,
    helpers.pyramid_mesh,
    helpers.wedge_mesh,
    helpers.add_point_data(helpers.tri_mesh, 1),
    helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
]

test_set_reduced = [
    helpers.tri_mesh,
    helpers.tri_mesh_2d,
    helpers.quad_mesh,
    helpers.tri_quad_mesh,
    helpers.tet_mesh,
    helpers.hex_mesh,
    helpers.add_point_data(helpers.tri_mesh, 1),
    helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
]


@pytest.mark.parametrize("mesh", test_set_full)
@pytest.mark.parametrize(
    "kwargs0",
    [
        {"data_format": "XML"},
        {"data_format": "Binary"},
        {"data_format": "HDF", "compression": None},
        {"data_format": "HDF", "compression": "gzip"},
    ],
)
def test_xdmf3(mesh, kwargs0, tmp_path):
    def write(*args, **kwargs):
        return meshioplusplus.xdmf.write(*args, **{**kwargs0, **kwargs})

    helpers.write_read(tmp_path, write, meshioplusplus.xdmf.read, mesh, 1.0e-14)


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.xdmf")
    # With additional, insignificant suffix:
    helpers.generic_io(tmp_path / "test.0.xdmf")


@pytest.mark.parametrize("version", ["2.0", "3.0"])
@pytest.mark.parametrize(
    "reference", ["XML", "/Xdmf/Domain/Grid/Attribute[2]/DataItem"]
)
def test_native_xpath_references(version, reference, tmp_path, monkeypatch):
    core = pytest.importorskip("meshioplusplus._core")
    path = tmp_path / "references.xdmf"
    xpath = "/Xdmf/Domain/Grid/Attribute[2]/DataItem"
    path.write_text(
        f'<Xdmf Version="{version}"><Domain><Grid GridType="Uniform">'
        '<Topology TopologyType="Triangle"><DataItem Dimensions="1 3" '
        'NumberType="Int" Precision="8" Format="XML">0 1 2</DataItem></Topology>'
        f'<Geometry GeometryType="XYZ"><DataItem Reference="{reference}">'
        f'{xpath if reference == "XML" else ""}</DataItem></Geometry>'
        '<Attribute Name="source" Center="Node"><DataItem Dimensions="3 3" '
        'NumberType="Float" Precision="8" Format="XML">0 0 0 1 0 0 0 1 0</DataItem></Attribute>'
        '<Attribute Name="alias" Center="Node"><DataItem '
        'Reference="/Xdmf/Domain/Grid/Attribute[1]/DataItem"/></Attribute>'
        "</Grid></Domain></Xdmf>"
    )
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    native = meshioplusplus.read(path)
    reference_mesh = meshioplusplus.xdmf.main.read(path)
    np.testing.assert_array_equal(native.points, reference_mesh.points)
    np.testing.assert_array_equal(native.cells[0].data, reference_mesh.cells[0].data)
    for key in ("source", "alias"):
        np.testing.assert_array_equal(
            native.point_data[key], reference_mesh.point_data[key]
        )
    meta = meshioplusplus.read_metadata(path)
    assert meta["num_points"] == 3
    assert meta["point_dim"] == 3
    assert core.xdmf_read(str(path), arrays=[]).point_data == {}


def test_native_xdmf2_information(tmp_path, monkeypatch):
    core = pytest.importorskip("meshioplusplus._core")
    path = tmp_path / "information.xdmf"
    path.write_text(
        '<Xdmf Version="2.0"><Domain><Grid>'
        '<Topology TopologyType="Triangle"><DataItem Dimensions="1 3" '
        'NumberType="Int" Format="XML">0 1 2</DataItem></Topology>'
        '<Geometry GeometryType="XYZ"><DataItem Dimensions="3 3" '
        'NumberType="Float" Format="XML">0 0 0 1 0 0 0 1 0</DataItem></Geometry>'
        '<Information><![CDATA[<main><map key="wall" dim="2">7</map></main>]]></Information>'
        "</Grid></Domain></Xdmf>"
    )
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    mesh = meshioplusplus.read(path)
    np.testing.assert_array_equal(mesh.field_data["wall"], [7, 2])
    assert core.xdmf_read(str(path), points_only=True).field_data == {}
    assert core.xdmf_read(str(path), arrays=[]).field_data == {}
    assert "wall" in meshioplusplus.read_metadata(path)["field_data_names"]


@pytest.mark.parametrize(
    "reference",
    ["relative/path", "/Xdmf/Domain/Missing", "/Xdmf/Domain/Grid", "///["],
)
def test_native_invalid_xpath_is_read_error(tmp_path, reference):
    core = pytest.importorskip("meshioplusplus._core")
    path = tmp_path / "bad-reference.xdmf"
    path.write_text(
        '<Xdmf Version="3.0"><Domain><Grid><Geometry>'
        f'<DataItem Reference="{reference}"/></Geometry></Grid></Domain></Xdmf>'
    )
    with pytest.raises(Exception, match="XDMF:.*reference"):
        core.xdmf_read(str(path))


def test_native_cyclic_xpath_is_read_error(tmp_path):
    core = pytest.importorskip("meshioplusplus._core")
    path = tmp_path / "cycle.xdmf"
    path.write_text(
        '<Xdmf Version="3.0"><Domain><Grid><Geometry>'
        '<DataItem Reference="/Xdmf/Domain/Grid/Geometry/DataItem"/>'
        "</Geometry></Grid></Domain></Xdmf>"
    )
    with pytest.raises(Exception, match="cyclic DataItem reference"):
        core.xdmf_read(str(path))


def test_time_series():
    # write the data
    filename = "out.xdmf"

    with meshioplusplus.xdmf.TimeSeriesWriter(filename) as writer:
        writer.write_points_cells(helpers.tri_mesh_2d.points, helpers.tri_mesh_2d.cells)
        n = helpers.tri_mesh_2d.points.shape[0]

        times = np.linspace(0.0, 1.0, 5)
        point_data = [
            {
                "phi": np.full(n, t),
                "u": np.full(helpers.tri_mesh_2d.points.shape, t),
            }
            for t in times
        ]
        for t, pd in zip(times, point_data):
            writer.write_data(
                t, point_data=pd, cell_data={"a": {"triangle": [3.0, 4.2]}}
            )

    # read it back in
    with meshioplusplus.xdmf.TimeSeriesReader(filename) as reader:
        points, cells = reader.read_points_cells()
        for k in range(reader.num_steps):
            t, pd, cd = reader.read_data(k)
            assert np.abs(times[k] - t) < 1.0e-12
            for key, value in pd.items():
                assert np.all(np.abs(value - point_data[k][key]) < 1.0e-12)


# def test_information_xdmf():
#     mesh_out = meshioplusplus.Mesh(
#         np.array(
#             [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
#         )
#         / 3,
#         [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
#         field_data={
#             "bottom": np.array([1, 1]),
#             "right": np.array([2, 1]),
#             "top": np.array([3, 1]),
#             "left": np.array([4, 1]),
#         },
#     )
#     # write the data
#     points, cells, field_data = mesh_out.points, mesh_out.cells, mesh_out.field_data
#
#     assert cells[0].type == "triangle"
#     meshioplusplus.write(
#         "mesh.xdmf",
#         meshioplusplus.Mesh(points=points, cells=[cells[0]], field_data=field_data),
#     )
#
#     # read it back in
#     mesh_in = meshioplusplus.read("mesh.xdmf")
#     assert len(mesh_in.field_data) == len(mesh_out.field_data)


if __name__ == "__main__":
    test_time_series()


# --- malformed-input / error-path coverage ---


def test_xdmf_wrong_root_tag_raises(tmp_path):
    p = tmp_path / "bad.xdmf"
    p.write_text('<?xml version="1.0"?><Nonsense/>')
    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.read(p, file_format="xdmf")


@pytest.mark.parametrize("engine", ["core", "python"])
def test_mixed_topology_with_vertices_round_trips(tmp_path, engine):
    # Both writers give a Polyvertex its node count (XDMF requires it); both
    # readers skipped the count only for Polyline, so every mixed topology
    # holding a vertex was misread.
    from meshioplusplus.xdmf.main import XdmfReader

    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [
            ("vertex", [[3]]),
            ("line", [[0, 1]]),
            ("triangle", [[0, 1, 2]]),
            ("tetra", [[0, 1, 2, 3]]),
        ],
    )
    path = tmp_path / "mixed.xdmf"
    meshioplusplus.write(path, mesh, data_format="XML")
    if engine == "core":
        from meshioplusplus import _core

        back = _core.xdmf_read(str(path))
        blocks = [(b.type, np.asarray(b.data).tolist()) for b in back.cells]
    else:
        back = XdmfReader(path).read()
        blocks = [(b.type, np.asarray(b.data).tolist()) for b in back.cells]
    assert blocks == [(b.type, b.data.tolist()) for b in mesh.cells]


def test_truncated_mixed_topology_is_a_read_error():
    from meshioplusplus.xdmf.common import translate_mixed_cells

    with pytest.raises(meshioplusplus.ReadError, match="ends inside"):
        translate_mixed_cells(np.array([6, 0, 1, 2]))


def test_python_reader_reads_quadratic_cells_in_a_mixed_topology(tmp_path):
    # The Python twin's node-count table lacked the quadratic types the C++
    # reader has, so a mixed topology holding one failed with a KeyError on
    # any platform whose core has no HDF5 (the Windows CI legs).
    from meshioplusplus.xdmf.main import XdmfReader

    mesh = meshioplusplus.Mesh(
        np.array(
            [
                [0, 0, 0],
                [2, 0, 0],
                [2, 2, 0],
                [0, 2, 0],
                [1, 0, 0],
                [2, 1, 0],
                [1, 2, 0],
                [0, 1, 0],
            ],
            dtype=float,
        ),
        [("triangle", [[0, 1, 2]]), ("quad8", [[0, 1, 2, 3, 4, 5, 6, 7]])],
    )
    path = tmp_path / "quadratic.xdmf"
    meshioplusplus.write(path, mesh, data_format="XML")
    back = XdmfReader(path).read()
    assert [(b.type, np.asarray(b.data).tolist()) for b in back.cells] == [
        (b.type, b.data.tolist()) for b in mesh.cells
    ]


@pytest.mark.parametrize("engine", ["core", "python"])
def test_hdf_dataset_larger_than_its_dimensions_is_a_read_error(tmp_path, engine):
    # The file's own extent is not bounded by the XML: a corrupt dataset that
    # declares millions of unwritten chunks made HDF5 1.10's H5Dread crawl (the
    # fuzzer's timeout, tests/fuzz/regressions/xdmf). A dataset larger than its
    # DataItem's Dimensions is refused before it is read, by both engines.
    pytest.importorskip("h5py")
    from meshioplusplus.xdmf.main import XdmfReader

    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [("triangle", [[0, 1, 2]])],
    )
    path = tmp_path / "big.xdmf"
    meshioplusplus.write(path, mesh, data_format="HDF")
    text = path.read_text()
    geometry = text.index("<Geometry")
    start = text.index('Dimensions="', geometry) + len('Dimensions="')
    end = text.index('"', start)
    path.write_text(text[:start] + "1 1" + text[end:])  # stored: 3 x 3
    if engine == "core":
        from meshioplusplus import _core

        if not getattr(_core, "__has_hdf5__", True):
            pytest.skip("build has no HDF5")
        with pytest.raises(meshioplusplus.ReadError):
            _core.xdmf_read(str(path))
    else:
        with pytest.raises(meshioplusplus.ReadError):
            XdmfReader(path).read()


_TOPOLOGY = (
    '<Topology TopologyType="Triangle" NumberOfElements="1">'
    '<DataItem Dimensions="1 3" Format="XML">0 1 2</DataItem></Topology>'
)


def _xdmf_document(geometry_item):
    return (
        '<?xml version="1.0"?>\n<Xdmf Version="3.0"><Domain><Grid>\n'
        f'{_TOPOLOGY}\n<Geometry GeometryType="XY">{geometry_item}</Geometry>\n'
        "</Grid></Domain></Xdmf>\n"
    )


def _read_with_engine(engine, path):
    from meshioplusplus import _core
    from meshioplusplus.xdmf.main import XdmfReader

    if engine == "core":
        return _core.xdmf_read(str(path))
    return XdmfReader(path).read()


@pytest.mark.parametrize("engine", ["core", "python"])
@pytest.mark.parametrize(
    "geometry",
    [
        # A negative extent: the C++ reader cast it to a huge count, the Python
        # reshape read it as "infer this one".
        '<DataItem Dimensions="-1 2" Format="XML">0 0 1 0 0 1</DataItem>',
        # A product the text cannot hold: the buffer was sized from it first.
        '<DataItem Dimensions="1000000000 1000000000" Format="XML">0 0</DataItem>',
    ],
    ids=["negative-dimensions", "product-beyond-text"],
)
def test_dimensions_the_data_cannot_back_are_read_errors(tmp_path, engine, geometry):
    # Both engines refuse the same declarations; the C++ core declines to the
    # Python reader on a ReadError, so the Python reader must refuse them too.
    path = tmp_path / "bad.xdmf"
    path.write_text(_xdmf_document(geometry))
    with pytest.raises(meshioplusplus.ReadError):
        _read_with_engine(engine, path)


def _unbundle(bundle, out_dir):
    # A MIOB bundle: the magic and version, a u16 entry count, then per entry a
    # u16 name length, the name, a u32 data length and the data.
    data = bundle.read_bytes()
    count = int.from_bytes(data[5:7], "little")
    pos = 7
    for _ in range(count):
        n = int.from_bytes(data[pos : pos + 2], "little")
        name = data[pos + 2 : pos + 2 + n].decode()
        pos += 2 + n
        size = int.from_bytes(data[pos : pos + 4], "little")
        (out_dir / name).write_bytes(data[pos + 4 : pos + 4 + size])
        pos += 4 + size


@pytest.mark.parametrize("engine", ["core", "python"])
def test_an_unfiltered_chunk_stored_in_fewer_bytes_is_a_read_error(tmp_path, engine):
    # The fuzzer's crash: an unfiltered chunk stored in 36 bytes where its 13
    # int64 elements take 104. libhdf5 copied past the stored bytes, so both
    # engines refuse the dataset before reading it.
    bundle = (
        pathlib.Path(__file__).parent.parent
        / "fuzz"
        / "regressions"
        / "xdmf"
        / "unfiltered-chunk-stored-short.bundle"
    )
    _unbundle(bundle, tmp_path)
    with pytest.raises(meshioplusplus.ReadError):
        _read_with_engine(engine, tmp_path / "seed.xdmf")


@pytest.mark.parametrize("engine", ["core", "python"])
def test_binary_payload_shorter_than_its_dimensions_is_a_read_error(tmp_path, engine):
    payload = tmp_path / "pts.bin"
    np.array([0.0, 1.0], dtype=np.float32).tofile(payload)
    path = tmp_path / "short.xdmf"
    path.write_text(
        _xdmf_document(
            f'<DataItem Dimensions="3 2" Precision="4" Format="Binary">{payload}</DataItem>'
        )
    )
    with pytest.raises(meshioplusplus.ReadError):
        _read_with_engine(engine, path)


@pytest.mark.parametrize("n", [64, 125, 1000])
@pytest.mark.parametrize("mixed", [False, True])
def test_native_high_order_hexahedra(n, mixed, tmp_path):
    core = pytest.importorskip("meshioplusplus._core")
    points = np.random.default_rng(n).random((n, 3))
    cells = [(f"hexahedron{n}", np.arange(n)[None, :])]
    if mixed:
        cells = [("hexahedron", np.arange(8)[None, :])] + cells
    mesh = meshioplusplus.Mesh(points, cells)
    path = tmp_path / "high_order.xdmf"
    core.xdmf_write(str(path), mesh, "XML")
    back = core.xdmf_read(str(path))
    assert [c.type for c in back.cells] == [c.type for c in mesh.cells]
    np.testing.assert_array_equal(back.cells[-1].data, mesh.cells[-1].data)


@pytest.mark.parametrize("engine", ["core", "python"])
def test_field_data_round_trips_as_an_information_element(engine, tmp_path):
    """Roadmap 3.3.7: (tag, dim) pairs are written, other field arrays skipped."""
    from meshioplusplus import _core
    from meshioplusplus.xdmf import main as xdmf_main

    mesh = meshioplusplus.Mesh(
        np.eye(3),
        [("triangle", [[0, 1, 2]])],
        field_data={"wall": np.array([7, 2]), "gravity": np.array([9.81])},
    )
    path = tmp_path / f"{engine}.xdmf"
    if engine == "core":
        _core.xdmf_write(str(path), mesh, "XML")
    else:
        xdmf_main.write(path, mesh, data_format="XML")
    for back in (meshioplusplus.read(path), xdmf_main.read(path)):
        assert sorted(back.field_data) == ["wall"]
        np.testing.assert_array_equal(back.field_data["wall"], [7, 2])
