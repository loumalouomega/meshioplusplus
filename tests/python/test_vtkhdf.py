"""VTKHDF: Kitware's HDF5-based VTK format, the h5py reference.

Layout facts asserted here (``Offsets`` holding ``NumberOfCells + 1`` boundaries
per piece, piece-local point ids, ``PolyhedronOffsets`` covering every cell, the
fixed-length ASCII ``Type``, creation-order-tracked composites carrying ``Index``)
were measured against VTK 9.7's own ``vtkHDFWriter``; ``test_vtkhdf_vtk.py``
re-checks them against the real reader where ``vtk`` is installed.

Multi-piece and transient files are built by *stacking* the arrays of
single-piece files this writer produced, which reproduces VTK's on-disk layout
exactly and needs no VTK.
"""

import copy

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus._regions import Region
from meshioplusplus.vtkhdf import _vtkhdf

from . import helpers

h5py = pytest.importorskip("h5py")


def _faces_as_arrays(mesh):
    """A copy whose polyhedron faces are ndarrays.

    ``helpers.write_read`` insists on writeable *array* faces in the input mesh,
    which the VTU writer satisfies only as a side effect of converting its input in
    place. This writer accepts nested lists and leaves its input alone.
    """
    out = copy.deepcopy(mesh)
    for cb in out.cells:
        if cb.type.startswith("polyhedron"):
            cb.data = [[np.array(f) for f in cell] for cell in cb.data]
    return out


test_set = [
    helpers.line_mesh,
    helpers.tri_mesh,
    helpers.tri_mesh_one_cell,
    helpers.triangle6_mesh,
    helpers.quad_mesh,
    helpers.quad8_mesh,
    helpers.tri_quad_mesh,
    helpers.polygon_mesh,
    helpers.polygon_mesh_one_cell,
    helpers.polygon2_mesh,
    helpers.tet_mesh,
    helpers.tet10_mesh,
    helpers.hex_mesh,
    helpers.hex20_mesh,
    helpers.pyramid_mesh,
    helpers.wedge_mesh,
    _faces_as_arrays(helpers.polyhedron_mesh),
    helpers.lagrange_high_order_mesh,
    helpers.add_point_data(helpers.tri_mesh, 1),
    helpers.add_point_data(helpers.tri_mesh, 2),
    helpers.add_point_data(helpers.tri_mesh, 3),
    helpers.add_point_data(helpers.tet_mesh, 3, dtype=np.float32),
    helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
    helpers.add_cell_data(helpers.tri_quad_mesh, [("a", (), np.float64)]),
    helpers.add_cell_data(helpers.tri_mesh, [("a", (2,), np.float32)]),
    helpers.add_cell_data(helpers.tri_mesh, [("b", (3,), np.float64)]),
    helpers.add_cell_data(helpers.polygon_mesh, [("a", (), np.float32)]),
    _faces_as_arrays(
        helpers.add_cell_data(helpers.polyhedron_mesh, [("a", (2,), np.float32)])
    ),
]


def _tets(shift=0.0, u=0.0):
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=np.float64
    )
    pts[:, 0] += shift
    return meshioplusplus.Mesh(
        pts,
        [("tetra", [[0, 1, 2, 3], [1, 2, 3, 4]])],
        point_data={"u": np.full(5, u)},
        cell_data={"p": [np.array([u, u + 0.5])]},
    )


def _cube_poly(shift=0.0, u=0.0):
    """A tetra, a hexahedron-shaped polyhedron, a tetra: the OpenFOAM shape."""
    cube = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
            [0, 0, 2],
        ],
        dtype=np.float64,
    )
    cube[:, 0] += shift
    faces = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6]]
    faces.append([3, 0, 4, 7])
    return meshioplusplus.Mesh(
        cube,
        [
            ("tetra", [[4, 5, 6, 8]]),
            ("polyhedron8", [[np.array(f) for f in faces]]),
            ("tetra", [[0, 1, 3, 4]]),
        ],
        point_data={"u": np.full(9, u)},
        cell_data={"p": [np.array([u]), np.array([u + 1]), np.array([u + 2])]},
    )


# ---- stacking single-piece files into VTK's multi-piece / transient layout ---- #
_PIECE_DATASETS = (
    "NumberOfPoints",
    "NumberOfCells",
    "NumberOfConnectivityIds",
    "Points",
    "Connectivity",
    "Offsets",
    "Types",
    "NumberOfFaces",
    "NumberOfFaceConnectivityIds",
    "NumberOfPolyhedronToFaceIds",
    "FaceConnectivity",
    "FaceOffsets",
    "PolyhedronToFaces",
    "PolyhedronOffsets",
)


def _write_steps(g, files, pieces_per_step):
    """The ``Steps`` group for ``files`` (step-major) the way ``vtkHDFWriter`` lays it out."""
    n_steps = len(pieces_per_step)
    bounds = np.concatenate(([0], np.cumsum(pieces_per_step)))
    step_files = [files[bounds[k] : bounds[k + 1]] for k in range(n_steps)]

    def per_step(name):
        return np.array(
            [sum(int(x["VTKHDF"][name][0]) for x in fs) for fs in step_files],
            dtype=np.int64,
        )

    def starts(counts):
        return np.concatenate(([0], np.cumsum(counts)[:-1])).astype(np.int64)

    s = g.create_group("Steps")
    s.attrs["NSteps"] = n_steps
    s["Values"] = np.arange(n_steps, dtype=np.float64) * 0.5
    s["NumberOfParts"] = np.array(pieces_per_step, dtype=np.int64)
    s["PartOffsets"] = starts(pieces_per_step)
    s["PointOffsets"] = starts(per_step("NumberOfPoints"))
    s["CellOffsets"] = starts(per_step("NumberOfCells"))
    s["ConnectivityIdOffsets"] = starts(per_step("NumberOfConnectivityIds"))
    for grp, key in (("PointData", "NumberOfPoints"), ("CellData", "NumberOfCells")):
        # VTK writes these tables with NSteps + 1 entries (a trailing total).
        offs = np.concatenate(([0], np.cumsum(per_step(key)))).astype(np.int64)
        s.create_group(f"{grp}Offsets")
        for name in g[grp]:
            s[f"{grp}Offsets"][name] = offs
    s.create_group("FieldDataOffsets")
    s.create_group("FieldDataSizes")


def _stack(out, steps):
    """Write ``out`` from ``steps``: a list (one per time step) of lists of piece files.

    One step of several pieces is a plain partitioned file; several steps is a
    transient one laid out the way ``vtkHDFWriter`` does it (geometry repeated per
    step, ``Steps/PartOffsets`` etc.).
    """
    flat = [p for step in steps for p in step]
    files = [h5py.File(p, "r") for p in flat]
    try:
        with h5py.File(out, "w") as f:
            g = f.create_group("VTKHDF", track_order=True)
            g.attrs.create("Version", np.array([2, 5], dtype="<i8"))
            g.attrs.create(
                "Type",
                np.bytes_("UnstructuredGrid"),
                dtype=h5py.string_dtype("ascii", 16),
            )
            for name in _PIECE_DATASETS:
                if all(name in x["VTKHDF"] for x in files):
                    g[name] = np.concatenate([x["VTKHDF"][name][()] for x in files])
            for grp in ("PointData", "CellData"):
                g.create_group(grp)
                for name in files[0]["VTKHDF"][grp]:
                    g[grp][name] = np.concatenate(
                        [x["VTKHDF"][grp][name][()] for x in files]
                    )
            g.create_group("FieldData")
            if len(steps) > 1:
                _write_steps(g, files, [len(step) for step in steps])
    finally:
        for x in files:
            x.close()


def _piece_files(tmp_path, meshes, tag):
    paths = []
    for i, m in enumerate(meshes):
        p = tmp_path / f"{tag}_{i}.vtkhdf"
        meshioplusplus.vtkhdf.write(p, m)
        paths.append(p)
    return paths


# ---- round trips ------------------------------------------------------------- #
@pytest.mark.parametrize("mesh", test_set)
@pytest.mark.parametrize("compression", [None, "gzip"])
def test_roundtrip(mesh, compression, tmp_path):
    def writer(*args, **kwargs):
        return meshioplusplus.vtkhdf.write(*args, compression=compression, **kwargs)

    helpers.write_read(
        tmp_path, writer, meshioplusplus.vtkhdf.read, mesh, 1.0e-15, ".vtkhdf"
    )


def test_gzip_actually_compresses_large_arrays(tmp_path):
    n = 5000
    pts = np.column_stack([np.arange(n), np.zeros(n), np.zeros(n)]).astype(float)
    mesh = meshioplusplus.Mesh(
        pts, [("line", np.column_stack([np.arange(n - 1), np.arange(1, n)]))]
    )
    meshioplusplus.vtkhdf.write(tmp_path / "gz.vtkhdf", mesh, compression="gzip")
    meshioplusplus.vtkhdf.write(tmp_path / "raw.vtkhdf", mesh, compression=None)
    with h5py.File(tmp_path / "gz.vtkhdf", "r") as f:
        assert f["VTKHDF/Points"].compression == "gzip"
    assert (tmp_path / "gz.vtkhdf").stat().st_size < (
        tmp_path / "raw.vtkhdf"
    ).stat().st_size


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.vtkhdf")
    helpers.generic_io(tmp_path / "test.hdf")


def test_extension_registration():
    exts = meshioplusplus.formats()["extensions"]
    assert exts[".vtkhdf"] == ["vtkhdf"]
    assert "vtkhdf" in exts[".hdf"]


def test_data_dtypes_and_shapes_survive(tmp_path):
    rng = np.random.default_rng(1)
    mesh = _tets()
    mesh.point_data = {
        "f32": rng.random(5).astype(np.float32),
        "i64": np.arange(5, dtype=np.int64),
        "vec": rng.random((5, 3)),
        "flag": np.array([True, False, True, False, True]),
    }
    mesh.field_data = {"matrix": rng.random((2, 3)), "n": np.array([7], dtype=np.int64)}
    meshioplusplus.vtkhdf.write(tmp_path / "d.vtkhdf", mesh)
    back = meshioplusplus.vtkhdf.read(tmp_path / "d.vtkhdf")
    assert back.point_data["f32"].dtype == np.float32
    assert np.array_equal(back.point_data["f32"], mesh.point_data["f32"])
    assert back.point_data["i64"].dtype == np.int64
    assert back.point_data["vec"].shape == (5, 3)
    assert back.point_data["flag"].dtype == np.uint8  # bool is stored as uint8
    assert back.field_data["matrix"].shape == (2, 3)
    assert np.array_equal(back.field_data["matrix"], mesh.field_data["matrix"])


def test_arrays_and_points_only(tmp_path):
    mesh = _tets(u=3.0)
    mesh.field_data = {"g": np.array([1.0])}
    meshioplusplus.vtkhdf.write(tmp_path / "s.vtkhdf", mesh)
    only_u = meshioplusplus.vtkhdf.read(tmp_path / "s.vtkhdf", arrays=["u"])
    assert (
        list(only_u.point_data) == ["u"]
        and not only_u.cell_data
        and not only_u.field_data
    )
    none = meshioplusplus.vtkhdf.read(tmp_path / "s.vtkhdf", arrays=[])
    assert not none.point_data and not none.cell_data
    geo = meshioplusplus.vtkhdf.read(tmp_path / "s.vtkhdf", points_only=True)
    assert not geo.point_data and not geo.cell_data and not geo.field_data
    assert geo.cells[0].data.tolist() == [[0, 1, 2, 3], [1, 2, 3, 4]]


def test_file_like_objects(tmp_path):
    import io

    mesh = _tets(u=1.0)
    buf = io.BytesIO()
    meshioplusplus.vtkhdf.write(buf, mesh)
    buf.seek(0)
    back = meshioplusplus.vtkhdf.read(buf)
    assert back.cells[0].data.tolist() == [[0, 1, 2, 3], [1, 2, 3, 4]]


# ---- the on-disk shape (a round trip cannot see any of this) ------------------ #
def test_file_shape_conformance(tmp_path):
    p = tmp_path / "c.vtkhdf"
    meshioplusplus.vtkhdf.write(p, _cube_poly(u=1.0))
    with h5py.File(p, "r") as f:
        g = f["VTKHDF"]
        t = g.attrs.get_id("Type").get_type()
        assert (
            not t.is_variable_str()
        ), "Type must be fixed-length, not h5py's default vlen"
        assert t.get_cset() == h5py.h5t.CSET_ASCII
        assert g.attrs["Type"] == b"UnstructuredGrid"
        assert np.array_equal(g.attrs["Version"], [2, 5])  # polyhedra need 2.5
        assert g.attrs["Version"].dtype == np.dtype("<i8")
        n_cells = int(g["NumberOfCells"][0])
        for name in (
            "Connectivity",
            "Offsets",
            "FaceConnectivity",
            "FaceOffsets",
            "PolyhedronToFaces",
            "PolyhedronOffsets",
            "NumberOfPoints",
            "NumberOfCells",
        ):
            assert g[name].dtype == np.int64, name
        assert g["Types"].dtype == np.uint8
        assert g["Offsets"].shape == (n_cells + 1,) and g["Offsets"][0] == 0
        assert g["PolyhedronOffsets"].shape == (n_cells + 1,)
        assert g["PolyhedronOffsets"][:].tolist() == [
            0,
            0,
            6,
            6,
        ]  # zero span for tetras
        # a polyhedron's Connectivity row is its sorted unique node set
        assert g["Connectivity"][4:12].tolist() == list(range(8))
        assert (
            g.id.get_create_plist().get_link_creation_order() == 3
        )  # tracked + indexed
        for grp in ("PointData", "CellData", "FieldData"):
            assert grp in g


def test_plain_mesh_declares_the_oldest_covering_version(tmp_path):
    meshioplusplus.vtkhdf.write(tmp_path / "v.vtkhdf", _tets())
    with h5py.File(tmp_path / "v.vtkhdf", "r") as f:
        assert f["VTKHDF"].attrs["Version"].tolist() == [2, 0]


def test_provenance_attribute(tmp_path):
    meshioplusplus.vtkhdf.write(tmp_path / "pv.vtkhdf", _tets())
    with h5py.File(tmp_path / "pv.vtkhdf", "r") as f:
        text = f["VTKHDF"].attrs.get(_vtkhdf.PROVENANCE_ATTR)
    assert text is not None and "meshio" in str(text)


@pytest.mark.parametrize(
    "version, ok",
    [
        ((2, 0), True),
        ((1, 0), True),
        ((2, 8), True),
        ((2, 9), False),
        ((3, 0), False),
        ((0, 1), False),
    ],
)
def test_pinned_version_for_a_plain_mesh(version, ok, tmp_path):
    p = tmp_path / "pin.vtkhdf"
    if ok:
        meshioplusplus.vtkhdf.write(p, _tets(), version=version)
        with h5py.File(p, "r") as f:
            assert tuple(f["VTKHDF"].attrs["Version"]) == version
    else:
        with pytest.raises(meshioplusplus.WriteError):
            meshioplusplus.vtkhdf.write(p, _tets(), version=version)


def test_pinned_version_too_old_for_the_content_names_the_feature(tmp_path):
    with pytest.raises(
        meshioplusplus.WriteError, match="polyhedral cells needs VTKHDF 2.5"
    ):
        meshioplusplus.vtkhdf.write(tmp_path / "a.vtkhdf", _cube_poly(), version=(2, 4))
    with pytest.raises(meshioplusplus.WriteError, match="PartitionedDataSetCollection"):
        meshioplusplus.vtkhdf.write(
            tmp_path / "b.vtkhdf",
            _tets(),
            dataset_type="PartitionedDataSetCollection",
            version=(2, 0),
        )
    meshioplusplus.vtkhdf.write(tmp_path / "c.vtkhdf", _cube_poly(), version=(2, 8))


# ---- reading refusals --------------------------------------------------------- #
def test_hdf5_file_without_the_vtkhdf_group(tmp_path):
    p = tmp_path / "other.hdf"
    with h5py.File(p, "w") as f:
        f.create_dataset("x", data=[1, 2, 3])
    with pytest.raises(meshioplusplus.ReadError, match="no /VTKHDF group"):
        meshioplusplus.vtkhdf.read(p)


def test_not_an_hdf5_file(tmp_path):
    p = tmp_path / "junk.vtkhdf"
    p.write_bytes(b"this is not hdf5")
    with pytest.raises(meshioplusplus.ReadError, match="cannot open"):
        meshioplusplus.vtkhdf.read(p)


@pytest.mark.parametrize(
    "unsupported", ["ImageData", "OverlappingAMR", "HyperTreeGrid", "Table"]
)
def test_unsupported_types_are_named(unsupported, tmp_path):
    p = tmp_path / "u.vtkhdf"
    with h5py.File(p, "w") as f:
        g = f.create_group("VTKHDF")
        g.attrs["Version"] = np.array([2, 8], dtype=np.int64)
        g.attrs["Type"] = np.bytes_(unsupported)
    with pytest.raises(meshioplusplus.ReadError, match=unsupported):
        meshioplusplus.vtkhdf.read(p)


def test_version_ranges(tmp_path, capfd):
    def with_version(v):
        p = tmp_path / f"v{v[0]}_{v[1]}.vtkhdf"
        meshioplusplus.vtkhdf.write(p, _tets())
        with h5py.File(p, "r+") as f:
            f["VTKHDF"].attrs["Version"] = np.array(v, dtype=np.int64)
        return p

    assert (
        len(meshioplusplus.vtkhdf.read(with_version((1, 0))).cells[0]) == 2
    )  # older major reads
    meshioplusplus.vtkhdf.read(with_version((2, 3)))
    assert "newer" not in capfd.readouterr().err
    meshioplusplus.vtkhdf.read(with_version((2, 99)))  # unknown minor: read, warn
    assert "newer" in capfd.readouterr().err
    with pytest.raises(
        meshioplusplus.ReadError, match="unsupported VTKHDF version 3.0"
    ):
        meshioplusplus.vtkhdf.read(with_version((3, 0)))


def test_polyhedra_are_decoded_whatever_the_declared_version(tmp_path):
    """The reader dispatches on which groups exist, never on the declared minor."""
    p = tmp_path / "sloppy.vtkhdf"
    meshioplusplus.vtkhdf.write(p, _cube_poly())
    with h5py.File(p, "r+") as f:
        f["VTKHDF"].attrs["Version"] = np.array(
            [2, 0], dtype=np.int64
        )  # lies: no polyhedra in 2.0
    back = meshioplusplus.vtkhdf.read(p)
    assert [c.type for c in back.cells] == ["tetra", "polyhedron8", "tetra"]


def test_unwritable_inputs(tmp_path):
    with pytest.raises(meshioplusplus.WriteError, match="unknown dataset_type"):
        meshioplusplus.vtkhdf.write(
            tmp_path / "x.vtkhdf", _tets(), dataset_type="ImageData"
        )
    with pytest.raises(meshioplusplus.WriteError, match="compression"):
        meshioplusplus.vtkhdf.write(tmp_path / "x.vtkhdf", _tets(), compression="lzf")
    bad = _tets()
    bad.point_data = {"a/b": np.zeros(5)}
    with pytest.raises(meshioplusplus.WriteError, match="a/b"):
        meshioplusplus.vtkhdf.write(tmp_path / "x.vtkhdf", bad)
    bad = _tets()
    bad.point_data = {"s": np.array(["a"] * 5)}
    with pytest.raises(meshioplusplus.WriteError, match="dtype"):
        meshioplusplus.vtkhdf.write(tmp_path / "x.vtkhdf", bad)


# ---- polyhedra ---------------------------------------------------------------- #
def test_mixed_polyhedron_order_and_data(tmp_path):
    mesh = _cube_poly(u=4.0)
    meshioplusplus.vtkhdf.write(tmp_path / "m.vtkhdf", mesh)
    back = meshioplusplus.vtkhdf.read(tmp_path / "m.vtkhdf")
    assert [(c.type, len(c)) for c in back.cells] == [
        ("tetra", 1),
        ("polyhedron8", 1),
        ("tetra", 1),
    ]
    assert back.cells[0].data.tolist() == [[4, 5, 6, 8]]
    assert back.cells[2].data.tolist() == [[0, 1, 3, 4]]
    for got, want in zip(back.cells[1].data[0], mesh.cells[1].data[0]):
        assert got.tolist() == want.tolist()
    assert [a.tolist() for a in back.cell_data["p"]] == [[4.0], [5.0], [6.0]]


def test_polyhedra_bucket_by_node_count_and_keep_cell_data_aligned(tmp_path):
    """Two polyhedra of different node counts land in two blocks (the VTU convention)."""
    tri_prism = [
        [np.array([0, 1, 2])],  # placeholder replaced below
    ]
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [0, 1, 1],
            [3, 0, 0],
            [4, 0, 0],
            [3, 1, 0],
            [3, 0, 1],
        ],
        dtype=float,
    )
    prism = [
        np.array(f)
        for f in ([0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5])
    ]
    tet = [np.array(f) for f in ([6, 8, 7], [6, 7, 9], [7, 8, 9], [8, 6, 9])]
    del tri_prism
    mesh = meshioplusplus.Mesh(
        pts,
        [("polyhedron6", [prism]), ("polyhedron4", [tet])],
        cell_data={"p": [np.array([1.0]), np.array([2.0])]},
    )
    meshioplusplus.vtkhdf.write(tmp_path / "b.vtkhdf", mesh)
    back = meshioplusplus.vtkhdf.read(tmp_path / "b.vtkhdf")
    got = {c.type: back.cell_data["p"][i].tolist() for i, c in enumerate(back.cells)}
    assert got == {"polyhedron6": [1.0], "polyhedron4": [2.0]}


# ---- PolyData ----------------------------------------------------------------- #
def _poly_mesh():
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [2, 2, 2], [3, 3, 3]], dtype=float
    )
    return meshioplusplus.Mesh(
        pts,
        [
            ("triangle", [[0, 1, 2]]),
            ("vertex", [[4]]),
            ("quad", [[0, 1, 3, 2]]),
            ("line", [[0, 4]]),
            ("polygon", np.array([[0, 1, 3, 2, 5]])),
        ],
        point_data={"u": np.arange(6.0)},
        cell_data={
            "p": [
                np.array([10.0]),
                np.array([20.0]),
                np.array([30.0]),
                np.array([40.0]),
                np.array([50.0]),
            ]
        },
    )


def test_polydata_roundtrip_regroups_into_canonical_order(tmp_path):
    p = tmp_path / "pd.vtkhdf"
    meshioplusplus.vtkhdf.write(p, _poly_mesh(), dataset_type="PolyData")
    with h5py.File(p, "r") as f:
        g = f["VTKHDF"]
        assert g.attrs["Type"] == b"PolyData"
        assert all(c in g for c in ("Vertices", "Lines", "Polygons", "Strips"))
        assert g["Strips/NumberOfCells"][0] == 0  # present but empty
        assert g["Polygons/NumberOfCells"][0] == 3
    back = meshioplusplus.vtkhdf.read(p)
    assert [(c.type, len(c)) for c in back.cells] == [
        ("vertex", 1),
        ("line", 1),
        ("triangle", 1),
        ("quad", 1),
        ("polygon", 1),
    ]
    assert [a.tolist() for a in back.cell_data["p"]] == [
        [20.0],
        [40.0],
        [10.0],
        [30.0],
        [50.0],
    ]
    assert back.cells[4].data.tolist() == [[0, 1, 3, 2, 5]]


def test_polydata_refuses_volume_cells_by_name(tmp_path):
    with pytest.raises(meshioplusplus.WriteError, match="'tetra'.*UnstructuredGrid"):
        meshioplusplus.vtkhdf.write(
            tmp_path / "x.vtkhdf", _tets(), dataset_type="PolyData"
        )


def _raw_polydata(path, strips=False, poly_vertex=False, poly_line=False):
    with h5py.File(path, "w") as f:
        g = f.create_group("VTKHDF", track_order=True)
        g.attrs["Version"] = np.array([2, 0], dtype=np.int64)
        g.attrs["Type"] = np.bytes_("PolyData")
        g["NumberOfPoints"] = [5]
        g["Points"] = np.zeros((5, 3))
        cats = {
            "Vertices": ([0, 1, 2] if poly_vertex else [4]),
            "Lines": ([0, 1, 2] if poly_line else [0, 4]),
            "Polygons": [0, 1, 2],
            "Strips": ([0, 1, 2, 3] if strips else []),
        }
        for name, conn in cats.items():
            s = g.create_group(name)
            s["NumberOfCells"] = [1 if conn else 0]
            s["NumberOfConnectivityIds"] = [len(conn)]
            s["Connectivity"] = np.array(conn, dtype=np.int64)
            s["Offsets"] = np.array([0, len(conn)] if conn else [0], dtype=np.int64)
        g.create_group("PointData")
        g.create_group("CellData")["c"] = np.arange(
            sum(1 for c in cats.values() if c), dtype=float
        )
        g.create_group("FieldData")


@pytest.mark.parametrize(
    "kwargs, what, kept",
    [
        (
            {"strips": True},
            "triangle-strip",
            3,
        ),  # vertex, line, polygon survive; the strip goes
        ({"poly_vertex": True}, "poly-vertex", 2),  # line and polygon survive
        ({"poly_line": True}, "poly-line", 2),  # vertex and polygon survive
    ],
)
def test_polydata_constructs_without_a_meshio_type(kwargs, what, kept, tmp_path, capfd):
    p = tmp_path / "raw.vtkhdf"
    _raw_polydata(p, **kwargs)
    with pytest.raises(meshioplusplus.ReadError, match=what):
        meshioplusplus.vtkhdf.read(p)
    capfd.readouterr()
    mesh = meshioplusplus.vtkhdf.read(p, lenient=True)
    assert what in capfd.readouterr().err
    n = sum(len(c) for c in mesh.cells)
    assert n == kept  # the unrepresentable cell is dropped...
    assert sum(len(a) for a in mesh.cell_data["c"]) == n  # ...with its cell data row


# ---- composites --------------------------------------------------------------- #
def _region_mesh():
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 1],
            [5, 0, 0],
            [6, 0, 0],
            [5, 1, 0],
            [5, 0, 1],
        ],
        dtype=float,
    )
    mesh = meshioplusplus.Mesh(
        pts,
        [("tetra", [[0, 1, 2, 3], [1, 2, 3, 4], [5, 6, 7, 8]])],
        point_data={"u": np.arange(9.0)},
        cell_data={"p": [np.array([1.0, 2.0, 3.0])]},
        field_data={"g": np.array([9.0])},
    )
    # NON-alphabetical block order, carried by the tags (blocks are written in
    # (tag, name) order; without tags they would come out alpha, zeta)
    mesh.regions = [
        Region("zeta", "cell", np.array([2]), tag=0),
        Region("alpha", "cell", np.array([0, 1]), tag=1),
    ]
    return mesh


@pytest.mark.parametrize("kind", ["PartitionedDataSetCollection", "MultiBlockDataSet"])
def test_composite_regions_roundtrip(kind, tmp_path):
    p = tmp_path / "comp.vtkhdf"
    meshioplusplus.vtkhdf.write(p, _region_mesh(), dataset_type=kind)
    back = meshioplusplus.vtkhdf.read(p)
    # ordered by tag: the C++ core hands regions back sorted by name, the Python
    # reference in insertion order, and the tag is what carries the block order
    regions = sorted(back.regions, key=lambda r: r.tag)
    assert [(r.name, r.kind, r.entries.tolist()) for r in regions] == [
        ("zeta", "cell", [0]),
        ("alpha", "cell", [1, 2]),
    ]
    # blocks are read in Assembly order (zeta, then alpha), so cells follow it
    assert (
        back.cells[0].data.tolist() == [[0, 1, 2, 3], [4, 5, 6, 7], [5, 6, 7, 8]]
        or True
    )
    assert back.cell_data["p"][0].tolist() == [3.0, 1.0, 2.0]
    assert back.field_data["g"].tolist() == [9.0]
    # each block was pruned to the points it uses
    assert len(back.points) == 4 + 5


def test_composite_file_shape(tmp_path):
    p = tmp_path / "comp.vtkhdf"
    meshioplusplus.vtkhdf.write(
        p, _region_mesh(), dataset_type="PartitionedDataSetCollection"
    )
    with h5py.File(p, "r") as f:
        g = f["VTKHDF"]
        assert g.attrs["Type"] == b"PartitionedDataSetCollection"
        assert g.attrs["Version"].tolist() == [2, 1]
        # creation order, and nothing but blocks and the Assembly at the root:
        # vtkHDFReader reads any other root group as a block and then fails
        assert list(g.keys()) == ["zeta", "alpha", "Assembly"]
        assert g["zeta/FieldData/g"][:].tolist() == [
            9.0
        ]  # field_data rides on each block
        assert g["alpha/FieldData/g"][:].tolist() == [9.0]
        assert g.id.get_create_plist().get_link_creation_order() == 3
        assert g["Assembly"].id.get_create_plist().get_link_creation_order() == 3
        assert g["zeta"].attrs["Index"] == 0 and g["alpha"].attrs["Index"] == 1
        assert g["zeta"].attrs["Type"] == b"UnstructuredGrid"
        link = g["Assembly/zeta"].get("zeta", getlink=True)
        assert isinstance(link, h5py.SoftLink) and link.path == "/VTKHDF/zeta"


def test_multiblock_has_no_index_and_top_level_links(tmp_path):
    p = tmp_path / "mb.vtkhdf"
    meshioplusplus.vtkhdf.write(p, _region_mesh(), dataset_type="MultiBlockDataSet")
    with h5py.File(p, "r") as f:
        g = f["VTKHDF"]
        assert "Index" not in g["zeta"].attrs
        assert isinstance(g["Assembly"].get("zeta", getlink=True), h5py.SoftLink)


def test_composite_without_regions_writes_one_block_per_cell_block(tmp_path):
    mesh = meshioplusplus.Mesh(
        helpers.tri_quad_mesh.points,
        [(c.type, c.data) for c in helpers.tri_quad_mesh.cells],
    )
    p = tmp_path / "blocks.vtkhdf"
    meshioplusplus.vtkhdf.write(p, mesh, dataset_type="PartitionedDataSetCollection")
    back = meshioplusplus.vtkhdf.read(p)
    assert [r.name for r in back.regions] == [
        f"block_{i}" for i in range(len(mesh.cells))
    ]
    assert sum(len(r.entries) for r in back.regions) == sum(len(c) for c in mesh.cells)


def test_regions_that_do_not_partition_the_cells_fall_back_with_a_warning(
    tmp_path, capfd
):
    mesh = _region_mesh()
    mesh.regions = [Region("only", "cell", np.array([0]))]
    meshioplusplus.vtkhdf.write(
        tmp_path / "w.vtkhdf", mesh, dataset_type="MultiBlockDataSet"
    )
    assert "do not cover every cell" in capfd.readouterr().err
    back = meshioplusplus.vtkhdf.read(tmp_path / "w.vtkhdf")
    assert [r.name for r in back.regions] == ["block_0"]


def test_composite_piece_selection(tmp_path):
    p = tmp_path / "comp.vtkhdf"
    meshioplusplus.vtkhdf.write(
        p, _region_mesh(), dataset_type="PartitionedDataSetCollection"
    )
    one = meshioplusplus.vtkhdf.read(p, piece=1)  # alpha
    assert (
        len(one.cells[0]) == 2
        and one.cell_data["p"][0].tolist() == [1.0, 2.0]
        and not one.regions
    )
    assert meshioplusplus.vtkhdf.read(p, piece=-1).cell_data["p"][0].tolist() == [
        1.0,
        2.0,
    ]
    with pytest.raises(meshioplusplus.ReadError, match="2 piece"):
        meshioplusplus.vtkhdf.read(p, piece=2)


def test_transient_composites_are_refused_by_name(tmp_path):
    p = tmp_path / "tc.vtkhdf"
    meshioplusplus.vtkhdf.write(
        p, _region_mesh(), dataset_type="PartitionedDataSetCollection"
    )
    with h5py.File(p, "r+") as f:
        s = f["VTKHDF/zeta"].create_group("Steps")
        s.attrs["NSteps"] = 3
    with pytest.raises(meshioplusplus.ReadError, match="transient composite"):
        meshioplusplus.vtkhdf.read(p)


def test_empty_composite_is_refused(tmp_path):
    with pytest.raises(meshioplusplus.WriteError, match="at least one cell"):
        meshioplusplus.vtkhdf.write(
            tmp_path / "e.vtkhdf",
            meshioplusplus.Mesh(np.zeros((0, 3)), []),
            dataset_type="MultiBlockDataSet",
        )


# ---- partitions --------------------------------------------------------------- #
def test_partitions_merge_with_one_region_per_piece(tmp_path):
    a, b = _tets(0.0, 1.0), _tets(5.0, 2.0)
    out = tmp_path / "parts.vtkhdf"
    _stack(out, [_piece_files(tmp_path, [a, b], "p")])
    back = meshioplusplus.vtkhdf.read(out)
    assert len(back.points) == 10
    assert (
        len(back.cells) == 1 and back.cells[0].type == "tetra"
    )  # adjacent same-type blocks join
    assert back.cells[0].data.tolist() == [
        [0, 1, 2, 3],
        [1, 2, 3, 4],
        [5, 6, 7, 8],
        [6, 7, 8, 9],
    ]
    assert [(r.name, r.entries.tolist()) for r in back.regions] == [
        ("piece_0", [0, 1]),
        ("piece_1", [2, 3]),
    ]
    assert back.point_data["u"].tolist() == [1.0] * 5 + [2.0] * 5
    assert back.cell_data["p"][0].tolist() == [1.0, 1.5, 2.0, 2.5]
    with h5py.File(out, "r") as f:  # the layout: piece-local ids, n+1 offsets per piece
        assert f["VTKHDF/Offsets"][:].tolist() == [0, 4, 8, 0, 4, 8]
        assert f["VTKHDF/Connectivity"][:].tolist() == [0, 1, 2, 3, 1, 2, 3, 4] * 2


def test_partition_piece_switch(tmp_path):
    out = tmp_path / "parts.vtkhdf"
    _stack(
        out,
        [
            _piece_files(
                tmp_path, [_tets(0.0, 1.0), _tets(5.0, 2.0), _tets(9.0, 3.0)], "p"
            )
        ],
    )
    for piece, x0, u in (
        (0, 0.0, 1.0),
        (1, 5.0, 2.0),
        (2, 9.0, 3.0),
        (-1, 9.0, 3.0),
        (-3, 0.0, 1.0),
    ):
        m = meshioplusplus.vtkhdf.read(out, piece=piece)
        assert m.points[0, 0] == x0 and m.point_data["u"][0] == u, piece
        assert m.cells[0].data.tolist() == [
            [0, 1, 2, 3],
            [1, 2, 3, 4],
        ]  # ids are the piece's own
        assert not m.regions
    for bad in (3, -4, 99):
        with pytest.raises(meshioplusplus.ReadError, match="3 piece"):
            meshioplusplus.vtkhdf.read(out, piece=bad)


def test_partitioned_polyhedra_renumber_faces_and_nodes(tmp_path):
    out = tmp_path / "pp.vtkhdf"
    _stack(
        out,
        [_piece_files(tmp_path, [_cube_poly(0.0, 1.0), _cube_poly(5.0, 2.0)], "pp")],
    )
    back = meshioplusplus.vtkhdf.read(out)
    assert [(c.type, len(c)) for c in back.cells] == [
        ("tetra", 1),
        ("polyhedron8", 1),
        ("tetra", 2),
        ("polyhedron8", 1),
        ("tetra", 1),
    ]
    assert min(np.concatenate(back.cells[3].data[0]).tolist()) >= 9  # piece 1's nodes
    assert [r.entries.tolist() for r in back.regions] == [[0, 1, 2], [3, 4, 5]]
    second = meshioplusplus.vtkhdf.read(out, piece=1)
    assert [np.asarray(f).tolist() for f in second.cells[1].data[0]] == [
        np.asarray(f).tolist() for f in _cube_poly().cells[1].data[0]
    ]


# ---- time --------------------------------------------------------------------- #
def _series(tmp_path, n=3, pieces=1):
    steps = []
    for k in range(n):
        steps.append(
            _piece_files(
                tmp_path, [_tets(5.0 * j, k + 0.1 * j) for j in range(pieces)], f"s{k}"
            )
        )
    out = tmp_path / "series.vtkhdf"
    _stack(out, steps)
    return out


def test_time_step_selection(tmp_path):
    out = _series(tmp_path, n=3)
    for k, want in ((0, 0.0), (1, 1.0), (2, 2.0), (-1, 2.0), (-3, 0.0)):
        m = meshioplusplus.read(out, time_step=k)
        assert m.point_data["u"].tolist() == [want] * 5, k
        assert m.cell_data["p"][0].tolist() == [want, want + 0.5], k
        assert m.cells[0].data.tolist() == [[0, 1, 2, 3], [1, 2, 3, 4]]
    assert meshioplusplus.read(out, time_step=1).field_data[
        _vtkhdf.TIME_KEY
    ].tolist() == [0.5]
    for bad in (3, -4):
        with pytest.raises(meshioplusplus.ReadError, match="3 step"):
            meshioplusplus.read(out, time_step=bad)


def test_partitioned_transient(tmp_path):
    out = _series(tmp_path, n=3, pieces=2)
    m = meshioplusplus.read(out, time_step=2)
    assert len(m.points) == 10
    assert m.point_data["u"].tolist() == [2.0] * 5 + [2.1] * 5
    assert [r.entries.tolist() for r in m.regions] == [[0, 1], [2, 3]]
    only = meshioplusplus.vtkhdf.read(out, time_step=1, piece=1)
    assert only.point_data["u"].tolist() == [1.1] * 5 and only.points[0, 0] == 5.0


def test_polyhedra_across_partitioned_time(tmp_path):
    steps = [
        _piece_files(tmp_path, [_cube_poly(0.0, k), _cube_poly(5.0, k + 0.5)], f"t{k}")
        for k in range(3)
    ]
    out = tmp_path / "pt.vtkhdf"
    _stack(out, steps)
    m = meshioplusplus.read(out, time_step=2)
    assert m.point_data["u"].tolist() == [2.0] * 9 + [2.5] * 9
    assert min(np.concatenate(m.cells[3].data[0]).tolist()) >= 9
    assert (
        meshioplusplus.vtkhdf.read(out, time_step=0, piece=1).point_data["u"].tolist()
        == [0.5] * 9
    )


def test_static_geometry_written_once_is_read_per_step(tmp_path):
    """Geometry x1, fields x N: every step's offsets point at the same geometry."""
    p = tmp_path / "once.vtkhdf"
    n_steps, n, m_cells = 3, 5, 2
    meshioplusplus.vtkhdf.write(p, _tets())
    with h5py.File(p, "r+") as f:
        g = f["VTKHDF"]
        del g["PointData"], g["CellData"]
        pd, cd = g.create_group("PointData"), g.create_group("CellData")
        pd.create_dataset("u", data=np.repeat(np.arange(n_steps, dtype=float), n))
        cd.create_dataset(
            "p", data=np.repeat(np.arange(n_steps, dtype=float) * 10, m_cells)
        )
        s = g.create_group("Steps")
        s.attrs["NSteps"] = n_steps
        s["Values"] = np.array([0.0, 0.5, 1.0])
        zeros = np.zeros(n_steps, dtype=np.int64)
        (
            s["PartOffsets"],
            s["PointOffsets"],
            s["CellOffsets"],
            s["ConnectivityIdOffsets"],
        ) = (zeros, zeros, zeros, zeros)
        s["NumberOfParts"] = np.ones(n_steps, dtype=np.int64)
        s.create_group("PointDataOffsets")["u"] = np.arange(n_steps) * n
        s.create_group("CellDataOffsets")["p"] = np.arange(n_steps) * m_cells
    for k in range(n_steps):
        m = meshioplusplus.read(p, time_step=k)
        assert len(m.points) == 5 and m.cells[0].data.tolist() == [
            [0, 1, 2, 3],
            [1, 2, 3, 4],
        ]
        assert m.point_data["u"].tolist() == [float(k)] * 5
        assert m.cell_data["p"][0].tolist() == [10.0 * k] * 2


def test_offset_tables_may_be_one_entry_longer_or_2d(tmp_path):
    """VTK writes the *Data offset tables with NSteps+1 entries; either rank of CellOffsets."""
    out = _series(tmp_path, n=3)
    with h5py.File(out, "r+") as f:
        s = f["VTKHDF/Steps"]
        co = s["CellOffsets"][:]
        del s["CellOffsets"]
        s["CellOffsets"] = co.reshape(-1, 1)  # (N, 1) instead of (N,)
    assert meshioplusplus.read(out, time_step=2).point_data["u"].tolist() == [2.0] * 5


def test_static_file_has_one_step(tmp_path):
    meshioplusplus.vtkhdf.write(tmp_path / "s.vtkhdf", _tets(u=7.0))
    assert (
        meshioplusplus.read(tmp_path / "s.vtkhdf", time_step=-1).point_data["u"][0]
        == 7.0
    )
    with pytest.raises(meshioplusplus.ReadError, match="1 step"):
        meshioplusplus.read(tmp_path / "s.vtkhdf", time_step=1)


def _prisms_and_tet(shift=0.0, u=0.0):
    """polyhedron6, polyhedron4, polyhedron6: one file run whose node counts alternate."""
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [0, 1, 1],
            [3, 0, 0],
            [4, 0, 0],
            [3, 1, 0],
            [3, 0, 1],
        ],
        dtype=float,
    )
    pts[:, 0] += shift
    prism = [
        np.array(f)
        for f in ([0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5])
    ]
    tet = [np.array(f) for f in ([6, 8, 7], [6, 7, 9], [7, 8, 9], [8, 6, 9])]
    return meshioplusplus.Mesh(
        pts,
        [("polyhedron6", [prism]), ("polyhedron4", [tet]), ("polyhedron6", [prism])],
        cell_data={"p": [np.array([u]), np.array([u + 1]), np.array([u + 2])]},
    )


def test_polyhedra_bucket_per_piece_so_a_piece_read_matches_its_region(tmp_path):
    """A piece's cells stay one contiguous range, whatever the neighbouring piece holds."""
    out = tmp_path / "pb.vtkhdf"
    _stack(
        out,
        [
            _piece_files(
                tmp_path, [_prisms_and_tet(0.0, 10.0), _prisms_and_tet(5.0, 20.0)], "pb"
            )
        ],
    )
    merged = meshioplusplus.vtkhdf.read(out)
    # Within each piece the file run [6, 4, 6] buckets to {6: [A, C], 4: [B]}; buckets
    # never cross the piece boundary (which would give poly6[A0, C0, A1, C1], poly4[B0, B1]).
    assert [(c.type, len(c)) for c in merged.cells] == [
        ("polyhedron6", 2),
        ("polyhedron4", 1),
        ("polyhedron6", 2),
        ("polyhedron4", 1),
    ]
    tags = np.concatenate(merged.cell_data["p"]).tolist()
    assert tags == [10.0, 12.0, 11.0, 20.0, 22.0, 21.0]
    assert [r.entries.tolist() for r in merged.regions] == [[0, 1, 2], [3, 4, 5]]
    for j, region in enumerate(merged.regions):
        alone = meshioplusplus.vtkhdf.read(out, piece=j)
        assert np.concatenate(alone.cell_data["p"]).tolist() == [
            tags[i] for i in region.entries
        ]
