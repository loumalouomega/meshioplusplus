"""VTK XML parallel indices `.pvtu` / `.pvtp` (v15.0.0).

An index that declares the arrays and names one `.vtu`/`.vtp` piece per part.
Every behavioural test runs against both engines: the C++ core (with strict-core
on, so a silent fallback to the Python twin would fail loudly) and the
pure-Python reference.
"""

import os

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import ReadError, WriteError, _core
from meshioplusplus._fallback import set_strict_core
from meshioplusplus.pvtp import _pvtp
from meshioplusplus.pvtu import _pvtu

from . import helpers

# The core-engine cases run the C++ path under strict-core with the default
# zlib codec, which a build without zlib (the Windows wheels) cannot honour --
# there the Python twin covers the format and these cases are skipped.
HAS_ZLIB = getattr(_core, "__has_zlib__", False)
requires_zlib = pytest.mark.skipif(not HAS_ZLIB, reason="build has no zlib")

INDEX_TAG = {"pvtu": "PUnstructuredGrid", "pvtp": "PPolyData"}


def _grid(n=6):
    """An n x n square split into triangles, with a point and a cell array."""
    xs, ys = np.meshgrid(np.arange(n + 1.0), np.arange(n + 1.0))
    points = np.c_[xs.ravel(), ys.ravel(), np.zeros((n + 1) ** 2)]
    tris = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            b, c = a + 1, a + n + 1
            tris += [[a, b, c + 1], [a, c + 1, c]]
    mesh = meshioplusplus.Mesh(points, [("triangle", np.array(tris))])
    mesh.point_data["u"] = points[:, 0]
    mesh.cell_data["c"] = [np.arange(len(tris), dtype=np.float64)]
    return mesh


def _num_cells(mesh):
    return sum(len(cb.data) for cb in mesh.cells)


def _assert_same_up_to_point_order(a, b):
    """Same point set and same cells, whatever the point numbering."""

    def canon(mesh):
        pts = np.asarray(mesh.points)
        order = np.lexsort(pts.T[::-1])
        rank = np.empty_like(order)
        rank[order] = np.arange(len(order))
        cells = {}
        for cb in mesh.cells:
            rows = sorted(map(tuple, np.sort(rank[np.asarray(cb.data)], axis=1)))
            cells.setdefault(cb.type, []).extend(rows)
        return pts[order], {k: sorted(v) for k, v in cells.items()}

    pa, ca = canon(a)
    pb, cb = canon(b)
    assert pa.shape == pb.shape and np.allclose(pa, pb)
    assert ca == cb


class _Engine:
    """One engine's `read`/`write`/`write_pieces` for one index kind."""

    def __init__(self, name, kind, mod):
        self.name, self.kind, self.mod = name, kind, mod
        self.read = mod.read
        self.write = mod.write
        self.write_pieces = mod.write_pieces


def _engines(kind):
    pkg = getattr(meshioplusplus, kind)
    twin = {"pvtu": _pvtu, "pvtp": _pvtp}[kind]
    return {"core": _Engine("core", kind, pkg), "python": _Engine("python", kind, twin)}


@pytest.fixture(params=["core", "python"])
def engine(request):
    """The pvtu engine under test; the core one runs with strict-core on."""
    if request.param == "core":
        if not HAS_ZLIB:
            pytest.skip("the core engine writes zlib and this build has none")
        if not hasattr(_core, "pvtu_read"):
            pytest.skip("this build has no pvtu core")
        set_strict_core(True)
        yield _engines("pvtu")["core"]
        set_strict_core(None)
    else:
        yield _engines("pvtu")["python"]


@pytest.fixture
def strict_core():
    set_strict_core(True)
    yield
    set_strict_core(None)


def _labelled(n=6, nparts=3):
    mesh = _grid(n)
    mesh.cell_data["partition:part"] = meshioplusplus.partition_labels(mesh, nparts)
    return mesh


@pytest.mark.parametrize(
    "binary, compression", [(False, None), (True, None), (True, "zlib")]
)
def test_roundtrip(engine, binary, compression, tmp_path):
    mesh = _labelled()
    path = tmp_path / "g.pvtu"
    engine.write(path, mesh, binary=binary, compression=compression)
    back = engine.read(path)

    assert _num_cells(back) == _num_cells(mesh)
    assert len(back.points) > len(mesh.points)  # interface points once per piece
    assert [r.name for r in back.regions] == ["piece_0", "piece_1", "piece_2"]
    assert all(r.kind == "cell" for r in back.regions)
    assert sum(len(r.entries) for r in back.regions) == _num_cells(mesh)
    assert sorted(back.point_data) == ["u"]
    assert sorted(back.cell_data) == ["c", "partition:part"]
    np.testing.assert_array_equal(
        np.sort(np.concatenate(back.cell_data["c"])),
        np.sort(np.concatenate(mesh.cell_data["c"])),
    )


def test_index_layout(engine, tmp_path):
    mesh = _labelled()
    engine.write(tmp_path / "g.pvtu", mesh)
    text = (tmp_path / "g.pvtu").read_text()
    assert '<VTKFile type="PUnstructuredGrid"' in text
    assert '<PUnstructuredGrid GhostLevel="0">' in text
    assert '<PDataArray type="Float64" Name="u"/>' in text
    assert '<PDataArray type="Float64" Name="c"/>' in text
    assert '<PDataArray type="Int64" Name="partition:part"/>' in text
    assert '<PDataArray type="Float64" Name="Points" NumberOfComponents="3"/>' in text
    for i in range(3):
        assert f'<Piece Source="g/g_{i:04d}.vtu"/>' in text
        assert (tmp_path / "g" / f"g_{i:04d}.vtu").is_file()


def test_pieces_are_independently_readable_vtu_files(engine, tmp_path):
    mesh = _labelled()
    engine.write(tmp_path / "g.pvtu", mesh)
    total = 0
    for i in range(3):
        piece = meshioplusplus.read(tmp_path / "g" / f"g_{i:04d}.vtu")
        assert set(piece.cell_data["partition:part"][0].tolist()) <= {i}
        total += _num_cells(piece)
    assert total == _num_cells(mesh)


def test_partition_then_pvtu_read_merge_is_the_original_up_to_point_order(
    engine, tmp_path
):
    """The roadmap's "done when"."""
    mesh = _grid(8)
    mesh.cell_data["partition:part"] = meshioplusplus.partition_labels(mesh, 4)
    engine.write(tmp_path / "p.pvtu", mesh)
    back = engine.read(tmp_path / "p.pvtu")
    welded = meshioplusplus.clean(back, weld=True)
    original = meshioplusplus.Mesh(mesh.points, mesh.cells)
    _assert_same_up_to_point_order(welded, original)


def test_partition_output_writes_directly_as_pieces(engine, tmp_path):
    mesh = _grid(8)
    pieces = meshioplusplus.partition(mesh, 4)
    engine.write_pieces(tmp_path / "p.pvtu", pieces)
    back = engine.read(tmp_path / "p.pvtu")
    assert [r.name for r in back.regions] == [f"piece_{i}" for i in range(4)]
    _assert_same_up_to_point_order(
        meshioplusplus.clean(back, weld=True),
        meshioplusplus.Mesh(mesh.points, mesh.cells),
    )


def test_no_part_array_writes_one_piece(engine, tmp_path):
    mesh = _grid(3)
    engine.write(tmp_path / "one.pvtu", mesh)
    assert (tmp_path / "one.pvtu").read_text().count("<Piece ") == 1
    back = engine.read(tmp_path / "one.pvtu")
    assert back.regions == [] or len(back.regions) == 0  # one piece: no region
    assert len(back.points) == len(mesh.points)
    assert _num_cells(back) == _num_cells(mesh)


def test_part_key_selects_the_carving_array(engine, tmp_path):
    mesh = _grid(4)
    mesh.cell_data["rank"] = [np.arange(_num_cells(mesh)) % 2]
    engine.write(tmp_path / "k.pvtu", mesh, part_key="rank")
    assert (tmp_path / "k.pvtu").read_text().count("<Piece ") == 2
    engine.write(tmp_path / "n.pvtu", mesh, part_key="partition:part")  # absent
    assert (tmp_path / "n.pvtu").read_text().count("<Piece ") == 1


def test_non_integer_part_array_is_refused(engine, tmp_path):
    mesh = _grid(3)
    mesh.cell_data["partition:part"] = [np.zeros(_num_cells(mesh), dtype=np.float64)]
    with pytest.raises(WriteError, match="integer part ids"):
        engine.write(tmp_path / "f.pvtu", mesh)


def test_parts_with_no_cells_are_kept(engine, tmp_path):
    """Part 1 owns nothing, but skipping it would renumber part 2."""
    mesh = _grid(3)
    labels = np.zeros(_num_cells(mesh), dtype=np.int64)
    labels[::2] = 2
    mesh.cell_data["partition:part"] = [labels]
    engine.write(tmp_path / "z.pvtu", mesh)
    assert (tmp_path / "z.pvtu").read_text().count("<Piece ") == 3
    back = engine.read(tmp_path / "z.pvtu")
    assert [(r.name, len(r.entries)) for r in back.regions][1] == ("piece_1", 0)
    assert _num_cells(back) == _num_cells(mesh)
    assert _num_cells(engine.read(tmp_path / "z.pvtu", piece=1)) == 0


def test_piece_selection(engine, tmp_path):
    mesh = _labelled()
    engine.write(tmp_path / "s.pvtu", mesh)
    sizes = [int(np.sum(mesh.cell_data["partition:part"][0] == i)) for i in range(3)]
    for k in range(3):
        one = engine.read(tmp_path / "s.pvtu", piece=k)
        assert _num_cells(one) == sizes[k]
        assert one.regions == [] or len(one.regions) == 0
    assert _num_cells(engine.read(tmp_path / "s.pvtu", piece=-1)) == sizes[2]
    with pytest.raises(ReadError, match=r"piece 3 is out of range.*3 pieces"):
        engine.read(tmp_path / "s.pvtu", piece=3)
    with pytest.raises(ReadError, match=r"piece -4 is out of range"):
        engine.read(tmp_path / "s.pvtu", piece=-4)


def test_piece_selection_through_generic_read(tmp_path):
    mesh = _labelled()
    meshioplusplus.write(tmp_path / "g.pvtu", mesh)
    one = meshioplusplus.read(tmp_path / "g.pvtu", piece=1)
    assert _num_cells(one) == int(np.sum(mesh.cell_data["partition:part"][0] == 1))


# --- declarations ------------------------------------------------------------


def _piece(dtype=np.float64, ncomp=1, name="u", cells=1):
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0.0]])
    tri = np.array([[0, 1, 2]] * cells)
    shape = (3,) if ncomp == 1 else (3, ncomp)
    return meshioplusplus.Mesh(
        pts,
        [("triangle", tri)],
        point_data={name: np.zeros(shape, dtype=dtype)},
        cell_data={"c": [np.zeros(cells)]},
    )


@pytest.mark.parametrize(
    "bad, match",
    [
        (
            _piece(np.float32),
            r"piece 1 declares point_data 'u' as Float32 with 1 component, but piece 0 declares it as Float64 with 1 component",
        ),
        (
            _piece(ncomp=3),
            r"point_data 'u' as Float64 with 3 components, but piece 0 declares it as Float64 with 1 component",
        ),
        (
            _piece(name="v"),
            r"piece 1 is missing point_data 'u', which piece 0 declares",
        ),
    ],
)
def test_declaration_mismatch_refuses_before_writing_anything(
    engine, bad, match, tmp_path
):
    with pytest.raises(WriteError, match=match) as info:
        engine.write_pieces(tmp_path / "bad.pvtu", [_piece(), bad])
    assert "every piece of a parallel index must declare identical arrays" in str(
        info.value
    )
    assert not (tmp_path / "bad").exists()
    assert not (tmp_path / "bad.pvtu").exists()


def test_extra_array_in_a_later_piece_is_refused(engine, tmp_path):
    extra = _piece()
    extra.point_data["w"] = np.zeros(3)
    with pytest.raises(
        WriteError, match="declares point_data 'w', which piece 0 does not"
    ):
        engine.write_pieces(tmp_path / "x.pvtu", [_piece(), extra])


def test_an_empty_piece_list_is_refused(engine, tmp_path):
    with pytest.raises(WriteError, match="at least one piece"):
        engine.write_pieces(tmp_path / "e.pvtu", [])


def test_a_cell_free_piece_may_omit_cell_arrays(engine, tmp_path):
    """An idle rank with no cell arrays at all is exempt from the cell_data check."""
    idle = meshioplusplus.Mesh(
        np.zeros((0, 3)),
        [("triangle", np.zeros((0, 3), dtype=np.int64))],
        point_data={"u": np.zeros(0)},
    )
    engine.write_pieces(tmp_path / "i.pvtu", [_piece(), idle, _piece()])
    back = engine.read(tmp_path / "i.pvtu")
    assert _num_cells(back) == 2
    assert [len(r.entries) for r in back.regions] == [1, 0, 1]


# --- ghosts ------------------------------------------------------------------


def _ghosted(nparts=3, layers=1):
    if not hasattr(_core, "partition"):
        pytest.skip("ghost layers need the C++ core")
    mesh = _grid(6)
    return mesh, meshioplusplus.partition(mesh, nparts, ghost_layers=layers)


def test_ghost_layers_become_vtkGhostType_and_ghost_level(engine, tmp_path):
    mesh, pieces = _ghosted()
    engine.write_pieces(tmp_path / "g.pvtu", pieces)
    text = (tmp_path / "g.pvtu").read_text()
    assert '<PUnstructuredGrid GhostLevel="1">' in text
    assert text.count('<PDataArray type="UInt8" Name="vtkGhostType"/>') == 2  # p + c

    piece = meshioplusplus.read(tmp_path / "g" / "g_0000.vtu")
    cell_gh = piece.cell_data["vtkGhostType"][0]
    point_gh = piece.point_data["vtkGhostType"]
    assert cell_gh.dtype == np.uint8 and point_gh.dtype == np.uint8
    layers = pieces[0].cell_data["partition:ghost"][0]
    np.testing.assert_array_equal(cell_gh, (layers > 0).astype(np.uint8))
    # a point is a duplicate exactly when no owned cell of the piece uses it
    owned = np.zeros(len(pieces[0].points), bool)
    owned[pieces[0].cells[0].data[layers == 0].ravel()] = True
    np.testing.assert_array_equal(point_gh, np.where(owned, 0, 1).astype(np.uint8))
    # the layer number itself survives: vtkGhostType alone would collapse it
    assert "partition:ghost" in piece.cell_data


def test_ghosts_are_kept_by_default_and_dropped_on_request(engine, tmp_path):
    mesh, pieces = _ghosted()
    engine.write_pieces(tmp_path / "g.pvtu", pieces)

    kept = engine.read(tmp_path / "g.pvtu")
    assert _num_cells(kept) == sum(_num_cells(p) for p in pieces) > _num_cells(mesh)
    assert "vtkGhostType" in kept.cell_data and "vtkGhostType" in kept.point_data

    dropped = engine.read(tmp_path / "g.pvtu", ghosts="drop")
    assert _num_cells(dropped) == _num_cells(mesh)
    assert "vtkGhostType" not in dropped.cell_data
    assert "vtkGhostType" not in dropped.point_data
    assert "partition:ghost" not in dropped.cell_data
    _assert_same_up_to_point_order(
        meshioplusplus.clean(dropped, weld=True),
        meshioplusplus.Mesh(mesh.points, mesh.cells),
    )


def test_ghost_policy_is_validated(engine, tmp_path):
    engine.write(tmp_path / "g.pvtu", _grid(3))
    with pytest.raises(ValueError, match="ghosts must be 'keep' or 'drop'"):
        engine.read(tmp_path / "g.pvtu", ghosts="maybe")


def test_piece_selection_drops_ghosts_too(engine, tmp_path):
    mesh, pieces = _ghosted()
    engine.write_pieces(tmp_path / "g.pvtu", pieces)
    one = engine.read(tmp_path / "g.pvtu", piece=0, ghosts="drop")
    owned = int(np.sum(pieces[0].cell_data["partition:ghost"][0] == 0))
    assert _num_cells(one) == owned


def test_a_supplied_vtkGhostType_is_passed_through(engine, tmp_path):
    piece = _piece(cells=2)
    piece.cell_data["vtkGhostType"] = [np.array([0, 8], dtype=np.uint8)]  # REFINEDCELL
    engine.write_pieces(tmp_path / "s.pvtu", [piece])
    back = engine.read(tmp_path / "s.pvtu")
    np.testing.assert_array_equal(back.cell_data["vtkGhostType"][0], [0, 8])
    assert _num_cells(engine.read(tmp_path / "s.pvtu", ghosts="drop")) == 1


def test_writing_does_not_mutate_the_callers_pieces(engine, tmp_path):
    _mesh, pieces = _ghosted()
    engine.write_pieces(tmp_path / "g.pvtu", pieces)
    assert all("vtkGhostType" not in p.cell_data for p in pieces)
    assert all("vtkGhostType" not in p.point_data for p in pieces)


# --- reading hand-made and foreign files -------------------------------------


def _write_vtu(path, points, tris):
    """A minimal ASCII .vtu, so a test controls exactly what a piece holds."""
    pts = " ".join(str(v) for row in points for v in row)
    conn = " ".join(str(v) for row in tris for v in row)
    offs = " ".join(str(3 * (i + 1)) for i in range(len(tris)))
    types = " ".join("5" for _ in tris)
    path.write_text(
        f'<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" '
        f'byte_order="LittleEndian">\n<UnstructuredGrid>\n'
        f'<Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(tris)}">\n'
        f"<Points>"
        f'<DataArray type="Float64" NumberOfComponents="3" format="ascii">{pts}'
        f"</DataArray></Points>\n<Cells>"
        f'<DataArray type="Int64" Name="connectivity" format="ascii">{conn}</DataArray>'
        f'<DataArray type="Int64" Name="offsets" format="ascii">{offs}</DataArray>'
        f'<DataArray type="UInt8" Name="types" format="ascii">{types}</DataArray>'
        "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n"
    )


def _write_index(path, sources, root="PUnstructuredGrid", extra=""):
    pieces = "".join(f"<Piece Source={s!r}/>".replace("'", '"') for s in sources)
    path.write_text(
        f'<?xml version="1.0"?>\n<VTKFile type="{root}" version="0.1" '
        f'byte_order="LittleEndian" header_type="UInt64">\n<{root} GhostLevel="0">\n'
        f'{extra}<PPoints><PDataArray type="Float64" NumberOfComponents="3"/></PPoints>\n'
        f"{pieces}\n</{root}>\n</VTKFile>\n"
    )


TRI = [[0, 1, 2]]
XY = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]


def test_piece_paths_with_spaces_and_entities(engine, tmp_path):
    d = tmp_path / "my pieces & more"
    d.mkdir()
    _write_vtu(d / "piece one.vtu", XY, TRI)
    _write_index(tmp_path / "i.pvtu", ["my pieces &amp; more/piece one.vtu"])
    back = engine.read(tmp_path / "i.pvtu")
    assert _num_cells(back) == 1


def test_an_index_with_a_space_in_its_own_name_round_trips(engine, tmp_path):
    mesh = _labelled()
    engine.write(tmp_path / "my case.pvtu", mesh)
    assert (tmp_path / "my case" / "my case_0000.vtu").is_file()
    assert _num_cells(engine.read(tmp_path / "my case.pvtu")) == _num_cells(mesh)


def test_header_type_and_byte_order_may_differ_between_index_and_pieces(
    engine, tmp_path
):
    """The index's header_type says UInt64; the piece is UInt32 (the default)."""
    _write_vtu(tmp_path / "a.vtu", XY, TRI)
    _write_index(tmp_path / "i.pvtu", ["a.vtu"])  # header_type="UInt64"
    assert _num_cells(engine.read(tmp_path / "i.pvtu")) == 1


def test_a_missing_relative_piece_names_the_attribute_and_the_file(engine, tmp_path):
    _write_index(tmp_path / "i.pvtu", ["gone/piece.vtu"])
    with pytest.raises(ReadError, match=r"gone/piece\.vtu.*Source=.*does not exist"):
        engine.read(tmp_path / "i.pvtu")


def test_an_absolute_path_from_another_machine_is_an_error_not_a_guess(
    engine, tmp_path
):
    """A sibling with the same basename must not be read instead."""
    _write_vtu(tmp_path / "piece.vtu", XY, TRI)
    _write_index(tmp_path / "i.pvtu", ["/scratch/run7/piece.vtu"])
    with pytest.raises(ReadError, match=r"/scratch/run7/piece\.vtu.*not searched"):
        engine.read(tmp_path / "i.pvtu")


def test_a_piece_that_is_not_vtu_or_vtp_is_refused_by_name(engine, tmp_path):
    (tmp_path / "m.msh").write_text("$MeshFormat\n")
    _write_index(tmp_path / "i.pvtu", ["m.msh"])
    with pytest.raises(ReadError, match=r"unsupported piece.*m\.msh.*\.vtu/\.vtp"):
        engine.read(tmp_path / "i.pvtu")


def test_an_index_is_not_a_piece(engine, tmp_path):
    """No recursion inside a parallel index, so no cycle is expressible."""
    _write_index(tmp_path / "self.pvtu", ["self.pvtu"])
    with pytest.raises(ReadError, match="unsupported piece"):
        engine.read(tmp_path / "self.pvtu")


def test_an_index_naming_no_piece_reads_as_an_empty_mesh(engine, tmp_path):
    _write_index(tmp_path / "e.pvtu", [])
    back = engine.read(tmp_path / "e.pvtu")
    assert len(back.points) == 0 and back.cells == []


def test_wrong_root_type_is_refused(engine, tmp_path):
    _write_index(tmp_path / "w.pvtu", [], root="PPolyData")
    with pytest.raises(ReadError, match="expected type PUnstructuredGrid"):
        engine.read(tmp_path / "w.pvtu")


def test_not_xml_is_refused(engine, tmp_path):
    (tmp_path / "x.pvtu").write_text("not xml at all")
    with pytest.raises(ReadError):
        engine.read(tmp_path / "x.pvtu")


# --- metadata, sniffing and generic dispatch ----------------------------------


def test_metadata_agrees_with_a_real_read(tmp_path):
    mesh = _labelled()
    path = tmp_path / "m.pvtu"
    meshioplusplus.pvtu.write(path, mesh)
    meta = meshioplusplus.read_metadata(path)
    back = meshioplusplus.read(path)
    assert meta["num_points"] == len(back.points)
    assert [b["type"] for b in meta["cell_blocks"]] == [c.type for c in back.cells]
    assert [b["num_cells"] for b in meta["cell_blocks"]] == [
        len(c.data) for c in back.cells
    ]
    assert meta["point_data_names"] == sorted(back.point_data)
    assert meta["cell_data_names"] == sorted(back.cell_data)
    assert meta["format"] == "pvtu"


def test_generic_io_and_sniffing(tmp_path):
    mesh = _labelled()
    path = tmp_path / "g.pvtu"
    meshioplusplus.write(path, mesh)
    assert meshioplusplus.sniff_format(path) == "pvtu"
    back = meshioplusplus.read(path)
    assert _num_cells(back) == _num_cells(mesh)
    # an explicit name works for a file whose extension says nothing
    os.rename(path, tmp_path / "g.dat")
    assert _num_cells(meshioplusplus.read(tmp_path / "g.dat", file_format="pvtu")) == (
        _num_cells(mesh)
    )


@pytest.mark.parametrize(
    "root, expected",
    [
        ("PUnstructuredGrid", "pvtu"),
        ("PPolyData", "pvtp"),
        ("Collection", "pvd"),
        # refused outright rather than mistaken for a serial twin
        ("PImageData", ""),
        ("PStructuredGrid", ""),
        ("PRectilinearGrid", ""),
        # the serial types still sniff as themselves
        ("UnstructuredGrid", "vtu"),
        ("PolyData", "vtp"),
        ("vtkMultiBlockDataSet", "vtm"),
    ],
)
@pytest.mark.parametrize("quote", ['"', "'"])
def test_sniff_index_types_before_their_serial_substrings(
    root, expected, quote, tmp_path
):
    """`PUnstructuredGrid` contains `UnstructuredGrid`; a bare `Collection`
    would also match `vtkPartitionedDataSetCollection`."""
    path = tmp_path / "x.bin"
    path.write_text(f"<?xml version='1.0'?>\n<VTKFile type={quote}{root}{quote}>\n")
    assert meshioplusplus.sniff_format(path) == expected
    from meshioplusplus._sniff import _sniff_format_py

    assert _sniff_format_py(path) == expected


def test_a_partitioned_dataset_collection_is_not_a_pvd(tmp_path):
    path = tmp_path / "x.bin"
    path.write_text('<VTKFile type="vtkPartitionedDataSetCollection">\n')
    assert meshioplusplus.sniff_format(path) != "pvd"


def test_buffers_are_refused_by_name(tmp_path):
    import io

    with pytest.raises(Exception, match="pvtu"):
        meshioplusplus.write(io.StringIO(), _grid(2), file_format="pvtu")


def test_the_index_records_provenance(tmp_path):
    meshioplusplus.pvtu.write(tmp_path / "p.pvtu", _labelled())
    meta = meshioplusplus.read_metadata(tmp_path / "p.pvtu")
    assert meta["provenance_recognised"]


# --- cross-compatibility between the engines ----------------------------------


@pytest.mark.skipif(not hasattr(_core, "pvtu_read"), reason="no pvtu core")
@pytest.mark.parametrize(
    "binary, compression",
    [(False, None), (True, None), pytest.param(True, "zlib", marks=requires_zlib)],
)
def test_cross_compat(binary, compression, tmp_path, strict_core):
    mesh = _labelled()
    (tmp_path / "core").mkdir()
    (tmp_path / "py").mkdir()
    a = tmp_path / "core" / "g.pvtu"
    meshioplusplus.pvtu.write(a, mesh, binary=binary, compression=compression)
    b = tmp_path / "py" / "g.pvtu"
    _pvtu.write(b, mesh, binary=binary, compression=compression)
    for path in (a, b):  # each engine reads the other's file
        from_core = meshioplusplus.pvtu.read(path)
        from_py = _pvtu.read(path)
        assert len(from_core.points) == len(from_py.points)
        assert [r.name for r in from_core.regions] == [r.name for r in from_py.regions]
        _assert_same_up_to_point_order(
            meshioplusplus.clean(from_core, weld=True),
            meshioplusplus.clean(from_py, weld=True),
        )

    def body(p):  # everything but the provenance note
        return [ln for ln in p.read_text().splitlines() if not ln.startswith("<!--")]

    assert body(a) == body(b)


# --- .pvtp -------------------------------------------------------------------


def _polydata(n=3):
    pts = np.array([[i, j, 0.0] for j in range(n + 1) for i in range(n + 1)])
    quads = [
        [
            j * (n + 1) + i,
            j * (n + 1) + i + 1,
            (j + 1) * (n + 1) + i + 1,
            (j + 1) * (n + 1) + i,
        ]
        for j in range(n)
        for i in range(n)
    ]
    mesh = meshioplusplus.Mesh(pts, [("quad", np.array(quads))])
    mesh.point_data["u"] = pts[:, 1]
    mesh.cell_data["partition:part"] = [np.arange(len(quads)) % 2]
    return mesh


@pytest.fixture(params=["core", "python"])
def pvtp_engine(request):
    if request.param == "core":
        if not HAS_ZLIB:
            pytest.skip("the core engine writes zlib and this build has none")
        if not hasattr(_core, "pvtp_read"):
            pytest.skip("this build has no pvtp core")
        set_strict_core(True)
        yield _engines("pvtp")["core"]
        set_strict_core(None)
    else:
        yield _engines("pvtp")["python"]


def test_pvtp_roundtrip(pvtp_engine, tmp_path):
    mesh = _polydata()
    path = tmp_path / "s.pvtp"
    pvtp_engine.write(path, mesh)
    text = path.read_text()
    assert '<VTKFile type="PPolyData"' in text and "<PPolyData GhostLevel" in text
    assert '<Piece Source="s/s_0000.vtp"/>' in text
    assert (tmp_path / "s" / "s_0001.vtp").is_file()
    back = pvtp_engine.read(path)
    assert _num_cells(back) == _num_cells(mesh)
    assert [r.name for r in back.regions] == ["piece_0", "piece_1"]
    assert meshioplusplus.sniff_format(path) == "pvtp"
    assert _num_cells(pvtp_engine.read(path, piece=1)) == _num_cells(mesh) // 2


def test_pvtp_declaration_check_is_shared(pvtp_engine, tmp_path):
    good, bad = _polydata(), _polydata()
    bad.point_data["u"] = bad.point_data["u"].astype(np.float32)
    with pytest.raises(WriteError, match="pvtp: piece 1 declares point_data 'u'"):
        pvtp_engine.write_pieces(tmp_path / "b.pvtp", [good, bad])


# --- the partition command ---------------------------------------------------


def _cli(*argv):
    return meshioplusplus._cli.main([str(a) for a in argv])


def test_partition_cli_writes_one_index_over_every_part_halo_included(tmp_path):
    if not hasattr(_core, "partition"):
        pytest.skip("ghost layers need the C++ core")
    mesh = _grid(6)
    meshioplusplus.vtu.write(tmp_path / "in.vtu", mesh)
    _cli(
        "partition",
        tmp_path / "in.vtu",
        tmp_path / "out.pvtu",
        "-n",
        3,
        "--method",
        "sfc",
        "--ghost-layers",
        1,
    )
    text = (tmp_path / "out.pvtu").read_text()
    assert text.count("<Piece ") == 3 and 'GhostLevel="1"' in text

    kept = meshioplusplus.read(tmp_path / "out.pvtu")
    assert _num_cells(kept) > _num_cells(mesh)
    dropped = meshioplusplus.pvtu.read(tmp_path / "out.pvtu", ghosts="drop")
    assert _num_cells(dropped) == _num_cells(mesh)


def test_partition_cli_pvtp_index_and_explicit_output_format(tmp_path):
    meshioplusplus.vtp.write(tmp_path / "in.vtp", _polydata())
    _cli(
        "partition",
        tmp_path / "in.vtp",
        tmp_path / "out.pvtp",
        "-n",
        2,
        "--method",
        "sfc",
    )
    assert (tmp_path / "out.pvtp").read_text().count("<Piece ") == 2
    # an explicit --output-format makes any path an index
    _cli(
        "partition",
        tmp_path / "in.vtp",
        tmp_path / "out.idx",
        "-n",
        2,
        "--method",
        "sfc",
        "-o",
        "pvtp",
    )
    assert (tmp_path / "out.idx").read_text().count("<Piece ") == 2


def test_partition_cli_token_still_writes_one_file_per_part(tmp_path):
    meshioplusplus.vtu.write(tmp_path / "in.vtu", _grid(4))
    _cli(
        "partition",
        tmp_path / "in.vtu",
        tmp_path / "p_{part}.vtu",
        "-n",
        2,
        "--method",
        "sfc",
    )
    assert (tmp_path / "p_0.vtu").is_file() and (tmp_path / "p_1.vtu").is_file()
    # the token wins over an index extension: each part is its own one-piece index
    _cli(
        "partition",
        tmp_path / "in.vtu",
        tmp_path / "q_{part}.pvtu",
        "-n",
        2,
        "--method",
        "sfc",
    )
    assert (tmp_path / "q_0.pvtu").read_text().count("<Piece ") == 1


def test_partition_cli_refusal_names_the_index_form(tmp_path):
    meshioplusplus.vtu.write(tmp_path / "in.vtu", _grid(3))
    with pytest.raises(ValueError, match=r"\.pvtu/\.pvtp index"):
        _cli(
            "partition",
            tmp_path / "in.vtu",
            tmp_path / "out.vtu",
            "-n",
            2,
            "--method",
            "sfc",
        )


# --- polyhedra: the core's job ----------------------------------------------


def _polyhedra(part):
    """The shared polyhedron fixture, with `part` (a function of the block's
    cell count) attached as an integer cell array."""
    src = helpers.polyhedron_mesh
    mesh = meshioplusplus.Mesh(src.points, src.cells)
    return mesh, [part(len(c.data)) for c in mesh.cells]


@requires_zlib
@pytest.mark.skipif(not hasattr(_core, "pvtu_read"), reason="no pvtu core")
def test_the_core_carves_ghosts_and_drops_polyhedron_blocks(tmp_path, strict_core):
    mesh, parts = _polyhedra(lambda n: np.arange(n) % 2)
    total = _num_cells(mesh)
    mesh.cell_data["partition:part"] = parts
    meshioplusplus.pvtu.write(tmp_path / "p.pvtu", mesh)
    back = meshioplusplus.pvtu.read(tmp_path / "p.pvtu")
    assert _num_cells(back) == total
    assert [r.name for r in back.regions] == ["piece_0", "piece_1"]

    ghosted, layers = _polyhedra(lambda n: np.arange(n) % 2)
    ghosted.cell_data["partition:ghost"] = layers
    meshioplusplus.pvtu.write_pieces(tmp_path / "g.pvtu", [ghosted])
    kept = meshioplusplus.pvtu.read(tmp_path / "g.pvtu")
    assert _num_cells(kept) == total and "vtkGhostType" in kept.cell_data
    dropped = meshioplusplus.pvtu.read(tmp_path / "g.pvtu", ghosts="drop")
    assert _num_cells(dropped) == sum(int(np.sum(x == 0)) for x in layers)


def test_the_python_reference_refuses_polyhedron_blocks_by_name(tmp_path):
    mesh, parts = _polyhedra(lambda n: np.arange(n) % 2)
    mesh.cell_data["partition:part"] = parts
    with pytest.raises(NotImplementedError, match=r"needs the C\+\+ core"):
        _pvtu.write(tmp_path / "x.pvtu", mesh)


# --- the ghost switch on the generic read surfaces (v15.0.0) ------------------


def _ghosted_index(tmp_path, name="g.pvtu"):
    if not hasattr(_core, "partition"):
        pytest.skip("ghost layers need the C++ core")
    mesh = _grid(6)
    pieces = meshioplusplus.partition(mesh, 3, ghost_layers=1)
    meshioplusplus.pvtu.write_pieces(tmp_path / name, pieces)
    return mesh, tmp_path / name


def test_the_generic_read_takes_ghosts(tmp_path):
    mesh, index = _ghosted_index(tmp_path)
    kept = meshioplusplus.read(index)
    dropped = meshioplusplus.read(index, ghosts="drop")
    assert _num_cells(kept) > _num_cells(mesh)
    assert _num_cells(dropped) == _num_cells(mesh)
    assert "vtkGhostType" not in dropped.cell_data
    # an explicit "keep" is the default
    assert _num_cells(meshioplusplus.read(index, ghosts="keep")) == _num_cells(kept)


def test_the_generic_read_refuses_a_bad_policy_by_name(tmp_path):
    _mesh, index = _ghosted_index(tmp_path)
    with pytest.raises(ValueError, match="ghosts must be 'keep' or 'drop'"):
        meshioplusplus.read(index, ghosts="maybe")


def test_every_other_reader_ignores_the_ghost_policy(tmp_path):
    """A file with no halo is already the answer "drop" asks for."""
    mesh = _grid(3)
    meshioplusplus.vtu.write(tmp_path / "a.vtu", mesh)
    assert _num_cells(meshioplusplus.read(tmp_path / "a.vtu", ghosts="drop")) == (
        _num_cells(mesh)
    )
    # ... and it does not swallow the options that reader *does* take
    assert sorted(
        meshioplusplus.read(tmp_path / "a.vtu", ghosts="drop", arrays=["u"]).point_data
    ) == ["u"]


def test_the_ghost_policy_and_piece_reach_a_reader_together(tmp_path):
    mesh, index = _ghosted_index(tmp_path)
    one = meshioplusplus.read(index, piece=0, ghosts="drop")
    assert 0 < _num_cells(one) < _num_cells(mesh)
    assert "vtkGhostType" not in one.cell_data


def test_read_sequence_forwards_the_ghost_policy(tmp_path):
    mesh, index = _ghosted_index(tmp_path)
    (tmp_path / "run.pvd").write_text(
        '<VTKFile type="Collection" version="0.1"><Collection>'
        f'<DataSet timestep="0" file="{index.name}"/></Collection></VTKFile>'
    )
    ((_t, back),) = meshioplusplus.read_sequence(tmp_path / "run.pvd", ghosts="drop")
    assert _num_cells(back) == _num_cells(mesh)


def test_convert_cli_drop_ghosts(tmp_path):
    mesh, index = _ghosted_index(tmp_path)
    _cli("convert", index, tmp_path / "kept.vtu")
    _cli("convert", index, tmp_path / "dropped.vtu", "--drop-ghosts")
    assert _num_cells(meshioplusplus.read(tmp_path / "kept.vtu")) > _num_cells(mesh)
    dropped = meshioplusplus.read(tmp_path / "dropped.vtu")
    assert _num_cells(dropped) == _num_cells(mesh)
    assert "vtkGhostType" not in dropped.cell_data


def test_mcp_tools_take_ghosts(tmp_path):
    from meshioplusplus.mcp import _tools

    mesh, index = _ghosted_index(tmp_path)
    kept = _tools.tool_stats(str(index))
    dropped = _tools.tool_stats(str(index), ghosts="drop")
    assert kept["num_cells"] > dropped["num_cells"] == _num_cells(mesh)
    out = _tools.tool_convert(str(index), str(tmp_path / "o.vtu"), ghosts="drop")
    assert "o.vtu" in str(out)
    assert _num_cells(meshioplusplus.read(tmp_path / "o.vtu")) == _num_cells(mesh)


# --- field data is dataset-global (v15.0.0) ------------------------------------


def _timed(cells=1):
    piece = _piece(cells=cells)
    piece.field_data["TimeValue"] = np.array([0.25])
    return piece


def test_field_data_is_the_union_across_pieces_not_namespaced(engine, tmp_path):
    """merge() would rename the copies `0:TimeValue`, `1:TimeValue`, ..."""
    a, b, c = _timed(), _timed(), _timed()
    b.field_data["only_in_b"] = np.arange(2.0)
    engine.write_pieces(tmp_path / "f.pvtu", [a, b, c])
    back = engine.read(tmp_path / "f.pvtu")
    assert sorted(back.field_data) == ["TimeValue", "only_in_b"]
    assert float(back.field_data["TimeValue"][0]) == 0.25
    meta = meshioplusplus.read_metadata(tmp_path / "f.pvtu")
    assert meta["field_data_names"] == sorted(back.field_data)


def test_a_single_write_carries_its_field_data_through_the_pieces(engine, tmp_path):
    mesh = _labelled()
    mesh.field_data["TimeValue"] = np.array([1.5])
    engine.write(tmp_path / "w.pvtu", mesh)
    back = engine.read(tmp_path / "w.pvtu")
    assert sorted(back.field_data) == ["TimeValue"]
    assert float(back.field_data["TimeValue"][0]) == 1.5


def test_vtkGhostType_is_always_uint8_on_disk_whatever_dtype_it_is_given(
    engine, tmp_path
):
    """ParaView takes vtkGhostType as the ghost array only when it is an unsigned
    char array (an Int64 one is ignored); the NATIVE and KRATOS backends hold every
    integer array as Int64, and a caller can hand any integer dtype to either
    engine, so the writers and the index declaration pin the on-disk type."""
    piece = _piece(cells=2)
    piece.cell_data["vtkGhostType"] = [np.array([0, 8], dtype=np.int64)]
    piece.point_data["vtkGhostType"] = np.array([0, 1, 0], dtype=np.int32)
    for binary in (False, True):
        engine.write_pieces(tmp_path / "g.pvtu", [piece], binary=binary)
        index = (tmp_path / "g.pvtu").read_text()
        assert index.count('<PDataArray type="UInt8" Name="vtkGhostType"/>') == 2
        text = (tmp_path / "g" / "g_0000.vtu").read_text()
        assert text.count('type="UInt8" Name="vtkGhostType"') == 2
        back = engine.read(tmp_path / "g.pvtu")
        assert back.cell_data["vtkGhostType"][0].tolist() == [0, 8]  # the bit survives
        assert back.point_data["vtkGhostType"].tolist() == [0, 1, 0]
    # ... and the caller's arrays keep their dtype
    assert piece.cell_data["vtkGhostType"][0].dtype == np.int64
    assert piece.point_data["vtkGhostType"].dtype == np.int32


@pytest.mark.parametrize("compression", [None, "zlib"])
def test_the_python_vtu_reader_reads_a_cell_free_piece(compression, tmp_path):
    """An empty piece has empty (block-less, when compressed) connectivity."""
    from meshioplusplus.vtu import _vtu

    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]), [("line", np.empty((0, 2), int))]
    )
    _vtu.write(tmp_path / "e.vtu", mesh, compression=compression)
    back = _vtu.read(tmp_path / "e.vtu")
    assert len(back.points) == 2 and _num_cells(back) == 0
