"""GiD postprocess format tests.

The reader is exercised against **hand-authored** fixtures rather than files
produced by meshio++'s own writer wherever the point is real-world tolerance.
That is deliberate: a fixture our writer produced makes the round trip a weak
oracle, and two of the variants covered here (a node table repeated in every
``MESH`` block, element ids restarting per block) are ones our writer
structurally cannot emit at all. Both were observed in a genuine
Kratos-produced GiD file.

There is no pure-Python GiD engine, so unlike most format tests here there is
no ``cpp``/``python`` parametrization -- the C++ core is the only engine.
"""

import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core

MESHES = pathlib.Path(__file__).resolve().parent / "meshes" / "gid"

# The hdf5 flavour needs an HDF5-enabled build; the binary flavour needs zlib.
_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)
_HAS_ZLIB = getattr(_core, "__has_zlib__", False)
_HAS_GIDPOST = getattr(_core, "__has_gidpost__", False)

FLAVOURS = [
    pytest.param(".post.msh", id="ascii"),
    pytest.param(
        ".post.bin",
        id="binary",
        marks=pytest.mark.skipif(not _HAS_ZLIB, reason="needs zlib"),
    ),
    pytest.param(
        ".post.h5",
        id="hdf5",
        marks=pytest.mark.skipif(not _HAS_HDF5, reason="needs HDF5"),
    ),
]

needs_writer = pytest.mark.skipif(
    not _HAS_GIDPOST, reason="writing GiD needs a gidpost-enabled build"
)


def _sample_mesh():
    points = np.array(
        [[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1]], dtype=float
    )
    return meshioplusplus.Mesh(
        points,
        [
            ("tetra", np.array([[0, 1, 2, 4]], dtype=np.int64)),
            ("triangle", np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int64)),
        ],
    )


# ---------------------------------------------------------------------------
# Round trips, all three flavours.


@needs_writer
@pytest.mark.parametrize("ext", FLAVOURS)
def test_roundtrip_geometry(tmp_path, ext):
    mesh = _sample_mesh()
    path = tmp_path / ("rt" + ext)
    meshioplusplus.write(path, mesh)
    back = meshioplusplus.read(path)

    # The binary flavour narrows every real to a 4-byte float on the way out
    # (gidpost's own CPostBinary_WriteDouble does), so its tolerance is
    # float32's, not float64's. That is a property of the format, not a defect.
    atol = 1.0e-6 if ext == ".post.bin" else 1.0e-8
    assert back.points.shape == mesh.points.shape
    assert np.allclose(back.points, mesh.points, atol=atol, rtol=0.0)
    assert [c.type for c in back.cells] == [c.type for c in mesh.cells]
    for got, want in zip(back.cells, mesh.cells):
        assert np.array_equal(got.data, want.data)


@needs_writer
@pytest.mark.parametrize("ext", FLAVOURS)
def test_roundtrip_data(tmp_path, ext):
    mesh = _sample_mesh()
    mesh.point_data["T"] = np.arange(len(mesh.points), dtype=float)
    mesh.cell_data["q"] = [np.arange(len(cb.data), dtype=float) for cb in mesh.cells]
    path = tmp_path / ("rt" + ext)
    meshioplusplus.write(path, mesh)
    back = meshioplusplus.read(path)

    atol = 1.0e-6 if ext == ".post.bin" else 1.0e-8
    assert np.allclose(back.point_data["T"], mesh.point_data["T"], atol=atol)
    # A cell_data array spanning several blocks is written as several Result
    # blocks sharing one name; the reader must MERGE them, not let the last
    # one win.
    assert "q" in back.cell_data
    for got, want in zip(back.cell_data["q"], mesh.cell_data["q"]):
        assert np.allclose(got, want, atol=atol)


@needs_writer
@pytest.mark.parametrize("ext", FLAVOURS)
def test_material_column_roundtrips_as_gmsh_physical(tmp_path, ext):
    mesh = _sample_mesh()
    mesh.cell_data["gmsh:physical"] = [
        np.full(len(cb.data), i + 3, dtype=np.int64) for i, cb in enumerate(mesh.cells)
    ]
    path = tmp_path / ("mat" + ext)
    meshioplusplus.write(path, mesh)
    back = meshioplusplus.read(path)
    for got, want in zip(
        back.cell_data["gmsh:physical"], mesh.cell_data["gmsh:physical"]
    ):
        assert np.array_equal(got, want)


@needs_writer
@pytest.mark.parametrize("ext", FLAVOURS)
def test_no_material_invents_nothing(tmp_path, ext):
    """The binary and HDF5 writers always emit a material column, ASCII does not.

    An all-zero column means "no materials", so surfacing it as a
    ``gmsh:physical`` array would invent data the source mesh never carried --
    and would make the three flavours disagree on a round trip.
    """
    path = tmp_path / ("nomat" + ext)
    meshioplusplus.write(path, _sample_mesh())
    assert "gmsh:physical" not in meshioplusplus.read(path).cell_data


# ---------------------------------------------------------------------------
# Real-world variants, from hand-authored fixtures.


def test_repeated_node_table_is_deduplicated():
    """Real files repeat the full node table in every MESH block.

    meshio++'s own writer emits it once and writes empty Coordinates pairs
    thereafter, so this shape is unreachable through our writer -- it comes
    from a genuine Kratos-produced file.
    """
    mesh = meshioplusplus.read(MESHES / "repeated_nodes.post.msh")
    assert mesh.points.shape == (5, 3)  # 5, not 10: de-duplicated by node id
    assert [c.type for c in mesh.cells] == ["tetra", "triangle"]
    assert mesh.cells[0].data.shape == (2, 4)
    assert mesh.cells[1].data.shape == (3, 3)


def test_element_ids_may_restart_per_block():
    """Element ids are NOT globally unique in real files.

    The fixture numbers both blocks from 1, as Kratos does; our own writer
    deliberately assigns a single running id across every block. A reader that
    kept one global id map would mis-associate results.
    """
    mesh = meshioplusplus.read(MESHES / "repeated_nodes.post.msh")
    # Both blocks read back at full length despite sharing ids 1 and 2.
    assert len(mesh.cells[0].data) == 2
    assert len(mesh.cells[1].data) == 3


def test_trailing_material_column_is_read():
    """`Nnode` is the only disambiguator for the optional trailing material id.

    There is no separator between connectivity and material, so a row is
    ``1 + Nnode`` or ``1 + Nnode + 1`` tokens wide and nothing else.
    """
    mesh = meshioplusplus.read(MESHES / "repeated_nodes.post.msh")
    assert np.array_equal(mesh.cell_data["gmsh:physical"][0], [7, 7])
    assert np.array_equal(mesh.cell_data["gmsh:physical"][1], [9, 9, 9])


def test_results_are_read():
    mesh = meshioplusplus.read(MESHES / "results.post.msh")
    assert np.array_equal(mesh.point_data["T"], [10, 20, 30, 40])
    # The Result header carries no ":N" suffix, so the component count comes
    # from the row width -- three here, despite the bare "Vector" declaration.
    assert mesh.point_data["flux"].shape == (4, 3)
    assert np.array_equal(mesh.cell_data["q"][0], [5, 6])


@pytest.mark.parametrize("step,expected", [(0, 10.0), (1, 100.0), (-1, 100.0)])
def test_time_step_selects_a_step(step, expected):
    """Multi-step results: mTimeStep cannot be emulated after the fact."""
    mesh = meshioplusplus.read(MESHES / "results.post.msh", time_step=step)
    assert mesh.point_data["T"][0] == expected


def test_out_of_range_time_step_raises():
    # Through gid.read directly, not meshioplusplus.read: `.post.msh` also
    # lists the shorter-suffix `.msh` candidates (ansys/gmsh/freefem), and
    # _read_file's candidate loop folds a specific ReadError into an aggregate
    # "as either of ..." message. The specific error is what we want to pin.
    with pytest.raises(Exception, match="out of range"):
        meshioplusplus.gid.read(MESHES / "results.post.msh", time_step=9)


def test_multi_gauss_point_result_is_dropped_not_guessed():
    """meshio++'s cell_data is (n,)/(n,k), never per-node-within-cell.

    The same structural limit MED's ELNO/ELGA documents. The rows still have to
    PARSE -- they exercise gidpost's id-suppression rule, where a repeated id is
    omitted and the row begins with whitespace -- and are then dropped rather
    than averaged, which would invent data.
    """
    mesh = meshioplusplus.read(MESHES / "results.post.msh")
    assert "twogp" not in mesh.cell_data
    assert "q" in mesh.cell_data  # the G=1 result alongside it is unaffected


# ---------------------------------------------------------------------------
# Sibling policy, dispatch and errors.


@needs_writer
def test_results_sibling_is_optional(tmp_path):
    """A geometry file with no results reads back as geometry only."""
    path = tmp_path / "geo.post.msh"
    meshioplusplus.write(path, _sample_mesh())
    (tmp_path / "geo.post.res").unlink()
    back = meshioplusplus.read(path)
    assert back.points.shape == (5, 3)
    assert not back.point_data


@needs_writer
def test_reading_via_the_res_sibling_finds_the_geometry(tmp_path):
    path = tmp_path / "pair.post.msh"
    mesh = _sample_mesh()
    mesh.point_data["T"] = np.arange(len(mesh.points), dtype=float)
    meshioplusplus.write(path, mesh)
    back = meshioplusplus.read(tmp_path / "pair.post.res")
    assert back.points.shape == (5, 3)
    assert np.allclose(back.point_data["T"], mesh.point_data["T"])


def test_missing_geometry_file_raises(tmp_path):
    with pytest.raises(Exception):
        meshioplusplus.read(tmp_path / "absent.post.msh")


def test_gid_is_now_readable():
    from meshioplusplus._helpers import reader_map

    assert "gid" in reader_map


@pytest.mark.parametrize(
    "name,expected",
    [
        ("x.post.msh", "gid"),
        ("x.post.res", "gid"),
        ("x.post.bin", "gid"),
        ("x.post.h5", "gid"),
        # Unchanged by the reader's arrival.
        ("x.post", "permas"),
        ("x.dato", "permas"),
        ("x.vol.gz", "netgen"),
        ("x.post.gz", "permas"),
        ("x.msh", "ansys"),
    ],
)
def test_extension_dispatch_is_unchanged(name, expected):
    from meshioplusplus._helpers import _filetypes_from_path

    assert _filetypes_from_path(pathlib.Path(name))[0] == expected


def test_unsupported_cell_type_raises_by_name(tmp_path):
    """Sphere/Circle rows carry a radius, not a node list -- a genuinely and
    permanently unsupported row shape, not merely an unverified ordering.

    hexahedron27/wedge15/pyramid13 used to be this test's fixture; they are
    now supported (Kratos-derived orderings) and read successfully -- see
    test_gid_cpp's GidOrdering suite and the C++ twin of this test.
    """
    path = tmp_path / "sphere.post.msh"
    path.write_text(
        'MESH "m" dimension 3 ElemType Sphere Nnode 1\n'
        "Coordinates\n1 0 0 0\nEnd Coordinates\n"
        "Elements\nEnd Elements\n"
    )
    # Direct, for the same reason as test_out_of_range_time_step_raises.
    with pytest.raises(Exception, match="Sphere"):
        meshioplusplus.gid.read(path)


@pytest.mark.parametrize(
    "cell_type,nnode",
    [("hexahedron27", 27), ("wedge15", 15), ("pyramid13", 13)],
)
def test_kratos_derived_orderings_round_trip(tmp_path, cell_type, nnode):
    """hexahedron27/wedge15/pyramid13, closing the last roadmap item for this
    format: orderings derived from Kratos's own geometry classes
    (kratos/geometries/hexahedra_3d_27.h, prism_3d_15.h, pyramid_3d_13.h),
    cross-checked against Kratos's Element-agnostic vtk_output.cpp/
    ensight_output.cpp conversion. The real oracle is the C++ GidOrdering
    bytes tests (which check the raw file against known geometry, not a round
    trip); this is the cheaper regression lock over the same path.
    """
    pts = np.arange(nnode * 3, dtype=float).reshape(nnode, 3)
    mesh = meshioplusplus.Mesh(pts, [(cell_type, np.arange(nnode).reshape(1, nnode))])
    path = tmp_path / f"{cell_type}.post.msh"
    meshioplusplus.gid.write(str(path), mesh)
    back = meshioplusplus.gid.read(str(path))
    assert back.cells[0].type == cell_type
    assert back.cells[0].data.shape == (1, nnode)


def test_buffer_is_refused():
    import io

    with pytest.raises(Exception, match="multiple files"):
        meshioplusplus.read(io.BytesIO(b""), file_format="gid")


# ---------------------------------------------------------------------------
# Grammar conformance against CIMNE's published specification.
#
# These three cases come from the *current* official grammar rather than from
# gidpost's behaviour, and none of them is reachable through our own writer --
# gidpost emits exactly one casing, always writes a mesh name, and always
# spells the 1-D type "Linear". So a round trip cannot exercise any of them,
# which is precisely why they are hand-authored.
# ---------------------------------------------------------------------------


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_keywords_are_case_insensitive(tmp_path):
    """CIMNE: "keywords ... are not case-sensitive".

    The specification's own worked example opens with ``Coordinates`` and
    closes with ``end coordinates``, so a case-sensitive reader rejects the
    grammar's own example file.
    """
    path = _write(
        tmp_path,
        "mixed.post.msh",
        'mesh "lower" DIMENSION 3 elemtype Triangle nnode 3\n'
        "coordinates\n1 0.0 0.0 0.0\n2 1.0 0.0 0.0\n3 0.0 1.0 0.0\nEND COORDINATES\n"
        "Elements\n1 1 2 3\nend elements\n",
    )
    mesh = meshioplusplus.gid.read(str(path))
    assert len(mesh.points) == 3
    assert mesh.cells[0].type == "triangle"


def test_mesh_name_is_optional(tmp_path):
    """A nameless header must parse. ``MESH dimension 3 ElemType Linear Nnode
    2`` is the specification's own example."""
    path = _write(
        tmp_path,
        "noname.post.msh",
        "MESH    dimension 3 ElemType Triangle  Nnode 3\n"
        "Coordinates\n1 0 0 0\n2 1 0 0\n3 0 1 0\nend coordinates\n"
        "Elements\n1 1 2 3\nend elements\n",
    )
    mesh = meshioplusplus.gid.read(str(path))
    assert len(mesh.points) == 3
    assert mesh.cells[0].type == "triangle"


def test_a_nameless_mesh_is_not_named_dimension(tmp_path):
    """The observable consequence of the optional-name rule.

    Parsing a nameless header alone cannot distinguish the two readings -- the
    ``ElemType``/``Nnode`` scan finds its keywords either way -- so this pins
    the one place the name is actually *used*: binding a ``GaussPoints`` set to
    a block by mesh name. Taking ``tok[1]`` unconditionally names the mesh
    "dimension", and a set declared ``OnMesh "dimension"`` then binds to it and
    attaches a result that does not belong to it.
    """
    _write(
        tmp_path,
        "nm.post.msh",
        "MESH dimension 3 ElemType Triangle Nnode 3\n"
        "Coordinates\n1 0 0 0\n2 1 0 0\n3 0 1 0\nend coordinates\n"
        "Elements\n1 1 2 3\nend elements\n",
    )
    _write(
        tmp_path,
        "nm.post.res",
        "GiD Post Results File 1.2\n"
        'GaussPoints "gp" ElemType Triangle "dimension"\n'
        "Number Of Gauss Points: 1\n"
        "Natural Coordinates: Internal\n"
        "End GaussPoints\n"
        'Result "q" "a" 1 Scalar OnGaussPoints "gp"\n'
        "Values\n1 5\nEnd Values\n",
    )
    mesh = meshioplusplus.gid.read(str(tmp_path / "nm.post.msh"))
    assert "q" not in mesh.cell_data


@pytest.mark.parametrize("spelling", ["Linear", "Line"])
def test_both_spellings_of_the_1d_type_are_read(tmp_path, spelling):
    """gidpost emits "Linear"; CIMNE's current grammar names it "Line"."""
    path = _write(
        tmp_path,
        f"{spelling}.post.msh",
        f'MESH "l" dimension 3 ElemType {spelling} Nnode 2\n'
        "Coordinates\n1 0 0 0\n2 1 0 0\nend coordinates\n"
        "Elements\n1 1 2\nend elements\n",
    )
    mesh = meshioplusplus.gid.read(str(path))
    assert mesh.cells[0].type == "line"
    assert len(mesh.points) == 2


# ---------------------------------------------------------------------------
# Result types (Matrix, Complex, and the rest of GiD_ResultType).
#
# The declaration rides a field_data key because the registry's (path, mesh)
# writers cannot carry a side-channel struct -- see gid.hpp. Component counts
# come from gidpost's own _ResultTypeInfo; component ORDER is quoted from
# CIMNE's Customization Manual and is pinned literally below, because the two
# complex families use opposite conventions and neither can be inferred from
# the other.
# ---------------------------------------------------------------------------

from meshioplusplus import gid as _gid  # noqa: E402

RT = _gid.ResultType
PREFIX = _gid.RESULT_TYPE_PREFIX

# (type, a legal component count) for every one of the nine types.
ALL_RESULT_TYPES = [
    (RT.SCALAR, 1),
    (RT.VECTOR, 2),
    (RT.VECTOR, 3),
    (RT.VECTOR, 4),
    (RT.MATRIX, 3),
    (RT.MATRIX, 6),
    (RT.PLAIN_DEFORMATION_MATRIX, 4),
    (RT.MAIN_MATRIX, 12),
    (RT.LOCAL_AXES, 3),
    (RT.COMPLEX_SCALAR, 2),
    (RT.COMPLEX_VECTOR, 4),
    (RT.COMPLEX_VECTOR, 6),
    (RT.COMPLEX_MATRIX, 6),
    (RT.COMPLEX_MATRIX, 12),
]


def _typed_mesh(rtype, k, name="f"):
    """A 2-triangle mesh carrying one declared point_data array.

    Values encode their own (row, component) position, so a permutation
    anywhere in the write/read path is visible rather than merely plausible.
    """
    pts = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    vals = np.arange(4 * k, dtype=float).reshape(4, k)
    arr = vals if k > 1 else vals.reshape(4)
    return meshioplusplus.Mesh(
        pts,
        [("triangle", np.array([[0, 1, 2], [1, 3, 2]]))],
        point_data={name: arr},
        field_data={PREFIX + name: np.array([int(rtype)])},
    )


@pytest.mark.parametrize("rtype,k", ALL_RESULT_TYPES)
def test_declared_result_type_roundtrips(tmp_path, rtype, k):
    path = tmp_path / "t.post.msh"
    mesh = _typed_mesh(rtype, k)
    meshioplusplus.gid.write(str(path), mesh)
    back = meshioplusplus.gid.read(str(path))

    expected = np.arange(4 * k, dtype=float).reshape(4, k)
    got = np.asarray(back.point_data["f"])
    assert got.reshape(4, k) == pytest.approx(expected)

    # The declaration survives -- except where it carries no information,
    # which is exactly the shapes the writer would have inferred anyway.
    inferred = {(RT.SCALAR, 1), (RT.VECTOR, 2), (RT.VECTOR, 3)}
    key = PREFIX + "f"
    if (rtype, k) in inferred:
        assert key not in back.field_data
    else:
        assert int(np.asarray(back.field_data[key]).ravel()[0]) == int(rtype)


@pytest.mark.parametrize("rtype,k", ALL_RESULT_TYPES)
def test_declared_result_type_reaches_the_file_header(tmp_path, rtype, k):
    """The written .post.res must name the declared type, not a substitute."""
    path = tmp_path / "h.post.msh"
    meshioplusplus.gid.write(str(path), _typed_mesh(rtype, k))
    text = (tmp_path / "h.post.res").read_text()
    names = {
        RT.SCALAR: "Scalar",
        RT.VECTOR: "Vector",
        RT.MATRIX: "Matrix",
        RT.PLAIN_DEFORMATION_MATRIX: "PlainDeformationMatrix",
        RT.MAIN_MATRIX: "MainMatrix",
        RT.LOCAL_AXES: "LocalAxes",
        RT.COMPLEX_SCALAR: "ComplexScalar",
        RT.COMPLEX_VECTOR: "ComplexVector",
        RT.COMPLEX_MATRIX: "ComplexMatrix",
    }
    assert f'Result "f" "meshio++" 1 {names[rtype]} OnNodes' in text


def _res_values(path):
    """The numeric rows of the single Result block in a .post.res."""
    rows, inside = [], False
    for line in pathlib.Path(path).read_text().splitlines():
        if line.strip().lower() == "values":
            inside = True
            continue
        if line.strip().lower() == "end values":
            break
        if inside:
            rows.append([float(t) for t in line.split()])
    return rows


def test_complex_vector_is_interleaved_and_complex_matrix_is_blocked(tmp_path):
    """The one trap in this feature, pinned literally.

    CIMNE's manual gives ``ComplexVector:6`` as ``x_real, x_imag, y_real,
    y_imag, z_real, z_imag`` -- **interleaved** -- but ``ComplexMatrix:12`` as
    ``Sxx_real .. Sxz_real, Sxx_imag .. Sxz_imag`` -- **blocked**. Same family,
    opposite conventions, so neither may be inferred from the other.

    meshio++ stores values verbatim in GiD's order rather than reinterpreting
    them, having no canonical complex layout of its own. This test therefore
    pins that the bytes come out in *input* order for both, which is what makes
    "verbatim" a checkable claim instead of a hope: if either path grew a
    re-pack, one of these two would fail.
    """
    for rtype, k in ((RT.COMPLEX_VECTOR, 6), (RT.COMPLEX_MATRIX, 12)):
        path = tmp_path / f"c{k}.post.msh"
        meshioplusplus.gid.write(str(path), _typed_mesh(rtype, k))
        rows = _res_values(tmp_path / f"c{k}.post.res")
        assert len(rows) == 4
        for r, row in enumerate(rows):
            assert row[0] == r + 1  # the node id
            assert row[1:] == pytest.approx(
                list(np.arange(r * k, (r + 1) * k, dtype=float))
            )


def test_matrix_needs_no_permutation(tmp_path):
    """GiD's ``Matrix:6`` is ``Sxx Syy Szz Sxy Syz Sxz`` -- already meshio/VTK's
    symmetric-tensor order, so a stress tensor is written straight through.

    This corrects a claim that stood in the docs until now (that the two orders
    differ), which was the stated reason for refusing ``GiD_Matrix`` at all.
    """
    path = tmp_path / "m.post.msh"
    meshioplusplus.gid.write(str(path), _typed_mesh(RT.MATRIX, 6))
    rows = _res_values(tmp_path / "m.post.res")
    assert rows[0][1:] == pytest.approx([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])


@pytest.mark.parametrize(
    "rtype,k,needle",
    [
        (RT.MATRIX, 5, "accepts 3, 6"),
        (RT.COMPLEX_SCALAR, 3, "accepts 2"),
        (RT.MAIN_MATRIX, 6, "accepts 12"),
    ],
)
def test_illegal_component_count_errors_by_name(tmp_path, rtype, k, needle):
    """Never a silent fallback to splitting -- write_options.hpp's rule."""
    with pytest.raises(Exception) as exc:
        meshioplusplus.gid.write(str(tmp_path / "bad.post.msh"), _typed_mesh(rtype, k))
    msg = str(exc.value)
    assert "'f'" in msg and needle in msg


def test_out_of_range_result_type_errors_by_name(tmp_path):
    mesh = _typed_mesh(RT.MATRIX, 6)
    mesh.field_data[PREFIX + "f"] = np.array([42])
    with pytest.raises(Exception) as exc:
        meshioplusplus.gid.write(str(tmp_path / "bad.post.msh"), mesh)
    assert "42" in str(exc.value) and "0..8" in str(exc.value)


def test_undeclared_arrays_are_unaffected(tmp_path):
    """A mesh with no declaration keeps the historical inference exactly.

    In particular a 6-component array still splits into six named scalars
    rather than silently becoming a Matrix -- ``(n, 6)`` is genuinely
    ambiguous, so inferring would pick one meaning without being asked.
    """
    pts = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    mesh = meshioplusplus.Mesh(
        pts,
        [("triangle", np.array([[0, 1, 2], [1, 3, 2]]))],
        point_data={"s6": np.arange(24, dtype=float).reshape(4, 6)},
    )
    path = tmp_path / "u.post.msh"
    meshioplusplus.gid.write(str(path), mesh)
    text = (tmp_path / "u.post.res").read_text()
    assert "Matrix" not in text
    for c in range(1, 7):
        assert f'Result "s6_{c}"' in text

    back = meshioplusplus.gid.read(str(path))
    assert not any(key.startswith(PREFIX) for key in back.field_data)
