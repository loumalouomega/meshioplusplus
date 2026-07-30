import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.exodus import _exodus as py_exodus

from . import helpers
from .exodus_fixture import (
    BLOCK_IDS,
    BLOCK_NAMES,
    NODE_SET_IDS,
    NODE_SETS,
    QA_RECORD,
    SIDE_SET_EXPECTED,
    SIDE_SET_ID,
    SIDE_SET_NAME,
    SPHERE_ATTRIBUTE_NAME,
    SPHERE_BLOCK_NODES,
    SPHERE_POINTS,
    SPHERE_RADII,
    TIME_VALUES,
    temperature_at,
    write_fixture,
    write_sphere_fixture,
)

netCDF4 = pytest.importorskip("netCDF4")

_HAS_CPP = getattr(_core, "__has_netcdf__", False)


@pytest.fixture
def seacas_file(tmp_path):
    """A SEACAS/Cubit-shaped file: qa_records, 2 same-type blocks, sets, 3 steps."""
    return write_fixture(tmp_path / "seacas.exo")


def _readers():
    """The C++ core and the pure-Python twin, which must agree.

    Both are exercised directly rather than through the shim, because the shim
    silently falls back on any core failure -- which is precisely what hid the
    `qa_records` throw for as long as it did.
    """
    readers = [pytest.param(py_exodus.read, id="python")]
    if _HAS_CPP:
        readers.append(
            pytest.param(
                lambda p, time_step=0: _core.exodus_read(str(p), time_step=time_step),
                id="cpp",
            )
        )
    return readers


@pytest.fixture
def sphere_file(tmp_path):
    """A peridynamics-shaped file: 2-D, NUL-padded SPHERE elem_type, attributes."""
    return write_sphere_fixture(tmp_path / "spheres.exo")


def _writers():
    """The C++ core writer and the pure-Python twin, which must agree."""
    writers = [pytest.param(py_exodus.write, id="python")]
    if _HAS_CPP:
        writers.append(
            pytest.param(lambda p, m: _core.exodus_write(str(p), m), id="cpp")
        )
    return writers


test_set = [
    helpers.empty_mesh,
    helpers.tri_mesh,
    helpers.tri_mesh_2d,
    helpers.triangle6_mesh,
    helpers.quad_mesh,
    helpers.quad8_mesh,
    helpers.tri_quad_mesh,
    helpers.tet_mesh,
    helpers.tet10_mesh,
    helpers.hex_mesh,
    helpers.hex20_mesh,
    helpers.add_point_data(helpers.tri_mesh, 1),
    helpers.add_point_data(helpers.tri_mesh, 2),
    helpers.add_point_data(helpers.tri_mesh, 3),
    helpers.add_point_sets(helpers.tri_mesh),
    helpers.add_point_sets(helpers.tet_mesh),
]


@pytest.mark.parametrize("mesh", test_set)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.exodus.write, meshioplusplus.exodus.read, mesh, 1.0e-15
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.e")
    # With additional, insignificant suffix:
    helpers.generic_io(tmp_path / "test.0.e")


# ---------------------------------------------------------------------------- #
# Real-world (SEACAS/Cubit-shaped) files. See tests/python/exodus_fixture.py --  #
# meshio++'s own writer cannot emit qa_records/eb_names/side sets, so a          #
# round-trip test could not have caught any of this.                            #
# ---------------------------------------------------------------------------- #


@pytest.mark.parametrize("reader", _readers())
def test_qa_records_do_not_fail_the_read(reader, seacas_file):
    """The regression: this used to throw, making the format unusable in WASM."""
    mesh = reader(seacas_file)
    assert len(mesh.points) == 12
    assert [cb.type for cb in mesh.cells] == ["hexahedron", "hexahedron"]
    # Provenance is preserved, not silently dropped -- it rides the ExodusInfo
    # side channel on the C++ path (NDArray has no string dtype).
    assert list(mesh.info) == QA_RECORD


@pytest.mark.parametrize("reader", _readers())
def test_element_blocks_become_named_cell_regions(reader, seacas_file):
    mesh = reader(seacas_file)
    cell_regions = {r.name: r for r in mesh.regions if r.kind == "cell"}
    assert sorted(cell_regions) == sorted(BLOCK_NAMES)
    # Both blocks are HEX8: they must NOT collapse together. Each holds exactly
    # its own cell, and carries its own Exodus block id as the tag.
    for k, name in enumerate(BLOCK_NAMES):
        region = cell_regions[name]
        assert region.entries.tolist() == [k], f"{name} must hold only its own cell"
        assert region.tag == BLOCK_IDS[k]
        assert region.dim == 3


@pytest.mark.parametrize("reader", _readers())
def test_node_sets_become_point_regions(reader, seacas_file):
    mesh = reader(seacas_file)
    point_regions = {r.name: r for r in mesh.regions if r.kind == "point"}
    assert sorted(point_regions) == sorted(NODE_SETS)
    for k, (name, nodes) in enumerate(NODE_SETS.items()):
        # Entries are canonicalized (sorted) by AddRegion, so compare as sets.
        assert sorted(point_regions[name].entries.tolist()) == sorted(nodes)
        assert point_regions[name].tag == NODE_SET_IDS[k]
    # point_sets is a compat view over the same regions.
    assert sorted(mesh.point_sets) == sorted(NODE_SETS)


@pytest.mark.parametrize("reader", _readers())
def test_side_sets_become_side_regions_with_remapped_facets(reader, seacas_file):
    mesh = reader(seacas_file)
    side_regions = [r for r in mesh.regions if r.kind == "side"]
    assert len(side_regions) == 1
    region = side_regions[0]
    assert region.name == SIDE_SET_NAME
    assert region.tag == SIDE_SET_ID
    # The two entries use DIFFERENT Exodus side numbers (4 and 2), so an
    # identity side->facet mapping would give a visibly wrong answer here.
    assert [tuple(row) for row in region.entries.tolist()] == SIDE_SET_EXPECTED


def test_side_facets_match_cell_faces(seacas_file):
    """Pin the exo_face_index table against the actual face topology.

    The mapping is only correct if the facet it names really is the face the
    Exodus side number described -- assert that by node set, not by trusting the
    transcription.
    """
    from meshioplusplus._regions import block_bases
    from meshioplusplus._skin import _CELL_FACES

    mesh = py_exodus.read(seacas_file)
    region = next(r for r in mesh.regions if r.kind == "side")
    bases = block_bases(mesh.cells)
    for cell, facet in region.entries.tolist():
        # Side entries are GLOBAL block-major cell indices; resolve through the
        # same bases the reader used rather than assuming one cell per block.
        b = int(np.searchsorted(bases, cell, side="right")) - 1
        row = mesh.cells[b].data[cell - bases[b]]
        # _CELL_FACES rows are (face_type, num_corners, local_node_indices).
        _, _, local = _CELL_FACES[mesh.cells[b].type][facet]
        face_nodes = {int(row[i]) for i in local}
        # Both side-set facets were chosen to be end faces, i.e. exactly the
        # nodes of the matching node set.
        expected = set(NODE_SETS["left"] if cell == 0 else NODE_SETS["right"])
        assert face_nodes == expected, f"cell {cell} facet {facet} is the wrong face"


@pytest.mark.parametrize("reader", _readers())
@pytest.mark.parametrize("step", [0, 1, 2, -1, -3])
def test_time_step_selects_the_requested_step(reader, step, seacas_file):
    mesh = reader(seacas_file, time_step=step)
    expected = temperature_at(step % len(TIME_VALUES))
    assert np.allclose(mesh.point_data["temperature"], expected)


@pytest.mark.parametrize("reader", _readers())
@pytest.mark.parametrize("step", [3, 7, -4])
def test_out_of_range_time_step_is_an_error_naming_the_count(reader, step, seacas_file):
    """Never a silent clamp: returning step 0 for a bad request is the failure
    mode this option exists to remove."""
    with pytest.raises(Exception, match=r"3 steps"):
        reader(seacas_file, time_step=step)


def test_cpp_matches_python(seacas_file):
    """The two paths must be interchangeable, since the shim swaps them freely."""
    if not _HAS_CPP:
        pytest.skip("core built without netCDF")
    for step in (0, 1, 2):
        cpp = _core.exodus_read(str(seacas_file), time_step=step)
        py = py_exodus.read(seacas_file, time_step=step)
        assert np.allclose(cpp.points, py.points)
        assert [cb.type for cb in cpp.cells] == [cb.type for cb in py.cells]
        for a, b in zip(cpp.cells, py.cells):
            assert np.array_equal(a.data, b.data)
        assert np.allclose(cpp.point_data["temperature"], py.point_data["temperature"])
        assert list(cpp.info) == list(py.info)

        # Regions are compared as a set: AddRegion canonicalizes entries, but
        # the two readers build the list in a different order (blocks first vs
        # sets first), which is not part of the contract.
        def key(regions):
            return sorted(
                (
                    r.name,
                    r.kind,
                    r.dim,
                    r.tag,
                    tuple(map(tuple, np.atleast_2d(r.entries))),
                )
                for r in regions
            )

        assert key(cpp.regions) == key(py.regions)


def test_generic_read_and_metadata_report_time_steps(seacas_file):
    mesh = meshioplusplus.read(seacas_file, time_step=-1)
    assert np.allclose(mesh.point_data["temperature"], temperature_at(2))

    meta = meshioplusplus.read_metadata(seacas_file)
    assert meta["format"] == "exodus"
    assert list(meta["time_values"]) == TIME_VALUES


def test_reader_supports_options():
    """Exodus is now an options-aware reader; before this change it was not."""
    if not _HAS_CPP:
        pytest.skip("core built without netCDF")
    assert _core.reader_supports_options("exodus")


# ---------------------------------------------------------------------------- #
# Peridynamics-shaped files: SPHERE elements, 2-D coordinates and per-element    #
# attributes. See tests/python/exodus_fixture.py's SPHERE_FIXTURE_FACTS.         #
# ---------------------------------------------------------------------------- #

ATTRIBUTE_KEY = "exodus:attr:" + SPHERE_ATTRIBUTE_NAME


@pytest.mark.parametrize("reader", _readers())
def test_nul_terminated_elem_type_is_read(reader, sphere_file):
    """The regression: an `elem_type` of "SPHERE\\0" must still be SPHERE.

    NetCDF.jl counts the terminating NUL as part of the attribute, so the C++
    reader compared 7 characters against its type table and failed with
    "unknown element type SPHERE" -- the NUL invisible in the message, because
    `what()` is a `const char*`. netCDF4-python strips it, so only the C++ path
    (and therefore only WASM, which has no Python fallback) ever failed.
    """
    mesh = reader(sphere_file)
    assert [cb.type for cb in mesh.cells] == ["vertex", "vertex"]
    for block, nodes in zip(mesh.cells, SPHERE_BLOCK_NODES):
        assert block.data.reshape(-1).tolist() == nodes


@pytest.mark.parametrize("reader", _readers())
def test_two_dimensional_points_pad_to_three_components(reader, sphere_file):
    """num_dim=2 with coordx/coordy and no coordz -- z fills with zeros."""
    mesh = reader(sphere_file)
    assert mesh.points.shape == (len(SPHERE_POINTS), 3)
    assert np.allclose(mesh.points[:, :2], SPHERE_POINTS)
    assert np.all(mesh.points[:, 2] == 0.0)


@pytest.mark.parametrize("reader", _readers())
def test_element_attributes_become_prefixed_cell_data(reader, sphere_file):
    """A SPHERE's radius lives in `attrib{k}`; it must reach the user as data.

    The prefix keeps it apart from a same-named element *variable*, which is a
    genuinely different concept (attributes are constant in time).
    """
    mesh = reader(sphere_file)
    assert ATTRIBUTE_KEY in mesh.cell_data
    values = mesh.cell_data[ATTRIBUTE_KEY]
    assert len(values) == len(mesh.cells), "one array per cell block"
    assert values[0].tolist() == SPHERE_RADII
    # Block 2 carries no attribute at all -- NaN, since a cell_data array must
    # cover every block and there is no value to report.
    assert np.all(np.isnan(values[1]))


@pytest.mark.parametrize("writer", _writers())
@pytest.mark.parametrize("reader", _readers())
def test_element_attributes_round_trip(reader, writer, sphere_file, tmp_path):
    """Both writers must put attributes back where both readers find them."""
    mesh = py_exodus.read(sphere_file)
    out = tmp_path / "out.exo"
    writer(out, mesh)

    back = reader(out)
    values = back.cell_data[ATTRIBUTE_KEY]
    assert values[0].tolist() == SPHERE_RADII
    # The all-NaN block must NOT gain an attribute of its own: that NaN means
    # "this block never had one", so writing it back would be inventing data.
    assert np.all(np.isnan(values[1]))
    with netCDF4.Dataset(out) as nc:
        assert "num_att_in_blk1" in nc.dimensions
        assert "num_att_in_blk2" not in nc.dimensions


@pytest.mark.parametrize("writer", _writers())
def test_non_scalar_element_attribute_is_an_error(writer, sphere_file, tmp_path):
    """An Exodus attribute is one value per element -- a vector cannot go there."""
    mesh = py_exodus.read(sphere_file)
    mesh.cell_data[ATTRIBUTE_KEY] = [np.zeros((len(cb.data), 2)) for cb in mesh.cells]
    with pytest.raises(Exception, match="scalar"):
        writer(tmp_path / "bad.exo", mesh)


def test_cpp_matches_python_on_spheres(sphere_file):
    """The shim swaps the two paths freely, so they must agree on this file too."""
    if not _HAS_CPP:
        pytest.skip("core built without netCDF")
    cpp = _core.exodus_read(str(sphere_file))
    py = py_exodus.read(sphere_file)
    assert np.allclose(cpp.points, py.points)
    assert [cb.type for cb in cpp.cells] == [cb.type for cb in py.cells]
    for a, b in zip(cpp.cells, py.cells):
        assert np.array_equal(a.data, b.data)
    assert sorted(cpp.cell_data) == sorted(py.cell_data)
    for name, arrays in cpp.cell_data.items():
        for a, b in zip(arrays, py.cell_data[name]):
            assert np.array_equal(a, b, equal_nan=True)


def test_generic_read_handles_spheres(sphere_file):
    """The public entry point, i.e. what the CLI and every consumer goes through."""
    mesh = meshioplusplus.read(sphere_file)
    assert [cb.type for cb in mesh.cells] == ["vertex", "vertex"]
    assert mesh.cell_data[ATTRIBUTE_KEY][0].tolist() == SPHERE_RADII


# ---------------------------------------------------------------------------- #
# The real PeriLab output from VSCode-MDPA-Preview#63. Everything above          #
# constructs the failing shape; this is the file that actually failed. See       #
# tests/python/meshes/exodus/README.md for provenance and licensing.             #
# ---------------------------------------------------------------------------- #

PERILAB_FILE = (
    pathlib.Path(__file__).resolve().parent
    / "meshes"
    / "exodus"
    / "DCBmodel_PD_solid.e"
)

#: 4 SPHERE element blocks, in file order.
PERILAB_BLOCK_SIZES = [216, 216, 36, 36]
PERILAB_POINT_DATA = [
    "Cauchy Stressxx",
    "Cauchy Stressxy",
    "Cauchy Stressyx",
    "Cauchy Stressyy",
    "Damage",
    "Displacementsx",
    "Displacementsy",
    "Forcesx",
    "Forcesy",
]
PERILAB_NUM_STEPS = 10


@pytest.fixture
def perilab_file():
    """The reference file, having checked it is not an unfetched LFS pointer.

    A hard failure with a readable message rather than a skip: the two CI jobs
    that run this suite (`test`, `coverage`) both check out with `lfs: true`, so
    a pointer here means the checkout is broken, and skipping would let that
    look green. The netCDF error you get from reading a 130-byte pointer says
    nothing about the cause.
    """
    with open(PERILAB_FILE, "rb") as f:
        if f.read(24).startswith(b"version https://git-lfs"):
            pytest.fail(f"{PERILAB_FILE.name} is a Git-LFS pointer; run `git lfs pull`")
    return PERILAB_FILE


def test_the_reference_file_still_carries_a_nul_terminated_elem_type(perilab_file):
    """Guard the fixture's defining property, not just the code that reads it.

    Everything below passes just as well against a file whose `elem_type` is a
    plain 6-character "SPHERE" — so if this file were ever re-fetched from an
    upstream that had switched netCDF writers, the regression suite would go on
    passing while testing nothing. The classic-format encoding of an attribute
    is `[nc_type][nelems][value]`, so assert the recorded length is 7, i.e. that
    the terminating NUL is genuinely counted in.
    """
    import struct

    blob = perilab_file.read_bytes()
    assert not blob.startswith(
        b"version https://git-lfs"
    ), "the fixture is an unfetched Git-LFS pointer; run `git lfs pull`"
    padded = struct.pack(">II", 2, 7) + b"SPHERE\x00"  # 2 = NC_CHAR
    plain = struct.pack(">II", 2, 6) + b"SPHERE"
    assert blob.count(padded) == len(PERILAB_BLOCK_SIZES)
    assert blob.count(plain) == 0, "no block may have the plain 6-character form"


@pytest.mark.parametrize("reader", _readers())
def test_perilab_reference_file(reader, perilab_file):
    """The end-to-end regression: this read used to fail on the C++ path.

    `Exodus: unknown element type SPHERE` — on a type both tables had mapped to
    `vertex` all along. Only the C++ reader compared the NUL too, so this was
    invisible from Python and fatal in WASM, which has no fallback.
    """
    mesh = reader(perilab_file)

    assert [cb.type for cb in mesh.cells] == ["vertex"] * len(PERILAB_BLOCK_SIZES)
    assert [len(cb.data) for cb in mesh.cells] == PERILAB_BLOCK_SIZES
    # num_dim = 2 (coordx/coordy, no coordz): points pad out to three columns.
    assert mesh.points.shape == (504, 3)
    assert np.all(mesh.points[:, 2] == 0.0)

    assert sorted(mesh.point_data) == PERILAB_POINT_DATA
    # None of these recombine: `categorize()` only pairs names ending in an
    # UPPERCASE X/Y/Z, and PeriLab writes "Displacementsx".
    assert mesh.point_data["Damage"].shape == (504,)

    cell_regions = {r.name: r for r in mesh.regions if r.kind == "cell"}
    assert sorted(cell_regions) == ["block_1", "block_2", "block_3", "block_4"]
    for name, size in zip(sorted(cell_regions), PERILAB_BLOCK_SIZES):
        # A SPHERE is topologically a point, so its blocks are dimension 0.
        assert len(cell_regions[name].entries) == size
        assert cell_regions[name].dim == 0
    point_regions = [r for r in mesh.regions if r.kind == "point"]
    assert len(point_regions) == 2
    assert all(len(r.entries) == 36 for r in point_regions)

    # qa_records/info_records: real provenance, the thing that used to throw.
    assert any("PeriLab" in line for line in mesh.info)


@pytest.mark.parametrize("reader", _readers())
def test_perilab_reference_file_time_steps(reader, perilab_file):
    """A genuine 10-step run: the damage field is what makes a wrong step visible."""
    first = reader(perilab_file, time_step=0)
    last = reader(perilab_file, time_step=-1)
    # Undamaged at t=0 and cracked at the end -- so a reader silently pinned to
    # step 0 (the pre-v8.6.0 behaviour) would fail this, not merely differ.
    assert np.nanmax(first.point_data["Damage"]) == 0.0
    assert np.nanmax(last.point_data["Damage"]) == pytest.approx(0.48148148, abs=1e-7)

    meta = meshioplusplus.read_metadata(perilab_file)
    assert meta["format"] == "exodus"
    assert len(meta["time_values"]) == PERILAB_NUM_STEPS
    assert meta["time_values"][0] == 0.0


def test_perilab_reference_file_cpp_matches_python(perilab_file):
    """The shim swaps the paths freely, so they must agree on the real file too."""
    if not _HAS_CPP:
        pytest.skip("core built without netCDF")
    cpp = _core.exodus_read(str(perilab_file))
    py = py_exodus.read(perilab_file)
    assert np.array_equal(cpp.points, py.points)
    assert [cb.type for cb in cpp.cells] == [cb.type for cb in py.cells]
    for a, b in zip(cpp.cells, py.cells):
        assert np.array_equal(a.data, b.data)
    assert sorted(cpp.point_data) == sorted(py.point_data)
    for name in cpp.point_data:
        assert np.array_equal(cpp.point_data[name], py.point_data[name])
    assert list(cpp.info) == list(py.info)

    def key(regions):
        return sorted(
            (r.name, r.kind, r.dim, r.tag, tuple(map(tuple, np.atleast_2d(r.entries))))
            for r in regions
        )

    assert key(cpp.regions) == key(py.regions)


# ---------------------------------------------------------------------------- #
# v9.9.0: the three things this writer used to drop -- ordinary cell_data,       #
# block names, and the time value. Both writers must agree, so every case is     #
# parametrized over both and read back with both readers.                       #
# ---------------------------------------------------------------------------- #


def _two_block_mesh():
    """Two same-type blocks, so a per-block element variable has real structure."""
    return meshioplusplus.Mesh(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [2.0, 0.0, 0.0],
            ]
        ),
        [
            ("triangle", np.array([[0, 1, 2]])),
            ("triangle", np.array([[0, 2, 3], [1, 4, 2]])),
        ],
    )


@pytest.mark.parametrize("writer", _writers())
@pytest.mark.parametrize("reader", _readers())
def test_cell_data_round_trips_as_element_variables(reader, writer, tmp_path):
    """The regression: this writer emitted no `vals_elem_var` at all."""
    mesh = _two_block_mesh()
    mesh.cell_data["material"] = [
        np.array([7.0]),
        np.array([8.0, 9.0]),
    ]
    out = tmp_path / "ev.exo"
    writer(out, mesh)

    back = reader(out)
    assert "material" in back.cell_data, "ordinary cell_data was dropped"
    values = back.cell_data["material"]
    assert len(values) == len(mesh.cells)
    assert values[0].ravel().tolist() == [7.0]
    assert values[1].ravel().tolist() == [8.0, 9.0]


@pytest.mark.parametrize("writer", _writers())
@pytest.mark.parametrize("reader", _readers())
def test_multi_component_cell_data_round_trips(reader, writer, tmp_path):
    """A vector cell field: trailing dims become extra netCDF dimensions."""
    mesh = _two_block_mesh()
    mesh.cell_data["stress"] = [
        np.array([[0.0, 1.0, 2.0]]),
        np.array([[3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]),
    ]
    out = tmp_path / "vec.exo"
    writer(out, mesh)

    back = reader(out)
    values = back.cell_data["stress"]
    assert values[0].shape == (1, 3)
    assert values[1].shape == (2, 3)
    assert values[0].tolist() == [[0.0, 1.0, 2.0]]
    assert values[1].tolist() == [[3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]


@pytest.mark.parametrize("writer", _writers())
@pytest.mark.parametrize("reader", _readers())
def test_block_names_round_trip_through_cell_regions(reader, writer, tmp_path):
    """A block name used to come back as the reader's synthetic "Block N"."""
    from meshioplusplus._regions import Region

    mesh = _two_block_mesh()
    mesh.regions = [
        Region("inner", "cell", np.array([0]), dim=2),
        Region("outer", "cell", np.array([1, 2]), dim=2),
    ]
    out = tmp_path / "named.exo"
    writer(out, mesh)

    back = reader(out)
    names = {r.name for r in back.regions if r.kind == "cell"}
    assert "inner" in names and "outer" in names, f"got {names}"


@pytest.mark.parametrize("writer", _writers())
@pytest.mark.parametrize("reader", _readers())
def test_time_value_comes_from_field_data(reader, writer, tmp_path):
    """`time_whole` was hard-coded 0, so a labelled frame lost its time."""
    mesh = _two_block_mesh()
    mesh.field_data["exodus:time"] = np.array([2.5])
    out = tmp_path / "t.exo"
    writer(out, mesh)

    with netCDF4.Dataset(out) as nc:
        assert float(nc.variables["time_whole"][0]) == pytest.approx(2.5)

    # And it comes back as `exodus:time`, so a frame survives a round trip.
    back = reader(out)
    assert back.field_data["exodus:time"].ravel()[0] == pytest.approx(2.5)


def test_a_region_less_mesh_still_writes_no_eb_names(tmp_path):
    """`eb_names` is written only when a block actually has a name."""
    py_exodus.write(tmp_path / "plain.exo", _two_block_mesh())
    with netCDF4.Dataset(tmp_path / "plain.exo") as nc:
        assert "eb_names" not in nc.variables
