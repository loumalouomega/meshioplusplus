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
    TIME_VALUES,
    temperature_at,
    write_fixture,
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
