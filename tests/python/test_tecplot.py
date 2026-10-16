import pathlib
from copy import deepcopy

import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        # helpers.empty_mesh,
        helpers.tri_mesh,
        helpers.quad_mesh,
        helpers.tet_mesh,
        helpers.hex_mesh,
    ],
)
def test(mesh, tmp_path):
    helpers.write_read(
        tmp_path,
        meshioplusplus.tecplot.write,
        meshioplusplus.tecplot.read,
        mesh,
        1.0e-15,
    )


@pytest.mark.parametrize(
    "filename", ["quad_zone_comma.tec", "quad_zone_space.tec", "quad_zone_multivar.tec"]
)
def test_comma_space(filename, tmp_path):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "tecplot" / filename
    mesh = meshioplusplus.read(filename)

    helpers.write_read(
        tmp_path,
        meshioplusplus.tecplot.write,
        meshioplusplus.tecplot.read,
        mesh,
        1.0e-15,
    )


def test_varlocation(tmp_path):
    # Test that VARLOCATION is correctly written and read depending on the
    # number of point and cell data.
    writer = meshioplusplus.tecplot.write
    reader = meshioplusplus.tecplot.read
    mesh = deepcopy(helpers.tri_mesh)
    num_points = len(mesh.points)
    num_cells = sum(len(c.data) for c in mesh.cells)

    # Add point data: no VARLOCATION
    mesh.point_data["one"] = np.ones(num_points)
    helpers.write_read(tmp_path, writer, reader, mesh, 1.0e-15)

    # Add cell data: VARLOCATION = ([5] = CELLCENTERED)
    mesh.cell_data["two"] = [np.ones(num_cells) * 2.0]
    helpers.write_read(tmp_path, writer, reader, mesh, 1.0e-15)

    # Add point data: VARLOCATION = ([6] = CELLCENTERED)
    mesh.point_data["three"] = np.ones(num_points) * 3.0
    helpers.write_read(tmp_path, writer, reader, mesh, 1.0e-15)

    # Add cell data: VARLOCATION = ([6-7] = CELLCENTERED)
    mesh.cell_data["four"] = [np.ones(num_cells) * 4.0]
    helpers.write_read(tmp_path, writer, reader, mesh, 1.0e-15)

    # Add point data: VARLOCATION = ([7-8] = CELLCENTERED)
    mesh.point_data["five"] = np.ones(num_points) * 5.0
    helpers.write_read(tmp_path, writer, reader, mesh, 1.0e-15)

    # Add cell data: VARLOCATION = ([7-9] = CELLCENTERED)
    mesh.cell_data["six"] = [np.ones(num_cells) * 6.0]
    helpers.write_read(tmp_path, writer, reader, mesh, 1.0e-15)


# --- malformed-input / error-path coverage (Python reference reader) ---
from meshioplusplus.tecplot._tecplot import read as _tecplot_py_read  # noqa: E402


def test_tecplot_missing_x_variable_raises(tmp_path):
    # VARIABLES list without the mandatory coordinate 'X'.
    p = tmp_path / "bad.dat"
    p.write_text('VARIABLES = "A" "B"\nZONE N=1 E=1 F=FEPOINT ET=TRIANGLE\n')
    with pytest.raises(meshioplusplus.ReadError):
        _tecplot_py_read(str(p))


def test_read_metadata_reports_both_zones_of_a_transient_file(tmp_path):
    """roadmap §1 tier B1: read_tecplot_metadata scans every ZONE header
    (SOLUTIONTIME/STRANDID) without decoding any data body, and
    ``time_step`` picks the one zone that SOLUTIONTIME/STRANDID timeline
    names -- two FEBLOCK zones sharing STRANDID=1, hand-written text since
    Tecplot ASCII needs no external tool."""
    p = tmp_path / "transient.dat"
    p.write_text(
        'VARIABLES = "X" "Y" "u"\n'
        "ZONE N=3 E=1 DATAPACKING=BLOCK ZONETYPE=FETRIANGLE SOLUTIONTIME=0.0 STRANDID=1\n"
        "0.0 1.0 0.0\n"
        "0.0 0.0 1.0\n"
        "10.0 20.0 30.0\n"
        "1 2 3\n"
        "ZONE N=3 E=1 DATAPACKING=BLOCK ZONETYPE=FETRIANGLE SOLUTIONTIME=2.5 STRANDID=1\n"
        "0.0 1.0 0.0\n"
        "0.0 0.0 1.0\n"
        "11.0 21.0 31.0\n"
        "1 2 3\n"
    )

    meta = meshioplusplus.read_metadata(p, "tecplot")
    assert meta["fell_back_to_full_read"] is False
    assert meta["format"] == "tecplot"
    assert meta["time_values"] == [0.0, 2.5]
    assert meta["num_points"] == 3

    mesh0 = meshioplusplus.tecplot.read(p, time_step=0)
    assert mesh0.point_data["u"][0] == 10.0
    mesh1 = meshioplusplus.tecplot.read(p, time_step=1)
    assert mesh1.point_data["u"][0] == 11.0
    mesh_last = meshioplusplus.tecplot.read(p, time_step=-1)
    assert mesh_last.point_data["u"][0] == 11.0

    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.tecplot.read(p, time_step=5)


_FUZZ_TECPLOT = (
    pathlib.Path(__file__).parent.parent / "fuzz" / "regressions" / "tecplot"
)


@pytest.mark.parametrize(
    "read", [meshioplusplus.tecplot.read, _tecplot_py_read], ids=["public", "python"]
)
@pytest.mark.parametrize(
    "name",
    ["face-node-count-overflows-sum", "plt-passive-count-exceeds-file"],
)
def test_counts_beyond_the_file_are_read_errors(read, name):
    # The fuzzer's crashes: a face count that overflowed the node list's running
    # sum (UBSan), and a passive variable sized from a zone's node count the
    # file cannot hold (ASan). Both engines refuse them as ReadError.
    with pytest.raises(meshioplusplus.ReadError):
        read(str(_FUZZ_TECPLOT / name))


@pytest.mark.parametrize(
    "read", [meshioplusplus.tecplot.read, _tecplot_py_read], ids=["public", "python"]
)
def test_connectivity_ids_outside_the_zone_are_read_errors(read, tmp_path):
    # The fuzzer's crash: a connectivity id near INT64_MAX, saturated by the C++
    # lenient parse, whose point offset overflowed (UBSan). An id past the zone's
    # node count is refused the same way, and the zone's last node is accepted.
    with pytest.raises(meshioplusplus.ReadError):
        read(str(_FUZZ_TECPLOT / "ascii-connectivity-id-overflows-offset"))
    header = (
        'VARIABLES = "X" "Y"\n'
        "ZONE N=3 E=1 DATAPACKING=POINT ZONETYPE=FETRIANGLE\n"
        "0.0 0.0\n1.0 0.0\n0.0 1.0\n"
    )
    last = tmp_path / "last-node.dat"
    last.write_text(header + "1 2 3\n")
    read(str(last))
    past = tmp_path / "past-node.dat"
    past.write_text(header + "1 2 4\n")
    with pytest.raises(meshioplusplus.ReadError):
        read(str(past))


def _plt_meshes():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1.0]])
    return {
        "fe": meshioplusplus.Mesh(
            pts,
            [("tetra", [[0, 1, 2, 3], [1, 2, 3, 4]]), ("triangle", [[0, 1, 2]])],
            point_data={"T": np.arange(5.0), "U": np.arange(15.0).reshape(5, 3)},
            cell_data={"p": [np.array([1.5, 2.5]), np.array([9.0])]},
        ),
        "polygon": meshioplusplus.Mesh(
            np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0.0]]),
            [meshioplusplus.CellBlock("polygon3", [[0, 1, 2], [0, 2, 3]])],
        ),
        "polyhedron": meshioplusplus.Mesh(
            np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1.0]]),
            [
                meshioplusplus.CellBlock(
                    "polyhedron4", [[[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]]]
                )
            ],
        ),
        "passive": meshioplusplus.Mesh(
            pts,
            [("tetra", [[0, 1, 2, 3]]), ("tetra", [[1, 2, 3, 4]])],
            cell_data={"c": [np.array([1.0]), np.array([2.0])]},
        ),
    }


@pytest.mark.parametrize("name", ["fe", "polygon", "polyhedron", "passive"])
def test_binary_plt_writer_is_byte_identical_across_engines_and_reads_back(
    name, tmp_path
):
    """Roadmap 3.3.2: a `.plt` name writes the binary #!TDV112 file."""
    from meshioplusplus import _core
    from meshioplusplus.tecplot import _tecplot

    mesh = _plt_meshes()[name]
    a, b = tmp_path / "core.plt", tmp_path / "twin.plt"
    _core.tecplot_write(str(a), mesh)
    _tecplot.write(b, mesh)
    assert a.read_bytes() == b.read_bytes()
    assert a.read_bytes()[:8] == b"#!TDV112"

    ascii_path = tmp_path / "same.dat"
    _tecplot.write(ascii_path, mesh)
    plt, dat = meshioplusplus.read(a), meshioplusplus.read(ascii_path)
    assert [c.type for c in plt.cells] == [c.type for c in dat.cells]
    np.testing.assert_allclose(plt.points, dat.points)
    for x, y in zip(plt.cells, dat.cells):
        assert len(x.data) == len(y.data)
        for r, q in zip(x.data, y.data):
            np.testing.assert_array_equal(
                (
                    np.concatenate([np.ravel(f) for f in np.atleast_1d(r)])
                    if x.type.startswith("polyhedron")
                    else r
                ),
                (
                    np.concatenate([np.ravel(f) for f in np.atleast_1d(q)])
                    if y.type.startswith("polyhedron")
                    else q
                ),
            )
    for key in dat.point_data:
        np.testing.assert_allclose(plt.point_data[key], dat.point_data[key])
    for key in dat.cell_data:
        for x, y in zip(plt.cell_data[key], dat.cell_data[key]):
            np.testing.assert_allclose(x, y)
