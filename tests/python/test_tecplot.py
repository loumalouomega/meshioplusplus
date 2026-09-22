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
