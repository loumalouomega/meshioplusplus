import pathlib

import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        # WKT/TIN is triangle-only.
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.tri_mesh_one_cell,
    ],
)
def test_wkt(mesh, tmp_path):
    helpers.write_read(
        tmp_path,
        meshioplusplus.wkt.write,
        meshioplusplus.wkt.read,
        mesh,
        1.0e-12,
        ".wkt",
    )


@pytest.mark.parametrize(
    "filename, ref_sum, ref_num_cells",
    [("simple.wkt", 4, 2), ("whitespaced.wkt", 3.2, 2)],
)
def test_reference_file(filename, ref_sum, ref_num_cells):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "wkt" / filename

    mesh = meshioplusplus.read(filename)
    tol = 1.0e-5
    s = np.sum(mesh.points)
    assert abs(s - ref_sum) < tol * abs(ref_sum)
    assert mesh.cells[0].type == "triangle"
    assert len(mesh.cells[0].data) == ref_num_cells
