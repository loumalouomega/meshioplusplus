import numpy as np
import pytest

import meshioplusplus

from . import helpers


def is_same_mesh(mesh0, mesh1, atol):
    if not np.allclose(mesh0.points, mesh1.points, atol=atol, rtol=0.0):
        return False
    for cells0, cells1 in zip(mesh0.cells, mesh1.cells):
        if cells0.type != cells1.type or not np.allclose(cells0.data, cells1.data):
            return False
    return True


def test_info(tmp_path):
    infile = tmp_path / "out.msh"
    meshioplusplus.write(infile, helpers.tri_mesh, file_format="gmsh")
    meshioplusplus._cli.main(["info", str(infile), "--input-format", "gmsh"])


def test_convert(tmp_path):
    input_mesh = helpers.tri_mesh

    infile = tmp_path / "in.msh"
    meshioplusplus.write(infile, helpers.tri_mesh, file_format="gmsh")

    outfile = tmp_path / "out.msh"
    meshioplusplus._cli.main(
        [
            "convert",
            str(infile),
            str(outfile),
            "--input-format",
            "gmsh",
            "--output-format",
            "vtk",
            "--sets-to-int-data",
        ]
    )

    mesh = meshioplusplus.read(outfile, file_format="vtk")

    atol = 1.0e-15
    assert np.allclose(input_mesh.points, mesh.points, atol=atol, rtol=0.0)

    for cells0, cells1 in zip(input_mesh.cells, mesh.cells):
        assert cells0.type == cells1.type
        assert np.allclose(cells0.data, cells1.data)


def test_compress(tmp_path):
    input_mesh = helpers.tri_mesh

    infile = tmp_path / "in.vtu"
    meshioplusplus.write(infile, input_mesh)

    meshioplusplus._cli.main(["decompress", str(infile)])
    mesh = meshioplusplus.read(infile)
    assert is_same_mesh(input_mesh, mesh, atol=1.0e-15)

    meshioplusplus._cli.main(["compress", str(infile)])
    mesh = meshioplusplus.read(infile)
    assert is_same_mesh(input_mesh, mesh, atol=1.0e-15)


def test_ascii_binary(tmp_path):
    input_mesh = helpers.tri_mesh

    infile = tmp_path / "in.vtu"
    meshioplusplus.write(infile, input_mesh)

    meshioplusplus._cli.main(["ascii", str(infile)])
    mesh = meshioplusplus.read(infile)
    assert is_same_mesh(input_mesh, mesh, atol=1.0e-12)

    meshioplusplus._cli.main(["binary", str(infile)])
    mesh = meshioplusplus.read(infile)
    assert is_same_mesh(input_mesh, mesh, atol=1.0e-12)


def test_version():
    # `--version` is an argparse action that prints and exits 0.
    with pytest.raises(SystemExit) as exc:
        meshioplusplus._cli.main(["--version"])
    assert exc.value.code == 0


def test_no_command():
    # A required subcommand is missing -> argparse exits with code 2.
    with pytest.raises(SystemExit) as exc:
        meshioplusplus._cli.main([])
    assert exc.value.code == 2


def test_convert_format_inference(tmp_path):
    # No --input-format/--output-format: the formats are inferred from the
    # extensions, and the `c` alias is used.
    infile = tmp_path / "in.vtu"
    outfile = tmp_path / "out.vtk"
    meshioplusplus.write(infile, helpers.tri_mesh)

    meshioplusplus._cli.main(["c", str(infile), str(outfile)])

    mesh = meshioplusplus.read(outfile)
    assert is_same_mesh(helpers.tri_mesh, mesh, atol=1.0e-12)


def test_info_inference(tmp_path):
    # `info` with extension inference, via the `i` alias.
    infile = tmp_path / "mesh.vtu"
    meshioplusplus.write(infile, helpers.tri_mesh)
    meshioplusplus._cli.main(["i", str(infile)])
