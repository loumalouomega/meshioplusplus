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


def test_decompress_standalone(tmp_path):
    # `decompress` on a compressed VTU produces a still-readable, larger file.
    infile = tmp_path / "in.vtu"
    meshioplusplus.write(infile, helpers.tri_mesh, compression="zlib")
    meshioplusplus._cli.main(["decompress", str(infile)])
    mesh = meshioplusplus.read(infile)
    assert is_same_mesh(helpers.tri_mesh, mesh, atol=1.0e-12)


def test_compress_max_flag(tmp_path):
    # The `--max` flag switches VTU to lzma; the result must round-trip.
    infile = tmp_path / "in.vtu"
    meshioplusplus.write(infile, helpers.tri_mesh, binary=True, compression=None)
    meshioplusplus._cli.main(["compress", str(infile), "--max"])
    mesh = meshioplusplus.read(infile)
    assert is_same_mesh(helpers.tri_mesh, mesh, atol=1.0e-12)


def test_convert_ascii_and_float_format(tmp_path):
    # Exercise the `-a`/`--ascii` and `-f`/`--float-format` flag plumbing.
    infile = tmp_path / "in.vtu"
    outfile = tmp_path / "out.vtk"
    meshioplusplus.write(infile, helpers.tri_mesh)
    meshioplusplus._cli.main(["convert", str(infile), str(outfile), "-a", "-f", ".8e"])
    mesh = meshioplusplus.read(outfile)
    assert is_same_mesh(helpers.tri_mesh, mesh, atol=1.0e-7)


def test_convert_int_data_to_sets(tmp_path):
    # The `-d`/`--int-data-to-sets` path (converse of `--sets-to-int-data`).
    # Regression test: `convert -d` used to raise "dictionary changed size
    # during iteration" because *_data_to_sets(key) mutates the dict being
    # iterated (fixed by snapshotting the keys in _cli/_convert.py).
    mesh = helpers.tri_mesh.copy()
    # Both loops (point_data and cell_data) must survive the in-place mutation.
    mesh.point_data = {"fixed-loose": np.array([0, 0, 1, 1])}
    mesh.cell_data = {"grain0-grain1": [np.array([0, 1])]}

    infile = tmp_path / "in.vtu"
    outfile = tmp_path / "out.vtu"
    meshioplusplus.write(infile, mesh)
    meshioplusplus._cli.main(["convert", str(infile), str(outfile), "-d"])
    # The command completes and the geometry survives.
    out = meshioplusplus.read(outfile)
    assert np.allclose(out.points, mesh.points)


def test_ascii_unsupported_format_returns_1(tmp_path):
    # `.off` is not in the ascii command's whitelist -> error + `return 1`.
    infile = tmp_path / "in.off"
    meshioplusplus.write(infile, helpers.tri_mesh)
    rc = meshioplusplus._cli.main(["ascii", str(infile)])
    assert rc == 1


@pytest.mark.parametrize("command", ["binary", "compress", "decompress"])
def test_unsupported_format_exits_1(command, tmp_path):
    # binary/compress/decompress use `exit(1)` (SystemExit), unlike ascii.
    infile = tmp_path / "in.off"
    meshioplusplus.write(infile, helpers.tri_mesh)
    with pytest.raises(SystemExit) as exc:
        meshioplusplus._cli.main([command, str(infile)])
    assert exc.value.code == 1


def test_info_inconsistent_points_warns(monkeypatch, capsys):
    # A mesh whose cells reference a nonexistent point triggers the
    # "Inconsistent mesh" warning branch of `info` (_info.py:28).
    bad = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        [("triangle", np.array([[0, 1, 9]]))],
    )
    monkeypatch.setattr(meshioplusplus._cli._info, "read", lambda *a, **k: bad)
    meshioplusplus._cli.main(["info", "ignored"])
    assert "Inconsistent" in capsys.readouterr().err


def test_info_unused_points_warns(monkeypatch, capsys):
    # A consistent mesh with an unused point triggers the second warning
    # branch of `info` (_info.py:37).
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [5.0, 5.0, 5.0]],
        [("triangle", np.array([[0, 1, 2]]))],
    )
    monkeypatch.setattr(meshioplusplus._cli._info, "read", lambda *a, **k: mesh)
    meshioplusplus._cli.main(["info", "ignored"])
    assert "not part of any cell" in capsys.readouterr().err
