import numpy as np
import pytest

import meshioplusplus
from meshioplusplus._sniff import _sniff_format_py


@pytest.mark.parametrize(
    "contents,expected",
    [
        (b"# vtk DataFile Version 3.0\n", "vtk"),
        (b"$MeshFormat\n2.2 0 8\n", "gmsh"),
        (b"ply\nformat ascii 1.0\n", "ply"),
        (b"OFF\n8 6 0\n", "off"),
        (b"solid mysolid\n facet normal 0 0 1\n", "stl"),
        (b'<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid">', "vtu"),
        (b'<?xml version="1.0"?>\n<VTKFile type="PolyData">', "vtp"),
        (b"*Heading\n test\n*Node\n", "abaqus"),
        (b"*KEYWORD\n*NODE\n", "lsdyna"),
        (b"$ a comment\n$ another\n*keyword long=y\n*NODE\n", "lsdyna"),
        (b"% written by Salome\nTITRE\n mesh\nFINSF\n", "code_aster"),
        (b"coor_3d NOM=INDEFINI\n N1 0 0 0\nFINSF\nFIN\n", "code_aster"),
        (b"    1C\n    1UCALCULIX\n    2C\n", "frd"),
        (b"    1C\r\n    2C                            20\r\n", "frd"),
        (b"GiD Post Results File 1.2\n", "gid"),
        (b"# Created by COMSOL\n0 1 \n1 # number of tags\n5 mesh1 \n", "mphtxt"),
        (
            b"\x00\x00\x00\x00\x01\x00\x00\x00\x01\x00\x00\x00\x05\x00\x00\x00m\x00\x00\x00",
            "mphbin",
        ),
        (b'MESH "m" dimension 3 ElemType Triangle Nnode 3\n', "gid"),
        (b'<?xml version="1.0"?>\n<febio_spec version="4.0">\n', "febio"),
        (b"BEF\x00\x00\x00\x00\x01", "xplt"),
    ],
)
def test_recognizes_signatures(tmp_path, contents, expected):
    f = tmp_path / "mesh.dat"
    f.write_bytes(contents)
    assert meshioplusplus.sniff_format(f) == expected
    # the pure-python twin agrees with the C++ core
    assert _sniff_format_py(f) == expected


@pytest.mark.parametrize(
    "contents",
    [
        b"\x89HDF\r\n\x1a\n----------",
        b"just some random text\n",
        # A bare "MESH " with no opening quote is exactly the generic English
        # token sniff_format's own contract refuses to claim (medit's .mesh,
        # FreeFem output and hand-written headers all start this way). Only
        # `MESH "` is unambiguous.
        b"MESH something else\n",
    ],
)
def test_ambiguous_returns_empty(tmp_path, contents):
    f = tmp_path / "mesh.dat"
    f.write_bytes(contents)
    assert meshioplusplus.sniff_format(f) == ""


def test_missing_file_returns_empty():
    assert meshioplusplus.sniff_format("/nonexistent/path/xyz.dat") == ""


def test_read_falls_back_to_sniff(tmp_path):
    # Write a real VTU, then read it through a file whose extension is unknown.
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    src = tmp_path / "mesh.vtu"
    meshioplusplus.write(src, mesh)
    unknown = tmp_path / "mesh.unknownext"
    unknown.write_bytes(src.read_bytes())
    back = meshioplusplus.read(unknown)  # no file_format -> extension unknown -> sniff
    assert len(back.points) == 4


def _touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("x\n")


@pytest.mark.parametrize(
    "files, sub, expected",
    [
        (["mesh.header"], "", "elmer"),
        (["mesh.header"], "mesh.header", "elmer"),
        (["partitioning.2/part.1.header"], "partitioning.2", "elmer"),
        (["partitioning.2/part.1.header"], "", "elmer"),
        (["constant/polyMesh/owner", "constant/polyMesh/faces"], "", "openfoam"),
        (
            ["constant/polyMesh/owner", "constant/polyMesh/faces"],
            "constant/polyMesh",
            "openfoam",
        ),
        (["processor0/constant/polyMesh/"], "", "openfoam"),
        (["constant/regionProperties"], "", "openfoam"),
        # both, neither, or half a polyMesh: not a guess
        (["mesh.header", "polyMesh/owner", "polyMesh/faces"], "", ""),
        ([], "", ""),
        (["polyMesh/owner"], "", ""),
    ],
)
def test_recognizes_directories(tmp_path, files, sub, expected):
    root = tmp_path / "case"
    root.mkdir()
    for name in files:
        if name.endswith("/"):
            (root / name).mkdir(parents=True)
        else:
            _touch(root / name)
    target = root / sub if sub else root
    assert meshioplusplus.sniff_format(target) == expected
    assert _sniff_format_py(target) == expected


def test_read_sniffs_an_openfoam_case_directory(tmp_path):
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    case = tmp_path / "case"
    try:
        meshioplusplus.write(case, mesh, file_format="openfoam")
    except meshioplusplus.WriteError:
        pytest.skip("openfoam writing needs the C++ core")
    back = meshioplusplus.read(case)  # no file_format, no extension -> directory sniff
    assert len(back.points) == 4
    # A write to an extension-less path still needs its format named.
    with pytest.raises((meshioplusplus.ReadError, meshioplusplus.WriteError)):
        meshioplusplus.write(tmp_path / "other", mesh)
