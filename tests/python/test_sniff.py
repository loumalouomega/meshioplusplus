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
        (b"libMesh-1.3.0\n2\t # number of elements\n", "libmesh"),
        (b"\x00\x00\x00\x0dlibMesh-1.8.0\x00\x00\x00", "libmesh"),
        (b"*I 19I 41921A6.23-1  A07-Nov-2A024     A16:50:01I 11I 18", "abaqus_fil"),
        (
            b"\x00\x10\x00\x00"
            + (9).to_bytes(8, "little")
            + (1921).to_bytes(8, "little"),
            "abaqus_fil",
        ),
        (
            b"\x00\x00\x10\x00" + (9).to_bytes(8, "big") + (1921).to_bytes(8, "big"),
            "abaqus_fil",
        ),
        (b"#RADIOSS STARTER\n/BEGIN\nrun\n", "radioss"),
        # Nastran OP2 header blocks: [3] date [7] tape code, 4- and 8-byte words
        (
            b"".join(
                n.to_bytes(4, "little") + p + n.to_bytes(4, "little")
                for p in (
                    (3).to_bytes(4, "little"),
                    b"\x09\x00\x00\x00\x18\x00\x00\x00\x1a\x00\x00\x00",
                    (7).to_bytes(4, "little"),
                    b"NASTRAN FORT TAPE ID CODE - ",
                )
                for n in [len(p)]
            ),
            "nastran_op2",
        ),
        (
            b"".join(
                n.to_bytes(4, "big") + p + n.to_bytes(4, "big")
                for p in (
                    (3).to_bytes(8, "big"),
                    (9).to_bytes(8, "big") * 3,
                    (7).to_bytes(8, "big"),
                    b"NAST    RAN     FORT    TAPE    ID C    ODE     -       ",
                )
                for n in [len(p)]
            ),
            "nastran_op2",
        ),
        # LS-DYNA d3plot control block (version 971, NDIM 4, one shell)
        (
            b" " * 40
            + b"".join((v).to_bytes(4, "little", signed=True) for v in [0, 1, 0, 0])
            + np.float32(971.0).tobytes()
            + b"".join(
                (v).to_bytes(4, "little", signed=True)
                for v in [4, 4, 6, 0, 0, 1, 0, 0] + [0] * 8 + [1, 1, 7] + [0] * 30
            ),
            "lsdyna_d3plot",
        ),
        (b"# a comment\n$ another\n/BEGIN\nrun\n      2019         0\n", "radioss"),
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
        # A starred line that is not a record of I items (an Abaqus deck's
        # `*Node` is claimed by the abaqus rule instead).
        b"*I am not a results file\n",
        b"/NODE\n1 0 0 0\n",
    ],
)
def test_ambiguous_returns_empty(tmp_path, contents):
    f = tmp_path / "mesh.dat"
    f.write_bytes(contents)
    assert meshioplusplus.sniff_format(f) == ""


def test_z88_is_recognised_by_its_file_name(tmp_path):
    f = tmp_path / "Z88I1.TXT"
    f.write_bytes(b"3 4 1 12 0\n")
    assert meshioplusplus.sniff_format(f) == "z88"
    assert _sniff_format_py(f) == "z88"
    other = tmp_path / "points.txt"
    other.write_bytes(b"3 4 1 12 0\n")
    assert meshioplusplus.sniff_format(other) != "z88"


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
