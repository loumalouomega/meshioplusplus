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
        (b"GiD Post Results File 1.2\n", "gid"),
        (b'MESH "m" dimension 3 ElemType Triangle Nnode 3\n', "gid"),
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
