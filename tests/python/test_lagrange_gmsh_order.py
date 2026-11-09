"""The fixed-order high cells as VTK Lagrange cells (roadmap 6.3.2).

meshio++ holds ``triangle10``, ``quad16``, ``tetra20``, ``wedge40``,
``hexahedron64`` and their higher orders in Gmsh's node order, which
``_lagrange.gmsh_lattice`` generates. The VTK writers write them as VTK
Lagrange cells through ``gmsh_to_vtk_order``. The digests below are Gmsh
4.15's own node lists (``tools/gen_gmsh_node_order_digests.py``, which runs
Gmsh); the generator is held to them, and to Gmsh itself where it is
installed. The writes are held to VTK's own measure of each cell: a node in
the wrong place would change it.
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _lagrange
from meshioplusplus.vtu import _vtu

#: Gmsh 4.15's node order of every element it has from order 3, as digests of
#: the lattice coordinates (tools/gen_gmsh_node_order_digests.py).
GMSH_DIGESTS = {
    "line:3": "1cd3ecd998db55f6",
    "line:4": "ed1a7725fa457278",
    "line:5": "e73b1cfba48d6da3",
    "line:6": "79135c949ef48cc6",
    "line:7": "62f3903809f2214a",
    "line:8": "151301039ae338f9",
    "line:9": "8136c4f8add1dc70",
    "line:10": "fe7e45d07adf3c60",
    "triangle:3": "49025b3bd153b585",
    "triangle:4": "d47f2a97da073dc1",
    "triangle:5": "c85552ccde4b423a",
    "triangle:6": "03487771713c9f6d",
    "triangle:7": "aa823258052a6f03",
    "triangle:8": "71a40a344308680c",
    "triangle:9": "eb92d25991c84e37",
    "triangle:10": "d1468bb4d6cfbb66",
    "quad:3": "c87f42a74956bceb",
    "quad:4": "28dd879c97a27037",
    "quad:5": "2680af44c87a53d3",
    "quad:6": "706489ca96cd20c1",
    "quad:7": "2e17f3ee46dc6fd1",
    "quad:8": "2cbaa2a7783db8c7",
    "quad:9": "51b11c5796334e7f",
    "quad:10": "de318d39ad4ff383",
    "tetra:3": "876b5245c2ea0fb7",
    "tetra:4": "9f4faa6385c05478",
    "tetra:5": "8f748646bfe57c48",
    "tetra:6": "636827126d432b1f",
    "tetra:7": "f2e1a361ccde538c",
    "tetra:8": "1be548da02892432",
    "tetra:9": "0e1fdecd0eb3b208",
    "tetra:10": "b0400f7645728d82",
    "hexahedron:3": "81aae557fc5a840d",
    "hexahedron:4": "25139bc464f11ec8",
    "hexahedron:5": "bb9efd55a3b7c5fe",
    "hexahedron:6": "969997b90d17b78f",
    "hexahedron:7": "2e63b02db701ceaf",
    "hexahedron:8": "b9a79641316b75ab",
    "hexahedron:9": "10ebfcd190ed28dc",
    "wedge:3": "fd3329b6d28fddcc",
    "wedge:4": "31be7a4bf9658af8",
    "wedge:5": "bf60bffb608f9dc4",
    "wedge:6": "6505c275165c8749",
    "wedge:7": "b9e8e11ed7eca620",
    "wedge:8": "031ca76b30d299c3",
    "wedge:9": "9cc5e5228f3b2522",
}


def _digest(shape, p):
    import importlib.util
    import pathlib

    tool = pathlib.Path(__file__).resolve().parents[2] / "tools"
    spec = importlib.util.spec_from_file_location(
        "gen_gmsh_node_order_digests", tool / "gen_gmsh_node_order_digests.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.digest(_lagrange.gmsh_lattice(shape, p))


@pytest.mark.parametrize("key", sorted(GMSH_DIGESTS))
def test_the_generator_is_gmshs_order(key):
    shape, p = key.split(":")
    assert _digest(shape, int(p)) == GMSH_DIGESTS[key]


def test_against_gmsh_itself():
    pytest.importorskip("gmsh")
    import importlib.util
    import pathlib

    tool = pathlib.Path(__file__).resolve().parents[2] / "tools"
    spec = importlib.util.spec_from_file_location(
        "gen_gmsh_node_order_digests", tool / "gen_gmsh_node_order_digests.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.gmsh_digests() == GMSH_DIGESTS


def test_the_permutation_is_a_permutation():
    for key in GMSH_DIGESTS:
        shape, p = key.split(":")
        perm = _lagrange.gmsh_to_vtk_order(shape, int(p))
        assert sorted(perm.tolist()) == list(range(len(perm)))


#: One cell of each family, its order, and its VTK measure (length, area or
#: volume of the reference cell). line4 is VTK's own cubic line, so the line
#: is taken at order 4.
_CELLS = [
    ("line5", "line", 4, 1.0),
    ("triangle10", "triangle", 3, 0.5),
    ("quad16", "quad", 3, 1.0),
    ("tetra20", "tetra", 3, 1 / 6),
    ("wedge40", "wedge", 3, 0.5),
    ("hexahedron64", "hexahedron", 3, 1.0),
    ("hexahedron125", "hexahedron", 4, 1.0),
    ("tetra35", "tetra", 4, 1 / 6),
]


def _reference_cell(name, shape, p=3):
    """The order-p reference cell, its nodes in Gmsh's order. A node VTK
    takes for another distorts its Lagrange map, which changes the measure."""
    points = np.array(_lagrange.gmsh_lattice(shape, p), dtype=float) / p
    return meshioplusplus.Mesh(points, [(name, [list(range(len(points)))])])


def _vtk_measure(path):
    vtk = pytest.importorskip("vtk")
    reader = (
        vtk.vtkXMLUnstructuredGridReader()
        if path.suffix == ".vtu"
        else vtk.vtkUnstructuredGridReader()
    )
    reader.SetFileName(str(path))
    reader.Update()
    sizes = vtk.vtkCellSizeFilter()
    sizes.SetInputData(reader.GetOutput())
    sizes.Update()
    data = sizes.GetOutput().GetCellData()
    for name in ("Volume", "Area", "Length"):
        value = data.GetArray(name).GetValue(0)
        if value:
            return value
    return 0.0


@pytest.mark.parametrize("name, shape, p, measure", _CELLS)
@pytest.mark.parametrize("engine", ["core", "python"])
@pytest.mark.parametrize("suffix", [".vtu", ".vtk"])
def test_written_as_vtk_lagrange_cells(
    tmp_path, name, shape, p, measure, engine, suffix
):
    """roadmap 6.3.2: the VTU writers raised on these types; VTK reads what
    they write now as the reference cell it is."""
    from meshioplusplus.vtk import _main as _vtk_py

    mesh = _reference_cell(name, shape, p)
    path = tmp_path / f"cell{suffix}"
    if engine == "python":
        (_vtu.write if suffix == ".vtu" else _vtk_py.write)(path, mesh)
    else:
        meshioplusplus.write(path, mesh)
    assert _vtk_measure(path) == pytest.approx(measure, rel=1e-9)
    back = meshioplusplus.read(path)
    assert back.cells[0].type.startswith("VTK_LAGRANGE_")
    # Read back in VTK's order: node j is where VTK's lattice puts it.
    expected = np.array(_lagrange.vtk_lattice(shape, p), dtype=float) / p
    got = back.points[back.cells[0].data[0]]
    np.testing.assert_allclose(got[:, : expected.shape[1]], expected, atol=1e-12)


def test_the_core_writes_them_itself(tmp_path, monkeypatch):
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    for name, shape, p, _ in _CELLS:
        path = tmp_path / f"{name}.vtu"
        meshioplusplus.write(path, _reference_cell(name, shape, p))
        assert meshioplusplus.read(path).cells[0].type.startswith("VTK_LAGRANGE_")
