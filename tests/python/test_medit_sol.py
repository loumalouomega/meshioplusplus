"""Medit solutions (``.sol``/``.solb``), roadmap 3.2.6: both engines."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.medit_sol import _medit_sol as py_sol

ENGINES = ["core", "python"]


def _read(engine, path):
    return _core.medit_sol_read(str(path)) if engine == "core" else py_sol.read(path)


def _write(engine, path, mesh):
    if engine == "core":
        _core.medit_sol_write(str(path), mesh)
    else:
        py_sol.write(path, mesh)


def _attach(engine, path, mesh):
    if engine == "core":
        return _core.medit_sol_read_into(str(path), mesh)
    dim, sections = py_sol.read_sections(path)
    return py_sol.attach(mesh, dim, sections)


def _square(dim=3):
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0.0]])[:, :dim]
    return meshioplusplus.Mesh(
        points,
        [
            ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
            ("line", np.array([[0, 1], [1, 2], [2, 3]])),
        ],
    )


def _with_fields(dim):
    mesh = _square(dim)
    rng = np.random.default_rng(3)
    mesh.point_data["a_scalar"] = rng.random(4)
    mesh.point_data["b_vector"] = rng.random((4, dim))
    mesh.point_data["c_metric"] = rng.random((4, dim * (dim + 1) // 2))
    mesh.point_data["d_matrix"] = rng.random((4, dim * dim))
    mesh.point_data["medit:ref"] = np.array([1, 2, 3, 4])
    return mesh


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("ext", [".sol", ".solb"])
@pytest.mark.parametrize("dim", [2, 3])
def test_vertex_fields_round_trip(engine, ext, dim, tmp_path):
    mesh = _with_fields(dim)
    path = tmp_path / f"out{ext}"
    _write(engine, path, mesh)
    sol = _read(engine, path)
    assert sol.points.shape == (4, 0)
    assert int(sol.field_data["medit:dimension"][0]) == dim
    # file order is the sorted name order; integer arrays are not solutions
    names = ["a_scalar", "b_vector", "c_metric", "d_matrix"]
    assert sorted(sol.point_data) == [f"medit:sol:{i}" for i in range(4)]
    for i, name in enumerate(names):
        got = sol.point_data[f"medit:sol:{i}"]
        assert got.shape == mesh.point_data[name].shape
        np.testing.assert_array_equal(got, mesh.point_data[name])
    # and back out of the geometry-less mesh, its dimension from field_data
    again = tmp_path / f"again{ext}"
    _write(engine, again, sol)
    assert again.read_bytes() == path.read_bytes()


@pytest.mark.parametrize("ext", [".sol", ".solb"])
def test_engines_write_the_same_bytes(ext, tmp_path):
    mesh = _with_fields(3)
    mesh.cell_data["e_cell"] = [np.array([0.5, 1.5]), np.array([1.0, 2.0, 3.0])]
    _write("core", tmp_path / f"c{ext}", mesh)
    _write("python", tmp_path / f"p{ext}", mesh)
    assert (tmp_path / f"c{ext}").read_bytes() == (tmp_path / f"p{ext}").read_bytes()


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("ext", [".sol", ".solb"])
def test_cell_fields_attach_to_their_mesh(engine, ext, tmp_path):
    mesh = _square()
    mesh.cell_data["size"] = [np.array([0.25, 0.5]), np.full(3, np.nan)]
    mesh.cell_data["grad"] = [np.ones((2, 3)), np.arange(9.0).reshape(3, 3)]
    path = tmp_path / f"out{ext}"
    _write(engine, path, mesh)
    got = _attach(engine, path, _square())
    # sorted names: grad (all types), then size (triangles only)
    np.testing.assert_array_equal(got.cell_data["medit:sol:0"][0], np.ones((2, 3)))
    np.testing.assert_array_equal(
        got.cell_data["medit:sol:0"][1], np.arange(9.0).reshape(3, 3)
    )
    np.testing.assert_array_equal(got.cell_data["medit:sol:1"][0], [0.25, 0.5])
    # the edges had no second field: the name's fill
    assert np.all(np.isnan(got.cell_data["medit:sol:1"][1]))
    assert got.cell_data["medit:sol:1"][1].shape == (3,)


@pytest.mark.parametrize("ext", [".mesh", ".meshb"])
def test_medit_companion_route(ext, tmp_path):
    mesh = _with_fields(3)
    mesh.cell_data["size"] = [np.array([0.25, 0.5]), np.array([1.0, 2.0, 3.0])]
    path = tmp_path / f"square{ext}"
    meshioplusplus.medit.write(path, mesh, sol=True)
    companion = tmp_path / ("square.solb" if ext.endswith("b") else "square.sol")
    assert companion.exists()
    got = meshioplusplus.medit.read(path, sol=True)
    np.testing.assert_array_equal(
        got.point_data["medit:sol:2"], mesh.point_data["c_metric"]
    )
    # the C++ ascii writer groups blocks by keyword: find the edges by type
    lines = next(k for k, cb in enumerate(got.cells) if cb.type == "line")
    np.testing.assert_array_equal(got.cell_data["medit:sol:0"][lines], [1.0, 2.0, 3.0])
    # a named path works too, and the mesh's own data stays
    got = meshioplusplus.medit.read(path, sol=companion)
    np.testing.assert_array_equal(got.point_data["medit:ref"], [1, 2, 3, 4])


def test_medit_sol_true_without_a_companion_raises(tmp_path):
    path = tmp_path / "lonely.mesh"
    meshioplusplus.medit.write(path, _square())
    with pytest.raises(meshioplusplus.ReadError, match="lonely.sol"):
        meshioplusplus.medit.read(path, sol=True)


_MMG_STYLE = """MeshVersionFormatted 2

Dimension 2
# a metric and a scalar, as MMG writes them
SolAtVertices
3
2 3 1
1.0 0.0 2.0   0.5
3.0 1.0 4.0   0.25
5.0 2.0 6.0   0.125

Time
0.5

End
"""


@pytest.mark.parametrize("engine", ENGINES)
def test_reads_an_mmg_style_file(engine, tmp_path):
    path = tmp_path / "metric.sol"
    path.write_text(_MMG_STYLE)
    sol = _read(engine, path)
    np.testing.assert_array_equal(
        sol.point_data["medit:sol:0"], [[1, 0, 2], [3, 1, 4], [5, 2, 6]]
    )
    np.testing.assert_array_equal(sol.point_data["medit:sol:1"], [0.5, 0.25, 0.125])
    assert int(sol.field_data["medit:dimension"][0]) == 2


def _big_endian_solb():
    """Version 2 (4-byte positions), 2-D, two scalar fields at two vertices."""
    out = [np.array([1, 2, 3, 20, 2], ">i4").tobytes()]
    body = np.array([2, 2, 1, 1], ">i4").tobytes()  # rows, fields, their types
    body += np.array([[1.0, 2.0], [3.0, 4.0]], ">f8").tobytes()
    out.append(np.array([62, 20 + 8 + len(body)], ">i4").tobytes())
    out.append(body)
    out.append(np.array([54], ">i4").tobytes())
    return b"".join(out)


@pytest.mark.parametrize("engine", ENGINES)
def test_reads_a_big_endian_solb(engine, tmp_path):
    path = tmp_path / "be.solb"
    path.write_bytes(_big_endian_solb())
    sol = _read(engine, path)
    np.testing.assert_array_equal(sol.point_data["medit:sol:0"], [1.0, 3.0])
    np.testing.assert_array_equal(sol.point_data["medit:sol:1"], [2.0, 4.0])


@pytest.mark.parametrize("engine", ENGINES)
def test_standalone_read_skips_element_sections(engine, tmp_path):
    mesh = _square()
    mesh.cell_data["size"] = [np.array([0.25, 0.5]), np.array([1.0, 2.0, 3.0])]
    path = tmp_path / "cells.sol"
    _write(engine, path, mesh)
    sol = _read(engine, path)
    assert sol.points.shape == (0, 0)
    assert not sol.point_data and not sol.cells


def test_generic_read_and_write_dispatch(tmp_path):
    mesh = _with_fields(3)
    meshioplusplus.write(tmp_path / "x.sol", mesh)
    meshioplusplus.write(tmp_path / "x.solb", mesh)
    for name in ("x.sol", "x.solb"):
        sol = meshioplusplus.read(tmp_path / name)
        np.testing.assert_array_equal(
            sol.point_data["medit:sol:1"], mesh.point_data["b_vector"]
        )


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize(
    "text, match",
    [
        ("Dimension 3\nSolAtVertices\n2\n1 1\n1.0\n", "truncated"),
        ("Dimension 3\nSolAtVertices\n1\n1 7\n1.0\n", "solution type"),
        ("Dimension 3\nSolAtVertices\n1\n0\n", "no field"),
        ("SolAtVertices\n1\n1 1\n1.0\n", "Dimension"),
        ("Dimension 4\n", "dimension"),
        ("Dimension 3\nVertices\n0\n", "unknown"),
    ],
)
def test_malformed_files_raise(engine, text, match, tmp_path):
    path = tmp_path / "bad.sol"
    path.write_text(text)
    with pytest.raises(meshioplusplus.ReadError, match=match):
        _read(engine, path)


@pytest.mark.parametrize("engine", ENGINES)
def test_truncated_binary_raises(engine, tmp_path):
    path = tmp_path / "bad.solb"
    path.write_bytes(_big_endian_solb()[:-12])
    with pytest.raises(meshioplusplus.ReadError):
        _read(engine, path)


@pytest.mark.parametrize("engine", ENGINES)
def test_attach_checks_the_mesh(engine, tmp_path):
    mesh = _with_fields(3)
    path = tmp_path / "out.sol"
    _write(engine, path, mesh)
    with pytest.raises(meshioplusplus.ReadError, match="2-D"):
        _attach(engine, path, _square(2))
    small = meshioplusplus.Mesh(np.zeros((3, 3)), [])
    with pytest.raises(meshioplusplus.ReadError, match="rows"):
        _attach(engine, path, small)


@pytest.mark.parametrize("engine", ENGINES)
def test_nothing_to_write_raises(engine, tmp_path):
    mesh = _square()
    mesh.point_data["labels"] = np.array([1, 2, 3, 4])
    with pytest.raises(meshioplusplus.WriteError, match="no floating-point"):
        _write(engine, tmp_path / "out.sol", mesh)


@pytest.mark.parametrize("engine", ENGINES)
def test_a_solution_inside_an_ascii_mesh_is_skipped(engine, tmp_path):
    from meshioplusplus.medit import _medit as py_medit

    path = tmp_path / "with_sol.mesh"
    path.write_text(
        "MeshVersionFormatted 2\nDimension 2\nVertices\n3\n0 0 1\n1 0 1\n0 1 1\n"
        "Triangles\n1\n1 2 3 1\nSolAtVertices\n3\n1 2\n1 2\n3 4\n5 6\nEnd\n"
    )
    mesh = (
        _core.medit_read_ascii(str(path)) if engine == "core" else py_medit.read(path)
    )
    assert len(mesh.points) == 3
    np.testing.assert_array_equal(mesh.cells[0].data, [[0, 1, 2]])
