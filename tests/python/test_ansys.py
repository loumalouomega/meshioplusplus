import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.ansys import _ansys

from . import helpers

_ROUND_TRIP = [
    helpers.tri_mesh,
    helpers.tri_mesh_2d,
    helpers.quad_mesh,
    helpers.tri_quad_mesh,
    helpers.tet_mesh,
    helpers.hex_mesh,
    helpers.pyramid_mesh,
    helpers.wedge_mesh,
]


def _cell_set(mesh, dim):
    """(type, sorted nodes) of every cell of dimension ``dim``, as points."""
    out = []
    for block in mesh.cells:
        if block.dim != dim:
            continue
        for row in np.asarray(block.data):
            out.append(
                (block.type, tuple(sorted(map(tuple, mesh.points[row].tolist()))))
            )
    return sorted(out)


@pytest.mark.parametrize("mesh", _ROUND_TRIP)
@pytest.mark.parametrize("binary", [False, True])
@pytest.mark.parametrize("engine", ["core", "python"])
def test(mesh, binary, engine, tmp_path):
    """A written file reads back as the same cells (rebuilt from their faces)
    plus the boundary facets; a planar mesh as 2-D points."""
    path = tmp_path / "out.msh"
    if engine == "core":
        meshioplusplus.ansys.write(path, mesh, binary=binary)
    else:
        _ansys.write(path, mesh, binary=binary)
    out = meshioplusplus.ansys.read(path)
    dim = max(b.dim for b in mesh.cells)
    assert out.points.shape[1] == dim
    points = np.asarray(mesh.points)[:, :dim]
    np.testing.assert_allclose(out.points, points, rtol=0, atol=1e-15)
    ref = meshioplusplus.Mesh(points, mesh.cells)
    assert _cell_set(out, dim) == _cell_set(ref, dim)
    # every boundary facet comes back as a (dim - 1) cell in a wall zone
    assert any(b.dim == dim - 1 for b in out.cells)
    assert "ansys:zone" in out.cell_data


def test_engines_write_the_same_bytes(tmp_path):
    for mesh in _ROUND_TRIP + [meshioplusplus.read(p) for p in sorted(_FLUENT_FILES())]:
        for binary in (False, True):
            a, b = tmp_path / "a.msh", tmp_path / "b.msh"
            meshioplusplus.ansys.write(a, mesh, binary=binary)
            _ansys.write(b, mesh, binary=binary)
            assert a.read_bytes() == b.read_bytes()


def test_an_empty_mesh_is_refused(tmp_path):
    with pytest.raises(meshioplusplus.WriteError):
        meshioplusplus.ansys.write(tmp_path / "e.msh", helpers.empty_mesh)


@pytest.mark.parametrize("engine", ["core", "python"])
def test_a_surface_off_the_xy_plane_is_refused(engine, tmp_path):
    """A 2-D Fluent mesh lies in z = 0; a 3-D surface is not dropped to it."""
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.5]], [("triangle", [[0, 1, 2]])]
    )
    write = meshioplusplus.ansys.write if engine == "core" else _ansys.write
    with pytest.raises(meshioplusplus.WriteError, match="z = 0 plane"):
        write(tmp_path / "s.msh", mesh)
    if engine == "core":
        with pytest.raises(meshioplusplus.WriteError, match="z = 0 plane"):
            _core.ansys_write(str(tmp_path / "s.msh"), mesh, True)


def _FLUENT_FILES():
    return (pathlib.Path(__file__).resolve().parent / "meshes" / "ansys").glob("*.msh")


# --- malformed-input / error-path coverage ---


def test_ansys_garbage_raises(tmp_path):
    p = tmp_path / "bad.msh"
    p.write_text("this is not an ansys mesh file\n")
    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.read(p, file_format="ansys")


# Fluent meshes whose cells exist only through their faces
# (tools/gen_feconv_quirk_fixtures.py).
_FLUENT = pathlib.Path(__file__).resolve().parent / "meshes" / "ansys"


def _core_read(path):
    from meshioplusplus import _core

    return _core.ansys_read(str(path))


_READERS = pytest.mark.parametrize(
    "reader", [_core_read, _ansys.read], ids=["core", "python"]
)


def _regions(mesh):
    return {r.name: (r.dim, r.tag, len(r.entries)) for r in mesh.regions}


@_READERS
def test_cells_from_faces_3d(reader):
    mesh = reader(_FLUENT / "cells3d.msh")
    types = [c.type for c in mesh.cells]
    assert types == [
        "tetra",
        "hexahedron",
        "polyhedron8",
        "wedge",
        "pyramid",
        "triangle",
        "quad",
    ]
    assert len(mesh.cells[0].data) == 2
    helpers.assert_well_formed(mesh)
    # The polyhedron keeps its seven faces, each wound outward.
    (poly,) = mesh.cells[2].data
    assert len(poly) == 7
    pts = mesh.points
    centre = pts[np.unique(np.concatenate(poly))].mean(axis=0)
    for f in poly:
        f = np.asarray(f)
        normal = np.cross(pts[f[1]] - pts[f[0]], pts[f[2]] - pts[f[0]])
        assert np.dot(normal, pts[f].mean(axis=0) - centre) > 0
    # Boundary faces are kept and wound outward; the interior face is not.
    assert sum(len(c.data) for c in mesh.cells[5:]) == 29
    owners = [set(np.asarray(c).tolist()) for b in mesh.cells[:2] for c in b.data]
    owners += [set(np.unique(np.concatenate(poly)).tolist())]
    owners += [set(np.asarray(c).tolist()) for b in mesh.cells[3:5] for c in b.data]
    for block in mesh.cells[5:]:
        for f in np.asarray(block.data):
            (cell,) = [c for c in owners if set(f.tolist()) <= c]
            normal = np.cross(pts[f[1]] - pts[f[0]], pts[f[2]] - pts[f[0]])
            inside = pts[sorted(cell)].mean(axis=0)
            assert np.dot(normal, pts[f].mean(axis=0) - inside) > 0
    zones = [np.asarray(z).tolist() for z in mesh.cell_data["ansys:zone"]]
    assert zones[:5] == [[4, 4], [4], [4], [4], [4]]
    assert _regions(mesh) == {"block-of-cells": (3, 4, 6), "outer-walls": (2, 3, 29)}


@_READERS
@pytest.mark.parametrize("name", ["tgrid2d", "gambit2d"])
def test_cells_from_edges_2d(reader, name):
    mesh = reader(_FLUENT / f"{name}.msh")
    assert mesh.points.shape == (5, 2)
    np.testing.assert_array_equal(mesh.points[3], [0.0, 1.0])
    assert [c.type for c in mesh.cells] == ["triangle", "quad", "line"]

    # counter-clockwise, from whichever node the cell's first edge starts at
    def rotations(ring):
        return [ring[k:] + ring[:k] for k in range(len(ring))]

    assert np.asarray(mesh.cells[0].data)[0].tolist() in rotations([1, 4, 2])
    assert np.asarray(mesh.cells[1].data)[0].tolist() in rotations([0, 1, 2, 3])
    assert len(mesh.cells[2].data) == 5
    regions = _regions(mesh)
    assert regions["plate"] == (2, 10, 2)


@_READERS
def test_binary_sections(reader):
    mesh = reader(_FLUENT / "binary3d.msh")
    assert [c.type for c in mesh.cells] == ["tetra", "triangle"]
    assert len(mesh.cells[0].data) == 2 and len(mesh.cells[1].data) == 6
    helpers.assert_well_formed(mesh)
    assert _regions(mesh)["tets"] == (3, 4, 2)


def test_engines_agree_on_the_fixtures():
    for path in sorted(_FLUENT.glob("*.msh")):
        a, b = _core_read(path), _ansys.read(path)
        np.testing.assert_array_equal(a.points, b.points)
        assert [c.type for c in a.cells] == [c.type for c in b.cells]
        for x, y in zip(a.cells, b.cells):
            if x.type.startswith("polyhedron"):
                assert [[list(f) for f in c] for c in x.data] == [
                    [list(f) for f in c] for c in y.data
                ]
            else:
                np.testing.assert_array_equal(x.data, y.data)
        assert _regions(a) == _regions(b)
