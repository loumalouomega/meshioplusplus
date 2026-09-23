import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus.ansys import _ansys

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.empty_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.quad_mesh,
        helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.hex_mesh,
        helpers.pyramid_mesh,
        helpers.wedge_mesh,
    ],
)
@pytest.mark.parametrize("binary", [False, True])
def test(mesh, binary, tmp_path):
    def writer(*args, **kwargs):
        return meshioplusplus.ansys.write(*args, binary=binary, **kwargs)

    helpers.write_read(tmp_path, writer, meshioplusplus.ansys.read, mesh, 1.0e-15)


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
