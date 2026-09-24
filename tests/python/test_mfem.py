"""MFEM mesh (``.mesh``) and grid functions (``.gf``): both engines, MFEM's own
sample meshes and fields checked against MFEM's evaluation frozen by
``tools/gen_mfem_fixtures.py`` (order 2 and arbitrary order), non-conforming
meshes, round trips, and the ``.mesh`` clash with Medit."""

import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.mfem import _mfem as py_mfem

from .test_node_order import _is_valid

MESHES = pathlib.Path(__file__).parent / "meshes" / "mfem"
REFERENCE = np.load(MESHES / "reference.npz")
LAGRANGE = np.load(MESHES / "reference_lagrange.npz")
CONFORMING = [
    "star-q2",
    "escher-p2",
    "fichera-q2",
    "fichera-mixed-p2",
    "compass",
    "tinyzoo-3d",
    "periodic-square",
]
# Order 3 and up: VTK Lagrange cells (MFEM's samples, and warped curved meshes).
HIGH_ORDER = [
    "escher-p3",
    "fichera-q3",
    "toroid-wedge",
    "rt-2d-p4-tri",
    "curved-tet-p4",
    "curved-hex-p3",
    "curved-prism-p4",
    "curved-tri-p5",
    "curved-quad-p3u",
    "curved-segment-p3",
]
_VTK_LAGRANGE = {
    68: "VTK_LAGRANGE_CURVE",
    69: "VTK_LAGRANGE_TRIANGLE",
    70: "VTK_LAGRANGE_QUADRILATERAL",
    71: "VTK_LAGRANGE_TETRAHEDRON",
    72: "VTK_LAGRANGE_HEXAHEDRON",
    73: "VTK_LAGRANGE_WEDGE",
}
GRID_FUNCTIONS = {
    "star-q2": ["u", "v"],
    "fichera-q2": ["w"],
    "compass": ["t", "q", "e"],
}

# meshio++ (VTK) non-corner nodes, by the corners they sit between.
_TRI_E = [(0, 1), (1, 2), (2, 0)]
_QUA_E = [(0, 1), (1, 2), (2, 3), (3, 0)]
_SLOTS = {
    "line3": [(0, 1)],
    "triangle6": _TRI_E,
    "quad9": _QUA_E + [(0, 1, 2, 3)],
    "tetra10": [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
    "wedge18": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]
    + [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    "hexahedron27": [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
    + [(0, 4), (1, 5), (2, 6), (3, 7)]
    + [(0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4), (3, 2, 6, 7), (0, 1, 2, 3)]
    + [(4, 5, 6, 7), tuple(range(8))],
}
_CORNERS = {
    "vertex": 1,
    "line": 2,
    "line3": 2,
    "triangle": 3,
    "triangle6": 3,
    "quad": 4,
    "quad9": 4,
    "tetra": 4,
    "tetra10": 4,
    "pyramid": 5,
    "wedge": 6,
    "wedge18": 6,
    "hexahedron": 8,
    "hexahedron27": 8,
}


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read/write signatures."""
    if request.param == "core":
        return meshioplusplus.mfem
    return py_mfem


def _gfs(name):
    return {g: str(MESHES / f"{name}.{g}.gf") for g in GRID_FUNCTIONS.get(name, [])}


def _blocks(mesh):
    return [(b.type, np.asarray(b.data)) for b in mesh.cells]


def _regions(mesh):
    return sorted(
        (r.kind, r.name, r.dim, r.tag, tuple(np.asarray(r.entries).ravel()))
        for r in mesh.regions
    )


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [t for t, _ in _blocks(a)] == [t for t, _ in _blocks(b)]
    for (_, x), (_, y) in zip(_blocks(a), _blocks(b)):
        np.testing.assert_array_equal(x, y)
    assert _regions(a) == _regions(b)
    assert sorted(a.point_data) == sorted(b.point_data)
    for k in a.point_data:
        np.testing.assert_array_equal(a.point_data[k], b.point_data[k])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for k in a.cell_data:
        for x, y in zip(a.cell_data[k], b.cell_data[k]):
            np.testing.assert_array_equal(x, y)


def _node_keys(mesh):
    """(point, sorted vertex key) of every node of every cell. The first points
    are MFEM's vertices, so a corner is its own key."""
    for cell_type, data in _blocks(mesh):
        nc = _CORNERS[cell_type]
        for row in data:
            corners = [int(v) for v in row[:nc]]
            for j, p in enumerate(corners):
                yield p, (p,)
            for k, slot in enumerate(_SLOTS.get(cell_type, [])):
                yield int(row[nc + k]), tuple(sorted(corners[q] for q in slot))


def _reference(name):
    keys = REFERENCE[f"{name}/keys"]
    lookup = {tuple(int(v) for v in k if v >= 0): i for i, k in enumerate(keys)}
    return lookup


@pytest.mark.parametrize("name", CONFORMING)
def test_engines_agree(name):
    gfs = _gfs(name)
    _same(
        meshioplusplus.mfem.read(MESHES / f"{name}.mesh", gfs),
        py_mfem.read(MESHES / f"{name}.mesh", gfs),
    )


@pytest.mark.parametrize(
    "name", sorted(GRID_FUNCTIONS) + ["escher-p2", "fichera-mixed-p2", "tinyzoo-3d"]
)
def test_points_and_fields_match_mfem(engine, name):
    """Every node sits where MFEM's own element transformation puts it, and every
    grid function has MFEM's value there (frozen in reference.npz)."""
    mesh = engine.read(MESHES / f"{name}.mesh", _gfs(name))
    lookup = _reference(name)
    xyz = REFERENCE[f"{name}/xyz"]
    pdim = mesh.points.shape[1]
    checked = 0
    for p, key in _node_keys(mesh):
        i = lookup[key]
        np.testing.assert_allclose(mesh.points[p], xyz[i, :pdim], atol=1e-12)
        for g in GRID_FUNCTIONS.get(name, []):
            if g in mesh.point_data:
                want = REFERENCE[f"{name}/{g}"][i]
                got = np.atleast_1d(mesh.point_data[g][p])
                np.testing.assert_allclose(got, want, rtol=1e-9, atol=1e-9)
        checked += 1
    assert checked > 0


def test_cells_are_valid_meshio_cells(engine):
    """Order-2 MFEM cells become valid, positively oriented meshio++ cells."""
    for name in ("escher-p2", "fichera-q2", "tinyzoo-3d"):
        mesh = engine.read(MESHES / f"{name}.mesh")
        for cell_type, data in _blocks(mesh):
            if cell_type in (
                "tetra10",
                "hexahedron27",
                "wedge",
                "tetra",
                "hexahedron",
                "pyramid",
            ):
                for row in data:
                    x = mesh.points[row]
                    if cell_type in ("tetra10", "hexahedron27"):
                        # curved cells: check the orientation of the corners only
                        base = "tetra" if cell_type == "tetra10" else "hexahedron"
                        nc = _CORNERS[base]
                        assert _is_valid(x[:nc], base), (name, cell_type)
                    else:
                        assert _is_valid(x, cell_type), (name, cell_type)


def test_cell_data_and_attributes(engine):
    mesh = engine.read(MESHES / "compass.mesh", _gfs("compass"))
    assert [t for t, _ in _blocks(mesh)] == ["triangle6", "quad9", "line3"]
    e = mesh.cell_data["e"]
    assert np.isnan(e[2]).all() and not np.isnan(e[0]).any()
    attrs = [list(a) for a in mesh.cell_data["mfem:attribute"]]
    assert attrs[2] == [1, 2, 3, 4, 5, 6, 7, 8]
    regions = {(r.name, r.dim): r for r in mesh.regions}
    assert ("attribute_9", 2) in regions and ("boundary_1", 1) in regions
    # attribute sets are regions of their own name
    np.testing.assert_array_equal(
        regions[("East", 2)].entries,
        np.sort(
            np.concatenate(
                [
                    regions[("attribute_16", 2)].entries,
                    regions[("attribute_17", 2)].entries,
                ]
            )
        ),
    )
    assert ("Eastern Boundary", 1) in regions


def test_discontinuous_nodes_give_each_element_its_own_points(engine):
    mesh = engine.read(MESHES / "periodic-square.mesh")
    assert [(t, len(d)) for t, d in _blocks(mesh)] == [("quad", 9)]
    assert len(mesh.points) == 36
    for row in mesh.cells[0].data:
        a, b = (
            mesh.points[row[1]] - mesh.points[row[0]],
            mesh.points[row[2]] - mesh.points[row[0]],
        )
        assert a[0] * b[1] - a[1] * b[0] > 0


def _high_order_gfs(name):
    gf = MESHES / f"{name}.u.gf"
    return {"u": str(gf)} if gf.exists() else None


@pytest.mark.parametrize("name", HIGH_ORDER)
def test_arbitrary_order_matches_mfem(engine, name):
    """Every node of every VTK Lagrange cell sits where MFEM's own high-order
    output puts it, in VTK order, and u has MFEM's value there (frozen in
    reference_lagrange.npz)."""
    mesh = engine.read(MESHES / f"{name}.mesh", _high_order_gfs(name))
    cells = {}
    for block in mesh.cells:
        if block.type.startswith("VTK_LAGRANGE"):
            cells.setdefault(block.type, []).extend(np.asarray(block.data))
    scale = max(1.0, np.abs(mesh.points).max())
    codes = sorted(
        {int(k.split(":")[1]) for k in LAGRANGE.files if k.startswith(name + ":")}
    )
    assert codes
    for code in codes:
        ref = LAGRANGE[f"{name}:{code}:points"]
        ref_u = (
            LAGRANGE[f"{name}:{code}:u"]
            if f"{name}:{code}:u" in LAGRANGE.files
            else None
        )
        mine = cells[_VTK_LAGRANGE[code]]
        assert len(mine) == len(ref)
        centres = np.array([mesh.points[row].mean(axis=0) for row in mine])
        for k, want in enumerate(ref):
            i = int(np.argmin(np.abs(centres - want.mean(axis=0)).max(axis=1)))
            got = mesh.points[mine[i]]
            if np.allclose(got, want, atol=1e-11 * scale):
                if ref_u is not None:
                    np.testing.assert_allclose(
                        mesh.point_data["u"][mine[i]], ref_u[k], atol=1e-11
                    )
                continue
            # MFEM re-marks triangles and tetrahedra on loading (reorders their
            # vertices), so its cell can list the same nodes in another order.
            assert code in (69, 71) and ref_u is None
            np.testing.assert_allclose(
                np.sort(got, axis=0), np.sort(want, axis=0), atol=1e-11 * scale
            )


@pytest.mark.parametrize("name", HIGH_ORDER)
def test_arbitrary_order_engines_agree(name):
    gfs = _high_order_gfs(name)
    a = meshioplusplus.mfem.read(MESHES / f"{name}.mesh", gfs)
    b = py_mfem.read(MESHES / f"{name}.mesh", gfs)
    np.testing.assert_allclose(a.points, b.points, rtol=0, atol=1e-13)
    assert [(t, x.tolist()) for t, x in _blocks(a)] == [
        (t, x.tolist()) for t, x in _blocks(b)
    ]
    if gfs:
        np.testing.assert_allclose(a.point_data["u"], b.point_data["u"], atol=1e-13)
    assert _regions(a) == _regions(b)


@pytest.mark.parametrize("writer", ["core", "python"])
@pytest.mark.parametrize(
    "name, order",
    # the hexahedra carry an order-5 field, so they are read (and written) at 5
    [
        ("curved-tet-p4", 4),
        ("curved-hex-p3", 5),
        ("curved-prism-p4", 4),
        ("fichera-q3", 3),
    ],
)
def test_lagrange_cells_write_as_h1_nodes(writer, name, order, tmp_path):
    """VTK Lagrange cells are written as order-p H1 (Gauss-Lobatto) nodes and
    read back to the same nodes; so is the field."""
    gfs = _high_order_gfs(name)
    mesh = meshioplusplus.mfem.read(MESHES / f"{name}.mesh", gfs)
    out = tmp_path / "out.mesh"
    if writer == "core":
        _core.mfem_write(str(out), mesh, True)
    else:
        py_mfem.write(out, mesh, grid_functions=True)
    assert f"FiniteElementCollection: H1_3D_P{order}" in out.read_text()
    back = meshioplusplus.mfem.read(
        out, {"u": str(tmp_path / "out.u.gf")} if gfs else None
    )
    np.testing.assert_allclose(back.points, mesh.points, atol=1e-12)
    for (t, x), (u, y) in zip(_blocks(mesh), _blocks(back)):
        assert t == u
        np.testing.assert_array_equal(x, y)
    if gfs:
        np.testing.assert_allclose(
            back.point_data["u"], mesh.point_data["u"], atol=1e-12
        )


def test_non_conforming_mesh_reads_as_its_leaves(engine, capfd):
    """An MFEM NC mesh: the leaves MFEM itself builds, with their attributes."""
    mesh = engine.read(MESHES / "amr-quad.mesh")
    assert "leaf element" in " ".join(capfd.readouterr().err.split())
    attrs = np.concatenate(mesh.cell_data["mfem:attribute"])
    got = {"elements": [], "boundary": []}
    g = 0
    for block in mesh.cells:
        part = "elements" if block.type == "quad" else "boundary"
        for row in np.asarray(block.data):
            got[part].append(
                (int(attrs[g]), sorted(map(tuple, mesh.points[row].tolist())))
            )
            g += 1
    for part in ("elements", "boundary"):
        want = [
            (int(a), [tuple(p) for p in c])
            for a, c in zip(
                LAGRANGE[f"amr-quad:{part}:attributes"],
                LAGRANGE[f"amr-quad:{part}:corners"].tolist(),
            )
        ]
        assert sorted(got[part]) == sorted(want)


def test_non_conforming_mesh_skips_grid_functions(engine, tmp_path, capfd):
    gf = tmp_path / "e.gf"
    gf.write_text(
        "FiniteElementSpace\nFiniteElementCollection: L2_2D_P0\nVDim: 1\nOrdering: 0\n\n1\n"
    )
    mesh = engine.read(MESHES / "amr-quad.mesh", {"e": str(gf)})
    assert "e" not in mesh.cell_data
    assert "space-filling-curve" in " ".join(capfd.readouterr().err.split())


@pytest.mark.parametrize("name", CONFORMING)
def test_engines_write_the_same_bytes(name, tmp_path):
    mesh = meshioplusplus.mfem.read(MESHES / f"{name}.mesh", _gfs(name))
    _core.mfem_write(str(tmp_path / "cpp.mesh"), mesh, True)
    py_mfem.write(tmp_path / "py.mesh", mesh, grid_functions=True)
    cpp = sorted(p.name for p in tmp_path.glob("cpp*"))
    py = sorted(p.name for p in tmp_path.glob("py*"))
    assert [n[3:] for n in cpp] == [n[2:] for n in py]
    for a, b in zip(cpp, py):
        assert (tmp_path / a).read_bytes() == (tmp_path / b).read_bytes(), a


@pytest.mark.parametrize(
    "name",
    ["star-q2", "escher-p2", "fichera-q2", "fichera-mixed-p2", "compass", "tinyzoo-3d"],
)
def test_round_trip(engine, name, tmp_path):
    gfs = _gfs(name)
    mesh = engine.read(MESHES / f"{name}.mesh", gfs)
    engine.write(tmp_path / "out.mesh", mesh, grid_functions=True)
    back = engine.read(
        tmp_path / "out.mesh",
        {g: str(tmp_path / f"out.{g}.gf") for g in gfs},
    )
    np.testing.assert_allclose(back.points, mesh.points, rtol=0, atol=0)
    for (t, x), (u, y) in zip(_blocks(mesh), _blocks(back)):
        assert t == u
        np.testing.assert_array_equal(x, y)
    assert _regions(back) == _regions(mesh)
    for g in gfs:
        if g in mesh.point_data:
            np.testing.assert_array_equal(back.point_data[g], mesh.point_data[g])
        else:
            for x, y in zip(mesh.cell_data[g], back.cell_data[g]):
                np.testing.assert_array_equal(x, y)


def test_writer_completes_serendipity_cells_and_drops_data(engine, tmp_path, capfd):
    """A hex20 from another format becomes an order-2 MFEM hex whose face and
    body centres follow the serendipity map; data is dropped with a warning."""
    corners = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )
    edges = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ]
    mids = np.array([(corners[a] + corners[b]) / 2 for a, b in edges])
    mids[0, 1] -= 0.1  # a curved bottom edge
    points = np.vstack([corners, mids])
    mesh = meshioplusplus.Mesh(points, [("hexahedron20", [list(range(20))])])
    mesh.point_data["T"] = np.arange(20.0)
    engine.write(tmp_path / "h.mesh", mesh)
    assert "dropped" in " ".join(capfd.readouterr().err.split())
    back = engine.read(tmp_path / "h.mesh")
    assert [t for t, _ in _blocks(back)] == ["hexahedron27"]
    row = back.cells[0].data[0]
    x = back.points[row]
    # the bottom face centre: -1/4 corners + 1/2 mid-edges
    want = -0.25 * corners[[0, 1, 2, 3]].sum(0) + 0.5 * mids[[0, 1, 2, 3]].sum(0)
    np.testing.assert_allclose(x[24], want)
    np.testing.assert_allclose(x[8], mids[0])


def test_side_regions_become_boundary_elements(engine, tmp_path):
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float
    )
    mesh = meshioplusplus.Mesh(points, [("tetra", [[0, 1, 2, 3], [1, 2, 3, 4]])])
    mesh.regions = [meshioplusplus.Region("floor", "side", [[0, 3]], dim=2, tag=5)]
    engine.write(tmp_path / "t.mesh", mesh)
    text = (tmp_path / "t.mesh").read_text()
    assert "MFEM mesh v1.3" in text and '"floor" 1 5' in text
    back = engine.read(tmp_path / "t.mesh")
    regions = {(r.name, r.dim): r for r in back.regions}
    assert ("boundary_5", 2) in regions and ("floor", 2) in regions
    assert [t for t, _ in _blocks(back)] == ["tetra", "triangle"]


def test_mesh_extension_is_shared_with_medit(tmp_path):
    assert meshioplusplus._helpers._filetypes_from_path(pathlib.Path("x.mesh")) == [
        "medit",
        "mfem",
    ]
    mesh = meshioplusplus.read(MESHES / "star-q2.mesh")
    assert mesh.cells[0].type == "quad9"
    assert meshioplusplus.sniff_format(MESHES / "star-q2.mesh") == "mfem"
    assert meshioplusplus.sniff_format(MESHES / "amr-quad.mesh") == "mfem"
    # The registry's content-aware default (native CLI, C, WASM ...): the
    # native metadata path resolves the format itself, and medit would refuse.
    meta = _core.read_metadata(str(MESHES / "star-q2.mesh"))
    assert "quad9" in str(meta)
    medit = tmp_path / "m.mesh"
    meshioplusplus.write(
        medit,
        meshioplusplus.Mesh(
            [[0.0, 0, 0], [1, 0, 0], [0, 1, 0]], [("triangle", [[0, 1, 2]])]
        ),
        file_format="medit",
    )
    assert "triangle" in str(_core.read_metadata(str(medit)))


@pytest.mark.parametrize(
    "body, match",
    [
        (
            "MFEM mesh v1.0\ndimension\n2\nelements\n1\n1 3 0 1 2 9\nboundary\n0\nvertices\n3\n2\n0 0\n1 0\n0 1\n",
            "out of range",
        ),
        ("MFEM mesh v1.0\ndimension\n2\nelements\n1\n1 9 0 1 2\n", "unknown geometry"),
        ("MFEM mesh v1.0\ndimension\n2\nelements\n1\n1 2 0 1\n", "ends"),
        ("MFEM mesh v7.0\n", "not an MFEM mesh"),
        (
            "MFEM NC mesh v1.0\ndimension\n2\nelements\n0\nboundary\n0\n",
            "no top-level coordinates",
        ),
        (
            "MFEM NC mesh v1.0\ndimension\n2\nelements\n0\nboundary\n0\nnodes\n",
            "curved non-conforming",
        ),
        ("MFEM NC mesh v2.0\n", "not supported"),
        ("MFEM NURBS mesh v1.0\n", "not supported"),
        (
            "MFEM mesh v1.0\ndimension\n2\nelements\n0\nboundary\n0\nvertices\n0\n\nnodes\nFiniteElementSpace\nFiniteElementCollection: RT_2D_P1\nVDim: 2\nOrdering: 0\n",
            "not supported",
        ),
        (
            "MFEM mesh v1.0\ndimension\n2\nelements\n0\nboundary\n0\nbogus\n",
            "unexpected 'bogus'",
        ),
    ],
)
def test_errors_name_the_culprit(engine, tmp_path, body, match):
    path = tmp_path / "bad.mesh"
    path.write_text(body)
    with pytest.raises(meshioplusplus.ReadError, match=match):
        engine.read(path)


def test_grid_function_mismatch_is_an_error(engine, tmp_path):
    bad = tmp_path / "bad.gf"
    bad.write_text(
        "FiniteElementSpace\nFiniteElementCollection: H1_2D_P1\nVDim: 1\nOrdering: 0\n\n1\n2\n"
    )
    with pytest.raises(meshioplusplus.ReadError, match="dofs"):
        engine.read(MESHES / "compass.mesh", {"x": str(bad)})
