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


NC_REF = np.load(MESHES / "reference_nc.npz")


@pytest.mark.parametrize(
    "name, cell_type", [("amr-quad", "quad"), ("amr-hex", "hexahedron")]
)
def test_non_conforming_mesh_follows_mfem_numbering(engine, name, cell_type):
    """The leaves come in MFEM's space-filling-curve order and the points in
    MFEM's vertex order, so MFEM's grid functions on the mesh apply: an H1
    field at every VTK Lagrange node and an L2 one per leaf (frozen in
    reference_nc.npz)."""
    gfs = {g: str(MESHES / f"{name}.{g}.gf") for g in ("u", "e")}
    mesh = engine.read(MESHES / f"{name}.mesh", gfs)
    order = int(NC_REF[f"{name}:order"])
    want_type = {2: f"{cell_type}{9 if cell_type == 'quad' else 27}"}.get(
        order, _VTK_LAGRANGE[70 if cell_type == "quad" else 72]
    )
    block = next(b for b in mesh.cells if b.type == want_type)
    cells = np.asarray(block.data)
    verts = NC_REF[f"{name}:element_vertices"]
    np.testing.assert_array_equal(cells[:, : verts.shape[1]], verts)
    np.testing.assert_array_equal(
        mesh.points[: len(NC_REF[f"{name}:vertices"])], NC_REF[f"{name}:vertices"]
    )
    np.testing.assert_allclose(
        mesh.point_data["u"][cells], NC_REF[f"{name}:u"], rtol=0, atol=1e-13
    )
    k = [b.type for b in mesh.cells].index(want_type)
    np.testing.assert_allclose(mesh.cell_data["e"][k], NC_REF[f"{name}:e"], atol=1e-14)


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
        ("MFEM NURBS NC-patch mesh v1.0\n", "not supported"),
        (
            "MFEM NURBS mesh v1.0\ndimension\n2\nelements\n0\nbogus\n",
            "expected 'boundary'",
        ),
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


# --- parallel meshes ------------------------------------------------------------------

PARALLEL = MESHES / "parallel"
PARALLEL_REFERENCE = np.load(PARALLEL / "reference.npz")
# case -> (ranks, merged point count of the conforming mesh)
PARALLEL_CASES = {"star-p2": (4, 361), "beam-tet": (3, 153)}
# MFEM's VTK cell type -> ours (it writes order-1 and order-2 cells as VTK
# Lagrange too)
_PARALLEL_TYPES = {
    10: ["tetra"],
    28: ["quad9"],
    70: ["quad9"],
    71: ["tetra"],
}


@pytest.mark.parametrize("layout", ["mesh", "pmesh"])
@pytest.mark.parametrize("name", sorted(PARALLEL_CASES))
def test_parallel_mesh_merges_its_ranks(engine, name, layout):
    """Every rank file of a parallel mesh, ParMesh::Save's (no communication
    groups: boundary vertices merged by position, the ranks' interface faces
    dropped) or ParPrint's (merged by group), reads as one conforming mesh whose
    cells are MFEM's own, with each cell's rank."""
    ranks, npoints = PARALLEL_CASES[name]
    mesh = engine.read(
        PARALLEL / f"{name}.{layout}.000001", {"u": str(PARALLEL / f"{name}.u.000000")}
    )
    assert len(mesh.points) == npoints
    labels = np.concatenate(mesh.cell_data["partition:part"])
    assert sorted(set(labels.tolist())) == list(range(ranks))
    cells = {}
    for block in mesh.cells:
        cells.setdefault(block.type, []).extend(np.asarray(block.data))
    codes = {
        int(k.split(":")[1])
        for k in PARALLEL_REFERENCE.files
        if k.startswith(name + ":")
    }
    for code in codes:
        ref = PARALLEL_REFERENCE[f"{name}:{code}:points"]
        ref_u = PARALLEL_REFERENCE[f"{name}:{code}:u"]
        mine = [row for t in _PARALLEL_TYPES[code] for row in cells[t]]
        assert len(mine) == len(ref)
        dim = mesh.points.shape[1]
        centres = np.array([mesh.points[row].mean(axis=0) for row in mine])
        for k, want in enumerate(ref):
            want = want[:, :dim]
            i = int(np.argmin(np.abs(centres - want.mean(axis=0)).max(axis=1)))
            np.testing.assert_allclose(mesh.points[mine[i]], want, atol=1e-12)
            np.testing.assert_allclose(
                mesh.point_data["u"][mine[i]], ref_u[k], atol=1e-12
            )
    # the interface faces a ParMesh::Save rank lists as boundary are gone
    boundary = sum(len(b.data) for b in mesh.cells if b.dim < mesh.cells[0].dim)
    # (serial MFEM: star 40, beam-tet 272 boundary elements once refined)
    serial_boundary = {"star-p2": 40, "beam-tet": 272}[name]
    assert boundary == serial_boundary


def test_parallel_engines_agree():
    for name in PARALLEL_CASES:
        a = meshioplusplus.mfem.read(PARALLEL / f"{name}.mesh.000000")
        b = py_mfem.read(PARALLEL / f"{name}.mesh.000000")
        np.testing.assert_allclose(a.points, b.points, atol=1e-13)
        assert [(t, x.tolist()) for t, x in _blocks(a)] == [
            (t, x.tolist()) for t, x in _blocks(b)
        ]


def test_parallel_piece_reads_one_rank(engine):
    whole = engine.read(PARALLEL / "star-p2.mesh.000000")
    piece = engine.read(PARALLEL / "star-p2.mesh.000000", piece=2)
    labels = np.concatenate(piece.cell_data["partition:part"])
    assert set(labels.tolist()) == {2}
    merged = np.concatenate(whole.cell_data["partition:part"])
    n2 = int((merged[: len(whole.cells[0].data)] == 2).sum())
    assert len(piece.cells[0].data) == n2
    with pytest.raises(meshioplusplus.ReadError, match="out of range"):
        engine.read(PARALLEL / "star-p2.mesh.000000", piece=4)
    # through the generic API too
    same = meshioplusplus.read(
        PARALLEL / "star-p2.mesh.000000", file_format="mfem", piece=2
    )
    assert len(same.cells[0].data) == n2


def test_parallel_grid_function_must_be_a_rank_file(engine, tmp_path):
    gf = tmp_path / "u.gf"
    gf.write_text(
        "FiniteElementSpace\nFiniteElementCollection: H1_2D_P2\nVDim: 1\nOrdering: 0\n\n"
    )
    with pytest.raises(meshioplusplus.ReadError, match="rank files"):
        engine.read(PARALLEL / "star-p2.mesh.000000", {"u": str(gf)})


# --- NURBS meshes ----------------------------------------------------------------------------

NURBS = [
    "ball-nurbs",
    "pipe-nurbs",
    "square-disc-nurbs-patch",
    "nurbs-segments2d-patches",
    "beam-quad-nurbs-sf",
    "cube-nurbs",
]
NURBS_REF = np.load(MESHES / "nurbs" / "reference_nurbs.npz")


def _nurbs(engine, name):
    path = MESHES / "nurbs" / f"{name}.mesh"
    return engine.read(path, {"u": str(MESHES / "nurbs" / f"{name}.u.gf")})


@pytest.mark.parametrize("name", NURBS)
def test_nurbs_matches_mfem(engine, name):
    """Every knot-span element is MFEM's, in MFEM's order; every node of its
    cell sits where MFEM's element transformation puts it and u has MFEM's
    value there (frozen in nurbs/reference_nurbs.npz)."""
    mesh = _nurbs(engine, name)
    q = int(NURBS_REF[f"{name}:order"])
    want_x, want_u = NURBS_REF[f"{name}:x"], NURBS_REF[f"{name}:u"]
    dim = {2: 1, 4: 2, 8: 3}[NURBS_REF[f"{name}:elements"].shape[1]]
    block = mesh.cells[0]
    cells = np.asarray(block.data)
    assert len(cells) == len(want_x)
    shape = {1: "line", 2: "quad", 3: "hexahedron"}[dim]
    lex = [
        i + (q + 1) * (j + (q + 1) * k)
        for i, j, k in py_mfem._lagrange.vtk_lattice(shape, q)
    ]
    np.testing.assert_allclose(
        mesh.points[cells], want_x[:, lex, :], rtol=0, atol=1e-13 * np.abs(want_x).max()
    )
    np.testing.assert_allclose(
        mesh.point_data["u"][cells], want_u[:, lex], rtol=0, atol=1e-12
    )
    # MFEM's knot-span vertices are the cells' corners (the first nodes)
    elements = NURBS_REF[f"{name}:elements"]
    corners = cells[:, : elements.shape[1]]
    ids = {}
    for mine, theirs in zip(corners.ravel(), elements.ravel()):
        assert ids.setdefault(int(theirs), int(mine)) == int(mine)
    boundary = NURBS_REF[f"{name}:boundary"]
    bdr = [b for b in mesh.cells[1:]]
    got = np.concatenate([np.asarray(b.data)[:, : boundary.shape[1]] for b in bdr])
    np.testing.assert_array_equal(got, np.vectorize(ids.get)(boundary))


@pytest.mark.parametrize("name", NURBS)
def test_nurbs_engines_agree(name):
    _same(_nurbs(meshioplusplus.mfem, name), _nurbs(py_mfem, name))


def test_nurbs_cells_and_attributes(engine):
    mesh = _nurbs(engine, "pipe-nurbs")
    # order-2 NURBS hexahedra sampled as complete quadratic cells; the boundary
    # MFEM builds for a file without one, attribute 1
    assert [(b.type, len(b.data)) for b in mesh.cells] == [
        ("hexahedron27", 8),
        ("quad9", 24),
    ]
    assert [r.name for r in mesh.regions] == [
        f"attribute_{a}" for a in (1, 2, 3, 4)
    ] + ["boundary_1"]
    ball = _nurbs(engine, "ball-nurbs")
    assert [b.type for b in ball.cells] == [
        "VTK_LAGRANGE_HEXAHEDRON",
        "VTK_LAGRANGE_QUADRILATERAL",
    ]


def test_nurbs_skips_other_grid_functions(engine, tmp_path, capfd):
    gf = tmp_path / "p.gf"
    gf.write_text(
        "FiniteElementSpace\nFiniteElementCollection: H1_2D_P1\nVDim: 1\nOrdering: 0\n\n1\n"
    )
    mesh = engine.read(MESHES / "nurbs" / "cube-nurbs.mesh", {"p": str(gf)})
    assert "p" not in mesh.point_data
    assert "not the NURBS mesh's own" in capfd.readouterr().err


# --- Bernstein and serendipity spaces ------------------------------------------------------

MODAL = [
    "pos-quad-p3",
    "pos-tri-p4",
    "pos-tet-p3",
    "pos-hex-p3",
    "pos-wedge-p3",
    "pos-seg-p4",
    "ser-quad-p3",
    "ser-quad-p5",
]
MODAL_REF = np.load(MESHES / "modal" / "reference_modal.npz")
_MODAL_TYPES = {
    1: ("line3", "VTK_LAGRANGE_CURVE"),
    2: ("triangle6", "VTK_LAGRANGE_TRIANGLE"),
    3: ("quad9", "VTK_LAGRANGE_QUADRILATERAL"),
    4: ("tetra10", "VTK_LAGRANGE_TETRAHEDRON"),
    5: ("hexahedron27", "VTK_LAGRANGE_HEXAHEDRON"),
    6: ("wedge18", "VTK_LAGRANGE_WEDGE"),
}


def _modal(engine, name):
    return engine.read(
        MESHES / "modal" / f"{name}.mesh",
        {"u": str(MESHES / "modal" / f"{name}.u.gf")},
    )


@pytest.mark.parametrize("name", MODAL)
def test_modal_spaces_match_mfem(engine, name):
    """Bernstein (H1Pos) and serendipity (H1Ser) nodes and fields are
    coefficients, not values: every VTK Lagrange node of every cell sits where
    MFEM's element transformation puts it, and u has MFEM's value there
    (frozen in modal/reference_modal.npz)."""
    mesh = _modal(engine, name)
    q = int(MODAL_REF[f"{name}:order"])
    geoms = sorted(
        {
            int(k.split(":")[1])
            for k in MODAL_REF.files
            if k.startswith(name + ":") and k.count(":") == 2
        }
    )
    for geom in geoms:
        want_x = MODAL_REF[f"{name}:{geom}:x"]
        want_u = MODAL_REF[f"{name}:{geom}:u"]
        cell_type = _MODAL_TYPES[geom][0 if q == 2 else 1]
        cells = np.asarray(next(b.data for b in mesh.cells if b.type == cell_type))
        assert len(cells) == len(want_x)
        np.testing.assert_allclose(
            mesh.points[cells], want_x, rtol=0, atol=1e-13 * np.abs(want_x).max()
        )
        np.testing.assert_allclose(
            mesh.point_data["u"][cells], want_u, rtol=0, atol=1e-12
        )


@pytest.mark.parametrize("name", MODAL)
def test_modal_spaces_engines_agree(name):
    a, b = _modal(meshioplusplus.mfem, name), _modal(py_mfem, name)
    np.testing.assert_allclose(a.points, b.points, rtol=0, atol=1e-13)
    assert [(t, x.tolist()) for t, x in _blocks(a)] == [
        (t, x.tolist()) for t, x in _blocks(b)
    ]
    np.testing.assert_allclose(a.point_data["u"], b.point_data["u"], atol=1e-13)


def test_serendipity_needs_quadrilaterals(engine, capfd):
    # a serendipity field on a mesh of other cells is skipped
    mesh = engine.read(
        MESHES / "modal" / "pos-tri-p4.mesh",
        {"u": str(MESHES / "modal" / "ser-quad-p3.u.gf")},
    )
    assert "u" not in mesh.point_data
    assert "quadrilateral meshes only" in capfd.readouterr().err
