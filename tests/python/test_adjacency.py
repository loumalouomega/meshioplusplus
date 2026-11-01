"""`adjacency` and `write_graph` (roadmap 7.2.2): the node, clique and cell
graphs in CSR form, their numpy twins, and the METIS / Scotch graph files."""

import copy
import warnings

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _adjacency

from . import helpers

try:
    from meshioplusplus import _core

    HAS_CORE = hasattr(_core, "adjacency")
except ImportError:  # pragma: no cover - the numpy-only install
    HAS_CORE = False

needs_core = pytest.mark.skipif(not HAS_CORE, reason="needs the C++ core")

MESHES = {
    "empty": helpers.empty_mesh,
    "line": helpers.line_mesh,
    "tri": helpers.tri_mesh,
    "tri_2d": helpers.tri_mesh_2d,
    "line_tri": helpers.line_tri_mesh,
    "triangle6": helpers.triangle6_mesh,
    "quad": helpers.quad_mesh,
    "tri_quad": helpers.tri_quad_mesh,
    "tet": helpers.tet_mesh,
    "tet10": helpers.tet10_mesh,
    "hex": helpers.hex_mesh,
    "wedge": helpers.wedge_mesh,
    "pyramid": helpers.pyramid_mesh,
    "hex20": helpers.hex20_mesh,
    "polygon": helpers.polygon_mesh,
    "polygon2": helpers.polygon2_mesh,
}


def _five_tets():
    pts = np.array(
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
    tets = np.array(
        [[0, 1, 3, 4], [1, 2, 3, 6], [1, 4, 5, 6], [3, 4, 6, 7], [1, 3, 4, 6]]
    )
    return meshioplusplus.Mesh(pts, [("tetra", tets)])


def _assert_csr(indptr, indices, n):
    assert indptr.dtype == np.int64 and indices.dtype == np.int64
    assert indptr.shape == (n + 1,)
    assert indptr[0] == 0 and indptr[-1] == len(indices)
    rows = np.repeat(np.arange(n), np.diff(indptr))
    assert not np.any(rows == indices)
    for i in range(n):
        row = indices[indptr[i] : indptr[i + 1]]
        assert np.all(np.diff(row) > 0)
    pairs = set(zip(rows.tolist(), indices.tolist()))
    assert pairs == {(b, a) for a, b in pairs}


def test_probe_node_graph_counts_each_edge_twice():
    indptr, indices = meshioplusplus.adjacency(_five_tets(), "node")
    _assert_csr(indptr, indices, 8)
    assert indptr[-1] == 2 * 18  # 12 cube edges + 6 face diagonals
    assert indices[indptr[0] : indptr[1]].tolist() == [1, 3, 4]


def test_cell_graph_of_five_tets():
    indptr, indices = meshioplusplus.adjacency(_five_tets(), "cell")
    _assert_csr(indptr, indices, 5)
    assert indices[indptr[4] : indptr[5]].tolist() == [0, 1, 2, 3]
    for c in range(4):
        assert indices[indptr[c] : indptr[c + 1]].tolist() == [4]


def test_clique_is_the_fem_pattern():
    mesh = helpers.hex_mesh
    indptr, indices = meshioplusplus.adjacency(mesh, "clique")
    n = len(mesh.points)
    _assert_csr(indptr, indices, n)
    # Every pair of the hexahedron's nodes: K8.
    assert indptr[-1] == 8 * 7
    node_ptr, node_idx = meshioplusplus.adjacency(mesh, "node")
    assert node_ptr[-1] == 2 * 12
    for i in range(n):
        assert set(node_idx[node_ptr[i] : node_ptr[i + 1]]) <= set(
            indices[indptr[i] : indptr[i + 1]]
        )


@pytest.mark.parametrize("name", sorted(MESHES))
def test_node_graph_matches_edge_index(name):
    """The CSR node graph is `edge_index`'s coordinate graph, row by row."""
    mesh = MESHES[name]
    indptr, indices = meshioplusplus.adjacency(mesh, "node")
    _assert_csr(indptr, indices, len(mesh.points))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # edge_index warns on untabulated types
        ei = meshioplusplus.edge_index(mesh, "node")
    rows = np.repeat(np.arange(len(mesh.points)), np.diff(indptr))
    np.testing.assert_array_equal(np.vstack([rows, indices]), ei)


@needs_core
@pytest.mark.parametrize("kind", ["node", "clique", "cell"])
@pytest.mark.parametrize("name", sorted(MESHES) + ["five_tets"])
def test_core_matches_numpy_twin(name, kind):
    mesh = _five_tets() if name == "five_tets" else MESHES[name]
    core = _core.adjacency(mesh, kind)
    twin = _adjacency._adjacency_py(mesh, kind)
    np.testing.assert_array_equal(core[0], twin[0])
    np.testing.assert_array_equal(core[1], twin[1])


@needs_core
def test_polyhedra_in_the_core():
    mesh = helpers.polyhedron_mesh
    for kind in ("node", "clique", "cell"):
        indptr, indices = meshioplusplus.adjacency(mesh, kind)
        n = len(mesh.points) if kind != "cell" else sum(len(b.data) for b in mesh.cells)
        _assert_csr(indptr, indices, n)
    # The node and clique twins handle polyhedra too; only the cell twin does not.
    for kind in ("node", "clique"):
        core = _core.adjacency(mesh, kind)
        twin = _adjacency._adjacency_py(mesh, kind)
        np.testing.assert_array_equal(core[0], twin[0])
        np.testing.assert_array_equal(core[1], twin[1])
    with pytest.raises(NotImplementedError, match="C\\+\\+ core"):
        _adjacency._adjacency_py(mesh, "cell")


def test_uniform_polygon_block_joins_by_rings():
    """A uniform n-gon block (rectangular storage) is rings, not skipped."""
    pts = np.array(
        [[0, 0], [1, 0], [1, 1], [0, 1], [-0.5, 0.5], [2, 0], [2.5, 0.5], [2, 1]],
        dtype=float,
    )
    mesh = meshioplusplus.Mesh(
        pts, [("polygon", np.array([[0, 1, 2, 3, 4], [1, 5, 6, 7, 2]]))]
    )
    indptr, indices = meshioplusplus.adjacency(mesh, "node")
    assert indptr[-1] == 2 * 9
    assert indices[indptr[1] : indptr[2]].tolist() == [0, 2, 5]
    cptr, cidx = meshioplusplus.adjacency(mesh, "cell")
    assert cptr.tolist() == [0, 1, 2] and cidx.tolist() == [1, 0]


def test_mixed_dimensions_isolate_lower_cells():
    mesh = helpers.line_tri_mesh
    indptr, _ = meshioplusplus.adjacency(mesh, "cell")
    nlines = len(mesh.cells[0].data)
    # Lines come first (block-major) and have no neighbours in the 2-D dual.
    assert np.all(np.diff(indptr)[:nlines] == 0)


def test_input_is_not_modified():
    mesh = copy.deepcopy(helpers.tet_mesh)
    before = [b.data.copy() for b in mesh.cells]
    for kind in ("node", "clique", "cell"):
        meshioplusplus.adjacency(mesh, kind)
    for a, b in zip(before, mesh.cells):
        np.testing.assert_array_equal(a, b.data)


def test_bad_kind():
    with pytest.raises(ValueError, match="'node', 'clique' or 'cell'"):
        meshioplusplus.adjacency(helpers.tri_mesh, "edge")


# --------------------------------------------------------------------------- #
# Graph files                                                                 #
# --------------------------------------------------------------------------- #
METIS_TRI = "4 5\n2 3 4\n1 3\n1 2 4\n1 3\n"
SCOTCH_TRI = "0\n4 10\n0 000\n3 1 2 3\n2 0 2\n3 0 1 3\n2 0 2\n"


def _two_triangles():
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    return meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))])


@pytest.mark.parametrize(
    "ext,fmt,text",
    [
        (".graph", None, METIS_TRI),
        (".metis", None, METIS_TRI),
        (".grf", None, SCOTCH_TRI),
        (".txt", "scotch", SCOTCH_TRI),
        (".grf", "metis", METIS_TRI),
    ],
)
def test_write_graph_files(tmp_path, ext, fmt, text):
    path = tmp_path / f"g{ext}"
    report = meshioplusplus.write_graph(path, _two_triangles(), file_format=fmt)
    assert path.read_text() == text
    assert report == {
        "format": fmt or ("scotch" if ext == ".grf" else "metis"),
        "num_vertices": 4,
        "num_edges": 5,
        "max_degree": 3,
    }


@pytest.mark.parametrize("fmt", ["metis", "scotch"])
@pytest.mark.parametrize("kind", ["node", "clique", "cell"])
def test_writers_byte_identical_between_engines(tmp_path, fmt, kind):
    mesh = _five_tets()
    indptr, indices = meshioplusplus.adjacency(mesh, kind)
    a, b = tmp_path / "core.out", tmp_path / "twin.out"
    rep_a = meshioplusplus.write_graph(a, (indptr, indices), file_format=fmt)
    rep_b = _adjacency._write_graph_py(str(b), indptr, indices, fmt)
    assert a.read_bytes() == b.read_bytes()
    assert rep_a == rep_b


def _parse_metis(text):
    lines = text.split("\n")
    n, m = map(int, lines[0].split())
    rows = [[int(v) - 1 for v in line.split()] for line in lines[1 : n + 1]]
    indptr = np.cumsum([0] + [len(r) for r in rows])
    return n, m, indptr, np.array([v for r in rows for v in r], dtype=np.int64)


def test_metis_round_trip(tmp_path):
    mesh = helpers.hex20_mesh
    for kind in ("clique", "node"):
        indptr, indices = meshioplusplus.adjacency(mesh, kind)
        path = tmp_path / f"{kind}.graph"
        meshioplusplus.write_graph(path, mesh, kind)
        n, m, ptr, idx = _parse_metis(path.read_text())
        assert n == len(mesh.points) and m == len(indices) // 2
        np.testing.assert_array_equal(ptr, indptr)
        np.testing.assert_array_equal(idx, indices)


def test_isolated_vertices(tmp_path):
    mesh = _two_triangles()
    mesh.points = np.vstack([mesh.points, [[5.0, 5.0, 5.0]]])
    meshioplusplus.write_graph(tmp_path / "g.graph", mesh)
    assert (tmp_path / "g.graph").read_text() == "5 5\n2 3 4\n1 3\n1 2 4\n1 3\n\n"
    meshioplusplus.write_graph(tmp_path / "g.grf", mesh)
    assert (tmp_path / "g.grf").read_text().endswith("2 0 2\n0\n")


@pytest.mark.parametrize("engine", ["core", "twin"])
@pytest.mark.parametrize(
    "indptr,indices,match",
    [
        ([0, 1, 1], [1], "not symmetric"),
        ([0, 1], [0], "self-loop"),
        ([0, 1, 2], [1, 5], "outside"),
        ([0, 2, 1], [1, 1], "indptr"),
        ([], [], "empty"),
        ([0, 2, 3], [1, 1, 0], "twice"),
    ],
)
def test_writer_refuses_malformed_graphs(tmp_path, engine, indptr, indices, match):
    if engine == "core" and not HAS_CORE:
        pytest.skip("needs the C++ core")
    ptr = np.asarray(indptr, dtype=np.int64)
    idx = np.asarray(indices, dtype=np.int64)
    path = str(tmp_path / "bad.graph")
    with pytest.raises(ValueError, match=match):
        if engine == "core":
            meshioplusplus.write_graph(path, (ptr, idx))
        else:
            _adjacency._write_graph_py(path, ptr, idx, "metis")


def test_write_graph_argument_errors(tmp_path):
    with pytest.raises(ValueError, match="cannot infer the format"):
        meshioplusplus.write_graph(tmp_path / "g.txt", _two_triangles())
    with pytest.raises(ValueError, match="unknown format 'chaco'"):
        meshioplusplus.write_graph(
            tmp_path / "g.graph", _two_triangles(), file_format="chaco"
        )
    with pytest.raises(TypeError, match="Mesh or an"):
        meshioplusplus.write_graph(tmp_path / "g.graph", 42)


def test_exported():
    assert {"adjacency", "write_graph"} <= set(meshioplusplus.__all__)
