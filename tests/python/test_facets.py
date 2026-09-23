"""The Python twin of detail/facet_index (tests/cpp/test_facet_index.cpp)."""

import numpy as np

import meshioplusplus
from meshioplusplus._facets import FacetIndex, facet_nodes


def test_two_tets_share_one_face():
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float),
        [("tetra", np.array([[0, 1, 2, 3], [1, 2, 3, 4]]))],
    )
    index = FacetIndex(mesh)
    assert len(index) == 7
    hit = index.find([3, 1, 2])
    assert hit.count == 2
    assert hit.first == (0, 1)
    assert hit.second[0] == 1
    assert index.find([0, 1, 3]).first == (0, 0)
    assert index.find([0, 1, 4]) is None


def test_surface_cells_index_their_edges_and_optionally_themselves():
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float),
        [("quad", np.array([[0, 1, 2, 3]]))],
    )
    edges = FacetIndex(mesh)
    assert edges.find([2, 1]).first == (0, 1)
    assert edges.find([3, 2, 1, 0]) is None
    self_ = FacetIndex(mesh, surface_edges=False, surface_self=True)
    assert self_.find([3, 2, 1, 0]).first == (0, 0)
    assert self_.find([2, 1]) is None


def test_facet_nodes_includes_mid_side_nodes():
    mesh = meshioplusplus.Mesh(
        np.zeros((10, 3)), [("tetra10", np.arange(10).reshape(1, 10))]
    )
    assert facet_nodes(mesh, 0, 0) == ("triangle6", [0, 1, 3, 4, 8, 7])
    assert facet_nodes(mesh, 0, 4) is None
    assert facet_nodes(mesh, 1, 0) is None
