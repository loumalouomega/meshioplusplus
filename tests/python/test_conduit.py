"""Conduit Mesh Blueprint in memory (roadmap 10.1).

Every test but the last runs on a plain nested dict of the Blueprint's paths,
which ``to_conduit``/``from_conduit`` accept as they accept a ``conduit.Node``;
the last checks a real node with ``conduit.blueprint.mesh.verify`` when Conduit's
Python module is installed (it is not on PyPI).
"""

import numpy as np
import pytest

import meshioplusplus as mio


def _mixed():
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [2.0, 0.0, 0.0],
        ]
    )
    return mio.Mesh(
        points,
        [
            ("tetra", np.array([[0, 1, 2, 4]])),
            ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
            ("line", np.array([[1, 5]])),
            ("vertex", np.array([[5]])),
        ],
        point_data={"T": np.arange(6.0), "u": np.arange(18.0).reshape(6, 3)},
        cell_data={
            "p": [
                np.array([1.0]),
                np.array([2.0, 3.0]),
                np.array([4.0]),
                np.array([5.0]),
            ],
            "id": [np.array([7]), np.array([8, 9]), np.array([10]), np.array([11])],
        },
        field_data={"meshio:time": np.array([0.25])},
    )


def _same(a, b):
    np.testing.assert_array_equal(np.asarray(a.points), np.asarray(b.points))
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        if x.type.startswith("poly"):
            assert [np.asarray(r).tolist() for r in x.data] == [
                np.asarray(r).tolist() for r in y.data
            ] or [[np.asarray(f).tolist() for f in r] for r in x.data] == [
                [np.asarray(f).tolist() for f in r] for r in y.data
            ]
        else:
            np.testing.assert_array_equal(np.asarray(x.data), np.asarray(y.data))
    assert sorted(a.point_data) == sorted(b.point_data)
    for k in a.point_data:
        np.testing.assert_array_equal(a.point_data[k], b.point_data[k])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for k in a.cell_data:
        for x, y in zip(a.cell_data[k], b.cell_data[k]):
            np.testing.assert_array_equal(x, y)


def test_mixed_shapes_round_trip_exactly():
    mesh = _mixed()
    node = mio.to_conduit(mesh, {})
    assert node["coordsets"]["coords"]["type"] == "explicit"
    elements = node["topologies"]["mesh"]["elements"]
    assert elements["shape"] == "mixed"
    assert set(elements["shape_map"]) == {"tet", "tri", "line", "point"}
    assert list(elements["shapes"]) == [10, 5, 5, 3, 1]
    assert list(elements["offsets"]) == [0, 4, 7, 10, 12]
    assert node["fields"]["u"]["association"] == "vertex"
    assert set(node["fields"]["u"]["values"]) == {"x", "y", "z"}
    assert node["state"]["time"] == 0.25
    back = mio.from_conduit(node)
    _same(mesh, back)
    assert back.field_data["meshio:time"].tolist() == [0.25]


def test_one_shape_is_not_mixed():
    mesh = mio.Mesh(np.eye(3), [("triangle", np.array([[0, 1, 2]]))])
    node = mio.to_conduit(mesh, {})
    assert node["topologies"]["mesh"]["elements"]["shape"] == "tri"
    _same(mesh, mio.from_conduit(node))


def test_polygons_and_polyhedra_round_trip():
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 1.5, 0], [0, 0, 1]], float
    )
    polygons = mio.Mesh(
        points, [("polygon", [np.array([0, 1, 2, 4, 3]), np.array([0, 2, 3])])]
    )
    node = mio.to_conduit(polygons, {})
    assert node["topologies"]["mesh"]["elements"]["shape"] == "polygonal"
    back = mio.from_conduit(node)
    assert [list(r) for r in back.cells[0].data] == [[0, 1, 2, 4, 3], [0, 2, 3]]

    tet_faces = [[0, 2, 1], [0, 1, 5], [1, 2, 5], [2, 0, 5]]
    poly = mio.Mesh(points, [("polyhedron4", [[np.array(f) for f in tet_faces]])])
    node = mio.to_conduit(poly, {})
    topo = node["topologies"]["mesh"]
    assert topo["elements"]["shape"] == "polyhedral"
    assert topo["subelements"]["shape"] == "polygonal"
    back = mio.from_conduit(node)
    assert back.cells[0].type.startswith("polyhedron")
    assert [[list(f) for f in c] for c in back.cells[0].data] == [tet_faces]


def test_blocks_without_a_shape_are_dropped_with_a_warning(capsys):
    mesh = mio.Mesh(
        np.zeros((6, 3)),
        [
            ("triangle", np.array([[0, 1, 2]])),
            ("triangle6", np.array([[0, 1, 2, 3, 4, 5]])),
        ],
        cell_data={"p": [np.array([1.0]), np.array([2.0])]},
    )
    node = mio.to_conduit(mesh, {})
    assert "triangle6" in capsys.readouterr().err
    back = mio.from_conduit(node)
    assert [c.type for c in back.cells] == ["triangle"]
    assert back.cell_data["p"][0].tolist() == [1.0]


def test_uniform_and_structured_grids_expand_to_cells():
    node = {
        "coordsets": {
            "g": {
                "type": "uniform",
                "dims": {"i": 3, "j": 2},
                "origin": {"x": 0.0, "y": 0.0},
                "spacing": {"dx": 0.5, "dy": 2.0},
            }
        },
        "topologies": {"t": {"type": "uniform", "coordset": "g"}},
        "fields": {
            "c": {
                "association": "element",
                "topology": "t",
                "values": np.array([1.0, 2.0]),
            }
        },
    }
    mesh = mio.from_conduit(node)
    np.testing.assert_allclose(
        mesh.points, [[0, 0], [0.5, 0], [1, 0], [0, 2], [0.5, 2], [1, 2]]
    )
    assert mesh.cells[0].type == "quad"
    assert mesh.cells[0].data.tolist() == [[0, 1, 4, 3], [1, 2, 5, 4]]
    assert mesh.cell_data["c"][0].tolist() == [1.0, 2.0]

    node = {
        "coordsets": {
            "c": {
                "type": "explicit",
                "values": {
                    "x": np.array([0.0, 1, 0, 1, 0, 1, 0, 1]),
                    "y": np.array([0.0, 0, 1, 1, 0, 0, 1, 1]),
                    "z": np.array([0.0, 0, 0, 0, 1, 1, 1, 1]),
                },
            }
        },
        "topologies": {
            "s": {
                "type": "structured",
                "coordset": "c",
                "elements": {"dims": {"i": 1, "j": 1, "k": 1}},
            }
        },
    }
    mesh = mio.from_conduit(node)
    assert mesh.cells[0].type == "hexahedron"
    assert mesh.cells[0].data.tolist() == [[0, 1, 3, 2, 4, 5, 7, 6]]


def test_a_multi_domain_node_reads_one_domain(capsys):
    one = mio.to_conduit(_mixed(), {})
    two = mio.to_conduit(mio.Mesh(np.eye(3), [("triangle", np.array([[0, 1, 2]]))]), {})
    node = {"domain_000000": one, "domain_000001": two}
    first = mio.from_conduit(node)
    assert "2 domains" in capsys.readouterr().err
    assert len(first.points) == 6
    assert len(mio.from_conduit(node, domain=1).points) == 3


def test_missing_conduit_is_a_named_error(monkeypatch):
    import sys

    monkeypatch.setitem(sys.modules, "conduit", None)
    assert not mio.has_conduit()
    with pytest.raises(ImportError, match="Conduit's Python module"):
        mio.to_conduit(_mixed())


def test_a_real_node_verifies():
    conduit = pytest.importorskip("conduit")
    blueprint = pytest.importorskip("conduit.blueprint")
    if not hasattr(conduit, "Node"):
        pytest.skip("the 'conduit' module installed is not LLNL Conduit")
    mesh = _mixed()
    node = mio.to_conduit(mesh)
    info = conduit.Node()
    assert blueprint.mesh.verify(node, info), info.to_yaml()
    _same(mesh, mio.from_conduit(node))
