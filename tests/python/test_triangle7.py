"""``triangle7`` (v16.0.0): a ``triangle6`` plus a centre node, as Code_Aster
(``TRIA7``), MED (``TR7``) and VTK (34) define it. Before it was a cell type the
VTK and MED readers named blocks ``triangle7`` that nothing else could size."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus._common import num_nodes_per_cell
from meshioplusplus.med import _med as py_med
from meshioplusplus.vtu import _vtu as py_vtu

POINTS = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [1 / 3, 1 / 3, 0.0],
    ]
)


def _mesh():
    return meshioplusplus.Mesh(POINTS, [("triangle7", [list(range(7))])])


def test_is_a_seven_node_cell():
    assert num_nodes_per_cell["triangle7"] == 7
    assert _mesh().cells[0].dim == 2


@pytest.mark.parametrize(
    "fmt, core, python",
    [
        ("vtu", meshioplusplus.vtu, py_vtu),
        ("med", meshioplusplus.med, py_med),
    ],
)
@pytest.mark.parametrize("engine", ["core", "python"])
def test_round_trips(tmp_path, fmt, core, python, engine):
    if fmt == "med":
        pytest.importorskip("h5py")
    module = core if engine == "core" else python
    path = tmp_path / f"t7.{fmt}"
    module.write(path, _mesh())
    back = meshioplusplus.read(path)
    assert [b.type for b in back.cells] == ["triangle7"]
    np.testing.assert_array_equal(back.cells[0].data, [list(range(7))])
    np.testing.assert_allclose(back.points[:, :2], POINTS[:, :2])


def test_linearizes_to_a_triangle():
    out = meshioplusplus.convert_cells(_mesh(), mode="linearize")
    assert [b.type for b in out.cells] == ["triangle"]
    assert len(out.points) == 3


def test_surface_is_three_quadratic_edges():
    surface = meshioplusplus.extract_surface(_mesh())
    assert [(b.type, len(b.data)) for b in surface.cells] == [("line3", 3)]
