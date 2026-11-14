"""ADCIRC fort.14 and fort.63/fort.64 (roadmap 2.4.3.3), files written here."""

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _fallback
from meshioplusplus.adcirc import _adcirc

GRID = """quarter annulus test grid ! comment
2 4   = NE NP
1 0.0 0.0 5.0
2 1.0 0.0 6.0
3 1.0 1.0 7.0
10 0.0 1.0 8.0
1 3 1 2 3
2 3 1 3 10
1 = Number of open boundaries
2 = Total number of open boundary nodes
2 = Number of nodes for open boundary 1
2
3
2 = Number of land boundaries
5 = Total number of land boundary nodes
3 20 = Number of nodes for land boundary 1
3
10
1
2 3 = external barrier
1 1.5 2.0
2 1.5 2.0
"""

FORT63 = """run fort.63
2 4 3600.0 1 1
3600.0 1
1 0.1
2 0.2
3 0.3
10 0.4
7200.0 2 2 -99999.0
2 1.2
10 1.4
"""

FORT64 = """run fort.64
2 4 3600.0 1 2
3600.0 1
1 1.0 0.0
2 1.0 0.5
3 0.0 1.0
10 0.0 0.0
7200.0 2
1 2.0 0.0
2 2.0 0.5
3 0.0 2.0
10 0.0 0.0
"""


@pytest.fixture(params=["core", "python"])
def reader(request, monkeypatch):
    if request.param == "core":
        monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
        monkeypatch.setattr(_fallback, "_strict", True)
        return mio.adcirc.read
    return _adcirc.read


@pytest.fixture
def case(tmp_path):
    (tmp_path / "fort.14").write_text(GRID)
    (tmp_path / "fort.63").write_text(FORT63)
    (tmp_path / "fort.64").write_text(FORT64)
    return tmp_path


def test_grid_boundaries_and_first_dataset(case, reader):
    mesh = reader(case / "fort.14")
    assert np.asarray(mesh.cells[0].data).tolist() == [[0, 1, 2], [0, 2, 3]]
    assert np.allclose(mesh.point_data["depth"], [5, 6, 7, 8])
    names = {
        r.name: (sorted(np.asarray(r.entries).tolist()), r.tag) for r in mesh.regions
    }
    assert names == {
        "open_boundary_0": ([1, 2], 0),
        "land_boundary_0": ([0, 2, 3], 20),
        "land_boundary_1": ([0, 1], 3),
    }
    assert np.asarray(mesh.field_data["adcirc:boundary:land_boundary_0"]).tolist() == [
        2,
        3,
        0,
    ]
    assert np.allclose(mesh.point_data["zeta"], [0.1, 0.2, 0.3, 0.4])
    assert np.allclose(mesh.point_data["velocity"][1], [1.0, 0.5])
    assert np.allclose(mesh.field_data["meshio:time"], [3600.0])


def test_sparse_last_dataset(case, reader):
    mesh = reader(case / "fort.14", time_step=-1)
    assert np.allclose(mesh.point_data["zeta"], [-99999.0, 1.2, -99999.0, 1.4])
    assert np.allclose(mesh.field_data["meshio:time"], [7200.0])
    with pytest.raises(mio.ReadError, match="out of range"):
        reader(case / "fort.14", time_step=2)


def test_explicit_results_and_no_siblings(case, reader):
    grid_only = reader(case / "fort.14", siblings=False)
    assert "zeta" not in grid_only.point_data
    other = case / "grid.14"
    other.write_text(GRID)
    mesh = reader(other, results=[case / "fort.64"])
    assert "velocity" in mesh.point_data and "zeta" not in mesh.point_data


def test_time_values(case):
    assert mio.adcirc.time_values(case / "fort.14") == [3600.0, 7200.0]
    assert _adcirc.time_values(case / "fort.14") == [3600.0, 7200.0]


def test_write_round_trips_and_engines_write_the_same_bytes(case, tmp_path):
    mesh = mio.adcirc.read(case / "fort.14", siblings=False)
    mesh.regions = [r for r in mesh.regions if r.name != "land_boundary_1"]
    a, b = tmp_path / "a.14", tmp_path / "b.14"
    mio.adcirc.write(a, mesh)
    _adcirc.write(b, mesh)
    assert a.read_bytes() == b.read_bytes()
    back = mio.read(a)
    assert np.allclose(back.points, mesh.points)
    assert np.allclose(back.point_data["depth"], mesh.point_data["depth"])
    assert np.asarray(back.field_data["adcirc:boundary:land_boundary_0"]).tolist() == [
        2,
        3,
        0,
    ]


def test_barrier_boundaries_and_quads_are_refused(case, tmp_path):
    mesh = mio.adcirc.read(case / "fort.14", siblings=False)
    for write in (mio.adcirc.write, _adcirc.write):
        with pytest.raises(mio.WriteError, match="barrier"):
            write(tmp_path / "x.14", mesh)
    quads = mio.Mesh(np.zeros((4, 3)), [("quad", np.array([[0, 1, 2, 3]]))])
    for write in (mio.adcirc.write, _adcirc.write):
        with pytest.raises(mio.WriteError, match="triangles"):
            write(tmp_path / "q.14", quads)


def test_a_truncated_grid_is_a_read_error(tmp_path, reader):
    path = tmp_path / "bad.14"
    path.write_text("title\n5 4\n1 0 0 0\n")
    with pytest.raises(mio.ReadError, match="adcirc"):
        reader(path)
