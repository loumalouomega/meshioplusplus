"""Descriptive statistics of data arrays per region (roadmap 15.5)."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import data_integrate, data_statistics
from meshioplusplus._regions import Region

POINTS = np.array([[0.0, 0, 0], [2, 0, 0], [2, 1, 0], [0, 1, 0], [4, 0, 0]])
TRIS = np.array([[0, 1, 2], [0, 2, 3], [1, 4, 2]])


def _mesh():
    mesh = mp.Mesh(POINTS.copy(), [("triangle", TRIS.copy())])
    mesh.point_data["f"] = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    mesh.point_data["v"] = np.column_stack(
        [mesh.point_data["f"], -mesh.point_data["f"]]
    )
    mesh.cell_data["c"] = [np.array([1.0, 3.0, 10.0])]
    return mesh


def test_point_statistics_match_numpy():
    st = data_statistics(_mesh(), ["f"], bins=4)[0]
    f = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    d = st["domain"]
    assert d["count"] == [5] and d["num_nan"] == [0]
    assert d["sum"][0] == f.sum()
    assert d["mean"][0] == pytest.approx(f.mean())
    assert d["variance"][0] == pytest.approx(f.var())
    assert d["std"][0] == pytest.approx(f.std())
    assert d["quantiles"]["0.5"][0] == pytest.approx(np.median(f))
    assert d["quantiles"]["0.05"][0] == pytest.approx(np.quantile(f, 0.05))
    counts, edges = np.histogram(f, bins=4)
    assert d["histogram"][0]["counts"] == counts.tolist()
    np.testing.assert_allclose(d["histogram"][0]["edges"], edges)


def test_ddof_and_components():
    st = data_statistics(_mesh(), ["v"], ddof=1)[0]
    assert st["num_components"] == 2
    f = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    assert st["domain"]["variance"][1] == pytest.approx(f.var(ddof=1))
    assert st["domain"]["max"][1] == -1.0


def test_regions_point_cell_and_side():
    mesh = _mesh()
    mesh.regions.append(Region("pts", "point", [0, 4]))
    mesh.regions.append(Region("left", "cell", [0, 1]))
    mesh.regions.append(Region("edge", "side", [[0, 0]]))
    st = data_statistics(mesh, ["f"])[0]
    by = {r["name"]: r for r in st["regions"]}
    assert set(by) == {"pts", "left"}  # the side region has no point image
    assert by["pts"]["mean"][0] == pytest.approx(8.5)
    # the cell region covers the nodes of cells 0 and 1: points 0..3
    assert by["left"]["count"] == [4]
    assert by["left"]["mean"][0] == pytest.approx(3.75)
    cell = {
        r["name"]: r
        for r in data_statistics(mesh, ["c"], location="cell")[0]["regions"]
    }
    assert set(cell) == {"left"} and cell["left"]["mean"] == [2.0]
    assert data_statistics(mesh, ["f"], per_region=False)[0]["regions"] == []


def test_weighted_cell_mean_equals_the_integral_mean():
    mesh = _mesh()
    mesh.regions.append(Region("left", "cell", [0, 1]))
    st = data_statistics(mesh, ["c"], location="cell", weighted=True)[0]
    ref = {a["name"]: a for a in data_integrate(mesh, ["c"])}["c"]
    assert st["domain"]["mean"][0] == pytest.approx(
        ref["domain"]["mean_per_component"][0]
    )
    # areas 1, 1, 1 -> the weighted mean is the plain one; then enlarge a cell
    mesh.points[4] = [6.0, 0.0, 0.0]
    st = data_statistics(mesh, ["c"], location="cell", weighted=True)[0]
    ref = {a["name"]: a for a in data_integrate(mesh, ["c"])}["c"]
    assert st["domain"]["mean"][0] == pytest.approx(
        ref["domain"]["mean_per_component"][0]
    )
    assert st["domain"]["mean"][0] != pytest.approx(np.mean([1.0, 3.0, 10.0]))


def test_nonfinite_values_are_counted_not_dropped():
    mesh = _mesh()
    mesh.point_data["f"] = np.array([1.0, np.nan, np.inf, 3.0, 5.0])
    d = data_statistics(mesh, ["f"])[0]["domain"]
    assert (d["count"], d["num_nan"], d["num_inf"]) == ([3], [1], [1])
    assert d["mean"][0] == pytest.approx(3.0)


@pytest.mark.parametrize(
    "kwargs,needle",
    [
        (dict(arrays=["nope"]), "no point_data array 'nope'"),
        (dict(location="field"), "field_data"),
        (dict(quantiles=(1.5,)), "quantiles"),
        (dict(bins=-1), "bins"),
        (dict(weighted=True), "weighted"),
    ],
)
def test_refusals(kwargs, needle):
    with pytest.raises(ValueError, match=needle):
        data_statistics(_mesh(), **kwargs)


def test_cli_verb(tmp_path, capsys):
    from meshioplusplus._cli._main import main

    path = str(tmp_path / "m.vtu")
    mp.write(path, _mesh())
    rc = main(["data", "stats", path, "--arrays", "f", "--bins", "2", "--json"])
    assert rc in (0, None)
    out = json.loads(capsys.readouterr().out)
    assert out[0]["name"] == "f" and out[0]["domain"]["count"] == [5]
    assert len(out[0]["domain"]["histogram"][0]["counts"]) == 2
