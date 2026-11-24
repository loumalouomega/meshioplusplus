"""Statistics and derivatives over the steps of a sequence (roadmap 15.5)."""

import json

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import sequence_derivative, sequence_statistics

POINTS = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
TRIS = np.array([[0, 1, 2], [0, 2, 3]])


def _step(f, cell=None):
    mesh = mp.Mesh(POINTS.copy(), [("triangle", TRIS.copy())])
    mesh.point_data["f"] = np.asarray(f, dtype=float)
    mesh.point_data["v"] = np.column_stack([f, 2 * np.asarray(f), np.zeros(4)])
    if cell is not None:
        mesh.cell_data["c"] = [np.asarray(cell, dtype=float)]
    return mesh


def _series(tmp_path, times, fn, cell=None):
    paths = []
    for i, t in enumerate(times):
        path = str(tmp_path / f"s{i:02d}.vtu")
        mp.write(path, _step(fn(t), None if cell is None else cell(t)))
        paths.append(path)
    return paths


def _stack(times, fn):
    return np.stack([np.asarray(fn(t), dtype=float) for t in times])


FN = lambda t: np.array([t, t * t, 1.0, -t])  # noqa: E731


def test_statistics_match_numpy(tmp_path):
    times = [0.0, 1.0, 2.5, 3.0, 5.0]
    paths = _series(tmp_path, times, FN)
    out, rep = sequence_statistics(
        paths,
        ["f", "v"],
        ("mean", "std", "min", "max", "sum"),
        time_from="index",
        return_report=True,
    )
    data = _stack(times, FN)
    np.testing.assert_allclose(out.point_data["f:mean"], data.mean(axis=0))
    np.testing.assert_allclose(out.point_data["f:std"], data.std(axis=0), atol=1e-12)
    np.testing.assert_allclose(out.point_data["f:min"], data.min(axis=0))
    np.testing.assert_allclose(out.point_data["f:max"], data.max(axis=0))
    np.testing.assert_allclose(out.point_data["f:sum"], data.sum(axis=0))
    assert out.point_data["v:mean"].shape == (4, 3)
    assert rep["num_steps"] == 5 and rep["arrays"] == ["f", "v"]
    # the data of the steps is replaced by the statistics, the geometry is kept
    assert set(out.point_data) == {
        f"{a}:{s}" for a in ("f", "v") for s in ("mean", "std", "min", "max", "sum")
    }
    np.testing.assert_array_equal(out.points, POINTS)


def test_ddof_gives_the_sample_standard_deviation(tmp_path):
    times = [0.0, 1.0, 2.0, 4.0]
    paths = _series(tmp_path, times, FN)
    out = sequence_statistics(paths, ["f"], "std", ddof=1, time_from="index")
    np.testing.assert_allclose(
        out.point_data["f:std"], _stack(times, FN).std(axis=0, ddof=1)
    )


def test_a_single_step_has_zero_std_and_the_step_itself_as_mean(tmp_path):
    paths = _series(tmp_path, [0.0], FN)
    out = sequence_statistics(paths, ["f"], ("mean", "std"), time_from="index")
    np.testing.assert_allclose(out.point_data["f:mean"], FN(0.0))
    np.testing.assert_allclose(out.point_data["f:std"], 0.0)


def test_streaming_equals_the_in_memory_reduction_even_for_a_poorly_scaled_series(
    tmp_path,
):
    # a large offset: the naive sum-of-squares variance loses it, Welford's does not
    times = list(np.arange(12.0))
    fn = lambda t: 1e8 + np.array([t, 2 * t, 0.5 * t, 1.0])  # noqa: E731
    paths = _series(tmp_path, times, fn)
    out = sequence_statistics(paths, ["f"], ("mean", "std"), time_from="index")
    data = _stack(times, fn)
    np.testing.assert_allclose(out.point_data["f:std"], data.std(axis=0), rtol=1e-6)
    np.testing.assert_allclose(out.point_data["f:mean"], data.mean(axis=0), rtol=1e-12)


def test_cell_data_statistics_keep_the_block_structure(tmp_path):
    times = [0.0, 1.0, 2.0]
    cell = lambda t: [t, 10.0 + t]  # noqa: E731
    paths = _series(tmp_path, times, FN, cell)
    out = sequence_statistics(
        paths, ["c"], ("mean", "max"), location="cell", time_from="index"
    )
    assert len(out.cell_data["c:mean"]) == 1
    np.testing.assert_allclose(out.cell_data["c:mean"][0], [1.0, 11.0])
    np.testing.assert_allclose(out.cell_data["c:max"][0], [2.0, 12.0])


def test_a_nan_propagates_and_is_not_dropped(tmp_path):
    paths = _series(
        tmp_path, [0.0, 1.0], lambda t: np.array([1.0, np.nan if t else 2.0, 3.0, 4.0])
    )
    out = sequence_statistics(paths, ["f"], ("mean", "sum"), time_from="index")
    assert np.isnan(out.point_data["f:mean"][1]) and np.isnan(
        out.point_data["f:sum"][1]
    )
    assert out.point_data["f:mean"][0] == 1.0


@pytest.mark.parametrize(
    "kwargs,needle",
    [
        (dict(arrays=["f"], stats=("median",)), "unknown statistic"),
        (dict(arrays=["nothing"]), "no numeric point_data array 'nothing'"),
        (dict(arrays=["f"], stats=()), "no statistic"),
        (dict(arrays=["f"], location="field"), "field_data"),
    ],
)
def test_statistics_refusals(tmp_path, kwargs, needle):
    paths = _series(tmp_path, [0.0, 1.0], FN)
    with pytest.raises(ValueError, match=needle):
        sequence_statistics(paths, time_from="index", **kwargs)


def test_a_step_with_another_topology_is_refused_by_name(tmp_path):
    paths = _series(tmp_path, [0.0, 1.0], FN)
    other = mp.Mesh(POINTS[:3].copy(), [("triangle", np.array([[0, 1, 2]]))])
    other.point_data["f"] = np.array([1.0, 2.0, 3.0])
    mp.write(str(tmp_path / "s99.vtu"), other)
    with pytest.raises(ValueError, match="different number of points or cells"):
        sequence_statistics(
            paths + [str(tmp_path / "s99.vtu")], ["f"], time_from="index"
        )


def test_a_derivative_of_a_linear_in_time_field_is_exact(tmp_path):
    times = [0.0, 0.5, 2.0, 2.5, 4.0]  # unequal steps
    fn = lambda t: np.array([3.0 * t + 1.0, -2.0 * t, 5.0, 0.5 * t])  # noqa: E731
    paths = _series(tmp_path, times, fn)
    for scheme in ("central", "forward"):
        steps = list(sequence_derivative(paths, "f", scheme=scheme, times=times))
        assert len(steps) == 5
        for (t, mesh), t_expect in zip(steps, times):
            assert t == t_expect
            np.testing.assert_allclose(
                mesh.point_data["f:dt"], [3.0, -2.0, 0.0, 0.5], atol=1e-12
            )


def test_the_central_difference_is_second_order_on_unequal_steps(tmp_path):
    times = [0.0, 0.4, 1.0, 1.3, 2.0]
    fn = lambda t: np.array([t**2, np.sin(t), 1.0, t**3])  # noqa: E731
    paths = _series(tmp_path, times, fn)
    steps = list(sequence_derivative(paths, "f", times=times))
    # the interior of a quadratic is exact with the three-point formula
    for t, mesh in steps[1:-1]:
        assert mesh.point_data["f:dt"][0] == pytest.approx(2.0 * t, abs=1e-12)
    # the ends are one-sided and first order: the forward difference of t^2
    assert steps[0][1].point_data["f:dt"][0] == pytest.approx((0.4**2 - 0.0) / 0.4)
    assert steps[-1][1].point_data["f:dt"][0] == pytest.approx((2.0**2 - 1.3**2) / 0.7)


def test_the_derivative_of_a_vector_and_a_cell_array(tmp_path):
    times = [0.0, 1.0, 2.0]
    paths = _series(tmp_path, times, FN, lambda t: [t, 3.0 * t])
    v = list(sequence_derivative(paths, "v", times=times))
    np.testing.assert_allclose(
        v[1][1].point_data["v:dt"][:, 1], 2.0 * np.array([1.0, 2.0 * 1.0, 0.0, -1.0])
    )
    c = list(sequence_derivative(paths, "c", location="cell", times=times, output="dc"))
    np.testing.assert_allclose(c[1][1].cell_data["dc"][0], [1.0, 3.0])


def test_the_derivative_streams_into_a_written_sequence(tmp_path):
    times = [0.0, 1.0, 2.0]
    paths = _series(tmp_path, times, lambda t: np.array([t, 2 * t, 0.0, 1.0]))
    steps = sequence_derivative(paths, "f", times=times)
    written = mp.write_sequence(str(tmp_path / "d_{step}.vtu"), steps)
    assert len(written) == 3
    back = mp.read(written[1])
    np.testing.assert_allclose(
        back.point_data["f:dt"], [1.0, 2.0, 0.0, 0.0], atol=1e-12
    )


@pytest.mark.parametrize(
    "paths_n,times,needle",
    [(1, [0.0], "at least two steps"), (3, [0.0, 1.0, 1.0], "must increase")],
)
def test_derivative_refusals(tmp_path, paths_n, times, needle):
    paths = _series(tmp_path, times, FN)
    with pytest.raises(ValueError, match=needle):
        list(sequence_derivative(paths, "f", times=times))


def test_an_unknown_scheme_and_array(tmp_path):
    paths = _series(tmp_path, [0.0, 1.0], FN)
    with pytest.raises(ValueError, match="unknown scheme"):
        list(sequence_derivative(paths, "f", scheme="spline", time_from="index"))
    with pytest.raises(ValueError, match="no numeric point_data array 'nothing'"):
        list(sequence_derivative(paths, "nothing", time_from="index"))


def test_the_cli_verbs(tmp_path, capsys):
    from meshioplusplus._cli._main import main

    times = [0.0, 1.0, 2.0]
    _series(tmp_path, times, lambda t: np.array([t, 2 * t, 0.0, 1.0]))
    pattern = str(tmp_path / "s*.vtu")
    (tmp_path / "out").mkdir()
    out = str(tmp_path / "out" / "stats.vtu")
    rc = main(
        [
            "sequence-stats",
            pattern,
            out,
            "--arrays",
            "f",
            "--stats",
            "mean,max",
            "--time-from",
            "index",
            "--json",
        ]
    )
    assert rc in (0, None)
    report = json.loads(capsys.readouterr().out)
    assert report["num_steps"] == 3 and report["arrays"] == ["f"]
    back = mp.read(out)
    np.testing.assert_allclose(back.point_data["f:mean"], [1.0, 2.0, 0.0, 1.0])
    rc = main(
        [
            "sequence-derivative",
            pattern,
            str(tmp_path / "d_{step}.vtu"),
            "--array",
            "f",
            "--time-from",
            "index",
            "--json",
        ]
    )
    assert rc in (0, None)
    assert json.loads(capsys.readouterr().out)["num_steps"] == 3
