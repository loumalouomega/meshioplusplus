"""Tests for sequence time resampling (``resample_sequence``, ``blend_steps``)
and the ``agglomerate`` coplanar merge / sphericity gate."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import agglomerate, blend_steps, resample_sequence
from meshioplusplus._sequence import _blend_py, parse_times, resample_plan

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")


def _quads(n=3):
    xs, ys = np.meshgrid(np.arange(n + 1.0), np.arange(n + 1.0))
    pts = np.column_stack([xs.ravel(), ys.ravel(), np.zeros(xs.size)])
    quads = np.array(
        [
            [
                j * (n + 1) + i,
                j * (n + 1) + i + 1,
                (j + 1) * (n + 1) + i + 1,
                (j + 1) * (n + 1) + i,
            ]
            for j in range(n)
            for i in range(n)
        ]
    )
    return pts, quads


def _hex_block(n=3, stretch=(1.0, 1.0, 1.0)):
    axes = [np.arange(n + 1.0) * s for s in stretch]
    xs, ys, zs = np.meshgrid(*axes, indexing="ij")
    pts = np.column_stack([xs.ravel(), ys.ravel(), zs.ravel()])

    def vid(i, j, k):
        return (i * (n + 1) + j) * (n + 1) + k

    hexes = np.array(
        [
            [
                vid(i, j, k),
                vid(i + 1, j, k),
                vid(i + 1, j + 1, k),
                vid(i, j + 1, k),
                vid(i, j, k + 1),
                vid(i + 1, j, k + 1),
                vid(i + 1, j + 1, k + 1),
                vid(i, j + 1, k + 1),
            ]
            for i in range(n)
            for j in range(n)
            for k in range(n)
        ]
    )
    return meshioplusplus.Mesh(pts, [("hexahedron", hexes)])


def _step(t, pts, quads):
    # u is linear in time, so linear resampling must reproduce it exactly
    return meshioplusplus.Mesh(
        pts,
        [("quad", quads)],
        point_data={"u": pts[:, 0] + 2.0 * t, "id": np.arange(len(pts))},
        cell_data={"c": [np.full(len(quads), 3.0 * t)]},
    )


@pytest.fixture
def series(tmp_path):
    pts, quads = _quads()
    times = [0.0, 1.0, 3.0]
    for k, t in enumerate(times):
        meshioplusplus.write(str(tmp_path / f"src_{k}.vtu"), _step(t, pts, quads))
    return tmp_path, pts, quads, times


def test_parse_times():
    assert parse_times("0:1:0.25") == [0.0, 0.25, 0.5, 0.75, 1.0]
    assert parse_times("0:1:0.1")[-1] == pytest.approx(1.0)
    assert len(parse_times("0:1:0.1")) == 11
    assert parse_times("0.5, 2") == [0.5, 2.0]
    with pytest.raises(ValueError, match="positive step"):
        parse_times("0:1:0")


@pytest.mark.parametrize("method", ["linear", "nearest", "previous"])
@pytest.mark.parametrize("extrapolate", ["error", "clamp"])
def test_plan_twin_matches_core(method, extrapolate):
    if _core is None:
        pytest.skip("needs the compiled core")
    rng = np.random.default_rng(3)
    source = np.cumsum(rng.uniform(0.1, 1.0, 12)).tolist()
    targets = rng.uniform(source[0], source[-1], 40).tolist() + source
    if extrapolate == "clamp":
        targets += [source[0] - 1.0, source[-1] + 1.0]
    py = resample_plan(source, targets, method, extrapolate)
    core = [
        tuple(s)
        for s in _core.resample_plan(source, targets, method, extrapolate == "clamp")
    ]
    assert len(py) == len(core)
    for (a, b, w), (c, d, x) in zip(py, core):
        assert (a, b) == (c, d)
        assert w == pytest.approx(x, abs=1e-15)


def test_plan_errors():
    with pytest.raises(ValueError, match="outside the source range"):
        resample_plan([0.0, 1.0], [2.0])
    with pytest.raises(ValueError, match="increase strictly"):
        resample_plan([0.0, 0.0], [0.0])
    with pytest.raises(ValueError, match="unknown method"):
        resample_plan([0.0, 1.0], [0.5], "cubic")
    assert resample_plan([0.0, 1.0], [0.4], "nearest") == [(0, 0, 0.0)]
    assert resample_plan([0.0, 1.0], [0.6], "nearest") == [(1, 1, 0.0)]
    assert resample_plan([0.0, 1.0], [0.9], "previous") == [(0, 0, 0.0)]


def test_blend_steps_core_and_twin_agree():
    pts, quads = _quads()
    a, b = _step(0.0, pts, quads), _step(1.0, pts, quads)
    b.points = b.points + [0.0, 0.0, 1.0]
    for blend_points in (False, True):
        ref = _blend_py(a, b, 0.25, blend_points)
        out = blend_steps(a, b, 0.25, blend_points=blend_points)
        np.testing.assert_allclose(out.points, ref.points)
        np.testing.assert_allclose(out.point_data["u"], pts[:, 0] + 0.5)
        np.testing.assert_allclose(out.point_data["u"], ref.point_data["u"])
        # integers come from the nearer step, never blended
        assert out.point_data["id"].dtype.kind == "i"
        np.testing.assert_array_equal(out.point_data["id"], np.arange(len(pts)))
        np.testing.assert_allclose(out.cell_data["c"][0], 0.75)
        np.testing.assert_array_equal(out.cells[0].data, quads)
    assert blend_steps(a, b, 0.25).points[:, 2].max() == 0.0


def test_blend_steps_rejects_another_topology():
    pts, quads = _quads()
    other_pts, other_quads = _quads(2)
    with pytest.raises(ValueError, match="points"):
        blend_steps(_step(0, pts, quads), _step(1, other_pts, other_quads), 0.5)


@pytest.mark.parametrize("strict", [False, True])
def test_resample_linear_is_exact_on_a_linear_field(series, monkeypatch, strict):
    if strict:
        if _core is None:
            pytest.skip("needs the compiled core")
        monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    d, pts, quads, _ = series
    # with time_from="index" the source times are 0, 1, 2 (u's times 0, 1, 3)
    targets = [0.0, 0.5, 1.0, 1.75, 2.0]
    resample_sequence(
        str(d / "src_*.vtu"),
        str(d / "out_{index}.vtu"),
        times=targets,
        time_from="index",
    )
    for k, t in enumerate(targets):
        m = meshioplusplus.read(str(d / f"out_{k}.vtu"))
        expected = pts[:, 0] + 2.0 * np.interp(t, [0.0, 1.0, 2.0], [0.0, 1.0, 3.0])
        np.testing.assert_allclose(m.point_data["u"], expected, rtol=0, atol=1e-12)
    with pytest.raises(ValueError, match="outside the source range"):
        resample_sequence(
            str(d / "src_*.vtu"),
            str(d / "bad_{index}.vtu"),
            times=[3.0],
            time_from="index",
        )


def test_resample_explicit_times_clamp_and_methods(series):
    d, pts, quads, times = series
    paths = [str(d / f"src_{k}.vtu") for k in range(3)]

    def run(out, targets, **kw):
        # explicit source times through the pipeline document
        return meshioplusplus.run_sequence_pipeline(
            {
                "Version": 1,
                "Input": {"Paths": paths, "Times": times},
                "Output": {"Path": str(d / out)},
                "Resample": {"Times": targets, **kw},
            }
        )

    run("lin_{index}.vtu", [2.0])
    m = meshioplusplus.read(str(d / "lin_0.vtu"))
    np.testing.assert_allclose(m.point_data["u"], pts[:, 0] + 4.0)
    np.testing.assert_allclose(m.cell_data["c"][0], 6.0)
    assert float(np.asarray(m.field_data["meshio:time"]).ravel()[0]) == 2.0

    run("near_{index}.vtu", [2.2], Method="nearest")
    np.testing.assert_allclose(
        meshioplusplus.read(str(d / "near_0.vtu")).point_data["u"], pts[:, 0] + 6.0
    )
    run("prev_{index}.vtu", [2.9], Method="previous")
    np.testing.assert_allclose(
        meshioplusplus.read(str(d / "prev_0.vtu")).point_data["u"], pts[:, 0] + 2.0
    )
    run("clamp_{index}.vtu", [-1.0, 9.0], Extrapolate="clamp")
    np.testing.assert_allclose(
        meshioplusplus.read(str(d / "clamp_1.vtu")).point_data["u"], pts[:, 0] + 6.0
    )
    run("range_{index}.vtu", {"Start": 0, "Stop": 3, "Step": 1.5})
    assert (d / "range_2.vtu").exists() and not (d / "range_3.vtu").exists()
    with pytest.raises(ValueError, match="exactly one of Times and TimesFrom"):
        run("x_{index}.vtu", [])
    with pytest.raises(ValueError, match="Extrapolate"):
        run("x_{index}.vtu", [1.0], Extrapolate="wrap")


def test_resample_times_from_another_sequence(series):
    d, pts, quads, _ = series
    for k in range(3):
        meshioplusplus.write(str(d / f"other_{k}.vtu"), _step(0.0, pts, quads))
    resample_sequence(
        str(d / "src_*.vtu"),
        str(d / "al_{index}.vtu"),
        times_from=str(d / "other_*.vtu"),
        time_from="index",
    )
    assert sorted(p.name for p in d.glob("al_*.vtu")) == [
        f"al_{k}.vtu" for k in range(3)
    ]


def test_resample_to_a_series_file(series):
    d, pts, quads, _ = series
    out = d / "series.xdmf"
    resample_sequence(
        str(d / "src_*.vtu"), str(out), times="0:2:0.5", time_from="index"
    )
    meta = meshioplusplus.read_metadata(str(out))
    np.testing.assert_allclose(meta["time_values"], [0.0, 0.5, 1.0, 1.5, 2.0])


@needs_core
def test_agglomerate_coplanar_merge_conserves_volume_and_fuses_faces():
    mesh = _hex_block()
    plain, rep0 = agglomerate(mesh, target_group_size=8, return_report=True)
    merged, rep = agglomerate(
        mesh, target_group_size=8, merge_coplanar_faces=True, return_report=True
    )
    assert rep0["num_faces_merged"] == 0
    assert rep["num_faces_merged"] > 0
    assert len(merged.cells[0].data) == len(plain.cells[0].data)
    vol = meshioplusplus.compute_stats(merged)["signed_volume"]
    assert vol == pytest.approx(27.0)

    def nfaces(m):
        return sum(len(cell) for cell in m.cells[0].data)

    assert nfaces(merged) < nfaces(plain)


@needs_core
def test_agglomerate_sphericity_gate_rejects_elongated_growth():
    mesh = _hex_block(n=4, stretch=(4.0, 1.0, 1.0))
    loose, rep0 = agglomerate(mesh, target_group_size=8, return_report=True)
    tight, rep = agglomerate(
        mesh, target_group_size=8, min_sphericity=0.8, return_report=True
    )
    assert rep0["num_rejected"] == 0
    assert rep["num_rejected"] > 0
    assert len(tight.cells[0].data) > len(loose.cells[0].data)
    assert meshioplusplus.compute_stats(tight)["signed_volume"] == pytest.approx(
        4.0**3 * 4.0
    )


@needs_core
def test_agglomerate_pipeline_keys(tmp_path):
    src = str(tmp_path / "hex.vtu")
    meshioplusplus.write(src, _hex_block())
    report = meshioplusplus.run_pipeline(
        {
            "Version": 1,
            "Input": {"Path": src},
            "Operations": [
                {
                    "Op": "Agglomerate",
                    "TargetGroupSize": 8,
                    "MergeCoplanarFaces": True,
                    "CoplanarAngle": 2.0,
                    "MinSphericity": 0.1,
                }
            ],
            "Output": {"Path": str(tmp_path / "out.vtu")},
        }
    )
    assert report["steps"][0]["op"] == "Agglomerate"
    assert report["steps"][0]["NumFacesMerged"] > 0
    assert report["steps"][0]["NumRejected"] >= 0
    with pytest.raises(ValueError, match="MinSphericityy|unknown key"):
        meshioplusplus.run_pipeline(
            {
                "Version": 1,
                "Input": {"Path": src},
                "Operations": [{"Op": "Agglomerate", "MinSphericityy": 0.1}],
                "Output": {"Path": str(tmp_path / "out2.vtu")},
            }
        )


def _poly_volume_area(points, faces):
    vol, area = 0.0, 0.0
    for face in faces:
        p = points[np.asarray(face)]
        c = p.mean(axis=0)
        a = 0.5 * sum(np.cross(p[k], p[(k + 1) % len(p)]) for k in range(len(p)))
        vol += np.dot(c, a) / 3.0
        area += np.linalg.norm(a)
    return vol, area


@needs_core
@pytest.mark.parametrize("threshold", [0.6, 0.7, 0.75])
def test_agglomerate_groups_meet_the_sphericity_threshold(threshold):
    mesh = _hex_block(n=6)
    out, rep = agglomerate(
        mesh, target_group_size=8, min_sphericity=threshold, return_report=True
    )
    assert rep["num_rejected"] > 0
    for cell in out.cells[0].data:
        v, a = _poly_volume_area(out.points, cell)
        assert v > 0
        assert np.pi ** (1 / 3) * (6 * v) ** (2 / 3) / a >= threshold - 1e-12
    assert meshioplusplus.compute_stats(out)["signed_volume"] == pytest.approx(216.0)
