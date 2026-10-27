"""project_to_cad (roadmap 6.10.1): a mesh's points moved back onto the CAD
faces and edges its cells name.

Through the Python wrapper: a smoothed sphere goes back onto the sphere, the
capped cylinder's rims stay on their circles, a point pushed into a hole is
reported outside its face's trims, the refusals, the CLI verb and the
pipeline step through both engines. The algorithm itself is tested in the
core (``tests/cpp/test_project_to_cad.cpp``).
"""

import json
import pathlib

import numpy as np
import pytest

import meshioplusplus as mio

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - the numpy-only install
    _core = None

needs_core = pytest.mark.skipif(
    not hasattr(_core, "project_to_cad"), reason="core predates project_to_cad"
)

CAD = pathlib.Path(__file__).parent / "meshes" / "cad_json"


def _read(name):
    return mio.read(CAD / name)


@needs_core
def test_a_smoothed_sphere_goes_back():
    model = _read("sphere.cad.json")
    mesh = mio.mesh_cad_surface(model, 0.2)
    smoothed = mio.smooth(mesh, iterations=10, method="laplacian")
    shrink = np.abs(np.linalg.norm(smoothed.points, axis=1) - 1.0).max()
    assert shrink > 1e-2
    out, report = mio.project_to_cad(smoothed, model, return_report=True)
    assert report["num_projected"] == len(mesh.points)
    assert report["num_ambiguous"] == report["num_failed"] == 0
    assert report["max_move"] == pytest.approx(shrink, rel=0.5)
    assert np.allclose(np.linalg.norm(out.points, axis=1), 1.0, atol=1e-12)
    # Cells, data and regions pass through.
    assert [c.type for c in out.cells] == [c.type for c in mesh.cells]
    assert set(out.cell_data) == set(mesh.cell_data)
    # The bare B-rep dict serves as well, and the result is deterministic.
    again = mio.project_to_cad(smoothed, model.info["brep"])
    assert np.array_equal(again.points, out.points)


@needs_core
def test_the_cylinder_keeps_its_rims():
    model = _read("capped_cylinder.cad.json")
    mesh = mio.mesh_cad_surface(model, 0.25)
    rng = np.random.default_rng(7)
    off = mio.Mesh(
        np.asarray(mesh.points) + rng.uniform(-0.02, 0.02, (len(mesh.points), 3)),
        mesh.cells,
        cell_data=mesh.cell_data,
    )
    out, report = mio.project_to_cad(off, model, return_report=True)
    assert report["num_on_edges"] > 0 and report["num_ambiguous"] == 0
    p, p0 = np.asarray(out.points), np.asarray(mesh.points)
    r, r0 = np.hypot(p[:, 0], p[:, 1]), np.hypot(p0[:, 0], p0[:, 1])
    rim = np.isclose(r0, 1.0) & (np.isclose(p0[:, 2], 0.0) | np.isclose(p0[:, 2], 2.0))
    assert rim.any()
    assert np.allclose(r[rim], 1.0, atol=1e-10)
    assert np.allclose(p[rim, 2], p0[rim, 2], atol=1e-10)
    on_side = np.isclose(r, 1.0, atol=1e-10)
    on_cap = np.isclose(p[:, 2], 0.0, atol=1e-10) | np.isclose(p[:, 2], 2.0, atol=1e-10)
    assert (on_side | on_cap).all()


@needs_core
def test_a_point_outside_the_trims_is_reported():
    """The 4 x 4 plate's hole is the circle of radius 1 about (2, 2): a plate
    point pushed into it lands on the plane there, outside the face's trims,
    and is moved, counted and named."""
    model = _read("plate_with_hole.cad.json")
    mesh = mio.mesh_cad_surface(model, 0.2)
    p = np.asarray(mesh.points).copy()
    on_line = np.unique(
        np.concatenate([c.data.ravel() for c in mesh.cells if c.type == "line"])
    )
    pick = int(np.setdiff1d(np.arange(len(p)), on_line)[0])
    p[pick] = [2.0, 2.0, 0.05]
    off = mio.Mesh(p, mesh.cells, cell_data=mesh.cell_data)
    out, report = mio.project_to_cad(off, model, return_report=True)
    assert report["num_outside_trim"] == 1
    assert report["outside_trim"] == [pick]
    assert np.allclose(np.asarray(out.points)[pick], [2.0, 2.0, 0.0], atol=1e-12)
    assert report["warnings"]
    _, far = mio.project_to_cad(off, model, max_distance=1e-3, return_report=True)
    assert far["num_too_far"] == 1 and far["num_outside_trim"] == 0


@needs_core
def test_points_on_a_straight_trim_are_inside(tmp_path):
    """A triangular face bounded by three straight trims: a mesh's own points on
    them lie on the trims to rounding, so none is outside; a point pushed well
    past the hypotenuse is still reported."""
    doc = {
        "breps": [
            {
                "brep_id": 1,
                "faces": [
                    {
                        "brep_id": 2,
                        "surface": {
                            "is_trimmed": True,
                            "degrees": [1, 1],
                            "knot_vectors": [[0, 2], [0, 2]],
                            "control_points": [
                                [1, [0, 0, 0, 1]],
                                [2, [2, 0, 0, 1]],
                                [3, [0, 2, 0, 1]],
                                [4, [2, 2, 0, 1]],
                            ],
                        },
                        "boundary_loops": [
                            {
                                "loop_type": "outer",
                                "trimming_curves": [
                                    {
                                        "trim_index": i,
                                        "curve_direction": True,
                                        "parameter_curve": {
                                            "degree": 1,
                                            "knot_vector": [0, 0, 1, 1],
                                            "active_range": [0, 1],
                                            "control_points": [
                                                [[a[0], a[1], 0, 1]],
                                                [[b[0], b[1], 0, 1]],
                                            ],
                                        },
                                    }
                                    for i, (a, b) in enumerate(
                                        [
                                            ((0, 0), (1.5, 0)),
                                            ((1.5, 0), (0, 1.5)),
                                            ((0, 1.5), (0, 0)),
                                        ]
                                    )
                                ],
                            }
                        ],
                    }
                ],
                "edges": [],
            }
        ]
    }
    path = tmp_path / "triangle.cad.json"
    path.write_text(json.dumps(doc))
    model = mio.read(path)
    mesh = mio.mesh_cad_surface(model, 0.1)
    _, report = mio.project_to_cad(mesh, model, return_report=True)
    assert report["num_outside_trim"] == 0 and report["outside_trim"] == []
    p = np.asarray(mesh.points).copy()
    on_hypotenuse = np.flatnonzero(
        np.isclose(p[:, 0] + p[:, 1], 1.5, atol=1e-9)
        & (p[:, 0] > 0.2)
        & (p[:, 0] < 1.3)
    )
    pick = int(on_hypotenuse[0])
    p[pick] += [0.05, 0.05, 0.0]
    off = mio.Mesh(p, mesh.cells, cell_data=mesh.cell_data)
    _, report = mio.project_to_cad(off, model, return_report=True)
    assert report["outside_trim"] == [pick]


@needs_core
def test_refusals_are_by_name():
    model = _read("sphere.cad.json")
    mesh = mio.mesh_cad_surface(model, 0.4)
    with pytest.raises(ValueError, match="cad:face"):
        mio.project_to_cad(model, model)
    with pytest.raises(ValueError, match="trim_samples"):
        mio.project_to_cad(mesh, model, trim_samples=1)
    with pytest.raises(ValueError, match="max_distance"):
        mio.project_to_cad(mesh, model, max_distance=-1.0)
    with pytest.raises(ValueError, match="B-rep"):
        mio.project_to_cad(mesh, mesh)


@needs_core
def test_the_cli_verb(tmp_path, capsys):
    model = CAD / "sphere.cad.json"
    meshed = tmp_path / "sphere.vtu"
    mio.write(
        meshed, mio.smooth(mio.mesh_cad_surface(mio.read(model), 0.3), iterations=5)
    )
    out = tmp_path / "back.vtu"
    code = mio._cli.main(
        ["project-to-cad", str(meshed), str(model), str(out), "--json"]
    )
    assert code in (0, None)
    report = json.loads(capsys.readouterr().out)
    assert report["num_projected"] > 0 and report["max_move"] > 0
    assert 0 <= report["num_on_curves"] <= report["num_on_edges"]
    back = mio.read(out)
    assert np.allclose(np.linalg.norm(back.points, axis=1), 1.0, atol=1e-12)


@needs_core
@pytest.mark.parametrize("version", [1, 2])
def test_the_pipeline_step(tmp_path, version):
    """MeshCadSurface, a Smooth that drops the model, then ProjectToCad on
    the model at Model: both engines put the sphere back."""
    doc = {
        "Version": version,
        "Input": {"Path": str(CAD / "sphere.cad.json")},
        "Operations": [
            {"Op": "MeshCadSurface", "TargetSize": 0.3},
            {"Op": "Smooth", "Iterations": 5},
            {"Op": "ProjectToCad", "Model": str(CAD / "sphere.cad.json")},
        ],
    }
    engines = [("py", mio.run_pipeline)]
    if hasattr(_core, "run_pipeline_json"):
        engines.append(("cpp", lambda d: _core.run_pipeline_json(json.dumps(d))))
    reports = []
    for name, run in engines:
        path = tmp_path / f"{name}.vtu"
        report = run(dict(doc, Output={"Path": str(path)}))
        back = mio.read(path)
        assert np.allclose(np.linalg.norm(back.points, axis=1), 1.0, atol=1e-12), name
        step = report["steps"][-1]
        assert step["NumProjected"] == len(back.points)
        reports.append(step)
    if len(reports) == 2:
        assert reports[0]["NumProjected"] == reports[1]["NumProjected"]
    # Without Model, after a step that changed the mesh: refused by name.
    bad = dict(doc, Operations=doc["Operations"][:2] + [{"Op": "ProjectToCad"}])
    for _, run in engines:
        with pytest.raises(ValueError, match="Model"):
            run(dict(bad, Output={"Path": str(tmp_path / "bad.vtu")}))
