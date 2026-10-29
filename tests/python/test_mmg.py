"""The optional MMG backend (roadmap 3.2): adapt and its helpers.

The probes run only in a core built with -DMESHIOPLUSPLUS_WITH_MMG=ON; the
named refusal runs only in one without. metric_from_hessian needs no MMG.
"""

import os
import subprocess
import sys

import numpy as np
import pytest

import meshioplusplus as mio

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - the pure-Python install
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")
needs_mmg = pytest.mark.skipif(
    not mio.has_mmg(), reason="needs a core built with -DMESHIOPLUSPLUS_WITH_MMG=ON"
)
without_mmg = pytest.mark.skipif(
    _core is None or mio.has_mmg(), reason="checks the build without MMG"
)


def _edge_lengths(mesh):
    pairs = {"triangle": [(0, 1), (1, 2), (2, 0)]}
    pairs["tetra"] = [(a, b) for a in range(4) for b in range(a + 1, 4)]
    edges = set()
    for block in mesh.cells:
        for a, b in pairs.get(block.type, []):
            for i, j in block.data[:, [a, b]]:
                edges.add((min(i, j), max(i, j)))
    edges = np.array(sorted(edges))
    return np.linalg.norm(mesh.points[edges[:, 0]] - mesh.points[edges[:, 1]], axis=1)


def _square():
    return mio.rectangle(dims=(4, 4), cell_type="triangle")


def test_capabilities_agree():
    assert mio.has_mmg() is bool(getattr(_core, "__has_mmg__", False))
    if not mio.has_mmg():
        assert not mio.has_mmg_elas()


@without_mmg
def test_adapt_without_mmg_names_the_flag():
    with pytest.raises(NotImplementedError, match="MESHIOPLUSPLUS_WITH_MMG"):
        mio.adapt(_square(), hsiz=0.1)


@needs_mmg
@pytest.mark.parametrize(
    "mesh, hsiz, library",
    [
        (_square(), 0.1, "mmg2d"),
        (mio.box(dims=(3, 3, 3), cell_type="tetra"), 0.2, "mmg3d"),
    ],
    ids=["planar", "volume"],
)
def test_a_constant_size_gives_edges_of_that_length(mesh, hsiz, library):
    # MMG aims at unit edges in the metric of hsiz; a mean within 20% of it
    # is the probe of roadmap 3.2.7.
    out, report = mio.adapt(mesh, hsiz=hsiz, return_report=True)
    assert report["library"] == library
    assert report["status"] == "success"
    assert abs(_edge_lengths(out).mean() - hsiz) < 0.2 * hsiz
    assert report["num_elements_out"] > report["num_elements_in"]
    assert 1.0 <= report["quality"]["min"] <= report["quality"]["mean"]


@needs_mmg
def test_a_curved_surface_goes_to_mmgs():
    sphere = mio.icosa_sphere(subdivisions=2)
    out, report = mio.adapt(sphere, hsiz=0.2, hausd=0.01, return_report=True)
    assert report["library"] == "mmgs"
    radius = np.linalg.norm(out.points, axis=1)
    assert np.all(np.abs(radius - 1.0) < 0.02)


@needs_mmg
def test_side_regions_come_back_on_their_edges():
    out = mio.adapt(_square(), hsiz=0.1)
    names = {r.name for r in out.regions}
    assert names == {"xmin", "xmax", "ymin", "ymax"}
    local_edges = {0: (0, 1), 1: (1, 2), 2: (2, 0)}
    triangles = out.cells[0].data
    for region in out.regions:
        axis, value = (0 if region.name[0] == "x" else 1), float(
            region.name.endswith("max")
        )
        assert len(region.entries) >= 8
        for cell, facet in np.asarray(region.entries):
            nodes = triangles[cell][list(local_edges[facet])]
            np.testing.assert_allclose(out.points[nodes][:, axis], value, atol=1e-12)


@needs_mmg
def test_frozen_points_survive_bit_for_bit():
    mesh = _square()
    frozen = [6, 12, 18]
    out = mio.adapt(mesh, hsiz=0.1, frozen=frozen)
    for p in mesh.points[frozen]:
        assert np.any(np.all(out.points == p, axis=1))


@needs_mmg
def test_a_required_point_region_survives_under_its_name():
    mesh = _square()
    mesh.regions = list(mesh.regions) + [mio.Region("pins", "point", np.array([6, 18]))]
    out = mio.adapt(mesh, hsiz=0.1, required=["pins"])
    pins = [r for r in out.regions if r.name == "pins"]
    assert len(pins) == 1
    np.testing.assert_array_equal(
        np.sort(out.points[pins[0].entries], axis=0),
        np.sort(mesh.points[[6, 18]], axis=0),
    )


@needs_mmg
def test_cell_regions_and_medit_refs_come_back():
    box = mio.box(dims=(2, 2, 2), cell_type="tetra")
    n = len(box.cells[0].data)
    box.cell_data["medit:ref"] = [np.where(np.arange(n) < n // 2, 3, 7)]
    box.regions = list(box.regions) + [mio.Region("half", "cell", np.arange(n // 2))]
    out = mio.adapt(box, hsiz=0.25)
    refs = out.cell_data["medit:ref"][0]
    assert set(np.unique(refs)) == {3, 7}
    half = [r for r in out.regions if r.name == "half"][0]
    np.testing.assert_array_equal(np.sort(half.entries), np.flatnonzero(refs == 3))


@needs_mmg
def test_an_anisotropic_metric_stretches_the_triangles():
    mesh = _square()
    metric = np.zeros((len(mesh.points), 3))
    metric[:, 0] = 1 / 0.2**2  # h_x = 0.2
    metric[:, 2] = 1 / 0.02**2  # h_y = 0.02
    out = mio.adapt(mesh, metric=metric)
    corners = out.points[out.cells[0].data]
    dx = np.ptp(corners[:, :, 0], axis=1)
    dy = np.ptp(corners[:, :, 1], axis=1)
    assert np.median(dx / dy) > 5


@needs_mmg
def test_a_size_field_grades_the_mesh():
    mesh = _square()
    size = 0.02 + 0.2 * mesh.points[:, 0]
    out = mio.adapt(mesh, size=size)
    corners = out.points[out.cells[0].data]
    centroid_x = corners[:, :, 0].mean(axis=1)
    area = 0.5 * np.abs(
        np.cross(corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0])
    )
    assert area[centroid_x < 0.2].mean() * 20 < area[centroid_x > 0.8].mean()


@needs_mmg
def test_a_local_parameter_refines_its_region():
    mesh = _square()
    out = mio.adapt(mesh, hsiz=0.2, local={"xmin": {"hmax": 0.02, "hmin": 0.01}})
    edge = [r for r in out.regions if r.name == "xmin"][0]
    far = [r for r in out.regions if r.name == "xmax"][0]
    assert len(edge.entries) > 3 * len(far.entries)


@needs_mmg
def test_refusals_are_by_name():
    with pytest.raises(ValueError, match="convert_cells"):
        mio.adapt(mio.box(cell_type="hexahedron"))
    with pytest.raises(ValueError, match="one of hsiz, size and metric"):
        mio.adapt(_square(), hsiz=0.1, size=np.ones(25))
    with pytest.raises(ValueError, match="no region"):
        mio.adapt(_square(), hsiz=0.1, required=["nope"])
    with pytest.raises(ValueError, match="6 components"):
        mio.adapt(mio.box(dims=(2, 2, 2), cell_type="tetra"), metric=np.ones((27, 3)))


@needs_mmg
def test_data_is_dropped_with_a_warning_and_field_data_kept():
    mesh = _square()
    mesh.point_data["u"] = mesh.points[:, 0].copy()
    mesh.field_data["note"] = np.array([42])
    out, report = mio.adapt(mesh, hsiz=0.2, return_report=True)
    assert "u" not in out.point_data
    assert any("'u' is dropped" in w for w in report["warnings"])
    np.testing.assert_array_equal(out.field_data["note"], [42])


@needs_mmg
def test_the_host_keeps_its_ctrl_c():
    # MMG installs exit() handlers for SIGINT and five others on every call;
    # the core puts the host's back, so Ctrl-C still raises afterwards.
    code = (
        "import os, signal, meshioplusplus as mio\n"
        "mio.adapt(mio.rectangle(dims=(2, 2), cell_type='triangle'), hsiz=0.3)\n"
        "try:\n"
        "    os.kill(os.getpid(), signal.SIGINT)\n"
        "    import time; time.sleep(1)\n"
        "except KeyboardInterrupt:\n"
        "    print('interrupted')\n"
    )
    done = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, env=os.environ
    )
    assert done.returncode == 0, done.stdout + done.stderr
    assert "interrupted" in done.stdout


@needs_mmg
def test_a_metric_from_a_medit_solution(tmp_path):
    # The file route: a metric in a .sol beside the mesh is the adapt input.
    mesh = _square()
    metric = np.zeros((len(mesh.points), 3))
    metric[:, 0] = metric[:, 2] = 1 / 0.1**2
    mesh.point_data["metric"] = metric
    mio.medit.write(tmp_path / "square.mesh", mesh, sol=True)
    back = mio.medit.read(tmp_path / "square.mesh", sol=True)
    out = mio.adapt(back, metric="medit:sol:0")
    assert abs(_edge_lengths(out).mean() - 0.1) < 0.02


@needs_core
def test_metric_from_hessian():
    mesh = _square()
    hessian = np.zeros((len(mesh.points), 9))
    hessian[:, 0], hessian[:, 4] = 2.0, 20.0  # f = x^2 + 10 y^2
    mesh.point_data["h"] = hessian
    out = mio.metric_from_hessian(mesh, "h", error=0.01, hmin=1e-3, hmax=10.0)
    c = 0.5 * (3 / 4) ** 2
    np.testing.assert_allclose(
        out.point_data["metric"][0], [c * 2 / 0.01, 0, c * 20 / 0.01, 0, 0, 1 / 100]
    )
    # 2-D Hessians give 2-D metrics; clamping bounds the sizes
    mesh.point_data["h2"] = np.tile([1e6, 0, 0, 0.0], (len(mesh.points), 1))
    out = mio.metric_from_hessian(
        mesh, "h2", error=1.0, hmin=0.1, hmax=1.0, output="m2"
    )
    np.testing.assert_allclose(out.point_data["m2"][0], [100, 0, 1])
    with pytest.raises(ValueError, match="positive"):
        mio.metric_from_hessian(mesh, "h", error=0.0)
    mesh.point_data["h3"] = np.ones((len(mesh.points), 3))
    with pytest.raises(ValueError, match="9 components"):
        mio.metric_from_hessian(mesh, "h3")


@needs_mmg
def test_the_adapt_pipeline_step_in_both_engines(tmp_path):
    import json

    src = tmp_path / "square.vtu"
    mio.write(str(src), _square(), compression=None)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "Adapt", "Hsiz": 0.1, "Required": ["xmin"]}],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    report = mio.run_pipeline(settings)
    step = report["steps"][0]
    assert step["op"] == "Adapt" and step["LowFailure"] == 0.0
    assert step["NumElementsOut"] > step["NumElementsIn"] == 32
    if hasattr(_core, "run_pipeline_json"):
        settings["Output"]["Path"] = str(tmp_path / "cpp.vtu")
        cpp = _core.run_pipeline_json(json.dumps(settings))["steps"][0]
        assert cpp["NumElementsOut"] == step["NumElementsOut"]
        assert cpp["NumPointsOut"] == step["NumPointsOut"]


def _interface_area(mesh):
    """The area between the cells of different medit:ref, by shared faces."""
    from collections import defaultdict

    tets, refs = mesh.cells[0].data, mesh.cell_data["medit:ref"][0]
    owners = defaultdict(list)
    for i, t in enumerate(tets):
        for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            owners[tuple(sorted(t[list(f)]))].append(i)
    area = 0.0
    for face, cells in owners.items():
        if len(cells) == 2 and refs[cells[0]] != refs[cells[1]]:
            p = mesh.points[list(face)]
            area += 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
    return area


@without_mmg
def test_the_level_set_and_domain_operations_name_the_flag_without_mmg():
    with pytest.raises(NotImplementedError, match="MESHIOPLUSPLUS_WITH_MMG"):
        mio.discretize_levelset(_square(), np.ones(25))
    lines = mio.Mesh(np.eye(3)[:, :2], [("line", np.array([[0, 1], [1, 2], [2, 0]]))])
    with pytest.raises(NotImplementedError, match="MESHIOPLUSPLUS_WITH_MMG"):
        mio.mesh_domain(lines, hsiz=0.1)


@needs_mmg
def test_a_sphere_level_set_gives_its_area():
    # Roadmap 3.2.7's probe: the interface of a sphere's signed distance has
    # area ~ 4 pi r^2 (2% under it here, the polygonal approximation).
    box = mio.box(
        dims=(10, 10, 10), size=(2, 2, 2), origin=(-1, -1, -1), cell_type="tetra"
    )
    ls = np.linalg.norm(box.points, axis=1) - 0.6
    out, report = mio.discretize_levelset(box, ls, hausd=0.005, return_report=True)
    assert report["library"] == "mmg3d" and report["status"] == "success"
    assert abs(_interface_area(out) - 4 * np.pi * 0.36) < 0.03 * 4 * np.pi * 0.36
    names = {r.name: r for r in out.regions}
    assert {"inside", "outside", "interface", "xmin", "zmax"} <= set(names)
    refs = out.cell_data["medit:ref"][0]
    assert set(np.unique(refs)) == {2, 3}
    np.testing.assert_array_equal(
        np.sort(names["inside"].entries), np.flatnonzero(refs == 3)
    )
    # the interface sides are on the inside elements
    assert np.all(refs[np.asarray(names["interface"].entries)[:, 0]] == 3)
    # the inside is the ball: its points are within the radius
    inside_points = out.points[np.unique(out.cells[0].data[refs == 3])]
    assert np.linalg.norm(inside_points, axis=1).max() < 0.6 + 0.01


@needs_mmg
def test_a_disc_level_set_in_the_plane():
    square = mio.rectangle(
        dims=(10, 10), size=(2, 2), origin=(-1, -1, 0), cell_type="triangle"
    )
    square.point_data["ls"] = np.linalg.norm(square.points[:, :2], axis=1) - 0.5
    out, report = mio.discretize_levelset(square, "ls", return_report=True)
    assert report["library"] == "mmg2d"
    names = {r.name for r in out.regions}
    assert {"inside", "outside", "interface", "xmin"} <= names
    # the input's level set itself is not carried: the points are new
    assert "ls" not in out.point_data


@needs_mmg
def test_level_set_materials_name_the_parts():
    square = mio.rectangle(
        dims=(8, 8), size=(2, 2), origin=(-1, -1, 0), cell_type="triangle"
    )
    n = len(square.cells[0].data)
    square.cell_data["medit:ref"] = [np.full(n, 7)]
    ls = np.linalg.norm(square.points[:, :2], axis=1) - 0.5
    out = mio.discretize_levelset(square, ls, materials=[(7, True, 11, 12)])
    assert set(np.unique(out.cell_data["medit:ref"][0])) == {11, 12}


def _annulus(n=32):
    th = np.linspace(0, 2 * np.pi, n + 1)[:-1]
    outer = np.c_[np.cos(th), np.sin(th)]
    lines = [[i, (i + 1) % n] for i in range(n)] + [
        [n + i, n + (i + 1) % n] for i in range(n)
    ]
    boundary = mio.Mesh(np.vstack([outer, 0.4 * outer]), [("line", np.array(lines))])
    boundary.regions = [
        mio.Region("outer", "cell", np.arange(n)),
        mio.Region("inner", "cell", np.arange(n, 2 * n)),
    ]
    return boundary


@needs_mmg
def test_mesh_domain_meshes_an_annulus_with_its_hole():
    # Without the hole point MMG meshes both regions the circles enclose.
    out, report = mio.mesh_domain(_annulus(), hsiz=0.1, return_report=True)
    assert report["num_subdomains"] == 2 and report["num_holes"] == 0
    out, report = mio.mesh_domain(
        _annulus(), hsiz=0.1, holes=[(0, 0)], return_report=True
    )
    assert report["num_subdomains"] == 1 and report["num_holes"] == 1
    triangles = out.cells[0].data
    corners = out.points[triangles]
    a, b = corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0]
    area = 0.5 * np.abs(a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]).sum()
    assert abs(area - np.pi * (1 - 0.16)) < 0.05
    assert np.linalg.norm(out.points, axis=1).min() > 0.39  # no point in the hole
    kinds = {(r.name, r.kind) for r in out.regions}
    assert {("inner", "side"), ("outer", "side"), ("inner", "cell")} <= kinds
    assert sum(r.name.startswith("subdomain:") for r in out.regions) == 1


@needs_mmg
def test_mesh_domain_refusals():
    with pytest.raises(ValueError, match="line cells only"):
        mio.mesh_domain(_square(), hsiz=0.1)
    with pytest.raises(ValueError, match="outside the domain"):
        mio.mesh_domain(_annulus(), hsiz=0.1, holes=[(5, 5)])


@needs_mmg
def test_the_level_set_and_domain_pipeline_steps_in_both_engines(tmp_path):
    import json

    square = mio.rectangle(
        dims=(8, 8), size=(2, 2), origin=(-1, -1, 0), cell_type="triangle"
    )
    square.point_data["ls"] = np.linalg.norm(square.points[:, :2], axis=1) - 0.5
    src = tmp_path / "square.vtu"
    mio.write(str(src), square, compression=None)
    ann = tmp_path / "annulus.vtu"
    mio.write(str(ann), _annulus(), compression=None)
    runs = [
        (
            src,
            {"Op": "DiscretizeLevelset", "LevelSet": "ls", "Materials": [0, 1, 3, 2]},
        ),
        (ann, {"Op": "MeshDomain", "Hsiz": 0.1, "Holes": [0, 0]}),
    ]
    for path, op in runs:
        settings = {
            "Version": 1,
            "Input": {"Path": str(path)},
            "Operations": [op],
            "Output": {"Path": str(tmp_path / "py.vtu")},
        }
        step = mio.run_pipeline(settings)["steps"][0]
        assert step["LowFailure"] == 0.0 and step["NumElementsOut"] > 0
        if op["Op"] == "MeshDomain":
            assert step["NumSubdomains"] == 1.0 and step["NumHoles"] == 1.0
        if hasattr(_core, "run_pipeline_json"):
            settings["Output"]["Path"] = str(tmp_path / "cpp.vtu")
            cpp = _core.run_pipeline_json(json.dumps(settings))["steps"][0]
            assert cpp["NumElementsOut"] == step["NumElementsOut"]
    bad = {"Op": "MeshDomain", "Hsiz": 0.1, "Holes": [0]}
    with pytest.raises(ValueError, match="x, y per hole"):
        mio.run_pipeline(
            {
                "Version": 1,
                "Input": {"Path": str(ann)},
                "Operations": [bad],
                "Output": {"Path": str(tmp_path / "bad.vtu")},
            }
        )


@needs_mmg
def test_the_level_set_and_domain_subcommands(tmp_path, capsys):
    from meshioplusplus._cli import main

    ann = str(tmp_path / "annulus.vtu")
    mio.write(ann, _annulus())
    out = str(tmp_path / "out.vtu")
    assert main(["mesh-2d", ann, out, "--hsiz", "0.1", "--hole", "0,0"]) == 0
    assert "1 subdomain(s), 1 hole(s)" in capsys.readouterr().out
    square = _square()
    square.point_data["ls"] = square.points[:, 0] - 0.5
    src = str(tmp_path / "square.vtu")
    mio.write(src, square)
    assert main(["discretize-levelset", src, out, "--level-set", "ls"]) == 0
    assert {"inside", "outside", "interface"} <= {r.name for r in mio.read(out).regions}


def _holed_square():
    """[0, 1]^2 around the hole [0.4, 0.6]^2, meshed; the hole's lines are "hole"."""
    outer = [[0, 0], [1, 0], [1, 1], [0, 1]]
    inner = [[0.4, 0.4], [0.6, 0.4], [0.6, 0.6], [0.4, 0.6]]
    lines = [[i, (i + 1) % 4] for i in range(4)] + [
        [4 + i, 4 + (i + 1) % 4] for i in range(4)
    ]
    boundary = mio.Mesh(np.array(outer + inner, float), [("line", np.array(lines))])
    boundary.regions = [
        mio.Region("outer", "cell", np.arange(4)),
        mio.Region("hole", "cell", np.arange(4, 8)),
    ]
    return mio.mesh_domain(boundary, hsiz=0.04, holes=[(0.5, 0.5)])


def _hole_points(mesh):
    x, y = mesh.points[:, 0], mesh.points[:, 1]
    return (x > 0.39) & (x < 0.61) & (y > 0.39) & (y < 0.61)


needs_elas = pytest.mark.skipif(
    not mio.has_mmg_elas(), reason="needs an MMG built with ELAS"
)


@pytest.mark.skipif(
    _core is None or (mio.has_mmg() and mio.has_mmg_elas()),
    reason="checks a build without MMG or ELAS",
)
def test_move_names_what_is_missing():
    box = mio.box(dims=(2, 2, 2), cell_type="tetra")
    word = "ELAS" if mio.has_mmg() else "MESHIOPLUSPLUS_WITH_MMG"
    with pytest.raises(NotImplementedError, match=word):
        mio.move(box, np.zeros((len(box.points), 3)), "xmax")


@needs_mmg
@needs_elas
def test_move_translates_a_hole_rigidly():
    # Roadmap 3.2.4's probe: the hole moves by +0.1 in x, every triangle
    # stays right side out, and the area is the square's less the hole's.
    mesh = _holed_square()
    d = np.zeros((len(mesh.points), 2))
    d[_hole_points(mesh), 0] = 0.1
    out, report = mio.move(mesh, d, "hole", return_report=True)
    assert report["library"] == "mmg2d" and report["status"] == "success"
    tri = out.cells_dict["triangle"]
    p = out.points[tri]
    a, b = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]
    signed = 0.5 * (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0])
    assert signed.min() > 0 and np.isclose(signed.sum(), 0.96)
    # the hole's boundary is now [0.5, 0.7] x [0.4, 0.6]
    lines = out.cells_dict["line"]
    edge_points = out.points[np.unique(lines)]
    inside = (edge_points[:, 0] > 0.05) & (edge_points[:, 0] < 0.95)
    inside &= (edge_points[:, 1] > 0.05) & (edge_points[:, 1] < 0.95)
    hole = edge_points[inside]
    np.testing.assert_allclose(hole[:, 0].min(), 0.5, atol=1e-9)
    np.testing.assert_allclose(hole[:, 0].max(), 0.7, atol=1e-9)
    names = {(r.name, r.kind) for r in out.regions}
    assert {("hole", "cell"), ("hole", "side"), ("outer", "side")} <= names


@needs_mmg
@needs_elas
def test_move_pushes_a_face_of_a_box():
    box = mio.box(dims=(6, 6, 6), cell_type="tetra")
    d = np.zeros((len(box.points), 3))
    d[np.isclose(box.points[:, 0], 1.0), 0] = 0.1
    out, report = mio.move(box, d, ["xmax"], return_report=True)
    assert report["library"] == "mmg3d"
    t = out.cells_dict["tetra"]
    p = out.points[t]
    vol = np.einsum(
        "ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0]
    )
    assert vol.min() > 0
    # the face's centre reaches 1.1; its rim, on the fixed faces, stays
    assert np.isclose(out.points[:, 0].max(), 1.1)
    assert 1.0 < vol.sum() / 6 < 1.1


@needs_mmg
@needs_elas
def test_move_refusals_and_its_pipeline_step(tmp_path):
    import json

    mesh = _holed_square()
    d = np.zeros((len(mesh.points), 2))
    with pytest.raises(ValueError, match="no cell or side region named 'nope'"):
        mio.move(mesh, d, "nope")
    with pytest.raises(ValueError, match="2 components"):
        mio.move(mesh, np.zeros(len(mesh.points)), "hole")
    with pytest.raises(ValueError, match="not a surface"):
        mio.move(mio.icosa_sphere(subdivisions=1), np.zeros((42, 3)), "x")
    d[_hole_points(mesh), 0] = 0.1
    mesh.point_data["d"] = d
    src = tmp_path / "holed.vtu"
    mio.write(str(src), mesh, compression=None)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "Move", "Displacement": "d", "Moving": ["hole"]}],
        "Output": {"Path": str(tmp_path / "py.vtu")},
    }
    step = mio.run_pipeline(settings)["steps"][0]
    assert step["op"] == "Move" and step["LowFailure"] == 0.0
    if hasattr(_core, "run_pipeline_json"):
        settings["Output"]["Path"] = str(tmp_path / "cpp.vtu")
        cpp = _core.run_pipeline_json(json.dumps(settings))["steps"][0]
        assert cpp["NumElementsOut"] == step["NumElementsOut"]


@needs_mmg
@needs_elas
def test_the_move_subcommand(tmp_path, capsys):
    from meshioplusplus._cli import main

    mesh = _holed_square()
    d = np.zeros((len(mesh.points), 2))
    d[_hole_points(mesh), 0] = 0.1
    mesh.point_data["d"] = d
    src, out = str(tmp_path / "holed.vtu"), str(tmp_path / "moved.vtu")
    mio.write(src, mesh)
    assert main(["move", src, out, "--displacement", "d", "--moving", "hole"]) == 0
    assert "moved with mmg2d (success)" in capsys.readouterr().out
