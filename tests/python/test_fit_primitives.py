"""Tests for ``fit_primitives`` (roadmap §6.3.2): the Python face, its numpy
twin against the compiled core, the CLI verb and the pipeline step. The C++
behaviour is pinned in ``tests/cpp/test_fit_primitives.cpp``."""

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _fit_primitives as fp

from .analytic_surfaces import (
    box_skin,
    cone,
    cylinder,
    jitter,
    sphere_band,
    torus,
    turned_part,
)

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")

_DEFAULTS = {
    "patch_data": "",
    "tolerance": 1e-3,
    "relative": True,
    "min_inlier_fraction": 0.9,
    "types": fp.types_mask(fp.TYPES[1:]),
    "max_iterations": 50,
}


def _whole(mesh, **kw):
    out = mio.fit_primitives(mesh, patch_data="", **kw)
    fd = out.field_data
    return (
        int(fd["fit:patch_type"][0]),
        fd["fit:patch_params"][0],
        fd["fit:patch_rms"][0],
    )


def test_box_gives_six_outward_planes():
    out, report = mio.fit_primitives(mio.segment(box_skin(3)), return_report=True)
    assert report["num_patches"] == 6 and report["num_plane"] == 6
    params = out.field_data["fit:patch_params"]
    normals, points = params[:, 3:6], params[:, 0:3]
    np.testing.assert_allclose(np.linalg.norm(normals, axis=1), 1.0, atol=1e-12)
    np.testing.assert_allclose(
        np.einsum("ij,ij->i", normals, points - 0.5), 0.5, atol=1e-12
    )
    assert np.all(out.cell_data["fit:type"][0] == 1)


@pytest.mark.parametrize(
    "mesh,kind,check",
    [
        (cylinder(1.25), 3, lambda p: abs(p[6] - 1.25) < 1e-9),
        (cone(20.0), 4, lambda p: abs(p[8] - 20.0) < 1e-7 and abs(p[:3]).max() < 1e-7),
        (sphere_band(1.5), 2, lambda p: abs(p[6] - 1.5) < 1e-9),
        (
            torus(3.0, 1.0),
            5,
            lambda p: abs(p[6] - 3.0) < 1e-9 and abs(p[7] - 1.0) < 1e-9,
        ),
    ],
    ids=["cylinder", "cone", "sphere", "torus"],
)
def test_recovers_each_primitive(mesh, kind, check):
    t, params, rms = _whole(mesh)
    assert t == kind
    assert check(params)
    assert rms < 1e-9


def test_noise_raises_the_residual_not_the_type():
    clean = _whole(cylinder())
    t, params, rms = _whole(jitter(cylinder(), 1e-3), tolerance=5e-3)
    assert t == 3 and rms > clean[2] and rms < 2e-3
    assert abs(params[6] - 1.0) < 2e-3


def _saddle():
    n = 8
    g = np.linspace(-1.0, 1.0, n + 1)
    x, y = np.meshgrid(g, g)
    pts = np.column_stack([x.ravel(), y.ravel(), 0.5 * (x * x - y * y).ravel()])
    quads = [
        [
            j * (n + 1) + i,
            j * (n + 1) + i + 1,
            (j + 1) * (n + 1) + i + 1,
            (j + 1) * (n + 1) + i,
        ]
        for j in range(n)
        for i in range(n)
    ]
    return mio.Mesh(pts, [("quad", np.array(quads))])


def test_a_saddle_is_reported_unfitted_with_its_residual():
    out, report = mio.fit_primitives(
        _saddle(), patch_data="", tolerance=1e-4, return_report=True
    )
    assert report["num_unfitted"] == 1
    assert out.field_data["fit:patch_type"][0] == 0
    assert out.field_data["fit:patch_rms"][0] > 1e-4
    assert np.all(np.isnan(out.field_data["fit:patch_params"][0]))


def test_types_restrict_the_candidates():
    assert _whole(cylinder(), types=["plane"])[0] == 0
    assert _whole(cylinder(), types=["plane", "cylinder"])[0] == 3


def test_bad_arguments_are_refused():
    with pytest.raises(ValueError, match="no cell data"):
        mio.fit_primitives(cylinder())
    with pytest.raises(ValueError, match="tolerance"):
        mio.fit_primitives(cylinder(), patch_data="", tolerance=0.0)
    with pytest.raises(ValueError, match="unknown primitive"):
        mio.fit_primitives(cylinder(), patch_data="", types=["blob"])
    with pytest.raises(ValueError, match="min_inlier_fraction"):
        mio.fit_primitives(cylinder(), patch_data="", min_inlier_fraction=2.0)


@needs_core
@pytest.mark.parametrize(
    "mesh",
    [
        mio.segment(box_skin(3)),
        cylinder(),
        cone(),
        sphere_band(),
        torus(),
        jitter(torus(), 2e-3),
        _saddle(),
    ],
    ids=["box", "cylinder", "cone", "sphere", "torus", "noisy-torus", "saddle"],
)
def test_twin_matches_the_core(mesh):
    name = "segment:patch" if "segment:patch" in mesh.cell_data else ""
    core, core_report = mio.fit_primitives(mesh, patch_data=name, return_report=True)
    twin, twin_report = fp._fit_primitives_py(mesh, dict(_DEFAULTS, patch_data=name))
    assert core_report == twin_report
    for key in ("fit:patch_id", "fit:patch_type"):
        np.testing.assert_array_equal(core.field_data[key], twin.field_data[key])
    for key in (
        "fit:patch_params",
        "fit:patch_rms",
        "fit:patch_max",
        "fit:patch_tolerance",
    ):
        np.testing.assert_allclose(
            core.field_data[key], twin.field_data[key], rtol=1e-6, atol=1e-9
        )
    np.testing.assert_array_equal(
        core.cell_data["fit:type"][0], twin.cell_data["fit:type"][0]
    )


def test_cli_verb(tmp_path, capsys):
    import json

    src = str(tmp_path / "cyl.vtu")
    mio.write(src, cylinder())
    out = str(tmp_path / "fit.vtu")
    args = ["fit-primitives", src, out, "--patch-data", "", "--types", "plane,cylinder"]
    assert mio._cli.main(args + ["--json"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["num_cylinder"] == 1
    assert "fit:patch_params" in mio.read(out).field_data


def test_a_noisy_turned_part_recovers_every_type():
    seg = mio.segment(
        turned_part(noise=1e-3), curvature_tolerance=1.0, min_patch_cells=20
    )
    out = mio.fit_primitives(seg, tolerance=2e-3)
    types = sorted(out.field_data["fit:patch_type"].tolist())
    assert types == [1, 1, 2, 3, 4, 5]  # base and shoulder planes, then one of each
    params = dict(
        zip(out.field_data["fit:patch_type"], out.field_data["fit:patch_params"])
    )
    assert abs(params[3][6] - 2.0) < 5e-3  # the wall's radius
    assert abs(params[5][6] - 1.5) < 5e-3 and abs(params[5][7] - 0.5) < 5e-3  # fillet
    assert abs(params[2][6] - 1.5) < 5e-2  # the shallow dome's radius
