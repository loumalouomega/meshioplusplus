"""fit_cad's freeform modes (v17.12.0): options, report keys, per-patch
records, the facet-area measure and the freeform fixtures.

The legacy default must leave every report key as it was; the other modes add
their keys and the ``fit_cad:face_table`` field data. The native backend runs
in every build, so most of these do not need OpenCASCADE.
"""

import json

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _fit_cad
from meshioplusplus._cli._main import main as cli_main

from .analytic_surfaces import (
    box_skin,
    cc_cage,
    double_torus,
    g1_loaf,
    gaussian_bump_box,
    gaussian_noise,
    hood_box,
    monkey_saddle_box,
    organic_torus,
    s_band_box,
    sphere_with_handle,
)

HAS_OCCT = _fit_cad.has_occt()
needs_occt = pytest.mark.skipif(
    not HAS_OCCT, reason="core built without -DMESHIOPLUSPLUS_WITH_OCCT=ON"
)
try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python install
    _core = None
needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")

# The report keys a default (legacy) call has always returned.
_LEGACY_KEYS = {"valid", "backend", "closed"}
_LEGACY_KEYS |= set(_fit_cad._REPORT_FLOATS) | set(_fit_cad._REPORT_INTS)
_MODE_KEYS = {"freeform_mode", "brep_checked", "faces"}
_MODE_KEYS |= set(_fit_cad._MODE_INTS) | set(_fit_cad._MODE_FLOATS)


def _box():
    return mio.fit_primitives(mio.segment(box_skin(3)))


def _closed_and_wound(mesh):
    half = []
    for block in mesh.cells:
        for row in block.data.tolist():
            k = len(row)
            half += [(row[j], row[(j + 1) % k]) for j in range(k)]
    edges = set(half)
    return len(edges) == len(half) and all((b, a) in edges for a, b in half)


@needs_core
def test_legacy_report_keys_are_unchanged():
    _, report = mio.fit_cad(_box(), backend="native", return_report=True)
    assert set(report) == _LEGACY_KEYS


@needs_core
def test_mode_adds_its_keys_and_the_face_table():
    out, report = mio.fit_cad(
        _box(), backend="native", freeform_mode="adaptive", return_report=True
    )
    assert set(report) == _LEGACY_KEYS | _MODE_KEYS
    assert report["freeform_mode"] == "adaptive"
    assert report["brep_checked"] is False
    assert report["facet_area_fraction"] == 0.0
    assert len(report["faces"]) == 6
    for row in report["faces"]:
        assert row["kind"] == 1
        assert row["face"] >= 0
        assert row["demotion_reason"] == "none"
        assert row["param_method"] == "none"
    table = np.asarray(out.field_data["fit_cad:face_table"])
    assert table.shape == (6, 13)


@needs_core
def test_native_modes_return_the_solid_as_nurbs(tmp_path):
    # Outside the legacy mode the native solid is trimmed NURBS faces too,
    # so Kratos CAD JSON needs no CAD kernel; the legacy mode leaves it empty.
    out, report = mio.fit_cad(
        _box(), backend="native", freeform_mode="adaptive", return_report=True
    )
    brep = out.info["brep"]
    assert len(brep["faces"]) == report["num_faces"] == 6
    assert all(f["trimmed"] for f in brep["faces"])
    assert len(brep["edges"]) == 12
    mio.write(tmp_path / "box.cad.json", out, file_format="cad_json")
    back = mio.read(tmp_path / "box.cad.json").info["brep"]
    assert len(back["faces"]) == 6
    legacy = mio.fit_cad(_box(), backend="native")
    assert not legacy.info["brep"] or not legacy.info["brep"].get("faces")


@needs_core
def test_native_network_brep_is_written_alike_by_both_engines(tmp_path):
    # The native patch network's NURBS B-rep: the core and the Python CAD
    # JSON writers give the same bytes, and a read and rewrite gives them again.
    from meshioplusplus.cad_json import _cad_json as ref

    out = mio.fit_cad(
        mio.fit_primitives(mio.segment(monkey_saddle_box(16)), tolerance=2e-3),
        backend="native",
        freeform_mode="network",
    )
    assert len(out.info["brep"]["faces"]) > 6
    mio.cad_json.write(tmp_path / "core.cad.json", out)
    ref.write(tmp_path / "py.cad.json", out)
    core = (tmp_path / "core.cad.json").read_bytes()
    assert core == (tmp_path / "py.cad.json").read_bytes()
    back = mio.cad_json.read(tmp_path / "core.cad.json")
    mio.cad_json.write(tmp_path / "again.cad.json", back)
    assert (tmp_path / "again.cad.json").read_bytes() == core


@needs_core
def test_demotion_reason_is_recorded():
    from .analytic_surfaces import jitter

    skin = mio.fit_primitives(mio.segment(jitter(box_skin(3), 1e-4)))
    _, report = mio.fit_cad(
        skin,
        backend="native",
        tolerance=1e-6,
        freeform_mode="adaptive",
        return_report=True,
    )
    assert report["num_demoted"] > 0
    reasons = {row["demotion_reason"] for row in report["faces"]}
    assert "hausdorff_exceeded" in reasons
    assert report["facet_area_fraction"] > 0.0


@pytest.mark.parametrize(
    "key, value",
    [
        ("freeform_mode", "fast"),
        ("freeform_parameterization", "conformal"),
        ("extraordinary_vertex", "g2"),
        ("local_refine", "tsplines"),
        ("freeform_g1_tolerance_deg", 0.0),
        ("freeform_max_control_points", 8),
        ("layout_target_patch_size", -0.1),
        ("layout_max_valence_deviation", 4),
        ("resegment_rounds", -1),
    ],
)
def test_options_are_validated(key, value):
    with pytest.raises(ValueError, match=key):
        mio.fit_cad(_box(), backend="native", **{key: value})


@needs_core
def test_cli_freeform_mode_json(tmp_path, capsys):
    src = tmp_path / "box.vtu"
    mio.write(src, _box())
    rc = cli_main(
        [
            "to-cad",
            str(src),
            str(tmp_path / "box.step"),
            "--backend",
            "native",
            "--freeform-mode",
            "adaptive",
            "--json",
        ]
    )
    assert rc == 0
    report = json.loads(capsys.readouterr().out)
    assert report["freeform_mode"] == "adaptive"
    assert report["facet_area_fraction"] == 0.0
    assert len(report["faces"]) == 6


@needs_core
def test_pipeline_freeform_mode_counters(tmp_path):
    src = tmp_path / "box.vtu"
    mio.write(src, _box())
    doc = {
        "Input": {"Path": str(src)},
        "Operations": [{"Op": "FitCad", "Backend": "native", "FreeformMode": "auto"}],
        "Output": {"Path": str(tmp_path / "out.vtu")},
    }
    step = mio.run_pipeline(doc)["steps"][0]
    assert step["FacetAreaFraction"] == 0.0
    assert step["NumNetworkPatches"] == 0
    doc["Operations"][0].pop("FreeformMode")
    step = mio.run_pipeline(doc)["steps"][0]
    assert "FacetAreaFraction" not in step


@needs_core
def test_freeform_fixtures_are_closed_skins():
    for mesh in (
        monkey_saddle_box(8),
        gaussian_bump_box(8),
        s_band_box(8),
        hood_box(8, warp=True),
        organic_torus(16, 8),
        cc_cage(3, 2),
        cc_cage(5, 2),
        cc_cage(6, 2),
        g1_loaf(8, 6),
    ):
        assert _closed_and_wound(mesh)


@needs_core
def test_implicit_fixtures_are_closed_skins():
    assert _closed_and_wound(sphere_with_handle(16))
    assert _closed_and_wound(double_torus(20))


def _fnv1a(data: bytes) -> int:
    h = 0xCBF29CE484222325
    for byte in data:
        h = ((h ^ byte) * 0x100000001B3) & ((1 << 64) - 1)
    return h


# The trig-free fixtures' points, bit for bit the C++ ones: the same digests
# are pinned in tests/cpp/test_fit_cad_modes.cpp (FitCadModes.FixturesMatchPython).
_POINT_DIGESTS = {
    "monkey_saddle_box": 0x06C59D1FC5A203C4,
    "s_band_box": 0x04808850AB51BF6D,
    "hood_box": 0x4B1FE1E9066260D6,
    "cc_cage3": 0x85F0B2DD4C5350E1,
    "noisy_hood": 0x6BE481CC0438D8C3,
}


def _fixture_points():
    return {
        "monkey_saddle_box": monkey_saddle_box(24),
        "s_band_box": s_band_box(24),
        "hood_box": hood_box(24, warp=True),
        "cc_cage3": cc_cage(3),
        "noisy_hood": gaussian_noise(hood_box(12), 1e-3, 11),
    }


def test_fixtures_match_the_cpp_ones():
    for name, mesh in _fixture_points().items():
        pts = np.ascontiguousarray(mesh.points, dtype="<f8")
        assert _fnv1a(pts.tobytes()) == _POINT_DIGESTS[name], name


def _fitted_abs(skin, tol, min_patch_cells=1):
    seg = mio.segment(skin, min_patch_cells=min_patch_cells)
    return mio.fit_primitives(seg, tolerance=tol, relative=False)


def _diag(mesh):
    pts = np.asarray(mesh.points, dtype=float)
    return float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))


@needs_core
@pytest.mark.parametrize("make", [hood_box, monkey_saddle_box, s_band_box])
def test_native_adaptive_writes_bspline_faces(tmp_path, make):
    skin = make(16)
    tol = 2e-3 * _diag(skin)
    fitted = _fitted_abs(skin, tol)
    _, legacy = mio.fit_cad(fitted, backend="native", tolerance=tol, return_report=True)
    assert legacy["num_bspline"] == 0
    path = tmp_path / "adaptive.step"
    out, report = mio.fit_cad(
        fitted,
        path,
        backend="native",
        tolerance=tol,
        freeform_mode="adaptive",
        return_report=True,
    )
    assert "B_SPLINE_SURFACE_WITH_KNOTS" in path.read_text()
    kinds = [row["kind"] for row in report["faces"]]
    assert 6 in kinds
    assert report["max_freeform_distance"] <= tol
    # The kernel-free writer's bytes are reproducible.
    (tmp_path / "again").mkdir()
    path2 = tmp_path / "again" / "adaptive.step"
    mio.fit_cad(
        fitted, path2, backend="native", tolerance=tol, freeform_mode="adaptive"
    )
    assert path.read_bytes() == path2.read_bytes()


@needs_core
def test_native_legacy_still_facets_the_freeform_top():
    skin = hood_box(16)
    tol = 2e-3 * _diag(skin)
    _, report = mio.fit_cad(
        _fitted_abs(skin, tol), backend="native", tolerance=tol, return_report=True
    )
    assert report["num_bspline"] == 0
    assert report["num_faceted_faces"] >= 1


def _merged_side(fitted):
    """The fitted box with its planar side at the smallest x merged into its
    freeform (type 0) patch, as a segmentation that missed the creases leaves
    it; the side's fit-table row dropped (tests/cpp/fit_cad_fixtures.hpp)."""
    pts = np.asarray(fitted.points)
    lo = pts[:, 0].min()
    ids = np.asarray(fitted.field_data["fit:patch_id"])
    types = np.asarray(fitted.field_data["fit:patch_type"])
    freeform = int(ids[types == 0][0])
    side = None
    labels = []
    for block, data in zip(fitted.cells, fitted.cell_data["segment:patch"]):
        lab = np.asarray(data, dtype=np.int64).copy()
        on = np.all(pts[block.data][:, :, 0] == lo, axis=1)
        if on.any():
            side = int(lab[on][0])
        lab[on] = freeform
        labels.append(lab)
    out = mio.Mesh(
        fitted.points,
        fitted.cells,
        cell_data={**fitted.cell_data, "segment:patch": labels},
        field_data=dict(fitted.field_data),
    )
    keep = ids != side
    for key in (
        "fit:patch_id",
        "fit:patch_type",
        "fit:patch_params",
        "fit:patch_tolerance",
    ):
        out.field_data[key] = np.asarray(fitted.field_data[key])[keep]
    return out


@needs_core
def test_resegmentation_cuts_a_missed_crease():
    skin = hood_box(16)
    tol = 2e-3 * _diag(skin)
    merged = _merged_side(_fitted_abs(skin, tol))
    kw = dict(
        backend="native", tolerance=tol, freeform_mode="adaptive", return_report=True
    )
    _, off = mio.fit_cad(merged, resegment_rounds=0, **kw)
    _, on = mio.fit_cad(merged, resegment_rounds=2, **kw)
    assert off["num_resegmented"] == 0
    assert on["num_resegmented"] > 0
    assert off["facet_area_fraction"] > 0.04
    assert on["facet_area_fraction"] < 0.01
    assert on["num_faces"] < off["num_faces"]
