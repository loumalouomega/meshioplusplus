"""Tests for the isogeometric layer: ``iga_evaluate``, ``iga_refine``, ``iga_sample``.

The defining property is that refinement is exact: knot insertion and degree
elevation change the control net and the knot vectors but not the geometry, the
parametrisation, the NURBS fields or a trimmed face's trimming loops. It is
checked on every NURBS fixture of the MFEM reader (rational and polynomial,
curves, surfaces and volumes, with a field) and every CAD JSON fixture (a
rational cylinder, a trimmed plate, two coupled patches), against closed-form
circles, and against central differences for the derivatives. The two engines
(the C++ core and the numpy twin) must give the same bits. The C++ behaviour is
pinned in ``tests/cpp/test_iga.cpp``.
"""

import pathlib

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _iga

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")

MESHES = pathlib.Path(__file__).parent / "meshes"
NURBS_DIR = MESHES / "mfem" / "nurbs"
NURBS = [
    "ball-nurbs",
    "pipe-nurbs",
    "square-disc-nurbs-patch",
    "nurbs-segments2d-patches",
    "beam-quad-nurbs-sf",
    "cube-nurbs",
]
CAD = sorted(p.name for p in (MESHES / "cad_json").glob("*.cad.json"))

# A model with MFEM's patch topology refines every direction alike (its patches
# share knot vectors along differently oriented directions); CAD faces and a
# hand-built model may refine each direction on its own.
VARIANTS = [
    dict(subdivisions=2),
    dict(elevate=1),
    dict(subdivisions=(2, 3, 2), elevate=(2, 1, 1)),
]
NURBS_VARIANTS = [
    dict(subdivisions=2),
    dict(elevate=1),
    dict(subdivisions=3, elevate=2),
]


def _nurbs(name):
    return mp.mfem.read(
        NURBS_DIR / f"{name}.mesh", {"u": str(NURBS_DIR / f"{name}.u.gf")}
    )


def _cad(name):
    return mp.read(MESHES / "cad_json" / name)


def _same_tree(a, b, where="model"):
    """Two models, key by key: equal values, dtypes and shapes."""
    if isinstance(a, dict):
        assert sorted(a) == sorted(b), where
        for k in a:
            _same_tree(a[k], b[k], f"{where}[{k!r}]")
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b), where
        for i, (x, y) in enumerate(zip(a, b)):
            _same_tree(x, y, f"{where}[{i}]")
    elif isinstance(a, np.ndarray):
        assert isinstance(b, np.ndarray) and a.dtype == b.dtype, where
        np.testing.assert_array_equal(a, b, err_msg=where)
    else:
        assert a == b, where


def _range(patch, d):
    p = int(patch["degrees"][d])
    k = patch["knots"][d]
    return float(k[p]), float(k[len(k) - p - 1])


def _params(patch, n, seed=1):
    rng = np.random.default_rng(seed)
    dim = len(patch["degrees"])
    cols = [rng.uniform(*_range(patch, d), size=n) for d in range(dim)]
    return np.stack(cols, axis=1)


def _engines(kind):
    """The two engines of a refinement, as ``(label, fn(model, sub, ele))``
    returning the refined model."""

    def py(model, sub, ele):
        if kind == "nurbs":
            return _iga._refine_nurbs(_iga._complete_nurbs(model), sub, ele)
        return _iga._refine_brep(model, sub, ele)

    out = [("python", py)]
    if _core is not None:

        def core(model, sub, ele):
            if kind == "nurbs":
                return _core.iga_refine_nurbs(model, sub, ele, 0).info["nurbs"]
            return _core.iga_refine_brep(model, sub, ele, 0, 16, True, True).info[
                "brep"
            ]

        out.append(("core", core))
    return out


# ---------------------------------------------------------------------------
# Refinement is exact
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("variant", NURBS_VARIANTS, ids=["h", "p", "hp"])
@pytest.mark.parametrize("name", NURBS)
def test_nurbs_refinement_keeps_geometry_and_fields(name, variant):
    info = _nurbs(name).info["nurbs"]
    sub = _iga._triple(variant.get("subdivisions", 1), "s", 1)
    ele = _iga._triple(variant.get("elevate", 0), "e", 0)
    for label, refine in _engines("nurbs"):
        refined = refine(info, sub, ele)
        assert len(refined["patches"]) == len(info["patches"]), label
        for pi, (p, q) in enumerate(zip(info["patches"], refined["patches"])):
            u = _params(p, 60, seed=pi)
            a = _iga._evaluate_patch(p, _iga._is_rational(p), u, 1, "p")
            b = _iga._evaluate_patch(q, _iga._is_rational(q), u, 1, "q")
            np.testing.assert_allclose(a["points"], b["points"], rtol=0, atol=1e-11)
            np.testing.assert_allclose(
                a["derivatives"], b["derivatives"], rtol=1e-10, atol=1e-10
            )
            for fa, fb in zip(a["fields"], b["fields"]):
                np.testing.assert_allclose(fa, fb, rtol=0, atol=1e-11)
            for d in range(len(p["degrees"])):
                assert q["degrees"][d] == p["degrees"][d] + ele[d]
                assert _range(q, d) == _range(p, d)


@pytest.mark.parametrize("variant", VARIANTS, ids=["h", "p", "hp"])
@pytest.mark.parametrize("name", CAD)
def test_brep_refinement_keeps_faces_and_trims(name, variant):
    brep = _cad(name).info["brep"]
    sub = _iga._triple(variant.get("subdivisions", 1), "s", 1)
    ele = _iga._triple(variant.get("elevate", 0), "e", 0)
    for label, refine in _engines("brep"):
        refined = refine(brep, sub, ele)
        # the trimming loops, edges and vertices are the parametrisation's: untouched
        for fa, fb in zip(brep["faces"], refined["faces"]):
            _same_tree(fa["loops"], fb["loops"])
            assert fb["control_point_ids"].size == 0
        _same_tree(brep["edges"], refined["edges"])
        _same_tree(brep["vertices"], refined["vertices"])
        for fi, face in enumerate(brep["faces"]):
            p = brep["patches"][face["patch"]]
            u = _params(p, 50, seed=fi)
            a = _iga._evaluate_brep_py(brep, fi, u, 1, 16)
            b = _iga._evaluate_brep_py(refined, fi, u, 1, 16)
            np.testing.assert_allclose(a["points"], b["points"], rtol=0, atol=1e-11)
            np.testing.assert_allclose(
                a["derivatives"], b["derivatives"], rtol=1e-10, atol=1e-10
            )
            np.testing.assert_array_equal(a["inside"], b["inside"])


# The numpy MFEM sampler is slow on the big volumes: compare whole meshes on the
# small fixtures, the refined models and every evaluation on all of them.
SMALL = [
    "square-disc-nurbs-patch",
    "nurbs-segments2d-patches",
    "beam-quad-nurbs-sf",
    "cube-nurbs",
]


@needs_core
@pytest.mark.parametrize("variant", NURBS_VARIANTS, ids=["h", "p", "hp"])
@pytest.mark.parametrize("name", NURBS)
def test_nurbs_engines_are_bit_identical(name, variant):
    mesh = _nurbs(name)
    info = mesh.info["nurbs"]
    sub = _iga._triple(variant.get("subdivisions", 1), "s", 1)
    ele = _iga._triple(variant.get("elevate", 0), "e", 0)
    core = _core.iga_refine_nurbs(info, sub, ele, 0)
    refined = _iga._refine_nurbs(info, sub, ele)
    # the core's side channel is re-gathered from MFEM's global numbering, where a
    # control point two patches share is one value: within rounding of the patches'
    assert core.info["nurbs"]["field_names"] == refined["field_names"]
    _same_tree(core.info["nurbs"]["topology"], refined["topology"])
    for cp, rp in zip(core.info["nurbs"]["patches"], refined["patches"]):
        np.testing.assert_array_equal(cp["degrees"], rp["degrees"])
        for kc, kr in zip(cp["knots"], rp["knots"]):
            np.testing.assert_allclose(kc, kr, rtol=0, atol=1e-15)
        np.testing.assert_allclose(
            cp["control_points"], rp["control_points"], rtol=0, atol=1e-14
        )
        np.testing.assert_allclose(cp["weights"], rp["weights"], rtol=0, atol=1e-14)
        for fc, fr in zip(cp["fields"], rp["fields"]):
            np.testing.assert_allclose(fc, fr, rtol=0, atol=1e-14)
    if name in SMALL:
        twin, mapped = _iga._sample_nurbs(refined, 0)
        _same_tree(core.info["nurbs"], mapped)
        np.testing.assert_array_equal(core.points, twin.points)
        assert [c.type for c in core.cells] == [c.type for c in twin.cells]
        for x, y in zip(core.cells, twin.cells):
            np.testing.assert_array_equal(x.data, y.data)
        for k in twin.point_data:
            np.testing.assert_array_equal(core.point_data[k], twin.point_data[k])
    # the generic sampler, both engines
    gen_core = _core.iga_sample_nurbs(refined, 1)
    gen_twin, gen_mapped = _iga._sample_nurbs(refined, 1)
    _same_tree(gen_core.info["nurbs"], gen_mapped)
    np.testing.assert_array_equal(gen_core.points, gen_twin.points)
    np.testing.assert_array_equal(gen_core.cells[0].data, gen_twin.cells[0].data)
    for k in gen_twin.point_data:
        np.testing.assert_array_equal(gen_core.point_data[k], gen_twin.point_data[k])
    for pi, p in enumerate(info["patches"]):
        u = _params(p, 30, seed=pi)
        a = _core.iga_evaluate_nurbs(info, pi, u, 1)
        b = _iga._evaluate_patch(p, _iga._is_rational(p), u, 1, "p")
        for k in ("points", "derivatives", "normals"):
            if k in a or k in b:
                np.testing.assert_array_equal(a[k], b[k], err_msg=k)
        for fa, fb in zip(a["fields"], b["fields"]):
            np.testing.assert_array_equal(fa, fb)


@needs_core
@pytest.mark.parametrize("variant", VARIANTS, ids=["h", "p", "hp"])
@pytest.mark.parametrize("name", CAD)
def test_brep_engines_are_bit_identical(name, variant):
    mesh = _cad(name)
    brep = mesh.info["brep"]
    sub = _iga._triple(variant.get("subdivisions", 1), "s", 1)
    ele = _iga._triple(variant.get("elevate", 0), "e", 0)
    core = _core.iga_refine_brep(brep, sub, ele, 0, 16, True, True)
    refined = _iga._refine_brep(brep, sub, ele)
    twin = _iga._sample_brep_py(refined, mesh, 0, 16, True, True)
    _same_tree(core.info["brep"], twin.info["brep"])
    np.testing.assert_array_equal(core.points, twin.points)
    for fi, face in enumerate(brep["faces"]):
        u = _params(brep["patches"][face["patch"]], 30, seed=fi)
        a = _core.iga_evaluate_brep(brep, fi, u, 1, 16)
        b = _iga._evaluate_brep_py(brep, fi, u, 1, 16)
        for k in ("points", "derivatives", "normals", "inside"):
            np.testing.assert_array_equal(a[k], b[k], err_msg=k)


def test_knot_multiplicities_follow_the_elevation():
    # one quadratic span pair with a C1 interior knot: elevated by 2 it is a
    # cubic-plus-one with the knot's multiplicity raised by 2 (C1 kept: degree
    # 4, multiplicity 3 -> continuity 4 - 3 = 1), then split in 3 by simple knots
    patch = {
        "degrees": [2],
        "knots": [[0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0]],
        "control_points": [[0.0, 0.0], [1.0, 2.0], [2.0, 0.0], [3.0, 2.0]],
    }
    model = {"patches": [patch]}
    mesh = _iga.iga_refine(model, subdivisions=3, elevate=2)
    q = mesh.info["nurbs"]["patches"][0]
    assert q["degrees"].tolist() == [4]
    np.testing.assert_allclose(
        q["knots"][0],
        [0, 0, 0, 0, 0, 1 / 3, 2 / 3, 1, 1, 1, 4 / 3, 5 / 3, 2, 2, 2, 2, 2],
        rtol=0,
        atol=1e-15,
    )
    assert len(q["control_points"]) == len(q["knots"][0]) - 4 - 1
    u = np.linspace(0, 2, 41)[:, None]
    a = _iga.iga_evaluate(model, u)["points"]
    b = _iga.iga_evaluate(mesh, u)["points"]
    np.testing.assert_allclose(a, b, atol=1e-13)


def test_circle_stays_exact_under_every_refinement():
    s = np.sqrt(0.5)
    arc = {
        "degrees": [2],
        "knots": [[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]],
        "control_points": [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        "weights": [1.0, s, 1.0],
    }
    u = np.linspace(0.0, 1.0, 37)
    for sub, ele in [(1, 0), (5, 0), (1, 3), (4, 2)]:
        out = _iga.iga_refine(
            {"patches": [arc]}, subdivisions=sub, elevate=ele, sample_subdivisions=1
        )
        r = np.hypot(*_iga.iga_evaluate(out, u[:, None])["points"].T)
        np.testing.assert_allclose(r, 1.0, atol=1e-14, rtol=0)
        # the sampled nodes are on the circle too
        np.testing.assert_allclose(np.hypot(*out.points[:, :2].T), 1.0, atol=1e-14)


def test_sphere_octant_normals_point_outwards():
    # a rational bi-quadratic sphere patch: unit-length position, normal parallel to it
    s = np.sqrt(0.5)
    a = [(1, 0, 0), (1, 1, 0), (0, 1, 0)]
    wa = [1.0, s, 1.0]
    # the octant is the product of the quarter circle in the xy plane with a
    # vertical quarter circle (radius, height) from the pole to the equator
    cps, w = [], []
    vert = [(0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]  # (radius, height) control points
    wv = [1.0, s, 1.0]
    for (rad, h), wvk in zip(vert, wv):
        for (x, y, _), wak in zip(a, wa):
            cps.append([rad * x, rad * y, h])
            w.append(wak * wvk)
    patch = {
        "degrees": [2, 2],
        "knots": [[0, 0, 0, 1, 1, 1.0], [0, 0, 0, 1, 1, 1.0]],
        "control_points": cps,
        "weights": w,
    }
    model = {"patches": [patch]}
    for sub, ele in [(1, 0), (3, 1)]:
        out = _iga.iga_refine(model, subdivisions=sub, elevate=ele)
        grid = np.stack(
            np.meshgrid(np.linspace(0, 1, 9), np.linspace(0, 1, 9)), -1
        ).reshape(-1, 2)
        res = _iga.iga_evaluate(out, grid, derivatives=1)
        np.testing.assert_allclose(
            np.linalg.norm(res["points"], axis=1), 1.0, atol=1e-13
        )
        mask = (
            np.linalg.norm(
                np.cross(res["derivatives"][:, 0], res["derivatives"][:, 1]), axis=1
            )
            > 1e-9
        )
        np.testing.assert_allclose(
            np.linalg.norm(res["normals"][mask], axis=1), 1.0, atol=1e-13
        )
        dots = np.abs(np.einsum("ij,ij->i", res["normals"][mask], res["points"][mask]))
        np.testing.assert_allclose(dots, 1.0, atol=1e-12)


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name", ["square-disc-nurbs-patch", "pipe-nurbs", "nurbs-segments2d-patches"]
)
def test_derivatives_match_central_differences(name):
    info = _nurbs(name).info["nurbs"]
    for pi, p in enumerate(info["patches"][:2]):
        dim = len(p["degrees"])
        lo = np.array([_range(p, d)[0] for d in range(dim)])
        hi = np.array([_range(p, d)[1] for d in range(dim)])
        u = lo + (hi - lo) * np.random.default_rng(pi).uniform(0.2, 0.8, size=(12, dim))
        got = _iga.iga_evaluate(info, u, patch=pi, derivatives=1)["derivatives"]
        h = 1e-6
        for d in range(dim):
            e = np.zeros(dim)
            e[d] = h * (hi[d] - lo[d])
            fwd = _iga.iga_evaluate(info, u + e, patch=pi)["points"]
            bwd = _iga.iga_evaluate(info, u - e, patch=pi)["points"]
            np.testing.assert_allclose(
                got[:, d], (fwd - bwd) / (2 * e[d]), rtol=1e-6, atol=1e-6
            )


def test_cylindrical_roof_is_exact_and_normals_are_unit():
    mesh = _cad("cylindrical_roof.cad.json")
    brep = mesh.info["brep"]
    u = _params(brep["patches"][0], 200)
    res = _iga.iga_evaluate(mesh, u, derivatives=1)
    # a cylinder of radius 25 about the y axis: the radial distance is constant
    r = np.hypot(res["points"][:, 0], res["points"][:, 2])
    np.testing.assert_allclose(r, 25.0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(res["normals"], axis=1), 1.0, atol=1e-13)
    assert res["inside"].all()  # an untrimmed face contains every point


def test_trimmed_plate_inside_matches_the_clipped_sampling():
    mesh = _cad("plate_with_hole.cad.json")
    brep = mesh.info["brep"]
    patch = brep["patches"][0]
    n = 12
    # the centres of an n x n grid of cells over the parameter domain
    (u0, u1), (v0, v1) = _range(patch, 0), _range(patch, 1)
    gu = u0 + (np.arange(n) + 0.5) * (u1 - u0) / n
    gv = v0 + (np.arange(n) + 0.5) * (v1 - v0) / n
    centres = np.array([[a, b] for b in gv for a in gu])
    inside = _iga.iga_evaluate(mesh, centres)["inside"]
    assert inside.dtype == bool and 0 < inside.sum() < len(inside)
    # the sampler keeps exactly the cells whose parametric centre is inside
    clipped = mp.cad_json.sample(brep, subdivisions=n, edges=False, trim=True)
    full = mp.cad_json.sample(brep, subdivisions=n, edges=False, trim=False)
    assert sum(len(c.data) for c in clipped.cells) == int(inside.sum())
    assert sum(len(c.data) for c in full.cells) == len(inside)


def test_fields_are_evaluated_with_the_geometry():
    mesh = _nurbs("ball-nurbs")
    info = mesh.info["nurbs"]
    # at the sampled nodes of the reader's mesh the field equals the grid function
    sampled = _iga.iga_sample(mesh, subdivisions=0)
    assert "u" in sampled.point_data
    res = _iga.iga_evaluate(info, _params(info["patches"][0], 5), patch=0)
    assert len(res["fields"]) == 1 and res["fields"][0].shape == (5, 1)


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", NURBS)
def test_default_sampling_is_the_mfem_readers(name):
    mesh = _nurbs(name)
    out = _iga.iga_sample(mesh)
    np.testing.assert_array_equal(out.points, mesh.points)
    assert [c.type for c in out.cells] == [c.type for c in mesh.cells]
    for x, y in zip(out.cells, mesh.cells):
        np.testing.assert_array_equal(x.data, y.data)
    for k in mesh.point_data:
        np.testing.assert_array_equal(out.point_data[k], mesh.point_data[k])
    _same_tree(out.info["nurbs"], mesh.info["nurbs"])
    # an identity refinement is the same mesh and the same model
    again = _iga.iga_refine(mesh)
    np.testing.assert_array_equal(again.points, mesh.points)
    _same_tree(again.info["nurbs"], mesh.info["nurbs"])


@pytest.mark.parametrize("name", NURBS)
def test_generic_sampling_types_maps_and_fields(name):
    mesh = _nurbs(name)
    out = _iga.iga_sample(mesh, subdivisions=1)
    info = out.info["nurbs"]
    dim = info["topology"]["dimension"]
    top = max(int(d) for p in info["patches"] for d in p["degrees"])
    family = {
        1: ("line", "line3"),
        2: ("quad", "quad9"),
        3: ("hexahedron", "hexahedron27"),
    }[dim]
    lagrange = "VTK_LAGRANGE_" + {1: "CURVE", 2: "QUADRILATERAL", 3: "HEXAHEDRON"}[dim]
    want = family[0] if top <= 1 else family[1] if top == 2 else lagrange
    assert [c.type for c in out.cells] == [want]
    n = len(out.cells[0].data)
    assert n == len(mesh.cells[0].data)  # one cell per knot span, as the reader
    assert info["cell_patch"].shape == (n,) and info["cell_span"].shape == (n, 3)
    np.testing.assert_array_equal(out.cell_data["iga:patch"][0], info["cell_patch"])
    np.testing.assert_array_equal(out.cell_data["iga:span"][0], info["cell_span"])
    assert info["boundary_cell_patch"].size == 0
    # the same geometry as the reader's, patch by patch (points are not shared
    # across patches here, so compare the sets of cell centres)
    got = np.sort(np.round(out.points[out.cells[0].data].mean(axis=1), 9), axis=0)
    want_c = np.sort(np.round(mesh.points[mesh.cells[0].data].mean(axis=1), 9), axis=0)
    np.testing.assert_allclose(got, want_c[:, : got.shape[1]], atol=1e-8)
    assert "u" in out.point_data


def test_sample_subdivisions_split_every_span():
    mesh = _nurbs("square-disc-nurbs-patch")
    one = _iga.iga_sample(mesh, subdivisions=1)
    three = _iga.iga_sample(mesh, subdivisions=3)
    assert len(three.cells[0].data) == 9 * len(one.cells[0].data)


def test_brep_sample_is_the_cad_json_sample():
    mesh = _cad("plate_with_hole.cad.json")
    a = _iga.iga_sample(mesh, subdivisions=6)
    b = mp.cad_json.sample(mesh.info["brep"], subdivisions=6)
    np.testing.assert_array_equal(a.points, b.points)


# ---------------------------------------------------------------------------
# Round trips through the writers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", NURBS)
def test_refined_nurbs_writes_and_reads_back(name, tmp_path):
    mesh = _nurbs(name)
    refined = _iga.iga_refine(mesh, subdivisions=2, elevate=1)
    path = tmp_path / "refined.mesh"
    mp.mfem.write(path, refined)
    assert path.read_text().startswith("MFEM NURBS mesh v1.0\n")
    back = mp.mfem.read(path)
    got, want = back.info["nurbs"], refined.info["nurbs"]
    assert len(got["patches"]) == len(want["patches"])
    for g, w in zip(got["patches"], want["patches"]):
        np.testing.assert_array_equal(g["degrees"], w["degrees"])
        for kg, kw in zip(g["knots"], w["knots"]):
            np.testing.assert_allclose(kg, kw, atol=1e-14)
        np.testing.assert_allclose(g["control_points"], w["control_points"], atol=1e-12)
        np.testing.assert_allclose(g["weights"], w["weights"], atol=1e-12)


@pytest.mark.parametrize("name", CAD)
def test_refined_brep_writes_and_reads_back(name, tmp_path):
    mesh = _cad(name)
    refined = _iga.iga_refine(mesh, subdivisions=2, elevate=1)
    path = tmp_path / "refined.cad.json"
    mp.write(path, refined)
    back = mp.read(path)
    for g, w in zip(back.info["brep"]["patches"], refined.info["brep"]["patches"]):
        np.testing.assert_array_equal(g["degrees"], w["degrees"])
        for kg, kw in zip(g["knots"], w["knots"]):
            np.testing.assert_array_equal(kg, kw)
        np.testing.assert_array_equal(g["control_points"], w["control_points"])
    # and the loops survived byte for byte
    for fa, fb in zip(back.info["brep"]["faces"], mesh.info["brep"]["faces"]):
        _same_tree(fa["loops"], fb["loops"])


# ---------------------------------------------------------------------------
# The API
# ---------------------------------------------------------------------------


def test_a_patch_topology_refuses_per_direction_refinement():
    mesh = _nurbs("square-disc-nurbs-patch")
    with pytest.raises(ValueError, match="every direction alike"):
        _iga.iga_refine(mesh, subdivisions=(2, 3))
    with pytest.raises(ValueError, match="every direction alike"):
        _iga.iga_refine(mesh, elevate=(0, 1))
    # a hand-built model has no neighbours to agree with
    patch = {
        "degrees": [1, 1],
        "knots": [[0, 0, 1, 1.0], [0, 0, 1, 1.0]],
        "control_points": [[0, 0], [1, 0], [0, 1], [1, 1.0]],
    }
    out = _iga.iga_refine({"patches": [patch]}, subdivisions=(2, 3), elevate=(0, 1))
    assert out.info["nurbs"]["patches"][0]["degrees"].tolist() == [1, 2]


def test_other_info_keys_are_kept_and_the_model_is_replaced():
    mesh = _nurbs("cube-nurbs")
    mesh.info["note"] = "kept"
    out = _iga.iga_refine(mesh, subdivisions=2)
    assert out.info["note"] == "kept"
    assert len(out.info["nurbs"]["patches"][0]["knots"][0]) > len(
        mesh.info["nurbs"]["patches"][0]["knots"][0]
    )


def test_a_bare_model_dict_is_accepted_and_filled_in():
    s = np.sqrt(0.5)
    arc = {
        "degrees": [2],
        "knots": [[0, 0, 0, 1, 1, 1]],
        "control_points": [[1, 0], [1, 1], [0, 1]],
        "weights": [1, s, 1],
    }
    out = _iga.iga_sample({"patches": [arc]}, subdivisions=4)
    assert out.cells[0].type == "line3" and len(out.cells[0].data) == 4
    assert out.points.shape[1] == 2


def test_public_names_are_exported():
    for name in ("iga_evaluate", "iga_refine", "iga_sample"):
        assert name in mp.__all__ and callable(getattr(mp, name))


@pytest.mark.parametrize(
    "call, match",
    [
        (lambda m: _iga.iga_refine(m, subdivisions=0), "subdivisions"),
        (lambda m: _iga.iga_refine(m, subdivisions=65), "subdivisions"),
        (lambda m: _iga.iga_refine(m, elevate=-1), "elevate"),
        (lambda m: _iga.iga_refine(m, elevate=40), "degree"),
        (lambda m: _iga.iga_refine(m, subdivisions=(1, 2, 3, 4)), "1 to 3"),
        (lambda m: _iga.iga_sample(m, subdivisions=-1), "subdivisions"),
        (lambda m: _iga.iga_evaluate(m, [[5.0, 5.0]]), "outside the knot range"),
        (lambda m: _iga.iga_evaluate(m, [[0.5]]), "parameters"),
        (lambda m: _iga.iga_evaluate(m, [[0.5, 0.5]], patch=9), "outside"),
        (lambda m: _iga.iga_evaluate(m, [[0.5, 0.5]], derivatives=2), "derivatives"),
        (lambda m: _iga.iga_refine(m, model="brep"), "no 'brep'"),
        (lambda m: _iga.iga_refine(m, model="stl"), "model must be"),
    ],
)
def test_bad_arguments_are_value_errors(call, match):
    mesh = _nurbs("square-disc-nurbs-patch")
    with pytest.raises(ValueError, match=match):
        call(mesh)


def test_a_mesh_without_a_model_is_refused():
    with pytest.raises(ValueError, match="no NURBS model"):
        _iga.iga_refine(
            mp.Mesh(np.zeros((3, 3)), [("triangle", np.array([[0, 1, 2]]))])
        )


def test_both_models_must_be_named():
    mesh = _cad("cylindrical_roof.cad.json")
    mesh.info["nurbs"] = _nurbs("cube-nurbs").info["nurbs"]
    with pytest.raises(ValueError, match="both"):
        _iga.iga_refine(mesh)
    assert _iga.iga_refine(mesh, model="brep").info["brep"]
    assert _iga.iga_refine(mesh, model="nurbs").info["nurbs"]


def test_faces_sharing_a_patch_must_agree_on_rational():
    mesh = _cad("two_patch_beam.cad.json")
    brep = mesh.info["brep"]
    brep["faces"][1]["patch"] = brep["faces"][0]["patch"]
    brep["faces"][1]["rational"] = True
    brep["faces"][0]["rational"] = False
    with pytest.raises(ValueError, match="disagree"):
        _iga.iga_refine(mesh)


def test_a_zero_weight_is_refused():
    arc = {
        "degrees": [2],
        "knots": [[0, 0, 0, 1, 1, 1]],
        "control_points": [[1, 0], [1, 1], [0, 1]],
        "weights": [1, 0.0, 1],
    }
    with pytest.raises(ValueError, match="weight"):
        _iga.iga_refine({"patches": [arc]}, subdivisions=2)


@needs_core
def test_strict_core_runs_every_operation_natively(monkeypatch):
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    from meshioplusplus import _fallback

    monkeypatch.setattr(_fallback, "_strict", None)
    try:
        mesh = _nurbs("pipe-nurbs")
        _iga.iga_refine(mesh, subdivisions=2, elevate=1)
        _iga.iga_sample(mesh)
        _iga.iga_evaluate(mesh, [[0.5, 0.5, 0.5]], derivatives=1)
        cad_mesh = _cad("plate_with_hole.cad.json")
        _iga.iga_refine(cad_mesh, subdivisions=2)
        _iga.iga_sample(cad_mesh)
        _iga.iga_evaluate(cad_mesh, [[0.5, 0.5]], derivatives=1)
    finally:
        monkeypatch.setattr(_fallback, "_strict", None)


def test_mixed_degree_faces_do_not_share_a_block():
    # faces of degree 3 and 4 are both VTK_LAGRANGE_QUADRILATERAL: the sampler
    # keeps them in blocks of their own node counts (a regression: they used to
    # be written into one block, whatever their node counts)
    mesh = _cad("two_patch_beam.cad.json")
    out = _iga.iga_refine(mesh, elevate=1)
    nodes = sorted(
        c.data.shape[1] for c in out.cells if c.type == "VTK_LAGRANGE_QUADRILATERAL"
    )
    assert nodes == [16, 25]
