"""Property-based tests (Hypothesis) over invariants the docs state.

Each property names the page that promises it. The meshes are small jittered
lattices -- a ``grid()`` hexahedral block, simplexified to tetrahedra, with its
*interior* points moved so the boundary stays the exact box (which is what
makes whole-mesh volume and mass comparisons exact rather than approximate).

Profiles: ``ci`` (the default; derandomized, 25 examples, so a run is
reproducible) and ``dev`` (``HYPOTHESIS_PROFILE=dev``, 300 random examples) for
hunting locally. A counterexample Hypothesis finds is shrunk and printed; pin
it as an ordinary regression test next to the fix.
"""

import os

import numpy as np
import pytest

import meshioplusplus

hypothesis = pytest.importorskip("hypothesis")

from hypothesis import HealthCheck, given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

settings.register_profile(
    "ci",
    max_examples=25,
    derandomize=True,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
    print_blob=True,
)
settings.register_profile("dev", max_examples=300, deadline=None)
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "ci"))


# ----------------------------------------------------------------------------
# Strategies
# ----------------------------------------------------------------------------


@st.composite
def boxes(draw):
    """(dims, origin, spacing) of a small lattice."""
    dims = tuple(draw(st.integers(1, 3)) for _ in range(3))
    origin = tuple(draw(st.floats(-5, 5, allow_nan=False)) for _ in range(3))
    spacing = tuple(draw(st.floats(0.5, 2.0, allow_nan=False)) for _ in range(3))
    return dims, origin, spacing


def _interior(points, lo, hi):
    tol = 1e-9 * max(1.0, float(np.max(np.abs(points))))
    return np.all((points > lo + tol) & (points < hi - tol), axis=1)


@st.composite
def tet_meshes(draw, jitter=True):
    """A simplexified lattice whose interior points are jittered."""
    dims, origin, spacing = draw(boxes())
    hexes = meshioplusplus.grid(dims, origin=origin, spacing=spacing)
    mesh = meshioplusplus.convert_cells(hexes, "simplexify")
    points = np.asarray(mesh.points, dtype=float).copy()
    if jitter:
        lo, hi = points.min(axis=0), points.max(axis=0)
        inner = _interior(points, lo, hi)
        amp = 0.15 * min(spacing)
        seed = draw(st.integers(0, 2**32 - 1))
        rng = np.random.default_rng(seed)
        points[inner] += rng.uniform(-amp, amp, size=(int(inner.sum()), 3))
    return meshioplusplus.Mesh(points, [(b.type, b.data) for b in mesh.cells])


@st.composite
def rigid_motions(draw):
    """A proper rotation and a translation."""
    q, r = np.linalg.qr(
        np.array(
            [
                [draw(st.floats(-1, 1, allow_nan=False)) for _ in range(3)]
                for _ in range(3)
            ]
        )
        + 3 * np.eye(3)
    )
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    t = np.array([draw(st.floats(-10, 10, allow_nan=False)) for _ in range(3)])
    return q, t


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def _tet_volumes(mesh):
    p = np.asarray(mesh.points, dtype=float)
    out = []
    for block in mesh.cells:
        assert block.type == "tetra", block.type
        c = p[np.asarray(block.data)]
        out.append(
            np.einsum(
                "ij,ij->i",
                np.cross(c[:, 1] - c[:, 0], c[:, 2] - c[:, 0]),
                c[:, 3] - c[:, 0],
            )
            / 6.0
        )
    return np.concatenate(out)


def _volume(mesh):
    return meshioplusplus.compute_stats(mesh)["unsigned_volume"]


def _area(mesh):
    return meshioplusplus.compute_stats(meshioplusplus.extract_surface(mesh))[
        "total_area"
    ]


def _strip(mesh):
    return meshioplusplus.Mesh(mesh.points, [(b.type, b.data) for b in mesh.cells])


# ----------------------------------------------------------------------------
# Volume conservation -- doc/refine.md "Guarantees", doc/agglomerate.md
# ----------------------------------------------------------------------------


@given(tet_meshes())
def test_refine_conserves_tetrahedral_volume(mesh):
    before = _volume(mesh)
    after = _volume(meshioplusplus.refine(mesh))
    assert after == pytest.approx(before, rel=1e-12)


@given(boxes(), rigid_motions())
def test_refine_conserves_volume_of_affine_hexahedra(box, motion):
    dims, origin, spacing = box
    q, t = motion
    hexes = meshioplusplus.grid(dims, origin=origin, spacing=spacing)
    hexes = meshioplusplus.Mesh(
        np.asarray(hexes.points) @ q.T + t, [(b.type, b.data) for b in hexes.cells]
    )
    assert _volume(meshioplusplus.refine(hexes)) == pytest.approx(
        _volume(hexes), rel=1e-12
    )


@given(tet_meshes(), st.integers(2, 12))
def test_agglomerate_keeps_signed_volume(mesh, group):
    before = meshioplusplus.compute_stats(mesh)["signed_volume"]
    after = meshioplusplus.compute_stats(meshioplusplus.agglomerate(mesh, group))[
        "signed_volume"
    ]
    assert after == pytest.approx(before, rel=1e-12)


# ----------------------------------------------------------------------------
# Conformity -- doc/refine.md "Closures": redgreen and propagate never leave a
# hanging node, so the boundary of the refined mesh is the old boundary.
# ----------------------------------------------------------------------------


@given(tet_meshes(), st.sampled_from(["redgreen", "propagate"]), st.data())
def test_adaptive_refinement_is_conforming(mesh, closure, data):
    n = sum(len(b.data) for b in mesh.cells)
    chosen = data.draw(
        st.lists(st.integers(0, n - 1), min_size=1, max_size=n, unique=True)
    )
    refined = meshioplusplus.refine(mesh, cells=sorted(chosen), closure=closure)
    # A hanging node splits a face on one side only: the unsplit face and its
    # split halves would then all be counted once, inflating the "surface".
    assert _area(refined) == pytest.approx(_area(mesh), rel=1e-11)
    assert _volume(refined) == pytest.approx(_volume(mesh), rel=1e-12)


# ----------------------------------------------------------------------------
# Mass conservation -- doc/conservative_interpolate.md
# ----------------------------------------------------------------------------


@given(tet_meshes(), tet_meshes(jitter=False), st.integers(0, 2**32 - 1))
def test_conservative_interpolation_conserves_mass(source, target, seed):
    # Same box for both meshes, so "the region the two share" is all of it.
    lo, hi = np.min(source.points, axis=0), np.max(source.points, axis=0)
    tp = np.asarray(target.points, dtype=float)
    tlo, thi = tp.min(axis=0), tp.max(axis=0)
    target = meshioplusplus.Mesh(
        lo + (tp - tlo) / (thi - tlo) * (hi - lo),
        [(b.type, b.data) for b in target.cells],
    )
    rng = np.random.default_rng(seed)
    n = sum(len(b.data) for b in source.cells)
    source.cell_data["m"] = [rng.uniform(0.1, 10.0, n)]
    out = meshioplusplus.conservative_interpolate(source, target, arrays=["m"])
    mass_in = float(np.dot(source.cell_data["m"][0], _tet_volumes(source)))
    mass_out = float(np.dot(np.concatenate(out.cell_data["m"]), _tet_volumes(out)))
    assert mass_out == pytest.approx(mass_in, rel=1e-9)


# ----------------------------------------------------------------------------
# Partition of unity -- doc/interpolate.md: barycentric weights sum to one and
# reproduce linear fields exactly.
# ----------------------------------------------------------------------------


@given(tet_meshes(), st.integers(0, 2**32 - 1), st.floats(-3, 3), st.floats(-3, 3))
def test_barycentric_interpolation_reproduces_linear_fields(mesh, seed, a, c):
    p = np.asarray(mesh.points, dtype=float)
    lo, hi = p.min(axis=0), p.max(axis=0)
    rng = np.random.default_rng(seed)
    samples = lo + rng.uniform(0.05, 0.95, size=(20, 3)) * (hi - lo)
    mesh.point_data["one"] = np.ones(len(p))
    mesh.point_data["lin"] = a * p[:, 0] - 2.0 * p[:, 1] + c * p[:, 2] + 1.0
    target = meshioplusplus.Mesh(samples, [("vertex", np.arange(20)[:, None])])
    out = meshioplusplus.interpolate(
        mesh, target, method="barycentric", arrays=["one", "lin"]
    )
    np.testing.assert_allclose(out.point_data["one"], 1.0, rtol=0, atol=1e-12)
    want = a * samples[:, 0] - 2.0 * samples[:, 1] + c * samples[:, 2] + 1.0
    np.testing.assert_allclose(out.point_data["lin"], want, rtol=1e-9, atol=1e-9)


# ----------------------------------------------------------------------------
# Determinism -- doc/refine.md, doc/interpolate.md "byte-identical"
# ----------------------------------------------------------------------------


def _same(a, b):
    assert np.array_equal(np.asarray(a.points), np.asarray(b.points))
    assert [x.type for x in a.cells] == [x.type for x in b.cells]
    for x, y in zip(a.cells, b.cells):
        assert np.array_equal(np.asarray(x.data), np.asarray(y.data))


@given(tet_meshes(), st.sampled_from(["redgreen", "propagate", "balanced"]))
def test_refine_is_deterministic(mesh, closure):
    cells = list(range(0, sum(len(b.data) for b in mesh.cells), 3))
    _same(
        meshioplusplus.refine(mesh, cells=cells, closure=closure),
        meshioplusplus.refine(mesh, cells=cells, closure=closure),
    )


@given(tet_meshes(), st.integers(0, 2**32 - 1))
def test_interpolation_is_byte_identical_across_the_core_boundary(mesh, seed):
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._interpolate import _interpolate_py

    p = np.asarray(mesh.points, dtype=float)
    rng = np.random.default_rng(seed)
    mesh.point_data["f"] = rng.normal(size=len(p))
    lo, hi = p.min(axis=0), p.max(axis=0)
    samples = lo + rng.uniform(-0.1, 1.1, size=(15, 3)) * (hi - lo)
    target = meshioplusplus.Mesh(samples, [("vertex", np.arange(15)[:, None])])
    got = core.interpolate(mesh, target, "barycentric", ["f"], True, -1.0, "error")
    want = _interpolate_py(
        mesh,
        target,
        method="barycentric",
        arrays=["f"],
        extrapolate=True,
        default_value=-1.0,
        on_conflict="error",
    )
    assert (
        np.asarray(got.point_data["f"]).tobytes()
        == np.asarray(want.point_data["f"]).tobytes()
    )


# ----------------------------------------------------------------------------
# Map composition -- doc/refine.md: `refine:parent_cell` names the *original*
# ancestor at any depth, so two single passes compose to one double pass.
# ----------------------------------------------------------------------------


@given(tet_meshes())
def test_refine_parent_maps_compose(mesh):
    once = meshioplusplus.refine(mesh, record_parent_ids=True)
    twice = meshioplusplus.refine(_strip(once), record_parent_ids=True)
    direct = meshioplusplus.refine(mesh, levels=2, record_parent_ids=True)
    _same(twice, direct)
    first = np.ravel(once.cell_data["refine:parent_cell"][0])
    second = np.ravel(twice.cell_data["refine:parent_cell"][0])
    np.testing.assert_array_equal(
        first[second], np.ravel(direct.cell_data["refine:parent_cell"][0])
    )


@given(boxes())
def test_simplexify_parent_map_is_total_and_in_range(box):
    dims, origin, spacing = box
    hexes = meshioplusplus.grid(dims, origin=origin, spacing=spacing)
    out = meshioplusplus.convert_cells(hexes, "simplexify", record_parent_ids=True)
    (key,) = [k for k in out.cell_data if k.endswith("parent_cell")]
    parents = np.concatenate(out.cell_data[key])
    n = len(hexes.cells[0].data)
    assert parents.min() >= 0 and parents.max() < n
    # Every hexahedron is covered by the same number of tetrahedra.
    counts = np.bincount(parents, minlength=n)
    assert np.all(counts == counts[0]) and counts[0] > 0
    assert _volume(out) == pytest.approx(_volume(hexes), rel=1e-12)
