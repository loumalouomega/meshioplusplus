"""Mesh <-> regular-grid transfer.

The fixture is deliberately **anisotropic** (extent 1 x 2 x 4). On a unit cube
every length scale is 1 and every axis is interchangeable, so a per-axis bug --
a transposed reshape, a spacing taken from the wrong axis -- is invisible; the
tests below were written against a cube first and several of them passed with the
z and x axes swapped.
"""

from __future__ import annotations

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import (
    GridArray,
    GridSpec,
    expand_grid,
    interpolate_grid,
    power_spectrum,
    resample_grid,
    sample_grid,
    scatter_grid,
    squeeze_grid,
)

EXTENT = (1.0, 2.0, 4.0)


def _linear(points):
    """A field that is linear in every coordinate, with distinct coefficients."""
    return 1.0 + 2.0 * points[:, 0] + 3.0 * points[:, 1] - 5.0 * points[:, 2]


def _source(n=3):
    """A tetrahedral mesh over the anisotropic box, carrying a scalar and a vector."""
    mesh = mio.grid((n, n, n), spacing=[e / n for e in EXTENT])
    mesh = mio.convert_cells(mesh, mode="simplexify")
    p = mesh.points
    mesh.point_data["u"] = _linear(p)
    mesh.point_data["v"] = np.column_stack([p[:, 0], 2.0 * p[:, 1], -p[:, 2]])
    return mesh


def _spec(dims=(4, 4, 4)):
    return GridSpec(
        origin=(0.0, 0.0, 0.0),
        spacing=[EXTENT[k] / dims[k] for k in range(3)],
        dims=dims,
    )


# --------------------------------------------------------------------------- #
# the lattice value                                                           #
# --------------------------------------------------------------------------- #
def test_points_are_bit_identical_to_grid():
    """`GridSpec.points()` must reproduce `grid()`'s own coordinates exactly.

    Not merely to a tolerance: both evaluate `origin + index * spacing`
    independently per point, and anything that accumulates instead would make the
    last plane depend on the cell count in the low bits.
    """
    spec = _spec((5, 3, 7))
    mesh = mio.grid(spec.dims, spec.origin, spec.spacing)
    assert np.array_equal(spec.points(), mesh.points)


def test_shape_is_tensor_order_and_dims_is_world_order():
    spec = _spec((2, 3, 4))
    assert tuple(spec.dims) == (2, 3, 4)  # nx, ny, nz cell counts
    assert spec.shape == (5, 4, 3)  # (D, H, W) = (nz+1, ny+1, nx+1) points
    assert spec.num_points == 5 * 4 * 3
    assert spec.num_cells == 2 * 3 * 4


def test_values_index_as_z_y_x():
    """`values[c, k, j, i]` is channel c at lattice point (i, j, k)."""
    spec = _spec()
    array = sample_grid(_source(), spec)
    field = array.channel("u")
    grid_points = spec.points().reshape(spec.shape + (3,))
    assert np.allclose(field, _linear(grid_points.reshape(-1, 3)).reshape(spec.shape))


def test_upscale_preserves_the_box_exactly():
    spec = _spec((4, 4, 4))
    fine = spec.upscale(2)
    assert tuple(fine.dims) == (8, 8, 8)
    assert spec.same_bounds(fine)
    lo_a, hi_a = spec.bounds
    lo_b, hi_b = fine.bounds
    assert np.array_equal(lo_a, lo_b) and np.allclose(hi_a, hi_b, rtol=0, atol=1e-15)
    assert spec.scaling_factor(fine) == (2, 2, 2)


def test_scaling_factor_is_none_when_not_a_whole_number():
    spec = _spec((4, 4, 4))
    other = GridSpec(origin=spec.origin, spacing=spec.spacing / 1.5, dims=(6, 6, 5))
    assert spec.scaling_factor(other) is None


def test_spec_round_trips_through_both_serializations():
    spec = _spec((2, 3, 4))
    assert GridSpec.from_dict(spec.to_dict()) == spec
    settings = spec.to_settings()
    assert set(settings) == {"Origin", "Spacing", "Dims"}  # PascalCase, hand-edited
    assert set(spec.to_dict()) >= {"origin", "spacing", "dims"}  # snake_case, machine
    assert spec.to_dict()["layout"] == "channels_first_zyx"


def test_spec_rejects_a_degenerate_axis():
    with pytest.raises(ValueError, match="dims must be positive"):
        GridSpec(origin=(0, 0, 0), spacing=(1, 1, 1), dims=(4, 0, 4))
    with pytest.raises(ValueError, match="spacing must be positive"):
        GridSpec(origin=(0, 0, 0), spacing=(1, 0, 1), dims=(4, 4, 4))


def test_from_lattice_mesh_recovers_the_spec_and_refuses_a_non_lattice():
    spec = _spec((3, 4, 5))
    assert GridSpec.from_lattice_mesh(spec.mesh()) == spec
    with pytest.raises(ValueError, match="not a dense regular lattice"):
        GridSpec.from_lattice_mesh(_source())


# --------------------------------------------------------------------------- #
# sampling                                                                    #
# --------------------------------------------------------------------------- #
def test_a_linear_field_is_exact_at_every_lattice_point():
    spec = _spec()
    array = sample_grid(_source(), spec)
    got = array.channel("u").reshape(-1)
    assert np.max(np.abs(got - _linear(spec.points()))) < 1e-12


def test_channel_order_expands_components_the_way_feature_matrix_does():
    array = sample_grid(_source(), _spec())
    assert array.channels == ("u", "v_0", "v_1", "v_2")
    assert array.values.shape == (4,) + _spec().shape
    assert array.schema["layout"] == "channels_first_zyx"
    assert array.schema["channels"] == list(array.channels)


def test_sampling_is_deterministic():
    spec = _spec()
    mesh = _source()
    a = sample_grid(mesh, spec).values
    b = sample_grid(mesh, spec).values
    assert np.array_equal(a, b)


def test_cell_data_is_refused_by_name():
    mesh = _source()
    mesh.cell_data["tag"] = [np.zeros(len(b)) for b in mesh.cells]
    with pytest.raises(ValueError, match="data to-point"):
        sample_grid(mesh, _spec(), fields=["tag"])


def test_an_unknown_field_names_what_is_available():
    with pytest.raises(ValueError, match="no point data named 'nope'"):
        sample_grid(_source(), _spec(), fields=["nope"])


def test_points_outside_the_source_take_the_fill_and_are_counted():
    """A grid larger than the mesh: the excess is fill, and coverage says so."""
    spec = GridSpec(origin=(-1.0, -1.0, -1.0), spacing=(0.5, 1.0, 2.0), dims=(4, 4, 4))
    array = sample_grid(_source(), spec, fill_value=np.nan)
    values = array.channel("u")
    assert np.isnan(values).any()
    assert 0.0 < array.coverage < 1.0
    # coverage counts exactly the points that are not fill
    assert array.coverage == pytest.approx(
        float(np.count_nonzero(~np.isnan(values))) / values.size
    )


def test_extrapolate_fills_from_the_nearest_source_value():
    spec = GridSpec(origin=(-1.0, -1.0, -1.0), spacing=(0.5, 1.0, 2.0), dims=(4, 4, 4))
    filled = sample_grid(_source(), spec, fill_value=np.nan)
    extrapolated = sample_grid(_source(), spec, extrapolate=True)
    assert np.isnan(filled.channel("u")).any()
    assert not np.isnan(extrapolated.channel("u")).any()
    # coverage still means "inside", and costs an extra pass to say so
    assert extrapolated.coverage == pytest.approx(filled.coverage)


def test_float32_is_an_explicit_opt_in():
    assert sample_grid(_source(), _spec()).values.dtype == np.float64
    assert sample_grid(_source(), _spec(), float32=True).values.dtype == np.float32


# --------------------------------------------------------------------------- #
# trilinear evaluation                                                        #
# --------------------------------------------------------------------------- #
def test_trilinear_is_exact_for_a_linear_field_at_random_points():
    spec = _spec()
    array = sample_grid(_source(), spec)
    rng = np.random.default_rng(20260904)
    lo, hi = spec.bounds
    query = lo + rng.random((500, 3)) * (hi - lo)
    got = interpolate_grid(array.channel("u"), spec, query)
    assert np.max(np.abs(got - _linear(query))) < 1e-12


def test_trilinear_is_not_the_barycentric_interpolant():
    """The measured difference, which is why this is its own implementation.

    On a unit cell carrying ``u = x*y``, trilinear reproduces the field exactly
    (``x*y`` *is* trilinear) while barycentric interpolation over the cell's
    simplex decomposition gives twice the value at the centre -- the chosen
    diagonal showing through. Both are exact for a *linear* field, so a test that
    only used one would see no difference at all.
    """
    cell = mio.grid((1, 1, 1))
    cell.point_data["u"] = cell.points[:, 0] * cell.points[:, 1]
    spec = GridSpec(origin=(0, 0, 0), spacing=(1, 1, 1), dims=(1, 1, 1))
    array = sample_grid(cell, spec)

    query = np.array([[0.5, 0.5, 0.5], [0.25, 0.75, 0.5], [0.75, 0.25, 0.5]])
    exact = query[:, 0] * query[:, 1]
    trilinear = interpolate_grid(array.channel("u"), spec, query)
    barycentric = mio.interpolate(
        cell, mio.Mesh(query, []), method="barycentric"
    ).point_data["u"]

    assert np.allclose(trilinear, exact, atol=1e-15)
    assert np.allclose(barycentric, [0.5, 0.25, 0.25], atol=1e-15)
    assert barycentric[0] == pytest.approx(2.0 * trilinear[0])


def test_queries_outside_the_box_are_clamped_not_extrapolated():
    spec = _spec()
    array = sample_grid(_source(), spec)
    lo, hi = spec.bounds
    far = np.array([hi + 1000.0, lo - 1000.0])
    got = interpolate_grid(array.channel("u"), spec, far)
    corners = interpolate_grid(array.channel("u"), spec, np.array([hi, lo]))
    assert np.allclose(got, corners)


# --------------------------------------------------------------------------- #
# resampling and the full chain                                               #
# --------------------------------------------------------------------------- #
def test_resampling_onto_the_same_spec_is_exactly_the_identity():
    spec = _spec()
    values = sample_grid(_source(), spec).values
    assert np.array_equal(resample_grid(values, spec, spec), values)


def test_the_whole_chain_reproduces_a_linear_field():
    """Mesh -> coarse grid -> upsample -> fine mesh, with the network replaced by
    the trilinear baseline. This is the production path; a linear field must
    survive all four steps.
    """
    fine_mesh = mio.grid((6, 6, 6), spacing=[e / 6 for e in EXTENT])
    # the shared-box rule: the box comes from the FINE mesh, so every fine node
    # is inside the coarse grid and the scatter back is interpolation everywhere
    coarse = GridSpec.from_mesh(fine_mesh, resolution=(4, 4, 4))
    array = sample_grid(_source(), coarse)

    fine = coarse.upscale(2)
    upsampled = resample_grid(array.values, coarse, fine)
    out = scatter_grid(
        GridArray(upsampled, fine, array.channels, dict(array.schema)), fine_mesh
    )

    assert np.max(np.abs(out.point_data["u"] - _linear(fine_mesh.points))) < 1e-12
    expected = np.column_stack(
        [fine_mesh.points[:, 0], 2.0 * fine_mesh.points[:, 1], -fine_mesh.points[:, 2]]
    )
    assert out.point_data["v"].shape == expected.shape
    assert np.max(np.abs(out.point_data["v"] - expected)) < 1e-12


# --------------------------------------------------------------------------- #
# scattering                                                                  #
# --------------------------------------------------------------------------- #
def test_scatter_rebuilds_multi_component_arrays_from_the_channel_contract():
    spec = _spec()
    array = sample_grid(_source(), spec)
    target = mio.grid((3, 3, 3), spacing=[e / 3 for e in EXTENT])
    out = scatter_grid(array, target)
    assert sorted(out.point_data) == ["u", "v"]
    assert out.point_data["u"].shape == (len(target.points),)
    assert out.point_data["v"].shape == (len(target.points), 3)


def test_scatter_does_not_mutate_the_target():
    array = sample_grid(_source(), _spec())
    target = mio.grid((2, 2, 2), spacing=[e / 2 for e in EXTENT])
    before = sorted(target.point_data)
    scatter_grid(array, target)
    assert sorted(target.point_data) == before


def test_scatter_conflict_policy():
    array = sample_grid(_source(), _spec())
    target = mio.grid((2, 2, 2), spacing=[e / 2 for e in EXTENT])
    target.point_data["u"] = np.zeros(len(target.points))

    with pytest.raises(ValueError, match="already has point data named 'u'"):
        scatter_grid(array, target)
    assert not np.array_equal(
        scatter_grid(array, target, on_conflict="overwrite").point_data["u"],
        target.point_data["u"],
    )
    suffixed = scatter_grid(array, target, on_conflict="suffix")
    assert "u_grid" in suffixed.point_data
    assert np.array_equal(suffixed.point_data["u"], target.point_data["u"])


# --------------------------------------------------------------------------- #
# GridArray                                                                   #
# --------------------------------------------------------------------------- #
def test_grid_array_round_trips_through_a_mesh_bit_identically():
    array = sample_grid(_source(), _spec())
    back = GridArray.from_mesh(array.to_mesh())
    assert back.channels == array.channels
    assert back.spec == array.spec
    assert np.array_equal(back.values, array.values)


def test_grid_array_round_trips_through_a_vti_file(tmp_path):
    """`.vti` is the format that keeps a grid a grid: it stores the lattice as
    origin/spacing/extent and regenerates the points arithmetically, so the spec
    is recovered exactly with no tolerance.
    """
    array = sample_grid(_source(), _spec())
    path = tmp_path / "grid.vti"
    mio.write(str(path), array.to_mesh())
    back = GridArray.from_mesh(mio.read(str(path)))
    assert back.spec == array.spec
    assert back.channels == array.channels
    assert np.allclose(back.values, array.values, rtol=0, atol=1e-15)


def test_grid_array_rejects_a_shape_that_disagrees_with_its_spec():
    spec = _spec()
    with pytest.raises(ValueError, match="but the spec describes"):
        GridArray(np.zeros((1, 2, 2, 2)), spec, ("u",))
    with pytest.raises(ValueError, match="channels in values but"):
        GridArray(np.zeros((1,) + spec.shape), spec, ("u", "v"))


def test_channel_lookup_names_what_is_available():
    array = sample_grid(_source(), _spec())
    with pytest.raises(ValueError, match="no channel named 'w'"):
        array.channel("w")


# --------------------------------------------------------------------------- #
# the thin-axis idiom                                                         #
# --------------------------------------------------------------------------- #
def test_squeeze_takes_a_world_axis_not_a_tensor_axis():
    values = np.arange(2 * 3 * 4 * 5, dtype=float).reshape(2, 3, 4, 5)  # C, z, y, x
    assert np.array_equal(squeeze_grid(values, 2, index=1), values[:, 1])  # z
    assert np.array_equal(squeeze_grid(values, 1, index=2), values[:, :, 2])  # y
    assert np.array_equal(squeeze_grid(values, 0, index=4), values[..., 4])  # x


def test_squeeze_requires_an_index_when_the_axis_is_not_a_single_plane():
    values = np.zeros((2, 3, 4, 5))
    with pytest.raises(ValueError, match="pass index="):
        squeeze_grid(values, 2)
    with pytest.raises(ValueError, match="out of range"):
        squeeze_grid(values, 2, index=9)


def test_expand_inverts_squeeze():
    values = np.arange(2 * 1 * 4 * 5, dtype=float).reshape(2, 1, 4, 5)
    squeezed = squeeze_grid(values, 2)
    assert squeezed.shape == (2, 4, 5)
    assert np.array_equal(expand_grid(squeezed, 2), values)


# --------------------------------------------------------------------------- #
# power spectrum                                                              #
# --------------------------------------------------------------------------- #
def _cubic(n=16):
    return GridSpec(origin=(0, 0, 0), spacing=(1.0 / n,) * 3, dims=(n - 1,) * 3)


def _index_grids(spec):
    d, h, w = spec.shape
    return np.meshgrid(np.arange(d), np.arange(h), np.arange(w), indexing="ij")


def test_power_sums_to_the_mean_square_exactly():
    """Parseval. It holds only because every mode is binned, including the ones
    past the isotropic Nyquist in the corners of the box -- truncating the tail
    would break it, which is why `counts` exists instead.
    """
    spec = _cubic()
    field = np.random.default_rng(11).standard_normal(spec.shape)
    ps = power_spectrum(field, spec)
    assert ps.power.sum() == pytest.approx(float((field * field).mean()), rel=1e-12)
    assert int(ps.counts.sum()) == spec.num_points


def test_a_constant_field_puts_all_its_power_in_bin_zero():
    spec = _cubic()
    ps = power_spectrum(np.full(spec.shape, 3.0), spec)
    assert ps.power[0] == pytest.approx(9.0)
    assert ps.power[1:].sum() == pytest.approx(0.0, abs=1e-20)


def test_a_single_mode_lands_in_a_single_bin():
    spec = _cubic()
    _, _, x = _index_grids(spec)
    ps = power_spectrum(np.cos(2 * np.pi * 3 * x / spec.shape[2]), spec)
    loud = np.nonzero(ps.power > 1e-12)[0]
    assert list(loud) == [3]
    assert ps.wavenumber[3] == pytest.approx(3.0)


def test_the_spectrum_is_rotation_invariant():
    """The property that discriminates an azimuthal average from per-axis binning.

    A 90-degree rotation of a cubic grid is an array transpose, and a genuine
    shell average cannot see it. Binning on one component of k instead --  a
    plausible wrong implementation -- disagrees here by about 7e-4, twelve orders
    of magnitude above what this asserts.
    """
    spec = _cubic()
    z, y, x = _index_grids(spec)
    field = np.exp(-((x - 7.3) ** 2 + (y - 4.1) ** 2 + (z - 9.7) ** 2) / 6.0) * np.sin(
        0.7 * x + 0.3 * y
    )
    a = power_spectrum(field, spec)
    b = power_spectrum(np.transpose(field, (2, 0, 1)), spec)
    assert np.max(np.abs(a.power - b.power)) < 1e-15


def test_a_vector_field_sums_the_power_of_its_components():
    spec = _cubic(8)
    rng = np.random.default_rng(5)
    field = rng.standard_normal((3,) + spec.shape)
    total = power_spectrum(field, spec)
    parts = sum(power_spectrum(field[c], spec).power for c in range(3))
    assert np.allclose(total.power, parts)


def test_an_anisotropic_lattice_is_refused_by_name():
    spec = GridSpec(origin=(0, 0, 0), spacing=(1.0, 1.0, 2.0), dims=(7, 7, 7))
    with pytest.raises(ValueError, match="needs an isotropic"):
        power_spectrum(np.zeros(spec.shape), spec)


def test_the_spectrum_separates_fields_a_pointwise_metric_cannot():
    """What the spectrum is *for*.

    Two single-mode fields of the same amplitude have exactly the same mean
    square, so RMS -- and any other pointwise reduction -- reports them as
    equally good. Their spectra are unambiguous about which one carries the fine
    structure, which is the whole reason a superresolution result is reported
    against a spectrum rather than an error norm.
    """
    spec = _cubic()
    _, _, x = _index_grids(spec)
    n = spec.shape[2]
    smooth = np.sin(2 * np.pi * 1 * x / n)
    rough = np.sin(2 * np.pi * 6 * x / n)

    assert (smooth * smooth).mean() == pytest.approx((rough * rough).mean())

    a = power_spectrum(smooth, spec)
    b = power_spectrum(rough, spec)
    assert int(np.argmax(a.power)) == 1
    assert int(np.argmax(b.power)) == 6
    assert a.power[6] < 1e-20 and b.power[1] < 1e-20
