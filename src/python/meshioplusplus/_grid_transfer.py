"""Mesh to regular-grid transfer (``meshioplusplus.sample_grid`` and friends).

A convolutional model does not take a mesh. It takes a dense rectangular array,
and getting a field from an unstructured mesh onto one -- and the prediction back
off it -- is the whole of the data path a grid-shaped architecture needs. That is
what this module is: a lattice as a value, the two transfers, and the one metric
that says whether a super-resolved field has the right small-scale content.

Everything here is pure numpy over machinery that already exists. Sampling is
:func:`meshioplusplus.interpolate` onto the lattice's points, so it inherits that
operation's exactness for a linear field and its own outside-the-domain warning;
the column contract is :func:`meshioplusplus.feature_matrix`'s, so a grid's
channels and a graph's node features are named by one rule rather than two.

Array layout
------------
A sampled grid is ``(C, D, H, W)``: channels first, then **z, y, x**. That is
`torch.nn.Conv3d`'s layout and VTK ImageData's, and it is what
:func:`meshioplusplus.grid`'s x-fastest point numbering already produces, so a
flat per-point array reshapes into it with no transpose::

    values[c, k, j, i]      is channel c at lattice point (i, j, k)

:class:`GridSpec` carries both spellings and they are deliberately different
words: ``spec.dims`` is **world** order ``(nx, ny, nz)`` cell counts, matching
``grid``/``voxelize``/``compute_sdf``, and ``spec.shape`` is **tensor** order
``(nz+1, ny+1, nx+1)`` point counts, matching the array. Anything reading a
checkpoint months later needs to be told which it is looking at, so the model
card records the layout as a literal string.

The shared-box rule
-------------------
Take the box from the **fine** mesh. Then every fine node lies inside the coarse
grid, and scattering the prediction back is interpolation everywhere rather than
extrapolation at the boundary. Sampling happens at lattice *points*, not cell
centres, for the same reason: the points reach the box's faces and the centres
do not.

Why trilinear is its own implementation
---------------------------------------
:func:`meshioplusplus.interpolate`'s ``"barycentric"`` mode would also evaluate a
field on a lattice, and both are exact for a *linear* field, but they are not the
same interpolant and the difference is not small. On a unit cell carrying
``u = x*y``, trilinear reproduces the field **exactly** (0.25 at the cell centre,
because ``x*y`` is trilinear) while barycentric-on-simplexified hexahedra gives
**0.5** -- a factor of two, measured, not estimated. The simplex interpolant
exposes the chosen diagonal; trilinear does not, it is what
``torch.nn.functional.interpolate`` does, and it is what every superresolution
baseline is reported against. Adding it as a third ``interpolate`` method would
also break that function's byte-identical C++/numpy twin contract for an
algorithm that only makes sense on a dense lattice.

Public API:

* :class:`GridSpec` -- a lattice as a serializable value.
* :class:`GridArray` -- values on one, plus the channel contract.
* :func:`sample_grid` -- mesh -> tensor.
* :func:`scatter_grid` -- tensor -> mesh.
* :func:`interpolate_grid` -- trilinear evaluation at arbitrary points.
* :func:`resample_grid` -- grid -> grid, the upsampling baseline.
* :func:`squeeze_grid` / :func:`expand_grid` -- the thin-axis idiom.
* :func:`power_spectrum` -- azimuthally averaged, Parseval-exact.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

from ._grid import DEFAULT_MAX_CELLS, _lattice_py, lattice_from_mesh
from ._mesh import Mesh

__all__ = [
    "GRID_SCHEMA_VERSION",
    "GridSpec",
    "GridArray",
    "PowerSpectrum",
    "sample_grid",
    "scatter_grid",
    "interpolate_grid",
    "resample_grid",
    "squeeze_grid",
    "expand_grid",
    "power_spectrum",
]

#: Bumped when the recorded channel/layout contract changes meaning.
GRID_SCHEMA_VERSION = 1

#: The layout, recorded verbatim in every schema. A checkpoint that does not say
#: which axis is which produces finite, plausible, transposed numbers.
GRID_LAYOUT = "channels_first_zyx"

_PREFIX = "meshio++: sample_grid: "

#: The constant channel ridden through the sampler to recover the inside mask.
#: Barycentric weights sum to one, so it comes back as 1 inside the source domain
#: and as the fill value outside -- an exact mask for the price of one extra
#: column in a call that is already being made.
_PROBE = "__meshioplusplus_coverage_probe__"


# --------------------------------------------------------------------------- #
# the lattice as a value                                                      #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class GridSpec:
    """A regular lattice: where it starts, how big its cells are, how many.

    The same triple ``detail::LatticeSpec`` carries, in the same order, so this
    is a *value* rather than a description -- two specs are equal when they
    describe the same lattice, and one serializes to a settings document and
    back without loss.
    """

    #: lo corner, ``(3,)`` float64.
    origin: np.ndarray
    #: cell size per axis, ``(3,)`` float64.
    spacing: np.ndarray
    #: **cell** counts ``(nx, ny, nz)``, world order, ``(3,)`` int64.
    dims: np.ndarray

    def __post_init__(self):
        origin = np.ascontiguousarray(
            np.asarray(self.origin, dtype=np.float64).reshape(-1)
        )
        spacing = np.ascontiguousarray(
            np.asarray(self.spacing, dtype=np.float64).reshape(-1)
        )
        dims = np.ascontiguousarray(np.asarray(self.dims, dtype=np.int64).reshape(-1))
        if origin.size != 3 or spacing.size != 3 or dims.size != 3:
            raise ValueError(
                "meshio++: GridSpec: origin, spacing and dims must each have three "
                "components"
            )
        for k in range(3):
            if dims[k] <= 0:
                raise ValueError(
                    f"meshio++: GridSpec: dims must be positive on every axis, got "
                    f"{int(dims[k])} on axis {k}"
                )
            if not spacing[k] > 0.0:
                raise ValueError(
                    f"meshio++: GridSpec: spacing must be positive on every axis, got "
                    f"{float(spacing[k])} on axis {k}"
                )
        origin.flags.writeable = False
        spacing.flags.writeable = False
        dims.flags.writeable = False
        object.__setattr__(self, "origin", origin)
        object.__setattr__(self, "spacing", spacing)
        object.__setattr__(self, "dims", dims)

    # -- shape -------------------------------------------------------------- #
    @property
    def shape(self) -> Tuple[int, int, int]:
        """Point counts in **tensor** order ``(D, H, W)`` = ``(z, y, x)``."""
        return (
            int(self.dims[2]) + 1,
            int(self.dims[1]) + 1,
            int(self.dims[0]) + 1,
        )

    @property
    def num_points(self) -> int:
        d, h, w = self.shape
        return d * h * w

    @property
    def num_cells(self) -> int:
        return int(self.dims[0]) * int(self.dims[1]) * int(self.dims[2])

    @property
    def bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """``(lo, hi)`` corners."""
        return (
            self.origin.copy(),
            self.origin + self.dims.astype(np.float64) * self.spacing,
        )

    @property
    def is_isotropic(self) -> bool:
        """Whether every axis has the same cell size (exactly)."""
        return bool(self.spacing[0] == self.spacing[1] == self.spacing[2])

    # -- geometry ----------------------------------------------------------- #
    def points(self) -> np.ndarray:
        """The lattice points, ``(N, 3)`` float64, x fastest.

        Bit-identical to ``grid(dims, origin, spacing).points`` -- the same
        ``origin + index * spacing`` evaluated the same way, never accumulated.
        """
        pts, _ = _lattice_py(self.dims, self.origin, self.spacing)
        return pts

    def mesh(self, max_cells: int = DEFAULT_MAX_CELLS) -> Mesh:
        """The lattice as an ordinary ``hexahedron`` mesh."""
        from ._grid import grid as _grid

        return _grid(self.dims, self.origin, self.spacing, max_cells=max_cells)

    # -- relations ---------------------------------------------------------- #
    def upscale(self, factor) -> "GridSpec":
        """The same box at ``factor`` times the resolution.

        The bounds are preserved exactly (the new spacing is recomputed from the
        box rather than divided, so ``lo + dims*spacing`` still lands on ``hi``),
        which is what makes a coarse/fine pair share a domain.
        """
        f = np.asarray(factor, dtype=np.int64).reshape(-1)
        if f.size == 1:
            f = np.repeat(f, 3)
        if f.size != 3 or np.any(f <= 0):
            raise ValueError(
                "meshio++: GridSpec.upscale: factor must be one positive integer, or "
                "three"
            )
        lo, hi = self.bounds
        dims = self.dims * f
        spacing = (hi - lo) / dims.astype(np.float64)
        return GridSpec(origin=lo, spacing=spacing, dims=dims)

    def same_bounds(self, other: "GridSpec", rtol: float = 1e-9) -> bool:
        """Whether two specs cover the same box, to a relative tolerance."""
        lo_a, hi_a = self.bounds
        lo_b, hi_b = other.bounds
        scale = np.maximum(np.abs(hi_a - lo_a), 1.0)
        return bool(
            np.all(np.abs(lo_a - lo_b) <= rtol * scale)
            and np.all(np.abs(hi_a - hi_b) <= rtol * scale)
        )

    def upscale_samples(self, factor) -> "GridSpec":
        """The same box with ``factor`` times as many **sample points**.

        This is :meth:`upscale`'s sibling, and the difference between them is an
        off-by-one that a superresolution model turns into a shape error deep
        inside a loss function.

        A lattice's samples are its *corners*, so ``dims`` cells carry
        ``dims + 1`` points. :meth:`upscale` multiplies the **cells**, which is
        what resampling wants -- every coarse point is then also a fine point,
        and the two grids nest. A convolutional upsampler multiplies the
        **samples**: ``SRResNet(scaling_factor=s)`` maps ``(B, C, D, H, W)`` to
        ``(B, C, sD, sH, sW)``. Measured on a 4x4x4-cell grid at ``s = 2``:
        ``upscale`` gives a 9x9x9-point grid (all 125 coarse points nest),
        while the model needs 10x10x10.

        So use this one to build the *fine* spec of a coarse/fine pair, and
        :meth:`upscale` to resample. The price is that only the box's eight
        corners are shared between the two lattices -- invisible to a model,
        which never sees a coordinate, and the bounds are identical either way.
        """
        f = np.asarray(factor, dtype=np.int64).reshape(-1)
        if f.size == 1:
            f = np.repeat(f, 3)
        if f.size != 3 or np.any(f <= 0):
            raise ValueError(
                "meshio++: GridSpec.upscale_samples: factor must be one positive "
                "integer, or three"
            )
        lo, hi = self.bounds
        dims = f * (self.dims + 1) - 1
        spacing = (hi - lo) / dims.astype(np.float64)
        return GridSpec(origin=lo, spacing=spacing, dims=dims)

    def scaling_factor(self, other: "GridSpec") -> Optional[Tuple[int, int, int]]:
        """``other.dims / self.dims`` when it is a whole number on every axis.

        The **cell** ratio, matching :meth:`upscale`. For the number an
        SRResNet's ``scaling_factor`` must equal, use
        :meth:`sample_scaling_factor`.

        ``None`` when the ratio is not whole -- that is a mismatched pair rather
        than a model hyperparameter.
        """
        ratio = other.dims.astype(np.float64) / self.dims.astype(np.float64)
        whole = np.rint(ratio).astype(np.int64)
        if np.any(whole < 1) or np.any(whole.astype(np.float64) != ratio):
            return None
        return (int(whole[0]), int(whole[1]), int(whole[2]))

    def sample_scaling_factor(self, other: "GridSpec") -> Optional[int]:
        """The single integer ``other.shape / self.shape``, or ``None``.

        The number a convolutional upsampler is parametrized by, and the one to
        pass as ``SRResNet(scaling_factor=...)``. It is a **scalar**, not a
        triple, because such a model applies one factor to every axis; a pair
        whose axes disagree is not upsamplable by one and returns ``None``.
        """
        mine = np.asarray(self.shape, dtype=np.float64)
        theirs = np.asarray(other.shape, dtype=np.float64)
        ratio = theirs / mine
        whole = np.rint(ratio).astype(np.int64)
        if np.any(whole < 1) or np.any(whole.astype(np.float64) != ratio):
            return None
        if not np.all(whole == whole[0]):
            return None
        return int(whole[0])

    # -- construction ------------------------------------------------------- #
    @classmethod
    def from_mesh(
        cls,
        mesh,
        *,
        resolution=None,
        cell_size=None,
        bounds=None,
        padding: float = 0.0,
        padding_relative: float = 0.0,
        max_cells: int = DEFAULT_MAX_CELLS,
    ) -> "GridSpec":
        """Resolve a lattice covering ``mesh``.

        The six fields are exactly ``voxelize``'s and ``compute_sdf``'s, resolved
        by the same helper, so "the grid around this mesh" means one thing across
        the library. Give exactly one of ``resolution`` (cell counts, world
        order) and ``cell_size``.
        """
        from ._voxelize import _resolve_lattice

        origin, spacing, dims = _resolve_lattice(
            mesh,
            resolution,
            cell_size,
            bounds,
            padding,
            padding_relative,
            max_cells,
            _PREFIX="meshio++: GridSpec.from_mesh: ",
        )
        return cls(origin=origin, spacing=spacing, dims=dims)

    @classmethod
    def from_lattice_mesh(cls, mesh) -> "GridSpec":
        """Recover the spec of a mesh that *is* a dense lattice.

        Exact rather than a fit, because ``grid`` writes every coordinate
        independently -- see :func:`meshioplusplus._grid.lattice_from_mesh`. A
        mesh that is not a dense lattice raises by name; in particular a *lossy*
        round trip (through ``.vtu``, say) is refused rather than accepted as a
        warped lattice, so cache grids as ``.vti``.
        """
        got = lattice_from_mesh(mesh)
        if got is None:
            raise ValueError(
                "meshio++: GridSpec.from_lattice_mesh: this mesh is not a dense "
                "regular lattice (it must be one hexahedron block whose points tile "
                "a box in x-fastest order); write and read grids as .vti, which "
                "stores the lattice as origin/spacing/extent and reproduces the "
                "points exactly"
            )
        dims, origin, spacing = got
        return cls(origin=origin, spacing=spacing, dims=dims)

    # -- serialization ------------------------------------------------------ #
    def to_dict(self) -> dict:
        """snake_case, for machine-written artefacts (model cards, field_data)."""
        return {
            "origin": [float(v) for v in self.origin],
            "spacing": [float(v) for v in self.spacing],
            "dims": [int(v) for v in self.dims],
            "shape": list(self.shape),
            "layout": GRID_LAYOUT,
        }

    @classmethod
    def from_dict(cls, doc: dict) -> "GridSpec":
        try:
            return cls(origin=doc["origin"], spacing=doc["spacing"], dims=doc["dims"])
        except KeyError as exc:
            raise ValueError(
                f"meshio++: GridSpec.from_dict: missing key {exc.args[0]!r} "
                "(origin, spacing and dims are all required)"
            ) from None

    def to_settings(self) -> dict:
        """PascalCase, for a hand-edited settings document.

        Deliberately a second spelling rather than the same one: a *resolved*
        lattice and the *rule* that resolves one are different objects, and the
        settings family (``DatasetManifest``, ``TrainSpec``) is PascalCase.
        """
        return {
            "Origin": [float(v) for v in self.origin],
            "Spacing": [float(v) for v in self.spacing],
            "Dims": [int(v) for v in self.dims],
        }

    def __eq__(self, other) -> bool:
        if not isinstance(other, GridSpec):
            return NotImplemented
        return bool(
            np.array_equal(self.origin, other.origin)
            and np.array_equal(self.spacing, other.spacing)
            and np.array_equal(self.dims, other.dims)
        )

    def __hash__(self):
        return hash(
            (self.origin.tobytes(), self.spacing.tobytes(), self.dims.tobytes())
        )

    def __repr__(self):
        return (
            f"GridSpec(origin={self.origin.tolist()}, spacing={self.spacing.tolist()}, "
            f"dims={self.dims.tolist()})"
        )


# --------------------------------------------------------------------------- #
# values on one                                                               #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class GridArray:
    """A dense field on a lattice, plus the contract that names its channels.

    ``channels`` is the recorded column order -- the same contract
    :class:`meshioplusplus.FeatureMatrix` carries, produced by the same rule, so
    a multi-component array expands into ``v_0``/``v_1``/``v_2`` here exactly as
    it does there. Record it at training time and assert it at inference time.
    """

    #: ``(C, D, H, W)`` float, channels first then z, y, x.
    values: np.ndarray
    #: the lattice the values live on.
    spec: GridSpec
    #: flat channel names, in ``values``' own order -- the contract.
    channels: Tuple[str, ...]
    #: JSON-serializable: schema version, layout, spec, sources, coverage.
    schema: dict = field(default_factory=dict)

    def __post_init__(self):
        values = np.asarray(self.values)
        if values.ndim != 4:
            raise ValueError(
                f"meshio++: GridArray: values must be (C, D, H, W), got "
                f"{values.ndim} dimensions"
            )
        if tuple(values.shape[1:]) != self.spec.shape:
            raise ValueError(
                f"meshio++: GridArray: values have grid shape {tuple(values.shape[1:])} "
                f"but the spec describes {self.spec.shape} (D, H, W)"
            )
        if len(self.channels) != values.shape[0]:
            raise ValueError(
                f"meshio++: GridArray: {values.shape[0]} channels in values but "
                f"{len(self.channels)} names"
            )
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "channels", tuple(self.channels))

    @property
    def num_channels(self) -> int:
        return int(self.values.shape[0])

    @property
    def coverage(self) -> Optional[float]:
        """Fraction of lattice points inside the sampled source, or ``None``.

        A concave domain gives a grid that is mostly fill, and a model trained on
        it learns the fill. This is the number that makes that visible.
        """
        return self.schema.get("coverage")

    def channel(self, name: str) -> np.ndarray:
        """One channel by name, ``(D, H, W)``."""
        try:
            return self.values[self.channels.index(name)]
        except ValueError:
            raise ValueError(
                f"meshio++: GridArray.channel: no channel named {name!r} "
                f"(have {list(self.channels)})"
            ) from None

    def to_mesh(self) -> Mesh:
        """The lattice as a ``hexahedron`` mesh carrying the values as point data.

        Every writer, ``view``, ``crop``, ``isosurface`` and the browser viewer
        then work on it with no new code -- which is the whole reason the grid is
        a mesh here rather than a bespoke container. Write it as ``.vti`` to keep
        the lattice recoverable.
        """
        mesh = self.spec.mesh()
        for name, values in _regroup(
            self.channels, self.values.reshape(self.num_channels, -1)
        ):
            mesh.point_data[name] = np.ascontiguousarray(values.T)
        return mesh

    @classmethod
    def from_mesh(cls, mesh, *, fields=None) -> "GridArray":
        """Read a lattice mesh's point data back into ``(C, D, H, W)``."""
        spec = GridSpec.from_lattice_mesh(mesh)
        names = sorted(mesh.point_data) if fields is None else list(fields)
        missing = [n for n in names if n not in mesh.point_data]
        if missing:
            raise ValueError(
                f"meshio++: GridArray.from_mesh: no point data named {missing} "
                f"(available: {sorted(mesh.point_data)})"
            )
        channels, columns = [], []
        for name in names:
            arr = np.asarray(mesh.point_data[name], dtype=np.float64)
            flat = arr.reshape(arr.shape[0], -1)
            for c in range(flat.shape[1]):
                channels.append(name if flat.shape[1] == 1 else f"{name}_{c}")
                columns.append(flat[:, c])
        values = np.stack(columns, axis=0).reshape((len(channels),) + spec.shape)
        return cls(
            values=np.ascontiguousarray(values),
            spec=spec,
            channels=tuple(channels),
            schema={
                "grid_schema_version": GRID_SCHEMA_VERSION,
                "layout": GRID_LAYOUT,
                "spec": spec.to_dict(),
                "channels": list(channels),
            },
        )


@dataclass(frozen=True)
class PowerSpectrum:
    """An azimuthally averaged power spectrum, and the bins it was built on."""

    #: bin centres, ``(B,)`` float64, in :attr:`units`.
    wavenumber: np.ndarray
    #: summed power per bin, ``(B,)`` float64. Sums to ``mean(field**2)``.
    power: np.ndarray
    #: modes per bin, ``(B,)`` int64 -- zero where the grid resolves none.
    counts: np.ndarray
    #: what ``wavenumber`` is measured in.
    units: str = "cycles per unit length"


def _regroup(channels, columns):
    """Invert ``feature_matrix``' suffix expansion: ``v_0``,``v_1`` -> ``v``.

    Runs of ``<base>_<k>`` with ``k`` counting from zero rebuild one
    multi-component array; anything else stays a scalar under its own name. The
    grouping is by *consecutive* run rather than by name, so two unrelated
    arrays that happen to share a prefix cannot be fused.
    """
    out = []
    i = 0
    n = len(channels)
    while i < n:
        name = channels[i]
        base, sep, tail = name.rpartition("_")
        if sep and tail == "0" and base:
            j = i + 1
            while j < n and channels[j] == f"{base}_{j - i}":
                j += 1
            if j - i > 1:
                out.append((base, np.stack(columns[i:j], axis=0)))
                i = j
                continue
        out.append((name, columns[i][None, :]))
        i += 1
    return out


# --------------------------------------------------------------------------- #
# mesh -> grid                                                                #
# --------------------------------------------------------------------------- #
def sample_grid(
    mesh,
    spec: GridSpec,
    *,
    fields=None,
    extrapolate: bool = False,
    fill_value: float = 0.0,
    coverage: bool = True,
    float32: bool = False,
) -> GridArray:
    """Sample a mesh's ``point_data`` onto a regular lattice.

    Barycentric interpolation at the lattice's own points, which is exact for a
    field that is linear over each source cell. The lattice is passed to
    :func:`meshioplusplus.interpolate` as a **cell-less** target, so the sampler
    never materializes hexahedron connectivity -- on a 128³ grid that is points
    only rather than points plus eight indices per cell.

    Parameters
    ----------
    mesh :
        the source (unmodified).
    spec :
        the lattice to sample onto. Take its box from the *fine* mesh when
        building a coarse/fine pair -- see the module docstring.
    fields :
        source ``point_data`` names, or ``None`` (default) for all of them in
        sorted order. Naming a ``cell_data`` array raises: a piecewise-constant
        field has no value at a point, so transfer it with
        ``cell_data_to_point_data`` (CLI ``data to-point``) first.
    extrapolate :
        give a lattice point outside the source domain its nearest source
        value instead of ``fill_value``.
    fill_value :
        what an outside point gets when ``extrapolate`` is off. ``float("nan")``
        makes the fill visible to every downstream reduction rather than
        plausible.
    coverage :
        record the fraction of lattice points inside the source domain. Free
        when ``extrapolate`` is off (a constant probe channel rides the same
        call); one extra pass when it is on, since "inside" then has to be
        established separately.
    float32 :
        cast the result to float32, the training convention.

    Returns
    -------
    GridArray
    """
    if not isinstance(spec, GridSpec):
        raise TypeError(f"{_PREFIX}spec must be a GridSpec")

    point_data = dict(getattr(mesh, "point_data", {}) or {})
    cell_data = dict(getattr(mesh, "cell_data", {}) or {})
    if fields is None:
        names = sorted(point_data)
    else:
        names = list(fields)
        for name in names:
            if name in point_data:
                continue
            if name in cell_data:
                raise ValueError(
                    f"{_PREFIX}{name!r} is cell data, and a piecewise-constant field "
                    "has no value at a point; convert it first with "
                    "cell_data_to_point_data (CLI: data to-point)"
                )
            raise ValueError(
                f"{_PREFIX}no point data named {name!r} "
                f"(available: {sorted(point_data)})"
            )
    if not names:
        raise ValueError(
            f"{_PREFIX}the source mesh carries no point data to sample "
            "(cell data must be converted with data to-point first)"
        )

    target = Mesh(spec.points(), [])
    probe = coverage and not extrapolate
    source = mesh
    if probe:
        # A shallow copy: the arrays and cell blocks are shared, only the dict is
        # new, so the caller's mesh is untouched and nothing large is duplicated.
        source = Mesh(mesh.points, mesh.cells, point_data=dict(point_data))
        source.point_data[_PROBE] = np.ones(len(mesh.points), dtype=np.float64)

    from ._interpolate import interpolate

    sampled = interpolate(
        source,
        target,
        method="barycentric",
        arrays=names + ([_PROBE] if probe else []),
        extrapolate=extrapolate,
        default_value=float(fill_value),
    )

    inside = None
    if probe:
        inside = np.asarray(sampled.point_data.pop(_PROBE), dtype=np.float64) > 0.5
    elif coverage:
        mask_target = Mesh(spec.points(), [])
        mask_source = Mesh(mesh.points, mesh.cells, point_data={})
        mask_source.point_data[_PROBE] = np.ones(len(mesh.points), dtype=np.float64)
        probed = interpolate(
            mask_source,
            mask_target,
            method="barycentric",
            arrays=[_PROBE],
            extrapolate=False,
            default_value=0.0,
        )
        inside = np.asarray(probed.point_data[_PROBE], dtype=np.float64) > 0.5

    from ._ml import feature_matrix

    fm = feature_matrix(sampled, "point", fields=names, coords=False, regions=False)
    values = np.ascontiguousarray(fm.matrix.T).reshape(
        (fm.matrix.shape[1],) + spec.shape
    )
    if float32:
        values = np.ascontiguousarray(values, dtype=np.float32)

    from .__about__ import __version__

    schema = {
        "grid_schema_version": GRID_SCHEMA_VERSION,
        "meshioplusplus_version": str(__version__),
        "layout": GRID_LAYOUT,
        "spec": spec.to_dict(),
        "channels": list(fm.columns),
        "sources": fm.schema["sources"],
        "fields": list(names),
        "extrapolate": bool(extrapolate),
        "fill_value": float(fill_value),
    }
    if inside is not None:
        schema["coverage"] = float(np.count_nonzero(inside)) / float(inside.size)
    return GridArray(
        values=values, spec=spec, channels=tuple(fm.columns), schema=schema
    )


# --------------------------------------------------------------------------- #
# grid -> anywhere                                                            #
# --------------------------------------------------------------------------- #
def interpolate_grid(values, spec: GridSpec, points) -> np.ndarray:
    """Evaluate a grid at arbitrary points by trilinear interpolation.

    Queries outside the box are **clamped** to it rather than extrapolated --
    a trilinear form extrapolates as a product of linears and diverges fast, and
    with the shared-box rule an outside query is a small numerical excursion at
    a face rather than a real request for the field beyond the domain.

    Parameters
    ----------
    values :
        ``(C, D, H, W)`` or ``(D, H, W)``.
    spec :
        the lattice ``values`` live on.
    points :
        ``(M, 3)`` query coordinates.

    Returns
    -------
    numpy.ndarray
        ``(C, M)`` float64, or ``(M,)`` when ``values`` was ``(D, H, W)``.
    """
    arr = np.asarray(values)
    bare = arr.ndim == 3
    if bare:
        arr = arr[None, ...]
    if arr.ndim != 4:
        raise ValueError(
            "meshio++: interpolate_grid: values must be (C, D, H, W) or (D, H, W)"
        )
    if tuple(arr.shape[1:]) != spec.shape:
        raise ValueError(
            f"meshio++: interpolate_grid: values have grid shape "
            f"{tuple(arr.shape[1:])} but the spec describes {spec.shape} (D, H, W)"
        )
    arr = np.asarray(arr, dtype=np.float64)

    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] < 2:
        raise ValueError("meshio++: interpolate_grid: points must be (M, 2) or (M, 3)")
    if pts.shape[1] == 2:
        pts = np.column_stack([pts, np.zeros(len(pts), dtype=np.float64)])
    pts = pts[:, :3]

    nx, ny, nz = (int(v) for v in spec.dims)
    # Continuous index in cell units, clamped into the box.
    t = (pts - spec.origin) / spec.spacing
    t[:, 0] = np.clip(t[:, 0], 0.0, float(nx))
    t[:, 1] = np.clip(t[:, 1], 0.0, float(ny))
    t[:, 2] = np.clip(t[:, 2], 0.0, float(nz))
    i0 = np.floor(t).astype(np.int64)
    i0[:, 0] = np.clip(i0[:, 0], 0, nx - 1)
    i0[:, 1] = np.clip(i0[:, 1], 0, ny - 1)
    i0[:, 2] = np.clip(i0[:, 2], 0, nz - 1)
    f = t - i0.astype(np.float64)
    ix, iy, iz = i0[:, 0], i0[:, 1], i0[:, 2]
    fx, fy, fz = f[:, 0], f[:, 1], f[:, 2]

    def corner(dx, dy, dz):
        return arr[:, iz + dz, iy + dy, ix + dx]

    # Fixed accumulation order: blend x, then y, then z.
    gx = 1.0 - fx
    c00 = corner(0, 0, 0) * gx + corner(1, 0, 0) * fx
    c10 = corner(0, 1, 0) * gx + corner(1, 1, 0) * fx
    c01 = corner(0, 0, 1) * gx + corner(1, 0, 1) * fx
    c11 = corner(0, 1, 1) * gx + corner(1, 1, 1) * fx
    gy = 1.0 - fy
    c0 = c00 * gy + c10 * fy
    c1 = c01 * gy + c11 * fy
    out = c0 * (1.0 - fz) + c1 * fz
    return out[0] if bare else out


def resample_grid(values, spec: GridSpec, target: GridSpec) -> np.ndarray:
    """Evaluate a grid on another lattice -- the upsampling baseline.

    This is what a superresolution model has to beat, and what stands in for the
    model in a test of the whole chain.

    Resampling a grid onto **its own** spec returns a copy without evaluating
    anything. That is not a shortcut around a numerical problem: trilinear
    interpolation at the sample points is the identity by definition, and
    evaluating it would only introduce the rounding of recovering an index from
    a coordinate.
    """
    arr = np.asarray(values)
    if spec == target:
        return arr.copy()
    out = interpolate_grid(arr, spec, target.points())
    if arr.ndim == 3:
        return np.ascontiguousarray(out.reshape(target.shape))
    return np.ascontiguousarray(out.reshape((arr.shape[0],) + target.shape))


def scatter_grid(
    array: GridArray,
    mesh,
    *,
    names=None,
    on_conflict: str = "error",
    suffix: str = "_grid",
) -> Mesh:
    """Write a grid's channels onto a mesh's points, by trilinear interpolation.

    The inverse of :func:`sample_grid`, and the step that turns a model's output
    back into something every format can hold.

    Parameters
    ----------
    array :
        the values and their channel contract.
    mesh :
        the target (unmodified); a copy is returned.
    names :
        override the array names written. By default the channel contract is
        inverted -- consecutive ``v_0``/``v_1``/``v_2`` runs rebuild one
        multi-component array ``v`` -- so a round trip through a grid preserves
        an array's shape as well as its values.
    on_conflict :
        ``"error"`` (default), ``"overwrite"``, or ``"suffix"``.
    suffix :
        what ``"suffix"`` appends.
    """
    if on_conflict not in ("error", "overwrite", "suffix"):
        raise ValueError(
            "meshio++: scatter_grid: on_conflict must be 'error', 'overwrite' or "
            f"'suffix', got {on_conflict!r}"
        )
    pts = np.asarray(mesh.points, dtype=np.float64)
    flat = interpolate_grid(array.values, array.spec, pts)

    if names is None:
        grouped = _regroup(array.channels, flat)
    else:
        names = list(names)
        if len(names) != array.num_channels:
            raise ValueError(
                f"meshio++: scatter_grid: {array.num_channels} channels but "
                f"{len(names)} names"
            )
        grouped = [(n, flat[i][None, :]) for i, n in enumerate(names)]

    out = _clone(mesh)
    for name, columns in grouped:
        key = name
        if key in out.point_data:
            if on_conflict == "error":
                raise ValueError(
                    f"meshio++: scatter_grid: the target already has point data named "
                    f"{key!r} (pass on_conflict='overwrite' or 'suffix')"
                )
            if on_conflict == "suffix":
                key = f"{name}{suffix}"
                if key in out.point_data:
                    raise ValueError(
                        f"meshio++: scatter_grid: both {name!r} and {key!r} already "
                        "exist on the target"
                    )
        values = columns.T
        out.point_data[key] = np.ascontiguousarray(
            values[:, 0] if values.shape[1] == 1 else values
        )
    return out


def _clone(mesh) -> Mesh:
    """A copy that shares nothing mutable at the top level."""
    out = Mesh(
        np.asarray(mesh.points),
        [(b.type, b.data) for b in mesh.cells],
        point_data=dict(getattr(mesh, "point_data", {}) or {}),
        cell_data={
            k: list(v) for k, v in (getattr(mesh, "cell_data", {}) or {}).items()
        },
        field_data=dict(getattr(mesh, "field_data", {}) or {}),
    )
    for region in getattr(mesh, "regions", []) or []:
        out.regions.append(region)
    return out


# --------------------------------------------------------------------------- #
# the thin-axis idiom                                                         #
# --------------------------------------------------------------------------- #
def _tensor_axis(world_axis: int) -> int:
    """World axis (0=x, 1=y, 2=z) -> tensor axis in a ``(C, D, H, W)`` array."""
    if world_axis not in (0, 1, 2):
        raise ValueError(
            f"meshio++: grid: axis must be 0 (x), 1 (y) or 2 (z), got {world_axis}"
        )
    return 3 - world_axis


def squeeze_grid(
    values, axis: int = 2, index: Optional[int] = None, reduce: Optional[str] = None
) -> np.ndarray:
    """Drop a **world** axis, turning a 3-D grid into a 2-D one.

    How a 2-D operator (FNO, AFNO, a 2-D U-Net) is fed a mesh that is planar in
    one direction. ``axis`` is a world axis (0=x, 1=y, 2=z), like every other
    three-vector here, not the tensor axis it maps to.

    A lattice always has at least two points on every axis -- one cell has two
    corners -- so the thin axis of a planar problem still arrives with a plane at
    each face. Say which of them you want:

    * ``index=k`` keeps plane ``k`` and discards the rest;
    * ``reduce="mean"`` averages over the axis, which is the right choice when
      the axis is thin but not exactly constant, and what a solver coupling
      normally does;
    * neither is required only when the axis is already a single plane.

    One of them is **required** whenever the axis is longer than one, because
    silently taking the first plane would discard the others without saying so.

    A free function on raw arrays rather than a :class:`GridArray` method,
    because the result no longer matches a three-dimensional spec -- that is the
    point of it, and pretending otherwise would let a squeezed array be scattered
    back through the wrong geometry.
    """
    arr = np.asarray(values)
    ax = _tensor_axis(axis)
    if arr.ndim != 4:
        raise ValueError("meshio++: squeeze_grid: values must be (C, D, H, W)")
    if index is not None and reduce is not None:
        raise ValueError("meshio++: squeeze_grid: give at most one of index and reduce")
    if reduce is not None:
        if reduce != "mean":
            raise ValueError(
                f"meshio++: squeeze_grid: reduce must be 'mean', got {reduce!r}"
            )
        return np.ascontiguousarray(arr.mean(axis=ax))
    n = arr.shape[ax]
    if index is None:
        if n != 1:
            raise ValueError(
                f"meshio++: squeeze_grid: world axis {axis} has {n} points, so there "
                "is no single plane to keep; pass index= to choose one, or "
                "reduce='mean' to average over them"
            )
        index = 0
    if not -n <= index < n:
        raise ValueError(
            f"meshio++: squeeze_grid: index {index} is out of range for world axis "
            f"{axis}, which has {n} points"
        )
    return np.ascontiguousarray(np.take(arr, index, axis=ax))


def expand_grid(values, axis: int = 2, size: int = 1) -> np.ndarray:
    """Reinsert a **world** axis -- :func:`squeeze_grid`'s inverse.

    ``size`` duplicates the plane that many times, which is how a 2-D model's
    prediction is written back onto a thin 3-D lattice: the model saw one plane,
    and every plane of the thin axis gets it.
    """
    arr = np.asarray(values)
    if arr.ndim != 3:
        raise ValueError("meshio++: expand_grid: values must be (C, H, W)")
    if size < 1:
        raise ValueError(f"meshio++: expand_grid: size must be at least 1, got {size}")
    ax = _tensor_axis(axis)
    out = np.expand_dims(arr, axis=ax)
    if size > 1:
        out = np.repeat(out, size, axis=ax)
    return np.ascontiguousarray(out)


# --------------------------------------------------------------------------- #
# spectra                                                                     #
# --------------------------------------------------------------------------- #
def power_spectrum(values, spec: GridSpec) -> PowerSpectrum:
    """The azimuthally averaged power spectrum of a field on a lattice.

    The honest measure of whether a super-resolved or generated field carries the
    right *small-scale content*, as opposed to a plausible-looking one with the
    high wavenumbers smoothed away -- a pointwise error does not see the
    difference, and a model can score well on it while producing a field that is
    physically wrong at every scale that matters.

    Parameters
    ----------
    values :
        ``(D, H, W)`` for a scalar field, or ``(C, D, H, W)`` for a vector one,
        in which case the power of the components is summed -- the energy
        spectrum, which is what a velocity field wants.
    spec :
        the lattice. An **anisotropic** one raises: averaging over shells is only
        meaningful when a shell means the same thing on every axis.

    Returns
    -------
    PowerSpectrum
        Bins cover every radius the grid resolves, including past the isotropic
        Nyquist where the corners of the box live, so the power sums exactly to
        ``mean(field**2)`` (Parseval). Truncating would break that; ``counts``
        is there so a caller can cut the tail themselves.
    """
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim == 3:
        arr = arr[None, ...]
    if arr.ndim != 4:
        raise ValueError(
            "meshio++: power_spectrum: values must be (C, D, H, W) or (D, H, W)"
        )
    if tuple(arr.shape[1:]) != spec.shape:
        raise ValueError(
            f"meshio++: power_spectrum: values have grid shape {tuple(arr.shape[1:])} "
            f"but the spec describes {spec.shape} (D, H, W)"
        )
    if not spec.is_isotropic:
        raise ValueError(
            "meshio++: power_spectrum: an azimuthal average needs an isotropic "
            f"lattice, but the spacing is {list(spec.spacing)}; resample onto an "
            "isotropic grid first"
        )

    nd, nh, nw = spec.shape
    total = float(nd * nh * nw)
    h = float(spec.spacing[0])

    # Power, normalized so that summing it gives mean(field**2) exactly.
    power = np.zeros((nd, nh, nw), dtype=np.float64)
    for c in range(arr.shape[0]):
        f = np.fft.fftn(arr[c])
        power += (f.real * f.real + f.imag * f.imag) / (total * total)

    # Physical wavenumber per mode, in cycles per unit length.
    kz = np.fft.fftfreq(nd, d=h)
    ky = np.fft.fftfreq(nh, d=h)
    kx = np.fft.fftfreq(nw, d=h)
    kkz, kky, kkx = np.meshgrid(kz, ky, kx, indexing="ij")
    kmag = np.sqrt(kkx * kkx + kky * kky + kkz * kkz)

    # Bin width: the finest fundamental of the three axes, so every axis' own
    # modes land on a bin boundary. A grid whose axes differ in length simply
    # leaves some bins empty, which `counts` reports rather than hides.
    dk = 1.0 / (float(max(nd, nh, nw)) * h)
    idx = np.rint(kmag / dk).astype(np.int64)
    nbins = int(idx.max()) + 1
    binned = np.bincount(idx.reshape(-1), weights=power.reshape(-1), minlength=nbins)
    counts = np.bincount(idx.reshape(-1), minlength=nbins).astype(np.int64)
    wavenumber = np.arange(nbins, dtype=np.float64) * dk
    return PowerSpectrum(
        wavenumber=wavenumber,
        power=binned,
        counts=counts,
        units="cycles per unit length",
    )
