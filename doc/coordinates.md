# Cylindrical and spherical coordinates

`convert_coordinates(mesh, system, axis=..., arrays=[...])` reads a mesh in cylindrical or spherical components about an axis: it records each point's `(r, θ, z)` or `(r, θ, φ)` as point data, and rewrites the named vector and tensor fields in the local basis of each point (a cell's at its centroid). The results of rotating machinery and pressure vessels are read this way: the radial, tangential and axial velocity of a rotor, the hoop and axial stress of a vessel. The points do not move; `inverse=True` turns local components back into Cartesian ones.

```python
import meshioplusplus as mp

rotor = mp.read("rotor.vtu")
out = mp.convert_coordinates(rotor, "cylindrical", axis=(0, 0, 1), arrays=["U", "stress"])
u_r, u_theta, u_z = out.point_data["U"].T
r, theta, z = out.point_data["coords:cylindrical"].T

vessel = mp.read("vessel.vtu")
out, report = mp.convert_coordinates(
    vessel, "spherical", origin=(0, 0, 1.5), arrays=["sigma"], return_report=True
)
report["num_arrays"], report["num_singular"]

back = mp.convert_coordinates(out, "spherical", origin=(0, 0, 1.5), arrays=["sigma"],
                              inverse=True, record_coordinates=False)
```

## The systems

With `a` the unit axis through the origin `o`, `d = p − o`, `z = d · a`, the radial vector `ρ = d − z a`, and `e₁` the `reference` direction projected off the axis (`(1, 0, 0)` by default; `(0, 1, 0)` stands in when the reference lies along the axis):

| System | `coords:<system>` | Local basis (rows of `R`) |
| --- | --- | --- |
| `cylindrical` | `(r, θ, z)`: `r = |ρ|`, `θ` the angle from `e₁` to `ρ` about `a`, in `(−π, π]` | `e_r = ρ / r`, `e_θ = a × e_r`, `e_z = a` |
| `spherical` | `(r, θ, φ)`: `r = |d|`, `θ` the polar angle from `a`, in `[0, π]`, `φ` the cylindrical `θ` | `e_r = d / r`, `e_θ = cos θ e_ρ − sin θ a`, `e_φ = a × e_ρ` |

The spherical system follows the physics convention (ISO 80000-2): `θ` is measured from the axis, `φ` about it. Both bases are right-handed and orthonormal.

A vector `v` becomes `R v` and a tensor `T` becomes `R T Rᵀ`; with `inverse=True` they become `Rᵀ v` and `Rᵀ T R`, so the two directions round trip to rounding. The trace and the norm are unchanged, which makes a quick check: the von Mises stress of a converted tensor equals the original's.

## Scope

- **Arrays** are point or cell data, by name; a name present as both converts both. A vector is `(n, 3)`; a tensor `(n, 9)` (row-major) or `(n, 3, 3)`. Anything else, or a missing name, is an error. A converted array becomes Float64; every other array passes through.
- **Cell data** is converted in the basis at the cell's centroid: the mean of its distinct nodes.
- **On the axis** (`r = 0` cylindrically, or the spherical centre) the radial direction is undefined: `e₁` stands in for it, the angle is 0, and the point counts in `num_singular`. A point on the axis but away from the spherical centre is singular too, because its `φ` is undefined, although its `e_r` (along the axis) is still exact.
- **Data-only**: points, connectivity, block order and regions are left bit-identical. 2-D points are read with `z = 0`.

Both engines, the C++ core and the numpy twin used when it declines, give the same values to rounding (tested to 1e-13).

## Report

`num_arrays` (the arrays converted, a point and a cell array of one name counting twice) and `num_singular` (the points on the axis or at the centre; 0 with `record_coordinates=False`).

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `convert_coordinates(mesh, system="cylindrical", *, axis=(0, 0, 1), origin=(0, 0, 0), reference=(1, 0, 0), arrays=(), inverse=False, record_coordinates=True, return_report=False)` |
| Pipeline | `{"Op": "ConvertCoordinates", "System": "cylindrical", "Axis": [0, 0, 1], "Origin": [0, 0, 0], "Reference": [1, 0, 0], "Arrays": ["U"], "Inverse": false, "RecordCoordinates": true}`; counters `NumArrays`, `NumSingular` |
| Both CLIs | `meshioplusplus coordinates IN OUT [--system cylindrical\|spherical] [--axis X,Y,Z] [--origin X,Y,Z] [--reference X,Y,Z] [--arrays A,B] [--inverse] [--no-record] [--json]` |
| C | `mio_coordinate_opts_init(&o); o.system = MIO_COORDS_SPHERICAL; o.arrays = names; o.num_arrays = 1; mio_mesh* out = mio_convert_coordinates(mesh, &o, &num_arrays, &num_singular);` |
| Fortran | `out = m%convert_coordinates(system="spherical", arrays=[character(8) :: "U"], num_singular=ns, stat=ierr)` |
| Julia | `r = convert_coordinates(mesh; system="cylindrical", arrays=["U"])`; `r.mesh`, `r.num_arrays`, `r.num_singular` |
| R | `r <- mio_convert_coordinates(mesh, "cylindrical", arrays = "U")`; `r$mesh`, `r$num_arrays`, `r$num_singular` |
| WASM | `m.convertCoordinates(mesh, {system, axis, origin, reference, arrays, inverse, recordCoordinates})` returns `{mesh, numArrays, numSingular}` |
| MCP | tool `convert_coordinates(input_path, output_path, system, axis, origin, reference, arrays, inverse, record_coordinates)` |
