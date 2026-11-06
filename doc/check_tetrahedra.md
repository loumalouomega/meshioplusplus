# Tetrahedron check

`check_tetrahedra(mesh)` reports whether the tetrahedra of a mesh are usable: how many are inverted, how many are degenerate, and how many pairs overlap across a shared face. It is the check behind the volume fix of [`octree_mesh`](octree_mesh.md#the-tetrahedral-output), as a standalone operation on any mesh, and it leaves the mesh unchanged.

```python
import meshioplusplus as mp

mesh = mp.read("model.vtu")

report = mp.check_tetrahedra(mesh)
report["valid"]           # False
report["num_inverted"]    # 3
report["num_overlaps"]    # 0
```

## What it reports

The `tetra` and `tetra10` cells are examined, a quadratic tetrahedron by its four corners; other cells are ignored, and a mesh without tetrahedra reports zeros with a warning.

| Defect | Rule | Counted as |
|---|---|---|
| inverted | the signed volume is below minus the tolerance | `num_inverted` |
| degenerate | the volume is within the tolerance, either sign | `num_degenerate` |
| overlap | two tetrahedra share a face, no other tetrahedron uses it, and their opposite nodes lie on the same side of it | `num_overlaps`, per pair |

The tolerance of a tetrahedron is `zero_volume_factor` (default `1e-6`) times the cube of its longest edge, plus the effect of six significant digits of written coordinates (`1e-6` times its largest coordinate times its longest edge squared), so a sliver of a fine mesh is not flagged because of a large coordinate offset. A face used by three or more tetrahedra is not examined. A conforming mesh of a valid solid has no overlap whatever its materials, since the tetrahedra on the two sides of an interface have their opposite nodes on opposite sides.

The report also has `num_tetrahedra`, `min_volume` and `max_volume` (signed, 0 without tetrahedra), `valid` (no inverted or degenerate tetrahedron and no overlap) and `warnings`.

## Where it is

| Where | Name |
|---|---|
| Python | `meshioplusplus.check_tetrahedra(mesh, zero_volume_factor=1e-6)`, with a numpy twin when the compiled core is missing |
| CLI | `check-tetrahedra INFILE [--zero-volume-factor X] [--json]`: exit 0 valid, 1 not valid, 2 could not check |
| Pipeline | the `CheckTetrahedra` step: `ZeroVolumeFactor` (1e-6), `Require` (true); a gate by default, counters `NumTetrahedra`, `NumInverted`, `NumDegenerate`, `NumOverlaps` |
| MCP | the `check_tetrahedra` tool |
| C++ | `check_tetrahedra(mesh, TetrahedraCheckOptions)` in `meshioplusplus/operations/check_tetrahedra.hpp` |
| C | `mio_check_tetrahedra` ([C API](c_api.md)) |
| Fortran, Julia, R, Octave / MATLAB, WASM | [Fortran](fortran.md), [Julia](julia.md), [R](r.md), [Octave](octave.md), [WASM](wasm.md) |
