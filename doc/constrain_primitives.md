# Constraining fitted primitives

`constrain_primitives(skin)` fits related primitives of a [`fit_primitives`](./fit_primitives.md) skin together, under control (roadmap 7.7.3): the pass [`fit_cad`](./fit_cad.md#constrained-primitives)'s `constrain_primitives` runs — coaxial groups, planes made parallel or perpendicular to an axis or to each other, equal radii — with a force and a forbid list of relations, angle and distance tolerances in place of its own, and a report of every relation applied with its residual.

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

constrained, report = meshioplusplus.constrain_primitives(
    skin,
    angle_tolerance=0.5,              # snap within half a degree
    forbid=[(3, 7, "perpendicular")],  # by segment:patch label
    force=[(1, 4, "coaxial")],
    return_report=True,
)
[(r["patch_a"], r["patch_b"], r["relation"], r["angle_deg"]) for r in report["relations"]]
solid = meshioplusplus.fit_cad(constrained, constrain_primitives=False)
```

| Option | Default | Meaning |
|---|---|---|
| `face_tolerance` | three times the largest `fit:patch_tolerance` | how far a face may lie from its patch before a writer demotes it (as `fit_cad`'s default) |
| `angle_tolerance` | the pass's own | the largest angle, in degrees, between two directions still related (the pass's own: what moves a patch by twice its tolerance over its extent) |
| `distance_tolerance` | the pass's own | the largest distance between two axes (or of a centre from an axis) still coaxial, and between two radii still equal (the pass's own: twice the patch tolerance; the smaller tolerance) |
| `min_fraction` | `fit_primitives`' 0.9 | the share of a patch's points that must stay within its tolerance |
| `force`, `forbid` | none | `(patch_a, patch_b, relation)`, relation `coaxial`, `parallel`, `perpendicular` or `equal_radius` |

A forbidden relation is never applied (a forbidden coaxial pair keeps either out of the other's group). A forced one is tried whatever the tolerances say, and kept only when the patches still pass the test they were accepted under, as every change is; a forced relation that fails is listed in `unmet`. The relations are reported as `(patch_a, patch_b, relation, angle_deg, distance, forced)`: `patch_a` the reference (a coaxial group's seed, the axis or earlier plane a plane took its normal from, the first of an equal-radius set), the residual measured before the change (the angle between the directions, for a perpendicular its distance from 90 degrees; a coaxial member's axis point from the seed's axis, two radii's difference). With the default options the result is the pass `fit_cad` runs, bit for bit (`tests/cpp/test_constrain_primitives.cpp`), and forbidding every relation of a plane leaves its fitted normal unchanged.

Run it before `fit_cad(constrain_primitives=False)`; `fit_cad`'s own pass keeps its defaults. The pass's limits stay as [recorded](./fit_cad.md#constrained-primitives): parallel but offset axes, angles other than right ones and patterns are not related.

## Surfaces

C++ `constrain_primitives`, `primitive_relation_from_name` and `primitive_relation_name` (`operations/constrain_primitives.hpp`); Python `meshioplusplus.constrain_primitives` (C++ core only); both CLIs (`constrain-primitives IN OUT [--face-tolerance T] [--angle-tolerance DEG] [--distance-tolerance D] [--min-fraction F] [--force A:B:REL] [--forbid A:B:REL] [--json]`, the relation flags repeatable); the `ConstrainPrimitives` [pipeline](./pipeline.md) step (`Force` and `Forbid` as `"A:B:REL"` texts separated by commas); C `mio_constrain_primitives` (`mio_constrain_primitives_opts`, 96 bytes, relations as `(a, b, mio_primitive_relation)` triples; `mio_constrain_primitives_report`, 64; `mio_constrain_primitives_relation`, `mio_constrain_primitives_unmet`); Fortran, Julia, R, Octave/MATLAB, WASM `constrainPrimitives`; the MCP tool `constrain_primitives`.
