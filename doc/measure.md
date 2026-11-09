# Measuring a scan: elements, dimensions and GD&T

`measure(mesh, spec)` turns a scanned surface into the numbers a first-article inspection report carries: the elements a drawing refers to (planes, bores, spheres, their axes and the lines and points built from them), the dimensions between them, and the GD&T form, orientation and location tolerances, each against a nominal or a limit with PASS or FAIL (roadmap 7.2). `align_datums` then lays the part on its datums, the 3-2-1 alignment a CMM starts from (roadmap 7.3.4), and `deviation` measures against any element (roadmap 7.4.2).

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

spec = {
    "elements": [
        {"name": "A", "fit": "plane", "region": "bottom"},
        {"name": "B", "fit": "plane", "seeds": [[40.0, 5.0, 10.0]]},
        {"name": "bore", "fit": "cylinder", "region": "hole"},
        {"name": "axis", "construct": "axis", "of": "bore"},
    ],
    "dimensions": [
        {"name": "bore_to_B", "type": "distance", "between": ["axis", "B"], "nominal": 25, "tolerance": 0.05},
    ],
    "tolerances": [
        {"name": "flat_A", "type": "flatness", "feature": "A", "limit": 0.02},
        {"name": "perp_B", "type": "perpendicularity", "feature": "B", "datum": "A", "limit": 0.05},
        {"name": "perp_bore", "type": "perpendicularity", "feature": "axis", "datum": "A", "limit": 0.05},
    ],
    "criterion": "minimum_zone",
}
measured, report = meshioplusplus.measure(scan, spec)
report["num_fail"], [t["value"] for t in report["tolerances"]]

aligned, matrix = meshioplusplus.align_datums(measured, "A", "B", origin="C")
dev = meshioplusplus.deviation(measured, measured, element="A")
```

## Design: the element table

The roadmap asked four questions before any code (7.2.1); this is how they are answered.

- **What an element is.** A point, a line, a plane, a sphere, a cylinder or a cone: a point and a unit direction (a plane's normal, an axis), a radius or a half-angle where they apply, the form error under the spec's criterion for a measured one, and for a measured one the points of its cells. Measured elements are fits; constructed ones are built from earlier elements; a literal point is given.
- **Where it lives.** As the typed `MeasureElement` list a run returns (the API), and on the mesh as the field data `measure:table`, the table's JSON text as UInt8 bytes, so a pipeline step, a later run or a file in a format that keeps field data (VTU, VTK, XDMF, ...) carries it; the fitted cells are the cell regions `measure:<name>`. A side channel like `NurbsInfo` would not survive a pipeline step or a write; the field-data copy holds no points, only the geometry.
- **How a step refers to one.** By name, as regions are: names are unique in a spec, and a construction, dimension or tolerance names the elements it takes.
- **Formula or value.** Both: a constructed element records its construction and its inputs (`source`, `inputs`) beside its value, and a run evaluates the whole spec in order, so refitting an element refits everything built on it.

## The spec

A JSON object (a dict in Python, its text, or the path of a JSON file).

### Elements

| Form | Keys | Result |
|---|---|---|
| fitted | `"fit"`: `plane`, `sphere`, `cylinder`, `cone` or `auto`; one of `"region"` (a cell region), `"seeds"` (`[[x, y, z], ...]`, grown as [fit primitives](./fit_primitives.md#a-feature-from-a-seed-or-a-region) grows them) or `"patch"` (a `segment:patch` label); optional `"tolerance"` (relative, `1e-3`), `"absolute"`, `"min_inlier_fraction"`, `"feature_angle"`, `"rim_angle"` | the fitted primitive; its cells become the region `measure:<name>` |
| constructed | `"construct"` and `"of"` (element names) | see below |
| literal | `"point"`: `[x, y, z]` | a point |

A fitted element runs `fit_primitives` on its cells (sigma trimming in the region and seed modes), under the spec's `"criterion"` (`least_squares` by default, or `minimum_zone`, `min_circumscribed`, `max_inscribed`; see [fit primitives](./fit_primitives.md#minimum-zone-and-envelope-fits)); a fit no primitive of the type explains within the tolerance is refused by name.

| Construction | Takes | Gives |
|---|---|---|
| `intersection` | two planes; a line (or axis) and a plane; two lines | a line; a point; the midpoint of the lines' closest points |
| `axis` | a cylinder, a cone or a line | its axis, a line (with the feature's extent along it) |
| `centre` | a sphere, a cylinder or a point | its centre (a cylinder's axis point nearest its points' mean) |
| `apex` | a cone | its apex |
| `mid_plane` | two planes | the plane midway, its normal the mean of theirs |
| `projection` | a point or a line, and a plane | the point or the line projected onto the plane |
| `plane_through` | three point-like elements | the plane through them |

A cylinder and a cone stand for their axes and a sphere for its centre wherever a line or a point is expected.

### Dimensions

`{"name", "type": "distance" | "angle", "between": [a, b], "nominal", "tolerance"}`. A distance is between points (their distance), a point and a plane or a line (perpendicular), two lines (the distance between them, for parallel lines the distance of one's point from the other), a line and a plane (0 when they meet, else the distance), or two planes (the distance of the second's point from the first, with a warning when they are not parallel). An angle is the acute angle between two lines or two planes' normals, or between a line and a plane (90 minus the angle to its normal), in degrees. With a `nominal`, the report gives `deviation = value - nominal`; with a `tolerance` as well, PASS when `|deviation| <= tolerance`.

### Tolerances

`{"name", "type", "feature", "datum", "limit", ...}`, PASS when the value is at most the `limit`. The definitions follow ISO 1101 and ASME Y14.5; the value is the width (or diameter) of the narrowest zone of the tolerance's shape, oriented as the tolerance says, that holds the feature.

| Tolerance | Feature (datum) | Value |
|---|---|---|
| `flatness`, `cylindricity`, `sphericity` | a measured plane, cylinder, sphere | the form error: the width of the zone about the fit under the spec's criterion (the minimum zone with `"criterion": "minimum_zone"`) |
| `circularity` | a measured cylinder; `"sections"` (10) | the largest radial spread about the fitted axis over that many sections along it |
| `parallelism`, `perpendicularity`, `angularity` (with `"basic_angle"` in degrees) | a measured plane (a plane or an axis) | the spread of the plane's points along the zone normal: the datum normal turned towards the feature's normal by the basic angle (0 parallel, 90 perpendicular) for a plane datum, by 90 minus it for an axis |
| the same | an axis (a plane or an axis) | between its ends (its points' extent along it): a cylindrical zone along the datum direction for an axis perpendicular to a plane or parallel to an axis (the ends' offset across it), else two planes at the basic angle (the ends' offset along their normal) |
| `coaxiality` | an axis (an axis or a point) | twice the larger distance of the axis' ends from the datum |
| `concentricity` | a centre (an axis or a point) | twice its distance from the datum |

A zone's orientation about the datum is the feature's own: for a plane perpendicular to a datum plane, the zone planes contain the datum normal and face the way the feature's fitted normal does (projected). Position, runout and datum reference frames with more than one datum are the next additions once this is used.

### The report

`elements` (name, kind, point, dir, radius, angle, form_error, measured, source, inputs, num_points), `dimensions` and `tolerances` (name, type, of, value, nominal, tolerance, deviation, pass), `num_fail` and `warnings`. The table on the mesh is the `elements` list without the points.

## Datum alignment

`align_datums(mesh, primary, secondary, origin=None)` reads the mesh's element table: the primary (a plane's normal, or an axis) sets Z, the secondary (a plane's normal or an axis, made perpendicular to Z) sets X, Y completes them, and the origin is a point-like element, the meeting point of the primary, the secondary and a third plane named as `origin`, an axis' meeting point with that plane, or the primary's point. It returns the moved mesh, its element table moved with it, and the row-major 4 × 4 matrix [`transform`](./transform.md) applies. A block fitted on three faces lands with those faces on the coordinate planes to `1e-12`.

## Against an element

`deviation(scan, reference, element="A")` measures the scan's points against the element `A` of the reference's table (often the scan itself, measured): a plane, a sphere, a cylinder or a cone, as [deviation against a primitive](./deviation.md) does. A point or a line has no surface and is refused.

## Surfaces

C++ `measure`, `measure_table`, `measure_element`, `measure_element_primitive` and `align_datums` (`operations/measure.hpp`); Python `meshioplusplus.measure`, `measure_table`, `align_datums` and `deviation(element=)` (C++ core only; `measure_table` reads the table in Python); both CLIs (`measure IN SPEC.json [OUT] [--json]`, exit 1 when a check fails; `align-datums IN OUT --primary A --secondary B [--origin C]`; `deviation --element NAME`); the `Measure` (`Spec`, a path or the JSON text, and `Require`), `AlignDatums` (`Primary`, `Secondary`, `Origin`) [pipeline](./pipeline.md) steps and `Deviation`'s `Element`; C `mio_measure`, `mio_measure_report_json`, `mio_align_datums`, `mio_measure_element_primitive`; Fortran, Julia, R, Octave/MATLAB and WASM (`measure`, `alignDatums`, `measureElementPrimitive`); the MCP tools `measure`, `align_datums` and `deviation`'s `element`.

The probes: a block whose side leans by 0.05 degrees reports a perpendicularity of `height · tan(0.05°)` to `1e-9` and FAIL against a 0.01 limit when 20 high; the alignment above (`tests/cpp/test_measure.cpp`, `tests/python/test_measure.py`).
