# DXF (`.dxf`)

[Drawing Exchange Format](https://help.autodesk.com/view/OARX/2024/ENU/?guid=GUID-235B22E0-A567-4CF6-92D3-38A2306D73F3) output for sketches and outlines: the format every CAD sketcher, laser cutter and CNC post-processor imports. A section of a scan chained by [`slice(chain=True)`](../slice.md#chained-sections) or `sections`, then fitted by [`fit_curves`](../fit_curves.md), crosses to CAD and CAM at 1:1 as a DXF drawing of lines, arcs and circles (roadmap 7.6.4). It is a drawing, not a mesh format: only line meshes and sketches are written, and only the curves come back.

| | |
|---|---|
| **Format name** | `dxf` |
| **Extensions** | `.dxf` |
| **Read / Write** | ✓ (a drawing: `LINE`, `ARC`, `CIRCLE`, `LWPOLYLINE`, `SPLINE`) / ✓ |
| **Extra dependencies** | — |

## Writing

`meshioplusplus.dxf.write(path, mesh)` takes no options.

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

profile = meshioplusplus.sections(scan, direction=(0, 0, 1), offsets=[2.0])
sketch = meshioplusplus.fit_curves(profile, tolerance=0.01)
meshioplusplus.write("profile.dxf", sketch)
```

The file is ASCII DXF R2000 (`$ACADVER AC1015`), the oldest release that has `LWPOLYLINE` and `SPLINE`, with what an R2000 reader expects around the entities: a HEADER (`$ACADVER`, `$HANDSEED`, `$INSUNITS`), the nine symbol tables (one layer, `0`; the `ByBlock`, `ByLayer` and `Continuous` line types; the `Standard` text and dimension styles; the `ACAD` application), the model and paper space block records and blocks, and an OBJECTS section with the root dictionary and the two layouts. Every object carries a handle and its owner. ezdxf's audit of each file the tests write finds no error and nothing to fix.

- **A sketch** (a mesh carrying `fit_curves`' `sketch:type`, `sketch:params` and `sketch:frame` field data) is written as its pieces. A line becomes a `LINE` between its ends in world coordinates. An arc or a full circle becomes an `ARC` or `CIRCLE` in the object coordinate system (OCS) of the sketch plane's normal, built by DXF's arbitrary-axis algorithm, with that normal as its extrusion direction (`210`/`220`/`230`); a clockwise arc is written with its start and end angles swapped, since DXF arcs run counterclockwise about the extrusion.
- **Any other line mesh** is written as one `LINE` per segment. Vertex cells are skipped. Surface and volume cells are refused by name (slice or section a surface first).
- **Units:** `$INSUNITS` is the code of `field_data["meshio:length_unit"]` (inches 1, feet 2, miles 3, millimetres 4, centimetres 5, metres 6, kilometres 7, yards 10, nanometres 12, micrometres 13, decimetres 14), or 0, unitless, for no unit or any other. Coordinates are written as they are, so a sketch in millimetres is 1:1 in a millimetre drawing.
- **Numbers** are written with 17 significant digits (`-0` as `0`), so every coordinate reads back to the same double.
- **Provenance:** the [provenance block](../provenance.md) leads the file as `999` comment pairs, which DXF readers skip.

`SPLINE` and `LWPOLYLINE` are read but not written: a `fit_curves` sketch has no B-spline pieces, and a polyline is written as its `LINE`s.

## Reading

`meshioplusplus.read("part.dxf")` returns the drawing's curves as a line mesh, in world coordinates:

| Entity | Read as |
|---|---|
| `LINE` | one segment |
| `ARC` | segments, 64 per full turn, through the OCS of its extrusion direction |
| `CIRCLE` | a closed chain of 64 segments |
| `LWPOLYLINE` | its vertices at the polyline's elevation, in its OCS; a bulge span sampled as an arc (64 segments per full turn); `70` bit 1 closes it |
| `SPLINE` | 8 samples per control point along the B-spline (de Boor, with weights when given) |

Points are shared where entities meet at exactly the same coordinates. The cell data `dxf:entity` holds each segment's entity index and `dxf:type` its kind (1 line, 2 arc, 3 circle, 4 polyline, 5 spline); `$INSUNITS` becomes `meshio:length_unit`. Any other entity in the `ENTITIES` section (`TEXT`, `INSERT`, `HATCH`, `POLYLINE`, `3DFACE`, …) is refused by name, as are inconsistent `SPLINE` knots and a non-positive radius: the `svg` and `tikz` rule, a reader limited to what the writer emits.

## Engines

The C++ core (`formats/dxf.cpp`) and the Python reference (`dxf/_dxf.py`) write the same bytes and read the same mesh, pinned by `tests/python/test_dxf.py`. The native reader returns a mesh or throws `ReadError` on any input, every number going through `detail::parse_double`; the fuzz seed corpus carries a line drawing.
