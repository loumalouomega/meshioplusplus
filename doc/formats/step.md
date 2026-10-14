# STEP (`.step`, `.stp`)

[ISO 10303-21](https://en.wikipedia.org/wiki/ISO_10303-21), read through the optional OpenCASCADE backend and written without it. The writer is a hand-off of a closed surface to CAD: a faceted solid in the AP214 (`AUTOMOTIVE_DESIGN`) schema that FreeCAD, OpenCASCADE-based tools and commercial CAD systems open as a solid body, made of one planar face per facet. It is the CAD end of the [skin-to-CAD workflow](/cad_fit): the baseline a CAD user can measure, section and boolean against, where an STL import gives only an unconnected triangle soup.

| | |
|---|---|
| **Format name** | `step` |
| **Extensions** | `.step`, `.stp` |
| **Read / Write** | ✓ (with OpenCASCADE) / ✓ |
| **Extra dependencies** | reading: OpenCASCADE 7.6 or later (`-DMESHIOPLUSPLUS_WITH_OCCT=ON`) |

## Reading

Reading STEP means evaluating its curves and surfaces, which is a CAD kernel's job: the reader runs on the optional OpenCASCADE backend that [`fit_cad`](/fit_cad) introduced (`-DMESHIOPLUSPLUS_WITH_OCCT=ON`, LGPL-2.1 with the OCCT exception, bring your own, absent from the wheels, the WASM build and the release CLI binaries; see the [OpenCASCADE decision](/spikes#the-opencascade-dependency)). Without it, reading raises `NotImplementedError` (C++ `Unsupported`) naming the flag, and the round trip on the [conformance page](/conformance#step) is checked only in a build with it. [IGES](/formats/iges) is read the same way.

```python
import meshioplusplus as mp

mesh = mp.read("part.step")                       # deflection 1e-3 of the bounding-box diagonal
mesh = mp.step.read("part.step", deflection=0.05, angular_deflection=0.3)
brep = mesh.info["brep"]                          # the faces, loops and edges as NURBS
mp.write("part.cad.json", mesh)                   # ... as Kratos CAD JSON
```

- **The shape.** `STEPControl_Reader` transfers every root into one shape, which `BRepBuilderAPI_NurbsConvert` turns into NURBS throughout. Converting each face's surface on its own would reparametrise a cylinder or a sphere and leave its trimming curves in the old parameters; converting the whole shape keeps every parameter curve on its converted surface.
- **The mesh.** `BRepMesh_IncrementalMesh` triangulates it at `deflection` (a length; at most 0 takes 1e-3 of the bounding-box diagonal) and `angular_deflection` (radians, 0.5 by default) into one `triangle` block, points shared between faces welded where they coincide. Cell data: `cad:kind` (0), `cad:index` and `cad:id` (the face), `cad:span` (the knot span of the triangle's centre in the face's parameters) and `cad:trimmed` (0), as for [CAD JSON](/formats/cad_json#the-mesh).
- **The side channel.** `mesh.info["brep"]` (C++ `BrepInfo`) holds every face's surface (degrees, knots, control net, weights; reversed faces flagged `swapped_normal`), its trimming loops (the outer one first, each trimming curve with its active range and direction, in the surface's own orientation) and every edge's 3-D curve and topology (the faces and trims it bounds). Faces have ids `1 .. F` and edges `F + 1 .. F + E`, control points ids from 1, so the B-rep writes as [Kratos CAD JSON](/formats/cad_json) that Kratos's `CadIoModeler` loads: a STEP model becomes Kratos IGA geometry with `mp.write(...)` or `meshio++ convert part.step part.cad.json`. Vertices are not carried.
- **Errors.** A file OpenCASCADE cannot read, or one that transfers no shape, is a `ReadError`; a non-finite deflection or a non-positive angular deflection a `ValueError`.

C++: `read_step(path)` and `read_step(path, info, options)` (`formats/step.hpp`, `CadReadOptions` in `brep.hpp`, `cad_read_has_occt()`). The C API's `mio_read_with_info` returns a `"step"` handle with the `mio_brep_info_*` accessors, which `mio_write_with_info` writes as CAD JSON; Fortran, Julia, R and the native CLI follow.

## Writing

Writing takes no options:

```python
import meshioplusplus as mp

mp.write("part.step", mp.read("part.stl"))       # a closed triangle skin
mp.write("block.stp", mp.read("block.msh"))       # a volume mesh: its skin is written
```

## What is written

- **The surface.** `triangle`, `quad` and `polygon` blocks (quadratic surface cells contribute their corners); when the mesh has volume cells, the skin of those, as the [`stl`](/formats/stl) writer takes it. Lines and vertices are dropped with a warning.
- **One face per facet.** Each facet is an `ADVANCED_FACE` on a `PLANE`, bounded by an `EDGE_LOOP` of `ORIENTED_EDGE`s over straight `EDGE_CURVE`s (`LINE`s) between `VERTEX_POINT`s. An edge two facets share is one `EDGE_CURVE`, used once in each direction, so the shell is topologically closed and needs no sewing. A facet of four or more corners that is not planar to `1e-9` of its size is fanned into triangles from its first corner, so every face lies exactly on its plane. Coplanar facets are not merged into one face.
- **One solid per piece.** Every edge-connected piece of the surface is a `CLOSED_SHELL` inside a `MANIFOLD_SOLID_BREP`; they are collected in one `ADVANCED_BREP_SHAPE_REPRESENTATION` under the usual product structure (`PRODUCT`, `PRODUCT_DEFINITION`, `SHAPE_DEFINITION_REPRESENTATION`). A piece wound inwards (negative signed volume) is written reversed, with a warning, so every solid has positive volume. A surface nested inside another (a void) becomes a second solid, not a cavity.
- **Units.** Millimetres, with a length uncertainty of `1e-7`. Coordinates are written as they are: scale with [`transform`](/transform) first if the mesh is in metres.

The roadmap item that queued this writer described the faces as `ADVANCED_FACE`s with `POLY_LOOP` bounds. That combination is not valid AP214, since an `ADVANCED_FACE` must be bounded by `EDGE_LOOP`s; the writer emits the edge-loop form, which is also what lets the faces share edges.

## Refusals

The writer refuses, with a `WriteError`, anything that is not the boundary of a solid:

- an **open**, **non-manifold** or **inconsistently wound** edge — the message gives the counts and names the first five as point-id pairs; run [`repair`](/repair) first;
- a **degenerate** facet (zero area), or an edge between two coincident points (run [`clean`](/clean) first);
- a non-finite coordinate, or a mesh with no surface facets at all.

## C++ vs Python

The C++ core and the Python reference write identical bytes: entity ids are sequential, every number is written with `%.17g` (with the decimal point STEP requires, `3.` and `1.E-05`), the header time stamp is fixed at `1970-01-01T00:00:00`, and the [provenance](/provenance) lines go in a `/* ... */` comment after `ISO-10303-21;`. `FILE_NAME` and `PRODUCT` carry the file's basename.

## Validation

Checked with OpenCASCADE 8.0.1 (the `cadquery-ocp` 8.0.1.1 wheels): `STEPControl_Reader` reads every file `tests/python/test_step.py` writes — a subdivided cube, an inward-wound cube, two disjoint cubes, a tetrahedron's skin, a box with one warped (fanned) quad, a polygon-block cube and a torus — as the expected number of solids, `BRepCheck_Analyzer` finds each valid, and `GProp` gives the mesh's own volume to `1e-9`. The test runs when `OCP` is importable (`pip install cadquery-ocp`) and skips otherwise. FreeCAD has not yet been checked directly; it reads STEP through the same OpenCASCADE reader. The reader is checked on the backend's own build (OpenCASCADE 7.6, the CI `occt` leg): a written box reads back as its faces, each trimmed by one loop and every edge bounding two of them, and its CAD JSON rewrite samples back onto the box (`test_step.py`, `tests/cpp/test_cad_read.cpp`).

## Not supported

- Reading without the OpenCASCADE backend (see above), and a reader's vertices (STEP's `VERTEX_POINT`s are carried as edge ends only).
- Analytic faces — this writer is faceted by design. [`fit_cad`](/fit_cad) writes the fitted cylinders, cones, spheres and tori, and B-spline faces for what no primitive explains, as CAD faces through the optional OpenCASCADE backend.
- Colours, layers, assemblies, and AP242 product manufacturing information.
