# STEP (`.step`, `.stp`)

[ISO 10303-21](https://en.wikipedia.org/wiki/ISO_10303-21), read through OpenCASCADE (included in the PyPI wheels and the release binaries) and written without it. The writer is a hand-off of a closed surface to CAD: a B-rep solid in the AP214 (`AUTOMOTIVE_DESIGN`) schema that FreeCAD, OpenCASCADE-based tools and commercial CAD systems open as a solid body. By default each flat region is one planar face with its holes as inner loops, and, on a skin [`segment`](/segment) and [`fit_primitives`](/fit_primitives) annotated, each fitted cylinder, cone, sphere and torus is one analytic face — the clean body a CAD user can select, measure, section and boolean against, where an STL import gives only an unconnected triangle soup. No CAD kernel is involved. It is the CAD end of the [skin-to-CAD workflow](/cad_fit).

| | |
|---|---|
| **Format name** | `step` |
| **Extensions** | `.step`, `.stp` |
| **Read / Write** | ✓ (with OpenCASCADE) / ✓ |
| **Extra dependencies** | reading: OpenCASCADE 7.6 or later (in the wheels and release binaries; on by default from source) |

## Reading

Reading STEP means evaluating its curves and surfaces, which is a CAD kernel's job: the reader runs on the OpenCASCADE backend that [`fit_cad`](/fit_cad) introduced (LGPL-2.1 with the OCCT exception, included in the PyPI wheels and the release CLI binaries, on by default and auto-detected in a build from source, absent from the lean WASM build; see [getting OpenCASCADE](/fit_cad#getting-opencascade) and the [OpenCASCADE decision](/spikes#the-opencascade-dependency)). The [full WASM artifact](/wasm#lean-and-full-builds) (v17.2.0) links it, so STEP and IGES read in the browser and in Node there, `readMeshSelective(path, {info: true})` returning the B-rep as CAD JSON's shape. Without it, reading raises `NotImplementedError` (C++ `Unsupported`) naming the flag, and the round trip on the [conformance page](/conformance#step) is checked only in a build with it. [IGES](/formats/iges) is read the same way.

```python
import meshioplusplus as mp

mesh = mp.read("part.step")                       # deflection 1e-3 of the bounding-box diagonal
mesh = mp.step.read("part.step", deflection=0.05, angular_deflection=0.3)
brep = mesh.info["brep"]                          # the faces, loops and edges as NURBS
mp.write("part.cad.json", mesh)                   # ... as Kratos CAD JSON
```

- **The shape.** `STEPControl_Reader` transfers every root into one shape, which `BRepBuilderAPI_NurbsConvert` turns into NURBS throughout. Converting each face's surface on its own would reparametrise a cylinder or a sphere and leave its trimming curves in the old parameters; converting the whole shape keeps every parameter curve on its converted surface. Where that conversion fails, or leaves a face without a B-spline surface or with trims outside its parameter range (OpenCASCADE 7.9 does all three on some shapes), the shape is kept as it is and each surface and curve that is not a B-spline is approximated by one in its own parameters, within `1e-6` of the bounding-box diagonal; planes and lines stay exact.
- **The mesh.** `BRepMesh_IncrementalMesh` triangulates it at `deflection` (a length; at most 0 takes 1e-3 of the bounding-box diagonal) and `angular_deflection` (radians, 0.5 by default) into one `triangle` block, points shared between faces welded where they coincide. Cell data: `cad:kind` (0), `cad:index` and `cad:id` (the face), `cad:span` (the knot span of the triangle's centre in the face's parameters) and `cad:trimmed` (0), as for [CAD JSON](/formats/cad_json#the-mesh).
- **The side channel.** `mesh.info["brep"]` (C++ `BrepInfo`) holds every face's surface (degrees, knots, control net, weights; reversed faces flagged `swapped_normal`), its trimming loops (the outer one first, each trimming curve with its active range and direction, in the surface's own orientation) and every edge's 3-D curve and topology (the faces and trims it bounds). Faces have ids `1 .. F` and edges `F + 1 .. F + E`, control points ids from 1, so the B-rep writes as [Kratos CAD JSON](/formats/cad_json) that Kratos's `CadIoModeler` loads: a STEP model becomes Kratos IGA geometry with `mp.write(...)` or `meshio++ convert part.step part.cad.json`. Vertices are not carried. Since v17.20.0 each edge's `extra` holds OpenCASCADE's tolerances, `{"tolerance": t, "end_tolerances": [t0, t1]}` (the edge's own and its two vertices', the approximation tolerance added where the edge was approximated), which a CAD JSON write keeps and [`mesh_cad_surface(model_tolerances=True)`](/mesh_cad) accepts as the gaps the model allows there.
- **Errors.** A file OpenCASCADE cannot read, or one that transfers no shape, is a `ReadError`; a non-finite deflection or a non-positive angular deflection a `ValueError`.

C++: `read_step(path)` and `read_step(path, info, options)` (`formats/step.hpp`, `CadReadOptions` in `brep.hpp`, `cad_read_has_occt()`). The C API's `mio_read_with_info` returns a `"step"` handle with the `mio_brep_info_*` accessors, which `mio_write_with_info` writes as CAD JSON; Fortran, Julia, R and the native CLI follow.

## Writing

```python
import numpy as np
import meshioplusplus as mp

mp.write("part.step", mp.read("part.stl"))       # one planar face per flat region
mp.write("block.stp", mp.read("block.msh"))       # a volume mesh: its skin is written

skin = mp.read("bracket.stl")
diag = np.linalg.norm(np.ptp(skin.points, axis=0))
seg = mp.segment(skin, method="primitive")        # tolerance 1e-3 of the diagonal
fit = mp.fit_primitives(seg, tolerance=1e-3 * diag, relative=False)
mp.write("bracket.step", fit)                     # planes, cylinders, cones ... as analytic faces
mp.write("bracket_facets.step", skin, faces="facet")
```

`faces` picks how facets become faces (C++ `StepWriteOptions::mFaces`, `StepFaces`):

| `faces` | Faces written |
| --- | --- |
| `"auto"` (default) | `"fitted"` when the mesh carries `segment:patch` and the `fit:patch_*` table, else `"planar"`. |
| `"planar"` | Each edge-connected group of coplanar facets (to `1e-6` of the bounding-box diagonal) is one `PLANE` face; on a face with several boundary loops the one enclosing the largest area is the `FACE_OUTER_BOUND` and the others are holes. |
| `"fitted"` | Each patch `fit_primitives` accepted becomes one face on its `PLANE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`, `SPHERICAL_SURFACE` or `TOROIDAL_SURFACE`, with the parameters of the fit table; the rest of the skin as with `"planar"`. |
| `"facet"` | One `PLANE` face per facet, the layout of the writer before v17.9.0 (byte-identical to it). |

`tolerance` (`"fitted"` only) is how far a patch's points may lie from its surface, as a length; at most 0, the default, takes three times the largest `fit:patch_tolerance` of an accepted patch, or `1e-3` of the bounding-box diagonal when there is none.

## What is written

- **The surface.** `triangle`, `quad` and `polygon` blocks (quadratic surface cells contribute their corners); when the mesh has volume cells, the skin of those, as the [`stl`](/formats/stl) writer takes it. Lines and vertices are dropped with a warning.
- **Faces.** Each face is an `ADVANCED_FACE` bounded by `EDGE_LOOP`s of `ORIENTED_EDGE`s between `VERTEX_POINT`s. An edge two faces share is one `EDGE_CURVE`, used once in each direction, so the shell is topologically closed and needs no sewing. With `"facet"`, a facet of four or more corners that is not planar to `1e-9` of its size is fanned into triangles from its first corner, so every face lies exactly on its plane; a merged planar face takes the plane of its first facet. An analytic face's `same_sense` follows the facets' normals.
- **Edges.** The mesh edges between two faces form chains from one corner (a point where three or more faces meet) to the next; a closed chain is split at its point farthest from its start. Between two planar faces a chain becomes the straight `LINE`s of its Douglas-Peucker corners, so the collinear points of the mesh disappear. Beside a curved analytic face a chain is a `CIRCLE` when it lies within the tolerance of one about the surface's axis (the rim of a bore, a cone's base), a `LINE` when it is straight (a cylinder's ruling), and otherwise a degree-1 `B_SPLINE_CURVE_WITH_KNOTS` through its points. A mesh point on analytic faces only is first moved onto them (alternating projections), when that moves it less than the tolerance.
- **Periodic faces.** A cylinder, cone, sphere or torus face has one bound per boundary loop and no seam edge, as ISO 10303-42 allows; a reader adds its seams. A whole sphere or torus, with no boundary at all, is split in two faces by a plane through its centre.
- **Demotion.** An accepted patch with a point further than the tolerance from its surface, or whose boundary does not close into loops, is written as planar faces instead and counted; nothing is silently wrong.
- **One solid per piece.** Every edge-connected piece of the surface is a `CLOSED_SHELL` inside a `MANIFOLD_SOLID_BREP`; they are collected in one `ADVANCED_BREP_SHAPE_REPRESENTATION` under the usual product structure (`PRODUCT`, `PRODUCT_DEFINITION`, `SHAPE_DEFINITION_REPRESENTATION`). A piece wound inwards (negative signed volume) is written reversed, with a warning, so every solid has positive volume. A surface nested inside another (a void) becomes a second solid, not a cavity.
- **Units.** Millimetres, with a length uncertainty of `1e-7`. Coordinates are written as they are: scale with [`transform`](/transform) first if the mesh is in metres.

The roadmap item that queued this writer described the faces as `ADVANCED_FACE`s with `POLY_LOOP` bounds. That combination is not valid AP214, since an `ADVANCED_FACE` must be bounded by `EDGE_LOOP`s; the writer emits the edge-loop form, which is also what lets the faces share edges.

## Refusals

The writer refuses, with a `WriteError`, anything that is not the boundary of a solid:

- an **open**, **non-manifold** or **inconsistently wound** edge — the message gives the counts and names the first five as point-id pairs; run [`repair`](/repair) first;
- a **degenerate** facet (zero area), or an edge between two coincident points (run [`clean`](/clean) first);
- a non-finite coordinate, or a mesh with no surface facets at all.

## C++ vs Python

The Python reference writes `faces="facet"` only: without the core, `"auto"` falls back to it with a warning and `"planar"` or `"fitted"` raise `NotImplementedError`. In that mode the C++ core and the Python reference write identical bytes: entity ids are sequential, every number is written with `%.17g` (with the decimal point STEP requires, `3.` and `1.E-05`), the header time stamp is fixed at `1970-01-01T00:00:00`, and the [provenance](/provenance) lines go in a `/* ... */` comment after `ISO-10303-21;`. `FILE_NAME` and `PRODUCT` carry the file's basename.

## Validation

Checked with OpenCASCADE 8.0.1 (the `cadquery-ocp` 8.0.1.1 wheels): `STEPControl_Reader` reads every file `tests/python/test_step.py` writes — a subdivided cube, an inward-wound cube, two disjoint cubes, a tetrahedron's skin, a box with one warped (fanned) quad, a polygon-block cube and a torus — as the expected number of solids, with `faces="facet"` and `"planar"` alike, `BRepCheck_Analyzer` finds each valid, and `GProp` gives the mesh's own volume to `1e-9`. With `"fitted"`, a CAD-style bracket reads back as one valid solid of twelve planes and four cylinders with its exact volume, and a filleted plate as eight faces. The test runs when `OCP` is importable (`pip install cadquery-ocp`) and skips otherwise. FreeCAD has not yet been checked directly; it reads STEP through the same OpenCASCADE reader. The reader is checked on the backend's own build (OpenCASCADE 7.6, the CI `occt` leg): a written box reads back as its faces, each trimmed by one loop and every edge bounding two of them, and its CAD JSON rewrite samples back onto the box (`test_step.py`, `tests/cpp/test_cad_read.cpp`).

## Not supported

- Reading without the OpenCASCADE backend (see above), and a reader's vertices (STEP's `VERTEX_POINT`s are carried as edge ends only).
- Trimmed freeform, fillet and thread faces: what no plane, sphere, cylinder, cone or torus explains stays planar faces here. [`fit_cad`](/fit_cad) builds B-spline, fillet and thread faces through OpenCASCADE, and its kernel-free `backend="native"` is this writer.
- Colours, layers, assemblies, and AP242 product manufacturing information.

The NURBS B-rep a STEP file reads into `mesh.info["brep"]` can be refined, evaluated and sampled exactly with the [isogeometric layer](../iga.md) (`iga_refine`, `iga_evaluate`, `iga_sample`); STEP files themselves are not written from a refined B-rep, but `mp.write("refined.cad.json", refined)` is. [`mesh_cad_surface`](../mesh_cad.md) triangulates the same B-rep at a target edge size, conforming across faces, in place of `BRepMesh`'s deflection-driven mesh.
