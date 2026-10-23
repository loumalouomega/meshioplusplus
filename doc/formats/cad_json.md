# Kratos CAD JSON (`.cad.json`)

The NURBS boundary representation that [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos)' isogeometric analysis (IGA) reads and writes: trimmed NURBS faces, the loops of trimming curves that bound them, edges with 3-D curves and the coupling topology between faces. Kratos's `CadJsonInput` (`kratos/input_output/cad_json_input.h`, behind the `CadIoModeler`) reads it, and `CadJsonOutput` writes it, in the schema of Teschemacher et al., [*Realization of CAD-integrated shell simulation based on isogeometric B-Rep analysis*](https://doi.org/10.1186/s40323-018-0109-4) (AMSES 2018). Kratos's `.mdpa` decks hold no NURBS (their `Begin Geometries` blocks take only Lagrange geometry names), so this file is where a Kratos IGA model's geometry lives; see [the MDPA sidecar](#the-mdpa-sidecar).

| | |
|---|---|
| **Format name** | `cad_json` |
| **Extensions** | `.cad.json` (a bare `.json` is recognised by its content) |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

```python
import meshioplusplus as mp

mesh = mp.read("geometry.cad.json")              # the faces sampled, the B-rep kept
brep = mesh.info["brep"]                          # faces, loops, edges, vertices
mp.write("copy.cad.json", mesh)                   # writes the B-rep back exactly
fine = mp.cad_json.read("geometry.cad.json", subdivisions=8)
```

## What is read

- **Faces.** Every `breps[].faces[]` entry: its `surface` (`degrees`, `knot_vectors`, `control_points` as `[id, [x, y, z, w]]` or `[[x, y, z, w]]`, `is_rational`, `is_trimmed`), its `boundary_loops` (`loop_type` `outer` or `inner`, each a list of `trimming_curves` with `trim_index`, `curve_direction` and a 2-D `parameter_curve` with its `active_range`), its `embedded_edges`, and `swapped_surface_normal`. Knot vectors may be in the full form (`n + p + 1` knots) or Kratos's reduced one (`n + p - 1`, without the outer knots), told apart from the control-point count the way Kratos's own `NurbsSurfaceGeometry` does (the reduced form first); they are held in the full form and written back in the form they came in.
- **Edges and vertices.** Every `edges[]` entry with its optional `3d_curve` and its `topology` (the faces and trim indices it is, with `relative_direction`), and every `vertices[]` entry with its topology and `local_coordinates`.
- **Identity.** `brep_id` or `brep_name` on every entity (an id of -1 in meshio++ means the entity is named instead).
- **The rest.** Members meshio++ does not model (Kratos's `tolerances` and `version_number`, a face's `embedded_loops` and `embedded_points`, anything a newer schema adds) are kept per object as compact JSON text (`extra`) and written back in place, so nothing a file carries is dropped.
- **The import dialect.** The `Surfaces`/`Lines` files Kratos's `ImportNurbsSbmModeler` reads (shifted-boundary IGA: `CPCoordinates`, `Weights`, `knotVector`/`knotVectors`, `pDegree`, `Layer`) are read too: surfaces as untrimmed faces, lines as rational edge curves (as that modeler builds them). That modeler reads no surface weights, and its files often give them as zeros: weights all zero or all one make a polynomial surface. Such a file is written back in the `breps` schema.

Parsing is strict JSON plus comments (Kratos's `Parameters` accept them). A malformed file, text that is not UTF-8, a control net that does not match its knots, a decreasing knot vector, a degree outside 1–32, or a non-positive weight on a rational entity is a `ReadError`.

## The mesh

The mesh is a sampling of the faces, so the usual operations, writers and viewers apply to it, while the exact geometry stays in the side channel.

- **Faces.** Every knot span of every face becomes one Lagrange quadrilateral of the face's higher degree: `quad`, `quad9`, or `VTK_LAGRANGE_QUADRILATERAL` from degree 3, its nodes the exact (rational) surface at the cell's lattice. `subdivisions` splits each span further. The default (`0`, automatic) keeps one cell per span, the IGA analysis elements, except on a trimmed face, which is refined to about 16 cells per direction so that its trim shows.
- **Trimming.** A trimmed face keeps the cells whose parametric centre lies inside its loops (even-odd over all of them, each trimming curve sampled at `trim_samples` points per knot span; `trim=False` keeps every span). The trimmed edge of the mesh is therefore staircase-accurate to a cell, while the trimming curves stay exact in the side channel. A cell the trim crosses has `cad:trimmed` = 1. For a mesh that follows the trims exactly and conforms across faces, triangulate the model with [`mesh_cad_surface`](/mesh_cad) (v17.13.0), which needs no CAD kernel.
- **Edges.** Every edge with a 3-D curve becomes Lagrange line cells (`line`, `line3`, `VTK_LAGRANGE_CURVE`) over its knot spans within its active range; `edges=False` skips them.
- **Cell data.** `cad:kind` (0 a face cell, 1 an edge cell), `cad:index` (into the side channel's faces or edges), `cad:id` (that entity's `brep_id`), `cad:span` (the knot span per direction, counted over the spans of nonzero length) and `cad:trimmed`.
- **Points** are shared within a face or an edge, not between them (adjacent faces of a Kratos model are coupled, not conforming).

## The side channel

`mesh.info["brep"]` (C++ `BrepInfo`, `brep.hpp`) holds the patches (`degrees`, full `knots`, the control net u-fastest, `weights`), the brep groups, the faces with their loops and trims, the edges, the vertices and the file's own extra members; see the module docstring of `meshioplusplus.cad_json._cad_json` for the dict layout. `mp.cad_json.sample(brep, ...)` samples it again, as the reader does.

`mp.iga_evaluate`, `mp.iga_refine` and `mp.iga_sample` act on this model directly (v16.43.0, [isogeometric layer](../iga.md)): a face is evaluated with derivatives and the trimming test, refined by knot insertion and degree elevation with its trimming loops untouched, and sampled again.

## What is written

- **With a side channel** (a mesh read from CAD JSON, STEP or IGES, or the result of [`fit_cad`](/fit_cad)), the B-rep is written exactly, in the layout of Kratos's own `PrettyPrintJsonString`. The keys Kratos reads are always written (`is_rational`, `swapped_surface_normal`, `boundary_loops` …), so a hand-written file gains its defaults on the first rewrite; a file meshio++ wrote reads back to the same side channel and rewrites to the same bytes. Kratos makes model-part nodes of surface and 3-D curve control points, which need ids: the writers warn when a net has none.
- **From a plain mesh**, every `quad` cell becomes an untrimmed bilinear face (degree 1, its corners the control points, ids the point ids + 1, so neighbouring faces share their nodes) and every `line` cell an edge; other cells are dropped with a warning, and a mesh with neither is a `WriteError`.

Numbers are written as integers when they are integral (below 2⁵³), else as the shortest of `%.15g`, `%.16g` and `%.17g` that reads back exactly; both engines write the same bytes.

## The MDPA sidecar

A Kratos IGA project keeps its geometry beside the deck: `mp.mdpa.read("model.mdpa", cad_json=True)` also loads it into `mesh.info["brep"]`, from `model.cad.json` beside the deck or else from the `geometry_file_name` of the `CadIoModeler` in a `ProjectParameters.json` beside it (with `.cad.json` appended when the name lacks it, as Kratos appends it); `cad_json="path"` names it. `mp.mdpa.write(..., cad_json=True)` writes it back as `<stem>.cad.json`. C++: `mdpa_cad_json_path`, and the `read_mdpa`/`write_mdpa` overloads taking a `BrepInfo` (`formats/mdpa.hpp`). See [MDPA](/formats/mdpa#the-kratos-iga-geometry-beside-a-deck).

## Other surfaces

The C API's `mio_read_with_info` returns a `"cad_json"` handle whose `mio_brep_info_*` accessors read the side channel section by section, and `mio_write_with_info` writes it back ([C API](/c_api#v9-1-0-additions)); Fortran (`info%brep_*`), Julia (`brep_info`), R (`mio_brep_info`) and WebAssembly (`readMeshSelective(path, {info: true}).info`, `CadJsonInfo`) mirror it. Both CLIs convert to and from it through the registry. STEP and IGES files read through the optional OpenCASCADE backend carry the same side channel, so `meshio++ convert part.step part.cad.json` turns a CAD model into Kratos IGA geometry.

## Verification

The fixtures in `tests/python/meshes/cad_json/` are generated by `tools/gen_cad_json_fixtures.py` from closed-form geometry (a plate with a circular hole, the Scordelis-Lo roof as an exact rational cylinder, a two-patch beam with a coupling edge, an import-dialect square, and three closed solids for [`mesh_cad_surface`](/mesh_cad): a unit cube, a capped cylinder and a sphere), not copied from Kratos. `tests/python/test_cad_json.py` checks the two engines against each other (mesh, side channel and bytes), the exactness of the rational geometry (radii to 1e-14), the trimming and the error paths; when `KratosMultiphysics` and its `IgaApplication` are installed (`pip install KratosMultiphysics KratosIgaApplication`), it loads every file meshio++ writes through Kratos's own `CadJsonInput` and compares `GlobalCoordinates` with meshio++'s evaluation to 1e-12. The cases of the Kratos [Examples](https://github.com/KratosMultiphysics/Examples/tree/master/iga) `iga/` folder were read and rewritten during development: each rewrite equals the original once the canonical defaults are filled in. `tests/cpp/test_cad_json.cpp` covers the native path. Notebook [35](https://github.com/loumalouomega/meshioplusplus/blob/main/example/python/35_kratos_iga_cad_json.ipynb) renders the fixtures, the MDPA sidecar and, with OpenCASCADE, a `fit_cad` solid with a freeform face written as CAD JSON.
