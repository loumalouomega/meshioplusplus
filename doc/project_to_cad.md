# Back onto the CAD — `project_to_cad`

`project_to_cad` moves a mesh's points back onto the CAD model it was made from (roadmap 6.10.1). A mesh [`mesh_cad_surface`](./mesh_cad.md) wrote names, cell by cell, the CAD face it lies on (`cad:face` on its triangles and quads) and the CAD edge it runs along (`cad:edge` on its lines). An operation that knows no CAD — smoothing, an adaptation, a remesh that keeps those arrays — moves points off the surfaces; `project_to_cad` puts them back.

```python
import meshioplusplus as mp

model = mp.read("part.cad.json")
surface = mp.mesh_cad_surface(model, 0.1)
smoothed = mp.smooth(surface, iterations=10, method="laplacian")  # shrinks a curved part
back, report = mp.project_to_cad(smoothed, model, return_report=True)
print(report["num_projected"], report["max_move"], report["outside_trim"])
```

## How a point is placed

Each point is classified by the cells that use it:

| the point lies on | it goes to |
|---|---|
| the lines of one edge | the nearest point of that edge: its first use's trim curve lifted through its face's surface, a scan of 64 steps, then Newton on `(X - p) · X' = 0` |
| the lines of two or more edges | a corner of the model: the nearest point of the lowest of those edges, which is the edge's end for a point near it (counted in `num_corners`) |
| cells that all name one face | the nearest point of that face's surface |
| cells naming two faces and no edge | nowhere: it stays, counted in `num_ambiguous` |
| no named cell | nowhere: it stays |

The nearest point of a face is found by Newton's method on the first-order condition `(S - p) · S_u = (S - p) · S_v = 0`, whose Hessian is `JᵀJ + (S - p) · S₂` with the surface's second derivatives; where that Hessian is not positive definite the step is Gauss–Newton's (the first fundamental form), and it is halved until the distance does not grow, the parameters staying in the knot range. Newton converges only from inside the right basin, so it starts from the nearest three of a grid of the patch's points, `(degree + 1)` per knot span in each direction (8 to 32), the subdivision start of Ma and Hewitt (2003), and the nearest result is kept. Points are independent and projected in parallel, so the result does not depend on threads.

A projection that lands outside the face's trimmed region (the trim loops sampled `trim_samples` times per span, even–odd) is still the nearest surface point: it is moved, counted in `num_outside_trim` and listed in `outside_trim`, never clamped to the trim. A point on a trim (within a relative 1e-9 of the loops' extent) counts as inside, since a point lying on a trim is only on it to rounding. A point farther than `max_distance` from where it would go stays where it is and is counted in `num_too_far`. Data, regions and property sets pass through; only coordinates change.

## Options and report

| option | meaning |
|---|---|
| `face_array`, `edge_array` | the cell data naming each cell's face and each line's edge (`cad:face`, `cad:edge`) |
| `max_distance` | a point farther than this from where it would go stays (0, the default, is no limit) |
| `trim_samples` | samples per knot span when testing a face's trims (16) |

The report has `num_projected` (points moved onto a face or an edge), `num_on_edges` and `num_corners` of those, `num_ambiguous`, `num_too_far`, `num_failed` (a face with no usable surface or an edge with no usable trim), `num_outside_trim`, `max_move`, `outside_trim` and `warnings`. A mesh with no `face_array` cell data, points that are not three-dimensional or an index outside the model is refused by name.

## Which route keeps a mesh on the CAD

Two routes keep an adapted mesh on the model (roadmap 6.10). **Re-mesh from the model** when the change can be said as a size field or a metric: [`mesh_cad_surface`](./mesh_cad.md) with a `background` mesh or a `metric_background` keeps every point on the model by construction, with conforming edges and every face's trims respected. **Project back** when an operation that knows no CAD has already made the mesh — a smoothing, an adaptation by a backend that keeps `cad:face` and `cad:edge`, a hand edit: `project_to_cad` restores the geometry without changing the connectivity. Projection moves points but cannot fix a triangle that straddles a trim or a fold the operation introduced; `outside_trim` names the points where that happened.

## Bindings

| surface | entry |
|---|---|
| Python | `meshioplusplus.project_to_cad(mesh, model, *, face_array="cad:face", edge_array="cad:edge", max_distance=0.0, trim_samples=16, return_report=False)`; `model` is a mesh read from a CAD JSON, STEP or IGES file or its `info["brep"]` dict; no numpy twin |
| C++ | `project_to_cad(const Mesh&, const BrepInfo&, const ProjectToCadOptions&)` → `ProjectToCadResult` (`operations/project_to_cad.hpp`) |
| C | `mio_project_to_cad(mesh, info, opts, report)` with `mio_project_to_cad_opts` (64 bytes, through `mio_project_to_cad_opts_init`) and `mio_project_to_cad_report` (96 bytes); `mio_project_to_cad_outside_trim(points, n)` for the thread's last call |
| Fortran | `mio_project_to_cad(mesh, info, face_array, edge_array, max_distance, trim_samples, report, stat, errmsg)`, `mio_project_to_cad_outside_trim()` (1-based) |
| Julia | `project_to_cad(mesh, info; face_array, edge_array, max_distance, trim_samples)` |
| R | `mio_project_to_cad(mesh, info, face_array, edge_array, max_distance, trim_samples)` |
| WASM | `projectToCad(mesh, info, {faceArray, edgeArray, maxDistance, trimSamples})` |
| CLI | `project-to-cad INFILE MODEL OUTFILE [--face-array NAME] [--edge-array NAME] [--max-distance D] [--trim-samples N] [--json]`, both CLIs |
| MCP | `project_to_cad(input_path, model_path, output_path, ...)` |
| Pipeline | `ProjectToCad` with `Model` (a CAD JSON, STEP or IGES path), or the B-rep a CAD Input carries while no step has changed the mesh since; the model stays with the run after it |

## References

- L. Ma and W. T. Hewitt, "Point inversion and projection for NURBS curve and surface: control polygon approach", *Computer Aided Geometric Design* 20 (2003).
- L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed., Springer (1997), section 6.1.
