# Mesh generation from CAD — `mesh_cad_surface`

`mesh_cad_surface` triangulates the faces of a CAD model into one conforming surface mesh, with edges of about one target size `h`, and needs no CAD kernel. The model is the B-rep a [Kratos CAD JSON](/formats/cad_json) read keeps in `mesh.info["brep"]` (the C++ `BrepInfo`), or a STEP or IGES read through OpenCASCADE. Every B-rep edge is sampled once and shared by the faces that use it, so the faces meet exactly and a closed solid comes out watertight. Each face is a constrained Delaunay triangulation of its parameter domain, lifted onto its surface. It is the first milestone of mesh generation in meshio++ ([roadmap §6](/roadmap#_6-mesh-generation)): one global size, no quality guarantee yet.

```python
import meshioplusplus as mp

model = mp.read("part.cad.json")  # or part.step / part.iges with OpenCASCADE
surface, report = mp.mesh_cad_surface(model, 0.1, return_report=True)
mp.write("part.vtu", surface)
report["num_failed"], report["num_open_edges"], report["min_angle_deg"]
```

```sh
meshioplusplus mesh-cad-surface part.cad.json part.vtu --size 0.1 --deflection 0.001 --json
```

![A unit cube, a capped cylinder and a sphere from CAD JSON, meshed watertight at h = 0.2 and coloured by source face](/images/mesh_cad_solids.png)

The planar triangulator underneath is public too: `triangulate_polygon_2d(outer, holes=(), target_size=0)` triangulates a polygon with holes in the plane z = 0.

## How it works

1. **Topology.** Every edge's uses (`BrepEdge.topology`: a face id or name, a trim index and a relative direction) are resolved to the face trims they name. A trim no edge names becomes its own edge; with `match_unreferenced_trims` (the default), unreferenced trims whose 3-D images coincide (both ends and the middle within `tolerance`) are paired, so models without edge topology, such as the untrimmed patches of the `Surfaces` dialect, still conform. An untrimmed face is its patch's parameter rectangle. A loop's trims are walked head to tail by their ends in the face's parameters, whatever order they are listed in.
2. **Corners from topology.** A trim's end and the next trim's start are one corner, and so are the ends of every use of one edge, joined by a union-find rather than by comparing coordinates. Two faces whose trims disagree by the model's tolerance still share their corner points exactly; the disagreement is measured and reported as `max_edge_gap`.
3. **Edges, sampled once.** Each edge is sampled on its first use, by arc length (Gauss–Legendre on the trim's knot spans), into `ceil(L / h)` equal pieces between the trim's corners (a polyline trim keeps its vertices), at least `min_points_per_edge` points inside. With `max_chordal_deviation`, a piece is halved while its sagitta, the distance of the edge from its chord, exceeds the cap (OpenCASCADE's "linear deflection"). Every other use, the neighbouring face or the same face across a seam, projects the same 3-D points onto its own trim, and the points' ids are shared: this is what makes the faces conform. An edge shorter than the tolerance, a pole, keeps its two ends, which are one corner.
4. **Faces, in their parameters.** Each face's loops become closed polylines in its parameters, mapped by one affine transformation (the Cholesky factor of the face's mean first fundamental form) so that a Delaunay triangulation there follows the surface's own metric on average rather than raw parameters, whose scales can differ by orders of magnitude (the Scordelis–Lo roof's differ fortyfold). The constrained Delaunay triangulation inserts the boundary points with Lawson's flips, recovers every boundary segment by re-triangulating the cavity it crosses (Anglada 1997), and keeps the triangles inside by the even-odd rule, so holes need no orientation. Interior points come from an equilateral lattice whose step is at most `h` on the surface, kept only farther than `h / 2` from the boundary points and the points already kept and `0.4 h` from the boundary segments. The points are lifted onto the surface, triangles a pole collapses are dropped, and the winding follows `dS/du × dS/dv`, reversed for a face whose normal is swapped.
5. **Exact predicates.** Every orientation and in-circle decision is Shewchuk's adaptive exact predicate, re-implemented from the paper. No coordinate is compared to another with a tolerance in the triangulation.

![The trimmed plate: the CAD JSON read's lattice clipped at cell centres beside mesh_cad_surface's triangulation following the hole](/images/mesh_cad_plate.png)

Notebook [42](https://github.com/loumalouomega/meshioplusplus/blob/main/example/python/42_mesh_cad_surface.ipynb) renders these meshes and the Scordelis–Lo roof, and checks their convergence.

## Options

| Python | C++ (`MeshCadSurfaceOptions`) | Meaning |
|---|---|---|
| `target_size` | `mTargetSize` | The target edge length `h`; required, > 0 |
| `max_chordal_deviation=0` | `mMaxChordalDeviation` | The sagitta cap of the boundary segments; 0 is off |
| `min_points_per_edge=0` | `mMinPointsPerEdge` | At least this many points inside every edge |
| `tolerance=0` | `mGeometricTolerance` | Edges shorter than this are degenerate, and unreferenced trims closer than this are paired; 0 is `1e-7` times the diagonal of the control nets' box |
| `record_face_ids=True` | `mRecordFaceIds` | `cad:face` cell data and a `face:<k>` cell region per face (failed faces keep an empty one) |
| `edges=True` | `mEdges` | A `line` block of the discretised edges with `cad:edge`, the `BrepEdge` index (-1 for an edge no `BrepEdge` names) |
| `match_unreferenced_trims=True` | `mMatchUnreferencedTrims` | Pair unreferenced trims by their images |
| `strict=False` | `mStrict` | Raise (`RuntimeError`; C++ `std::runtime_error`) when any face fails, instead of reporting it |
| `max_points=50000000` | `mMaxPoints` | The most points in the result; a face or edge that would exceed it fails |

## Output and report

The mesh has the points and one `triangle` block (and the `line` block of the edges). It writes to VTU, Gmsh, STL and every other surface format unchanged. Gmsh gets the faces as the `cad:face` data: with every face in one triangle block, its writer cannot give each face its own physical group.

The report (`return_report=True`; `MeshCadSurfaceResult` in C++) has the counts (`num_points`, `num_triangles`, `num_edge_segments`, `num_faces`, `num_failed`, `num_edges`), the conformity of the result (`num_open_edges`, triangle edges used once: zero on a closed solid, the model's own boundary on an open one; `num_non_manifold_edges`; `num_misoriented_edges`, used twice the same way round), `num_degenerate_dropped`, the extreme angles and edge lengths, `tolerance`, `max_edge_gap`, one record per face (`status`, `reason`, counts, extremes) and the warnings.

**A face that cannot be meshed is left out and reported, never meshed wrong.** Its status names why: `unsupported` (an invalid surface or trim), `open_loop` (its loops do not close, or a corner the topology joins it to lies elsewhere), `self_intersecting` (its boundary crosses itself in its parameters at this size: try a smaller one), `degenerate`, `projection_failed` (an edge's uses lie farther apart than `h / 4`, such as a closed edge whose uses start at different points), `budget_exceeded`, `empty`. The face's edges stay in the `line` block, and `num_open_edges` shows the hole.

## Limitations

This first milestone deliberately leaves out:

- **Quality guarantees.** The lattice gives near-equilateral triangles inside and edges within about `[0.5 h, 1.5 h]` on planar faces, but there is no Delaunay refinement (Ruppert, Chew): near sharp boundary corners and poles angles can be small, and the report says how small.
- **Sizing.** One global size, plus the chordal cap on the boundary: no size field, no curvature-adaptive sizing inside a face.
- **Metrics.** One affine map per face: a face whose parameterisation stretches very differently from place to place is meshed by its average.
- **Poles.** A pole is a fan apex, ungraded.
- **Other topologies.** Embedded edges are not meshed (a warning names them), and a closed edge whose uses start at different points is reported rather than re-aligned.
- **Other elements.** No quadrilaterals, no volume meshes, no boundary layers.

[Roadmap §6](/roadmap#_6-mesh-generation) lists what comes next. Everything is serial within a face; faces and edges run in parallel and the output is byte-identical across backends and thread counts.

## Bindings

| Surface | Entry points |
|---|---|
| Python | `meshioplusplus.mesh_cad_surface(mesh_or_brep, target_size, …, return_report=False)`, `triangulate_polygon_2d(outer, holes=(), target_size=0)` (core only: without the compiled core both raise `NotImplementedError`) |
| C++ | `operations/mesh_cad.hpp`: `mesh_cad_surface(const BrepInfo&, const MeshCadSurfaceOptions&)`, `triangulate_polygon_2d` |
| C | `mio_mesh_cad_surface` on a [`mio_format_info`](/c_api), `mio_mesh_cad_surface_opts` / `_report`, per-face `mio_mesh_cad_surface_face_status` / `_reason`, `mio_mesh_cad_surface_warnings` |
| Fortran | `info%mesh_cad_surface(target_size, …)`, `mio_mesh_cad_surface_report`, `mio_mesh_cad_face_status` / `_reason` |
| Julia | `mesh_cad_surface(info; target_size, …)` |
| R | `mio_mesh_cad_surface(info, target_size, …)` |
| WebAssembly | `meshCadSurface(info, {targetSize, …})` on `readMeshSelective(path, {info: true}).info`, in both artifacts |
| CLIs | `mesh-cad-surface IN OUT --size H [--deflection D] [--min-points-per-edge N] [--tolerance T] [--no-face-ids] [--no-edges] [--no-match] [--strict] [--json]` in both CLIs |
| MCP | `mesh_cad_surface` |

The pipeline carries no CAD model, so there is no pipeline step; `MeshCadSurface` is refused by name.

## Verification

The predicates are checked against 128-bit integer determinants, including collinear, cocircular, coplanar and cospherical configurations and Kettner et al.'s classroom grid; the triangulation on squares, L-shapes, holes, collinear and near-duplicate points, slivers and random polygons for exact area, segment recovery, Euler characteristic and the Delaunay property of every unconstrained edge; the edge sampler on rational arcs for arc length to `1e-12` and on reparametrised, reversed uses. The meshes of boxes, a capped cylinder through a seam, a sphere through its poles, an untrimmed rational patch, the trimmed plate and the CAD JSON fixtures are checked for watertightness, area and enclosed volume converging at second order, every point on the exact surface, the Hausdorff distance to a dense reference, determinism over thread counts and VTU and Gmsh round trips; with OpenCASCADE, a STEP solid against BRepMesh's own triangulation (`tests/cpp/test_predicates.cpp`, `test_cdt2d.cpp`, `test_edge_sampler.cpp`, `test_mesh_cad.cpp`, `tests/python/test_mesh_cad.py`).

## References

- J. R. Shewchuk, *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates*, Discrete & Computational Geometry 18 (1997).
- J. R. Shewchuk, *Delaunay Refinement Algorithms for Triangular Mesh Generation*, Computational Geometry 22 (2002).
- C. L. Lawson, *Software for C1 surface interpolation*, Mathematical Software III (1977).
- M. V. Anglada, *An improved incremental algorithm for constructing restricted Delaunay triangulations*, Computers & Graphics 21 (1997).
- J. Ruppert, *A Delaunay refinement algorithm for quality 2-dimensional mesh generation*, Journal of Algorithms 18 (1995).
- L. P. Chew, *Guaranteed-quality triangular meshes*, Cornell TR 89-983 (1989).
- C. Geuzaine and J.-F. Remacle, *Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities*, IJNME 79 (2009).
- J. Schöberl, *NETGEN: an advancing front 2D/3D-mesh generator based on abstract rules*, Computing and Visualization in Science 1 (1997).
- H. Si, *TetGen, a Delaunay-based quality tetrahedral mesh generator*, ACM TOMS 41 (2015).
- S.-W. Cheng, T. K. Dey and J. R. Shewchuk, *Delaunay Mesh Generation*, CRC Press (2012).
- S. J. Owen, *A survey of unstructured mesh generation technology*, 7th International Meshing Roundtable (1998).
- Open CASCADE Technology, *Mesh* user's guide (BRepMesh, linear and angular deflection).

No code from Triangle, TetGen, Gmsh, CGAL or the CDT library is read or vendored; see [`THIRD_PARTY_NOTICES.md`](https://github.com/loumalouomega/meshioplusplus/blob/main/THIRD_PARTY_NOTICES.md).
