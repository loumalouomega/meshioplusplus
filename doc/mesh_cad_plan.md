# Mesh generation from CAD: milestone 1 plan

This is the plan for `mesh_cad_surface`, a new operation that triangulates a CAD B-rep's faces into one conforming surface mesh at a single global target edge size, without a CAD kernel. It is milestone 1 of a longer programme (refinement, size fields, anisotropic metrics, volume meshing; see the [roadmap](./roadmap.md)). Every existing operation keeps its behaviour and its output byte for byte.

**Status.** Shipped in v17.13.0; the contract is [mesh generation from CAD](./mesh_cad.md), and what followed it is on that page (the roadmap's mesh-generation section closed in October 2026; its parked notes are [measured and parked](./mesh_cad.md#measured-and-parked)). This page is kept as the record of the plan; where the two differ, the contract wins (the per-face metric map and the loop chaining, for instance, were added while meshing the real fixtures).

## What exists

- **The model.** CAD JSON (natively) and STEP and IGES (through OpenCASCADE) are read into the `BrepInfo` side channel (`brep.hpp`, `mesh.info["brep"]` in Python). Every face is a `NurbsPatch`; analytic surfaces arrive as their exact rational NURBS. Every trim of a face's loops carries a 2-D curve in the face's parameters (`BrepTrim::mCurve`), and edges are shared topology: `BrepEdge::mTopology` lists the `(face id or name, trim index, relative direction)` of each use. A 3-D edge curve is optional. Seams are explicit trims (one edge used twice in a face, a period apart) and poles are degenerate trims whose edge has no curve, so a periodic face is already cut open in its parameters. Vertices are filled by CAD JSON only, and there is no tolerance field.
- **Evaluation.** `detail/iga_kernels.hpp` evaluates rational curves and surfaces with first derivatives; `detail/brep_trim.hpp` samples loops and has an even-odd inside test.
- **Sampling, not meshing.** `iga_sample` and `cad_json`'s sampler put a quad lattice on each face and keep the cells whose centre is inside the trim: the boundary is a staircase and faces share no points. The only triangulator is OpenCASCADE's `BRepMesh`, in builds that link it.
- **No predicates.** Nothing in the core decides an orientation or an in-circle test exactly; the roadmap dropped a Delaunay item for that reason and recorded native triangulation as a non-goal, which this milestone narrows.

## Algorithm

1. **Topology.** Every edge use is resolved to a `(face, trim)`, and the reverse map built. A trim no edge names becomes its own edge; unreferenced trims of different faces whose 3-D images coincide (ends and midpoint within the tolerance) are paired. An untrimmed face gets the patch's parameter rectangle as its outer loop.
2. **Vertices from topology.** A union-find joins the end of each trim to the start of the next in its loop and the ends of each edge across its uses, so corners merge exactly even where the trims of two faces disagree by the model's tolerance; the distance is only measured and reported.
3. **Edges, sampled once.** Each edge's first use (lowest face, then trim) is its master: its 3-D image `S(c(t))` is sampled into `max(1, ceil(L / h))` segments of equal arc length, split further while a segment's sagitta exceeds the optional chordal cap. Every other use (the neighbouring face, or the same face across a seam) projects those same 3-D points onto its own trim curve, so it gets its own `(u, v)` while sharing the point ids: this is what makes the faces conform. A zero-length edge (a pole) keeps its two ends, which merge into one vertex.
4. **Faces, triangulated in their parameters.** The loops become closed polylines of `(u, v)` samples, triangulated by a constrained Delaunay triangulation: incremental insertion with Lawson flips, segment recovery by re-triangulating the cavity a missing segment crosses (Anglada 1997), and an even-odd flood fill across the constrained edges to drop the outside and the holes. Interior points come from an equilateral lattice whose spacing in each parameter direction is `h` over the largest `|∂S/∂u|` and `|∂S/∂v|` on the face, inserted only inside the domain and kept only if farther than `h / 2` in 3-D from every boundary sample and every point already kept. The points are lifted with the surface evaluator, triangles a pole collapses are dropped, and the winding follows `∂S/∂u × ∂S/∂v`, reversed for a face with a swapped normal.
5. **Determinism.** Point ids run vertices, then edge interiors in edge order, then face interiors in face order; faces and edges run in parallel into their own buffers and are assembled in index order, so the output is byte-identical across backends and thread counts.
6. **Robust predicates.** Shewchuk's public-domain adaptive `orient2d`, `incircle`, `orient3d` and `insphere`, ported to a core-private header, decide every orientation and in-circle test. Floating-point contraction is switched off inside them, since a fused multiply-add breaks their error-free transformations.

A face that fails (a loop that does not close, crossing boundary segments at this size, a projection that does not converge) is left out of the triangles and reported by index and reason, never meshed wrong; the open-edge count in the report makes the gap visible, and `strict` turns the report into an exception.

## Commits

0. This plan.
1. `detail/predicates.hpp` and its tests: signs against exact rational evaluation, near-degenerate grids, permutation antisymmetry.
2. `detail/cdt2d` and the planar tests: square, L-shape, square with a hole, collinear and near-duplicate points, a sliver, random domains; exact area, every segment an edge, positive counter-clockwise triangles, Euler characteristic, the Delaunay test on unconstrained edges.
3. `detail/edge_sampler` and its tests: lines, rational arcs, equal spacing, the chordal cap, reversed and seam uses sharing point ids, zero-length edges.
4. `operations/mesh_cad` and the CAD tests: a planar square, the trimmed plate, a cube, a capped cylinder and a sphere; watertightness and orientation, area and enclosed volume against the analytic values with second-order convergence, edge lengths and angles, Hausdorff distance, determinism, refusals, and OpenCASCADE's `BRepMesh` as an oracle where it is linked.
5. Python, both CLIs (`mesh-cad-surface`), MCP; generated CAD JSON fixtures for the closed solids.
6. C, Fortran, Julia, R and WebAssembly.
7. Documentation, roadmap, notebook, single header.
8. The release.

## Not in this milestone

No quality guarantee (Ruppert or Chew refinement), no smoothing, no size field or curvature sizing beyond the chordal cap, no anisotropic metric, no quadrilaterals, tetrahedra or boundary layers, no new dependency, no numpy twin (the core-only precedent of `remesh_volume` and `octree_mesh`: a numpy triangulator could not stay bit-identical), and no pipeline step (the pipeline carries no model, the `iga_*` precedent). OpenCASCADE's mesher is a test oracle only. The triangulator follows the published algorithms; no code is taken from Triangle, TetGen, Gmsh, CGAL or the CDT library.

## References

- J. R. Shewchuk, *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates*, Discrete & Computational Geometry 18 (1997).
- J. R. Shewchuk, *Delaunay Refinement Algorithms for Triangular Mesh Generation*, Computational Geometry 22 (2002).
- M. V. Anglada, *An improved incremental algorithm for constructing restricted Delaunay triangulations*, Computers & Graphics 21 (1997).
- C. L. Lawson, *Software for C1 surface interpolation*, Mathematical Software III (1977).
- J. Ruppert, *A Delaunay refinement algorithm for quality 2-dimensional mesh generation*, Journal of Algorithms 18 (1995).
- C. Geuzaine and J.-F. Remacle, *Gmsh: a three-dimensional finite element mesh generator*, IJNME 79 (2009), for the parametric-domain approach.
