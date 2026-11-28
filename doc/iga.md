# Isogeometric layer: evaluate, refine, sample

NURBS models reach meshio++ as a sampling: the [MFEM](/formats/mfem#the-patches-themselves-the-nurbs-side-channel) reader samples each patch into Lagrange cells, and [Kratos CAD JSON](/formats/cad_json), [STEP](/formats/step#reading) and [IGES](/formats/iges) sample trimmed faces. What the sampling drops (knot vectors, weights, control nets, trimming loops) stays beside the mesh as a side channel: `mesh.info["nurbs"]` (the C++ `NurbsInfo`) or `mesh.info["brep"]` (`BrepInfo`). The isogeometric layer acts on that model directly: **`iga_evaluate`** evaluates a patch with derivatives, **`iga_refine`** inserts knots and elevates degrees exactly, and **`iga_sample`** turns a model into Lagrange cells again.

It lives beside `Mesh` (the model is the side channel, not a second entity kind inside the mesh), so no mesh backend, `Mesh` layout or binding of the mesh changes. A mesh that is not NURBS has no model, and the operations say so.

```python
import numpy as np
import meshioplusplus as mp

mesh = mp.read("pipe-nurbs.mesh", file_format="mfem")      # the sampled cells + mesh.info["nurbs"]
fine = mp.iga_refine(mesh, subdivisions=2, elevate=1)      # exact: same geometry, finer net, degree + 1
mp.write("fine.mesh", fine, file_format="mfem")            # an MFEM NURBS mesh again, patches intact

ev = mp.iga_evaluate(fine, np.array([[0.3, 0.7, 0.4]]), patch=0, derivatives=1)
ev["points"], ev["derivatives"], ev["fields"]               # position, d x / d u_k, NURBS fields

face = mp.read("plate_with_hole.cad.json")
mp.iga_evaluate(face, [[0.5, 0.5], [2.0, 2.0]])["inside"]   # [True, False]: the trimming loops' test
mp.iga_sample(face, subdivisions=8)                         # the face as Lagrange quads, trim clipped
```

## Operations

| Operation | What it does |
|---|---|
| `iga_evaluate(mesh, params, patch=0, derivatives=0)` | Points, NURBS fields and, with `derivatives=1`, first derivatives `(N, d, space dim)` and unit normals of a surface in 3-D, and with `derivatives=2` also the symmetric second derivatives `(N, d, d, space dim)`, at the parametric coordinates `params` `(N, d)`. For a B-rep face, `inside`: whether each point lies in the face's trimming loops (even-odd over all of them). Rational, in 1, 2 or 3 parametric dimensions. |
| `iga_refine(mesh, subdivisions=1, elevate=0, sample_subdivisions=0)` | Degree elevation by `elevate`, then knot insertion splitting every nonzero knot span in `subdivisions` parts: the refined model, sampled again, with the refined model in `info`. |
| `iga_sample(mesh, subdivisions=0)` | The model as Lagrange cells. |
| `iga_remove_knots(mesh, direction=0, knots=(), tolerance=0.0, patch=0)` | Knot removal on one patch ([below](#knot-removal)). |
| `iga_closest_point(mesh, points, patch=0)` | The closest point of a curve or surface patch to each query point ([below](#closest-point)). |
| `iga_bezier_extract(mesh, patch=0)` | One patch of a NURBS model as VTK Bezier cells ([below](#bezier-extraction)). |

`mesh` is a mesh whose `info` holds a model; `model="nurbs"` or `"brep"` picks one when it holds both, and a bare model dict is accepted too.

## Refinement is exact

Both refinements are exact rearrangements of the same NURBS: the geometry, the parametrisation and the NURBS fields do not change (to rounding), and a trimmed face's trimming curves, which live in the face's parameter space, stay valid and untouched.

- **h-refinement** is knot insertion (The NURBS Book, A5.4): each interior knot of every span is added `subdivisions − 1` times, equally spaced; continuity at the existing knots is kept.
- **p-refinement** is degree elevation (A5.9): the degree rises by `elevate` and every knot's multiplicity rises with it, so the continuity across it stays as it was. Elevating first and then inserting is k-refinement, the highest continuity for the degree.
- A control net is held *homogeneous* while it is refined, `(w x, w f…, w)`, so the weights and every NURBS field refine with the geometry and one rounding. A net whose weights are all exactly 1 stays polynomial: its weights are written back as 1.
- Refining is uniform over every patch of a model. An MFEM NURBS mesh shares knot vectors between patches through its `edges` table and orients each patch's directions along its own edges, so it refines every direction alike (`subdivisions=2`, not `(2, 3)`); a CAD face, and a model you build yourself, may refine each direction on its own: `subdivisions=(2, 3)`, `elevate=(0, 1)`.
- A face's Kratos control-point ids (`control_point_ids`) have no meaning on a refined net and are cleared; edge curves and vertices are untouched.

The refined model writes back through the format it came from: `mp.write("fine.mesh", fine, file_format="mfem")` writes an `MFEM NURBS mesh` from `fine.info["nurbs"]`, `mp.write("fine.cad.json", fine)` the refined B-rep. Any other format gets the sampled cells.

## Bezier extraction

`iga_bezier_extract(mesh, patch=0, return_operators=False)` (C++ `iga_bezier_extract(info, patch)` and `iga_bezier_extract_operators(info, patch)`, `meshioplusplus iga-bezier`, the MCP tool of the same name and the flat bindings) takes one patch of a NURBS model or one face of a B-rep model to the piecewise-Bezier form finite element codes assemble on. Every interior knot is raised to multiplicity `p` in each direction by exact knot insertion on the homogeneous net (a rational patch stays exact), and the refined net is cut at its nonzero knot spans: one `VTK_BEZIER_CURVE`, `VTK_BEZIER_QUADRILATERAL` or `VTK_BEZIER_HEXAHEDRON` cell per span product, with `(p+1)^d` points in the node order of the same-shape Lagrange cell, in Cartesian coordinates (the homogeneous values divided by the weight). Points are not shared between cells, as VTK's Bezier cells expect.

The cells carry the cell data `HigherOrderDegrees` (Float64 `(cells, 3)`, an unused direction 0), `iga:patch` and `iga:span` (as `iga_sample` writes them), the point data `RationalWeights` when the patch is rational (any weight other than 1) and one point array per NURBS field. A VTU file keeps them as Bezier cells and a reader or ParaView evaluates them exactly. A patch whose directions have different degrees is refused by name when emitting VTK Bezier cells, as VTK's Bezier cells take one order: raise them with `iga_refine(..., elevate=...)` first. Passing `return_operators=True` (C++ `iga_bezier_extract_operators`) also returns the per-element Borden–Scott–Evans–Hughes (2011) extraction matrices `C^e` of shape `(num_elements, n_local, n_local)` together with `ien`, `degrees`, and `element_spans` so an IGA solver can assemble directly on Bernstein shape functions, including mixed-degree patches and B-rep faces.

## Knot removal

`iga_remove_knots(mesh, direction=0, knots=(), tolerance=0.0, patch=0)` (C++ `iga_remove_knots(info, patch, direction, knots, tolerance, report)`, `meshioplusplus iga-remove-knots`, the MCP tool of the same name and the flat bindings) is the inverse of knot insertion, The NURBS Book's A5.8 run on the homogeneous net. Interior knots of one direction of one patch (or B-rep face) are removed while the geometry, and the NURBS fields, move by at most `tolerance`: the book's bound `d * wmin / (1 + Pmax)` on the homogeneous net (`wmin` the smallest weight, `Pmax` the largest control-point norm) keeps the curve or surface within `d` of the original. A tolerance of 0 removes only what is exact up to rounding (`1e-12` of the net's size), which is what undoes a knot insertion exactly, a rational patch included.

`knots` lists the values to remove, one copy per entry (a value listed twice asks for two); empty asks for every copy of every interior knot. A knot of a surface or volume goes only if it goes from every line of the net, as the net must stay a tensor product, so the operation finds the number of copies every line can lose first and removes exactly that many. On a multi-patch NURBS model with shared knot vectors (joined patches), the operation finds the common removable multiplicity across every patch sharing that knot vector and removes the same copies from all of them together. The report gives, per distinct knot value, the copies asked for and removed. The result is the model with the modified patch(es) replaced and sampled again (the cell maps cleared, as `iga_refine` does). A knot that is not an interior knot of the direction, or that is listed more often than its multiplicity, is refused by name.

## Closest point

`iga_closest_point(mesh, points, patch=0)` (C++ `iga_closest_point(info, patch, points)`, `meshioplusplus iga-closest-point`, the MCP tool of the same name and the flat bindings) finds, for each query point, the point of one patch (or B-rep face) nearest to it: the parameters, the point and the distance. It takes a curve (1 parametric direction), surface (2 directions), or volume (3 directions), rational or not. A lattice of samples (`degree + 1` per knot span and direction, capped for surfaces and volumes) seeds a Newton iteration on the squared distance with the exact second derivatives of `iga_evaluate(derivatives=2)`, clamped to the knot range (a point beyond an edge or face lands on the boundary), falling back on Gauss-Newton where the Hessian is not positive definite and halving the step until the distance drops. The three nearest samples seed three runs and the closest result wins, so a point near a fold of the patch usually finds the global minimum, though no global guarantee is given, and a point equidistant from several places of the patch gets one of them. Each point is independent, so the result does not depend on the thread count.

On a surface in 3-D it agrees with the projector behind [`project_to_cad`](./project_to_cad.md) to the projector's own tolerance, and it is the more accurate where a point lands on an edge.

## Sampling

`iga_sample(mesh, subdivisions)` puts `subdivisions` cells in every nonzero knot span and direction (`0`, the default, is automatic).

- **An MFEM model with its patch topology** samples as the MFEM reader does (`subdivisions=0`): the same cells, nodes shared across patches, boundary cells, so a refined model is the mesh `read_mfem` returns for the refined file. With an explicit `subdivisions`, or a model you built yourself, the generic sampler runs instead.
- **The generic NURBS sampler** makes one block of the highest degree's Lagrange cells (`line`, `line3`, `VTK_LAGRANGE_CURVE`; `quad`, `quad9`, `VTK_LAGRANGE_QUADRILATERAL`; `hexahedron`, `hexahedron27`, `VTK_LAGRANGE_HEXAHEDRON`, by parametric dimension), the nodes the exact rational geometry at each cell's lattice, the NURBS fields as point data, and `iga:patch` and `iga:span` cell data. Points are shared within a patch, not between patches; the model's `cell_patch` and `cell_span` maps are filled and `boundary_cell_patch` emptied.
- **A B-rep** samples as [`cad_json.sample`](/formats/cad_json#the-mesh) does (`trim_samples`, `edges` and `trim` as there).

## Surfaces

| Surface | Entry points |
|---|---|
| Python | `meshioplusplus.iga_evaluate`, `iga_refine`, `iga_sample` (core with a bit-identical numpy twin; `MESHIOPLUSPLUS_STRICT_CORE=1` proves the core handles a call) |
| C++ | `operations/iga.hpp`: `iga_evaluate`, `iga_refine`, `iga_sample` on `NurbsInfo` and `BrepInfo`, and the patch-level `iga_insert_knots` and `iga_elevate_degree` (one direction at a time, arbitrary patches) |
| C | `mio_iga_refine`, `mio_iga_sample` and `mio_iga_evaluate` on the [`mio_format_info`](/c_api) handle of `mio_read_with_info`; `mio_iga_opts` |
| Fortran | `info%iga_refine`, `info%iga_sample`, `info%iga_evaluate` and the `mio_iga_evaluation` result type |
| Julia | `iga_refine(info)`, `iga_sample(info)`, `iga_evaluate(info, patch, params)` |
| R | `mio_iga_refine()`, `mio_iga_sample()`, `mio_iga_evaluate()` |
| WebAssembly | `igaRefine(info, options)`, `igaSample`, `igaEvaluate` on `readMeshSelective(path, {info: true}).info` |
| CLIs | `iga-refine`, `iga-sample`, `iga-evaluate` (both CLIs); `meshioplusplus iga-refine pipe.mesh fine.mesh -o mfem -s 2 -e 1` |
| MCP | `iga_refine`, `iga_sample`, `iga_evaluate` |

The [pipeline](/pipeline#the-rules) carries the model a CAD JSON, STEP, IGES or MFEM NURBS Input reads to its `IgaSample` and `IgaRefine` steps (v17.17.0), the refined model on to the steps after `IgaRefine`, and drops it at the first step that changes the mesh; `iga_evaluate`, which returns arrays rather than a mesh, has no step.

## What this is not

- **Not a second entity kind in `Mesh`.** Operations on a mesh (`refine`, `smooth`, `partition`, …) still see only the sampling; the model changes only through these three operations. Putting patches and trims inside `Mesh` itself (every backend, binding and the ABI) is a recorded [non-goal](/roadmap#non-goals-and-decisions-taken) that a consumer needing it can reopen.
- **Not a mesher.** `iga_sample` keeps the cells of a trimmed face whose centre its loops contain, so its trimmed edges are staircases; the loops themselves stay exact. A conforming triangulation that follows the trims is [`mesh_cad_surface`](/mesh_cad), on the same side channel.
- **Not yet:** refinement of one patch of an MFEM model (the patch-level `iga_insert_knots` and `iga_elevate_degree` act on a single `NurbsPatch`, leaving the neighbours to the caller).

## Verification

### Second derivatives

`derivatives=2` adds `second_derivatives`, the full symmetric `d2 x / d u_k d u_l` of shape `(N, d, d, space dim)`, for curves, surfaces and volumes, rational or not (rational ones by the quotient rule `R_kl = (((A_kl - W_kl R) - W_k R_l) - W_l R_k) / W`). The first derivatives of order 2 come from the second-order basis recurrence and may differ from those of `derivatives=1` in the last bits. C: `MIO_IGA_SECOND_DERIVATIVES` of `mio_iga_array` (`derivatives` 0 to 2). The Python twin repeats the arithmetic order of the core, so both engines are bit-identical.

The tests refine all six NURBS fixtures of the [MFEM](/formats/mfem) reader and every [CAD JSON](/formats/cad_json) fixture by h, p and both, on dense parameter grids: points, first derivatives, NURBS fields and the trimming test agree with the unrefined model to 1e-11, rational circles and spheres stay on their analytic surface to 1e-14, derivatives match central differences, a refined MFEM model writes and reads back, and the C++ core and the numpy twin return the same bits (`tests/python/test_iga.py`, `tests/cpp/test_iga.cpp`).
