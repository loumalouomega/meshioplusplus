# Fit CAD

`fit_cad(mesh, path)` turns the fitted patches of a closed triangle skin into a CAD solid: every patch that [`fit_primitives`](/fit_primitives) accepted becomes a plane, sphere, cylinder, cone or torus face trimmed by the patch boundary, a rolling-ball fillet or blend (of one radius or a radius varying along it) or a helical thread face its own [feature face](#fillets-blends-and-threads), any other patch no primitive explains a [freeform B-spline face](#freeform-faces), the faces are sewn into a solid that OpenCASCADE checks, and the solid is written as STEP. It is the last step of [skin to CAD](/cad_fit), after [`segment`](/segment) and `fit_primitives`, and it needs the optional [OpenCASCADE backend](#getting-opencascade).

```python
import meshioplusplus as mp

skin = mp.repair(mp.read("part.stl"))
seg = mp.segment(skin, curvature_tolerance=1.0, min_patch_cells=20)
fit = mp.fit_primitives(seg, tolerance=2e-3)

solid, report = mp.fit_cad(fit, "part.step", return_report=True)
report["valid"], report["num_cylinder"], report["num_demoted"]
report["volume_ratio"], report["max_deviation"]   # against the skin
```

![The turned part of skin to CAD as a CAD solid: two planes, a cylinder, a torus fillet, a cone and a spherical dome, re-meshed and coloured by face type](/images/fit_cad_solid.png)

Unlike the faceted [`step`](/formats/step) writer, whose solid is one flat face per triangle, the STEP file `fit_cad` writes carries `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`, `SPHERICAL_SURFACE` and `TOROIDAL_SURFACE` faces a CAD tool can select, dimension and edit. A skin without `segment:patch` or the `fit:patch_*` table is segmented and fitted first with the default options; a noisy scan wants a curvature test and a `min_patch_cells` in the tens, so run those two yourself, or pass them to the [CLI](#cli).

## How the solid is built

- **Faces.** A patch whose fit was accepted, and whose type is in `types`, becomes one face on its fitted surface. A patch no primitive explains is tried, in order, as a [rolling-ball fillet](#fillets-blends-and-threads) between two analytic faces (`fillets`), as a [helical thread face](#fillets-blends-and-threads) (`threads`) and, with `freeform` (the default), as a [B-spline face](#freeform-faces); once the freeform faces exist, a patch between two fitted faces, one of them freeform, is tried as a [rolling-ball blend](#fillets-blends-and-threads) against them (`fillets`). Every other surface cell (a type left out, a freeform fit that failed, a cell with no patch) becomes planar facets, one per cell, as the faceted writer writes them.
- **Edges.** Where two patches meet, the mesh edges between them form chains from one junction (a point where three or more patches meet) to the next; a closed chain is split in two. Between two analytic faces the edge is the branch of the two surfaces' exact intersection that lies within the tolerance of the chain, when `exact_edges` is set and one does (a plane meeting a cylinder gives a circle). Otherwise, notably where a fillet meets its wall tangentially and no intersection is stable, the edge is a B-spline through the chain, smoothed and moved onto both surfaces, kept only if it stays within the tolerance of both. Beside a facet, the edge is the mesh's own straight edge.
- **Vertices.** A junction among analytic faces only is moved onto all of their surfaces; one beside a facet keeps the mesh position.
- **Solid.** Each face gets parameter-space curves that are continuous around its loops (so a boundary running along a cylinder's seam or a torus's equator does not fall apart), OpenCASCADE adds the seams and repairs the faces, sews them and closes each shell into a solid, one per edge-connected piece of the skin. A patch with no boundary at all, a whole sphere or torus, is its surface's natural face.

## Freeform faces

A patch `fit_primitives` left unexplained (no accepted primitive, or one with non-finite parameters) gets a B-spline surface fitted to it, so a smooth blend, a dent or a styled panel becomes one CAD face instead of hundreds of facets:

1. **A parametrisation.** The patch is flattened over its best-fit plane (principal components) or, when that folds, unrolled around a cylinder: the axis is the direction its normals are most perpendicular to (the smallest eigenvector of Σ nnᵀ), the radius a least-squares circle fit of the points around it, the angle cut at the widest gap. A flattening whose triangles do not all keep one orientation is not a parametrisation and is not used. A patch that folds over both (a doubly curled one, or a band closed around an axis) is, with `freeform_split` (the default), [split](#splitting-curled-patches) into pieces that one of them does parametrise.
2. **A height field.** The offset from the plane or the cylinder is resampled on a grid over the flattened bounding box, interpolated inside the patch's triangles and continued from the nearest boundary point outside them, so the surface extends past the patch's trim.
3. **The fit.** `GeomAPI_PointsToBSplineSurface` fits the grid (iso-parametric, C2, degrees up to `freeform_degree`, default 3, valid 3–8). The grid starts near √n points a side and doubles up to `freeform_max_control` (default 64, valid 4–1024) until every patch vertex lies within half the tolerance and every facet centre within the tolerance; a patch that never gets there, or a skin too coarse for any smooth surface to pass within it, stays faceted.
4. **The face.** The surface is oriented by the facet normals and goes through the same trimming, sewing, check and [demotion](#checks-and-demotion) as an analytic face; its `fit_cad:type` is 6 (`kFitCadBsplineType`), and the report's `num_bspline` counts it.

`freeform=False` keeps every unexplained patch faceted.

### Splitting curled patches

A patch that folds over both flattenings is halved at the area-weighted median of `n · e`, where `e` is the direction its normals spread along most (the principal axis of the area-weighted Σ nnᵀ): a cut along an iso-line of the normal, a smooth curve on a smooth surface. Each half is split into its edge-connected components, and a component that still folds is halved again, up to four times. Each piece that now flattens becomes its own freeform face under a new patch label, numbered after the largest existing one in the order of its lowest facet, so `fit_cad:patch` tells the pieces apart; a piece under 16 triangles, or one still folding at the fourth halving, stays faceted. The edges between pieces are approximated B-splines like any other tangent junction. The split is deterministic and needs no OpenCASCADE (`detail/patch_flatten`); `num_split_patches` counts the pieces, and `freeform_split=False` keeps such a patch faceted, as before.

![A capped band whose radius waves along its axis: the wall, closed around the axis, is split into two freeform B-spline faces between two planar caps](/images/fit_cad_split.png)

## Fillets, blends and threads

Features that no primitive explains are recognised from the patches around them, before the freeform fallback or, for a blend, after it. The fits need no OpenCASCADE (`detail/feature_fit`); the faces they give are then trimmed, sewn, checked and demoted like any other.

- **Rolling-ball fillets** (`fillets`, type 7). A patch that runs between two analytic faces (the two neighbours sharing the most boundary, together at least 60% of it) is tried as a constant-radius fillet between them: each interior vertex, moved by `r` against or along its normal, should be a ball centre at `r` from both primitives. The radius minimises that misfit over a log-spaced scan of `[tolerance / 2, diagonal]` refined by golden section, for a convex (centres inside) and a concave (centres outside) ball. The face is the pipe of that radius around the spine where the two surfaces, each offset by `r` towards the centres, meet (`GeomFill_Pipe` around the `GeomAPI_IntSS` branch nearest the centres, kept whole and periodic when the centres go most of the way round it). It covers a straight rounded edge, a swept fillet along a curved edge (the foot of a leaning boss, along an ellipse) and a ring. A constant-radius fillet between coaxial surfaces is a torus, which [`fit_primitives`](/fit_primitives) already fits; this stage only sees what it left unexplained.
- **Fillets of varying radius** (`fillets`, type 7, counted in `num_variable_fillet`). When no one radius explains the patch, each interior vertex gets its own ball radius, and the radii become a law along the rail, the patch's boundary on the first surface: the parameter of a vertex is where its own ball touches that surface, a Chebyshev series (a Fourier series round a closed rail) of the lowest degree up to 8 that explains the radii within a quarter of the tolerance. The face is a B-spline through cross-sections of the ball family: along the rail, the ball touching both surfaces at the law's radius is solved for in the plane through the rail point, and each section is the arc of the ball's characteristic circle from one foot to the other, widened by 15% each way, the sections running a tenth of the rail past each open end so the face trims inside its surface.
- **Blends beside freeform faces** (`fillets`, type 7, counted in `num_blend`). Once the [freeform faces](#freeform-faces) exist, a patch between two fitted faces, at least one of them freeform (plane–B-spline or B-spline–B-spline), is tried as a rolling ball against those surfaces, with the radius found per vertex and fitted as a law exactly as above, constant or varying; a freeform patch can become a blend this way. The neighbours are taken as the freeform stage left them, so the result does not depend on the order the patches are tried in. A ball resting on a nearly flat stretch of a fitted face touches it where its normal happens to point, which wanders far more than the face does, so the balls' centres and radii, not their feet, set the face: smooth laws are fitted through the solved centres (judged across the spine, as a wavering rail stretches the arc length along it) and the arcs ride on them. A blend stays tangent to its fitted neighbours, which lie within half the tolerance of the skin, so a vertex on the blend's own boundary may lie the whole tolerance from it.
- **Helical thread faces** (`threads`, type 8). A patch beside an accepted cylinder or cone is first tried as a straight-flank thread about its axis: the ruled helicoid `z = c + k θ + t ρ` (`k` the pitch over 2π, `t` the flank's slope), one least-squares solve once each vertex's angle is unwrapped over the patch's edges; the face is a B-spline through a grid on it. A patch that is not such a flank, or that has no cylinder or cone beside it, is tried as a **helical surface of any profile**. Without an axis, the screw comes from the patch itself, as the screw motion whose velocity field is most nearly tangent to it (Pottmann and Randrup, *Rotational and helical surface approximation for reverse engineering*, Computing 60, 1998): the smallest eigenvector of a normal-weighted 6 × 6 moment, reduced to 3 × 3, refused when a plane, cylinder or sphere leaves more than one screw nearly free, or when the advance vanishes (a surface of revolution). The axis and the advance are then refined together by Levenberg-Marquardt on the residual of the profile: each vertex, with the screw's advance taken off its height, lands in the plane through the axis, and the profile is a Chebyshev series of the lowest degree up to 8, `ζ(ρ)` or `ρ(ζ)`, within 0.4 of the tolerance. A straight profile is a flank; any other (a rounded root or crest) is counted in `num_rounded_thread`. A trapezoidal or ISO thread then comes back as cylinders (root and crest) and two helical flanks, and a thread with rounded roots and crests, with no cylinder anywhere, as four helical faces.

Every patch must span at least half a turn to be a thread (less of a helical surface is as well something else). A candidate face must keep every patch vertex within half the tolerance and every facet centre within it, as a freeform face must. Edges beside a fillet, a blend or a thread face skip the exact intersection, which is slow and unstable between tangent or long B-spline surfaces, and go straight to the approximated B-spline. `num_fillet` and `num_helical` count the faces, and of them `num_variable_fillet`, `num_blend` and `num_rounded_thread` the new kinds; `fillets=False` and `threads=False` leave such patches to the freeform fit.

![A boss leaning on a round plate, its foot rounded by a fillet swept along an ellipse: the fillet is one rolling-ball face between the plate's plane and the boss's cylinder](/images/fit_cad_fillet.png)

![A block whose edge is rounded by a ball shrinking from 0.3 to 0.2 along it: one fillet face of varying radius between two planes](/images/fit_cad_variable_fillet.png)

![A bar whose wavy top and leaning front are freeform B-spline faces, the edge between them one rolling-ball blend](/images/fit_cad_blend.png)

![A rod with a trapezoidal thread: root and crest are cylinders, both flanks helical faces, the shallow end caps facets](/images/fit_cad_thread.png)

![A rod whose thread has rounded roots and crests and no cylinder anywhere: four helical faces, found from the screw each patch fits](/images/fit_cad_rounded_thread.png)

## Checks and demotion

A face is not kept when OpenCASCADE's B-rep check rejects it, when sewing leaves an edge of it free, when its boundary cannot be built, or when it strays: the solid is re-meshed and a face with a point further than `tolerance` from the skin, or a patch with a skin point further than that from the solid, fails. A failed patch is **demoted** to facets and the solid is rebuilt; after `max_demotion_rounds` rounds (default 8) every remaining patch is faceted, which is always a valid solid for a closed skin. So a feature the fitter cannot explain, or a face the kernel cannot trim, comes back faceted, never silently wrong, and `num_demoted` counts the features dropped. `MESHIOPLUSPLUS_LOG_LEVEL=debug` logs why each patch was demoted.

`tolerance` defaults to three times the largest `fit:patch_tolerance` of an accepted patch, so it follows the scale `fit_primitives` worked at; with no accepted patch it is `1e-3` of the bounding-box diagonal. The re-meshing uses a chordal deflection of a quarter of the tolerance unless `mesh_deflection` says otherwise.

## Output and report

The return value is the solid re-meshed by OpenCASCADE as one `triangle` block, with points welded where faces share them exactly, and Int64 cell data `fit_cad:face` (the B-rep face), `fit_cad:type` (0 a facet, else the `fit:type` value: 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus; 6 a freeform B-spline, 7 a rolling-ball fillet or blend, 8 a helical thread face) and `fit_cad:patch` (the `segment:patch` label, -1 for a cell with no patch). With `return_report=True` the report holds:

| Key | Meaning |
| --- | --- |
| `valid` | the solid passes OpenCASCADE's B-rep check |
| `tolerance` | the tolerance applied |
| `volume`, `skin_volume`, `volume_ratio` | the solid's volume, the skin's enclosed volume and their ratio |
| `max_deviation` | the symmetric Hausdorff distance between the skin and the re-meshed solid |
| `num_faces`, `num_plane` … `num_torus`, `num_bspline`, `num_fillet`, `num_helical`, `num_faceted_faces` | B-rep faces in all, by type (freeform B-splines, fillets and thread flanks included), and planar facets |
| `num_variable_fillet`, `num_blend` | of the fillet faces, those whose radius varies along them, and the blends beside a freeform face |
| `num_rounded_thread` | of the helical faces, those whose profile is not a straight line (a rounded root or crest) |
| `num_split_patches` | pieces that freeform patches folding over every flattening were split into |
| `num_patches`, `num_demoted` | patches in the fit table, and accepted ones demoted to facets |
| `num_exact_edges`, `num_approx_edges` | edges between analytic faces built as exact intersections and as B-splines |
| `num_solids` | solids written, one per piece of the skin |

The STEP file is AP214 in millimetres. Its header timestamp is pinned and its product is named `meshio++ fit_cad`, so the same input writes the same bytes; the bytes do depend on the OpenCASCADE version. The work is serial, so nothing depends on the thread count.

## The solid as NURBS

`solid.info["brep"]` (C++ `FitCadResult::mBrep`, a `BrepInfo`) holds the solid with every face converted to NURBS (`BRepBuilderAPI_NurbsConvert`): the surfaces, their trimming loops and the edges with the faces they bound, as a [STEP read](/formats/step#reading) gives them. Every face of the solid is there, trimmed, and every trimming curve lies in its surface's parameter range. When the exact conversion does not deliver that (OpenCASCADE 7.9 throws on a trimmed cylinder or one wound several turns round, leaves a face that already was a B-spline without a surface, and can stop a converted surface short of its trims), the solid is kept as it is and each surface and edge that is not a B-spline is approximated by one in its own parameters, within `1e-6` of the bounding-box diagonal, so the trimming curves stay valid; planes and lines stay exact. Writing it as [Kratos CAD JSON](/formats/cad_json) (`mp.write("part.cad.json", solid)`, `to-cad --cad-json`) makes a scanned part Kratos IGA geometry.

## Getting OpenCASCADE

[Open CASCADE Technology](https://github.com/Open-Cascade-SAS/OCCT) is licensed under the GNU LGPL 2.1 with the additional exception in its `OCCT_LGPL_EXCEPTION.txt`. meshio++ stays MIT: the backend is **off by default and bring-your-own** (never vendored or downloaded), linked dynamically, and absent from the PyPI wheels, the WASM build and the release CLI binaries. A build that enables it links an LGPL library; distributing that build carries the LGPL's obligations for OpenCASCADE, not for meshio++. OpenCASCADE 7.6 or later is supported (7.6 is Ubuntu 24.04's package; the CI leg tests it).

```bash
sudo apt-get install libocct-foundation-dev libocct-modeling-data-dev \
    libocct-modeling-algorithms-dev libocct-data-exchange-dev libtbb-dev
CMAKE_ARGS="-DMESHIOPLUSPLUS_WITH_OCCT=ON" pip install meshioplusplus --no-binary meshioplusplus
```

From source, `build/configure.sh --with-occt` (or `-DMESHIOPLUSPLUS_WITH_OCCT=ON`) finds OpenCASCADE through its CMake config; point `OpenCASCADE_DIR` or `CMAKE_PREFIX_PATH` at an install outside the default prefixes, and make its libraries findable at run time as for any shared library. Debian's packages export toolkit targets that link TBB, hence `libtbb-dev`.

Without the backend, `fit_cad` raises `NotImplementedError` (C++ `Unsupported`) naming `-DMESHIOPLUSPLUS_WITH_OCCT=ON`; there is no pure-Python fallback, because trimming, sewing and checking surfaces is a CAD kernel's work. `meshioplusplus._fit_cad.has_occt()`, `fit_cad_has_occt()` in C++, `mio_has_occt()` in C and the language surfaces below say in advance whether a build has it.

## CLI

```sh
meshioplusplus to-cad IN OUT.step [--curvature-tolerance C] [--min-patch-cells N] \
    [--fit-tolerance T] [--tolerance T] [--no-exact-edges] [--types plane,cylinder,...] \
    [--mesh-deflection D] [--max-demotion-rounds N] [--no-freeform] [--freeform-degree P] \
    [--freeform-max-control N] [--no-freeform-split] [--no-fillets] [--no-threads] [--mesh-out FILE] [--cad-json FILE] [--json] [--quiet]
```

A skin without `segment:patch` is segmented with `--curvature-tolerance` and `--min-patch-cells`, and fitted with `--fit-tolerance` (relative, as in `fit-primitives`), before the solid is built. `--mesh-out` also writes the re-meshed solid, `--cad-json` the solid as Kratos CAD JSON; `--no-freeform`, `--freeform-degree`, `--freeform-max-control` and `--no-freeform-split` are the [freeform](#freeform-faces) options, `--no-fillets` and `--no-threads` the [feature](#fillets-blends-and-threads) ones. Both CLIs exit with 2 when the solid fails the B-rep check, as `check` does for a failed gate.

```sh
meshioplusplus to-cad part.stl part.step --curvature-tolerance 1 --min-patch-cells 20 --fit-tolerance 2e-3
```

## Pipeline

```json
{"Op": "FitCad", "Path": "part.step", "Tolerance": 0.0, "ExactEdges": true,
 "Types": ["plane", "cylinder"], "MeshDeflection": 0.0, "MaxDemotionRounds": 8,
 "Freeform": true, "FreeformDegree": 3, "FreeformMaxControl": 64,
 "FreeformSplit": true, "Fillets": true, "Threads": true}
```

The step writes the solid to `Path` when given, and the mesh continues as the re-meshed solid; it reports `Valid`, `VolumeRatio`, `MaxDeviation`, `NumFaces`, `NumBspline`, `NumFillet`, `NumHelical`, `NumVariableFillet`, `NumBlend`, `NumRoundedThread`, `NumSplitPatches`, `NumFacetedFaces` and `NumDemoted`. A `Tolerance` or `MeshDeflection` of 0 is automatic; an empty `Types` facets everything.

## Other languages

```c
mio_fit_cad_opts opts;
mio_fit_cad_opts_init(&opts);
opts.types = (1u << MIO_FIT_PLANE) | (1u << MIO_FIT_CYLINDER);
opts.freeform = 1;             /* the default; freeform_degree, freeform_max_control */
mio_fit_cad_report report;
if (mio_has_occt()) {
    mio_mesh* solid = mio_fit_cad(fitted, "part.step", &opts, &report);
}
```

```fortran
solid = fit%fit_cad(step_path='part.step', valid=ok, num_demoted=ndropped, stat=ierr)
```

```julia
has_occt() && (r = fit_cad(fit; step_path="part.step"); r.valid)
```

```r
if (mio_has_occt()) cad <- mio_fit_cad(fit, step_path = "part.step")
```

The C report carries `num_bspline`, `num_fillet`, `num_helical`, `num_variable_fillet`, `num_blend`, `num_rounded_thread` and `num_split_patches`, as do Fortran's, Julia's and R's results (`freeform`, `freeform_degree`, `freeform_max_control`, `freeform_split`, `fillets` and `threads` are their options; a zero-filled C struct keeps the last three off, `mio_fit_cad_opts_init` turns them on). The MCP server's `fit_cad` tool takes the same options plus `curvature_tolerance`, `min_patch_cells`, `fit_tolerance` and `cad_json_output_path`, and reports `occt_available`; without the backend it answers with the named error and `occt_available: false`. There is no WASM binding.

## Not supported

- **Vertex blends, chamfers and blends of a section other than a circle**: a fillet or blend is a rolling ball between exactly two faces; where three or more blends meet (a setback vertex), or the section is a flat chamfer or a conic, the patch comes back freeform or faceted. They are recorded as a [non-goal](/roadmap#non-goals-and-decisions-taken) until a consumer asks.
- **Threads that are not screw-invariant**: a tapered (pipe) thread or one of varying pitch is no helical surface; it comes back freeform or faceted, as does a thread patch spanning less than half a turn.
- **A whole thread as one patch**: a thread surface whose profile spans the full pitch (one patch over root, flanks and crest) has no single trimmed face to be; [`segment`](/segment) it at its creases or with a curvature test first.
- **Fitting a STEP or IGES model**: those are [read](/formats/step#reading) as their own NURBS B-rep, not refitted.
