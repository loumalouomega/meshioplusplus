# Fit CAD

`fit_cad(mesh, path)` turns the fitted patches of a closed triangle skin into a CAD solid: every patch that [`fit_primitives`](/fit_primitives) accepted becomes a plane, sphere, cylinder, cone or torus face trimmed by the patch boundary, a rolling-ball fillet or blend (of one radius or a radius varying along it) or a helical thread face its own [feature face](#fillets-blends-and-threads), a chamfer, a conic blend, a tapered or variable-pitch thread or a corner where blends meet its own [face too](#chamfers-conic-blends-vertex-blends-and-tapered-threads), any other patch no primitive explains a [freeform B-spline face](#freeform-faces), the faces are sewn into a solid that OpenCASCADE checks, and the solid is written as STEP. It is the last step of [skin to CAD](/cad_fit), after [`segment`](/segment) and `fit_primitives`. It runs on [OpenCASCADE](#getting-opencascade), which the PyPI wheels and the release binaries include, and falls back, in a build without it, to a [kernel-free native backend](#the-native-backend) that writes the fitted primitives as analytic STEP faces.

```python
import meshioplusplus as mp

skin = mp.repair(mp.read("part.stl"))
seg = mp.segment(skin, curvature_tolerance=1.0, min_patch_cells=20)
fit = mp.fit_primitives(seg, tolerance=2e-3)

solid, report = mp.fit_cad(fit, "part.step", return_report=True)
report["valid"], report["num_cylinder"], report["num_demoted"]
report["volume_ratio"], report["max_deviation"]   # against the skin
report["backend"]                                 # "occt" or "native"

mp.fit_cad(mp.read("bracket.stl"), "bracket.step")  # segments and fits it first
```

![The turned part of skin to CAD as a CAD solid: two planes, a cylinder, a torus fillet, a cone and a spherical dome, re-meshed and coloured by face type](/images/fit_cad_solid.png)

Unlike the faceted [`step`](/formats/step) writer, whose solid is one flat face per triangle, the STEP file `fit_cad` writes carries `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`, `SPHERICAL_SURFACE` and `TOROIDAL_SURFACE` faces a CAD tool can select, dimension and edit. A skin without `segment:patch` is segmented first with [`method="primitive"`](/segment#primitive-method), which ends a patch where the plane, cylinder, cone, sphere or torus its facets lie on changes (so a rounded end comes apart from the flat it is tangent to), and fitted at the same absolute tolerance, `1e-3` of the bounding-box diagonal; a skin with patches but no `fit:patch_*` table is fitted with the default options. A CAD-exported STL, whose points lie exactly on its surfaces, comes out as one face per surface, not a patchwork. A noisy scan wants a curvature test and a `min_patch_cells` in the tens, so segment it yourself, or pass those to the [CLI](#cli) with `--segment-method angle`.

## How the solid is built

- **Faces.** A patch whose fit was accepted, and whose type is in `types`, becomes one face on its fitted surface. A patch no primitive explains is tried, in order, as a [rolling-ball fillet](#fillets-blends-and-threads) between two analytic faces (`fillets`), as a [helical thread face](#fillets-blends-and-threads) (`threads`; then [of varying pitch or on a taper](#chamfers-conic-blends-vertex-blends-and-tapered-threads)), as a [chamfer or conic blend](#chamfers-conic-blends-vertex-blends-and-tapered-threads) between two analytic faces (`fillets`) and, with `freeform` (the default), as a [B-spline face](#freeform-faces); once the freeform faces exist, a patch between two fitted faces, one of them freeform, is tried as a [rolling-ball blend](#fillets-blends-and-threads) against them, and a corner patch among three or more fitted faces as a [vertex blend](#chamfers-conic-blends-vertex-blends-and-tapered-threads) (`fillets`). Every other surface cell (a type left out, a freeform fit that failed, a cell with no patch) becomes planar faces, each edge-connected group of coplanar facets one face with its holes as inner loops.
- **Edges.** Where two patches meet, the mesh edges between them form chains from one junction (a point where three or more patches meet) to the next; a closed chain is split in two. Between two analytic faces the edge is the branch of the two surfaces' exact intersection that lies within the tolerance of the chain, when `exact_edges` is set and one does (a plane meeting a cylinder gives a circle). Otherwise, notably where a fillet meets its wall tangentially and no intersection is stable, the edge is a B-spline through the chain, smoothed and moved onto both surfaces, kept only if it stays within the tolerance of both. Beside a facet, the edge is the mesh's own straight edge.
- **Vertices.** A junction among analytic faces only is moved onto all of their surfaces; one beside a facet keeps the mesh position.
- **Solid.** Each face gets parameter-space curves that are continuous around its loops (so a boundary running along a cylinder's seam or a torus's equator does not fall apart), OpenCASCADE adds the seams and repairs the faces, sews them and closes each shell into a solid, one per edge-connected piece of the skin. A patch with no boundary at all, a whole sphere or torus, is its surface's natural face.

## Refining the patches first

With `refine_segmentation=True` (C++ `mRefineSegmentation`, v17.11.0; off by default, so the output is what it was before the option existed) the patches are refined by [`refine_patches`](/refine_patches) before the solid is built, after any segmentation and fit `fit_cad` runs itself: a variational shape approximation (Cohen-Steiner, Alliez and Desbrun 2004, with the primitive proxies of Wu and Kobbelt 2005) moves patch boundaries to the surfaces, splits a primitive out of a patch no primitive explains by seeded RANSAC (Schnabel, Wahl and Klein 2007) and merges neighbours one primitive explains, so a fillet strip an angle segmentation gave to the plane beside it, or a cylinder cut in two at a smooth edge, becomes the analytic face it is. The refinement's options are `max_refine_iterations` (10), `refine_epsilon` (`1e-3`, the least gain in explained area worth another iteration), `normal_weight` (1), `allow_merge` and `allow_split` (on), `ransac_seed` (a fixed default, so the output is reproducible) and `refine_tolerance` with `relative_tolerance` (`1e-3` of the bounding-box diagonal). The report then holds `refined`, `refine_iterations`, `num_refine_splits`, `num_refine_merges` and the area share primitives explain before and after (`initial_explained_area_fraction`, `explained_area_fraction`), and the returned mesh carries the per-patch `refine:*` field data. On an angle-segmented turned part the refinement raises the explained area from 38% to 100% and every surface becomes an analytic face of a valid solid.

The refinement knows primitives only: a patch it finds no primitive for within its tolerance stays unexplained here (tried as a freeform face) even when the fit table, by `fit_primitives`' looser measure, names one, and a blend or a corner the [feature fits](#fillets-blends-and-threads) explain whole may come out cut into strips. So the OpenCASCADE backend builds the solid from the refined patches and from the patches as given, and keeps the refined one unless the other is better: valid where it is not, then fewer facet faces, fewer demotions, fewer faces. `refinement_kept` in the report says which; the refinement's counts describe what it did either way, and its `refine:*` field data comes along only when it was kept. Refinement therefore never makes the OpenCASCADE result worse by that measure, and costs a second build. The native backend, which has no feature fits, always keeps the refinement.

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

### General freeform faces

With `freeform_general=True` (C++ `mFreeformGeneral`, v17.11.0; off by default) a patch the plane and unrolled-cylinder fit leaves unexplained gets a second, general attempt before it is split or faceted, so every face the earlier fit built stays as it was:

1. **A disc parametrisation.** A patch with one boundary loop is mapped to the plane by least-squares conformal maps (Lévy, Petitjean, Ray and Maillot, ACM TOG 2002, [doi:10.1145/566654.566590](https://doi.org/10.1145/566654.566590)), its two farthest boundary vertices pinned, or, if that flips a triangle, by mean-value coordinates with the boundary on a circle (Floater, CAGD 2003), which is bijective for a disc; the map is turned so its bounding box is the smallest. Both are solved by OpenCASCADE-free sparse iterative solvers whose results do not depend on the thread count.
2. **The fit.** A cubic B-spline on uniform knots over that box is fitted by least squares to the vertices and, at a quarter of their weight, the facet centroids, with a thin-plate term near zero where the samples are and strong where none are, so the surface runs on past the patch as the least-bending extension. Each iteration then corrects the samples' parameters by projecting them onto the surface (Hoschek's parameter correction) and fits again; when that gains less than 10% it takes a damped tangent-distance step (Wang, Pottmann and Liu, ACM TOG 2006, [doi:10.1145/1138450.1138453](https://doi.org/10.1145/1138450.1138453)), and when that stalls too it inserts a knot at the middle of every knot span whose samples are still too far, as long as the span holds enough samples. It stops at the target (`freeform_tolerance`, a quarter of the tolerance by default, at the vertices; the tolerance at the centroids), at `max_freeform_poles` per direction (40) or after `max_freeform_iterations` (30), keeping the best state it reached.
3. **The face.** OpenCASCADE gets only the poles and knots (`Geom_BSplineSurface`, never `GeomAPI_PointsToBSplineSurface`), and the face is accepted by the test every freeform face passes, half the tolerance at a vertex and the tolerance at a facet centre, then trimmed, sewn, checked and, if it fails, [demoted](#checks-and-demotion) like any other.

A patch that is no disc, or that the fit cannot explain, is halved along its normals' spread, up to `max_freeform_splits` times (2), and each half tried in turn; whatever is left is facets. `freeform_smoothing` (`1e-12`, relative to the sample count) is the thin-plate weight over the spans that hold samples. The report counts `num_freeform_general` (of the B-spline faces, those this path built), `num_freeform_refinements` (their knot insertion rounds) and `avg_freeform_poles`. A closed wavy band whose halves the grid fit cannot reach at four grid points a side comes back as two general B-spline faces instead of nearly six thousand facets.


## Freeform modes

`freeform_mode` (C++ `mFreeformMode`, v17.12.0) chooses how a patch no primitive, fillet, blend, chamfer or thread explains becomes a face. `"legacy"`, the default, is everything above, byte for byte. The other modes build on the [general freeform faces](#general-freeform-faces) (they imply `freeform_general` and `extend_surfaces`), climb an escalation ladder before a patch is demoted, and give every demotion a reason:

| Mode | What it adds |
| --- | --- |
| `"adaptive"` | each freeform patch as one adaptive B-spline face, or pieces of one; with the native backend too |
| `"network"` | the [patch network](#the-patch-network): a quad layout of each freeform region, fitted as B-spline patches that share their edge curves exactly, G1 across them; a patch it cannot fit falls back to an adaptive face (native backend; with OpenCASCADE, adaptive faces for now) |
| `"auto"` | adaptive faces first, the network for what they leave |

### The adaptive face

1. **Discs first.** A patch that is not a disc is cut along its topology before anything is halved (`detail/cut_graph`): an annulus by two shortest paths between its loops, a sphere with more holes by a cycle of shortest paths through every loop, a closed genus-0 patch along the median level of the geodesic distance difference to its two farthest vertices. A patch of higher genus is halved across its longest extent.
2. **The parametrisation ladder** (`freeform_parameterization`: `"auto"`, `"projection"`, `"lscm"`, `"harmonic"`): the best-fit plane, then the cylinder unrolled around the axis the normals are most perpendicular to, each kept only when it folds no triangle; then least-squares conformal maps; then mean-value coordinates with the boundary on a circle. Every map is turned so its bounding box, the spline's domain, is the smallest.
3. **The fit**: the general fit's cubic B-spline, with the poles of a face capped at `freeform_max_control_points` (default 400, valid 16 to 1048576) on top of `max_freeform_poles` per direction, and a domain five per cent wider than the patch's, so the face runs on past its trim as the least-bending extension.
4. **The rungs**: the given budgets, then twice the poles per direction and twice the iterations (refine knots), then, under `"auto"`, conformal maps and mean-value coordinates explicitly (re-parametrise), then the topology cut and the halving (split). `escalation_steps` counts the rungs a patch climbed.
5. **Acceptance.** Every vertex within half the tolerance and every facet centre within the tolerance, by projection onto the surface, and, with `validate_by_surface_sampling` (the default), the surface's Greville points and a grid per knot span inside the trim within the tolerance of the patch. A clean patch is never held looser. Only on a skin that shows noise (each vertex's residual from the quadric through its two-ring, above a twentieth of the tolerance), whose strict fit fails with residuals that do not correlate across the patch's edges, does a rung try again at `min(4.5 sigma, 3 tol)`: a scanned skin is no more accurate than its noise, as a coarse tessellation is no more accurate than its chordal sag, and a structured misfit is shape, not noise.

The native backend writes the accepted faces as `B_SPLINE_SURFACE_WITH_KNOTS` faces bounded by polylines on the skin's points, so `adaptive` needs no OpenCASCADE; `brep_checked` is then false. With OpenCASCADE the faces are trimmed, sewn and checked like any other, and a failing round demotes the freeform faces involved before an analytic neighbour, so no mode leaves a face `legacy` would have kept.

### Report and records

Outside `legacy` the report adds:

| Key | Meaning |
| --- | --- |
| `freeform_mode` | the mode that ran |
| `brep_checked` | OpenCASCADE's B-rep check ran (false for the native backend) |
| `facet_area_fraction` | the share of the skin's area left as facets, the measure the modes lower |
| `num_freeform_adaptive` | adaptive B-spline faces |
| `num_network_patches`, `num_network_edges`, `num_extraordinary_vertices`, `max_g1_angle_deg` | the patch network's patches, shared edge curves, vertices where other than four patches meet, and the largest angle between neighbouring normals along a shared edge |
| `num_resegmented`, `num_escalations` | regions re-segmented on the ladder, and rungs climbed in all |
| `max_freeform_distance`, `rms_freeform_distance` | the largest and the RMS vertex distance from a freeform face |
| `faces` | one record per patch tried as a surface |

A record holds `patch`, `face` (the B-rep face, -1 when demoted), `kind` (the `fit_cad:type`), `param_method` (`none`, `plane`, `cylinder`, `lscm`, `mean_value`, `square`), `demotion_reason` (`none`, `param_not_injective`, `hausdorff_exceeded`, `trim_failed`, `sew_failed`, `brep_check_failed`, `too_small`, `non_disc`, `kernel_mesh_coarse`, `budget_exceeded` or `other`), `num_ctrl`, `iterations`, `refinements`, `escalation_steps`, `max_dist`, `rms_dist`, `reverse_max_dist` and `max_g1_angle_deg`. The same table is the `fit_cad:face_table` field data, an `(n, 13)` float64 array with its columns in that order and the enums as their integer values. In `legacy` the C++ `FitCadResult` carries all of it, and the Python report and the field data are as before.

### The faceted area

`facet_area_fraction` on the freeform fixtures (`tests/cpp/fit_cad_fixtures.hpp`, tolerance `2e-3` of the diagonal; `MESHIOPLUSPLUS_FIT_CAD_KPI=1` runs the full set):

| Fixture | Native, legacy | Native, adaptive | Native, network | Native, auto | OpenCASCADE 7.6.3, legacy | OpenCASCADE 7.6.3, adaptive |
| --- | --- | --- | --- | --- | --- | --- |
| monkey saddle | 0.182 | 0 | 0 (16 network patches) | 0 | 0.661 | 0.490 |
| car hood | 0.375 | 0 | 0.005 (25 network patches, 10 adaptive faces) | 0 | 0.375 | 0.375 |
| Catmull-Clark valence 3 | 1 | 0 | 0 (40 network patches, 8 extraordinary vertices) | 0 | 1 | 1 |

The native backend writes every face the fit accepts; OpenCASCADE must also trim, sew and check it, and the hood's top, segmented together with a side across a 25-degree crease, still fails its corners there.

### The patch network

`network` (and `auto`, for what the adaptive faces leave) groups the unexplained patches into edge-connected regions and builds each as one network, in the native backend:

1. **The cross field and the layout** (`detail/cross_field`, `detail/quad_layout`). The smoothest cross field of the region (Knöppel, Crane, Pinkall and Schröder, "Globally optimal direction fields", SIGGRAPH 2013) is aligned with its boundary and its creases, and its singularities' indices sum to the Euler characteristic. The layout cuts along the creases (`layout_snap_to_features`) and the topology (discs first), splits a piece at its one singularity along the separatrices traced from it, turns a three-cornered piece into three quads round its centre and a disc without real corners into an O-grid, and refines to `layout_target_patch_size` (relative to the diagonal, default 0.1; at least about 48 triangles a patch) in strips that conform across shared sides. Valences stay within `layout_max_valence_deviation` (2) of four. The cube cage's limit surface becomes 24 patches with its eight valence-3 vertices and no T-junction; regions of genus one and two get valid layouts, with T-junctions where their topology cuts meet.
2. **Shared edges** (`detail/patch_network`, `detail/bspline_curve`). Every layout chain is fitted once as a cubic B-spline. A crease or the region's boundary is fitted through its vertices; a cut of the layout's own, a zigzag of mesh edges, is smoothed and moved onto free fits of the patches beside it, so it lies on the smooth surface, not on the facets. Each patch's four boundary rows of poles are exact pieces of those curves (knot insertion), so neighbours meet exactly (G0) whatever their interiors.
3. **Interiors and G1.** The interior poles are fitted by least squares to the patch's vertices, facet centres and edge middles on its unit-square map, more strongly smoothed where little data holds them. Knots are inserted at the worst spans until each patch passes. G1 sweeps then pull each patch's cross-boundary derivative into its neighbours' tangent planes (block Gauss-Seidel, the weight doubling to eight times while the largest angle exceeds `freeform_g1_tolerance_deg`), and a patch a sweep pushes out of the tolerance keeps its accurate surface.
4. **Extraordinary vertices** (`extraordinary_vertex`). `g0` leaves them alone (the G1 measure skips their first span); `g1` puts the second pole of every edge curve leaving one, and every corner's twist pole, in a common tangent plane (after Peters, "Patching Catmull-Clark meshes", SIGGRAPH 2000); `g1x` (the default) also smooths the poles round it; `g1xx` adds a double knot in the first span towards it (after Loop and Schaefer, ACM TOG 2008).
5. **Acceptance.** A patch is held to the tolerance against the network (a vertex just across a smooth edge may lie on the neighbour): every vertex within half of it and every facet centre within it, no fold (the surface's normal on a grid agrees with the skin's), and the surface's samples within it of its own and its neighbours' triangles. A patch that fails becomes an adaptive face or, failing that, facets with its reason.

The native STEP writer gives each accepted patch a `B_SPLINE_SURFACE_WITH_KNOTS` face and writes every layout edge between two network faces once, as the cubic `B_SPLINE_CURVE_WITH_KNOTS` both faces took their boundary from (its end poles on the edge's vertices), so the file's faces meet on one shared curve and not on a polyline; an edge between a network face and any other face stays a polyline on skin points. The file is reproducible byte for byte.

On the cube cage with one more subdivision level, every patch passes in `g0`, `g1` and `g1x` (largest G1 angle 1.3 degrees away from the extraordinary vertices, 11.5 degrees with them), and 17 of 24 in `g1xx`. On the monkey saddle the network is 16 patches whose largest G1 angle is 1.2 degrees. A torus's eight patches from its topology cuts are too long and twisted to fit and fall back to adaptive faces. `resegment_rounds` (2) and `local_refine` (`"knots"`, `"thb"`) are validated and reach every surface, for the rungs still to come: today no region is re-segmented, and `"thb"` warns and refines by knot insertion.

### Provenance

The algorithms are reimplemented from the papers cited, in MIT-licensed code; no GPL code (QuadWild, libigl's GPL parts, CGAL's GPL packages) is used or ported, and nothing is vendored. Instant Meshes (BSD-3) and QuadriFlow (MIT) were read as references only.

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

## Chamfers, conic blends, vertex blends and tapered threads

Since v17.2.0 (roadmap 6.3) four more kinds of patch get their own faces. Each is tried only where every fit above has failed for that patch, so a skin they do not concern builds the same solid, to the byte, as before; their fits are in `detail/feature_fit` too.

- **Chamfers** (`fillets`, type 10, `num_chamfer`) and **conic blends** (`fillets`, type 7, `num_conic_blend`). A patch between two analytic faces (the two neighbours sharing the most boundary, together at least 60% of it) that no rolling ball and no thread explains is tried as a blend of the edge between them with a cross-section other than the rolling ball's circle. Its rails are the patch's boundary edges across from each neighbour; its cross-sections lie in the planes normal to the edge the two surfaces would meet along (`n₁ × n₂`, through each point of the first rail), from the first rail to where that plane cuts the second. When every interior vertex lies within 0.4 of the tolerance of its section's chord, the patch is a **chamfer**, of constant or varying width (a chamfer between two planes of varying width is a twisted ruled surface, a constant one a plane, and a chamfer coaxial with a cylinder a cone: those two stay the primitives `fit_primitives` fits). Otherwise each section is taken as a rational quadratic from rail to rail through the apex where the two surfaces' tangent lines meet, the vertex's weight `w = b₁ / (2 √(b₀ b₂))` read off its barycentric coordinates in that triangle, and a weight law of the lowest degree (Chebyshev, Fourier round a closed rail) whose conics pass within 0.4 of the tolerance of every vertex makes it a **conic blend**: an elliptic, parabolic or hyperbolic section tangent to both faces, of legs that may differ and change along it. The face is a B-spline interpolating the sections, each widened by 15% and carried a tenth of the rail past an open end.
- **Tapered and variable-pitch threads** (`threads`, type 8, `num_tapered_thread`, `num_variable_pitch_thread`). A straight flank beside a cone is the same ruled helicoid as a flank of a straight thread (with an apparent pitch), so it is built as before and reported tapered. A flank whose straight helicoid fails is refitted with one more column, `z = c + k θ + q (θ − θ₀)² + t ρ`, a pitch changing linearly along it, and kept when the pitch changes by more than the tolerance across the patch. A patch no constant screw explains is refitted as a helical surface of any profile on a generalised screw, its advance `g(θ) = k θ + q (θ − θ₀)²` and its profile moved out by `τ g(θ)` as it advances, a thread cut on a cone whose radius grows by `τ` per unit length: `τ` starts at a neighbouring cone's half-angle and is refined with the axis, the advance and `q` by the same Levenberg-Marquardt profile fit. A trapezoidal pipe thread then comes back as cones and two tapered flanks, a thread whose pitch grows as cylinders and two variable-pitch flanks, and a rounded one on a taper as four helical faces, the rounded root and crest tapered.
- **Vertex blends** (`fillets`, type 9, `num_vertex_blend`). After the freeform and blend stages, a patch that is unexplained or freeform, whose boundary is one closed loop around three or more fitted faces (each met along one run) and creased nowhere (no facet across it turns more than 20°), is tried as the corner where blends meet: a plate surface (`GeomPlate_BuildPlateSurface`) through its vertices, tangent (G1) to each neighbour's surface along the loop, approximated by a B-spline. It is kept when it passes the freeform face's test and its normal stays within 3° of every neighbour's along the loop. This covers a setback corner, three blends of different radii or sections meeting at one vertex; a corner of three equal rolling balls is a sphere, which `fit_primitives` already fits.

Segmentation is the condition: a conic blend is tangent to both its faces and its curvature varies across it, so no curvature tolerance of [`segment`](/segment) isolates it alone; label such a skin's patches yourself (`segment:patch`, from the CAD faces it was tessellated from, say), as the fixtures behind the tests do. A chamfer meets its faces at creases and segments by angle.

![A block whose edge is cut by a chamfer of growing width: one chamfer face, a twisted ruled surface, between two planes](/images/fit_cad_chamfer.png)

![A rounded box scaled unequally: twelve conic edge blends and eight vertex blends at its corners, between six planes](/images/fit_cad_vertex_blend.png)

![A pipe thread: a trapezoidal thread cut on a cone, its root and crest cones and both flanks tapered helical faces](/images/fit_cad_tapered_thread.png)

## Checks and demotion

A face is not kept when OpenCASCADE's B-rep check rejects it, when sewing leaves an edge of it free, when its boundary cannot be built, or when it strays: the solid is re-meshed and a face with a point further than `tolerance` from the skin, or a patch with a skin point further than that from the solid, fails. A failed patch is **demoted** to facets and the solid is rebuilt. A failure is blamed on the faces the B-rep check names or that a free edge after sewing touches (a face the kernel cannot map back is matched to the patch of its nearest facet), never on every patch at once; after `max_demotion_rounds` rounds (default 8) each further round demotes the blamed patches together with their neighbours, so the loop ends while the patches nothing blames stay analytic, and an all-facet solid is always valid for a closed skin. A group of coplanar facets whose merged face fails is split back into single facets. A face's distance allowance is `tolerance` or, when larger, the skin's own chordal sag over its patch (1.05 times the largest distance of a facet centre from the fitted surface), so a coarse tessellation of a small bore is not demoted for its own chords. So a feature the fitter cannot explain, or a face the kernel cannot trim, comes back faceted, never silently wrong, and `num_demoted` counts the features dropped. `MESHIOPLUSPLUS_LOG_LEVEL=debug` logs why each patch was demoted.

### Edges, sewing and validation options

Four options, all off by default so the output is what it was without them (OpenCASCADE backend, v17.11.0):

- `sew_tolerance` (`0`): the tolerance the faces are sewn at; `0` sews at `tolerance` itself. The edges are shared already, so sewing only closes what the face repair re-made, and a tighter value makes a gap the repair left fail the check instead of being bridged.
- `extend_surfaces`: each freeform B-spline surface gets an extended copy, past its patch on each open side by a tenth of the patch's size (`GeomLib::ExtendSurfByLength`, tangent-continuous). Where the intersection of two faces' surfaces has no branch along their shared chain, the copies are intersected instead, so an edge that runs just outside the fitted area is still an exact curve; the face keeps its own surface, and the edge passes the same checks. `num_extended_edges` counts the edges found that way. An extended branch is kept only where it stays within `tolerance` of both faces' own surfaces, since the faces are trimmed on those. It is opt-in for a reason: on the generated threaded rod with `freeform_general` it finds 9 exact edges but leaves 21 facet faces and 2 demotions against 20 and 1 without it, at two and a half times the time, and swapping the extended surfaces in for the faces themselves (tried first) did worse still.
- `snap_junctions`: a B-rep vertex where analytic faces meet is placed by a Gauss-Newton least-squares solve on the surfaces' tangent planes (OpenCASCADE-free, `detail::junction_closest`): each step minimises the summed squared distances to the tangent planes at the current point's feet, with a small Levenberg damping that leaves a direction no plane constrains where it was. At a transversal corner it lands in a few steps, where the averaged alternating projection creeps. The point is accepted on the same test as before (within `tolerance` of every surface and of the mesh vertex), and the projection is the fallback; `num_snapped_junctions` counts the vertices it placed.
- `face_deviation`: the two-sided distance between each face and the skin is measured, from the face's re-mesh (vertices and triangle centres) to the skin and from the skin (vertices and facet centres, each charged to the face of the nearest re-mesh triangle) to the solid. The field data `fit_cad:face_hausdorff` holds per face (in `fit_cad:face` order) both maxima, an `(n, 2)` float64 array, and `fit_cad:face_rms` the root mean square of both sides' distances; the report adds the global `hausdorff_max` and `hausdorff_rms`. `max_deviation` keeps its meaning.

Every run also reports `closed`: every shell's edges border two faces each (a seam counts for its face). The native backend reports its writer's closed-loop check.

**Merging.** With `merge_faces` (the default), once the solid is built, faces on one surface are merged and collinear edges joined (`ShapeUpgrade_UnifySameDomain`); the result is kept only when it stays valid with the same volume (to `1e-9`), and the `fit_cad:face` labels follow the merged faces.

`tolerance` defaults to three times the largest `fit:patch_tolerance` of an accepted patch, so it follows the scale `fit_primitives` worked at; with no accepted patch it is `1e-3` of the bounding-box diagonal. The re-meshing uses a chordal deflection of a quarter of the tolerance unless `mesh_deflection` says otherwise.

## Output and report

The return value is the solid re-meshed by OpenCASCADE as one `triangle` block, with points welded where faces share them exactly, and Int64 cell data `fit_cad:face` (the B-rep face), `fit_cad:type` (0 a facet, else the `fit:type` value: 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus; 6 a freeform B-spline, 7 a rolling-ball fillet, blend or conic blend, 8 a helical thread face, 9 a vertex blend, 10 a chamfer) and `fit_cad:patch` (the `segment:patch` label, -1 for a cell with no patch). With `return_report=True` the report holds:

| Key | Meaning |
| --- | --- |
| `backend` | `"occt"` or `"native"`, the engine that built the solid |
| `valid` | the solid passes OpenCASCADE's B-rep check (native: the writer's closed-shell check) |
| `tolerance` | the tolerance applied |
| `volume`, `skin_volume`, `volume_ratio` | the solid's volume, the skin's enclosed volume and their ratio |
| `max_deviation` | the symmetric Hausdorff distance between the skin and the re-meshed solid |
| `num_faces`, `num_plane` … `num_torus`, `num_bspline`, `num_fillet`, `num_helical`, `num_faceted_faces` | B-rep faces in all, by type (freeform B-splines, fillets and thread flanks included), and planar facets |
| `num_variable_fillet`, `num_blend` | of the fillet faces, those whose radius varies along them, and the blends beside a freeform face |
| `num_rounded_thread` | of the helical faces, those whose profile is not a straight line (a rounded root or crest) |
| `num_vertex_blend`, `num_chamfer` | vertex blends and chamfers (v17.2.0) |
| `num_conic_blend` | of the fillet faces, the conic blends (v17.2.0) |
| `num_tapered_thread`, `num_variable_pitch_thread` | of the helical faces, those cut on a taper and those whose pitch varies (v17.2.0) |
| `num_split_patches` | pieces that freeform patches folding over every flattening were split into |
| `num_patches`, `num_demoted` | patches in the fit table, and accepted ones demoted to facets |
| `num_exact_edges`, `num_approx_edges` | edges between analytic faces built as exact intersections and as B-splines |
| `num_solids` | solids written, one per piece of the skin |
| `refined`, `refinement_kept`, `refine_iterations`, `num_refine_splits`, `num_refine_merges`, `initial_explained_area_fraction`, `explained_area_fraction` | with `refine_segmentation` only: the [refinement](#refining-the-patches-first), and whether the solid was built from it (v17.11.0) |
| `num_freeform_general`, `num_freeform_refinements`, `avg_freeform_poles` | with `freeform_general` only: the [general freeform faces](#general-freeform-faces) (v17.11.0) |
| `closed` | every shell is closed: no edge borders only one face (v17.11.0) |
| `num_extended_edges` | with `extend_surfaces` only: exact edges found on extended surfaces (v17.11.0) |
| `num_snapped_junctions` | with `snap_junctions` only: vertices placed by the least-squares junction solve (v17.11.0) |
| `hausdorff_max`, `hausdorff_rms` | with `face_deviation` only: the two-sided distance between the solid and the skin over both sides' samples, with the per-face `fit_cad:face_hausdorff` and `fit_cad:face_rms` field data (v17.11.0) |

The STEP file is AP214 in millimetres. Its header timestamp is pinned and its product is named `meshio++ fit_cad`, so the same input writes the same bytes; the bytes do depend on the OpenCASCADE version. The work is serial, so nothing depends on the thread count.

## The solid as NURBS

`solid.info["brep"]` (C++ `FitCadResult::mBrep`, a `BrepInfo`) holds the solid with every face converted to NURBS (`BRepBuilderAPI_NurbsConvert`): the surfaces, their trimming loops and the edges with the faces they bound, as a [STEP read](/formats/step#reading) gives them. Every face of the solid is there, trimmed, and every trimming curve lies in its surface's parameter range. When the exact conversion does not deliver that (OpenCASCADE 7.9 throws on a trimmed cylinder or one wound several turns round, leaves a face that already was a B-spline without a surface, and can stop a converted surface short of its trims), the solid is kept as it is and each surface and edge that is not a B-spline is approximated by one in its own parameters, within `1e-6` of the bounding-box diagonal, so the trimming curves stay valid; planes and lines stay exact. Writing it as [Kratos CAD JSON](/formats/cad_json) (`mp.write("part.cad.json", solid)`, `to-cad --cad-json`) makes a scanned part Kratos IGA geometry.

## Getting OpenCASCADE

`fit_cad` assumes OpenCASCADE for the best result, and **the PyPI wheels and the release CLI binaries include it**: `pip install meshioplusplus` is all it takes. [Open CASCADE Technology](https://github.com/Open-Cascade-SAS/OCCT) is licensed under the GNU LGPL 2.1 with the additional exception in its `OCCT_LGPL_EXCEPTION.txt`; those binaries link the pinned OpenCASCADE 7.6.3 statically, so they are `MIT AND (LGPL-2.1-only WITH OCCT-exception-1.0)`, with [`THIRD_PARTY_NOTICES.md`](https://github.com/loumalouomega/meshioplusplus/blob/master/THIRD_PARTY_NOTICES.md), OpenCASCADE's licence texts and the relink recipe shipped beside them. The meshio++ sources stay MIT. The [full WASM artifact](/wasm#lean-and-full-builds) links the same commit; the lean one does not.

A build from source turns `MESHIOPLUSPLUS_WITH_OCCT` on by default and links an OpenCASCADE 7.6 or later when it finds one through its CMake config (`OpenCASCADE_ROOT`, `OpenCASCADE_DIR` or `CMAKE_PREFIX_PATH` for an install outside the default prefixes): a system package, linked dynamically, or the static build `build/build-occt.sh PREFIX` makes from the pinned commit. `-DMESHIOPLUSPLUS_REQUIRE_OCCT=ON` makes a missing kernel a configure error; `-DMESHIOPLUSPLUS_WITH_OCCT=OFF` (`build/configure.sh --without-occt`) leaves it out.

```bash
sudo apt-get install libocct-foundation-dev libocct-modeling-data-dev \
    libocct-modeling-algorithms-dev libocct-data-exchange-dev libtbb-dev   # or:
build/build-occt.sh ~/occt && export OpenCASCADE_ROOT=~/occt
pip install --no-binary meshioplusplus meshioplusplus                       # finds it
```

Debian's packages export toolkit targets that link TBB, hence `libtbb-dev`. OpenCASCADE 7.6 or later is supported (7.6 is Ubuntu 24.04's package and the binaries carry 7.6.3); CI runs the `fit_cad` suites on Ubuntu's 7.6 and on 7.9.3, which `MIO_OCCT_PIN=7.9 build/build-occt.sh PREFIX` builds from the same recipe. The STEP bytes depend on the kernel's version, so the golden digests that pin `fit_cad`'s output (`tests/cpp/test_fit_cad_goldens.cpp`) are kept per version; `fit_cad_occt_version()` in C++ and `_core.__occt_version__` in Python name the version a build links (`""` without one).

Without the kernel, `fit_cad` falls back to the [native backend](#the-native-backend) and says so once (a warning naming what is lost); `backend="occt"` raises `NotImplementedError` (C++ `Unsupported`) naming `-DMESHIOPLUSPLUS_WITH_OCCT=ON`. There is no pure-Python fallback of either: without the compiled core `fit_cad` raises. `meshioplusplus._fit_cad.has_occt()`, `fit_cad_has_occt()` in C++, `mio_has_occt()` in C and the language surfaces below say in advance whether a build has it.

## The native backend

`backend="native"` (C++ `FitCadBackend::Native`), and `"auto"` (the default) in a build that found no OpenCASCADE, needs no CAD kernel: it segments and fits the skin as above when it has no patches, then writes the [`step`](/formats/step) writer's `faces="fitted"` B-rep. Every accepted plane, sphere, cylinder, cone and torus patch whose type is in `types` is one analytic face on its fitted surface, bounded by `LINE`, `CIRCLE` and polyline edges, its corners moved onto the surfaces they lie on; every other cell is part of a coplanar planar face. A patch with a point further than `tolerance` from its surface, or whose boundary does not close into loops, is written as planar faces and counted in `num_demoted`. On a CAD-exported bracket it gives the same twelve planes and four cylinders as the OpenCASCADE backend, which OpenCASCADE reads back as one valid solid of the exact volume.

In the [freeform modes](#freeform-modes) it also writes adaptive B-spline faces (`B_SPLINE_SURFACE_WITH_KNOTS`). What it does not do: legacy freeform, fillet, blend and thread faces (`freeform`, `fillets`, `threads`, `exact_edges`, `mesh_deflection`, `max_demotion_rounds` and `merge_faces` are ignored), sewing and the B-rep check, and re-meshing. The returned mesh is the input skin with the `fit_cad:*` cell data, `info["brep"]` is empty, `volume` and `volume_ratio` are NaN, `max_deviation` is the largest distance of a skin point from its analytic face, and `valid` says the shell was written closed, every edge used once each way. Its STEP bytes do not depend on any library version.

## CLI

```sh
meshioplusplus to-cad IN OUT.step [--curvature-tolerance C] [--min-patch-cells N] \
    [--fit-tolerance T] [--tolerance T] [--no-exact-edges] [--types plane,cylinder,...] \
    [--mesh-deflection D] [--max-demotion-rounds N] [--no-freeform] [--freeform-degree P] \
    [--freeform-max-control N] [--no-freeform-split] [--no-fillets] [--no-threads] \
    [--segment-method primitive|angle] [--segment-tolerance T] [--no-merge-faces] \
    [--backend auto|occt|native] [--refine] [--max-refine-iterations N] [--refine-epsilon E] \
    [--normal-weight W] [--no-refine-merge] [--no-refine-split] [--ransac-seed S] \
    [--refine-tolerance T] [--freeform-mode legacy|adaptive|network|auto] \
    [--freeform-parameterization auto|projection|lscm|harmonic] [--g1-tolerance DEG] \
    [--extraordinary-vertex g0|g1|g1x|g1xx] [--local-refine knots|thb] \
    [--freeform-max-control-points N] [--layout-target-patch-size S] \
    [--layout-max-valence-deviation N] [--resegment-rounds N] [--no-layout-snap-to-features] \
    [--no-validate-by-surface-sampling] [--mesh-out FILE] [--cad-json FILE] [--json] [--quiet]
```

A skin without `segment:patch` is segmented by `--segment-method` (default `primitive`, at `--segment-tolerance` of the diagonal, default `1e-3`; `angle` uses `--curvature-tolerance` and `--min-patch-cells`) and fitted, at the segmentation's own tolerance as a length after the primitive method, or with `--fit-tolerance` (relative, as in `fit-primitives`) when given, before the solid is built. `--backend` picks the engine and `--no-merge-faces` keeps faces on one surface apart; the summary and `--json` name the backend. `--mesh-out` also writes the re-meshed solid, `--cad-json` the solid as Kratos CAD JSON; `--no-freeform`, `--freeform-degree`, `--freeform-max-control` and `--no-freeform-split` are the [freeform](#freeform-faces) options, `--no-fillets` and `--no-threads` the [feature](#fillets-blends-and-threads) ones. Both CLIs exit with 2 when the solid fails the B-rep check, as `check` does for a failed gate.

`--refine` (v17.11.0) first refines the patches by variational shape approximation (`refine_patches`, `refine_segmentation=True`), so more of the skin becomes analytic faces: `--max-refine-iterations` (default `10`) caps the iterations, `--refine-epsilon` (default `1e-3`) stops them when one adds less than that to the area share primitives explain, `--normal-weight` (default `1`) weighs the normal term of a cell's cost, `--no-refine-merge` and `--no-refine-split` turn off merging neighbours one primitive explains and splitting a primitive out of a freeform patch, `--ransac-seed` (default `0x6d657368696f2b2b`; the native CLI also takes `0x` hexadecimal) seeds the split, and `--refine-tolerance` (default `1e-3` of the diagonal) is how far a point may lie from its refined patch's primitive. Only with `--refine` does the summary add a `refined:` line and `--json` add `refined`, `refine_iterations`, `num_refine_splits`, `num_refine_merges`, `initial_explained_area_fraction` and `explained_area_fraction`. The `refine-patches` verb of both CLIs runs the refinement alone ([CLI reference](/cli#meshioplusplus-refine-patches)).

`--freeform-mode` and the flags after it are the [freeform modes](#freeform-modes)' options (v17.12.0); outside `legacy` the summary adds a `freeform:` line with the mode, the faceted area and the adaptive and network faces, and `--json` adds the mode's report keys and `faces`.

```sh
meshioplusplus to-cad bracket.stl bracket.step
meshioplusplus to-cad scan.stl scan.step --segment-method angle --curvature-tolerance 1 --min-patch-cells 20 --fit-tolerance 2e-3
meshioplusplus to-cad scan.stl scan.step --refine --ransac-seed 7 --json
```

## Pipeline

```json
{"Op": "FitCad", "Path": "part.step", "Tolerance": 0.0, "ExactEdges": true,
 "Types": ["plane", "cylinder"], "MeshDeflection": 0.0, "MaxDemotionRounds": 8,
 "Freeform": true, "FreeformDegree": 3, "FreeformMaxControl": 64,
 "FreeformSplit": true, "Fillets": true, "Threads": true, "MergeFaces": true,
 "Backend": "auto"}
```

The [freeform modes](#freeform-modes) take `FreeformMode`, `FreeformParameterization`, `FreeformG1ToleranceDeg`, `FreeformMaxControlPoints`, `LayoutTargetPatchSize`, `LayoutMaxValenceDeviation`, `LayoutSnapToFeatures`, `ExtraordinaryVertex`, `ResegmentRounds`, `LocalRefine` and `ValidateBySurfaceSampling`, and outside `legacy` the step also reports `BrepChecked`, `FacetAreaFraction`, `NumFreeformAdaptive`, `NumNetworkPatches`, `NumNetworkEdges`, `NumExtraordinaryVertices`, `NumResegmented`, `NumEscalations` and `MaxG1AngleDeg`.

The step writes the solid to `Path` when given, and the mesh continues as the re-meshed solid; it reports `Native` (whether the native backend built it), `Valid`, `VolumeRatio`, `MaxDeviation`, `NumFaces`, `NumBspline`, `NumFillet`, `NumHelical`, `NumVariableFillet`, `NumBlend`, `NumRoundedThread`, `NumVertexBlend`, `NumChamfer`, `NumConicBlend`, `NumTaperedThread`, `NumVariablePitchThread`, `NumSplitPatches`, `NumFacetedFaces` and `NumDemoted`. A `Tolerance` or `MeshDeflection` of 0 is automatic; an empty `Types` facets everything.

## Other languages

```c
mio_fit_cad_opts opts;
mio_fit_cad_opts_init(&opts);
opts.types = (1u << MIO_FIT_PLANE) | (1u << MIO_FIT_CYLINDER);
opts.freeform = 1;             /* the default; freeform_degree, freeform_max_control */
opts.backend = MIO_FIT_CAD_AUTO;  /* or MIO_FIT_CAD_OCCT, MIO_FIT_CAD_NATIVE */
mio_fit_cad_report report;
mio_fit_cad_report_ex report_ex;  /* v17.2.0: num_vertex_blend, ...; backend */
mio_mesh* solid = mio_fit_cad_ex(fitted, "part.step", &opts, &report, &report_ex);
```

```fortran
solid = fit%fit_cad(step_path='part.step', valid=ok, num_demoted=ndropped, used_native=native, stat=ierr)
```

```julia
r = fit_cad(fit; step_path="part.step"); r.valid, r.backend
```

```r
cad <- mio_fit_cad(fit, step_path = "part.step"); cad$backend
```

The freeform modes are `mio_fit_cad3` (v17.12.0) with `mio_fit_cad_opts_ex2` (`freeform_mode`, `freeform_parameterization`, `extraordinary_vertex`, `local_refine` as the `MIO_FREEFORM_*`, `MIO_PARAM_*`, `MIO_EV_*` and `MIO_LOCAL_REFINE_*` constants, and the layout options; a zero in `freeform_g1_tolerance_deg`, `freeform_max_control_points`, `layout_target_patch_size` or `layout_max_valence_deviation` reads as its default, the enums' zero is their first value and the flags' zero is off, so `mio_fit_cad_opts_ex2_init` is the way to get every default) and `mio_fit_cad_report_ex3` (the mode's counts and measures, and `num_face_records`; the records themselves are the `fit_cad:face_table` field data, `demotion_reason` as `MIO_DEMOTION_*`); `mio_fit_cad2` forwards to it with the defaults. Fortran's `fit_cad` takes them as optional arguments (the enums as strings), Julia's and R's as keyword arguments, and WASM's `fitCad` in camelCase, returning `faces` outside `legacy`.

The C report carries `num_bspline`, `num_fillet`, `num_helical`, `num_variable_fillet`, `num_blend`, `num_rounded_thread` and `num_split_patches`, and `mio_fit_cad_ex`'s extended report (`mio_fit_cad_report_ex`, v17.2.0, as `mio_fit_cad_report` has no room left) `num_vertex_blend`, `num_chamfer`, `num_conic_blend`, `num_tapered_thread` and `num_variable_pitch_thread`; Fortran's, Julia's and R's results carry all of them (`freeform`, `freeform_degree`, `freeform_max_control`, `freeform_split`, `fillets`, `threads`, `merge_faces` and `backend` are their options; a zero-filled C struct keeps `freeform_split`, `fillets`, `threads` and `merge_faces` off, `mio_fit_cad_opts_init` turns them on, and `backend` 0 is auto), and `report_ex.backend` (Julia's and R's `backend`, Fortran's `used_native`) names the engine. The MCP server's `fit_cad` tool takes the same options plus `segment_method`, `segment_tolerance`, `curvature_tolerance`, `min_patch_cells`, `fit_tolerance` and `cad_json_output_path`, and reports `backend` and `occt_available`. In the browser, `fitCad(mesh, {stepPath, mergeFaces, backend, ...})` runs the native backend in the lean artifact, and the [full WASM artifact](/wasm#lean-and-full-builds) links OpenCASCADE, with the same options and report in camelCase and the solid's B-rep as `info`; `backend: 'occt'` throws by name in the lean one.

## Not supported

- **Blends of a section that is no conic**, a chamfer or conic blend beside a freeform face (only analytic neighbours are tried), and a corner patch whose loop meets a neighbour twice or at a crease: they come back freeform or faceted.
- **Threads whose depth varies along them**, or whose pitch varies other than linearly, and a tapered or variable-pitch thread with no cylinder or cone beside it that is also too far from one screw for the axis to be found from the patch itself: freeform or faceted, as is a thread patch spanning less than half a turn.
- **A whole thread as one patch**: a thread surface whose profile spans the full pitch (one patch over root, flanks and crest) has no single trimmed face to be; [`segment`](/segment) it at its creases or with a curvature test first.
- **Fitting a STEP or IGES model**: those are [read](/formats/step#reading) as their own NURBS B-rep, not refitted.
- **Feature faces without OpenCASCADE**: the [native backend](#the-native-backend) writes planes, spheres, cylinders, cones and tori, and, in the [freeform modes](#freeform-modes), adaptive B-spline faces only.
