# Volume meshes on an octree — `octree_mesh`

`octree_mesh` meshes the inside of a skin on a balanced octree. The skin is any surface (`triangle`, `quad` and `polygon` blocks, read from STL, OBJ, PLY or any other format; the boundary of a volume mesh is taken), closed or slightly dirty: inside and outside come from the fast generalised winding number, so holes, self-intersections and non-manifold edges are tolerated. It is the volume-meshing sibling of [`octree`](octree.md), which builds the tree around a skin without meshing its inside; `octree` itself is unchanged.

```python
import meshioplusplus as mp

skin = mp.read("bracket.stl")
cells, report = mp.octree_mesh(skin, output="castellated", min_size=0.5, max_size=4.0, return_report=True)
mp.write("bracket_cells.vtu", cells)
```

Three outputs are planned ([plan](octree_mesh_3d_plan.md)); the castellated and tetrahedral ones are here, and `output="hex"` raises `NotImplementedError` by name until its phase lands.

| `output` | what it is |
|---|---|
| `"castellated"` | the leaves inside the skin as `hexahedron` cells, with surrogate-boundary tags, for embedded and Shifted Boundary Method solvers |
| `"tet"` (default) | conforming graded tetrahedra fitted to the skin, sharp features preserved |
| `"hex"` | all-hexahedra by the octree dual |

## The size field

Leaves are no larger than the size field `h(x) = clamp(min(h_max, h_curv, h_feat, h_gap, h_user), min_size, max_size)`, evaluated on every skin triangle as the smallest size its vertices ask for:

| option | term |
|---|---|
| `max_size` | the size away from the skin (default: an eighth of the skin's bounding-box diagonal) |
| `curvature_angle` | `h_curv = 2 sin(θ/2) / κ_max`: a leaf edge subtends at most `θ` degrees on the skin's largest principal curvature (default 15; 0 disables) |
| `feature_angle`, `feature_size` | `h_feat` along edges sharper than `feature_angle` degrees (default 30) and non-manifold edges (default: `min_size`) |
| `gap_cells` | `h_gap = t / gap_cells`, `t` the distance an inward ray from a skin vertex travels to the far wall (default 2; 0 disables); `gap_outside` also resolves the gaps between separate parts |
| `size_array` | `h_user`, a `point_data` array of target sizes on the skin |
| `min_size` | the floor (default `max_size / 32`) |

`min_size` and `max_size` only clamp the field: a sphere meshed with `min_size=0.01` is still as coarse as its curvature asks for (about 0.26 at the default 15° on a unit sphere). To refine, tighten a criterion (`curvature_angle`, `feature_size`, `size_array`) or, under uniform grading, set `surface_size`.

The root lattice has `L = ceil(log2(max_size / min_size))` levels (at most 19; `max_level` caps it) with root cells of edge `min_size · 2^L`, centred on the padded box (`padding`, `padding_relative`, or `bounds`). A leaf the skin cuts (exact triangle–box overlap) is split until it is no larger than the size its triangles ask for, every other leaf until it is no larger than `max_size`; then `grading_layers` rings of leaves around every cut leaf are kept at its level, and the tree is 2:1-balanced across faces, edges and corners. `skin_grading="uniform"` puts every cut leaf at one level instead (`surface_size`, by default the smallest size the field asks for on the skin), with two rings around it and at least three skin sizes of padding.

## Inside and outside

A point is inside where its generalised winding number exceeds `winding_threshold` (0.5). The winding number is evaluated with the tree approximation of Barill et al. (2018): far triangles are replaced by their area-weighted dipoles and first and second moments when they are more than `winding_beta` (2) times their extent away, which gives about three correct digits; values within 0.1 of the threshold are evaluated again at a quarter of the error. `winding_beta=0` sums every triangle exactly. A skin wound inside out is detected by its negative volume, its winding numbers negated, and a warning reported.

## The castellated output

The kept leaves (`fill`: `"inside"`, the default, those whose centre is inside; `"cover"`, those and every leaf the skin cuts; `"all"`) are one `hexahedron` block in Morton order, with:

| name | kind | meaning |
|---|---|---|
| `inside` | cell region | the cells whose centre is inside |
| `boundary_cells` | cell region | the cells the skin passes through |
| `surrogate_boundary` | side region | the faces of inside cells with no inside cell across them: the surrogate boundary of the Shifted Boundary Method |
| `sdf:distance` | point data | the signed distance to the skin (negative inside, 0 on it) |
| `sdf:closest_point` | point data, 3 components | the nearest point on the skin |
| `octree:hanging` | point data | 1 for a node in the interior of a face or edge of a cell |
| `octree:level`, `octree:cut` | cell data | as for [`octree`](octree.md) |

![Castellated cells of a sphere, cut open and coloured by octree level](/images/octree_mesh_castellated.png)

An inside leaf beside a finer leaf that is not inside is split until every face of the surrogate boundary is a whole face of an inside cell, so the side region describes it exactly; the surrogate faces enclose the inside cells' volume (the tests check the divergence theorem). `attach_distance=False` skips the two `sdf:*` arrays.

## The tetrahedral output

The default output is a conforming, graded tetrahedral mesh of the skin's inside, in three stages.

**The lattice.** The balanced octree is tetrahedralised as a graded body-centred cubic (BCC) lattice on the doubled integer lattice, so every vertex is exact: each leaf is a cone from its centre over fans of its face portions (a face shared with four finer leaves is four portions), and across a face between two leaves of the same size the two cones become the BCC tetrahedra around the segment between the centres. Six shapes occur, with dihedral angles between 35.26° and 125.26°; a region of equal leaves is all BCC tetrahedra (60°–90°). This is our own construction after Labelle and Shewchuk's lattice, and the tests pin its angle table.

**Stuffing.** Isosurface stuffing (Labelle and Shewchuk, SIGGRAPH 2007) cuts the lattice to the skin: a lattice edge whose ends lie on either side gets a cut point where it crosses the skin (exactly, through the winding-number tree, or by bisection where it passes through a hole); a lattice vertex closer to a cut point than `alpha_short` (edges between a leaf centre and another vertex) or `alpha_long` (the others) times the edge is *warped* onto it; then every tetrahedron the skin crosses is cut by a stencil. Quadrilaterals take the diagonal through their lower-numbered lattice vertex on lattice faces and their shorter diagonal on the skin, so both sides of every face agree. With the paper's `alpha_long = 0.24799` and `alpha_short = 0.40173`, every tetrahedron of a uniform lattice has dihedral angles between 10.7° and 164.8°; the tests check this on spheres, tori and random ellipsoids. The lattice's own boundary is never warped, so `keep_outside=True`, which also keeps the tetrahedra outside the skin as the `outside` region, fills the root lattice exactly.

![Graded tetrahedra of a sphere, cut open and coloured by the level of their leaf](/images/octree_mesh_tets.png)

**Features.** With `preserve_features` (the default), the skin's edges sharper than `feature_angle` form a graph of corners (where it branches, ends, or turns by more than `corner_angle`) and chains between them. Each corner moves the nearest boundary node within reach onto itself, exactly; each chain is traced on the boundary by a shortest path in a tube around it, whose nodes move onto the chain. Every move is refused if it would invert a tetrahedron or create a dihedral angle below 2° that was not already there; a corner or chain that cannot be recovered is counted (`num_corners_snapped`, `num_feature_chains_missed`) and reported in `warnings`. The output carries the `feature_corners` and `feature_edges` point regions.

![A cube meshed with its corners and edges kept, coloured by each tetrahedron's smallest dihedral angle](/images/octree_mesh_cube.png)

`conformity="cleaving"` meshes nested shells (a skin whose winding number reaches 2 or more, such as a sphere inside a sphere): the interface at each winding level `winding_threshold + k` is stuffed in turn through the tetrahedra between the previous interfaces, with the earlier interfaces held fixed, so the materials conform. It follows the idea of lattice cleaving (Bronson, Levine and Whitaker, IEEE TVCG 2014), one interface at a time rather than with its multi-material stencils. The materials are the `octree_mesh:material` cell data and the `material_<k>` cell regions.

| name | kind | meaning |
|---|---|---|
| `inside` | cell region | the tetrahedra inside the skin |
| `outside` | cell region | `keep_outside`: the tetrahedra outside it |
| `boundary` | side region | the boundary faces |
| `feature_corners`, `feature_edges` | point regions | the nodes on the skin's corners and feature chains |
| `material_<k>` | cell regions | cleaving: each material |
| `octree:level` | cell data | the level of the leaf a tetrahedron came from |
| `octree_mesh:shape` | cell data | the lattice shape (1–6), 0 for a tetrahedron stuffing cut |
| `octree_mesh:material` | cell data | the material, when there is more than one |
| `sdf:distance` | point data | the signed distance to the skin |

**Optimisation.** `fit_iterations` rounds (default 5) then improve the mesh: an optimal-Delaunay (ODT) relocation of the interior nodes (the pass of [`smooth`](smooth.md)), a slide of every skin node towards the mean of its skin neighbours projected back onto the skin, and the 2-3 and 3-2 flips of [`optimize_volume`](optimize_volume.md). With `dihedral_floor_guard` (the default) a node is put back when it would lower the mesh's smallest dihedral angle or raise its largest, and a round of flips is kept only when it does neither, so optimisation never makes the worst angles worse; feature nodes, the box and interfaces between materials never move, and flips are made only for one material. `num_optimized_moves` and `num_flips` count what changed.

![Histograms of the dihedral angles of a sphere's tetrahedra before and after optimisation, graded and uniform](/images/octree_mesh_dihedrals.png)

`guarantee_applies` in the report is true when the 10.7°–164.8° bound covers the whole mesh: uniform skin grading, no warp refused, and no node moved onto a feature. Otherwise the actual angles are in `min_dihedral` and `max_dihedral`; snapping to a sharp corner can lower them locally, and graded skins cut transition shapes.

## Report

`return_report=True` returns the frame (`origin`, `root_size`, `root_cells`, `max_level`, `finest_level`, `surface_level`, `min_size`, `max_size`), the tree (`num_leaves`, `leaves_per_level`), the counts (`num_points`, `num_cells`, `num_hexahedra`, `num_tetrahedra`, `num_inside`, `num_boundary_cells`, `num_surrogate_faces`, `num_hanging`, `num_cut_points`, `num_warped`, `num_warp_refused`), the features (`num_corners`, `num_corners_snapped`, `num_feature_chains`, `num_feature_chains_missed`, `num_feature_moves`), the optimisation (`num_optimized_moves`, `num_flips`), `guarantee_applies`, the quality (`min_scaled_jacobian`, `min_dihedral`, `max_dihedral`, `num_inverted`, `hausdorff_to_skin` (the two-sided distance between the mesh's boundary and the skin), `num_boundary_edges_open` and `num_non_manifold_boundary_edges` (of the boundary, both 0 for a watertight one), and `quality`, the [`compute_quality`](mesh_quality.md) summary with its histograms), `num_winding_rechecks`, `timings` (seconds per stage) and `warnings`. Anything that misses its target is reported in `warnings`, never hidden.

## Determinism

Every leaf and point is an integer coordinate on the finest lattice, deduplicated by sorting, and every parallel loop writes disjoint slots: the output is byte-identical across parallel backends and thread counts.

## Bindings

Python `mp.octree_mesh`, the pipeline step `OctreeMesh`, the CLI verb `octree-mesh` (both CLIs) and the MCP tool `octree_mesh` take the same options; see [pipelines](pipeline.md), [CLI](cli.md) and [MCP](mcp.md). The C API has `mio_octree_mesh` with `mio_octree_mesh_opts` and `mio_octree_mesh_report` ([C API](c_api.md)), and the [Fortran](fortran.md), [Julia](julia.md), [R](r.md) and [WASM](wasm.md) (`octreeMesh`) bindings wrap it. The notebook `example/python/41_octree_mesh.ipynb` walks through every output and draws the figures on this page.
