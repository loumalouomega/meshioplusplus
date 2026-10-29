# 3-D octree volume mesher: plan

This is the plan for `octree_mesh`, a new operation that meshes the volume a skin encloses on an octree: a closed triangle surface, or a slightly dirty one, from any readable surface format (STL, OBJ, PLY, …). It extends the [`octree`](./octree.md) of v17.8.0 into a complete 3-D volume mesher. It is a plan, not a feature page: each phase below narrows it as it lands, and the user-facing contract will live in `doc/octree_mesh.md`.

The operation has three outputs:

- **(A) Castellated hexahedra.** Surrogate-boundary tags and the signed distance and closest skin point at every node, for Kratos embedded and Shifted Boundary Method users.
- **(B) Conforming graded tetrahedra (the default).** Fitted to the skin, with its sharp features preserved.
- **(C) All-hexahedra.** Octree dual schemes, then snapping, pillowing and untangling.

The 2-D quadtree path and every existing operation keep their behaviour and their output byte for byte.

**Status (v17.12.0).** Phases 1–8 have landed and the contract is [`doc/octree_mesh.md`](./octree_mesh.md). Three things differ from the plan below:
- **Phase 5's `volume_fit`.** It became the operation's own optimisation stage: ODT relocation, sliding along the skin and flips under a guard on the smallest and largest dihedral angle.
- **Phase 7 (C1) landed later (v17.23.0).** `hex_grading="graded"` is the dual of the paired octree with the transition schemes ported from CinoLib; the default stays phase 7's uniform grid (C0), made manifold, pillowed, snapped and untangled.
- **Feature recovery.** A cylinder's rims at a small `min_size` were reported missed; since the snap guard counts only the tetrahedra the output keeps (v17.23.0), they are recovered.

Nothing of this plan remains open.

## What exists: `octree` (v17.8.0, PR #168)

`operations/octree.hpp` / `src/cpp/src/operations/octree.cpp` build a `hexahedron` octree around a surface, or a `quad` quadtree around a closed planar polyline. They reach Python (`_octree.py`, with a numpy twin of the hanging-node tree), C (`mio_octree`), Fortran, Julia, R, WASM, both CLIs' `octree` verb, the `Octree` pipeline step and the MCP `octree` tool.

**How the tree is built.**
- Leaves are a flat vector with integer lo corners on the finest lattice.
- Neighbour lookup uses one hash map per level, keyed by Morton code and rebuilt on every balance pass.
- Refinement culls each leaf's skin candidates top-down with the exact triangle–box overlap of `detail/tri_box.hpp`. The criteria are `sample_distance`, `surface_level`, `curvature_threshold`, `feature_angle` and `size_array`.
- `balance` is `none`, `face` or `full` (2:1), rippled to a fixed point.
- `fill` (`all`, `surface` or `inside`) is decided by the sign of the distance at each leaf centre: pseudonormal, or an O(T)-per-query winding number.

**What it emits.**
- Transitions are `hanging` (flagged in `octree:hanging`), `polyhedra` (split-hex polyhedra) or `templates` (cones to the centroid: pyramids and tetrahedra).
- `snap` projects the outer nodes through `shrinkwrap` and reverts any leaf with a non-positive corner Jacobian.
- The outputs are `octree:level`, `octree:cut`, `octree:leaf` and optionally `sdf:distance`.

**Tests.** `tests/cpp/test_octree.cpp`, `tests/python/test_octree.py` (with byte parity against the numpy twin), `tests/python/test_octree_transitions.py` (face sharing and volume over 30 neighbour patterns), the binding smoke tests and the `OctreeOptions` layout pin.

**What it lacks for volume meshing.**
- **Outputs.**
  - It produces no regions: no inside or cut-cell tags and no surrogate-boundary sides.
  - It produces no closest points.
  - It produces no conforming tetrahedral or hexahedral mesh fitted to the skin. `snap` moves only the hexahedra's outer nodes, and the transitions are not fitted.
- **Robustness to dirty skins.** A skin with holes, self-intersections or non-manifold edges is classified at leaf centres by a pseudonormal sign or an unaccelerated winding number, with no fast winding number.
- **Sizing.**
  - The tree is driven by levels, not by sizes: it has no `min_size`/`max_size` and no gap or thickness term.
  - It does not keep a uniform band next to the skin, which stuffing's guarantees need.
- **Reuse and coverage.**
  - There is no reusable linear octree: neighbour queries are per-level hash maps rebuilt on every balance pass.
  - There is no golden-digest test of its output and no benchmark row.

## Reusable pieces

| piece | where | used for |
|---|---|---|
| Morton keys, 21 bits per axis | `detail/space_filling.hpp` (`sfc_morton_key`, `sfc_stable_argsort`) | the linear octree, lattice vertex keys, BVH build order |
| Triangle–box SAT | `detail/tri_box.hpp` (`tri_box_overlap`, 13 axes) | candidate culling in refinement |
| Triangle soup, distance, closest point, feature-aware projection | `detail/surface_distance.hpp` (`build_triangle_soup`, `DistanceQuery`, `query_distances`, `query_closest_points`, `query_surface_projections`), `detail/point_triangle.hpp`, `detail/bucket_table.hpp` | `sdf:distance`, `sdf:closest_point`, cut-point root finding, patch projection |
| Deterministic numbering | `src/cpp/src/detail/slot_runs.hpp` | first-seen ids, deduplication |
| Feature edges | `src/cpp/src/detail/crease_edges.hpp`, `operations/feature_edges.hpp` | the feature graph |
| Curvature | `operations/curvature.hpp` (`mRecordPrincipal`) | the curvature size term |
| ODT smoothing | `src/cpp/src/operations/smooth_odt.hpp` (`smooth_odt_pass`) | the optimisation stage |
| Flips | `operations/optimize_volume.hpp` (2-3 / 3-2) | the optimisation stage |
| Quality | `operations/quality.hpp` (`compute_quality`, `attach_quality`, histograms), `operations/quality_gate.hpp` | the report and the quality gate |
| Two-sided distance | `operations/hausdorff.hpp` | `hausdorff_to_skin` |
| Watertightness | `skin.hpp` (`extract_skin`) | non-manifold boundary edge count |
| Parallel primitives | `parallel.hpp` (`parallel_for`, `parallel_sort` with a total order, `parallel_reduce`, `parallel_exclusive_scan`) | every hot loop, deterministically |

Two existing pieces are left untouched:
- `remesh_volume`'s lattice is a cube-centre fan of 12 tetrahedra, with dihedrals of 45°–120°, on a uniform grid. It is not Labelle and Shewchuk's BCC lattice, so their guarantee does not carry over. The new stuffing module is separate.
- `compute_sdf`'s adaptive octree is a sequence of `refine` calls on a hexahedral mesh; it shares nothing with this design.

## Algorithms

- **Linear octree.** Sorted 64-bit Morton keys (21 bits per axis) with a level per leaf, and a forest of cubic roots. Locating a cell is an `upper_bound`. 2:1 balance runs across faces, edges and corners, finest level first, by required parent-neighbour octants (Sundar, Sampath and Biros, SIAM J. Sci. Comput. 30(5), 2008; p4est, Burstedde, Wilcox and Ghattas, SIAM J. Sci. Comput. 33(3), 2011). The minimal balanced refinement of a tree is unique, so the result is deterministic. Vertex coordinates are doubled so that cell centres are integers, which limits `octree_mesh` to `root_cells << max_level < 2^20` and `max_level ≤ 19`, one bit below `octree`.
- **Inside and outside.** The fast generalised winding number (Jacobson, Kavan and Sorkine-Hornung, SIGGRAPH 2013; Barill, Dickson, Schmidt, Levin and Jacobson, SIGGRAPH 2018), thresholded at 0.5. It is evaluated on a BVH whose triangles are sorted by the Morton key of their centroid and split at the median count, so the tree is a function of the triangle count alone. The far field uses the dipole and second-order moments (β = 2), the near field the exact solid angle (Van Oosterom and Strackee 1983), and a value within 0.1 of the threshold is re-evaluated at β = 8. Holes, self-intersections and non-manifold edges are tolerated by construction. A globally inverted skin is negated with a warning, and a point on the skin gets sign 0.
- **Size field.** `h(x) = clamp(min(h_max, h_curv, h_feat, h_gap, h_user), h_min, h_max)`, with:
  - `h_curv = 2 sin(θ/2) / κ_max`, θ being `curvature_angle`;
  - `h_feat`, the size along feature edges;
  - `h_gap = t / gap_cells`, where t is the local thickness from an inward ray cast on the BVH;
  - `h_user` from `size_array`.

  The sizes map to levels in log space: `L = ceil(log2(h_max/h_min))` and root size `R = h_min · 2^L`. A leaf splits while its edge exceeds the smallest size of the skin primitives that reach it, the rule `size_array` already applies. Candidates are culled with the exact triangle–box SAT (Akenine-Möller, JGT 2001).
- **Tetrahedral background mesh.** A conforming graded BCC tetrahedralisation of the fully balanced octree. Each leaf is a cone, from its centre, of fans over its face portions: a face shared with four finer leaves contributes four squares, and an edge splits where its midpoint is a vertex. Across a face between two leaves of the same size, the face centre is removed and its star retriangulated around the axis joining the two cell centres, which gives exactly the BCC tetrahedra. Vertices are ranked by Morton key, with no hashing and no tolerance. The six resulting shapes (T1–T6) have dihedral angles between 35.26° and 125.26°; a uniform region is all T1 (60°–90°). The construction is our own derivation and its angle table is pinned by tests, as `remesh_volume`'s lattice is.
- **Conforming to the skin.**
  - Isosurface stuffing (Labelle and Shewchuk, SIGGRAPH 2007) with α_long = 0.24799 and α_short = 0.40173 per edge class, the pair behind the 10.7°–164.8° dihedral bound. Cut points are found by Illinois root-finding on the signed field and keyed by lattice edge, so both sides of a face share them. Vertices are warped in a deterministic order, then the 1-, 2- and 3-inside stencils are applied with a shared-face diagonal rule; that rule is chosen by an exhaustive quality sweep in the tests.
  - `conformity="cleaving"` follows lattice cleaving (Bronson, Levine and Whitaker, IEEE TVCG 20(2), 2014) for skins with several regions, such as nested shells.
- **Sharp features.** Corners and chains come from the crease edges at `feature_angle`. Each corner snaps to the nearest boundary vertex that keeps every incident tetrahedron valid, and is never moved again. Each chain follows a shortest boundary path inside a tube around its polyline and is projected onto it monotonically. The other boundary vertices project onto their own patch. Every move is guarded: positive volume, and no drop in the local minimum dihedral.
- **Fit and untangle.**
  - Volume fitting with automatic feature preservation, after Protais, Cherchi and Livesu (2026, hal-05574616): boundary vertices are drawn onto the skin and its features while a foldover-free barrier untangler (Garanzha, Kaporin, Kudryavtseva, Protais, Ray and Sokolov, ACM TOG 40(4), 2021) keeps every element valid. The reference implementation (fprotais/Mesh_optimization) is studied only once its licence is checked, and nothing is copied from it.
  - Then ODT smoothing alternates with 2-3 / 3-2 flips, both under a dihedral-floor guard.
- **All-hexahedra.** The octree dual schemes of Livesu, Pitzalis and Cherchi (ACM TOG 41(2), 2022) on the balanced octree, taking the pairing rule and transition templates from the paper's text (CinoLib is MIT and may be consulted). Generalised pairing follows Pitzalis, Livesu, Cherchi, Gobbetti and Scateni (ACM TOG 40(6), 2021). The hexahedra inside the skin are selected by winding number and repaired to a manifold boundary, then projected, pillowed and untangled (Owen et al., Sculpt; Gao, Shen and Panozzo, CGF 38(5), 2019). snappyHexMesh is GPL and is an algorithmic reference only.

Every module states in its header comment which publications it follows and that no code was taken from an implementation with an incompatible licence.

## API

C++: a new installed header `operations/octree_mesh.hpp`.

```cpp
enum class OctreeMeshOutput { Castellated, Tetra, Hexa };       // default Tetra
enum class OctreeMeshConformity { Stuffing, Cleaving };         // default Stuffing
enum class OctreeMeshSkinGrading { Graded, Uniform };           // default Graded
enum class OctreeMeshFill { All, Inside, Cover };               // castellated output; default Inside

struct OctreeMeshOptions {
    OctreeMeshOutput mOutput = OctreeMeshOutput::Tetra;
    double mMinSize = 0.0, mMaxSize = 0.0;                      // 0: derived from the skin's bounding box
    std::int64_t mMaxLevel = -1;                                // -1: from the sizes
    std::optional<std::array<double, 6>> mBounds;
    double mPadding = 0.0, mPaddingRelative = 0.1;
    double mCurvatureAngle = 15.0;                              // degrees; 0 disables
    double mFeatureAngle = 30.0, mCornerAngle = 30.0;           // degrees
    double mFeatureSize = 0.0;                                  // 0: min size
    double mGapCells = 2.0;                                     // 0 disables
    bool mGapOutside = false;
    std::string mSizeArray, mSurfaceRegion;
    std::int64_t mGradingLayers = 1;
    OctreeMeshSkinGrading mSkinGrading = OctreeMeshSkinGrading::Graded;
    double mSurfaceSize = 0.0;                                  // uniform skin grading; 0: the smallest size asked for on the skin
    double mWindingThreshold = 0.5, mWindingBeta = 2.0;
    OctreeMeshConformity mConformity = OctreeMeshConformity::Stuffing;
    double mAlphaLong = 0.24799, mAlphaShort = 0.40173;
    bool mPreserveFeatures = true;
    std::int64_t mFitIterations = 5;
    bool mDihedralFloorGuard = true;
    OctreeMeshFill mFill = OctreeMeshFill::Inside;
    bool mKeepOutside = false, mAttachDistance = true;
    bool mPillow = true;
    std::int64_t mUntangleIterations = 50;
    std::int64_t mMaxCells = 20000000, mMaxTets = 100000000;
};

MESHIOPLUSPLUS_API OctreeMeshResult octree_mesh(const Mesh& rSkin, const OctreeMeshOptions& rOptions = {});
```

`OctreeMeshResult` carries:
- the mesh;
- the frame (origin, root size, root counts, finest and skin levels);
- leaves per level and cells per type;
- hanging, warped and cut counts;
- minimum and maximum dihedral, minimum scaled Jacobian and inverted cells;
- `hausdorff_to_skin` and `num_non_manifold_boundary_edges`;
- corners and chains snapped and missed;
- `guarantee_applies` (uniform skin grading and nothing snapped to features);
- winding-number re-checks;
- the `QualityReport` with its histograms;
- the time spent in each stage;
- the warnings.

A stage that misses its quality target says so in the warnings and the report. Nothing degrades silently.

**Outputs.**
- Point data: `sdf:distance` and `sdf:closest_point`.
- Cell regions: `inside` and `boundary_cells`.
- Side regions: `surrogate_boundary` for (A), `boundary` for (B) and (C).
- Point regions: `feature_corners` and `feature_edges`.
- Cell data: `octree:level`, `octree:cut`, `octree_mesh:shape` (the lattice shape class) and `quality:*`.

In (A), hanging nodes are kept and flagged in `octree:hanging`, and an inside leaf beside a finer one that is not inside is split until every surrogate face is a whole cell face.

**Every other surface mirrors the C++ options:**

| surface | entry point |
|---|---|
| Python | `meshioplusplus.octree_mesh(skin, min_size=…, max_size=…, output="tet", return_report=True)`, core only; without the compiled core it raises `NotImplementedError` by name |
| Both CLIs | `octree-mesh skin.stl out.vtu --min-size 0.01 --max-size 0.2 --output tet --json` |
| Pipeline | the `OctreeMesh` step |
| MCP | the `octree_mesh` tool |
| C | `mio_octree_mesh`, `mio_octree_mesh_opts` and `mio_octree_mesh_report` (with reserved tails), `mio_octree_mesh_opts_init` |
| Fortran | `mesh%octree_mesh` |
| Julia | `octree_mesh` |
| R | `mio_octree_mesh` |
| WASM | `octreeMesh` |

## ABI impact

**C++.**
- `operations/octree_mesh.hpp` is a new header of new types and a new free function: tier C, additive.
- `OctreeMeshOptions` is pinned in `tests/cpp/test_abi_layout.cpp` from its first release.
- No existing declaration changes. `octree.hpp`, `remesh_volume.hpp` and every other installed header are untouched, and the new modules are core-private (`src/cpp/src/detail/`).
- The release that ships it (17.12.0) already carries ABI 33 from the `fit_cad` freeform work. `tools/check-abi-version.sh` decides whether the new header also needs a row in [ABI reviews](./abi_reviews.md).

**C.** New structs and functions only.

## Phases

One branch, one or more commits per phase, and each phase green in CI before the next starts.

| phase | delivers | tests |
|---|---|---|
| 1. Octree core | `src/cpp/src/detail/linear_octree.{hpp,cpp}`; golden digests of the existing `octree` output (2-D and 3-D, every transition, snap) before anything else | keys ascending; leaves aligned; no overlap; full coverage of the roots; balance per adjacency against brute force on seeded random trees; idempotence; refinement only; equality with `octree(balance="full")`; identical across SEQ/OpenMP/TBB and 1/N threads |
| 2. Skin refinement and classification, output (A) | `tri_bvh` (fast winding number, ray cast), `octree_size_field`, castellated `octree_mesh` with its regions, `sdf:distance` and `sdf:closest_point`; Python, both CLIs, MCP | winding number against brute force on a sphere, torus, open hemisphere, overlapping spheres, inverted sphere and a triangle soup; curvature size on a sphere; gap cells across a thin plate; regions against brute force; closed surrogate surface; hanging count |
| 3. Conforming tetrahedra on the balanced octree | `bcc_octree_tets` | exact integer volume sum; every interior face shared by exactly two tetrahedra with opposite orientation; positive orientations; shape classes T1–T6 and their dihedrals; a uniform tree gives only T1; no vertex inside an edge |
| 4. Stuffing and feature snapping, output (B) | `isosurface_stuffing` (with cleaving), `feature_snap` | every sign pattern on every shape gives positive volumes that add up to the clipped volume; the stencil quality sweep stays within 10.7°–164.8° on T1; uniform sphere, torus and ellipsoid ≥ 10.7°; nested spheres give two regions; cube, fichera, cylinder and a fandisk-like part keep their corners (within 1e-12) and edges (feature Hausdorff ≤ 0.05 h) |
| 5. Fit, untangle, optimise, quality gate | `volume_fit`, the report as JSON, benchmark rows (`octree` first as a baseline, then `octree_mesh_castellated`, `octree_mesh_tet`) | no inverted cell; uniform skin grading ≥ 10°; graded minimum dihedral reported; a closed 2-manifold `extract_skin`; Hausdorff to the skin ≤ c·h_min; volume close to the analytic one |
| 6. Bindings and docs | C, Fortran, Julia, R, WASM, both CLIs, pipeline and MCP in full; `doc/octree_mesh.md` with diagrams; the binding pages; README; `example/python/41_octree_mesh.ipynb`; CHANGELOG | each binding's tests and layout checks; MCP parity; CLI JSON parity |
| 7. All-hexahedra, output (C) | uniform dual first (C0), then the graded dual with transition templates (C1) | only hexahedra; each quad shared by at most two cells; closed manifold boundary; minimum scaled Jacobian > 0 (target ≥ 0.2, reported); Hausdorff bound; determinism |
| 8. Release 17.12.0 | the ten version files, `find_package` pins, `BASELINE_HASHES`, CHANGELOG, single header | the version-bump checks of `AGENTS.md`; the Spack recipes follow the tag |

**Test skins.** Generated in the tests: unit sphere, cube (sharp edges and corners), torus (genus 1), cylinder, a fandisk-like part with creases, a thin plate (thinner than 2·h_min), two nested spheres, a slightly holed sphere and a self-intersecting triangle soup. The repository's bracket fixture is added where available.

**Every phase keeps these:**
- 2:1 balance.
- No inverted cell.
- A watertight boundary.
- Byte-identical output across SEQ, OpenMP and TBB and across thread counts (`tools/bench_ops.sh --hash`).
- Unchanged `octree` goldens.
- No O(N²) pass, and about a million tetrahedra in well under a minute on a laptop with OpenMP.

## Skin grading and the angle guarantee

Labelle and Shewchuk's 10.7°–164.8° bound holds only where the skin cuts uniform BCC tetrahedra. By default the size field also grades along the skin, so the tetrahedra cut where the level changes are transition shapes, and the report states the actual minimum dihedral.

`skin_grading="uniform"` instead places every leaf the skin cuts at one level. It keeps two rings of leaves around them at that level and pads the domain by at least three skin sizes, so every tetrahedron that stuffing touches is a T1 and the bound applies. The tests assert at least 10° in that mode, before and after optimisation. Snapping to features voids the guarantee near the features, and the report's `guarantee_applies` says so.

## Roadmap

[`doc/roadmap.md`](./roadmap.md) admitted lattice generation and listed hexahedral meshing from a CAD-like boundary as a non-goal. The octree dual of this plan is lattice generation (index arithmetic on a balanced octree, then projection and untangling), so the non-goal is narrowed to exclude it. Hex-dominant meshing and general hexahedral meshing from a boundary stay out. The work has landed; see [`octree_mesh`](./octree_mesh.md).

## Risks

| risk | mitigation |
|---|---|
| Labelle and Shewchuk's exact shared-face diagonal (parity) rule | the rule is pluggable and chosen by the exhaustive stencil quality sweep |
| The dihedral angles of the T5 and T6 transition shapes are derived by hand | the tests recompute them; conformity does not depend on them |
| Stuffing cost: root finding needs many field evaluations per cut edge | an exact segment–triangle hit on the BVH can give the first bracket |
| A feature chain on a coarse graded skin may not snap | it is reported as missed, never dropped silently |
| All-hexahedra: the template set and pairing rule, concave features and thin regions | taken from the papers' text; the minimum scaled Jacobian is always reported |
| Licences of reference implementations | checked before reading; nothing copied |
