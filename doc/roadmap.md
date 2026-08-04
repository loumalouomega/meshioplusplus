# meshio++ roadmap

Status at time of writing: **v9.19.0** — 41 formats, twenty mesh operations + five data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers, an MCP server, a settings-driven pipeline engine, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 5).

This document lists what is *not* built. Items are grouped by theme, each with an effort estimate and the reason it matters. Nothing here duplicates shipped functionality; where a feature partially exists, the gap is stated explicitly.

Effort key: **S** = days, **M** = a couple of weeks, **L** = a month or more, **XL** = a project in its own right.

Multi-file / transient workflows (glob input, fan-in/fan-out, per-step pipeline
execution, and Python's `TimeSeries` for random-access "hold a series as one
value") shipped in full in v9.12.0 across every language surface including
WASM — see [`doc/sequences.md`](sequences.md) — and so no longer appears here.

MDPA's arbitrary/gapped node ids (v9.13.0, read side) and original-id
preservation on write, including a fixed `SubModelPart` stale-reference bug
found along the way (v9.14.0, write side) shipped in full — see
[`doc/formats/mdpa.md`](formats/mdpa.md#original-ids-preserved-on-write-v9-14-0)
— and so no longer appears here.

---

## 1. Polyhedral meshes

**The gap.** Ragged polyhedron blocks exist in all three backends (`AddPolyhedronBlock`, CSR / nested storage, `CellView::NumFaces`/`Face`), and MED, EnSight `nfaced` and OpenFOAM read them. Since v9.15.0 they also cross the flat C ABI in both directions, so Fortran / Julia / R / WASM can build and read one — see [`doc/polyhedra.md`](polyhedra.md). Since v9.16.0 there is also a geometric kernel for them (`detail/polyhedron.hpp`: volume, centroid, area, and the winding repair all of those depend on), and `stats`, `gradient`, `quality`, `smooth`'s inversion guard, `extract_surface`/`extract_skin` and `data_average`'s measure weighting all go through it. Since v9.17.0 `convert_cells(simplexify)` decomposes them into tetrahedra with the same fan, which also makes `slice`, `isosurface` and `interpolate --barycentric` work. **What remains is the writers** — and the fixtures, which depend on being able to read a real polyhedral file at all. Since OpenFOAM — the most-used open CFD code — is natively polyhedral, this is a real ceiling.

- **OpenFOAM writer** — the format is still read-only, and it is the last polyhedral writer missing (VTU type 42, MED `POE` and EnSight `nfaced`/`nsided` all landed in v9.19.0). It is also the hardest: unlike the other three it needs a **globally deduplicated face table** with owner/neighbour, and an ordering contract that silently produces a mesh `checkMesh` rejects if broken (internal faces first sorted by owner then neighbour, `owner < neighbour`, normal pointing owner→neighbour, boundary faces contiguous per patch). It would be the only writer taking a *directory* path, and round-tripping a patch's `type` (`wall` vs `patch`, which changes solver behaviour) needs a new `OpenFoamInfo` member — a **Tier A ABI change**, so this release bumps `MESHIOPLUSPLUS_ABI_VERSION`. Round-tripping an OpenFOAM case is the headline outcome. **M–L**
- **CGNS `NGON_n`/`NFACE_n` on the write side**: `write_cgns` still refuses a ragged block. It needs the *same* global face dedup as OpenFOAM (an `NFACE_n` cell is a list of signed face ids into a shared `NGON_n`), so the two belong in one change over a shared `detail/face_mesh.hpp` rather than building that primitive twice. **S–M**
- **Test fixtures**: there is no genuinely polyhedral file fixture anywhere in the tree — `med/voronoi_hex.med` is 2-D `POG` polygons, `ensight/simple.geo` has `nsided` but no `nfaced`, the `vtu/*.vtu` files are tetra, and no OpenFOAM case is checked in at all. The only polyhedral fixtures are synthetic and in-code (`tests/python/helpers.py`'s `polyhedron_mesh`, plus `AddPolyhedronBlock` literals in the gtests), and `helpers.py`'s is **inconsistently wound** — both of its `polyhedron5` cells traverse an edge in the same direction from two incident faces — which is exactly why the geometric kernel above has to repair winding rather than assume it. Generate fixtures by round-tripping through meshio++'s own writers once they exist (the repo's stated preference, and how the `.cgns`/`.ex2` fixtures came about) rather than importing third-party files; anything third-party follows `tests/python/meshes/exodus/` to the letter (`.license` sidecar + `LICENSE.<SPDX>` + README + `CITATION.cff` credit), and the rule against GPL-sourced test data still applies. **S**

  *Note:* the CGNS project's examples page ([cgns.org/current/examples.html](https://cgns.org/current/examples.html)) became usable in v9.18.0, on a build with the optional [cgnslib backend](formats/cgns.md) — which reads both `NGON_n`/`NFACE_n` and the ADF container those files often use. Without that flag they remain unreadable.


*Transferable since v9.15.0 (C ABI), measurable since v9.16.0 (the kernel), operable since v9.17.0 (the tetrahedra decomposition). Recommended next: the cgnslib backend, which is what makes real polyhedral fixtures readable at all, and then the writers.*

---

## 2. Signed distance fields for skin meshes (octree)

**The gap.** meshio++ has surface extraction (`extract_surface`/`extract_skin`), slicing and isosurfacing, but no way to answer "how far is this point from the surface" — the primitive collision detection, offsetting, and voxelization-style ML preprocessing all need. A closed skin mesh (STL and friends: watertight triangle soup, no volume topology) is exactly the input shape this needs. Two spatial structures matter here and both should come out the other end as an ordinary, **exportable** `Mesh`, not a bespoke in-memory-only object: a sparse **octree** (adaptive, coarse away from the surface, refined near it — cheap on RAM, the acceleration structure for the nearest-triangle queries themselves) and a dense **regular voxel grid** (uniform cell size, the shape most voxel/ML tooling and simple boolean pipelines actually expect). They should be two generation modes of the same feature, not two separate ones.

- **Spike: where does the grid live in the data model?** The natural encoding is a `hexahedron` (or `custom` for an octree with T-junction/hanging-node leaves) `CellBlock` per cell — a voxel grid is trivially one uniform block, an octree needs the ragged/`custom` path or a balancing pass (`refine`'s `Balanced` closure is the existing precedent for 2:1-balanced hanging nodes). Getting this right is what makes "exportable" free: once it *is* a `Mesh`, every existing writer (VTU, VTK, gmsh, …) already handles it, `view`/`screenshot` already render it, and no new file format needs writing. Decide the exact cell layout and any non-mesh metadata's (bounds, depth, cell size) home (a `SdfInfo` side-channel, the `MedInfo`/`GmshInfo` precedent) before committing, following the NURBS-spike precedent (§10). **S**
- **Voxel grid generation**: `voxelize(mesh, VoxelOptions{resolution|cell_size, bounds, band})` — a dense uniform `hexahedron` block over the mesh's (optionally padded) bounding box. The simpler of the two structures; a good first deliverable and the base case the octree generalizes. **M**
- **Octree construction over a closed triangle skin**: adaptive subdivision to a max depth / target cell size, refined near the surface (leaf-triangle bucketing reusing `detail/spatial_hash.hpp`'s bucket-grid idiom), coarse away from it, and balanced/exported the same way the voxel grid is. **M**
- **Signed distance evaluation**, shared by both structures, per cell (corner or center): nearest-triangle distance (the same bucket search) plus a sign — fast winding number, or angle-weighted pseudo-normals for a possibly-imperfect/non-watertight STL. Stored as ordinary `point_data`/`cell_data` (e.g. `sdf:distance`), so it is exported and colored exactly like any other field, with no new data convention needed. **M**
- **`compute_sdf(mesh, SdfOptions{structure: octree|voxel, resolution|max_depth, band, watertight_check})`** as a new operation, exposed like the others (pybind / C-ABI / Fortran / WASM / a `sdf` CLI verb) returning the grid as a `Mesh` plus a point-query API (`sample(x, y, z)`) for callers who just want values, not the mesh. **L**
- **Downstream writers/consumers**: since the result is a real `Mesh`, `write(grid, "out.vtu")` and friends already work; add only what a plain mesh writer cannot express — a dense binary voxel dump (NRRD/raw-array style) for ML tooling that wants a tensor, not a mesh. Once the primitive exists, offsetting and inside/outside queries for `crop`/`merge` follow cheaply. **S–M**

*Recommended entry point: the data-model spike (get the cell-layout choice right so both structures are exportable for free), then voxel generation (the simpler case), then the octree and shared SDF evaluation on a single closed STL skin.*

---

## 3. Machine-learning data handling

**The gap.** v8.2.0 gave Arrow/Parquet export of `point_data`/`cell_data`, which is the right primitive but only the first step. ML pipelines want *datasets* (many meshes), tabular frames, batched tensors and stable feature layouts — none of which exist.

- **pandas / polars frames** directly (`to_pandas(mesh, location=...)`), not only via pyarrow. Trivially thin over the existing table payload, and it is what people actually reach for. **S**
- **Dataset-level export**: a *directory* of meshes → one partitioned Parquet dataset with a `mesh_id` column, using the sequence machinery (see `doc/sequences.md`). This is the format an ML training loop wants. **M**
- **Feature extraction helpers** — assemble a canonical per-node or per-cell feature matrix (coordinates, selected fields, derived quantities like `quality:*` or `|∇f|`, region one-hots) with a **stable, recorded column order**, so training and inference cannot silently disagree. The column-order contract is the whole value; make it explicit and versioned. **M**
- **Graph export for GNNs**: node/edge index arrays in the layout PyTorch Geometric and DGL expect (`edge_index` as (2, E)), from the existing `node_adjacency` / cell-dual machinery. The cheapest genuinely ML-shaped feature in this list. **S–M**
- **PyTorch / JAX tensor handoff** via the DLPack path already built for GPU (v-GPU work), so a mesh becomes a batch of device tensors without a file round-trip. Mostly already there — needs the framework-facing convenience layer and docs. **S**
- **HDF5/Zarr chunked dataset writer** for datasets too large to hold in memory, with a documented on-disk layout. Only worth it once the dataset-level export exists. **L**

*Recommended entry point: pandas + `edge_index` + the feature-matrix contract — a week of work that makes meshio++ directly usable from a training script.*

---

## 4. NVIDIA PhysicsNeMo integration

**The gap.** PhysicsNeMo (github.com/NVIDIA/physicsnemo) is the mainstream open Physics-ML framework, and its data ingestion is where most users write bespoke glue. meshio++ already has 41 readers, GPU/DLPack handoff, and the operations (`interpolate`, `partition`, `gradient`, `decimate`) that a training pipeline needs for preprocessing. A thin, well-documented bridge would let people train on simulation output without integrating their solver at all — which is exactly the friction PhysicsNeMo users hit.

- **Reconnaissance first, and treat it as a real deliverable.** PhysicsNeMo's dataset/datapipe contracts, its mesh and point-cloud conventions, and its dependency weight (CUDA-specific, container-oriented) all need checking against the repo's "optional, gated, never in `[all]`" rule. Write the findings down before writing code — the CuPy packaging finding is the precedent for how this repo handles such constraints. **S**
- **A `physicsnemo` optional extra + dataset adapter**: a meshio++-backed dataset class yielding the tensors PhysicsNeMo's datapipes expect, built on the §3 feature-matrix contract and the existing DLPack handoff. Pure Python, lazily imported, named install error. **M**
- **Preprocessing recipes as pipeline documents** — sampling, normalisation, surface extraction, decimation, partitioning into training patches — expressed as v9.11.0 `settings.json` files so they are reproducible and reviewable rather than notebook cells. A strong fit for the pipeline engine, and cheap once the adapter exists. **S–M**
- **A worked end-to-end example**: simulation output → meshio++ preprocessing → PhysicsNeMo training → inference results read back as a mesh and rendered. The example *is* the feature; without it the adapter will not be adopted. **M**
- **CI reality check**: PhysicsNeMo needs a GPU, which public runners do not have. Follow the precedent set for the GPU work — test the pure adapter logic without the framework, gate the rest, and state plainly that the integration path is not covered by public CI. **S**

*Recommended entry point: the reconnaissance note, then the dataset adapter plus one worked example. Do not build the adapter before writing down what PhysicsNeMo actually expects.*

---

## 5. Remaining refinement and coarsening gaps

`refine` is adaptive (v9.5.0) and `decimate` exists, but the pair still has holes.

- **Volume decimation** — `decimate` is surface-only by documented design; tet-collapse validity is the hard part. **L**
- **Refinement hierarchy across passes** — `refine:level` exists per pass; a persistent parent/child hierarchy is what multigrid and green-element undo need. **M**
- **Error-estimator helpers** — now that `gradient` exists, a gradient-jump or recovery-based indicator that feeds `refine`'s selection directly closes the adaptive loop end to end. **M**

---

## 6. Field capability beyond derivatives

- **Conservative (mass-preserving) interpolation** — `interpolate`'s barycentric mode is pointwise; CFD remapping needs conservation. **L**
- **Field integration** — total, mean, and per-region reductions over cells as a `data` verb; the natural companion to `gradient`. **S**
- **Second derivatives / Hessian**, for curvature-based adaptivity. **M**

---

## 7. Scale

The benchmark is a ~52k-node bracket; nothing addresses meshes that do not fit in RAM.

- **A large-mesh benchmark tier** (10M+ cells) — cheap, and it would show whether the parallel paths actually hold. Do this before the two below, since it decides whether they matter. **S**
- **Streaming / chunked writes**, the counterpart to the selective-read work. **L**
- **Out-of-core operations** for the ops that are already block-local. **XL**

---

## 8. Ecosystem reach

- **Blender add-on** — Blender ships Python and reads almost no FEA formats; unusually high visibility per line of code. **S–M**
- **Rust bindings** over the C API — the next language by scientific adoption after Julia/R, and the ABI/`SOVERSION` work makes it cheap. **M**
- **Registration and distribution** — conda-forge, CRAN, Julia General, a proper ParaView reader plugin. All deferred at binding time; all pure logistics, and all blocking real adoption. **M**

---

## 9. Quality of implementation

- **Fuzzing the readers** (libFuzzer / AFL, OSS-Fuzz if it will take the project). 41 mostly hand-rolled parsers, reachable from a C ABI, a browser and an MCP server — untrusted input reaches them by design. The highest-value non-feature item in this document. **M**
- **A format conformance matrix** — one canonical mesh written to and read back from every format, with declared per-format lossiness, generalising the region round-trip test into executable documentation of what survives what. **M**
- **Property-based testing** (Hypothesis) over the invariants already articulated in the docs: partition-of-unity, volume conservation, conformity, byte-identical determinism. **M**

---

## 10. NURBS and higher-order geometry (long run)

**The gap.** The data model is strictly linear/Lagrange polytopes: a `CellBlock` is a cell-type string plus a node-index array. NURBS is a genuinely different object — control points, weights, knot vectors, and a parametric mapping — and CAD/IGA formats (STEP, IGES, Rhino 3dm, `.iga`) express geometry that no current cell type can hold. This is the most architecturally invasive item on the list and should be approached as a research spike, not a feature.

- **Spike: how far can the current model stretch?** Higher-order Lagrange cells already exist (`hexahedron27`, VTK-Lagrange types); a rational Bézier/NURBS patch needs *weights* and a *knot vector*, which have nowhere to live. Determine whether a side-channel struct (the `MedInfo`/`GmshInfo` precedent) suffices, or whether the `Mesh` needs a genuine second entity kind. Write the finding up before committing. **M**
- **Read-only CAD ingestion first**: a NURBS surface tessellated to a triangle mesh at a requested tolerance, with the parametric data carried out-of-band. This delivers most of the practical value (getting CAD into the mesh world) without touching the data model, and is the natural first release. **L**
- **A real IGA data model** — patches, control nets, weights, knots, trimming curves — plus formats and evaluation. This is XL, likely a separate library or a major version, and should only be attempted if the spike shows real demand.
- **Dependency reality**: robust STEP/IGES parsing effectively means OpenCASCADE, which is a heavyweight LGPL dependency. If ingestion goes ahead, it must follow the KaHIP/Polyscope pattern — strictly optional, off by default, never in the core, licence implications documented. **Findings before code.**

*Recommended posture: spike and document; do not schedule implementation until the spike says what shape it takes.*

---

## 11. Mesh generation

**The gap.** Every operation transforms a mesh you already have; nothing creates one. This is the only empty category in the operations layer.

- **Primitive constructors** — `box`, `grid(nx,ny,nz)`, `sphere`, `cylinder`, `disk`. Trivial, dependency-free, and it removes the fixture-file dependency from tests, docs, notebooks, the browser demo and the MCP server. Highest leverage per line of code in this document. **S**
- **`extrude`** — 2D → 3D sweep (triangle→wedge, quad→hexahedron), `nlayers`, per-layer offsets. The most-requested generation primitive; repeatedly deferred. **M**
- **`revolve`** — extrude's rotational sibling, sweeping around an axis. **M**
- **Delaunay / constrained 2D meshing** — genuinely useful, but robust geometric predicates are where dependency-free stops paying. Better as an optional Triangle or Gmsh backend, following the KaHIP pattern. **L**

---

## Suggested sequencing

1. **Primitive constructors (§11, first item)** — a few days, and it improves testing, docs and every demo surface at once.
2. **ML data handling (§3)** — pandas, `edge_index`, and the feature-matrix contract; this is also the prerequisite for §4.
3. **PhysicsNeMo reconnaissance (§4, first item)** — a written findings note before any code.
4. **Polyhedral C-ABI exposure + geometric kernel (§1)** — lifts a real ceiling and unblocks OpenFOAM round-tripping.
5. **SDF/octree spike (§2, first item)** — a documented data-model decision before the octree/SDF work is scheduled.
6. **Fuzzing (§9)** — should start in parallel with all of the above; it is not a feature and does not compete for the same attention.
7. **NURBS spike (§10)** — a documented investigation, scheduled independently of the rest.
