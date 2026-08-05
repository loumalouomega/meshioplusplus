# meshio++ roadmap

Status at time of writing: **v9.24.0** — 41 formats, twenty-three mesh operations + five data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers, an MCP server, a settings-driven pipeline engine, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 6).

This document lists what is *not* built. Items are grouped by theme, each with an effort estimate and the reason it matters. Nothing here duplicates shipped functionality; where a feature partially exists, the gap is stated explicitly.

Effort key: **S** = days, **M** = a couple of weeks, **L** = a month or more, **XL** = a project in its own right.

Multi-file / transient workflows (glob input, fan-in/fan-out, per-step pipeline
execution, and Python's `TimeSeries` for random-access "hold a series as one
value") shipped in full in v9.12.0 across every language surface including
WASM — see [`doc/sequences.md`](sequences.md) — and so no longer appears here.

Polyhedral meshes are first-class end to end as of v9.22.0: ragged blocks cross
the flat C ABI in both directions, a geometric kernel measures them, they
decompose into tetrahedra, and MED, EnSight, VTU, OpenFOAM and CGNS all read
*and* write them — see [`doc/polyhedra.md`](polyhedra.md) — and so no longer
appears here. Two leftovers live where they belong rather than here: OpenFOAM's
writer is ASCII-only (a per-format encoding gap, recorded under
[Quirks & limitations](formats/openfoam.md#quirks-limitations)), and `refine`
and `decimate` still raise on a polyhedron, which is a refinement gap and is
listed as one below.

MDPA's arbitrary/gapped node ids (v9.13.0, read side) and original-id
preservation on write, including a fixed `SubModelPart` stale-reference bug
found along the way (v9.14.0, write side) shipped in full — see
[`doc/formats/mdpa.md`](formats/mdpa.md#original-ids-preserved-on-write-v9-14-0)
— and so no longer appears here.

---

## 1. Signed distance fields for skin meshes (octree)

**Mostly shipped in v9.24.0.** `grid`, `voxelize` (fills `all`/`surface`/`inside`),
`sample_distance`, `distance_to_surface` and `surface_watertight_check` are live on
every surface — see [`doc/voxelize.md`](voxelize.md) and [`doc/sdf.md`](sdf.md).
The data-model question the spike existed to answer is settled: both structures are
an ordinary `hexahedron` `CellBlock`, so export, rendering and every existing
operation come for free. The roadmap's `sample(x, y, z)` point-query API is
satisfied by `sample_distance`, which takes an array of points rather than one at
a time — a batch call, deliberately, so no opaque sampler handle is needed on five
flat bindings.

What remains:

- **`compute_sdf(mesh, SdfOptions{...})`** — generate the grid *and* fill it in one
  call. Declared with its option and result layouts final (so adding it is a pure
  `.cpp` change) but currently throws by name; compose `voxelize` with
  `distance_to_surface` in the meantime. **S**
- **Octree construction** — adaptive subdivision refined near the surface, coarse
  away from it. The mechanism is settled and unblocked: iterated selective
  `refine(closure="balanced")` over a root lattice, which already gives 2:1
  balance, `refine:level` and `refine:hanging` (and whose cross-pass tearing bug was
  fixed in v9.23.0 precisely for this). Two traps are recorded in the plan: compute
  the field only after the last pass (`refine` interpolates `point_data`), and
  rebuild the cell selection from each pass's own output. **M**
- **A dense binary voxel dump** for ML tooling that wants a tensor rather than a
  mesh. `.vti` (VTK XML ImageData) is the recommendation over NRRD/raw: it reuses
  `detail/vtk_xml.hpp`'s DataArray codec and the existing block codecs verbatim, and
  its `Origin`/`Spacing`/`WholeExtent` attributes *are* the grid header — which
  matters because **no format persists arbitrary `field_data`**, so a
  written-and-reread grid otherwise loses its lattice metadata and has to recover it
  from the geometry. This is a *format*, with the whole format checklist. **S–M**
- **Offsetting and inside/outside predicates for `crop`/`merge`**, which the
  primitive now makes cheap. **S–M**

*Recommended entry point: `compute_sdf` (small, and it makes the octree a flag
rather than a new entry point), then the octree.*

## 2. Machine-learning data handling

**The gap.** v8.2.0 gave Arrow/Parquet export of `point_data`/`cell_data`, which is the right primitive but only the first step. ML pipelines want *datasets* (many meshes), tabular frames, batched tensors and stable feature layouts — none of which exist.

- **pandas / polars frames** directly (`to_pandas(mesh, location=...)`), not only via pyarrow. Trivially thin over the existing table payload, and it is what people actually reach for. **S**
- **Dataset-level export**: a *directory* of meshes → one partitioned Parquet dataset with a `mesh_id` column, using the sequence machinery (see `doc/sequences.md`). This is the format an ML training loop wants. **M**
- **Feature extraction helpers** — assemble a canonical per-node or per-cell feature matrix (coordinates, selected fields, derived quantities like `quality:*` or `|∇f|`, region one-hots) with a **stable, recorded column order**, so training and inference cannot silently disagree. The column-order contract is the whole value; make it explicit and versioned. **M**
- **Graph export for GNNs**: node/edge index arrays in the layout PyTorch Geometric and DGL expect (`edge_index` as (2, E)), from the existing `node_adjacency` / cell-dual machinery. The cheapest genuinely ML-shaped feature in this list. **S–M**
- **PyTorch / JAX tensor handoff** via the DLPack path already built for GPU (v-GPU work), so a mesh becomes a batch of device tensors without a file round-trip. Mostly already there — needs the framework-facing convenience layer and docs. **S**
- **HDF5/Zarr chunked dataset writer** for datasets too large to hold in memory, with a documented on-disk layout. Only worth it once the dataset-level export exists. **L**

*Recommended entry point: pandas + `edge_index` + the feature-matrix contract — a week of work that makes meshio++ directly usable from a training script.*

---

## 3. NVIDIA PhysicsNeMo integration

**The gap.** PhysicsNeMo (github.com/NVIDIA/physicsnemo) is the mainstream open Physics-ML framework, and its data ingestion is where most users write bespoke glue. meshio++ already has 41 readers, GPU/DLPack handoff, and the operations (`interpolate`, `partition`, `gradient`, `decimate`) that a training pipeline needs for preprocessing. A thin, well-documented bridge would let people train on simulation output without integrating their solver at all — which is exactly the friction PhysicsNeMo users hit.

- **Reconnaissance first, and treat it as a real deliverable.** PhysicsNeMo's dataset/datapipe contracts, its mesh and point-cloud conventions, and its dependency weight (CUDA-specific, container-oriented) all need checking against the repo's "optional, gated, never in `[all]`" rule. Write the findings down before writing code — the CuPy packaging finding is the precedent for how this repo handles such constraints. **S**
- **A `physicsnemo` optional extra + dataset adapter**: a meshio++-backed dataset class yielding the tensors PhysicsNeMo's datapipes expect, built on the §2 feature-matrix contract and the existing DLPack handoff. Pure Python, lazily imported, named install error. **M**
- **Preprocessing recipes as pipeline documents** — sampling, normalisation, surface extraction, decimation, partitioning into training patches — expressed as v9.11.0 `settings.json` files so they are reproducible and reviewable rather than notebook cells. A strong fit for the pipeline engine, and cheap once the adapter exists. **S–M**
- **A worked end-to-end example**: simulation output → meshio++ preprocessing → PhysicsNeMo training → inference results read back as a mesh and rendered. The example *is* the feature; without it the adapter will not be adopted. **M**
- **CI reality check**: PhysicsNeMo needs a GPU, which public runners do not have. Follow the precedent set for the GPU work — test the pure adapter logic without the framework, gate the rest, and state plainly that the integration path is not covered by public CI. **S**

*Recommended entry point: the reconnaissance note, then the dataset adapter plus one worked example. Do not build the adapter before writing down what PhysicsNeMo actually expects.*

---

## 4. Remaining refinement and coarsening gaps

`refine` is adaptive (v9.5.0) and `decimate` exists, but the pair still has holes.

- **Volume decimation** — `decimate` is surface-only by documented design; tet-collapse validity is the hard part. **L**
- **Polyhedral refinement and coarsening** — both `refine` and `decimate` raise by name on a polyhedron, pointing at `convert_cells(simplexify)`. Both are built on fixed subdivision templates and an arbitrary polyhedron has none, so closing this means polyhedral agglomeration — a genuinely different algorithm, not another template table. **L**
- **Refinement hierarchy across passes** — `refine:level` exists per pass; a persistent parent/child hierarchy is what multigrid and green-element undo need. **M**
- **Error-estimator helpers** — now that `gradient` exists, a gradient-jump or recovery-based indicator that feeds `refine`'s selection directly closes the adaptive loop end to end. **M**

---

## 5. Field capability beyond derivatives

- **Conservative (mass-preserving) interpolation** — `interpolate`'s barycentric mode is pointwise; CFD remapping needs conservation. **L**
- **Field integration** — total, mean, and per-region reductions over cells as a `data` verb; the natural companion to `gradient`. **S**
- **Second derivatives / Hessian**, for curvature-based adaptivity. **M**

---

## 6. Scale

The benchmark is a ~52k-node bracket; nothing addresses meshes that do not fit in RAM.

- **A large-mesh benchmark tier** (10M+ cells) — cheap, and it would show whether the parallel paths actually hold. Do this before the two below, since it decides whether they matter. **S**
- **Streaming / chunked writes**, the counterpart to the selective-read work. **L**
- **Out-of-core operations** for the ops that are already block-local. **XL**

---

## 7. Ecosystem reach

- **Blender add-on** — Blender ships Python and reads almost no FEA formats; unusually high visibility per line of code. **S–M**
- **Rust bindings** over the C API — the next language by scientific adoption after Julia/R, and the ABI/`SOVERSION` work makes it cheap. **M**
- **Registration and distribution** — conda-forge, CRAN, Julia General, a proper ParaView reader plugin. All deferred at binding time; all pure logistics, and all blocking real adoption. **M**

---

## 8. Quality of implementation

- **Fuzzing the readers** (libFuzzer / AFL, OSS-Fuzz if it will take the project). 41 mostly hand-rolled parsers, reachable from a C ABI, a browser and an MCP server — untrusted input reaches them by design. The highest-value non-feature item in this document. **M**
- **A format conformance matrix** — one canonical mesh written to and read back from every format, with declared per-format lossiness, generalising the region round-trip test into executable documentation of what survives what. **M**
- **Property-based testing** (Hypothesis) over the invariants already articulated in the docs: partition-of-unity, volume conservation, conformity, byte-identical determinism. **M**

---

## 9. NURBS and higher-order geometry (long run)

**The gap.** The data model is strictly linear/Lagrange polytopes: a `CellBlock` is a cell-type string plus a node-index array. NURBS is a genuinely different object — control points, weights, knot vectors, and a parametric mapping — and CAD/IGA formats (STEP, IGES, Rhino 3dm, `.iga`) express geometry that no current cell type can hold. This is the most architecturally invasive item on the list and should be approached as a research spike, not a feature.

- **Spike: how far can the current model stretch?** Higher-order Lagrange cells already exist (`hexahedron27`, VTK-Lagrange types); a rational Bézier/NURBS patch needs *weights* and a *knot vector*, which have nowhere to live. Determine whether a side-channel struct (the `MedInfo`/`GmshInfo` precedent) suffices, or whether the `Mesh` needs a genuine second entity kind. Write the finding up before committing. **M**
- **Read-only CAD ingestion first**: a NURBS surface tessellated to a triangle mesh at a requested tolerance, with the parametric data carried out-of-band. This delivers most of the practical value (getting CAD into the mesh world) without touching the data model, and is the natural first release. **L**
- **A real IGA data model** — patches, control nets, weights, knots, trimming curves — plus formats and evaluation. This is XL, likely a separate library or a major version, and should only be attempted if the spike shows real demand.
- **Dependency reality**: robust STEP/IGES parsing effectively means OpenCASCADE, which is a heavyweight LGPL dependency. If ingestion goes ahead, it must follow the KaHIP/Polyscope pattern — strictly optional, off by default, never in the core, licence implications documented. **Findings before code.**

*Recommended posture: spike and document; do not schedule implementation until the spike says what shape it takes.*

---

## 10. Mesh generation

**The gap.** Every operation transforms a mesh you already have; nothing creates one. This is the only empty category in the operations layer.

- **Primitive constructors** — `box`, `sphere`, `cylinder`, `disk`. (`grid(nx,ny,nz)` shipped in v9.24.0 as part of §1's lattice work, over the same `detail/grid_lattice.hpp`; the rest follow the same shape.) Trivial, dependency-free, and it removes the fixture-file dependency from tests, docs, notebooks, the browser demo and the MCP server. Highest leverage per line of code in this document. **S**
- **`extrude`** — 2D → 3D sweep (triangle→wedge, quad→hexahedron), `nlayers`, per-layer offsets. The most-requested generation primitive; repeatedly deferred. **M**
- **`revolve`** — extrude's rotational sibling, sweeping around an axis. **M**
- **Delaunay / constrained 2D meshing** — genuinely useful, but robust geometric predicates are where dependency-free stops paying. Better as an optional Triangle or Gmsh backend, following the KaHIP pattern. **L**

---

## Suggested sequencing

1. **Primitive constructors (§10, first item)** — a few days, and it improves testing, docs and every demo surface at once.
2. **ML data handling (§2)** — pandas, `edge_index`, and the feature-matrix contract; this is also the prerequisite for §3.
3. **PhysicsNeMo reconnaissance (§3, first item)** — a written findings note before any code.
4. **SDF/octree spike (§1, first item)** — a documented data-model decision before the octree/SDF work is scheduled.
5. **Fuzzing (§8)** — should start in parallel with all of the above; it is not a feature and does not compete for the same attention.
6. **NURBS spike (§9)** — a documented investigation, scheduled independently of the rest.
