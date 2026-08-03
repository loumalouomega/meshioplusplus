# meshio++ roadmap

Status at time of writing: **v9.14.0** — 41 formats, twenty mesh operations + five data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers, an MCP server, a settings-driven pipeline engine, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 5).

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

**The gap.** Ragged polyhedron blocks exist in all three backends (`AddPolyhedronBlock`, CSR / nested storage, `CellView::NumFaces`/`Face`), and MED, EnSight `nfaced` and OpenFOAM read them. But they are second-class almost everywhere else: the C ABI reports `is_ragged` and then **cannot expose the connectivity at all**, most operations raise on them, and the geometric kernel (`cell_faces.hpp`) is a fixed table of canonical types with no polyhedral entry. Since OpenFOAM — the most-used open CFD code — is natively polyhedral, this is a real ceiling.

- **Polyhedral connectivity across the C ABI** (and therefore Fortran / Julia / R / WASM). The documented v1 gap; needs a face-offsets + face-nodes CSR pair on the flat surface. Everything else in this section depends on it. **M**
- **Geometric kernel for arbitrary polyhedra**: volume and centroid by face-fan decomposition about an interior point, face normal/area for arbitrary planar-ish polygons (Newell), and a documented answer for non-planar faces. Unblocks `stats`, `quality` (a reduced metric set), `smooth`'s inversion guard, `gradient` (Green–Gauss is *naturally* polyhedral — it only needs faces), and `data_average`'s measure weighting. **M**
- **Operations that should accept polyhedra**: `gradient`, `stats`, `crop`/`slice` (cut a polyhedron → still a polyhedron), `merge`, `partition` (the dual graph is face-based already), `interpolate` (nearest works; barycentric needs the tet decomposition below). **L**
- **Polyhedron → tetrahedra decomposition** in `convert_cells simplexify` (face-fan about the centroid, adding one point per cell), which is the escape hatch that makes every remaining op work. Cheaper than making each op polyhedral, and honest. **M**
- **Polyhedral writers**: OpenFOAM (currently read-only), VTU `polyhedron` cells, MED `POG`. Round-tripping an OpenFOAM case is the headline outcome. **M–L**

*Recommended entry point: the C-ABI exposure plus the geometric kernel — they convert polyhedra from "stored" to "usable" and unblock everything else.*

---

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

- **Primitive constructors** — `box`, `grid(nx,ny,nz)`, `sphere`, `cylinder`, `disk`. Trivial, dependency-free, and it removes the fixture-file dependency from tests, docs, notebooks, the browser demo and the MCP server. Highest leverage per line of code in this document. **S**
- **`extrude`** — 2D → 3D sweep (triangle→wedge, quad→hexahedron), `nlayers`, per-layer offsets. The most-requested generation primitive; repeatedly deferred. **M**
- **`revolve`** — extrude's rotational sibling, sweeping around an axis. **M**
- **Delaunay / constrained 2D meshing** — genuinely useful, but robust geometric predicates are where dependency-free stops paying. Better as an optional Triangle or Gmsh backend, following the KaHIP pattern. **L**

---

## Suggested sequencing

1. **Primitive constructors (§10, first item)** — a few days, and it improves testing, docs and every demo surface at once.
2. **ML data handling (§2)** — pandas, `edge_index`, and the feature-matrix contract; this is also the prerequisite for §3.
3. **PhysicsNeMo reconnaissance (§3, first item)** — a written findings note before any code.
4. **Polyhedral C-ABI exposure + geometric kernel (§1)** — lifts a real ceiling and unblocks OpenFOAM round-tripping.
5. **Fuzzing (§8)** — should start in parallel with all of the above; it is not a feature and does not compete for the same attention.
6. **NURBS spike (§9)** — a documented investigation, scheduled independently of the rest.
