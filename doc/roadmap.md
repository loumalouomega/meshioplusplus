# meshio++ roadmap

Status at time of writing: **v9.28.0** — 42 formats, twenty-three mesh operations + five data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers, an MCP server, a settings-driven pipeline engine, a dataset-manifest layer with a PhysicsNeMo adapter, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 6).

This document lists what is *not* built. Items are grouped by theme, each with an effort estimate and the reason it matters. Nothing here duplicates shipped functionality; where a feature partially exists, the gap is stated explicitly.

Effort key: **S** = days, **M** = a couple of weeks, **L** = a month or more, **XL** = a project in its own right.

---

## 1. NVIDIA PhysicsNeMo integration

**What shipped in v9.28.0** (everything except the UI): the reconnaissance note (`doc/physicsnemo.md` — which settled *against* the pip extra this section originally proposed: `nvidia-physicsnemo` hard-depends on torch, the exact wheel the repo's no-`[torch]`-extra precedent refuses to pin), the adapter (`meshioplusplus.physicsnemo`: `graph_sample`/`field_stats`/`edge_stats`/`make_dataset`/`make_reader`), the dataset manager (`DatasetManifest` + the `dataset` CLI group + three MCP tools, `doc/datasets.md`), preprocessing recipes as settings documents, the GPU-executed worked example (`example/physicsnemo/`), and the CI honesty statement. What remains:

- **Dataset-manager UI**: a browser page — a fourth `src/viewer/`-family app, or a mode of the existing one, reusing its MEMFS staging and WASM worker rather than a new render stack — to build and curate a `DatasetManifest` visually: point it at a set of local files/directories, preview each solution (and, for a time series, scrub through its steps) through the existing viewer pipeline, assign splits/tags/group, and edit each case's notes — writing straight back to the same `DatasetManifest` JSON (the format shipped in v9.28.0, `doc/datasets.md`) rather than separate UI-only state, so a session can freely mix hand edits and UI edits without either clobbering the other. This is what makes the dataset manager usable by someone who is not scripting the manifest by hand, and is the natural place to surface per-entry `data_info`/`quality` summaries so a bad case is visible before it corrupts a training split. **M**
- **Adapter follow-ups, recorded with reasons in `doc/physicsnemo.md`**: autoregressive t→t+1 target pairing (`target_fields` is same-step in v1), and a `physicsnemo.mesh.Mesh` bridge — deferred while that type stays simplicial-only with a self-declared-unstable `.pmsh` format; if it stabilizes, an upstream `io_meshio.py` mirroring `io_pyvista`'s shape is the natural form. **S–M**

---

## 2. Remaining refinement and coarsening gaps

`refine` is adaptive (v9.5.0) and `decimate` exists, but the pair still has holes.

- **Volume decimation** — `decimate` is surface-only by documented design; tet-collapse validity is the hard part. **L**
- **Polyhedral refinement and coarsening** — both `refine` and `decimate` raise by name on a polyhedron, pointing at `convert_cells(simplexify)`. Both are built on fixed subdivision templates and an arbitrary polyhedron has none, so closing this means polyhedral agglomeration — a genuinely different algorithm, not another template table. **L**
- **Refinement hierarchy across passes** — `refine:level` exists per pass; a persistent parent/child hierarchy is what multigrid and green-element undo need. **M**
- **Error-estimator helpers** — now that `gradient` exists, a gradient-jump or recovery-based indicator that feeds `refine`'s selection directly closes the adaptive loop end to end. **M**

---

## 3. Field capability beyond derivatives

- **Conservative (mass-preserving) interpolation** — `interpolate`'s barycentric mode is pointwise; CFD remapping needs conservation. **L**
- **Field integration** — total, mean, and per-region reductions over cells as a `data` verb; the natural companion to `gradient`. **S**
- **Second derivatives / Hessian**, for curvature-based adaptivity. **M**

---

## 4. Scale

The benchmark is a ~52k-node bracket; nothing addresses meshes that do not fit in RAM.

- **A large-mesh benchmark tier** (10M+ cells) — cheap, and it would show whether the parallel paths actually hold. Do this before the two below, since it decides whether they matter. **S**
- **Streaming / chunked writes**, the counterpart to the selective-read work. **L**
- **Out-of-core operations** for the ops that are already block-local. **XL**

---

## 5. Ecosystem reach

- **Blender add-on** — Blender ships Python and reads almost no FEA formats; unusually high visibility per line of code. **S–M**
- **Rust bindings** over the C API — the next language by scientific adoption after Julia/R, and the ABI/`SOVERSION` work makes it cheap. **M**
- **Registration and distribution** — conda-forge, CRAN, Julia General, a proper ParaView reader plugin. All deferred at binding time; all pure logistics, and all blocking real adoption. **M**

---

## 6. Quality of implementation

- **Fuzzing the readers** (libFuzzer / AFL, OSS-Fuzz if it will take the project). 42 mostly hand-rolled parsers, reachable from a C ABI, a browser and an MCP server — untrusted input reaches them by design. The highest-value non-feature item in this document. **M**
- **A format conformance matrix** — one canonical mesh written to and read back from every format, with declared per-format lossiness, generalising the region round-trip test into executable documentation of what survives what. **M**
- **Property-based testing** (Hypothesis) over the invariants already articulated in the docs: partition-of-unity, volume conservation, conformity, byte-identical determinism. **M**

---

## 7. NURBS and higher-order geometry (long run)

**The gap.** The data model is strictly linear/Lagrange polytopes: a `CellBlock` is a cell-type string plus a node-index array. NURBS is a genuinely different object — control points, weights, knot vectors, and a parametric mapping — and CAD/IGA formats (STEP, IGES, Rhino 3dm, `.iga`) express geometry that no current cell type can hold. This is the most architecturally invasive item on the list and should be approached as a research spike, not a feature.

- **Spike: how far can the current model stretch?** Higher-order Lagrange cells already exist (`hexahedron27`, VTK-Lagrange types); a rational Bézier/NURBS patch needs *weights* and a *knot vector*, which have nowhere to live. Determine whether a side-channel struct (the `MedInfo`/`GmshInfo` precedent) suffices, or whether the `Mesh` needs a genuine second entity kind. Write the finding up before committing. **M**
- **Read-only CAD ingestion first**: a NURBS surface tessellated to a triangle mesh at a requested tolerance, with the parametric data carried out-of-band. This delivers most of the practical value (getting CAD into the mesh world) without touching the data model, and is the natural first release. **L**
- **A real IGA data model** — patches, control nets, weights, knots, trimming curves — plus formats and evaluation. This is XL, likely a separate library or a major version, and should only be attempted if the spike shows real demand.
- **Dependency reality**: robust STEP/IGES parsing effectively means OpenCASCADE, which is a heavyweight LGPL dependency. If ingestion goes ahead, it must follow the KaHIP/Polyscope pattern — strictly optional, off by default, never in the core, licence implications documented. **Findings before code.**

*Recommended posture: spike and document; do not schedule implementation until the spike says what shape it takes.*

---

## 8. Mesh generation

**The gap.** Every operation transforms a mesh you already have; nothing creates one. This is the only empty category in the operations layer.

- **Primitive constructors** — `box`, `sphere`, `cylinder`, `disk`. (`grid(nx,ny,nz)` shipped in v9.24.0 as part of the signed-distance work, over the same `detail/grid_lattice.hpp`; the rest follow the same shape.) Trivial, dependency-free, and it removes the fixture-file dependency from tests, docs, notebooks, the browser demo and the MCP server. Highest leverage per line of code in this document. **S**
- **`extrude`** — 2D → 3D sweep (triangle→wedge, quad→hexahedron), `nlayers`, per-layer offsets. The most-requested generation primitive; repeatedly deferred. **M**
- **`revolve`** — extrude's rotational sibling, sweeping around an axis. **M**
- **Delaunay / constrained 2D meshing** — genuinely useful, but robust geometric predicates are where dependency-free stops paying. Better as an optional Triangle or Gmsh backend, following the KaHIP pattern. **L**

---

## Suggested sequencing

1. **Primitive constructors (§8, first item)** — a few days, and it improves testing, docs and every demo surface at once. `grid` already shipped over `detail/grid_lattice.hpp`; `box`/`sphere`/`cylinder`/`disk` follow the same shape.
2. **PhysicsNeMo integration (§1)** — shipped in v9.28.0 through the worked example (recon note, adapter, dataset manager, recipes); only the dataset-manager UI remains, once the manifest format has seen some real use.
3. **Fuzzing (§6)** — should start in parallel with all of the above; it is not a feature and does not compete for the same attention.
4. **NURBS spike (§7)** — a documented investigation, scheduled independently of the rest.
