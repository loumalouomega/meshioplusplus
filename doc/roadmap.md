# meshio++ roadmap

Status at time of writing: **v10.11.0** — 42 formats, thirty-two mesh operations + five data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers plus a browser dataset manager, an MCP server, a settings-driven pipeline engine, a dataset-manifest layer with a PhysicsNeMo adapter, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 8).

This document lists what is *not* built. Items are grouped by theme, each with an effort estimate and the reason it matters. Nothing here duplicates shipped functionality; where a feature partially exists, the gap is stated explicitly.

Effort key: **S** = days, **M** = a couple of weeks, **L** = a month or more, **XL** = a project in its own right.

---

## 1. Remeshing (approximated centroidal Voronoi diagrams)

**The gap.** Every resolution-changing operation shipped so far works *on the input's own triangulation*: `refine` subdivides its cells in place, `decimate`/`decimate_volume` collapse its edges, `subdivide`/`agglomerate` restructure it into polyhedra, `smooth` moves its points and touches nothing else. None of them can *retriangulate* — produce a new, near-uniformly-sized, well-shaped triangulation of the same surface at a caller-chosen vertex count. That is a different capability, not a tuning of decimation: QEM edge collapse inherits the input's element shapes and can only remove elements, so a badly-shaped input stays badly-shaped at every target count and there is no way to *raise* the quality of a surface mesh at all. The reference method is Valette's ACVD (`github.com/valette/ACVD`) — discrete approximated centroidal Voronoi clustering over the mesh's own elements, then a dual triangulation of the clusters — which is the standard answer for isotropic surface remeshing, and whose MIT reimplementation `pyvista/pyacvd` confirms the algorithm ports cleanly out of VTK.

**Closed in v10.10.0**: the isotropic CVD core and the ACVDQ feature-preserving metric — `remesh(mesh, RemeshOptions{...})` (`operations/remesh.hpp`, [`doc/remesh.md`](doc/remesh.md)) is a clean-room isotropic clustering engine derived from the MIT `pyvista/pyacvd`, plus a `metric="quadric"` mode (no pyacvd counterpart, synthesized from meshio++'s own pre-existing Garland-Heckbert quadric machinery) that adds a compactness-stabilized quadric term to the clustering objective and places each dual vertex at its cluster's quadric-optimal point — shipped across every binding surface (Python, C API, Fortran, Julia, R, WASM, both CLIs, the settings pipeline, MCP) in one release, unlike `decimate_volume`'s deferred-bindings precedent.

**Closed in v10.11.0**: curvature gradation and boundaries/output manifoldness, both extending `remesh` rather than adding a new operation. Curvature gradation (`RemeshOptions::mGradation`, the exponent `γ` in the item weight `area · κ^γ`) is a genuine new surface-curvature estimator — a local osculating-paraboloid fit over each vertex's 1-ring, principal curvatures from the eigenvalues of the resulting 2×2 shape-operator matrix — reusing the clustering's own `detail::NodeAdjacency` graph rather than adding new neighbourhood machinery; `γ = 0` (default) skips curvature entirely and reproduces plain area weighting byte-for-byte. Boundaries/manifoldness closed both of its stated halves: an open surface's boundary vertices are now seeded before the interior BFS and pinned by name (`RemeshOptions::mPreserveBoundary`, default true), with a second, optional `line` dual cell block emitted along the boundary; and non-manifold "bowtie" output vertices are now detected and reported separately from disconnected clusters (`RemeshResult::mNumNonManifoldVertices`, distinct from `mNumIsolatedClusters`), both repaired by the same regrow-and-reminimise loop. This is a clean-room design achieving the roadmap's stated goals, not a reproduction of ACVD's own boundary-fixing algorithm. Topology preservation remains explicitly out of scope and documented as such (thin features can still merge; the genus can still change). See `doc/remesh.md`.

- **Anisotropic metric** — clusters shaped by a local curvature tensor rather than isotropic distance, so elongated features are meshed with elongated elements at a fraction of the vertex count. The largest quality win and the most delicate numerically; now unblocked, since curvature gradation's density-weighting and curvature-estimation machinery (v10.11.0) is exactly what it was sequenced to share. **L**
- **Volumetric CVD/ODT** — the tetrahedral counterpart (optimal Delaunay triangulation / CVD in 3D), which would do for `decimate_volume` what everything above does for `decimate`. A genuinely separate algorithm needing robust Delaunay predicates, i.e. exactly where dependency-free stops paying — the same trade §7's last bullet already records. Not a follow-on task; a project in its own right. **XL**

*Recommended posture: the two remaining bullets are independent of each other and of the shipped core — the anisotropic metric is a quality refinement layered on the existing `RemeshMetric` dispatch and curvature machinery, and volumetric CVD/ODT is a separate project in its own right with no shared code to reuse beyond the pattern.*

---

## 2. Dataset dashboard and training integration

**The gap.** The browser dataset manager (`src/viewer/`, `dataset.html`, v9.29.0) curates a single `DatasetManifest` at a time — add/list/split/tag entries, stage one entry into MEMFS and preview it (per-entry quality summaries landed in v9.30.0). It has no aggregate view: there is no way to see many datasets, or many entries within one, at a glance, and once a dataset is curated there is no path from "manifest is ready" to "a PhysicsNeMo run is training against it" without leaving the browser entirely for the CLI/Python recipe in `example/physicsnemo/`. This is a UI/workflow gap, not a numerical one — the underlying capability (`edge_index`, `feature_matrix`, `write_dataset`, dataset manifests, the `physicsnemo.mesh` bridge) is complete as of v9.30.0; nothing surfaces it as one connected experience.

- **Multi-dataset overview** — a landing dashboard listing every known manifest as a summary card (entry count, split sizes, tags/groups, last-modified, a thumbnail from one representative entry) instead of today's single-manifest-at-a-time view. **M**
- **Per-dataset drill-down** — clicking a card opens the existing entry list/preview flow scoped to that one manifest, so the aggregate and detail views become two depths of one page rather than two disconnected surfaces. **S** (mostly re-plumbing `dataset.html`'s existing state around a selected-manifest concept)
- **Launch PhysicsNeMo training from the dashboard** — a "train" action against a manifest/split that kicks off training (generalizing the `example/physicsnemo/` recipe). This is the one item here that crosses an architectural line the rest of the viewer never has: everything else runs client-side in WASM with no server, and real GPU training cannot. It needs a small companion process — a local service exposing job start/status/logs, a natural extension of the existing `meshioplusplus-mcp` server rather than a new protocol — that the page talks to; **that server-side dependency should be designed and documented explicitly, not treated as a detail.** **L**
- **Training monitoring** — once a job exists, live loss/metric curves (train/val), progress/ETA, and a stop control, polling or streaming from the companion process above. **M**
- **Run history and comparison** — a sortable/filterable table of past runs (hyperparameters, dataset/split, final metrics), so a new run can be compared against prior ones without leaving the page. **M**
- **Log tailing** — raw stdout/stderr from a running or finished job, for the moment a metric alone doesn't explain a failure. **S**
- **Checkpoint browser** — list, download, and mark-as-best the checkpoints a run produced. **S**
- **Prediction preview in-viewer** — run a checkpoint over a held-out entry and show the predicted field (and its error against ground truth) through the *existing* mesh viewer's colour-by machinery — the PhysicsNeMo example's `T_pred`/`T_error` write-back, made interactive instead of a one-off script. **M**
- **Dataset health summaries** — per-dataset checks on the card/drill-down: split balance, fields missing across entries, degenerate/inverted-cell counts (reusing `compute_quality`) — so a bad dataset is visible before a run wastes GPU time on it. **S**
- **Manifest diffing/versioning** — manifests are hand-editable JSON (`doc/datasets.md`); a lightweight diff view between two versions of the same file (or two git revisions of it) would make manual edits auditable. **S**
- **Run-completion notifications** — a browser notification, or a webhook the companion process posts to, when a launched run finishes or fails, so the dashboard need not stay the active tab. **S**
- **Visual/UX design pass** — once the layout above is functional, loop it through Claude Design (Claude in a design-iteration capacity — layout, spacing, colour, motion) rather than shipping the first working arrangement as the final one; the rest of the viewer already has a deliberate icon set and colour system (`doc/icons/`, the `dataviz` skill's palette) and this should read as one piece with it, not a bolted-on admin panel. **S**

*Recommended posture: the dashboard/drill-down/health-summary items are ordinary viewer work and can proceed independently and incrementally. Training launch and monitoring is the one item requiring a server-side companion process and should get its own short design pass — what talks to what, auth, where jobs actually run — before implementation, the same way the NURBS spike (§6) is scoped before its own implementation. Once the pieces work, run a design-polish loop (Claude Design) over the whole dashboard before calling it done — functional and pleasant are two different bars.*

---

## 3. Scale

The benchmark is a ~52k-node bracket; nothing addresses meshes that do not fit in RAM.

- **A large-mesh benchmark tier** (10M+ cells) — cheap, and it would show whether the parallel paths actually hold. Do this before the two below, since it decides whether they matter. **S**
- **Streaming / chunked writes**, the counterpart to the selective-read work. **L**
- **Out-of-core operations** for the ops that are already block-local. **XL**

---

## 4. Ecosystem reach

- **Blender add-on** — Blender ships Python and reads almost no FEA formats; unusually high visibility per line of code. **S–M**
- **Rust bindings** over the C API — the next language by scientific adoption after Julia/R, and the ABI/`SOVERSION` work makes it cheap. **M**
- **Registration and distribution** — conda-forge, CRAN, Julia General, a proper ParaView reader plugin. All deferred at binding time; all pure logistics, and all blocking real adoption. **M**

---

## 5. Quality of implementation

- **Fuzzing the readers** (libFuzzer / AFL, OSS-Fuzz if it will take the project). 42 mostly hand-rolled parsers, reachable from a C ABI, a browser and an MCP server — untrusted input reaches them by design. The highest-value non-feature item in this document. **M**
- **A format conformance matrix** — one canonical mesh written to and read back from every format, with declared per-format lossiness, generalising the region round-trip test into executable documentation of what survives what. **M**
- **Property-based testing** (Hypothesis) over the invariants already articulated in the docs: partition-of-unity, volume conservation, conformity, byte-identical determinism. **M**

---

## 6. NURBS and higher-order geometry (long run)

**The gap.** The data model is strictly linear/Lagrange polytopes: a `CellBlock` is a cell-type string plus a node-index array. NURBS is a genuinely different object — control points, weights, knot vectors, and a parametric mapping — and CAD/IGA formats (STEP, IGES, Rhino 3dm, `.iga`) express geometry that no current cell type can hold. This is the most architecturally invasive item on the list and should be approached as a research spike, not a feature.

- **Spike: how far can the current model stretch?** Higher-order Lagrange cells already exist (`hexahedron27`, VTK-Lagrange types); a rational Bézier/NURBS patch needs *weights* and a *knot vector*, which have nowhere to live. Determine whether a side-channel struct (the `MedInfo`/`GmshInfo` precedent) suffices, or whether the `Mesh` needs a genuine second entity kind. Write the finding up before committing. **M**
- **Read-only CAD ingestion first**: a NURBS surface tessellated to a triangle mesh at a requested tolerance, with the parametric data carried out-of-band. This delivers most of the practical value (getting CAD into the mesh world) without touching the data model, and is the natural first release. **L**
- **A real IGA data model** — patches, control nets, weights, knots, trimming curves — plus formats and evaluation. This is XL, likely a separate library or a major version, and should only be attempted if the spike shows real demand.
- **Dependency reality**: robust STEP/IGES parsing effectively means OpenCASCADE, which is a heavyweight LGPL dependency. If ingestion goes ahead, it must follow the KaHIP/Polyscope pattern — strictly optional, off by default, never in the core, licence implications documented. **Findings before code.**

*Recommended posture: spike and document; do not schedule implementation until the spike says what shape it takes.*

---

## 7. Mesh generation

**The gap.** Every operation transforms a mesh you already have; nothing creates one. This is the only empty category in the operations layer.

- **Primitive constructors** — `box`, `sphere`, `cylinder`, `disk`. (`grid(nx,ny,nz)` shipped in v9.24.0 as part of the signed-distance work, over the same `detail/grid_lattice.hpp`; the rest follow the same shape.) Trivial, dependency-free, and it removes the fixture-file dependency from tests, docs, notebooks, the browser demo and the MCP server. Highest leverage per line of code in this document. **S**
- **`extrude`** — 2D → 3D sweep (triangle→wedge, quad→hexahedron), `nlayers`, per-layer offsets. The most-requested generation primitive; repeatedly deferred. **M**
- **`revolve`** — extrude's rotational sibling, sweeping around an axis. **M**
- **Delaunay / constrained 2D meshing** — genuinely useful, but robust geometric predicates are where dependency-free stops paying. Better as an optional Triangle or Gmsh backend, following the KaHIP pattern. **L**

---

## 8. CLI chatbot / conversational assistant (MCP-driven)

**The gap.** The MCP server (`src/python/meshioplusplus/mcp/`, `doc/mcp.md`) exposes the whole Python surface to an AI *agent* — but only to one already speaking MCP over stdio (Claude Desktop, an IDE, a custom host). There is no way to have a natural-language conversation about a mesh **from the terminal itself**: a user with an LLM API key on hand cannot ask `meshioplusplus` "why does this file fail to convert" or "clean this mesh and tell me what changed" and get a tool-calling assistant that drives the existing operations for them. Every other surface (Python API, CLI verbs, MCP tools) is imperative-only; this is the one conversational entry point missing.

- **Minimal REPL verb**: a `meshioplusplus chat` command (Python CLI only, `_cli/_chat.py`, alongside `data`/`dataset` as a third nested concern) that starts a terminal chat loop. Gated cleanly on an LLM API key being present (`ANTHROPIC_API_KEY` first-class, matching the `claude-api` skill's model; provider-agnostic wiring behind the same interface is a stretch goal, not a v1 requirement) — absent, it fails by name (`pip install meshioplusplus[chat]` / "set ANTHROPIC_API_KEY"), the `meshioplusplus-mcp` entry-point precedent, never a bare traceback. **S–M**
- **Reuse, not reimplementation**: the assistant's tool-calling loop drives `mcp/_tools.py`'s existing `TOOL_REGISTRY` directly (in-process, no stdio round trip needed for a same-process CLI) — never a second copy of the tool dispatch, sandbox (`_resolve`), or `_json_safe` result-sanitizing logic MCP already owns. A new tool added to `TOOL_REGISTRY` is automatically available to the chatbot for free, the same parity guarantee `test_every_operation_has_a_tool` already gives the MCP surface. **M** (mostly wiring an agentic loop over an existing, stable tool registry — see the `claude-api` skill for the Messages API tool-use shape)
- **Scope for v1**: a stateless-per-turn loop (send message + conversation history + tool schemas from `TOOL_REGISTRY`, execute any tool calls the model requests, feed results back, repeat until a plain-text reply) operating on files under the CWD or an explicit `--root`, mirroring the MCP server's own path sandbox rather than inventing a second one. Multi-turn context lives only in the terminal session; no persistence, no server. **S**
- **New optional extra**: `chat = ["anthropic>=0.40,<1"]` (the SDK, not a hard dependency of `meshioplusplus` itself — the `mcp` extra's precedent), kept out of `[all]`. `_tools.py` stays untouched (no SDK import); only the new `_cli/_chat.py` / a `mcp/_chat.py` module imports it, so the default CI matrix and every other surface are unaffected. **S**
- **Docs**: a `doc/chat.md` page (modelled on `doc/mcp.md`) and a CLAUDE.md entry once shipped, per the "Keeping docs in sync" rule at the top of this file.

*Recommended posture: this is a thin client over work that already exists (the MCP tool registry) rather than new mesh functionality — the honest estimate for a usable v1 is small (S–M), with provider-agnostic support and a persisted chat history as documented follow-ups rather than v1 requirements.*

---

## Suggested sequencing

1. **Isotropic CVD remeshing, the ACVDQ feature-preserving metric, curvature gradation and boundaries/manifoldness (§1, first four items)** — shipped in full and removed from this document: v10.10.0 (`remesh`, a clean-room isotropic clustering engine derived from the MIT `pyvista/pyacvd` plus an original `metric="quadric"` synthesis over meshio++'s own pre-existing Garland-Heckbert quadric machinery, across every binding surface in one release); v10.11.0 (curvature gradation via a new osculating-paraboloid vertex-curvature estimator, plus boundary pinning/dual-edge insertion and separately-reported non-manifold-vertex detection, both extending `remesh` and its full binding surface). The anisotropic metric (now unblocked) and the volumetric counterpart remain open in §1. See `doc/remesh.md`.
2. **Primitive constructors (§7, first item)** — a few days, and it improves testing, docs and every demo surface at once. `grid` already shipped over `detail/grid_lattice.hpp`; `box`/`sphere`/`cylinder`/`disk` follow the same shape.
3. **PhysicsNeMo integration** — shipped in full and removed from this document: v9.28.0 (recon note, adapter, dataset manager, recipes, GPU-executed example), v9.29.0 (dataset-manager UI), v9.30.0 (t→t+1 target pairing, the `physicsnemo.mesh.Mesh` bridge, persisted directory handles, per-entry quality summaries). See `doc/physicsnemo.md` and `doc/datasets.md`.
4. **Remaining refinement and coarsening gaps** — shipped in full and removed from this document: v10.2.0 (error-estimator helpers — `estimate_error`), v10.3.0 (polyhedral refinement — `subdivide`), v10.4.0 (polyhedral coarsening — `agglomerate`), v10.5.0 (green-element undo — `undo_green`), v10.6.0 (volume decimation — `decimate_volume`, the section's last open item). See `doc/error.md`, `doc/subdivide.md`, `doc/agglomerate.md`, `doc/undo_green.md` and `doc/decimate_volume.md`.
5. **Field capability beyond derivatives** — shipped in full and removed from this document: v10.8.0 (field integration — `data_integrate`, a cell-measure-weighted total/mean over cells, whole-mesh and per named Cell region), v10.9.0 (second derivatives — `hessian`, a composition of two `gradient` calls, exact for a linear field everywhere and for a quadratic field away from a structured mesh's own boundary, closing the section's last open item). See `doc/field_integration.md` and `doc/hessian.md`.
6. **Fuzzing (§5)** — should start in parallel with all of the above; it is not a feature and does not compete for the same attention.
7. **NURBS spike (§6)** — a documented investigation, scheduled independently of the rest.
8. **Dataset dashboard (§2)** — the non-training pieces (multi-dataset overview, drill-down, health summaries, diffing) can proceed any time; training launch/monitoring waits on its own design pass (server-side companion process) before implementation.
9. **CLI chatbot (§8)** — small and self-contained (a thin client over the existing MCP tool registry); can proceed independently whenever a maintainer wants it, no sequencing dependency on anything above.
