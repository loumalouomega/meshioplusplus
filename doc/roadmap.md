# meshio++ roadmap

Status at time of writing: **v10.5.0** — 42 formats, twenty-seven mesh operations + five data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers plus a browser dataset manager, an MCP server, a settings-driven pipeline engine, a dataset-manifest layer with a PhysicsNeMo adapter, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 7).

This document lists what is *not* built. Items are grouped by theme, each with an effort estimate and the reason it matters. Nothing here duplicates shipped functionality; where a feature partially exists, the gap is stated explicitly.

Effort key: **S** = days, **M** = a couple of weeks, **L** = a month or more, **XL** = a project in its own right.

---

## 1. Remaining refinement and coarsening gaps

`refine` is adaptive (v9.5.0) and `decimate` exists, but the pair still has holes.

- **Volume decimation** — `decimate` is surface-only by documented design; tet-collapse validity is the hard part. **L**

**Closed in v10.2.0**: the error-estimator-helpers gap — `estimate_error` (`operations/error.hpp`, [`doc/error.md`](doc/error.md)) is the Zienkiewicz-Zhu recovery-based indicator, composed from `gradient` and the existing point↔cell averaging round trip, plus absolute/fraction/Dörfler marking into a `error:marked` array `refine`'s own `--where` selector consumes directly with no change to `refine` — closing the adaptive loop (estimate → mark → refine) end to end across every binding surface.

**Closed in v10.3.0**: the polyhedral-*refinement* half of the former "polyhedral coarsening" gap — `subdivide` (`operations/subdivide.hpp`, [`doc/subdivide.md`](doc/subdivide.md)) splits every eligible 3D cell (tabulated or an existing polyhedron block, handled uniformly through `detail::cell_rings`/`orient_rings` — no per-type template table at all) into one polyhedral child per face, connected to a new interior point. Automatically conforming (a shared face is never touched, so there is no closure/hanging-node bookkeeping to do), and shipped across every binding surface.

**Closed in v10.4.0**: the polyhedral-*coarsening* (agglomeration) half of that same gap, now closed in full — `agglomerate` (`operations/agglomerate.hpp`, [`doc/agglomerate.md`](doc/agglomerate.md)) merges groups of cells into single larger polyhedral cells via greedy seed-and-grow over the mesh's shared-face dual (`detail::build_global_faces`'s owner/neighbour pairing), absorbing face-adjacent neighbours by accumulated shared-face area until a target group size. Each group emits one polyhedron whose faces are exactly its external boundary, conserving volume exactly (an identity of surviving faces, not a divergence-theorem coincidence). Regions carry through `CellMapKind::Global`, a single flat input-cell→output-cell map — simpler than `merge`'s own per-input-mesh usage of the same map kind, since agglomerate has exactly one input mesh. Shipped across every binding surface.

**Closed in v10.5.0**: green-element undo — `undo_green` (`operations/undo_green.hpp`, [`doc/undo_green.md`](doc/undo_green.md)) restores `refine`'s transitional (green) cells back to their coarse parent before a new selective pass, closing the quality-degradation issue this section had documented since v10.1.0. The originally-planned mechanism — inverting `refine_templates.cpp`'s subdivision tables against a sibling run — turned out to be unnecessary: since `refine()`'s point map is always the identity, a green parent's exact connectivity and cell_data are already sitting, byte-for-byte, in the caller-supplied **coarse** mesh (the mesh the `refine()` call that produced `fine` was run on), resolved via `refine:parent_id` against `coarse`'s own `refine:cell_id` (or its implicit id when absent). `undo_green(coarse, fine)` is therefore a lookup-and-substitution, not a reconstruction — no template inversion, no winding repair, no discrete sign branch — which is also what gives it a full numpy twin, unlike `subdivide`/`agglomerate`. Cell regions carry through the first genuinely non-injective `CellMapKind::Direct` use in the C++ core (several fine cells collapsing onto one output row), deduplicated by `Region`'s existing sort+unique. Shipped across every binding surface, as a two-mesh operation mirroring `interpolate`'s own shape (module-level Fortran, excluded from the settings pipeline like `Merge`/`Interpolate`/`Split`/`Diff`). `decimate` (surface QEM) and volume decimation, above, remain the one still-open piece of this document's refinement/coarsening section.

---

## 2. Field capability beyond derivatives

- **Conservative (mass-preserving) interpolation** — `interpolate`'s barycentric mode is pointwise; CFD remapping needs conservation. **L**
- **Field integration** — total, mean, and per-region reductions over cells as a `data` verb; the natural companion to `gradient`. **S**
- **Second derivatives / Hessian**, for curvature-based adaptivity. **M**

---

## 3. Dataset dashboard and training integration

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

*Recommended posture: the dashboard/drill-down/health-summary items are ordinary viewer work and can proceed independently and incrementally. Training launch and monitoring is the one item requiring a server-side companion process and should get its own short design pass — what talks to what, auth, where jobs actually run — before implementation, the same way the NURBS spike (§7) is scoped before its own implementation. Once the pieces work, run a design-polish loop (Claude Design) over the whole dashboard before calling it done — functional and pleasant are two different bars.*

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

## 9. CLI chatbot / conversational assistant (MCP-driven)

**The gap.** The MCP server (`src/python/meshioplusplus/mcp/`, `doc/mcp.md`) exposes the whole Python surface to an AI *agent* — but only to one already speaking MCP over stdio (Claude Desktop, an IDE, a custom host). There is no way to have a natural-language conversation about a mesh **from the terminal itself**: a user with an LLM API key on hand cannot ask `meshioplusplus` "why does this file fail to convert" or "clean this mesh and tell me what changed" and get a tool-calling assistant that drives the existing operations for them. Every other surface (Python API, CLI verbs, MCP tools) is imperative-only; this is the one conversational entry point missing.

- **Minimal REPL verb**: a `meshioplusplus chat` command (Python CLI only, `_cli/_chat.py`, alongside `data`/`dataset` as a third nested concern) that starts a terminal chat loop. Gated cleanly on an LLM API key being present (`ANTHROPIC_API_KEY` first-class, matching the `claude-api` skill's model; provider-agnostic wiring behind the same interface is a stretch goal, not a v1 requirement) — absent, it fails by name (`pip install meshioplusplus[chat]` / "set ANTHROPIC_API_KEY"), the `meshioplusplus-mcp` entry-point precedent, never a bare traceback. **S–M**
- **Reuse, not reimplementation**: the assistant's tool-calling loop drives `mcp/_tools.py`'s existing `TOOL_REGISTRY` directly (in-process, no stdio round trip needed for a same-process CLI) — never a second copy of the tool dispatch, sandbox (`_resolve`), or `_json_safe` result-sanitizing logic MCP already owns. A new tool added to `TOOL_REGISTRY` is automatically available to the chatbot for free, the same parity guarantee `test_every_operation_has_a_tool` already gives the MCP surface. **M** (mostly wiring an agentic loop over an existing, stable tool registry — see the `claude-api` skill for the Messages API tool-use shape)
- **Scope for v1**: a stateless-per-turn loop (send message + conversation history + tool schemas from `TOOL_REGISTRY`, execute any tool calls the model requests, feed results back, repeat until a plain-text reply) operating on files under the CWD or an explicit `--root`, mirroring the MCP server's own path sandbox rather than inventing a second one. Multi-turn context lives only in the terminal session; no persistence, no server. **S**
- **New optional extra**: `chat = ["anthropic>=0.40,<1"]` (the SDK, not a hard dependency of `meshioplusplus` itself — the `mcp` extra's precedent), kept out of `[all]`. `_tools.py` stays untouched (no SDK import); only the new `_cli/_chat.py` / a `mcp/_chat.py` module imports it, so the default CI matrix and every other surface are unaffected. **S**
- **Docs**: a `doc/chat.md` page (modelled on `doc/mcp.md`) and a CLAUDE.md entry once shipped, per the "Keeping docs in sync" rule at the top of this file.

*Recommended posture: this is a thin client over work that already exists (the MCP tool registry) rather than new mesh functionality — the honest estimate for a usable v1 is small (S–M), with provider-agnostic support and a persisted chat history as documented follow-ups rather than v1 requirements.*

---

## Suggested sequencing

1. **Primitive constructors (§8, first item)** — a few days, and it improves testing, docs and every demo surface at once. `grid` already shipped over `detail/grid_lattice.hpp`; `box`/`sphere`/`cylinder`/`disk` follow the same shape.
2. **PhysicsNeMo integration** — shipped in full and removed from this document: v9.28.0 (recon note, adapter, dataset manager, recipes, GPU-executed example), v9.29.0 (dataset-manager UI), v9.30.0 (t→t+1 target pairing, the `physicsnemo.mesh.Mesh` bridge, persisted directory handles, per-entry quality summaries). See `doc/physicsnemo.md` and `doc/datasets.md`.
3. **Fuzzing (§6)** — should start in parallel with all of the above; it is not a feature and does not compete for the same attention.
4. **NURBS spike (§7)** — a documented investigation, scheduled independently of the rest.
5. **Dataset dashboard (§3)** — the non-training pieces (multi-dataset overview, drill-down, health summaries, diffing) can proceed any time; training launch/monitoring waits on its own design pass (server-side companion process) before implementation.
6. **CLI chatbot (§9)** — small and self-contained (a thin client over the existing MCP tool registry); can proceed independently whenever a maintainer wants it, no sequencing dependency on anything above.
