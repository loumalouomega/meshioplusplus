# meshio++ roadmap

Status at time of writing: **v10.21.0** — 43 readable and 46 writable formats, 34 mesh operations and 5 data operations, six language surfaces (Python / C / Fortran / Julia / R / WebAssembly), two viewers plus a browser dataset manager, a Blender add-on, a ParaView reader plugin, an MCP server with 57 tools, a settings-driven pipeline engine with a transient sequence driver, a dataset-manifest layer with a PhysicsNeMo adapter, provenance in written files, and a versioned C++ ABI (`MESHIOPLUSPLUS_ABI_VERSION` 11).

This document lists what is *not* built. Items are grouped by theme, each with an effort estimate and the reason it matters. Nothing here duplicates shipped functionality; where a feature partially exists, the gap is stated explicitly, and when a change closes or narrows an item the entry moves to the [recently closed](#recently-closed) table in the same change.

Effort key: **S** = days, **M** = a couple of weeks, **L** = a month or more, **XL** = a project in its own right.

## The map

![The roadmap at a glance: open items grouped by theme, shaded by effort, with dependency arrows and the items that need a design pass or a research spike first](/diagrams/roadmap_map.svg)

| Section | Open items | Effort | Depends on | Posture |
| --- | --- | --- | --- | --- |
| [§1 Dataset dashboard and training](#1-dataset-dashboard-and-training-integration) | 12 | S–L | a server-side companion process for anything that trains | proceed on the client-side pieces; design pass first for training |
| [§2 Scale](#2-scale) | 3 | S–XL | the benchmark tier decides whether the other two matter | benchmark first |
| [§3 Ecosystem reach](#3-ecosystem-reach) | 2 | M | nothing technical; registration is logistics | proceed |
| [§4 Quality of implementation](#4-quality-of-implementation) | 3 | M | nothing | proceed, in parallel with everything |
| [§5 NURBS and higher-order geometry](#5-nurbs-and-higher-order-geometry-long-run) | 3 | M–XL | the spike's findings | research spike first |
| [§6 Mesh generation](#6-mesh-generation) | 4 | S–L | primitives before extrude, extrude before revolve | proceed, primitives first |
| [§7 CLI chatbot](#7-cli-chatbot--conversational-assistant-mcp-driven) | 1 | S–M | the existing MCP tool registry | proceed whenever wanted |
| [§8 Binding-surface parity](#8-binding-surface-parity) | 2 | S–M | nothing | proceed |
| [§9 Format completeness](#9-format-completeness) | 5 | S–M | nothing | proceed |
| [§10 Follow-ups recorded elsewhere](#10-follow-ups-recorded-in-other-pages) | 6 | S–M | see each host page | pointers only |

---

## 1. Dataset dashboard and training integration

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
- **Visual/UX design pass** — once the layout above is functional, loop it through a design-iteration pass (layout, spacing, colour, motion) rather than shipping the first working arrangement as the final one; the rest of the viewer already has a deliberate icon set and colour system (`doc/icons/`, the palette `doc/diagrams/` also uses) and this should read as one piece with it, not a bolted-on admin panel. **S**

*Recommended posture: the dashboard/drill-down/health-summary items are ordinary viewer work and can proceed independently and incrementally. Training launch and monitoring is the one item requiring a server-side companion process and should get its own short design pass — what talks to what, auth, where jobs actually run — before implementation, the same way the NURBS spike (§5) is scoped before its own implementation.*

---

## 2. Scale

The benchmark is a ~52k-node bracket (293k cells; the backend benchmark's largest input is 257k tetrahedra); nothing addresses meshes that do not fit in RAM.

- **A large-mesh benchmark tier** (10M+ cells) — cheap, and it would show whether the parallel paths actually hold. Do this before the two below, since it decides whether they matter. **S**
- **Streaming / chunked writes of one mesh**, the counterpart to the selective-read work. The [sequence layer](./sequences.md) already streams at *file* granularity — a fan-in, fan-out or per-step run keeps at most one mesh alive, and that invariant is pinned by a test — but nothing streams the cells of a single mesh that does not fit in memory, and every writer takes a whole `Mesh`. **L**
- **Out-of-core operations** for the ops that are already block-local. **XL**

---

## 3. Ecosystem reach

- **Rust bindings** over the C API — the next language by scientific adoption after Julia/R, and the ABI/`SOVERSION` work makes it cheap. **M**
- **Registration and distribution** — conda-forge, CRAN, Julia General (plus a JLL for the binary), the Blender Extensions Platform listing, and a listed ParaView plugin. All deferred at binding time; all pure logistics, and all blocking real adoption. The artefacts themselves exist: the Blender extension zips are attached to every `v*` release (`doc/blender.md`), and the ParaView reader plugin ships as `tools/paraview-meshioplusplus-plugin.py` with its own [install page](./paraview_plugin.md), so what remains on each is the *listing*, not the code. Spack recipes are already upstream. **M**

---

## 4. Quality of implementation

- **Fuzzing the readers** (libFuzzer / AFL, OSS-Fuzz if it will take the project). 43 mostly hand-rolled parsers, reachable from a C ABI, a browser and an MCP server — untrusted input reaches them by design. The highest-value non-feature item in this document. **M**
- **A format conformance matrix** — one canonical mesh written to and read back from every format, with declared per-format lossiness, generalising the region round-trip test into executable documentation of what survives what. **M**
- **Property-based testing** (Hypothesis) over the invariants already articulated in the docs: partition-of-unity, volume conservation, conformity, byte-identical determinism. **M**

---

## 5. NURBS and higher-order geometry (long run)

**The gap.** The data model is strictly linear/Lagrange polytopes: a `CellBlock` is a cell-type string plus a node-index array. NURBS is a genuinely different object — control points, weights, knot vectors, and a parametric mapping — and CAD/IGA formats (STEP, IGES, Rhino 3dm, `.iga`) express geometry that no current cell type can hold. This is the most architecturally invasive item on the list and should be approached as a research spike, not a feature.

- **Spike: how far can the current model stretch?** Higher-order Lagrange cells already exist (`hexahedron27`, VTK-Lagrange types); a rational Bézier/NURBS patch needs *weights* and a *knot vector*, which have nowhere to live. Determine whether a side-channel struct (the `MedInfo`/`GmshInfo` precedent) suffices, or whether the `Mesh` needs a genuine second entity kind. Write the finding up before committing. **M**
- **Read-only CAD ingestion first**: a NURBS surface tessellated to a triangle mesh at a requested tolerance, with the parametric data carried out-of-band. This delivers most of the practical value (getting CAD into the mesh world) without touching the data model, and is the natural first release. **L**
- **A real IGA data model** — patches, control nets, weights, knots, trimming curves — plus formats and evaluation. This is XL, likely a separate library or a major version, and should only be attempted if the spike shows real demand.
- **Dependency reality**: robust STEP/IGES parsing effectively means OpenCASCADE, which is a heavyweight LGPL dependency. If ingestion goes ahead, it must follow the KaHIP/Polyscope pattern — strictly optional, off by default, never in the core, licence implications documented. **Findings before code.**

*Recommended posture: spike and document; do not schedule implementation until the spike says what shape it takes.*

---

## 6. Mesh generation

**The gap.** Every operation transforms a mesh you already have; almost nothing creates one. Two generators exist and neither closes a bullet below: [`grid`](./voxelize.md) (v9.24.0) builds a regular hexahedron lattice from nothing and ships across Python, C++, C, Fortran, Julia, R, WebAssembly and the MCP server, and [`remesh_volume`](./remesh_volume.md) (v10.13.0) accepts a closed surface and generates a genuinely new tetrahedral volume mesh by isosurface stuffing. There are still no primitive constructors, no extrude/revolve, and no 2-D triangulator.

- **Primitive constructors** — `box`, `sphere`, `cylinder`, `disk`, over the same `detail/grid_lattice.hpp` shape `grid` established. Trivial, dependency-free, and it removes the fixture-file dependency from tests, docs, notebooks, the browser demo and the MCP server. Highest leverage per line of code in this document. **S**
- **`extrude`** — 2D → 3D sweep (triangle→wedge, quad→hexahedron), `nlayers`, per-layer offsets. The most-requested generation primitive; repeatedly deferred. **M**
- **`revolve`** — extrude's rotational sibling, sweeping around an axis. **M**
- **Delaunay / constrained 2D meshing** — genuinely useful, but robust geometric predicates are where dependency-free stops paying. Better as an optional Triangle or Gmsh backend, following the KaHIP pattern. **L**

---

## 7. CLI chatbot / conversational assistant (MCP-driven)

**The gap.** The MCP server (`src/python/meshioplusplus/mcp/`, `doc/mcp.md`) exposes the whole Python surface to an AI *agent* — but only to one already speaking MCP over stdio (Claude Desktop, an IDE, a custom host). There is no way to have a natural-language conversation about a mesh **from the terminal itself**: a user with an LLM API key on hand cannot ask `meshioplusplus` "why does this file fail to convert" or "clean this mesh and tell me what changed" and get a tool-calling assistant that drives the existing operations for them. Every other surface (Python API, CLI verbs, MCP tools) is imperative-only; this is the one conversational entry point missing.

- **Minimal REPL verb**: a `meshioplusplus chat` command (Python CLI only, `_cli/_chat.py`, alongside `data`/`dataset` as a third nested concern) that starts a terminal chat loop. Gated cleanly on an LLM API key being present (`ANTHROPIC_API_KEY` first-class; provider-agnostic wiring behind the same interface is a stretch goal, not a v1 requirement) — absent, it fails by name (`pip install meshioplusplus[chat]` / "set ANTHROPIC_API_KEY"), the `meshioplusplus-mcp` entry-point precedent, never a bare traceback. **S–M**
- **Reuse, not reimplementation**: the assistant's tool-calling loop drives `mcp/_tools.py`'s existing `TOOL_REGISTRY` directly (in-process, no stdio round trip needed for a same-process CLI) — never a second copy of the tool dispatch, sandbox (`_resolve`), or `_json_safe` result-sanitizing logic MCP already owns. A new tool added to `TOOL_REGISTRY` is automatically available to the chatbot for free, the same parity guarantee `test_every_operation_has_a_tool` already gives the MCP surface. **M** (mostly wiring an agentic loop over an existing, stable tool registry)
- **Scope for v1**: a stateless-per-turn loop (send message + conversation history + tool schemas from `TOOL_REGISTRY`, execute any tool calls the model requests, feed results back, repeat until a plain-text reply) operating on files under the CWD or an explicit `--root`, mirroring the MCP server's own path sandbox rather than inventing a second one. Multi-turn context lives only in the terminal session; no persistence, no server. **S**
- **New optional extra**: `chat = ["anthropic>=0.40,<1"]` (the SDK, not a hard dependency of `meshioplusplus` itself — the `mcp` extra's precedent), kept out of `[all]`. `_tools.py` stays untouched (no SDK import); only the new `_cli/_chat.py` / a `mcp/_chat.py` module imports it, so the default CI matrix and every other surface are unaffected. **S**
- **Docs**: a `doc/chat.md` page (modelled on `doc/mcp.md`) and a CLAUDE.md entry once shipped, per the "Keeping docs in sync" rule.

*Recommended posture: this is a thin client over work that already exists (the MCP tool registry) rather than new mesh functionality — the honest estimate for a usable v1 is small (S–M), with provider-agnostic support and a persisted chat history as documented follow-ups rather than v1 requirements.*

---

## 8. Binding-surface parity

**The gap.** Every operation is meant to reach every surface, and almost all do. The exceptions are recorded per page as "documented gaps"; gathering them here is what makes them visible as work rather than as footnotes.

- **`decimate_volume` on the flat surfaces** — it shipped in v10.6.0 on Python, C++, the C API, both CLIs, the pipeline and MCP, with the Fortran, Julia, R and WebAssembly bindings recorded as a follow-up that has not landed (no `decimate_volume` symbol exists under `bindings/fortran`, `bindings/julia`, `bindings/r` or `bindings/wasm`). Each is a transcription of the existing `decimate` binding with the tet-flavoured counters. **S** per surface
- **The documented flat-ABI gaps** — each deferred because the C ABI cannot carry the shape cheaply, and each closable with a small additive entry point: the `frozen` node mask on `smooth`, `decimate`, `decimate_volume` and `optimize_volume` (S); named `point_sets`/`cell_sets` in `diff`, `merge` and `split`, which today only the Python shim compares and remaps (M); the combined `data_manage` call (the three primitives exist) (S); per-cell-type counts in `mio_stats` (S); the SVG/TikZ colouring parameters, which the registry's `(path, mesh)` writer lambdas cannot carry (S); the `MedInfo`/`ExodusInfo`/`GmshInfo` side channels the registry drops (M); R's data setters always writing `Float64`, which is what stops `mio_split(by = "region")` on a tag built in R (S); a structured pipeline report over the C ABI instead of JSON text (S); and a native-CLI `data export` (Parquet is written through pyarrow, so this needs an Arrow writer or a documented "Python only") (S–M).

---

## 9. Format completeness

**The gap.** A format is "supported" when it reads and writes; several still lose a specific construct on one side, and the losses are recorded in tests so they cannot be forgotten.

- **Exodus writer** — the reader maps element blocks, node sets and side sets to regions, but the writer emits neither node sets nor side sets, so exodus is a region *source* and not a round-trip target (`tests/python/test_region_roundtrip.py`, `READ_ONLY_REGIONS`); and there is no multi-step Exodus writer, which is a stateful object of `XdmfTimeSeriesWriter`'s shape rather than a `(path, mesh)` writer. **M** each
- **Regions, Phase 2** — UNV groups and Ansys components (absorbing `UnvInfo`/`AnsysInfo`), OpenFOAM boundary patches (face groups, so side regions), XDMF sets and a VTU convention are the formats whose native group concept still does not map onto `Region` (`PHASE_2` in the same test file; see [regions](./regions.md#the-format-matrix)). **M**
- **A MED metadata reader** — `read_metadata` on a MED file is a full read that still reports one time step, so MED is invisible to the [sequence layer](./sequences.md)'s step probing even though its reader honours `time_step`; a `read_med_metadata` filling `mTimeValues` closes that for free. **S**
- **Gmsh in the C++ reader** — `$Periodic` still throws (the Python shim falls back, every other surface cannot), format 4.0 declines, and gmsh's 14-node pyramid is passed through in gmsh's edge-lexicographic order while every consumer in the core (`detail/cell_faces.hpp`, the CGNS table) treats `pyramid14` as `pyramid13` plus a base-centre node — the permutation the reader applies to `pyramid13` is missing for type 14, in both the C++ and the Python readers. **M**
- **Provenance assumption notes** — two sites are wired as the reference pattern (`warn_regions_dropped`, the OFF cell-type skip); the sweep over the remaining `warn()` calls that are genuine conversion assumptions, each checked rather than transcribed, is the open half of [provenance](./provenance.md). **M**

---

## 10. Follow-ups recorded in other pages

Items that a page already owns in detail are listed here as pointers only; the host page states the constraints and the reasoning, and this table exists so the roadmap stays a complete list.

| Host page | Item | Effort |
| --- | --- | --- |
| [pipeline](./pipeline.md#follow-ups-recorded-not-implemented) | multi-mesh steps (`Inputs:` and `Output.Pattern`), a structured report accessor on the C ABI, Conan/vcpkg shipping the JSON parser | S–M |
| [GPU handoff](./gpu.md#phase-2-reading-directly-into-pinned-memory) | the `pinned_reads()` context manager and `to_cupy(pinned=True)` fast path over the C++ buffer-allocator hook that already ships | S |
| [interoperability](./interop.md#phase-2-open3d-and-dolfinx) | Open3D and DOLFINx targets, whose constraints (a copying API and a 400 MB wheel; a single cell type, MPI and a basix permutation) are recorded there | M |
| [Julia](./julia.md) | registration in the General registry and a JLL for the binary | S |
| [MCP server](./mcp.md) | the SDK is pinned below 2.0 because 2.0 removed `mcp.server.fastmcp`; porting `_server.py` to the 2.x API is a separate piece of work | S |
| [WebAssembly](./wasm.md) | `Parallel` in a sequence document is accepted and ignored with a warning there, since a wasm module has no process pool | S |

---

## Recently closed

Shipped items leave the sections above and land here, so the list stays an accurate "not built" list without losing the history of what closed and when.

| Version | What closed | Where |
| --- | --- | --- |
| v10.21.0 | The Blender add-on: `to_blender`/`from_blender`, the 4.2+ extension, per-platform zips on every release | [blender](./blender.md) |
| v10.18.0 – v10.20.0 | GiD postprocess: all three write flavours, the hand-rolled reader, Gauss points, `ResultGroup`, multi-step series, and `gid` in every release binary and wheel via a vendored static zlib | [gid](./formats/gid.md) |
| v10.15.0 – v10.17.0 | Provenance in written files: one engine-identical credit line, the opt-in record (source, operations, assumptions, timestamp), read-back through `read_metadata`, and the Fortran/Julia/R/WASM bindings | [provenance](./provenance.md) |
| v10.10.0 – v10.14.0 | Remeshing: `remesh` (isotropic, quadric and anisotropic metrics, curvature gradation, boundaries), `remesh_volume` (isosurface stuffing), ODT smoothing and `optimize_volume` (flips) | [remesh](./remesh.md), [remesh_volume](./remesh_volume.md), [optimize_volume](./optimize_volume.md) |
| v10.8.0 – v10.9.0 | Field integration (`data_integrate`) and second derivatives (`hessian`) | [field integration](./field_integration.md), [hessian](./hessian.md) |
| v10.2.0 – v10.6.0 | Error estimation (`estimate_error`), polyhedral refinement (`subdivide`) and coarsening (`agglomerate`), green-element undo (`undo_green`), volume decimation (`decimate_volume`) | [error](./error.md), [subdivide](./subdivide.md), [agglomerate](./agglomerate.md), [undo_green](./undo_green.md), [decimate_volume](./decimate_volume.md) |
| v9.28.0 – v9.30.0 | PhysicsNeMo: the adapter, dataset manifests, the browser dataset manager, t→t+1 pairing, the `physicsnemo.mesh.Mesh` bridge | [physicsnemo](./physicsnemo.md), [datasets](./datasets.md) |
| v9.24.0 – v9.25.0 | `grid`, `voxelize`, signed distance (`compute_sdf`, the octree), the `.vti` format, `crop_predicate` | [voxelize](./voxelize.md), [sdf](./sdf.md) |

---

## Suggested sequencing

1. **Primitive constructors (§6)** — a few days, and it improves testing, docs and every demo surface at once.
2. **The large-mesh benchmark tier (§2)** — cheap, and it decides whether streaming and out-of-core work matter.
3. **Fuzzing (§4)** — start in parallel with everything else; it is not a feature and does not compete for the same attention.
4. **Binding parity (§8)** — small, mechanical, any time; the `decimate_volume` bindings first, since they are the only operation missing from a surface outright.
5. **Format follow-ups (§9)** — the MED metadata reader and the gmsh `pyramid14` permutation are days each; the Exodus writer and Phase-2 regions are the real work.
6. **NURBS spike (§5)** — a documented investigation, scheduled independently of the rest; no implementation until it reports.
7. **Dataset dashboard (§1)** — the non-training pieces can proceed any time; training launch and monitoring wait on their own design pass for the server-side companion process.
8. **CLI chatbot (§7)** — small and self-contained; can proceed whenever a maintainer wants it.
9. **Rust bindings and the registrations (§3)** — logistics, scheduled around releases rather than features.
