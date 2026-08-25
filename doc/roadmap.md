# meshio++ roadmap

Status at time of writing: **v10.18.0** — 43 formats (`gid` write-only), thirty-four mesh operations + five data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers plus a browser dataset manager, an MCP server, a settings-driven pipeline engine, a dataset-manifest layer with a PhysicsNeMo adapter, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 11).

This document lists what is *not* built. Items are grouped by theme, each with an effort estimate and the reason it matters. Nothing here duplicates shipped functionality; where a feature partially exists, the gap is stated explicitly.

Effort key: **S** = days, **M** = a couple of weeks, **L** = a month or more, **XL** = a project in its own right.

---


## 1. GiD postprocess format: the reader (write path shipped in v10.18.0)

**What shipped in v10.18.0.** GiD is a widely used pre/postprocessor in the same FE community meshio++ serves; its postprocess format is now written (not yet read) through a vendored hardcopy of CIMNE's [gidpost 2.14](https://www.gidsimulation.com/downloads/gidpost-2-14-library-to-write-postprocess-results-for-gid-in-ascii-binary-or-hdf5-format/) (`src/cpp/third_party/gidpost/`, BSD-2-Clause-Views), a hardcopy rather than a submodule (the `pugixml` precedent). All three on-disk flavours (ascii `.post.msh`/`.post.res` pair, binary `.post.bin`, hdf5 `.post.h5`) are covered, hard-gated on zlib (`gidpostInt.h` includes `<zlib.h>` unconditionally); reachable from every registry-driven surface (WASM, C API, Fortran, native CLI) with no per-binding code. Node ordering (`vertex`/`line`/`line3`/`triangle`/`triangle6`/`quad`/`quad8`/`quad9`/`tetra`/`tetra10`/`hexahedron`/`hexahedron20`/`wedge`/`pyramid`, all identity) was independently cross-checked against Kratos Multiphysics's own production GiD writer and pinned by a written-bytes ordering oracle (`tests/cpp/test_gid.cpp`'s `GidOrdering` suite), never a round trip — gidpost has no reader to round-trip through. See `doc/formats/gid.md` and `CHANGELOG.md`'s v10.18.0 entry for the full write-side design (material ids, the Gauss-point-set idiom for `cell_data`, the result component-count rule, provenance).

**The gap that remains.** gidpost's own public API has **zero read functions**, so the reader is ordinary meshio++ format-module work against the GiD postprocess-data-files specification, entirely independent of the vendored library.

- **ASCII reader** — `.post.msh`/`.post.res` is human-readable text; a hand-rolled parser (mesh sections, coordinate/element blocks, Gauss-point-set declarations, result blocks) is the natural first cut, mirroring the shape of any other hand-rolled text-format reader in this repo. **M**
- **Compressed-binary and HDF5 readers** — separate, later steps with their own cost (binary needs the inverse of gidpost's own record layout; HDF5 needs walking the group structure the writer already produces). Can land after the ASCII reader closes the format's most-requested direction (round-tripping a GiD case). **M–L**
- **Node orderings still unverified** — `hexahedron27`, `wedge15`, `pyramid13` throw `WriteError` by name on write (v10.18.0's deliberate scope cut); a reader for these needs the same independent verification the shipped types received before they can be added to either direction.
- **Tensor-aware results** — `GiD_Matrix`/`MainMatrix`/`ComplexVector` result types have no meshio++ counterpart today (a 6-component array splits into six scalars on write rather than risking a wrong component-order guess); an opt-in `GidInfo` side channel (the `MedInfo`/`GmshInfo` shape) carrying the declared result type/component order is the natural home, needed by both directions once attempted.
- **Regions → GiD mesh groups** — `GiD_fBeginMeshGroup`/`EndMeshGroup` is the natural target for named regions, not yet attempted in either direction; interacts with results in ways not yet verified.
- **Multi-step results** — the writer emits exactly one step per call today; a reader needs to handle a file with several, and a multi-step writer is a separate, stateful-object-shaped addition (the `XdmfTimeSeriesWriter` precedent).
- **`GiD_PostAsciiZipped`** — a fourth gidpost flavour with no extension of its own, not exposed by `GidMode` today.
- **`gid` in the release CLI binaries** — the statically-linked release CLI and the Windows wheels build `-DMESHIOPLUSPLUS_WITH_ZLIB=OFF` and so ship without `gid` today; flipping that default needs its own cross-platform static-link validation, deliberately deferred rather than bundled into the write-path release.

*Recommended posture: the ASCII reader is the natural next slice — self-contained, and it closes the format's most-requested direction (a GiD case round-tripping through meshio++) without needing the binary/HDF5 read flavours at all.*

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

**The gap.** Every operation transforms a mesh you already have; nothing creates one. This is the only empty category in the operations layer — **partly closed as a side effect of `remesh_volume`'s own delivery** (v10.13.0), which accepts a closed surface directly and generates a genuinely new tetrahedral volume mesh (isosurface stuffing over a BCC lattice) rather than transforming the input's own cells; a surface-in/volume-out capability this section previously lacked entirely. It does not close any bullet below, though: no primitive constructors, no extrude/revolve, and no 2D output (it is a 3D volume mesher, not a 2D triangulator) — those gaps stand as stated.

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

1. **Remeshing — shipped in full**, removed from this document entirely: v10.10.0 (`remesh`, a clean-room isotropic clustering engine derived from the MIT `pyvista/pyacvd` plus an original `metric="quadric"` synthesis over meshio++'s own pre-existing Garland-Heckbert quadric machinery, across every binding surface in one release); v10.11.0 (curvature gradation via a new osculating-paraboloid vertex-curvature estimator, plus boundary pinning/dual-edge insertion and separately-reported non-manifold-vertex detection, both extending `remesh` and its full binding surface); v10.12.0 (`metric="anisotropic"`, packing a curvature tensor into the existing quadric accumulator, closing the section's L estimate faster than expected once gradation had already supplied the curvature machinery it needed — `MESHIOPLUSPLUS_ABI_VERSION` 8→9); v10.13.0 (`remesh_volume`, isosurface stuffing over a BCC lattice, closing the volumetric bullet's *capability* by a predicate-free algorithm rather than the literal Delaunay/CVD method named, plus `SmoothMethod::Odt` closing its "ODT" half's *smoothing* piece on fixed connectivity — `MESHIOPLUSPLUS_ABI_VERSION` 9→10, not because either operation's own new option/result struct needed it, but because `SmoothMethod` gaining an explicit `: std::uint8_t` underlying type for the first time is itself a Tier A layout change under `doc/abi.md`'s own rule); v10.14.0 (`optimize_volume`, closing the "ODT" bullet's *remeshing* piece — predicate-free 2-3/3-2 topological flips alternated with the ODT vertex relocation, genuinely changing connectivity, Tier C additive with no ABI bump, full binding surface in one release). See `doc/remesh.md`, `doc/remesh_volume.md` and `doc/optimize_volume.md`.
2. **Provenance in written files** — shipped in full and removed from this document: v10.15.0 (the audit-and-normalize bullet — one canonical, engine-identical credit line, plus the per-format comment-syntax table in `doc/formats.md#provenance`), v10.16.0 (the opt-in record — source/target, operation chain, conversion assumptions, timestamp, over a thread-local scope with a Python↔C++ bridge, and the `mio_provenance_*` C ABI), v10.17.0 (read-back through `read_metadata`/`info`, plus the Fortran/Julia/R/WASM bindings — `MESHIOPLUSPLUS_ABI_VERSION` 10→11, since `MeshMetadata` grew two fields). Widening conversion-assumption coverage past the wired formats is a documented follow-up in `doc/provenance.md`, not a roadmap gap. See `doc/provenance.md`.
3. **GiD postprocess write path — shipped in v10.18.0**: vendoring, build integration, all three write flavours (ascii/binary/hdf5), extension dispatch, and reaching every registry-driven surface. The hand-rolled reader (§1) remains open and is self-contained, blocking nothing else in this document.
4. **Primitive constructors (§7, first item)** — a few days, and it improves testing, docs and every demo surface at once. `grid` already shipped over `detail/grid_lattice.hpp`; `box`/`sphere`/`cylinder`/`disk` follow the same shape.
5. **PhysicsNeMo integration** — shipped in full and removed from this document: v9.28.0 (recon note, adapter, dataset manager, recipes, GPU-executed example), v9.29.0 (dataset-manager UI), v9.30.0 (t→t+1 target pairing, the `physicsnemo.mesh.Mesh` bridge, persisted directory handles, per-entry quality summaries). See `doc/physicsnemo.md` and `doc/datasets.md`.
6. **Remaining refinement and coarsening gaps** — shipped in full and removed from this document: v10.2.0 (error-estimator helpers — `estimate_error`), v10.3.0 (polyhedral refinement — `subdivide`), v10.4.0 (polyhedral coarsening — `agglomerate`), v10.5.0 (green-element undo — `undo_green`), v10.6.0 (volume decimation — `decimate_volume`, the section's last open item). See `doc/error.md`, `doc/subdivide.md`, `doc/agglomerate.md`, `doc/undo_green.md` and `doc/decimate_volume.md`.
7. **Field capability beyond derivatives** — shipped in full and removed from this document: v10.8.0 (field integration — `data_integrate`, a cell-measure-weighted total/mean over cells, whole-mesh and per named Cell region), v10.9.0 (second derivatives — `hessian`, a composition of two `gradient` calls, exact for a linear field everywhere and for a quadratic field away from a structured mesh's own boundary, closing the section's last open item). See `doc/field_integration.md` and `doc/hessian.md`.
8. **Fuzzing (§5)** — should start in parallel with all of the above; it is not a feature and does not compete for the same attention.
9. **NURBS spike (§6)** — a documented investigation, scheduled independently of the rest.
10. **Dataset dashboard (§2)** — the non-training pieces (multi-dataset overview, drill-down, health summaries, diffing) can proceed any time; training launch/monitoring waits on its own design pass (server-side companion process) before implementation.
11. **CLI chatbot (§8)** — small and self-contained (a thin client over the existing MCP tool registry); can proceed independently whenever a maintainer wants it, no sequencing dependency on anything above.
