# meshio++ roadmap

Status at time of writing: **v10.16.0** — 42 formats, thirty-four mesh operations + five data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers plus a browser dataset manager, an MCP server, a settings-driven pipeline engine, a dataset-manifest layer with a PhysicsNeMo adapter, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 10).

This document lists what is *not* built. Items are grouped by theme, each with an effort estimate and the reason it matters. Nothing here duplicates shipped functionality; where a feature partially exists, the gap is stated explicitly.

Effort key: **S** = days, **M** = a couple of weeks, **L** = a month or more, **XL** = a project in its own right.

---

## 1. Provenance and traceability in written files

**The gap.** A meshio++ output file used to say almost nothing about where it came from. **v10.15.0 closed the audit-and-normalize bullet**: every writer with a free-text header slot emits one canonical, character-identical, deterministic line (`"Written by meshio++ v<release>"`) from a single source on each side — see [`doc/formats.md`](formats.md#provenance) for the per-format audit table this produced. **v10.16.0 closed the substantive half** (the remaining bullets below, in full): a caller can now opt into a richer block recording the source path/format, the target format/encoding/codec/float-format actually used, the operation chain, and the conversion assumptions accepted on the way (`Side` regions dropped by name, wedge node orderings permuted, ragged blocks tessellated, dtypes promoted) — see [`doc/provenance.md`](provenance.md) for the design note and the mechanism. **What remains is bullet 8 alone**: reading a file's own embedded block back into a report (`info`/`read_metadata`), so the feature is legible on read as well as write.

- **Surfacing it on read** — report a file's provenance block back through `info`/`read_metadata`, so the feature is legible rather than write-only. Round-trip safety (closed in v10.16.0, see `doc/provenance.md#round-trip-safety`) already guarantees no reader mistakes the block for data; this bullet is purely additive: parse the block a supporting reader finds, expose it in the metadata/info report, never re-emit it. **S**
- **Extending conversion-assumption coverage** — v10.16.0 wired the mechanism and two reference call sites (`detail::warn_regions_dropped`'s 9-operation choke point, and the OFF writer's cell-type-skip warning in both engines); the remaining ~60-80 `log::warn`/`warn()` sites across the format writers each need one more `provenance_note(category, detail)` call beside the existing warning. Mechanical, not a design question — the pattern and the cross-engine parity test (`test_provenance.py::test_off_writer_records_dropped_cell_types`) are the template. **M**
- **Fortran/Julia/R/WASM bindings for the scope API** — the C ABI (`mio_provenance_scope_begin`/`_end`, `mio_provenance_note`, `mio_provenance_set_source`, `mio_provenance_set_target`) exists and is what these bindings would ride, following the `decimate_volume` staged-rollout precedent; none has been written yet. **S**

*Recommended posture: bullet 8 (read-back) is the one design-worthy piece left — it wants to decide whether the parsed block lives in `MeshMetadata` directly or a side-channel struct, following the `MedInfo`/`GmshInfo` precedent, before implementation starts. The conversion-assumption coverage extension and the additional-language bindings are both mechanical follow-ups with an established pattern to copy, not new design.*

---

## 2. GiD postprocess format, via a vendored gidpost

**The gap.** GiD is a widely used pre/postprocessor in the same FE community meshio++ serves, and its postprocess format is absent from the 42. CIMNE ships [gidpost 2.14](https://www.gidsimulation.com/downloads/gidpost-2-14-library-to-write-postprocess-results-for-gid-in-ascii-binary-or-hdf5-format/), a small C library that writes GiD postprocess files in ASCII, compressed binary and HDF5, under a **BSD-2-Clause-Views** licence (CIMNE, 2015–2024) — permissive, and so vendorable directly into this MIT repo provided the notice travels with it. The decision recorded here is to **hardcopy** it (not a submodule) and build a full round-trip `gid` format on top. One fact shapes the whole item: **gidpost is write-only** — its public API contains zero read functions, because it exists to emit postprocess output — so the library covers the write path and the reader has to be meshio++'s own code.

- **Vendor the library** — a hardcopy under `src/cpp/third_party/gidpost/` (the bundled-`pugixml` precedent; `eigen`/`json`/`polyscope` are submodules and are the wrong model, since the point here is that the source travels with the repo). Copy `source/`'s C files and headers plus `LICENSE`; **exclude** `binaries/` (96 MB of prebuilt libraries), `doc/gidpost.pdf` (7 MB, reference material rather than a vendored asset), `examples/`, `gidpost-swig/`, `win/`, `fortran_module/`, and — this one is a licence obligation, not a preference — `cfortran/` together with `source/gidpostforAPI.c`, since the gidpost licence explicitly carves out `cfortran.h` as third-party code under its own distribution policy. meshio++ has its own Fortran binding over its own C API and needs none of that layer. The vendored surface is roughly 370 KB. Licence obligations: the `LICENSE` copied alongside, a `CITATION.cff` entry in the shape of the existing pyacvd and PeriLab notes, and a line in the docs. It is C rather than C++, so it needs its own compile settings and must not leak warnings into the core build. **S**
- **Build integration** — `MESHIOPLUSPLUS_WITH_GIDPOST`, default ON behind an `EXISTS` probe (the `MESHIOPLUSPLUS_WITH_JSON` shape) → `MESHIOPLUSPLUS_HAS_GIDPOST`. **Hard-gated on zlib**: `gidpostInt.h` includes `zlib.h` unconditionally and binary output is always compressed, so this is not an optional flavour but a build prerequisite; the HDF5 flavour rides `MESHIOPLUSPLUS_WITH_HDF5` on top. Compiled out, the entry point still exists and throws naming the flag (`partition_kahip_parts`' contract), never a silent downgrade to another format. Four surfaces have to keep building and each needs checking rather than assuming: Emscripten/WASM (does the vendored C compile against the zlib port?), the amalgamation (a vendored C library inside a single C++ header is a real question — the honest answer may be that `gid` is simply absent there, as HDF5-dependent formats already are), conan/vcpkg, and the statically-linked release CLI binaries. **S**
- **Write path** — map `Mesh` onto gidpost's API in all three flavours, emitting the `.post.msh` + `.post.res` pair. GiD element types and their node orderings must be derived from the specification and then pinned by a test that inspects the **written bytes**, not by a round trip through our own reader: the CGNS work established why, since a self-inverse permutation makes `read(write(m)) == m` hold even when the ordering table is wrong. **M**
- **Extension dispatch is a real obstacle, not a detail** — `registry.cpp` already maps `.post` to `permas` and `.msh` to `gmsh`, while GiD wants the *double* extensions `.post.msh` / `.post.res` / `.post.bin` / `.post.h5` and a two-file pair. `resolve_format`'s single-extension map cannot express either today. The multi-file precedents to follow are `ensight` (`.case`/`.geo`) and `tetgen`/`triangle` (`.node`/`.ele`), along with `_helpers.py`'s multi-file buffer guards, which will need extending. **S**
- **Read path, hand-rolled** — gidpost supplies nothing here, so the ASCII `.post.msh`/`.post.res` reader is ordinary meshio++ format-module work against the GiD postprocess-data-files specification. Compressed-binary and HDF5 reading are separate later steps with their own cost, and the format can ship write-only in between (the `svg`/`tikz`/`gmsh22` precedent for a writer with no matching reader key). **M–L**
- **Scope boundaries, stated rather than discovered** — GiD results are natively defined on Gauss points, and meshio++'s `cell_data` is `(n,)`/`(n,k)`, never per-node-within-cell 3-D: this is the same structural limit MED's ELNO/ELGA already documents, a boundary rather than a defect to fix. Result range tables, mesh groups and multi-step results need a decision too, most likely a side-channel struct in the `MedInfo`/`GmshInfo`/`OpenFoamInfo` shape. **S**
- **Reaching every surface** — registering the format in `registry.cpp` is what makes it reachable from WASM, the C API, Fortran, Julia, R and the native CLI at once; the WASM leg specifically depends on the Emscripten question above. **S**
- **Docs and tests** — a `doc/formats/gid.md` page, the `doc/formats.md` table row and the README list, round-trip tests in the `helpers.write_read()` shape, and — where a maintainer has GiD itself — an env-gated external-validation test that skips with an actionable reason, following the `cgnscheck` and `checkMesh` precedent. **S**

*Recommended posture: vendoring plus the ASCII write path is a well-bounded first release and can land on its own. The hand-rolled reader and the binary/HDF5 read flavours are where the real cost sits and should be sequenced after it, not bundled into the same change.*

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

**The gap.** Every operation transforms a mesh you already have; nothing creates one. This is the only empty category in the operations layer — **partly closed as a side effect of `remesh_volume`'s own delivery** (v10.13.0), which accepts a closed surface directly and generates a genuinely new tetrahedral volume mesh (isosurface stuffing over a BCC lattice) rather than transforming the input's own cells; a surface-in/volume-out capability this section previously lacked entirely. It does not close any bullet below, though: no primitive constructors, no extrude/revolve, and no 2D output (it is a 3D volume mesher, not a 2D triangulator) — those gaps stand as stated.

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

1. **Remeshing — shipped in full**, removed from this document entirely: v10.10.0 (`remesh`, a clean-room isotropic clustering engine derived from the MIT `pyvista/pyacvd` plus an original `metric="quadric"` synthesis over meshio++'s own pre-existing Garland-Heckbert quadric machinery, across every binding surface in one release); v10.11.0 (curvature gradation via a new osculating-paraboloid vertex-curvature estimator, plus boundary pinning/dual-edge insertion and separately-reported non-manifold-vertex detection, both extending `remesh` and its full binding surface); v10.12.0 (`metric="anisotropic"`, packing a curvature tensor into the existing quadric accumulator, closing the section's L estimate faster than expected once gradation had already supplied the curvature machinery it needed — `MESHIOPLUSPLUS_ABI_VERSION` 8→9); v10.13.0 (`remesh_volume`, isosurface stuffing over a BCC lattice, closing the volumetric bullet's *capability* by a predicate-free algorithm rather than the literal Delaunay/CVD method named, plus `SmoothMethod::Odt` closing its "ODT" half's *smoothing* piece on fixed connectivity — `MESHIOPLUSPLUS_ABI_VERSION` 9→10, not because either operation's own new option/result struct needed it, but because `SmoothMethod` gaining an explicit `: std::uint8_t` underlying type for the first time is itself a Tier A layout change under `doc/abi.md`'s own rule); v10.14.0 (`optimize_volume`, closing the "ODT" bullet's *remeshing* piece — predicate-free 2-3/3-2 topological flips alternated with the ODT vertex relocation, genuinely changing connectivity, Tier C additive with no ABI bump, full binding surface in one release). See `doc/remesh.md`, `doc/remesh_volume.md` and `doc/optimize_volume.md`.
2. **Provenance in written files (§1)** — shipped in full except read-back: the audit-and-normalize bullet in v10.15.0 (see `doc/formats.md#provenance`), the opt-in richer record (source/target, operation chain, conversion assumptions, timestamp, the thread-local scope mechanism and its Python<->C++ bridge) in v10.16.0 (see `doc/provenance.md`). Surfacing the block on read is the one bullet left, and wants its own short design note (where the parsed block lives in `MeshMetadata`) before implementation.
3. **GiD postprocess via vendored gidpost (§2)** — self-contained and blocking nothing: vendoring plus the ASCII write path is one bounded change, and the hand-rolled reader and binary/HDF5 read flavours follow as separate ones.
4. **Primitive constructors (§8, first item)** — a few days, and it improves testing, docs and every demo surface at once. `grid` already shipped over `detail/grid_lattice.hpp`; `box`/`sphere`/`cylinder`/`disk` follow the same shape.
5. **PhysicsNeMo integration** — shipped in full and removed from this document: v9.28.0 (recon note, adapter, dataset manager, recipes, GPU-executed example), v9.29.0 (dataset-manager UI), v9.30.0 (t→t+1 target pairing, the `physicsnemo.mesh.Mesh` bridge, persisted directory handles, per-entry quality summaries). See `doc/physicsnemo.md` and `doc/datasets.md`.
6. **Remaining refinement and coarsening gaps** — shipped in full and removed from this document: v10.2.0 (error-estimator helpers — `estimate_error`), v10.3.0 (polyhedral refinement — `subdivide`), v10.4.0 (polyhedral coarsening — `agglomerate`), v10.5.0 (green-element undo — `undo_green`), v10.6.0 (volume decimation — `decimate_volume`, the section's last open item). See `doc/error.md`, `doc/subdivide.md`, `doc/agglomerate.md`, `doc/undo_green.md` and `doc/decimate_volume.md`.
7. **Field capability beyond derivatives** — shipped in full and removed from this document: v10.8.0 (field integration — `data_integrate`, a cell-measure-weighted total/mean over cells, whole-mesh and per named Cell region), v10.9.0 (second derivatives — `hessian`, a composition of two `gradient` calls, exact for a linear field everywhere and for a quadratic field away from a structured mesh's own boundary, closing the section's last open item). See `doc/field_integration.md` and `doc/hessian.md`.
8. **Fuzzing (§6)** — should start in parallel with all of the above; it is not a feature and does not compete for the same attention.
9. **NURBS spike (§7)** — a documented investigation, scheduled independently of the rest.
10. **Dataset dashboard (§3)** — the non-training pieces (multi-dataset overview, drill-down, health summaries, diffing) can proceed any time; training launch/monitoring waits on its own design pass (server-side companion process) before implementation.
11. **CLI chatbot (§9)** — small and self-contained (a thin client over the existing MCP tool registry); can proceed independently whenever a maintainer wants it, no sequencing dependency on anything above.
