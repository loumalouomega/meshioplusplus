# Provenance

Record where a written file came from, what ran on it, and what conversion assumptions were accepted on the way — as an **opt-in** block a writer renders alongside the unconditional one-line credit every writer has emitted since v10.15.0.

```python
import meshioplusplus as mp
from meshioplusplus import _provenance

with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
    _provenance.set_source("bracket.msh", "gmsh")
    mesh = mp.read("bracket.msh")
    mesh = mp.clean(mesh, weld=True)
    _provenance.add_operation("clean(weld=true)")
    mp.write("bracket.vtu", mesh)
    print(s.get())  # everything recorded, whether or not the file could hold it
```

```
<!--
Written by meshio++ v10.16.0
Converted from bracket.msh (gmsh)
Operation: clean(weld=true)
Timestamp: 2026-08-24T15:00:00Z
-->
```

With no scope open, output is byte-for-byte what v10.15.0 wrote — the default has not moved. This page is the design note `doc/roadmap.md` section 1 asked for before implementing bullets 2 through 7; it also is the reference for the mechanism itself.

## The gap this closes, and the one it does not

Section 1's [audit-and-normalize bullet](./roadmap.md) (v10.15.0) fixed the drift in the one line every writer already emitted. It said nothing about what a file was converted *from*, what ran on it, or which of this library's many documented per-format compromises — `Side` regions with no equivalent, a permuted node order, a 2D point silently padded to 3D — were accepted while writing it. That is the gap this page closes: bullets 2 ("where the record lives"), 3 ("what it contains"), 4 ("the operation chain"), 5 ("determinism and the timestamp"), 6 ("round-trip safety") and 7 ("formats that admit no comment").

**Bullet 8 (surfacing provenance on read) is out of scope here.** The record this page describes is write-side only; reading a file's own embedded block back into a report — so `info`/`read_metadata` can show it — is real, useful, separate work, deferred to keep this change bounded.

## Two decisions the roadmap's own wording left in tension

The roadmap bullets, read together, ask for two things that cannot both hold.

**Bullet 3 asks the file to record "which surface actually wrote it (C++ core, numpy fallback, or a named binding)."** Bullet 5 requires the C++ core and its Python fallback to emit **character-identical** bytes — the same invariant that drove the whole v10.15.0 change, since a `(C++ core)`-vs-`v{version}` split in the old per-writer strings was exactly what made the fallback boundary visible in output bytes. An engine marker is that same mistake with a different name. This page resolves the tension in favour of byte-identity: which engine ran is a fact about the *process*, not the *mesh*, and it costs nothing to report through `current_record()`/`scope.get()` instead — see [Reading back what was recorded](#reading-back-what-was-recorded) below. Nothing in the file distinguishes the two engines.

**The record defaults to opt-in, not on.** The roadmap does not settle this explicitly, but v10.15.0's tag is unconditional and this page's richer block is not: opening no scope reproduces v10.15.0's bytes exactly, which is what keeps `tests/python/test_io_baseline.py`'s pinned hashes and every existing cross-engine parity claim (`test_tikz.py`, `test_svg.py`, `test_gmsh.py`, `test_mdpa.py`) untouched by this change rather than needing new exemptions. It also means round-trip safety (bullet 6, below) only has to hold for callers who asked for the block in the first place.

## Where the record lives

The roadmap names two candidate homes and reconnaissance ruled out the more obvious one.

`registry_write_ex` — `write_options.hpp`'s documented "single owner of write-this-format-with-these-parameters" — is reached by only **four** of the write paths in this library: the native CLI, `mio_write_ex`, the settings pipeline, and the sequence driver. Python's own path (`_helpers.write` → the per-format shim → `_core.<fmt>_write` → the format's free function) never touches it; neither does the WASM path (`registry_writers()`'s raw lambdas) nor plain `mio_write`. A record grown onto `WriteOptions` would therefore be invisible from the primary user surface — the one most callers actually use. It is also structurally the wrong shape: growing `WriteOptions` is a **Tier A** ABI change (it is pinned at 48 bytes in `tests/cpp/test_abi_layout.cpp`), and the C mirror `mio_write_opts` only preserves `sizeof` across an appended field where `sizeof(void*) == 8` — the exact hazard the roadmap bullet already flagged.

The record instead lives in a **thread-local, RAII-scoped context**, the shape `meshioplusplus::set_buffer_allocator` already established for a similar problem (`ndarray.hpp`), with one deliberate difference: **thread-local, not process-global**. Writes are not serialized against each other, so a global record would let one thread's provenance leak into a concurrent write on another. Every writer this page touches reads the active record exactly where it already read `kProvenanceTag` — inside its own body — so no function signature anywhere had to change to add this.

```cpp
// C++
meshioplusplus::detail::ProvenanceScope scope(
    meshioplusplus::detail::ProvenanceMode::BestEffort);
meshioplusplus::detail::provenance_set_source(path, format);
meshioplusplus::write_vtu(out_path, mesh, /*binary=*/true, /*zlib=*/true);
```

```python
# Python -- the exact twin
with _provenance.scope(_provenance.Mode.BEST_EFFORT):
    _provenance.set_source(path, fmt)
    mp.vtu.write(out_path, mesh, binary=True)
```

### The Python ↔ C++ bridge

Python's `_provenance.scope` does not only manage a Python-side stack — when the compiled extension is present, opening it **also** pushes a matching `ProvenanceScope` on the C++ side (`bindings/python/_core.cpp`'s `provenance_scope_push`/`_pop`), and every `note()`/`set_source()`/`set_target()`/`add_operation()` call mirrors into both. Without this, a scope opened from Python would be invisible to the ~40 of 44 formats whose write goes through the compiled C++ writer, and the feature would silently do nothing for the surface most callers use. `lines()` — what the pure-Python fallback writers call to render — reads back through the same bridge, so a Python writer running under a scope agrees with a C++ writer byte for byte: one rendering engine, not two independently-maintained ones. The bridge degrades to a pure-Python thread-local stack when the extension is absent, the same fallback posture every other part of this library takes.

## What the record contains

```
ProvenanceRecord {
    source_path, source_format          -- where the mesh came from, if a single call spans read+write
    target_format, encoding, codec, float_format   -- what was actually used, not what was requested
    operations: [str]                    -- the chain, in order, each rendered "op(key=value, ...)"
    notes: [(category, detail)]          -- conversion assumptions accepted; duplicates collapsed
    timestamp                            -- ISO-8601 UTC, or empty when suppressed
}
```

`source_path`/`source_format` and the four `target_*` fields are set explicitly (`provenance_set_source`/`provenance_set_target`) by whichever call spans both the read and the write — the pipeline and sequence drivers do this automatically; a bare `read()` then `write()` two calls apart has nothing to set them from unless the caller does so themselves. `operations` accumulates one entry per step a driver ran. `notes` is the conversion-assumptions list described next.

## Sourcing the operation chain

Reused, not reinvented: the settings pipeline already models a chain as `PipelineStep{op, params}` objects run in order (`operations/pipeline.hpp`), so `run_pipeline` opens nothing itself — a caller wraps the call in a scope — but every step it runs calls `provenance_add_operation` with the step rendered as `Op(Key=value, ...)`, keys sorted. `_pipeline.py`, documented elsewhere as a genuinely separate pure-Python engine (not a wrapper over the C++ one), renders the same shape from the same step dictionaries and is pinned to agree with the C++ engine byte for byte for the common bool/int/string parameter cases; a step carrying an unusual float value can render with a different textual form between `str()` and `ostream`'s default formatting, a known, minor, documented gap rather than a design goal.

**The honest limit is the in-memory API path.** A caller doing `mesh = mp.read(...)`, then `mesh = mp.refine(mesh, ...)`, then `mp.write(...)` as three separate calls has given nothing a way to know the three belong together — there is no mesh-carried history the data model holds, and inventing one would mean tracking every mutation of every `Mesh` regardless of whether anyone asked. The escape hatch is that the scope itself is public API: a caller who wants the chain recorded opens one and calls `add_operation` themselves around their own steps, exactly as the pipeline driver does internally.

## Conversion assumptions

The category of information the roadmap calls out as having "no substitute elsewhere" — node-order permutations, data dropped for want of a target concept, 2D point padding, dtype promotions, ragged-to-tessellated fallbacks — is captured by one call, `provenance_note(category, detail)` (`note(category, detail)` in Python), placed beside an *existing* `log::warn`/`warn()` call at the site the assumption is made. It is a **no-op outside a scope**, so it costs nothing to call unconditionally, and identical `(category, detail)` pairs collapse — a warning that would otherwise fire once per cell must not produce one record per cell.

Two call sites are wired as the reference pattern for the rest: `detail::warn_regions_dropped` (`detail/region_remap.hpp`), the single choke point nine operations already share for "the output cells have no correspondence with the input's," now also records a `"regions-dropped"` note; and the OFF writer's cell-type-skip warning, wired in **both** engines with matching wording so the rendered note is character-identical whichever one produced the file (`tests/python/test_provenance.py::test_off_writer_records_dropped_cell_types` pins this). Extending coverage to the ~60-80 remaining `log::warn`/`warn()` sites across the format writers is mechanical from here — the same one-line addition at each site — and is recorded as follow-up work rather than attempted exhaustively in this change; the mechanism, the cross-engine parity guarantee, and the two worked examples are what this page commits to.

## Determinism and the timestamp

`ProvenanceScope`'s constructor stamps `mTimestamp` once, at open, honouring `SOURCE_DATE_EPOCH` (the reproducible-builds convention — `provenance_timestamp()`/`timestamp()` read it directly) so a byte-comparable build can pin an exact instant, and `MESHIOPLUSPLUS_PROVENANCE_TIMESTAMP=off` suppresses it to an empty string (dropped from the rendered block entirely) when neither reproducibility nor a wall-clock record is wanted. `SOURCE_DATE_EPOCH` wins when both are set — a caller exporting both wants the pinned epoch, not silence. Host and user are never recorded; these files get shared.

Because the block is opt-in, this needed no exemption anywhere: `test_io_baseline.py`'s pinned hashes, and the C++/Python byte-parity suites for tikz/svg/gmsh/mdpa, all exercise the *default* (no scope) path, which is unconditionally what v10.15.0 wrote.

## Round-trip safety

No reader in this library retains any slot a provenance-bearing writer uses. The legacy VTK reader discards its title line outright (`vtk_read.cpp`); Exodus never reads its own `title` netCDF attribute back; OBJ/OFF/PLY/FLAC3D/MDPA all skip a comment line generically by its marker, never inspecting its content. "Replace, never append" is therefore structural rather than a rule a writer has to remember: every writer renders the block fresh from whatever the *live* record holds, never by reading and appending to what the input file already carried, so converting a file through N provenance-opted-in steps produces one block, not N.

## Slot tiers, and the `Mode`/`Required` interaction

`doc/formats.md`'s [Provenance table](./formats.md#provenance) records, per format, the comment syntax admitted and where a header may sit; `SlotTier` is that same classification made machine-readable, at the granularity a *writer's own call site* actually has (which can be narrower than what the format's grammar in general allows — Ansys admits a genuine `(0 "...")` multi-line comment section, for instance, but that writer stays at `SingleLine` here because populating it with arbitrary rendered text risks corrupting the reader's own paren-balance parsing on a note containing literal parentheses):

| Tier | Meaning | Formats |
|---|---|---|
| `Block` | Arbitrarily many lines | abaqus, ansys*, avsucd, exodus, flac3d, flux, mphtxt, nastran, netgen, obj, off, permas, ply, tetgen, triangle, vti, vtp, vtu |
| `SingleLine` | Exactly one line | ansys, tecplot, vtk/vtk42/vtk51 |
| `Bounded` | One line, hard byte cap | ensight, openfoam, stl (binary only) |
| `None` | No slot at all | every format `doc/formats.md`'s table marks "—" |

\* ansys's own comment grammar would support `Block`; the writer is pinned to `SingleLine` for the reason in the parenthetical above.

`Mode::Off` (the default) always renders just the tag, at every tier. `Mode::BestEffort` renders the full block at `Block` and degrades silently to the tag everywhere else — a `SingleLine`/`Bounded` slot cannot structurally hold more, and a `None` slot cannot hold even the tag (which is why no writer at that tier calls this machinery at all; nothing changed there since v10.15.0). `Mode::Required` is identical to `BestEffort` at every tier **except** `None`, where it raises `WriteError` naming the format. This is the resolution of `write_options.hpp`'s standing rule — "an option a format cannot honour is an error, never silently ignored" — that keeps `Required` from being unusable: degrading to the honest maximum a structurally smaller slot can hold is not a failure to report; having nothing to write at all is.

## Reading back what was recorded

A scope's `Get()`/`get()` returns the live record at any point while it is open — including the one fact this page deliberately keeps out of every file, which engine rendered it:

```python
with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
    mp.write("out.vtu", mesh)
    print(s.get())  # Record(source_path=..., operations=[...], notes=[...], ...)
```

## Language surfaces

| Surface | Entry points |
|---|---|
| C++ | `meshioplusplus::detail::ProvenanceScope`, `provenance_note`, `provenance_set_source`, `provenance_set_target`, `provenance_add_operation`, `current_provenance()` (`detail/provenance.hpp`) |
| Python | `meshioplusplus._provenance.scope`, `.note`, `.set_source`, `.set_target`, `.add_operation`, `.current_record()` |
| C API | `mio_provenance_scope_begin`/`_end`, `mio_provenance_note`, `mio_provenance_set_source`, `mio_provenance_set_target` |

The C API mirrors the C++/Python shape exactly, with the one adjustment every C ABI in this library makes for a stateful, multi-call object: `mio_provenance_scope_begin`/`_end` are two separate calls rather than RAII, backed by a thread-local stack on the binding side (the identical shape the Python bridge already uses) so a scope survives between them. Fortran, Julia, R and WASM bindings are a recorded follow-up, not yet shipped as of v10.16.0 — the same staged-rollout precedent `decimate_volume`'s bindings followed.

## What this page does not cover

- **Reading a file's own provenance block back on read** (roadmap bullet 8) — genuinely separate work, deferred.
- **Instrumenting the remaining ~60-80 conversion-assumption call sites** — mechanical, one `provenance_note` call per existing warning, not attempted exhaustively here (see [Conversion assumptions](#conversion-assumptions) above).
- **Fortran/Julia/R/WASM bindings for the scope API** — the C ABI exists and is what those bindings would ride; none has been written yet.
