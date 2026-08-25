# Provenance

Record where a written file came from, what ran on it, and what conversion assumptions were accepted on the way — as a block a writer renders alongside the one-line credit every writer has emitted since v10.15.0. **On by default** since v10.17.0: an ordinary `write()` records the assumptions raised while it ran, with no scope needed.

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
Written by meshio++ v10.17.0
Converted from bracket.msh (gmsh)
Operation: clean(weld=true)
Timestamp: 2026-08-24T15:00:00Z
-->
```

With no scope open a writer still renders any conversion assumptions raised while it ran. The richer fields — source, target, operation chain, timestamp — need a scope or a driver to set them, so a plain `write()` of a mesh that loses nothing is byte-for-byte what v10.15.0 wrote. This page is the design note `doc/roadmap.md` section 1 asked for, and the reference for the mechanism itself. Section 1 is now closed in full: v10.15.0 normalized the credit line, v10.16.0 added the record, v10.17.0 added read-back and the remaining language bindings.

## The gap this closes, and the one it does not

Section 1's [audit-and-normalize bullet](./roadmap.md) (v10.15.0) fixed the drift in the one line every writer already emitted. It said nothing about what a file was converted *from*, what ran on it, or which of this library's many documented per-format compromises — `Side` regions with no equivalent, a permuted node order, a 2D point silently padded to 3D — were accepted while writing it. That is the gap this page closes: bullets 2 ("where the record lives"), 3 ("what it contains"), 4 ("the operation chain"), 5 ("determinism and the timestamp"), 6 ("round-trip safety") and 7 ("formats that admit no comment").

**v10.17.0 closed the read side too** (bullet 8), so §1 is now complete and the section has been removed from the roadmap. See [Reading a block back](#reading-a-block-back) below.

## Two decisions the roadmap's own wording left in tension

The roadmap bullets, read together, ask for two things that cannot both hold.

**Bullet 3 asks the file to record "which surface actually wrote it (C++ core, numpy fallback, or a named binding)."** Bullet 5 requires the C++ core and its Python fallback to emit **character-identical** bytes — the same invariant that drove the whole v10.15.0 change, since a `(C++ core)`-vs-`v{version}` split in the old per-writer strings was exactly what made the fallback boundary visible in output bytes. An engine marker is that same mistake with a different name. This page resolves the tension in favour of byte-identity: which engine ran is a fact about the *process*, not the *mesh*, and it costs nothing to report through `current_record()`/`scope.get()` instead — see [Reading the live record](#reading-the-live-record-mid-write) below. Nothing in the file distinguishes the two engines.

**The record is on by default** (v10.17.0; it shipped opt-in in v10.16.0). What is on is the part that needs no caller cooperation: the conversion assumptions a writer raises while it runs. The fields that require someone to supply them — source, target, operation chain — stay absent without a scope, and **no timestamp is added**, so a write that loses nothing still reproduces v10.15.0's bytes exactly and every byte-pinned test (`test_io_baseline.py`, `test_tikz.py`, `test_svg.py`, `test_gmsh.py`, `test_mdpa.py`) is untouched. Turn it off with `MESHIOPLUSPLUS_PROVENANCE=off` or `set_default_provenance_mode(Off)`.

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

## How long a note lives

A note has no natural lifetime without a scope, and getting this wrong writes a *false* statement into a file. Before the bounding below existed, `extract_surface(A)` dropping a region left its note in place, and the next `write()` — of an unrelated mesh B — put `Note [regions-dropped]: extract_surface: ...` into B's header.

So the **public write entry points bound it**: `meshioplusplus.write()`, `registry_write_ex` and `mio_write` each call `provenance_begin_write()`, which resets the scope-less record so only notes raised by *this* write can be rendered. It is a no-op while a scope is open — the caller's scope owns the lifetime then, and spanning more than one write is the whole point of opening one.

Two consequences worth stating plainly:

- **A note raised before the write begins is dropped, not misattributed.** That covers the operation-side sites (`warn_regions_dropped` and friends). To capture those, open a scope spanning the operations *and* the write — which is exactly the example at the top of this page.
- **The low-level C++ writers do not bound notes for themselves.** Calling `meshioplusplus::write_vtu(...)` directly skips the public path, so notes from earlier work are still in scope-less storage. That is the same contract the rest of that layer has: it is the API where the caller manages state. Go through `registry_write_ex`, or open a scope.

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

## Reading the live record, mid-write

A scope's `Get()`/`get()` returns the live record at any point while it is open — including the one fact this page deliberately keeps out of every file, which engine rendered it:

```python
with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
    mp.write("out.vtu", mesh)
    print(s.get())  # Record(source_path=..., operations=[...], notes=[...], ...)
```

## Reading a block back

A file's own block comes back through the ordinary summary call — no second entry point, no second file open:

```python
meta = mp.read_metadata("bracket.vtu")
meta["provenance"]              # ['Written by meshio++ v10.17.0', 'Converted from ...', ...]
meta["provenance_recognised"]   # True when the first line is our own tag format
```

```
$ meshioplusplus info --fast bracket.vtu
  ...
  Provenance:
    Written by meshio++ v10.17.0
    Converted from bracket.msh (gmsh)
    Operation: Clean(Weld=true)
    Timestamp: 2026-08-24T15:00:00Z
```

### One scanner, not forty-four parsers

The block's *content* lines are format-independent by construction — only the comment punctuation wrapping them differs — so recovery is a single pass over the file's head bytes rather than a parser per format: strip whatever marker the format uses (`#`, `!`, `*`, `$`, `%`, `//`, PLY's `comment `, an XML `<!--`), then keep the lines matching the six shapes the renderer emits. `detail::read_provenance_lines` lives beside the renderer it inverts, which is what stops the two drifting.

Three wrapping shapes need more than a leading-marker strip, and each earns its place: a keyword slot that *wraps* the text (Tecplot's `TITLE = "…"`, Ansys's `(1 "…")`), where the closing punctuation is removed only because that opener put it there — never unconditionally, since `Converted from x (fmt)` and `Operation: Clean(Weld=true)` both end in a legitimate `)`; a box-drawn banner (OpenFOAM), whose credit sits in an interior `|`-delimited cell rather than at the line's start; and a NUL-padded fixed-width binary slot (STL's 80-byte header, EnSight's `str80`), handled by treating NUL as a line break so the same line-oriented scan finds it.

**Exodus is the one format the scanner does not serve**: its block rides a netCDF `title` attribute rather than the head bytes, so `read_exodus_metadata` reads it there — it already has the file open for `mTimeValues` — and hands it to the same scanner, so the two paths cannot disagree about what counts as a block.

### Raw lines, not a re-parsed record

Read-back deliberately returns the block's lines as found, plus a `recognised` flag, rather than re-parsing them into a `ProvenanceRecord`. A block can have been hand-edited, truncated by a `Bounded` slot, or written by a later release carrying fields this build has never heard of; handing back what is actually there — and saying separately whether it *starts* like ours — cannot silently drop or mis-attribute any of those. `recognised` is false both for a file with no block and for one carrying a comment that merely looks like a header, which is the honest distinction between "meshio++ wrote this" and "something left a comment here".

### Never re-emitted

Nothing carries a read block into a write: writers render from the live record only. That is what makes [round-trip safety](#round-trip-safety) structural rather than a rule each writer must remember — converting a file through N provenance-opted-in steps leaves one block, not N, and the first file's source never leaks into the second.

## Language surfaces

| Surface | Entry points |
|---|---|
| C++ | `meshioplusplus::detail::ProvenanceScope`, `provenance_note`, `provenance_set_source`, `provenance_set_target`, `provenance_add_operation`, `current_provenance()` (`detail/provenance.hpp`) |
| Python | `meshioplusplus._provenance.scope`, `.note`, `.set_source`, `.set_target`, `.add_operation`, `.current_record()` |
| C API | `mio_provenance_scope_begin`/`_end`, `mio_provenance_note`, `mio_provenance_set_source`, `mio_provenance_set_target`; read-back via `mio_read_metadata_num_provenance_lines`/`_provenance_line`/`_provenance_recognised` |
| Fortran | `mio_provenance_begin`/`_end`, `mio_provenance_note`, `mio_provenance_source`, `mio_provenance_target`, and the `MIO_PROVENANCE_*` mode constants |
| Julia | `provenance(f, mode)` (a `do`-block form; the scope closes even if `f` throws), `provenance_note`, `provenance_source`, `provenance_target` |
| R | `mio_with_provenance(expr, mode)` (paired with `on.exit`), `mio_provenance_note`, `mio_provenance_source`, `mio_provenance_target` |
| WASM | `withProvenance(mode, fn)` (a callback form; the scope closes even if `fn` throws), plus `provenanceNote`/`provenanceSetSource`/`provenanceSetTarget` and `readProvenance(path)` |

A scope cannot cross an ABI boundary as RAII, so every non-C++ surface has to give the begin/end pair a lifetime of its own. Each does it in its own idiom rather than exposing the raw pair: Julia a `do`-block, R an `on.exit`-paired wrapper, WASM a callback — all three close the scope even when the body throws, which is the whole reason they exist. Fortran keeps the explicit pair, matching this module's own explicit `m%free()` convention. The C ABI itself is backed by a thread-local stack of owned scopes on the binding side, the identical shape the Python bridge uses.

## What this page does not cover

- **Instrumenting every conversion-assumption call site.** The mechanism and the category vocabulary are settled, and the sites where *both* engines share a lossy path are wired (OFF, PLY, STL, UNV, CGNS) along with the C++-only ones (MDPA, MED, OpenFOAM). What is not done is a format-by-format sweep of the rest.

  Two findings from doing this much are worth carrying into that work. First, most `warn()` calls are **not** conversion assumptions — the majority are reader-side diagnostics or user-error messages, and sweeping them in indiscriminately would fill a block with things that are not assumptions at all. Second, some are simply **stale**: `avsucd`'s "can only write one cell data array" fires while *both* engines demonstrably write every array, so wiring it would have recorded a loss that does not happen. Each site needs checking, not translating.

  `tests/python/test_provenance.py::test_engines_record_the_same_notes` is the guard for that work: it asserts that wherever both engines write a given mesh, they record the same notes, so a one-sided addition fails at the format that caused it.
