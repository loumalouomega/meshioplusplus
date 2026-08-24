# Selective reads and `read_metadata`

Reading a mesh file is normally all-or-nothing: every data array is decoded whether or not you want it. Two options narrow that, a third picks *which time step* to decode, and a fourth call skips the arrays entirely.

All are **additive**. Default `read()` behaviour is unchanged.

![Identical geometry either way; the data array is loaded only on the left](/images/selective_read.png)

## Reading part of a file

```python
import meshioplusplus

mesh = meshioplusplus.read("big.vtu")                       # everything (default)
mesh = meshioplusplus.read("big.vtu", points_only=True)     # geometry, no data arrays
mesh = meshioplusplus.read("big.vtu", arrays=["u", "p"])    # only these data arrays
```

- `points_only=True` keeps points **and connectivity** — it narrows data, not topology, so the result is still a usable mesh you can convert or render.
- `arrays=None` (the default) reads every array; `arrays=[]` reads **none**. That distinction is deliberate and is preserved down to the C ABI.
- Names not present in the file are ignored rather than raising, so the same call works across a directory of meshes that do not all carry the same fields.

## Picking a time step

Formats carrying a time series hold every array once *per step*, but a `Mesh` holds one. Before `time_step` existed, every such reader silently took the first — correct for the overwhelmingly common single-step file, and quietly wrong for the rest.

```python
mesh = meshioplusplus.read("run.exo")                 # step 0, as always
mesh = meshioplusplus.read("run.exo", time_step=2)    # the third step
mesh = meshioplusplus.read("run.exo", time_step=-1)   # the last step

meshioplusplus.read_metadata("run.exo")["time_values"]   # [0.0, 0.5, 1.0]
```

- `0` (the default) is the first step, so existing behaviour is bit-for-bit unchanged.
- Negative values count from the end. `-1` means "the final state of the simulation", which is otherwise unsayable without knowing the step count up front.
- Out of range is an **error naming the available count**, never a silent clamp — handing back step 0 for a request of step 7 is exactly the failure this option removes.
- Unlike `points_only`/`arrays`, this is **not** a narrowing option: no caller-side filter can recover a step that was never read. So a format whose reader has no time concept **raises** rather than quietly returning the first step.
- `read_metadata(...)["time_values"]` reports the recorded times, so a request is checkable before it is issued. It is always present — empty for a format with no time concept — so `len(meta["time_values"])` needs no key test.

Currently honoured by **`exodus`**, **`xdmf`** (temporal collections — the counterpart to `XdmfTimeSeriesWriter`, see [XDMF time series](xdmf_time_series.md); the C++ reader resolves the collection structurally rather than running an XInclude/XPointer pass, and `read_metadata`'s `time_values` come off each step's `<Time Value>` attribute without touching a payload) and, since v9.9.0, **`med`** (a `CHA` field's `(NDT, NOR)` step subgroups, whose zero-padded group names sort into step order; `read_metadata`'s `time_values` are the steps' `PDT` attributes). A multi-step MED field used to fail the read outright unless there was a Python fallback to defer to. CGNS also has a time concept and is the natural next adopter; it still takes the first step today.

## Summarizing without loading

```python
meta = meshioplusplus.read_metadata("big.vtu")
meta["num_points"]              # 1_048_576
meta["cell_blocks"]             # [{"type": "tetra", "num_cells": 5_000_000, ...}]
meta["point_data_names"]        # ["p", "u"]
meta["fell_back_to_full_read"]  # False -> the summary really was cheap
```

`bbox_min`/`bbox_max` are **absent** rather than `None` when no bounding box was computed, so "not computed" cannot be mistaken for a real box at the origin.

## Which formats are actually faster

| Format | Selective reads | `read_metadata` | `time_step` |
| --- | --- | --- | --- |
| XDMF | native | native, and genuinely O(1) | ✗ (takes the first step) |
| VTU / VTP | native | native | n/a |
| Gmsh 4.1 | native | native | n/a |
| Gmsh 2.2 | native | falls back to a full read | n/a |
| Exodus | read whole, then filtered | falls back to a full read, but reports `time_values` | ✅ |
| MED | read whole, then filtered | falls back to a full read | ✅ |
| everything else | read whole, then filtered | falls back to a full read | ✗ |

`reader_supports_options(fmt)` reports whether a format has a native options-aware path at all. Note that "options-aware" is not one capability but several: Exodus is on that list for `time_step`, not because it narrows arrays natively, and MED is on it for `time_step` and `lenient` (see [MED](formats/med.md#lenient-reads)).

**A fallback is correct, just not fast**, and it always says so via `fell_back_to_full_read`. A partial read that silently wasn't partial would be worse than no feature at all, so the flag is exposed on every binding surface.

### How much faster, honestly

- **XDMF** is the best case: every `<DataItem>` declares its shape in a `Dimensions` attribute, so counts are exact without touching any payload — and on the HDF path without opening the sibling `.h5` at all.
- **VTU/VTP**: the file is still read and XML-parsed, because pugixml always materializes PCDATA. What is skipped is base64 decoding (which is sequential, and often the larger half), decompression, allocation and byte-swapping. Expect a solid multiple, **not** an asymptotic change. Truly O(1) VTU metadata would need the *appended* data format with `offset=` attributes, which this reader does not accept — recorded as future work.
- **Gmsh 4.1** groups elements into typed blocks whose headers carry the type and count, so the summary walks block headers and skips each payload — by exact byte arithmetic for binary, line counts for ascii. Measured ~4× (ascii) and ~8× (binary) against a full read. Unwanted `$NodeData`/`$ElementData` bodies are skipped wholesale.
- **Gmsh 2.2** stores a type on *every element*, so there is no cheap summary to be had; it declines and falls back rather than pretending. That is the honest answer for the format, not an omission.

## Command line

```bash
meshioplusplus info --fast big.vtu               # summarize from the header
meshioplusplus convert --points-only in.vtu out.vtu
meshioplusplus convert --arrays u,p in.vtu out.vtu
meshioplusplus convert --time-step=-1 run.exo last.vtu   # negatives need the = form
```

`info --fast` prints `(no header-only path for this format; the file was read in full)` when the summary was not actually cheap, and `Time steps: N [...]` when the file records more than one (a single-step file gives you nothing to choose). `--points-only`/`--arrays` cannot be combined with `-s`/`-d`, which operate on exactly the arrays that were skipped.

Both CLIs — the Python one and the native binary — behave identically.

## Other languages

```c
mio_read_opts opts;
mio_read_opts_init(&opts);          /* always: fields added later default sensibly */
opts.points_only = 1;
opts.time_step = -1;  /* the last step; 0 (the default) is the first */
mio_mesh* mesh = mio_read_ex("big.vtu", "vtu", &opts);

mio_read_metadata* meta = mio_read_metadata_create("big.vtu", NULL);
int64_t n = mio_read_metadata_num_points(meta);
mio_read_metadata_free(meta);
```

`mio_read` is unchanged. `mio_read_opts` is permanent ABI: it carries reserved capacity and may only grow additively.

```fortran
call m%read('big.msh', points_only=.true.)
call m%read('run.exo', time_step=-1)   ! the last step
meta = mio_read_metadata('big.msh')
```

```javascript
const mesh = m.readMeshSelective('big.vtu', { arrays: ['u'] });
const last = m.readMeshSelective('run.exo', { format: 'exodus', timeStep: -1 });
m.readMetadata('run.exo', 'exodus').timeValues;  // [0, 0.5, 1]
const meta = m.readMetadata('big.vtu');
```

## `lenient`: skipping what a reader cannot represent

A reader that meets a construct it has no way to express throws `ReadError` naming it, so a caller with a Python fallback can take it and a caller without one at least learns what was in the file. That is the right default — but it makes formats whose *optional* sections are commonplace unreadable in practice for a caller that only wants the mesh.

`ReadOptions::mLenient` (C `mio_read_opts.lenient`, Fortran/Julia/R/WASM `lenient`, `--lenient` on the native CLI) turns those specific rejections into a warning plus a skip.

It is deliberately **not** "ignore all errors". It applies only where a reader can skip a construct and still return a *correct* mesh; a malformed file, a truncated block, a bad node reference or an unknown element type still throw, because continuing past those would hand back a mesh that is quietly wrong rather than merely incomplete.

**Currently honoured by `mdpa` only** — its `Table`, `Geometries`, `Mesh`, `Constraints` and non-empty `SubModelPart*` blocks. Every other reader ignores the flag. `MdpaInfo::mSkippedConstructs` records what was skipped, so "lenient" never means "silently lossy". See [`doc/formats/mdpa.md`](formats/mdpa.md).

The Python `read()` deliberately does **not** take this parameter: mdpa's Python path is the pure-Python reference reader, which already accepts every construct the flag covers, so it would be a dead argument.
