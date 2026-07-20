# Selective reads and `read_metadata`

Reading a mesh file is normally all-or-nothing: every data array is decoded whether or not
you want it. Two options narrow that, and a third call skips the arrays entirely.

Both are **additive**. Default `read()` behaviour is unchanged.

## Reading part of a file

```python
import meshioplusplus

mesh = meshioplusplus.read("big.vtu")                       # everything (default)
mesh = meshioplusplus.read("big.vtu", points_only=True)     # geometry, no data arrays
mesh = meshioplusplus.read("big.vtu", arrays=["u", "p"])    # only these data arrays
```

- `points_only=True` keeps points **and connectivity** — it narrows data, not topology, so
  the result is still a usable mesh you can convert or render.
- `arrays=None` (the default) reads every array; `arrays=[]` reads **none**. That distinction
  is deliberate and is preserved down to the C ABI.
- Names not present in the file are ignored rather than raising, so the same call works
  across a directory of meshes that do not all carry the same fields.

## Summarizing without loading

```python
meta = meshioplusplus.read_metadata("big.vtu")
meta["num_points"]              # 1_048_576
meta["cell_blocks"]             # [{"type": "tetra", "num_cells": 5_000_000, ...}]
meta["point_data_names"]        # ["p", "u"]
meta["fell_back_to_full_read"]  # False -> the summary really was cheap
```

`bbox_min`/`bbox_max` are **absent** rather than `None` when no bounding box was computed, so
"not computed" cannot be mistaken for a real box at the origin.

## Which formats are actually faster

| Format | Selective reads | `read_metadata` |
| --- | --- | --- |
| XDMF | native | native, and genuinely O(1) |
| VTU / VTP | native | native |
| Gmsh 4.1 | native | native |
| Gmsh 2.2 | native | falls back to a full read |
| everything else | read whole, then filtered | falls back to a full read |

**A fallback is correct, just not fast**, and it always says so via
`fell_back_to_full_read`. A partial read that silently wasn't partial would be worse than no
feature at all, so the flag is exposed on every binding surface.

### How much faster, honestly

- **XDMF** is the best case: every `<DataItem>` declares its shape in a `Dimensions`
  attribute, so counts are exact without touching any payload — and on the HDF path without
  opening the sibling `.h5` at all.
- **VTU/VTP**: the file is still read and XML-parsed, because pugixml always materializes
  PCDATA. What is skipped is base64 decoding (which is sequential, and often the larger half),
  decompression, allocation and byte-swapping. Expect a solid multiple, **not** an asymptotic
  change. Truly O(1) VTU metadata would need the *appended* data format with `offset=`
  attributes, which this reader does not accept — recorded as future work.
- **Gmsh 4.1** groups elements into typed blocks whose headers carry the type and count, so the
  summary walks block headers and skips each payload — by exact byte arithmetic for binary, line
  counts for ascii. Measured ~4× (ascii) and ~8× (binary) against a full read. Unwanted
  `$NodeData`/`$ElementData` bodies are skipped wholesale.
- **Gmsh 2.2** stores a type on *every element*, so there is no cheap summary to be had; it
  declines and falls back rather than pretending. That is the honest answer for the format, not
  an omission.

## Command line

```bash
meshioplusplus info --fast big.vtu               # summarize from the header
meshioplusplus convert --points-only in.vtu out.vtu
meshioplusplus convert --arrays u,p in.vtu out.vtu
```

`info --fast` prints `(no header-only path for this format; the file was read in full)` when
the summary was not actually cheap. `--points-only`/`--arrays` cannot be combined with
`-s`/`-d`, which operate on exactly the arrays that were skipped.

Both CLIs — the Python one and the native binary — behave identically.

## Other languages

```c
mio_read_opts opts;
mio_read_opts_init(&opts);          /* always: fields added later default sensibly */
opts.points_only = 1;
mio_mesh* mesh = mio_read_ex("big.vtu", "vtu", &opts);

mio_read_metadata* meta = mio_read_metadata_create("big.vtu", NULL);
int64_t n = mio_read_metadata_num_points(meta);
mio_read_metadata_free(meta);
```

`mio_read` is unchanged. `mio_read_opts` is permanent ABI: it carries reserved capacity and
may only grow additively.

```fortran
call m%read('big.msh', points_only=.true.)
meta = mio_read_metadata('big.msh')
```

```javascript
const mesh = m.readMeshSelective('big.vtu', { arrays: ['u'] });
const meta = m.readMetadata('big.vtu');
```
