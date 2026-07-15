# Modulef Formatted Mesh — MFM (`.mfm`)

The Modulef Formatted Mesh is a compact ASCII format used by
[FEconv](https://github.com/victorsndvg/FEconv) (a simplified NOPO/Modulef
mesh). A file holds a **single element type** (non-hybrid): an 8-integer
header, the connectivity, per-entity reference arrays, the vertex coordinates
and a per-element subdomain array.

| | |
|---|---|
| **Format name** | `mfm` |
| **Extensions** | `.mfm` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.mfm")
meshioplusplus.mfm.write("out.mfm", mesh, float_fmt=".16e")
```

- **`float_fmt`** (default `".16e"`) — coordinate format string.

## File structure

```
nel nnod nver dim lnn lnv lne lnf   # header, 8 ints
```

`nel`/`nver` are the element/vertex counts; `lnn`/`lnv`/`lne`/`lnf` are the
*local* (per-element) node/vertex/edge/face counts, which together identify
the element type. After the header, the file is a flat whitespace-token
stream (no line-structure requirement — Python reads it with `f.read().split()`
and the C++ reader does the equivalent) laid out as:

1. **`mm`** — connectivity, `nel × lnv` integers, 1-based, element-major.
2. **`nrc`** — face reference array, `nel × lnf` integers — **only present if
   `dim == 3`** — read but discarded.
3. **`nra`** — edge reference array, `nel × lne` integers — **only present if
   `dim >= 2`** — discarded.
4. **`nrv`** — vertex reference array, `nel × lnv` integers — always present —
   discarded.
5. **`z`** — vertex coordinates, `nver × dim` floats, vertex-major.
6. **`nsd`** — per-element subdomain/material integer, `nel` values →
   `cell_data["mfm:ref"]`.

The writer emits the same sections in the same order, with `nrc`/`nra`/`nrv`
written as all-zero placeholders (they carry no meshio++-side data).

## Cell types

The element type is recovered from `(lnv, lne, lnf)` plus `lnn == lnv`:

| type | lnv | lne | lnf |
|---|---|---|---|
| `line` | 2 | 1 | 0 |
| `triangle` | 3 | 3 | 1 |
| `quad` | 4 | 4 | 1 |
| `tetra` | 4 | 6 | 4 |
| `hexahedron` | 8 | 12 | 6 |
| `wedge` | 6 | 9 | 5 |

Linear elements only. Because MFM only stores **vertex** coordinates (no
mid-edge node positions), second-order (P2) elements would necessarily be
straight-sided if forced through this representation — meshio++ rejects them
outright (requires `lnn == lnv`) rather than silently discarding curvature
information.

## Data mapping

- `cell_data["mfm:ref"]` — per-element subdomain/material reference (a single
  int array; defaults to all-ones on write if absent).

## Quirks & limitations

- **Single element type per file**: the writer requires exactly one cell type
  across the whole mesh, raising `WriteError` for a mixed-type mesh — MFM is
  fundamentally single-type ("non-hybrid").
- `nrc`/`nra`/`nrv` (face/edge/vertex reference arrays) are read and discarded
  entirely — there is no meshio++-side representation for them, and they're
  always written as zeros.
- Requires `nnod == nver` (no separate "node" vs "vertex" numbering) in
  addition to `lnn == lnv` — both checks reject P2 elements.

## Notes

- Fully handled by the C++ core (Python fallback only for buffers or a
  non-default `float_fmt`).
- No reference fixture exists under `tests/meshes/mfm/`; tests round-trip
  every supported linear type and explicitly verify `WriteError` on a
  mixed-type mesh.
