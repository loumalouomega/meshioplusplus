# PERMAS (`.post`, `.dato`)

The [PERMAS](https://www.intes.de) data-file format: `$`-delimited keyword sections in plain text, optionally gzip-compressed.

| | |
|---|---|
| **Format name** | `permas` |
| **Extensions** | `.post`, `.post.gz`, `.dato`, `.dato.gz` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.post")
meshioplusplus.permas.write("out.post", mesh)
```

`write` takes no keyword arguments.

## File structure

Lines starting `!` are comments and are skipped. A keyword line is `$KEYWORD[=...]` (leading `$` stripped, upper-cased for matching).

- `$COOR...` → node block: lines `<gid> <x> <y> <z...>` (split on single spaces — multiple consecutive spaces create empty tokens, which are not filtered here), read until a `$`/`!` line, building a `gid → running index` map.
- `$ELEMENT TYPE=<permas type>` → element block: lines are `<running_index> <node gid> <node gid> ...`, resolved through the node-gid map. PERMAS uses a **trailing `!` as a line-continuation marker**: node ids accumulate across lines until one does *not* end in `!`, at which point the accumulated ids are emitted as one completed cell.
- `$NSET`/`$ESET` blocks are parsed (including `GENERATE`, expanding exactly 3 ids `start,end,step` via `np.arange` with an **exclusive** stop, i.e. literal Python `range` semantics rather than an inclusive PERMAS/Abaqus- style generate) but **the parsed sets are never attached to the returned Mesh** — this is a currently-dead read path.
- Everything else is silently ignored ("too many PERMAS keywords to explicitly skip").

Write emits: `!PERMAS DataFile Version 18.0`, a `! Written by meshio++ v<version>` provenance line (character-identical between the C++ and Python writers — see [`doc/formats.md`](../formats.md#provenance)), `$ENTER COMPONENT NAME=DFLT_COMP`, `$STRUCTURE`, `$COOR`, then node rows (sequential 1-based index, not the original PERMAS gid — none is tracked on write), then per element type a `!` separator + `$ELEMENT TYPE=...` + rows (`<running_id> <n1+1> ...`, id counting continuously across all cell blocks), then `$END STRUCTURE` / `$EXIT COMPONENT` / `$FIN`.

## Cell types

| PERMAS | meshio++ | PERMAS | meshio++ |
|---|---|---|---|
| `PLOT1` | `vertex` | `HEXFO8` | `hexahedron` |
| `FSCPIPE2` (and 10 other beam/rod names) | `line` | `HEXE20` | `hexahedron20` |
| `PLOTL3` | `line3` | `HEXE27` | `hexahedron27` |
| `TRIMS3` (and 6 other names) | `triangle` | `TET4` | `tetra` |
| `TRIMS6` | `triangle6` | `TET10` | `tetra10` |
| `SHELL4` (and 5 other names) | `quad` | `PYRA5` | `pyramid` |
| `QUAMS8` | `quad8` | `PENTA6` | `wedge` |
| `QUAMS9` | `quad9` | `PENTA15` | `wedge15` |

(Several PERMAS names collapse onto one meshio++ type on read; the write-side reverse map picks one canonical name per meshio++ type as shown, determined by Python dict-insertion-order "last wins" — the C++ port hardcodes the resulting map directly to guarantee it matches.)

Write-only node-order permutation for second-order elements (no analogous reorder exists on read — see quirks):

| type | permutation |
|---|---|
| `triangle6` | `[0, 3, 1, 4, 2, 5]` |
| `tetra10` | `[0, 4, 1, 5, 2, 6, 7, 8, 9, 3]` |
| `quad9` | `[0, 4, 1, 7, 8, 5, 3, 6, 2]` |
| `wedge15` | `[0, 6, 1, 7, 2, 8, 9, 10, 11, 3, 12, 4, 13, 5, 14]` |

## Data mapping

None — PERMAS produces no `point_data`/`cell_data`/`field_data` at all; `$NSET`/`$ESET` are parsed but discarded (see above).

## Quirks & limitations

- **Asymmetric quadratic-element round-trip**: the four second-order node reorder tables above are applied **only on write** — there is no inverse permutation on read. A file written by meshio++ and read back by meshio++'s own reader would therefore not restore the original meshio++ node order for `triangle6`/`tetra10`/`quad9`/`wedge15` without external correction.
- `$NSET`/`$ESET` sets are read into local dicts but never surface on the returned `Mesh` — effectively dead code in the current reader.
- The `!`-continuation parsing is essential to get right: a standalone `!` separator line between element blocks must yield **no** cell (not an empty one) — getting this wrong previously caused an out-of-bounds crash on a multi-block mesh during development of the C++ port, now fixed and covered by tests.
- `GENERATE` set expansion uses exclusive-stop `np.arange(start, end, step)` semantics, not an inclusive range — a possible off-by-one relative to true PERMAS/Abaqus-style `GENERATE` (moot in practice since sets are discarded regardless).

## Notes

- The C++ core handles the plain-text `.post`/`.dato` files. The gzip `.gz` containers fall back to the Python implementation.
- No reference fixture exists under `tests/python/meshes/permas/`; tests round-trip synthetic meshes and check `.post`/`.post.gz`/`.0.post`/`.0.post.gz` extension/suffix dispatch.
