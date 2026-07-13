# Abaqus (`.inp`)

The [Abaqus](http://abaqus.software.polimi.it/v6.14/index.html) input-deck
format: keyword-driven ASCII (`*NODE`, `*ELEMENT`, `*NSET`, `*ELSET`,
`*INCLUDE`, …).

| | |
|---|---|
| **Format name** | `abaqus` |
| **Extensions** | `.inp` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("model.inp")
meshio.abaqus.write("out.inp", mesh,
    float_fmt=".16e",
    translate_cell_names=True,
)
```

- **`float_fmt`** — coordinate format.
- **`translate_cell_names`** — if `True` (default), meshio cell types are
  mapped to an Abaqus element name via the lookup table below; if `False`,
  the meshio type string is written verbatim as `TYPE=`.

## File structure

Comment lines start `**`; keyword lines start `*`. A keyword is extracted as
`line.partition(",")[0].strip().replace("*","").upper()`. Recognized:
`NODE`, `ELEMENT`, `NSET`, `ELSET`, `INCLUDE`; anything else is simply
skipped line-by-line (not block-skipped) until the next `*` line.

- `*NODE`: comma-separated `id, x, y, [z]` rows.
- `*ELEMENT, TYPE=<abaqus_type>[, ELSET=<name>]`: comma-separated integer
  rows flattened into one token stream, then split into fixed-width records
  of `(node count for TYPE) + 1` (element id + node ids). An inline
  `ELSET=` also registers a `cell_sets` entry spanning that whole block.
- `*NSET, NSET=<name>[, GENERATE]` / `*ELSET, ELSET=<name>[, GENERATE]`:
  member ids, or (with `GENERATE`) exactly 3 ids `start, end, step` expanded
  via `np.arange`. `*ELSET` may also reference **other set names**
  (non-numeric first token), resolved recursively.
- `*INCLUDE, INPUT=<path>`: recursively reads and merges another `.inp` file
  (path resolved relative to the current file's directory if not found
  as-is).

## Cell types

The Abaqus element-type table is large (trusses, beams, shells, solids); a
representative excerpt:

| Abaqus | meshio | Abaqus | meshio |
|---|---|---|---|
| `T2D2`, `T3D2`, `B21`, `B31` | `line` | `C3D8`, `C3D8R`, `S4`, `CPS4` | `hexahedron`* / `quad` |
| `T2D3`, `T3D3`, `B22`, `B32` | `line3` | `C3D20`, `C3D20R` | `hexahedron20` |
| `CPS3`, `STRI3`, `S3` | `triangle` | `C3D4` | `tetra` |
| `STRI65`, `CPE6` | `triangle6` | `C3D4H` | `tetra4`** |
| `S8R`, `S8R5` | `quad8` | `C3D10`, `C3D10M` | `tetra10` |
| `S9R5` | `quad9` | `C3D6` | `wedge` |
| | | `C3D15` | `wedge15` |

(*`C3D8*` → `hexahedron`, `S4`/`CPS4`/etc. → `quad`; both map to distinct
meshio types depending on whether the card is a solid or shell element.
**`C3D4H` maps to the type string `"tetra4"`, not `"tetra"` — this is an
asymmetric entry relative to `C3D4`→`tetra` and doesn't round-trip through
meshio's standard type vocabulary; noted here as a known table quirk rather
than a deliberate feature.)

The reverse map (meshio → Abaqus) is lossy: several Abaqus names collapse to
one meshio type, so the writer always emits whichever Abaqus name happens to
be *last* in the internal table for that meshio type — the originating
keyword is not preserved through a read→write round trip.

## Data mapping

- `point_sets` / `cell_sets` — from `*NSET`/`*ELSET` (including inline
  `ELSET=` on `*ELEMENT`), keyed by set name.
- No point_data/cell_data/field_data is produced — the Abaqus reader never
  populates them.

## Quirks & limitations

- `translate_cell_names=False` bypasses the lookup table entirely and writes
  the meshio type string as-is; useful for pass-through of unsupported types,
  but has no C++ equivalent (the C++ writer always looks up the table and
  throws if the type isn't found).
- `GENERATE` sets require exactly 3 numbers (`start,end,step`) or raise
  `ReadError`.
- `*ELSET` can reference other set names transitively, including elsets
  implicitly created by an inline `ELEMENT ... ELSET=`.
- The C++ reader explicitly refuses `NSET`/`ELSET`/`INCLUDE` keywords with a
  hard error, deferring the entire file to the Python reader whenever any of
  them appear.
- The C++ writer is only attempted when `float_fmt == ".16e"`,
  `translate_cell_names == True`, and the mesh has **no** `point_sets`/
  `cell_sets`.
- The C++ reader's Abaqus-type lookup tries the upper-cased type string
  first, then falls back to the as-written (case-sensitive) string — a
  leniency the Python reader (a plain case-sensitive dict lookup) doesn't
  have.

## Notes

- `tests/meshes/abaqus/UUea.inp` — point-sum 4950.0, 50 cells, 10 cell_sets.
- `nle1xf3c.inp` (MIT Abaqus course material) — point-sum ≈32.215275528, 12
  cells, 3 cell_sets.
- `element_elset.inp` — exercises inline `ELSET=` on `*ELEMENT`.
- `wInclude_main.inp` + `wInclude_bulk.inp` — exercises `*INCLUDE`.
