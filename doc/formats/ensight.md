# EnSight Gold (`.case` / `.geo`)

The [EnSight Gold](https://vis.lbl.gov/archive/NERSC/Software/ensight/doc/OnlineHelp/UM-C11.pdf) format: a small `.case` index file plus a Gold geometry file (conventionally `.geo`), the de-facto exchange format for EnSight/ParaView post-processing. meshio++ handles the mesh-geometry subset in both ASCII and C-binary form.

| | |
|---|---|
| **Format name** | `ensight` |
| **Extensions** | `.case`, `.geo` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("out.case")          # or the .geo directly
meshioplusplus.ensight.write("out.case", mesh, binary=True)
```

- **`binary`** — write C-binary geometry (default `True`); `False` writes ASCII (`%12.5e` floats — EnSight caps usable precision at 6 significant digits).

Either path (`.case` or `.geo`) selects the sibling pair on write; both files are always written. On read, a `.case` path is parsed for its `GEOMETRY`/`model:` entry (resolved relative to the case file's directory), while any other path is treated as a Gold geometry file directly.

## File structure

**`.case`**: only the `FORMAT` (must declare `type: ensight gold`) and `GEOMETRY` (`model: [ts] [fs] <file>.geo`) sections are consumed; `VARIABLE`/`TIME` sections are ignored (mesh-only scope). Transient wildcard geometry names (`model: name.****.geo`) are rejected.

**`.geo`** (Gold): 2 description records, `node id <off|given|assign|ignore>`, `element id <...>`, an optional `extents` block, then per part: `part`, part number, description, `coordinates`, node count, optional node-id array, and the X, Y, Z coordinate arrays (blocked), followed by element sections (`tria3`, `tetra4`, ..., each with a count, an optional element-id array, and connectivity). Binary files start with an 80-char `"C Binary"` record; all strings are 80-char records and all numbers 32-bit ints/floats in the writer's native byte order. `"Fortran Binary"` files are rejected.

Per the Gold specification, **connectivity is positional** (1-based index into the part's coordinate list); `given`/`ignore` id arrays are present in the file but skipped — the same interpretation VTK's EnSight readers use.

The writer emits a single part (id 1, description `Mesh`) with `node id assign` / `element id assign` (no id arrays) and always writes three coordinate components (z padded with 0 for 2D meshes).

## Cell types

| EnSight keyword | meshio++ type |
|---|---|
| `point` | `vertex` |
| `bar2` / `bar3` | `line` / `line3` |
| `tria3` / `tria6` | `triangle` / `triangle6` |
| `quad4` / `quad8` | `quad` / `quad8` |
| `tetra4` / `tetra10` | `tetra` / `tetra10` |
| `pyramid5` / `pyramid13` | `pyramid` / `pyramid13` |
| `penta6` / `penta15` | `wedge` / `wedge15` |
| `hexa8` / `hexa20` | `hexahedron` / `hexahedron20` |
| `nsided` | `polygon` (read only) |
| `nfaced` | `polyhedron<N>` (read only) |

Node ordering matches meshio for every type except `penta15`, which differs from `wedge15` by the involution `{0,2,1,3,5,4, 8,7,6, 11,10,9, 12,14,13}` (the same map VTK's EnSight readers apply; note that Kratos's penta15/hexa20 permutations fix Kratos-specific ordering and do **not** apply here). Ghost-cell sections (`g_*`) are read as their base type.

## Data mapping

- Multi-part files are concatenated into one point array; every per-part element section becomes its own cell block.
- `cell_data["ensight:part"]` — the owning part number per cell (Int64), emitted only when the file has **two or more** parts, so single-part round-trips stay clean. The writer ignores it (always one part).
- The writer drops `point_data`/`cell_data`/`field_data` with a warning — EnSight variables live in separate per-variable files, which are out of scope (v1).

## Quirks & limitations

- The format spans two files and **cannot** be read from or written to a buffer.
- Foreign-endian binaries are auto-detected: the reader checks the plausibility of the part-number/count records (preferring the smaller of two plausible interpretations at the tiny part-number record) and byte-swaps all subsequent reads.
- Writing `nsided`/`nfaced` (polygon/polyhedron) blocks raises `WriteError` (v1); reading them is supported (`nfaced` cells are grouped by node count into `polyhedron<N>` blocks, the openfoam convention).
- Transient/moving geometry (`change_coords_only`, filename wildcards) and EnSight 5/6 files are not supported.
- The `.geo` extension is also used by Gmsh *script* files — those are not meshes and not claimed by meshio++'s gmsh reader, but a stray Gmsh `.geo` passed to `ensight` will simply fail to parse.

## Notes

- `tests/meshes/ensight/simple.case`/`simple.geo` — hand-authored two-part ASCII example (non-sequential `node id given` ids, `tetra4` + `tria3` + `nsided`) exercising part concatenation, positional connectivity, and `ensight:part` tagging.
- Fully handled by the C++ core (ASCII + binary, both endiannesses); the pure-Python reference implements the identical feature set for fallback platforms.
