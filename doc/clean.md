# Clean (weld / prune / de-dup)

`meshioplusplus.clean(mesh, ...)` cleans a mesh in one pass: **weld** coincident points, **drop degenerate** cells, **drop exact-duplicate** cells, and **remove orphaned** (unused) points — each step individually toggleable — remapping connectivity and data throughout. It is a mesh **operation** (like [merge](/merge) and [reordering](/reorder)), not a file format, and uses only standard C++/numpy, so it runs under every mesh backend.

![Coincident duplicate points welded back down to the original point count](/images/clean_weld.png)

```python
import meshioplusplus

mesh = meshioplusplus.read("messy.vtu")

out = meshioplusplus.clean(mesh)                        # prune + de-dup, no weld
out = meshioplusplus.clean(mesh, weld=True, atol=1e-8)  # also fuse coincident pts

# get a report of what changed
out, report = meshioplusplus.clean(mesh, weld=True, return_report=True)
print(report)   # points_welded / points_removed_orphan / cells_dropped_*
```

## Options

| option | default | effect |
|---|---|---|
| `weld` | `False` | fuse points within `atol` (spatial-hash bucket grid, never O(N²)) |
| `atol` | `1e-8` | absolute coincidence tolerance for welding |
| `remove_orphans` | `True` | drop points referenced by no surviving cell |
| `drop_degenerate` | `True` | drop cells with a repeated corner node or a near-zero measure |
| `drop_duplicate_cells` | `True` | drop exact-duplicate cells (keep-first, per block) |

The pass runs **weld → drop degenerate → drop duplicate → remove orphans**. On a weld, the surviving point keeps the value of the *first* contributing point (keep-first tie-break).

## What changes

- **Points** are welded/pruned; connectivity is remapped.
- **`point_data`** follows the surviving points (keep-first on weld); **`cell_data`** follows the surviving cells; **`field_data`** is preserved.
- **`point_sets` / `cell_sets`** are remapped; entries pointing at removed points/cells are dropped (done in the Python layer).

The numpy fallback handles rectangular cell blocks; ragged/polyhedron blocks are handled by the C++ core only.

## CLI

```bash
meshioplusplus clean in.vtu out.vtu                     # default set (no weld)
meshioplusplus clean in.vtu out.vtu --weld --atol 1e-6
meshioplusplus clean in.vtu out.vtu --remove-orphans --drop-degenerate
```

With no step flags the default set runs (remove-orphans + drop-degenerate + drop-duplicates, **no** weld). Passing any step flag selects exactly the requested steps. A removal summary is printed. See the [CLI reference](/cli).

## Other languages

- **C API** — `mio_clean(mesh, weld, atol, remove_orphans, drop_degenerate, drop_duplicate_cells, &n_welded, &n_orphan, &n_degen, &n_dup)`. See the [C API reference](/c_api).
- **Fortran** — `mesh%clean(weld=..., atol=..., points_welded=..., ...)`. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `clean(mesh, weld, atol, removeOrphans, dropDegenerate, dropDuplicateCells)` returns `{mesh, pointsWelded, ...}`. See the [WebAssembly reference](/wasm).
