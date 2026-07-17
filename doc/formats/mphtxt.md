# COMSOL text mesh (`.mphtxt`)

The [COMSOL Multiphysics](https://www.comsol.com) `.mphtxt` text mesh format (as handled by [FEconv](https://github.com/victorsndvg/FEconv)). ASCII, storing a version, tag/type name tables and one or more mesh objects; comments run from `#` to end of line.

| | |
|---|---|
| **Format name** | `mphtxt` |
| **Extensions** | `.mphtxt` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.mphtxt")
meshioplusplus.mphtxt.write("out.mphtxt", mesh)
```

`write` takes no keyword arguments.

## File structure

The whole file (after stripping `#`-comments) is read as a flat token stream:

1. Version major, version minor (2 ints, discarded).
2. `n_tags` length-prefixed tag-name strings (discarded).
3. `n_types` length-prefixed type-name strings — this count is also the number of "object" records that follow.
4. For each object (**only the first is actually parsed** — the reader `break`s after it):
   - 3 ints (object type indices, discarded), a class-name string (expected `"Mesh"`), an object-version int, `sdim` (spatial dimension), `n_points`, `lowest` (the file's lowest node index — normally 1, but not assumed to be).
   - `n_points * sdim` floats — node coordinates, row-major `(n_points, sdim)`.
   - `n_eltypes` element-type blocks, each:
     - a COMSOL type-name string (`"tet"`, `"tri2"`, …);
     - `n_nodes`, `n_elem`;
     - `n_elem * n_nodes` ints — connectivity, shifted by `-lowest` to become 0-based;
     - `n_par_per`, `n_par`, then `n_par * n_par_per` discarded parameter tokens;
     - `n_geom`, then `n_geom` ints — one **geometric entity index** per element → `cell_data["mphtxt:geom"]`;
     - `n_ud`, then `n_ud * 2` discarded up/down-pair ints.

## Cell types

Hybrid (multi-type) meshes are supported. COMSOL ↔ meshio++ type map:

| COMSOL | meshio++ | COMSOL | meshio++ |
|---|---|---|---|
| `vtx` | `vertex` | | |
| `edg` | `line` | `edg2` | `line3` |
| `tri` | `triangle` | `tri2` | `triangle6` |
| `quad` | `quad` | `quad2` | `quad9` |
| `tet` | `tetra` | `tet2` | `tetra10` |
| `prism` | `wedge` | `prism2` | `wedge18` |
| `pyr` | `pyramid` | | |
| `hex` | `hexahedron` | `hex2` | `hexahedron27` |

Node-order permutation (COMSOL ↔ meshio++, self-inverse swaps — identity for every other type):

| type | permutation |
|---|---|
| `quad` | `[0, 1, 3, 2]` |
| `hexahedron` | `[0, 1, 3, 2, 4, 5, 7, 6]` |

## Data mapping

- `cell_data["mphtxt:geom"]` — per-element geometric entity index, one array per cell block.

## Quirks & limitations

- Only the **first mesh object** in the file is read; subsequent objects are ignored entirely — a hard limitation for multi-mesh `.mphtxt` files.
- Element parameter values (`n_par`) and up/down topology-linking pairs (`n_ud`) are always discarded on read and always written as empty/zero on write — any parametrization or entity-linking data COMSOL stores there does not round-trip.
- The connectivity shift uses the file's actual `lowest` value (not a hardcoded `1`), so files with an unusual lowest-index convention are handled correctly.
- Unsupported cell types on write: the Python writer `warn`s and skips them; the C++ writer raises `WriteError`, which forces the Python fallback.

## Notes

- Fully handled by the C++ core.
- No reference fixture exists under `tests/meshes/mphtxt/`; tests round-trip linear, second-order, and hybrid (`tri_quad_mesh`) meshes.
