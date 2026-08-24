# Merge / combine (merge & weld)

`meshioplusplus.merge([mesh_a, mesh_b, ...], ...)` combines an ordered list of two or more meshes into a **single** mesh. It is a mesh **operation** (like [reordering](/reorder), [comparison](/diff), [surface extraction](/extract_surface) and [quality metrics](/mesh_quality)), not a file format, and uses only standard C++/numpy — no third-party dependencies — so it runs under every mesh backend.

![Two meshes combined, coloured by which input each cell came from](/images/merge_sources.png)

```python
import meshioplusplus

a = meshioplusplus.read("block_a.vtu")
b = meshioplusplus.read("block_b.vtu")

# concatenate (default)
combined = meshioplusplus.merge([a, b])

# weld coincident nodes on the shared interface
welded = meshioplusplus.merge([a, b], weld=True, atol=1e-8)
```

## Concatenate vs. weld

**Concatenation** (the default) is a pure, O(N) append:

- all points are appended; each input's connectivity is offset by the running point count so indices stay valid;
- cell blocks of the same [cell type](/formats) are merged into one block, in input order;
- `point_data` / `cell_data` are combined per the [data policy](#data-policy);
- a per-cell `source_mesh_id` cell_data records each cell's originating input (see [provenance](#provenance-tag)).

**Welding** (`weld=True`) additionally detects points that are coincident **across all inputs** within an absolute tolerance and merges them into one, remapping connectivity. Coincidence is found with a **spatial-hash bucket grid** (coordinates quantized to `atol`-sized cells, then a 3×3×3 neighbourhood search) — never an O(N²) scan. This is the standard way to stitch adjacent blocks that share an interface into a watertight mesh.

```python
# two adjacent bricks sharing a face -> one mesh, interface nodes fused
mesh = meshioplusplus.merge([left, right], weld=True, atol=1e-9)
```

## Tolerance

`atol` (default `1e-8`) is the **absolute** coincidence distance: two points are welded when the Euclidean distance between them is `<= atol`. Points farther apart than `atol` are kept distinct. Pick `atol` well below your smallest real feature size and above floating-point/round-trip noise.

## Data policy

A `point_data` / `cell_data` key present in only **some** inputs is handled by `data_policy`:

| `data_policy` | behaviour |
|---|---|
| `"intersection"` (default) | keep only keys present in **every** input (with matching column count); partial keys are dropped (with a warning) |
| `"fill"` | keep every key; rows contributed by inputs lacking it are filled with **NaN** (the array is emitted as `float64`) |

A key present in all inputs is concatenated into one array. When inputs disagree on a key's dtype the result is promoted to `float64`.

### Weld tie-break

When welding fuses points that carry `point_data`, the surviving point keeps the value of the **first** contributing point — the earliest input, lowest original index (a fixed, deterministic *keep-first* rule).

## Provenance tag

By default `merge` attaches an `int64` `source_mesh_id` cell_data whose value is the index (0-based, in input order) of the mesh each output cell came from. Turn it off with `source_tag=False` (CLI `--no-source-tag`).

## Dropping duplicate cells

With `drop_duplicate_cells=True` (only meaningful together with `weld`), cells that become **identical** after welding — same connectivity as an earlier cell in the same output block — are removed, keeping the first. Useful when overlapping inputs share not just an interface but whole coincident elements.

## Sets, regions and field data

Overlapping **names** from different inputs are **namespaced by source id**: on a collision the name is prefixed with the input index (`"0:boundary"`, `"1:boundary"`); a name unique to one input is kept as-is. This applies to `point_sets`, `cell_sets` and `field_data` keys, so nothing is silently overwritten. `point_sets` / `cell_sets` live only on the Python `Mesh`, so they are remapped in the Python layer (using the point/cell index maps the core returns) and are **not** carried by the C / Fortran / WebAssembly bindings.

## What is carried through

Points, cells (rectangular and ragged/polyhedron blocks), `point_data`, `cell_data`, `field_data`, `point_sets` and `cell_sets`. `mesh.info` and `gmsh_periodic` are not.

## CLI

The last positional argument is the **output**; everything before it is a variadic list of inputs:

```bash
meshioplusplus merge a.vtu b.vtu c.vtu out.vtu               # concatenate
meshioplusplus merge a.vtu b.vtu out.vtu --weld --atol 1e-8  # weld coincident nodes
meshioplusplus merge a.vtu b.vtu out.vtu --data-policy fill  # keep partial data (NaN fill)
meshioplusplus merge a.vtu b.vtu out.vtu --no-source-tag     # skip provenance cell_data
```

Unless `--quiet`, a short summary is printed (total points/cells in, points welded, points/cells out). See the [CLI reference](/cli).

## Other languages

The operation is exposed across every binding surface:

- **C API** — `mio_merge(meshes, count, weld, atol, source_tag, data_policy, drop_duplicate_cells)` takes an array of mesh handles and returns a new mesh handle (`data_policy` `0` = intersection, `1` = fill; free with `mio_mesh_free`). See the [C API reference](/c_api).
- **Fortran** — `mio_merge(meshes, weld=..., atol=..., source_tag=..., data_policy=..., drop_duplicate_cells=...)` takes an array of `type(mio_mesh)` and returns a new `type(mio_mesh)`. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `merge([meshA, meshB, ...], weld, atol, sourceTag, dataPolicy, dropDuplicateCells)` returns a new mesh object. See the [WebAssembly reference](/wasm).
