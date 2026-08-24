# Split (partition into multiple meshes)

`meshioplusplus.split(mesh, by=...)` partitions one mesh into several submeshes — by **cell type**, by **connected component**, by **region** (a named `cell_sets` group, or an integer `cell_data` tag), or by **named Cell regions** (`by="regions"`, plural) — each pruned to its own points with connectivity and data remapped. It is a mesh **operation** (like [crop](/crop) and [merge](/merge)), not a file format, and uses only standard C++/numpy, so it runs under every mesh backend.

![A mesh split into regions, coloured by piece](/images/split_regions.png)

```python
import meshioplusplus

mesh = meshioplusplus.read("mixed.vtu")

pieces = meshioplusplus.split(mesh, by="type")        # {"triangle": ..., "quad": ...}
pieces = meshioplusplus.split(mesh, by="component")   # {"0": ..., "1": ...}
pieces = meshioplusplus.split(mesh, by="region")      # by cell_sets / tag
pieces = meshioplusplus.split(mesh, by="regions")     # one piece per named Cell region

for key, piece in pieces.items():
    meshioplusplus.write(f"part_{key}.vtu", piece)
```

`split` returns an ordered `dict` `{key: Mesh}`.

## Criteria

| `by` | partitions by | keys |
|---|---|---|
| `"type"` | cell type (one piece per type) | the meshio type name (`"triangle"`, …) |
| `"component"` | connected components (flood-fill over node-sharing adjacency) | `"0"`, `"1"`, … |
| `"region"` | named `cell_sets` if present, else an integer `cell_data` tag | set name / tag value |
| `"regions"` | named **Cell** [regions](/regions) (cross-binding — runs in the C++ core, not just the Python shim) | region name |

For `by="region"`, if the mesh carries `cell_sets` they define the pieces (one per set, keyed by name); otherwise the split falls back to an integer `cell_data` tag — `tag="name"` selects it, or the first integer cell_data is auto-detected. Connected components use a union-find over cells sharing a node (never O(N²)).

`by="regions"` (plural) is the odd one out: **it is not a partition**. Each named `Cell` region is processed independently, so a cell belonging to two regions lands in both output pieces, and a cell in no region lands in none. `Point`/`Side` regions produce no piece at all — there is no sound default for turning "these facets" or "these points" alone into a whole submesh without first deciding which surrounding cells to pull in, so `split` does not guess. Unlike `by="region"`, this criterion runs in the C++ core and so reaches every binding uniformly (see [Other languages](#other-languages)), not just Python; `--tag` is unused by it. Use [`meshioplusplus regions`](/cli#regions) first to see what a mesh's regions are before splitting by them.

## What changes

Each piece is **pruned to only its own points**, with connectivity remapped and `point_data` / `cell_data` subset to it (`field_data` carried through). `point_sets` / `cell_sets` are remapped into each piece (done in the Python layer). The numpy fallback handles rectangular cell blocks; ragged/polyhedron blocks are handled by the C++ core only.

## CLI

```bash
meshioplusplus split in.vtu 'out_{key}.vtu' --by type
meshioplusplus split in.vtu 'out_{key}.vtu' --by component
meshioplusplus split in.vtu 'out_{key}.vtu' --by region --tag gmsh:physical
meshioplusplus split in.vtu 'out_{key}.vtu' --by regions
```

The output pattern must contain `{key}`, which is replaced by each piece's key. The number of pieces and their sizes are printed. See the [CLI reference](/cli).

## Other languages

The C++ core returns the pieces (the file-writing pattern is a caller concern). `by="region"` there means the integer `cell_data` tag (`cell_sets` are Python-only); `by="regions"` (plural) means named Cell regions, and reaches every binding below identically, since `split`'s `by` string is resolved once by the shared `split_by_from_name`:

- **C API** — `mio_split(mesh, by, tag_name)` returns a `mio_split_result` (`mio_split_result_count`, `..._key`, `..._mesh` / `..._take_mesh`, `..._free`). See the [C API reference](/c_api).
- **Fortran** — `mesh%split(by, tag_name=..., keys=...)` returns an array of meshes. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `split(mesh, by, tagName)` returns an array of `{key, mesh}`. See the [WebAssembly reference](/wasm).
