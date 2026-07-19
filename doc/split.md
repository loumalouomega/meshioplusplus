# Split (partition into multiple meshes)

`meshioplusplus.split(mesh, by=...)` partitions one mesh into several submeshes —
by **cell type**, by **connected component**, or by **region** (a named
`cell_sets` group, or an integer `cell_data` tag) — each pruned to its own points
with connectivity and data remapped. It is a mesh **operation** (like
[crop](/crop) and [merge](/merge)), not a file format, and uses only standard
C++/numpy, so it runs under every mesh backend.

```python
import meshioplusplus

mesh = meshioplusplus.read("mixed.vtu")

pieces = meshioplusplus.split(mesh, by="type")        # {"triangle": ..., "quad": ...}
pieces = meshioplusplus.split(mesh, by="component")   # {"0": ..., "1": ...}
pieces = meshioplusplus.split(mesh, by="region")      # by cell_sets / tag

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

For `by="region"`, if the mesh carries `cell_sets` they define the pieces (one
per set, keyed by name); otherwise the split falls back to an integer `cell_data`
tag — `tag="name"` selects it, or the first integer cell_data is auto-detected.
Connected components use a union-find over cells sharing a node (never O(N²)).

## What changes

Each piece is **pruned to only its own points**, with connectivity remapped and
`point_data` / `cell_data` subset to it (`field_data` carried through).
`point_sets` / `cell_sets` are remapped into each piece (done in the Python
layer). The numpy fallback handles rectangular cell blocks; ragged/polyhedron
blocks are handled by the C++ core only.

## CLI

```bash
meshioplusplus split in.vtu 'out_{key}.vtu' --by type
meshioplusplus split in.vtu 'out_{key}.vtu' --by component
meshioplusplus split in.vtu 'out_{key}.vtu' --by region --tag gmsh:physical
```

The output pattern must contain `{key}`, which is replaced by each piece's key.
The number of pieces and their sizes are printed. See the [CLI reference](/cli).

## Other languages

The C++ core returns the pieces (the file-writing pattern is a caller concern),
and `by="region"` there means the integer `cell_data` tag (`cell_sets` are
Python-only):

- **C API** — `mio_split(mesh, by, tag_name)` returns a `mio_split_result`
  (`mio_split_result_count`, `..._key`, `..._mesh` / `..._take_mesh`,
  `..._free`). See the [C API reference](/c_api).
- **Fortran** — `mesh%split(by, tag_name=..., keys=...)` returns an array of
  meshes. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `split(mesh, by, tagName)` returns an array of
  `{key, mesh}`. See the [WebAssembly reference](/wasm).
