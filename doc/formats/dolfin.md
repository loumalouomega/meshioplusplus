# DOLFIN XML (`.xml`)

The legacy [DOLFIN/FEniCS](https://manpages.ubuntu.com/manpages/jammy/en/man1/dolfin-convert.1.html)
XML mesh format. A file holds one mesh (triangle or tetrahedron only); each
cell-data array associated with that mesh lives in its own sibling file.

| | |
|---|---|
| **Format name** | `dolfin-xml` |
| **Extensions** | `.xml` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.xml")
meshioplusplus.dolfin.write("out.xml", mesh)
```

Both `read(filename)` and `write(filename, mesh)` take no keyword arguments.

## File structure

The main file:

```xml
<dolfin nsmap="{'dolfin': 'https://fenicsproject.org/'}">
  <mesh celltype="triangle|tetrahedron" dim="2|3">
    <vertices size="N">
      <vertex index="0" x="..." y="..." [z="..."] />
      ...
    </vertices>
    <cells size="M">
      <triangle index="0" v0="..." v1="..." v2="..." />
      <!-- or <tetrahedron index="0" v0=".." v1=".." v2=".." v3=".." /> -->
      ...
    </cells>
  </mesh>
</dolfin>
```

Vertices and cells are placed **by their `index` attribute**, not by document
order — the Python reader streams the file with `ElementTree.iterparse` and
clears each element after processing (memory-light, no full-DOM parse).

Each `cell_data` array is stored in a **separate sibling file**, not inline:
for a mesh file `mesh.xml` and a cell_data key `"a"`, the file is
`mesh_a.xml`, matched by the reader via the regex `"{stem}_([^.]+)\.xml"`.
Its content:

```xml
<dolfin>
  <mesh_function type="int|uint|float" dim="D" size="N">
    <entity index="0" value="..." />
    ...
  </mesh_function>
</dolfin>
```

`type` is derived from the numpy dtype family (`int`, `uint`, or `float`); the
C++ writer only distinguishes float vs. integer (no separate `uint`), a minor
naming difference from the Python writer that doesn't affect numeric
round-trip. The `dim` attribute is **not** the cell's topological dimension —
it's `2` if the mesh is 2D or all point z-coordinates are exactly zero
(checked with `np.allclose(..., atol=1e-14)`), else `3`.

## Cell types

`triangle` and `tetra` only — no node reordering (DOLFIN's node order matches
meshio++'s).

## Data mapping

- `cell_data["<name>"]` — one sibling `<stem>_<name>.xml` file per key; the
  reader scans the mesh file's directory for matches.
- No point_data, no field_data.

## Quirks & limitations

- If a mesh has both `triangle` and `tetra` cells, the writer prefers `tetra`
  and discards everything else with a warning — DOLFIN XML stores exactly one
  cell type per mesh.
- Writing always emits the warning `"DOLFIN XML is a legacy format. Consider
  using XDMF instead."`
- Each cell-data XML file supports exactly one `<mesh_function>`; a file with
  more than one raises `ReadError`.
- The `dim` heuristic for cell-data files (2D-or-all-z-zero → 2, else 3) is a
  z-flatness check, not a request for the actual topological dimension of the
  data.

## Notes

- Fully handled by the C++ core (via the vendored pugixml + `std::filesystem`
  for the directory scan) — no Python fallback path is needed for this format.
- No reference fixture exists under `tests/meshes/dolfin/`; tests round-trip
  synthetic meshes (`tri_mesh`, `tri_mesh_2d`, `tet_mesh`) plus a cell-data
  variant.
