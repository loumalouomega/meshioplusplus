# DOLFIN XML (`.xml`)

The legacy [DOLFIN/FEniCS](https://manpages.ubuntu.com/manpages/jammy/en/man1/dolfin-convert.1.html) XML mesh format. A file holds one mesh (triangle or tetrahedron only); each data array associated with that mesh lives in its own sibling file.

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

Vertices and cells are placed **by their `index` attribute**, not by document order — the Python reader streams the file with `ElementTree.iterparse` and clears each element after processing (memory-light, no full-DOM parse).

Each `point_data` or `cell_data` array is stored in a **separate sibling file**, not inline: for a mesh file `mesh.xml` and a key `"a"`, the file is `mesh_a.xml`, matched by the reader via the regex `"{stem}_([^.]+)\.xml"`. Its content:

```xml
<dolfin>
  <mesh_function type="int|uint|float" dim="D" size="N">
    <entity index="0" value="..." />
    ...
  </mesh_function>
</dolfin>
```

`type` is derived from the numpy dtype family (`int`, `uint`, or `float`); the C++ writer only distinguishes float vs. integer (no separate `uint`), a minor naming difference from the Python writer that doesn't affect numeric round-trip.

The `dim` attribute is the topological dimension of the entities the function is defined on, and **that is the whole point/cell discriminator** (v9.9.0): `dim="0"` means *vertices*, so a `dim="0"` file is read as `point_data` and anything else as `cell_data`. Point data therefore needs no new file convention — it uses the format's own notion.

For a **cell** function the value written is not the actual topological dimension but a z-flatness heuristic: `2` if the mesh is 2-D or all point z-coordinates are exactly zero (checked with `np.allclose(..., atol=1e-14)`), else `3`. That is a pre-existing quirk of both writers, kept for compatibility; it never collides with `0`, which is why the point/cell split is unambiguous.

## Cell types

`triangle` and `tetra` only — no node reordering (DOLFIN's node order matches meshio++'s).

## Data mapping

- `point_data["<name>"]` — one sibling `<stem>_<name>.xml` file per key, written with `dim="0"` (v9.9.0; previously dropped on write).
- `cell_data["<name>"]` — one sibling `<stem>_<name>.xml` file per key, written with the z-flatness `dim` above. The reader scans the mesh file's directory for matches and routes each by its `dim`.
- No field_data.

## Quirks & limitations

- If a mesh has both `triangle` and `tetra` cells, the writer prefers `tetra` and discards everything else with a warning — DOLFIN XML stores exactly one cell type per mesh.
- Writing always emits the warning `"DOLFIN XML is a legacy format. Consider using XDMF instead."`
- Each cell-data XML file supports exactly one `<mesh_function>`; a file with more than one raises `ReadError`.
- The `dim` heuristic for cell-data files (2D-or-all-z-zero → 2, else 3) is a z-flatness check, not a request for the actual topological dimension of the data. It never yields `0`, so it cannot be confused with a point (`dim="0"`) function.
- A name present in **both** `point_data` and `cell_data` wants the same sibling file. Cell data has always owned it, so the cell array is written and the point array is skipped with a warning rather than one clobbering the other.
- A **non-scalar** `point_data` array is skipped with a warning: a `mesh_function` holds one value per entity, so there is nowhere to put the extra components.
- **Triangles and tetrahedra only is correct by format**, not a meshio++ limitation — DOLFIN XML is simplicial. Both writers raise a named `WriteError` rather than silently dropping the unsupported blocks.

## Notes

- Fully handled by the C++ core (via the vendored pugixml + `std::filesystem` for the directory scan) — no Python fallback path is needed for this format.
- No reference fixture exists under `tests/python/meshes/dolfin/`; tests round-trip synthetic meshes (`tri_mesh`, `tri_mesh_2d`, `tet_mesh`) plus cell-data and point-data variants. `Dolfin.PointDataRoundTripsAsADimZeroMeshFunction` asserts the written `dim="0"` in the raw file rather than only reading our own output back — the latter would pass even with a wrong `dim`, since the reader would simply put the array in `cell_data` and the values would still be there.
