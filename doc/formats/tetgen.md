# TetGen (`.node` / `.ele`, `.smesh`)

The [TetGen](https://wias-berlin.de/software/tetgen/fformats.html) mesh format: a pair of sibling files sharing a stem — `<name>.node` (points, attributes, boundary markers) and `<name>.ele` (tetrahedra, region attributes) — and the single-file `.smesh` surface TetGen meshes (its points and polygonal facets).

| | |
|---|---|
| **Format name** | `tetgen` |
| **Extensions** | `.node`, `.ele`, `.smesh` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.node")   # reads the .node/.ele pair
meshioplusplus.tetgen.write("out.node", mesh, float_fmt=".16e")
```

- **`float_fmt`** — coordinate format (default `".16e"`; the C++ fast path is only used at this default).

Either path (`.node` or `.ele`) selects the sibling pair. A `.smesh` path reads or writes the surface file instead (below).

## The `.smesh` surface

A `.smesh` (TetGen's [surface mesh](https://wias-berlin.de/software/tetgen/fformats.smesh.html) file, a `.poly` whose facets are plain polygons) holds a surface for TetGen to mesh: part 1 the nodes (a `.node` section, header and rows in one file), part 2 the facets, one line each, `n v1 .. vn [marker]` under a header `nfacets markers`, then the holes and the regions. It is what [`tetrahedralize(backend="tetgen")`](../tetrahedralize.md#external-meshers-tetgen-and-gmsh) hands TetGen (since roadmap 7.7.2), and it reads and writes in both engines:

- **Writing** takes `triangle`, `quad` and rectangular `polygon` blocks, in block order, and the cell data `tetgen:ref` as each facet's boundary marker (0 for a block without it); any other cell type is refused by name. Points are written at `%.16e`, so they read back bit-equal. The hole and region lists are written empty.
- **Reading** gives one cell block per run of facets with the same corner count (`triangle`, `quad`, or `polygon`), in file order, and the markers as `cell_data["tetgen:ref"]`; node attributes and point markers in part 1 are skipped, and the holes and regions (which describe a volume to mesh, not the surface) are not kept.

## File structure

**`.node`**: header `npoints dim nattrs nbmarkers` (`dim` must be 3), then rows `idx x y z attr1..attrN marker1..markerM`. The **node index base** (0 or 1) is auto-detected from the first row's `idx` value, and all indices are then required to be **exactly consecutive** from that base (`ReadError` otherwise).

**`.ele`**: header `ntets 4 nattrs`, then rows `idx n0 n1 n2 n3 attr1..attrK`. Connectivity is shifted by the `.node` file's detected index base (not hardcoded), so files using either 0- or 1-based numbering read correctly.

Write: for `.node`, attribute/marker keys are partitioned into at most one "ref" key (the first `point_data` key containing the substring `":ref"`, or else the first key present) plus the remaining keys as plain attributes; header comments record the attribute/marker names. Ref columns use plain `{}` (`str()`) formatting in the Python writer, not the float format string — the C++ writer instead special-cases exact-integer values to print as plain integers, falling back to `%.16e` otherwise, a refinement that can format float-valued refs slightly differently between the two backends. For `.ele`, the same "first `:ref`-containing key floats to the front" rule applies to `cell_data` keys, and **each `tetra` cell block gets its own header line with its own index counting restarting at 0** — a mesh with multiple `tetra` blocks would therefore produce duplicate element ids across blocks (TetGen expects globally unique ids).

## Cell types

`tetra` only in a `.node`/`.ele` pair — TetGen only ever represents tetrahedra by construction; `triangle`, `quad` and `polygon` facets in a `.smesh`.

## Data mapping

- `point_data["tetgen:attr{k}"]` — the k-th node attribute column (`tetgen:attr1`, `tetgen:attr2`, ...).
- `point_data["tetgen:ref"]` / `"tetgen:ref2"` / ... — boundary marker columns.
- `cell_data["tetgen:ref"]` / `"tetgen:ref2"` / ... — region attribute columns (a single-element list, since TetGen only has one cell block).

## Quirks & limitations

- The format spans two files and **cannot** be read from or written to a buffer.
- Node index base auto-detection plus the consecutive-numbering check means a `.node` file with gaps in its index sequence is rejected outright, not silently accepted with holes.
- Writing multiple `tetra` blocks resets the element-id counter per block — a genuine round-trip risk for multi-block meshes (uncommon in practice, since TetGen conventionally produces exactly one tetrahedron block).

## Notes

- `tests/python/meshes/tetgen/mesh.node` (89 points, 1 attribute column named `"moje_data"` + 1 boundary-marker column named `"medit:ref"`) and `mesh.ele` (304 tetrahedra, 1 region attribute `"medit:ref"`) — checked via `mesh.point_data["tetgen:ref"].sum() == 12` and `mesh.cell_data["tetgen:ref"][0].sum() == 373`.
- Fully handled by the C++ core (Python fallback only for buffers or a non-default `float_fmt`).
