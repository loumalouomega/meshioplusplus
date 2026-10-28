# Medit solution (`.sol`, `.solb`)

The solution files of the [Medit / GMF](./medit.md) family: fields stored beside a `.mesh`/`.meshb`, one row per vertex or per element, under libMeshb's `SolAt*` keywords. They are the size maps and metrics [MMG](https://www.mmgtools.org/) reads and writes (`mmg3d_O3 mesh.mesh -sol mesh.sol`), so this format is the file route to MMG with no link to it. ASCII (`.sol`) or binary (`.solb`), told apart by whether the name ends in `b`, as for the mesh files.

| | |
|---|---|
| **Format name** | `medit_sol` |
| **Extensions** | `.sol`, `.solb` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

A solution needs the mesh it belongs to for its element fields, so the usual route is the Medit reader's `sol` argument:

```python
import meshioplusplus

mesh = meshioplusplus.medit.read("mesh.mesh", sol=True)       # mesh.sol (else mesh.solb) beside it
mesh = meshioplusplus.medit.read("mesh.meshb", sol="metric.solb")
meshioplusplus.medit.write("out.mesh", mesh, sol=True)       # out.mesh and out.sol
meshioplusplus.medit_sol.attach(mesh, "other.sol")            # onto a mesh already in memory
```

- **Vertex fields** (`SolAtVertices`) become `point_data["medit:sol:<i>"]`, numbered in file order.
- **Element fields** (`SolAtEdges`, `SolAtTriangles`, `SolAtQuadrilaterals`, `SolAtTetrahedra`, `SolAtPrisms`, `SolAtHexahedra`, `SolAtPyramids`) become `cell_data["medit:sol:<i>"]`, `<i>` numbered within each section, on the blocks of that cell type in block order. The blocks of other types hold `NaN`. When two types give the same name different shapes, the second keeps it with a `:<cell type>` suffix.
- **`medit.write(..., sol=True)`** writes `<stem>.sol` beside an ASCII mesh and `<stem>.solb` beside a binary one.

Read on its own, a solution is a geometry-less mesh, as for [`mff`](./mff.md): one zero-width point per row of `SolAtVertices`, the vertex fields as point data, and the dimension in `field_data["medit:dimension"]`. Element sections are skipped with a warning naming the `sol` route. Writing that mesh back gives the same file.

```python
sol = meshioplusplus.read("metric.sol")        # sol.points.shape == (n, 0)
meshioplusplus.write("copy.solb", sol)
```

Both engines read and write both encodings and write the same bytes. The C++ core has `read_medit_sol`, `read_medit_sol_into` and `write_medit_sol` (`formats/medit.hpp`), registered as `medit_sol`, so the native CLI, the C API and WASM convert solutions too.

## File structure

```
MeshVersionFormatted 2        # 0/1 -> float32 values, 2 -> float64
Dimension 3
SolAtVertices
<n>
<k> <type_1> ... <type_k>     # the fields of each row and their types
v v v ...                     # n rows, the k fields concatenated
Time                          # optional, read and ignored (as is Iterations)
0.5
End
```

The binary file is the libMeshb record stream of [`.meshb`](./medit.md#binary-meshb-libmeshb-gmf): the magic code (byte-swapped for the other endianness), the version that sets the widths, `GmfDimension`, then records `(code, next position, row count, field count, field types, values)`. The field count and types are always 4-byte integers. The codes are `GmfSolAtVertices` (62), `GmfSolAtEdges` (63), `GmfSolAtTriangles` (64), `GmfSolAtQuadrilaterals` (65), `GmfSolAtTetrahedra` (66), `GmfSolAtPrisms` (67), `GmfSolAtHexahedra` (68) and `GmfSolAtPyramids` (26). Other records (`GmfTime`, `GmfIterations`, the high-order `HOSolAt*`) are skipped through their next-section position, with a warning.

## Data mapping

| type | meaning | components in `d`-D | array shape |
|---|---|---|---|
| 1 | scalar | 1 | `(n,)` |
| 2 | vector | `d` | `(n, d)` |
| 3 | symmetric matrix | `d(d+1)/2` | `(n, 3)` in 2-D, `(n, 6)` in 3-D |
| 4 | full matrix | `d²` | `(n, 4)` in 2-D, `(n, 9)` in 3-D |

A symmetric matrix keeps the file's order, the lower triangle row by row: `m11 m21 m22` in 2-D and `m11 m21 m22 m31 m32 m33` in 3-D. Since the matrix is symmetric that is also `m11 m12 m22 m13 m23 m33`. This is not the `xx yy zz xy yz zx` order of a six-component VTK array ([data model](../mesh_data_model.md)), and readers keep the order their file uses, so convert it before handing a VTK tensor to MMG. A full matrix is row-major, the layout of the 9-component output of [`hessian`](../hessian.md).

On write, each floating-point point or cell array becomes the field whose type its component count names in the mesh's dimension. They are written in sorted-name order, vertex fields in one `SolAtVertices` section and element fields in one section per cell type. Integer arrays (a solution holds reals; `medit:ref` silently), arrays of another width and cell types with no `SolAt*` keyword are skipped with a warning. So is a type whose values for a field are all `NaN`, the fill another type's section left. Values are written as `float64` (`MeshVersionFormatted 2` in ASCII, version 3 in binary, version 4 for a section of more than 2³¹ rows). A mesh with nothing to write is refused.

## Quirks & limitations

- The dimension comes from the mesh's points (2-D or 3-D), else from `field_data["medit:dimension"]`; a geometry-less mesh with neither is refused, because the dimension decides what a type-2 or type-3 field holds.
- Attaching checks the solution against the mesh: the same dimension, as many `SolAtVertices` rows as points and as many rows of each element section as cells of that type, else `ReadError`.
- A `SolAt*` section inside a `.mesh` is skipped with a warning, as the binary reader already skipped one inside a `.meshb`; keep it in a `.sol`.
- Field names do not survive: a solution has no names, so a read numbers them.
