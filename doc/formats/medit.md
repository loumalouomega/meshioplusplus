# Medit / GMF (`.mesh`, `.meshb`)

The [Medit](https://people.sc.fsu.edu/~jburkardt/data/medit/medit.html) mesh format, also known as the INRIA `libMeshb` GMF (Gamma Mesh Format): keyword sections in ASCII (`.mesh`) or a binary, position-indexed record stream (`.meshb`). Dispatch is purely by whether the filename ends in `"b"`.

| | |
|---|---|
| **Format name** | `medit` |
| **Extensions** | `.mesh`, `.meshb` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.mesh")
meshioplusplus.medit.write("out.mesh", mesh, float_fmt=".16e")
```

- **`float_fmt`** (default `".16e"`) — coordinate format for the ASCII writer (the binary path's float width is dictated by the file's `MeshVersionFormatted` version, not by this kwarg).

## File structure

### ASCII (`.mesh`)

Whitespace/`#`-comment tokenized, keyword-driven:

```
MeshVersionFormatted <v>      # 0/1 -> float32 coords, 2 -> float64
Dimension <d>
Vertices
<n>
x1 x2 [x3] ref                # n rows
Triangles
<n>
v0 v1 v2 ref                  # n rows, 1-based
...
End
```

Recognized element keywords: `Edges`, `Triangles`, `Quadrilaterals`, `Tetrahedra`, `Prisms`, `Pyramids`, `Hexahedra` (and the alternate spelling `Hexaedra`). Every keyword listed above is followed by a count then that many rows of `<node ids> <ref>`. Several other keywords are recognized but discarded on read: `Corners`, `Normals`, `NormalAtVertices`, `SubDomainFromMesh`, `VertexOnGeometricVertex`/`Edge`, `EdgeOnGeometricEdge`, `Identifier`, `Geometry`, `RequiredVertices`, `TangentAtVertices`, `Tangents`, `Ridges`.

### Binary (`.meshb`, libMeshb GMF)

A leading `int32` magic code — `1` for native byte order, or its byte-swap (`16777216`) to signal the opposite endianness, which flips every subsequent typed read. Then an `int32` version (1-4), which fixes the integer/float/ position field widths: v1 → 4-byte int, 4-byte float; v2 → 4-byte int, 8-byte float; v3 → 4-byte int, 8-byte float, 8-byte positions; v4 → 8-byte int, 8-byte float, 8-byte positions. The rest of the file is a sequence of `(keyword_code, [position], [count], payload)` records, each keyed by a numeric GMF field code (from a large internal table of ~200 known libMeshb field codes — only a handful are meaningful to meshio++: `GmfVertices` (4), `GmfEdges` (5), `GmfTriangles` (6), `GmfQuadrilaterals` (7), `GmfTetrahedra` (8), `GmfPrisms` (9), `GmfPyramids` (49), plus `GmfMeshVersionFormatted` (1), `GmfDimension` (3) and `GmfEnd` (54)). Records with any other code are skipped with a warning.

## Cell types

| keyword | meshio++ type | nodes |
|---|---|---|
| `Edges` | `line` | 2 |
| `Triangles` | `triangle` | 3 |
| `Quadrilaterals` | `quad` | 4 |
| `Tetrahedra` | `tetra` | 4 |
| `Prisms` | `wedge` | 6 |
| `Pyramids` | `pyramid` | 5 |
| `Hexahedra` / `Hexaedra` | `hexahedron` | 8 |

## Data mapping

- `point_data["medit:ref"]` — the per-vertex trailing reference integer.
- `cell_data["medit:ref"]` — the per-element trailing reference integer, one array per cell block.

## Quirks & limitations

- Only **one** integer point-data array and **one** integer cell-data array can be written (Medit's single-`ref`-column limitation): if more than one candidate exists, the first is used and the rest silently dropped (a warning is emitted for the dropped ones).
- Coordinate dtype is version-driven, not user-selectable in binary mode: ASCII `MeshVersionFormatted` 0 or 1 → `float32`, 2 → `float64`; binary version 1 → `float32`, versions 2-4 → `float64`. The writer auto-upgrades to binary version 4 if any cell block's connectivity needs 8-byte integers.
- Binary record headers embed **absolute byte offsets** (`pos` fields) that the writer must track precisely while emitting records — a strict requirement for libMeshb compatibility, not merely a convenience field.
- `Corners`/`Normals`/etc. sections are read-but-discarded in the ASCII path; no equivalent handling exists in the C++ reader beyond simple token-skipping.
- The C++ core implements **only the ASCII `.mesh` variant**. The binary `.meshb` GMF format (including its little/big-endian variants) always falls back to the Python implementation.

## Notes

- `tests/meshes/medit/cube86.mesh` (ascii) — 39 points, 72 triangles, 86 tets, boundary tag counts `{1:14, 2:14, 3:14, 4:8, 5:14, 6:8}`.
- `sphere_mixed.1.meshb` (binary) — 3270 points, 864 triangles, 3024 wedges, 9072 tets, tag counts `{1:432, 2:216, 3:216}`.
- `hch_strct.4.meshb` / `hch_strct.4.be.meshb` (binary, little/big-endian pair, exercising the endian-swap path) — 306 points, 12 triangles, 178 quads, 96 wedges, 144 hexahedra, tag counts `{1:15, 2:15, 3:160}`.
- All four originate from UGRID files converted with UGC (simcenter.msstate.edu).
