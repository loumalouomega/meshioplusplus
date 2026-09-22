# PCD (`.pcd`)

The [Point Cloud Library format](https://pointclouds.org/documentation/tutorials/pcd_file_format.html) (v0.7): a short text header (`VERSION`, `FIELDS`, `SIZE`, `TYPE`, `COUNT`, `WIDTH`, `HEIGHT`, `VIEWPOINT`, `POINTS`, `DATA`) followed by the points as ASCII rows (`DATA ascii`), packed little-endian records (`DATA binary`) or a field-by-field (struct-of-arrays) block inside an LZF stream (`DATA binary_compressed`, two `uint32` sizes then the stream).

| | |
|---|---|
| **Format name** | `pcd` |
| **Extensions** | `.pcd` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("cloud.pcd")
meshioplusplus.pcd.write("out.pcd", mesh, data="binary")
```

- **`data`** — `"ascii"`, `"binary"` (the default) or `"binary_compressed"`; when omitted, `binary=True/False` picks between the first two.
- **`point_dtype`** — the precision of `x y z` (and of normals, curvature and intensity): `"float32"` (the default: PCL's typed `PointXYZ` loaders reject any other field type), `"float64"`, or `"keep"` to follow the mesh's own points (what the in-place `ascii`/`binary`/`compress`/`decompress` verbs use). A lossy narrowing is recorded in the provenance.
- **`drop_invalid`** (read) — drop the points whose x, y or z is not finite (an organised cloud's invalid returns) and, when any is dropped, the organisation with them.

## Mesh mapping

The points plus exactly one `vertex` block — what [`subsample_points`](../point_budgets.md) produces. Point data is keyed by field name, with three conventions: `normal_x/y/z` → `"normals"` (n, 3); PCL's `rgb`/`rgba`, a `uint32` `0x00RRGGBB`/`0xAARRGGBB` that lives in a `float32` slot (unpacked by bit-cast, never by value) → `"rgb"` (n, 3) / `"rgba"` (n, 4) `uint8`; any other field keeps its own name and dtype, with `COUNT > 1` giving (n, count). `_` padding fields are skipped. `x`/`y`/`z` keep the file's precision (all `F4` → float32, otherwise float64).

## File structure

```
# .PCD v0.7 - Point Cloud Data file format
VERSION 0.7
FIELDS x y z [...]
SIZE 4 4 4 [...]
TYPE F F F [...]
COUNT 1 1 1 [...]
WIDTH <n>
HEIGHT 1
VIEWPOINT <tx ty tz qw qx qy qz>
POINTS <n>
DATA ascii | binary | binary_compressed
<points>
```

## Quirks & limitations

- **Organised clouds** (`HEIGHT > 1`) are kept whole, NaN rows included, with `WIDTH`/`HEIGHT` in `field_data["pcd:width"]`/`["pcd:height"]` (restored on write when `WIDTH * HEIGHT` equals the point count).
- **`VIEWPOINT` is recorded, never applied**: a non-identity viewpoint goes to `field_data["pcd:viewpoint"]`.
- Only `vertex` cells are writable; anything else, all cell data and unrelated field data are dropped with a warning and a provenance note.
- Point-data names are sanitised to the header's word vocabulary (spaces become `_`, collisions get a numeric suffix).

## Notes

- `tests/python/meshes/pcd/` — PCL's own test corpus (BSD-3, see `LICENSE.PCL` there): `bun0.pcd` (ascii, normals + curvature), `colored_cloud.pcd` (binary, organised, uint `rgb`), `pcl_logo.pcd` (binary_compressed, float-slot `rgb`, non-identity viewpoint), `milk_color.pcd` (binary_compressed, `rgba`).
- The C++ core handles all three `DATA` modes. Its LZF codec is an independent implementation of the liblzf stream format (no liblzf code is copied), and the Python reference writes byte-identical files. The CLI verbs `ascii`, `binary`, `compress` (`binary_compressed`) and `decompress` rewrite a `.pcd` in place, and the MCP `convert` tool takes `mode` and `compression: lzf`. The flat bindings (C, Fortran, Julia, R, WASM) read every mode and write `ascii`/`binary`; `binary_compressed` writing is on the [roadmap](../roadmap.md).

## Web output

A PCD cloud is a `vertex` block (or no cells at all), which the [glTF writer](./gltf.md) turns into a `POINTS` primitive, keeping the `normals` point data as `NORMAL` and every other point array as a raw `_NAME` attribute: `meshioplusplus convert cloud.pcd cloud.glb` puts a scanned cloud on the web. Run [`compute_normals`](../normals.md) first when the cloud has none.
