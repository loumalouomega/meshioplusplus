# XDMF (`.xdmf`, `.xmf`)

The [XDMF](https://xdmf.org/index.php/XDMF_Model_and_Format) format: an XML "light data" description of topology/geometry/attributes, whose "heavy data" lives inline in the XML, in external raw binary files, or in a companion HDF5 file. Both XDMF2 and XDMF3 variants exist in the wild; meshio++ reads both, but the C++ core only handles XDMF3.

| | |
|---|---|
| **Format name** | `xdmf` |
| **Extensions** | `.xdmf`, `.xmf` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` (for `data_format="HDF"`) |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.xdmf")
meshioplusplus.xdmf.write("out.xdmf", mesh,
    data_format="HDF",   # "HDF", "XML", or "Binary"
    compression="gzip",  # HDF only
    compression_opts=4,
)
```

- **`data_format`** — where DataItem payloads are stored: `"XML"` (inline text), `"Binary"` (external `.bin` files, one per array), or `"HDF"` (a companion `.h5` file).
- **`compression`** / **`compression_opts`** — gzip filter for HDF data.

## File structure

```xml
<Xdmf Version="2.x"|"3.x">
  <Domain>
    <Grid Name="Grid" [GridType="Uniform"]>
      <Topology Type|TopologyType="<xdmf type>|Mixed" NumberOfElements=".." [NodesPerElement=".."]>
        <DataItem DataType|NumberType="Int|UInt|Float" Precision="1|2|4|8"
                  Dimensions=".. .." Format="XML|Binary|HDF">...</DataItem>
      </Topology>
      <Geometry Type|GeometryType="X|XY|XYZ"><DataItem .../></Geometry>
      <Attribute Name=".." AttributeType="Scalar|Vector|Tensor|Tensor6|Matrix" Center="Node|Cell|Grid">
        <DataItem .../>
      </Attribute>
    </Grid>
  </Domain>
</Xdmf>
```

`Format="HDF"` DataItem text is `"<file.h5>:/path/to/dataset"` (the HDF5 file path is resolved relative to the `.xdmf` file, not the CWD). `Format="Binary"` text is a raw-binary sibling file path. `Format="XML"` is whitespace- separated inline numbers.

`Mixed` topology encodes a flat array of `(xdmf_type_index, node0, node1, ...)` tuples concatenated across all cells; the numeric type-index table is shared with the per-type `TopologyType` names (below). A `line` (`Polyline`) entry in a Mixed array carries an extra "number of points in this polyline" field immediately after its type index, which meshio++ requires to equal exactly `2`.

## Cell types

meshio++ ↔ XDMF `TopologyType` names (each accepts one or more spellings on read; the first is used on write):

| meshio++ | XDMF | meshio++ | XDMF |
|---|---|---|---|
| `vertex` | `Polyvertex` | `tetra` | `Tetrahedron` |
| `line` | `Polyline` | `tetra10` | `Tetrahedron_10` / `Tet_10` |
| `line3` | `Edge_3` | `wedge` | `Wedge` |
| `quad` | `Quadrilateral` | `wedge15` | `Wedge_15` |
| `quad8` | `Quadrilateral_8` / `Quad_8` | `wedge18` | `Wedge_18` |
| `quad9` | `Quadrilateral_9` / `Quad_9` | `hexahedron` | `Hexahedron` |
| `pyramid` | `Pyramid` | `hexahedron20` | `Hexahedron_20` / `Hex_20` |
| `pyramid13` | `Pyramid_13` | `hexahedron27` | `Hexahedron_27` / `Hex_27` |
| `triangle` | `Triangle` | `hexahedron64..1331` | `Hexahedron_{64..1331}` / `Hex_{64..1331}` (Python-only, high-order) |
| `triangle6` | `Triangle_6` / `Tri_6` | | |

Numeric Mixed-topology type indices (a subset shared with the per-type names): `0x1`=vertex, `0x2`=line, `0x4`=triangle, `0x5`=quad, `0x6`=tetra, `0x7`=pyramid, `0x8`=wedge, `0x9`=hexahedron, `0x22`=line3, `0x23`=quad9, `0x24`=triangle6, `0x25`=quad8, `0x26`=tetra10, `0x27`=pyramid13, `0x28`=wedge15, `0x29`=wedge18, `0x30`=hexahedron20, `0x31`=hexahedron24, `0x32`=hexahedron27 (and further codes through `0x40` for the higher-order hexahedra, Python-only).

## Data mapping

`Attribute Name="..."` maps generically to `point_data`/`cell_data`; XDMF2 also supports `field_data` via an `Information` element holding `[num_tag, dim]` per key (XDMF2-read-only; never emitted on write, since an earlier attempt hit XML `CDATA` serialization bugs and was abandoned).

## Quirks & limitations

- **XDMF2 vs XDMF3**: dispatched by the major version digit in the root `Version` attribute. **The C++ core only implements version 3** — any XDMF2 file (`Version="2.x"`) always falls back to Python.
- XDMF2 uses `TopologyType`/`GeometryType`; XDMF3 accepts either that or the shorter `Type`, but errors if both are given on the same element simultaneously.
- A `Reference="XML"` / `Reference="<xpath>"` attribute on a `DataItem` supports XInclude-like references to another `DataItem` elsewhere in the document (only absolute `/`-rooted XPaths are resolved) — not implemented in the C++ core at all; any file using it would throw inside the C++ path and transparently fall back to Python.
- The **only** supported node count for a Mixed-topology `line` (`Polyline`) entry is exactly 2 — anything else raises `ReadError`.
- **The C++ core's type table is a strict subset** of the Python one — it covers up through `hexahedron27` but omits `hexahedron64` through `hexahedron1331`; files using those higher-order types fall back to Python.
- `data_format="HDF"` is handled by the C++ core only when built with `MESHIO_WITH_HDF5` and `compression in (None, "gzip")`; otherwise Python handles it via `h5py`.
- Points are restricted to dimension ≤3 on write (`WriteError` otherwise).

## Time series

Temporal XDMF is written/read with the `TimeSeriesWriter`/`TimeSeriesReader` classes — see [XDMF time series](../xdmf_time_series.md). These remain pure Python (stateful, HDF5-backed) regardless of the C++ core's availability.

## Notes

- No reference fixture exists under `tests/meshes/xdmf/`; tests use synthetic meshes across all three `data_format` values.
