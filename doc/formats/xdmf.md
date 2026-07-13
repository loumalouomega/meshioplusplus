# XDMF (`.xdmf`, `.xmf`)

The [XDMF](https://xdmf.org/index.php/XDMF_Model_and_Format) format: an XML
"light data" file describing topology, geometry and attributes, whose "heavy
data" lives inline in the XML, in external raw binary files, or in a companion
HDF5 file.

| | |
|---|---|
| **Format name** | `xdmf` |
| **Extensions** | `.xdmf`, `.xmf` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` (for `data_format="HDF"`) |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.xdmf")
meshio.xdmf.write("out.xdmf", mesh,
    data_format="HDF",   # "HDF", "XML", or "Binary"
    compression="gzip",  # HDF only
    compression_opts=4,
)
```

- **`data_format`** — where the arrays are stored: `"XML"` (inline text),
  `"Binary"` (external `.bin` files), or `"HDF"` (a companion `.h5`).
- **`compression`** / **`compression_opts`** — gzip filter for HDF data.

## File structure

`<Xdmf><Domain><Grid>` with `<Geometry>` (XY/XYZ), `<Topology>` (single type or
`Mixed`), and `<Attribute>` elements (`Center="Node"`/`"Cell"`). Each holds a
`<DataItem>` whose `Format` is `XML`, `Binary` or `HDF`.

## Data mapping

- Node/cell `<Attribute>` → `point_data` / `cell_data`.
- Mixed topology is expanded into per-type cell blocks.

## Time series

Temporal XDMF is written/read with the `TimeSeriesWriter` / `TimeSeriesReader`
classes — see [XDMF time series](../xdmf_time_series.md).

## Notes

- The C++ core handles the **XML** and **Binary** data paths always, and the
  **HDF** path when built with `MESHIO_WITH_HDF5` (otherwise `h5py`). Single and
  Mixed topology, XY/XYZ geometry and node/cell attributes are supported. The
  time-series classes remain Python.
