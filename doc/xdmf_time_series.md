# XDMF Time Series

XDMF is the only format in meshio with built-in support for temporal (time series) data. The mesh topology is written once; field data is written per time step.

Requires `h5py` when using the default `data_format="HDF"`.

---

## Writing a time series

```python
import meshio

with meshio.xdmf.TimeSeriesWriter("simulation.xdmf") as writer:
    writer.write_points_cells(points, cells)
    for t, phi in time_steps:
        writer.write_data(t, point_data={"phi": phi})
```

### `TimeSeriesWriter(filename, data_format="HDF")`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `filename` | — | Path to the `.xdmf` file |
| `data_format` | `"HDF"` | `"HDF"` (companion `.h5`), `"XML"` (inline), or `"Binary"` (separate `.bin` files) |

Must be used as a context manager (`with` statement). The `.xdmf` file is written on `__exit__`.

### `writer.write_points_cells(points, cells)`

Write the shared mesh topology. Must be called before `write_data`.

### `writer.write_data(t, point_data=None, cell_data=None)`

Write field data for one time step `t` (a float). Both `point_data` and `cell_data` are dicts of `str -> numpy array`.

---

## Reading a time series

```python
with meshio.xdmf.TimeSeriesReader("simulation.xdmf") as reader:
    points, cells = reader.read_points_cells()
    for k in range(reader.num_steps):
        t, point_data, cell_data = reader.read_data(k)
```

### `TimeSeriesReader(filename)`

Parses the XDMF file on construction. Only XDMF version 3 is supported for time series.

### `reader.num_steps`

Total number of time steps stored in the file.

### `reader.read_points_cells()`

Returns `(points, cells)` — the shared mesh topology as a numpy array and a list of `CellBlock`.

### `reader.read_data(k)`

Returns `(t, point_data, cell_data)` for time step index `k`.

---

## Notes

- The mesh topology is stored once in the XDMF file and referenced by each time step using XInclude.
- With `data_format="HDF"`, all numerical data goes into a companion `<stem>.h5` file. Both files must be present to read.
- `data_format="XML"` embeds all data directly into the XML, which avoids external files but produces large `.xdmf` files.
- `data_format="Binary"` writes one `.bin` file per data array; useful when HDF5 is not available.
