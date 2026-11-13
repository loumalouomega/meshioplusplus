# UGRID NetCDF — unstructured grids in netCDF (`.nc`)

The [UGRID conventions](https://ugrid-conventions.github.io/ugrid-conventions/) for unstructured meshes in netCDF: what Delft3D-FM, D-Flow FM, SCHISM, FVCOM and xugrid write, and QGIS's mesh layer reads (roadmap 2.4.3.2). Both engines read and write it (`formats/ugrid_nc.cpp` on the shared [netCDF helpers](../common_utilities.md); `meshioplusplus.ugrid_nc` through `netCDF4`). The AFLR3 [`ugrid`](./ugrid.md) format keeps its own name.

| | |
|---|---|
| **Format name** | `ugrid_nc` |
| **Extensions** | `.nc` |
| **Read / Write** | ✓ / ✓ (1-D and 2-D topologies) |
| **Extra dependencies** | netCDF (`MESHIOPLUSPLUS_WITH_NETCDF`) in the core; `netCDF4` for the Python reference |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("model_map.nc")                       # the first step
last = meshioplusplus.ugrid_nc.read("model_map.nc", time_step=-1)
times = meshioplusplus.ugrid_nc.time_values("model_map.nc")
```

- The first variable with `cf_role = "mesh_topology"` is the mesh; a volume (`topology_dimension = 3`) topology is refused by name, and an Exodus file (`num_nodes`, `num_dim`) is told to read as `exodus`.
- `node_coordinates` gives x, y and an optional third coordinate (z is 0 otherwise). Spherical `lon`/`lat` coordinates are kept as they are.
- `face_node_connectivity`: `start_index` (0 or 1) is honoured, the `_FillValue` padding of mixed polygons is dropped, and an array stored `(max_nodes, faces)` is read through the topology's `face_dimension`. Faces become a `triangle`, a `quad` and one `polygon` block per larger node count, in that order, each in the file's face order.
- `edge_node_connectivity` becomes a trailing `line` block.
- Data variables naming the mesh (`mesh`) with `location` `node`, `face` or `edge` become point and cell data (a trailing dimension becomes components). A face array is NaN (0 for an integer array) on the line block, and an edge array on the face blocks; a variable carrying `meshioplusplus_name` is read under that name, which is how a written `<name>_edge` merges back.
- A leading time dimension is the step axis: `time_step` (0 the first, negatives from the end; `ReadOptions::mTimeStep` in C++) picks one, its value is `field_data["meshio:time"]`, and `time_values` (the metadata reader in the core) lists them all. Variables laid out otherwise are skipped with a warning.

## Writing

```python
meshioplusplus.write("out.nc", mesh)
```

A netCDF-4 file with one topology `mesh` (2-D, or 1-D when the mesh has only lines): `mesh_node_x`, `mesh_node_y` (and `mesh_node_z` when any z is non-zero), `mesh_face_nodes` from the triangle, quad and polygon blocks (`_FillValue` -1, `start_index` 0) and `mesh_edge_nodes` from the line blocks. Point data are `node` variables, cell data `face` variables, and the line-block part of a cell data array an `edge` variable `<name>_edge` carrying `meshioplusplus_name`. Names are made netCDF-safe (the original kept in `meshioplusplus_name`). Volume and other cells are refused by name. One step is written; a series writer is not there yet.

## Verification

xugrid (0.15) opens a file meshio++ wrote, and meshio++ reads the file xugrid writes back (`tests/python/test_ugrid_nc.py`); the mixed, transposed, 1-based and time-dependent layouts are built in the tests with netCDF4. No Delft3D-FM `_map.nc` has been read yet.
