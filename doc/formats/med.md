# MED / Salome (`.med`)

The [MED](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/med-file.html)
format (Salome/Code-Aster), stored in HDF5.

| | |
|---|---|
| **Format name** | `med` |
| **Extensions** | `.med` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` (or a C++ build with HDF5) |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.med")
meshio.med.write("out.med", mesh)
```

`write` takes no keyword arguments (MED does not compress).

## File structure

`ENS_MAA/<mesh>` with a time-step group holding `NOE/COO` (Fortran-order node
coordinates) and `MAI/<MED type>/NOD` (1-based, Fortran-order connectivity);
`FAS` holds named families (point/cell tags); `CHA` holds fields (`NOEU` nodal,
`ELEM` per-cell, `ELNO` per-node-per-cell).

## Cell types

The MED ↔ meshio type table including second-order elements.

## Data mapping

- Point/cell family tags → `point_data["point_tags"]` / `cell_data["cell_tags"]`
  with named families exposed as `mesh.point_tags` / `mesh.cell_tags`.
- Field names → `field_data["med:nom"]`; fields → `point_data` / `cell_data`.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_HDF5`, otherwise
  through `h5py`. `ELGA` (Gauss-point) fields and non-default profiles fall back
  to the Python reader.
