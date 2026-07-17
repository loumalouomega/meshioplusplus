# HMF (`.hmf`)

An experimental HDF5 mesh container specific to meshio++, reusing the XDMF topology-name vocabulary. **The format may change at any time** — writing always emits a warning to that effect.

| | |
|---|---|
| **Format name** | `hmf` |
| **Extensions** | `.hmf` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` (or a C++ build with HDF5) |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.hmf")
meshioplusplus.hmf.write("out.hmf", mesh, compression="gzip", compression_opts=4)
```

- **`compression`** / **`compression_opts`** — HDF5 gzip filter and level.

## File structure

```
(file attrs) type="hmf", version="0.1-alpha"
domain/
  grid/
    Geometry              # dataset = mesh.points; attr GeometryType = "X"|"XY"|"XYZ"
    Topology{k}            # one dataset per cell block; attr TopologyType = XDMF type name
    NodeAttributes/<name>   # one dataset per point-data key
    CellAttributes/<name>   # one dataset per cell-data key (blocks concatenated)
```

Only one `domain` and one `grid` are supported. The design rationale (from the source): HDF5 doesn't allow multiple same-named "Attribute" entries the way XDMF's XML does, so HMF uses separate `NodeAttributes`/`CellAttributes` subgroups keyed by name, and separate numbered `Topology{k}` datasets rather than one grouped-by-type "Topology".

## Cell types

Reuses [XDMF](./xdmf.md)'s `meshio_to_xdmf_type`/`xdmf_to_meshio_type` tables for the `TopologyType` attribute — see the XDMF page for the full mapping.

## Data mapping

Fully generic: any `point_data`/`cell_data` key name is preserved verbatim as an HDF5 dataset name under `NodeAttributes`/`CellAttributes`.

## Quirks & limitations

- If two `Topology{k}` datasets happen to resolve to the **same** meshio++ type, the Python reader's dict-based accumulation means the **later** one silently **replaces** the earlier one (since `cells` is keyed by meshio++ type name, not by dataset index). The C++ reader deliberately replicates this exact "later entry wins" semantics rather than merging or erroring.
- `GeometryType` is asserted to be one of `"X"`/`"XY"`/`"XYZ"` but is otherwise unused after validation.
- The C++ reader correctly round-trips **multi-block** cell data (several cell blocks with data under the same `CellAttributes` name); the reference Python/`h5py` reader has a known issue handling this case correctly for the same input, making the C++ path strictly more correct here — one of the few formats where this is true.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_HDF5`, otherwise through `h5py`.
- No reference fixture exists under `tests/meshes/hmf/`; tests round-trip synthetic meshes only.
