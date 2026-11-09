# Removing islands

`remove_islands(mesh)` drops a scan's crumbs: the connected components whose bounding-box diagonal, area or cell count is below a threshold (roadmap 7.8.1). A scanner leaves them where it caught dust, a fixture's edge or a reflection; they spoil a fit, a deviation map and a watertight check.

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

clean, report = meshioplusplus.remove_islands(scan, return_report=True)
report["num_removed"], report["num_cells_removed"]

# What would go, without removing it: islands:component and islands:removed.
labelled = meshioplusplus.remove_islands(scan, dry_run=True)
```

| Option | Default | Meaning |
|---|---|---|
| `connectivity` | `"edge"` | `"edge"`: cells sharing a facet (an edge of a surface, a face of a volume) are one component; `"point"`: cells sharing a point are |
| `min_diagonal` | `0.01` | a component whose bounding box's diagonal is below this goes |
| `min_area` | `0` (off) | a component whose area (of its triangles, quads and polygons) is below this goes |
| `min_cells` | `0` (off) | a component with fewer cells goes |
| `relative` | `true` | `min_diagonal` and `min_area` are fractions of the whole mesh's diagonal and area |
| `dry_run` | `false` | remove nothing; return the mesh with the Int64 cell data `islands:component` (its component) and `islands:removed` (1 where it would go) |

A component goes when any threshold says so. Components are numbered in the order of each one's lowest cell (block-major), the edge-connected ones by the facet-sharing dual graph (`detail/cell_components.hpp`'s `cell_components_by_edge`, beside the point-connected `cell_components` that `split(by="component")` uses). The kept cells keep their data and regions, and the points no kept cell uses are dropped. The report gives `num_components`, `num_removed`, `num_cells_removed`, `num_points_removed`, and per component `component_cells`, `component_diagonal`, `component_area` and `component_removed`.

A sphere plus ten crumbs a hundredth of its size keeps the sphere only (`tests/cpp/test_remove_islands.cpp`, `tests/python/test_scan_edit.py`).

## Surfaces

C++ `remove_islands` (`operations/remove_islands.hpp`), Python `meshioplusplus.remove_islands` (C++ core only), both CLIs (`remove-islands IN OUT [--connectivity edge|point] [--min-diagonal F] [--min-area F] [--min-cells N] [--absolute] [--dry-run] [--json]`), the `RemoveIslands` [pipeline](./pipeline.md) step, C `mio_remove_islands`, Fortran, Julia, R, Octave/MATLAB, WASM `removeIslands` and the MCP tool `remove_islands`.
