# Mass properties

`meshioplusplus.mass_properties(mesh, density=...)` integrates a density over a body and reports its **mass**, **centre of mass**, **inertia tensor**, **principal moments and axes** and **radii of gyration**, for the whole mesh and for every named `Cell` region. [`compute_stats`](/stats) reports the mean of the points as its `centroid`, and [`measure`](/measure) is GD&T; neither is the centre of mass of a body.

```python
import meshioplusplus as mio

mesh = mio.read("part.stl")
r = mio.mass_properties(mesh, density=7850.0)
print(r["total"]["mass"], r["total"]["centre"], r["total"]["principal_moments"])
```

```sh
meshioplusplus mass part.stl --density 7850 --json
```

## What is integrated

| Cells | Body | Notes |
| --- | --- | --- |
| Volume cells (every tabulated 3D type with a face table, polyhedra) | the cells | each cell is the corner-average fan of the [polyhedron kernel](/polyhedra), the decomposition `poly_measure` and OpenFOAM use, integrated exactly per tetrahedron for a constant density |
| Triangles and polygons forming a **closed, consistently wound** shell | the body they enclose | the divergence theorem, exact for a polyhedron; an inward-wound shell is turned outward |
| Surface cells with `thickness` above zero | a **shell** | the mass is `density * thickness * area`; open surfaces are allowed |

An open or inconsistently wound surface without `thickness` throws by name, pointing at [`repair`](/repair). A curved (higher-order) cell is measured on its corners, as everywhere else in the library, and a cell with no face table is skipped with a warning. Cells of a lower dimension than the highest present are ignored, with a warning.

## Density

`density` is the density of every cell. `density_array` names a scalar `cell_data` array that replaces it per cell (volume cells and shells; a closed surface has no cells of its own to carry a density, so it is refused there). A non-finite density throws.

## The report

`mass_properties` returns `{"kind", "total", "regions", "warnings"}`. `kind` is `"solid"` or `"shell"`. Each entry (`total`, and one per named cell region in sorted name order) holds:

| Key | Meaning |
| --- | --- |
| `name` | the region's name; empty for the total |
| `num_cells` | cells that took part |
| `measure` | volume, or area for a shell |
| `mass` | `integral(density)` |
| `centre` | centre of mass |
| `inertia` | `xx yy zz xy yz zx`, with `Ixy = -integral(rho x y)` |
| `principal_moments` | ascending |
| `principal_axes` | a row-major 3x3, the axis of moment *k* in column *k* (the convention of the [tensor invariants](/tensor_invariants)) |
| `gyration_radii` | `sqrt(principal moment / mass)` |

`about="centroid"` (default) takes the inertia tensor about the centre of mass; `about="origin"` takes it about the coordinate origin (the parallel-axis theorem relates the two). `region=` restricts the total to one cell region and lists no sub-regions.

## Determinism

The mesh is never modified. Per-cell contributions are computed in parallel and summed serially in cell order about the mean of the points (so a body far from the origin keeps its accuracy); the report is identical on every backend and thread count. There is no numpy fallback: the function raises `NotImplementedError` without the compiled core.

## Surfaces

C++ `mass_properties` (`operations/mass_properties.hpp`), Python, both CLIs (`meshioplusplus mass`) and the [MCP server](/mcp) tool `mass_properties`.
