# Surface reconstruction

`reconstruct_surface(cloud)` turns an oriented point cloud into a triangle surface. A scan read from PLY, PCD or XYZ has positions, often normals, and no cells; meshing it, measuring it against its design, or handing it to a solver needs a surface through it. The normals come from the cloud's `normals` array (the PCD and XYZ normal columns), or the `nx`, `ny` and `nz` scalars a PLY scan reads as, or are estimated first with [`estimate_normals`](./point_normals.md).

```python
import meshioplusplus as mp

scan = mp.read("part_scan.ply")
surface = mp.reconstruct_surface(scan)                          # normals read or estimated, k = 8
surface = mp.reconstruct_surface(points, normals=normals)      # arrays: (n, 3) points and normals
surface, report = mp.reconstruct_surface(scan, k=1, resolution=128, return_report=True)
report["num_components"], report["num_undefined"]
mp.write("part_scan.stl", surface)
```

## The field

The surface is the zero level set of the signed distance to the samples' tangent planes (Hoppe, DeRose, Duchamp, McDonald and Stuetzle, "Surface reconstruction from unorganized points", SIGGRAPH 1992). At a lattice node `x` the field blends the distances `(x − pᵢ) · nᵢ` to the planes of its `k` nearest samples, weighted by `exp(−|x − pᵢ|² / s²)`, where `s` is the larger of the sample spacing and the cell size. With `k=1` this is Hoppe's field exactly: piecewise planar, so the surface shows facets the size of the sample spacing. The default `k=8` smooths them away and averages out noise along the normals. The nearest samples come from an exact search (ties to the lower index), so the field is the same on every backend and thread count.

The samples are the points with a finite, non-zero normal; the rest are skipped. Normals need not be unit length, but they must be oriented consistently: a sample whose normal points the wrong way pulls the surface through itself. `estimate_normals`'s default orientation is consistent and outward for a closed scan.

## Trimming: open scans stay open

A node is **undefined** when its projection onto the nearest sample's tangent plane, `x − ((x − p₀) · n₀) n₀`, lies farther than `max_distance` from that sample: the node is beside the data rather than over it (Hoppe's test). Cells with an undefined corner are not cut. A one-sided scan therefore ends near its rim instead of being closed over the side nobody scanned, and a hole in the sampling wider than about twice `max_distance` stays a hole. The default is twice the larger of the spacing and the cell size; raising it bridges sparser sampling and lets the surface run further past a rim.

## The lattice

The cell size is `cell_size`, or the cloud's longest extent over `resolution`, or by default the **sample spacing**, the median distance from a sample to its fourth-nearest other sample: about the nearest-neighbour distance on a regular sampling, and the size of the gaps on a random one, whose nearest neighbours clump. Only the cells within reach of a sample are evaluated (a band of about two spacings either side of the samples), and only those the surface crosses are cut, by the marching-tetrahedra cutter [`isosurface`](./isosurface.md) uses, after each cell is split into six tetrahedra along its main diagonal. The cost therefore follows the scan's area, not its bounding box's volume. Crossings shared by neighbouring cells are one point, quads are split along their first diagonal and crossings that coincide are welded, so the result is one `triangle` block, watertight wherever the field is defined and wound the way the normals point: outward for outward normals.

`max_cells` caps the band (20 million cells by default); a cell size small enough to exceed it raises, naming the size used.

## Results

On Fibonacci samples of a unit sphere with exact normals, at the default cell size (the spacing), the result is one closed piece, wound outward, whose [Hausdorff distance](./hausdorff.md) to a finely tessellated sphere is a few hundredths of the spacing:

| Samples | Spacing | `k` | Triangles | Hausdorff / spacing | Time (4 cores) |
| --- | --- | --- | --- | --- | --- |
| 4 000 | 0.058 | 1 | 33 520 | 0.019 | 0.05 s |
| 4 000 | 0.058 | 8 | 33 608 | 0.026 | 0.06 s |
| 100 000 | 0.0123 | 1 | 746 956 | 0.022 | 1.7 s |
| 100 000 | 0.0123 | 8 | 747 072 | 0.025 | 1.8 s |
| 1 000 000 | 0.0038 | 8 | 7 988 616 | | 24 s |

Random samples clump, which is why the spacing is measured to the fourth-nearest sample: 20 000 uniformly random points on the sphere give one closed piece within 0.04 of a spacing, as do random scans of lumpy bodies with estimated normals. With 100 000 Fibonacci samples displaced by Gaussian noise of 0.002 (a sixth of the spacing) and no normals, the normals are estimated and the surface lies within 0.0044 of the sphere, in 1.9 s. The default cell size gives about eight triangles per sample; a `resolution` or a larger `cell_size` gives fewer. The field and the band are evaluated in parallel and the cut runs serially; the output is the same at every thread count (50 000 noisy samples: 2.9 s on one thread, 1.1 s on four, identical bytes).

The output has no data arrays; it keeps the input's [length unit](/formats#length-units). There is no numpy twin (the `estimate_normals` precedent): the nearest-sample ties and the cutter's serial numbering would have to match bit for bit. The Python function raises `NotImplementedError` without the compiled core.

What it does not do: Poisson reconstruction (a global solve, smoother on noisy scans and always watertight), adaptive (octree) cells, or sharp features, which the blended planes round off at the scale of `s`.

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `reconstruct_surface(points, *, normals="normals", k=8, cell_size=0.0, resolution=0, max_distance=0.0, estimate_normals=True, normals_k=16, max_cells=20_000_000, return_report=False)`; `normals` may also be an `(n, 3)` array |
| Pipeline | `{"Op": "ReconstructSurface", "Normals": "normals", "K": 8, "CellSize": 0, "Resolution": 0, "MaxDistance": 0, "EstimateNormals": true, "NormalsK": 16, "MaxCells": 20000000}`; counters `NumSamples`, `Spacing`, `CellSize`, `MaxDistance`, `NumLatticePoints`, `NumUndefined`, `NumTriangles`, `NumComponents`, `EstimatedNormals` |
| Both CLIs | `meshioplusplus reconstruct IN OUT [--normals NAME] [--k N] [--cell-size H] [--resolution N] [--max-distance D] [--no-estimate] [--normals-k N] [--max-cells N] [--json]` |
| C | `mio_reconstruct_surface_opts_init(&o); o.k = 1; mio_mesh* out = mio_reconstruct_surface(mesh, &o, &report);` |
| Fortran | `out = m%reconstruct_surface(k=1, resolution=128, num_triangles=nt, num_components=nc, stat=ierr)` |
| Julia | `reconstruct_surface(mesh; normals="normals", k=8, cell_size=0.0, resolution=0, max_distance=0.0, estimate_normals=true, normals_k=16, max_cells=20_000_000)` |
| R | `mio_reconstruct_surface(mesh, normals = "normals", k = 8, cell_size = 0, resolution = 0, max_distance = 0, estimate_normals = TRUE, normals_k = 16, max_cells = 2e7)` |
| WASM | `m.reconstructSurface(mesh, {normals, k, cellSize, resolution, maxDistance, estimateNormals, normalsK, maxCells})` returns `{mesh, numSamples, spacing, cellSize, maxDistance, numLatticePoints, numUndefined, numTriangles, numComponents, estimatedNormals}` |
| MCP | tool `reconstruct_surface(input_path, output_path, normals, k, cell_size, resolution, max_distance, estimate_normals, normals_k, max_cells)` |
