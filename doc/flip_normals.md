# Flipping normals

`flip_normals(mesh, region=None)` reverses the winding of every cell, or of the cells of a cell region, with no point moved (roadmap 7.8.2). [`repair`](./repair.md) orients a surface consistently and outward; this is the operation for the opposite orientation, or for one part of a surface flipped on purpose.

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

inward = meshioplusplus.flip_normals(scan)
patched, report = meshioplusplus.flip_normals(scan, "lid", return_report=True)
report["num_flipped"]
```

Each cell's nodes are permuted by the same table [`mirror`](./mirror.md) uses (`detail/cell_flip.hpp`'s `block_flip`): a triangle or a polygon reverses its ring after node 0, a higher-order cell keeps its mid-edge and mid-face nodes in their slots, a volume cell turns inside out. A `Side` region entry `(cell, facet)` moves with its facet, as in `mirror`. Points, the block order, the cell data and the other regions are untouched. Every permutation is an involution, so flipping twice gives the input back, bit for bit (`tests/cpp/test_flip_normals.cpp`). An unknown region is refused by name.

## Surfaces

C++ `flip_normals` (`operations/flip_normals.hpp`), Python `meshioplusplus.flip_normals` (C++ core only), both CLIs (`flip-normals IN OUT [--region NAME] [--json]`), the `FlipNormals` [pipeline](./pipeline.md) step (`Region`), C `mio_flip_normals(mesh, region, &num_flipped)`, Fortran, Julia, R, Octave/MATLAB, WASM `flipNormals` and the MCP tool `flip_normals`.
