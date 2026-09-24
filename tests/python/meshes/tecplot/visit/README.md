# Tecplot files of older binary versions

Binary Tecplot files from VisIt's test data ([`data/tecplot_test_data.tar.xz`](https://github.com/visit-dav/visit/raw/develop/data/tecplot_test_data.tar.xz), BSD-3-Clause, see `LICENSE`). They are Tecplot's own example datasets, written by the Tecplot versions of their magic numbers; the `.tec` files are the ASCII twins shipped beside them. They check the reader's handling of the `#!TDV71`–`111` layouts, which follows VisIt's `TecplotFile.C`, against files neither meshio++ nor TecIO wrote.

| File | Version | Content | Checked against |
| --- | --- | --- | --- |
| `simpscat.plt` | `#!TDV71` | big-endian, one point-packed IJ-ordered zone, two text records | `simpscat.tec` |
| `fluid.plt` | `#!TDV71` | big-endian, block-packed IJK-ordered zone | `jetflow.plt` (the same `RHO-V` and `E`, cell by cell) |
| `jetflow.plt` | `#!TDV75` | block-packed IJK-ordered zone | `fluid.plt` |
| `eddy.plt` | `#!TDV75` | block-packed FE tetrahedra, 1-based connectivity | `fetebk.tec` |
| `polarplot.plt` | `#!TDV106` | double-precision ordered zones (a line and a surface) | both engines agree |
| `fetebk.plt`, `fetetpt.plt` | `#!TDV108` | FE tetrahedra, block and point packing | `fetebk.tec` |
| `movie.plt` | `#!TDV108` | two point-packed IJ-ordered zones | `movie.tec` (comma-separated values) |
