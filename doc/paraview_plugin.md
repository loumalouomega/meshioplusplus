# ParaView Plugin

meshio++ ships a ParaView Python plugin (`tools/paraview-meshioplusplus-plugin.py`) with a reader and a writer, so ParaView can open any file meshio++ reads and save to any format it writes, without converting first.

## Installation

The plugin runs inside ParaView's own Python, so meshio++ must be installed for that interpreter, and its compiled core must match that interpreter's version.

1. Find the Python that your ParaView uses:

   ```sh
   pvpython -c "import sys; print(sys.version)"
   ```

2. Install meshio++ for that Python (a conda-forge ParaView: into the same environment; a system ParaView: into the system Python, or a virtual environment on `PYTHONPATH`):

   ```sh
   pip install "meshioplusplus[all]"
   ```

3. The wheel installs the plugin as data, at a path with no ParaView version in it:

   ```
   <prefix>/share/meshioplusplus/paraview/paraview-meshioplusplus-plugin.py
   ```

   where `<prefix>` is the environment's root (`python -c "import sys; print(sys.prefix)"`). A source checkout has the same file at `tools/paraview-meshioplusplus-plugin.py`.

4. Load it in ParaView, either:

   - **Tools → Manage Plugins → Load New**, browse to the file, and tick **Auto Load** to load it on every start; or
   - point `PV_PLUGIN_PATH` at the directory before starting ParaView, which loads every plugin in it:

     ```sh
     export PV_PLUGIN_PATH="$(python -c 'import sys; print(sys.prefix)')/share/meshioplusplus/paraview"
     ```

From `pvpython` or `pvbatch`, `paraview.simple.LoadPlugin(path, ns=globals())` loads it and defines `meshioreader` and `meshioWriter`.

## Proxy names

The reader and the writer are the server-manager proxies `MeshioPlusPlusReader` (group `sources`) and `MeshioPlusPlusWriter` (group `writers`): stable identifiers that state files and `servermanager.ProxyManager()` lookups use. `paraview.simple` names its functions after their labels, `meshio++ reader` and `meshio++ Writer`, so scripts keep calling `meshioreader(...)` and `meshioWriter(...)`. The names the plugin used before, `meshio++ reader` and `meshio++ Writer`, stay registered as aliases (labelled "(old name)", without a file-dialog entry), so a state file saved with them still loads.

The reader registers every extension meshio++ knows, `.vtu` and the other formats ParaView reads itself included. In the GUI, ParaView then asks which reader to open such a file with and can remember the answer as the default for that extension; `paraview.simple.OpenDataFile` keeps picking ParaView's own reader for those extensions (CI checks this for `.vtu`), so call `meshioreader(FileName=...)` to read one through meshio++.

## Usage

After loading, every extension meshio++ knows appears in ParaView's file dialogs. The reader's **FileFormat** property picks the format explicitly (`automatic` infers it from the file), and the file is exposed as a `vtkUnstructuredGrid` with its point, cell and field data. The writer saves an unstructured grid through `meshioplusplus.write`, choosing the format from the file name.

## What crosses into ParaView

The reader hands ParaView the mesh's points (padded to three columns when the mesh is 1-D or 2-D), its cells, and its point, cell and field data, in one vectorised pass with no per-cell Python loop:

- Every cell keeps meshio++'s node order, which is VTK's: a linear `wedge` (Gmsh's prism) is VTK's wedge as it stands, so every wedge ParaView shows has a positive volume and a wedge made or edited in ParaView is written the right way out (meshio++'s VTK writers used to permute it, which turned it inside out; see the [changelog](https://github.com/loumalouomega/meshioplusplus/blob/main/CHANGELOG.md)).
- A block of uniform n-gons (`polygon3`, `polygon5`, ..., the boundary blocks the OpenFOAM, Fluent and CGNS readers produce) is a `VTK_POLYGON` block.
- The writer keeps the grid's cell order: each run of consecutive cells of one type (and, for polygons, one size) becomes a cell block, so a mixed grid can come back as several blocks of one type rather than one block per type.
- Polyhedra cross in both directions. A `polyhedron<n>` block becomes `VTK_POLYHEDRON` cells, whose points are the distinct nodes of their faces, with the faces set through `SetPolyhedralCells` on VTK 9.4 and later or the face stream `SetCells` took before (the same `set_cells` serves `to_pyvista`); the writer reads either back, one `polyhedron<n>` block per run of equal point count.
- Arbitrary-order VTK Lagrange and Bézier cells keep their degrees and weights: a `HigherOrderDegrees` cell array and a `RationalWeights` point array (the names the [VTU format](./formats/vtu.md) uses) become the active attributes of those names, which is where VTK's Lagrange and Bézier cells read them, and come back as arrays of the same names.
- A point or cell array meshio++ cannot hold (a `vtkStringArray`, say) is skipped with a warning naming it; a string array of field data (the QA records of ParaView's Exodus reader) is kept as a string array.

## Architecture: the payload layer

The plugin file declares the ParaView proxies and nothing else. The conversions live in `meshioplusplus._paraview`, which is the [Blender add-on](./blender.md#architecture-the-pure-payload-layer)'s pattern: the mesh-to-VTK direction is the payload `to_pyvista` uses (`_interop._to_vtk_payload`), the VTK-to-mesh direction is the VTU reader's `_vtk_common.vtk_cells_from_data`, and the two VTK adapters over them (`grid_from_mesh`, `mesh_from_grid`) import VTK inside each function. They need VTK but not ParaView, so `tests/python/test_paraview_payload.py` runs them with the PyPI `vtk` wheel in the default test matrix.

## Limitations

- Optional format dependencies (`h5py`, `netCDF4`, ...) must be installed in ParaView's Python as well.

## Testing

`tests/python/test_paraview_payload.py` tests the conversions with the PyPI `vtk` wheel and no ParaView. `tests/python/test_paraview_plugin.py` runs `pvpython` in a subprocess: it loads the plugin, reads files through the reader (wedges must have positive volumes in ParaView's `CellSize`, `polygon5` blocks and 1-D meshes must keep their cells), writes them back through the writer, and reads the result with meshio++. It is skipped when `pvpython` is absent or runs a different Python than the build under test; CI runs it in its `paraview-plugin` job against conda-forge ParaView, with HDF5 and netCDF on so XDMF, Exodus and CGNS files go through the plugin too, together with `tests/python/test_pvd_paraview.py`, the ParaView-side checks of the `.pvd` and `.pvtu` files meshio++ writes.
