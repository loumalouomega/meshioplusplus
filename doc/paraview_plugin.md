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

## Usage

After loading, every extension meshio++ knows appears in ParaView's file dialogs. The reader's **FileFormat** property picks the format explicitly (`automatic` infers it from the file), and the file is exposed as a `vtkUnstructuredGrid` with its point, cell and field data. The writer saves an unstructured grid through `meshioplusplus.write`, choosing the format from the file name.

## Limitations

- Polyhedron blocks are skipped by the reader (the legacy cell layout the plugin builds has no face stream); a message on ParaView's output names each skipped block.
- The writer groups cells by VTK type and node count, so a mesh's cell order is not preserved across a write.
- Optional format dependencies (`h5py`, `netCDF4`, ...) must be installed in ParaView's Python as well.

## Testing

`tests/python/test_paraview_plugin.py` runs `pvpython` in a subprocess: it loads the plugin, reads a file through the reader, writes it back through the writer, and reads the result with meshio++. It is skipped when `pvpython` is absent or runs a different Python than the build under test; CI runs it in its `paraview-plugin` job against conda-forge ParaView.
