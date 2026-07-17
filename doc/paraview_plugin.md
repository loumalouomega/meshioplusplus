# ParaView Plugin

meshio++ ships a ParaView plugin (`tools/paraview-meshioplusplus-plugin.py`) that lets you open any meshio++-supported file directly in ParaView without converting it first.

## Installation

1. Find the Python version that your ParaView uses:

   ```sh
   pvpython --version
   ```

2. Install meshio++ for that Python:

   ```sh
   pip install meshioplusplus[all]
   ```

3. Open ParaView and navigate to **Tools → Manage Plugins → Load New**.

4. Browse to the plugin file. On Linux it is typically installed at:

   ```
   ~/.local/share/paraview-5.9/plugins/paraview-meshioplusplus-plugin.py
   ```

You can also point directly to the file in the meshio++ source tree at `tools/paraview-meshioplusplus-plugin.py`.

5. *(Optional)* Tick **Auto Load** so the plugin is active every time ParaView starts.

## Usage

After loading the plugin, any meshio++-supported file extension appears in the ParaView file open dialog. ParaView will use meshio++ to read the file and expose it as a `vtkUnstructuredGrid`.

The plugin provides both a reader and a writer. The writer allows exporting ParaView data to any format meshio++ can write.

## Notes

- The plugin requires the same Python environment that ParaView's `pvpython` uses. If they differ (e.g. system Python vs. Conda), the plugin will not find meshioplusplus.
- All optional format dependencies (`h5py`, `netCDF4`) must be installed in that same environment.
