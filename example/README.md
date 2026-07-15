# Examples

Jupyter notebooks demonstrating meshio++ on the bundled `example.msh` (a Gmsh
4.1 mesh of a mechanical bracket, ~52k nodes / ~298k elements).

| Notebook | What it shows |
|----------|---------------|
| [`01_read_and_visualize.ipynb`](01_read_and_visualize.ipynb) | Read the Gmsh file with meshio++, inspect the mesh, and render it with PyVista (full view + a clipped interior view). |
| [`02_convert_and_inspect.ipynb`](02_convert_and_inspect.ipynb) | Convert the geometry to VTU / VTK / XDMF / Gmsh / PLY, compare file sizes, and verify the round trip. |

## Running them

The notebooks are committed **with their outputs** so they render on GitHub
without any setup. To re-run, install the notebook/rendering extras and execute
head-lessly:

```sh
uv pip install --python ../.venv pyvista matplotlib jupyter nbconvert ipykernel
PYVISTA_OFF_SCREEN=true \
  ../.venv/bin/jupyter nbconvert --to notebook --execute --inplace *.ipynb
```

PyVista renders off-screen through VTK's EGL backend (VTK ≥ 9.5), so no display
or `xvfb` is required; if GL is unavailable the notebook falls back to a
matplotlib surface plot.

## Other files

`example.stp` / `example.stp.geo` / `example.stp.mesh.json` are the CAD source
and meshing options `example.msh` was generated from.
