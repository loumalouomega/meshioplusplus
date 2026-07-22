# Examples (Python)

Jupyter notebooks demonstrating meshio++'s Python bindings on the bundled
[`../example.msh`](../example.msh) (a Gmsh 4.1 mesh of a mechanical bracket,
~52k nodes / ~298k elements).

| Notebook | What it shows |
|----------|---------------|
| [`01_read_and_visualize.ipynb`](01_read_and_visualize.ipynb) | Read the Gmsh file with meshio++, inspect the mesh, and render it with PyVista (full view + a clipped interior view). |
| [`02_convert_and_inspect.ipynb`](02_convert_and_inspect.ipynb) | Convert the geometry to VTU / VTK / XDMF / Gmsh / PLY, compare file sizes, and verify the round trip. |
| [`03_mesh_operations.ipynb`](03_mesh_operations.ipynb) | Tour of the operations layer -- surface/skin extraction, quality, reorder, diff, transform, clean, crop, merge, split, stats, convert_cells, refine, partition, smooth, interpolate, the five data operations, and selective reads -- each rendered with PyVista. |

See [`../cpp/`](../cpp/) for the same tour written directly against the **C++ core** (no Python), running on the [xeus-cpp](https://github.com/compiler-research/xeus-cpp) Jupyter kernel.

## Running them

The notebooks are committed **with their outputs** so they render on GitHub
without any setup. To re-run, install the notebook/rendering extras and execute
head-lessly:

```sh
uv pip install --python ../../.venv pyvista matplotlib jupyter nbconvert ipykernel
PYVISTA_OFF_SCREEN=true \
  ../../.venv/bin/jupyter nbconvert --to notebook --execute --inplace *.ipynb
```

PyVista renders off-screen through VTK's EGL backend (VTK ≥ 9.5), so no display
or `xvfb` is required; if GL is unavailable the notebook falls back to a
matplotlib surface plot.
