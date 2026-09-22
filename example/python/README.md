# Examples (Python)

Jupyter notebooks demonstrating meshio++'s Python bindings on the bundled
[`../example.msh`](../example.msh) (a Gmsh 4.1 mesh of a mechanical bracket,
~52k nodes / ~298k elements).

| Notebook | What it shows |
|----------|---------------|
| [`01_read_and_visualize.ipynb`](01_read_and_visualize.ipynb) | Read the Gmsh file with meshio++, inspect the mesh, and render it with PyVista (full view + a clipped interior view). |
| [`02_convert_and_inspect.ipynb`](02_convert_and_inspect.ipynb) | Convert the geometry to VTU / VTK / XDMF / Gmsh / PLY, compare file sizes, and verify the round trip. |
| [`03_mesh_operations.ipynb`](03_mesh_operations.ipynb) | Tour of the operations layer -- surface/skin extraction, quality, reorder, diff, transform, clean, crop (bounding box / half-space / data predicate), merge, split, stats, convert_cells, refine, partition, smooth, interpolate, regular grids and signed distance fields (`grid`, `voxelize`, `compute_sdf` voxel/octree, `.vti`), slice, isosurface, gradient, the settings pipeline, sequences, the five data operations, and selective reads -- each rendered with PyVista. |
| [`04_interop.ipynb`](04_interop.ipynb) | Hand the mesh to PyVista, trimesh and Arrow/Parquet with no file round-trip -- shared buffers checked with `np.shares_memory`, `zero_copy_only`, regions through a PyVista round-trip, and a Parquet export read back with pandas. Needs `meshioplusplus[interop]`. |
| [`05_point_clouds.ipynb`](05_point_clouds.ipynb) | Read PCL's own `.pcd` files (ASCII, binary, `binary_compressed`, organised), render the colour-coded cloud, write every `DATA` mode and compare sizes, then `subsample_points` → `proximity_graph`, and the column rules and chemistry-XYZ refusal of `.xyz`. |
| [`06_lsdyna_decks.ipynb`](06_lsdyna_decks.ipynb) | Read a multi-file LS-DYNA keyword deck (long-format nodes, I10 solids, free-format shells, `*INCLUDE_PATH`), render it one colour per `*PART`, see the tetra/pyramid/wedge collapse with a positive-volume check, the four card formats reading to one mesh, `*SET_*` as point/cell/side regions, and a write → read round trip. |
| [`07_calculix_frd.ipynb`](07_calculix_frd.ipynb) | Read CalculiX `.frd` results written by `ccx`: the increments as sequence steps, the deformed cantilever coloured by von Mises (with a beam-theory check), the first four bending modes of the same beam, `read_sequence` → `{step}` files, shells and beams that `ccx` expands into solids, the he20/pe15/be3 node permutations checked geometrically, the short layout and the cgx shell types, and derived fields against NumPy. |
| [`08_gltf_export.ipynb`](08_gltf_export.ipynb) | Write a mesh as glTF 2.0 and look at what a viewer sees: an L-bracket with two cell regions and a boundary patch, the skin's point normals smooth versus split at a crease angle (`compute_normals`), a temperature field baked into an unlit `COLOR_0`, the `.glb` decoded by hand (root node rotation, scale, recentring; linear colours displayed as sRGB), the source coordinates recovered exactly, and a point cloud exported as `POINTS` with normals. |

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
