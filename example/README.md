# Examples

Jupyter notebooks demonstrating meshio++ on the bundled `example.msh` (a Gmsh
4.1 mesh of a mechanical bracket, ~52k nodes / ~298k elements), in two
languages against the same geometry:

- [`python/`](python/) -- the Python bindings + PyVista.
- [`cpp/`](cpp/) -- the **C++ core** directly, no Python, on the
  [xeus-cpp](https://github.com/compiler-research/xeus-cpp) Jupyter kernel;
  see [`cpp/README.md`](cpp/README.md) for kernel setup.

Each has its own three-notebook tour (read/visualize, convert/inspect,
operations); see the READMEs in each subfolder for the per-notebook
breakdown and how to re-run them.

## Shared files (this folder)

`example.msh` is the bracket geometry both languages' notebooks read.
`example.stp` / `example.stp.geo` / `example.stp.mesh.json` are the CAD source
and meshing options it was generated from.

`Bunny.stl` is ["Stanford Bunny -- Digitized!"](https://www.thingiverse.com/thing:88208)
by MakerBot, licensed under the
[Creative Commons - Attribution](https://creativecommons.org/licenses/by/3.0/)
license (~112k triangles, binary STL). It feeds the project logo
(`doc/logo/gen_logo_tikz.py` renders it through meshio++'s 3D TikZ path) and makes
a handy real-world surface mesh for trying the STL/PLY/SVG/TikZ writers.
