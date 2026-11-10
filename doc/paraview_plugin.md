# ParaView Plugin

meshio++ ships a ParaView Python plugin (`tools/paraview-meshioplusplus-plugin.py`): a reader and a writer, so ParaView can open any file meshio++ reads and save to any format it writes without converting first, and filters that bring the meshing and analysis operations ParaView has no equivalent for, plus one that runs any [settings pipeline](./pipeline.md).

## Installation

The plugin runs inside ParaView's own Python, so meshio++ must be importable by that interpreter, and its compiled core must match that interpreter's version. There are two ways to get there.

### The self-contained bundle

Every release attaches one zip per platform and ParaView Python to its GitHub release: `meshioplusplus-paraview-<version>-<platform>-cp<XY>.zip`, for Linux x64, Windows x64 and macOS (x64 and arm64), and Python 3.10 to 3.13. It holds the plugin beside a `site-packages` folder with meshio++, the OpenCASCADE, Netgen and MMG add-ons the platform has, and `rich`, already unpacked (numpy is left out: ParaView ships it). Nothing is installed: the plugin puts that folder on `sys.path` when it finds it beside itself, so it works with Kitware's binaries, which have their own Python and no pip.

1. Find the Python of your ParaView: `pvpython -c "import sys; print(sys.version)"` (3.12 is `cp312`).
2. Download the matching zip and unzip it anywhere; keep the folder whole.
3. Load `paraview-meshioplusplus-plugin.py` from it (*Tools → Manage Plugins → Load New*, *Auto Load*), or point `PV_PLUGIN_PATH` at the folder.

`src/paraview/build_bundle.py` builds it (`--platform`, `--python-tag`, `--local-wheel`, `--addon-wheel`), as `src/blender/build_extension.py` builds the Blender extension; CI's `paraview-bundle` job loads a bundle in the newest Kitware Linux binary of each ParaView series it lists, with no meshio++ installed.

### Installing into ParaView's Python

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

The reader registers every extension meshio++ knows, `.vtu` and the other formats ParaView reads itself included. In the GUI, ParaView then asks which reader to open such a file with and can remember the answer as the default for that extension; `paraview.simple.OpenDataFile` may pick meshio++'s reader for those extensions too once the plugin is loaded, so call the reader you want by name: `meshioreader(FileName=...)` reads through meshio++, and ParaView's own reader proxy (for example `XMLUnstructuredGridReader`) reads through ParaView.

## Usage

After loading, every extension meshio++ knows appears in ParaView's file dialogs. The reader exposes a file as a `vtkUnstructuredGrid` with its point, cell and field data (or, with **RegionsAsBlocks**, a collection of its regions), and the writer saves any dataset through `meshioplusplus.write`. Their properties are below; the filters are under *Filters → meshio++*.

## Reader properties

| Property | What it does |
| --- | --- |
| **FileFormat** | The format, or `automatic` to infer it from the file name. |
| **PointArrays**, **CellArrays** | Checkboxes filled from [`read_metadata`](./selective_read.md); an unticked array is not read (`read(arrays=...)`, which formats with a selective path never materialise). Only a changed selection is passed on. |
| **PointsOnly** | Geometry only: points and cells, no data arrays (`read(points_only=True)`). |
| **Piece** | One piece of a partitioned file (`pvtu`, `pvtp`, `pvd`, VTKHDF, MFEM, Elmer); `-1`, the default, reads every piece. |
| **Ghosts** | `keep` or `drop` the ghost cells of a partitioned file ([selective reads](./selective_read.md)). |
| **Deflection**, **AngularDeflection** | The tessellation of a [STEP](./formats/step.md#reading) or IGES file: a length (0 takes 1e-3 of the bounding-box diagonal) and an angle in radians (0.5). |
| **RegionsAsBlocks** | Output one block per named region instead of one grid (see [named regions](#named-regions)). |
| **AddOns** | Information only: whether the OpenCASCADE, Netgen and MMG add-ons are loaded (`occt: yes`, ...), so a STEP read that fails for want of OpenCASCADE explains itself. |
| **Refresh** | Plans the time steps again (see below). |

## Writer properties

The writer takes any dataset ParaView makes, not only unstructured grids: poly data, image data and structured grids become unstructured through `vtkAppendFilter`, and a composite dataset (a multiblock, a partitioned collection) has its leaves merged with [`merge`](./merge.md), as the [Blender add-on](./blender.md)'s export merges objects, each leaf's cells a cell region named after its block (whole cell blocks of their own, so a format such as Gmsh 4.1 keeps one physical group per block).

| Property | What it does |
| --- | --- |
| **FileFormat** | The format, or `automatic` to infer it from the file name. |
| **Encoding** | `default`, `ascii`, `binary` or `raw_appended`. |
| **Codec** | `none`, `zlib`, `lz4`, `zstd` (the VTK XML formats' block codecs) or `lzf` (PCD). |
| **FloatFormat** | A printf-style float format for the text formats, e.g. `.6e`. |
| **Provenance** | On by default: the file's [provenance](./provenance.md) block names ParaView and its version as the source, with the target format and the options above. |
| **WriteTimeSteps** | Every time step rather than the current one (see below). |

Encoding, codec and float format take the [settings pipeline](./pipeline.md)'s `Output` mapping, so each format's spelling has one owner and a combination the format cannot take (an ASCII file with a codec, say) is refused by name.

## Time steps

A file that holds several steps (every format [`read_metadata`](./sequences.md) reports `time_values` for: XDMF, Exodus, VTKHDF, CalculiX `frd`, ...) shows them on ParaView's time toolbar, and each time reads its step: the last step at or before the time asked for. Selecting a file series in the file dialog (`out_0.vtu`, `out_1.vtu`, ...) opens it as one reader with a step per file, timed as [sequences](./sequences.md) time a series: a number in the file name, else the file's own `meshio:time`, else its position. The plan of steps is cached per file and its size and modification time, so *File → Reload Files* on a rewritten file plans it again, and the reader's **Refresh** button drops every cached plan.

The writer's **WriteTimeSteps** property writes every time step of its input rather than the current one: into one file when the format holds a series (XDMF, VTKHDF, `frd`, PVD and the other [`write_sequence`](./sequences.md) targets, which hold every step until the last arrives), else one file per step, `name_0000.vtu`, `name_0001.vtu`, ..., each with its time, written as it arrives.

## Named regions

Regions cross as arrays by default, so they survive a round trip through ParaView and can be thresholded or coloured by: a point or cell region is a 0/1 `region:<name>` point or cell array (the masks [`to_pyvista`](./interop.md) writes), with each region's `dim` and `tag` in a `meshioplusplus:regions` JSON field array; a side region, which VTK has no `(cell, facet)` concept for, is the [VTU convention](./regions.md)'s `region:side:<name>` field array of `(cell, facet)` pairs with its `region-meta:side:<name>`. The writer turns all of them back into regions, and so it does the `region:<kind>:<name>` field arrays of a meshio++ `.vtu` that ParaView opened with its own reader.

With **RegionsAsBlocks** the reader outputs a `vtkPartitionedDataSetCollection` instead, one block per region in a `vtkDataAssembly` (groups `cell_regions`, `side_regions` and `point_regions`), as ParaView's own Exodus reader shows its blocks and sets, so a region is a block to tick in the *Multiblock Inspector*: a cell region's piece of [`split(by="regions")`](./regions.md) (regions may overlap, and a cell can be in several blocks), a side region's facets, a point region's points as vertices. A mesh with no cell region is one more block, `mesh`.

## The meshio++ Pipeline filter

*Filters → meshio++ → meshio++ Pipeline* (`meshioPipeline` in `pvpython`) applies the [settings pipeline](./pipeline.md)'s `Operations` list, given as JSON in its **Operations** property, to its input through [`apply_operations`](./pipeline.md#in-memory-apply-operations): every operation the pipeline has, with its validation, in one filter. Any dataset or composite input is accepted, as by the writer; the output is an unstructured grid, and the run's report (`{"steps": [...], "warnings": [...]}`) is its `meshioplusplus:report` field array, with the warnings on ParaView's output window. An unknown operation or key fails the filter with the pipeline's error, which names it.

```json
[{"Op": "Refine"}, {"Op": "Smooth", "Iterations": 5}]
```

The model steps (`MeshCadSurface`, `FitCad`, `ProjectToCad`, the IGA steps) are reachable only within one run of the filter, since a CAD model cannot travel through VTK; `Split` and a fanned-out `Partition`, which write a file per piece, need `run_pipeline` instead.

## Dedicated filters

The operations ParaView has no equivalent for are filters of their own under *Filters → meshio++*, each a thin proxy over the [Pipeline filter](#the-meshio-pipeline-filter) with typed properties: one per key of its [pipeline step](./pipeline.md), with the kind and default `tools/pipeline_steps.json` gives it (`_pipeline_table.KEY_DEFAULTS`), and an **ExtraKeys** JSON object for the keys that take a list (`Refine`'s `Cells`, `Tetrahedralize`'s sizes, MMG's options, ...). Only the properties that differ from their defaults reach the step.

| Filter | Operation |
| --- | --- |
| Remesh (surface), Remesh Volume, Tetrahedralize, Octree Mesh, Boundary Layer, Optimize Volume, Decimate Volume | [`remesh`](./remesh.md), [`remesh_volume`](./remesh_volume.md), [`tetrahedralize`](./tetrahedralize.md), [`octree_mesh`](./octree_mesh.md), [`boundary_layer`](./boundary_layer.md), [`optimize_volume`](./optimize_volume.md), [`decimate_volume`](./decimate_volume.md) |
| Refine, Agglomerate, Tetrahedra to Polyhedra, Polyhedra to Tetrahedra | `refine` with its selectors, `agglomerate`, `tet_to_poly`, `poly_to_tet` |
| Hessian, Estimate Error, Error Size, Adapt (MMG) | the adaptation chain; Adapt needs the MMG add-on |
| Repair, Find Intersections, Thickness, Check Conformity, Select | [`repair`](./repair.md), [`find_intersections`](./find_intersections.md), [`thickness`](./thickness.md), `check_conformity`, [`select`](./select.md) |
| Region Adjacency, Find Interface, Reorder | the interface tools and `reorder` |

Five more call their function directly, since they are not one in-memory pipeline step:

| Filter | Operation |
| --- | --- |
| Partition | `partition` with **GhostLayers**: one partition per piece in a `vtkPartitionedDataSet` (the pipeline step writes a file per piece). |
| Conservative Interpolate | `conservative_interpolate`, the **Source** input's data onto the **Input** (two inputs), with **OnConflict**. |
| Shrinkwrap | `shrinkwrap`, the **Input**'s points onto the **Target** surface. |
| Undo Green | `undo_green`, the **Input** (a `refine` output) with its green cells restored from **Coarse**. |
| Match Periodic Nodes, Contact Pairs | `match_periodic_nodes` as a `periodic:master` point array (-1 off the slave region); `contact_pairs`, the slave points as vertices carrying `contact:gap`, `contact:normal`, `contact:master_cell` and the rest. |

Not wrapped, because ParaView has them natively (check per release): slice, contour, gradient, clip, threshold, extract surface, clean, merge blocks, transform, reflect, extrusion, surface decimation, smoothing, normals, curvature, mesh quality, integrate variables, resampling, calculator, point and cell data conversion, tessellation, subdivision, Delaunay 3D, redistribution, feature edges, Hausdorff distance and connectivity. The [Pipeline filter](#the-meshio-pipeline-filter) still reaches every pipeline operation, these included.

## What crosses into ParaView

The reader hands ParaView the mesh's points (padded to three columns when the mesh is 1-D or 2-D), its cells, and its point, cell and field data, in one vectorised pass with no per-cell Python loop:

- Every cell keeps meshio++'s node order, which is VTK's: a linear `wedge` (Gmsh's prism) is VTK's wedge as it stands, so every wedge ParaView shows has a positive volume and a wedge made or edited in ParaView is written the right way out (meshio++'s VTK writers used to permute it, which turned it inside out; see the [changelog](https://github.com/loumalouomega/meshioplusplus/blob/main/CHANGELOG.md)).
- A block of uniform n-gons (`polygon3`, `polygon5`, ..., the boundary blocks the OpenFOAM, Fluent and CGNS readers produce) is a `VTK_POLYGON` block.
- The writer keeps the grid's cell order: each run of consecutive cells of one type (and, for polygons, one size) becomes a cell block, so a mixed grid can come back as several blocks of one type rather than one block per type.
- Polyhedra cross in both directions. A `polyhedron<n>` block becomes `VTK_POLYHEDRON` cells, whose points are the distinct nodes of their faces, with the faces set through `SetPolyhedralCells` on VTK 9.4 and later or the face stream `SetCells` took before (the same `set_cells` serves `to_pyvista`); the writer reads either back, one `polyhedron<n>` block per run of equal point count.
- The fixed-order high cells VTK has no id of its own for (`triangle10`, `quad16`, `tetra20`, `wedge40`, `hexahedron64`, ...) are VTK Lagrange cells, their Gmsh node order permuted into VTK's, as the [VTU writers](./formats/vtu.md#cell-types) write them.
- Arbitrary-order VTK Lagrange and Bézier cells keep their degrees and weights: a `HigherOrderDegrees` cell array and a `RationalWeights` point array (the names the [VTU format](./formats/vtu.md) uses) become the active attributes of those names, which is where VTK's Lagrange and Bézier cells read them, and come back as arrays of the same names.
- A point or cell array meshio++ cannot hold (a `vtkStringArray`, say) is skipped with a warning naming it; a string array of field data (the QA records of ParaView's Exodus reader) is kept as a string array.

## Architecture: the payload layer

The plugin file declares the ParaView proxies and nothing else. The conversions live in `meshioplusplus._paraview`, which is the [Blender add-on](./blender.md#architecture-the-pure-payload-layer)'s pattern: the mesh-to-VTK direction is the payload `to_pyvista` uses (`_interop._to_vtk_payload`), the VTK-to-mesh direction is the VTU reader's `_vtk_common.vtk_cells_from_data`, and the two VTK adapters over them (`grid_from_mesh`, `mesh_from_grid`) import VTK inside each function. They need VTK but not ParaView, so `tests/python/test_paraview_payload.py` runs them with the PyPI `vtk` wheel in the default test matrix.

## Limitations

- Optional format dependencies (`h5py`, `netCDF4`, ...) must be installed in ParaView's Python as well; the bundle carries meshio++'s compiled HDF5 and netCDF support where its wheel has it, not the Python packages.
- A `pyramid14` block and a `hexahedron1331` block are dropped with a warning: VTK has no 14-node pyramid, and an order-10 hexahedron has no Gmsh element to check its node order against. A `custom` block is dropped likewise.
- A file series is timed as [sequences](./sequences.md) time one, so a number in the file names (`out_0003.vtu`) wins over the time stored in each file.
- Side regions, which VTK has no `(cell, facet)` concept for, cross as field arrays: ParaView carries them through a pipeline but cannot show them; *RegionsAsBlocks* shows them as blocks of facets.
- The writer holds every step of a series-format output until the last one arrives (one-file-per-step formats are written as each step comes).

## Testing

`tests/python/test_paraview_payload.py` tests the conversions with the PyPI `vtk` wheel and no ParaView, and `tests/python/test_paraview_plugin_algorithms.py` runs the plugin's reader, writer and filters through a VTK pipeline on that wheel, with stand-ins for ParaView's `smproxy` decorators. `tests/python/test_paraview_plugin.py` runs `pvpython` in a subprocess: it loads the plugin, reads files through the reader (wedges must have positive volumes in ParaView's `CellSize`, `polygon5` blocks and 1-D meshes must keep their cells), writes them back through the writer, and reads the result with meshio++. It is skipped when `pvpython` is absent or runs a different Python than the build under test; CI runs it in its `paraview-plugin` job against conda-forge ParaView, with HDF5 and netCDF on so XDMF, Exodus and CGNS files go through the plugin too, together with `tests/python/test_pvd_paraview.py`, the ParaView-side checks of the `.pvd` and `.pvtu` files meshio++ writes.
