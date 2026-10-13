# Long-run spikes: findings

[Roadmap §6](./roadmap.md#_6-long-run-spike-first) admits work whose shape is unknown until an investigation writes it down: findings before code. This page is that write-up for the three groups of §6 — scale, NURBS, and skin to CAD — as of October 2026. Each finding ends in a decision, and the roadmap was narrowed to match.

## Scale at 10M cells

**The question.** Memory-mapped reads halve the peak of a large read, and the XDMF series appender and Python's chunked `write_dataset` write a series or a dataset without holding it, but nothing writes one mesh larger than memory and no operation streams. The roadmap asked for the operations harness to be run at `--tier XL` (a Kuhn cube of 10,368,000 tetrahedra on 1,771,561 points) before deciding whether streaming writes (6.1.1) or out-of-core operations (6.1.2) matter.

**The harness itself did not fit.** `meshioplusplus_bench_ops --tier XL` was killed by the out-of-memory killer at 13.6 GB on a 15 GiB machine, during setup, before its first row: it built every derived input of every row up front (the jittered, doubled, quadratic, green-refined, narrowed and ragged copies of the cube) and held them all. The harness now builds each derived input on first use, so a run selected with `--ops` holds only what its rows read; the change is digest-neutral (every row's digest at tier S is identical before and after). That is a finding about the harness, not the library — but it is the first sign of what a 10M-cell workload costs when several copies of a mesh coexist.

**Operations.** Each row was then run alone, `SEQ` backend, one thread, one timed run after the warmup, on a 4-core machine with 15 GiB and no swap; the peak is the process's resident high-water mark. Building the input cube and its skin alone peaks at 6.4 GB (the generator's own temporaries; the cube itself is 0.4 GB), so a peak equal to that means the operation stayed under it.

Of 76 rows, 70 completed, 54 of them in under 30 seconds on one thread; the slowest that finished were `optimize_volume` (823 s), `pipeline_smooth_surface` (345 s), `smooth_volume` (317 s) and `interpolate_cells` (268 s). The heaviest that fit added about 5 GB over the setup peak (`feature_edges` on the volume, `agglomerate`, the skin-extracting renders): roughly a dozen times the 0.4 GB cube. Six did not finish. Three were killed by the out-of-memory killer at 13.3 GB, and all three read a copy of the cube that the harness builds for them: `reorder_ragged` and `simplexify_ragged` the cube as 10M polyhedra, `undo_green` a green-refined cube with its refinement hierarchy recorded. `simplexify_ragged` was rerun alone on an idle machine and died the same way. Three exceeded the 15-minute cap on a warmup plus one timed run (so more than about 7.5 minutes each) without using memory beyond the setup peak: `remesh_volume`, and the two rows that sample the distance from every volume point to the skin, `sample_distance` and `distance_to_surface`. Times are one thread on the `SEQ` backend; the four-thread OpenMP times of the render rows in the [terminal-rendering page](./tui.md#design-decisions) are about six times shorter.

| Row | Time (s) | Peak (GB) | |
| --- | ---: | ---: | --- |
| `extract_surface` | 29.38 | 6.3 |  |
| `extract_skin` | 28.39 | 6.3 |  |
| `smooth` | 0.98 | 6.3 |  |
| `smooth_volume` | 317.41 | 7.3 | +1.1 GB over setup |
| `pipeline_smooth_surface` | 344.85 | 8.1 | +1.8 GB over setup |
| `refine` | 16.40 | 6.3 |  |
| `elevate` | 8.33 | 6.3 |  |
| `merge` | 9.07 | 6.3 |  |
| `clean` | 11.35 | 6.3 |  |
| `clean_weld` | 25.93 | 6.3 |  |
| `compute_sdf` | 21.73 | 6.3 |  |
| `decimate` | 24.04 | 6.3 |  |
| `partition` | 13.24 | 6.3 |  |
| `reorder` | 6.77 | 6.3 |  |
| `reorder_hilbert` | 1.83 | 6.3 |  |
| `linearize_ragged` | 3.10 | 6.3 |  |
| `subdivide_ragged` | 36.96 | 10.0 | +3.7 GB over setup |
| `quality` | 9.27 | 6.3 |  |
| `quality_ragged` | 36.37 | 6.3 |  |
| `data_average` | 1.03 | 6.3 |  |
| `data_average_ragged` | 6.07 | 6.4 |  |
| `interpolate_cells` | 268.24 | 7.8 | +1.5 GB over setup |
| `partition_ragged` | 33.87 | 6.7 | +0.4 GB over setup |
| `feature_edges_ragged` | 37.52 | 9.9 | +3.6 GB over setup |
| `read_ragged_ensight` | 3.75 | 6.3 |  |
| `read_ragged_tecplot` | 30.67 | 10.9 | +4.6 GB over setup |
| `read_stl` | 0.36 | 6.3 |  |
| `read_obj` | 0.13 | 6.3 |  |
| `write_obj` | 0.32 | 6.3 |  |
| `read_tetgen` | 5.69 | 6.3 |  |
| `read_triangle` | 6.06 | 6.3 |  |
| `read_freefem` | 2.84 | 6.3 |  |
| `read_ugrid` | 1.93 | 6.3 |  |
| `read_openfoam` | 28.47 | 10.5 | +4.2 GB over setup |
| `read_elmer` | 7.34 | 8.7 | +2.4 GB over setup |
| `optimize_volume` | 822.85 | 10.5 | +4.2 GB over setup |
| `segment` | 0.12 | 6.3 |  |
| `fit_primitives` | 0.07 | 6.3 |  |
| `reorder_ragged` | — | 13.3 | killed (OOM) |
| `agglomerate` | 51.66 | 11.0 | +4.8 GB over setup |
| `region_adjacency` | 39.25 | 10.3 | +4.1 GB over setup |
| `find_interface` | 36.20 | 7.4 | +1.1 GB over setup |
| `split` | 2.67 | 6.3 |  |
| `undo_green` | — | 13.3 | killed (OOM) |
| `hessian` | 21.42 | 7.9 | +1.7 GB over setup |
| `compute_normals` | 0.10 | 6.3 |  |
| `compute_curvature` | 0.09 | 6.3 |  |
| `feature_edges` | 56.85 | 11.3 | +5.0 GB over setup |
| `render_frame_160x96` | 19.83 | 10.2 | +3.9 GB over setup |
| `render_ss4_800x480` | 19.92 | 10.2 | +3.9 GB over setup |
| `render_streamlines_320x192` | 26.60 | 10.5 | +4.3 GB over setup |
| `render_synthwave_320x192` | 20.31 | 10.5 | +4.3 GB over setup |
| `hausdorff` | 118.76 | 6.3 |  |
| `shrinkwrap` | 6.36 | 6.3 |  |
| `voxelize` | 0.02 | 6.3 |  |
| `remesh` | 1.08 | 6.3 |  |
| `remesh_volume` | > 450 | 10.0 | timeout |
| `sample_distance` | > 450 | 6.3 | timeout |
| `distance_to_surface` | > 450 | 6.3 | timeout |
| `distance_build` | 0.12 | 6.3 |  |
| `isosurface` | 3.09 | 6.3 |  |
| `slice` | 2.12 | 6.3 |  |
| `linearize` | 0.98 | 7.3 | +1.1 GB over setup |
| `interpolate` | 59.72 | 9.9 | +3.7 GB over setup |
| `repair` | 0.33 | 6.3 |  |
| `sobolev_deform` | 0.65 | 6.3 |  |
| `gradient` | 30.48 | 7.2 | +1.0 GB over setup |
| `gradient_gg` | 9.27 | 6.3 |  |
| `gradient_narrow` | 8.00 | 6.3 |  |
| `gradient_surface` | 0.08 | 6.3 |  |
| `diff` | 1.83 | 6.3 |  |
| `diff_narrow` | 2.62 | 6.3 |  |
| `refine_narrow` | 24.91 | 9.7 | +3.4 GB over setup |
| `linearize_narrow` | 1.78 | 7.3 | +1.1 GB over setup |
| `elevate_narrow` | 8.32 | 7.6 | +1.3 GB over setup |
| `simplexify_ragged` | — | 13.3 | killed (OOM) |

**Writes.** One XL cube (357 MB of points and connectivity) written through the Python API with the native writers where the build has them, the resident set sampled during the write:

| Format | Extra memory while writing | Write time | File |
| --- | --- | --- | --- |
| `vtu` (binary, zlib) | +1,024 MB (2.9× the mesh) | 17.7 s | 111 MB |
| `vtk` (legacy binary) | +399 MB (1.1×) | 2.2 s | 476 MB |
| `gmsh` (4.1 binary) | +413 MB (1.2×) | 1.7 s | 450 MB |
| `su2` | +334 MB (0.9×) | 7.3 s | 438 MB |
| `exodus` (netCDF4, Python) | +333 MB (0.9×) | 3.8 s | 357 MB |
| `xdmf` (HDF5 through h5py) | +29 MB | 5.4 s | 61 MB |
| `medit` | +12 MB | 4.5 s | 432 MB |

**Decision on 6.1.1 (streaming writes of one mesh).** No writer needs more than three times the mesh it writes, and most need about one: a mesh that fits in memory writes. A mesh that does *not* fit has no producer in meshio++ to begin with — every reader, operation and binding builds a whole `Mesh` — so a general streaming writer would serve only callers that generate a mesh piecewise, and none has asked. The item is narrowed to the one real multiple: the VTU writer holds every array's compressed form beside the raw one; compressing and writing block by block brings it to the others' footprint. **Kept, narrowed, S–M.**

**Decision on 6.1.2 (out-of-core operations).** Not needed for 10M cells: every operation that failed did so holding at most a few copies of a mesh that itself fits comfortably, and the ones that ran long did so in memory. Streaming the block-local operations through a mesh that does not fit would not fix any row above. The item is replaced by what the measurement did find: the **peak memory of the copying paths** — the ragged (polyhedral) paths of `reorder` and `convert_cells` and `undo_green` hold their input, their output and full-size scratch at once, more than 7 GB over the setup peak for a 0.4 GB tetrahedral cube at XL, and are the first to fail as meshes grow. **Replaced, S–M**, with `--tier XL --ops` as the probe. The two distance rows' run time is a performance question, recorded here rather than reopened: the roadmap's performance section is closed and its [landing rules](./benchmarks.md#landing-a-performance-change) apply to any change.

### What landed

Both items were then fixed with every output byte unchanged: the VTU, VTP, VTR, VTS and VTI files written before and after are identical, and every operation digest of the harness at tiers S and M agrees with the earlier build on the `SEQ` backend and with `OPENMP` at one and four threads.

**6.1.1, the VTK XML writers.** Each binary array now goes through a streaming encoder that asks the writer for its bytes 8 MiB at a time, compresses that batch's 32 KiB blocks in parallel into one buffer of compressed bytes, and writes the base64 text as it goes; the VTU writer generates its padded points and, for fixed-size cell blocks, the connectivity, offsets and types a window at a time instead of building them whole. Only the compressed bytes are ever held (the block header, written first, needs every block's size). Writing the XL cube as binary zlib VTU now needs +178 MB, half the mesh, where it needed +1,024 MB, and takes 11.2 s instead of 17.7 s. Ragged and polyhedral blocks and ASCII output still build their cell arrays before writing them, through the same encoder.

**6.1.2, the copying paths.** Each row's largest scratch is released or never built:

- `reorder` on a ragged mesh: the node adjacency was a vector per node fed every ordered pair of a cell's node *list*, so a polyhedron, which lists each corner once per face, added each edge nine times and the rows rounded up to 8 KiB per node (14 GB at XL). It is now one flat array of raw entries over each cell's distinct nodes, 32-bit when the ids fit, counted and then filled; the rows are the same sets, so the CSR and the Reverse Cuthill-McKee order are too.
- `convert_cells(simplexify)`: the tetrahedra and triangles are written straight into their final arrays, where they were staged in a vector and copied (4 GB at XL); the parent of each child is built only when cell data or `record_parent_ids` reads it; the sources of the new points are freed before the blocks are added.
- `undo_green`: the classification's scratch (copies of the three `refine:*` columns, the coarse id lookup, the sibling order) is released before the output is built, roles are one byte, and the sibling sort is skipped when the parent ids arrive in order, as `refine` writes them.

| Row | XL before | XL after | XL time after | L peak before → after |
| --- | --- | --- | --- | --- |
| `reorder_ragged` | killed at 13.3 GB | 6.3 GB (the setup peak) | 42.5 s | 1,209 → 470 MB |
| `simplexify_ragged` | killed at 13.3 GB | 9.7 GB | 61.7 s | 1,111 → 726 MB |
| `undo_green` | killed at 13.3 GB | killed at 13.3 GB | — | 1,294 → 1,039 → 994 MB |

`undo_green` is still killed at XL, and by its input more than by itself: the harness's green-refined cube (every seventh cell refined and closed red-green, about 83M tetrahedra with their hierarchy) leaves the process at 11.1 GB once built, and the 3.6 GB output does not fit beside it. At tier L its peak above that input fell from 0.48 GB to 0.23 GB. Reading the `refine:*` columns in place and narrower group ids were to bring it to its floor of input plus output, about 16 GB here. They landed next (the third L figure above): each `refine:*` column is now read through a view per block, which is the array itself for the Int64 columns `refine` writes (a column of another dtype is converted one block at a time), the group ids are 32-bit and the per-group sizes are gone, and with a coarse mesh that records no hierarchy, as the harness's, the coarse id table is never built. At XL that removes about 2 GB of column copies, 0.33 GB of group ids and 0.25 GB of coarse ids from the classification; at tier L the row also runs in 0.24 s instead of 1.34 s, and every digest at tiers S, M and L is unchanged. What is left is the floor itself, the input and the output side by side, which does not fit in this machine's 15 GiB; [roadmap 6.1.1](./roadmap.md#_6-long-run-spike-first) keeps only the probe on a larger one.

## How far the data model stretches

**The question (formerly roadmap 6.2.1).** A `CellBlock` is a cell-type string and a node-index array; a rational Bézier or NURBS patch also needs weights and knot vectors, which have nowhere to live. Does a side-channel struct suffice, or does the `Mesh` need a second entity kind?

**What exists.** The [MFEM](./formats/mfem.md) reader is the only NURBS reader: it samples each patch into its knot-span elements as Lagrange cells (`quad9`, `hexahedron27`, the VTK Lagrange types for higher orders), evaluates NURBS grid functions at the nodes, and discards the knots, weights and control points. VTK's Bézier cell types (ids 75–81) exist only in the name tables, and nothing reads `RationalWeights` or `HigherOrderDegrees`, so a rational Bézier VTU loses its weights. The precedent for format data the `Mesh` cannot hold is the side-channel struct — `MedInfo`, `GmshInfo`, `ExodusInfo`, `MdpaInfo` — filled by a reader overload, consumed by a writer overload, carried explicitly by the Python binding and through `mio_format_info` in the C API, and pinned by the ABI layout tests.

**Finding.** A `NurbsInfo` side channel — per patch: the parametric dimension, the degrees, one knot vector per direction, the weights and the control net as an `NDArray`, plus the map from the sampled Lagrange cells back to their patch and knot span — carries everything read-only ingestion needs, and gives MFEM's patches a round trip, without touching `Mesh`, the mesh backends or any existing layout. It does not make NURBS a first-class citizen: no operation could act on a patch, and a trimmed surface's trimming curves would be one more side table. A real IGA model (patches, control nets and trimming curves as entities beside cells, with evaluation) needs a second entity kind in the `Mesh`, every backend, every binding and the ABI: a major version. **Decision: the side channel first, filled by MFEM (roadmap 6.2.1, M); the IGA model stays XL and unscheduled (6.2.3).**

**What landed.** The side channel, as found (v16.40.0): `NurbsInfo` (`nurbs.hpp`) holds per patch the degrees, the knot vector of each direction, the control net, the weights and the control values of every NURBS field read with it, then the file's patch topology and the cell-to-patch-and-span map. The MFEM reader fills it in both engines and the MFEM writer writes it back as an `MFEM NURBS mesh` in MFEM's global form, so all six NURBS fixtures make an exact round trip of mesh, field and side channel, and the two engines write the same bytes; it reaches C, Fortran, Julia, R, WASM and the MCP server through the format side-channel paths ([MFEM](./formats/mfem.md#the-patches-themselves-the-nurbs-side-channel)). The Bézier half took a different shape than the roadmap first sketched: the VTK Bézier types became cell types (ragged, in the Lagrange node order), and the VTU `RationalWeights` and `HigherOrderDegrees` markers name a point array and a cell array of those names on every surface, so a rational Bézier VTU needs no side channel; a [VTU](./formats/vtu.md#bezier-cells-and-rational-weights) written with them is read by VTK with its weights. Checking them against VTK turned up an older defect, fixed with them: VTK 9.1 swapped the nodes of edges 10 and 11 of its higher-order hexahedra and marks the new order by file version (XML 2.2, legacy 5.1), so the `VTK_LAGRANGE_HEXAHEDRON` cells meshio++ wrote as VTU 0.1 or legacy 4.2 were read by VTK with those nodes exchanged. The writers now declare XML 2.2 when such cells are present and convert on a 4.2 write, and the readers convert an older file. The IGA model (now 6.2.2) stays XL and unscheduled.

**Trimmed B-reps.** The `NurbsInfo` finding held for trimmed geometry too, as a second side channel rather than a wider first one. Kratos keeps its IGA geometry in CAD JSON beside the `.mdpa` deck, whose `Geometries` blocks take only Lagrange names; `BrepInfo` (`brep.hpp`) carries that file's faces, trimming loops, edges, vertices and the members meshio++ does not model, exactly, and STEP, IGES and `fit_cad` fill the same struct, so a CAD model becomes Kratos geometry by a format conversion ([Kratos CAD JSON](./formats/cad_json.md)). The sampling follows the non-goal on constrained triangulation: each knot span is a Lagrange cell, and a trimmed face keeps the cells whose centre its loops contain. Two small core-private kernels came with it, a JSON DOM (the core had only the optional nlohmann submodule) and a NURBS evaluator whose fixed arithmetic order its Python twin repeats, so both engines sample bit for bit.

## The OpenCASCADE dependency

**The question (formerly roadmap 6.2.4 and 6.3.5).** Reading STEP or IGES robustly, and writing fitted analytic faces as a valid solid, both mean a CAD kernel; the only open one of production quality is [Open CASCADE Technology](https://github.com/Open-Cascade-SAS/OCCT).

**Facts checked.** OCCT is distributed under the GNU LGPL 2.1 with the additional exception in its `OCCT_LGPL_EXCEPTION.txt`; the current major version is 8.0. The Python route stl2curves uses, the `cadquery-ocp` wheels, installs OCCT 8.0.1 (`libTKernel.so.8.0.1`, `cadquery-ocp` 8.0.1.1) and with its VTK dependency took 1.2 GB here. It is what validated the [`step`](./formats/step.md) writer: `STEPControl_Reader`, `BRepCheck_Analyzer` and `GProp` read every file the tests write as a valid solid of the mesh's volume. [stl2curves](https://github.com/Duvdivan/stl2curves), the reference for the whole skin-to-CAD pipeline, is MIT-licensed and built on `cadquery-ocp`, NumPy and SciPy; meshio++ read its README only.

**Finding.** Everything that needs no kernel has shipped: [`segment`](./segment.md), [`fit_primitives`](./fit_primitives.md) and the faceted [`step`](./formats/step.md) writer, on every surface (see [skin to CAD](./cad_fit.md)). The STEP writer needed no kernel because a faceted solid is only planes and lines; the roadmap's description of it — `ADVANCED_FACE`s bounded by `POLY_LOOP`s — is not valid AP214 (an advanced face takes `EDGE_LOOP`s), and the writer emits the edge-loop form, which also lets neighbouring faces share their edges. What remains — CAD faces on the fitted cylinders, cones, spheres and tori, B-spline faces for the rest, fillets and threads, and reading STEP — is trimming, sewing and validating surfaces, which is the kernel's work and not something to re-implement.

**Decision.** An OpenCASCADE backend, if it comes, follows the KaHIP / ADIOS2 / MMG rule: `MESHIOPLUSPLUS_WITH_OCCT`, off by default, bring-your-own through `find_package(OpenCASCADE CONFIG)`, never vendored or downloaded, linked dynamically, absent from the wheels, the WASM build and the native CLI binaries, its licence implications stated on its page, and `Unsupported` naming the flag wherever it is absent (Python raises `NotImplementedError`; there is no numpy twin of a kernel). One backend serves reading (now 6.2.1) and writing (6.3.3). **Not scheduled**: the fits are already useful as measurements and the faceted solid already opens in CAD; the backend is worth its weight when a consumer asks for *editable* CAD faces rather than a measured skin.

**Taken.** The backend landed under exactly that rule, with OpenCASCADE 7.6 as the floor (Ubuntu 24.04's packages, which the CI leg installs): [`fit_cad`](./fit_cad.md) builds the analytic faces, sews and checks the solid and writes STEP, demoting any face the kernel rejects to facets. Two things the kernel did not do by itself had to be done around it: parameter-space curves continuous around each loop, since OpenCASCADE's own projection splits a boundary that runs along a seam and its missing-seam repair then drops the wire; and a sphere's pole axis chosen away from the patch, so a cap or a band closes in parameter space. Reading (6.2.1) and the freeform fallback (6.3.2) then landed on the same backend. [STEP](./formats/step.md#reading) and [IGES](./formats/iges.md) read as a `BRepMesh` triangulation plus the shape as a NURBS B-rep; converting the whole shape with `BRepBuilderAPI_NurbsConvert`, not each face's surface on its own, was what kept every trimming curve in its converted surface's parameters. `fit_cad` fits a B-spline face to a patch no primitive explains, over its plane or an unrolled cylinder, and holds the fit to half the tolerance at the vertices and the tolerance at the facet centres, since the solid's own check measures the face against the skin between vertices, where a smooth surface lies a facet's chordal error from it.
