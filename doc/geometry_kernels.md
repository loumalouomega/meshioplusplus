# Shared geometry and FEM kernels

Many operations evaluate the same small formulas per cell: a tetrahedron's signed volume, the scaled Jacobian at a corner, a centroid, a cross product. Before November 2026 most of them carried private copies. This page records where the shared versions live, the rules for using and extending them, and the design question of how far to take them: a layer of plain kernels, or a geometry class hierarchy in the style of [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos). The cleaning work still open is [roadmap item 3.4](./roadmap.md#quality-of-implementation), ordered by SonarCloud's duplication report ([static analysis](./static_analysis.md)).

## Where the shared code lives

This page covers the geometry and FEM kernels and the binding helpers. Every other shared utility, from number parsing to row gathers, region remaps and the Python shims, is listed by role in [common utilities](./common_utilities.md).

| Header | Contents |
| --- | --- |
| `include/meshioplusplus/detail/geometry.hpp` (installed) | `Vec3`, `vec3_sub`, `vec3_add`, `vec3_scale`, `vec3_dot`, `vec3_cross`, `vec3_norm_sq`, `vec3_norm`, `vec3_normalize`, `triple_product`, `det3`; `read_point`, `read_corner_coords`, `cell_corner_count`. |
| `src/cpp/src/detail/fem_kernels.hpp` (private) | `tet_volume6` and `tet_volume`; `det_unit3`, the per-corner factor of the scaled Jacobian that `quality` reports; `corner_jacobian_ratio`, the same factor as one division; `kHexCornerEdges` and `kWedgeCornerEdges`, the corner frames of a hexahedron and a wedge; `tet_min_dihedral` (moved verbatim from `optimize_volume`, whose digests pin it), `tet_dihedral_range`, `tet_mean_ratio` and `tet_amips` (its reciprocal), the tetrahedron quality measures `optimize_volume` and `adapt_volume` share. |
| `src/cpp/src/detail/cell_centroids.hpp` (private) | `cell_centroids`, per-cell centroids of a whole mesh in global cell order (`partition`, `interpolate`). |
| `include/meshioplusplus/detail/polyhedron.hpp` | Polygon and polyhedron measures: `poly_measure`, `polygon_area_vector`, `polygon_centroid`, `cell_volume_from_corners`. |
| `include/meshioplusplus/detail/cell_faces.hpp`, `cell_edges.hpp`, `cell_subdivision.hpp` | Topology tables per cell type: outward faces, edges, refinement edges and quad faces. |
| `src/cpp/src/detail/predicates.hpp` | Exact `orient2d`, `orient3d`, `incircle`, `insphere`. |
| `src/cpp/src/detail/math_utils.hpp` (private) | `MathUtils`: `Sub`, `Add`, `Scale`, `Dot`, `Cross`, `NormSq`, `Norm`, `Distance` over `Vec3` and flat `double` triples (the `Vec3` ones forward to `geometry.hpp`); `UnitByDivision`, the divide-by-length unit vector with a zero guard; `AngleDeg`, `atan2(‖a × b‖, a · b)` in degrees; `SolveSymmetric3`, the cofactor solve of a symmetric 3×3 system with the `1e-12 s³` singularity bound. |
| `src/cpp/src/detail/geometry_utils.hpp` (private) | `GeometryUtils`: `TriangleNormal` (unnormalised, `(b − a) × (c − a)`), `TriangleTwiceArea`, `TriangleArea`, `TriangleAnglesDeg`, over `Vec3` and flat triples; `SignedTwiceArea2D` in a parameter plane. |
| `src/cpp/src/detail/bounding_box_utils.hpp` (private) | `BoundingBox3` and `BoundingBoxUtils`: `FromPoints` (seeded with the first point), `FromCoords` (coordinates gathered on the fly), `Expand`, `Extent`, `MaxExtent`, `HalfSurfaceArea`, `Diagonal`. |
| `src/cpp/src/detail/facet_slots.hpp` (private) | `FacetSlotUtils`: the global facet-slot numbering of `extract_surface` and `smooth`'s boundary pass (a polyhedron block's face prefix sum, its inverse, the block owning a slot). |
| `src/cpp/src/detail/dtype_table.hpp`, `ndarray_utils.hpp` (private) | `DTypeUtils::Pick`, one value per `DType` from a table in enumerator order (the netCDF, HDF5, XDMF, VTK and PLY type names); `NDArrayUtils::Int64Vector` and `Float64Vector`, owning 1-D arrays from `std::vector`s. |

The private headers carry no ABI cost, so a kernel starts there. It moves to an installed header only when an installed inline function needs it, and that move follows the [ABI policy](./abi.md).

## Rules

- **Do not add a file-local copy.** Use the kernels above, and add the missing formula to them instead of to the operation.
- **Utility classes, Kratos-style.** A group of related helpers is a class of static inline members (`MathUtils::Cross(a, b)`, `GeometryUtils::TriangleArea(a, b, c)`, `BoundingBoxUtils::FromPoints(points)`), as Kratos' `MathUtils` and `GeometryUtils` are: no state, called in place of the block it replaced. Each starts in a private header under `src/cpp/src/detail/`, included directly and by relative path from every `.cpp` that uses it (the amalgamation emits only the private headers a source names), and has its own unit test (`tests/cpp/test_*_utils.cpp`) pinning known values and bit equality with the expression it replaced.
- **Keep the arithmetic.** A shared kernel replaces a copy only when it computes the same operations in the same order. `a·(b×c)` and a cofactor expansion of the same determinant are equal in exact arithmetic but round differently, and the operations' outputs are pinned bit for bit by the [determinism check](./benchmarks.md#determinism-check). A copy whose arithmetic differs stays until its own output is shown unchanged, or until a change in its results is intended and recorded.
- **Keep deliberate copies.** A few copies say why they exist. `smooth_facefan_volume` is a flat-buffer variant of `quality`'s face-fan volume, kept out of a hot loop's allocation. `grad_ring_area_vector` is a recentred Newell sum, better conditioned than the plain one. Those comments stay with the code.
- **Prove it.** Run `tools/bench_ops.sh` with `--hash` before and after the change, with `BASELINE` pointing at the first run. Run the [instruction-count gate](./benchmarks.md#instruction-count-gate) too: an inline kernel must cost what the copy cost.

## What the first pass shared

These kernels replaced copies with identical arithmetic:

- the tetrahedron's six-times signed volume (`quality`, `optimize_volume`, `tetrahedralize`, `isosurface_stuffing`, `feature_snap`, `conservative_interpolate`);
- the unit-edge corner determinant (`quality`, `optimize_volume`);
- the hexahedron and wedge corner tables (`quality`, `hex_fit`, `layer_sweep`);
- the centroid gather (`partition`, `interpolate`);
- 70 file-local sub, dot, cross, norm, add, scale and distance helpers across 27 sources, now `geometry.hpp`'s.

The third pass (November 2026) added the utility classes and moved onto them, with the same arithmetic:

- the unnormalised triangle normal and its length (`external_mesher`, `mesh_boolean`, `hull`, `cut_graph`, `remesh` and seven sites of `mesh_cad`), and `mesh_cad`'s corner angles (`quality`'s `tri_area` lambda moved too and was put back: GCC compiled the rest of `quality` 1.8% longer, [roadmap 3.4.11](./roadmap.md#quality-of-implementation));
- the first-point bounding boxes of `point_normals`, `registration`, `reconstruct_surface`, `external_mesher` and `interpolate`;
- the symmetric 3×3 cofactor solve of `decimate`'s optimal placement and `gradient`'s least squares (`decimate` negates the solution rather than the right-hand side, which keeps the sign of a zero);
- the division-form unit vectors (`ql_unit`, `ff_unit`, `cf_unit`, `bl_unit` and the external mesher's STL normals) and the `double*` dot, cross and difference helpers of ten sources.

The unit vectors without a zero guard, the boxes seeded with ±infinity and the Cramer-rule solves (`quality_solve3`, `remesh_solve_sym3`) round or behave differently and stay where they are. What remains is listed under roadmap 3.4.

## What the binding pass shared

The second duplication pass (November 2026, [roadmap §5](./roadmap.md#code-duplication)) gave the binding layers and a few operations their own helpers, each file-private and each a call in place of the block it replaced:

| Where | Helper | Replaces |
| --- | --- | --- |
| `src/cpp/src/detail/data_plumbing.hpp` | `append_rows(a, rows, new, comps, fill)`: a copy of an array with rows appended, the new values written by the caller's `fill` through one `DoubleSink` commit | the point and point_data growth of `convert_cells`' simplexify and elevate, `refine` and `tet_to_poly` (each caller keeps its mean, midpoint or key average) |
| `src/cpp/src/detail/decimate_shared.hpp` | `decim_emit_blocks<K>`, the kept cells per block and their cell maps; `decim_remap_regions` | the emission tails of `decimate` (first surviving child of `convert:parent_cell`, `FirstChild`) and `decimate_volume` (identity, `Direct`) |
| `src/cpp/src/operations/primitives.cpp` | `prim_soup(map, radius, center)`, `prim_cell_type_is(shape, cell_type, alt, default, choices)` | 19 soup set-ups and 16 two-way cell-type checks |
| `src/cpp/src/registry.cpp` | `RegistryRow` and `registry_rows()`: one row per format with its reader, writer, options reader and metadata reader | the four tables `registry_readers()`, `registry_writers()`, `registry_readers_ex()` and `registry_metadata_readers()` are built from |
| `src/cpp/cli/main.cpp` | `CliVerb` and `kCliVerbs[]`: name, alias, usage block, handler | `print_usage`'s verb list and `cli_run_verb`'s if chain |
| `bindings/c/c_api.cpp` | `capi_result_mesh`, `capi_result_take_mesh`, `capi_result_point_map`, `capi_result_num_cell_maps`, `capi_result_cell_map`; `capi_decimate_options`; `capi_indexed_reason`, `capi_copy_indexed_edges` | the bodies of 29 result accessors, the four decimate entry points' option fills, `mesh_cad_surface`'s seven indexed getters |
| `bindings/python/_core.cpp` | `core_def_mesh_fn`, `core_def_parallel_xml`, `core_def_decimate`, `core_def_step_writer`, `core_def_finalizing_context`, `core_def_time_series_read`, `core_named_arrays` | the near-copy `m.def`s of the primitives, pvtu/pvtp, the decimate pair, the step writers, the finalizing contexts and the time-series readers |
| `src/python/meshioplusplus/mcp/_server.py` | `_forward(target)`: a tool whose body only forwards its own arguments | 170 forwarding bodies (the tools that remap or render stay hand-written) |
| `bindings/r/meshioplusplus/src/util.c` | `mio_r_report_list` (a named list whose `NULL` entries take the next of an array of doubles), `mio_r_wrap_handle`, `mio_r_handle` | 47 `PROTECT(Rf_ScalarReal(…))` report chains and the five handle types' wrap and check |

Its second half added:

| Where | Helper | Replaces |
| --- | --- | --- |
| `src/cpp/src/detail/always_inline.hpp` | `MIO_ALWAYS_INLINE` | nothing by itself: it lets `stage_flat_or_polygon`/`emit_flat_or_polygon` serve `subdivide` and `ccells_finish_identity` end `convert_cells`' linearize and elevate without changing their callers' code generation ([roadmap 3.4.11](./roadmap.md#quality-of-implementation)) |
| `src/cpp/src/detail/nastran_write_common.hpp` | `NastranWriteBlock`, `NastranWriteRow`, `nastran_card_index`, `nastran_property_values`, `nastran_named_arrays`, `nastran_fill_block`, `nastran_read_components`, `nastran_read_field`, `nastran_compact_block`, `nastran_rows_by_card` | the copies of the Nastran HDF5 and OP2 writers' model building |
| `bindings/c/c_api.cpp` | `capi_take_mesh_maps`, `capi_grid_extent`, `capi_grid_bounds`, `capi_checked_enum` | the result moves of six operations, the lattice extents of four converters, 22 enum range checks |
| `bindings/python/_core.cpp` | `core_distance_options`, `core_metric_summary_dict`, `core_quality_dict` | four distance-option fills, three quality report conversions |
| `src/cpp/cli/main.cpp` | `cli_opt_double`, `cli_opt_int64`, `cli_opt_int`, `cli_opt_string`, `cli_print_warnings`, `cli_grid_bounds`, `cli_decimate` | 195 option reads, 27 warning loops and arrays, five lattice blocks, the decimate pair's bodies |
| `src/python/meshioplusplus/mcp/_tools.py` | `_run_op(op, locals(), ...)` | the bodies of 54 tools that only forward their parameters |
| `src/python/meshioplusplus/_shims.py` | `plain_shims`, `vtk_xml_shims` | the read and write shims of 20 format packages |
| `tools/pipeline_steps.json`, `tools/gen_pipeline_table.py` | the step tables of both pipeline engines, generated | `_pipeline.py`'s `_OP_TABLE` and `pipeline.cpp`'s `pipe_op_table()`, kept by hand |
| `bindings/r/meshioplusplus/src/util.c`, `bindings/octave/meshioplusplus/src/util.c` | `mio_r_getbuffer`, `mio_r_cstrings`; `mio_o_getbuffer`, `mio_o_row3`, `mio_o_surface_quality`, `mio_o_ids_0based` | the per-operation getter adapters, three character-vector conversions, the Octave files' own copies |


The installed-header pass (November 2026, ABI 42, [roadmap §5](./roadmap.md#code-duplication)) closed the last three items:

| Where | Helper | Replaces |
| --- | --- | --- |
| `src/cpp/src/detail/subset_gather.hpp` | `subset_gather_points`, forced inline into `build_cell_subset` and `clean` | the points-and-point_data gather of `build_cell_subset` and `clean` ([roadmap 3.4.11](./roadmap.md#quality-of-implementation)); `build_cell_subset` compiles as it did with its own copy |
| `src/cpp/src/detail/fit_cad_native.hpp` | `native_local_triangles`, `native_adaptive_options`, `native_resegment_pieces`, `native_refine_proxy` | the OpenCASCADE provider's copies of the native `fit_cad`'s local triangulation, adaptive-fit options, re-segmentation rung and refinement verdicts |
| `src/cpp/src/detail/step_brep.hpp` | `step_boundary_chains`, `StepCurveIndex`, `step_ring_next`, `step_fan_next` | the boundary-chain walk, the network-curve lookup and the facet-fan turn of the native STEP writer and the OpenCASCADE provider |
| `src/cpp/providers/occt/` | `fc_occt_knots`, `fc_occt_bspline_curve`; `occt_weld_nodes`, `occt_welded_triangles` (`occt_brep.hpp`) | two pole and knot builds and a third for curves; the node welding of `fc_triangulate` and `occt_tessellate` |
| `detail/spatial_hash.hpp` | `grid_for_each_in_shell`, `grid_for_each_in_box` | the shell and box walks of `SpatialGrid` and `BucketTable` |
| `detail/link_check.hpp` | `MESHIOPLUSPLUS_PP_CAT`, `MESHIOPLUSPLUS_PP_STR`, `MESHIOPLUSPLUS_LINK_CHECK_REF` | the paste, stringize and `gnu::used` reference of the ABI and backend sentinels |
| `formats/time_series_common.hpp` | `SeriesMode`, `SeriesNamedArray`, `MESHIOPLUSPLUS_SERIES_WRITER_SPECIAL_MEMBERS` | the XDMF and VTKHDF series modes and raw arrays; the special members of 19 series writers |
| `formats/drawing_options.hpp` | `DrawingView`, `DrawingColoring`, in `SvgWriteOptions` and `TikzWriteOptions` | the camera and colouring parameters `write_svg` and `write_tikz` spelled apart |
| `brep.hpp`, `operations/sdf_names.hpp` | `BrepFaceCounts`; `kSdfDistanceName`, `kSdfClosestPointName` | the nine face counts of `FitCadResult` and `StepBrepResult`; the distance names `octree_mesh` and `compute_sdf` both spelled |

The exported names, error messages, status codes and Python signatures are unchanged: the C++ and C API ones are pinned by the determinism digests, the gtest suite and byte diffs of the CLI help and completion scripts, the registry's key sets and the `_core` docstrings taken before the change; the MCP tool schemas by a diff of `list_tools()`.

## Kernels or a geometry hierarchy?

Kratos models each element shape as a class: `Geometry<TPointType>` is the base, and `Triangle2D3`, `Tetrahedra3D4`, `Hexahedra3D8` and the rest derive from it. Each instance holds its points and answers `ShapeFunctionsValues`, `ShapeFunctionsLocalGradients`, `Jacobian`, `DeterminantOfJacobian`, `IntegrationPoints` and `ShapeFunctionsIntegrationPointsGradients`, with the shape-function values at the integration points computed once per geometry type and shared. Two designs are possible for meshio++.

### Option A: a static kernel layer (recommended)

Free inline functions over gathered corners, plus `constexpr` per-type tables, dispatched by the `switch (CellType)` the operations already use:

- **Measures and frames**: what `fem_kernels.hpp` holds today, plus barycentric coordinates, dihedral angles and the P1 gradient of a simplex.
- **Reference elements**: per cell type, its corner coordinates, its corner frames and its faces and edges (reusing `cell_faces` and `cell_edges` rather than adding a second table), as `constexpr` data.
- **Interpolation**: shape functions `N` and their local derivatives `dN/dξ` per cell type and order, then `jacobian(corners, dN)`, its determinant and inverse, and the physical gradients.
- **Quadrature**: Gauss–Legendre points and weights on the line, and the standard rules on the triangle, tetrahedron, quadrilateral, hexahedron, wedge and pyramid.

This layer fits how operations already work. They gather corners through the uniform mesh API (`CellView`, `read_corner_coords`) and loop over cells in `parallel_for`. Kernels on the stack cost nothing over the copies they replace, need no allocation per cell, and work unchanged on the MESHIO, NATIVE and KRATOS mesh backends. Private kernels are no ABI surface.

### Option B: a Kratos-style geometry hierarchy

A `Geometry` base class with virtual shape functions, Jacobians and integration rules, one derived class per cell type, constructed per cell. In its favour:

- FEM users already know its vocabulary.
- One object can cache shape-function data per integration point across many evaluations.
- An operation can be written against the base class with no per-type switch.

Against it, for this library:

- **Cost in the hot loops.** A virtual call or an object per cell inside the per-cell loops would move the instruction counts the gate pins. Meshes store connectivity, not geometry objects, so each cell would have to be wrapped first.
- **A second hierarchy.** The KRATOS mesh backend already wraps Kratos' own geometries, and a second hierarchy would shadow them.
- **ABI and API surface.** An installed class hierarchy is layout and ABI surface (every virtual function counts), and every binding would want it exposed.
- **Little need today.** Few operations integrate: `data_integrate` uses a one-point measure, `cut` uses the external `cutcells` rules, and `surface_fem` assembles P1 triangles in closed form. None needs cached integration-point data.

### A middle path

A non-virtual `GeometryView` can give Kratos-like call sites without virtual dispatch: a `CellView` plus the points, whose methods (`Volume()`, `Jacobian(xi)`, `ShapeFunctions(xi)`) forward to the Option A kernels through a `switch` on the cell type or a template on it. It can be added on top of Option A whenever an operation reads better with it, and costs nothing until then.

### Decision rule

Build Option A first. Take Option B only when an operation needs integration-point data cached across many evaluations of the same cells, for example an assembly-heavy operation or a solver-facing interpolation. Even then, it belongs in a private header first, built on Option A's kernels, so the arithmetic stays in one place.
