# Common utilities

Before writing a helper, look here. Most small pieces an operation, a format or a binding needs (a number parser, a row gather, a facet table, a region remap, a bounding box, a CLI option reader) already exist once, shared, with tests and with the arithmetic the determinism digests pin. This page lists them by role, with the header or module to include and what it holds. The geometry and FEM kernels and the binding helpers have their own page, [shared geometry and FEM kernels](./geometry_kernels.md), linked below rather than repeated.

## Adding or extending a utility

1. **Search first.** Grep for the formula or loop you are about to write, read the tables below, and run `python3 tools/dup_report.py --lang all` on your change ([static analysis](./static_analysis.md#local-duplication-report)): a copy it reports is a helper you should have called.
2. **Extend the closest utility.** A missing variant (another dtype, a flat-buffer overload, an extra out-parameter) belongs beside the existing one, in its header, with its test.
3. **Otherwise add a shared one.** A helper two places need goes in a header both include, never into a second file's anonymous namespace. Start it core-private, in `src/cpp/src/detail/`, where it carries no ABI cost; move it to `include/meshioplusplus/detail/` only when an installed inline function needs it, following the [ABI policy](./abi.md). A group of related helpers is a Kratos-style class of static members (`MathUtils`, `GeometryUtils`, `BoundingBoxUtils`) with a `tests/cpp/test_*_utils.cpp`; a loop body shared by callers is a template taking the caller's part as a functor (`step_boundary_chains`, `occt_welded_triangles`). Anonymous-namespace helpers still take a prefix unique to their `.cpp`, because the [single header](./cpp_api.md) puts every source in one translation unit.
4. **Keep the arithmetic.** A share replaces a copy only when it computes the same operations in the same order; prove it with `tools/bench_ops.sh --hash` against a `BASELINE` and with the [instruction-count gate](./benchmarks.md#instruction-count-gate). When a share moves a caller's code generation, force it inline with `MIO_ALWAYS_INLINE` (`detail/always_inline.hpp`), and check the inlining of the caller's own loops too: an extra out-of-line function in a translation unit can cost its other functions their inlining (`detail/subset_gather.hpp` records a case).
5. **List it here** in the same change, in the table it belongs to.

## C++: installed `detail/` headers

Under `src/cpp/include/meshioplusplus/detail/`. Installed, so an inline body or a layout change there follows the [ABI policy](./abi.md).

### Numbers, text and bytes

| Header | Contents |
| --- | --- |
| `fast_number.hpp` | Locale-independent number parsing and formatting: `parse_double`, `snprintf_c`. Every float read or written goes through it. |
| `classic_stream.hpp` | Streams pinned to the classic locale: `make_classic_ifstream`, never a bare `std::ifstream`. |
| `parse_guard.hpp` | The checks a reader makes before trusting a count, an offset or a token from the file. |
| `read_guard.hpp` | Makes every reader entry point fail with `ReadError`, whatever the parser threw. |
| `binary_stream.hpp`, `byteswap.hpp` | A bounds-checked cursor over a byte buffer in a chosen byte order; the swap primitives. |
| `file_source.hpp`, `compressed_input.hpp`, `zlib_inflate.hpp` | Whole-file read access (memory-mapped where it helps); compressed input (`.gz`, `.bz2`, `.xz`, `.zst`, `.zip`) read as if it were not; one deflate stream of unknown size. |
| `fortran_records.hpp`, `keyword_card.hpp` | Fortran sequential records split and framed; the fixed-width card tokenizer of keyword decks. |
| `format_compat.hpp`, `source_location_compat.hpp` | Stand-ins for `std::format` and `std::source_location` where the toolchain lacks them. |

### Arrays, data and regions

| Header | Contents |
| --- | --- |
| `value_io.hpp` | Scalars read out of an `NDArray` whatever its runtime dtype (`read_int`, `read_double`, …). |
| `data_ops.hpp` | The data operations' shared vocabulary: array location lookup, owned copies (`data_owned_copy`), component checks. |
| `named_arrays.hpp`, `map_order.hpp` | Insertion-ordered `name -> NDArray` containers; the sorted-name order every writer iterates in (`sorted_keys`). |
| `ragged_csr.hpp` | The CSR form of ragged (polygon, polyhedron) blocks, its conversions and checks. |
| `subset.hpp` | `build_cell_subset`, a submesh from kept cells (points pruned, data gathered, regions remapped), and `subset_gather_rows`. |
| `region_remap.hpp` | Carrying regions through an operation: `Direct`, `FirstChild` and `Global` maps, side facets by containment ([regions](./regions.md)). |
| `mesh_carve.hpp` | A mesh split into named cell-index parts, for the multi-part writers. |
| `provenance.hpp`, `length_unit.hpp` | The provenance record every writer renders; the length-unit field-data key and its accessors. |

### Topology

| Header | Contents |
| --- | --- |
| `cell_faces.hpp`, `cell_edges.hpp`, `cell_subdivision.hpp` | Per-cell-type face, edge, refinement-edge and quad-face tables. |
| `cell_index.hpp` | The global (block-major) cell index and its inverse. |
| `cell_adjacency.hpp`, `node_adjacency.hpp`, `facet_index.hpp` | Cell-to-cell adjacency by shared node; the node graph in CSR; a cell facet found by its nodes. |
| `node_order.hpp`, `vtk_cells.hpp` | Node-ordering permutations between formats and the VTK layout; cell blocks from VTK connectivity, offsets and types. |
| `refine_templates.hpp`, `refine_hierarchy.hpp` | `refine`'s subdivision templates; the `refine:cell_id`/`refine:parent_id` hierarchy reader. |
| `face_mesh.hpp`, `planar_skin.hpp`, `degenerate_solid.hpp` | Volume cells as a deduplicated face list; the 2-D skin of `voxelize` and `octree`; a brick with repeated nodes collapsed to its solid. |

### Geometry and search

| Header | Contents |
| --- | --- |
| `geometry.hpp`, `polyhedron.hpp` | `Vec3` arithmetic, `triple_product`, `det3`, corner readers; polygon and polyhedron measures ([kernels](./geometry_kernels.md)). |
| `spatial_hash.hpp`, `bucket_table.hpp` | The growable bucket grid `SpatialGrid`, its frozen twin `BucketTable`, and the shell and box walks both run (`grid_for_each_in_shell`, `grid_for_each_in_box`). |
| `point_triangle.hpp`, `tri_box.hpp`, `surface_distance.hpp`, `surface_normals.hpp` | Closest point on a triangle and its feature; exact triangle/box overlap; the signed-distance kernel, with `query_surface_projections_filtered` for a nearest-triangle search restricted per query (a `TriangleFilter`: `deviation`'s far-side rule, `thickness`'s sphere method); vertex normals of a soup. |
| `space_filling.hpp`, `projection.hpp`, `sym3_eigen.hpp`, `marching.hpp` | Morton and Hilbert keys; the orthographic camera and painter's order of the drawing writers; symmetric 3×3 invariants; the marching-tetrahedra cutter. |

### Format families and link checks

| Header | Contents |
| --- | --- |
| `vtk_xml.hpp`, `vtu_binary.hpp`, `xdmf_common.hpp`, `hdf5_util.hpp` | VTK XML `<DataArray>` helpers and codecs; XDMF type names and raw cell data; HDF5 helpers. |
| `nastran_model.hpp`, `ansys_model.hpp`, `abaqus_types.hpp`, `colormap.hpp`, `face_color.hpp` | The models the Nastran, Ansys and Abaqus readers share; the drawing writers' colormaps and per-face colours. |
| `link_check.hpp`, `abi_version_check.hpp`, `mesh_backend_check.hpp` | The two link-time sentinels and the paste, stringize and `gnu::used` reference they share. |

Beside them, `formats/time_series_common.hpp` and `formats/drawing_options.hpp` what the format headers share (`SeriesMode`, `SeriesNamedArray`, `MESHIOPLUSPLUS_SERIES_WRITER_SPECIAL_MEMBERS`; `DrawingView`, `DrawingColoring`).

## C++: private `src/cpp/src/detail/` headers

Core-private: no installed header names them, so they change freely. The general-purpose ones:

| Header | Contents |
| --- | --- |
| `math_utils.hpp`, `geometry_utils.hpp`, `bounding_box_utils.hpp` | The Kratos-style utility classes `MathUtils`, `GeometryUtils`, `BoundingBoxUtils` ([kernels](./geometry_kernels.md)). |
| `fem_kernels.hpp`, `cell_centroids.hpp`, `predicates.hpp`, `expansion.hpp` | Tet volume, scaled-Jacobian corner factors and corner tables; per-cell centroids; exact predicates and expansion arithmetic. |
| `size_field_utils.hpp` | `SizeFieldUtils`: metric edge length of a geometrically graded isotropic size field, the metric midpoint of an edge, size interpolation and the Dijkstra gradation `h(v) <= h(u) + (hgrad - 1) |uv|` (`adapt_volume`). |
| `tet_flips.hpp` | `optimize_volume`'s worst-cell moves over a flat xyz buffer and a vector of tetrahedra, moved verbatim so `adapt_volume` runs them too: the 2-3 and 3-2 flips, edge removal and multi-face removal, the guarded ODT relocation and the worst-cell smoothing, each touching only faces and edges interior to the vector it is given. |
| `tet_remesh.hpp`, `tet_scope.hpp` | `TetRemeshStore`, the editable tetrahedral store of `adapt_volume`: tombstoned tetrahedra, per-vertex incident tetrahedra, the tagged boundary and interface faces, ridge chains parametrised by arc length, the corner / ridge / surface / interior classification and the invariant checker; `require_tet_only_blocks`, the tet-only scope check `optimize_volume` and `adapt_volume` share. |
| `facet_slots.hpp`, `slot_runs.hpp`, `dtype_table.hpp`, `ndarray_utils.hpp` | `FacetSlotUtils`; the sort-based facet and edge table; `DTypeUtils::Pick`; `NDArrayUtils` owning arrays from vectors. |
| `typed_view.hpp` | An array's elements as `int64_t` or `double`, the dtype switch hoisted (`DoubleView`, `DoubleSink`). |
| `vtk_write_cell.hpp`, `vtk_hex_order.hpp` | The VTK cell type and node order a VTK writer writes a block as (`vtk_type_for_write`, `vtk_order_for_write`: VTK's own types, and the fixed-order high types as VTK Lagrange cells); VTK 9.1's edge order of higher-order hexahedra across file versions. |
| `data_plumbing.hpp`, `subset_gather.hpp`, `staged_block.hpp` | Row plumbing of the topology operations (`compact_points_by_mask`, `replicate_parent_rows`, `append_rows`; `stack_padded_rows`, a multi-zone reader's points stacked with narrower zones zero-padded); the points-and-point_data gather `build_cell_subset` and `clean` share; an unchanged block staged and emitted. |
| `gzip_deflate.hpp` | `gzip_compress`: one gzip member byte-identical to Python's `gzip.compress(data, mtime=0)`, so a native writer and its twin agree (libMesh `.gz`, OpenFOAM `writeCompression on`); `gzip_compress_available`. The reading side is `compressed_input.hpp`'s `decompress_bytes`. |
| `always_inline.hpp` | `MIO_ALWAYS_INLINE` and `MIO_ALWAYS_INLINE_LAMBDA`, for a share the instruction gate shows moved its caller's code generation. |
| `text_cursor.hpp`, `chunked_parse.hpp`, `row_writer.hpp`, `open_source.hpp`, `cow_source.hpp` | Text readers' lines, tokens and numbers; numbers parsed in parallel; ASCII rows formatted in parallel and streamed in order; a reader's file through `FileSource`; copy-on-write file bytes. |
| `exact_weld.hpp`, `weld.hpp`, `node_table.hpp`, `utf8.hpp` | First-seen numbering of identical rows (STL, glTF, WKT); keep-first welding (`clean`, `merge`); the TetGen/Triangle `.node` table; UTF-8 encoding. |
| `point_knn.hpp`, `tri_bvh.hpp`, `tri_pairs.hpp`, `tri_tri.hpp`, `bucket_build.hpp` | k-nearest points over a frozen grid; a triangle BVH; overlapping triangle pairs; exact triangle/triangle intersection; the grid build. |
| `primitive_constraints.hpp` | `fit_cad`'s constraint pass (coaxial groups, parallel and perpendicular planes, equal radii) with `PrimitiveConstraintControl`'s tolerances and force/forbid lists and a record of every relation (`fit_cad`, `constrain_primitives`). |
| `hole_fill.hpp` | A faired patch over one boundary loop (Liepa 2003): minimum-area triangulation, refinement to the rim's density with Delaunay flips, membrane or thin-plate fairing through `sparse_solve` (`repair`'s `fill_style`). |
| `hausdorff_samples.hpp` | Every sample distance `hausdorff_distance` reduces, kept (`hausdorff_sample_distances`, and `HausdorffSampler` with the reference's query built once), and `coverage_quantile` (`decimate`'s `max_deviation`). |
| `offset_normals.hpp` | Each point's offset direction over a triangle and quad surface (Aubry and Löhner's most normal normal, smoothed) and its stretch, one over the least cosine to its faces (`boundary_layer`, `offset_surface`, `thicken`, `hollow`). |
| `cell_flip.hpp` | The node permutation that turns a cell inside out and where it sends its facets; `BlockFlip` and `block_flip(cells, op)` per block, `block_flip_facet` for a `Side` entry (`mirror`, `flip_normals`). |
| `vertex_cell_csr.hpp`, `cell_components.hpp`, `dual_graph.hpp`, `crease_edges.hpp`, `closed_facets.hpp`, `face_skin.hpp`, `cell_region_mask.hpp` | Cells around each vertex; node-connected components and their edge-connected sibling `cell_components_by_edge` (a flood fill over the facet dual graph: `remove_islands`, `deviation_regions`); the facet dual graph; the one crease test; a surface as a closed facet shell; the skin of a selection; a region as a cell mask. |
| `dense_chol.hpp`, `sparse_solve.hpp`, `jacobi_eigen.hpp`, `weighted_matching.hpp` | Dense Cholesky and LDLᵀ; sparse CSR and conjugate gradients; small symmetric eigen-decomposition; maximum-weight matching. |
| `nurbs_eval.hpp`, `bspline_ndu.hpp`, `iga_kernels.hpp` | B-spline and NURBS evaluation; the shared Cox–de Boor triangle; knot insertion and degree elevation. |
| `decimate_common.hpp`, `decimate_shared.hpp`, `nastran_write_common.hpp`, `transfer_names.hpp`, `mesh_from_soup.hpp`, `polyhedron_groups.hpp`, `sfc_quantize.hpp` | The two decimators' quadric machinery and emission; the Nastran writers' model; the transfer operations' store and diagnostics; a one-block mesh from a soup; polyhedra grouped by node count; space-filling-curve quantization. |
| `fit_cad_native.hpp`, `step_brep.hpp` | The native `fit_cad` helpers the OpenCASCADE provider builds on (`native_local_triangles`, `native_adaptive_options`, `native_resegment_pieces`, `native_refine_proxy`, …); the kernel-free STEP writer and the chain walks both backends share (`step_boundary_chains`, `StepCurveIndex`, `step_fan_next`). |
| `primitive_grow.hpp` | Regions grown from a seed facet while one primitive explains them: the surface builder (rings, normals, centroids, smooth edges, neighbours), the smooth areas and `primitive_grow_from_seed`, which `segment`'s primitive method and `fit_primitives`' seed mode share; `nearest_centroid_face`. |
| `polyline_chain.hpp` | Segments walked into polylines in a canonical order (open chains from their lowest-id end, then loops from their lowest point, each walk to the lowest-id neighbour): `polyline_chains`, which `slice(chain=True)` and `sections` use; the numpy twin `_slice._polyline_chains` replicates it. |
| `sketch_pieces.hpp` | A `fit_curves` sketch read back from the `sketch:*` field data as pieces in world coordinates (`origin + u e1 + v e2`): `has_sketch`, `sketch_world_pieces`, which the DXF and STEP writers draw from. |
| `dense_lp.hpp` | A small dense linear-programme solver: the two-phase simplex with Bland's rule, `dense_lp_equality` (`min c.y`, `A y = b`, `y >= 0`, with the rows' duals) and `dense_lp_inequality` (`min c.x`, `A x <= b`, `x` free, solved on its dual: one row per variable, a column per constraint). Deterministic. |
| `zone_fit.hpp` | Minimum-zone, minimum-circumscribed and maximum-inscribed fits of a plane, sphere or cylinder by a sequence of `dense_lp` programmes in a trust region: `zone_fit`, which `fit_primitives`' `criterion` calls. |
| `planar_union.hpp` | The region a soup of planar triangles covers minus the region another covers, exactly: candidate boundary edges (on one triangle, or on two on the same side), split where they cross (`cdt2d_find_crossings`), constrained in a `Cdt2d`, coverage decided once per component of triangles joined across unconstrained edges (`orient2d` over a bucket grid), and the boundary walked into loops (region on the left, touching loops kept apart): `planar_region`, `planar_loop_twice_area`. `silhouette` unions a projected surface with it and offsets the outline through it. |
| `lsq_fit.hpp` | Small dense least squares: `lsq_sym3_eigen`, `lsq_solve`, `lsq_frame`, Levenberg–Marquardt, and the circles `lsq_circle_kasa` (Kåsa's algebraic circle, the arithmetic `fit_cad`'s cylinder unroll has always used, lifted bit for bit) and `lsq_circle_geometric` (Gauss–Newton on the geometric distance), which `fit_curves` fits its arcs with. |
| `primitive_geom.hpp`, `hex_color.hpp` | Queries against a fitted primitive's ten `fit:patch_params`: distance, signed distance (`primitive_signed_distance`), projection, normal, and the form error and extremes of points under it (`primitive_form_metrics`); `parse_hex_color`, the `#rrggbb` every `nan_color` option takes. |
| `json_value.hpp`, `subprocess.hpp`, `dynlib.hpp`, `library_lock.hpp`, `parallel_workers.hpp` | A small JSON model; running an external program; dynamic loading; the lock around non-thread-safe libraries; the backend's worker count. |

Headers for one algorithm (`cdt2d`, `hex_fit`, `layer_sweep`, `isosurface_stuffing`, `mesh_boolean`, `quad_layout`, `patch_network`, `cross_field` and the like) are modules of their operation rather than utilities; reuse them when you need that algorithm, but do not grow them into general helpers. The OpenCASCADE provider's own shared pieces live in `src/cpp/providers/occt/occt_brep.hpp` (`occt_weld_nodes`, `occt_welded_triangles`, `occt_tessellate`, `OcctQuiet`).

## Bindings and CLIs

The C API's `capi_result_*`, `capi_take_mesh_maps`, `capi_grid_extent` and `capi_checked_enum`; `_core.cpp`'s `core_def_*`, `core_named_arrays` and `core_distance_options`; the native CLI's `kCliVerbs` and `cli_opt_*`; the R binding's `mio_r_report_list`, `mio_r_wrap_handle` and `mio_r_getbuffer`; the Octave binding's `mio_o_row3` and `mio_o_getbuffer`: see [what the binding pass shared](./geometry_kernels.md#what-the-binding-pass-shared).

## Python

Under `src/python/meshioplusplus/`.

| Module | Contents |
| --- | --- |
| `_helpers.py` | Format dispatch and registration (`register_format`, `read`, `write`, extension lookup). |
| `_fallback.py` | Whether a declined C++ fast path may fall back to its Python twin: `core_declined`, `core_op_declined`; `MESHIOPLUSPLUS_STRICT_CORE=1` makes every fallback an error. |
| `_shims.py` | `plain_shims` and `vtk_xml_shims`, the read and write shims of the format packages that only try one `_core` entry point. |
| `_common.py`, `_vtk_common.py`, `vtk/_legacy_common.py`, `vti/_vti.py` | Cell-type tables, warnings and XML helpers; the VTK cell conversions (with `FIXED_LAGRANGE`, the fixed-order high types written as VTK Lagrange cells); the legacy VTK sections and writer pieces; the lattice reader and writer. |
| `_lagrange.py` (C++ twin `formats/lagrange_common.hpp`) | Arbitrary-order Lagrange cells: VTK's node lattice (`vtk_lattice`), Gmsh's (`gmsh_lattice`) and the permutation between them (`gmsh_to_vtk_order`), lattice weights and nodal interpolation. |
| `_set_remap.py`, `_side_carry.py`, `_regions.py` | The shim-only `point_sets`/`cell_sets` carried through an operation; side regions by containment (the twin of `region_remap.cpp`); region helpers. |
| `_data_common.py` | The data operations' shared vocabulary. |
| `_ragged.py`, `_files.py`, `_compressed.py` | Ragged blocks that carry their CSR (`CsrRows`, `csr_of`), a polyhedron block's CSR whatever its storage (`polyhedron_csr`) and the gather of segments laid end to end (`ragged_gather`); buffer and path handling; compressed input. |
| `_cli/_common.py` | The CLI's argument groups and value parsers (`add_box_args`, `add_transfer_io_args`, `parse_bounds`, `parse_ints`, `split_csv`, …). |
| `mcp/_server.py`, `mcp/_tools.py` | The MCP server's `_forward` decorator and the tools' `_run_op`. |
| `_pipeline_table.py`, `_format_table.py` | Generated by `tools/gen_pipeline_table.py` and `tools/gen_format_table.py`; edit their sources, not them. `_pipeline_table.KEY_DEFAULTS` holds each step's typed keys with their kinds and defaults, which the ParaView plugin's dedicated filters expose. |
