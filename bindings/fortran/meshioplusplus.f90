!  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
! ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
!  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
!  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
!  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
!  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
!  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
! ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
!
!
!  License:         MIT License
!                   meshio++ default license: LICENSE
!
!  Main authors:    Vicente Mataix Ferrandiz
!
!
! The meshio++ Fortran interface: a modern OO Fortran 2008 module layered on
! the C API (bindings/c/include/meshioplusplus/meshioplusplus.h) via
! ISO_C_BINDING, in the HDF5/PETSc style:
!
!     use meshioplusplus
!     type(mio_mesh) :: m
!     call m%read("bracket.msh")
!     print *, m%num_points()
!     call m%write("bracket.vtu")
!     call m%free()
!
! Conventions (differences from the C API):
!  - Arrays are Fortran-shaped: points are `points(dim, num_points)` and
!    connectivity `conn(nodes_per_cell, num_cells)`. Because Fortran is
!    column-major and the C core row-major, this is the SAME memory -- no
!    transpose happens anywhere in this module.
!  - Connectivity is 1-based here; the +-1 shift happens inside the copying
!    setters/getters (where a copy is made anyway). Zero-copy borrows
!    (points_ptr) carry no indices, so nothing zero-copy ever needs shifting.
!  - All indices (cell blocks, data names) are 1-based.
!  - Every fallible procedure takes `optional` `stat` (integer, 0 = success)
!    and `errmsg` (deferred-length character) arguments. If `stat` is absent
!    and the call fails, the message is printed and the program error stops
!    -- pass `stat` to handle errors yourself (the stdlib pattern).
!  - Handles are freed explicitly with `call m%free()` (no finalizer).
!
! Compiled .mod files are compiler-(major-version-)specific; this source file
! is installed next to the .mod so consumers on a different compiler can
! simply recompile the module (the HDF5 approach).
module meshioplusplus
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64, error_unit
    implicit none
    private

    public :: mio_mesh
    public :: mio_xdmf_series
    public :: mio_exodus_series
    public :: mio_format_info
    public :: mio_iga_evaluation
    public :: mio_mesh_cad_surface_report
    public :: mio_mesh_cad_face_status, mio_mesh_cad_face_reason
    public :: mio_mesh_cad_face_size_limits
    public :: mio_mesh_cad_num_reparametrised, mio_mesh_cad_get_reparametrised
    public :: mio_mesh_cad_num_merged_faces, mio_mesh_cad_get_merged_face
    public :: mio_mesh_cad_num_tolerance_uses, mio_mesh_cad_get_tolerance_use
    public :: mio_mesh_cad_num_lopped_corners, mio_mesh_cad_get_lopped_corner
    public :: mio_mesh_cad_lopped_corner_reason
    public :: MIO_MESH_CAD_METRIC_MEAN, MIO_MESH_CAD_METRIC_VARYING
    public :: MIO_MESH_CAD_TOLERANCE_CLOSED_START, MIO_MESH_CAD_TOLERANCE_CORNER
    public :: MIO_MESH_CAD_TOLERANCE_PROJECTION
    public :: MIO_MESH_CAD_OK, MIO_MESH_CAD_EMPTY, MIO_MESH_CAD_OPEN_LOOP
    public :: MIO_MESH_CAD_SELF_INTERSECTING, MIO_MESH_CAD_DEGENERATE
    public :: MIO_MESH_CAD_PROJECTION_FAILED, MIO_MESH_CAD_BUDGET_EXCEEDED
    public :: MIO_MESH_CAD_UNSUPPORTED, MIO_MESH_CAD_OTHER
    public :: mio_mesh_cad_corner, mio_mesh_cad_angle_histogram
    public :: mio_mesh_cad_num_corners, mio_mesh_cad_get_corner
    public :: MIO_MESH_CAD_SMOOTH_LAPLACIAN, MIO_MESH_CAD_SMOOTH_ODT, MIO_MESH_CAD_ANGLE_BINS
    public :: mio_gltf_options
    public :: MIO_MDPA_PROPERTIES, MIO_MDPA_ENTITY_NAMES, MIO_MDPA_SKIPPED
    public :: MIO_MDPA_MODEL_PART_DATA, MIO_MDPA_TABLES, MIO_MDPA_GEOMETRIES
    public :: MIO_MDPA_MESH_BLOCKS, MIO_MDPA_SUBMODELPARTS, MIO_MDPA_RAW_BLOCKS
    public :: MIO_MDPA_VALUE_NUMBER, MIO_MDPA_VALUE_TEXT, MIO_MDPA_VALUE_TABLE
    public :: mio_stats_report
    public :: mio_octree_mesh_report
    public :: mio_surface_quality
    public :: mio_find_interface_result, mio_contact_pairs_result
    public :: mio_grid
    public :: mio_rectangle, mio_box, mio_disk, mio_uv_sphere, mio_icosa_sphere
    public :: mio_cylinder, mio_reference_cells
    public :: mio_reference_triangle, mio_ngon, mio_annulus, mio_lshape
    public :: mio_plate_with_hole, mio_honeycomb
    public :: mio_frame, mio_render_settings, mio_text_settings, mio_rgba
    public :: mio_detect_color_depth
    public :: mio_data_array_info
    public :: mio_field_integral_info
    public :: mio_convert, mio_version, mio_mesh_backend, mio_error_message
    public :: mio_pipeline_run_file, mio_pipeline_run_json, mio_pipeline_has_json
    public :: mio_has_occt
    public :: mio_pipeline_run_file_report, mio_pipeline_run_json_report
    public :: mio_sequence_pipeline_run_file_report, mio_sequence_pipeline_run_json_report
    public :: mio_med_mesh_names, mio_med_read_named, mio_med_write_multi
    ! `mio_sequence`'s fan-in is the type-bound `%to_timeseries`; only the
    ! fan-out (which starts from a path, not a handle) is module-level.
    public :: mio_sequence, mio_timeseries_to_sequence
    public :: mio_sequence_pipeline_run_file, mio_sequence_pipeline_run_json
    public :: mio_format_readable, mio_format_writable
    public :: mio_sniff_format
    public :: mio_read_metadata, mio_metadata, mio_cell_block_info
    public :: mio_merge
    public :: mio_provenance_begin, mio_provenance_end
    public :: mio_provenance_note, mio_provenance_source, mio_provenance_target
    public :: MIO_PROVENANCE_OFF, MIO_PROVENANCE_BEST_EFFORT, MIO_PROVENANCE_REQUIRED
    public :: mio_interpolate
    public :: mio_conservative_interpolate
    public :: mio_undo_green
    public :: mio_shrinkwrap
    ! Length of the fixed-width string buffers the `keys` out-arguments of
    ! `split` and `data_info` use; consumers need it to declare those arrays.
    public :: STRBUF_LEN
    ! Data-operation enumerations (must match the C enums in meshioplusplus.h).
    public :: MIO_DATA_POINT, MIO_DATA_CELL, MIO_DATA_FIELD
    public :: MIO_WEIGHT_UNIFORM, MIO_WEIGHT_MEASURE
    public :: MIO_COND_CLAMP, MIO_COND_NORMALIZE, MIO_COND_STANDARDIZE
    public :: MIO_SCOPE_COMPONENT, MIO_SCOPE_MAGNITUDE
    public :: MIO_NAN_IGNORE, MIO_NAN_REPLACE, MIO_NAN_FAIL
    ! Named regions (see doc/regions.md).
    public :: MIO_REGION_POINT, MIO_REGION_CELL, MIO_REGION_SIDE
    public :: mio_blend_steps
    public :: mio_region_info
    public :: MIO_TINV_MISES, MIO_TINV_PRINCIPAL, MIO_TINV_HYDROSTATIC, MIO_TINV_DEVIATORIC
    public :: MIO_TINV_ALL

    ! Mirror of mio_write_opts; keep its size and reserved tail unchanged.
    type, bind(c) :: mio_write_opts_c
        integer(c_int) :: encoding = 0, codec = 0
        type(c_ptr) :: float_format = c_null_ptr
        integer(c_int64_t) :: reserved(5) = 0
    end type

    ! Public owning options; strings never have to be passed as C pointers.
    type :: mio_gltf_options
        integer :: container = 0, up_axis = 0, normal_weight = 0
        logical :: normals = .true., fields = .true., recenter = .true.
        logical :: by_region = .true., unlit = .true.
        integer :: component = 1 ! Fortran component indices are 1-based.
        logical :: component_set = .false., vmin_set = .false., vmax_set = .false.
        real(c_double) :: split_angle = 30, scale = 1, vmin = 0, vmax = 0
        character(:), allocatable :: color_by, cmap, nan_color
    end type

    type, bind(c) :: mio_gltf_opts_c
        integer(c_int) :: container = 0, up_axis = 0, normal_weight = 0
        integer(c_int) :: normals = 1, fields = 1, recenter = 1, by_region = 1, unlit = 1
        integer(c_int) :: component = 0, component_set = 0, vmin_set = 0, vmax_set = 0
        real(c_double) :: split_angle = 30, scale = 1, vmin = 0, vmax = 0
        type(c_ptr) :: color_by = c_null_ptr, cmap = c_null_ptr, nan_color = c_null_ptr
        integer(c_int64_t) :: reserved(8) = 0
    end type

    ! What is wrong with a surface (bind(c); layout must match
    ! mio_surface_quality in meshioplusplus.h). The four counts are separate
    ! because they need different fixes: "12 boundary edges" is actionable,
    ! "not watertight" is not.
    type, bind(c) :: mio_surface_quality
        integer(c_int64_t) :: boundary_edges = 0
        integer(c_int64_t) :: non_manifold_edges = 0
        integer(c_int64_t) :: inconsistent_pairs = 0
        integer(c_int64_t) :: degenerate_triangles = 0
        integer(c_int32_t) :: watertight = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_curvature_opts`. Field order and types are ABI
    !> and must match bindings/c/include/meshioplusplus/meshioplusplus.h
    !> exactly; `reserved` is padding for additive growth and must stay zero.
    !> The defaults here mirror `mio_curvature_opts_init`, which is called
    !> anyway -- an all-zero struct is NOT the default, since `mean` and
    !> `gaussian` are on.
    type, bind(c) :: mio_curvature_opts_t
        type(c_ptr) :: region = c_null_ptr
        integer(c_int32_t) :: mean = 1
        integer(c_int32_t) :: gaussian = 1
        integer(c_int32_t) :: dual_area = 0
        integer(c_int32_t) :: include_boundary = 0
        integer(c_int32_t) :: record_area = 0
        integer(c_int32_t) :: record_principal = 0
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_curvature_report`. Field order/types are ABI.
    type, bind(c) :: mio_curvature_report_t
        type(mio_surface_quality) :: quality
        integer(c_int64_t) :: num_boundary = 0
        integer(c_int64_t) :: num_isolated = 0
        integer(c_int64_t) :: num_degenerate = 0
        real(c_double) :: total_angle_defect = 0.0_c_double
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_normals_opts`. Field order and types are ABI
    !> and must match bindings/c/include/meshioplusplus/meshioplusplus.h
    !> exactly; `reserved` is padding for additive growth and must stay zero.
    !> The defaults mirror `mio_normals_opts_init`, which is called anyway --
    !> an all-zero struct is NOT the default.
    type, bind(c) :: mio_normals_opts_t
        type(c_ptr) :: region = c_null_ptr
        integer(c_int32_t) :: point_normals = 1
        integer(c_int32_t) :: cell_normals = 0
        integer(c_int32_t) :: weight = 0
        integer(c_int32_t) :: split = 0
        real(c_double) :: split_angle = 30.0_c_double
        integer(c_int32_t) :: record_parent_ids = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_normals_report`. Field order/types are ABI.
    type, bind(c) :: mio_normals_report_t
        type(mio_surface_quality) :: quality
        integer(c_int64_t) :: num_isolated = 0
        integer(c_int64_t) :: num_undefined = 0
        integer(c_int64_t) :: num_degenerate = 0
        integer(c_int64_t) :: num_split_points = 0
        integer(c_int64_t) :: num_added_points = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_feature_edges_opts`. Field order/types are ABI;
    !> the defaults mirror `mio_feature_edges_opts_init` (called anyway).
    type, bind(c) :: mio_feature_edges_opts_t
        type(c_ptr) :: region = c_null_ptr
        real(c_double) :: feature_angle = 30.0_c_double
        integer(c_int32_t) :: feature = 1
        integer(c_int32_t) :: boundary = 1
        integer(c_int32_t) :: non_manifold = 1
        integer(c_int32_t) :: inconsistent = 1
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_feature_edges_report`.
    type, bind(c) :: mio_feature_edges_report_t
        integer(c_int64_t) :: num_feature = 0
        integer(c_int64_t) :: num_boundary = 0
        integer(c_int64_t) :: num_non_manifold = 0
        integer(c_int64_t) :: num_inconsistent = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_segment_opts`. Field order/types are ABI; the
    !> defaults mirror `mio_segment_opts_init` (called anyway).
    type, bind(c) :: mio_segment_opts_t
        type(c_ptr) :: region = c_null_ptr
        type(c_ptr) :: prefix = c_null_ptr
        real(c_double) :: feature_angle = 30.0_c_double
        real(c_double) :: normal_tolerance = 0.0_c_double
        real(c_double) :: curvature_tolerance = 0.0_c_double
        integer(c_int64_t) :: min_patch_cells = 1
        integer(c_int32_t) :: method = 0
        integer(c_int32_t) :: reserved32 = 0
        real(c_double) :: tolerance = 1.0e-3_c_double
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_segment_report`.
    type, bind(c) :: mio_segment_report_t
        integer(c_int64_t) :: num_patches = 0
        integer(c_int64_t) :: num_merged = 0
        integer(c_int64_t) :: num_feature_edges = 0
        integer(c_int64_t) :: num_boundary_edges = 0
        integer(c_int64_t) :: num_non_manifold_edges = 0
        integer(c_int64_t) :: reserved(3) = 0
    end type

    !> Interop mirror of C `mio_fit_primitives_opts` (`types` is the C
    !> `uint32_t` bit mask, `1 << type`; 62 selects every primitive).
    type, bind(c) :: mio_fit_primitives_opts_t
        type(c_ptr) :: patch_data = c_null_ptr
        real(c_double) :: tolerance = 1.0e-3_c_double
        real(c_double) :: min_inlier_fraction = 0.9_c_double
        integer(c_int32_t) :: relative = 1
        integer(c_int32_t) :: types = 62
        integer(c_int64_t) :: max_iterations = 50
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_fit_primitives_report`.
    type, bind(c) :: mio_fit_primitives_report_t
        integer(c_int64_t) :: num_patches = 0
        integer(c_int64_t) :: num_plane = 0
        integer(c_int64_t) :: num_sphere = 0
        integer(c_int64_t) :: num_cylinder = 0
        integer(c_int64_t) :: num_cone = 0
        integer(c_int64_t) :: num_torus = 0
        integer(c_int64_t) :: num_unfitted = 0
        integer(c_int64_t) :: reserved(5) = 0
    end type

    !> Interop mirror of C `mio_fit_cad_opts` (`types` as for fit_primitives;
    !> 0 facets every patch).
    type, bind(c) :: mio_fit_cad_opts_t
        real(c_double) :: tolerance = 0.0_c_double
        real(c_double) :: mesh_deflection = 0.0_c_double
        integer(c_int64_t) :: max_demotion_rounds = 8
        integer(c_int32_t) :: exact_edges = 1
        integer(c_int32_t) :: types = 62
        integer(c_int32_t) :: freeform = 1
        integer(c_int32_t) :: freeform_degree = 3
        integer(c_int64_t) :: freeform_max_control = 64
        integer(c_int32_t) :: freeform_split = 1
        integer(c_int32_t) :: fillets = 1
        integer(c_int32_t) :: threads = 1
        integer(c_int32_t) :: merge_faces = 1
        integer(c_int32_t) :: backend = 0
        integer(c_int32_t) :: reserved32 = 0
        integer(c_int64_t) :: reserved(1) = 0
    end type

    !> Interop mirror of C `mio_fit_cad_report`.
    type, bind(c) :: mio_fit_cad_report_t
        real(c_double) :: tolerance = 0.0_c_double
        real(c_double) :: volume = 0.0_c_double
        real(c_double) :: skin_volume = 0.0_c_double
        real(c_double) :: volume_ratio = 0.0_c_double
        real(c_double) :: max_deviation = 0.0_c_double
        integer(c_int64_t) :: num_faces = 0
        integer(c_int64_t) :: num_plane = 0
        integer(c_int64_t) :: num_sphere = 0
        integer(c_int64_t) :: num_cylinder = 0
        integer(c_int64_t) :: num_cone = 0
        integer(c_int64_t) :: num_torus = 0
        integer(c_int64_t) :: num_faceted_faces = 0
        integer(c_int64_t) :: num_patches = 0
        integer(c_int64_t) :: num_demoted = 0
        integer(c_int64_t) :: num_exact_edges = 0
        integer(c_int64_t) :: num_approx_edges = 0
        integer(c_int64_t) :: num_solids = 0
        integer(c_int32_t) :: valid = 0
        integer(c_int32_t) :: num_rounded_thread = 0
        integer(c_int64_t) :: num_bspline = 0
        integer(c_int64_t) :: num_split_patches = 0
        integer(c_int64_t) :: num_fillet = 0
        integer(c_int64_t) :: num_helical = 0
        integer(c_int64_t) :: num_variable_fillet = 0
        integer(c_int64_t) :: num_blend = 0
    end type

    !> Interop mirror of C `mio_fit_cad_report_ex` (v17.2.0).
    type, bind(c) :: mio_fit_cad_report_ex_t
        integer(c_int64_t) :: num_vertex_blend = 0
        integer(c_int64_t) :: num_chamfer = 0
        integer(c_int64_t) :: num_conic_blend = 0
        integer(c_int64_t) :: num_tapered_thread = 0
        integer(c_int64_t) :: num_variable_pitch_thread = 0
        integer(c_int64_t) :: backend = 0
        integer(c_int64_t) :: reserved(10) = 0
    end type

    !> Interop mirror of C `mio_refine_patches_opts` (v17.11.0, 144 bytes;
    !> `types` is the C `uint32_t` bit mask and `ransac_seed` the C `uint64_t`
    !> seed, held here as a signed integer of the same bits). The defaults
    !> mirror `mio_refine_patches_opts_init` (called anyway).
    type, bind(c) :: mio_refine_patches_opts_t
        real(c_double) :: tolerance = 1.0e-3_c_double
        real(c_double) :: normal_weight = 1.0_c_double
        real(c_double) :: epsilon = 1.0e-3_c_double
        real(c_double) :: feature_angle = 30.0_c_double
        integer(c_int64_t) :: max_iterations = 10
        integer(c_int64_t) :: min_patch_cells = 3
        integer(c_int64_t) :: ransac_seed = 7882833615064607531_c_int64_t
        integer(c_int32_t) :: types = 62
        integer(c_int32_t) :: relative = 1
        integer(c_int32_t) :: allow_split = 1
        integer(c_int32_t) :: allow_merge = 1
        type(c_ptr) :: prefix = c_null_ptr
        integer(c_int64_t) :: reserved(8) = 0
    end type

    !> Interop mirror of C `mio_refine_patches_report` (v17.11.0, 192 bytes).
    !> `stop` is a `mio_refine_stop`: 0 unchanged, 1 small gain, 2 iteration
    !> cap, 3 energy rose.
    type, bind(c) :: mio_refine_patches_report_t
        integer(c_int64_t) :: num_patches = 0
        integer(c_int64_t) :: num_iterations = 0
        integer(c_int64_t) :: num_splits = 0
        integer(c_int64_t) :: num_merges = 0
        integer(c_int64_t) :: num_absorbed = 0
        integer(c_int64_t) :: num_primitive = 0
        integer(c_int64_t) :: num_freeform = 0
        integer(c_int64_t) :: num_facets = 0
        real(c_double) :: tolerance = 0.0_c_double
        real(c_double) :: explained_area_fraction = 0.0_c_double
        real(c_double) :: initial_explained_area_fraction = 0.0_c_double
        real(c_double) :: hausdorff_max = 0.0_c_double
        real(c_double) :: hausdorff_rms = 0.0_c_double
        real(c_double) :: initial_energy = 0.0_c_double
        real(c_double) :: energy = 0.0_c_double
        integer(c_int32_t) :: stop = 0
        integer(c_int32_t) :: reserved32 = 0
        integer(c_int64_t) :: reserved(8) = 0
    end type

    !> Interop mirror of C `mio_fit_cad_opts_ex` (v17.11.0, 280 bytes): the
    !> fit_cad options `mio_fit_cad_opts` has no room for.
    type, bind(c) :: mio_fit_cad_opts_ex_t
        integer(c_int32_t) :: refine_segmentation = 0
        integer(c_int32_t) :: freeform_general = 0
        type(mio_refine_patches_opts_t) :: refine
        real(c_double) :: freeform_tolerance = 0.0_c_double
        real(c_double) :: freeform_smoothing = 1.0e-12_c_double
        integer(c_int64_t) :: max_freeform_poles = 40
        integer(c_int64_t) :: max_freeform_iterations = 30
        integer(c_int64_t) :: max_freeform_splits = 2
        real(c_double) :: sew_tolerance = 0.0_c_double
        integer(c_int32_t) :: extend_surfaces = 0
        integer(c_int32_t) :: snap_junctions = 0
        integer(c_int32_t) :: face_deviation = 0
        integer(c_int32_t) :: reserved32b = 0
        integer(c_int64_t) :: reserved(8) = 0
    end type

    !> Interop mirror of C `mio_fit_cad_report_ex2` (v17.11.0, 176 bytes).
    type, bind(c) :: mio_fit_cad_report_ex2_t
        integer(c_int64_t) :: refined = 0
        integer(c_int64_t) :: refine_iterations = 0
        integer(c_int64_t) :: num_refine_splits = 0
        integer(c_int64_t) :: num_refine_merges = 0
        real(c_double) :: initial_explained_area_fraction = 0.0_c_double
        real(c_double) :: explained_area_fraction = 0.0_c_double
        integer(c_int64_t) :: num_freeform_general = 0
        integer(c_int64_t) :: num_freeform_refinements = 0
        real(c_double) :: avg_freeform_poles = 0.0_c_double
        integer(c_int64_t) :: num_extended_edges = 0
        integer(c_int64_t) :: num_snapped_junctions = 0
        real(c_double) :: hausdorff_max = 0.0_c_double
        real(c_double) :: hausdorff_rms = 0.0_c_double
        integer(c_int64_t) :: closed = 0
        integer(c_int64_t) :: refinement_kept = 0
        integer(c_int64_t) :: reserved(7) = 0
    end type

    !> Interop mirror of C `mio_fit_cad_opts_ex2` (v17.12.0, 128 bytes): the
    !> freeform-mode options.
    type, bind(c) :: mio_fit_cad_opts_ex2_t
        integer(c_int32_t) :: freeform_mode = 0
        integer(c_int32_t) :: freeform_parameterization = 0
        integer(c_int32_t) :: extraordinary_vertex = 2
        integer(c_int32_t) :: local_refine = 0
        real(c_double) :: freeform_g1_tolerance_deg = 1.0_c_double
        integer(c_int64_t) :: freeform_max_control_points = 400
        real(c_double) :: layout_target_patch_size = 0.1_c_double
        integer(c_int64_t) :: layout_max_valence_deviation = 2
        integer(c_int64_t) :: resegment_rounds = 2
        integer(c_int32_t) :: layout_snap_to_features = 1
        integer(c_int32_t) :: validate_by_surface_sampling = 1
        integer(c_int64_t) :: reserved(8) = 0
    end type

    !> Interop mirror of C `mio_fit_cad_report_ex3` (v17.12.0, 184 bytes).
    type, bind(c) :: mio_fit_cad_report_ex3_t
        integer(c_int64_t) :: brep_checked = 0
        integer(c_int32_t) :: freeform_mode = 0
        integer(c_int32_t) :: reserved32 = 0
        integer(c_int64_t) :: num_freeform_adaptive = 0
        integer(c_int64_t) :: num_network_patches = 0
        integer(c_int64_t) :: num_network_edges = 0
        integer(c_int64_t) :: num_extraordinary_vertices = 0
        integer(c_int64_t) :: num_resegmented = 0
        integer(c_int64_t) :: num_escalations = 0
        real(c_double) :: facet_area_fraction = 0.0_c_double
        real(c_double) :: max_g1_angle_deg = 0.0_c_double
        real(c_double) :: max_freeform_distance = 0.0_c_double
        real(c_double) :: rms_freeform_distance = 0.0_c_double
        integer(c_int64_t) :: num_face_records = 0
        integer(c_int64_t) :: reserved(10) = 0
    end type

    !> Interop mirror of C `mio_iga_opts` (the subdivisions default to 1 and
    !> `trim_samples` to 16: an all-zero struct is not the default).
    type, bind(c) :: mio_iga_opts_t
        integer(c_int32_t) :: subdivisions(3) = 1
        integer(c_int32_t) :: elevate(3) = 0
        integer(c_int32_t) :: sample_subdivisions = 0
        integer(c_int32_t) :: trim_samples = 16
        integer(c_int32_t) :: edges = 1
        integer(c_int32_t) :: trim = 1
        integer(c_int32_t) :: reserved32(2) = 0
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Why `info%mesh_cad_surface` left a face out (C `mio_mesh_cad_face_status`).
    integer, parameter :: MIO_MESH_CAD_OK = 0, MIO_MESH_CAD_EMPTY = 1, &
                          MIO_MESH_CAD_OPEN_LOOP = 2, MIO_MESH_CAD_SELF_INTERSECTING = 3, &
                          MIO_MESH_CAD_DEGENERATE = 4, MIO_MESH_CAD_PROJECTION_FAILED = 5, &
                          MIO_MESH_CAD_BUDGET_EXCEEDED = 6, MIO_MESH_CAD_UNSUPPORTED = 7, &
                          MIO_MESH_CAD_OTHER = 8

    !> How `info%mesh_cad_surface` smooths (C `mio_mesh_cad_smoothing`), and the
    !> angle histogram's bins (10 degrees each, 0 to 180).
    integer, parameter :: MIO_MESH_CAD_SMOOTH_LAPLACIAN = 0, MIO_MESH_CAD_SMOOTH_ODT = 1
    integer, parameter :: MIO_MESH_CAD_ANGLE_BINS = 18
    integer, parameter :: MIO_MESH_CAD_METRIC_MEAN = 0, MIO_MESH_CAD_METRIC_VARYING = 1
    !> Where a declared tolerance accepted a gap (v17.20.0, mio_mesh_cad_tolerance_where).
    integer, parameter :: MIO_MESH_CAD_TOLERANCE_CLOSED_START = 0, MIO_MESH_CAD_TOLERANCE_CORNER = 1, &
                          MIO_MESH_CAD_TOLERANCE_PROJECTION = 2

    !> Interop mirror of C `mio_mesh_cad_surface_opts` (104 bytes, checked once
    !> per process by `check_mesh_cad_layout`). Filled by
    !> `mio_mesh_cad_surface_opts_init` and the optional arguments.
    type, bind(c) :: mio_mesh_cad_surface_opts_t
        real(c_double) :: target_size = 0.0_c_double
        real(c_double) :: max_chordal_deviation = 0.0_c_double
        real(c_double) :: tolerance = 0.0_c_double
        integer(c_int64_t) :: min_points_per_edge = 0
        integer(c_int64_t) :: max_points = 50000000_c_int64_t
        integer(c_int32_t) :: record_face_ids = 1
        integer(c_int32_t) :: edges = 1
        integer(c_int32_t) :: match_unreferenced_trims = 1
        integer(c_int32_t) :: strict = 0
        real(c_double) :: min_angle_deg = 0.0_c_double
        integer(c_int32_t) :: smooth_iterations = 0
        integer(c_int32_t) :: smooth_method = 0
        integer(c_int32_t) :: order = 1
        integer(c_int32_t) :: pole_rings = 0
        real(c_double) :: grading = 0.0_c_double
        real(c_double) :: curvature_deviation = 0.0_c_double
        real(c_double) :: cells_across_gaps = 0.0_c_double
    end type

    !> Interop mirror of C `mio_mesh_cad_surface_opts_ex` (v17.19.0; 136 bytes
    !> with 8-byte pointers, checked by `check_mesh_cad_layout`): the metric,
    !> virtual edges, given sizes and the background, borrowed for one call;
    !> since v17.20.0 the reparametrisation threshold, model tolerances,
    !> virtual faces and corner lopping from its reserved tail.
    type, bind(c) :: mio_mesh_cad_surface_opts_ex_t
        integer(c_int32_t) :: metric = 0
        integer(c_int32_t) :: virtual_edges = 0
        type(c_ptr) :: face_size_faces = c_null_ptr
        type(c_ptr) :: face_sizes = c_null_ptr
        integer(c_int64_t) :: num_face_sizes = 0
        type(c_ptr) :: edge_size_edges = c_null_ptr
        type(c_ptr) :: edge_sizes = c_null_ptr
        integer(c_int64_t) :: num_edge_sizes = 0
        type(c_ptr) :: background = c_null_ptr
        type(c_ptr) :: background_array = c_null_ptr
        real(c_double) :: reparam_condition_number = 0
        integer(c_int32_t) :: model_tolerances = 0
        integer(c_int32_t) :: virtual_faces = 0
        integer(c_int32_t) :: lop_corners = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: reserved(5) = 0
    end type

    !> `mio_mesh_cad_surface_report` (200 bytes; field order and types are ABI,
    !> see doc/mesh_cad.md for every field). Counts are counts, not indices.
    type, bind(c) :: mio_mesh_cad_surface_report
        integer(c_int64_t) :: num_points = 0
        integer(c_int64_t) :: num_triangles = 0
        integer(c_int64_t) :: num_edge_segments = 0
        integer(c_int64_t) :: num_faces = 0
        integer(c_int64_t) :: num_failed = 0
        integer(c_int64_t) :: num_edges = 0
        integer(c_int64_t) :: num_open_edges = 0
        integer(c_int64_t) :: num_non_manifold_edges = 0
        integer(c_int64_t) :: num_misoriented_edges = 0
        integer(c_int64_t) :: num_degenerate_dropped = 0
        real(c_double) :: tolerance = 0.0_c_double
        real(c_double) :: min_angle_deg = 0.0_c_double
        real(c_double) :: max_angle_deg = 0.0_c_double
        real(c_double) :: min_edge = 0.0_c_double
        real(c_double) :: max_edge = 0.0_c_double
        real(c_double) :: max_edge_gap = 0.0_c_double
        integer(c_int64_t) :: num_warnings = 0
        integer(c_int64_t) :: num_refinement_points = 0
        integer(c_int64_t) :: num_faces_below_bound = 0
        integer(c_int64_t) :: num_smoothed = 0
        real(c_double) :: min_input_angle_deg = 0.0_c_double
        integer(c_int64_t) :: num_curved_straightened = 0
        integer(c_int64_t) :: reserved(3) = 0
    end type

    !> A corner that keeps triangles under `min_angle` (C `mio_mesh_cad_corner`,
    !> 120 bytes): two trims meeting at less than 60 degrees in the face's
    !> plane, or a pole. From `mio_mesh_cad_get_corner`, `face` and `point`
    !> are 1-based.
    type, bind(c) :: mio_mesh_cad_corner
        integer(c_int64_t) :: face = 0
        integer(c_int64_t) :: point = 0
        real(c_double) :: u = 0.0_c_double
        real(c_double) :: v = 0.0_c_double
        real(c_double) :: x(3) = 0.0_c_double
        real(c_double) :: angle_deg = 0.0_c_double
        real(c_double) :: min_forced_deg = 0.0_c_double
        real(c_double) :: reach = 0.0_c_double
        integer(c_int32_t) :: pole = 0
        integer(c_int32_t) :: reserved0 = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_hausdorff_opts`.
    type, bind(c) :: mio_hausdorff_opts_t
        type(c_ptr) :: region_a = c_null_ptr
        type(c_ptr) :: region_b = c_null_ptr
        integer(c_int64_t) :: face_samples = 0
        real(c_double) :: grid_cell_size = 0.0_c_double
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_hausdorff_report`.
    type, bind(c) :: mio_hausdorff_report_t
        real(c_double) :: distance = 0.0_c_double
        real(c_double) :: a_to_b = 0.0_c_double
        real(c_double) :: b_to_a = 0.0_c_double
        real(c_double) :: mean_a_to_b = 0.0_c_double
        real(c_double) :: rms_a_to_b = 0.0_c_double
        real(c_double) :: mean_b_to_a = 0.0_c_double
        real(c_double) :: rms_b_to_a = 0.0_c_double
        integer(c_int64_t) :: num_samples_a = 0
        integer(c_int64_t) :: num_samples_b = 0
        real(c_double) :: worst_point_a(3) = 0.0_c_double
        real(c_double) :: worst_point_b(3) = 0.0_c_double
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_region_selector` (`MIO_REGION_ANY` = -2).
    type, bind(c) :: mio_region_selector_t
        type(c_ptr) :: name = c_null_ptr
        integer(c_int32_t) :: kind = -1
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: dim = -2
        integer(c_int64_t) :: tag = -2
        integer(c_int64_t) :: reserved(2) = 0
    end type

    !> C-layout mirror of mio_find_interface_opts.
    type, bind(c) :: mio_find_interface_opts_t
        integer(c_int32_t) :: mode = 0
        integer(c_int32_t) :: master = 0
        real(c_double) :: gap_tolerance = 0.0_c_double
        real(c_double) :: angle_tolerance = 30.0_c_double
        real(c_double) :: overlap_tolerance = 0.0_c_double
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> C-layout mirror of mio_contact_pairs_opts.
    type, bind(c) :: mio_contact_pairs_opts_t
        real(c_double) :: tolerance = 0.0_c_double
        integer(c_int32_t) :: require_complete = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> C-layout mirror of mio_split_interface_opts.
    type, bind(c) :: mio_split_interface_opts_t
        integer(c_int32_t) :: add_cohesive = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_periodic_opts`; `mio_periodic_opts_init` sets
    !> the identity matrix, atol 1e-8 and require_complete.
    type, bind(c) :: mio_periodic_opts_t
        real(c_double) :: matrix(16) = 0.0_c_double
        real(c_double) :: atol = 1.0e-8_c_double
        integer(c_int32_t) :: require_complete = 1
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_agglomerate_opts` (v16.25.0).
    type, bind(c) :: mio_agglomerate_opts_t
        integer(c_int64_t) :: target_group_size = 8
        integer(c_int32_t) :: merge_coplanar_faces = 0
        integer(c_int32_t) :: reserved_pad = 0
        real(c_double) :: coplanar_angle = 1.0_c_double
        real(c_double) :: min_sphericity = 0.0_c_double
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_resample_opts` (v16.25.0).
    type, bind(c) :: mio_resample_opts_t
        type(c_ptr) :: times = c_null_ptr
        integer(c_int64_t) :: num_times = 0
        integer(c_int32_t) :: method = 0
        integer(c_int32_t) :: extrapolate = 0
        integer(c_int32_t) :: blend_points = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_quality_gate_report` (v16.24.0).
    type, bind(c) :: mio_quality_gate_report_t
        integer(c_int32_t) :: passed = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: num_checks = 0
        integer(c_int64_t) :: num_failed = 0
        integer(c_int64_t) :: num_cells = 0
        integer(c_int64_t) :: num_inverted = 0
        integer(c_int64_t) :: num_degenerate = 0
        integer(c_int64_t) :: summary_length = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_repair_opts`. Field order/types are ABI; the
    !> defaults mirror `mio_repair_opts_init` (every pass on), which is called
    !> anyway -- an all-zero struct is NOT the default.
    type, bind(c) :: mio_repair_opts_t
        integer(c_int32_t) :: fix_orientation = 1
        integer(c_int32_t) :: orient_outward = 1
        integer(c_int32_t) :: fill_holes = 1
        integer(c_int32_t) :: split_non_manifold = 1
        integer(c_int32_t) :: record_provenance = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: max_hole_edges = 10
        real(c_double) :: weld_tolerance = 0.0_c_double
        integer(c_int64_t) :: reserved(5) = 0
    end type

    !> Interop mirror of C `mio_repair_report`. Field order/types are ABI.
    type, bind(c) :: mio_repair_report_t
        type(mio_surface_quality) :: quality_before
        type(mio_surface_quality) :: quality_after
        integer(c_int64_t) :: num_flipped = 0
        integer(c_int64_t) :: num_components = 0
        integer(c_int64_t) :: largest_component = 0
        integer(c_int64_t) :: num_oriented_outward = 0
        integer(c_int64_t) :: num_unorientable = 0
        integer(c_int64_t) :: num_vertices_split = 0
        integer(c_int64_t) :: num_holes_detected = 0
        integer(c_int64_t) :: num_holes_filled = 0
        integer(c_int64_t) :: num_holes_skipped = 0
        integer(c_int64_t) :: num_faces_added = 0
        integer(c_int64_t) :: num_points_added = 0
        integer(c_int64_t) :: points_welded = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_shrinkwrap_opts`. Field order/types are ABI.
    type, bind(c) :: mio_shrinkwrap_opts_t
        type(c_ptr) :: weights = c_null_ptr
        type(c_ptr) :: target_region = c_null_ptr
        real(c_double) :: offset = 0.0_c_double
        real(c_double) :: max_distance = 0.0_c_double
        real(c_double) :: grid_cell_size = 0.0_c_double
        integer(c_int32_t) :: normal_weight = 0
        integer(c_int32_t) :: record_distance = 0
        integer(c_int32_t) :: record_closest_cell = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: reserved(5) = 0
    end type

    !> Interop mirror of C `mio_shrinkwrap_report`. Field order/types are ABI.
    type, bind(c) :: mio_shrinkwrap_report_t
        type(mio_surface_quality) :: quality
        integer(c_int64_t) :: num_projected = 0
        integer(c_int64_t) :: num_missed = 0
        integer(c_int64_t) :: num_skipped = 0
        real(c_double) :: max_displacement = 0.0_c_double
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_sobolev_opts`. Field order/types are ABI.
    type, bind(c) :: mio_sobolev_opts_t
        type(c_ptr) :: array = c_null_ptr
        type(c_ptr) :: fixed_points_array = c_null_ptr
        real(c_double) :: length_scale = 0.0_c_double
        real(c_double) :: tolerance = 1.0e-10_c_double
        integer(c_int32_t) :: max_iterations = 128
        integer(c_int32_t) :: fix_boundary = 0
        integer(c_int32_t) :: record_filtered = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: reserved(5) = 0
    end type

    !> Interop mirror of C `mio_sobolev_report`. Field order/types are ABI.
    type, bind(c) :: mio_sobolev_report_t
        integer(c_int64_t) :: num_iterations = 0
        integer(c_int64_t) :: num_fixed = 0
        integer(c_int64_t) :: num_isolated = 0
        real(c_double) :: residual = 0.0_c_double
        real(c_double) :: max_displacement = 0.0_c_double
        integer(c_int32_t) :: converged = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_sdf_opts`. Field order and types are ABI and
    !> must match bindings/c/include/meshioplusplus/meshioplusplus.h exactly;
    !> `reserved` is padding for additive growth and must stay zero.
    type, bind(c) :: mio_sdf_opts_t
        type(c_ptr) :: surface_region = c_null_ptr
        real(c_double) :: band = 0.0_c_double
        real(c_double) :: grid_cell_size = 0.0_c_double
        real(c_double) :: max_winding_work = 2.0e9_c_double
        integer(c_int32_t) :: sign = 1
        integer(c_int32_t) :: weight = 0
        integer(c_int32_t) :: location = 0
        integer(c_int32_t) :: watertight_check = 1
        integer(c_int32_t) :: record_closest_cell = 0
        integer(c_int32_t) :: record_inside = 0
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_voxel_opts`. Same ABI rules as above. Exactly
    !> one of `resolution` and `cell_size` may be given.
    type, bind(c) :: mio_voxel_opts_t
        type(c_ptr) :: resolution = c_null_ptr
        type(c_ptr) :: bounds = c_null_ptr
        real(c_double) :: cell_size = 0.0_c_double
        real(c_double) :: padding = 0.0_c_double
        real(c_double) :: padding_relative = 0.0_c_double
        integer(c_int64_t) :: max_cells = 20000000
        integer(c_int32_t) :: fill = 0
        integer(c_int32_t) :: attach_occupancy = 0
        integer(c_int32_t) :: sign = 1
        integer(c_int32_t) :: watertight_check = 1
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_gen_opts`, the options of the parametric
    !> generators. Same ABI rules as above: 160 bytes, checked once per process
    !> (`check_gen_opts_layout`).
    type, bind(c) :: mio_gen_opts_t
        integer(c_int64_t) :: counts(3) = 8
        real(c_double) :: lengths(3) = 1.0_c_double
        real(c_double) :: origin(3) = 0.0_c_double
        real(c_double) :: radius = 1.0_c_double
        real(c_double) :: height = 1.0_c_double
        integer(c_int64_t) :: max_cells = 20000000
        integer(c_int32_t) :: cell_type = 0
        integer(c_int32_t) :: split = 0
        integer(c_int32_t) :: order = 1
        integer(c_int32_t) :: pad = 0
        real(c_double) :: radius2 = 0.5_c_double
        real(c_double) :: angle = 6.283185307179586_c_double
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_octree_opts`. Same ABI rules as above.
    type, bind(c) :: mio_octree_opts_t
        type(c_ptr) :: root_cells = c_null_ptr
        type(c_ptr) :: bounds = c_null_ptr
        real(c_double) :: padding = 0.0_c_double
        real(c_double) :: padding_relative = 0.0_c_double
        real(c_double) :: sample_distance = 0.0_c_double
        integer(c_int64_t) :: min_level = 0
        integer(c_int64_t) :: max_level = 4
        integer(c_int64_t) :: dimension = 0
        integer(c_int64_t) :: max_cells = 20000000
        integer(c_int32_t) :: balance = 1
        integer(c_int32_t) :: fill = 0
        integer(c_int32_t) :: attach_distance = 0
        integer(c_int32_t) :: even_odd = 0
        integer(c_int32_t) :: sign = 1
        integer(c_int32_t) :: watertight_check = 1
        type(c_ptr) :: size_array = c_null_ptr
        real(c_double) :: curvature_threshold = 0.0_c_double
        real(c_double) :: feature_angle = 0.0_c_double
        integer(c_int64_t) :: surface_level = -1
        integer(c_int32_t) :: transition = 0
        integer(c_int32_t) :: snap = 0
        integer(c_int32_t) :: plane = 0
        integer(c_int32_t) :: reserved32 = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_octree_report`.
    type, bind(c) :: mio_octree_report_t
        integer(c_int64_t) :: dimension = 0
        integer(c_int64_t) :: num_leaves = 0
        integer(c_int64_t) :: num_hanging = 0
        integer(c_int64_t) :: finest_level = 0
        real(c_double) :: origin(3) = 0.0_c_double
        real(c_double) :: root_size = 0.0_c_double
        integer(c_int64_t) :: num_snapped = 0
        integer(c_int64_t) :: num_refused = 0
        integer(c_int64_t) :: reserved(2) = 0
    end type

    !> Interop mirror of C `mio_octree_mesh_opts`. Same ABI rules as above: 264
    !> bytes, checked once per process (`check_octree_mesh_layout`). Always
    !> filled by `mio_octree_mesh_opts_init`, so the defaults come from C.
    type, bind(c) :: mio_octree_mesh_opts_t
        type(c_ptr) :: bounds = c_null_ptr
        type(c_ptr) :: size_array = c_null_ptr
        type(c_ptr) :: surface_region = c_null_ptr
        real(c_double) :: min_size = 0.0_c_double
        real(c_double) :: max_size = 0.0_c_double
        real(c_double) :: padding = 0.0_c_double
        real(c_double) :: padding_relative = 0.1_c_double
        real(c_double) :: curvature_angle = 15.0_c_double
        real(c_double) :: feature_angle = 30.0_c_double
        real(c_double) :: corner_angle = 30.0_c_double
        real(c_double) :: feature_size = 0.0_c_double
        real(c_double) :: gap_cells = 2.0_c_double
        real(c_double) :: surface_size = 0.0_c_double
        real(c_double) :: winding_threshold = 0.5_c_double
        real(c_double) :: winding_beta = 2.0_c_double
        real(c_double) :: alpha_long = 0.24799_c_double
        real(c_double) :: alpha_short = 0.40173_c_double
        integer(c_int64_t) :: max_level = -1
        integer(c_int64_t) :: grading_layers = 1
        integer(c_int64_t) :: fit_iterations = 5
        integer(c_int64_t) :: untangle_iterations = 50
        integer(c_int64_t) :: max_cells = 20000000
        integer(c_int64_t) :: max_tets = 100000000
        integer(c_int32_t) :: output = 1
        integer(c_int32_t) :: conformity = 0
        integer(c_int32_t) :: skin_grading = 0
        integer(c_int32_t) :: fill = 1
        integer(c_int32_t) :: gap_outside = 0
        integer(c_int32_t) :: preserve_features = 1
        integer(c_int32_t) :: dihedral_floor_guard = 1
        integer(c_int32_t) :: keep_outside = 0
        integer(c_int32_t) :: attach_distance = 1
        integer(c_int32_t) :: attach_quality = 0
        integer(c_int32_t) :: pillow = 1
        integer(c_int32_t) :: reserved32 = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> What `m%octree_mesh` reports: an interop mirror of C
    !> `mio_octree_mesh_report` (392 bytes; field order and types are ABI, see
    !> doc/octree_mesh.md for every field). `output` is 0 (castellated), 1 (tet)
    !> or 2 (hex); `guarantee_applies` is nonzero when the 10.7-164.8 degree
    !> bound covers every tetrahedron. Counts are counts, not indices: nothing
    !> here is shifted to 1-based.
    type, bind(c) :: mio_octree_mesh_report
        integer(c_int32_t) :: output = 0
        integer(c_int32_t) :: guarantee_applies = 0
        real(c_double) :: origin(3) = 0.0_c_double
        real(c_double) :: root_size = 0.0_c_double
        integer(c_int64_t) :: root_cells(3) = 0
        integer(c_int64_t) :: max_level = 0
        integer(c_int64_t) :: finest_level = 0
        integer(c_int64_t) :: surface_level = -1
        real(c_double) :: min_size = 0.0_c_double
        real(c_double) :: max_size = 0.0_c_double
        integer(c_int64_t) :: num_leaves = 0
        integer(c_int64_t) :: num_points = 0
        integer(c_int64_t) :: num_cells = 0
        integer(c_int64_t) :: num_hexahedra = 0
        integer(c_int64_t) :: num_tetrahedra = 0
        integer(c_int64_t) :: num_inside = 0
        integer(c_int64_t) :: num_boundary_cells = 0
        integer(c_int64_t) :: num_surrogate_faces = 0
        integer(c_int64_t) :: num_hanging = 0
        integer(c_int64_t) :: num_cut_points = 0
        integer(c_int64_t) :: num_warped = 0
        integer(c_int64_t) :: num_warp_refused = 0
        real(c_double) :: min_dihedral = 0.0_c_double
        real(c_double) :: max_dihedral = 0.0_c_double
        real(c_double) :: min_scaled_jacobian = 0.0_c_double
        real(c_double) :: hausdorff_to_skin = -1.0_c_double
        integer(c_int64_t) :: num_inverted = 0
        integer(c_int64_t) :: num_boundary_edges_open = 0
        integer(c_int64_t) :: num_non_manifold_boundary_edges = 0
        integer(c_int64_t) :: num_corners = 0
        integer(c_int64_t) :: num_corners_snapped = 0
        integer(c_int64_t) :: num_feature_chains = 0
        integer(c_int64_t) :: num_feature_chains_missed = 0
        integer(c_int64_t) :: num_feature_moves = 0
        integer(c_int64_t) :: num_optimized_moves = 0
        integer(c_int64_t) :: num_flips = 0
        integer(c_int64_t) :: num_winding_rechecks = 0
        integer(c_int64_t) :: num_warnings = 0
        integer(c_int64_t) :: reserved(8) = 0
    end type

    !> Interop mirror of C `mio_compute_sdf_opts`. Same ABI rules as above. For
    !> `structure = 0` (voxel) exactly one of `resolution` and `cell_size` may be
    !> given; for `structure = 1` (octree) NEITHER may be -- its finest cell is
    !> root/2**depth and is therefore already determined.
    type, bind(c) :: mio_compute_sdf_opts_t
        type(c_ptr) :: resolution = c_null_ptr
        type(c_ptr) :: bounds = c_null_ptr
        real(c_double) :: cell_size = 0.0_c_double
        real(c_double) :: padding = 0.0_c_double
        real(c_double) :: padding_relative = 0.1_c_double
        real(c_double) :: band_cells = 1.0_c_double
        integer(c_int64_t) :: max_cells = 20000000
        integer(c_int64_t) :: root_resolution = 8
        integer(c_int64_t) :: max_depth = 4
        integer(c_int32_t) :: structure = 0
        integer(c_int32_t) :: record_levels = 1
        integer(c_int64_t) :: reserved(6) = 0
        type(mio_sdf_opts_t) :: distance
    end type

    ! Geometric statistics (bind(c); layout must match mio_stats_report in
    ! meshioplusplus.h). Per-cell-type counts are not carried across the C ABI.
    type, bind(c) :: mio_stats_report
        integer(c_int64_t) :: num_points
        integer(c_int64_t) :: num_cells
        real(c_double) :: bbox_min(3)
        real(c_double) :: bbox_max(3)
        real(c_double) :: extent(3)
        real(c_double) :: centroid(3)
        real(c_double) :: total_area
        real(c_double) :: signed_volume
        real(c_double) :: unsigned_volume
        integer(c_int64_t) :: num_inverted
    end type

    ! Summary of one data array (bind(c); layout must match
    ! mio_data_array_info in meshioplusplus.h). Per-component statistics are
    ! retrieved separately via the `data_info` procedure's optional arguments.
    type, bind(c) :: mio_data_array_info
        integer(c_int) :: location            !< a MIO_DATA_* value
        integer(c_int) :: dtype               !< a MIO_* dtype, as stored
        integer(c_int64_t) :: num_blocks      !< cell_data: cell blocks; else 1
        integer(c_int64_t) :: num_entries     !< rows (points / cells / length)
        integer(c_int64_t) :: num_components  !< product of trailing dimensions
        integer(c_int64_t) :: num_values      !< num_entries * num_components
        real(c_double) :: min                 !< over finite values, else NaN
        real(c_double) :: max                 !< over finite values, else NaN
        real(c_double) :: mean                !< over finite values, else NaN
        integer(c_int64_t) :: num_nan
        integer(c_int64_t) :: num_inf
        integer(c_int64_t) :: num_finite
        integer(c_int) :: inconsistent_blocks !< nonzero if blocks disagree
    end type

    ! Whole-mesh (or one named Cell region's) reduction of one array (bind(c);
    ! layout must match mio_field_integral_info in meshioplusplus.h). Per-
    ! component total/mean/domain_measure/num_nan are retrieved separately via
    ! the `data_integrate`/`data_integrate_region` procedures' optional
    ! arguments.
    type, bind(c) :: mio_field_integral_info
        integer(c_int64_t) :: num_components !< product of trailing dimensions
        integer(c_int64_t) :: num_cells      !< cells with a computable measure
        integer(c_int64_t) :: num_skipped    !< cells excluded: unmeasurable geometry
    end type

    ! mio_dtype values (must match the C enum in meshioplusplus.h).
    integer(c_int), parameter :: MIO_FLOAT32 = 0, MIO_FLOAT64 = 1
    integer(c_int), parameter :: MIO_INT32 = 4, MIO_INT64 = 5

    ! Data-operation enumerations (must match the C enums in meshioplusplus.h).
    integer(c_int), parameter :: MIO_DATA_POINT = 0, MIO_DATA_CELL = 1, MIO_DATA_FIELD = 2
    integer(c_int), parameter :: MIO_WEIGHT_UNIFORM = 0, MIO_WEIGHT_MEASURE = 1
    integer(c_int), parameter :: MIO_COND_CLAMP = 0, MIO_COND_NORMALIZE = 1
    !> Provenance modes (see doc/provenance.md). OFF is the default
    !> everywhere: a writer emits only the one-line credit.
    integer(c_int), parameter :: MIO_PROVENANCE_OFF = 0
    integer(c_int), parameter :: MIO_PROVENANCE_BEST_EFFORT = 1
    integer(c_int), parameter :: MIO_PROVENANCE_REQUIRED = 2
    integer(c_int), parameter :: MIO_COND_STANDARDIZE = 2
    integer(c_int), parameter :: MIO_SCOPE_COMPONENT = 0, MIO_SCOPE_MAGNITUDE = 1
    integer(c_int), parameter :: MIO_NAN_IGNORE = 0, MIO_NAN_REPLACE = 1, MIO_NAN_FAIL = 2
    integer(c_int), parameter :: MIO_TINV_MISES = 1, MIO_TINV_PRINCIPAL = 2
    integer(c_int), parameter :: MIO_TINV_HYDROSTATIC = 4, MIO_TINV_DEVIATORIC = 8
    integer(c_int), parameter :: MIO_TINV_ALL = 15

    ! Region kinds (must match the C enum mio_region_kind).
    integer(c_int), parameter :: MIO_REGION_POINT = 0, MIO_REGION_CELL = 1
    integer(c_int), parameter :: MIO_REGION_SIDE = 2

    ! MDPA side-channel sections and value kinds (must match the C enums
    ! mio_mdpa_section / mio_mdpa_value_kind; see mio_format_info).
    integer(c_int32_t), parameter, public :: MIO_GMSH_BOUNDING_ENTITIES = 0, MIO_GMSH_PERIODIC = 1
    !> NURBS side-channel sections (an 'mfem' handle; see meshioplusplus.h).
    integer(c_int32_t), parameter, public :: MIO_NURBS_PATCHES = 0, MIO_NURBS_FIELDS = 1, &
                                             MIO_NURBS_TOPOLOGY = 2, MIO_NURBS_CELL_MAP = 3
    !> B-rep side-channel sections (a 'cad_json' handle, or 'step'/'iges' with
    !> OpenCASCADE; see meshioplusplus.h). Curve fields are 10 and up.
    integer(c_int32_t), parameter, public :: MIO_BREP_FILE = 0, MIO_BREP_BREPS = 1, &
                                             MIO_BREP_PATCHES = 2, MIO_BREP_FACES = 3, &
                                             MIO_BREP_LOOPS = 4, MIO_BREP_TRIMS = 5, &
                                             MIO_BREP_EDGES = 6, MIO_BREP_VERTICES = 7, &
                                             MIO_BREP_REFS = 8
    !> What `mio_iga_evaluation` holds (the C enum mio_iga_array).
    integer(c_int32_t), parameter, public :: MIO_IGA_POINTS = 0, MIO_IGA_DERIVATIVES = 1, &
                                             MIO_IGA_NORMALS = 2, MIO_IGA_INSIDE = 3, &
                                             MIO_IGA_FIELD = 4
    integer(c_int), parameter :: MIO_MDPA_PROPERTIES = 0, MIO_MDPA_ENTITY_NAMES = 1
    integer(c_int), parameter :: MIO_MDPA_SKIPPED = 2, MIO_MDPA_MODEL_PART_DATA = 3
    integer(c_int), parameter :: MIO_MDPA_TABLES = 4, MIO_MDPA_GEOMETRIES = 5
    integer(c_int), parameter :: MIO_MDPA_MESH_BLOCKS = 6, MIO_MDPA_SUBMODELPARTS = 7
    integer(c_int), parameter :: MIO_MDPA_RAW_BLOCKS = 8
    integer(c_int), parameter :: MIO_MDPA_VALUE_NUMBER = 0, MIO_MDPA_VALUE_TEXT = 1
    integer(c_int), parameter :: MIO_MDPA_VALUE_TABLE = 2

    ! Description of one named region (bind(c); layout must match
    ! mio_region_info in meshioplusplus.h). The entries themselves come back
    ! from the `regions` procedure as a separate allocatable array.
    type, bind(c) :: mio_region_info
        integer(c_int) :: kind               !< a MIO_REGION_* value
        integer(c_int) :: dim                !< topological dimension, or -1
        integer(c_int64_t) :: tag            !< format-native id, or -1
        integer(c_int64_t) :: num_entries    !< grouped entities (rows)
        integer(c_int64_t) :: stride         !< 2 for SIDE, else 1
    end type

    ! One cell block's full shape (bind(c); layout must match
    ! mio_cell_block_info in meshioplusplus.h). Deliberately NOT named
    ! mio_cell_block_info -- that name is already this module's own metadata
    ! summary type, which is a different thing (it carries the type NAME and
    ! comes from a file summary, not from a live mesh). Private, like every
    ! other c_mio_* interop detail; the public surface is the scalar
    ! cell_block_* queries below.
    type, bind(c) :: c_mio_cell_block_info
        integer(c_int64_t) :: num_cells
        integer(c_int64_t) :: nodes_per_cell   !< 0 when ragged
        integer(c_int32_t) :: is_ragged
        integer(c_int32_t) :: is_polyhedron    !< implies is_ragged
        integer(c_int64_t) :: num_faces
        integer(c_int64_t) :: num_nodes
        integer(c_int64_t) :: reserved(6)
    end type

    ! One ragged block's CSR shape (bind(c); matches mio_poly_conn_shape).
    type, bind(c) :: c_mio_poly_conn_shape
        integer(c_int32_t) :: is_polyhedron
        integer(c_int32_t) :: reserved0
        integer(c_int64_t) :: num_cells
        integer(c_int64_t) :: num_faces
        integer(c_int64_t) :: num_nodes
        integer(c_int64_t) :: reserved(4)
    end type

    integer, parameter :: MIO_MAX_NDIM = 8
    integer, parameter :: STRBUF_LEN = 4096

    !> Interop mirror of C `mio_read_opts`. Field order and types are ABI and
    !> must match bindings/c/include/meshioplusplus/meshioplusplus.h exactly;
    !> `reserved` is padding for additive growth and must stay zero.
    type, bind(c) :: mio_read_opts_t
        integer(c_int) :: points_only = 0
        integer(c_int) :: metadata_only = 0
        type(c_ptr) :: arrays = c_null_ptr
        integer(c_int64_t) :: num_arrays = 0
        integer(c_int) :: mmap_mode = 0
        !> Which step of a multi-step file to materialize: 0 (default) is the
        !> first, negative counts from the end. Takes one of the former
        !> `reserved` slots, so the struct size and every preceding field's
        !> offset are unchanged.
        integer(c_int64_t) :: time_step = 0
        !> Nonzero downgrades "this reader cannot represent construct X" errors
        !> to a warning plus a skip (currently mdpa's Table/Geometries/Mesh/
        !> Constraints blocks). Not "ignore all errors": a malformed file still
        !> fails. Takes a second former `reserved` slot; size unchanged.
        integer(c_int64_t) :: lenient = 0
        !> Which piece of a partitioned file to keep (VTKHDF partitions or composite
        !> blocks); meaningful only when `piece_set` is nonzero. 0 is the first,
        !> negative counts from the end. Takes a third former `reserved` slot; size
        !> unchanged.
        integer(c_int64_t) :: piece = 0
        !> Nonzero means `piece` was chosen; 0 (the default) merges every piece into
        !> one mesh with one cell region per piece. Takes a fourth former `reserved`
        !> slot; size unchanged.
        integer(c_int64_t) :: piece_set = 0
        !> Nonzero removes the ghost cells (halo) of a partitioned .pvtu/.pvtp/.pvd
        !> (every cell with a vtkGhostType bit set, and the points only they used);
        !> 0 keeps them. Every other reader ignores it. Takes a fifth former
        !> `reserved` slot; size unchanged.
        integer(c_int64_t) :: drop_ghosts = 0
        integer(c_int64_t) :: reserved(1) = 0
    end type

    !> Interop mirror of C `mio_refine_opts`. Field order and types are ABI and
    !> must match bindings/c/include/meshioplusplus/meshioplusplus.h exactly;
    !> `reserved` is padding for additive growth and must stay zero.
    !>
    !> All zero (plus `levels = 1`) means "refine every cell once". At most ONE
    !> of `cells`, `region` and `predicate_array` may be given.
    type, bind(c) :: mio_refine_opts_t
        type(c_ptr) :: cells = c_null_ptr
        integer(c_int64_t) :: num_cells = 0
        type(c_ptr) :: region = c_null_ptr
        type(c_ptr) :: predicate_array = c_null_ptr
        real(c_double) :: predicate_value = 0.0_c_double
        integer(c_int32_t) :: levels = 1
        integer(c_int32_t) :: record_parent_ids = 0
        integer(c_int32_t) :: record_levels = 0
        integer(c_int32_t) :: closure = 0
        integer(c_int32_t) :: predicate_op = 0
        integer(c_int32_t) :: reserved_pad = 0
        !> Nonzero to attach refine:cell_id/refine:parent_id -- the persistent
        !> parent/child hierarchy a multigrid caller resolves across the
        !> sequence of meshes it keeps. Also forces refine:entity to be
        !> attached even when the closure leaves no hanging node. Takes one of
        !> the former `reserved` slots; size unchanged.
        integer(c_int64_t) :: record_hierarchy = 0
        integer(c_int64_t) :: reserved(5) = 0
    end type

    !> Interop mirror of C `mio_tessellate_opts`. Field order/types are ABI; the
    !> defaults mirror `mio_tessellate_opts_init` (called anyway).
    type, bind(c) :: mio_tessellate_opts_t
        integer(c_int64_t) :: levels = 2
        integer(c_int32_t) :: curved = 1
        integer(c_int32_t) :: fields = 1
        integer(c_int32_t) :: record_stencil = 0
        integer(c_int32_t) :: reserved0 = 0
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_tessellate_report`.
    type, bind(c) :: mio_tessellate_report_t
        integer(c_int64_t) :: num_source_points = 0
        integer(c_int64_t) :: num_source_cells = 0
        integer(c_int64_t) :: num_curved_source_cells = 0
        integer(c_int64_t) :: num_pass_through_source_cells = 0
        integer(c_int64_t) :: num_points = 0
        integer(c_int64_t) :: num_cells = 0
        integer(c_int64_t) :: num_boundary_facets = 0
        integer(c_int64_t) :: num_interior_facets = 0
        integer(c_int64_t) :: num_facets_with_bad_count = 0
        integer(c_int64_t) :: reserved(5) = 0
    end type

    !> Interop mirror of C `mio_smooth_opts`. Field order and types are ABI
    !> and must match bindings/c/include/meshioplusplus/meshioplusplus.h
    !> exactly; `reserved_pad0`/`reserved_pad1`/`reserved` are padding for
    !> additive growth and must stay zero. `method`/`frozen` are a C string
    !> and an array, so they need the same TARGET-buffer + c_loc idiom
    !> `region`/`cells` use in mio_refine_opts_t -- see mesh_smooth below.
    type, bind(c) :: mio_smooth_opts_t
        type(c_ptr) :: method = c_null_ptr
        integer(c_int32_t) :: iterations = 10
        integer(c_int32_t) :: reserved_pad0 = 0
        real(c_double) :: lambda = -1.0_c_double
        real(c_double) :: mu = -0.34_c_double
        integer(c_int32_t) :: fix_boundary = 1
        integer(c_int32_t) :: preserve_features = 1
        real(c_double) :: feature_angle = 30.0_c_double
        integer(c_int32_t) :: guard_inversion = 1
        integer(c_int32_t) :: reserved_pad1 = 0
        type(c_ptr) :: frozen = c_null_ptr
        integer(c_int64_t) :: num_frozen = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_decimate_opts`. Field order and types are ABI
    !> and must match bindings/c/include/meshioplusplus/meshioplusplus.h
    !> exactly; `reserved` is padding for additive growth and must stay zero.
    type, bind(c) :: mio_decimate_opts_t
        real(c_double) :: target_ratio = -1.0_c_double
        integer(c_int64_t) :: target_faces = -1_c_int64_t
        real(c_double) :: max_error = -1.0_c_double
        type(c_ptr) :: placement = c_null_ptr
        integer(c_int32_t) :: preserve_boundary = 1
        integer(c_int32_t) :: preserve_features = 1
        real(c_double) :: feature_angle = 30.0_c_double
        type(c_ptr) :: frozen = c_null_ptr
        integer(c_int64_t) :: num_frozen = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_remesh_opts`. Field order and types are ABI
    !> and must match bindings/c/include/meshioplusplus/meshioplusplus.h
    !> exactly; `reserved`/`reserved_d` are padding for additive growth and
    !> must stay zero. `metric` is a C string, so it needs the same
    !> TARGET-buffer + c_loc idiom `region`/`predicate_array` use in
    !> mio_refine_opts_t -- see mesh_remesh below.
    type, bind(c) :: mio_remesh_opts_t
        integer(c_int64_t) :: num_clusters = 0
        integer(c_int32_t) :: subdivide = -1
        integer(c_int32_t) :: max_subdivide = 4
        real(c_double) :: subsample_ratio = 10.0_c_double
        integer(c_int32_t) :: max_iterations = 100
        integer(c_int32_t) :: max_repair_passes = 10
        type(c_ptr) :: metric = c_null_ptr
        real(c_double) :: gradation = 0.0_c_double
        integer(c_int32_t) :: preserve_boundary = 1
        integer(c_int32_t) :: reserved_pad = 0
        real(c_double) :: max_anisotropy = 4.0_c_double
        real(c_double) :: reserved_d(3) = 0.0_c_double
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_remesh_report`. Field order/types are ABI.
    type, bind(c) :: mio_remesh_report_t
        integer(c_int64_t) :: num_clusters = 0
        integer(c_int64_t) :: num_iterations = 0
        integer(c_int32_t) :: subdivide_applied = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: num_isolated_clusters = 0
        integer(c_int64_t) :: num_non_manifold_vertices = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_remesh_volume_opts`. Field order and types are
    !> ABI and must match bindings/c/include/meshioplusplus/meshioplusplus.h
    !> exactly; `reserved` is padding for additive growth and must stay zero.
    !> `distance` is embedded by value, the `mio_compute_sdf_opts_t` precedent
    !> above -- its own `reserved` tail absorbs future growth without shifting
    !> anything here.
    type, bind(c) :: mio_remesh_volume_opts_t
        type(c_ptr) :: resolution = c_null_ptr
        type(c_ptr) :: bounds = c_null_ptr
        real(c_double) :: cell_size = 0.0_c_double
        real(c_double) :: padding = 0.0_c_double
        real(c_double) :: padding_relative = 0.1_c_double
        integer(c_int64_t) :: max_cells = 20000000
        integer(c_int64_t) :: max_tets = 20000000
        real(c_double) :: warp_fraction = 0.35_c_double
        integer(c_int64_t) :: reserved(6) = 0
        type(mio_sdf_opts_t) :: distance
    end type

    !> Interop mirror of C `mio_remesh_volume_report`. Field order/types ABI.
    type, bind(c) :: mio_remesh_volume_report_t
        type(mio_surface_quality) :: input_quality
        integer(c_int64_t) :: num_tets = 0
        integer(c_int64_t) :: num_vertices_warped = 0
        integer(c_int64_t) :: num_tets_rejected = 0
        integer(c_int64_t) :: num_non_manifold_edges = 0
        integer(c_int64_t) :: reserved(4) = 0
    end type

    !> Interop mirror of C `mio_xdmf_series_opts`. Field order and types are ABI
    !> and must match bindings/c/include/meshioplusplus/meshioplusplus.h exactly;
    !> `reserved` is padding for additive growth and must stay zero.
    type, bind(c) :: mio_xdmf_series_opts_t
        type(c_ptr) :: data_format = c_null_ptr
        integer(c_int32_t) :: gzip_level = -1
        integer(c_int32_t) :: mode = 0
        integer(c_int32_t) :: auto_flush = 0
        integer(c_int32_t) :: reserved_pad = 0
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> One named region's shape, without its entries (see `mio_metadata%regions`).
    type :: mio_region_summary
        character(len=STRBUF_LEN) :: name = ''
        integer(c_int) :: kind = MIO_REGION_POINT
        integer(c_int) :: dim = -1
        integer(c_int64_t) :: tag = -1
        integer(c_int64_t) :: num_entries = 0
    end type

    !> One cell block's shape, without its connectivity.
    type :: mio_cell_block_info
        character(len=STRBUF_LEN) :: cell_type = ''
        integer(c_int64_t) :: num_cells = 0
        integer(c_int64_t) :: nodes_per_cell = 0
        logical :: ragged = .false.
    end type

    !> A file summary produced without loading the heavy arrays.
    !> `has_bbox` is false for a native summary, which never decodes the point
    !> coordinates; `fell_back_to_full_read` says whether the answer was cheap.
    type :: mio_metadata
        integer(c_int64_t) :: num_points = 0
        integer(c_int64_t) :: point_dim = 0
        integer(c_int64_t) :: num_cells = 0
        type(mio_cell_block_info), allocatable :: cell_blocks(:)
        character(len=STRBUF_LEN), allocatable :: point_data_names(:)
        character(len=STRBUF_LEN), allocatable :: cell_data_names(:)
        character(len=STRBUF_LEN), allocatable :: field_data_names(:)
        logical :: has_bbox = .false.
        logical :: fell_back_to_full_read = .false.
        !> The file's recorded time-series values; size 0 for a format with no
        !> time concept. This is the count `time_step` may name.
        real(c_double), allocatable :: time_values(:)
        !> The file's named regions, without their entries. Populated whenever
        !> the summary came from an already-read mesh; size 0 on a native
        !> metadata path, since none of those formats currently map regions.
        type(mio_region_summary), allocatable :: regions(:)
    end type

    !> Interop mirror of C `mio_render_opts`. Field order/types are ABI; use the
    !> pure-Fortran `mio_render_settings` instead of filling this by hand.
    type, bind(c) :: mio_render_opts_t
        type(c_ptr) :: view = c_null_ptr
        type(c_ptr) :: color_by = c_null_ptr
        type(c_ptr) :: cmap = c_null_ptr
        type(c_ptr) :: reduce = c_null_ptr
        type(c_ptr) :: expr = c_null_ptr
        type(c_ptr) :: vectors = c_null_ptr
        type(c_ptr) :: warp = c_null_ptr
        type(c_ptr) :: quality_metric = c_null_ptr
        type(c_ptr) :: iso_levels = c_null_ptr
        integer(c_int64_t) :: num_iso_levels = 0
        real(c_double) :: pixel_aspect = 1.0_c_double
        real(c_double) :: azimuth = 45.0_c_double
        real(c_double) :: elevation = 35.264389682754654_c_double
        real(c_double) :: roll = 0.0_c_double
        real(c_double) :: fov_deg = 30.0_c_double
        real(c_double) :: zoom = 1.0_c_double
        real(c_double) :: pan_x = 0.0_c_double
        real(c_double) :: pan_y = 0.0_c_double
        real(c_double) :: ambient = 0.25_c_double
        real(c_double) :: light_x = 0.0_c_double
        real(c_double) :: light_y = 0.0_c_double
        real(c_double) :: light_z = 1.0_c_double
        real(c_double) :: split_angle = 30.0_c_double
        real(c_double) :: feature_angle = 30.0_c_double
        real(c_double) :: point_radius = 1.5_c_double
        real(c_double) :: vmin = 0.0_c_double
        real(c_double) :: vmax = 0.0_c_double
        real(c_double) :: clip_low = 0.0_c_double
        real(c_double) :: clip_high = 0.0_c_double
        real(c_double) :: scale_threshold = 1.0_c_double
        real(c_double) :: vector_length = 0.0_c_double
        real(c_double) :: warp_scale = 1.0_c_double
        integer(c_int32_t) :: width = 320
        integer(c_int32_t) :: height = 240
        integer(c_int32_t) :: supersample = 1
        integer(c_int32_t) :: perspective = 0
        integer(c_int32_t) :: shading = 1
        integer(c_int32_t) :: two_sided = 1
        integer(c_int32_t) :: edges = 0
        integer(c_int32_t) :: component = -1
        integer(c_int32_t) :: has_vmin = 0
        integer(c_int32_t) :: has_vmax = 0
        integer(c_int32_t) :: has_clip_low = 0
        integer(c_int32_t) :: has_clip_high = 0
        integer(c_int32_t) :: colorbar = 0
        integer(c_int32_t) :: axes = 0
        integer(c_int32_t) :: scale_bar = 0
        integer(c_int32_t) :: symmetric = 0
        integer(c_int32_t) :: scale = 0
        integer(c_int32_t) :: categorical = 0
        integer(c_int32_t) :: color_regions = 0
        integer(c_int32_t) :: category_edges = 0
        integer(c_int32_t) :: isolines = 0
        integer(c_int32_t) :: vector_count = 200
        integer(c_int32_t) :: warp_outline = 0
        integer(c_int32_t) :: diagnostic = 0
        integer(c_int32_t) :: edge_color = 255
        integer(c_int32_t) :: fill_color = -926564865
        integer(c_int32_t) :: line_color = 33023
        integer(c_int32_t) :: background = 0
        integer(c_int32_t) :: nan_color = -2139062017
        integer(c_int32_t) :: iso_color = 505290495
        integer(c_int32_t) :: vector_color = -600689921
        integer(c_int32_t) :: outline_color = -1768515841
        type(c_ptr) :: cutaways = c_null_ptr
        integer(c_int32_t) :: num_cutaways = 0
        integer(c_int32_t) :: cutaway_tint = -394251009
        type(c_ptr) :: streamlines = c_null_ptr
        integer(c_int32_t) :: stream_seeds = 40
        integer(c_int32_t) :: stream_color = -263151361
        real(c_double) :: stream_length = 0.5_c_double
        integer(c_int64_t) :: reserved(1) = 0
    end type

    !> Interop mirror of C `mio_text_opts`.
    type, bind(c) :: mio_text_opts_t
        real(c_double) :: cell_aspect = 2.0_c_double
        integer(c_int32_t) :: encoding = 0
        integer(c_int32_t) :: color_depth = 0
        integer(c_int32_t) :: format = 0
        integer(c_int32_t) :: cols = 80
        integer(c_int32_t) :: rows = 24
        integer(c_int32_t) :: cell_pixel_width = 8
        integer(c_int32_t) :: cell_pixel_height = 16
        integer(c_int32_t) :: tmux = 0
        integer(c_int32_t) :: notes = 1
        integer(c_int32_t) :: flags = 0
        integer(c_int64_t) :: reserved(6) = 0
    end type

    !> Interop mirror of C `mio_snapshot_opts`.
    type, bind(c) :: mio_snapshot_opts_t
        real(c_double) :: cast_fps = 12.0_c_double
        real(c_double) :: cast_degrees = 360.0_c_double
        integer(c_int32_t) :: png_compress = 0
        integer(c_int32_t) :: cast_frames = 36
        integer(c_int64_t) :: reserved(6) = 0
    end type

    type :: mio_mesh
        private
        type(c_ptr) :: handle = c_null_ptr
    contains
        procedure :: create => mesh_create
        procedure :: free => mesh_free
        procedure :: is_valid => mesh_is_valid
        procedure :: read => mesh_read
        procedure :: write => mesh_write
        procedure :: write_gltf => mesh_write_gltf
        procedure :: read_with_info => mesh_read_with_info
        procedure :: write_with_info => mesh_write_with_info
        ! -- operations --
        procedure :: extract_surface => mesh_extract_surface
        procedure :: extract_skin => mesh_extract_skin
        procedure :: attach_quality => mesh_attach_quality
        procedure :: quality_counts => mesh_quality_counts
        procedure :: reorder => mesh_reorder
        procedure :: transform => mesh_transform
        procedure :: clean => mesh_clean
        procedure :: smooth => mesh_smooth
        procedure :: optimize_volume => mesh_optimize_volume
        procedure :: crop_bbox => mesh_crop_bbox
        procedure :: crop_plane => mesh_crop_plane
        procedure :: crop_predicate => mesh_crop_predicate
        procedure :: slice => mesh_slice
        procedure :: isosurface => mesh_isosurface
        procedure :: gradient => mesh_gradient
        procedure :: hessian => mesh_hessian
        procedure :: estimate_error => mesh_estimate_error
        procedure :: remesh => mesh_remesh
        procedure :: curvature => mesh_curvature
        procedure :: normals => mesh_normals
        procedure :: feature_edges => mesh_feature_edges
        procedure :: segment => mesh_segment
        procedure :: fit_primitives => mesh_fit_primitives
        procedure :: refine_patches => mesh_refine_patches
        procedure :: fit_cad => mesh_fit_cad
        procedure :: render => mesh_render
        procedure :: render_text => mesh_render_text
        procedure :: write_snapshot => mesh_write_snapshot
        procedure :: hausdorff_distance => mesh_hausdorff_distance
        procedure :: edit_regions => mesh_edit_regions
        procedure :: remove_region => mesh_remove_region
        procedure :: match_periodic_nodes => mesh_match_periodic_nodes
        procedure :: check_quality => mesh_check_quality
        procedure :: repair => mesh_repair
        procedure :: sobolev_deform => mesh_sobolev_deform
        procedure :: remesh_volume => mesh_remesh_volume
        procedure :: split => mesh_split
        procedure :: convert_cells => mesh_convert_cells
        procedure :: subdivide => mesh_subdivide
        procedure :: agglomerate => mesh_agglomerate
        procedure :: refine => mesh_refine
        procedure :: tessellate => mesh_tessellate
        procedure :: decimate => mesh_decimate
        procedure :: partition => mesh_partition
        procedure :: partition_labels => mesh_partition_labels
        procedure :: stats => mesh_stats
        procedure :: voxelize => mesh_voxelize
        procedure :: octree => mesh_octree
        procedure :: octree_mesh => mesh_octree_mesh
        procedure :: compute_sdf => mesh_compute_sdf
        procedure :: watertight_check => mesh_watertight_check
        procedure :: sample_distance => mesh_sample_distance
        procedure :: distance_to_surface => mesh_distance_to_surface
        !> Named regions (doc/regions.md): the groups a set-capable format
        !> carries. `regions` returns one mio_region_info per group, with the
        !> names in `keys` and the flat int64 entries in `entries`.
        procedure :: regions => mesh_regions
        procedure :: region_adjacency => mesh_region_adjacency
        procedure :: find_interface => mesh_find_interface
        procedure :: contact_pairs => mesh_contact_pairs
        procedure :: split_interface => mesh_split_interface
        procedure :: add_region => mesh_add_region
        procedure :: compute_bandwidth => mesh_compute_bandwidth
        procedure :: equals => mesh_equals
        procedure :: diff => mesh_diff
        ! -- data operations (act on data arrays; geometry is never modified) --
        procedure :: data_drop => mesh_data_drop
        procedure :: data_keep => mesh_data_keep
        procedure :: data_rename => mesh_data_rename
        procedure :: sets_to_data => mesh_sets_to_data
        procedure :: data_to_sets => mesh_data_to_sets
        procedure :: data_point_to_cell => mesh_data_point_to_cell
        procedure :: data_cell_to_point => mesh_data_cell_to_point
        procedure :: data_calc => mesh_data_calc
        procedure :: data_condition => mesh_data_condition
        procedure :: tensor_invariants => mesh_tensor_invariants
        procedure :: data_info => mesh_data_info
        procedure :: data_integrate => mesh_data_integrate
        procedure :: data_integrate_region => mesh_data_integrate_region
        ! -- building --
        procedure :: set_points => mesh_set_points
        procedure, private :: mesh_add_cell_block_i32
        procedure, private :: mesh_add_cell_block_i64
        generic :: add_cell_block => mesh_add_cell_block_i32, mesh_add_cell_block_i64
        procedure :: add_polygon_block => mesh_add_polygon_block
        procedure :: add_polyhedron_block => mesh_add_polyhedron_block
        procedure, private :: mesh_add_point_data_r1
        procedure, private :: mesh_add_point_data_r2
        generic :: add_point_data => mesh_add_point_data_r1, mesh_add_point_data_r2
        procedure :: add_cell_data => mesh_add_cell_data_r1
        procedure :: add_field_data => mesh_add_field_data_r1
        ! -- inspecting --
        procedure :: num_points => mesh_num_points
        procedure :: point_dim => mesh_point_dim
        procedure :: num_cell_blocks => mesh_num_cell_blocks
        procedure :: get_points => mesh_get_points
        procedure :: points_ptr => mesh_points_ptr
        procedure :: cell_block_type => mesh_cell_block_type
        procedure :: cell_block_num_cells => mesh_cell_block_num_cells
        procedure :: cell_block_nodes_per_cell => mesh_cell_block_nodes_per_cell
        procedure :: cell_block_is_ragged => mesh_cell_block_is_ragged
        procedure :: cell_block_is_polyhedron => mesh_cell_block_is_polyhedron
        procedure :: get_cell_block => mesh_get_cell_block
        procedure :: get_polygon_block => mesh_get_polygon_block
        procedure :: get_polyhedron_block => mesh_get_polyhedron_block
        procedure :: num_point_data => mesh_num_point_data
        procedure :: point_data_name => mesh_point_data_name
        procedure, private :: mesh_get_point_data_r1
        procedure, private :: mesh_get_point_data_r2
        generic :: get_point_data => mesh_get_point_data_r1, mesh_get_point_data_r2
        procedure :: num_cell_data => mesh_num_cell_data
        procedure :: cell_data_name => mesh_cell_data_name
        procedure :: cell_data_num_blocks => mesh_cell_data_num_blocks
        procedure :: get_cell_data => mesh_get_cell_data_r1
        procedure :: num_field_data => mesh_num_field_data
        procedure :: field_data_name => mesh_field_data_name
        procedure :: get_field_data => mesh_get_field_data_r1
    end type mio_mesh

    integer, parameter, public :: MIO_SHADING_NONE = 0, MIO_SHADING_FLAT = 1, MIO_SHADING_SMOOTH = 2
    integer, parameter, public :: MIO_EDGES_NONE = 0, MIO_EDGES_ALL = 1, MIO_EDGES_FEATURE = 2
    integer, parameter, public :: MIO_SCALE_LINEAR = 0, MIO_SCALE_LOG = 1, MIO_SCALE_SYMLOG = 2
    integer, parameter, public :: MIO_DIAGNOSTIC_NONE = 0, MIO_DIAGNOSTIC_QUALITY = 1
    integer, parameter, public :: MIO_DIAGNOSTIC_INVERTED = 2, MIO_DIAGNOSTIC_DEGENERATE = 3
    integer, parameter, public :: MIO_DIAGNOSTIC_ORIENTATION = 4, MIO_DIAGNOSTIC_FREE_EDGES = 5
    integer, parameter, public :: MIO_DIAGNOSTIC_EDGE_LENGTH = 6
    integer, parameter, public :: MIO_ENCODING_HALFBLOCK = 0, MIO_ENCODING_QUADRANT = 1
    integer, parameter, public :: MIO_ENCODING_SEXTANT = 2, MIO_ENCODING_BRAILLE = 3
    integer, parameter, public :: MIO_ENCODING_ASCII = 4, MIO_ENCODING_KITTY = 5
    integer, parameter, public :: MIO_ENCODING_ITERM2 = 6, MIO_ENCODING_SIXEL = 7
    integer, parameter, public :: MIO_COLOR_TRUECOLOR = 0, MIO_COLOR_256 = 1
    integer, parameter, public :: MIO_COLOR_16 = 2, MIO_COLOR_MONO = 3
    integer, parameter, public :: MIO_TEXT_ANSI = 0, MIO_TEXT_PLAIN = 1, MIO_TEXT_HTML = 2

    !> Options of `mio_mesh%render` and the functions built on it, in plain Fortran:
    !> allocatable strings, `.true.`/`.false.` flags, and an allocated `vmin`,
    !> `vmax`, `clip_low` or `clip_high` meaning "set". Colours are `integer(int64)`
    !> `0xRRGGBBAA` values (see `mio_rgba`); alpha 0 is transparent. Every default
    !> matches `mio_render_opts_init`. See doc/tui.md for the meaning of each field.
    type :: mio_render_settings
        character(:), allocatable :: view         !< "iso", "+x" ... "-z"
        character(:), allocatable :: color_by     !< array to colour by
        character(:), allocatable :: cmap         !< colormap (default viridis)
        character(:), allocatable :: reduce       !< "mises", "hydrostatic", "principal"
        character(:), allocatable :: expr         !< a data_calc expression to colour by
        character(:), allocatable :: vectors      !< vector point array drawn as arrows
        character(:), allocatable :: warp         !< displacement point array
        character(:), allocatable :: quality_metric !< metric of MIO_DIAGNOSTIC_QUALITY
        real(real64), allocatable :: iso_levels(:) !< explicit contour levels
        real(real64), allocatable :: vmin, vmax    !< explicit range ends
        real(real64), allocatable :: clip_low, clip_high !< percentiles bounding the range
        integer :: width = 320, height = 240, supersample = 1
        integer :: shading = 1, edges = 0, scale = 0, diagnostic = 0
        integer :: component = -1                  !< -1 takes the magnitude
        integer :: isolines = 0, vector_count = 200
        logical :: perspective = .false., two_sided = .true.
        logical :: colorbar = .false., axes = .false., scale_bar = .false.
        logical :: symmetric = .false., categorical = .false.
        logical :: color_regions = .false., category_edges = .false.
        logical :: warp_outline = .false.
        real(real64) :: pixel_aspect = 1.0_real64
        real(real64) :: azimuth = 45.0_real64
        real(real64) :: elevation = 35.264389682754654_real64
        real(real64) :: roll = 0.0_real64
        real(real64) :: fov_deg = 30.0_real64
        real(real64) :: zoom = 1.0_real64
        real(real64) :: pan_x = 0.0_real64
        real(real64) :: pan_y = 0.0_real64
        real(real64) :: ambient = 0.25_real64
        real(real64) :: light_x = 0.0_real64
        real(real64) :: light_y = 0.0_real64
        real(real64) :: light_z = 1.0_real64
        real(real64) :: split_angle = 30.0_real64
        real(real64) :: feature_angle = 30.0_real64
        real(real64) :: point_radius = 1.5_real64
        real(real64) :: scale_threshold = 1.0_real64
        real(real64) :: vector_length = 0.0_real64
        real(real64) :: warp_scale = 1.0_real64
        integer(int64) :: edge_color = 255_int64
        integer(int64) :: fill_color = 3368402431_int64
        integer(int64) :: line_color = 33023_int64
        integer(int64) :: background = 0_int64
        integer(int64) :: nan_color = 2155905279_int64
        integer(int64) :: iso_color = 505290495_int64
        integer(int64) :: vector_color = 3694277375_int64
        integer(int64) :: outline_color = 2526451455_int64
        real(real64), allocatable :: cutaways(:) !< up to two cut-away planes, six numbers each: a point, then the normal of the side kept
        integer(int64) :: cutaway_tint = 3900716287_int64 !< back faces seen through a cut
        character(:), allocatable :: streamlines !< vector point array to follow with streamlines
        integer :: stream_seeds = 40               !< about this many seeds
        real(real64) :: stream_length = 0.5_real64 !< longest line each way, in model diagonals
        integer(int64) :: stream_color = 4031815935_int64 !< colour of the lines
    end type

    !> Options of the text encodings. `encoding`, `color_depth` and `format` take
    !> the MIO_ENCODING_*, MIO_COLOR_* and MIO_TEXT_* constants.
    type :: mio_text_settings
        integer :: encoding = 0, color_depth = 0, format = 0
        integer :: cols = 80, rows = 24
        integer :: cell_pixel_width = 8, cell_pixel_height = 16
        logical :: tmux = .false., notes = .true.
        real(real64) :: cell_aspect = 2.0_real64
    end type

    !> An owned rendered image. Free it with `call frame%free()`.
    type :: mio_frame
        private
        type(c_ptr) :: handle = c_null_ptr
    contains
        procedure :: free => frame_free
        procedure :: is_valid => frame_is_valid
        procedure :: width => frame_width
        procedure :: height => frame_height
        procedure :: get_rgba => frame_get_rgba
        procedure :: get_cell_ids => frame_get_cell_ids
        procedure :: range => frame_range
        procedure :: num_notes => frame_num_notes
        procedure :: note => frame_note
        procedure :: text => frame_text
        procedure :: png => frame_png
    end type

    !> Owned interface mesh, report values and per-source Side entries.
    type :: mio_find_interface_result
        type(mio_mesh) :: mesh
        integer(c_int64_t) :: num_pairs = 0
        real(c_double) :: area = 0.0_c_double
        real(c_double) :: max_gap = 0.0_c_double
        integer(c_int64_t) :: unmatched_a = 0
        integer(c_int64_t) :: unmatched_b = 0
        integer(c_int64_t), allocatable :: side_a(:, :)
        integer(c_int64_t), allocatable :: side_b(:, :)
    end type

    !> Contact projection arrays; per-point coordinate arrays have shape (3,n).
    type :: mio_contact_pairs_result
        integer(c_int64_t), allocatable :: slave_point(:)
        integer(c_int64_t), allocatable :: master_cell(:)
        integer(c_int64_t), allocatable :: master_facet(:)
        integer(c_int64_t), allocatable :: master_subfacet(:)
        real(c_double), allocatable :: local_coordinates(:, :)
        real(c_double), allocatable :: closest_point(:, :)
        real(c_double), allocatable :: gap(:)
        real(c_double), allocatable :: normal(:, :)
        integer(c_int64_t), allocatable :: unmatched(:)
    end type

    !> A multi-file / transient dataset: an ordered PLAN over a set of files
    !> (or the steps inside one multi-step file), read ONE STEP AT A TIME.
    !>
    !> That is the whole point -- a 500-step dataset must be traversable
    !> without materializing it -- so the handle stores paths, per-file step
    !> indices and time values, and never a mesh. `read_step` therefore hands
    !> back an OWNED `mio_mesh` the caller frees, rather than a borrow the
    !> sequence would have to cache.
    !>
    !> Ordering is natural-numeric, so `out_9.vtu` precedes `out_10.vtu`.
    !> Handles are freed explicitly, exactly like `mio_mesh` -- no finalizer.
    !>
    !> ```fortran
    !> type(mio_sequence) :: seq
    !> type(mio_mesh) :: step
    !> integer :: i
    !> call seq%open('out_*.vtu')
    !> do i = 1, seq%count()
    !>     step = seq%read_step(i)      ! 1-based, Fortran style
    !>     ! ... use step, whose time is seq%time(i) ...
    !>     call step%free()
    !> end do
    !> call seq%free()
    !> ```
    !>
    !> See doc/sequences.md.
    type :: mio_sequence
        private
        type(c_ptr) :: handle = c_null_ptr
    contains
        procedure :: open => sequence_open
        procedure :: open_list => sequence_open_list
        procedure :: free => sequence_free
        procedure :: is_valid => sequence_is_valid
        procedure :: count => sequence_count
        procedure :: path => sequence_path
        procedure :: step => sequence_step
        procedure :: time => sequence_time
        procedure :: time_source => sequence_time_source
        procedure :: read_step => sequence_read_step
        procedure :: to_timeseries => sequence_to_timeseries
        procedure :: resample => sequence_resample
    end type mio_sequence

    !> A transient (time-series) XDMF writer: the write half of what the
    !> `time_step` read option and `mio_metadata%time_values` expose on the read
    !> side. A solver writes the mesh ONCE and then one cheap step per solve, so
    !> this is a stateful handle rather than a `write` call -- the one writer
    !> here that `m%write()` cannot express.
    !>
    !>     type(mio_xdmf_series) :: series
    !>     call series%create('simulation.xdmf')        ! "HDF" by default
    !>     call series%write_points_cells(m)            ! the static grid, once
    !>     do k = 1, nsteps
    !>         call solve(m)
    !>         call series%write_data(k*dt, m)          ! point_data/cell_data only
    !>     end do
    !>     call series%finalize()                       ! free() would do it too
    !>     call series%free()
    !>
    !> The `.xdmf` light data is buffered and written at finalize, so by default
    !> a series is only readable after `finalize()` (or `free()`, which finalizes
    !> first). Call `flush()` to make it readable *now* -- what keeps a run that
    !> is killed or still going from leaving nothing openable -- or pass
    !> `auto_flush=.true.` to `create()` to do that after every `write_data`.
    !> `mode='append'` continues a series already at the path instead of
    !> overwriting it; a path with no file yet is simply a fresh series, so a
    !> restartable solver can pass it unconditionally.
    !> Heavy data for `"HDF"` goes to a `<path minus extension>.h5`
    !> SIBLING of the `.xdmf`. Handles are freed explicitly, exactly like
    !> `mio_mesh` -- there is no finalizer.
    type :: mio_xdmf_series
        private
        type(c_ptr) :: handle = c_null_ptr
    contains
        procedure :: create => xdmf_series_create
        procedure :: free => xdmf_series_free
        procedure :: is_valid => xdmf_series_is_valid
        procedure :: write_points_cells => xdmf_series_write_points_cells
        procedure :: write_data => xdmf_series_write_data
        procedure :: flush => xdmf_series_flush
        procedure :: finalize => xdmf_series_finalize
        procedure :: finalized => xdmf_series_finalized
        procedure :: num_steps => xdmf_series_num_steps
    end type mio_xdmf_series

    !> Fixed Exodus grid/sets/attributes with time-dependent fields. Requires netCDF.
    type :: mio_exodus_series
        private
        type(c_ptr) :: handle = c_null_ptr
    contains
        procedure :: create => exodus_series_create
        procedure :: free => exodus_series_free
        procedure :: is_valid => exodus_series_is_valid
        procedure :: write_points_cells => exodus_series_write_points_cells
        procedure :: write_data => exodus_series_write_data
        procedure :: flush => exodus_series_flush
        procedure :: finalize => exodus_series_finalize
        procedure :: finalized => exodus_series_finalized
        procedure :: num_steps => exodus_series_num_steps
    end type mio_exodus_series

    !> A format's side channel: what `m%read_with_info` kept that a mesh
    !> cannot hold (MDPA's tables, geometries, Mesh blocks, constraints, ...),
    !> for `m%write_with_info` to put back. Freed explicitly, like `mio_mesh`.
    !>
    !>     type(mio_mesh) :: m
    !>     type(mio_format_info) :: info
    !>     call m%read_with_info('model.mdpa', info)
    !>     print *, info%mdpa_count(MIO_MDPA_GEOMETRIES)
    !>     call m%write_with_info('out.mdpa', info)
    !>     call info%free(); call m%free()
    !>
    !> Items are 1-based. Point references (geometry connectivity, Mesh-block
    !> nodes) come back 1-based; file ids (geometry, element, condition and
    !> table ids) come back as the file spelled them.
    type :: mio_format_info
        private
        type(c_ptr) :: handle = c_null_ptr
    contains
        procedure :: free => format_info_free
        procedure :: is_valid => format_info_is_valid
        procedure :: format => format_info_format
        procedure :: mdpa_count => format_info_mdpa_count
        procedure :: gmsh_count => format_info_gmsh_count
        procedure :: gmsh_tags => format_info_gmsh_tags
        procedure :: gmsh_affine => format_info_gmsh_affine
        procedure :: gmsh_pairs => format_info_gmsh_pairs
        procedure :: nurbs_count => format_info_nurbs_count
        procedure :: nurbs_int => format_info_nurbs_int
        procedure :: nurbs_degrees => format_info_nurbs_degrees
        procedure :: nurbs_knots => format_info_nurbs_knots
        procedure :: nurbs_control_points => format_info_nurbs_control_points
        procedure :: nurbs_weights => format_info_nurbs_weights
        procedure :: nurbs_field => format_info_nurbs_field
        procedure :: nurbs_field_name => format_info_nurbs_field_name
        procedure :: nurbs_ids => format_info_nurbs_ids
        procedure :: brep_count => format_info_brep_count
        procedure :: brep_int => format_info_brep_int
        procedure :: brep_values => format_info_brep_values
        procedure :: brep_points => format_info_brep_points
        procedure :: brep_ids => format_info_brep_ids
        procedure :: brep_string => format_info_brep_string
        procedure :: mdpa_string => format_info_mdpa_string
        procedure :: mdpa_int => format_info_mdpa_int
        procedure :: mdpa_ids => format_info_mdpa_ids
        procedure :: mdpa_table => format_info_mdpa_table
        procedure :: mdpa_geometry => format_info_mdpa_geometry
        procedure :: mdpa_data_count => format_info_mdpa_data_count
        procedure :: mdpa_data_kind => format_info_mdpa_data_kind
        procedure :: mdpa_data_key => format_info_mdpa_data_key
        procedure :: mdpa_data_text => format_info_mdpa_data_text
        procedure :: mdpa_data_values => format_info_mdpa_data_values
        procedure :: iga_refine => format_info_iga_refine
        procedure :: iga_sample => format_info_iga_sample
        procedure :: iga_evaluate => format_info_iga_evaluate
        procedure :: mesh_cad_surface => format_info_mesh_cad_surface
    end type mio_format_info

    !> The result of `info%iga_evaluate`: points (and, with `derivatives`, their
    !> first derivatives and the surface normals) of a NURBS patch or B-rep face
    !> at parametric coordinates, with its NURBS fields and a trimmed face's
    !> inside mask. Freed explicitly. Arrays are laid out column-major, so the
    !> points are `(space dimension, N)`, the derivatives `(space dimension,
    !> parametric dimension, N)`, the normals `(3, N)`, a field `(components, N)`.
    type :: mio_iga_evaluation
        private
        type(c_ptr) :: handle = c_null_ptr
    contains
        procedure :: free => iga_evaluation_free
        procedure :: is_valid => iga_evaluation_is_valid
        procedure :: num_fields => iga_evaluation_num_fields
        procedure :: points => iga_evaluation_points
        procedure :: derivatives => iga_evaluation_derivatives
        procedure :: normals => iga_evaluation_normals
        procedure :: inside => iga_evaluation_inside
        procedure :: field => iga_evaluation_field
    end type mio_iga_evaluation

    ! ------------------------------------------------------------------
    ! Raw bind(c) interfaces to libmeshioplusplus (private; the OO layer
    ! above is the public surface). Data pointers cross as type(c_ptr) via
    ! c_loc() so one C symbol serves every Fortran type/kind.
    ! ------------------------------------------------------------------
    interface
        function c_mio_version() bind(c, name="mio_version") result(p)
            import :: c_ptr
            type(c_ptr) :: p
        end function

        function c_mio_mesh_backend() bind(c, name="mio_mesh_backend") result(p)
            import :: c_ptr
            type(c_ptr) :: p
        end function

        function c_mio_last_error() bind(c, name="mio_last_error") result(p)
            import :: c_ptr
            type(c_ptr) :: p
        end function

        function c_mio_format_readable(format) bind(c, name="mio_format_readable") result(r)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: format
            integer(c_int) :: r
        end function

        function c_mio_format_writable(format) bind(c, name="mio_format_writable") result(r)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: format
            integer(c_int) :: r
        end function

        function c_mio_mesh_create() bind(c, name="mio_mesh_create") result(h)
            import :: c_ptr
            type(c_ptr) :: h
        end function

        subroutine c_mio_mesh_free(h) bind(c, name="mio_mesh_free")
            import :: c_ptr
            type(c_ptr), value :: h
        end subroutine

        function c_mio_read(path, format) bind(c, name="mio_read") result(h)
            import :: c_ptr, c_char
            character(kind=c_char), dimension(*), intent(in) :: path, format
            type(c_ptr) :: h
        end function

        function c_mio_read_ex(path, format, opts) bind(c, name="mio_read_ex") result(h)
            import :: c_ptr, c_char, mio_read_opts_t
            character(kind=c_char), dimension(*), intent(in) :: path, format
            type(mio_read_opts_t), intent(in) :: opts
            type(c_ptr) :: h
        end function

        subroutine c_mio_read_opts_init(opts) bind(c, name="mio_read_opts_init")
            import :: mio_read_opts_t
            type(mio_read_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_read_metadata_create(path, format) &
                bind(c, name="mio_read_metadata_create") result(h)
            import :: c_ptr, c_char
            character(kind=c_char), dimension(*), intent(in) :: path, format
            type(c_ptr) :: h
        end function

        function c_mio_read_metadata_num_points(h) &
                bind(c, name="mio_read_metadata_num_points") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_read_metadata_point_dim(h) &
                bind(c, name="mio_read_metadata_point_dim") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_read_metadata_num_cells(h) &
                bind(c, name="mio_read_metadata_num_cells") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_read_metadata_num_cell_blocks(h) &
                bind(c, name="mio_read_metadata_num_cell_blocks") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_read_metadata_cell_block(h, index, num_cells, npc, ragged) &
                bind(c, name="mio_read_metadata_cell_block") result(s)
            import :: c_ptr, c_int64_t, c_int
            type(c_ptr), value :: h
            integer(c_int64_t), value :: index
            integer(c_int64_t), intent(out) :: num_cells, npc
            integer(c_int), intent(out) :: ragged
            integer(c_int) :: s
        end function

        function c_mio_read_metadata_cell_block_type(h, index, buf, buflen) &
                bind(c, name="mio_read_metadata_cell_block_type") result(n)
            import :: c_ptr, c_int64_t, c_char
            type(c_ptr), value :: h
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(out) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_read_metadata_num_names(h, location) &
                bind(c, name="mio_read_metadata_num_names") result(n)
            import :: c_ptr, c_int64_t, c_int
            type(c_ptr), value :: h
            integer(c_int), value :: location
            integer(c_int64_t) :: n
        end function

        function c_mio_read_metadata_num_time_values(h) &
                bind(c, name="mio_read_metadata_num_time_values") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_read_metadata_time_values(h, out, count) &
                bind(c, name="mio_read_metadata_time_values") result(n)
            import :: c_ptr, c_int64_t, c_double
            type(c_ptr), value :: h
            real(c_double), intent(out) :: out(*)
            integer(c_int64_t), value :: count
            integer(c_int64_t) :: n
        end function

        function c_mio_read_metadata_num_regions(h) &
                bind(c, name="mio_read_metadata_num_regions") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_read_metadata_region_name(h, index, buf, buflen) &
                bind(c, name="mio_read_metadata_region_name") result(n)
            import :: c_ptr, c_int64_t, c_char
            type(c_ptr), value :: h
            integer(c_int64_t), value :: index
            character(kind=c_char), intent(out) :: buf(*)
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        function c_mio_read_metadata_region_info(h, index, out) &
                bind(c, name="mio_read_metadata_region_info") result(s)
            import :: c_ptr, c_int64_t, c_int, mio_region_info
            type(c_ptr), value :: h
            integer(c_int64_t), value :: index
            type(mio_region_info), intent(out) :: out
            integer(c_int) :: s
        end function

        function c_mio_read_metadata_name(h, location, index, buf, buflen) &
                bind(c, name="mio_read_metadata_name") result(n)
            import :: c_ptr, c_int64_t, c_int, c_char
            type(c_ptr), value :: h
            integer(c_int), value :: location
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(out) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_read_metadata_fell_back(h) &
                bind(c, name="mio_read_metadata_fell_back") result(n)
            import :: c_ptr, c_int
            type(c_ptr), value :: h
            integer(c_int) :: n
        end function

        subroutine c_mio_read_metadata_free(h) bind(c, name="mio_read_metadata_free")
            import :: c_ptr
            type(c_ptr), value :: h
        end subroutine

        function c_mio_write(path, h, format) bind(c, name="mio_write") result(s)
            import :: c_ptr, c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: path, format
            type(c_ptr), value :: h
            integer(c_int) :: s
        end function

        function c_mio_write_ex(path, h, format, opts) bind(c, name="mio_write_ex") result(s)
            import :: c_ptr, c_char, c_int, mio_write_opts_c
            character(kind=c_char), dimension(*), intent(in) :: path, format
            type(c_ptr), value :: h
            type(mio_write_opts_c), intent(in) :: opts
            integer(c_int) :: s
        end function

        function c_mio_write_gltf(path, h, opts) bind(c, name="mio_write_gltf") result(s)
            import :: c_ptr, c_char, c_int, mio_gltf_opts_c
            character(kind=c_char), dimension(*), intent(in) :: path
            type(c_ptr), value :: h
            type(mio_gltf_opts_c), intent(in) :: opts
            integer(c_int) :: s
        end function

        function c_mio_convert(in_path, in_format, out_path, out_format) &
                bind(c, name="mio_convert") result(s)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: in_path, in_format
            character(kind=c_char), dimension(*), intent(in) :: out_path, out_format
            integer(c_int) :: s
        end function

        function c_mio_med_mesh_count(path) bind(c, name="mio_med_mesh_count") result(n)
            import :: c_char, c_int64_t
            character(c_char), intent(in) :: path(*)
            integer(c_int64_t) :: n
        end function

        function c_mio_med_mesh_name(path, index, buf, buflen) bind(c, name="mio_med_mesh_name") result(n)
            import :: c_char, c_int64_t
            character(c_char), intent(in) :: path(*)
            integer(c_int64_t), value :: index, buflen
            character(c_char), intent(out) :: buf(*)
            integer(c_int64_t) :: n
        end function

        function c_mio_med_read_named(path, name, opts) bind(c, name="mio_med_read_named") result(h)
            import :: c_char, c_ptr, mio_read_opts_t
            character(c_char), intent(in) :: path(*), name(*)
            type(mio_read_opts_t), intent(in) :: opts
            type(c_ptr) :: h
        end function

        function c_mio_med_write_multi(path, meshes, names, count, version) &
                bind(c, name="mio_med_write_multi") result(s)
            import :: c_char, c_ptr, c_int64_t, c_int
            character(c_char), intent(in) :: path(*), version(*)
            type(c_ptr), intent(in) :: meshes(*), names(*)
            integer(c_int64_t), value :: count
            integer(c_int) :: s
        end function

        function c_mio_pipeline_run_file(settings_path) &
                bind(c, name="mio_pipeline_run_file") result(s)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: settings_path
            integer(c_int) :: s
        end function

        function c_mio_pipeline_run_file_report(path) bind(c, name="mio_pipeline_run_file_report") result(h)
            import :: c_char, c_ptr
            character(c_char), intent(in) :: path(*)
            type(c_ptr) :: h
        end function

        function c_mio_pipeline_run_json_report(text) bind(c, name="mio_pipeline_run_json_report") result(h)
            import :: c_char, c_ptr
            character(c_char), intent(in) :: text(*)
            type(c_ptr) :: h
        end function

        function c_mio_sequence_pipeline_run_file_report(path) &
                bind(c, name="mio_sequence_pipeline_run_file_report") result(h)
            import :: c_char, c_ptr
            character(c_char), intent(in) :: path(*)
            type(c_ptr) :: h
        end function

        function c_mio_sequence_pipeline_run_json_report(text) &
                bind(c, name="mio_sequence_pipeline_run_json_report") result(h)
            import :: c_char, c_ptr
            character(c_char), intent(in) :: text(*)
            type(c_ptr) :: h
        end function

        function c_mio_pipeline_report_json(h, buf, buflen) bind(c, name="mio_pipeline_report_json") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: h
            character(c_char), intent(out) :: buf(*)
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        subroutine c_mio_pipeline_report_free(h) bind(c, name="mio_pipeline_report_free")
            import :: c_ptr
            type(c_ptr), value :: h
        end subroutine

        function c_mio_pipeline_run_json(json_text) &
                bind(c, name="mio_pipeline_run_json") result(s)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: json_text
            integer(c_int) :: s
        end function

        function c_mio_sequence_open(pattern) &
                bind(c, name="mio_sequence_open") result(p)
            import :: c_char, c_ptr
            character(kind=c_char), dimension(*), intent(in) :: pattern
            type(c_ptr) :: p
        end function

        function c_mio_sequence_open_list(paths, num_paths) &
                bind(c, name="mio_sequence_open_list") result(p)
            import :: c_ptr, c_int64_t
            type(c_ptr), dimension(*), intent(in) :: paths
            integer(c_int64_t), value :: num_paths
            type(c_ptr) :: p
        end function

        function c_mio_sequence_count(seq) &
                bind(c, name="mio_sequence_count") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: seq
            integer(c_int64_t) :: n
        end function

        function c_mio_sequence_path(seq, index, buf, buflen) &
                bind(c, name="mio_sequence_path") result(n)
            import :: c_ptr, c_int64_t, c_char
            type(c_ptr), value :: seq
            integer(c_int64_t), value :: index
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        function c_mio_sequence_step(seq, index) &
                bind(c, name="mio_sequence_step") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: seq
            integer(c_int64_t), value :: index
            integer(c_int64_t) :: n
        end function

        function c_mio_sequence_time(seq, index, out_time) &
                bind(c, name="mio_sequence_time") result(s)
            import :: c_ptr, c_int64_t, c_double, c_int
            type(c_ptr), value :: seq
            integer(c_int64_t), value :: index
            real(c_double), intent(out) :: out_time
            integer(c_int) :: s
        end function

        function c_mio_sequence_time_source(seq, index) &
                bind(c, name="mio_sequence_time_source") result(r)
            import :: c_ptr, c_int64_t, c_int32_t
            type(c_ptr), value :: seq
            integer(c_int64_t), value :: index
            integer(c_int32_t) :: r
        end function

        function c_mio_sequence_read(seq, index) &
                bind(c, name="mio_sequence_read") result(p)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: seq
            integer(c_int64_t), value :: index
            type(c_ptr) :: p
        end function

        subroutine c_mio_sequence_free(seq) bind(c, name="mio_sequence_free")
            import :: c_ptr
            type(c_ptr), value :: seq
        end subroutine

        function c_mio_sequence_resample(seq, out_path, out_format, opts) &
                bind(c, name="mio_sequence_resample") result(s)
            import :: c_ptr, c_char, c_int, mio_resample_opts_t
            type(c_ptr), value :: seq
            character(kind=c_char), dimension(*), intent(in) :: out_path
            character(kind=c_char), dimension(*), intent(in) :: out_format
            type(mio_resample_opts_t), intent(in) :: opts
            integer(c_int) :: s
        end function

        function c_mio_sequence_to_timeseries(seq, out_path, out_format) &
                bind(c, name="mio_sequence_to_timeseries") result(s)
            import :: c_ptr, c_char, c_int
            type(c_ptr), value :: seq
            character(kind=c_char), dimension(*), intent(in) :: out_path
            character(kind=c_char), dimension(*), intent(in) :: out_format
            integer(c_int) :: s
        end function

        function c_mio_timeseries_to_sequence(in_path, in_format, out_pattern, out_format) &
                bind(c, name="mio_timeseries_to_sequence") result(s)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: in_path
            character(kind=c_char), dimension(*), intent(in) :: in_format
            character(kind=c_char), dimension(*), intent(in) :: out_pattern
            character(kind=c_char), dimension(*), intent(in) :: out_format
            integer(c_int) :: s
        end function

        function c_mio_sequence_pipeline_run_file(settings_path) &
                bind(c, name="mio_sequence_pipeline_run_file") result(s)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: settings_path
            integer(c_int) :: s
        end function

        function c_mio_sequence_pipeline_run_json(json_text) &
                bind(c, name="mio_sequence_pipeline_run_json") result(s)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: json_text
            integer(c_int) :: s
        end function

        function c_mio_pipeline_has_json() &
                bind(c, name="mio_pipeline_has_json") result(r)
            import :: c_int32_t
            integer(c_int32_t) :: r
        end function

        function c_mio_extract_surface(h, record_parent_ids) &
                bind(c, name="mio_extract_surface") result(r)
            import :: c_ptr, c_int
            type(c_ptr), value :: h
            integer(c_int), value :: record_parent_ids
            type(c_ptr) :: r
        end function

        function c_mio_extract_skin(h, linearize) bind(c, name="mio_extract_skin") result(r)
            import :: c_ptr, c_int
            type(c_ptr), value :: h
            integer(c_int), value :: linearize
            type(c_ptr) :: r
        end function

        function c_mio_attach_quality(h) bind(c, name="mio_attach_quality") result(r)
            import :: c_ptr
            type(c_ptr), value :: h
            type(c_ptr) :: r
        end function

        function c_mio_quality_counts(h, nc, ninv, ndeg) &
                bind(c, name="mio_quality_counts") result(s)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), intent(out) :: nc, ninv, ndeg
            integer(c_int) :: s
        end function

        function c_mio_merge(meshes, count, weld, atol, source_tag, data_policy, drop_dup) &
                bind(c, name="mio_merge") result(r)
            import :: c_ptr, c_int, c_int64_t, c_double
            type(c_ptr), intent(in) :: meshes(*)
            integer(c_int64_t), value :: count
            integer(c_int), value :: weld, source_tag, data_policy, drop_dup
            real(c_double), value :: atol
            type(c_ptr) :: r
        end function

        function c_mio_interpolate(source, target, method, arrays, count, extrapolate, &
                                   default_value, on_conflict) &
                bind(c, name="mio_interpolate") result(r)
            import :: c_ptr, c_char, c_int, c_int64_t, c_double
            type(c_ptr), value :: source, target
            character(kind=c_char), dimension(*), intent(in) :: method
            type(c_ptr), value :: arrays
            integer(c_int64_t), value :: count
            integer(c_int), value :: extrapolate
            real(c_double), value :: default_value
            character(kind=c_char), dimension(*), intent(in) :: on_conflict
            type(c_ptr) :: r
        end function

        function c_mio_conservative_interpolate(source, target, arrays, count, &
                                                default_value, on_conflict) &
                bind(c, name="mio_conservative_interpolate") result(r)
            import :: c_ptr, c_char, c_int64_t, c_double
            type(c_ptr), value :: source, target
            type(c_ptr), value :: arrays
            integer(c_int64_t), value :: count
            real(c_double), value :: default_value
            character(kind=c_char), dimension(*), intent(in) :: on_conflict
            type(c_ptr) :: r
        end function

        function c_mio_undo_green(coarse, fine, ngroups, nremoved) &
                bind(c, name="mio_undo_green") result(r)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: coarse, fine
            integer(c_int64_t), intent(out) :: ngroups, nremoved
            type(c_ptr) :: r
        end function

        function c_mio_reorder(h, method) bind(c, name="mio_reorder") result(r)
            import :: c_ptr, c_char
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: method
            type(c_ptr) :: r
        end function

        function c_mio_transform(h, matrix, rotate_data) bind(c, name="mio_transform") result(r)
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: h
            real(c_double), intent(in) :: matrix(*)
            integer(c_int), value :: rotate_data
            type(c_ptr) :: r
        end function

        function c_mio_clean(h, weld, atol, rmorph, ddeg, ddup, nweld, norph, ndeg, ndup) &
                bind(c, name="mio_clean") result(r)
            import :: c_ptr, c_int, c_int64_t, c_double
            type(c_ptr), value :: h
            integer(c_int), value :: weld, rmorph, ddeg, ddup
            real(c_double), value :: atol
            integer(c_int64_t), intent(out) :: nweld, norph, ndeg, ndup
            type(c_ptr) :: r
        end function

        function c_mio_smooth(h, method, iters, lambda, mu, fixb, feat, fangle, guard, &
                              nmoved, maxdisp, nskip) bind(c, name="mio_smooth") result(r)
            import :: c_ptr, c_char, c_int, c_int64_t, c_double
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: method
            integer(c_int), value :: iters, fixb, feat, guard
            real(c_double), value :: lambda, mu, fangle
            integer(c_int64_t), intent(out) :: nmoved, nskip
            real(c_double), intent(out) :: maxdisp
            type(c_ptr) :: r
        end function

        function c_mio_optimize_volume(h, maxiter, relocate, flip, presb, minimp, nflips, n23, &
                                       n32, nmoved, ntets, qbefore, qafter) &
                bind(c, name="mio_optimize_volume") result(r)
            import :: c_ptr, c_int, c_int64_t, c_double
            type(c_ptr), value :: h
            integer(c_int), value :: maxiter, relocate, flip, presb
            real(c_double), value :: minimp
            integer(c_int64_t), intent(out) :: nflips, n23, n32, nmoved, ntets
            real(c_double), intent(out) :: qbefore, qafter
            type(c_ptr) :: r
        end function

        function c_mio_crop_bbox(h, lo, hi, mode, record_ids) &
                bind(c, name="mio_crop_bbox") result(r)
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: h
            real(c_double), intent(in) :: lo(*), hi(*)
            integer(c_int), value :: mode, record_ids
            type(c_ptr) :: r
        end function

        function c_mio_crop_plane(h, point, normal, mode, record_ids) &
                bind(c, name="mio_crop_plane") result(r)
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: h
            real(c_double), intent(in) :: point(*), normal(*)
            integer(c_int), value :: mode, record_ids
            type(c_ptr) :: r
        end function

        function c_mio_slice(h, origin, normal, record_parent_ids) &
                bind(c, name="mio_slice") result(r)
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: h
            real(c_double), intent(in) :: origin(*), normal(*)
            integer(c_int), value :: record_parent_ids
            type(c_ptr) :: r
        end function

        function c_mio_isosurface(h, array_name, isovalues, n_isovalues, component, &
                                  record_parent_ids) &
                bind(c, name="mio_isosurface") result(r)
            import :: c_ptr, c_int, c_double, c_char
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: array_name
            real(c_double), intent(in) :: isovalues(*)
            integer(c_int), value :: n_isovalues, component, record_parent_ids
            type(c_ptr) :: r
        end function

        function c_mio_gradient(h, array_name, op, method, location, output_name, &
                               component, overwrite, n_skipped, n_fallback) &
                bind(c, name="mio_gradient") result(r)
            import :: c_ptr, c_int, c_int64_t, c_char
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: array_name, op, method
            character(kind=c_char), dimension(*), intent(in) :: location, output_name
            integer(c_int), value :: component, overwrite
            integer(c_int64_t), intent(out) :: n_skipped, n_fallback
            type(c_ptr) :: r
        end function

        function c_mio_hessian(h, array_name, method, location, output_name, &
                              overwrite, n_skipped, n_fallback) &
                bind(c, name="mio_hessian") result(r)
            import :: c_ptr, c_int, c_int64_t, c_char
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: array_name, method, location
            character(kind=c_char), dimension(*), intent(in) :: output_name
            integer(c_int), value :: overwrite
            integer(c_int64_t), intent(out) :: n_skipped, n_fallback
            type(c_ptr) :: r
        end function

        function c_mio_estimate_error(h, array_name, method, marking, marking_value, &
                                      output_name, marked_name, overwrite, global_error, &
                                      n_skipped, n_marked) &
                bind(c, name="mio_estimate_error") result(r)
            import :: c_ptr, c_int, c_int64_t, c_char, c_double
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: array_name, method, marking
            character(kind=c_char), dimension(*), intent(in) :: output_name, marked_name
            real(c_double), value :: marking_value
            integer(c_int), value :: overwrite
            real(c_double), intent(out) :: global_error
            integer(c_int64_t), intent(out) :: n_skipped, n_marked
            type(c_ptr) :: r
        end function
        function c_mio_remesh(h, num_clusters, subdivide, subsample_ratio, max_subdivide, &
                              max_iterations, max_repair_passes, metric, gradation, &
                              preserve_boundary, num_clusters_out, num_iterations, &
                              subdivide_applied, num_isolated_clusters, &
                              num_non_manifold_vertices) &
                bind(c, name="mio_remesh") result(r)
            import :: c_ptr, c_int, c_int64_t, c_char, c_double
            type(c_ptr), value :: h
            integer(c_int64_t), value :: num_clusters
            integer(c_int), value :: subdivide, max_subdivide, max_iterations, max_repair_passes
            real(c_double), value :: subsample_ratio
            character(kind=c_char), dimension(*), intent(in) :: metric
            real(c_double), value :: gradation
            integer(c_int), value :: preserve_boundary
            integer(c_int64_t), intent(out) :: num_clusters_out, num_iterations
            integer(c_int), intent(out) :: subdivide_applied
            integer(c_int64_t), intent(out) :: num_isolated_clusters
            integer(c_int64_t), intent(out) :: num_non_manifold_vertices
            type(c_ptr) :: r
        end function

        subroutine c_mio_remesh_opts_init(opts) bind(c, name="mio_remesh_opts_init")
            import :: mio_remesh_opts_t
            type(mio_remesh_opts_t), intent(out) :: opts
        end subroutine

        subroutine c_mio_curvature_opts_init(opts) bind(c, name="mio_curvature_opts_init")
            import :: mio_curvature_opts_t
            type(mio_curvature_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_compute_curvature(h, opts, report) &
                bind(c, name="mio_compute_curvature") result(r)
            import :: c_ptr, mio_curvature_opts_t, mio_curvature_report_t
            type(c_ptr), value :: h
            type(mio_curvature_opts_t), intent(in) :: opts
            type(mio_curvature_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        subroutine c_mio_normals_opts_init(opts) bind(c, name="mio_normals_opts_init")
            import :: mio_normals_opts_t
            type(mio_normals_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_compute_normals(h, opts, report) &
                bind(c, name="mio_compute_normals") result(r)
            import :: c_ptr, mio_normals_opts_t, mio_normals_report_t
            type(c_ptr), value :: h
            type(mio_normals_opts_t), intent(in) :: opts
            type(mio_normals_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        subroutine c_mio_feature_edges_opts_init(opts) &
                bind(c, name="mio_feature_edges_opts_init")
            import :: mio_feature_edges_opts_t
            type(mio_feature_edges_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_feature_edges(h, opts, report) &
                bind(c, name="mio_feature_edges") result(r)
            import :: c_ptr, mio_feature_edges_opts_t, mio_feature_edges_report_t
            type(c_ptr), value :: h
            type(mio_feature_edges_opts_t), intent(in) :: opts
            type(mio_feature_edges_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        subroutine c_mio_segment_opts_init(opts) bind(c, name="mio_segment_opts_init")
            import :: mio_segment_opts_t
            type(mio_segment_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_segment(h, opts, report) bind(c, name="mio_segment") result(r)
            import :: c_ptr, mio_segment_opts_t, mio_segment_report_t
            type(c_ptr), value :: h
            type(mio_segment_opts_t), intent(in) :: opts
            type(mio_segment_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        subroutine c_mio_fit_primitives_opts_init(opts) &
                bind(c, name="mio_fit_primitives_opts_init")
            import :: mio_fit_primitives_opts_t
            type(mio_fit_primitives_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_fit_primitives(h, opts, report) &
                bind(c, name="mio_fit_primitives") result(r)
            import :: c_ptr, mio_fit_primitives_opts_t, mio_fit_primitives_report_t
            type(c_ptr), value :: h
            type(mio_fit_primitives_opts_t), intent(in) :: opts
            type(mio_fit_primitives_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        subroutine c_mio_fit_cad_opts_init(opts) &
                bind(c, name="mio_fit_cad_opts_init")
            import :: mio_fit_cad_opts_t
            type(mio_fit_cad_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_fit_cad(h, path, opts, report) &
                bind(c, name="mio_fit_cad") result(r)
            import :: c_ptr, mio_fit_cad_opts_t, mio_fit_cad_report_t
            type(c_ptr), value :: h, path
            type(mio_fit_cad_opts_t), intent(in) :: opts
            type(mio_fit_cad_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        function c_mio_fit_cad_ex(h, path, opts, report, report_ex) &
                bind(c, name="mio_fit_cad_ex") result(r)
            import :: c_ptr, mio_fit_cad_opts_t, mio_fit_cad_report_t, mio_fit_cad_report_ex_t
            type(c_ptr), value :: h, path
            type(mio_fit_cad_opts_t), intent(in) :: opts
            type(mio_fit_cad_report_t), intent(out) :: report
            type(mio_fit_cad_report_ex_t), intent(out) :: report_ex
            type(c_ptr) :: r
        end function

        subroutine c_mio_refine_patches_opts_init(opts) &
                bind(c, name="mio_refine_patches_opts_init")
            import :: mio_refine_patches_opts_t
            type(mio_refine_patches_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_refine_patches(h, opts, report) &
                bind(c, name="mio_refine_patches") result(r)
            import :: c_ptr, mio_refine_patches_opts_t, mio_refine_patches_report_t
            type(c_ptr), value :: h
            type(mio_refine_patches_opts_t), intent(in) :: opts
            type(mio_refine_patches_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        subroutine c_mio_fit_cad_opts_ex_init(opts) &
                bind(c, name="mio_fit_cad_opts_ex_init")
            import :: mio_fit_cad_opts_ex_t
            type(mio_fit_cad_opts_ex_t), intent(out) :: opts
        end subroutine

        function c_mio_fit_cad2(h, path, opts, opts_ex, report, report_ex, report_ex2) &
                bind(c, name="mio_fit_cad2") result(r)
            import :: c_ptr, mio_fit_cad_opts_t, mio_fit_cad_opts_ex_t, mio_fit_cad_report_t, &
                      mio_fit_cad_report_ex_t, mio_fit_cad_report_ex2_t
            type(c_ptr), value :: h, path
            type(mio_fit_cad_opts_t), intent(in) :: opts
            type(mio_fit_cad_opts_ex_t), intent(in) :: opts_ex
            type(mio_fit_cad_report_t), intent(out) :: report
            type(mio_fit_cad_report_ex_t), intent(out) :: report_ex
            type(mio_fit_cad_report_ex2_t), intent(out) :: report_ex2
            type(c_ptr) :: r
        end function

        subroutine c_mio_fit_cad_opts_ex2_init(opts) &
                bind(c, name="mio_fit_cad_opts_ex2_init")
            import :: mio_fit_cad_opts_ex2_t
            type(mio_fit_cad_opts_ex2_t), intent(out) :: opts
        end subroutine

        function c_mio_fit_cad3(h, path, opts, opts_ex, opts_ex2, report, report_ex, &
                                report_ex2, report_ex3) bind(c, name="mio_fit_cad3") result(r)
            import :: c_ptr, mio_fit_cad_opts_t, mio_fit_cad_opts_ex_t, mio_fit_cad_opts_ex2_t, &
                      mio_fit_cad_report_t, mio_fit_cad_report_ex_t, mio_fit_cad_report_ex2_t, &
                      mio_fit_cad_report_ex3_t
            type(c_ptr), value :: h, path
            type(mio_fit_cad_opts_t), intent(in) :: opts
            type(mio_fit_cad_opts_ex_t), intent(in) :: opts_ex
            type(mio_fit_cad_opts_ex2_t), intent(in) :: opts_ex2
            type(mio_fit_cad_report_t), intent(out) :: report
            type(mio_fit_cad_report_ex_t), intent(out) :: report_ex
            type(mio_fit_cad_report_ex2_t), intent(out) :: report_ex2
            type(mio_fit_cad_report_ex3_t), intent(out) :: report_ex3
            type(c_ptr) :: r
        end function

        function c_mio_has_occt() bind(c, name="mio_has_occt") result(r)
            import :: c_int
            integer(c_int) :: r
        end function

        function c_mio_hausdorff_distance(a, b, opts, report) &
                bind(c, name="mio_hausdorff_distance") result(r)
            import :: c_ptr, c_int, mio_hausdorff_opts_t, mio_hausdorff_report_t
            type(c_ptr), value :: a, b
            type(mio_hausdorff_opts_t), intent(in) :: opts
            type(mio_hausdorff_report_t), intent(out) :: report
            integer(c_int) :: r
        end function

        function c_mio_edit_regions(h, op, inputs, num_inputs, output, dim, tag, &
                                    keep_inputs) bind(c, name="mio_edit_regions") result(r)
            import :: c_ptr, c_int32_t, c_int64_t, mio_region_selector_t
            type(c_ptr), value :: h
            integer(c_int32_t), value :: op
            type(mio_region_selector_t), intent(in) :: inputs(*)
            integer(c_int64_t), value :: num_inputs
            type(c_ptr), value :: output
            integer(c_int64_t), value :: dim, tag
            integer(c_int32_t), value :: keep_inputs
            type(c_ptr) :: r
        end function

        function c_mio_region_adjacency(h, regions, num_regions) &
                bind(c, name="mio_region_adjacency") result(r)
            import :: c_ptr, c_int64_t, mio_region_selector_t
            type(c_ptr), value :: h
            type(mio_region_selector_t), intent(in) :: regions(*)
            integer(c_int64_t), value :: num_regions
            type(c_ptr) :: r
        end function

        subroutine c_mio_find_interface_opts_init(opts) &
                bind(c, name="mio_find_interface_opts_init")
            import :: mio_find_interface_opts_t
            type(mio_find_interface_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_find_interface(a, region_a, b, region_b, opts) &
                bind(c, name="mio_find_interface") result(r)
            import :: c_ptr, mio_region_selector_t, mio_find_interface_opts_t
            type(c_ptr), value :: a, b
            type(mio_region_selector_t), intent(in) :: region_a, region_b
            type(mio_find_interface_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_find_interface_result_take_mesh(r) &
                bind(c, name="mio_find_interface_result_take_mesh") result(h)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: h
        end function

        function c_mio_find_interface_result_report(r, pairs, area, gap, unmatched_a, &
                                                    unmatched_b) &
                bind(c, name="mio_find_interface_result_report") result(status)
            import :: c_ptr, c_int, c_int64_t, c_double
            type(c_ptr), value :: r
            integer(c_int64_t), intent(out) :: pairs, unmatched_a, unmatched_b
            real(c_double), intent(out) :: area, gap
            integer(c_int) :: status
        end function

        function c_mio_find_interface_result_side_a(r, count) &
                bind(c, name="mio_find_interface_result_side_a") result(data)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t), intent(out) :: count
            type(c_ptr) :: data
        end function

        function c_mio_find_interface_result_side_b(r, count) &
                bind(c, name="mio_find_interface_result_side_b") result(data)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t), intent(out) :: count
            type(c_ptr) :: data
        end function

        subroutine c_mio_find_interface_result_free(r) &
                bind(c, name="mio_find_interface_result_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        subroutine c_mio_contact_pairs_opts_init(opts) &
                bind(c, name="mio_contact_pairs_opts_init")
            import :: mio_contact_pairs_opts_t
            type(mio_contact_pairs_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_contact_pairs(slave, slave_points, master, master_cells, opts) &
                bind(c, name="mio_contact_pairs") result(r)
            import :: c_ptr, mio_region_selector_t, mio_contact_pairs_opts_t
            type(c_ptr), value :: slave, master
            type(mio_region_selector_t), intent(in) :: slave_points, master_cells
            type(mio_contact_pairs_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_contact_pairs_result_info(r, count, unmatched) &
                bind(c, name="mio_contact_pairs_result_info") result(status)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t), intent(out) :: count, unmatched
            integer(c_int) :: status
        end function

        function c_mio_contact_pairs_slave_point(r) &
                bind(c, name="mio_contact_pairs_slave_point") result(data)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: data
        end function

        function c_mio_contact_pairs_master_cell(r) &
                bind(c, name="mio_contact_pairs_master_cell") result(data)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: data
        end function

        function c_mio_contact_pairs_master_facet(r) &
                bind(c, name="mio_contact_pairs_master_facet") result(data)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: data
        end function

        function c_mio_contact_pairs_master_subfacet(r) &
                bind(c, name="mio_contact_pairs_master_subfacet") result(data)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: data
        end function

        function c_mio_contact_pairs_local_coordinates(r) &
                bind(c, name="mio_contact_pairs_local_coordinates") result(data)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: data
        end function

        function c_mio_contact_pairs_closest_point(r) &
                bind(c, name="mio_contact_pairs_closest_point") result(data)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: data
        end function

        function c_mio_contact_pairs_gap(r) bind(c, name="mio_contact_pairs_gap") result(data)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: data
        end function

        function c_mio_contact_pairs_normal(r) &
                bind(c, name="mio_contact_pairs_normal") result(data)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: data
        end function

        function c_mio_contact_pairs_unmatched(r) &
                bind(c, name="mio_contact_pairs_unmatched") result(data)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: data
        end function

        subroutine c_mio_contact_pairs_result_free(r) &
                bind(c, name="mio_contact_pairs_result_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        subroutine c_mio_split_interface_opts_init(opts) &
                bind(c, name="mio_split_interface_opts_init")
            import :: mio_split_interface_opts_t
            type(mio_split_interface_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_split_interface(h, side, opts, duplicated, cohesive) &
                bind(c, name="mio_split_interface") result(r)
            import :: c_ptr, c_int64_t, mio_region_selector_t, mio_split_interface_opts_t
            type(c_ptr), value :: h
            type(mio_region_selector_t), intent(in) :: side
            type(mio_split_interface_opts_t), intent(in) :: opts
            integer(c_int64_t), intent(out) :: duplicated, cohesive
            type(c_ptr) :: r
        end function

        function c_mio_mesh_remove_region(h, index) &
                bind(c, name="mio_mesh_remove_region") result(r)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: index
            integer(c_int) :: r
        end function

        subroutine c_mio_periodic_opts_init(opts) bind(c, name="mio_periodic_opts_init")
            import :: mio_periodic_opts_t
            type(mio_periodic_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_match_periodic_nodes(h, slave, master, opts) &
                bind(c, name="mio_match_periodic_nodes") result(r)
            import :: c_ptr, mio_region_selector_t, mio_periodic_opts_t
            type(c_ptr), value :: h
            type(mio_region_selector_t), intent(in) :: slave, master
            type(mio_periodic_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_periodic_pairs_info(p, num_pairs, num_unmatched, num_fixed, &
                                           max_residual) &
                bind(c, name="mio_periodic_pairs_info") result(r)
            import :: c_ptr, c_int, c_int64_t, c_double
            type(c_ptr), value :: p
            integer(c_int64_t), intent(out) :: num_pairs, num_unmatched, num_fixed
            real(c_double), intent(out) :: max_residual
            integer(c_int) :: r
        end function

        function c_mio_periodic_pairs_slave(p, count) &
                bind(c, name="mio_periodic_pairs_slave") result(r)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: p
            integer(c_int64_t), intent(out) :: count
            type(c_ptr) :: r
        end function

        function c_mio_periodic_pairs_master(p, count) &
                bind(c, name="mio_periodic_pairs_master") result(r)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: p
            integer(c_int64_t), intent(out) :: count
            type(c_ptr) :: r
        end function

        function c_mio_periodic_pairs_unmatched(p, count) &
                bind(c, name="mio_periodic_pairs_unmatched") result(r)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: p
            integer(c_int64_t), intent(out) :: count
            type(c_ptr) :: r
        end function

        subroutine c_mio_periodic_pairs_free(p) bind(c, name="mio_periodic_pairs_free")
            import :: c_ptr
            type(c_ptr), value :: p
        end subroutine

        function c_mio_check_quality(h, spec, max_inverted, max_degenerate, report, summary, &
                                     summary_len) bind(c, name="mio_check_quality") result(r)
            import :: c_ptr, c_int, c_int64_t, c_char, mio_quality_gate_report_t
            type(c_ptr), value :: h
            character(kind=c_char), intent(in) :: spec(*)
            integer(c_int64_t), value :: max_inverted, max_degenerate
            type(mio_quality_gate_report_t), intent(out) :: report
            character(kind=c_char), intent(out) :: summary(*)
            integer(c_int64_t), value :: summary_len
            integer(c_int) :: r
        end function

        subroutine c_mio_repair_opts_init(opts) bind(c, name="mio_repair_opts_init")
            import :: mio_repair_opts_t
            type(mio_repair_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_repair(h, opts, report) bind(c, name="mio_repair") result(r)
            import :: c_ptr, mio_repair_opts_t, mio_repair_report_t
            type(c_ptr), value :: h
            type(mio_repair_opts_t), intent(in) :: opts
            type(mio_repair_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        subroutine c_mio_shrinkwrap_opts_init(opts) bind(c, name="mio_shrinkwrap_opts_init")
            import :: mio_shrinkwrap_opts_t
            type(mio_shrinkwrap_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_shrinkwrap(h, target, opts, report) bind(c, name="mio_shrinkwrap") result(r)
            import :: c_ptr, mio_shrinkwrap_opts_t, mio_shrinkwrap_report_t
            type(c_ptr), value :: h, target
            type(mio_shrinkwrap_opts_t), intent(in) :: opts
            type(mio_shrinkwrap_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        subroutine c_mio_sobolev_opts_init(opts) bind(c, name="mio_sobolev_opts_init")
            import :: mio_sobolev_opts_t
            type(mio_sobolev_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_sobolev_deform(h, opts, report) bind(c, name="mio_sobolev_deform") result(r)
            import :: c_ptr, mio_sobolev_opts_t, mio_sobolev_report_t
            type(c_ptr), value :: h
            type(mio_sobolev_opts_t), intent(in) :: opts
            type(mio_sobolev_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        function c_mio_remesh_ex(h, opts, report) bind(c, name="mio_remesh_ex") result(r)
            import :: c_ptr, mio_remesh_opts_t, mio_remesh_report_t
            type(c_ptr), value :: h
            type(mio_remesh_opts_t), intent(in) :: opts
            type(mio_remesh_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        subroutine c_mio_remesh_volume_opts_init(opts) &
                bind(c, name="mio_remesh_volume_opts_init")
            import :: mio_remesh_volume_opts_t
            type(mio_remesh_volume_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_remesh_volume_ex(h, opts, report) &
                bind(c, name="mio_remesh_volume_ex") result(r)
            import :: c_ptr, mio_remesh_volume_opts_t, mio_remesh_volume_report_t
            type(c_ptr), value :: h
            type(mio_remesh_volume_opts_t), intent(in) :: opts
            type(mio_remesh_volume_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        function c_mio_split(h, by, tag_name) bind(c, name="mio_split") result(r)
            import :: c_ptr, c_char
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: by, tag_name
            type(c_ptr) :: r
        end function

        function c_mio_split_result_count(r) bind(c, name="mio_split_result_count") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t) :: n
        end function

        function c_mio_split_result_key(r, index, buf, buflen) &
                bind(c, name="mio_split_result_key") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(out) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_split_result_take_mesh(r, index) &
                bind(c, name="mio_split_result_take_mesh") result(m)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index
            type(c_ptr) :: m
        end function

        subroutine c_mio_split_result_free(r) bind(c, name="mio_split_result_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_convert_cells(h, mode, record_parent_ids) &
                bind(c, name="mio_convert_cells") result(r)
            import :: c_ptr, c_char, c_int
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: mode
            integer(c_int), value :: record_parent_ids
            type(c_ptr) :: r
        end function

        function c_mio_convert_cells_result_take_mesh(r) &
                bind(c, name="mio_convert_cells_result_take_mesh") result(m)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: m
        end function

        function c_mio_convert_cells_result_point_map(r, data, dtype, n) &
                bind(c, name="mio_convert_cells_result_point_map") result(s)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: r
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int64_t), intent(out) :: n
            integer(c_int) :: s
        end function

        subroutine c_mio_convert_cells_result_free(r) &
                bind(c, name="mio_convert_cells_result_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_subdivide(h, record_parent_ids) &
                bind(c, name="mio_subdivide") result(r)
            import :: c_ptr, c_int
            type(c_ptr), value :: h
            integer(c_int), value :: record_parent_ids
            type(c_ptr) :: r
        end function

        function c_mio_subdivide_result_take_mesh(r) &
                bind(c, name="mio_subdivide_result_take_mesh") result(m)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: m
        end function

        subroutine c_mio_subdivide_result_free(r) &
                bind(c, name="mio_subdivide_result_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_agglomerate_ex(h, opts, num_faces_merged, num_rejected) &
                bind(c, name="mio_agglomerate_ex") result(r)
            import :: c_ptr, c_int64_t, mio_agglomerate_opts_t
            type(c_ptr), value :: h
            type(mio_agglomerate_opts_t), intent(in) :: opts
            integer(c_int64_t), intent(out) :: num_faces_merged, num_rejected
            type(c_ptr) :: r
        end function

        function c_mio_blend_steps(a, b, w, blend_points) bind(c, name="mio_blend_steps") result(r)
            import :: c_ptr, c_double, c_int32_t
            type(c_ptr), value :: a, b
            real(c_double), value :: w
            integer(c_int32_t), value :: blend_points
            type(c_ptr) :: r
        end function

        function c_mio_agglomerate(h, target_group_size) &
                bind(c, name="mio_agglomerate") result(r)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: target_group_size
            type(c_ptr) :: r
        end function

        function c_mio_agglomerate_result_take_mesh(r) &
                bind(c, name="mio_agglomerate_result_take_mesh") result(m)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: m
        end function

        subroutine c_mio_agglomerate_result_free(r) &
                bind(c, name="mio_agglomerate_result_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_refine(h, levels, record_parent_ids) &
                bind(c, name="mio_refine") result(r)
            import :: c_ptr, c_int
            type(c_ptr), value :: h
            integer(c_int), value :: levels
            integer(c_int), value :: record_parent_ids
            type(c_ptr) :: r
        end function

        subroutine c_mio_sdf_opts_init(opts) bind(c, name="mio_sdf_opts_init")
            import :: mio_sdf_opts_t
            type(mio_sdf_opts_t), intent(out) :: opts
        end subroutine

        subroutine c_mio_voxel_opts_init(opts) bind(c, name="mio_voxel_opts_init")
            import :: mio_voxel_opts_t
            type(mio_voxel_opts_t), intent(out) :: opts
        end subroutine

        subroutine c_mio_octree_opts_init(opts) bind(c, name="mio_octree_opts_init")
            import :: mio_octree_opts_t
            type(mio_octree_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_octree(h, opts, report) bind(c, name="mio_octree") result(r)
            import :: c_ptr, mio_octree_opts_t, mio_octree_report_t
            type(c_ptr), value :: h
            type(mio_octree_opts_t), intent(in) :: opts
            type(mio_octree_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        subroutine c_mio_octree_mesh_opts_init(opts) bind(c, name="mio_octree_mesh_opts_init")
            import :: mio_octree_mesh_opts_t
            type(mio_octree_mesh_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_octree_mesh(h, opts, report) bind(c, name="mio_octree_mesh") result(r)
            import :: c_ptr, mio_octree_mesh_opts_t, mio_octree_mesh_report
            type(c_ptr), value :: h
            type(mio_octree_mesh_opts_t), intent(in) :: opts
            type(mio_octree_mesh_report), intent(out) :: report
            type(c_ptr) :: r
        end function

        function c_mio_octree_mesh_warnings(buf, buflen) &
                bind(c, name="mio_octree_mesh_warnings") result(n)
            import :: c_char, c_int64_t
            character(c_char), intent(out) :: buf(*)
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        function c_mio_grid(dims, origin, spacing, max_cells) bind(c, name="mio_grid") result(r)
            import :: c_ptr, c_int64_t, c_double
            integer(c_int64_t), intent(in) :: dims(*)
            real(c_double), intent(in) :: origin(*), spacing(*)
            integer(c_int64_t), value :: max_cells
            type(c_ptr) :: r
        end function

        function c_mio_rectangle(opts) bind(c, name="mio_rectangle") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_box(opts) bind(c, name="mio_box") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_disk(opts) bind(c, name="mio_disk") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_uv_sphere(opts) bind(c, name="mio_uv_sphere") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_icosa_sphere(opts) bind(c, name="mio_icosa_sphere") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_cylinder(opts) bind(c, name="mio_cylinder") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_reference_triangle(opts) bind(c, name="mio_reference_triangle") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_ngon(opts) bind(c, name="mio_ngon") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_annulus(opts) bind(c, name="mio_annulus") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_lshape(opts) bind(c, name="mio_lshape") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_plate_with_hole(opts) bind(c, name="mio_plate_with_hole") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_honeycomb(opts) bind(c, name="mio_honeycomb") result(r)
            import :: c_ptr, mio_gen_opts_t
            type(mio_gen_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_reference_cells(types, pairs) bind(c, name="mio_reference_cells") result(r)
            import :: c_ptr, c_char, c_int32_t
            character(kind=c_char), intent(in) :: types(*)
            integer(c_int32_t), value :: pairs
            type(c_ptr) :: r
        end function

        function c_mio_voxelize(h, opts, dims, origin, spacing, occupied) &
                bind(c, name="mio_voxelize") result(r)
            import :: c_ptr, c_int64_t, c_double, mio_voxel_opts_t
            type(c_ptr), value :: h
            type(mio_voxel_opts_t), intent(in) :: opts
            integer(c_int64_t), intent(out) :: dims(*)
            real(c_double), intent(out) :: origin(*), spacing(*)
            integer(c_int64_t), intent(out) :: occupied
            type(c_ptr) :: r
        end function

        subroutine c_mio_compute_sdf_opts_init(opts) &
                bind(c, name="mio_compute_sdf_opts_init")
            import :: mio_compute_sdf_opts_t
            type(mio_compute_sdf_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_compute_sdf(h, opts, dims, origin, spacing, max_depth, banded, q) &
                bind(c, name="mio_compute_sdf") result(r)
            import :: c_ptr, c_int64_t, c_double, mio_compute_sdf_opts_t, mio_surface_quality
            type(c_ptr), value :: h
            type(mio_compute_sdf_opts_t), intent(in) :: opts
            integer(c_int64_t), intent(out) :: dims(*)
            real(c_double), intent(out) :: origin(*), spacing(*)
            integer(c_int64_t), intent(out) :: max_depth, banded
            type(mio_surface_quality), intent(out) :: q
            type(c_ptr) :: r
        end function

        function c_mio_crop_predicate(h, array, compare, value, record_ids) &
                bind(c, name="mio_crop_predicate") result(r)
            import :: c_ptr, c_char, c_int, c_double
            type(c_ptr), value :: h
            character(kind=c_char), intent(in) :: array(*)
            integer(c_int), value :: compare, record_ids
            real(c_double), value :: value
            type(c_ptr) :: r
        end function

        function c_mio_surface_watertight_check(h, out) &
                bind(c, name="mio_surface_watertight_check") result(r)
            import :: c_ptr, c_int, mio_surface_quality
            type(c_ptr), value :: h
            type(mio_surface_quality), intent(out) :: out
            integer(c_int) :: r
        end function

        function c_mio_sample_distance(h, points, n_points, opts, out) &
                bind(c, name="mio_sample_distance") result(r)
            import :: c_ptr, c_int, c_int64_t, c_double, mio_sdf_opts_t
            type(c_ptr), value :: h
            real(c_double), intent(in) :: points(*)
            integer(c_int64_t), value :: n_points
            type(mio_sdf_opts_t), intent(in) :: opts
            real(c_double), intent(out) :: out(*)
            integer(c_int) :: r
        end function

        function c_mio_distance_to_surface(q, s, opts, banded, quality) &
                bind(c, name="mio_distance_to_surface") result(r)
            import :: c_ptr, c_int64_t, mio_sdf_opts_t, mio_surface_quality
            type(c_ptr), value :: q, s
            type(mio_sdf_opts_t), intent(in) :: opts
            integer(c_int64_t), intent(out) :: banded
            type(mio_surface_quality), intent(out) :: quality
            type(c_ptr) :: r
        end function

        subroutine c_mio_refine_opts_init(opts) bind(c, name="mio_refine_opts_init")
            import :: mio_refine_opts_t
            type(mio_refine_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_refine_ex(h, opts) bind(c, name="mio_refine_ex") result(r)
            import :: c_ptr, mio_refine_opts_t
            type(c_ptr), value :: h
            type(mio_refine_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        subroutine c_mio_tessellate_opts_init(opts) bind(c, name="mio_tessellate_opts_init")
            import :: mio_tessellate_opts_t
            type(mio_tessellate_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_tessellate(h, opts, report) bind(c, name="mio_tessellate") result(r)
            import :: c_ptr, mio_tessellate_opts_t, mio_tessellate_report_t
            type(c_ptr), value :: h
            type(mio_tessellate_opts_t), intent(in) :: opts
            type(mio_tessellate_report_t), intent(out) :: report
            type(c_ptr) :: r
        end function

        subroutine c_mio_smooth_opts_init(opts) bind(c, name="mio_smooth_opts_init")
            import :: mio_smooth_opts_t
            type(mio_smooth_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_smooth_ex(h, opts, nodes_moved, max_displacement, skipped_inversion) &
                bind(c, name="mio_smooth_ex") result(r)
            import :: c_ptr, c_int64_t, c_double, mio_smooth_opts_t
            type(c_ptr), value :: h
            type(mio_smooth_opts_t), intent(in) :: opts
            integer(c_int64_t), intent(out) :: nodes_moved, skipped_inversion
            real(c_double), intent(out) :: max_displacement
            type(c_ptr) :: r
        end function

        subroutine c_mio_decimate_opts_init(opts) bind(c, name="mio_decimate_opts_init")
            import :: mio_decimate_opts_t
            type(mio_decimate_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_decimate_ex(h, opts) bind(c, name="mio_decimate_ex") result(r)
            import :: c_ptr, mio_decimate_opts_t
            type(c_ptr), value :: h
            type(mio_decimate_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        function c_mio_refine_result_take_mesh(r) &
                bind(c, name="mio_refine_result_take_mesh") result(m)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: m
        end function

        function c_mio_refine_result_point_map(r, data, dtype, n) &
                bind(c, name="mio_refine_result_point_map") result(s)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: r
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int64_t), intent(out) :: n
            integer(c_int) :: s
        end function

        subroutine c_mio_refine_result_free(r) &
                bind(c, name="mio_refine_result_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_decimate(h, ratio, faces, max_error, placement, preserve_boundary, &
                                preserve_features, feature_angle) &
                bind(c, name="mio_decimate") result(r)
            import :: c_ptr, c_char, c_int, c_int64_t, c_double
            type(c_ptr), value :: h
            real(c_double), value :: ratio, max_error, feature_angle
            integer(c_int64_t), value :: faces
            character(kind=c_char), dimension(*), intent(in) :: placement
            integer(c_int), value :: preserve_boundary, preserve_features
            type(c_ptr) :: r
        end function

        function c_mio_decimate_result_take_mesh(r) &
                bind(c, name="mio_decimate_result_take_mesh") result(m)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: m
        end function

        function c_mio_decimate_result_point_map(r, data, dtype, n) &
                bind(c, name="mio_decimate_result_point_map") result(s)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: r
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int64_t), intent(out) :: n
            integer(c_int) :: s
        end function

        function c_mio_decimate_result_faces_removed(r) &
                bind(c, name="mio_decimate_result_faces_removed") result(v)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t) :: v
        end function

        function c_mio_decimate_result_points_removed(r) &
                bind(c, name="mio_decimate_result_points_removed") result(v)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t) :: v
        end function

        function c_mio_decimate_result_collapses_rejected(r) &
                bind(c, name="mio_decimate_result_collapses_rejected") result(v)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t) :: v
        end function

        function c_mio_decimate_result_max_error_applied(r) &
                bind(c, name="mio_decimate_result_max_error_applied") result(v)
            import :: c_ptr, c_double
            type(c_ptr), value :: r
            real(c_double) :: v
        end function

        subroutine c_mio_decimate_result_free(r) &
                bind(c, name="mio_decimate_result_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_partition(h, nparts, method, imbalance, mode, seed, record_ids, &
                                 ghost_layers, weights_key) &
                bind(c, name="mio_partition") result(r)
            import :: c_ptr, c_char, c_int, c_double
            type(c_ptr), value :: h
            integer(c_int), value :: nparts
            character(kind=c_char), dimension(*), intent(in) :: method
            real(c_double), value :: imbalance
            character(kind=c_char), dimension(*), intent(in) :: mode
            integer(c_int), value :: seed, record_ids, ghost_layers
            character(kind=c_char), dimension(*), intent(in) :: weights_key
            type(c_ptr) :: r
        end function

        function c_mio_partition_result_num_pieces(r) &
                bind(c, name="mio_partition_result_num_pieces") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t) :: n
        end function

        function c_mio_partition_result_take_mesh(r, index) &
                bind(c, name="mio_partition_result_take_mesh") result(m)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index
            type(c_ptr) :: m
        end function

        subroutine c_mio_partition_result_free(r) &
                bind(c, name="mio_partition_result_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_partition_labels(h, nparts, method, imbalance, mode, seed, &
                                        weights_key, labels, labels_size) &
                bind(c, name="mio_partition_labels") result(s)
            import :: c_ptr, c_char, c_int, c_int64_t, c_double
            type(c_ptr), value :: h
            integer(c_int), value :: nparts
            character(kind=c_char), dimension(*), intent(in) :: method
            real(c_double), value :: imbalance
            character(kind=c_char), dimension(*), intent(in) :: mode
            integer(c_int), value :: seed
            character(kind=c_char), dimension(*), intent(in) :: weights_key
            integer(c_int64_t), dimension(*), intent(out) :: labels
            integer(c_int64_t), value :: labels_size
            integer(c_int) :: s
        end function

        function c_mio_stats(h, out) bind(c, name="mio_stats") result(s)
            import :: c_ptr, c_int, mio_stats_report
            type(c_ptr), value :: h
            type(mio_stats_report), intent(out) :: out
            integer(c_int) :: s
        end function

        ! -- data operations --
        ! Name lists cross as a C array of char* plus an explicit count; the
        ! `names` argument is a c_ptr to that array (built by c_str_array).

        function c_mio_data_drop(h, location, names, count, ignore_missing) &
                bind(c, name="mio_data_drop") result(m)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: h, names
            integer(c_int), value :: location, ignore_missing
            integer(c_int64_t), value :: count
            type(c_ptr) :: m
        end function

        function c_mio_data_keep(h, location, names, count, ignore_missing) &
                bind(c, name="mio_data_keep") result(m)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: h, names
            integer(c_int), value :: location, ignore_missing
            integer(c_int64_t), value :: count
            type(c_ptr) :: m
        end function

        function c_mio_data_rename(h, location, from_name, to_name) &
                bind(c, name="mio_data_rename") result(m)
            import :: c_ptr, c_int, c_char
            type(c_ptr), value :: h
            integer(c_int), value :: location
            character(kind=c_char), dimension(*), intent(in) :: from_name, to_name
            type(c_ptr) :: m
        end function

        function c_mio_sets_to_data(h, location, data_name, join_char, names, count) &
                bind(c, name="mio_sets_to_data") result(m)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: h, data_name, join_char, names
            integer(c_int), value :: location
            integer(c_int64_t), value :: count
            type(c_ptr) :: m
        end function

        function c_mio_data_to_sets(h, location, key) bind(c, name="mio_data_to_sets") result(m)
            import :: c_ptr, c_int, c_char
            type(c_ptr), value :: h
            integer(c_int), value :: location
            character(kind=c_char), dimension(*), intent(in) :: key
            type(c_ptr) :: m
        end function

        function c_mio_data_point_to_cell(h, names, count, suffix) &
                bind(c, name="mio_data_point_to_cell") result(m)
            import :: c_ptr, c_int64_t, c_char
            type(c_ptr), value :: h, names
            integer(c_int64_t), value :: count
            character(kind=c_char), dimension(*), intent(in) :: suffix
            type(c_ptr) :: m
        end function

        function c_mio_data_cell_to_point(h, names, count, weight, suffix) &
                bind(c, name="mio_data_cell_to_point") result(m)
            import :: c_ptr, c_int, c_int64_t, c_char
            type(c_ptr), value :: h, names
            integer(c_int64_t), value :: count
            integer(c_int), value :: weight
            character(kind=c_char), dimension(*), intent(in) :: suffix
            type(c_ptr) :: m
        end function

        function c_mio_data_calc(h, expression, location, output_name, overwrite) &
                bind(c, name="mio_data_calc") result(m)
            import :: c_ptr, c_int, c_char
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: expression, output_name
            integer(c_int), value :: location, overwrite
            type(c_ptr) :: m
        end function

        function c_mio_data_condition(h, location, names, count, mode, lo, hi, scope, &
                                      nan_policy, nan_replacement, suffix) &
                bind(c, name="mio_data_condition") result(m)
            import :: c_ptr, c_int, c_int64_t, c_double, c_char
            type(c_ptr), value :: h, names
            integer(c_int), value :: location, mode, scope, nan_policy
            integer(c_int64_t), value :: count
            real(c_double), value :: lo, hi, nan_replacement
            character(kind=c_char), dimension(*), intent(in) :: suffix
            type(c_ptr) :: m
        end function

        function c_mio_tensor_invariants(h, location, names, count, outputs, prefix, suffix, &
                                         overwrite) &
                bind(c, name="mio_tensor_invariants") result(m)
            import :: c_ptr, c_int, c_int64_t, c_char
            type(c_ptr), value :: h, names
            integer(c_int), value :: location, outputs, overwrite
            integer(c_int64_t), value :: count
            character(kind=c_char), dimension(*), intent(in) :: prefix, suffix
            type(c_ptr) :: m
        end function

        function c_mio_data_info_create(h) bind(c, name="mio_data_info_create") result(r)
            import :: c_ptr
            type(c_ptr), value :: h
            type(c_ptr) :: r
        end function

        function c_mio_data_info_count(r) bind(c, name="mio_data_info_count") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t) :: n
        end function

        function c_mio_data_info_name(r, index, buf, buflen) &
                bind(c, name="mio_data_info_name") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(out) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_data_info_entry(r, index, out) &
                bind(c, name="mio_data_info_entry") result(s)
            import :: c_ptr, c_int, c_int64_t, mio_data_array_info
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index
            type(mio_data_array_info), intent(out) :: out
            integer(c_int) :: s
        end function

        function c_mio_regions_create(h) bind(c, name="mio_regions_create") result(r)
            import :: c_ptr
            type(c_ptr), value :: h
            type(c_ptr) :: r
        end function

        function c_mio_regions_count(r) bind(c, name="mio_regions_count") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t) :: n
        end function

        function c_mio_regions_name(r, index, buf, buflen) &
                bind(c, name="mio_regions_name") result(n)
            import :: c_ptr, c_int64_t, c_char
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), intent(inout) :: buf(*)
            integer(c_int64_t) :: n
        end function

        function c_mio_regions_info(r, index, out) &
                bind(c, name="mio_regions_info") result(s)
            import :: c_ptr, c_int64_t, c_int, mio_region_info
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index
            type(mio_region_info), intent(out) :: out
            integer(c_int) :: s
        end function

        function c_mio_regions_entries(r, index, count) &
                bind(c, name="mio_regions_entries") result(p)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index
            integer(c_int64_t), intent(out) :: count
            type(c_ptr) :: p
        end function

        subroutine c_mio_regions_free(r) bind(c, name="mio_regions_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_mesh_add_region(h, name, kind, dim, tag, entries, count) &
                bind(c, name="mio_mesh_add_region") result(s)
            import :: c_ptr, c_char, c_int, c_int64_t
            type(c_ptr), value :: h
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), value :: kind, dim
            integer(c_int64_t), value :: tag, count
            integer(c_int64_t), intent(in) :: entries(*)
            integer(c_int) :: s
        end function

        subroutine c_mio_data_info_free(r) bind(c, name="mio_data_info_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_data_integrate_create(h, names, count) &
                bind(c, name="mio_data_integrate_create") result(r)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h, names
            integer(c_int64_t), value :: count
            type(c_ptr) :: r
        end function

        function c_mio_data_integrate_count(r) &
                bind(c, name="mio_data_integrate_count") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t) :: n
        end function

        function c_mio_data_integrate_name(r, index, buf, buflen) &
                bind(c, name="mio_data_integrate_name") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(out) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_data_integrate_entry(r, index, out) &
                bind(c, name="mio_data_integrate_entry") result(s)
            import :: c_ptr, c_int, c_int64_t, mio_field_integral_info
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index
            type(mio_field_integral_info), intent(out) :: out
            integer(c_int) :: s
        end function

        function c_mio_data_integrate_component(r, index, comp, total, mean, domain_measure, &
                num_nan) bind(c, name="mio_data_integrate_component") result(s)
            import :: c_ptr, c_int, c_int64_t, c_double
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index, comp
            real(c_double), intent(out) :: total, mean, domain_measure
            integer(c_int64_t), intent(out) :: num_nan
            integer(c_int) :: s
        end function

        function c_mio_data_integrate_region_count(r, index) &
                bind(c, name="mio_data_integrate_region_count") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index
            integer(c_int64_t) :: n
        end function

        function c_mio_data_integrate_region_name(r, index, region, buf, buflen) &
                bind(c, name="mio_data_integrate_region_name") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index, region, buflen
            character(kind=c_char), dimension(*), intent(out) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_data_integrate_region_entry(r, index, region, out) &
                bind(c, name="mio_data_integrate_region_entry") result(s)
            import :: c_ptr, c_int, c_int64_t, mio_field_integral_info
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index, region
            type(mio_field_integral_info), intent(out) :: out
            integer(c_int) :: s
        end function

        function c_mio_data_integrate_region_component(r, index, region, comp, total, mean, &
                domain_measure, num_nan) &
                bind(c, name="mio_data_integrate_region_component") result(s)
            import :: c_ptr, c_int, c_int64_t, c_double
            type(c_ptr), value :: r
            integer(c_int64_t), value :: index, region, comp
            real(c_double), intent(out) :: total, mean, domain_measure
            integer(c_int64_t), intent(out) :: num_nan
            integer(c_int) :: s
        end function

        subroutine c_mio_data_integrate_free(r) bind(c, name="mio_data_integrate_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_reorder_result_take_mesh(r) &
                bind(c, name="mio_reorder_result_take_mesh") result(m)
            import :: c_ptr
            type(c_ptr), value :: r
            type(c_ptr) :: m
        end function

        function c_mio_reorder_result_node_perm(r, data, dtype, n) &
                bind(c, name="mio_reorder_result_node_perm") result(s)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: r
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int64_t), intent(out) :: n
            integer(c_int) :: s
        end function

        subroutine c_mio_reorder_result_free(r) bind(c, name="mio_reorder_result_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_compute_bandwidth(h) bind(c, name="mio_compute_bandwidth") result(bw)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: bw
        end function

        function c_mio_meshes_equal(a, b, atol, rtol, unordered, out_equal) &
                bind(c, name="mio_meshes_equal") result(s)
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: a, b
            real(c_double), value :: atol, rtol
            integer(c_int), value :: unordered
            integer(c_int), intent(out) :: out_equal
            integer(c_int) :: s
        end function

        function c_mio_diff(a, b, atol, rtol, unordered, out) &
                bind(c, name="mio_diff") result(s)
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: a, b
            real(c_double), value :: atol, rtol
            integer(c_int), value :: unordered
            type(c_ptr), intent(out) :: out
            integer(c_int) :: s
        end function

        function c_mio_diff_result_verdict(r) &
                bind(c, name="mio_diff_result_verdict") result(v)
            import :: c_ptr, c_int
            type(c_ptr), value :: r
            integer(c_int) :: v
        end function

        subroutine c_mio_diff_result_free(r) bind(c, name="mio_diff_result_free")
            import :: c_ptr
            type(c_ptr), value :: r
        end subroutine

        function c_mio_sniff_format(path, buf, buflen) &
                bind(c, name="mio_sniff_format") result(n)
            import :: c_char, c_int64_t
            character(kind=c_char), dimension(*), intent(in) :: path
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_set_points(h, dtype, num_points, dim, xyz) &
                bind(c, name="mio_mesh_set_points") result(s)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: h, xyz
            integer(c_int), value :: dtype
            integer(c_int64_t), value :: num_points, dim
            integer(c_int) :: s
        end function

        function c_mio_mesh_add_cell_block(h, cell_type, num_cells, nodes_per_cell, dtype, conn) &
                bind(c, name="mio_mesh_add_cell_block") result(s)
            import :: c_ptr, c_char, c_int, c_int64_t
            type(c_ptr), value :: h, conn
            character(kind=c_char), dimension(*), intent(in) :: cell_type
            integer(c_int64_t), value :: num_cells, nodes_per_cell
            integer(c_int), value :: dtype
            integer(c_int) :: s
        end function

        function c_mio_mesh_add_point_data(h, name, dtype, ndim, shape, data) &
                bind(c, name="mio_mesh_add_point_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h, data
            character(kind=c_char), dimension(*), intent(in) :: name
            integer(c_int), value :: dtype
            integer(c_int32_t), value :: ndim
            integer(c_int64_t), dimension(*), intent(in) :: shape
            integer(c_int) :: s
        end function

        function c_mio_mesh_append_cell_data(h, name, dtype, ndim, shape, data) &
                bind(c, name="mio_mesh_append_cell_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h, data
            character(kind=c_char), dimension(*), intent(in) :: name
            integer(c_int), value :: dtype
            integer(c_int32_t), value :: ndim
            integer(c_int64_t), dimension(*), intent(in) :: shape
            integer(c_int) :: s
        end function

        function c_mio_mesh_add_field_data(h, name, dtype, ndim, shape, data) &
                bind(c, name="mio_mesh_add_field_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h, data
            character(kind=c_char), dimension(*), intent(in) :: name
            integer(c_int), value :: dtype
            integer(c_int32_t), value :: ndim
            integer(c_int64_t), dimension(*), intent(in) :: shape
            integer(c_int) :: s
        end function

        function c_mio_mesh_num_points(h) bind(c, name="mio_mesh_num_points") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_point_dim(h) bind(c, name="mio_mesh_point_dim") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_get_points(h, data, dtype) &
                bind(c, name="mio_mesh_get_points") result(s)
            import :: c_ptr, c_int
            type(c_ptr), value :: h
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int) :: s
        end function

        function c_mio_mesh_num_cell_blocks(h) bind(c, name="mio_mesh_num_cell_blocks") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_cell_block_info(h, block, num_cells, nodes_per_cell, is_ragged) &
                bind(c, name="mio_mesh_cell_block_info") result(s)
            import :: c_ptr, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: block
            integer(c_int64_t), intent(out) :: num_cells, nodes_per_cell
            integer(c_int32_t), intent(out) :: is_ragged
            integer(c_int) :: s
        end function

        function c_mio_mesh_cell_block_info_ex(h, block, out) &
                bind(c, name="mio_mesh_cell_block_info_ex") result(s)
            import :: c_ptr, c_int, c_int64_t, c_mio_cell_block_info
            type(c_ptr), value :: h
            integer(c_int64_t), value :: block
            type(c_mio_cell_block_info), intent(out) :: out
            integer(c_int) :: s
        end function

        function c_mio_mesh_add_polygon_block(h, ctype, num_cells, row_offsets, nodes, &
                                              num_nodes) &
                bind(c, name="mio_mesh_add_polygon_block") result(s)
            import :: c_ptr, c_char, c_int, c_int64_t
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: ctype
            integer(c_int64_t), value :: num_cells, num_nodes
            integer(c_int64_t), dimension(*), intent(in) :: row_offsets, nodes
            integer(c_int) :: s
        end function

        function c_mio_mesh_add_polyhedron_block(h, ctype, num_cells, cell_offsets, num_faces, &
                                                 face_offsets, nodes, num_nodes) &
                bind(c, name="mio_mesh_add_polyhedron_block") result(s)
            import :: c_ptr, c_char, c_int, c_int64_t
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: ctype
            integer(c_int64_t), value :: num_cells, num_faces, num_nodes
            integer(c_int64_t), dimension(*), intent(in) :: cell_offsets, face_offsets, nodes
            integer(c_int) :: s
        end function

        function c_mio_poly_conn_create(h, block) bind(c, name="mio_poly_conn_create") result(p)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: block
            type(c_ptr) :: p
        end function

        function c_mio_poly_conn_get_shape(p, out) &
                bind(c, name="mio_poly_conn_get_shape") result(s)
            import :: c_ptr, c_int, c_mio_poly_conn_shape
            type(c_ptr), value :: p
            type(c_mio_poly_conn_shape), intent(out) :: out
            integer(c_int) :: s
        end function

        function c_mio_poly_conn_nodes(p, count) bind(c, name="mio_poly_conn_nodes") result(q)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: p
            integer(c_int64_t), intent(out) :: count
            type(c_ptr) :: q
        end function

        function c_mio_poly_conn_face_offsets(p, count) &
                bind(c, name="mio_poly_conn_face_offsets") result(q)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: p
            integer(c_int64_t), intent(out) :: count
            type(c_ptr) :: q
        end function

        function c_mio_poly_conn_cell_offsets(p, count) &
                bind(c, name="mio_poly_conn_cell_offsets") result(q)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: p
            integer(c_int64_t), intent(out) :: count
            type(c_ptr) :: q
        end function

        subroutine c_mio_poly_conn_free(p) bind(c, name="mio_poly_conn_free")
            import :: c_ptr
            type(c_ptr), value :: p
        end subroutine

        function c_mio_mesh_cell_block_type(h, block, buf, buflen) &
                bind(c, name="mio_mesh_cell_block_type") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: block, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_cell_block_conn(h, block, conn, dtype) &
                bind(c, name="mio_mesh_cell_block_conn") result(s)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: block
            type(c_ptr), intent(out) :: conn
            integer(c_int), intent(out) :: dtype
            integer(c_int) :: s
        end function

        function c_mio_mesh_num_point_data(h) bind(c, name="mio_mesh_num_point_data") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_point_data_name(h, index, buf, buflen) &
                bind(c, name="mio_mesh_point_data_name") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_get_point_data(h, name, data, dtype, ndim, shape) &
                bind(c, name="mio_mesh_get_point_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: name
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), dimension(*), intent(out) :: shape
            integer(c_int) :: s
        end function

        function c_mio_mesh_num_cell_data(h) bind(c, name="mio_mesh_num_cell_data") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_cell_data_name(h, index, buf, buflen) &
                bind(c, name="mio_mesh_cell_data_name") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_cell_data_num_blocks(h, name) &
                bind(c, name="mio_mesh_cell_data_num_blocks") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: name
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_get_cell_data(h, name, block, data, dtype, ndim, shape) &
                bind(c, name="mio_mesh_get_cell_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: name
            integer(c_int64_t), value :: block
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), dimension(*), intent(out) :: shape
            integer(c_int) :: s
        end function

        function c_mio_mesh_num_field_data(h) bind(c, name="mio_mesh_num_field_data") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_field_data_name(h, index, buf, buflen) &
                bind(c, name="mio_mesh_field_data_name") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_get_field_data(h, name, data, dtype, ndim, shape) &
                bind(c, name="mio_mesh_get_field_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: name
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), dimension(*), intent(out) :: shape
            integer(c_int) :: s
        end function

        ! -- transient XDMF series --

        function c_mio_exodus_series_create(path) bind(c, name="mio_exodus_series_create") result(h)
            import :: c_char, c_ptr
            character(c_char), intent(in) :: path(*)
            type(c_ptr) :: h
        end function
        function c_mio_exodus_series_write_points_cells(s, mesh) &
                bind(c, name="mio_exodus_series_write_points_cells") result(st)
            import :: c_ptr, c_int
            type(c_ptr), value :: s, mesh
            integer(c_int) :: st
        end function
        function c_mio_exodus_series_write_data(s, time, mesh) &
                bind(c, name="mio_exodus_series_write_data") result(st)
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: s, mesh
            real(c_double), value :: time
            integer(c_int) :: st
        end function
        function c_mio_exodus_series_flush(s) bind(c, name="mio_exodus_series_flush") result(st)
            import :: c_ptr, c_int
            type(c_ptr), value :: s
            integer(c_int) :: st
        end function
        function c_mio_exodus_series_finalize(s) bind(c, name="mio_exodus_series_finalize") result(st)
            import :: c_ptr, c_int
            type(c_ptr), value :: s
            integer(c_int) :: st
        end function
        function c_mio_exodus_series_finalized(s) bind(c, name="mio_exodus_series_finalized") result(f)
            import :: c_ptr, c_int32_t
            type(c_ptr), value :: s
            integer(c_int32_t) :: f
        end function
        function c_mio_exodus_series_num_steps(s) bind(c, name="mio_exodus_series_num_steps") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: s
            integer(c_int64_t) :: n
        end function
        subroutine c_mio_exodus_series_free(s) bind(c, name="mio_exodus_series_free")
            import :: c_ptr
            type(c_ptr), value :: s
        end subroutine

        function c_mio_xdmf_series_create(path, data_format, gzip_level) &
                bind(c, name="mio_xdmf_series_create") result(h)
            import :: c_ptr, c_char, c_int32_t
            character(kind=c_char), dimension(*), intent(in) :: path, data_format
            integer(c_int32_t), value :: gzip_level
            type(c_ptr) :: h
        end function

        function c_mio_xdmf_series_write_points_cells(s, mesh) &
                bind(c, name="mio_xdmf_series_write_points_cells") result(st)
            import :: c_ptr, c_int
            type(c_ptr), value :: s, mesh
            integer(c_int) :: st
        end function

        function c_mio_xdmf_series_write_data(s, time, mesh) &
                bind(c, name="mio_xdmf_series_write_data") result(st)
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: s, mesh
            real(c_double), value :: time
            integer(c_int) :: st
        end function

        function c_mio_xdmf_series_finalize(s) &
                bind(c, name="mio_xdmf_series_finalize") result(st)
            import :: c_ptr, c_int
            type(c_ptr), value :: s
            integer(c_int) :: st
        end function

        function c_mio_xdmf_series_flush(s) &
                bind(c, name="mio_xdmf_series_flush") result(st)
            import :: c_ptr, c_int
            type(c_ptr), value :: s
            integer(c_int) :: st
        end function

        function c_mio_xdmf_series_finalized(s) &
                bind(c, name="mio_xdmf_series_finalized") result(f)
            import :: c_ptr, c_int32_t
            type(c_ptr), value :: s
            integer(c_int32_t) :: f
        end function

        function c_mio_xdmf_series_create_ex(path, opts) &
                bind(c, name="mio_xdmf_series_create_ex") result(h)
            import :: c_ptr, c_char, mio_xdmf_series_opts_t
            character(kind=c_char), dimension(*), intent(in) :: path
            type(mio_xdmf_series_opts_t), intent(in) :: opts
            type(c_ptr) :: h
        end function

        function c_mio_xdmf_series_num_steps(s) &
                bind(c, name="mio_xdmf_series_num_steps") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: s
            integer(c_int64_t) :: n
        end function

        subroutine c_mio_xdmf_series_free(s) bind(c, name="mio_xdmf_series_free")
            import :: c_ptr
            type(c_ptr), value :: s
        end subroutine

        function c_strlen(s) bind(c, name="strlen") result(n)
            import :: c_ptr, c_size_t
            type(c_ptr), value :: s
            integer(c_size_t) :: n
        end function

        ! Provenance (see doc/provenance.md). A scope cannot cross the ABI as
        ! RAII, so the C surface is a begin/end pair; Fortran keeps that shape
        ! explicitly, matching `mio_mesh`'s own explicit-free convention.
        function c_mio_provenance_scope_begin(mode) &
            bind(c, name="mio_provenance_scope_begin") result(r)
            import :: c_int
            integer(c_int), value :: mode
            integer(c_int) :: r
        end function

        subroutine c_mio_provenance_scope_end() bind(c, name="mio_provenance_scope_end")
        end subroutine

        subroutine c_mio_provenance_note(cat, det) bind(c, name="mio_provenance_note")
            import :: c_char
            character(kind=c_char), dimension(*), intent(in) :: cat, det
        end subroutine

        subroutine c_mio_provenance_set_source(path, fmt) &
            bind(c, name="mio_provenance_set_source")
            import :: c_char
            character(kind=c_char), dimension(*), intent(in) :: path, fmt
        end subroutine

        subroutine c_mio_provenance_set_target(fmt, enc, codec, ff) &
            bind(c, name="mio_provenance_set_target")
            import :: c_char
            character(kind=c_char), dimension(*), intent(in) :: fmt, enc, codec, ff
        end subroutine
        ! -- format side channels (mio_format_info) --
        function c_mio_read_with_info(path, format, opts, info) &
                bind(c, name="mio_read_with_info") result(h)
            import :: c_ptr, c_char
            character(kind=c_char), dimension(*), intent(in) :: path, format
            type(c_ptr), value :: opts
            type(c_ptr), intent(out) :: info
            type(c_ptr) :: h
        end function

        function c_mio_write_with_info(path, mesh, format, info) &
                bind(c, name="mio_write_with_info") result(st)
            import :: c_ptr, c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: path, format
            type(c_ptr), value :: mesh, info
            integer(c_int) :: st
        end function

        function c_mio_format_info_format(info, buf, buflen) &
                bind(c, name="mio_format_info_format") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: info
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        subroutine c_mio_format_info_free(info) bind(c, name="mio_format_info_free")
            import :: c_ptr
            type(c_ptr), value :: info
        end subroutine

        function c_mio_mdpa_info_count(info, section) &
                bind(c, name="mio_mdpa_info_count") result(n)
            import :: c_ptr, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section
            integer(c_int64_t) :: n
        end function

        function c_mio_gmsh_info_count(info, section) &
                bind(c, name="mio_gmsh_info_count") result(n)
            import :: c_ptr, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section
            integer(c_int64_t) :: n
        end function

        function c_mio_gmsh_info_array(info, section, index, field, data, dtype, ndim, shape) &
                bind(c, name="mio_gmsh_info_array") result(st)
            import :: c_ptr, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section, field
            integer(c_int64_t), value :: index
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), intent(out) :: shape(*)
            integer(c_int) :: st
        end function

        function c_mio_nurbs_info_count(info, section) &
                bind(c, name="mio_nurbs_info_count") result(n)
            import :: c_ptr, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section
            integer(c_int64_t) :: n
        end function

        function c_mio_nurbs_info_int(info, section, index, field, value) &
                bind(c, name="mio_nurbs_info_int") result(st)
            import :: c_ptr, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section, field
            integer(c_int64_t), value :: index
            integer(c_int64_t), intent(out) :: value
            integer(c_int) :: st
        end function

        function c_mio_nurbs_info_array(info, section, index, field, data, dtype, ndim, shape) &
                bind(c, name="mio_nurbs_info_array") result(st)
            import :: c_ptr, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section, field
            integer(c_int64_t), value :: index
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), intent(out) :: shape(*)
            integer(c_int) :: st
        end function

        function c_mio_nurbs_info_string(info, section, index, field, buf, buflen) &
                bind(c, name="mio_nurbs_info_string") result(n)
            import :: c_ptr, c_char, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section, field
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        subroutine c_mio_iga_opts_init(opts) bind(c, name="mio_iga_opts_init")
            import :: mio_iga_opts_t
            type(mio_iga_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_iga_refine(info, opts, refined) bind(c, name="mio_iga_refine") result(h)
            import :: c_ptr, mio_iga_opts_t
            type(c_ptr), value :: info
            type(mio_iga_opts_t), intent(in) :: opts
            type(c_ptr), intent(out) :: refined
            type(c_ptr) :: h
        end function

        subroutine c_mio_mesh_cad_surface_opts_init(opts) &
                bind(c, name="mio_mesh_cad_surface_opts_init")
            import :: mio_mesh_cad_surface_opts_t
            type(mio_mesh_cad_surface_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_mesh_cad_surface(info, opts, report) &
                bind(c, name="mio_mesh_cad_surface") result(h)
            import :: c_ptr, mio_mesh_cad_surface_opts_t, mio_mesh_cad_surface_report
            type(c_ptr), value :: info
            type(mio_mesh_cad_surface_opts_t), intent(in) :: opts
            type(mio_mesh_cad_surface_report), intent(out) :: report
            type(c_ptr) :: h
        end function

        function c_mio_mesh_cad_surface2(info, opts, opts_ex, report) &
                bind(c, name="mio_mesh_cad_surface2") result(h)
            import :: c_ptr, mio_mesh_cad_surface_opts_t, mio_mesh_cad_surface_opts_ex_t, &
                      mio_mesh_cad_surface_report
            type(c_ptr), value :: info
            type(mio_mesh_cad_surface_opts_t), intent(in) :: opts
            type(mio_mesh_cad_surface_opts_ex_t), intent(in) :: opts_ex
            type(mio_mesh_cad_surface_report), intent(out) :: report
            type(c_ptr) :: h
        end function

        function c_mio_mesh_cad_surface_num_reparametrised() &
                bind(c, name="mio_mesh_cad_surface_num_reparametrised") result(m)
            import :: c_int64_t
            integer(c_int64_t) :: m
        end function

        function c_mio_mesh_cad_surface_reparametrised(idx, face, cond, cond_disc, method) &
                bind(c, name="mio_mesh_cad_surface_reparametrised") result(st)
            import :: c_int64_t, c_int32_t, c_double, c_int
            integer(c_int64_t), value :: idx
            integer(c_int64_t), intent(out) :: face
            real(c_double), intent(out) :: cond, cond_disc
            integer(c_int32_t), intent(out) :: method
            integer(c_int) :: st
        end function

        function c_mio_mesh_cad_surface_num_merged_faces() &
                bind(c, name="mio_mesh_cad_surface_num_merged_faces") result(m)
            import :: c_int64_t
            integer(c_int64_t) :: m
        end function

        function c_mio_mesh_cad_surface_merged_face(idx, sliver, host, width, merged) &
                bind(c, name="mio_mesh_cad_surface_merged_face") result(st)
            import :: c_int64_t, c_int32_t, c_double, c_int
            integer(c_int64_t), value :: idx
            integer(c_int64_t), intent(out) :: sliver, host
            real(c_double), intent(out) :: width
            integer(c_int32_t), intent(out) :: merged
            integer(c_int) :: st
        end function

        function c_mio_mesh_cad_surface_num_tolerance_uses() &
                bind(c, name="mio_mesh_cad_surface_num_tolerance_uses") result(m)
            import :: c_int64_t
            integer(c_int64_t) :: m
        end function

        function c_mio_mesh_cad_surface_tolerance_use(idx, edge, declared, gap, where) &
                bind(c, name="mio_mesh_cad_surface_tolerance_use") result(st)
            import :: c_int64_t, c_int32_t, c_double, c_int
            integer(c_int64_t), value :: idx
            integer(c_int64_t), intent(out) :: edge
            real(c_double), intent(out) :: declared, gap
            integer(c_int32_t), intent(out) :: where
            integer(c_int) :: st
        end function

        function c_mio_mesh_cad_surface_num_lopped_corners() &
                bind(c, name="mio_mesh_cad_surface_num_lopped_corners") result(m)
            import :: c_int64_t
            integer(c_int64_t) :: m
        end function

        function c_mio_mesh_cad_surface_lopped_corner(idx, face, point, angle, lfs0, radius, &
                                                      turning, bound, forced, min_forced, &
                                                      within, guaranteed) &
                bind(c, name="mio_mesh_cad_surface_lopped_corner") result(st)
            import :: c_int64_t, c_int32_t, c_double, c_int
            integer(c_int64_t), value :: idx
            integer(c_int64_t), intent(out) :: face, point, forced
            real(c_double), intent(out) :: angle, lfs0, radius, turning, bound, min_forced
            integer(c_int32_t), intent(out) :: within, guaranteed
            integer(c_int) :: st
        end function

        function c_mio_mesh_cad_surface_lopped_corner_reason(idx, buf, buflen) &
                bind(c, name="mio_mesh_cad_surface_lopped_corner_reason") result(m)
            import :: c_int64_t, c_char
            integer(c_int64_t), value :: idx
            character(kind=c_char), intent(out) :: buf(*)
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: m
        end function

        function c_mio_mesh_cad_surface_face_size_limits(face, edges, n) &
                bind(c, name="mio_mesh_cad_surface_face_size_limits") result(m)
            import :: c_int64_t
            integer(c_int64_t), value :: face
            integer(c_int64_t), intent(out) :: edges(*)
            integer(c_int64_t), value :: n
            integer(c_int64_t) :: m
        end function

        function c_mio_mesh_cad_surface_face_status(face) &
                bind(c, name="mio_mesh_cad_surface_face_status") result(st)
            import :: c_int32_t, c_int64_t
            integer(c_int64_t), value :: face
            integer(c_int32_t) :: st
        end function

        function c_mio_mesh_cad_surface_face_reason(face, buf, buflen) &
                bind(c, name="mio_mesh_cad_surface_face_reason") result(n)
            import :: c_char, c_int64_t
            integer(c_int64_t), value :: face
            character(c_char), intent(out) :: buf(*)
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_cad_surface_warnings(buf, buflen) &
                bind(c, name="mio_mesh_cad_surface_warnings") result(n)
            import :: c_char, c_int64_t
            character(c_char), intent(out) :: buf(*)
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_cad_surface_angle_histogram(bins, n) &
                bind(c, name="mio_mesh_cad_surface_angle_histogram") result(m)
            import :: c_int64_t
            integer(c_int64_t), intent(out) :: bins(*)
            integer(c_int64_t), value :: n
            integer(c_int64_t) :: m
        end function

        function c_mio_mesh_cad_surface_num_corners() &
                bind(c, name="mio_mesh_cad_surface_num_corners") result(n)
            import :: c_int64_t
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_cad_surface_corner(idx, corner) &
                bind(c, name="mio_mesh_cad_surface_corner") result(st)
            import :: c_int, c_int64_t, mio_mesh_cad_corner
            integer(c_int64_t), value :: idx
            type(mio_mesh_cad_corner), intent(out) :: corner
            integer(c_int) :: st
        end function

        function c_mio_iga_sample(info, opts, mapped) bind(c, name="mio_iga_sample") result(h)
            import :: c_ptr, mio_iga_opts_t
            type(c_ptr), value :: info
            type(mio_iga_opts_t), intent(in) :: opts
            type(c_ptr), intent(out) :: mapped
            type(c_ptr) :: h
        end function

        function c_mio_iga_evaluate(info, patch, params, num_points, derivatives, trim_samples) &
                bind(c, name="mio_iga_evaluate") result(h)
            import :: c_ptr, c_double, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int64_t), value :: patch, num_points
            real(c_double), intent(in) :: params(*)
            integer(c_int32_t), value :: derivatives, trim_samples
            type(c_ptr) :: h
        end function

        function c_mio_iga_evaluation_num_fields(ev) &
                bind(c, name="mio_iga_evaluation_num_fields") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: ev
            integer(c_int64_t) :: n
        end function

        function c_mio_iga_evaluation_array(ev, which, index, data, dtype, ndim, shape) &
                bind(c, name="mio_iga_evaluation_array") result(st)
            import :: c_ptr, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: ev
            integer(c_int32_t), value :: which
            integer(c_int64_t), value :: index
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), intent(out) :: shape(*)
            integer(c_int) :: st
        end function

        subroutine c_mio_iga_evaluation_free(ev) bind(c, name="mio_iga_evaluation_free")
            import :: c_ptr
            type(c_ptr), value :: ev
        end subroutine

        function c_mio_brep_info_count(info, section) &
                bind(c, name="mio_brep_info_count") result(n)
            import :: c_ptr, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section
            integer(c_int64_t) :: n
        end function

        function c_mio_brep_info_int(info, section, index, field, value) &
                bind(c, name="mio_brep_info_int") result(st)
            import :: c_ptr, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section, field
            integer(c_int64_t), value :: index
            integer(c_int64_t), intent(out) :: value
            integer(c_int) :: st
        end function

        function c_mio_brep_info_array(info, section, index, field, data, dtype, ndim, shape) &
                bind(c, name="mio_brep_info_array") result(st)
            import :: c_ptr, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section, field
            integer(c_int64_t), value :: index
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), intent(out) :: shape(*)
            integer(c_int) :: st
        end function

        function c_mio_brep_info_string(info, section, index, field, buf, buflen) &
                bind(c, name="mio_brep_info_string") result(n)
            import :: c_ptr, c_char, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section, field
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_mdpa_info_string(info, section, index, field, buf, buflen) &
                bind(c, name="mio_mdpa_info_string") result(n)
            import :: c_ptr, c_char, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section, field
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_mdpa_info_int(info, section, index, field, value) &
                bind(c, name="mio_mdpa_info_int") result(st)
            import :: c_ptr, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section, field
            integer(c_int64_t), value :: index
            integer(c_int64_t), intent(out) :: value
            integer(c_int) :: st
        end function

        function c_mio_mdpa_info_array(info, section, index, field, data, dtype, ndim, shape) &
                bind(c, name="mio_mdpa_info_array") result(st)
            import :: c_ptr, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section, field
            integer(c_int64_t), value :: index
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), intent(out) :: shape(*)
            integer(c_int) :: st
        end function

        function c_mio_mdpa_info_data_count(info, section, index) &
                bind(c, name="mio_mdpa_info_data_count") result(n)
            import :: c_ptr, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section
            integer(c_int64_t), value :: index
            integer(c_int64_t) :: n
        end function

        function c_mio_mdpa_info_data_kind(info, section, index, entry) &
                bind(c, name="mio_mdpa_info_data_kind") result(k)
            import :: c_ptr, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section
            integer(c_int64_t), value :: index, entry
            integer(c_int32_t) :: k
        end function

        function c_mio_mdpa_info_data_string(info, section, index, entry, field, buf, buflen) &
                bind(c, name="mio_mdpa_info_data_string") result(n)
            import :: c_ptr, c_char, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section, field
            integer(c_int64_t), value :: index, entry, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_mdpa_info_data_array(info, section, index, entry, data, dtype, ndim, &
                                            shape) bind(c, name="mio_mdpa_info_data_array") &
                result(st)
            import :: c_ptr, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: info
            integer(c_int32_t), value :: section
            integer(c_int64_t), value :: index, entry
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), intent(out) :: shape(*)
            integer(c_int) :: st
        end function
        subroutine c_mio_render_opts_init(opts) bind(c, name="mio_render_opts_init")
            import :: mio_render_opts_t
            type(mio_render_opts_t), intent(out) :: opts
        end subroutine

        function c_mio_render(h, opts) bind(c, name="mio_render") result(r)
            import :: c_ptr, mio_render_opts_t
            type(c_ptr), value :: h
            type(mio_render_opts_t), intent(in) :: opts
            type(c_ptr) :: r
        end function

        subroutine c_mio_frame_free(f) bind(c, name="mio_frame_free")
            import :: c_ptr
            type(c_ptr), value :: f
        end subroutine

        function c_mio_frame_width(f) bind(c, name="mio_frame_width") result(n)
            import :: c_ptr, c_int32_t
            type(c_ptr), value :: f
            integer(c_int32_t) :: n
        end function

        function c_mio_frame_height(f) bind(c, name="mio_frame_height") result(n)
            import :: c_ptr, c_int32_t
            type(c_ptr), value :: f
            integer(c_int32_t) :: n
        end function

        function c_mio_frame_rgba(f) bind(c, name="mio_frame_rgba") result(p)
            import :: c_ptr
            type(c_ptr), value :: f
            type(c_ptr) :: p
        end function

        function c_mio_frame_cell_ids(f) bind(c, name="mio_frame_cell_ids") result(p)
            import :: c_ptr
            type(c_ptr), value :: f
            type(c_ptr) :: p
        end function

        function c_mio_frame_range(f, vmin, vmax) bind(c, name="mio_frame_range") result(n)
            import :: c_ptr, c_int32_t, c_double
            type(c_ptr), value :: f
            real(c_double), intent(out) :: vmin, vmax
            integer(c_int32_t) :: n
        end function

        function c_mio_frame_num_notes(f) bind(c, name="mio_frame_num_notes") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: f
            integer(c_int64_t) :: n
        end function

        function c_mio_frame_note(f, index, buf, buflen) bind(c, name="mio_frame_note") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: f
            integer(c_int64_t), value :: index, buflen
            character(c_char), intent(out) :: buf(*)
            integer(c_int64_t) :: n
        end function

        function c_mio_frame_text(f, opts, buf, buflen) bind(c, name="mio_frame_text") result(n)
            import :: c_ptr, c_char, c_int64_t, mio_text_opts_t
            type(c_ptr), value :: f
            type(mio_text_opts_t), intent(in) :: opts
            character(c_char), intent(out) :: buf(*)
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        function c_mio_frame_png(f, compress, buf, buflen) bind(c, name="mio_frame_png") result(n)
            import :: c_ptr, c_int32_t, c_int64_t, c_int8_t
            type(c_ptr), value :: f
            integer(c_int32_t), value :: compress
            integer(c_int8_t), intent(out) :: buf(*)
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        function c_mio_render_text(h, ropts, topts, buf, buflen) bind(c, name="mio_render_text") result(n)
            import :: c_ptr, c_char, c_int64_t, mio_render_opts_t, mio_text_opts_t
            type(c_ptr), value :: h
            type(mio_render_opts_t), intent(in) :: ropts
            type(mio_text_opts_t), intent(in) :: topts
            character(c_char), intent(out) :: buf(*)
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        function c_mio_write_snapshot(path, h, ropts, topts, sopts) bind(c, name="mio_write_snapshot") result(s)
            import :: c_ptr, c_char, c_int, mio_render_opts_t, mio_text_opts_t, mio_snapshot_opts_t
            character(c_char), intent(in) :: path(*)
            type(c_ptr), value :: h
            type(mio_render_opts_t), intent(in) :: ropts
            type(mio_text_opts_t), intent(in) :: topts
            type(mio_snapshot_opts_t), intent(in) :: sopts
            integer(c_int) :: s
        end function

        function c_mio_detect_color_depth(no_color, color_term, term) &
                bind(c, name="mio_detect_color_depth") result(n)
            import :: c_ptr, c_int32_t
            type(c_ptr), value :: no_color, color_term, term
            integer(c_int32_t) :: n
        end function

    end interface

contains

    ! ------------------------------------------------------------------
    ! String / error helpers
    ! ------------------------------------------------------------------

    !> NUL-terminated copy for passing to C ("" stays "" + NUL, meaning
    !> "infer format" on the C side).
    pure function c_str(f) result(c)
        character(*), intent(in) :: f
        character(kind=c_char, len=:), allocatable :: c
        c = trim(f)//c_null_char
    end function

    !> Fortran string from a NUL-terminated C pointer (static storage).
    function c_ptr_to_string(p) result(s)
        type(c_ptr), intent(in) :: p
        character(:), allocatable :: s
        character(kind=c_char), pointer :: chars(:)
        integer :: n, i
        s = ''
        if (.not. c_associated(p)) return
        n = int(c_strlen(p))
        if (n <= 0) return
        call c_f_pointer(p, chars, [n])
        s = repeat(' ', n)  ! reallocation on assignment (s starts as '')
        do i = 1, n
            s(i:i) = chars(i)
        end do
    end function

    !> Build the `char**` a C name-list parameter expects.
    !>
    !> `storage` holds NUL-terminated copies of every name (one column each) and
    !> `cptrs` holds a c_ptr to each column. BOTH must stay in scope for the
    !> duration of the C call -- they are the actual backing memory, so the
    !> caller keeps them as local variables until the call returns.
    !> `names` may be zero-sized, in which case the returned pointer is NULL and
    !> the count is 0 ("every array at that location").
    subroutine c_str_array(names, storage, cptrs, arr, count)
        character(*), intent(in) :: names(:)
        character(kind=c_char), allocatable, target, intent(out) :: storage(:, :)
        type(c_ptr), allocatable, target, intent(out) :: cptrs(:)
        type(c_ptr), intent(out) :: arr
        integer(c_int64_t), intent(out) :: count
        integer :: i, j, n, width
        n = size(names)
        count = int(n, c_int64_t)
        if (n == 0) then
            allocate (storage(1, 1))
            allocate (cptrs(1))
            arr = c_null_ptr
            return
        end if
        width = len(names) + 1
        allocate (storage(width, n))
        allocate (cptrs(n))
        do i = 1, n
            do j = 1, len_trim(names(i))
                storage(j, i) = names(i) (j:j)
            end do
            storage(len_trim(names(i)) + 1, i) = c_null_char
            cptrs(i) = c_loc(storage(1, i))
        end do
        arr = c_loc(cptrs(1))
    end subroutine

    !> Fortran string from the first `n` chars of a C char buffer.
    function from_c_buf(buf, n) result(s)
        character(kind=c_char), intent(in) :: buf(*)
        integer, intent(in) :: n
        character(:), allocatable :: s
        integer :: i
        allocate (character(max(n, 0)) :: s)
        do i = 1, n
            s(i:i) = buf(i)
        end do
    end function

    !> The failure message of the most recent failed meshio++ call on this
    !> thread ('' if none).
    function mio_error_message() result(msg)
        character(:), allocatable :: msg
        msg = c_ptr_to_string(c_mio_last_error())
    end function

    !> Map a C status to the optional stat/errmsg pair; no stat + failure =
    !> print and error stop.
    subroutine handle_status(status, what, stat, errmsg)
        integer(c_int), intent(in) :: status
        character(*), intent(in) :: what
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (present(stat)) stat = int(status)
        if (present(errmsg)) errmsg = ''
        if (status /= 0_c_int) then
            if (present(errmsg)) errmsg = mio_error_message()
            if (.not. present(stat)) then
                write (error_unit, '(a)') 'meshio++ ('//what//'): '//mio_error_message()
                error stop 1
            end if
        end if
    end subroutine

    !> Report a Fortran-side failure (bad handle, shape mismatch, ...)
    !> through the same stat/errmsg protocol.
    subroutine handle_failure(what, msg, stat, errmsg)
        character(*), intent(in) :: what, msg
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (present(stat)) stat = 1
        if (present(errmsg)) errmsg = msg
        if (.not. present(stat)) then
            write (error_unit, '(a)') 'meshio++ ('//what//'): '//msg
            error stop 1
        end if
    end subroutine

    subroutine clear_status(stat, errmsg)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (present(stat)) stat = 0
        if (present(errmsg)) errmsg = ''
    end subroutine

    ! ------------------------------------------------------------------
    ! Module-level procedures
    ! ------------------------------------------------------------------

    !> The meshio++ version string, e.g. "6.1.0".
    function mio_version() result(v)
        character(:), allocatable :: v
        v = c_ptr_to_string(c_mio_version())
    end function

    !> The compile-time mesh backend: "meshio", "native", or "kratos".
    function mio_mesh_backend() result(b)
        character(:), allocatable :: b
        b = c_ptr_to_string(c_mio_mesh_backend())
    end function

    !> .true. if `format` (e.g. "gmsh", "vtu", "med") is readable in this build.
    function mio_format_readable(format) result(r)
        character(*), intent(in) :: format
        logical :: r
        r = c_mio_format_readable(c_str(format)) /= 0_c_int
    end function

    !> .true. if `format` is writable in this build.
    function mio_format_writable(format) result(r)
        character(*), intent(in) :: format
        logical :: r
        r = c_mio_format_writable(c_str(format)) /= 0_c_int
    end function

    !> Read `in_path` and immediately write it to `out_path` (the CLI's
    !> `convert`). Formats are inferred from the extensions unless given.
    subroutine mio_convert(in_path, out_path, in_format, out_format, stat, errmsg)
        character(*), intent(in) :: in_path, out_path
        character(*), intent(in), optional :: in_format, out_format
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: ifmt, ofmt
        ifmt = ''; if (present(in_format)) ifmt = in_format
        ofmt = ''; if (present(out_format)) ofmt = out_format
        call handle_status(c_mio_convert(c_str(in_path), c_str(ifmt), c_str(out_path), &
                                         c_str(ofmt)), 'convert', stat, errmsg)
    end subroutine

    !> Enumerate MED mesh names without geometry. Names are dynamically sized.
    subroutine mio_med_mesh_names(path, names, stat, errmsg)
        character(*), intent(in) :: path
        character(:), allocatable, intent(out) :: names(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(c_char), allocatable :: buf(:)
        integer(c_int64_t) :: n, len_name, longest, i
        n = c_mio_med_mesh_count(c_str(path))
        if (n < 0) then
            call handle_failure('med_mesh_names', mio_error_message(), stat, errmsg)
            return
        end if
        longest = 0
        allocate(buf(1))
        do i = 0, n - 1
            len_name = c_mio_med_mesh_name(c_str(path), i, buf, 0_c_int64_t)
            if (len_name < 0) then
                call handle_failure('med_mesh_names', mio_error_message(), stat, errmsg)
                return
            end if
            longest = max(longest, len_name)
        end do
        deallocate(buf)
        allocate(character(len=int(longest)) :: names(int(n)))
        allocate(buf(int(longest) + 1))
        do i = 0, n - 1
            len_name = c_mio_med_mesh_name(c_str(path), i, buf, longest + 1)
            if (len_name < 0) then
                call handle_failure('med_mesh_names', mio_error_message(), stat, errmsg)
                return
            end if
            names(int(i) + 1) = from_c_buf(buf, int(len_name))
        end do
        call clear_status(stat, errmsg)
    end subroutine

    !> Read a named MED mesh, replacing output only after a successful read.
    subroutine mio_med_read_named(path, name, mesh, time_step, lenient, stat, errmsg)
        character(*), intent(in) :: path, name
        type(mio_mesh), intent(inout) :: mesh
        integer, intent(in), optional :: time_step
        logical, intent(in), optional :: lenient
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_read_opts_t) :: opts
        type(c_ptr) :: h
        call c_mio_read_opts_init(opts)
        if (present(time_step)) opts%time_step = int(time_step, c_int64_t)
        if (present(lenient)) then
            if (lenient) opts%lenient = 1
        end if
        h = c_mio_med_read_named(c_str(path), c_str(name), opts)
        if (.not. c_associated(h)) then
            call handle_failure('med_read_named', mio_error_message(), stat, errmsg)
            return
        end if
        call mesh%free()
        mesh%handle = h
        call clear_status(stat, errmsg)
    end subroutine

    !> Write several named meshes to one MED file; no input handle is modified.
    subroutine mio_med_write_multi(path, meshes, names, version, stat, errmsg)
        character(*), intent(in) :: path, names(:)
        type(mio_mesh), intent(in) :: meshes(:)
        character(*), intent(in), optional :: version
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: ver
        character(kind=c_char, len=:), allocatable, target :: bufs(:)
        type(c_ptr), allocatable :: ptrs(:), handles(:)
        integer :: i
        if (size(meshes) /= size(names) .or. size(meshes) == 0) then
            call handle_failure('med_write_multi', 'provide one name per mesh', stat, errmsg)
            return
        end if
        ver = '4.1.0'; if (present(version)) ver = version
        allocate(character(kind=c_char, len=len(names) + 1) :: bufs(size(names)))
        allocate(ptrs(size(names)), handles(size(names)))
        do i = 1, size(names)
            bufs(i) = trim(names(i))//c_null_char
            ptrs(i) = c_loc(bufs(i)(1:1))
            handles(i) = meshes(i)%handle
        end do
        call handle_status(c_mio_med_write_multi(c_str(path), handles, ptrs, &
            int(size(meshes), c_int64_t), c_str(ver)), 'med_write_multi', stat, errmsg)
    end subroutine

    !> Run a whole settings.json pipeline (read -> operation chain -> write;
    !> PascalCase vocabulary, see doc/pipeline.md). Needs a build with the
    !> JSON parser (-DMESHIOPLUSPLUS_WITH_JSON=ON); otherwise the error names
    !> the flag. mio_pipeline_has_json() reports which build this is.
    subroutine mio_pipeline_run_file(settings_path, stat, errmsg)
        character(*), intent(in) :: settings_path
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call handle_status(c_mio_pipeline_run_file(c_str(settings_path)), &
                           'pipeline_run_file', stat, errmsg)
    end subroutine

    !> Run a settings pipeline given as JSON text (see mio_pipeline_run_file).
    subroutine mio_pipeline_run_json(json_text, stat, errmsg)
        character(*), intent(in) :: json_text
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call handle_status(c_mio_pipeline_run_json(c_str(json_text)), &
                           'pipeline_run_json', stat, errmsg)
    end subroutine

    !> Copy an owning report's JSON into an allocatable string, then release it.
    subroutine pipeline_copy_report(h, report, stat, errmsg)
        type(c_ptr), value :: h
        character(:), allocatable, intent(out) :: report
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(c_char), allocatable :: buf(:)
        integer(c_int64_t) :: n
        report = ''
        if (.not. c_associated(h)) then
            call handle_failure('pipeline_report', mio_error_message(), stat, errmsg)
            return
        end if
        allocate(buf(1))
        n = c_mio_pipeline_report_json(h, buf, 0_c_int64_t)
        if (n < 0 .or. n >= int(huge(0) - 1, c_int64_t)) then
            call c_mio_pipeline_report_free(h)
            call handle_failure('pipeline_report', 'invalid report length', stat, errmsg)
            return
        end if
        deallocate(buf)
        allocate(buf(int(n) + 1))
        n = c_mio_pipeline_report_json(h, buf, n + 1)
        if (n < 0) then
            call c_mio_pipeline_report_free(h)
            call handle_failure('pipeline_report', mio_error_message(), stat, errmsg)
            return
        end if
        report = from_c_buf(buf, int(n))
        call c_mio_pipeline_report_free(h)
        call clear_status(stat, errmsg)
    end subroutine

    subroutine mio_pipeline_run_file_report(path, report, stat, errmsg)
        character(*), intent(in) :: path
        character(:), allocatable, intent(out) :: report
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call pipeline_copy_report(c_mio_pipeline_run_file_report(c_str(path)), report, stat, errmsg)
    end subroutine

    subroutine mio_pipeline_run_json_report(text, report, stat, errmsg)
        character(*), intent(in) :: text
        character(:), allocatable, intent(out) :: report
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call pipeline_copy_report(c_mio_pipeline_run_json_report(c_str(text)), report, stat, errmsg)
    end subroutine

    subroutine mio_sequence_pipeline_run_file_report(path, report, stat, errmsg)
        character(*), intent(in) :: path
        character(:), allocatable, intent(out) :: report
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call pipeline_copy_report(c_mio_sequence_pipeline_run_file_report(c_str(path)), report, stat, errmsg)
    end subroutine

    subroutine mio_sequence_pipeline_run_json_report(text, report, stat, errmsg)
        character(*), intent(in) :: text
        character(:), allocatable, intent(out) :: report
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call pipeline_copy_report(c_mio_sequence_pipeline_run_json_report(c_str(text)), report, stat, errmsg)
    end subroutine

    !> Whether this build carries the JSON pipeline parser.
    !> Whether the library carries the OpenCASCADE backend fit_cad needs.
    function mio_has_occt() result(has)
        logical :: has
        has = c_mio_has_occt() /= 0
    end function

    function mio_pipeline_has_json() result(has)
        logical :: has
        has = c_mio_pipeline_has_json() /= 0
    end function

    !> Guess a mesh file's format from its contents ('' if undetermined).
    function mio_sniff_format(path) result(fmt)
        character(*), intent(in) :: path
        character(:), allocatable :: fmt
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: n
        fmt = ''
        n = c_mio_sniff_format(c_str(path), buf, int(STRBUF_LEN, c_int64_t))
        if (n > 0) fmt = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
    end function

    ! ------------------------------------------------------------------
    ! mio_mesh: lifecycle & file I/O
    ! ------------------------------------------------------------------

    !> Allocate an empty mesh (read() does this implicitly).
    subroutine mesh_create(self, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call mesh_free(self)
        self%handle = c_mio_mesh_create()
        if (.not. c_associated(self%handle)) then
            call handle_failure('create', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end subroutine

    !> Release the mesh and every pointer borrowed from it. Idempotent.
    subroutine mesh_free(self)
        class(mio_mesh), intent(inout) :: self
        if (c_associated(self%handle)) call c_mio_mesh_free(self%handle)
        self%handle = c_null_ptr
    end subroutine

    !> .true. between a successful create()/read() and free().
    logical function mesh_is_valid(self)
        class(mio_mesh), intent(in) :: self
        mesh_is_valid = c_associated(self%handle)
    end function

    !> Read a mesh file, replacing any previous content of this handle.
    !>
    !> `points_only` skips every data array; `arrays` keeps only the named ones
    !> (a zero-sized `arrays` keeps none, which is distinct from omitting the
    !> argument entirely -- that reads everything). Formats without a native
    !> selective path are read whole and filtered, so the result is the same
    !> either way; only the cost differs.
    subroutine mesh_read(self, path, format, points_only, arrays, time_step, lenient, &
                         stat, errmsg, piece, drop_ghosts)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: path
        character(*), intent(in), optional :: format
        logical, intent(in), optional :: points_only
        character(*), intent(in), optional :: arrays(:)
        !> Which step of a multi-step file to read: 0 (default) is the first,
        !> negative counts from the end. Out of range fails, never clamps.
        integer, intent(in), optional :: time_step
        !> Downgrade "this reader cannot represent construct X" to a warning and
        !> a skip. Not "ignore all errors": a malformed file still fails.
        logical, intent(in), optional :: lenient
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        !> Keep only this piece of a partitioned file (VTKHDF partitions or
        !> composite blocks): 0 is the first, negative counts from the end.
        !> Omitted merges every piece into one mesh with one cell region per
        !> piece. Out of range fails, never clamps. Last in the argument list so
        !> a caller passing `stat`/`errmsg` positionally keeps working.
        integer, intent(in), optional :: piece
        !> Remove the ghost cells (halo) of a partitioned .pvtu/.pvtp/.pvd; every
        !> other reader ignores it. Omitted keeps them. Last in the argument list
        !> for the same reason as `piece`.
        logical, intent(in), optional :: drop_ghosts
        character(:), allocatable :: fmt
        type(c_ptr) :: h
        type(mio_read_opts_t) :: opts
        ! NUL-terminated copies must outlive the call, so they are held here
        ! rather than built inline; c_loc needs them contiguous and TARGET.
        character(kind=c_char, len=STRBUF_LEN), allocatable, target :: bufs(:)
        type(c_ptr), allocatable, target :: ptrs(:)
        integer :: i, n

        fmt = ''; if (present(format)) fmt = format

        if (.not. present(points_only) .and. .not. present(arrays) &
            .and. .not. present(time_step) .and. .not. present(lenient) &
            .and. .not. present(piece) .and. .not. present(drop_ghosts)) then
            h = c_mio_read(c_str(path), c_str(fmt))
        else
            call c_mio_read_opts_init(opts)
            if (present(points_only)) then
                if (points_only) opts%points_only = 1
            end if
            if (present(time_step)) opts%time_step = int(time_step, c_int64_t)
            if (present(lenient)) then
                if (lenient) opts%lenient = 1
            end if
            if (present(piece)) then
                opts%piece = int(piece, c_int64_t)
                opts%piece_set = 1
            end if
            if (present(drop_ghosts)) then
                if (drop_ghosts) opts%drop_ghosts = 1
            end if
            if (present(arrays)) then
                n = size(arrays)
                allocate (bufs(max(n, 1)))
                allocate (ptrs(max(n, 1)))
                do i = 1, n
                    bufs(i) = trim(arrays(i))//c_null_char
                    ptrs(i) = c_loc(bufs(i)(1:1))
                end do
                ! Non-null pointer with count 0 means "no arrays" -- the C side
                ! distinguishes that from a null pointer ("every array").
                opts%arrays = c_loc(ptrs(1))
                opts%num_arrays = int(n, c_int64_t)
            end if
            h = c_mio_read_ex(c_str(path), c_str(fmt), opts)
        end if

        if (.not. c_associated(h)) then
            call handle_failure('read', mio_error_message(), stat, errmsg)
            return
        end if
        call mesh_free(self)
        self%handle = h
        call clear_status(stat, errmsg)
    end subroutine

    !> Summarize a mesh file without loading its heavy arrays.
    function mio_read_metadata(path, format, stat, errmsg) result(meta)
        character(*), intent(in) :: path
        character(*), intent(in), optional :: format
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_metadata) :: meta
        character(:), allocatable :: fmt
        type(c_ptr) :: h
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: nblocks, i, n, ncells, npc, nsteps, nregions
        integer(c_int) :: ragged, s
        type(mio_region_info) :: rinfo

        fmt = ''; if (present(format)) fmt = format
        h = c_mio_read_metadata_create(c_str(path), c_str(fmt))
        if (.not. c_associated(h)) then
            call handle_failure('read_metadata', mio_error_message(), stat, errmsg)
            return
        end if

        meta%num_points = c_mio_read_metadata_num_points(h)
        meta%point_dim = c_mio_read_metadata_point_dim(h)
        meta%num_cells = c_mio_read_metadata_num_cells(h)
        meta%fell_back_to_full_read = (c_mio_read_metadata_fell_back(h) == 1)

        nblocks = c_mio_read_metadata_num_cell_blocks(h)
        allocate (meta%cell_blocks(max(nblocks, 0_c_int64_t)))
        do i = 1, nblocks
            s = c_mio_read_metadata_cell_block(h, i - 1_c_int64_t, ncells, npc, ragged)
            if (s /= 0) cycle
            meta%cell_blocks(i)%num_cells = ncells
            meta%cell_blocks(i)%nodes_per_cell = npc
            meta%cell_blocks(i)%ragged = (ragged /= 0)
            n = c_mio_read_metadata_cell_block_type(h, i - 1_c_int64_t, buf, &
                                                    int(STRBUF_LEN, c_int64_t))
            if (n > 0) meta%cell_blocks(i)%cell_type = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
        end do

        call metadata_read_names(h, MIO_DATA_POINT, meta%point_data_names)
        call metadata_read_names(h, MIO_DATA_CELL, meta%cell_data_names)
        call metadata_read_names(h, MIO_DATA_FIELD, meta%field_data_names)

        nsteps = c_mio_read_metadata_num_time_values(h)
        if (nsteps < 0) nsteps = 0
        allocate (meta%time_values(nsteps))
        if (nsteps > 0) n = c_mio_read_metadata_time_values(h, meta%time_values, nsteps)

        nregions = c_mio_read_metadata_num_regions(h)
        if (nregions < 0) nregions = 0
        allocate (meta%regions(nregions))
        do i = 1, nregions
            s = c_mio_read_metadata_region_info(h, i - 1_c_int64_t, rinfo)
            if (s /= 0) cycle
            meta%regions(i)%kind = rinfo%kind
            meta%regions(i)%dim = rinfo%dim
            meta%regions(i)%tag = rinfo%tag
            meta%regions(i)%num_entries = rinfo%num_entries
            n = c_mio_read_metadata_region_name(h, i - 1_c_int64_t, buf, &
                                                int(STRBUF_LEN, c_int64_t))
            if (n > 0) meta%regions(i)%name = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
        end do

        call c_mio_read_metadata_free(h)
        call clear_status(stat, errmsg)
    end function

    !> Fill `names` with the data-array names at `location`.
    subroutine metadata_read_names(h, location, names)
        type(c_ptr), intent(in) :: h
        integer(c_int), intent(in) :: location
        character(len=STRBUF_LEN), allocatable, intent(out) :: names(:)
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: count, i, n
        count = c_mio_read_metadata_num_names(h, location)
        if (count < 0) count = 0
        allocate (names(count))
        do i = 1, count
            names(i) = ''
            n = c_mio_read_metadata_name(h, location, i - 1_c_int64_t, buf, &
                                          int(STRBUF_LEN, c_int64_t))
            if (n > 0) names(i) = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
        end do
    end subroutine

    !> Write the mesh to a file.
    subroutine mesh_write(self, path, format, stat, errmsg, encoding, codec, float_format)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: path
        character(*), intent(in), optional :: format
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(*), intent(in), optional :: encoding, codec, float_format
        character(:), allocatable :: fmt
        character(kind=c_char), allocatable, target :: ff(:)
        type(mio_write_opts_c) :: opts
        fmt = ''; if (present(format)) fmt = format
        if (present(encoding)) then
            select case (encoding)
            case ('default'); opts%encoding = 0
            case ('ascii'); opts%encoding = 1
            case ('binary'); opts%encoding = 2
            case ('raw_appended'); opts%encoding = 3
            case default
                call handle_failure('write', 'invalid write encoding: '//encoding, stat, errmsg)
                return
            end select
        end if
        if (present(codec)) then
            select case (codec)
            case ('default'); opts%codec = 0
            case ('none'); opts%codec = 1
            case ('zlib'); opts%codec = 2
            case ('lz4'); opts%codec = 3
            case ('zstd'); opts%codec = 4
            case ('lzf'); opts%codec = 5
            case default
                call handle_failure('write', 'invalid write codec: '//codec, stat, errmsg)
                return
            end select
        end if
        if (present(float_format)) then
            ff = c_str(float_format)
            opts%float_format = c_loc(ff(1))
        end if
        call handle_status(c_mio_write_ex(c_str(path), self%handle, c_str(fmt), opts), 'write', &
                           stat, errmsg)
    end subroutine

    subroutine mesh_write_gltf(self, path, options, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: path
        type(mio_gltf_options), intent(in), optional :: options
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_gltf_opts_c) :: opts
        character(kind=c_char), allocatable, target :: color(:), cmap(:), nan_color(:)
        if (present(options)) then
            opts%container = options%container
            opts%up_axis = options%up_axis
            opts%normal_weight = options%normal_weight
            opts%normals = merge(1_c_int, 0_c_int, options%normals)
            opts%fields = merge(1_c_int, 0_c_int, options%fields)
            opts%recenter = merge(1_c_int, 0_c_int, options%recenter)
            opts%by_region = merge(1_c_int, 0_c_int, options%by_region)
            opts%unlit = merge(1_c_int, 0_c_int, options%unlit)
            opts%component = options%component - 1
            opts%component_set = merge(1_c_int, 0_c_int, options%component_set)
            opts%vmin_set = merge(1_c_int, 0_c_int, options%vmin_set)
            opts%vmax_set = merge(1_c_int, 0_c_int, options%vmax_set)
            opts%split_angle = options%split_angle
            opts%scale = options%scale
            opts%vmin = options%vmin
            opts%vmax = options%vmax
            if (allocated(options%color_by)) then
                color = c_str(options%color_by)
                opts%color_by = c_loc(color(1))
            end if
            if (allocated(options%cmap)) then
                cmap = c_str(options%cmap)
                opts%cmap = c_loc(cmap(1))
            end if
            if (allocated(options%nan_color)) then
                nan_color = c_str(options%nan_color)
                opts%nan_color = c_loc(nan_color(1))
            end if
        end if
        if (c_sizeof(opts) /= 168_c_size_t) then
            call handle_failure('write_gltf', 'glTF option layout mismatch', stat, errmsg)
            return
        end if
        call handle_status(c_mio_write_gltf(c_str(path), self%handle, opts), 'write_gltf', stat, errmsg)
    end subroutine

    !> Extract the boundary of the mesh's highest-dimension cells as a new mesh
    !> (volume -> faces, surface -> edges).
    function mesh_extract_surface(self, record_parent_ids, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        logical, intent(in), optional :: record_parent_ids
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: rec
        rec = 0
        if (present(record_parent_ids)) then
            if (record_parent_ids) rec = 1
        end if
        out%handle = c_mio_extract_surface(self%handle, rec)
        if (.not. c_associated(out%handle)) then
            call handle_failure('extract_surface', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Extract the boundary skin of a volume mesh as a new surface mesh.
    function mesh_extract_skin(self, linearize, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        logical, intent(in), optional :: linearize
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: lin
        lin = 0
        if (present(linearize)) then
            if (linearize) lin = 1
        end if
        out%handle = c_mio_extract_skin(self%handle, lin)
        if (.not. c_associated(out%handle)) then
            call handle_failure('extract_skin', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> A copy of the mesh with per-cell quality metrics attached as cell_data.
    function mesh_attach_quality(self, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        out%handle = c_mio_attach_quality(self%handle)
        if (.not. c_associated(out%handle)) then
            call handle_failure('attach_quality', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Aggregate quality counts (total / inverted / degenerate cells).
    subroutine mesh_quality_counts(self, num_cells, num_inverted, num_degenerate, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        integer(int64), intent(out), optional :: num_cells, num_inverted, num_degenerate
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t) :: nc, ninv, ndeg
        if (c_mio_quality_counts(self%handle, nc, ninv, ndeg) /= 0_c_int) then
            call handle_failure('quality_counts', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(num_cells)) num_cells = nc
        if (present(num_inverted)) num_inverted = ninv
        if (present(num_degenerate)) num_degenerate = ndeg
        call clear_status(stat, errmsg)
    end subroutine

    !> Merge two or more meshes into one (concatenate, with optional welding of
    !> coincident nodes within `atol`). `data_policy` is "intersection" (default)
    !> or "fill". point_sets/cell_sets are not carried across the C ABI.
    function mio_merge(meshes, weld, atol, source_tag, data_policy, &
                       drop_duplicate_cells, stat, errmsg) result(out)
        type(mio_mesh), intent(in) :: meshes(:)
        logical, intent(in), optional :: weld, source_tag, drop_duplicate_cells
        real(real64), intent(in), optional :: atol
        character(*), intent(in), optional :: data_policy
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr), allocatable :: handles(:)
        integer(c_int) :: cweld, csrc, cpolicy, cdrop
        real(c_double) :: catol
        integer :: i, n
        n = size(meshes)
        allocate (handles(max(n, 1)))
        do i = 1, n
            handles(i) = meshes(i)%handle
        end do
        cweld = 0
        if (present(weld)) then
            if (weld) cweld = 1
        end if
        csrc = 1
        if (present(source_tag)) then
            if (.not. source_tag) csrc = 0
        end if
        cdrop = 0
        if (present(drop_duplicate_cells)) then
            if (drop_duplicate_cells) cdrop = 1
        end if
        cpolicy = 0
        if (present(data_policy)) then
            if (data_policy == 'fill') cpolicy = 1
        end if
        catol = 1.0e-8_c_double
        if (present(atol)) catol = real(atol, c_double)
        out%handle = c_mio_merge(handles, int(n, c_int64_t), cweld, catol, csrc, cpolicy, cdrop)
        if (.not. c_associated(out%handle)) then
            call handle_failure('merge', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Sample data arrays from `source` onto `target` (cross-mesh field
    !> transfer): source point_data at the target's points, source cell_data by
    !> nearest source-cell centroid regardless of the method. `method` is
    !> 'nearest' (default) or 'barycentric' (linear in a simplexified source;
    !> a target point outside the source domain receives `default_value`
    !> unless `extrapolate`). `arrays` omitted or zero-sized means every source
    !> point_data array; cell_data transfers only when named. `on_conflict` is
    !> 'error' (default), 'overwrite' or 'suffix' (name + '_interp').
    function mio_interpolate(source, target, method, arrays, extrapolate, &
                             default_value, on_conflict, stat, errmsg) result(out)
        type(mio_mesh), intent(in) :: source, target
        character(*), intent(in), optional :: method, on_conflict
        character(*), intent(in), optional :: arrays(:)
        logical, intent(in), optional :: extrapolate
        real(real64), intent(in), optional :: default_value
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr
        integer(c_int64_t) :: count
        integer(c_int) :: cextrap
        real(c_double) :: cdefault
        character(:), allocatable :: cmethod, cconflict
        cmethod = 'nearest'
        if (present(method)) cmethod = trim(method)
        cconflict = 'error'
        if (present(on_conflict)) cconflict = trim(on_conflict)
        cextrap = 0
        if (present(extrapolate)) then
            if (extrapolate) cextrap = 1
        end if
        cdefault = 0.0_c_double
        if (present(default_value)) cdefault = real(default_value, c_double)
        if (present(arrays)) then
            call c_str_array(arrays, storage, cptrs, arr, count)
        else
            arr = c_null_ptr
            count = 0_c_int64_t
        end if
        out%handle = c_mio_interpolate(source%handle, target%handle, c_str(cmethod), arr, &
                                       count, cextrap, cdefault, c_str(cconflict))
        if (.not. c_associated(out%handle)) then
            call handle_failure('interpolate', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Mass-preserving cross-mesh field transfer: an exact overlap-measure
    !> weighted remap, so sum(target value * target measure) equals
    !> sum(source value * source measure) over the shared region -- the
    !> property mio_interpolate's 'barycentric' mode does not have. Both
    !> meshes are simplexified (accepting ragged/polyhedron blocks for free).
    !> `arrays` omitted or zero-sized means every source point_data AND
    !> cell_data array (one algorithm regardless of location, unlike
    !> mio_interpolate). `on_conflict` is 'error' (default), 'overwrite' or
    !> 'suffix' (name + '_interp'). Output arrays are always Float64.
    function mio_conservative_interpolate(source, target, arrays, default_value, &
                                          on_conflict, stat, errmsg) result(out)
        type(mio_mesh), intent(in) :: source, target
        character(*), intent(in), optional :: on_conflict
        character(*), intent(in), optional :: arrays(:)
        real(real64), intent(in), optional :: default_value
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr
        integer(c_int64_t) :: count
        real(c_double) :: cdefault
        character(:), allocatable :: cconflict
        cconflict = 'error'
        if (present(on_conflict)) cconflict = trim(on_conflict)
        cdefault = 0.0_c_double
        if (present(default_value)) cdefault = real(default_value, c_double)
        if (present(arrays)) then
            call c_str_array(arrays, storage, cptrs, arr, count)
        else
            arr = c_null_ptr
            count = 0_c_int64_t
        end if
        out%handle = c_mio_conservative_interpolate(source%handle, target%handle, arr, &
                                                     count, cdefault, c_str(cconflict))
        if (.not. c_associated(out%handle)) then
            call handle_failure('conservative_interpolate', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Green-element undo: restore `fine`'s transitional (closure-only) cells
    !> back to their original parent, read verbatim from `coarse` -- a lookup,
    !> not a reconstruction, since refine() never renumbers or prunes points.
    !> `fine` must carry refine:cell_id/refine:parent_id/refine:level (i.e. it
    !> must come from a refine() call with record_hierarchy=.true.,
    !> record_levels=.true.); `coarse` must be the mesh that call was run on.
    !> The optional num_groups_undone/num_cells_removed report how many green
    !> sibling groups were substituted and how many cells that removed.
    function mio_undo_green(coarse, fine, num_groups_undone, num_cells_removed, &
                            stat, errmsg) result(out)
        type(mio_mesh), intent(in) :: coarse, fine
        integer(int64), intent(out), optional :: num_groups_undone, num_cells_removed
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int64_t) :: ngroups, nremoved
        out%handle = c_mio_undo_green(coarse%handle, fine%handle, ngroups, nremoved)
        if (.not. c_associated(out%handle)) then
            call handle_failure('undo_green', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(num_groups_undone)) num_groups_undone = int(ngroups, int64)
        if (present(num_cells_removed)) num_cells_removed = int(nremoved, int64)
        call clear_status(stat, errmsg)
    end function

    !> Renumber the mesh (method: "rcm", "morton", or "hilbert") as a pure
    !> permutation. Returns the renumbered mesh; the optional `node_perm`
    !> receives the node permutation (1-based old index -> 1-based new index).
    function mesh_reorder(self, method, node_perm, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: method
        integer(int64), allocatable, intent(out), optional :: node_perm(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res, cdata
        integer(c_int) :: s, dt
        integer(c_int64_t) :: nlen
        integer(c_int64_t), pointer :: fp(:)
        res = c_mio_reorder(self%handle, c_str(method))
        if (.not. c_associated(res)) then
            call handle_failure('reorder', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(node_perm)) then
            s = c_mio_reorder_result_node_perm(res, cdata, dt, nlen)
            if (s /= 0_c_int) then
                call c_mio_reorder_result_free(res)
                call handle_failure('reorder', mio_error_message(), stat, errmsg)
                return
            end if
            allocate (node_perm(nlen))
            if (nlen > 0) then
                call c_f_pointer(cdata, fp, [nlen])
                node_perm = int(fp, int64) + 1_int64  ! 0-based new id -> 1-based
            end if
        end if
        out%handle = c_mio_reorder_result_take_mesh(res)
        call c_mio_reorder_result_free(res)
        if (.not. c_associated(out%handle)) then
            call handle_failure('reorder', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Apply an affine transform to the mesh's point coordinates. `matrix` is a
    !> row-major 4x4 affine matrix flattened to 16 doubles (point p maps to
    !> M * [p, 1]). `rotate_vector_data` (default .false.) rotates vector/tensor
    !> point_data by the transform's linear part. Returns the transformed mesh.
    function mesh_transform(self, matrix, rotate_vector_data, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        real(real64), intent(in) :: matrix(16)
        logical, intent(in), optional :: rotate_vector_data
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        real(c_double) :: cmat(16)
        integer(c_int) :: crot
        cmat = real(matrix, c_double)
        crot = 0
        if (present(rotate_vector_data)) then
            if (rotate_vector_data) crot = 1
        end if
        out%handle = c_mio_transform(self%handle, cmat, crot)
        if (.not. c_associated(out%handle)) then
            call handle_failure('transform', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Clean the mesh in one pass (weld / prune / de-dup). Each step is optional;
    !> defaults: no weld, remove orphans, drop degenerate, drop duplicate cells.
    !> The optional integer out-args receive the removal counts.
    function mesh_clean(self, weld, atol, remove_orphans, drop_degenerate, &
                        drop_duplicate_cells, points_welded, points_removed_orphan, &
                        cells_dropped_degenerate, cells_dropped_duplicate, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        logical, intent(in), optional :: weld, remove_orphans, drop_degenerate, &
                                         drop_duplicate_cells
        real(real64), intent(in), optional :: atol
        integer(int64), intent(out), optional :: points_welded, points_removed_orphan, &
                                                 cells_dropped_degenerate, cells_dropped_duplicate
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: cweld, crm, cdeg, cdup
        real(c_double) :: catol
        integer(c_int64_t) :: nweld, norph, ndeg, ndup
        cweld = 0
        if (present(weld)) then
            if (weld) cweld = 1
        end if
        crm = 1
        if (present(remove_orphans)) then
            if (.not. remove_orphans) crm = 0
        end if
        cdeg = 1
        if (present(drop_degenerate)) then
            if (.not. drop_degenerate) cdeg = 0
        end if
        cdup = 1
        if (present(drop_duplicate_cells)) then
            if (.not. drop_duplicate_cells) cdup = 0
        end if
        catol = 1.0e-8_c_double
        if (present(atol)) catol = real(atol, c_double)
        out%handle = c_mio_clean(self%handle, cweld, catol, crm, cdeg, cdup, &
                                 nweld, norph, ndeg, ndup)
        if (.not. c_associated(out%handle)) then
            call handle_failure('clean', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(points_welded)) points_welded = int(nweld, int64)
        if (present(points_removed_orphan)) points_removed_orphan = int(norph, int64)
        if (present(cells_dropped_degenerate)) cells_dropped_degenerate = int(ndeg, int64)
        if (present(cells_dropped_duplicate)) cells_dropped_duplicate = int(ndup, int64)
        call clear_status(stat, errmsg)
    end function

    !> One-time ABI layout guard for mio_smooth_opts_t, the same
    !> check_refine_opts_layout precedent applied here. Runs once per
    !> process (a SAVE'd flag), not on every smooth() call.
    subroutine check_smooth_opts_layout()
        logical, save :: checked = .false.
        type(mio_smooth_opts_t) :: probe
        if (checked) return
        checked = .true.
        if (c_sizeof(probe) /= 104_c_size_t) then
            write (error_unit, '(a,i0,a)') &
                'meshio++: mio_smooth_opts_t layout mismatch (', c_sizeof(probe), ' bytes)'
            error stop 1
        end if
    end subroutine

    !> Smooth the mesh's point coordinates, leaving topology and data intact.
    !> `method` is "laplacian" or "taubin"; `iterations` is the pass count (for
    !> taubin one iteration is two passes). A negative `lambda` — the default —
    !> means "this method's own default" (0.5 laplacian, 0.33 taubin). The
    !> optional out-args receive the run summary.
    function mesh_smooth(self, method, iterations, lambda, mu, fix_boundary, preserve_features, &
                         feature_angle, guard_inversion, nodes_moved, max_displacement, &
                         skipped_inversion, stat, errmsg, frozen) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: method
        integer, intent(in) :: iterations
        real(real64), intent(in), optional :: lambda, mu, feature_angle
        logical, intent(in), optional :: fix_boundary, preserve_features, guard_inversion
        integer(int64), intent(out), optional :: nodes_moved, skipped_inversion
        real(real64), intent(out), optional :: max_displacement
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        !> 1-based point ids to pin outright, unioned with the boundary/
        !> feature pins. Shifted to the C API's 0-based numbering inside.
        integer(int64), intent(in), optional :: frozen(:)
        type(mio_mesh) :: out
        type(c_ptr) :: res
        integer(c_int64_t) :: nmoved, nskip
        real(c_double) :: maxdisp
        type(mio_smooth_opts_t) :: opts
        ! NUL-terminated copy and the frozen-id buffer must outlive the call;
        ! c_loc needs them contiguous and TARGET, the mesh_refine idiom.
        character(kind=c_char, len=STRBUF_LEN), target :: method_buf
        integer(c_int64_t), allocatable, target :: frozen_ids(:)

        call check_smooth_opts_layout()
        call c_mio_smooth_opts_init(opts)
        method_buf = trim(method)//c_null_char
        opts%method = c_loc(method_buf(1:1))
        opts%iterations = int(iterations, c_int32_t)
        if (present(lambda)) opts%lambda = real(lambda, c_double)
        if (present(mu)) opts%mu = real(mu, c_double)
        if (present(feature_angle)) opts%feature_angle = real(feature_angle, c_double)
        if (present(fix_boundary)) then
            opts%fix_boundary = 0
            if (fix_boundary) opts%fix_boundary = 1
        end if
        if (present(preserve_features)) then
            opts%preserve_features = 0
            if (preserve_features) opts%preserve_features = 1
        end if
        if (present(guard_inversion)) then
            opts%guard_inversion = 0
            if (guard_inversion) opts%guard_inversion = 1
        end if
        if (present(frozen)) then
            allocate (frozen_ids(max(size(frozen), 1)))
            frozen_ids = 0_c_int64_t
            if (size(frozen) > 0) frozen_ids(1:size(frozen)) = int(frozen, c_int64_t) - 1_c_int64_t
            opts%frozen = c_loc(frozen_ids(1))
            opts%num_frozen = int(size(frozen), c_int64_t)
        end if

        res = c_mio_smooth_ex(self%handle, opts, nmoved, maxdisp, nskip)
        if (.not. c_associated(res)) then
            call handle_failure('smooth', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(nodes_moved)) nodes_moved = int(nmoved, int64)
        if (present(max_displacement)) max_displacement = real(maxdisp, real64)
        if (present(skipped_inversion)) skipped_inversion = int(nskip, int64)
        call clear_status(stat, errmsg)
    end function

    !> ODT-remesh a tetrahedral mesh: relocate vertices AND flip connectivity
    !> (2-3/3-2, predicate-free). Tet-only. The point set is invariant, so
    !> point_data and named Point regions carry; cell_data and Cell/Side
    !> regions are dropped. See doc/optimize_volume.md.
    function mesh_optimize_volume(self, max_iterations, relocate, flip, preserve_boundary, &
                                  min_improvement, num_flips, num_23_flips, num_32_flips, &
                                  num_vertices_moved, num_tets, min_quality_before, &
                                  min_quality_after, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer, intent(in), optional :: max_iterations
        logical, intent(in), optional :: relocate, flip, preserve_boundary
        real(real64), intent(in), optional :: min_improvement
        integer(int64), intent(out), optional :: num_flips, num_23_flips, num_32_flips, &
                                                 num_vertices_moved, num_tets
        real(real64), intent(out), optional :: min_quality_before, min_quality_after
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: cmaxiter, crel, cflip, cpresb
        real(c_double) :: cminimp
        integer(c_int64_t) :: nflips, n23, n32, nmoved, ntets
        real(c_double) :: qbefore, qafter
        cmaxiter = 10
        if (present(max_iterations)) cmaxiter = int(max_iterations, c_int)
        crel = 1
        if (present(relocate)) then
            if (.not. relocate) crel = 0
        end if
        cflip = 1
        if (present(flip)) then
            if (.not. flip) cflip = 0
        end if
        cpresb = 1
        if (present(preserve_boundary)) then
            if (.not. preserve_boundary) cpresb = 0
        end if
        cminimp = 1.0e-6_c_double
        if (present(min_improvement)) cminimp = real(min_improvement, c_double)
        out%handle = c_mio_optimize_volume(self%handle, cmaxiter, crel, cflip, cpresb, cminimp, &
                                           nflips, n23, n32, nmoved, ntets, qbefore, qafter)
        if (.not. c_associated(out%handle)) then
            call handle_failure('optimize_volume', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(num_flips)) num_flips = int(nflips, int64)
        if (present(num_23_flips)) num_23_flips = int(n23, int64)
        if (present(num_32_flips)) num_32_flips = int(n32, int64)
        if (present(num_vertices_moved)) num_vertices_moved = int(nmoved, int64)
        if (present(num_tets)) num_tets = int(ntets, int64)
        if (present(min_quality_before)) min_quality_before = real(qbefore, real64)
        if (present(min_quality_after)) min_quality_after = real(qafter, real64)
        call clear_status(stat, errmsg)
    end function

    !> Crop the mesh to an axis-aligned bounding box. `lo`/`hi` are the box
    !> corners (3 each). `mode` is "all" (default) or "any". Returns the crop.
    function mesh_crop_bbox(self, lo, hi, mode, record_ids, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        real(real64), intent(in) :: lo(3), hi(3)
        character(*), intent(in), optional :: mode
        logical, intent(in), optional :: record_ids
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: cmode, crec
        cmode = 0
        if (present(mode)) then
            if (mode == 'any') cmode = 1
        end if
        crec = 0
        if (present(record_ids)) then
            if (record_ids) crec = 1
        end if
        out%handle = c_mio_crop_bbox(self%handle, real(lo, c_double), real(hi, c_double), &
                                     cmode, crec)
        if (.not. c_associated(out%handle)) then
            call handle_failure('crop_bbox', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Crop the mesh to the half-space (p - point) . normal >= 0. `mode` is
    !> "all" (default) or "any". Returns the crop.
    function mesh_crop_plane(self, point, normal, mode, record_ids, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        real(real64), intent(in) :: point(3), normal(3)
        character(*), intent(in), optional :: mode
        logical, intent(in), optional :: record_ids
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: cmode, crec
        cmode = 0
        if (present(mode)) then
            if (mode == 'any') cmode = 1
        end if
        crec = 0
        if (present(record_ids)) then
            if (record_ids) crec = 1
        end if
        out%handle = c_mio_crop_plane(self%handle, real(point, c_double), &
                                      real(normal, c_double), cmode, crec)
        if (.not. c_associated(out%handle)) then
            call handle_failure('crop_plane', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Planar cross-section of the mesh: the intersection with the plane through
    !> `origin` with `normal`, one topological dimension below the cut cells (a
    !> volume mesh -> a triangle/quad surface, a 2D surface -> a line mesh). The
    !> input is simplexified first (marching tetrahedra); shared cut points are
    !> deduped so the section is watertight. `record_parent_ids` (default
    !> .false.) attaches a slice:parent_cell cell_data array.
    function mesh_slice(self, origin, normal, record_parent_ids, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        real(real64), intent(in) :: origin(3), normal(3)
        logical, intent(in), optional :: record_parent_ids
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: crec
        crec = 0
        if (present(record_parent_ids)) then
            if (record_parent_ids) crec = 1
        end if
        out%handle = c_mio_slice(self%handle, real(origin, c_double), &
                                 real(normal, c_double), crec)
        if (.not. c_associated(out%handle)) then
            call handle_failure('slice', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Isosurfaces / contours: the level set of the scalar point_data array
    !> `array` at each of `isovalues`, one topological dimension below the cut
    !> cells (a volume mesh -> a triangle/quad surface, a 2D surface -> a line
    !> mesh). The data-driven sibling of `slice`, sharing its marching-tetrahedra
    !> cutter, so contours are watertight; faces are wound toward increasing
    !> field. `array` must be point_data -- cell_data is piecewise constant and
    !> has no level set (convert it with data_cell_to_point first). All contours
    !> land in the one result, tagged per cell with a Float64 iso:value and an
    !> Int64 iso:index (the ordinal, which is the integer tag `split` needs).
    !> `component` (default -1) picks a component of a multi-component array,
    !> negative meaning the row magnitude. `record_parent_ids` (default .false.)
    !> attaches an iso:parent_cell cell_data array.
    function mesh_isosurface(self, array, isovalues, component, record_parent_ids, &
                             stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: array
        real(real64), intent(in) :: isovalues(:)
        integer, intent(in), optional :: component
        logical, intent(in), optional :: record_parent_ids
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: ccomp, crec
        ccomp = -1_c_int
        if (present(component)) ccomp = int(component, c_int)
        crec = 0
        if (present(record_parent_ids)) then
            if (record_parent_ids) crec = 1
        end if
        out%handle = c_mio_isosurface(self%handle, c_str(array), &
                                      real(isovalues, c_double), &
                                      int(size(isovalues), c_int), ccomp, crec)
        if (.not. c_associated(out%handle)) then
            call handle_failure('isosurface', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Differentiate a point_data field: its gradient, divergence or curl.
    !>
    !> `op` is "gradient" (default), "divergence" or "curl"; `method` is
    !> "green-gauss" (default) or "least-squares"; `location` is "cell"
    !> (default) or "point". The result is named `<array>:<op>` unless `output`
    !> overrides it. A gradient of an nc-component field has 3*nc components
    !> laid out [component][derivative]; divergence gives 1 and curl 3, both
    !> needing a 2- or 3-component field. `component` (default -1 = every
    !> component -- note this is the OPPOSITE of `isosurface`, where negative
    !> means the row magnitude) selects one component of the gradient.
    !>
    !> Cells that cannot be differentiated yield NaN and are reported in
    !> `num_skipped`; least-squares cells with a degenerate neighbourhood fall
    !> back to Green-Gauss and are reported in `num_fallback`.
    function mesh_gradient(self, array, op, method, location, output, component, &
                           overwrite, num_skipped, num_fallback, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: array
        character(*), intent(in), optional :: op, method, location, output
        integer, intent(in), optional :: component
        logical, intent(in), optional :: overwrite
        integer(int64), intent(out), optional :: num_skipped, num_fallback
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: ccomp, cover
        integer(c_int64_t) :: nskip, nfall
        character(:), allocatable :: cop, cmethod, cloc, cout
        cop = ''
        if (present(op)) cop = op
        cmethod = ''
        if (present(method)) cmethod = method
        cloc = 'cell'
        if (present(location)) cloc = location
        cout = ''
        if (present(output)) cout = output
        ccomp = -1_c_int
        if (present(component)) ccomp = int(component, c_int)
        cover = 0
        if (present(overwrite)) then
            if (overwrite) cover = 1
        end if
        nskip = 0
        nfall = 0
        out%handle = c_mio_gradient(self%handle, c_str(array), c_str(cop), c_str(cmethod), &
                                    c_str(cloc), c_str(cout), ccomp, cover, nskip, nfall)
        if (.not. c_associated(out%handle)) then
            call handle_failure('gradient', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(num_skipped)) num_skipped = int(nskip, int64)
        if (present(num_fallback)) num_fallback = int(nfall, int64)
        call clear_status(stat, errmsg)
    end function

    !> The Hessian (second derivative) of a scalar point_data field --
    !> `gradient`'s companion one order further, for curvature-based adaptive
    !> refinement.
    !>
    !> A composition of TWO `gradient` calls, not a new numerical kernel: the
    !> field is differentiated once (point location), then that (n,3)
    !> gradient is differentiated again with the default "gradient" operator,
    !> producing (n,9) -- the flattened row-major 3x3 Hessian, H(i,j) at
    !> index i*3+j. `method` (default "green-gauss") is forwarded to BOTH
    !> internal passes. `output` overrides the default "<array>:hessian"
    !> name.
    !>
    !> A field that is at most LINEAR has an exactly zero Hessian everywhere
    !> -- the one mesh-shape-independent guarantee. For a genuinely
    !> quadratic field the composition is exact on a structured/symmetric
    !> mesh away from its own boundary and a good, standard, but genuinely
    !> approximate curvature estimate on an irregular mesh (see
    !> doc/hessian.md). Input must have exactly one component -- a vector
    !> field's Hessian is a separate quantity per component.
    !>
    !> Cells that cannot be evaluated yield NaN and are reported in
    !> `num_skipped`; least-squares cells with a degenerate neighbourhood in
    !> either internal pass fall back to Green-Gauss and are reported in
    !> `num_fallback` (summed over both passes).
    function mesh_hessian(self, array, method, location, output, overwrite, num_skipped, &
                          num_fallback, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: array
        character(*), intent(in), optional :: method, location, output
        logical, intent(in), optional :: overwrite
        integer(int64), intent(out), optional :: num_skipped, num_fallback
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: cover
        integer(c_int64_t) :: nskip, nfall
        character(:), allocatable :: cmethod, cloc, cout
        cmethod = ''
        if (present(method)) cmethod = method
        cloc = 'cell'
        if (present(location)) cloc = location
        cout = ''
        if (present(output)) cout = output
        cover = 0
        if (present(overwrite)) then
            if (overwrite) cover = 1
        end if
        nskip = 0
        nfall = 0
        out%handle = c_mio_hessian(self%handle, c_str(array), c_str(cmethod), c_str(cloc), &
                                   c_str(cout), cover, nskip, nfall)
        if (.not. c_associated(out%handle)) then
            call handle_failure('hessian', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(num_skipped)) num_skipped = int(nskip, int64)
        if (present(num_fallback)) num_fallback = int(nfall, int64)
        call clear_status(stat, errmsg)
    end function

    !> Estimate the per-cell recovered-gradient (Zienkiewicz-Zhu) error of a
    !> point_data field, and optionally mark cells for refinement.
    !>
    !> A composition of `gradient` (Green-Gauss, cell location) with the
    !> measure-weighted point<->cell averaging round trip: the indicator is
    !> `sqrt(|measure| * sum((recovered - raw)^2))` per cell, attached as
    !> `output` (default "error:zz"). `marking` is "none" (default), "absolute",
    !> "fraction", or "dorfler"; when not "none" a second Int64 0/1 array
    !> `marked` (default "error:marked") is attached too, so `refine`'s own
    !> `where` selector needs no change at all -- the intended use is
    !> `refine(where="error:marked > 0.5")`. `marking_value`'s meaning depends
    !> on `marking`: an absolute indicator threshold, a fraction in (0, 1] of
    !> cells, or the Doerfler bulk fraction theta in (0, 1].
    !>
    !> Cells that cannot be evaluated read NaN in the indicator array and 0
    !> (never NaN) in the marking array, and are reported in `num_skipped`
    !> (excluded from `global_error` and from `num_marked`).
    function mesh_estimate_error(self, array, method, marking, marking_value, output, &
                                 marked, overwrite, global_error, num_skipped, num_marked, &
                                 stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: array
        character(*), intent(in), optional :: method, marking, output, marked
        real(real64), intent(in), optional :: marking_value
        logical, intent(in), optional :: overwrite
        real(real64), intent(out), optional :: global_error
        integer(int64), intent(out), optional :: num_skipped, num_marked
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: cover
        real(c_double) :: cvalue, cerror
        integer(c_int64_t) :: nskip, nmark
        character(:), allocatable :: cmethod, cmarking, cout, cmarked
        cmethod = ''
        if (present(method)) cmethod = method
        cmarking = ''
        if (present(marking)) cmarking = marking
        cvalue = 0.0_c_double
        if (present(marking_value)) cvalue = real(marking_value, c_double)
        cout = ''
        if (present(output)) cout = output
        cmarked = ''
        if (present(marked)) cmarked = marked
        cover = 0
        if (present(overwrite)) then
            if (overwrite) cover = 1
        end if
        cerror = 0.0_c_double
        nskip = 0
        nmark = 0
        out%handle = c_mio_estimate_error(self%handle, c_str(array), c_str(cmethod), &
                                          c_str(cmarking), cvalue, c_str(cout), c_str(cmarked), &
                                          cover, cerror, nskip, nmark)
        if (.not. c_associated(out%handle)) then
            call handle_failure('estimate_error', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(global_error)) global_error = real(cerror, real64)
        if (present(num_skipped)) num_skipped = int(nskip, int64)
        if (present(num_marked)) num_marked = int(nmark, int64)
        call clear_status(stat, errmsg)
    end function

    !> Remesh a surface by approximated centroidal Voronoi diagram (ACVD)
    !> clustering: replace its triangulation with a new, near-uniformly-sized,
    !> well-shaped one at `num_clusters` vertices. Unlike every other
    !> resolution-changing procedure on this type, the output has NEW points
    !> and NEW connectivity with no correspondence to `self` -- point/cell
    !> data and named regions are dropped, field data is carried.
    !>
    !> `metric` is "isotropic" (default; area-weighted centroidal distance),
    !> "quadric" (Garland-Heckbert quadric error, preserves sharp edges and
    !> corners), or "anisotropic" (clusters shaped by a local curvature
    !> tensor -- see `max_anisotropy`). `subdivide` defaults to automatic
    !> (the smallest count reaching `subsample_ratio` items per cluster,
    !> capped at `max_subdivide`); pass 0 to disable subdivision.
    !>
    !> Goes through mio_remesh_ex/mio_remesh_opts rather than the flat
    !> mio_remesh -- the mio_refine_ex precedent, needed because mio_remesh
    !> is a flat C function with no room to grow (it already changed once).
    !> Per-vertex mean and Gaussian curvature of this surface, by the angle
    !> defect (K) and the cotangent Laplace-Beltrami operator (H) -- the
    !> signed distance's companion as a node feature.
    !>
    !> Writes `curvature:mean` and `curvature:gaussian` as point data,
    !> optionally `curvature:area` and the (n, 2) `curvature:principal`.
    !> Geometry, connectivity and existing data are carried through unchanged.
    !>
    !> `total_angle_defect` is the oracle: on a CLOSED surface it is
    !> `2*pi*chi` exactly -- `4*pi` for anything sphere-like -- whatever the
    !> tessellation and whichever `dual_area`, so a value that is not that
    !> means the input is not closed or the result is not sane.
    !>
    !> `H` is orientation-dependent and `K` is not, so check
    !> `inconsistent_pairs` before trusting a sign: a nonzero count means
    !> facets disagree about which side is out. This never repairs its input.
    function mesh_curvature(self, mean, gaussian, dual_area, include_boundary, &
                            record_area, record_principal, region, num_boundary, &
                            num_isolated, num_degenerate, total_angle_defect, &
                            boundary_edges, non_manifold_edges, inconsistent_pairs, &
                            degenerate_triangles, watertight, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        logical, intent(in), optional :: mean, gaussian, include_boundary
        logical, intent(in), optional :: record_area, record_principal
        character(*), intent(in), optional :: dual_area, region
        integer(int64), intent(out), optional :: num_boundary, num_isolated, num_degenerate
        real(real64), intent(out), optional :: total_angle_defect
        integer(int64), intent(out), optional :: boundary_edges, non_manifold_edges
        integer(int64), intent(out), optional :: inconsistent_pairs, degenerate_triangles
        logical, intent(out), optional :: watertight
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_curvature_opts_t) :: opts
        type(mio_curvature_report_t) :: report
        ! NUL-terminated copy must outlive the call, so it is held here rather
        ! than built inline; c_loc needs it contiguous and TARGET -- the
        ! region_buf idiom mesh_refine and mesh_remesh already use.
        character(kind=c_char, len=STRBUF_LEN), target :: region_buf

        call c_mio_curvature_opts_init(opts)
        if (present(mean)) then
            if (.not. mean) opts%mean = 0
        end if
        if (present(gaussian)) then
            if (.not. gaussian) opts%gaussian = 0
        end if
        if (present(dual_area)) then
            if (trim(dual_area) == 'barycentric') then
                opts%dual_area = 1
            else if (trim(dual_area) == 'mixed-voronoi') then
                opts%dual_area = 0
            else
                call handle_failure('curvature', &
                    "meshio++: curvature: unknown dual area '"//trim(dual_area)// &
                    "' (expected 'mixed-voronoi' or 'barycentric')", stat, errmsg)
                return
            end if
        end if
        if (present(include_boundary)) then
            if (include_boundary) opts%include_boundary = 1
        end if
        if (present(record_area)) then
            if (record_area) opts%record_area = 1
        end if
        if (present(record_principal)) then
            if (record_principal) opts%record_principal = 1
        end if
        if (present(region)) then
            region_buf = trim(region)//c_null_char
            opts%region = c_loc(region_buf(1:1))
        end if

        res = c_mio_compute_curvature(self%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('curvature', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_boundary)) num_boundary = int(report%num_boundary, int64)
        if (present(num_isolated)) num_isolated = int(report%num_isolated, int64)
        if (present(num_degenerate)) num_degenerate = int(report%num_degenerate, int64)
        if (present(total_angle_defect)) &
            total_angle_defect = real(report%total_angle_defect, real64)
        if (present(boundary_edges)) &
            boundary_edges = int(report%quality%boundary_edges, int64)
        if (present(non_manifold_edges)) &
            non_manifold_edges = int(report%quality%non_manifold_edges, int64)
        if (present(inconsistent_pairs)) &
            inconsistent_pairs = int(report%quality%inconsistent_pairs, int64)
        if (present(degenerate_triangles)) &
            degenerate_triangles = int(report%quality%degenerate_triangles, int64)
        if (present(watertight)) watertight = (report%quality%watertight /= 0)
        call clear_status(stat, errmsg)
    end function

    !> The sharp, open, non-manifold and inconsistently wound edges of this
    !> surface (or of a volume mesh's skin), as a mesh of `line` cells over the
    !> input's points, with cell data `feature:kind` (1 feature, 2 boundary,
    !> 3 non-manifold, 4 inconsistent) and `feature:angle` (degrees). Each
    !> category is reported unless its flag is `.false.`; `feature_angle`
    !> defaults to 30. See doc/feature_edges.md.
    function mesh_feature_edges(self, feature_angle, feature, boundary, non_manifold, &
                                inconsistent, region, num_feature, num_boundary, &
                                num_non_manifold, num_inconsistent, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        real(real64), intent(in), optional :: feature_angle
        logical, intent(in), optional :: feature, boundary, non_manifold, inconsistent
        character(*), intent(in), optional :: region
        integer(int64), intent(out), optional :: num_feature, num_boundary
        integer(int64), intent(out), optional :: num_non_manifold, num_inconsistent
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_feature_edges_opts_t) :: opts
        type(mio_feature_edges_report_t) :: report
        character(kind=c_char, len=STRBUF_LEN), target :: region_buf

        call c_mio_feature_edges_opts_init(opts)
        if (present(feature_angle)) opts%feature_angle = real(feature_angle, c_double)
        if (present(feature)) opts%feature = merge(1_c_int32_t, 0_c_int32_t, feature)
        if (present(boundary)) opts%boundary = merge(1_c_int32_t, 0_c_int32_t, boundary)
        if (present(non_manifold)) &
            opts%non_manifold = merge(1_c_int32_t, 0_c_int32_t, non_manifold)
        if (present(inconsistent)) &
            opts%inconsistent = merge(1_c_int32_t, 0_c_int32_t, inconsistent)
        if (present(region)) then
            region_buf = trim(region)//c_null_char
            opts%region = c_loc(region_buf(1:1))
        end if
        res = c_mio_feature_edges(self%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('feature_edges', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_feature)) num_feature = int(report%num_feature, int64)
        if (present(num_boundary)) num_boundary = int(report%num_boundary, int64)
        if (present(num_non_manifold)) num_non_manifold = int(report%num_non_manifold, int64)
        if (present(num_inconsistent)) num_inconsistent = int(report%num_inconsistent, int64)
        call clear_status(stat, errmsg)
    end function

    !> Split this surface into smooth patches by region growing, bounded by its
    !> feature edges: the input plus cell data `segment:patch`, field data
    !> `segment:adjacency` and one cell region per patch (`<prefix>:<k>`).
    !> `feature_angle` defaults to 30; the tolerances default to off.
    !> `method = 'primitive'` (v17.9.0; default 'angle') ends a patch where the
    !> surface its facets lie on changes, at `tolerance` (default 1e-3 of the
    !> diagonal), so a fillet comes apart from its tangent planes. See
    !> doc/segment.md.
    function mesh_segment(self, feature_angle, normal_tolerance, curvature_tolerance, &
                          min_patch_cells, region, prefix, num_patches, num_merged, &
                          num_feature_edges, num_boundary_edges, stat, errmsg, &
                          method, tolerance) result(out)
        class(mio_mesh), intent(in) :: self
        real(real64), intent(in), optional :: feature_angle, normal_tolerance
        real(real64), intent(in), optional :: curvature_tolerance
        integer(int64), intent(in), optional :: min_patch_cells
        character(*), intent(in), optional :: region, prefix
        integer(int64), intent(out), optional :: num_patches, num_merged
        integer(int64), intent(out), optional :: num_feature_edges, num_boundary_edges
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(*), intent(in), optional :: method
        real(real64), intent(in), optional :: tolerance
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_segment_opts_t) :: opts
        type(mio_segment_report_t) :: report
        character(kind=c_char, len=STRBUF_LEN), target :: region_buf, prefix_buf

        call c_mio_segment_opts_init(opts)
        if (present(feature_angle)) opts%feature_angle = real(feature_angle, c_double)
        if (present(normal_tolerance)) &
            opts%normal_tolerance = real(normal_tolerance, c_double)
        if (present(curvature_tolerance)) &
            opts%curvature_tolerance = real(curvature_tolerance, c_double)
        if (present(min_patch_cells)) opts%min_patch_cells = int(min_patch_cells, c_int64_t)
        if (present(region)) then
            region_buf = trim(region)//c_null_char
            opts%region = c_loc(region_buf(1:1))
        end if
        if (present(prefix)) then
            prefix_buf = trim(prefix)//c_null_char
            opts%prefix = c_loc(prefix_buf(1:1))
        end if
        if (present(method)) then
            select case (trim(method))
            case ('angle')
                opts%method = 0
            case ('primitive')
                opts%method = 1
            case default
                call handle_failure('segment', "method must be 'angle' or 'primitive'", &
                                    stat, errmsg)
                return
            end select
        end if
        if (present(tolerance)) opts%tolerance = real(tolerance, c_double)
        res = c_mio_segment(self%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('segment', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_patches)) num_patches = int(report%num_patches, int64)
        if (present(num_merged)) num_merged = int(report%num_merged, int64)
        if (present(num_feature_edges)) num_feature_edges = int(report%num_feature_edges, int64)
        if (present(num_boundary_edges)) &
            num_boundary_edges = int(report%num_boundary_edges, int64)
        call clear_status(stat, errmsg)
    end function

    !> Fit a plane, sphere, cylinder, cone or torus to each patch of this
    !> surface: the input plus cell data `fit:type`, `fit:rms`, `fit:max` and
    !> per-patch field data `fit:patch_*`. `patch_data` names the label array
    !> (default `segment:patch`; '' fits the whole surface); `types` is the
    !> bit mask `2**type` (default 62, every primitive). See doc/fit_primitives.md.
    function mesh_fit_primitives(self, patch_data, tolerance, relative, &
                                 min_inlier_fraction, types, max_iterations, num_patches, &
                                 num_plane, num_sphere, num_cylinder, num_cone, num_torus, &
                                 num_unfitted, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in), optional :: patch_data
        real(real64), intent(in), optional :: tolerance, min_inlier_fraction
        logical, intent(in), optional :: relative
        integer, intent(in), optional :: types
        integer(int64), intent(in), optional :: max_iterations
        integer(int64), intent(out), optional :: num_patches, num_plane, num_sphere
        integer(int64), intent(out), optional :: num_cylinder, num_cone, num_torus, num_unfitted
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_fit_primitives_opts_t) :: opts
        type(mio_fit_primitives_report_t) :: report
        character(kind=c_char, len=STRBUF_LEN), target :: patch_buf

        call c_mio_fit_primitives_opts_init(opts)
        if (present(patch_data)) then
            patch_buf = trim(patch_data)//c_null_char
            opts%patch_data = c_loc(patch_buf(1:1))
        end if
        if (present(tolerance)) opts%tolerance = real(tolerance, c_double)
        if (present(min_inlier_fraction)) &
            opts%min_inlier_fraction = real(min_inlier_fraction, c_double)
        if (present(relative)) opts%relative = merge(1_c_int32_t, 0_c_int32_t, relative)
        if (present(types)) opts%types = int(types, c_int32_t)
        if (present(max_iterations)) opts%max_iterations = int(max_iterations, c_int64_t)
        res = c_mio_fit_primitives(self%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('fit_primitives', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_patches)) num_patches = int(report%num_patches, int64)
        if (present(num_plane)) num_plane = int(report%num_plane, int64)
        if (present(num_sphere)) num_sphere = int(report%num_sphere, int64)
        if (present(num_cylinder)) num_cylinder = int(report%num_cylinder, int64)
        if (present(num_cone)) num_cone = int(report%num_cone, int64)
        if (present(num_torus)) num_torus = int(report%num_torus, int64)
        if (present(num_unfitted)) num_unfitted = int(report%num_unfitted, int64)
        call clear_status(stat, errmsg)
    end function

    !> Build a CAD solid from the fitted patches of this closed skin: with
    !> OpenCASCADE (a library built with -DMESHIOPLUSPLUS_WITH_OCCT=ON;
    !> mio_has_occt() tells), accepted fits become plane, sphere, cylinder,
    !> cone and torus faces, the rest planar facets, and a face that fails
    !> OpenCASCADE's check or strays beyond `tolerance` is demoted to facets;
    !> without it, the native backend below runs instead (since v17.9.0). Writes AP214 STEP to
    !> `step_path` when given; returns the solid re-meshed as triangles with
    !> cell data `fit_cad:face`, `fit_cad:type` and `fit_cad:patch`. A patch no
    !> primitive explains becomes, with `fillets` (default on), a rolling-ball
    !> fillet face between two analytic faces (of one radius, or a radius
    !> varying along it: `num_variable_fillet`), with `threads` (default on), a
    !> helical thread face about a neighbouring cylinder's or cone's axis or
    !> the screw it fits (a profile not straight: `num_rounded_thread`), and
    !> otherwise, with `freeform` (default on), a B-spline face of degree up to
    !> `freeform_degree` on at most `freeform_max_control` grid points per
    !> direction; with `freeform_split` (default on), a freeform patch that
    !> folds over both flattenings is split into pieces. With `fillets`, a
    !> patch between two fitted faces, one of them freeform, then becomes a
    !> rolling-ball blend (`num_blend`). Since v17.2.0 a patch no rolling ball
    !> or thread explains between two analytic faces becomes a chamfer
    !> (`num_chamfer`) or a conic blend (`num_conic_blend`), a thread flank of
    !> varying pitch or cut on a taper is recognised (`num_variable_pitch_thread`,
    !> `num_tapered_thread`), and a corner where three or more blends meet
    !> becomes a vertex blend (`num_vertex_blend`). Since v17.9.0, `backend`
    !> ('auto' by default, 'occt' or 'native') picks the engine: 'native' needs
    !> no CAD kernel and writes the fitted patches as analytic STEP faces
    !> (`used_native` says which one ran), and `merge_faces` (default on)
    !> merges faces on one surface. Since v17.11.0, `refine_segmentation`
    !> (default off) first refines the patches as `refine_patches` does, with
    !> `refine_tolerance` (default 1e-3, of the diagonal when
    !> `relative_tolerance`, the default), `normal_weight`, `refine_epsilon`,
    !> `max_refine_iterations`, `allow_split`, `allow_merge` and `ransac_seed`;
    !> `refined`, `num_refine_splits`, `num_refine_merges` and
    !> `explained_area_fraction` report it, and `refinement_kept` whether the
    !> solid was built from the refined patches (OpenCASCADE also builds it from
    !> the patches as given and keeps the better). Also since v17.11.0,
    !> `freeform_general` (default off) tries a patch the plane and cylinder
    !> freeform fit leaves unexplained as a general B-spline face
    !> (`freeform_tolerance`, `freeform_smoothing`, `max_freeform_poles`,
    !> `max_freeform_iterations`, `max_freeform_splits`; reported by
    !> `num_freeform_general`, `num_freeform_refinements`, `avg_freeform_poles`),
    !> and, OpenCASCADE only, `sew_tolerance` (0: the tolerance) sets the sewing
    !> tolerance, `extend_surfaces` extends freeform surfaces past their patches
    !> (`num_extended_edges`), `snap_junctions` places the corners where
    !> analytic faces meet by least squares (`num_snapped_junctions`) and
    !> `face_deviation` measures each face against the skin (the
    !> `fit_cad:face_hausdorff` and `fit_cad:face_rms` field data and
    !> `hausdorff_max` / `hausdorff_rms`, NaN when not measured). `closed` says
    !> whether every shell is closed. See doc/fit_cad.md.
    function mesh_fit_cad(self, step_path, tolerance, exact_edges, types, &
                          mesh_deflection, max_demotion_rounds, valid, volume_ratio, &
                          max_deviation, num_faces, num_faceted_faces, num_demoted, &
                          stat, errmsg, freeform, freeform_degree, freeform_max_control, &
                          num_bspline, freeform_split, fillets, threads, &
                          num_split_patches, num_fillet, num_helical, &
                          num_variable_fillet, num_blend, num_rounded_thread, &
                          num_vertex_blend, num_chamfer, num_conic_blend, &
                          num_tapered_thread, num_variable_pitch_thread, &
                          merge_faces, backend, used_native, refine_segmentation, &
                          max_refine_iterations, refine_epsilon, normal_weight, &
                          allow_merge, allow_split, ransac_seed, refine_tolerance, &
                          relative_tolerance, refined, num_refine_splits, &
                          num_refine_merges, explained_area_fraction, freeform_general, &
                          freeform_tolerance, freeform_smoothing, max_freeform_poles, &
                          max_freeform_iterations, max_freeform_splits, sew_tolerance, &
                          extend_surfaces, snap_junctions, face_deviation, closed, &
                          num_freeform_general, num_freeform_refinements, &
                          avg_freeform_poles, num_extended_edges, &
                          num_snapped_junctions, hausdorff_max, hausdorff_rms, &
                          refinement_kept, freeform_mode, freeform_parameterization, &
                          extraordinary_vertex, local_refine, freeform_g1_tolerance_deg, &
                          freeform_max_control_points, layout_target_patch_size, &
                          layout_max_valence_deviation, resegment_rounds, &
                          layout_snap_to_features, validate_by_surface_sampling, &
                          brep_checked, facet_area_fraction, num_freeform_adaptive, &
                          num_network_patches, num_network_edges, &
                          num_extraordinary_vertices, num_resegmented, num_escalations, &
                          max_g1_angle_deg, max_freeform_distance, rms_freeform_distance, &
                          num_face_records) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in), optional :: step_path
        real(real64), intent(in), optional :: tolerance, mesh_deflection
        logical, intent(in), optional :: exact_edges
        integer, intent(in), optional :: types
        integer(int64), intent(in), optional :: max_demotion_rounds
        logical, intent(out), optional :: valid
        real(real64), intent(out), optional :: volume_ratio, max_deviation
        integer(int64), intent(out), optional :: num_faces, num_faceted_faces, num_demoted
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        logical, intent(in), optional :: freeform
        integer, intent(in), optional :: freeform_degree
        integer(int64), intent(in), optional :: freeform_max_control
        integer(int64), intent(out), optional :: num_bspline
        logical, intent(in), optional :: freeform_split, fillets, threads
        integer(int64), intent(out), optional :: num_split_patches, num_fillet, num_helical
        integer(int64), intent(out), optional :: num_variable_fillet, num_blend, num_rounded_thread
        integer(int64), intent(out), optional :: num_vertex_blend, num_chamfer, num_conic_blend
        integer(int64), intent(out), optional :: num_tapered_thread, num_variable_pitch_thread
        logical, intent(in), optional :: merge_faces
        character(*), intent(in), optional :: backend
        logical, intent(out), optional :: used_native
        logical, intent(in), optional :: refine_segmentation
        integer(int64), intent(in), optional :: max_refine_iterations
        real(real64), intent(in), optional :: refine_epsilon, normal_weight
        logical, intent(in), optional :: allow_merge, allow_split
        integer(int64), intent(in), optional :: ransac_seed
        real(real64), intent(in), optional :: refine_tolerance
        logical, intent(in), optional :: relative_tolerance
        logical, intent(out), optional :: refined
        integer(int64), intent(out), optional :: num_refine_splits, num_refine_merges
        real(real64), intent(out), optional :: explained_area_fraction
        logical, intent(in), optional :: freeform_general
        real(real64), intent(in), optional :: freeform_tolerance, freeform_smoothing
        integer(int64), intent(in), optional :: max_freeform_poles, max_freeform_iterations
        integer(int64), intent(in), optional :: max_freeform_splits
        real(real64), intent(in), optional :: sew_tolerance
        logical, intent(in), optional :: extend_surfaces, snap_junctions, face_deviation
        logical, intent(out), optional :: closed
        integer(int64), intent(out), optional :: num_freeform_general, num_freeform_refinements
        real(real64), intent(out), optional :: avg_freeform_poles
        integer(int64), intent(out), optional :: num_extended_edges, num_snapped_junctions
        real(real64), intent(out), optional :: hausdorff_max, hausdorff_rms
        logical, intent(out), optional :: refinement_kept
        character(*), intent(in), optional :: freeform_mode, freeform_parameterization
        character(*), intent(in), optional :: extraordinary_vertex, local_refine
        real(real64), intent(in), optional :: freeform_g1_tolerance_deg, layout_target_patch_size
        integer(int64), intent(in), optional :: freeform_max_control_points
        integer(int64), intent(in), optional :: layout_max_valence_deviation, resegment_rounds
        logical, intent(in), optional :: layout_snap_to_features, validate_by_surface_sampling
        logical, intent(out), optional :: brep_checked
        real(real64), intent(out), optional :: facet_area_fraction, max_g1_angle_deg
        real(real64), intent(out), optional :: max_freeform_distance, rms_freeform_distance
        integer(int64), intent(out), optional :: num_freeform_adaptive, num_network_patches
        integer(int64), intent(out), optional :: num_network_edges, num_extraordinary_vertices
        integer(int64), intent(out), optional :: num_resegmented, num_escalations
        integer(int64), intent(out), optional :: num_face_records
        type(mio_mesh) :: out
        type(c_ptr) :: res, path_ptr
        type(mio_fit_cad_opts_t) :: opts
        type(mio_fit_cad_opts_ex_t) :: opts_ex
        type(mio_fit_cad_opts_ex2_t) :: opts_ex2
        type(mio_fit_cad_report_ex3_t) :: report_ex3
        integer :: code
        type(mio_fit_cad_report_t) :: report
        type(mio_fit_cad_report_ex_t) :: report_ex
        type(mio_fit_cad_report_ex2_t) :: report_ex2
        character(kind=c_char, len=STRBUF_LEN), target :: path_buf

        call check_refine_patches_layout()
        call c_mio_fit_cad_opts_init(opts)
        call c_mio_fit_cad_opts_ex_init(opts_ex)
        call c_mio_fit_cad_opts_ex2_init(opts_ex2)
        if (present(freeform_mode)) then
            code = keyword_index(freeform_mode, [character(8) :: 'legacy', 'adaptive', &
                                                 'network', 'auto'])
            if (code < 0) then
                call handle_failure('fit_cad', &
                    "freeform_mode must be 'legacy', 'adaptive', 'network' or 'auto'", &
                    stat, errmsg)
                return
            end if
            opts_ex2%freeform_mode = int(code, c_int32_t)
        end if
        if (present(freeform_parameterization)) then
            code = keyword_index(freeform_parameterization, [character(10) :: 'auto', &
                                 'projection', 'lscm', 'harmonic'])
            if (code < 0) then
                call handle_failure('fit_cad', "freeform_parameterization must be 'auto', " // &
                    "'projection', 'lscm' or 'harmonic'", stat, errmsg)
                return
            end if
            opts_ex2%freeform_parameterization = int(code, c_int32_t)
        end if
        if (present(extraordinary_vertex)) then
            code = keyword_index(extraordinary_vertex, [character(4) :: 'g0', 'g1', 'g1x', 'g1xx'])
            if (code < 0) then
                call handle_failure('fit_cad', &
                    "extraordinary_vertex must be 'g0', 'g1', 'g1x' or 'g1xx'", stat, errmsg)
                return
            end if
            opts_ex2%extraordinary_vertex = int(code, c_int32_t)
        end if
        if (present(local_refine)) then
            code = keyword_index(local_refine, [character(5) :: 'knots', 'thb'])
            if (code < 0) then
                call handle_failure('fit_cad', "local_refine must be 'knots' or 'thb'", &
                                    stat, errmsg)
                return
            end if
            opts_ex2%local_refine = int(code, c_int32_t)
        end if
        if (present(freeform_g1_tolerance_deg)) &
            opts_ex2%freeform_g1_tolerance_deg = real(freeform_g1_tolerance_deg, c_double)
        if (present(freeform_max_control_points)) &
            opts_ex2%freeform_max_control_points = int(freeform_max_control_points, c_int64_t)
        if (present(layout_target_patch_size)) &
            opts_ex2%layout_target_patch_size = real(layout_target_patch_size, c_double)
        if (present(layout_max_valence_deviation)) &
            opts_ex2%layout_max_valence_deviation = int(layout_max_valence_deviation, c_int64_t)
        if (present(resegment_rounds)) &
            opts_ex2%resegment_rounds = int(resegment_rounds, c_int64_t)
        if (present(layout_snap_to_features)) opts_ex2%layout_snap_to_features = &
            merge(1_c_int32_t, 0_c_int32_t, layout_snap_to_features)
        if (present(validate_by_surface_sampling)) opts_ex2%validate_by_surface_sampling = &
            merge(1_c_int32_t, 0_c_int32_t, validate_by_surface_sampling)
        if (present(refine_segmentation)) &
            opts_ex%refine_segmentation = merge(1_c_int32_t, 0_c_int32_t, refine_segmentation)
        if (present(max_refine_iterations)) &
            opts_ex%refine%max_iterations = int(max_refine_iterations, c_int64_t)
        if (present(refine_epsilon)) opts_ex%refine%epsilon = real(refine_epsilon, c_double)
        if (present(normal_weight)) opts_ex%refine%normal_weight = real(normal_weight, c_double)
        if (present(allow_merge)) &
            opts_ex%refine%allow_merge = merge(1_c_int32_t, 0_c_int32_t, allow_merge)
        if (present(allow_split)) &
            opts_ex%refine%allow_split = merge(1_c_int32_t, 0_c_int32_t, allow_split)
        if (present(ransac_seed)) opts_ex%refine%ransac_seed = int(ransac_seed, c_int64_t)
        if (present(refine_tolerance)) &
            opts_ex%refine%tolerance = real(refine_tolerance, c_double)
        if (present(relative_tolerance)) &
            opts_ex%refine%relative = merge(1_c_int32_t, 0_c_int32_t, relative_tolerance)
        if (present(freeform_general)) &
            opts_ex%freeform_general = merge(1_c_int32_t, 0_c_int32_t, freeform_general)
        if (present(freeform_tolerance)) &
            opts_ex%freeform_tolerance = real(freeform_tolerance, c_double)
        if (present(freeform_smoothing)) &
            opts_ex%freeform_smoothing = real(freeform_smoothing, c_double)
        if (present(max_freeform_poles)) &
            opts_ex%max_freeform_poles = int(max_freeform_poles, c_int64_t)
        if (present(max_freeform_iterations)) &
            opts_ex%max_freeform_iterations = int(max_freeform_iterations, c_int64_t)
        if (present(max_freeform_splits)) &
            opts_ex%max_freeform_splits = int(max_freeform_splits, c_int64_t)
        if (present(sew_tolerance)) opts_ex%sew_tolerance = real(sew_tolerance, c_double)
        if (present(extend_surfaces)) &
            opts_ex%extend_surfaces = merge(1_c_int32_t, 0_c_int32_t, extend_surfaces)
        if (present(snap_junctions)) &
            opts_ex%snap_junctions = merge(1_c_int32_t, 0_c_int32_t, snap_junctions)
        if (present(face_deviation)) &
            opts_ex%face_deviation = merge(1_c_int32_t, 0_c_int32_t, face_deviation)
        path_ptr = c_null_ptr
        if (present(step_path)) then
            path_buf = trim(step_path)//c_null_char
            path_ptr = c_loc(path_buf(1:1))
        end if
        if (present(tolerance)) opts%tolerance = real(tolerance, c_double)
        if (present(mesh_deflection)) opts%mesh_deflection = real(mesh_deflection, c_double)
        if (present(exact_edges)) opts%exact_edges = merge(1_c_int32_t, 0_c_int32_t, exact_edges)
        if (present(types)) opts%types = int(types, c_int32_t)
        if (present(max_demotion_rounds)) &
            opts%max_demotion_rounds = int(max_demotion_rounds, c_int64_t)
        if (present(freeform)) opts%freeform = merge(1_c_int32_t, 0_c_int32_t, freeform)
        if (present(freeform_degree)) opts%freeform_degree = int(freeform_degree, c_int32_t)
        if (present(freeform_max_control)) &
            opts%freeform_max_control = int(freeform_max_control, c_int64_t)
        if (present(freeform_split)) &
            opts%freeform_split = merge(1_c_int32_t, 0_c_int32_t, freeform_split)
        if (present(fillets)) opts%fillets = merge(1_c_int32_t, 0_c_int32_t, fillets)
        if (present(threads)) opts%threads = merge(1_c_int32_t, 0_c_int32_t, threads)
        if (present(merge_faces)) &
            opts%merge_faces = merge(1_c_int32_t, 0_c_int32_t, merge_faces)
        if (present(backend)) then
            select case (trim(backend))
            case ('auto')
                opts%backend = 0
            case ('occt')
                opts%backend = 1
            case ('native')
                opts%backend = 2
            case default
                call handle_failure('fit_cad', "backend must be 'auto', 'occt' or 'native'", &
                                    stat, errmsg)
                return
            end select
        end if
        res = c_mio_fit_cad3(self%handle, path_ptr, opts, opts_ex, opts_ex2, report, report_ex, &
                             report_ex2, report_ex3)
        if (.not. c_associated(res)) then
            call handle_failure('fit_cad', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(valid)) valid = report%valid /= 0
        if (present(volume_ratio)) volume_ratio = real(report%volume_ratio, real64)
        if (present(max_deviation)) max_deviation = real(report%max_deviation, real64)
        if (present(num_faces)) num_faces = int(report%num_faces, int64)
        if (present(num_faceted_faces)) num_faceted_faces = int(report%num_faceted_faces, int64)
        if (present(num_demoted)) num_demoted = int(report%num_demoted, int64)
        if (present(num_bspline)) num_bspline = int(report%num_bspline, int64)
        if (present(num_split_patches)) num_split_patches = int(report%num_split_patches, int64)
        if (present(num_fillet)) num_fillet = int(report%num_fillet, int64)
        if (present(num_helical)) num_helical = int(report%num_helical, int64)
        if (present(num_variable_fillet)) &
            num_variable_fillet = int(report%num_variable_fillet, int64)
        if (present(num_blend)) num_blend = int(report%num_blend, int64)
        if (present(num_rounded_thread)) num_rounded_thread = int(report%num_rounded_thread, int64)
        if (present(num_vertex_blend)) num_vertex_blend = int(report_ex%num_vertex_blend, int64)
        if (present(num_chamfer)) num_chamfer = int(report_ex%num_chamfer, int64)
        if (present(num_conic_blend)) num_conic_blend = int(report_ex%num_conic_blend, int64)
        if (present(num_tapered_thread)) &
            num_tapered_thread = int(report_ex%num_tapered_thread, int64)
        if (present(num_variable_pitch_thread)) &
            num_variable_pitch_thread = int(report_ex%num_variable_pitch_thread, int64)
        if (present(used_native)) used_native = report_ex%backend == 2
        if (present(refined)) refined = report_ex2%refined /= 0
        if (present(num_refine_splits)) &
            num_refine_splits = int(report_ex2%num_refine_splits, int64)
        if (present(num_refine_merges)) &
            num_refine_merges = int(report_ex2%num_refine_merges, int64)
        if (present(explained_area_fraction)) &
            explained_area_fraction = real(report_ex2%explained_area_fraction, real64)
        if (present(closed)) closed = report_ex2%closed /= 0
        if (present(refinement_kept)) refinement_kept = report_ex2%refinement_kept /= 0
        if (present(num_freeform_general)) &
            num_freeform_general = int(report_ex2%num_freeform_general, int64)
        if (present(num_freeform_refinements)) &
            num_freeform_refinements = int(report_ex2%num_freeform_refinements, int64)
        if (present(avg_freeform_poles)) &
            avg_freeform_poles = real(report_ex2%avg_freeform_poles, real64)
        if (present(num_extended_edges)) &
            num_extended_edges = int(report_ex2%num_extended_edges, int64)
        if (present(num_snapped_junctions)) &
            num_snapped_junctions = int(report_ex2%num_snapped_junctions, int64)
        if (present(hausdorff_max)) hausdorff_max = real(report_ex2%hausdorff_max, real64)
        if (present(hausdorff_rms)) hausdorff_rms = real(report_ex2%hausdorff_rms, real64)
        if (present(brep_checked)) brep_checked = report_ex3%brep_checked /= 0
        if (present(facet_area_fraction)) &
            facet_area_fraction = real(report_ex3%facet_area_fraction, real64)
        if (present(num_freeform_adaptive)) &
            num_freeform_adaptive = int(report_ex3%num_freeform_adaptive, int64)
        if (present(num_network_patches)) &
            num_network_patches = int(report_ex3%num_network_patches, int64)
        if (present(num_network_edges)) num_network_edges = int(report_ex3%num_network_edges, int64)
        if (present(num_extraordinary_vertices)) &
            num_extraordinary_vertices = int(report_ex3%num_extraordinary_vertices, int64)
        if (present(num_resegmented)) num_resegmented = int(report_ex3%num_resegmented, int64)
        if (present(num_escalations)) num_escalations = int(report_ex3%num_escalations, int64)
        if (present(max_g1_angle_deg)) max_g1_angle_deg = real(report_ex3%max_g1_angle_deg, real64)
        if (present(max_freeform_distance)) &
            max_freeform_distance = real(report_ex3%max_freeform_distance, real64)
        if (present(rms_freeform_distance)) &
            rms_freeform_distance = real(report_ex3%rms_freeform_distance, real64)
        if (present(num_face_records)) num_face_records = int(report_ex3%num_face_records, int64)
        call clear_status(stat, errmsg)
    end function

    !> The 0-based position of `word` (trimmed) in `words`, or -1.
    pure function keyword_index(word, words) result(idx)
        character(*), intent(in) :: word
        character(*), intent(in) :: words(:)
        integer :: idx, i
        idx = -1
        do i = 1, size(words)
            if (trim(word) == trim(words(i))) then
                idx = i - 1
                return
            end if
        end do
    end function

    !> The layouts of the refine_patches and fit_cad extension mirrors against
    !> the C structs (144, 192, 280 and 176 bytes), checked once per process.
    subroutine check_refine_patches_layout()
        logical, save :: checked = .false.
        type(mio_refine_patches_opts_t) :: opts
        type(mio_refine_patches_report_t) :: report
        type(mio_fit_cad_opts_ex_t) :: opts_ex
        type(mio_fit_cad_report_ex2_t) :: report_ex2
        type(mio_fit_cad_opts_ex2_t) :: opts_ex2
        type(mio_fit_cad_report_ex3_t) :: report_ex3
        if (checked) return
        checked = .true.
        if (c_sizeof(opts) /= 144_c_size_t .or. c_sizeof(report) /= 192_c_size_t .or. &
            c_sizeof(opts_ex) /= 280_c_size_t .or. c_sizeof(report_ex2) /= 176_c_size_t .or. &
            c_sizeof(opts_ex2) /= 128_c_size_t .or. c_sizeof(report_ex3) /= 184_c_size_t) then
            write (error_unit, '(a,6(i0,1x),a)') &
                'meshio++: refine_patches layout mismatch (', c_sizeof(opts), c_sizeof(report), &
                c_sizeof(opts_ex), c_sizeof(report_ex2), c_sizeof(opts_ex2), &
                c_sizeof(report_ex3), 'bytes)'
            error stop 1
        end if
    end subroutine

    !> Refine the `segment:patch` labels of this surface by variational shape
    !> approximation (v17.11.0): patch boundaries move, a freeform patch is
    !> split where seeded RANSAC finds a primitive inside it, and neighbouring
    !> patches one primitive explains merge, until each patch is one plane,
    !> sphere, cylinder, cone or torus, or a freeform proxy. Returns the input
    !> with `segment:patch`, `segment:adjacency` and the `<prefix>:<k>`
    !> regions rewritten, the `fit:*` data of `fit_primitives` on the refined
    !> patches and per-patch `refine:*` field data. `tolerance` (default 1e-3)
    !> is a fraction of the diagonal when `relative` (default); `normal_weight`
    !> defaults to 1, `epsilon` to 1e-3, `feature_angle` to 30, `max_iterations`
    !> to 10, `min_patch_cells` to 3, `types` to 62 (every primitive), `prefix`
    !> to 'patch' and `ransac_seed` to the bytes of "meshio++"; `allow_split`
    !> and `allow_merge` are on. See doc/refine_patches.md.
    function mesh_refine_patches(self, tolerance, relative, normal_weight, epsilon, &
                                 feature_angle, max_iterations, min_patch_cells, &
                                 ransac_seed, allow_split, allow_merge, num_patches, &
                                 num_splits, num_merges, explained_area_fraction, stat, &
                                 errmsg, types, prefix, num_iterations, num_primitive) &
                                 result(out)
        class(mio_mesh), intent(in) :: self
        real(real64), intent(in), optional :: tolerance
        logical, intent(in), optional :: relative
        real(real64), intent(in), optional :: normal_weight, epsilon, feature_angle
        integer(int64), intent(in), optional :: max_iterations, min_patch_cells
        integer(c_int64_t), intent(in), optional :: ransac_seed
        logical, intent(in), optional :: allow_split, allow_merge
        integer(int64), intent(out), optional :: num_patches, num_splits, num_merges
        real(real64), intent(out), optional :: explained_area_fraction
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer, intent(in), optional :: types
        character(*), intent(in), optional :: prefix
        integer(int64), intent(out), optional :: num_iterations, num_primitive
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_refine_patches_opts_t) :: opts
        type(mio_refine_patches_report_t) :: report
        character(kind=c_char, len=STRBUF_LEN), target :: prefix_buf

        call check_refine_patches_layout()
        call c_mio_refine_patches_opts_init(opts)
        if (present(tolerance)) opts%tolerance = real(tolerance, c_double)
        if (present(relative)) opts%relative = merge(1_c_int32_t, 0_c_int32_t, relative)
        if (present(normal_weight)) opts%normal_weight = real(normal_weight, c_double)
        if (present(epsilon)) opts%epsilon = real(epsilon, c_double)
        if (present(feature_angle)) opts%feature_angle = real(feature_angle, c_double)
        if (present(max_iterations)) opts%max_iterations = int(max_iterations, c_int64_t)
        if (present(min_patch_cells)) opts%min_patch_cells = int(min_patch_cells, c_int64_t)
        if (present(ransac_seed)) opts%ransac_seed = ransac_seed
        if (present(allow_split)) opts%allow_split = merge(1_c_int32_t, 0_c_int32_t, allow_split)
        if (present(allow_merge)) opts%allow_merge = merge(1_c_int32_t, 0_c_int32_t, allow_merge)
        if (present(types)) opts%types = int(types, c_int32_t)
        if (present(prefix)) then
            prefix_buf = trim(prefix)//c_null_char
            opts%prefix = c_loc(prefix_buf(1:1))
        end if
        res = c_mio_refine_patches(self%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('refine_patches', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_patches)) num_patches = int(report%num_patches, int64)
        if (present(num_splits)) num_splits = int(report%num_splits, int64)
        if (present(num_merges)) num_merges = int(report%num_merges, int64)
        if (present(explained_area_fraction)) &
            explained_area_fraction = real(report%explained_area_fraction, real64)
        if (present(num_iterations)) num_iterations = int(report%num_iterations, int64)
        if (present(num_primitive)) num_primitive = int(report%num_primitive, int64)
        call clear_status(stat, errmsg)
    end function

    !> The (sampled) Hausdorff distance between this surface and `other`; a
    !> volume mesh contributes its skin. `face_samples = s > 0` also samples
    !> the centroids of the s*s sub-triangles of every triangle (vertices only
    !> is a lower bound). See doc/hausdorff.md.
    function mesh_hausdorff_distance(self, other, face_samples, region_a, region_b, &
                                     a_to_b, b_to_a, mean_a_to_b, mean_b_to_a, &
                                     rms_a_to_b, rms_b_to_a, stat, errmsg) result(d)
        class(mio_mesh), intent(in) :: self
        type(mio_mesh), intent(in) :: other
        integer(int64), intent(in), optional :: face_samples
        character(*), intent(in), optional :: region_a, region_b
        real(real64), intent(out), optional :: a_to_b, b_to_a, mean_a_to_b, mean_b_to_a
        real(real64), intent(out), optional :: rms_a_to_b, rms_b_to_a
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64) :: d
        type(mio_hausdorff_opts_t) :: opts
        type(mio_hausdorff_report_t) :: report
        character(kind=c_char, len=STRBUF_LEN), target :: buf_a, buf_b
        integer(c_int) :: s

        d = -1.0_real64
        if (present(face_samples)) opts%face_samples = int(face_samples, c_int64_t)
        if (present(region_a)) then
            buf_a = trim(region_a)//c_null_char
            opts%region_a = c_loc(buf_a(1:1))
        end if
        if (present(region_b)) then
            buf_b = trim(region_b)//c_null_char
            opts%region_b = c_loc(buf_b(1:1))
        end if
        s = c_mio_hausdorff_distance(self%handle, other%handle, opts, report)
        if (s /= 0_c_int) then
            call handle_failure('hausdorff_distance', mio_error_message(), stat, errmsg)
            return
        end if
        d = real(report%distance, real64)
        if (present(a_to_b)) a_to_b = real(report%a_to_b, real64)
        if (present(b_to_a)) b_to_a = real(report%b_to_a, real64)
        if (present(mean_a_to_b)) mean_a_to_b = real(report%mean_a_to_b, real64)
        if (present(mean_b_to_a)) mean_b_to_a = real(report%mean_b_to_a, real64)
        if (present(rms_a_to_b)) rms_a_to_b = real(report%rms_a_to_b, real64)
        if (present(rms_b_to_a)) rms_b_to_a = real(report%rms_b_to_a, real64)
        call clear_status(stat, errmsg)
    end function

    !> The C value of a region kind name ('point', 'cell', 'side'), or -1 for
    !> any other string.
    pure function region_kind_code(kind) result(k)
        character(*), intent(in) :: kind
        integer(c_int32_t) :: k
        select case (trim(kind))
        case ('point')
            k = 0
        case ('cell')
            k = 1
        case ('side')
            k = 2
        case default
            k = -1
        end select
    end function

    !> Apply one region edit to a copy of this mesh: `op` is 'union',
    !> 'intersection', 'difference' (two or more `inputs`, one kind, result
    !> named `output`), 'rename' (one input, to `output`), 'retag' (one input,
    !> new `tag` and/or `dim`) or 'delete'. `kind` pins the inputs' kind
    !> ('point', 'cell', 'side'). Points, cells and data are untouched. See
    !> doc/regions.md.
    function mesh_edit_regions(self, op, inputs, output, kind, dim, tag, keep_inputs, &
                               stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: op
        character(*), intent(in) :: inputs(:)
        character(*), intent(in), optional :: output, kind
        integer, intent(in), optional :: dim
        integer(int64), intent(in), optional :: tag
        logical, intent(in), optional :: keep_inputs
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr, res, out_ptr
        integer(c_int64_t) :: n, c_dim, c_tag
        type(mio_region_selector_t), allocatable :: sel(:)
        character(kind=c_char, len=STRBUF_LEN), target :: out_buf
        integer(c_int32_t) :: c_op, keep
        integer :: i

        select case (trim(op))
        case ('union')
            c_op = 0
        case ('intersection', 'intersect')
            c_op = 1
        case ('difference')
            c_op = 2
        case ('rename')
            c_op = 3
        case ('retag')
            c_op = 4
        case ('delete')
            c_op = 5
        case default
            call handle_failure('edit_regions', "meshio++: edit_regions: unknown operation '"// &
                                trim(op)//"'", stat, errmsg)
            return
        end select
        call c_str_array(inputs, storage, cptrs, arr, n)
        allocate (sel(max(size(inputs), 1)))
        do i = 1, size(inputs)
            sel(i)%name = cptrs(i)
            if (present(kind)) sel(i)%kind = region_kind_code(kind)
        end do
        out_ptr = c_null_ptr
        if (present(output)) then
            out_buf = trim(output)//c_null_char
            out_ptr = c_loc(out_buf(1:1))
        end if
        c_dim = -2_c_int64_t
        if (present(dim)) c_dim = int(dim, c_int64_t)
        c_tag = -2_c_int64_t
        if (present(tag)) c_tag = int(tag, c_int64_t)
        keep = 1
        if (present(keep_inputs)) keep = merge(1_c_int32_t, 0_c_int32_t, keep_inputs)
        res = c_mio_edit_regions(self%handle, c_op, sel, n, out_ptr, c_dim, c_tag, keep)
        if (.not. c_associated(res)) then
            call handle_failure('edit_regions', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        call clear_status(stat, errmsg)
    end function

    !> Return the conforming facets shared by named Cell regions. With no
    !> names, all Cell regions are used; when there are fewer than two, cell
    !> blocks become the groups. See doc/region_adjacency.md.
    function mesh_region_adjacency(self, regions, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in), optional :: regions(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr, out_ptr
        type(mio_region_selector_t), allocatable :: selectors(:)
        integer(c_int64_t) :: n
        integer :: i

        n = 0_c_int64_t
        if (present(regions)) then
            call c_str_array(regions, storage, cptrs, arr, n)
        else
            allocate (selectors(1))
        end if
        if (.not. allocated(selectors)) allocate (selectors(max(int(n), 1)))
        do i = 1, int(n)
            selectors(i)%name = cptrs(i)
            selectors(i)%kind = 1_c_int32_t
            selectors(i)%dim = -2_c_int64_t
            selectors(i)%tag = -2_c_int64_t
        end do
        out_ptr = c_mio_region_adjacency(self%handle, selectors, n)
        if (.not. c_associated(out_ptr)) then
            call handle_failure('region_adjacency', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = out_ptr
        call clear_status(stat, errmsg)
    end function

    !> Find conforming/proximity interfaces between Cell regions. `side_a` and
    !> `side_b` are (2,n) arrays of (global cell, local facet) ids; see doc/region_adjacency.md.
    function mesh_find_interface(self, region_a, region_b, mesh_b, mode, master, &
                                 gap_tolerance, angle_tolerance, overlap_tolerance, &
                                 stat, errmsg) result(result)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: region_a, region_b
        type(mio_mesh), intent(in), optional :: mesh_b
        character(*), intent(in), optional :: mode, master
        real(real64), intent(in), optional :: gap_tolerance, angle_tolerance, overlap_tolerance
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_find_interface_result) :: result
        type(mio_region_selector_t) :: selector_a, selector_b
        type(mio_find_interface_opts_t) :: opts
        character(kind=c_char, len=STRBUF_LEN), target :: buf_a, buf_b
        type(c_ptr) :: mesh_b_handle, result_handle, entries_ptr
        integer(c_int64_t), pointer :: entries(:)
        integer(c_int64_t) :: side_count
        integer(c_int) :: status

        if (.not. c_associated(self%handle)) then
            call handle_failure('find_interface', 'mesh is not open', stat, errmsg)
            return
        end if
        buf_a = trim(region_a)//c_null_char
        buf_b = trim(region_b)//c_null_char
        selector_a%name = c_loc(buf_a(1:1))
        selector_a%kind = 1_c_int32_t
        selector_b%name = c_loc(buf_b(1:1))
        selector_b%kind = 1_c_int32_t
        mesh_b_handle = c_null_ptr
        if (present(mesh_b)) mesh_b_handle = mesh_b%handle

        call c_mio_find_interface_opts_init(opts)
        if (present(mode)) then
            select case (trim(mode))
            case ('conforming')
                opts%mode = 0_c_int32_t
            case ('proximity')
                opts%mode = 1_c_int32_t
            case default
                call handle_failure('find_interface', &
                                    'mode must be "conforming" or "proximity"', stat, errmsg)
                return
            end select
        end if
        if (present(master)) then
            select case (trim(master))
            case ('a')
                opts%master = 0_c_int32_t
            case ('b')
                opts%master = 1_c_int32_t
            case default
                call handle_failure('find_interface', 'master must be "a" or "b"', stat, errmsg)
                return
            end select
        end if
        if (present(gap_tolerance)) opts%gap_tolerance = real(gap_tolerance, c_double)
        if (present(angle_tolerance)) opts%angle_tolerance = real(angle_tolerance, c_double)
        if (present(overlap_tolerance)) &
            opts%overlap_tolerance = real(overlap_tolerance, c_double)
        result_handle = c_mio_find_interface(self%handle, selector_a, mesh_b_handle, &
                                             selector_b, opts)
        if (.not. c_associated(result_handle)) then
            call handle_failure('find_interface', mio_error_message(), stat, errmsg)
            return
        end if
        result%mesh%handle = c_mio_find_interface_result_take_mesh(result_handle)
        if (.not. c_associated(result%mesh%handle)) then
            call c_mio_find_interface_result_free(result_handle)
            call handle_failure('find_interface', mio_error_message(), stat, errmsg)
            return
        end if
        status = c_mio_find_interface_result_report(result_handle, result%num_pairs, &
            result%area, result%max_gap, result%unmatched_a, result%unmatched_b)
        if (status /= 0_c_int) then
            call result%mesh%free()
            call c_mio_find_interface_result_free(result_handle)
            call handle_failure('find_interface', mio_error_message(), stat, errmsg)
            return
        end if

        entries_ptr = c_mio_find_interface_result_side_a(result_handle, side_count)
        allocate (result%side_a(2, int(side_count)))
        if (side_count > 0_c_int64_t) then
            call c_f_pointer(entries_ptr, entries, [2 * int(side_count)])
            result%side_a = reshape(entries, [2, int(side_count)])
        end if
        entries_ptr = c_mio_find_interface_result_side_b(result_handle, side_count)
        allocate (result%side_b(2, int(side_count)))
        if (side_count > 0_c_int64_t) then
            call c_f_pointer(entries_ptr, entries, [2 * int(side_count)])
            result%side_b = reshape(entries, [2, int(side_count)])
        end if
        call c_mio_find_interface_result_free(result_handle)
        call clear_status(stat, errmsg)
    end function

    !> Project a Point region to Cell-region facets; returned id arrays are
    !> 0-based and coordinate matrices have shape (3,n).
    function mesh_contact_pairs(self, slave_points, master_cells, master_mesh, tolerance, &
                                require_complete, stat, errmsg) result(result)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: slave_points, master_cells
        type(mio_mesh), intent(in), optional :: master_mesh
        real(real64), intent(in), optional :: tolerance
        logical, intent(in), optional :: require_complete
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_contact_pairs_result) :: result
        type(mio_region_selector_t) :: slave_selector, master_selector
        type(mio_contact_pairs_opts_t) :: opts
        character(kind=c_char, len=STRBUF_LEN), target :: slave_buf, master_buf
        type(c_ptr) :: master_handle, result_handle, array_ptr
        integer(c_int64_t), pointer :: values_i(:)
        real(c_double), pointer :: values_d(:)
        integer(c_int64_t) :: count, unmatched_count
        integer(c_int) :: status
        integer :: n

        if (.not. c_associated(self%handle)) then
            call handle_failure('contact_pairs', 'mesh is not open', stat, errmsg)
            return
        end if
        slave_buf = trim(slave_points)//c_null_char
        master_buf = trim(master_cells)//c_null_char
        slave_selector%name = c_loc(slave_buf(1:1))
        slave_selector%kind = 0_c_int32_t
        master_selector%name = c_loc(master_buf(1:1))
        master_selector%kind = 1_c_int32_t
        master_handle = c_null_ptr
        if (present(master_mesh)) master_handle = master_mesh%handle
        call c_mio_contact_pairs_opts_init(opts)
        if (present(tolerance)) opts%tolerance = real(tolerance, c_double)
        if (present(require_complete)) opts%require_complete = merge(1_c_int32_t, 0_c_int32_t, &
                                                                      require_complete)
        result_handle = c_mio_contact_pairs(self%handle, slave_selector, master_handle, &
                                            master_selector, opts)
        if (.not. c_associated(result_handle)) then
            call handle_failure('contact_pairs', mio_error_message(), stat, errmsg)
            return
        end if
        status = c_mio_contact_pairs_result_info(result_handle, count, unmatched_count)
        if (status /= 0_c_int) then
            call c_mio_contact_pairs_result_free(result_handle)
            call handle_failure('contact_pairs', mio_error_message(), stat, errmsg)
            return
        end if
        n = int(count)
        allocate (result%slave_point(n), result%master_cell(n), result%master_facet(n), &
                  result%master_subfacet(n), result%gap(n), result%local_coordinates(3, n), &
                  result%closest_point(3, n), result%normal(3, n), &
                  result%unmatched(int(unmatched_count)))

        if (n > 0) then
            array_ptr = c_mio_contact_pairs_slave_point(result_handle)
            call c_f_pointer(array_ptr, values_i, [n])
            result%slave_point = values_i
            array_ptr = c_mio_contact_pairs_master_cell(result_handle)
            call c_f_pointer(array_ptr, values_i, [n])
            result%master_cell = values_i
            array_ptr = c_mio_contact_pairs_master_facet(result_handle)
            call c_f_pointer(array_ptr, values_i, [n])
            result%master_facet = values_i
            array_ptr = c_mio_contact_pairs_master_subfacet(result_handle)
            call c_f_pointer(array_ptr, values_i, [n])
            result%master_subfacet = values_i
            array_ptr = c_mio_contact_pairs_gap(result_handle)
            call c_f_pointer(array_ptr, values_d, [n])
            result%gap = values_d
        end if
        if (unmatched_count > 0_c_int64_t) then
            array_ptr = c_mio_contact_pairs_unmatched(result_handle)
            call c_f_pointer(array_ptr, values_i, [int(unmatched_count)])
            result%unmatched = values_i
        end if
        if (n > 0) then
            array_ptr = c_mio_contact_pairs_local_coordinates(result_handle)
            call c_f_pointer(array_ptr, values_d, [3 * n])
            result%local_coordinates = reshape(values_d, [3, n])
            array_ptr = c_mio_contact_pairs_closest_point(result_handle)
            call c_f_pointer(array_ptr, values_d, [3 * n])
            result%closest_point = reshape(values_d, [3, n])
            array_ptr = c_mio_contact_pairs_normal(result_handle)
            call c_f_pointer(array_ptr, values_d, [3 * n])
            result%normal = reshape(values_d, [3, n])
        end if
        call c_mio_contact_pairs_result_free(result_handle)
        call clear_status(stat, errmsg)
    end function

    !> Split point fans along a Side region. Polyhedron input is rejected.
    function mesh_split_interface(self, side, add_cohesive, num_duplicated_points, &
                                  num_cohesive_cells, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: side
        logical, intent(in), optional :: add_cohesive
        integer(int64), intent(out), optional :: num_duplicated_points, num_cohesive_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_region_selector_t) :: selector
        type(mio_split_interface_opts_t) :: opts
        character(kind=c_char, len=STRBUF_LEN), target :: side_buf
        type(c_ptr) :: result_handle
        integer(c_int64_t) :: duplicated, cohesive

        side_buf = trim(side)//c_null_char
        selector%name = c_loc(side_buf(1:1))
        selector%kind = 2_c_int32_t
        call c_mio_split_interface_opts_init(opts)
        if (present(add_cohesive)) opts%add_cohesive = &
            merge(1_c_int32_t, 0_c_int32_t, add_cohesive)
        result_handle = c_mio_split_interface(self%handle, selector, opts, duplicated, cohesive)
        if (.not. c_associated(result_handle)) then
            call handle_failure('split_interface', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = result_handle
        if (present(num_duplicated_points)) num_duplicated_points = int(duplicated, int64)
        if (present(num_cohesive_cells)) num_cohesive_cells = int(cohesive, int64)
        call clear_status(stat, errmsg)
    end function

    !> Remove the `index`-th region (1-based, `regions()` order).
    subroutine mesh_remove_region(self, index, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (c_mio_mesh_remove_region(self%handle, int(index - 1, c_int64_t)) /= 0_c_int) then
            call handle_failure('remove_region', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end subroutine

    !> The master node each node of region `slave` maps onto, under the affine
    !> transform `translate` (3) or `matrix` (16, row-major 4x4), within `atol`.
    !> `slave_ids`/`master_ids` receive the pairs (1-based, ascending in the
    !> slave id). With `require_complete = .false.` unmatched slave nodes are
    !> returned in `unmatched` instead of failing. See doc/periodic.md.
    subroutine mesh_match_periodic_nodes(self, slave, master, slave_ids, master_ids, &
                                         translate, matrix, atol, require_complete, &
                                         unmatched, num_fixed, max_residual, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: slave, master
        integer(int64), allocatable, intent(out) :: slave_ids(:), master_ids(:)
        real(real64), intent(in), optional :: translate(3), matrix(16)
        real(real64), intent(in), optional :: atol
        logical, intent(in), optional :: require_complete
        integer(int64), allocatable, intent(out), optional :: unmatched(:)
        integer(int64), intent(out), optional :: num_fixed
        real(real64), intent(out), optional :: max_residual
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_periodic_opts_t) :: opts
        type(mio_region_selector_t) :: s_sel, m_sel
        character(kind=c_char, len=STRBUF_LEN), target :: s_buf, m_buf
        type(c_ptr) :: pairs, p
        integer(c_int64_t) :: np, nu, nf, cnt
        real(c_double) :: res
        integer(c_int64_t), pointer :: ids(:)

        call c_mio_periodic_opts_init(opts)
        if (present(matrix)) then
            opts%matrix = real(matrix, c_double)
        else if (present(translate)) then
            opts%matrix(4) = real(translate(1), c_double)
            opts%matrix(8) = real(translate(2), c_double)
            opts%matrix(12) = real(translate(3), c_double)
        end if
        if (present(atol)) opts%atol = real(atol, c_double)
        if (present(require_complete)) &
            opts%require_complete = merge(1_c_int32_t, 0_c_int32_t, require_complete)
        s_buf = trim(slave)//c_null_char
        m_buf = trim(master)//c_null_char
        s_sel%name = c_loc(s_buf(1:1))
        m_sel%name = c_loc(m_buf(1:1))
        pairs = c_mio_match_periodic_nodes(self%handle, s_sel, m_sel, opts)
        if (.not. c_associated(pairs)) then
            allocate (slave_ids(0), master_ids(0))
            call handle_failure('match_periodic_nodes', mio_error_message(), stat, errmsg)
            return
        end if
        if (c_mio_periodic_pairs_info(pairs, np, nu, nf, res) /= 0_c_int) then
            allocate (slave_ids(0), master_ids(0))
            call c_mio_periodic_pairs_free(pairs)
            call handle_failure('match_periodic_nodes', mio_error_message(), stat, errmsg)
            return
        end if
        allocate (slave_ids(np), master_ids(np))
        if (np > 0) then
            p = c_mio_periodic_pairs_slave(pairs, cnt)
            call c_f_pointer(p, ids, [cnt])
            slave_ids = int(ids, int64) + 1_int64
            p = c_mio_periodic_pairs_master(pairs, cnt)
            call c_f_pointer(p, ids, [cnt])
            master_ids = int(ids, int64) + 1_int64
        end if
        if (present(unmatched)) then
            allocate (unmatched(nu))
            if (nu > 0) then
                p = c_mio_periodic_pairs_unmatched(pairs, cnt)
                call c_f_pointer(p, ids, [cnt])
                unmatched = int(ids, int64) + 1_int64
            end if
        end if
        if (present(num_fixed)) num_fixed = int(nf, int64)
        if (present(max_residual)) max_residual = real(res, real64)
        call c_mio_periodic_pairs_free(pairs)
        call clear_status(stat, errmsg)
    end subroutine

    !> Quality gate: score the cells and test `spec` thresholds
    !> ("scaled_jacobian >= 0.2; aspect_ratio <= 5 @ 1%"), plus the inverted
    !> and degenerate counts (`max_inverted`/`max_degenerate`, default 0,
    !> negative disables). Returns whether every check passed; `summary`
    !> receives the text both CLIs print. See doc/mesh_quality.md.
    function mesh_check_quality(self, spec, max_inverted, max_degenerate, num_checks, &
                                num_failed, summary, stat, errmsg) result(passed)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in), optional :: spec
        integer(int64), intent(in), optional :: max_inverted, max_degenerate
        integer(int64), intent(out), optional :: num_checks, num_failed
        character(:), allocatable, intent(out), optional :: summary
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        logical :: passed
        type(mio_quality_gate_report_t) :: report
        character(kind=c_char), allocatable :: buf(:)
        integer(c_int64_t) :: mi, md
        integer :: s, i

        passed = .false.
        mi = 0_c_int64_t
        md = 0_c_int64_t
        if (present(max_inverted)) mi = int(max_inverted, c_int64_t)
        if (present(max_degenerate)) md = int(max_degenerate, c_int64_t)
        allocate (buf(4096))
        if (present(spec)) then
            s = c_mio_check_quality(self%handle, c_str(spec), mi, md, report, buf, 4096_c_int64_t)
        else
            s = c_mio_check_quality(self%handle, c_str(''), mi, md, report, buf, 4096_c_int64_t)
        end if
        if (s /= 0) then
            call handle_failure('check_quality', mio_error_message(), stat, errmsg)
            return
        end if
        passed = report%passed /= 0
        if (present(num_checks)) num_checks = int(report%num_checks, int64)
        if (present(num_failed)) num_failed = int(report%num_failed, int64)
        if (present(summary)) then
            allocate (character(int(min(report%summary_length, 4095_c_int64_t))) :: summary)
            do i = 1, len(summary)
                summary(i:i) = buf(i)
            end do
        end if
        call clear_status(stat, errmsg)
    end function

    !> Linear interpolation between two steps of one mesh (v16.25.0): `a` with
    !> every floating-point data array replaced by (1 - w) a + w b. The steps
    !> must share a topology; `blend_points` also blends the coordinates. See
    !> doc/sequences.md.
    function mio_blend_steps(a, b, w, blend_points, stat, errmsg) result(out)
        type(mio_mesh), intent(in) :: a, b
        real(real64), intent(in) :: w
        logical, intent(in), optional :: blend_points
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int32_t) :: bp
        bp = 0
        if (present(blend_points)) bp = merge(1_c_int32_t, 0_c_int32_t, blend_points)
        out%handle = c_mio_blend_steps(a%handle, b%handle, real(w, c_double), bp)
        if (.not. c_associated(out%handle)) then
            call handle_failure('blend_steps', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Point and/or cell normals of this surface, optionally splitting
    !> vertices at creases so every point carries exactly one normal.
    !>
    !> Writes `normals` as point data (n, 3) and, with `cell_normals`, as
    !> cell data. With `split_angle` present the surface is split wherever it
    !> creases by more than that many degrees (in [0, 180]): the copies are
    !> appended after the original points, cells keep their numbering, and
    !> `record_parent_ids` adds `normals:parent_point`. `weight` is 'angle'
    !> (default) or 'area'.
    !>
    !> Never reorients: a nonzero `inconsistent_pairs` means some normals
    !> average faces that disagree about which side is out. A volume block is
    !> refused by name pointing at `extract_surface`.
    function mesh_normals(self, point_normals, cell_normals, weight, split_angle, &
                          record_parent_ids, region, num_isolated, num_undefined, &
                          num_degenerate, num_split_points, num_added_points, &
                          boundary_edges, non_manifold_edges, inconsistent_pairs, &
                          degenerate_triangles, watertight, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        logical, intent(in), optional :: point_normals, cell_normals, record_parent_ids
        character(*), intent(in), optional :: weight, region
        real(real64), intent(in), optional :: split_angle
        integer(int64), intent(out), optional :: num_isolated, num_undefined
        integer(int64), intent(out), optional :: num_degenerate, num_split_points
        integer(int64), intent(out), optional :: num_added_points
        integer(int64), intent(out), optional :: boundary_edges, non_manifold_edges
        integer(int64), intent(out), optional :: inconsistent_pairs, degenerate_triangles
        logical, intent(out), optional :: watertight
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_normals_opts_t) :: opts
        type(mio_normals_report_t) :: report
        character(kind=c_char, len=STRBUF_LEN), target :: region_buf

        call c_mio_normals_opts_init(opts)
        if (present(point_normals)) then
            if (.not. point_normals) opts%point_normals = 0
        end if
        if (present(cell_normals)) then
            if (cell_normals) opts%cell_normals = 1
        end if
        if (present(weight)) then
            if (trim(weight) == 'area') then
                opts%weight = 1
            else if (trim(weight) == 'angle') then
                opts%weight = 0
            else
                call handle_failure('normals', &
                    "meshio++: normals: unknown weight '"//trim(weight)// &
                    "' (expected 'angle' or 'area')", stat, errmsg)
                return
            end if
        end if
        if (present(split_angle)) then
            opts%split = 1
            opts%split_angle = real(split_angle, c_double)
        end if
        if (present(record_parent_ids)) then
            if (record_parent_ids) opts%record_parent_ids = 1
        end if
        if (present(region)) then
            region_buf = trim(region)//c_null_char
            opts%region = c_loc(region_buf(1:1))
        end if

        res = c_mio_compute_normals(self%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('normals', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_isolated)) num_isolated = int(report%num_isolated, int64)
        if (present(num_undefined)) num_undefined = int(report%num_undefined, int64)
        if (present(num_degenerate)) num_degenerate = int(report%num_degenerate, int64)
        if (present(num_split_points)) &
            num_split_points = int(report%num_split_points, int64)
        if (present(num_added_points)) &
            num_added_points = int(report%num_added_points, int64)
        if (present(boundary_edges)) &
            boundary_edges = int(report%quality%boundary_edges, int64)
        if (present(non_manifold_edges)) &
            non_manifold_edges = int(report%quality%non_manifold_edges, int64)
        if (present(inconsistent_pairs)) &
            inconsistent_pairs = int(report%quality%inconsistent_pairs, int64)
        if (present(degenerate_triangles)) &
            degenerate_triangles = int(report%quality%degenerate_triangles, int64)
        if (present(watertight)) watertight = (report%quality%watertight /= 0)
        call clear_status(stat, errmsg)
    end function

    !> Repair this surface's orientation, holes and pinched vertices: weld
    !> (opt-in) -> triangulate (blocks 1:1) -> split bowties -> orient by the
    !> topological half-edge rule per component -> fan-fill boundary loops of
    !> at most `max_hole_edges` edges -> orient closed components outward.
    !> Lower-dimensional blocks ride along; fill triangles land in one
    !> trailing `triangle` block. Every counter is an optional out-argument.
    function mesh_repair(self, fix_orientation, orient_outward, fill_holes, &
                         split_non_manifold, max_hole_edges, weld_tolerance, &
                         record_provenance, num_flipped, num_components, &
                         num_holes_detected, num_holes_filled, num_holes_skipped, &
                         num_faces_added, num_points_added, num_vertices_split, &
                         num_oriented_outward, num_unorientable, points_welded, &
                         watertight_before, watertight_after, inconsistent_pairs_after, &
                         boundary_edges_after, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        logical, intent(in), optional :: fix_orientation, orient_outward, fill_holes
        logical, intent(in), optional :: split_non_manifold, record_provenance
        integer(int64), intent(in), optional :: max_hole_edges
        real(real64), intent(in), optional :: weld_tolerance
        integer(int64), intent(out), optional :: num_flipped, num_components
        integer(int64), intent(out), optional :: num_holes_detected, num_holes_filled
        integer(int64), intent(out), optional :: num_holes_skipped, num_faces_added
        integer(int64), intent(out), optional :: num_points_added, num_vertices_split
        integer(int64), intent(out), optional :: num_oriented_outward, num_unorientable
        integer(int64), intent(out), optional :: points_welded
        logical, intent(out), optional :: watertight_before, watertight_after
        integer(int64), intent(out), optional :: inconsistent_pairs_after, boundary_edges_after
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_repair_opts_t) :: opts
        type(mio_repair_report_t) :: report

        call c_mio_repair_opts_init(opts)
        if (present(fix_orientation)) then
            if (.not. fix_orientation) opts%fix_orientation = 0
        end if
        if (present(orient_outward)) then
            if (.not. orient_outward) opts%orient_outward = 0
        end if
        if (present(fill_holes)) then
            if (.not. fill_holes) opts%fill_holes = 0
        end if
        if (present(split_non_manifold)) then
            if (.not. split_non_manifold) opts%split_non_manifold = 0
        end if
        if (present(record_provenance)) then
            if (record_provenance) opts%record_provenance = 1
        end if
        if (present(max_hole_edges)) opts%max_hole_edges = int(max_hole_edges, c_int64_t)
        if (present(weld_tolerance)) opts%weld_tolerance = real(weld_tolerance, c_double)

        res = c_mio_repair(self%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('repair', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_flipped)) num_flipped = int(report%num_flipped, int64)
        if (present(num_components)) num_components = int(report%num_components, int64)
        if (present(num_holes_detected)) num_holes_detected = int(report%num_holes_detected, int64)
        if (present(num_holes_filled)) num_holes_filled = int(report%num_holes_filled, int64)
        if (present(num_holes_skipped)) num_holes_skipped = int(report%num_holes_skipped, int64)
        if (present(num_faces_added)) num_faces_added = int(report%num_faces_added, int64)
        if (present(num_points_added)) num_points_added = int(report%num_points_added, int64)
        if (present(num_vertices_split)) num_vertices_split = int(report%num_vertices_split, int64)
        if (present(num_oriented_outward)) &
            num_oriented_outward = int(report%num_oriented_outward, int64)
        if (present(num_unorientable)) num_unorientable = int(report%num_unorientable, int64)
        if (present(points_welded)) points_welded = int(report%points_welded, int64)
        if (present(watertight_before)) watertight_before = (report%quality_before%watertight /= 0)
        if (present(watertight_after)) watertight_after = (report%quality_after%watertight /= 0)
        if (present(inconsistent_pairs_after)) &
            inconsistent_pairs_after = int(report%quality_after%inconsistent_pairs, int64)
        if (present(boundary_edges_after)) &
            boundary_edges_after = int(report%quality_after%boundary_edges, int64)
        call clear_status(stat, errmsg)
    end function

    !> Sobolev (Helmholtz-filtered) deformation: solve `(M + l^2 K) u = M d`
    !> over this mesh's own P1 operators and move the points by `u`. `array`
    !> names the `(n, dim)` point-data displacement; `fixed_points_array` an
    !> optional point-data array whose nonzero entries pin their point. Every
    !> top-dimensional block must be a linear simplex. Non-convergence is
    !> reported through `converged`, never raised.
    function mesh_sobolev_deform(self, array, length_scale, fixed_points_array, &
                                 fix_boundary, record_filtered, max_iterations, &
                                 tolerance, num_iterations, residual, converged, &
                                 num_fixed, num_isolated, max_displacement, &
                                 stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: array
        real(real64), intent(in) :: length_scale
        character(*), intent(in), optional :: fixed_points_array
        logical, intent(in), optional :: fix_boundary, record_filtered
        integer, intent(in), optional :: max_iterations
        real(real64), intent(in), optional :: tolerance
        integer(int64), intent(out), optional :: num_iterations, num_fixed, num_isolated
        real(real64), intent(out), optional :: residual, max_displacement
        logical, intent(out), optional :: converged
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_sobolev_opts_t) :: opts
        type(mio_sobolev_report_t) :: report
        character(kind=c_char, len=STRBUF_LEN), target :: array_buf, fixed_buf

        call c_mio_sobolev_opts_init(opts)
        array_buf = trim(array)//c_null_char
        opts%array = c_loc(array_buf(1:1))
        opts%length_scale = real(length_scale, c_double)
        if (present(fixed_points_array)) then
            fixed_buf = trim(fixed_points_array)//c_null_char
            opts%fixed_points_array = c_loc(fixed_buf(1:1))
        end if
        if (present(fix_boundary)) then
            if (fix_boundary) opts%fix_boundary = 1
        end if
        if (present(record_filtered)) then
            if (record_filtered) opts%record_filtered = 1
        end if
        if (present(max_iterations)) opts%max_iterations = int(max_iterations, c_int32_t)
        if (present(tolerance)) opts%tolerance = real(tolerance, c_double)

        res = c_mio_sobolev_deform(self%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('sobolev_deform', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_iterations)) num_iterations = int(report%num_iterations, int64)
        if (present(num_fixed)) num_fixed = int(report%num_fixed, int64)
        if (present(num_isolated)) num_isolated = int(report%num_isolated, int64)
        if (present(residual)) residual = real(report%residual, real64)
        if (present(max_displacement)) max_displacement = real(report%max_displacement, real64)
        if (present(converged)) converged = (report%converged /= 0)
        call clear_status(stat, errmsg)
    end function

    !> Project every (selected) point of `mesh` onto the surface of `target`:
    !> `x' = x + w (p + offset n - x)`, one projection, no iteration. Every
    !> point of the source moves whatever cells it carries; only the target
    !> must be a surface. `weights` names a point-data array on the source
    !> (nonzero selects, a float value blends); `normal_weight` is 'angle'
    !> (default) or 'area'. Module-level like `mio_interpolate`: two meshes in.
    function mio_shrinkwrap(mesh, target, offset, max_distance, weights, target_region, &
                            normal_weight, record_distance, record_closest_cell, &
                            num_projected, num_missed, num_skipped, max_displacement, &
                            target_watertight, stat, errmsg) result(out)
        type(mio_mesh), intent(in) :: mesh, target
        real(real64), intent(in), optional :: offset, max_distance
        character(*), intent(in), optional :: weights, target_region, normal_weight
        logical, intent(in), optional :: record_distance, record_closest_cell
        integer(int64), intent(out), optional :: num_projected, num_missed, num_skipped
        real(real64), intent(out), optional :: max_displacement
        logical, intent(out), optional :: target_watertight
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_shrinkwrap_opts_t) :: opts
        type(mio_shrinkwrap_report_t) :: report
        character(kind=c_char, len=STRBUF_LEN), target :: weights_buf, region_buf

        call c_mio_shrinkwrap_opts_init(opts)
        if (present(offset)) opts%offset = real(offset, c_double)
        if (present(max_distance)) opts%max_distance = real(max_distance, c_double)
        if (present(weights)) then
            weights_buf = trim(weights)//c_null_char
            opts%weights = c_loc(weights_buf(1:1))
        end if
        if (present(target_region)) then
            region_buf = trim(target_region)//c_null_char
            opts%target_region = c_loc(region_buf(1:1))
        end if
        if (present(normal_weight)) then
            if (trim(normal_weight) == 'area') then
                opts%normal_weight = 1
            else if (trim(normal_weight) == 'angle') then
                opts%normal_weight = 0
            else
                call handle_failure('shrinkwrap', &
                    "meshio++: shrinkwrap: unknown normal weight '"//trim(normal_weight)// &
                    "' (expected 'angle' or 'area')", stat, errmsg)
                return
            end if
        end if
        if (present(record_distance)) then
            if (record_distance) opts%record_distance = 1
        end if
        if (present(record_closest_cell)) then
            if (record_closest_cell) opts%record_closest_cell = 1
        end if

        res = c_mio_shrinkwrap(mesh%handle, target%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('shrinkwrap', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_projected)) num_projected = int(report%num_projected, int64)
        if (present(num_missed)) num_missed = int(report%num_missed, int64)
        if (present(num_skipped)) num_skipped = int(report%num_skipped, int64)
        if (present(max_displacement)) max_displacement = real(report%max_displacement, real64)
        if (present(target_watertight)) target_watertight = (report%quality%watertight /= 0)
        call clear_status(stat, errmsg)
    end function

    function mesh_remesh(self, num_clusters, subdivide, subsample_ratio, max_subdivide, &
                         max_iterations, max_repair_passes, metric, gradation, &
                         preserve_boundary, max_anisotropy, num_iterations, subdivide_applied, &
                         num_isolated_clusters, num_non_manifold_vertices, stat, errmsg) &
            result(out)
        class(mio_mesh), intent(in) :: self
        integer(int64), intent(in) :: num_clusters
        integer, intent(in), optional :: subdivide, max_subdivide, max_iterations
        integer, intent(in), optional :: max_repair_passes
        real(real64), intent(in), optional :: subsample_ratio
        character(*), intent(in), optional :: metric
        real(real64), intent(in), optional :: gradation
        logical, intent(in), optional :: preserve_boundary
        real(real64), intent(in), optional :: max_anisotropy
        integer(int64), intent(out), optional :: num_iterations, num_isolated_clusters
        integer, intent(out), optional :: subdivide_applied
        integer(int64), intent(out), optional :: num_non_manifold_vertices
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_remesh_opts_t) :: opts
        type(mio_remesh_report_t) :: report
        ! NUL-terminated copy must outlive the call, so it is held here
        ! rather than built inline; c_loc needs it contiguous and TARGET --
        ! the region_buf/array_buf idiom mesh_refine already uses.
        character(kind=c_char, len=STRBUF_LEN), target :: metric_buf

        call c_mio_remesh_opts_init(opts)
        opts%num_clusters = int(num_clusters, c_int64_t)
        if (present(subdivide)) opts%subdivide = int(subdivide, c_int32_t)
        if (present(subsample_ratio)) opts%subsample_ratio = real(subsample_ratio, c_double)
        if (present(max_subdivide)) opts%max_subdivide = int(max_subdivide, c_int32_t)
        if (present(max_iterations)) opts%max_iterations = int(max_iterations, c_int32_t)
        if (present(max_repair_passes)) opts%max_repair_passes = int(max_repair_passes, c_int32_t)
        if (present(metric)) then
            metric_buf = trim(metric)//c_null_char
            opts%metric = c_loc(metric_buf(1:1))
        end if
        if (present(gradation)) opts%gradation = real(gradation, c_double)
        if (present(preserve_boundary)) then
            if (.not. preserve_boundary) opts%preserve_boundary = 0
        end if
        if (present(max_anisotropy)) opts%max_anisotropy = real(max_anisotropy, c_double)

        res = c_mio_remesh_ex(self%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('remesh', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_iterations)) num_iterations = int(report%num_iterations, int64)
        if (present(subdivide_applied)) subdivide_applied = int(report%subdivide_applied)
        if (present(num_isolated_clusters)) &
            num_isolated_clusters = int(report%num_isolated_clusters, int64)
        if (present(num_non_manifold_vertices)) &
            num_non_manifold_vertices = int(report%num_non_manifold_vertices, int64)
        call clear_status(stat, errmsg)
    end function

    !> Retetrahedralize this mesh (a volume mesh, or a closed surface) at a
    !> caller-chosen resolution by isosurface stuffing -- remesh's volumetric
    !> sibling. Unlike `remesh`, `self` may be a VOLUME mesh directly (its
    !> boundary is extracted internally) as well as a closed surface. Same
    !> no-correspondence-with-`self` output contract: new points, new
    !> connectivity, point/cell data and named regions dropped, field data
    !> carried.
    !>
    !> Exactly one of `resolution`/`cell_size` must be given. `warp_fraction`
    !> (default 0.35) trades a small, measured chance of non-manifold
    !> boundary edges (reported in `num_non_manifold_edges`) for
    !> substantially better boundary tet quality; 0 disables warping. See
    !> doc/remesh_volume.md for the measured tradeoff.
    function mesh_remesh_volume(self, resolution, cell_size, bounds, padding, &
                                padding_relative, max_cells, max_tets, warp_fraction, sign, &
                                watertight_check, num_vertices_warped, num_tets_rejected, &
                                num_non_manifold_edges, quality, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer, intent(in), optional :: resolution(3)
        real(real64), intent(in), optional :: cell_size, bounds(6)
        real(real64), intent(in), optional :: padding, padding_relative, warp_fraction
        integer(int64), intent(in), optional :: max_cells, max_tets
        character(*), intent(in), optional :: sign, watertight_check
        integer(int64), intent(out), optional :: num_vertices_warped, num_tets_rejected
        integer(int64), intent(out), optional :: num_non_manifold_edges
        type(mio_surface_quality), intent(out), optional :: quality
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_remesh_volume_opts_t) :: opts
        type(mio_remesh_volume_report_t) :: report
        ! The buffers the option pointers reference must outlive the call,
        ! the mesh_compute_sdf idiom above.
        integer(c_int64_t), target :: res_buf(3)
        real(c_double), target :: bounds_buf(6)
        integer(c_int32_t) :: code

        call c_mio_remesh_volume_opts_init(opts)
        if (present(resolution)) then
            res_buf = int(resolution, c_int64_t)
            opts%resolution = c_loc(res_buf(1))
        end if
        if (present(cell_size)) opts%cell_size = real(cell_size, c_double)
        if (present(bounds)) then
            bounds_buf = real(bounds, c_double)
            opts%bounds = c_loc(bounds_buf(1))
        end if
        if (present(padding)) opts%padding = real(padding, c_double)
        if (present(padding_relative)) opts%padding_relative = real(padding_relative, c_double)
        if (present(max_cells)) opts%max_cells = int(max_cells, c_int64_t)
        if (present(max_tets)) opts%max_tets = int(max_tets, c_int64_t)
        if (present(warp_fraction)) opts%warp_fraction = real(warp_fraction, c_double)
        if (present(sign)) then
            code = sdf_sign_code(sign)
            if (code < 0) then
                call handle_failure('remesh_volume', 'unknown sign '//trim(sign), stat, errmsg)
                return
            end if
            opts%distance%sign = code
        end if
        if (present(watertight_check)) then
            select case (trim(watertight_check))
            case ('off'); opts%distance%watertight_check = 0
            case ('', 'warn'); opts%distance%watertight_check = 1
            case ('error'); opts%distance%watertight_check = 2
            case default
                call handle_failure('remesh_volume', &
                                    'unknown watertight check '//trim(watertight_check), &
                                    stat, errmsg)
                return
            end select
        end if

        res = c_mio_remesh_volume_ex(self%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('remesh_volume', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_vertices_warped)) &
            num_vertices_warped = int(report%num_vertices_warped, int64)
        if (present(num_tets_rejected)) num_tets_rejected = int(report%num_tets_rejected, int64)
        if (present(num_non_manifold_edges)) &
            num_non_manifold_edges = int(report%num_non_manifold_edges, int64)
        if (present(quality)) quality = report%input_quality
        call clear_status(stat, errmsg)
    end function

    !> Convert the element representation of the mesh. `mode` is "linearize"
    !> (higher-order cells -> their linear base, pruning the orphaned nodes),
    !> "simplexify" (decompose into same-dimension simplices), or "elevate"
    !> (linear -> serendipity quadratic, adding one node per unique edge).
    !> `record_parent_ids` (default .false.) attaches a convert:parent_cell
    !> cell_data array. The optional `point_map` receives the input-point ->
    !> output-point mapping, 1-based, or 0 where the point was pruned.
    function mesh_convert_cells(self, mode, record_parent_ids, point_map, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: mode
        logical, intent(in), optional :: record_parent_ids
        integer(int64), allocatable, intent(out), optional :: point_map(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res, cdata
        integer(c_int) :: crec, s, dt
        integer(c_int64_t) :: nlen
        integer(c_int64_t), pointer :: fp(:)
        crec = 0
        if (present(record_parent_ids)) then
            if (record_parent_ids) crec = 1
        end if
        res = c_mio_convert_cells(self%handle, c_str(mode), crec)
        if (.not. c_associated(res)) then
            call handle_failure('convert_cells', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(point_map)) then
            s = c_mio_convert_cells_result_point_map(res, cdata, dt, nlen)
            if (s /= 0_c_int) then
                call c_mio_convert_cells_result_free(res)
                call handle_failure('convert_cells', mio_error_message(), stat, errmsg)
                return
            end if
            allocate (point_map(nlen))
            if (nlen > 0) then
                call c_f_pointer(cdata, fp, [nlen])
                point_map = int(fp, int64) + 1_int64  ! 0-based (-1 = pruned) -> 1-based (0 = pruned)
            end if
        end if
        out%handle = c_mio_convert_cells_result_take_mesh(res)
        call c_mio_convert_cells_result_free(res)
        if (.not. c_associated(out%handle)) then
            call handle_failure('convert_cells', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Polyhedrally refine the mesh: one polyhedral child per face of every
    !> eligible 3D cell, connected to a new interior point. Needs no per-type
    !> template table -- tabulated types (reduced to corners for a quadratic
    !> variant) and existing polyhedron blocks are handled uniformly.
    !> Automatically conforming, unlike `refine`. Non-3D blocks and the
    !> full-Lagrange family (no face table) pass through unchanged.
    !> `record_parent_ids` (default .false.) attaches a
    !> subdivide:parent_cell cell_data array. Unlike `convert_cells`, there is
    !> no point map to request -- subdivide never prunes or renumbers an
    !> original point.
    function mesh_subdivide(self, record_parent_ids, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        logical, intent(in), optional :: record_parent_ids
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res
        integer(c_int) :: crec
        crec = 0
        if (present(record_parent_ids)) then
            if (record_parent_ids) crec = 1
        end if
        res = c_mio_subdivide(self%handle, crec)
        if (.not. c_associated(res)) then
            call handle_failure('subdivide', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = c_mio_subdivide_result_take_mesh(res)
        call c_mio_subdivide_result_free(res)
        if (.not. c_associated(out%handle)) then
            call handle_failure('subdivide', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Polyhedrally coarsen the mesh: merge groups of cells into single
    !> larger polyhedral cells via greedy seed-and-grow over the shared-face
    !> dual, absorbing face-adjacent neighbours until a group reaches
    !> `target_group_size` (default 8; a short group at a mesh boundary or
    !> pocket is expected, not an error). Non-volume blocks pass through
    !> unchanged; points are never pruned or renumbered (`clean` with
    !> `remove_orphans=.true.` is the follow-up for a minimal point set).
    !> `target_group_size=1` groups every cell by itself.
    function mesh_agglomerate(self, target_group_size, stat, errmsg, merge_coplanar_faces, &
                              coplanar_angle, min_sphericity, num_faces_merged, &
                              num_rejected) result(out)
        class(mio_mesh), intent(in) :: self
        integer(int64), intent(in), optional :: target_group_size
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        logical, intent(in), optional :: merge_coplanar_faces
        real(real64), intent(in), optional :: coplanar_angle, min_sphericity
        integer(int64), intent(out), optional :: num_faces_merged, num_rejected
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_agglomerate_opts_t) :: opts
        integer(c_int64_t) :: nmerged, nrejected
        if (present(target_group_size)) opts%target_group_size = int(target_group_size, c_int64_t)
        if (present(merge_coplanar_faces)) &
            opts%merge_coplanar_faces = merge(1_c_int32_t, 0_c_int32_t, merge_coplanar_faces)
        if (present(coplanar_angle)) opts%coplanar_angle = real(coplanar_angle, c_double)
        if (present(min_sphericity)) opts%min_sphericity = real(min_sphericity, c_double)
        res = c_mio_agglomerate_ex(self%handle, opts, nmerged, nrejected)
        if (present(num_faces_merged)) num_faces_merged = int(nmerged, int64)
        if (present(num_rejected)) num_rejected = int(nrejected, int64)
        if (.not. c_associated(res)) then
            call handle_failure('agglomerate', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = c_mio_agglomerate_result_take_mesh(res)
        call c_mio_agglomerate_result_free(res)
        if (.not. c_associated(out%handle)) then
            call handle_failure('agglomerate', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Uniformly refine the mesh, subdividing every cell into same-type children
    !> (line -> 2, triangle -> 4, quad -> 4, tetra -> 8, wedge -> 8,
    !> hexahedron -> 8). New nodes sit at edge / quad-face / body midpoints and
    !> are shared between neighbouring cells, so the result has no hanging
    !> nodes. `levels` (default 1) applies the templates repeatedly; 0 or less
    !> returns an unchanged copy. `record_parent_ids` (default .false.) attaches
    !> a refine:parent_cell cell_data array naming each output cell's ORIGINAL
    !> ancestor. The optional `point_map` receives the input-point ->
    !> output-point mapping, 1-based (the identity: refinement never prunes).
    !> Higher-order cells, pyramids and ragged blocks have no same-type
    !> subdivision and fail.
    !> The `mio_refine_compare` code for an operator spelling, or -1 if unknown.
    function refine_compare_code(op) result(code)
        character(*), intent(in) :: op
        integer(c_int32_t) :: code
        select case (trim(op))
        case ('<', 'lt'); code = 0
        case ('<=', 'le'); code = 1
        case ('>', 'gt'); code = 2
        case ('>=', 'ge'); code = 3
        case ('==', '=', 'eq'); code = 4
        case ('/=', '!=', 'ne'); code = 5
        case default; code = -1
        end select
    end function

    !> The `mio_sdf_sign` code for a sign name, or -1 if unknown.
    function sdf_sign_code(name) result(code)
        character(*), intent(in) :: name
        integer(c_int32_t) :: code
        select case (trim(name))
        case ('unsigned'); code = 0
        case ('', 'pseudonormal'); code = 1
        case ('winding-number', 'winding_number'); code = 2
        case default; code = -1
        end select
    end function

    !> The `mio_voxel_fill` code for a fill name, or -1 if unknown.
    function voxel_fill_code(name) result(code)
        character(*), intent(in) :: name
        integer(c_int32_t) :: code
        select case (trim(name))
        case ('', 'all'); code = 0
        case ('surface'); code = 1
        case ('inside'); code = 2
        case default; code = -1
        end select
    end function

    !> Build a regular hexahedron lattice from nothing -- the one entry point
    !> here that takes no input mesh. `dims` is the cell count per axis;
    !> `origin` (default 0) and `spacing` (default 1) place it. A zero count on
    !> any axis gives an empty mesh, which is a legal request rather than an
    !> error. `max_cells` (default 20000000) refuses a grid larger than that.
    function mio_grid(dims, origin, spacing, max_cells, stat, errmsg) result(out)
        integer, intent(in) :: dims(3)
        real(real64), intent(in), optional :: origin(3), spacing(3)
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int64_t) :: cdims(3), cmax
        real(c_double) :: corigin(3), cspacing(3)
        cdims = int(dims, c_int64_t)
        corigin = 0.0_c_double
        cspacing = 1.0_c_double
        if (present(origin)) corigin = real(origin, c_double)
        if (present(spacing)) cspacing = real(spacing, c_double)
        cmax = 20000000_c_int64_t
        if (present(max_cells)) cmax = int(max_cells, c_int64_t)
        out%handle = c_mio_grid(cdims, corigin, cspacing, cmax)
        if (.not. c_associated(out%handle)) then
            call handle_failure('grid', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> The layout of `mio_gen_opts_t` against the C struct (160 bytes), checked
    !> once per process.
    subroutine check_gen_opts_layout()
        logical, save :: checked = .false.
        type(mio_gen_opts_t) :: probe
        if (checked) return
        checked = .true.
        if (c_sizeof(probe) /= 160_c_size_t) then
            write (error_unit, '(a,i0,a)') &
                'meshio++: mio_gen_opts_t layout mismatch (', c_sizeof(probe), ' bytes)'
            error stop 1
        end if
    end subroutine

    !> The `mio_gen_cell` code of a cell type name, or -1 for an unknown name.
    !> An empty name is the shape's default.
    pure function gen_cell_code(name) result(code)
        character(*), intent(in) :: name
        integer(c_int32_t) :: code
        select case (trim(name))
        case (''); code = 0
        case ('quad'); code = 1
        case ('triangle'); code = 2
        case ('hexahedron'); code = 3
        case ('tetra'); code = 4
        case ('wedge'); code = 5
        case default; code = -1
        end select
    end function

    !> Fill the options shared by the generators and run the one named by `which`.
    function gen_run(which, opts, cell_type, split, order, max_cells, stat, errmsg) result(out)
        character(*), intent(in) :: which
        type(mio_gen_opts_t), intent(inout) :: opts
        character(*), intent(in), optional :: cell_type, split
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int32_t) :: code
        call check_gen_opts_layout()
        if (present(order)) opts%order = int(order, c_int32_t)
        if (present(max_cells)) opts%max_cells = int(max_cells, c_int64_t)
        if (present(cell_type)) then
            code = gen_cell_code(cell_type)
            if (code < 0) then
                call handle_failure(which, 'unknown cell type '//trim(cell_type), stat, errmsg)
                return
            end if
            opts%cell_type = code
        end if
        if (present(split)) then
            select case (trim(split))
            case ('', 'zigzag'); opts%split = 0
            case ('up'); opts%split = 1
            case ('down'); opts%split = 2
            case ('center'); opts%split = 3
            case default
                call handle_failure(which, 'unknown split '//trim(split), stat, errmsg)
                return
            end select
        end if
        select case (which)
        case ('rectangle'); out%handle = c_mio_rectangle(opts)
        case ('box'); out%handle = c_mio_box(opts)
        case ('disk'); out%handle = c_mio_disk(opts)
        case ('uv_sphere'); out%handle = c_mio_uv_sphere(opts)
        case ('icosa_sphere'); out%handle = c_mio_icosa_sphere(opts)
        case ('reference_triangle'); out%handle = c_mio_reference_triangle(opts)
        case ('ngon'); out%handle = c_mio_ngon(opts)
        case ('annulus'); out%handle = c_mio_annulus(opts)
        case ('lshape'); out%handle = c_mio_lshape(opts)
        case ('plate_with_hole'); out%handle = c_mio_plate_with_hole(opts)
        case ('honeycomb'); out%handle = c_mio_honeycomb(opts)
        case default; out%handle = c_mio_cylinder(opts)
        end select
        if (.not. c_associated(out%handle)) then
            call handle_failure(which, mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> A regular 2-D grid of `dims` = (nx, ny) quads or triangles at
    !> z = origin(3). `cell_type` is 'quad' (default) or 'triangle'; `split`
    !> ('zigzag', 'up', 'down', 'center') picks the triangle layout; `order` 2
    !> gives quadratic cells. Side regions xmin, xmax, ymin, ymax.
    function mio_rectangle(dims, lengths, origin, cell_type, split, order, max_cells, &
                           stat, errmsg) result(out)
        integer, intent(in), optional :: dims(2)
        real(real64), intent(in), optional :: lengths(2), origin(3)
        character(*), intent(in), optional :: cell_type, split
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        opts%counts = 1
        if (present(dims)) opts%counts(1:2) = int(dims, c_int64_t)
        if (present(lengths)) opts%lengths(1:2) = real(lengths, c_double)
        if (present(origin)) opts%origin = real(origin, c_double)
        out = gen_run('rectangle', opts, cell_type, split, order, max_cells, stat, errmsg)
    end function

    !> A regular 3-D grid of `dims` = (nx, ny, nz) hexahedra (default), 'tetra'
    !> (six per cell) or 'wedge' cells. Side regions xmin ... zmax.
    function mio_box(dims, lengths, origin, cell_type, order, max_cells, stat, errmsg) &
        result(out)
        integer, intent(in), optional :: dims(3)
        real(real64), intent(in), optional :: lengths(3), origin(3)
        character(*), intent(in), optional :: cell_type
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        opts%counts = 1
        if (present(dims)) opts%counts = int(dims, c_int64_t)
        if (present(lengths)) opts%lengths = real(lengths, c_double)
        if (present(origin)) opts%origin = real(origin, c_double)
        out = gen_run('box', opts, cell_type, order=order, max_cells=max_cells, &
                      stat=stat, errmsg=errmsg)
    end function

    !> A disk of `n` x `n` quads (default) or `n` rings of triangles; side
    !> region rim. `center` places it.
    function mio_disk(radius, n, center, cell_type, order, max_cells, stat, errmsg) &
        result(out)
        real(real64), intent(in), optional :: radius, center(3)
        integer, intent(in), optional :: n
        character(*), intent(in), optional :: cell_type
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        if (present(radius)) opts%radius = real(radius, c_double)
        if (present(n)) opts%counts(1) = int(n, c_int64_t)
        if (present(center)) opts%origin = real(center, c_double)
        out = gen_run('disk', opts, cell_type, order=order, max_cells=max_cells, &
                      stat=stat, errmsg=errmsg)
    end function

    !> A latitude / longitude sphere surface with `nu` points per ring and `nv`
    !> bands, oriented outward. `cell_type` is 'triangle' (default) or 'quad'.
    function mio_uv_sphere(radius, nu, nv, center, cell_type, order, max_cells, stat, errmsg) &
        result(out)
        real(real64), intent(in), optional :: radius, center(3)
        integer, intent(in), optional :: nu, nv
        character(*), intent(in), optional :: cell_type
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        opts%counts(1) = 16
        opts%counts(2) = 8
        if (present(radius)) opts%radius = real(radius, c_double)
        if (present(nu)) opts%counts(1) = int(nu, c_int64_t)
        if (present(nv)) opts%counts(2) = int(nv, c_int64_t)
        if (present(center)) opts%origin = real(center, c_double)
        out = gen_run('uv_sphere', opts, cell_type, order=order, max_cells=max_cells, &
                      stat=stat, errmsg=errmsg)
    end function

    !> An icosahedron-based sphere surface: 20 * 4**subdivisions triangles,
    !> oriented outward.
    function mio_icosa_sphere(radius, subdivisions, center, order, max_cells, stat, errmsg) &
        result(out)
        real(real64), intent(in), optional :: radius, center(3)
        integer, intent(in), optional :: subdivisions
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        opts%counts(1) = 2
        if (present(radius)) opts%radius = real(radius, c_double)
        if (present(subdivisions)) opts%counts(1) = int(subdivisions, c_int64_t)
        if (present(center)) opts%origin = real(center, c_double)
        out = gen_run('icosa_sphere', opts, order=order, max_cells=max_cells, &
                      stat=stat, errmsg=errmsg)
    end function

    !> A solid cylinder along z of `nxy` cells across and `nz` layers:
    !> 'hexahedron' (default) or 'wedge' cells; side regions side, top, bottom.
    function mio_cylinder(radius, height, nxy, nz, center, cell_type, order, max_cells, &
                          stat, errmsg) result(out)
        real(real64), intent(in), optional :: radius, height, center(3)
        integer, intent(in), optional :: nxy, nz
        character(*), intent(in), optional :: cell_type
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        if (present(radius)) opts%radius = real(radius, c_double)
        if (present(height)) opts%height = real(height, c_double)
        if (present(nxy)) opts%counts(1) = int(nxy, c_int64_t)
        if (present(nz)) opts%counts(2) = int(nz, c_int64_t)
        if (present(center)) opts%origin = real(center, c_double)
        out = gen_run('cylinder', opts, cell_type, order=order, max_cells=max_cells, &
                      stat=stat, errmsg=errmsg)
    end function

    !> The reference triangle (0,0) (1,0) (0,1) scaled by `lengths` (default 1, 1)
    !> and refined into `n` * `n` triangles at z = origin(3). Side regions ymin,
    !> xmin, hypotenuse.
    function mio_reference_triangle(n, lengths, origin, order, max_cells, stat, errmsg) &
        result(out)
        integer, intent(in), optional :: n
        real(real64), intent(in), optional :: lengths(2), origin(3)
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        opts%counts(1) = 4
        if (present(n)) opts%counts(1) = int(n, c_int64_t)
        if (present(lengths)) opts%lengths(1:2) = real(lengths, c_double)
        if (present(origin)) opts%origin = real(origin, c_double)
        out = gen_run('reference_triangle', opts, order=order, max_cells=max_cells, &
                      stat=stat, errmsg=errmsg)
    end function

    !> A regular polygon of `sides` sides fanned from its centre (`rings` = 1) or
    !> refined into `rings` rings of triangles; side region rim.
    function mio_ngon(sides, radius, rings, center, order, max_cells, stat, errmsg) &
        result(out)
        integer, intent(in), optional :: sides, rings
        real(real64), intent(in), optional :: radius, center(3)
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        opts%counts(1) = 6
        opts%counts(2) = 1
        if (present(sides)) opts%counts(1) = int(sides, c_int64_t)
        if (present(rings)) opts%counts(2) = int(rings, c_int64_t)
        if (present(radius)) opts%radius = real(radius, c_double)
        if (present(center)) opts%origin = real(center, c_double)
        out = gen_run('ngon', opts, order=order, max_cells=max_cells, stat=stat, errmsg=errmsg)
    end function

    !> An annulus (or a sector of `angle` radians) of `ntheta` x `nr` quads
    !> (default) or triangles. Side regions inner, outer and, for a sector, start
    !> and end.
    function mio_annulus(inner_radius, outer_radius, ntheta, nr, angle, center, &
                         cell_type, order, max_cells, stat, errmsg) result(out)
        real(real64), intent(in), optional :: inner_radius, outer_radius, angle, center(3)
        integer, intent(in), optional :: ntheta, nr
        character(*), intent(in), optional :: cell_type
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        opts%counts(1) = 16
        opts%counts(2) = 4
        if (present(inner_radius)) opts%radius2 = real(inner_radius, c_double)
        if (present(outer_radius)) opts%radius = real(outer_radius, c_double)
        if (present(ntheta)) opts%counts(1) = int(ntheta, c_int64_t)
        if (present(nr)) opts%counts(2) = int(nr, c_int64_t)
        if (present(angle)) opts%angle = real(angle, c_double)
        if (present(center)) opts%origin = real(center, c_double)
        out = gen_run('annulus', opts, cell_type, order=order, max_cells=max_cells, &
                      stat=stat, errmsg=errmsg)
    end function

    !> The square of side `size` (default 2) around `center` minus its fourth
    !> quadrant (x > 0, y < 0); `n` cells along half a side. Side regions xmin,
    !> xmax, ymin, ymax, notch_x, notch_y.
    function mio_lshape(n, size, center, cell_type, order, max_cells, stat, errmsg) &
        result(out)
        integer, intent(in), optional :: n
        real(real64), intent(in), optional :: size, center(3)
        character(*), intent(in), optional :: cell_type
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        opts%counts(1) = 4
        opts%lengths(1) = 2.0_c_double
        if (present(n)) opts%counts(1) = int(n, c_int64_t)
        if (present(size)) opts%lengths(1) = real(size, c_double)
        if (present(center)) opts%origin = real(center, c_double)
        out = gen_run('lshape', opts, cell_type, order=order, max_cells=max_cells, &
                      stat=stat, errmsg=errmsg)
    end function

    !> A square plate of side `size` around a circular hole of `radius`, in four
    !> quad blocks of `nt` x `nr` cells. Side regions hole, xmin, xmax, ymin, ymax.
    function mio_plate_with_hole(radius, size, nt, nr, center, cell_type, order, &
                                 max_cells, stat, errmsg) result(out)
        real(real64), intent(in), optional :: radius, size, center(3)
        integer, intent(in), optional :: nt, nr
        character(*), intent(in), optional :: cell_type
        integer, intent(in), optional :: order
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        opts%counts(1) = 8
        opts%counts(2) = 4
        opts%radius = 0.25_c_double
        if (present(radius)) opts%radius = real(radius, c_double)
        if (present(size)) opts%lengths(1) = real(size, c_double)
        if (present(nt)) opts%counts(1) = int(nt, c_int64_t)
        if (present(nr)) opts%counts(2) = int(nr, c_int64_t)
        if (present(center)) opts%origin = real(center, c_double)
        out = gen_run('plate_with_hole', opts, cell_type, order=order, max_cells=max_cells, &
                      stat=stat, errmsg=errmsg)
    end function

    !> A tiling of `dims` = (nx, ny) pointy-top hexagons of side `edge` as one
    !> ragged polygon block; side region boundary.
    function mio_honeycomb(dims, edge, origin, max_cells, stat, errmsg) result(out)
        integer, intent(in), optional :: dims(2)
        real(real64), intent(in), optional :: edge, origin(3)
        integer(int64), intent(in), optional :: max_cells
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_gen_opts_t) :: opts
        opts%counts(1:2) = 4
        if (present(dims)) opts%counts(1:2) = int(dims, c_int64_t)
        if (present(edge)) opts%lengths(1) = real(edge, c_double)
        if (present(origin)) opts%origin = real(origin, c_double)
        out = gen_run('honeycomb', opts, max_cells=max_cells, stat=stat, errmsg=errmsg)
    end function

    !> One cell of each type named in `types` (comma-separated meshio names; every
    !> supported type when omitted) in its reference coordinates, side by side.
    !> `pairs` gives two cells sharing a facet (linear types only).
    function mio_reference_cells(types, pairs, stat, errmsg) result(out)
        character(*), intent(in), optional :: types
        logical, intent(in), optional :: pairs
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int32_t) :: cpairs
        cpairs = 0
        if (present(pairs)) then
            if (pairs) cpairs = 1
        end if
        if (present(types)) then
            out%handle = c_mio_reference_cells(trim(types)//c_null_char, cpairs)
        else
            out%handle = c_mio_reference_cells(c_null_char, cpairs)
        end if
        if (.not. c_associated(out%handle)) then
            call handle_failure('reference_cells', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Build a regular grid around the mesh. Exactly one of `resolution` and
    !> `cell_size` must be given. `fill` is 'all' (the whole bounding box),
    !> 'surface' (only cells a triangle passes through) or 'inside' (only cells
    !> whose centre is inside the surface). The lattice geometry comes back
    !> through the optional `dims`/`origin`/`spacing`/`num_occupied` arguments.
    function mesh_voxelize(self, resolution, cell_size, bounds, padding, &
                           padding_relative, fill, sign, watertight_check, &
                           attach_occupancy, max_cells, &
                           dims, origin, spacing, num_occupied, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer, intent(in), optional :: resolution(3)
        real(real64), intent(in), optional :: cell_size, bounds(6)
        real(real64), intent(in), optional :: padding, padding_relative
        character(*), intent(in), optional :: fill, sign, watertight_check
        logical, intent(in), optional :: attach_occupancy
        integer(int64), intent(in), optional :: max_cells
        integer(int64), intent(out), optional :: dims(3), num_occupied
        real(real64), intent(out), optional :: origin(3), spacing(3)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_voxel_opts_t) :: opts
        ! The buffers the option pointers reference must outlive the call, so
        ! they are locals here rather than temporaries.
        integer(c_int64_t), target :: res_buf(3)
        real(c_double), target :: bounds_buf(6)
        integer(c_int64_t) :: cdims(3), coccupied
        real(c_double) :: corigin(3), cspacing(3)
        integer(c_int32_t) :: code

        call c_mio_voxel_opts_init(opts)
        if (present(resolution)) then
            res_buf = int(resolution, c_int64_t)
            opts%resolution = c_loc(res_buf(1))
        end if
        if (present(cell_size)) opts%cell_size = real(cell_size, c_double)
        if (present(bounds)) then
            bounds_buf = real(bounds, c_double)
            opts%bounds = c_loc(bounds_buf(1))
        end if
        if (present(padding)) opts%padding = real(padding, c_double)
        if (present(padding_relative)) opts%padding_relative = real(padding_relative, c_double)
        if (present(attach_occupancy)) then
            if (attach_occupancy) opts%attach_occupancy = 1
        end if
        if (present(max_cells)) opts%max_cells = int(max_cells, c_int64_t)
        if (present(fill)) then
            code = voxel_fill_code(fill)
            if (code < 0) then
                call handle_failure('voxelize', 'unknown fill '//trim(fill), stat, errmsg)
                return
            end if
            opts%fill = code
        end if
        if (present(sign)) then
            code = sdf_sign_code(sign)
            if (code < 0) then
                call handle_failure('voxelize', 'unknown sign '//trim(sign), stat, errmsg)
                return
            end if
            opts%sign = code
        end if
        if (present(watertight_check)) then
            select case (trim(watertight_check))
            case ('off'); opts%watertight_check = 0
            case ('', 'warn'); opts%watertight_check = 1
            case ('error'); opts%watertight_check = 2
            case default
                call handle_failure('voxelize', &
                                    'unknown watertight check '//trim(watertight_check), &
                                    stat, errmsg)
                return
            end select
        end if

        out%handle = c_mio_voxelize(self%handle, opts, cdims, corigin, cspacing, coccupied)
        if (.not. c_associated(out%handle)) then
            call handle_failure('voxelize', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(dims)) dims = int(cdims, int64)
        if (present(origin)) origin = real(corigin, real64)
        if (present(spacing)) spacing = real(cspacing, real64)
        if (present(num_occupied)) num_occupied = int(coccupied, int64)
        call clear_status(stat, errmsg)
    end function

    !> Build an octree (a quadtree for a planar `line` skin) around the mesh:
    !> leaves are fine only where the skin is. `root_cells` (default 1,1,1),
    !> `min_level` (default 0) and `max_level` (default 4) size the tree;
    !> `balance` is 'none', 'face' (default) or 'full'; `fill` is 'all',
    !> 'surface' or 'inside' as for voxelize. `dimension` is 2, 3 or 0 (decide
    !> from the cell types). The result carries `octree:level` and `octree:cut`
    !> cell data and `octree:hanging` point data. The tree's summary comes back
    !> through the optional `num_leaves`/`num_hanging`/`finest_level`/`root_size`
    !> arguments.
    function mesh_octree(self, root_cells, min_level, max_level, bounds, padding, &
                         padding_relative, dimension, sample_distance, balance, fill, &
                         sign, watertight_check, attach_distance, even_odd, max_cells, &
                         transition, snap, surface_level, curvature_threshold, &
                         feature_angle, size_array, plane, &
                         num_leaves, num_hanging, finest_level, root_size, &
                         num_snapped, num_refused, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer, intent(in), optional :: root_cells(3)
        integer(int64), intent(in), optional :: min_level, max_level, dimension, max_cells
        real(real64), intent(in), optional :: bounds(6)
        real(real64), intent(in), optional :: padding, padding_relative, sample_distance
        character(*), intent(in), optional :: balance, fill, sign, watertight_check
        character(*), intent(in), optional :: transition, size_array, plane
        logical, intent(in), optional :: attach_distance, even_odd, snap
        integer(int64), intent(in), optional :: surface_level
        real(real64), intent(in), optional :: curvature_threshold, feature_angle
        integer(int64), intent(out), optional :: num_leaves, num_hanging, finest_level
        integer(int64), intent(out), optional :: num_snapped, num_refused
        real(real64), intent(out), optional :: root_size
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_octree_opts_t) :: opts
        type(mio_octree_report_t) :: report
        integer(c_int64_t), target :: roots_buf(3)
        real(c_double), target :: bounds_buf(6)
        character(kind=c_char), allocatable, target :: size_buf(:)
        integer(c_int32_t) :: code
        integer :: i

        call c_mio_octree_opts_init(opts)
        if (present(transition)) then
            select case (trim(transition))
            case ('', 'hanging'); opts%transition = 0
            case ('polyhedra'); opts%transition = 1
            case ('templates'); opts%transition = 2
            case default
                call handle_failure('octree', 'unknown transition '//trim(transition), stat, errmsg)
                return
            end select
        end if
        if (present(plane)) then
            select case (trim(plane))
            case ('', 'auto'); opts%plane = 0
            case ('xy'); opts%plane = 1
            case ('xz'); opts%plane = 2
            case ('yz'); opts%plane = 3
            case default
                call handle_failure('octree', 'unknown plane '//trim(plane), stat, errmsg)
                return
            end select
        end if
        if (present(snap)) then
            if (snap) opts%snap = 1
        end if
        if (present(surface_level)) opts%surface_level = int(surface_level, c_int64_t)
        if (present(curvature_threshold)) &
            opts%curvature_threshold = real(curvature_threshold, c_double)
        if (present(feature_angle)) opts%feature_angle = real(feature_angle, c_double)
        if (present(size_array)) then
            if (len_trim(size_array) > 0) then
                allocate (size_buf(len_trim(size_array) + 1))
                do i = 1, len_trim(size_array)
                    size_buf(i) = size_array(i:i)
                end do
                size_buf(len_trim(size_array) + 1) = c_null_char
                opts%size_array = c_loc(size_buf(1))
            end if
        end if
        if (present(root_cells)) then
            roots_buf = int(root_cells, c_int64_t)
            opts%root_cells = c_loc(roots_buf(1))
        end if
        if (present(bounds)) then
            bounds_buf = real(bounds, c_double)
            opts%bounds = c_loc(bounds_buf(1))
        end if
        if (present(min_level)) opts%min_level = int(min_level, c_int64_t)
        if (present(max_level)) opts%max_level = int(max_level, c_int64_t)
        if (present(dimension)) opts%dimension = int(dimension, c_int64_t)
        if (present(max_cells)) opts%max_cells = int(max_cells, c_int64_t)
        if (present(padding)) opts%padding = real(padding, c_double)
        if (present(padding_relative)) opts%padding_relative = real(padding_relative, c_double)
        if (present(sample_distance)) opts%sample_distance = real(sample_distance, c_double)
        if (present(attach_distance)) then
            if (attach_distance) opts%attach_distance = 1
        end if
        if (present(even_odd)) then
            if (even_odd) opts%even_odd = 1
        end if
        if (present(balance)) then
            select case (trim(balance))
            case ('none'); opts%balance = 0
            case ('', 'face'); opts%balance = 1
            case ('full'); opts%balance = 2
            case default
                call handle_failure('octree', 'unknown balance '//trim(balance), stat, errmsg)
                return
            end select
        end if
        if (present(fill)) then
            code = voxel_fill_code(fill)
            if (code < 0) then
                call handle_failure('octree', 'unknown fill '//trim(fill), stat, errmsg)
                return
            end if
            opts%fill = code
        end if
        if (present(sign)) then
            code = sdf_sign_code(sign)
            if (code < 0) then
                call handle_failure('octree', 'unknown sign '//trim(sign), stat, errmsg)
                return
            end if
            opts%sign = code
        end if
        if (present(watertight_check)) then
            select case (trim(watertight_check))
            case ('off'); opts%watertight_check = 0
            case ('', 'warn'); opts%watertight_check = 1
            case ('error'); opts%watertight_check = 2
            case default
                call handle_failure('octree', &
                                    'unknown watertight check '//trim(watertight_check), &
                                    stat, errmsg)
                return
            end select
        end if

        out%handle = c_mio_octree(self%handle, opts, report)
        if (.not. c_associated(out%handle)) then
            call handle_failure('octree', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(num_leaves)) num_leaves = int(report%num_leaves, int64)
        if (present(num_hanging)) num_hanging = int(report%num_hanging, int64)
        if (present(finest_level)) finest_level = int(report%finest_level, int64)
        if (present(root_size)) root_size = real(report%root_size, real64)
        if (present(num_snapped)) num_snapped = int(report%num_snapped, int64)
        if (present(num_refused)) num_refused = int(report%num_refused, int64)
        call clear_status(stat, errmsg)
    end function

    !> Mesh the inside of a skin on a balanced octree (see
    !> doc/octree_mesh.md): conforming graded tetrahedra by default, or the
    !> castellated leaves with surrogate-boundary regions and `sdf:distance` /
    !> `sdf:closest_point`. Every option is optional and takes its default
    !> from `mio_octree_mesh_opts_init`, i.e. from the C side. The string
    !> options are `output` ('tet' default, 'castellated' or 'hex'),
    !> `conformity` ('stuffing' default or 'cleaving'), `skin_grading`
    !> ('graded' default or 'uniform') and `fill` ('inside' default, 'all' or
    !> 'cover'); an unknown one fails by name. The whole summary comes back
    !> through the optional `report` (a `mio_octree_mesh_report`), and the
    !> run's warnings, one per line, through the optional `warnings`.
    function mesh_octree_mesh(self, output, min_size, max_size, max_level, bounds, padding, &
                              padding_relative, curvature_angle, feature_angle, &
                              corner_angle, feature_size, gap_cells, gap_outside, &
                              size_array, surface_region, grading_layers, skin_grading, &
                              surface_size, winding_threshold, winding_beta, conformity, &
                              alpha_long, alpha_short, preserve_features, fit_iterations, &
                              dihedral_floor_guard, fill, keep_outside, attach_distance, &
                              attach_quality, pillow, untangle_iterations, max_cells, &
                              max_tets, report, warnings, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in), optional :: output, skin_grading, conformity, fill
        character(*), intent(in), optional :: size_array, surface_region
        real(real64), intent(in), optional :: min_size, max_size, bounds(6)
        real(real64), intent(in), optional :: padding, padding_relative
        real(real64), intent(in), optional :: curvature_angle, feature_angle, corner_angle
        real(real64), intent(in), optional :: feature_size, gap_cells, surface_size
        real(real64), intent(in), optional :: winding_threshold, winding_beta
        real(real64), intent(in), optional :: alpha_long, alpha_short
        integer(int64), intent(in), optional :: max_level, grading_layers, fit_iterations
        integer(int64), intent(in), optional :: untangle_iterations, max_cells, max_tets
        logical, intent(in), optional :: gap_outside, preserve_features
        logical, intent(in), optional :: dihedral_floor_guard, keep_outside
        logical, intent(in), optional :: attach_distance, attach_quality, pillow
        type(mio_octree_mesh_report), intent(out), optional :: report
        character(:), allocatable, intent(out), optional :: warnings
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_octree_mesh_opts_t) :: opts
        type(mio_octree_mesh_report) :: crep
        ! The buffers the option pointers reference must outlive the call.
        real(c_double), target :: bounds_buf(6)
        character(kind=c_char), allocatable, target :: size_buf(:), region_buf(:)
        character(kind=c_char), allocatable :: wbuf(:)
        integer(c_int64_t) :: n
        integer :: i

        call check_octree_mesh_layout()
        if (present(warnings)) warnings = ''
        call c_mio_octree_mesh_opts_init(opts)
        if (present(output)) then
            select case (trim(output))
            case ('castellated'); opts%output = 0
            case ('', 'tet'); opts%output = 1
            case ('hex'); opts%output = 2
            case default
                call handle_failure('octree_mesh', 'unknown output '//trim(output), stat, errmsg)
                return
            end select
        end if
        if (present(conformity)) then
            select case (trim(conformity))
            case ('', 'stuffing'); opts%conformity = 0
            case ('cleaving'); opts%conformity = 1
            case default
                call handle_failure('octree_mesh', 'unknown conformity '//trim(conformity), &
                                    stat, errmsg)
                return
            end select
        end if
        if (present(skin_grading)) then
            select case (trim(skin_grading))
            case ('', 'graded'); opts%skin_grading = 0
            case ('uniform'); opts%skin_grading = 1
            case default
                call handle_failure('octree_mesh', 'unknown skin grading '//trim(skin_grading), &
                                    stat, errmsg)
                return
            end select
        end if
        if (present(fill)) then
            select case (trim(fill))
            case ('all'); opts%fill = 0
            case ('', 'inside'); opts%fill = 1
            case ('cover'); opts%fill = 2
            case default
                call handle_failure('octree_mesh', 'unknown fill '//trim(fill), stat, errmsg)
                return
            end select
        end if
        if (present(bounds)) then
            bounds_buf = real(bounds, c_double)
            opts%bounds = c_loc(bounds_buf(1))
        end if
        if (present(size_array)) then
            if (len_trim(size_array) > 0) then
                allocate (size_buf(len_trim(size_array) + 1))
                do i = 1, len_trim(size_array)
                    size_buf(i) = size_array(i:i)
                end do
                size_buf(len_trim(size_array) + 1) = c_null_char
                opts%size_array = c_loc(size_buf(1))
            end if
        end if
        if (present(surface_region)) then
            if (len_trim(surface_region) > 0) then
                allocate (region_buf(len_trim(surface_region) + 1))
                do i = 1, len_trim(surface_region)
                    region_buf(i) = surface_region(i:i)
                end do
                region_buf(len_trim(surface_region) + 1) = c_null_char
                opts%surface_region = c_loc(region_buf(1))
            end if
        end if
        if (present(min_size)) opts%min_size = real(min_size, c_double)
        if (present(max_size)) opts%max_size = real(max_size, c_double)
        if (present(padding)) opts%padding = real(padding, c_double)
        if (present(padding_relative)) opts%padding_relative = real(padding_relative, c_double)
        if (present(curvature_angle)) opts%curvature_angle = real(curvature_angle, c_double)
        if (present(feature_angle)) opts%feature_angle = real(feature_angle, c_double)
        if (present(corner_angle)) opts%corner_angle = real(corner_angle, c_double)
        if (present(feature_size)) opts%feature_size = real(feature_size, c_double)
        if (present(gap_cells)) opts%gap_cells = real(gap_cells, c_double)
        if (present(surface_size)) opts%surface_size = real(surface_size, c_double)
        if (present(winding_threshold)) &
            opts%winding_threshold = real(winding_threshold, c_double)
        if (present(winding_beta)) opts%winding_beta = real(winding_beta, c_double)
        if (present(alpha_long)) opts%alpha_long = real(alpha_long, c_double)
        if (present(alpha_short)) opts%alpha_short = real(alpha_short, c_double)
        if (present(max_level)) opts%max_level = int(max_level, c_int64_t)
        if (present(grading_layers)) opts%grading_layers = int(grading_layers, c_int64_t)
        if (present(fit_iterations)) opts%fit_iterations = int(fit_iterations, c_int64_t)
        if (present(untangle_iterations)) &
            opts%untangle_iterations = int(untangle_iterations, c_int64_t)
        if (present(max_cells)) opts%max_cells = int(max_cells, c_int64_t)
        if (present(max_tets)) opts%max_tets = int(max_tets, c_int64_t)
        if (present(gap_outside)) opts%gap_outside = merge(1_c_int32_t, 0_c_int32_t, gap_outside)
        if (present(preserve_features)) &
            opts%preserve_features = merge(1_c_int32_t, 0_c_int32_t, preserve_features)
        if (present(dihedral_floor_guard)) &
            opts%dihedral_floor_guard = merge(1_c_int32_t, 0_c_int32_t, dihedral_floor_guard)
        if (present(keep_outside)) &
            opts%keep_outside = merge(1_c_int32_t, 0_c_int32_t, keep_outside)
        if (present(attach_distance)) &
            opts%attach_distance = merge(1_c_int32_t, 0_c_int32_t, attach_distance)
        if (present(attach_quality)) &
            opts%attach_quality = merge(1_c_int32_t, 0_c_int32_t, attach_quality)
        if (present(pillow)) opts%pillow = merge(1_c_int32_t, 0_c_int32_t, pillow)

        out%handle = c_mio_octree_mesh(self%handle, opts, crep)
        if (.not. c_associated(out%handle)) then
            call handle_failure('octree_mesh', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(report)) report = crep
        if (present(warnings)) then
            ! String rule 5: query the length with an empty buffer, then copy.
            allocate (wbuf(1))
            n = c_mio_octree_mesh_warnings(wbuf, 0_c_int64_t)
            if (n > 0 .and. n < int(huge(0) - 1, c_int64_t)) then
                deallocate (wbuf)
                allocate (wbuf(int(n) + 1))
                n = c_mio_octree_mesh_warnings(wbuf, n + 1)
                if (n > 0) warnings = from_c_buf(wbuf, int(n))
            end if
        end if
        call clear_status(stat, errmsg)
    end function

    !> The layouts of `mio_octree_mesh_opts_t` and `mio_octree_mesh_report`
    !> against the C structs (264 and 392 bytes), checked once per process.
    subroutine check_octree_mesh_layout()
        logical, save :: checked = .false.
        type(mio_octree_mesh_opts_t) :: opts
        type(mio_octree_mesh_report) :: report
        if (checked) return
        checked = .true.
        if (c_sizeof(opts) /= 264_c_size_t .or. c_sizeof(report) /= 392_c_size_t) then
            write (error_unit, '(a,2(i0,1x),a)') &
                'meshio++: octree_mesh layout mismatch (', c_sizeof(opts), c_sizeof(report), &
                'bytes)'
            error stop 1
        end if
    end subroutine

    !> Crop to the cells whose value in a scalar cell_data array satisfies a
    !> comparison. `compare` is one of '<', '<=', '>', '>=', '==', '!=' -- the
    !> SAME vocabulary `refine`'s `where_op` uses, evaluated by the same C++
    !> function, so the two cannot drift on the boundary cases. A NON-FINITE
    !> cell value never matches, whatever the comparison.
    !>
    !> There is deliberately no `mode`: `crop_bbox`/`crop_plane` test points and
    !> then need an all/any rule, whereas a cell_data predicate is already one
    !> value per cell and has nothing to reduce.
    !>
    !> Inside/outside a surface composes rather than being its own procedure:
    !> `distance_to_surface(..., location='center')` then this on 'sdf:distance'.
    function mesh_crop_predicate(self, array, compare, value, record_ids, stat, errmsg) &
            result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: array
        character(*), intent(in), optional :: compare
        real(real64), intent(in), optional :: value
        logical, intent(in), optional :: record_ids
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: ccmp, crec
        real(c_double) :: cval
        ccmp = 0
        if (present(compare)) then
            ccmp = refine_compare_code(compare)
            if (ccmp < 0) then
                call handle_failure('crop_predicate', 'unknown comparison '//trim(compare), &
                                    stat, errmsg)
                return
            end if
        end if
        cval = 0.0_c_double
        if (present(value)) cval = real(value, c_double)
        crec = 0
        if (present(record_ids)) then
            if (record_ids) crec = 1
        end if
        out%handle = c_mio_crop_predicate(self%handle, c_str(array), ccmp, cval, crec)
        if (.not. c_associated(out%handle)) then
            call handle_failure('crop_predicate', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Generate a grid over this surface and fill it with signed distances.
    !>
    !> `structure` is 'voxel' (a dense lattice, the default) or 'octree' (refined
    !> near the surface, and therefore 1-irregular -- it has hanging nodes; see
    !> `refine`'s 'balanced' closure). `resolution`/`cell_size` size a voxel grid
    !> and are an ERROR with 'octree', whose finest cell is
    !> root_resolution/2**max_depth and is therefore already determined.
    !>
    !> `spacing` reports the FINEST cell size, `dims` the ROOT cell counts. The
    !> generated mesh also carries the numeric `sdf:*` field_data header; no file
    !> format persists arbitrary field_data, so write it as `.vti`, whose
    !> Origin/Spacing/WholeExtent attributes are the same information.
    function mesh_compute_sdf(self, structure, resolution, cell_size, bounds, padding, &
                              padding_relative, root_resolution, max_depth, band_cells, &
                              record_levels, max_cells, sign, weight, location, band, &
                              watertight_check, dims, origin, spacing, depth, num_banded, &
                              quality, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in), optional :: structure, sign, weight, location
        character(*), intent(in), optional :: watertight_check
        integer, intent(in), optional :: resolution(3)
        real(real64), intent(in), optional :: cell_size, bounds(6)
        real(real64), intent(in), optional :: padding, padding_relative, band_cells, band
        integer(int64), intent(in), optional :: root_resolution, max_depth, max_cells
        logical, intent(in), optional :: record_levels
        integer(int64), intent(out), optional :: dims(3), depth, num_banded
        real(real64), intent(out), optional :: origin(3), spacing(3)
        type(mio_surface_quality), intent(out), optional :: quality
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_compute_sdf_opts_t) :: opts
        ! The buffers the option pointers reference must outlive the call, so
        ! they are locals here rather than temporaries.
        integer(c_int64_t), target :: res_buf(3)
        real(c_double), target :: bounds_buf(6)
        integer(c_int64_t) :: cdims(3), cdepth, cbanded
        real(c_double) :: corigin(3), cspacing(3)
        type(mio_surface_quality) :: q
        integer(c_int32_t) :: code

        call c_mio_compute_sdf_opts_init(opts)
        if (present(structure)) then
            select case (trim(structure))
            case ('', 'voxel'); opts%structure = 0
            case ('octree'); opts%structure = 1
            case default
                call handle_failure('compute_sdf', 'unknown structure '//trim(structure), &
                                    stat, errmsg)
                return
            end select
        end if
        if (present(resolution)) then
            res_buf = int(resolution, c_int64_t)
            opts%resolution = c_loc(res_buf(1))
        end if
        if (present(cell_size)) opts%cell_size = real(cell_size, c_double)
        if (present(bounds)) then
            bounds_buf = real(bounds, c_double)
            opts%bounds = c_loc(bounds_buf(1))
        end if
        if (present(padding)) opts%padding = real(padding, c_double)
        if (present(padding_relative)) opts%padding_relative = real(padding_relative, c_double)
        if (present(root_resolution)) opts%root_resolution = int(root_resolution, c_int64_t)
        if (present(max_depth)) opts%max_depth = int(max_depth, c_int64_t)
        if (present(band_cells)) opts%band_cells = real(band_cells, c_double)
        if (present(record_levels)) then
            opts%record_levels = 0
            if (record_levels) opts%record_levels = 1
        end if
        if (present(max_cells)) opts%max_cells = int(max_cells, c_int64_t)
        if (present(band)) opts%distance%band = real(band, c_double)
        if (present(sign)) then
            code = sdf_sign_code(sign)
            if (code < 0) then
                call handle_failure('compute_sdf', 'unknown sign '//trim(sign), stat, errmsg)
                return
            end if
            opts%distance%sign = code
        end if
        if (present(weight)) then
            select case (trim(weight))
            case ('', 'angle'); opts%distance%weight = 0
            case ('area'); opts%distance%weight = 1
            case default
                call handle_failure('compute_sdf', 'unknown weight '//trim(weight), stat, errmsg)
                return
            end select
        end if
        if (present(location)) then
            select case (trim(location))
            case ('', 'corner', 'point'); opts%distance%location = 0
            case ('center', 'centre', 'cell'); opts%distance%location = 1
            case default
                call handle_failure('compute_sdf', 'unknown location '//trim(location), &
                                    stat, errmsg)
                return
            end select
        end if
        if (present(watertight_check)) then
            select case (trim(watertight_check))
            case ('off'); opts%distance%watertight_check = 0
            case ('', 'warn'); opts%distance%watertight_check = 1
            case ('error'); opts%distance%watertight_check = 2
            case default
                call handle_failure('compute_sdf', &
                                    'unknown watertight check '//trim(watertight_check), &
                                    stat, errmsg)
                return
            end select
        end if

        out%handle = c_mio_compute_sdf(self%handle, opts, cdims, corigin, cspacing, &
                                       cdepth, cbanded, q)
        if (.not. c_associated(out%handle)) then
            call handle_failure('compute_sdf', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(dims)) dims = int(cdims, int64)
        if (present(origin)) origin = real(corigin, real64)
        if (present(spacing)) spacing = real(cspacing, real64)
        if (present(depth)) depth = int(cdepth, int64)
        if (present(num_banded)) num_banded = int(cbanded, int64)
        if (present(quality)) quality = q
        call clear_status(stat, errmsg)
    end function

    !> What is wrong with this surface, in numbers: boundary edges,
    !> non-manifold edges, inconsistently wound pairs and degenerate triangles.
    function mesh_watertight_check(self, stat, errmsg) result(q)
        class(mio_mesh), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_surface_quality) :: q
        integer(c_int) :: s
        s = c_mio_surface_watertight_check(self%handle, q)
        if (s /= 0_c_int) then
            call handle_failure('surface_watertight_check', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Signed distances from `points` (shaped (3, n), the column-major memory
    !> the C ABI reads as n rows of three) to this surface. Negative is inside.
    function mesh_sample_distance(self, points, sign, band, watertight_check, &
                                  stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        real(real64), intent(in) :: points(:, :)
        character(*), intent(in), optional :: sign, watertight_check
        real(real64), intent(in), optional :: band
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: out(:)
        type(mio_sdf_opts_t) :: opts
        real(c_double), allocatable :: buf(:), cpoints(:)
        integer(c_int) :: s
        integer(c_int32_t) :: code
        integer(int64) :: n

        n = size(points, 2)
        allocate (out(max(n, 1_int64)))
        out = 0.0_real64
        if (size(points, 1) /= 3) then
            call handle_failure('sample_distance', 'points must be shaped (3, n)', stat, errmsg)
            return
        end if

        call c_mio_sdf_opts_init(opts)
        if (present(band)) opts%band = real(band, c_double)
        if (present(sign)) then
            code = sdf_sign_code(sign)
            if (code < 0) then
                call handle_failure('sample_distance', 'unknown sign '//trim(sign), stat, errmsg)
                return
            end if
            opts%sign = code
        end if
        if (present(watertight_check)) then
            select case (trim(watertight_check))
            case ('off'); opts%watertight_check = 0
            case ('', 'warn'); opts%watertight_check = 1
            case ('error'); opts%watertight_check = 2
            case default
                call handle_failure('sample_distance', &
                                    'unknown watertight check '//trim(watertight_check), &
                                    stat, errmsg)
                return
            end select
        end if

        allocate (cpoints(max(3_int64*n, 1_int64)), buf(max(n, 1_int64)))
        cpoints = 0.0_c_double
        buf = 0.0_c_double
        if (n > 0) cpoints(1:3*n) = real(reshape(points, [3*n]), c_double)
        s = c_mio_sample_distance(self%handle, cpoints, int(n, c_int64_t), opts, buf)
        if (s /= 0_c_int) then
            call handle_failure('sample_distance', mio_error_message(), stat, errmsg)
            return
        end if
        if (n > 0) out(1:n) = real(buf(1:n), real64)
        call clear_status(stat, errmsg)
    end function

    !> Attach the signed distance from this mesh's points (or cell centres) to
    !> `surface`, as the `sdf:distance` array.
    function mesh_distance_to_surface(self, surface, sign, location, band, &
                                      record_inside, watertight_check, num_banded, &
                                      quality, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        type(mio_mesh), intent(in) :: surface
        character(*), intent(in), optional :: sign, location, watertight_check
        real(real64), intent(in), optional :: band
        logical, intent(in), optional :: record_inside
        integer(int64), intent(out), optional :: num_banded
        type(mio_surface_quality), intent(out), optional :: quality
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_sdf_opts_t) :: opts
        type(mio_surface_quality) :: q
        integer(c_int64_t) :: banded
        integer(c_int32_t) :: code

        call c_mio_sdf_opts_init(opts)
        if (present(band)) opts%band = real(band, c_double)
        if (present(record_inside)) then
            if (record_inside) opts%record_inside = 1
        end if
        if (present(sign)) then
            code = sdf_sign_code(sign)
            if (code < 0) then
                call handle_failure('distance_to_surface', 'unknown sign '//trim(sign), &
                                    stat, errmsg)
                return
            end if
            opts%sign = code
        end if
        if (present(location)) then
            select case (trim(location))
            case ('', 'corner', 'point'); opts%location = 0
            case ('center', 'centre', 'cell'); opts%location = 1
            case default
                call handle_failure('distance_to_surface', &
                                    'unknown location '//trim(location), stat, errmsg)
                return
            end select
        end if
        if (present(watertight_check)) then
            select case (trim(watertight_check))
            case ('off'); opts%watertight_check = 0
            case ('', 'warn'); opts%watertight_check = 1
            case ('error'); opts%watertight_check = 2
            case default
                call handle_failure('distance_to_surface', &
                                    'unknown watertight check '//trim(watertight_check), &
                                    stat, errmsg)
                return
            end select
        end if

        out%handle = c_mio_distance_to_surface(self%handle, surface%handle, opts, banded, q)
        if (.not. c_associated(out%handle)) then
            call handle_failure('distance_to_surface', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(num_banded)) num_banded = int(banded, int64)
        if (present(quality)) quality = q
        call clear_status(stat, errmsg)
    end function

    !> The `mio_refine_closure` code for a closure name, or -1 if unknown.
    function refine_closure_code(name) result(code)
        character(*), intent(in) :: name
        integer(c_int32_t) :: code
        select case (trim(name))
        case ('', 'redgreen', 'red-green', 'green'); code = 0
        case ('propagate', 'red'); code = 1
        case ('balanced', '2:1'); code = 2
        case default; code = -1
        end select
    end function

    !> One-time ABI layout guard for mio_refine_opts_t, the Julia
    !> `_check_abi_layout()` precedent applied to Fortran. The C header pins
    !> `sizeof(mio_refine_opts) == 112` with a `static_assert`, and Julia
    !> checks its own mirror at load time, but Fortran has no compile-time
    !> equivalent -- without this, a field added to one side and not the
    !> other silently corrupts every call rather than failing loudly. Runs
    !> once per process (a SAVE'd flag), not on every refine() call.
    subroutine check_refine_opts_layout()
        logical, save :: checked = .false.
        type(mio_refine_opts_t) :: probe
        if (checked) return
        checked = .true.
        if (c_sizeof(probe) /= 112_c_size_t) then
            write (error_unit, '(a,i0,a)') &
                'meshio++: mio_refine_opts_t layout mismatch (', c_sizeof(probe), ' bytes)'
            error stop 1
        end if
    end subroutine

    function mesh_refine(self, levels, record_parent_ids, point_map, stat, errmsg, &
                         cells, region, where_array, where_op, where_value, closure, &
                         record_levels, record_hierarchy) result(out)
        class(mio_mesh), intent(in) :: self
        integer, intent(in), optional :: levels
        logical, intent(in), optional :: record_parent_ids
        integer(int64), allocatable, intent(out), optional :: point_map(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        !> Global (block-major) 1-BASED indices of the cells to refine. Shifted
        !> to the C API's 0-based numbering inside.
        integer(int64), intent(in), optional :: cells(:)
        !> Name of a region to refine. A cell region selects its own cells; a
        !> point region selects every cell with ANY node in it; a side region is
        !> an error. At most one of cells/region/where_array may be given.
        character(*), intent(in), optional :: region
        !> Name of a scalar cell_data array to threshold, with `where_op` one of
        !> '<', '<=', '>', '>=', '==', '!='. A non-finite value never matches.
        character(*), intent(in), optional :: where_array
        character(*), intent(in), optional :: where_op
        real(real64), intent(in), optional :: where_value
        !> 'redgreen' (default, local, conforming), 'propagate' (conforming but
        !> spreads to the whole edge-connected component) or 'balanced' (keeps
        !> the hanging nodes and only enforces 2:1 balance -- NOT conforming).
        character(*), intent(in), optional :: closure
        !> Attach the refine:level cell_data array.
        logical, intent(in), optional :: record_levels
        !> Attach refine:cell_id/refine:parent_id -- the persistent parent/
        !> child hierarchy a multigrid caller resolves across the sequence of
        !> meshes it keeps. An input already carrying refine:cell_id is
        !> updated whatever this says. Also forces refine:entity to be
        !> attached even when the closure leaves no hanging node.
        logical, intent(in), optional :: record_hierarchy
        type(mio_mesh) :: out
        type(c_ptr) :: res, cdata
        integer(c_int) :: s, dt
        integer(c_int64_t) :: nlen
        integer(c_int64_t), pointer :: fp(:)
        type(mio_refine_opts_t) :: opts
        ! NUL-terminated copies must outlive the call, so they are held here
        ! rather than built inline; c_loc needs them contiguous and TARGET.
        character(kind=c_char, len=STRBUF_LEN), target :: region_buf, array_buf
        integer(c_int64_t), allocatable, target :: cell_ids(:)

        call check_refine_opts_layout()
        call c_mio_refine_opts_init(opts)
        if (present(levels)) opts%levels = int(levels, c_int32_t)
        if (present(record_parent_ids)) then
            if (record_parent_ids) opts%record_parent_ids = 1
        end if
        if (present(record_levels)) then
            if (record_levels) opts%record_levels = 1
        end if
        if (present(record_hierarchy)) then
            if (record_hierarchy) opts%record_hierarchy = 1_c_int64_t
        end if
        if (present(cells)) then
            allocate (cell_ids(max(size(cells), 1)))
            cell_ids = 0_c_int64_t
            if (size(cells) > 0) cell_ids(1:size(cells)) = int(cells, c_int64_t) - 1_c_int64_t
            opts%cells = c_loc(cell_ids(1))
            opts%num_cells = int(size(cells), c_int64_t)
        end if
        if (present(region)) then
            region_buf = trim(region)//c_null_char
            opts%region = c_loc(region_buf(1:1))
        end if
        if (present(where_array)) then
            array_buf = trim(where_array)//c_null_char
            opts%predicate_array = c_loc(array_buf(1:1))
        end if
        if (present(where_value)) opts%predicate_value = real(where_value, c_double)
        if (present(where_op)) then
            opts%predicate_op = refine_compare_code(where_op)
            if (opts%predicate_op < 0) then
                call handle_failure('refine', "unknown comparison '"//trim(where_op)//"'", &
                                    stat, errmsg)
                return
            end if
        end if
        if (present(closure)) then
            opts%closure = refine_closure_code(closure)
            if (opts%closure < 0) then
                call handle_failure('refine', "unknown closure '"//trim(closure)//"'", &
                                    stat, errmsg)
                return
            end if
        end if

        res = c_mio_refine_ex(self%handle, opts)
        if (.not. c_associated(res)) then
            call handle_failure('refine', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(point_map)) then
            s = c_mio_refine_result_point_map(res, cdata, dt, nlen)
            if (s /= 0_c_int) then
                call c_mio_refine_result_free(res)
                call handle_failure('refine', mio_error_message(), stat, errmsg)
                return
            end if
            allocate (point_map(nlen))
            if (nlen > 0) then
                call c_f_pointer(cdata, fp, [nlen])
                point_map = int(fp, int64) + 1_int64  ! 0-based -> 1-based
            end if
        end if
        out%handle = c_mio_refine_result_take_mesh(res)
        call c_mio_refine_result_free(res)
        if (.not. c_associated(out%handle)) then
            call handle_failure('refine', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Tessellate this mesh's curved (quadratic) cells isoparametrically:
    !> triangle6/quad8/quad9 become triangles and tetra10/hexahedron27 tetras
    !> on a reference lattice; every other cell passes through. The output
    !> carries point data `tessellate:source_point` (0-based source point, or
    !> -1) and cell data `tessellate:source_cell` / `tessellate:sub_index`;
    !> `record_stencil` adds `tessellate:stencil` / `tessellate:weights`.
    !> `levels` defaults to 2, `curved` and `fields` to true. See
    !> doc/tessellation.md.
    function mesh_tessellate(self, levels, curved, fields, record_stencil, &
                             num_source_points, num_source_cells, num_curved_source_cells, &
                             num_pass_through_source_cells, num_points, num_cells, &
                             num_boundary_facets, num_interior_facets, &
                             num_facets_with_bad_count, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer(int64), intent(in), optional :: levels
        logical, intent(in), optional :: curved, fields, record_stencil
        integer(int64), intent(out), optional :: num_source_points, num_source_cells
        integer(int64), intent(out), optional :: num_curved_source_cells
        integer(int64), intent(out), optional :: num_pass_through_source_cells
        integer(int64), intent(out), optional :: num_points, num_cells
        integer(int64), intent(out), optional :: num_boundary_facets, num_interior_facets
        integer(int64), intent(out), optional :: num_facets_with_bad_count
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(c_ptr) :: res
        type(mio_tessellate_opts_t) :: opts
        type(mio_tessellate_report_t) :: report

        call c_mio_tessellate_opts_init(opts)
        if (present(levels)) opts%levels = int(levels, c_int64_t)
        if (present(curved)) opts%curved = merge(1_c_int32_t, 0_c_int32_t, curved)
        if (present(fields)) opts%fields = merge(1_c_int32_t, 0_c_int32_t, fields)
        if (present(record_stencil)) &
            opts%record_stencil = merge(1_c_int32_t, 0_c_int32_t, record_stencil)
        res = c_mio_tessellate(self%handle, opts, report)
        if (.not. c_associated(res)) then
            call handle_failure('tessellate', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(num_source_points)) &
            num_source_points = int(report%num_source_points, int64)
        if (present(num_source_cells)) num_source_cells = int(report%num_source_cells, int64)
        if (present(num_curved_source_cells)) &
            num_curved_source_cells = int(report%num_curved_source_cells, int64)
        if (present(num_pass_through_source_cells)) &
            num_pass_through_source_cells = int(report%num_pass_through_source_cells, int64)
        if (present(num_points)) num_points = int(report%num_points, int64)
        if (present(num_cells)) num_cells = int(report%num_cells, int64)
        if (present(num_boundary_facets)) &
            num_boundary_facets = int(report%num_boundary_facets, int64)
        if (present(num_interior_facets)) &
            num_interior_facets = int(report%num_interior_facets, int64)
        if (present(num_facets_with_bad_count)) &
            num_facets_with_bad_count = int(report%num_facets_with_bad_count, int64)
        call clear_status(stat, errmsg)
    end function

    !> Decimate a SURFACE mesh by quadric-error-metric edge collapse — the
    !> resolution-reducing inverse of refine. Exactly one of `ratio`,
    !> `target_faces`, `max_error` must be given. The output is all-triangle;
    !> boundary and feature vertices are pinned by default. The optional
    !> `point_map` receives, 1-based, each input point's surviving output index
    !> (0 when the survivor itself was pruned).
    !> One-time ABI layout guard for mio_decimate_opts_t.
    subroutine check_decimate_opts_layout()
        logical, save :: checked = .false.
        type(mio_decimate_opts_t) :: probe
        if (checked) return
        checked = .true.
        if (c_sizeof(probe) /= 96_c_size_t) then
            write (error_unit, '(a,i0,a)') &
                'meshio++: mio_decimate_opts_t layout mismatch (', c_sizeof(probe), ' bytes)'
            error stop 1
        end if
    end subroutine

    function mesh_decimate(self, ratio, target_faces, max_error, placement, &
                           preserve_boundary, preserve_features, feature_angle, &
                           faces_removed, points_removed, collapses_rejected, &
                           max_error_applied, point_map, stat, errmsg, frozen) result(out)
        class(mio_mesh), intent(in) :: self
        real(real64), intent(in), optional :: ratio
        integer(int64), intent(in), optional :: target_faces
        real(real64), intent(in), optional :: max_error
        character(*), intent(in), optional :: placement
        logical, intent(in), optional :: preserve_boundary
        logical, intent(in), optional :: preserve_features
        real(real64), intent(in), optional :: feature_angle
        integer(int64), intent(out), optional :: faces_removed
        integer(int64), intent(out), optional :: points_removed
        integer(int64), intent(out), optional :: collapses_rejected
        real(real64), intent(out), optional :: max_error_applied
        integer(int64), allocatable, intent(out), optional :: point_map(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        !> 1-based point ids to pin outright. Shifted to the C API's 0-based
        !> numbering inside.
        integer(int64), intent(in), optional :: frozen(:)
        type(mio_mesh) :: out
        type(c_ptr) :: res, cdata
        integer(c_int64_t) :: nlen
        integer(c_int) :: s, dt
        integer(c_int64_t), pointer :: fp(:)
        type(mio_decimate_opts_t) :: opts
        character(kind=c_char, len=STRBUF_LEN), target :: placement_buf
        integer(c_int64_t), allocatable, target :: frozen_ids(:)

        call check_decimate_opts_layout()
        call c_mio_decimate_opts_init(opts)
        if (present(ratio)) opts%target_ratio = real(ratio, c_double)
        if (present(target_faces)) opts%target_faces = int(target_faces, c_int64_t)
        if (present(max_error)) opts%max_error = real(max_error, c_double)
        if (present(placement)) then
            placement_buf = trim(placement)//c_null_char
            opts%placement = c_loc(placement_buf(1:1))
        end if
        if (present(preserve_boundary)) then
            opts%preserve_boundary = 0
            if (preserve_boundary) opts%preserve_boundary = 1
        end if
        if (present(preserve_features)) then
            opts%preserve_features = 0
            if (preserve_features) opts%preserve_features = 1
        end if
        if (present(feature_angle)) opts%feature_angle = real(feature_angle, c_double)
        if (present(frozen)) then
            allocate (frozen_ids(max(size(frozen), 1)))
            frozen_ids = 0_c_int64_t
            if (size(frozen) > 0) frozen_ids(1:size(frozen)) = int(frozen, c_int64_t) - 1_c_int64_t
            opts%frozen = c_loc(frozen_ids(1))
            opts%num_frozen = int(size(frozen), c_int64_t)
        end if

        res = c_mio_decimate_ex(self%handle, opts)
        if (.not. c_associated(res)) then
            call handle_failure('decimate', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(faces_removed)) faces_removed = int(c_mio_decimate_result_faces_removed(res), int64)
        if (present(points_removed)) points_removed = int(c_mio_decimate_result_points_removed(res), int64)
        if (present(collapses_rejected)) then
            collapses_rejected = int(c_mio_decimate_result_collapses_rejected(res), int64)
        end if
        if (present(max_error_applied)) then
            max_error_applied = real(c_mio_decimate_result_max_error_applied(res), real64)
        end if
        if (present(point_map)) then
            s = c_mio_decimate_result_point_map(res, cdata, dt, nlen)
            if (s /= 0_c_int) then
                call c_mio_decimate_result_free(res)
                call handle_failure('decimate', mio_error_message(), stat, errmsg)
                return
            end if
            allocate (point_map(nlen))
            if (nlen > 0) then
                call c_f_pointer(cdata, fp, [nlen])
                point_map = int(fp, int64) + 1_int64  ! 0-based -> 1-based
            end if
        end if
        out%handle = c_mio_decimate_result_take_mesh(res)
        call c_mio_decimate_result_free(res)
        if (.not. c_associated(out%handle)) then
            call handle_failure('decimate', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Split the mesh into pieces by "type", "component", or "region"/"tag".
    !> Returns an array of meshes (one per piece); the optional `keys` array
    !> receives each piece's key. For "region"/"tag", `tag_name` selects the
    !> integer cell_data to split on (default: first integer cell_data).
    function mesh_split(self, by, tag_name, keys, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: by
        character(*), intent(in), optional :: tag_name
        character(len=STRBUF_LEN), allocatable, intent(out), optional :: keys(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh), allocatable :: out(:)
        type(c_ptr) :: res
        integer(c_int64_t) :: count, i, n
        character(kind=c_char) :: buf(STRBUF_LEN)
        character(len=:), allocatable :: tag
        tag = ''
        if (present(tag_name)) tag = tag_name
        res = c_mio_split(self%handle, c_str(by), c_str(tag))
        if (.not. c_associated(res)) then
            call handle_failure('split', mio_error_message(), stat, errmsg)
            allocate (out(0))
            return
        end if
        count = c_mio_split_result_count(res)
        allocate (out(count))
        if (present(keys)) allocate (keys(count))
        do i = 1, count
            out(i)%handle = c_mio_split_result_take_mesh(res, i - 1_c_int64_t)
            if (present(keys)) then
                n = c_mio_split_result_key(res, i - 1_c_int64_t, buf, int(STRBUF_LEN, c_int64_t))
                keys(i) = ''
                if (n > 0) keys(i) = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
            end if
        end do
        call c_mio_split_result_free(res)
        call clear_status(stat, errmsg)
    end function

    !> Decompose the mesh into exactly `nparts` balanced pieces (the
    !> count-driven complement to `split`). method: "sfc" (Hilbert curve cut,
    !> always available), "kahip" (needs a KaHIP-enabled build; fails by name
    !> otherwise), or "auto" (default: kahip when built, else sfc). Every piece
    !> keeps the input's cell-block structure 1:1, so the pieces recombine into
    !> the input. `weights_key` names a scalar cell_data array of per-cell
    !> weights. `ghost_layers` > 0 grows each piece by that many shared-node
    !> BFS layers of other parts' cells (a halo), tagged `partition:ghost`.
    function mesh_partition(self, nparts, method, imbalance, mode, seed, record_ids, &
                            ghost_layers, weights_key, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: nparts
        character(*), intent(in), optional :: method
        real(real64), intent(in), optional :: imbalance
        character(*), intent(in), optional :: mode
        integer, intent(in), optional :: seed
        logical, intent(in), optional :: record_ids
        integer, intent(in), optional :: ghost_layers
        character(*), intent(in), optional :: weights_key
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh), allocatable :: out(:)
        type(c_ptr) :: res
        integer(c_int64_t) :: count, i
        character(len=:), allocatable :: cmethod, cmode, cweights
        real(c_double) :: cimb
        integer(c_int) :: cseed, crec, cghost
        cmethod = 'auto'
        if (present(method)) cmethod = method
        cimb = 0.03_c_double
        if (present(imbalance)) cimb = real(imbalance, c_double)
        cmode = 'eco'
        if (present(mode)) cmode = mode
        cseed = 0_c_int
        if (present(seed)) cseed = int(seed, c_int)
        crec = 0_c_int
        if (present(record_ids)) then
            if (record_ids) crec = 1_c_int
        end if
        cghost = 0_c_int
        if (present(ghost_layers)) cghost = int(ghost_layers, c_int)
        cweights = ''
        if (present(weights_key)) cweights = weights_key
        res = c_mio_partition(self%handle, int(nparts, c_int), c_str(cmethod), cimb, &
                              c_str(cmode), cseed, crec, cghost, c_str(cweights))
        if (.not. c_associated(res)) then
            call handle_failure('partition', mio_error_message(), stat, errmsg)
            allocate (out(0))
            return
        end if
        count = c_mio_partition_result_num_pieces(res)
        allocate (out(count))
        do i = 1, count
            out(i)%handle = c_mio_partition_result_take_mesh(res, i - 1_c_int64_t)
        end do
        call c_mio_partition_result_free(res)
        call clear_status(stat, errmsg)
    end function

    !> The per-cell part assignment only: a flat block-major Int64 array (the
    !> mesh's cell blocks concatenated in order), values in [0, nparts). The
    !> values are part ids, not indices -- no 1-based shift applies.
    function mesh_partition_labels(self, nparts, method, imbalance, mode, seed, &
                                   weights_key, stat, errmsg) result(labels)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: nparts
        character(*), intent(in), optional :: method
        real(real64), intent(in), optional :: imbalance
        character(*), intent(in), optional :: mode
        integer, intent(in), optional :: seed
        character(*), intent(in), optional :: weights_key
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int64), allocatable :: labels(:)
        character(len=:), allocatable :: cmethod, cmode, cweights
        real(c_double) :: cimb
        integer(c_int) :: cseed, s
        integer(c_int64_t) :: total, b
        integer(c_int64_t), allocatable :: buf(:)
        cmethod = 'auto'
        if (present(method)) cmethod = method
        cimb = 0.03_c_double
        if (present(imbalance)) cimb = real(imbalance, c_double)
        cmode = 'eco'
        if (present(mode)) cmode = mode
        cseed = 0_c_int
        if (present(seed)) cseed = int(seed, c_int)
        cweights = ''
        if (present(weights_key)) cweights = weights_key
        total = 0_c_int64_t
        do b = 1, self%num_cell_blocks()
            total = total + self%cell_block_num_cells(int(b))
        end do
        allocate (buf(max(total, 1_c_int64_t)))
        s = c_mio_partition_labels(self%handle, int(nparts, c_int), c_str(cmethod), cimb, &
                                   c_str(cmode), cseed, c_str(cweights), buf, total)
        if (s /= 0_c_int) then
            call handle_failure('partition_labels', mio_error_message(), stat, errmsg)
            allocate (labels(0))
            return
        end if
        allocate (labels(total))
        if (total > 0) labels = int(buf(1:total), int64)
        call clear_status(stat, errmsg)
    end function

    !> Geometric statistics of the mesh (bbox / centroid / area / volume /
    !> inverted count). Per-cell-type counts are not carried across the C ABI.
    function mesh_stats(self, stat, errmsg) result(rep)
        class(mio_mesh), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_stats_report) :: rep
        integer(c_int) :: s
        s = c_mio_stats(self%handle, rep)
        if (s /= 0_c_int) then
            call handle_failure('stats', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    ! ------------------------------------------------------------------
    ! Data operations
    !
    ! These act on the mesh's point/cell/field data arrays; the geometry is
    ! never modified. Each returns a NEW mesh which the caller must free.
    ! `location` is one of MIO_DATA_POINT / _CELL / _FIELD.
    ! ------------------------------------------------------------------

    !> Convert point/cell sets into scalar labels (labels remain zero-based).
    function mesh_sets_to_data(self, location, data_name, join_char, order, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer(c_int), intent(in) :: location
        character(*), intent(in), optional :: data_name, join_char, order(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        character(kind=c_char), allocatable, target :: cname(:), cjoin(:), storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: pname, pjoin, arr
        integer(c_int64_t) :: count
        pname = c_null_ptr
        pjoin = c_null_ptr
        arr = c_null_ptr
        count = 0_c_int64_t
        if (present(data_name)) then
            cname = c_str(data_name)
            pname = c_loc(cname(1))
        end if
        if (present(join_char)) then
            cjoin = c_str(join_char)
            pjoin = c_loc(cjoin(1))
        end if
        if (present(order)) call c_str_array(order, storage, cptrs, arr, count)
        out%handle = c_mio_sets_to_data(self%handle, location, pname, pjoin, arr, count)
        if (.not. c_associated(out%handle)) then
            call handle_failure('sets_to_data', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Convert one scalar integer point/cell field to sets on a new mesh.
    function mesh_data_to_sets(self, location, key, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer(c_int), intent(in) :: location
        character(*), intent(in) :: key
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        out%handle = c_mio_data_to_sets(self%handle, location, c_str(key))
        if (.not. c_associated(out%handle)) then
            call handle_failure('data_to_sets', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Drop the named data arrays at `location`.
    function mesh_data_drop(self, location, names, ignore_missing, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer(c_int), intent(in) :: location
        character(*), intent(in) :: names(:)
        logical, intent(in), optional :: ignore_missing
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr
        integer(c_int64_t) :: count
        integer(c_int) :: cignore
        cignore = 0
        if (present(ignore_missing)) then
            if (ignore_missing) cignore = 1
        end if
        call c_str_array(names, storage, cptrs, arr, count)
        out%handle = c_mio_data_drop(self%handle, location, arr, count, cignore)
        if (.not. c_associated(out%handle)) then
            call handle_failure('data_drop', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Keep only the named data arrays at `location`, dropping the rest there.
    !> The other two locations are untouched; a zero-sized `names` drops all.
    function mesh_data_keep(self, location, names, ignore_missing, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer(c_int), intent(in) :: location
        character(*), intent(in) :: names(:)
        logical, intent(in), optional :: ignore_missing
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr
        integer(c_int64_t) :: count
        integer(c_int) :: cignore
        cignore = 0
        if (present(ignore_missing)) then
            if (ignore_missing) cignore = 1
        end if
        call c_str_array(names, storage, cptrs, arr, count)
        out%handle = c_mio_data_keep(self%handle, location, arr, count, cignore)
        if (.not. c_associated(out%handle)) then
            call handle_failure('data_keep', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Rename one data array, preserving its values, dtype and shape.
    function mesh_data_rename(self, location, from_name, to_name, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer(c_int), intent(in) :: location
        character(*), intent(in) :: from_name, to_name
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        out%handle = c_mio_data_rename(self%handle, location, c_str(from_name), c_str(to_name))
        if (.not. c_associated(out%handle)) then
            call handle_failure('data_rename', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Average point_data onto the cells (mean over each cell's own nodes).
    !> The output is always double precision. Omit `names` to convert all.
    function mesh_data_point_to_cell(self, names, suffix, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in), optional :: names(:)
        character(*), intent(in), optional :: suffix
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr
        integer(c_int64_t) :: count
        character(:), allocatable :: sfx
        character(len=1) :: empty(0)
        sfx = ''
        if (present(suffix)) sfx = suffix
        if (present(names)) then
            call c_str_array(names, storage, cptrs, arr, count)
        else
            call c_str_array(empty, storage, cptrs, arr, count)
        end if
        out%handle = c_mio_data_point_to_cell(self%handle, arr, count, c_str(sfx))
        if (.not. c_associated(out%handle)) then
            call handle_failure('data_point_to_cell', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Average cell_data onto the points. `weight` is MIO_WEIGHT_UNIFORM
    !> (default) or MIO_WEIGHT_MEASURE. Points touched by no cell get NaN.
    function mesh_data_cell_to_point(self, names, weight, suffix, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in), optional :: names(:)
        integer(c_int), intent(in), optional :: weight
        character(*), intent(in), optional :: suffix
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr
        integer(c_int64_t) :: count
        integer(c_int) :: cweight
        character(:), allocatable :: sfx
        character(len=1) :: empty(0)
        cweight = MIO_WEIGHT_UNIFORM
        if (present(weight)) cweight = weight
        sfx = ''
        if (present(suffix)) sfx = suffix
        if (present(names)) then
            call c_str_array(names, storage, cptrs, arr, count)
        else
            call c_str_array(empty, storage, cptrs, arr, count)
        end if
        out%handle = c_mio_data_cell_to_point(self%handle, arr, count, cweight, c_str(sfx))
        if (.not. c_associated(out%handle)) then
            call handle_failure('data_cell_to_point', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Evaluate an elementwise expression over the arrays at `location` and
    !> store the result there as `output_name`. The grammar accepts + - * /,
    !> unary minus, parentheses, numbers, array names, and the functions abs,
    !> sqrt, min, max and norm -- nothing else is evaluated.
    function mesh_data_calc(self, expression, output_name, location, overwrite, &
                            stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: expression, output_name
        integer(c_int), intent(in), optional :: location
        logical, intent(in), optional :: overwrite
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: cloc, cover
        cloc = MIO_DATA_POINT
        if (present(location)) cloc = location
        cover = 0
        if (present(overwrite)) then
            if (overwrite) cover = 1
        end if
        out%handle = c_mio_data_calc(self%handle, c_str(expression), cloc, &
                                     c_str(output_name), cover)
        if (.not. c_associated(out%handle)) then
            call handle_failure('data_calc', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Condition the values of the selected arrays. `mode` is MIO_COND_CLAMP
    !> (default), _NORMALIZE or _STANDARDIZE; `scope` is MIO_SCOPE_COMPONENT
    !> (default) or _MAGNITUDE. For cell_data the statistics are computed
    !> jointly across all cell blocks.
    function mesh_data_condition(self, location, names, mode, lo, hi, scope, &
                                 nan_policy, nan_replacement, suffix, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer(c_int), intent(in) :: location
        character(*), intent(in), optional :: names(:)
        integer(c_int), intent(in), optional :: mode, scope, nan_policy
        real(real64), intent(in), optional :: lo, hi, nan_replacement
        character(*), intent(in), optional :: suffix
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr
        integer(c_int64_t) :: count
        integer(c_int) :: cmode, cscope, cnan
        real(c_double) :: clo, chi, crep
        character(:), allocatable :: sfx
        character(len=1) :: empty(0)
        cmode = MIO_COND_CLAMP
        if (present(mode)) cmode = mode
        cscope = MIO_SCOPE_COMPONENT
        if (present(scope)) cscope = scope
        cnan = MIO_NAN_IGNORE
        if (present(nan_policy)) cnan = nan_policy
        clo = 0.0_c_double
        if (present(lo)) clo = real(lo, c_double)
        chi = 1.0_c_double
        if (present(hi)) chi = real(hi, c_double)
        crep = 0.0_c_double
        if (present(nan_replacement)) crep = real(nan_replacement, c_double)
        sfx = ''
        if (present(suffix)) sfx = suffix
        if (present(names)) then
            call c_str_array(names, storage, cptrs, arr, count)
        else
            call c_str_array(empty, storage, cptrs, arr, count)
        end if
        out%handle = c_mio_data_condition(self%handle, location, arr, count, cmode, clo, chi, &
                                          cscope, cnan, crep, c_str(sfx))
        if (.not. c_associated(out%handle)) then
            call handle_failure('data_condition', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> von Mises / principal / hydrostatic / deviatoric fields of a symmetric
    !> (6-component) or general 3x3 (9-component) tensor array. `names` empty
    !> processes every 6- or 9-component array at `location` (field_data is
    !> rejected). `outputs` is a bitwise OR of MIO_TINV_* (default MIO_TINV_ALL).
    function mesh_tensor_invariants(self, location, names, outputs, prefix, suffix, &
                                    overwrite, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer(c_int), intent(in) :: location
        character(*), intent(in), optional :: names(:)
        integer(c_int), intent(in), optional :: outputs
        character(*), intent(in), optional :: prefix, suffix
        logical, intent(in), optional :: overwrite
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr
        integer(c_int64_t) :: count
        integer(c_int) :: couts, cover
        character(:), allocatable :: pfx, sfx
        character(len=1) :: empty(0)
        couts = MIO_TINV_ALL
        if (present(outputs)) couts = outputs
        cover = 1
        if (present(overwrite)) then
            if (.not. overwrite) cover = 0
        end if
        pfx = ''
        if (present(prefix)) pfx = prefix
        sfx = ''
        if (present(suffix)) sfx = suffix
        if (present(names)) then
            call c_str_array(names, storage, cptrs, arr, count)
        else
            call c_str_array(empty, storage, cptrs, arr, count)
        end if
        out%handle = c_mio_tensor_invariants(self%handle, location, arr, count, couts, &
                                             c_str(pfx), c_str(sfx), cover)
        if (.not. c_associated(out%handle)) then
            call handle_failure('tensor_invariants', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Summarize every data array the mesh carries (read-only). Returns one
    !> mio_data_array_info per array; pass `keys` to also receive their names.
    function mesh_data_info(self, keys, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(len=STRBUF_LEN), allocatable, intent(out), optional :: keys(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_data_array_info), allocatable :: out(:)
        type(c_ptr) :: res
        integer(c_int64_t) :: count, i, n
        integer(c_int) :: s
        character(kind=c_char) :: buf(STRBUF_LEN)
        res = c_mio_data_info_create(self%handle)
        if (.not. c_associated(res)) then
            call handle_failure('data_info', mio_error_message(), stat, errmsg)
            allocate (out(0))
            return
        end if
        count = c_mio_data_info_count(res)
        if (count < 0) count = 0
        allocate (out(count))
        if (present(keys)) allocate (keys(count))
        do i = 1, count
            s = c_mio_data_info_entry(res, i - 1_c_int64_t, out(i))
            if (s /= 0_c_int) then
                call c_mio_data_info_free(res)
                call handle_failure('data_info', mio_error_message(), stat, errmsg)
                return
            end if
            if (present(keys)) then
                n = c_mio_data_info_name(res, i - 1_c_int64_t, buf, int(STRBUF_LEN, c_int64_t))
                keys(i) = ''
                if (n > 0) keys(i) = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
            end if
        end do
        call c_mio_data_info_free(res)
        call clear_status(stat, errmsg)
    end function

    !> Cell-measure-weighted total/mean of one or more cell_data arrays over
    !> the whole mesh -- gradient's integration counterpart. `arrays` absent
    !> means every cell_data array (sorted name order); a point_data-only name
    !> fails, naming data_point_to_cell as the fix. Returns one
    !> mio_field_integral_info per array; pass `keys` for their names and
    !> `totals`/`means`/`domain_measures`/`num_nans` for the flat,
    !> array-major, per-component buffers (length sum(out%num_components)).
    !> Per-region breakdown is a separate call: `data_integrate_region`.
    function mesh_data_integrate(self, arrays, keys, totals, means, domain_measures, num_nans, &
            stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in), optional :: arrays(:)
        character(len=STRBUF_LEN), allocatable, intent(out), optional :: keys(:)
        real(real64), allocatable, intent(out), optional :: totals(:), means(:), &
            domain_measures(:)
        integer(int64), allocatable, intent(out), optional :: num_nans(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_field_integral_info), allocatable :: out(:)
        type(c_ptr) :: res
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr
        integer(c_int64_t) :: count, i, n, k, comp, off, total_comp
        integer(c_int) :: s
        character(kind=c_char) :: buf(STRBUF_LEN)
        real(c_double) :: ctotal, cmean, cdomain
        integer(c_int64_t) :: cnan

        if (present(arrays)) then
            call c_str_array(arrays, storage, cptrs, arr, count)
        else
            arr = c_null_ptr
            count = 0_c_int64_t
        end if
        res = c_mio_data_integrate_create(self%handle, arr, count)
        if (.not. c_associated(res)) then
            call handle_failure('data_integrate', mio_error_message(), stat, errmsg)
            allocate (out(0))
            if (present(keys)) allocate (keys(0))
            if (present(totals)) allocate (totals(0))
            if (present(means)) allocate (means(0))
            if (present(domain_measures)) allocate (domain_measures(0))
            if (present(num_nans)) allocate (num_nans(0))
            return
        end if

        n = c_mio_data_integrate_count(res)
        if (n < 0) n = 0
        allocate (out(n))
        if (present(keys)) allocate (keys(n))
        do i = 1, n
            s = c_mio_data_integrate_entry(res, i - 1_c_int64_t, out(i))
            if (s /= 0_c_int) then
                call c_mio_data_integrate_free(res)
                call handle_failure('data_integrate', mio_error_message(), stat, errmsg)
                return
            end if
            if (present(keys)) then
                k = c_mio_data_integrate_name(res, i - 1_c_int64_t, buf, int(STRBUF_LEN, c_int64_t))
                keys(i) = ''
                if (k > 0) keys(i) = from_c_buf(buf, min(int(k), STRBUF_LEN - 1))
            end if
        end do

        if (present(totals) .or. present(means) .or. present(domain_measures) .or. &
                present(num_nans)) then
            total_comp = 0
            do i = 1, n
                total_comp = total_comp + out(i)%num_components
            end do
            if (present(totals)) allocate (totals(total_comp))
            if (present(means)) allocate (means(total_comp))
            if (present(domain_measures)) allocate (domain_measures(total_comp))
            if (present(num_nans)) allocate (num_nans(total_comp))
            off = 0
            do i = 1, n
                do comp = 0, out(i)%num_components - 1
                    s = c_mio_data_integrate_component(res, i - 1_c_int64_t, comp, ctotal, cmean, &
                                                        cdomain, cnan)
                    if (s /= 0_c_int) then
                        call c_mio_data_integrate_free(res)
                        call handle_failure('data_integrate', mio_error_message(), stat, errmsg)
                        return
                    end if
                    off = off + 1
                    if (present(totals)) totals(off) = real(ctotal, real64)
                    if (present(means)) means(off) = real(cmean, real64)
                    if (present(domain_measures)) domain_measures(off) = real(cdomain, real64)
                    if (present(num_nans)) num_nans(off) = int(cnan, int64)
                end do
            end do
        end if

        call c_mio_data_integrate_free(res)
        call clear_status(stat, errmsg)
    end function

    !> Per-named-Cell-region breakdown of one array's field integral (see
    !> `data_integrate`, which this complements): one mio_field_integral_info
    !> per region present on the mesh, plus the same optional flat,
    !> region-major, per-component buffers.
    function mesh_data_integrate_region(self, array_name, keys, totals, means, domain_measures, &
            num_nans, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: array_name
        character(len=STRBUF_LEN), allocatable, intent(out), optional :: keys(:)
        real(real64), allocatable, intent(out), optional :: totals(:), means(:), &
            domain_measures(:)
        integer(int64), allocatable, intent(out), optional :: num_nans(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_field_integral_info), allocatable :: out(:)
        type(c_ptr) :: res
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr
        integer(c_int64_t) :: count, n, k, comp, off, total_comp
        integer(c_int) :: s
        character(kind=c_char) :: buf(STRBUF_LEN)
        real(c_double) :: ctotal, cmean, cdomain
        integer(c_int64_t) :: cnan

        call c_str_array([array_name], storage, cptrs, arr, count)
        res = c_mio_data_integrate_create(self%handle, arr, count)
        if (.not. c_associated(res)) then
            call handle_failure('data_integrate_region', mio_error_message(), stat, errmsg)
            allocate (out(0))
            if (present(keys)) allocate (keys(0))
            if (present(totals)) allocate (totals(0))
            if (present(means)) allocate (means(0))
            if (present(domain_measures)) allocate (domain_measures(0))
            if (present(num_nans)) allocate (num_nans(0))
            return
        end if

        n = c_mio_data_integrate_region_count(res, 0_c_int64_t)
        if (n < 0) n = 0
        allocate (out(n))
        if (present(keys)) allocate (keys(n))
        do k = 1, n
            s = c_mio_data_integrate_region_entry(res, 0_c_int64_t, k - 1_c_int64_t, out(k))
            if (s /= 0_c_int) then
                call c_mio_data_integrate_free(res)
                call handle_failure('data_integrate_region', mio_error_message(), stat, errmsg)
                return
            end if
            if (present(keys)) then
                off = c_mio_data_integrate_region_name(res, 0_c_int64_t, k - 1_c_int64_t, buf, &
                                                        int(STRBUF_LEN, c_int64_t))
                keys(k) = ''
                if (off > 0) keys(k) = from_c_buf(buf, min(int(off), STRBUF_LEN - 1))
            end if
        end do

        if (present(totals) .or. present(means) .or. present(domain_measures) .or. &
                present(num_nans)) then
            total_comp = 0
            do k = 1, n
                total_comp = total_comp + out(k)%num_components
            end do
            if (present(totals)) allocate (totals(total_comp))
            if (present(means)) allocate (means(total_comp))
            if (present(domain_measures)) allocate (domain_measures(total_comp))
            if (present(num_nans)) allocate (num_nans(total_comp))
            off = 0
            do k = 1, n
                do comp = 0, out(k)%num_components - 1
                    s = c_mio_data_integrate_region_component(res, 0_c_int64_t, k - 1_c_int64_t, &
                                                               comp, ctotal, cmean, cdomain, cnan)
                    if (s /= 0_c_int) then
                        call c_mio_data_integrate_free(res)
                        call handle_failure('data_integrate_region', mio_error_message(), stat, &
                                             errmsg)
                        return
                    end if
                    off = off + 1
                    if (present(totals)) totals(off) = real(ctotal, real64)
                    if (present(means)) means(off) = real(cmean, real64)
                    if (present(domain_measures)) domain_measures(off) = real(cdomain, real64)
                    if (present(num_nans)) num_nans(off) = int(cnan, int64)
                end do
            end do
        end if

        call c_mio_data_integrate_free(res)
        call clear_status(stat, errmsg)
    end function

    !> Named regions: a group of points, cells or cell facets carried by a
    !> set-capable format (gmsh physical groups, Abaqus NSET / ELSET / SURFACE,
    !> ...). See doc/regions.md.
    !>
    !> `entries` is the flat int64 buffer, `sum(out%num_entries * out%stride)`
    !> long, with the groups laid out back to back in `out` order. Point and
    !> cell indices are shifted to Fortran's 1-based convention on the way out;
    !> the facet column of a side region is **not** shifted, matching the
    !> `partition_labels` rule that a value which is not an index stays as it is.
    function mesh_regions(self, keys, entries, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        character(len=STRBUF_LEN), allocatable, intent(out), optional :: keys(:)
        integer(int64), allocatable, intent(out), optional :: entries(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_region_info), allocatable :: out(:)
        type(c_ptr) :: res, p
        integer(c_int64_t) :: count, i, n, nvals, total, off, k
        integer(c_int) :: s
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t), pointer :: vals(:)

        res = c_mio_regions_create(self%handle)
        if (.not. c_associated(res)) then
            call handle_failure('regions', mio_error_message(), stat, errmsg)
            allocate (out(0))
            if (present(keys)) allocate (keys(0))
            if (present(entries)) allocate (entries(0))
            return
        end if
        count = c_mio_regions_count(res)
        if (count < 0) count = 0
        allocate (out(count))
        if (present(keys)) allocate (keys(count))

        total = 0
        do i = 1, count
            s = c_mio_regions_info(res, i - 1_c_int64_t, out(i))
            if (s /= 0_c_int) then
                call c_mio_regions_free(res)
                call handle_failure('regions', mio_error_message(), stat, errmsg)
                return
            end if
            total = total + out(i)%num_entries * out(i)%stride
            if (present(keys)) then
                n = c_mio_regions_name(res, i - 1_c_int64_t, buf, int(STRBUF_LEN, c_int64_t))
                keys(i) = ''
                if (n > 0) keys(i) = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
            end if
        end do

        if (present(entries)) then
            allocate (entries(total))
            off = 0
            do i = 1, count
                p = c_mio_regions_entries(res, i - 1_c_int64_t, nvals)
                if (.not. c_associated(p) .or. nvals <= 0) cycle
                call c_f_pointer(p, vals, [nvals])
                if (out(i)%stride == 2) then
                    ! (cell, facet) pairs: shift only the cell column.
                    do k = 1, nvals, 2
                        entries(off + k) = vals(k) + 1_c_int64_t
                        entries(off + k + 1) = vals(k + 1)
                    end do
                else
                    do k = 1, nvals
                        entries(off + k) = vals(k) + 1_c_int64_t
                    end do
                end if
                off = off + nvals
            end do
        end if

        call c_mio_regions_free(res)
        call clear_status(stat, errmsg)
    end function

    !> Add a named region, replacing any with the same (kind, name, dim, tag).
    !> `entries` is 1-based for point/cell indices; a side region takes flat
    !> (cell, facet) pairs whose facet column is passed through unshifted.
    subroutine mesh_add_region(self, name, kind, entries, dim, tag, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer(c_int), intent(in) :: kind
        integer(int64), intent(in) :: entries(:)
        integer, intent(in), optional :: dim
        integer(int64), intent(in), optional :: tag
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t), allocatable :: buf(:)
        integer(c_int) :: s, c_dim
        integer(c_int64_t) :: c_tag, k

        c_dim = -1_c_int
        if (present(dim)) c_dim = int(dim, c_int)
        c_tag = -1_c_int64_t
        if (present(tag)) c_tag = int(tag, c_int64_t)

        allocate (buf(size(entries)))
        if (kind == MIO_REGION_SIDE) then
            do k = 1, int(size(entries), c_int64_t), 2
                buf(k) = entries(k) - 1_c_int64_t
                if (k + 1 <= size(entries)) buf(k + 1) = entries(k + 1)
            end do
        else
            do k = 1, int(size(entries), c_int64_t)
                buf(k) = entries(k) - 1_c_int64_t
            end do
        end if

        s = c_mio_mesh_add_region(self%handle, c_str(name), kind, c_dim, c_tag, buf, &
                                  int(size(entries), c_int64_t))
        if (s /= 0_c_int) then
            call handle_failure('add_region', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end subroutine

    !> Connectivity bandwidth: max over cells of (max - min) node index.
    function mesh_compute_bandwidth(self, stat, errmsg) result(bw)
        class(mio_mesh), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int64) :: bw
        integer(c_int64_t) :: cbw
        cbw = c_mio_compute_bandwidth(self%handle)
        if (cbw < 0_c_int64_t) then
            call handle_failure('compute_bandwidth', mio_error_message(), stat, errmsg)
            bw = -1_int64
            return
        end if
        bw = int(cbw, int64)
        call clear_status(stat, errmsg)
    end function

    !> Are two meshes equal within tolerance (abs_err <= atol + rtol*|expected|)?
    !> `unordered` (default .false.) matches points by spatial proximity. Named
    !> point_sets/cell_sets are not compared.
    function mesh_equals(self, other, atol, rtol, unordered, stat, errmsg) result(eq)
        class(mio_mesh), intent(in) :: self
        type(mio_mesh), intent(in) :: other
        real(real64), intent(in), optional :: atol, rtol
        logical, intent(in), optional :: unordered
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        logical :: eq
        real(c_double) :: catol, crtol
        integer(c_int) :: cuno, ceq, s
        eq = .false.
        catol = 1.0e-12_c_double
        crtol = 1.0e-9_c_double
        cuno = 0_c_int
        if (present(atol)) catol = real(atol, c_double)
        if (present(rtol)) crtol = real(rtol, c_double)
        if (present(unordered)) then
            if (unordered) cuno = 1_c_int
        end if
        s = c_mio_meshes_equal(self%handle, other%handle, catol, crtol, cuno, ceq)
        if (s /= 0_c_int) then
            call handle_failure('equals', mio_error_message(), stat, errmsg)
            return
        end if
        eq = ceq /= 0_c_int
        call clear_status(stat, errmsg)
    end function

    !> Compare two meshes; returns the verdict (0 = identical, 1 = equal within
    !> tolerance, 2 = different). `unordered` matches points by proximity.
    function mesh_diff(self, other, atol, rtol, unordered, stat, errmsg) result(verdict)
        class(mio_mesh), intent(in) :: self
        type(mio_mesh), intent(in) :: other
        real(real64), intent(in), optional :: atol, rtol
        logical, intent(in), optional :: unordered
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer :: verdict
        real(c_double) :: catol, crtol
        integer(c_int) :: cuno, s
        type(c_ptr) :: res
        verdict = 2
        catol = 1.0e-12_c_double
        crtol = 1.0e-9_c_double
        cuno = 0_c_int
        if (present(atol)) catol = real(atol, c_double)
        if (present(rtol)) crtol = real(rtol, c_double)
        if (present(unordered)) then
            if (unordered) cuno = 1_c_int
        end if
        s = c_mio_diff(self%handle, other%handle, catol, crtol, cuno, res)
        if (s /= 0_c_int) then
            call handle_failure('diff', mio_error_message(), stat, errmsg)
            return
        end if
        verdict = int(c_mio_diff_result_verdict(res))
        call c_mio_diff_result_free(res)
        call clear_status(stat, errmsg)
    end function

    ! ------------------------------------------------------------------
    ! mio_mesh: building (setters copy; see module header for layout rules)
    ! ------------------------------------------------------------------

    !> Assign the point coordinates from a `points(dim, num_points)` array.
    subroutine mesh_set_points(self, points, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        real(real64), intent(in), contiguous, target :: points(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        call handle_status(c_mio_mesh_set_points(self%handle, MIO_FLOAT64, &
                                                 int(size(points, 2), c_int64_t), &
                                                 int(size(points, 1), c_int64_t), &
                                                 c_loc(points)), 'set_points', stat, errmsg)
    end subroutine

    !> Append one cell block from a 1-based `conn(nodes_per_cell, num_cells)`
    !> array (default-integer version).
    subroutine mesh_add_cell_block_i32(self, cell_type, conn, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: cell_type
        integer(int32), intent(in) :: conn(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t), allocatable, target :: shifted(:, :)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shifted = int(conn, c_int64_t) - 1_c_int64_t
        call handle_status(c_mio_mesh_add_cell_block(self%handle, c_str(cell_type), &
                                                     int(size(conn, 2), c_int64_t), &
                                                     int(size(conn, 1), c_int64_t), MIO_INT64, &
                                                     c_loc(shifted)), 'add_cell_block', &
                           stat, errmsg)
    end subroutine

    !> As above, for `integer(int64)` connectivity.
    subroutine mesh_add_cell_block_i64(self, cell_type, conn, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: cell_type
        integer(int64), intent(in) :: conn(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t), allocatable, target :: shifted(:, :)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shifted = int(conn, c_int64_t) - 1_c_int64_t
        call handle_status(c_mio_mesh_add_cell_block(self%handle, c_str(cell_type), &
                                                     int(size(conn, 2), c_int64_t), &
                                                     int(size(conn, 1), c_int64_t), MIO_INT64, &
                                                     c_loc(shifted)), 'add_cell_block', &
                           stat, errmsg)
    end subroutine

    !> Append a 1-level ragged (jagged polygon) cell block.
    !>
    !> `row_offsets(num_cells + 1)` and `nodes(:)` are BOTH 1-based on this
    !> side: `nodes(row_offsets(c) : row_offsets(c + 1) - 1)` is cell `c`'s node
    !> list. The offsets are shifted as well as the node ids because they index
    !> a Fortran array -- the same reasoning that makes MED's own INN 1-based.
    subroutine mesh_add_polygon_block(self, cell_type, row_offsets, nodes, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: cell_type
        integer(int64), intent(in) :: row_offsets(:), nodes(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t), allocatable :: c_off(:), c_nodes(:)
        if (size(row_offsets) < 1) then
            call handle_failure('add_polygon_block', &
                                'row_offsets must have num_cells + 1 entries', stat, errmsg)
            return
        end if
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        c_off = int(row_offsets, c_int64_t) - 1_c_int64_t
        c_nodes = int(nodes, c_int64_t) - 1_c_int64_t
        call handle_status(c_mio_mesh_add_polygon_block(self%handle, c_str(cell_type), &
                                                        int(size(row_offsets) - 1, c_int64_t), &
                                                        c_off, c_nodes, &
                                                        int(size(nodes), c_int64_t)), &
                           'add_polygon_block', stat, errmsg)
    end subroutine

    !> Append a 2-level ragged (polyhedron) cell block: each cell is a list of
    !> faces, each face a list of node ids.
    !>
    !> `cell_offsets(num_cells + 1)` indexes the FACE list,
    !> `face_offsets(num_faces + 1)` indexes `nodes`, and all three are 1-based
    !> here, so face `f` of cell `c` is
    !> `nodes(face_offsets(j) : face_offsets(j + 1) - 1)` with
    !> `j = cell_offsets(c) + f - 1`.
    subroutine mesh_add_polyhedron_block(self, cell_type, cell_offsets, face_offsets, nodes, &
                                         stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: cell_type
        integer(int64), intent(in) :: cell_offsets(:), face_offsets(:), nodes(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t), allocatable :: c_cells(:), c_faces(:), c_nodes(:)
        if (size(cell_offsets) < 1 .or. size(face_offsets) < 1) then
            call handle_failure('add_polyhedron_block', &
                                'cell_offsets and face_offsets must each have one more entry '// &
                                'than the count they describe', stat, errmsg)
            return
        end if
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        c_cells = int(cell_offsets, c_int64_t) - 1_c_int64_t
        c_faces = int(face_offsets, c_int64_t) - 1_c_int64_t
        c_nodes = int(nodes, c_int64_t) - 1_c_int64_t
        call handle_status(c_mio_mesh_add_polyhedron_block( &
                           self%handle, c_str(cell_type), &
                           int(size(cell_offsets) - 1, c_int64_t), c_cells, &
                           int(size(face_offsets) - 1, c_int64_t), c_faces, c_nodes, &
                           int(size(nodes), c_int64_t)), 'add_polyhedron_block', stat, errmsg)
    end subroutine

    !> Attach a scalar per-point field: `data(num_points)`.
    subroutine mesh_add_point_data_r1(self, name, data, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: name
        real(real64), intent(in), contiguous, target :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t) :: shape(1)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shape(1) = int(size(data), c_int64_t)
        call handle_status(c_mio_mesh_add_point_data(self%handle, c_str(name), MIO_FLOAT64, &
                                                     1_c_int32_t, shape, c_loc(data)), &
                           'add_point_data', stat, errmsg)
    end subroutine

    !> Attach a vector per-point field: `data(num_components, num_points)`
    !> (same memory as the C API's row-major `(num_points, num_components)`).
    subroutine mesh_add_point_data_r2(self, name, data, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: name
        real(real64), intent(in), contiguous, target :: data(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t) :: shape(2)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shape(1) = int(size(data, 2), c_int64_t)  ! C shape is the reverse of the Fortran one
        shape(2) = int(size(data, 1), c_int64_t)
        call handle_status(c_mio_mesh_add_point_data(self%handle, c_str(name), MIO_FLOAT64, &
                                                     2_c_int32_t, shape, c_loc(data)), &
                           'add_point_data', stat, errmsg)
    end subroutine

    !> Append the named cell-data field's array for the next cell block (call
    !> once per block, in block order): `data(num_cells_in_block)`.
    subroutine mesh_add_cell_data_r1(self, name, data, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: name
        real(real64), intent(in), contiguous, target :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t) :: shape(1)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shape(1) = int(size(data), c_int64_t)
        call handle_status(c_mio_mesh_append_cell_data(self%handle, c_str(name), MIO_FLOAT64, &
                                                       1_c_int32_t, shape, c_loc(data)), &
                           'add_cell_data', stat, errmsg)
    end subroutine

    !> Attach a named mesh-level (field) data array.
    subroutine mesh_add_field_data_r1(self, name, data, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: name
        real(real64), intent(in), contiguous, target :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t) :: shape(1)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shape(1) = int(size(data), c_int64_t)
        call handle_status(c_mio_mesh_add_field_data(self%handle, c_str(name), MIO_FLOAT64, &
                                                     1_c_int32_t, shape, c_loc(data)), &
                           'add_field_data', stat, errmsg)
    end subroutine

    ! ------------------------------------------------------------------
    ! mio_mesh: inspection
    ! ------------------------------------------------------------------

    integer(int64) function mesh_num_points(self)
        class(mio_mesh), intent(in) :: self
        mesh_num_points = c_mio_mesh_num_points(self%handle)
    end function

    integer(int64) function mesh_point_dim(self)
        class(mio_mesh), intent(in) :: self
        mesh_point_dim = c_mio_mesh_point_dim(self%handle)
    end function

    integer(int64) function mesh_num_cell_blocks(self)
        class(mio_mesh), intent(in) :: self
        mesh_num_cell_blocks = c_mio_mesh_num_cell_blocks(self%handle)
    end function

    !> Copy the coordinates into `points(point_dim, num_points)` (allocated
    !> here), converting to real64 if the mesh stores single precision.
    subroutine mesh_get_points(self, points, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        real(real64), allocatable, intent(out) :: points(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: p
        integer(c_int) :: dtype, status
        integer(int64) :: n, dim
        real(c_double), pointer :: d64(:)
        real(c_float), pointer :: d32(:)
        status = c_mio_mesh_get_points(self%handle, p, dtype)
        if (status /= 0_c_int) then
            call handle_status(status, 'get_points', stat, errmsg)
            return
        end if
        n = c_mio_mesh_num_points(self%handle)
        dim = c_mio_mesh_point_dim(self%handle)
        allocate (points(dim, n))
        if (n*dim == 0) then
            call clear_status(stat, errmsg)
            return
        end if
        select case (dtype)
        case (MIO_FLOAT64)
            call c_f_pointer(p, d64, [n*dim])
            points = reshape(d64, [dim, n])
        case (MIO_FLOAT32)
            call c_f_pointer(p, d32, [n*dim])
            points = reshape(real(d32, real64), [dim, n])
        case default
            call handle_failure('get_points', 'unexpected points dtype', stat, errmsg)
            return
        end select
        call clear_status(stat, errmsg)
    end subroutine

    !> Zero-copy borrow of the coordinates as `ptr(point_dim, num_points)`.
    !> Valid until the next mutating call on this mesh or free(); fails if
    !> the mesh stores anything but real64 (then use get_points).
    subroutine mesh_points_ptr(self, ptr, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        real(real64), pointer, intent(out) :: ptr(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: p
        integer(c_int) :: dtype, status
        ptr => null()
        status = c_mio_mesh_get_points(self%handle, p, dtype)
        if (status /= 0_c_int) then
            call handle_status(status, 'points_ptr', stat, errmsg)
            return
        end if
        if (dtype /= MIO_FLOAT64) then
            call handle_failure('points_ptr', 'points are not real64; use get_points', &
                                stat, errmsg)
            return
        end if
        call c_f_pointer(p, ptr, [c_mio_mesh_point_dim(self%handle), &
                                  c_mio_mesh_num_points(self%handle)])
        call clear_status(stat, errmsg)
    end subroutine

    !> The meshio type name of 1-based cell block `block` (e.g. "tetra10").
    function mesh_cell_block_type(self, block, stat, errmsg) result(type_name)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: type_name
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: n
        type_name = ''
        n = c_mio_mesh_cell_block_type(self%handle, int(block - 1, c_int64_t), buf, &
                                       int(STRBUF_LEN, c_int64_t))
        if (n < 0) then
            call handle_failure('cell_block_type', mio_error_message(), stat, errmsg)
            return
        end if
        type_name = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
        call clear_status(stat, errmsg)
    end function

    integer(int64) function mesh_cell_block_num_cells(self, block)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer(c_int64_t) :: nc, npc
        integer(c_int32_t) :: ragged
        mesh_cell_block_num_cells = -1
        if (c_mio_mesh_cell_block_info(self%handle, int(block - 1, c_int64_t), nc, npc, &
                                       ragged) == 0_c_int) mesh_cell_block_num_cells = nc
    end function

    integer(int64) function mesh_cell_block_nodes_per_cell(self, block)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer(c_int64_t) :: nc, npc
        integer(c_int32_t) :: ragged
        mesh_cell_block_nodes_per_cell = -1
        if (c_mio_mesh_cell_block_info(self%handle, int(block - 1, c_int64_t), nc, npc, &
                                       ragged) == 0_c_int) mesh_cell_block_nodes_per_cell = npc
    end function

    logical function mesh_cell_block_is_ragged(self, block)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer(c_int64_t) :: nc, npc
        integer(c_int32_t) :: ragged
        mesh_cell_block_is_ragged = .false.
        if (c_mio_mesh_cell_block_info(self%handle, int(block - 1, c_int64_t), nc, npc, &
                                       ragged) == 0_c_int) mesh_cell_block_is_ragged = ragged /= 0
    end function

    !> Whether 1-based cell block `block` is 2-level ragged (a polyhedron block:
    !> each cell a list of faces). Implies `cell_block_is_ragged`.
    logical function mesh_cell_block_is_polyhedron(self, block)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        type(c_mio_cell_block_info) :: info
        mesh_cell_block_is_polyhedron = .false.
        if (c_mio_mesh_cell_block_info_ex(self%handle, int(block - 1, c_int64_t), info) &
            == 0_c_int) mesh_cell_block_is_polyhedron = info%is_polyhedron /= 0
    end function

    !> Copy 1-based ragged cell block `block` into a 1-based CSR pair,
    !> `row_offsets(num_cells + 1)` and `nodes(:)`, both allocated here.
    !> Cell `c`'s node list is `nodes(row_offsets(c) : row_offsets(c + 1) - 1)`.
    !>
    !> Fails on a rectangular block (use `get_cell_block`) and on a polyhedron
    !> block (use `get_polyhedron_block`).
    subroutine mesh_get_polygon_block(self, block, row_offsets, nodes, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer(int64), allocatable, intent(out) :: row_offsets(:), nodes(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: poly, p
        type(c_mio_poly_conn_shape) :: shape
        integer(c_int64_t) :: n
        integer(c_int64_t), pointer :: vals(:)

        allocate (row_offsets(0), nodes(0))
        poly = c_mio_poly_conn_create(self%handle, int(block - 1, c_int64_t))
        if (.not. c_associated(poly)) then
            call handle_failure('get_polygon_block', mio_error_message(), stat, errmsg)
            return
        end if
        if (c_mio_poly_conn_get_shape(poly, shape) /= 0_c_int) then
            call c_mio_poly_conn_free(poly)
            call handle_failure('get_polygon_block', mio_error_message(), stat, errmsg)
            return
        end if
        if (shape%is_polyhedron /= 0) then
            call c_mio_poly_conn_free(poly)
            call handle_failure('get_polygon_block', &
                                'cell block is a polyhedron block; use get_polyhedron_block', &
                                stat, errmsg)
            return
        end if
        p = c_mio_poly_conn_face_offsets(poly, n)
        deallocate (row_offsets)
        allocate (row_offsets(max(n, 0_c_int64_t)))
        if (c_associated(p) .and. n > 0) then
            call c_f_pointer(p, vals, [n])
            row_offsets = vals + 1_int64
        end if
        p = c_mio_poly_conn_nodes(poly, n)
        deallocate (nodes)
        allocate (nodes(max(n, 0_c_int64_t)))
        if (c_associated(p) .and. n > 0) then
            call c_f_pointer(p, vals, [n])
            nodes = vals + 1_int64
        end if
        call c_mio_poly_conn_free(poly)
        call clear_status(stat, errmsg)
    end subroutine

    !> Copy 1-based polyhedron cell block `block` into a 1-based CSR triple,
    !> all allocated here: `cell_offsets(num_cells + 1)` indexes the FACE list,
    !> `face_offsets(num_faces + 1)` indexes `nodes`. Face `f` of cell `c` is
    !> `nodes(face_offsets(j) : face_offsets(j + 1) - 1)` with
    !> `j = cell_offsets(c) + f - 1`.
    subroutine mesh_get_polyhedron_block(self, block, cell_offsets, face_offsets, nodes, &
                                         stat, errmsg)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer(int64), allocatable, intent(out) :: cell_offsets(:), face_offsets(:), nodes(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: poly, p
        type(c_mio_poly_conn_shape) :: shape
        integer(c_int64_t) :: n
        integer(c_int64_t), pointer :: vals(:)

        allocate (cell_offsets(0), face_offsets(0), nodes(0))
        poly = c_mio_poly_conn_create(self%handle, int(block - 1, c_int64_t))
        if (.not. c_associated(poly)) then
            call handle_failure('get_polyhedron_block', mio_error_message(), stat, errmsg)
            return
        end if
        if (c_mio_poly_conn_get_shape(poly, shape) /= 0_c_int) then
            call c_mio_poly_conn_free(poly)
            call handle_failure('get_polyhedron_block', mio_error_message(), stat, errmsg)
            return
        end if
        if (shape%is_polyhedron == 0) then
            call c_mio_poly_conn_free(poly)
            call handle_failure('get_polyhedron_block', &
                                'cell block is a 1-level polygon block; use get_polygon_block', &
                                stat, errmsg)
            return
        end if
        p = c_mio_poly_conn_cell_offsets(poly, n)
        deallocate (cell_offsets)
        allocate (cell_offsets(max(n, 0_c_int64_t)))
        if (c_associated(p) .and. n > 0) then
            call c_f_pointer(p, vals, [n])
            cell_offsets = vals + 1_int64
        end if
        p = c_mio_poly_conn_face_offsets(poly, n)
        deallocate (face_offsets)
        allocate (face_offsets(max(n, 0_c_int64_t)))
        if (c_associated(p) .and. n > 0) then
            call c_f_pointer(p, vals, [n])
            face_offsets = vals + 1_int64
        end if
        p = c_mio_poly_conn_nodes(poly, n)
        deallocate (nodes)
        allocate (nodes(max(n, 0_c_int64_t)))
        if (c_associated(p) .and. n > 0) then
            call c_f_pointer(p, vals, [n])
            nodes = vals + 1_int64
        end if
        call c_mio_poly_conn_free(poly)
        call clear_status(stat, errmsg)
    end subroutine

    !> Copy 1-based cell block `block` into `conn(nodes_per_cell, num_cells)`
    !> (allocated here), shifting the 0-based core indices to 1-based.
    subroutine mesh_get_cell_block(self, block, conn, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer(int64), allocatable, intent(out) :: conn(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: p
        integer(c_int) :: dtype, status
        integer(c_int64_t) :: nc, npc
        integer(c_int32_t) :: ragged
        integer(c_int64_t), pointer :: i64(:)
        integer(c_int32_t), pointer :: i32(:)
        status = c_mio_mesh_cell_block_info(self%handle, int(block - 1, c_int64_t), nc, npc, &
                                            ragged)
        if (status == 0_c_int) status = c_mio_mesh_cell_block_conn(self%handle, &
                                                                   int(block - 1, c_int64_t), &
                                                                   p, dtype)
        if (status /= 0_c_int) then
            call handle_status(status, 'get_cell_block', stat, errmsg)
            return
        end if
        allocate (conn(npc, nc))
        if (nc*npc == 0) then
            call clear_status(stat, errmsg)
            return
        end if
        select case (dtype)
        case (MIO_INT64)
            call c_f_pointer(p, i64, [nc*npc])
            conn = reshape(i64, [npc, nc]) + 1_int64
        case (MIO_INT32)
            call c_f_pointer(p, i32, [nc*npc])
            conn = reshape(int(i32, int64), [npc, nc]) + 1_int64
        case default
            call handle_failure('get_cell_block', 'unexpected connectivity dtype', stat, errmsg)
            return
        end select
        call clear_status(stat, errmsg)
    end subroutine

    integer(int64) function mesh_num_point_data(self)
        class(mio_mesh), intent(in) :: self
        mesh_num_point_data = c_mio_mesh_num_point_data(self%handle)
    end function

    !> The 1-based `index`-th point-data name (ascending lexicographic order).
    function mesh_point_data_name(self, index, stat, errmsg) result(name)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: name
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: n
        name = ''
        n = c_mio_mesh_point_data_name(self%handle, int(index - 1, c_int64_t), buf, &
                                       int(STRBUF_LEN, c_int64_t))
        if (n < 0) then
            call handle_failure('point_data_name', mio_error_message(), stat, errmsg)
            return
        end if
        name = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
        call clear_status(stat, errmsg)
    end function

    !> Copy the named scalar point-data array into `data(:)` (allocated here).
    subroutine mesh_get_point_data_r1(self, name, data, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: name
        real(real64), allocatable, intent(out) :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call get_named_r1(self, 'point_data', name, -1_int64, data, stat, errmsg)
    end subroutine

    !> Copy the named vector point-data array into
    !> `data(num_components, num_points)` (allocated here).
    subroutine mesh_get_point_data_r2(self, name, data, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: name
        real(real64), allocatable, intent(out) :: data(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: p
        integer(c_int) :: dtype, status
        integer(c_int32_t) :: ndim
        integer(c_int64_t) :: shape(MIO_MAX_NDIM)
        status = c_mio_mesh_get_point_data(self%handle, c_str(name), p, dtype, ndim, shape)
        if (status /= 0_c_int) then
            call handle_status(status, 'get_point_data', stat, errmsg)
            return
        end if
        if (ndim /= 2_c_int32_t) then
            call handle_failure('get_point_data', 'point_data "'//trim(name)// &
                                '" is not rank-2', stat, errmsg)
            return
        end if
        call copy_out_r2(p, dtype, shape(1), shape(2), data, 'get_point_data', stat, errmsg)
    end subroutine

    integer(int64) function mesh_num_cell_data(self)
        class(mio_mesh), intent(in) :: self
        mesh_num_cell_data = c_mio_mesh_num_cell_data(self%handle)
    end function

    !> The 1-based `index`-th cell-data name (sorted, as above).
    function mesh_cell_data_name(self, index, stat, errmsg) result(name)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: name
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: n
        name = ''
        n = c_mio_mesh_cell_data_name(self%handle, int(index - 1, c_int64_t), buf, &
                                      int(STRBUF_LEN, c_int64_t))
        if (n < 0) then
            call handle_failure('cell_data_name', mio_error_message(), stat, errmsg)
            return
        end if
        name = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
        call clear_status(stat, errmsg)
    end function

    integer(int64) function mesh_cell_data_num_blocks(self, name)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: name
        mesh_cell_data_num_blocks = c_mio_mesh_cell_data_num_blocks(self%handle, c_str(name))
    end function

    !> Copy the named cell-data field's array for 1-based cell block `block`
    !> into `data(:)` (allocated here).
    subroutine mesh_get_cell_data_r1(self, name, block, data, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: name
        integer, intent(in) :: block
        real(real64), allocatable, intent(out) :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call get_named_r1(self, 'cell_data', name, int(block - 1, int64), data, stat, errmsg)
    end subroutine

    integer(int64) function mesh_num_field_data(self)
        class(mio_mesh), intent(in) :: self
        mesh_num_field_data = c_mio_mesh_num_field_data(self%handle)
    end function

    !> The 1-based `index`-th field-data name (sorted, as above).
    function mesh_field_data_name(self, index, stat, errmsg) result(name)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: name
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: n
        name = ''
        n = c_mio_mesh_field_data_name(self%handle, int(index - 1, c_int64_t), buf, &
                                       int(STRBUF_LEN, c_int64_t))
        if (n < 0) then
            call handle_failure('field_data_name', mio_error_message(), stat, errmsg)
            return
        end if
        name = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
        call clear_status(stat, errmsg)
    end function

    !> Copy the named field-data array into `data(:)` (allocated here).
    subroutine mesh_get_field_data_r1(self, name, data, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: name
        real(real64), allocatable, intent(out) :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call get_named_r1(self, 'field_data', name, -2_int64, data, stat, errmsg)
    end subroutine

    ! ------------------------------------------------------------------
    ! Private implementation helpers
    ! ------------------------------------------------------------------

    !> Setters on a never-created handle allocate one on the fly, so
    !> `type(mio_mesh) :: m; call m%set_points(...)` just works.
    subroutine ensure_handle(self, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (.not. c_associated(self%handle)) call mesh_create(self, stat, errmsg)
    end subroutine

    !> Shared rank-1 copy-getter over the three named-data families
    !> (block >= 0: cell_data block; -1: point_data; -2: field_data).
    subroutine get_named_r1(self, family, name, block, data, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: family, name
        integer(int64), intent(in) :: block
        real(real64), allocatable, intent(out) :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: p
        integer(c_int) :: dtype, status
        integer(c_int32_t) :: ndim
        integer(c_int64_t) :: shape(MIO_MAX_NDIM)
        select case (family)
        case ('point_data')
            status = c_mio_mesh_get_point_data(self%handle, c_str(name), p, dtype, ndim, shape)
        case ('cell_data')
            status = c_mio_mesh_get_cell_data(self%handle, c_str(name), block, p, dtype, ndim, &
                                              shape)
        case default
            status = c_mio_mesh_get_field_data(self%handle, c_str(name), p, dtype, ndim, shape)
        end select
        if (status /= 0_c_int) then
            call handle_status(status, 'get_'//family, stat, errmsg)
            return
        end if
        if (ndim /= 1_c_int32_t) then
            call handle_failure('get_'//family, family//' "'//trim(name)// &
                                '" is not rank-1 (use the rank-2 getter)', stat, errmsg)
            return
        end if
        call copy_out_r1(p, dtype, shape(1), data, 'get_'//family, stat, errmsg)
    end subroutine

    !> Copy a C array of `n` elements of any supported dtype into a real64
    !> rank-1 allocatable.
    subroutine copy_out_r1(p, dtype, n, data, what, stat, errmsg)
        type(c_ptr), intent(in) :: p
        integer(c_int), intent(in) :: dtype
        integer(c_int64_t), intent(in) :: n
        real(real64), allocatable, intent(out) :: data(:)
        character(*), intent(in) :: what
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(c_double), pointer :: d64(:)
        real(c_float), pointer :: d32(:)
        integer(c_int64_t), pointer :: i64(:)
        integer(c_int32_t), pointer :: i32(:)
        allocate (data(n))
        if (n == 0) then
            call clear_status(stat, errmsg)
            return
        end if
        select case (dtype)
        case (MIO_FLOAT64)
            call c_f_pointer(p, d64, [n]); data = d64
        case (MIO_FLOAT32)
            call c_f_pointer(p, d32, [n]); data = real(d32, real64)
        case (MIO_INT64)
            call c_f_pointer(p, i64, [n]); data = real(i64, real64)
        case (MIO_INT32)
            call c_f_pointer(p, i32, [n]); data = real(i32, real64)
        case default
            call handle_failure(what, 'unsupported dtype for real64 copy', stat, errmsg)
            return
        end select
        call clear_status(stat, errmsg)
    end subroutine

    !> As copy_out_r1 for a C row-major `(n0, n1)` array, delivered as the
    !> Fortran `(n1, n0)` view of the same memory order.
    subroutine copy_out_r2(p, dtype, n0, n1, data, what, stat, errmsg)
        type(c_ptr), intent(in) :: p
        integer(c_int), intent(in) :: dtype
        integer(c_int64_t), intent(in) :: n0, n1
        real(real64), allocatable, intent(out) :: data(:, :)
        character(*), intent(in) :: what
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(c_double), pointer :: d64(:)
        real(c_float), pointer :: d32(:)
        integer(c_int64_t), pointer :: i64(:)
        integer(c_int32_t), pointer :: i32(:)
        allocate (data(n1, n0))
        if (n0*n1 == 0) then
            call clear_status(stat, errmsg)
            return
        end if
        select case (dtype)
        case (MIO_FLOAT64)
            call c_f_pointer(p, d64, [n0*n1]); data = reshape(d64, [n1, n0])
        case (MIO_FLOAT32)
            call c_f_pointer(p, d32, [n0*n1]); data = reshape(real(d32, real64), [n1, n0])
        case (MIO_INT64)
            call c_f_pointer(p, i64, [n0*n1]); data = reshape(real(i64, real64), [n1, n0])
        case (MIO_INT32)
            call c_f_pointer(p, i32, [n0*n1]); data = reshape(real(i32, real64), [n1, n0])
        case default
            call handle_failure(what, 'unsupported dtype for real64 copy', stat, errmsg)
            return
        end select
        call clear_status(stat, errmsg)
    end subroutine

    ! ------------------------------------------------------------------
    ! mio_xdmf_series: transient (time-series) XDMF writing
    !
    ! Same lifecycle rules as mio_mesh: create() replaces whatever the handle
    ! held, free() is idempotent, and nothing is freed implicitly.
    ! ------------------------------------------------------------------

    !> Open a transient XDMF series for writing. Nothing is written yet.
    !>
    !> `data_format` is `"HDF"` (the default; needs an HDF5-enabled build),
    !> `"XML"` (inline) or `"Binary"`. `gzip_level` applies to `"HDF"` datasets
    !> only; negative (the default) means no compression. An unknown format, or
    !> `"HDF"` without HDF5 support, fails through `stat`/`errmsg`.
    subroutine exodus_series_create(self, path, stat, errmsg)
        class(mio_exodus_series), intent(inout) :: self
        character(*), intent(in) :: path
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call self%free()
        self%handle = c_mio_exodus_series_create(c_str(path))
        if (.not. c_associated(self%handle)) then
            call handle_failure('exodus_series create', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end subroutine

    subroutine exodus_series_free(self)
        class(mio_exodus_series), intent(inout) :: self
        if (c_associated(self%handle)) call c_mio_exodus_series_free(self%handle)
        self%handle = c_null_ptr
    end subroutine

    logical function exodus_series_is_valid(self)
        class(mio_exodus_series), intent(in) :: self
        exodus_series_is_valid = c_associated(self%handle)
    end function

    subroutine exodus_series_write_points_cells(self, mesh, stat, errmsg)
        class(mio_exodus_series), intent(in) :: self
        class(mio_mesh), intent(in) :: mesh
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call handle_status(c_mio_exodus_series_write_points_cells(self%handle, mesh%handle), &
                           'exodus_series write_points_cells', stat, errmsg)
    end subroutine

    subroutine exodus_series_write_data(self, time, mesh, stat, errmsg)
        class(mio_exodus_series), intent(in) :: self
        real(real64), intent(in) :: time
        class(mio_mesh), intent(in) :: mesh
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call handle_status(c_mio_exodus_series_write_data(self%handle, real(time, c_double), mesh%handle), &
                           'exodus_series write_data', stat, errmsg)
    end subroutine

    subroutine exodus_series_flush(self, stat, errmsg)
        class(mio_exodus_series), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call handle_status(c_mio_exodus_series_flush(self%handle), 'exodus_series flush', stat, errmsg)
    end subroutine

    subroutine exodus_series_finalize(self, stat, errmsg)
        class(mio_exodus_series), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call handle_status(c_mio_exodus_series_finalize(self%handle), 'exodus_series finalize', stat, errmsg)
    end subroutine

    function exodus_series_finalized(self) result(f)
        class(mio_exodus_series), intent(in) :: self
        logical :: f
        f = .true.
        if (c_associated(self%handle)) f = c_mio_exodus_series_finalized(self%handle) == 1_c_int32_t
    end function

    function exodus_series_num_steps(self) result(n)
        class(mio_exodus_series), intent(in) :: self
        integer(int64) :: n
        n = 0_int64
        if (c_associated(self%handle)) n = int(c_mio_exodus_series_num_steps(self%handle), int64)
    end function

    subroutine xdmf_series_create(self, path, data_format, gzip_level, mode, auto_flush, &
                                  stat, errmsg)
        class(mio_xdmf_series), intent(inout) :: self
        character(*), intent(in) :: path
        character(*), intent(in), optional :: data_format
        integer, intent(in), optional :: gzip_level
        !> 'truncate' (default) or 'append'.
        character(*), intent(in), optional :: mode
        !> Flush the light data after every write_data (default .false.).
        logical, intent(in), optional :: auto_flush
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: fmt, md
        character(len=:, kind=c_char), allocatable, target :: fmt_buf
        type(mio_xdmf_series_opts_t) :: opts
        integer(c_int32_t) :: level
        call xdmf_series_free(self)
        fmt = 'HDF'; if (present(data_format)) fmt = data_format
        md = 'truncate'; if (present(mode)) md = mode
        level = -1_c_int32_t
        if (present(gzip_level)) level = int(gzip_level, c_int32_t)

        if (md /= 'truncate' .and. md /= 'append') then
            call handle_failure('xdmf_series create', &
                                "mode must be 'truncate' or 'append', got '"//md//"'", stat, errmsg)
            return
        end if

        if (.not. present(mode) .and. .not. present(auto_flush)) then
            self%handle = c_mio_xdmf_series_create(c_str(path), c_str(fmt), level)
        else
            ! The NUL-terminated format string must outlive the call, so it is
            ! held in a TARGET local rather than built inline.
            fmt_buf = fmt//c_null_char
            opts%data_format = c_loc(fmt_buf(1:1))
            opts%gzip_level = level
            if (md == 'append') opts%mode = 1_c_int32_t
            if (present(auto_flush)) then
                if (auto_flush) opts%auto_flush = 1_c_int32_t
            end if
            self%handle = c_mio_xdmf_series_create_ex(c_str(path), opts)
        end if
        if (.not. c_associated(self%handle)) then
            call handle_failure('xdmf_series create', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end subroutine

    !> Write the `.xdmf` as it currently stands, without finalizing, so a run
    !> that is killed or still going leaves a readable file covering every
    !> flushed step. Safe to call repeatedly; a no-op once finalized.
    subroutine xdmf_series_flush(self, stat, errmsg)
        class(mio_xdmf_series), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (.not. c_associated(self%handle)) then
            call handle_failure('xdmf_series flush', 'series handle is not open', stat, errmsg)
            return
        end if
        call handle_status(c_mio_xdmf_series_flush(self%handle), 'xdmf_series flush', stat, errmsg)
    end subroutine

    !> .true. once finalize() has run for this series (also on a closed handle).
    function xdmf_series_finalized(self) result(f)
        class(mio_xdmf_series), intent(in) :: self
        logical :: f
        f = .true.
        if (.not. c_associated(self%handle)) return
        f = c_mio_xdmf_series_finalized(self%handle) == 1_c_int32_t
    end function

    !> Finalize (if it has not happened yet) and release the series. Idempotent.
    !> A write failure during the implicit finalize is swallowed here -- call
    !> `finalize()` explicitly to see it.
    subroutine xdmf_series_free(self)
        class(mio_xdmf_series), intent(inout) :: self
        if (c_associated(self%handle)) call c_mio_xdmf_series_free(self%handle)
        self%handle = c_null_ptr
    end subroutine

    !> .true. between a successful create() and free().
    logical function xdmf_series_is_valid(self)
        class(mio_xdmf_series), intent(in) :: self
        xdmf_series_is_valid = c_associated(self%handle)
    end function

    !> Write the static grid every step shares. Call once, before the first
    !> `write_data`. Uses the mesh's points, cells and fixed named regions.
    subroutine xdmf_series_write_points_cells(self, mesh, stat, errmsg)
        class(mio_xdmf_series), intent(in) :: self
        class(mio_mesh), intent(in) :: mesh
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (.not. c_associated(self%handle)) then
            call handle_failure('xdmf_series write_points_cells', &
                                'series handle is not open', stat, errmsg)
            return
        end if
        call handle_status(c_mio_xdmf_series_write_points_cells(self%handle, mesh%handle), &
                           'xdmf_series write_points_cells', stat, errmsg)
    end subroutine

    !> Write one time step's point_data and cell_data at time `time`. The mesh's
    !> geometry is ignored, so a solver can pass the same object it updates in
    !> place; its cell blocks must match those of the static grid.
    subroutine xdmf_series_write_data(self, time, mesh, stat, errmsg)
        class(mio_xdmf_series), intent(in) :: self
        real(real64), intent(in) :: time
        class(mio_mesh), intent(in) :: mesh
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (.not. c_associated(self%handle)) then
            call handle_failure('xdmf_series write_data', 'series handle is not open', &
                                stat, errmsg)
            return
        end if
        call handle_status(c_mio_xdmf_series_write_data(self%handle, &
                                                        real(time, c_double), mesh%handle), &
                           'xdmf_series write_data', stat, errmsg)
    end subroutine

    !> Write the `.xdmf` and close the heavy-data container. Idempotent, and
    !> `free()` does it too -- call this explicitly so a write failure surfaces
    !> through `stat` instead of being swallowed.
    subroutine xdmf_series_finalize(self, stat, errmsg)
        class(mio_xdmf_series), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (.not. c_associated(self%handle)) then
            call handle_failure('xdmf_series finalize', 'series handle is not open', &
                                stat, errmsg)
            return
        end if
        call handle_status(c_mio_xdmf_series_finalize(self%handle), &
                           'xdmf_series finalize', stat, errmsg)
    end subroutine

    !> Number of steps written so far (0 on a closed handle).
    function xdmf_series_num_steps(self) result(n)
        class(mio_xdmf_series), intent(in) :: self
        integer(int64) :: n
        n = 0_int64
        if (.not. c_associated(self%handle)) return
        n = int(c_mio_xdmf_series_num_steps(self%handle), int64)
        if (n < 0_int64) n = 0_int64
    end function

    ! ----------------------------------------------------------------------
    ! Sequences (multi-file / transient datasets). See doc/sequences.md.
    !
    ! Indices are 1-based on this side, like every other Fortran accessor;
    ! the C ABI's 0-based index is produced by the `- 1` inside each wrapper.
    ! ----------------------------------------------------------------------

    !> Open a sequence from a glob pattern (`*` and `?` only; the directory
    !> part is taken literally). A pattern matching nothing is an error, never
    !> an empty sequence. Replaces any sequence this handle already held.
    subroutine sequence_open(self, pattern, stat, errmsg)
        class(mio_sequence), intent(inout) :: self
        character(*), intent(in) :: pattern
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call sequence_free(self)
        self%handle = c_mio_sequence_open(c_str(pattern))
        if (.not. c_associated(self%handle)) then
            call handle_failure('sequence open', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end subroutine

    !> Open a sequence from an explicit, ordered path list. The order is the
    !> caller's and is kept -- a stated order is not second-guessed.
    subroutine sequence_open_list(self, paths, stat, errmsg)
        class(mio_sequence), intent(inout) :: self
        character(*), intent(in) :: paths(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        ! `storage`/`cptrs` are the actual backing memory for the C string
        ! array and must stay in scope for the whole call -- see c_str_array.
        character(kind=c_char), allocatable, target :: storage(:, :)
        type(c_ptr), allocatable, target :: cptrs(:)
        type(c_ptr) :: arr
        integer(c_int64_t) :: count
        call sequence_free(self)
        if (size(paths) <= 0) then
            call handle_failure('sequence open_list', 'the path list is empty', stat, errmsg)
            return
        end if
        call c_str_array(paths, storage, cptrs, arr, count)
        self%handle = c_mio_sequence_open_list(cptrs, count)
        if (.not. c_associated(self%handle)) then
            call handle_failure('sequence open_list', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end subroutine

    !> Release the sequence handle. Meshes it produced are unaffected.
    subroutine sequence_free(self)
        class(mio_sequence), intent(inout) :: self
        if (c_associated(self%handle)) call c_mio_sequence_free(self%handle)
        self%handle = c_null_ptr
    end subroutine

    !> Whether this handle currently holds a sequence.
    function sequence_is_valid(self) result(v)
        class(mio_sequence), intent(in) :: self
        logical :: v
        v = c_associated(self%handle)
    end function

    !> Number of steps (0 on a closed handle).
    function sequence_count(self) result(n)
        class(mio_sequence), intent(in) :: self
        integer(int64) :: n
        n = 0_int64
        if (.not. c_associated(self%handle)) return
        n = int(c_mio_sequence_count(self%handle), int64)
        if (n < 0_int64) n = 0_int64
    end function

    !> Entry `index`'s file path (1-based; empty on error).
    function sequence_path(self, index) result(p)
        class(mio_sequence), intent(in) :: self
        integer, intent(in) :: index
        character(:), allocatable :: p
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: n
        p = ''
        if (.not. c_associated(self%handle)) return
        n = c_mio_sequence_path(self%handle, int(index - 1, c_int64_t), buf, &
                                int(STRBUF_LEN, c_int64_t))
        if (n > 0) p = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
    end function

    !> Entry `index`'s step WITHIN its own file (0 for a single-step file).
    function sequence_step(self, index) result(k)
        class(mio_sequence), intent(in) :: self
        integer, intent(in) :: index
        integer(int64) :: k
        k = -1_int64
        if (.not. c_associated(self%handle)) return
        k = int(c_mio_sequence_step(self%handle, int(index - 1, c_int64_t)), int64)
    end function

    !> Entry `index`'s time value.
    function sequence_time(self, index) result(t)
        class(mio_sequence), intent(in) :: self
        integer, intent(in) :: index
        real(real64) :: t
        real(c_double) :: v
        t = 0.0_real64
        if (.not. c_associated(self%handle)) return
        if (c_mio_sequence_time(self%handle, int(index - 1, c_int64_t), v) == 0) &
            t = real(v, real64)
    end function

    !> Where entry `index`'s time came from: 0 explicit, 1 file, 2 filename,
    !> 3 index (the fallback). Reported so a caller can tell "the file said
    !> 0.25" from "nothing said anything, so this is position 3".
    function sequence_time_source(self, index) result(src)
        class(mio_sequence), intent(in) :: self
        integer, intent(in) :: index
        integer :: src
        src = -1
        if (.not. c_associated(self%handle)) return
        src = int(c_mio_sequence_time_source(self%handle, int(index - 1, c_int64_t)))
    end function

    !> Read entry `index` (1-based). The result is OWNED: free it with
    !> `%free()`. The sequence caches nothing, which is what keeps a long
    !> dataset traversable.
    function sequence_read_step(self, index, stat, errmsg) result(m)
        class(mio_sequence), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: m
        if (.not. c_associated(self%handle)) then
            call handle_failure('sequence read_step', 'sequence handle is not open', &
                                stat, errmsg)
            return
        end if
        m%handle = c_mio_sequence_read(self%handle, int(index - 1, c_int64_t))
        if (.not. c_associated(m%handle)) then
            call handle_failure('sequence read_step', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Fan-in: write every step into one multi-step file. Streams -- one mesh
    !> alive at a time. A format that cannot hold a series fails by name.
    subroutine sequence_to_timeseries(self, out_path, out_format, stat, errmsg)
        class(mio_sequence), intent(in) :: self
        character(*), intent(in) :: out_path
        character(*), intent(in), optional :: out_format
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: fmt
        fmt = ''
        if (present(out_format)) fmt = out_format
        if (.not. c_associated(self%handle)) then
            call handle_failure('sequence to_timeseries', 'sequence handle is not open', &
                                stat, errmsg)
            return
        end if
        call handle_status(c_mio_sequence_to_timeseries(self%handle, c_str(out_path), &
                                                        c_str(fmt)), &
                           'sequence to_timeseries', stat, errmsg)
    end subroutine

    !> Resample the sequence onto `times` and write it (v16.25.0): one file per
    !> target time when `out_path` carries '{step}'/'{index}', else a series.
    !> `method` is 'linear' (default), 'nearest' or 'previous'; `clamp` takes
    !> the end steps for a time outside the source range (else it fails).
    subroutine sequence_resample(self, out_path, times, method, clamp, blend_points, &
                                 out_format, stat, errmsg)
        class(mio_sequence), intent(in) :: self
        character(*), intent(in) :: out_path
        real(real64), intent(in), target :: times(:)
        character(*), intent(in), optional :: method, out_format
        logical, intent(in), optional :: clamp, blend_points
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_resample_opts_t) :: opts
        real(c_double), allocatable, target :: t(:)
        character(:), allocatable :: fmt
        fmt = ''
        if (present(out_format)) fmt = out_format
        if (present(method)) then
            select case (trim(method))
            case ('linear')
                opts%method = 0
            case ('nearest')
                opts%method = 1
            case ('previous')
                opts%method = 2
            case default
                call handle_failure('sequence resample', "meshio++: resample: unknown method '"// &
                                    trim(method)//"'", stat, errmsg)
                return
            end select
        end if
        if (present(clamp)) opts%extrapolate = merge(1_c_int32_t, 0_c_int32_t, clamp)
        if (present(blend_points)) opts%blend_points = merge(1_c_int32_t, 0_c_int32_t, blend_points)
        t = real(times, c_double)
        opts%num_times = int(size(t), c_int64_t)
        if (size(t) > 0) opts%times = c_loc(t(1))
        call handle_status(c_mio_sequence_resample(self%handle, c_str(out_path), c_str(fmt), opts), &
                           'sequence resample', stat, errmsg)
    end subroutine

    !> Fan-out: write each step of a multi-step file to `out_pattern`, which
    !> must contain '{step}' or '{index}'.
    subroutine mio_timeseries_to_sequence(in_path, out_pattern, in_format, out_format, &
                                          stat, errmsg)
        character(*), intent(in) :: in_path
        character(*), intent(in) :: out_pattern
        character(*), intent(in), optional :: in_format
        character(*), intent(in), optional :: out_format
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: ifmt, ofmt
        ifmt = ''
        ofmt = ''
        if (present(in_format)) ifmt = in_format
        if (present(out_format)) ofmt = out_format
        call handle_status(c_mio_timeseries_to_sequence(c_str(in_path), c_str(ifmt), &
                                                        c_str(out_pattern), c_str(ofmt)), &
                           'timeseries_to_sequence', stat, errmsg)
    end subroutine

    !> Run a sequence settings document (the pipeline schema plus Mode /
    !> Input.Pattern / Input.Paths / Input.Times / Input.TimeFrom). A document
    !> using none of those behaves exactly as `mio_pipeline_run_file`.
    subroutine mio_sequence_pipeline_run_file(settings_path, stat, errmsg)
        character(*), intent(in) :: settings_path
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call handle_status(c_mio_sequence_pipeline_run_file(c_str(settings_path)), &
                           'sequence_pipeline_run_file', stat, errmsg)
    end subroutine

    !> `mio_sequence_pipeline_run_file` over JSON text.
    subroutine mio_sequence_pipeline_run_json(json_text, stat, errmsg)
        character(*), intent(in) :: json_text
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call handle_status(c_mio_sequence_pipeline_run_json(c_str(json_text)), &
                           'sequence_pipeline_run_json', stat, errmsg)
    end subroutine

    !> Open a provenance scope for this thread. Every mio write made before the
    !> matching mio_provenance_end records the richer block; with no scope open
    !> a writer emits only the one-line credit. Explicit begin/end rather than a
    !> block construct, matching this module's own explicit m%free() convention.
    subroutine mio_provenance_begin(mode, stat, errmsg)
        integer, intent(in), optional :: mode
        integer, intent(out), optional :: stat
        character(len=:), allocatable, intent(out), optional :: errmsg
        integer(c_int) :: m, rc
        m = MIO_PROVENANCE_BEST_EFFORT
        if (present(mode)) m = int(mode, c_int)
        rc = c_mio_provenance_scope_begin(m)
        call handle_status(rc, "mio_provenance_begin", stat, errmsg)
    end subroutine

    !> Close the most recently opened provenance scope. A stray call with none
    !> open is a no-op.
    subroutine mio_provenance_end()
        call c_mio_provenance_scope_end()
    end subroutine

    !> Record one conversion assumption. A no-op outside a scope, so it is safe
    !> to call unconditionally.
    subroutine mio_provenance_note(category, detail)
        character(*), intent(in) :: category, detail
        call c_mio_provenance_note(c_str(category), c_str(detail))
    end subroutine

    !> Record where the mesh came from. A no-op outside a scope.
    subroutine mio_provenance_source(path, format)
        character(*), intent(in) :: path, format
        call c_mio_provenance_set_source(c_str(path), c_str(format))
    end subroutine

    !> Record what was actually written (not what was requested). A no-op
    !> outside a scope.
    subroutine mio_provenance_target(format, encoding, codec, float_format)
        character(*), intent(in) :: format
        character(*), intent(in), optional :: encoding, codec, float_format
        character(len=:), allocatable :: e, c, f
        e = ""; c = ""; f = ""
        if (present(encoding)) e = encoding
        if (present(codec)) c = codec
        if (present(float_format)) f = float_format
        call c_mio_provenance_set_target(c_str(format), c_str(e), c_str(c), c_str(f))
    end subroutine

    ! ------------------------------------------------------------------
    ! mio_format_info: format side channels (v16.27.0)
    ! ------------------------------------------------------------------

    !> Read a mesh keeping its format's side channel in `info` (currently MDPA's
    !> tables, geometries, Mesh blocks, sub-model-part data and raw blocks).
    !> A format without one reads as `m%read` and leaves `info` invalid.
    subroutine mesh_read_with_info(self, path, info, format, lenient, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: path
        class(mio_format_info), intent(inout) :: info
        character(*), intent(in), optional :: format
        !> Skip (with a warning) what not even the side channel holds.
        logical, intent(in), optional :: lenient
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: fmt
        type(mio_read_opts_t), target :: opts
        type(c_ptr) :: h, ih
        fmt = ''; if (present(format)) fmt = format
        call c_mio_read_opts_init(opts)
        if (present(lenient)) then
            if (lenient) opts%lenient = 1
        end if
        h = c_mio_read_with_info(c_str(path), c_str(fmt), c_loc(opts), ih)
        if (.not. c_associated(h)) then
            call handle_failure('read_with_info', mio_error_message(), stat, errmsg)
            return
        end if
        call mesh_free(self)
        self%handle = h
        call format_info_free(info)
        info%handle = ih
        call clear_status(stat, errmsg)
    end subroutine

    !> Write a mesh restoring the side channel `info` came back with. An
    !> invalid (empty) `info` writes exactly as `m%write`.
    subroutine mesh_write_with_info(self, path, info, format, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: path
        class(mio_format_info), intent(in) :: info
        character(*), intent(in), optional :: format
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: fmt
        fmt = ''; if (present(format)) fmt = format
        call handle_status(c_mio_write_with_info(c_str(path), self%handle, c_str(fmt), &
                                                 info%handle), 'write_with_info', stat, errmsg)
    end subroutine

    !> Release the side channel. Idempotent.
    subroutine format_info_free(self)
        class(mio_format_info), intent(inout) :: self
        if (c_associated(self%handle)) call c_mio_format_info_free(self%handle)
        self%handle = c_null_ptr
    end subroutine

    !> .true. when a read left a side channel here.
    logical function format_info_is_valid(self)
        class(mio_format_info), intent(in) :: self
        format_info_is_valid = c_associated(self%handle)
    end function

    !> The side channel's format name ('mdpa'), or '' when invalid.
    function format_info_format(self) result(name)
        class(mio_format_info), intent(in) :: self
        character(:), allocatable :: name
        character(kind=c_char) :: buf(64)
        integer(c_int64_t) :: n
        name = ''
        if (.not. c_associated(self%handle)) return
        n = c_mio_format_info_format(self%handle, buf, 64_c_int64_t)
        if (n > 0) name = from_c_buf(buf, int(min(n, 63_c_int64_t)))
    end function

    !> Number of items in an MDPA section (a MIO_MDPA_* value); -1 on error.
    !> Number of Gmsh blocks/links in a MIO_GMSH_* section; -1 on error.
    integer function format_info_gmsh_count(self, section)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section
        format_info_gmsh_count = int(c_mio_gmsh_info_count(self%handle, int(section, c_int32_t)))
    end function

    !> Signed bounding-entity tags, or periodic (dimension, slave tag, master tag).
    !> Item index is 1-based; tags are raw file ids, not point indices.
    function format_info_gmsh_tags(self, section, index, stat, errmsg) result(tags)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int32), allocatable :: tags(:)
        type(c_ptr) :: p
        integer(c_int) :: dt, st
        integer(c_int32_t) :: nd
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        integer(c_int32_t), pointer :: vals(:)
        allocate (tags(0))
        st = c_mio_gmsh_info_array(self%handle, int(section, c_int32_t), &
                                  int(index - 1, c_int64_t), 0_c_int32_t, p, dt, nd, shp)
        call handle_status(st, 'gmsh_tags', stat, errmsg)
        if (st /= 0_c_int) return
        if (shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(1)])
        tags = int(vals, int32)
    end function

    !> A periodic link's affine coefficients (empty or 16 doubles), copied.
    function format_info_gmsh_affine(self, index, stat, errmsg) result(affine)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: affine(:)
        type(c_ptr) :: p
        integer(c_int) :: dt, st
        integer(c_int32_t) :: nd
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        real(c_double), pointer :: vals(:)
        allocate (affine(0))
        st = c_mio_gmsh_info_array(self%handle, MIO_GMSH_PERIODIC, &
                                  int(index - 1, c_int64_t), 1_c_int32_t, p, dt, nd, shp)
        call handle_status(st, 'gmsh_affine', stat, errmsg)
        if (st /= 0_c_int) return
        if (shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(1)])
        affine = real(vals, real64)
    end function

    !> Periodic slave/master pairs (2,N), copied as 1-based point rows.
    !> Ordering/duplicates survive I/O; mesh operations do not remap this info.
    function format_info_gmsh_pairs(self, index, stat, errmsg) result(pairs)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int64), allocatable :: pairs(:, :)
        type(c_ptr) :: p
        integer(c_int) :: dt, st
        integer(c_int32_t) :: nd
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        integer(c_int64_t), pointer :: vals(:, :)
        allocate (pairs(2, 0))
        st = c_mio_gmsh_info_array(self%handle, MIO_GMSH_PERIODIC, &
                                  int(index - 1, c_int64_t), 2_c_int32_t, p, dt, nd, shp)
        call handle_status(st, 'gmsh_pairs', stat, errmsg)
        if (st /= 0_c_int) return
        if (shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(2), shp(1)])
        pairs = int(vals, int64) + 1
    end function

    !> Number of items in a MIO_NURBS_* section; -1 on error (not an 'mfem'
    !> handle). A handle of a mesh that is not NURBS has no patches.
    integer function format_info_nurbs_count(self, section)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section
        format_info_nurbs_count = int(c_mio_nurbs_info_count(self%handle, int(section, c_int32_t)))
    end function

    !> Integer field `field` (0-based, see meshioplusplus.h) of item `index`
    !> (1-based): a patch's parametric (0) or space (1) dimension, the
    !> topology's dimension (0) or vertex count (1).
    function format_info_nurbs_int(self, section, index, field, stat, errmsg) result(v)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, field
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int64) :: v
        integer(c_int64_t) :: cv
        cv = 0
        call handle_status(c_mio_nurbs_info_int(self%handle, int(section, c_int32_t), &
                                                int(index - 1, c_int64_t), int(field, c_int32_t), cv), &
                           'nurbs_int', stat, errmsg)
        v = int(cv, int64)
    end function

    !> One borrowed array of the NURBS side channel, as a C pointer and shape.
    subroutine nurbs_array(self, section, index, field, what, p, shp, ok, stat, errmsg)
        class(mio_format_info), intent(in) :: self
        integer(c_int32_t), intent(in) :: section
        integer, intent(in) :: index, field
        character(*), intent(in) :: what
        type(c_ptr), intent(out) :: p
        integer(c_int64_t), intent(out) :: shp(MIO_MAX_NDIM)
        logical, intent(out) :: ok
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int) :: dt, st
        integer(c_int32_t) :: nd
        shp = 0
        st = c_mio_nurbs_info_array(self%handle, section, int(index - 1, c_int64_t), &
                                    int(field, c_int32_t), p, dt, nd, shp)
        call handle_status(st, what, stat, errmsg)
        ok = st == 0_c_int
    end subroutine

    ! ------------------------------------------------------------------
    ! Isogeometric layer (v16.43.0; see doc/iga.md)
    ! ------------------------------------------------------------------

    !> The options an `info%iga_*` call takes, from its optional arguments:
    !> `subdivisions` and `elevate` are one value (every direction) or one per
    !> direction.
    subroutine iga_opts_from(opts, subdivisions, elevate, sample_subdivisions, trim_samples, &
                             edges, trim)
        type(mio_iga_opts_t), intent(out) :: opts
        integer, intent(in), optional :: subdivisions(:), elevate(:)
        integer, intent(in), optional :: sample_subdivisions, trim_samples
        logical, intent(in), optional :: edges, trim
        integer :: d
        call c_mio_iga_opts_init(opts)
        if (present(subdivisions)) then
            do d = 1, 3
                opts%subdivisions(d) = int(subdivisions(min(d, size(subdivisions))), c_int32_t)
            end do
        end if
        if (present(elevate)) then
            do d = 1, 3
                opts%elevate(d) = int(elevate(min(d, size(elevate))), c_int32_t)
            end do
        end if
        if (present(sample_subdivisions)) opts%sample_subdivisions = int(sample_subdivisions, c_int32_t)
        if (present(trim_samples)) opts%trim_samples = int(trim_samples, c_int32_t)
        if (present(edges)) opts%edges = merge(1_c_int32_t, 0_c_int32_t, edges)
        if (present(trim)) opts%trim = merge(1_c_int32_t, 0_c_int32_t, trim)
    end subroutine

    !> Refine every patch (every face surface of a B-rep) of this NURBS model
    !> exactly: degree elevation by `elevate`, then knot insertion splitting each
    !> knot span in `subdivisions` parts (an integer, or one per direction; an
    !> 'mfem' model refines every direction alike). The geometry, fields and
    !> trimming loops do not change. Returns the refined model sampled
    !> (`sample_subdivisions` cells per span, 0 automatic) and, in `refined`,
    !> the refined model for `m%write_with_info`. `info` must come from
    !> `m%read_with_info` of an 'mfem' NURBS mesh, 'cad_json', 'step' or 'iges'.
    function format_info_iga_refine(self, subdivisions, elevate, sample_subdivisions, &
                                    trim_samples, edges, trim, refined, stat, errmsg) result(out)
        class(mio_format_info), intent(in) :: self
        integer, intent(in), optional :: subdivisions(:), elevate(:)
        integer, intent(in), optional :: sample_subdivisions, trim_samples
        logical, intent(in), optional :: edges, trim
        class(mio_format_info), intent(inout), optional :: refined
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_iga_opts_t) :: opts
        type(c_ptr) :: res, ih
        call iga_opts_from(opts, subdivisions, elevate, sample_subdivisions, trim_samples, &
                           edges, trim)
        res = c_mio_iga_refine(self%handle, opts, ih)
        if (.not. c_associated(res)) then
            call handle_failure('iga_refine', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(refined)) then
            call format_info_free(refined)
            refined%handle = ih
        else
            call c_mio_format_info_free(ih)
        end if
        call clear_status(stat, errmsg)
    end function

    !> Sample this NURBS model into Lagrange cells: `sample_subdivisions` cells
    !> per knot span and direction (0 automatic). `mapped` receives the model
    !> with the sampling's cell-to-patch maps.
    function format_info_iga_sample(self, sample_subdivisions, trim_samples, edges, trim, &
                                    mapped, stat, errmsg) result(out)
        class(mio_format_info), intent(in) :: self
        integer, intent(in), optional :: sample_subdivisions, trim_samples
        logical, intent(in), optional :: edges, trim
        class(mio_format_info), intent(inout), optional :: mapped
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_iga_opts_t) :: opts
        type(c_ptr) :: res, ih
        call iga_opts_from(opts, sample_subdivisions=sample_subdivisions, &
                           trim_samples=trim_samples, edges=edges, trim=trim)
        res = c_mio_iga_sample(self%handle, opts, ih)
        if (.not. c_associated(res)) then
            call handle_failure('iga_sample', mio_error_message(), stat, errmsg)
            return
        end if
        out%handle = res
        if (present(mapped)) then
            call format_info_free(mapped)
            mapped%handle = ih
        else
            call c_mio_format_info_free(ih)
        end if
        call clear_status(stat, errmsg)
    end function

    !> Triangulate the faces of this B-rep (`info` from `m%read_with_info` of
    !> 'cad_json', 'step' or 'iges') into one conforming surface mesh with edges
    !> of about `target_size` (doc/mesh_cad.md): every edge sampled once and
    !> shared by its faces, each face a constrained Delaunay triangulation of
    !> its parameters. A face that cannot be meshed is left out and recorded
    !> (`mio_mesh_cad_face_status` / `mio_mesh_cad_face_reason`), never meshed
    !> wrong; with `strict` the call fails instead. The summary comes back
    !> through the optional `report`, the warnings one per line in `warnings`.
    !> `min_angle` (degrees, under 60) refines every face until its triangles'
    !> smallest angles on the surface reach it, corners under 60 degrees and
    !> poles reported (`mio_mesh_cad_get_corner`); `smooth_iterations` passes
    !> of guarded smoothing (`smooth_method`, `MIO_MESH_CAD_SMOOTH_*`) follow.
    !> `order` 2 or 3 writes curved triangle6 or triangle10 elements, their
    !> nodes on the surface and edges (`report%num_curved_straightened`).
    !> `curvature_deviation` and `cells_across_gaps` size the mesh by the
    !> surface's curvature and its local feature size, `grading` limiting the
    !> size's growth; `pole_rings` grades the poles (v17.18.0).
    !> Sizes given per face (`face_size_faces`, 1-based, with `face_sizes`), per
    !> edge (`edge_size_edges`, `edge_sizes`) or by a `background` mesh's
    !> one-component point data `background_array` (default 'size') join the
    !> same field; `metric` (`MIO_MESH_CAD_METRIC_*`) and `virtual_edges` as in
    !> doc/mesh_cad.md; `mio_mesh_cad_face_size_limits` reads which edges
    !> limited a face given a size (v17.19.0). `reparam_condition_number`
    !> flattens a face whose metric's condition number passes it,
    !> `model_tolerances` accepts the gaps the model declares and
    !> `virtual_faces` meshes slivers with a neighbour and `lop_corners` lops
    !> acute corners in place of concentric shells (v17.20.0; read back by
    !> `mio_mesh_cad_get_reparametrised`, `_merged_face`, `_tolerance_use`,
    !> `_lopped_corner`).
    function format_info_mesh_cad_surface(self, target_size, max_chordal_deviation, &
                                          min_points_per_edge, tolerance, record_face_ids, &
                                          edges, match_unreferenced_trims, strict, max_points, &
                                          report, warnings, stat, errmsg, min_angle, &
                                          smooth_iterations, smooth_method, order, grading, &
                                          curvature_deviation, cells_across_gaps, pole_rings, &
                                          face_size_faces, face_sizes, edge_size_edges, &
                                          edge_sizes, background, background_array, metric, &
                                          virtual_edges, reparam_condition_number, &
                                          model_tolerances, virtual_faces, lop_corners) &
        result(out)
        class(mio_format_info), intent(in) :: self
        real(real64), intent(in) :: target_size
        real(real64), intent(in), optional :: max_chordal_deviation, tolerance, min_angle
        integer, intent(in), optional :: smooth_iterations, smooth_method, order, pole_rings
        real(real64), intent(in), optional :: grading, curvature_deviation, cells_across_gaps
        integer(int64), intent(in), optional :: face_size_faces(:), edge_size_edges(:)
        real(real64), intent(in), optional :: face_sizes(:), edge_sizes(:)
        class(mio_mesh), intent(in), optional :: background
        character(*), intent(in), optional :: background_array
        integer, intent(in), optional :: metric
        logical, intent(in), optional :: virtual_edges
        real(real64), intent(in), optional :: reparam_condition_number
        logical, intent(in), optional :: model_tolerances, virtual_faces, lop_corners
        integer(int64), intent(in), optional :: min_points_per_edge, max_points
        logical, intent(in), optional :: record_face_ids, edges, match_unreferenced_trims, strict
        type(mio_mesh_cad_surface_report), intent(out), optional :: report
        character(:), allocatable, intent(out), optional :: warnings
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        type(mio_mesh_cad_surface_opts_t) :: opts
        type(mio_mesh_cad_surface_report) :: crep
        type(mio_mesh_cad_surface_opts_ex_t) :: ex
        integer(c_int64_t), allocatable, target :: fidx(:), eidx(:)
        real(c_double), allocatable, target :: fsz(:), esz(:)
        character(kind=c_char, len=:), allocatable, target :: bname
        character(kind=c_char), allocatable :: wbuf(:)
        integer(c_int64_t) :: n

        call check_mesh_cad_layout()
        if (present(warnings)) warnings = ''
        call c_mio_mesh_cad_surface_opts_init(opts)
        opts%target_size = real(target_size, c_double)
        if (present(max_chordal_deviation)) &
            opts%max_chordal_deviation = real(max_chordal_deviation, c_double)
        if (present(tolerance)) opts%tolerance = real(tolerance, c_double)
        if (present(min_points_per_edge)) &
            opts%min_points_per_edge = int(min_points_per_edge, c_int64_t)
        if (present(max_points)) opts%max_points = int(max_points, c_int64_t)
        if (present(record_face_ids)) &
            opts%record_face_ids = merge(1_c_int32_t, 0_c_int32_t, record_face_ids)
        if (present(edges)) opts%edges = merge(1_c_int32_t, 0_c_int32_t, edges)
        if (present(match_unreferenced_trims)) &
            opts%match_unreferenced_trims = merge(1_c_int32_t, 0_c_int32_t, &
                                                  match_unreferenced_trims)
        if (present(strict)) opts%strict = merge(1_c_int32_t, 0_c_int32_t, strict)
        if (present(min_angle)) opts%min_angle_deg = real(min_angle, c_double)
        if (present(smooth_iterations)) &
            opts%smooth_iterations = int(smooth_iterations, c_int32_t)
        if (present(smooth_method)) opts%smooth_method = int(smooth_method, c_int32_t)
        if (present(order)) opts%order = int(order, c_int32_t)
        if (present(pole_rings)) opts%pole_rings = int(pole_rings, c_int32_t)
        if (present(grading)) opts%grading = real(grading, c_double)
        if (present(curvature_deviation)) &
            opts%curvature_deviation = real(curvature_deviation, c_double)
        if (present(cells_across_gaps)) opts%cells_across_gaps = real(cells_across_gaps, c_double)
        if (present(metric)) ex%metric = int(metric, c_int32_t)
        if (present(virtual_edges)) ex%virtual_edges = merge(1_c_int32_t, 0_c_int32_t, virtual_edges)
        if (present(reparam_condition_number)) &
            ex%reparam_condition_number = real(reparam_condition_number, c_double)
        if (present(model_tolerances)) &
            ex%model_tolerances = merge(1_c_int32_t, 0_c_int32_t, model_tolerances)
        if (present(virtual_faces)) ex%virtual_faces = merge(1_c_int32_t, 0_c_int32_t, virtual_faces)
        if (present(lop_corners)) ex%lop_corners = merge(1_c_int32_t, 0_c_int32_t, lop_corners)
        if (present(face_size_faces) .and. present(face_sizes)) then
            if (size(face_size_faces) /= size(face_sizes)) then
                call handle_failure('mesh_cad_surface', &
                                    'face_size_faces and face_sizes differ in length', stat, errmsg)
                return
            end if
            ! 1-based face indices, as the face status queries take them.
            fidx = int(face_size_faces - 1, c_int64_t)
            fsz = real(face_sizes, c_double)
            ex%num_face_sizes = size(fidx, kind=c_int64_t)
            if (ex%num_face_sizes > 0) then
                ex%face_size_faces = c_loc(fidx)
                ex%face_sizes = c_loc(fsz)
            end if
        end if
        if (present(edge_size_edges) .and. present(edge_sizes)) then
            if (size(edge_size_edges) /= size(edge_sizes)) then
                call handle_failure('mesh_cad_surface', &
                                    'edge_size_edges and edge_sizes differ in length', stat, errmsg)
                return
            end if
            eidx = int(edge_size_edges - 1, c_int64_t)
            esz = real(edge_sizes, c_double)
            ex%num_edge_sizes = size(eidx, kind=c_int64_t)
            if (ex%num_edge_sizes > 0) then
                ex%edge_size_edges = c_loc(eidx)
                ex%edge_sizes = c_loc(esz)
            end if
        end if
        if (present(background)) ex%background = background%handle
        if (present(background_array)) then
            bname = c_str(background_array)
            ex%background_array = c_loc(bname(1:1))
        end if

        out%handle = c_mio_mesh_cad_surface2(self%handle, opts, ex, crep)
        if (.not. c_associated(out%handle)) then
            call handle_failure('mesh_cad_surface', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(report)) report = crep
        if (present(warnings)) then
            ! String rule 5: query the length with an empty buffer, then copy.
            allocate (wbuf(1))
            n = c_mio_mesh_cad_surface_warnings(wbuf, 0_c_int64_t)
            if (n > 0 .and. n < int(huge(0) - 1, c_int64_t)) then
                deallocate (wbuf)
                allocate (wbuf(int(n) + 1))
                n = c_mio_mesh_cad_surface_warnings(wbuf, n + 1)
                if (n > 0) warnings = from_c_buf(wbuf, int(n))
            end if
        end if
        call clear_status(stat, errmsg)
    end function

    !> The edges (1-based) whose sampling this thread's last
    !> `info%mesh_cad_surface` call found finer than face `face`'s (1-based)
    !> given size allows (v17.19.0); empty for a face given no size or out of
    !> range.
    function mio_mesh_cad_face_size_limits(face) result(edges)
        integer, intent(in) :: face
        integer(int64), allocatable :: edges(:)
        integer(c_int64_t) :: buf(1), n
        integer(c_int64_t), allocatable :: all(:)
        allocate (edges(0))
        n = c_mio_mesh_cad_surface_face_size_limits(int(face - 1, c_int64_t), buf, 0_c_int64_t)
        if (n <= 0) return
        allocate (all(n))
        n = c_mio_mesh_cad_surface_face_size_limits(int(face - 1, c_int64_t), all, n)
        edges = int(all + 1, int64)
    end function

    !> How many faces this thread's last `info%mesh_cad_surface` call found
    !> over `reparam_condition_number` (v17.20.0), flattened or declined.
    integer function mio_mesh_cad_num_reparametrised() result(n)
        n = int(c_mio_mesh_cad_surface_num_reparametrised())
    end function

    !> Reparametrised face record `idx` (1-based): its face (1-based), its
    !> condition number and its flattening's, and `method` (1 conformal, 2
    !> mean-value, 0 declined); `stat` nonzero when out of range.
    subroutine mio_mesh_cad_get_reparametrised(idx, face, condition_number, &
                                               condition_number_disc, method, stat)
        integer, intent(in) :: idx
        integer(int64), intent(out) :: face
        real(real64), intent(out) :: condition_number, condition_number_disc
        integer, intent(out) :: method
        integer, intent(out), optional :: stat
        integer(c_int64_t) :: f
        real(c_double) :: c, d
        integer(c_int32_t) :: m
        integer(c_int) :: st
        st = c_mio_mesh_cad_surface_reparametrised(int(idx - 1, c_int64_t), f, c, d, m)
        face = int(f + 1, int64)
        condition_number = real(c, real64)
        condition_number_disc = real(d, real64)
        method = int(m)
        if (present(stat)) stat = int(st)
    end subroutine

    !> How many sliver faces this thread's last call found (v17.20.0).
    integer function mio_mesh_cad_num_merged_faces() result(n)
        n = int(c_mio_mesh_cad_surface_num_merged_faces())
    end function

    !> Sliver record `idx` (1-based): the sliver's and its host's faces
    !> (1-based; 0 for no host), its width and whether it was merged.
    subroutine mio_mesh_cad_get_merged_face(idx, sliver, host, width, merged, stat)
        integer, intent(in) :: idx
        integer(int64), intent(out) :: sliver, host
        real(real64), intent(out) :: width
        logical, intent(out) :: merged
        integer, intent(out), optional :: stat
        integer(c_int64_t) :: s, h
        real(c_double) :: w
        integer(c_int32_t) :: m
        integer(c_int) :: st
        st = c_mio_mesh_cad_surface_merged_face(int(idx - 1, c_int64_t), s, h, w, m)
        sliver = int(s + 1, int64)
        host = int(h + 1, int64)
        width = real(w, real64)
        merged = m /= 0
        if (present(stat)) stat = int(st)
    end subroutine

    !> How many places a declared tolerance accepted a gap (v17.20.0).
    integer function mio_mesh_cad_num_tolerance_uses() result(n)
        n = int(c_mio_mesh_cad_surface_num_tolerance_uses())
    end function

    !> Tolerance use `idx` (1-based): the edge (1-based), the tolerance it
    !> declares, the gap accepted and where (`MIO_MESH_CAD_TOLERANCE_*`).
    subroutine mio_mesh_cad_get_tolerance_use(idx, edge, declared, gap, where, stat)
        integer, intent(in) :: idx
        integer(int64), intent(out) :: edge
        real(real64), intent(out) :: declared, gap
        integer, intent(out) :: where
        integer, intent(out), optional :: stat
        integer(c_int64_t) :: e
        real(c_double) :: d, g
        integer(c_int32_t) :: w
        integer(c_int) :: st
        st = c_mio_mesh_cad_surface_tolerance_use(int(idx - 1, c_int64_t), e, d, g, w)
        edge = int(e + 1, int64)
        declared = real(d, real64)
        gap = real(g, real64)
        where = int(w)
        if (present(stat)) stat = int(st)
    end subroutine

    !> How many corners this thread's last call lopped (v17.20.0).
    integer function mio_mesh_cad_num_lopped_corners() result(n)
        n = int(c_mio_mesh_cad_surface_num_lopped_corners())
    end function

    !> Lopped corner `idx` (1-based): its face and mesh point (1-based), its
    !> wedge in the face's plane, the distance to the face's nearest other
    !> corner, its ball's radius, the face's turning, Corollary 3's bound (0
    !> on a curved face), the triangles left in its ball and their smallest
    !> angle, whether they lie within twice the radius, and whether the bound
    !> is guaranteed (degrees throughout).
    subroutine mio_mesh_cad_get_lopped_corner(idx, face, point, angle, lfs0, radius, turning, &
                                              bound, num_forced, min_forced, within_two_radii, &
                                              guaranteed, stat)
        integer, intent(in) :: idx
        integer(int64), intent(out) :: face, point, num_forced
        real(real64), intent(out) :: angle, lfs0, radius, turning, bound, min_forced
        logical, intent(out) :: within_two_radii, guaranteed
        integer, intent(out), optional :: stat
        integer(c_int64_t) :: f, p, nf
        real(c_double) :: a, l, r, t, b, m
        integer(c_int32_t) :: w, g
        integer(c_int) :: st
        st = c_mio_mesh_cad_surface_lopped_corner(int(idx - 1, c_int64_t), f, p, a, l, r, t, b, &
                                                  nf, m, w, g)
        face = int(f + 1, int64)
        point = int(p + 1, int64)
        num_forced = int(nf, int64)
        angle = real(a, real64)
        lfs0 = real(l, real64)
        radius = real(r, real64)
        turning = real(t, real64)
        bound = real(b, real64)
        min_forced = real(m, real64)
        within_two_radii = w /= 0
        guaranteed = g /= 0
        if (present(stat)) stat = int(st)
    end subroutine

    !> Why lopped corner `idx`'s (1-based) bound is not guaranteed; empty when
    !> it is or out of range.
    function mio_mesh_cad_lopped_corner_reason(idx) result(reason)
        integer, intent(in) :: idx
        character(:), allocatable :: reason
        character(kind=c_char), allocatable :: buf(:)
        integer(c_int64_t) :: n
        reason = ''
        allocate (buf(1))
        n = c_mio_mesh_cad_surface_lopped_corner_reason(int(idx - 1, c_int64_t), buf, 0_c_int64_t)
        if (n > 0 .and. n < int(huge(0) - 1, c_int64_t)) then
            deallocate (buf)
            allocate (buf(int(n) + 1))
            n = c_mio_mesh_cad_surface_lopped_corner_reason(int(idx - 1, c_int64_t), buf, n + 1)
            if (n > 0) reason = from_c_buf(buf, int(n))
        end if
    end function

    !> The status (`MIO_MESH_CAD_*`) of face `face` (1-based) in this thread's
    !> last `info%mesh_cad_surface` call, or -1 when out of range.
    integer function mio_mesh_cad_face_status(face) result(st)
        integer, intent(in) :: face
        st = int(c_mio_mesh_cad_surface_face_status(int(face - 1, c_int64_t)))
    end function

    !> Why face `face` (1-based) of this thread's last `info%mesh_cad_surface`
    !> call was left out; empty when it was meshed or is out of range.
    function mio_mesh_cad_face_reason(face) result(reason)
        integer, intent(in) :: face
        character(:), allocatable :: reason
        character(kind=c_char), allocatable :: buf(:)
        integer(c_int64_t) :: n
        reason = ''
        allocate (buf(1))
        n = c_mio_mesh_cad_surface_face_reason(int(face - 1, c_int64_t), buf, 0_c_int64_t)
        if (n > 0 .and. n < int(huge(0) - 1, c_int64_t)) then
            deallocate (buf)
            allocate (buf(int(n) + 1))
            n = c_mio_mesh_cad_surface_face_reason(int(face - 1, c_int64_t), buf, n + 1)
            if (n > 0) reason = from_c_buf(buf, int(n))
        end if
    end function

    !> The angle histogram of this thread's last `info%mesh_cad_surface` call:
    !> every angle of every triangle on the surface, in `MIO_MESH_CAD_ANGLE_BINS`
    !> bins of 10 degrees.
    function mio_mesh_cad_angle_histogram() result(bins)
        integer(int64) :: bins(MIO_MESH_CAD_ANGLE_BINS)
        integer(c_int64_t) :: cbins(MIO_MESH_CAD_ANGLE_BINS)
        integer(c_int64_t) :: n
        cbins = 0
        n = c_mio_mesh_cad_surface_angle_histogram(cbins, &
                                                   int(MIO_MESH_CAD_ANGLE_BINS, c_int64_t))
        bins = int(cbins, int64)
    end function

    !> How many corners this thread's last `info%mesh_cad_surface` call
    !> reported, every face's.
    integer function mio_mesh_cad_num_corners() result(n)
        n = int(c_mio_mesh_cad_surface_num_corners())
    end function

    !> Corner `k` (1-based) of this thread's last call, its `face` and `point`
    !> 1-based; `stat` nonzero when out of range.
    function mio_mesh_cad_get_corner(k, stat, errmsg) result(corner)
        integer, intent(in) :: k
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh_cad_corner) :: corner
        call check_mesh_cad_layout()
        if (c_mio_mesh_cad_surface_corner(int(k - 1, c_int64_t), corner) /= 0) then
            call handle_failure('mesh_cad_get_corner', mio_error_message(), stat, errmsg)
            return
        end if
        corner%face = corner%face + 1
        corner%point = corner%point + 1
        call clear_status(stat, errmsg)
    end function

    !> The layouts of `mio_mesh_cad_surface_opts_t`, `mio_mesh_cad_surface_report`
    !> and `mio_mesh_cad_corner` against the C structs (104, 200 and 120 bytes),
    !> checked once per process.
    subroutine check_mesh_cad_layout()
        logical, save :: checked = .false.
        type(mio_mesh_cad_surface_opts_t) :: opts
        type(mio_mesh_cad_surface_report) :: report
        type(mio_mesh_cad_corner) :: corner
        type(mio_mesh_cad_surface_opts_ex_t) :: ex
        if (checked) return
        checked = .true.
        if (c_sizeof(opts) /= 104_c_size_t .or. c_sizeof(report) /= 200_c_size_t .or. &
            c_sizeof(corner) /= 120_c_size_t .or. &
            (c_sizeof(c_null_ptr) == 8_c_size_t .and. c_sizeof(ex) /= 136_c_size_t)) then
            write (error_unit, '(a,3(i0,1x),a)') &
                'meshio++: mesh_cad_surface layout mismatch (', c_sizeof(opts), &
                c_sizeof(report), c_sizeof(corner), 'bytes)'
            error stop 1
        end if
    end subroutine

    !> Evaluate patch `patch` (1-based; the face of a B-rep) at the parametric
    !> coordinates `params(d, N)` (one column per point, `d` the parametric
    !> dimension), each inside the patch's knot range. `derivatives` 0 or 1.
    function format_info_iga_evaluate(self, patch, params, derivatives, trim_samples, stat, &
                                      errmsg) result(ev)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: patch
        real(real64), intent(in) :: params(:, :)
        integer, intent(in), optional :: derivatives, trim_samples
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_iga_evaluation) :: ev
        type(c_ptr) :: res
        integer(c_int32_t) :: der, ts
        real(c_double), allocatable :: flat(:)
        der = 0; if (present(derivatives)) der = int(derivatives, c_int32_t)
        ts = 16; if (present(trim_samples)) ts = int(trim_samples, c_int32_t)
        flat = reshape(real(params, c_double), [size(params)])
        res = c_mio_iga_evaluate(self%handle, int(patch - 1, c_int64_t), flat, &
                                 int(size(params, 2), c_int64_t), der, ts)
        if (.not. c_associated(res)) then
            call handle_failure('iga_evaluate', mio_error_message(), stat, errmsg)
            return
        end if
        ev%handle = res
        call clear_status(stat, errmsg)
    end function

    subroutine iga_evaluation_free(self)
        class(mio_iga_evaluation), intent(inout) :: self
        if (c_associated(self%handle)) call c_mio_iga_evaluation_free(self%handle)
        self%handle = c_null_ptr
    end subroutine

    logical function iga_evaluation_is_valid(self)
        class(mio_iga_evaluation), intent(in) :: self
        iga_evaluation_is_valid = c_associated(self%handle)
    end function

    !> The number of NURBS fields the evaluation carries.
    integer function iga_evaluation_num_fields(self)
        class(mio_iga_evaluation), intent(in) :: self
        iga_evaluation_num_fields = int(c_mio_iga_evaluation_num_fields(self%handle))
    end function

    subroutine iga_array(self, which, index, what, p, shp, ok, stat, errmsg)
        class(mio_iga_evaluation), intent(in) :: self
        integer(c_int32_t), intent(in) :: which
        integer, intent(in) :: index
        character(*), intent(in) :: what
        type(c_ptr), intent(out) :: p
        integer(c_int64_t), intent(out) :: shp(MIO_MAX_NDIM)
        logical, intent(out) :: ok
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int) :: dt, st
        integer(c_int32_t) :: nd
        shp = 0
        st = c_mio_iga_evaluation_array(self%handle, which, int(index, c_int64_t), p, dt, nd, shp)
        call handle_status(st, what, stat, errmsg)
        ok = st == 0_c_int
    end subroutine

    !> The points, `(space dimension, N)`.
    function iga_evaluation_points(self, stat, errmsg) result(xyz)
        class(mio_iga_evaluation), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: xyz(:, :)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        real(c_double), pointer :: vals(:, :)
        allocate (xyz(0, 0))
        call iga_array(self, MIO_IGA_POINTS, 0, 'iga_points', p, shp, ok, stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(2), shp(1)])
        xyz = real(vals, real64)
    end function

    !> The first derivatives, `(space dimension, parametric dimension, N)`
    !> (`derivatives = 1` only).
    function iga_evaluation_derivatives(self, stat, errmsg) result(d)
        class(mio_iga_evaluation), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: d(:, :, :)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        real(c_double), pointer :: vals(:, :, :)
        allocate (d(0, 0, 0))
        call iga_array(self, MIO_IGA_DERIVATIVES, 0, 'iga_derivatives', p, shp, ok, stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(3), shp(2), shp(1)])
        d = real(vals, real64)
    end function

    !> The unit normals of a surface in 3-D, `(3, N)` (`derivatives = 1` only).
    function iga_evaluation_normals(self, stat, errmsg) result(nrm)
        class(mio_iga_evaluation), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: nrm(:, :)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        real(c_double), pointer :: vals(:, :)
        allocate (nrm(0, 0))
        call iga_array(self, MIO_IGA_NORMALS, 0, 'iga_normals', p, shp, ok, stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(2), shp(1)])
        nrm = real(vals, real64)
    end function

    !> Whether each point of a B-rep face lies inside its trimming loops.
    function iga_evaluation_inside(self, stat, errmsg) result(inside)
        class(mio_iga_evaluation), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        logical, allocatable :: inside(:)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        integer(c_int8_t), pointer :: vals(:)
        allocate (inside(0))
        call iga_array(self, MIO_IGA_INSIDE, 0, 'iga_inside', p, shp, ok, stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(1)])
        inside = vals /= 0_c_int8_t
    end function

    !> NURBS field `k` (1-based) at the points, `(components, N)`.
    function iga_evaluation_field(self, k, stat, errmsg) result(values)
        class(mio_iga_evaluation), intent(in) :: self
        integer, intent(in) :: k
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: values(:, :)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        real(c_double), pointer :: vals(:, :)
        allocate (values(0, 0))
        call iga_array(self, MIO_IGA_FIELD, k - 1, 'iga_field', p, shp, ok, stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(2), shp(1)])
        values = real(vals, real64)
    end function

    !> A patch's degree per parametric direction (patch index 1-based).
    function format_info_nurbs_degrees(self, patch, stat, errmsg) result(degrees)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: patch
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int32), allocatable :: degrees(:)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        integer(c_int32_t), pointer :: vals(:)
        allocate (degrees(0))
        call nurbs_array(self, MIO_NURBS_PATCHES, patch, 0, 'nurbs_degrees', p, shp, ok, stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(1)])
        degrees = int(vals, int32)
    end function

    !> A patch's knot vector along direction `dir` (both 1-based), copied.
    function format_info_nurbs_knots(self, patch, dir, stat, errmsg) result(knots)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: patch, dir
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: knots(:)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        real(c_double), pointer :: vals(:)
        allocate (knots(0))
        call nurbs_array(self, MIO_NURBS_PATCHES, patch, 2 + dir, 'nurbs_knots', p, shp, ok, &
                         stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(1)])
        knots = real(vals, real64)
    end function

    !> A patch's control net (space dimension, points), direction 0 fastest.
    function format_info_nurbs_control_points(self, patch, stat, errmsg) result(xyz)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: patch
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: xyz(:, :)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        real(c_double), pointer :: vals(:, :)
        allocate (xyz(0, 0))
        call nurbs_array(self, MIO_NURBS_PATCHES, patch, 1, 'nurbs_control_points', p, shp, ok, &
                         stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(2), shp(1)])
        xyz = real(vals, real64)
    end function

    !> A patch's control-point weights.
    function format_info_nurbs_weights(self, patch, stat, errmsg) result(weights)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: patch
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: weights(:)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        real(c_double), pointer :: vals(:)
        allocate (weights(0))
        call nurbs_array(self, MIO_NURBS_PATCHES, patch, 2, 'nurbs_weights', p, shp, ok, stat, &
                         errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(1)])
        weights = real(vals, real64)
    end function

    !> The control values (components, points) of field `k` on a patch (1-based).
    function format_info_nurbs_field(self, patch, k, stat, errmsg) result(values)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: patch, k
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: values(:, :)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        real(c_double), pointer :: vals(:, :)
        allocate (values(0, 0))
        call nurbs_array(self, MIO_NURBS_PATCHES, patch, 5 + k, 'nurbs_field', p, shp, ok, &
                         stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(2), shp(1)])
        values = real(vals, real64)
    end function

    !> The point-data name of field `k` (1-based); '' on error.
    function format_info_nurbs_field_name(self, k) result(name)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: k
        character(:), allocatable :: name
        character(kind=c_char), allocatable :: buf(:)
        integer(c_int64_t) :: n
        name = ''
        allocate (buf(1))
        n = c_mio_nurbs_info_string(self%handle, MIO_NURBS_FIELDS, int(k - 1, c_int64_t), &
                                    0_c_int32_t, buf, 0_c_int64_t)
        if (n <= 0) return
        deallocate (buf)
        allocate (buf(n + 1))
        n = c_mio_nurbs_info_string(self%handle, MIO_NURBS_FIELDS, int(k - 1, c_int64_t), &
                                    0_c_int32_t, buf, n + 1)
        name = from_c_buf(buf, int(n))
    end function

    !> An index array of the TOPOLOGY or CELL_MAP section, (columns, rows),
    !> a one-dimensional one as (1, n). Vertices, patches and knot spans come
    !> back 1-based; attributes and the edges' knot-vector ids as the file
    !> spelled them.
    function format_info_nurbs_ids(self, section, field, stat, errmsg) result(ids)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, field
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int64), allocatable :: ids(:, :)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        integer(c_int64_t), pointer :: vals(:, :)
        integer(c_int64_t) :: cols
        allocate (ids(0, 0))
        call nurbs_array(self, int(section, c_int32_t), 1, field, 'nurbs_ids', p, shp, ok, stat, &
                         errmsg)
        if (.not. ok .or. shp(1) == 0) return
        cols = max(shp(2), 1_c_int64_t)
        call c_f_pointer(p, vals, [cols, shp(1)])
        ids = int(vals, int64)
        if (section == MIO_NURBS_TOPOLOGY .and. (field == 1 .or. field == 3)) return
        if (section == MIO_NURBS_TOPOLOGY .and. field == 4) then
            ids(2:3, :) = ids(2:3, :) + 1
            return
        end if
        ids = ids + 1
    end function

    !> Number of items in a MIO_BREP_* section; -1 on error (not a B-rep handle).
    integer function format_info_brep_count(self, section)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section
        format_info_brep_count = int(c_mio_brep_info_count(self%handle, int(section, c_int32_t)))
    end function

    !> Whether int field `field` of a B-rep section is an index into another
    !> section (returned 1-based here, 0 for none).
    logical function brep_int_is_index(section, field)
        integer, intent(in) :: section, field
        select case (section)
        case (MIO_BREP_FACES)
            brep_int_is_index = field == 1 .or. field == 6 .or. field == 8
        case (MIO_BREP_LOOPS)
            brep_int_is_index = field == 0 .or. field == 1
        case (MIO_BREP_TRIMS)
            brep_int_is_index = field == 0 .or. field == 1
        case (MIO_BREP_EDGES)
            brep_int_is_index = field == 2
        case (MIO_BREP_VERTICES)
            brep_int_is_index = field == 1
        case default
            brep_int_is_index = .false.
        end select
    end function

    !> Integer field `field` (0-based, see meshioplusplus.h) of item `index`
    !> (1-based). Fields that index another section (a face's patch, first
    !> loop or trim; a loop's or trim's face; a trim's loop; an edge's or
    !> vertex's first reference) come back 1-based, 0 for none; ids, flags,
    !> counts, degrees and trim indices as the file spells them.
    function format_info_brep_int(self, section, index, field, stat, errmsg) result(v)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, field
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int64) :: v
        integer(c_int64_t) :: cv
        cv = 0
        call handle_status(c_mio_brep_info_int(self%handle, int(section, c_int32_t), &
                                               int(index - 1, c_int64_t), int(field, c_int32_t), cv), &
                           'brep_int', stat, errmsg)
        v = int(cv, int64)
        if (brep_int_is_index(section, field)) v = v + 1
    end function

    !> One borrowed array of the B-rep side channel, as a C pointer and shape.
    subroutine brep_array(self, section, index, field, what, p, shp, ok, stat, errmsg)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, field
        character(*), intent(in) :: what
        type(c_ptr), intent(out) :: p
        integer(c_int64_t), intent(out) :: shp(MIO_MAX_NDIM)
        logical, intent(out) :: ok
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int) :: dt, st
        integer(c_int32_t) :: nd
        shp = 0
        st = c_mio_brep_info_array(self%handle, int(section, c_int32_t), int(index - 1, c_int64_t), &
                                   int(field, c_int32_t), p, dt, nd, shp)
        call handle_status(st, what, stat, errmsg)
        ok = st == 0_c_int
    end subroutine

    !> A Float64 vector field (weights, knots, an active range, local
    !> coordinates) of item `index` (1-based), copied.
    function format_info_brep_values(self, section, index, field, stat, errmsg) result(values)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, field
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: values(:)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        real(c_double), pointer :: vals(:)
        allocate (values(0))
        call brep_array(self, section, index, field, 'brep_values', p, shp, ok, stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(1)])
        values = real(vals, real64)
    end function

    !> A control net (3, points) of item `index` (1-based): a patch's (field
    !> 0) or a curve's (field 10), copied.
    function format_info_brep_points(self, section, index, field, stat, errmsg) result(xyz)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, field
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: xyz(:, :)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        real(c_double), pointer :: vals(:, :)
        allocate (xyz(0, 0))
        call brep_array(self, section, index, field, 'brep_points', p, shp, ok, stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(2), shp(1)])
        xyz = real(vals, real64)
    end function

    !> An Int64 field of item `index` (1-based): control-point ids as the file
    !> spells them, or a brep's face/edge/vertex indices (fields 0/1/2 of
    !> MIO_BREP_BREPS), 1-based.
    function format_info_brep_ids(self, section, index, field, stat, errmsg) result(ids)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, field
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int64), allocatable :: ids(:)
        type(c_ptr) :: p
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        logical :: ok
        integer(c_int64_t), pointer :: vals(:)
        allocate (ids(0))
        call brep_array(self, section, index, field, 'brep_ids', p, shp, ok, stat, errmsg)
        if (.not. ok .or. shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(1)])
        ids = int(vals, int64)
        if (section == MIO_BREP_BREPS) ids = ids + 1
    end function

    !> String field `field` of item `index` (1-based); '' on error.
    function format_info_brep_string(self, section, index, field) result(text)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, field
        character(:), allocatable :: text
        character(kind=c_char), allocatable :: buf(:)
        integer(c_int64_t) :: n
        text = ''
        allocate (buf(1))
        n = c_mio_brep_info_string(self%handle, int(section, c_int32_t), int(index - 1, c_int64_t), &
                                   int(field, c_int32_t), buf, 0_c_int64_t)
        if (n <= 0) return
        deallocate (buf)
        allocate (buf(n + 1))
        n = c_mio_brep_info_string(self%handle, int(section, c_int32_t), int(index - 1, c_int64_t), &
                                   int(field, c_int32_t), buf, n + 1)
        text = from_c_buf(buf, int(n))
    end function

    integer function format_info_mdpa_count(self, section)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section
        format_info_mdpa_count = int(c_mio_mdpa_info_count(self%handle, int(section, c_int32_t)))
    end function

    !> String field `field` (0-based, per section: see meshioplusplus.h) of
    !> item `index` (1-based). '' on error.
    function format_info_mdpa_string(self, section, index, field) result(str)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, field
        character(:), allocatable :: str
        character(kind=c_char), allocatable :: buf(:)
        integer(c_int64_t) :: n
        str = ''
        allocate (buf(1))
        n = c_mio_mdpa_info_string(self%handle, int(section, c_int32_t), &
                                   int(index - 1, c_int64_t), int(field, c_int32_t), buf, &
                                   0_c_int64_t)
        if (n <= 0) return
        deallocate (buf)
        allocate (buf(n + 1))
        n = c_mio_mdpa_info_string(self%handle, int(section, c_int32_t), &
                                   int(index - 1, c_int64_t), int(field, c_int32_t), buf, n + 1)
        str = from_c_buf(buf, int(n))
    end function

    !> Integer field `field` of item `index` (1-based): a Properties or Mesh
    !> id, or 1 for a Conditions entity-name entry.
    function format_info_mdpa_int(self, section, index, field, stat, errmsg) result(v)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, field
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int64) :: v
        integer(c_int64_t) :: cv
        cv = 0
        call handle_status(c_mio_mdpa_info_int(self%handle, int(section, c_int32_t), &
                                               int(index - 1, c_int64_t), int(field, c_int32_t), &
                                               cv), 'mdpa_int', stat, errmsg)
        v = int(cv, int64)
    end function

    !> An id list of item `index` (1-based): geometry ids (GEOMETRIES field 1),
    !> Mesh-block nodes (MESH_BLOCKS field 0, 1-based points) / element ids (1)
    !> / condition ids (2), or table ids (SUBMODELPARTS field 0).
    !> SUBMODELPARTS fields: 0 tables, 1 geometry ids, 2 constraint ids (raw file ids).
    function format_info_mdpa_ids(self, section, index, field, stat, errmsg) result(ids)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, field
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int64), allocatable :: ids(:)
        type(c_ptr) :: p
        integer(c_int) :: dt
        integer(c_int32_t) :: nd
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        integer(c_int64_t), pointer :: vals(:)
        integer(c_int) :: st
        allocate (ids(0))
        st = c_mio_mdpa_info_array(self%handle, int(section, c_int32_t), &
                                   int(index - 1, c_int64_t), int(field, c_int32_t), p, dt, nd, shp)
        if (st == 0_c_int .and. (nd /= 1 .or. dt /= MIO_INT64)) then
            call handle_failure('mdpa_ids', 'not an id list (use mdpa_table/mdpa_geometry)', &
                                stat, errmsg)
            return
        end if
        call handle_status(st, 'mdpa_ids', stat, errmsg)
        if (st /= 0_c_int) return
        if (shp(1) == 0) return
        call c_f_pointer(p, vals, [shp(1)])
        ids = int(vals, int64)
        if (section == MIO_MDPA_MESH_BLOCKS .and. field == 0) ids = ids + 1
    end function

    !> A top-level table's values, `(rows, columns)`, of item `index` (1-based).
    function format_info_mdpa_table(self, index, stat, errmsg) result(t)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: t(:, :)
        type(c_ptr) :: p
        integer(c_int) :: dt, st
        integer(c_int32_t) :: nd
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        real(c_double), pointer :: vals(:, :)
        allocate (t(0, 0))
        st = c_mio_mdpa_info_array(self%handle, MIO_MDPA_TABLES, int(index - 1, c_int64_t), &
                                   0_c_int32_t, p, dt, nd, shp)
        call handle_status(st, 'mdpa_table', stat, errmsg)
        if (st /= 0_c_int) return
        if (nd /= 2 .or. shp(1) * shp(2) == 0) return
        ! C is row-major (rows, cols): the same bytes are Fortran (cols, rows).
        call c_f_pointer(p, vals, [shp(2), shp(1)])
        t = transpose(vals)
    end function

    !> A geometry run's connectivity, `conn(nodes_per_geometry, num_geometries)`,
    !> as 1-based points -- the same layout as a cell block's.
    function format_info_mdpa_geometry(self, index, stat, errmsg) result(conn)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(int64), allocatable :: conn(:, :)
        type(c_ptr) :: p
        integer(c_int) :: dt, st
        integer(c_int32_t) :: nd
        integer(c_int64_t) :: shp(MIO_MAX_NDIM)
        integer(c_int64_t), pointer :: vals(:, :)
        allocate (conn(0, 0))
        st = c_mio_mdpa_info_array(self%handle, MIO_MDPA_GEOMETRIES, int(index - 1, c_int64_t), &
                                   0_c_int32_t, p, dt, nd, shp)
        call handle_status(st, 'mdpa_geometry', stat, errmsg)
        if (st /= 0_c_int) return
        if (nd /= 2 .or. shp(1) * shp(2) == 0) return
        call c_f_pointer(p, vals, [shp(2), shp(1)])
        conn = int(vals, int64) + 1
    end function

    !> Number of key/value entries of item `index` (1-based) of PROPERTIES,
    !> MODEL_PART_DATA, MESH_BLOCKS or SUBMODELPARTS; -1 on error.
    integer function format_info_mdpa_data_count(self, section, index)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index
        format_info_mdpa_data_count = int(c_mio_mdpa_info_data_count( &
                                          self%handle, int(section, c_int32_t), &
                                          int(index - 1, c_int64_t)))
    end function

    !> A MIO_MDPA_VALUE_* kind for entry `entry` (1-based); -1 on error.
    integer function format_info_mdpa_data_kind(self, section, index, entry)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, entry
        format_info_mdpa_data_kind = int(c_mio_mdpa_info_data_kind( &
                                         self%handle, int(section, c_int32_t), &
                                         int(index - 1, c_int64_t), int(entry - 1, c_int64_t)))
    end function

    !> Entry `entry`'s key (a table's header arguments for a TABLE entry).
    function format_info_mdpa_data_key(self, section, index, entry) result(str)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, entry
        character(:), allocatable :: str
        str = format_info_data_string(self, section, index, entry, 0)
    end function

    !> Entry `entry`'s text value ('' unless its kind is MIO_MDPA_VALUE_TEXT).
    function format_info_mdpa_data_text(self, section, index, entry) result(str)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, entry
        character(:), allocatable :: str
        str = format_info_data_string(self, section, index, entry, 1)
    end function

    function format_info_data_string(self, section, index, entry, field) result(str)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, entry, field
        character(:), allocatable :: str
        character(kind=c_char), allocatable :: buf(:)
        integer(c_int64_t) :: n
        str = ''
        allocate (buf(1))
        n = c_mio_mdpa_info_data_string(self%handle, int(section, c_int32_t), &
                                        int(index - 1, c_int64_t), int(entry - 1, c_int64_t), &
                                        int(field, c_int32_t), buf, 0_c_int64_t)
        if (n <= 0) return
        deallocate (buf)
        allocate (buf(n + 1))
        n = c_mio_mdpa_info_data_string(self%handle, int(section, c_int32_t), &
                                        int(index - 1, c_int64_t), int(entry - 1, c_int64_t), &
                                        int(field, c_int32_t), buf, n + 1)
        str = from_c_buf(buf, int(n))
    end function

    !> Entry `entry`'s numbers, flattened row-major: one value for a number,
    !> rows*columns for a table (see mdpa_data_count/kind), none for text.
    function format_info_mdpa_data_values(self, section, index, entry, stat, errmsg) result(v)
        class(mio_format_info), intent(in) :: self
        integer, intent(in) :: section, index, entry
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(real64), allocatable :: v(:)
        type(c_ptr) :: p
        integer(c_int) :: dt, st
        integer(c_int32_t) :: nd, k
        integer(c_int64_t) :: shp(MIO_MAX_NDIM), total
        real(c_double), pointer :: vals(:)
        allocate (v(0))
        st = c_mio_mdpa_info_data_array(self%handle, int(section, c_int32_t), &
                                        int(index - 1, c_int64_t), int(entry - 1, c_int64_t), &
                                        p, dt, nd, shp)
        call handle_status(st, 'mdpa_data_values', stat, errmsg)
        if (st /= 0_c_int) return
        total = 1
        do k = 1, nd
            total = total * shp(k)
        end do
        if (total == 0 .or. .not. c_associated(p)) return
        call c_f_pointer(p, vals, [total])
        v = real(vals, real64)
    end function

    !> `0xRRGGBBAA` from four channel values in 0..255, for the colour fields of
    !> `mio_render_settings`.
    pure function mio_rgba(r, g, b, a) result(c)
        integer, intent(in) :: r, g, b
        integer, intent(in), optional :: a
        integer(int64) :: c
        integer :: alpha
        alpha = 255
        if (present(a)) alpha = a
        c = int(r, int64) * 16777216_int64 + int(g, int64) * 65536_int64 &
            + int(b, int64) * 256_int64 + int(alpha, int64)
    end function

    !> The bit pattern of a `0xRRGGBBAA` value as the C `uint32_t` field.
    pure function rgba_to_c(c) result(i)
        integer(int64), intent(in) :: c
        integer(c_int32_t) :: i
        if (c >= 2147483648_int64) then
            i = int(c - 4294967296_int64, c_int32_t)
        else
            i = int(c, c_int32_t)
        end if
    end function

    !> A Fortran string as a NUL-terminated C buffer, or NULL when unset/empty.
    subroutine set_cstr(src, buf, p)
        character(:), allocatable, intent(in) :: src
        character(kind=c_char, len=STRBUF_LEN), target, intent(out) :: buf
        type(c_ptr), intent(out) :: p
        p = c_null_ptr
        buf = c_null_char
        if (.not. allocated(src)) return
        if (len_trim(src) == 0) return
        if (len_trim(src) >= STRBUF_LEN) return
        buf = trim(src)//c_null_char
        p = c_loc(buf(1:1))
    end subroutine

    !> Fill the C option struct from the settings. The string buffers are owned by
    !> the caller so they outlive the C call.
    subroutine render_opts_to_c(s, o, b, levels)
        type(mio_render_settings), intent(in) :: s
        type(mio_render_opts_t), intent(inout) :: o
        character(kind=c_char, len=STRBUF_LEN), target, intent(out) :: b(9)
        real(c_double), target, intent(inout) :: levels(:)
        real(c_double) :: cuts(12)

        ! `levels` is the caller's scratch: the contour levels fill all but its last
        ! twelve slots, which carry the cut-away planes (two of six numbers).
        call c_mio_render_opts_init(o)
        call set_cstr(s%view, b(1), o%view)
        call set_cstr(s%color_by, b(2), o%color_by)
        call set_cstr(s%cmap, b(3), o%cmap)
        call set_cstr(s%reduce, b(4), o%reduce)
        call set_cstr(s%expr, b(5), o%expr)
        call set_cstr(s%vectors, b(6), o%vectors)
        call set_cstr(s%warp, b(7), o%warp)
        call set_cstr(s%quality_metric, b(8), o%quality_metric)
        call set_cstr(s%streamlines, b(9), o%streamlines)
        o%stream_seeds = int(s%stream_seeds, c_int32_t)
        o%stream_length = real(s%stream_length, c_double)
        o%stream_color = rgba_to_c(s%stream_color)
        if (allocated(s%iso_levels)) then
            if (size(s%iso_levels) > 0) then
                levels(1:size(s%iso_levels)) = real(s%iso_levels, c_double)
                o%iso_levels = c_loc(levels(1))
                o%num_iso_levels = int(size(s%iso_levels), c_int64_t)
            end if
        end if
        o%pixel_aspect = real(s%pixel_aspect, c_double)
        o%azimuth = real(s%azimuth, c_double)
        o%elevation = real(s%elevation, c_double)
        o%roll = real(s%roll, c_double)
        o%fov_deg = real(s%fov_deg, c_double)
        o%zoom = real(s%zoom, c_double)
        o%pan_x = real(s%pan_x, c_double)
        o%pan_y = real(s%pan_y, c_double)
        o%ambient = real(s%ambient, c_double)
        o%light_x = real(s%light_x, c_double)
        o%light_y = real(s%light_y, c_double)
        o%light_z = real(s%light_z, c_double)
        o%split_angle = real(s%split_angle, c_double)
        o%feature_angle = real(s%feature_angle, c_double)
        o%point_radius = real(s%point_radius, c_double)
        o%scale_threshold = real(s%scale_threshold, c_double)
        o%vector_length = real(s%vector_length, c_double)
        o%warp_scale = real(s%warp_scale, c_double)
        o%width = int(s%width, c_int32_t)
        o%height = int(s%height, c_int32_t)
        o%supersample = int(s%supersample, c_int32_t)
        o%shading = int(s%shading, c_int32_t)
        o%edges = int(s%edges, c_int32_t)
        o%scale = int(s%scale, c_int32_t)
        o%diagnostic = int(s%diagnostic, c_int32_t)
        o%isolines = int(s%isolines, c_int32_t)
        o%vector_count = int(s%vector_count, c_int32_t)
        o%component = int(s%component, c_int32_t)
        o%perspective = merge(1_c_int32_t, 0_c_int32_t, s%perspective)
        o%two_sided = merge(1_c_int32_t, 0_c_int32_t, s%two_sided)
        o%colorbar = merge(1_c_int32_t, 0_c_int32_t, s%colorbar)
        o%axes = merge(1_c_int32_t, 0_c_int32_t, s%axes)
        o%scale_bar = merge(1_c_int32_t, 0_c_int32_t, s%scale_bar)
        o%symmetric = merge(1_c_int32_t, 0_c_int32_t, s%symmetric)
        o%categorical = merge(1_c_int32_t, 0_c_int32_t, s%categorical)
        o%color_regions = merge(1_c_int32_t, 0_c_int32_t, s%color_regions)
        o%category_edges = merge(1_c_int32_t, 0_c_int32_t, s%category_edges)
        o%warp_outline = merge(1_c_int32_t, 0_c_int32_t, s%warp_outline)
        if (allocated(s%vmin)) then
            o%vmin = real(s%vmin, c_double)
            o%has_vmin = 1
        end if
        if (allocated(s%vmax)) then
            o%vmax = real(s%vmax, c_double)
            o%has_vmax = 1
        end if
        if (allocated(s%clip_low)) then
            o%clip_low = real(s%clip_low, c_double)
            o%has_clip_low = 1
        end if
        if (allocated(s%clip_high)) then
            o%clip_high = real(s%clip_high, c_double)
            o%has_clip_high = 1
        end if
        o%edge_color = rgba_to_c(s%edge_color)
        o%fill_color = rgba_to_c(s%fill_color)
        o%line_color = rgba_to_c(s%line_color)
        o%background = rgba_to_c(s%background)
        o%nan_color = rgba_to_c(s%nan_color)
        o%iso_color = rgba_to_c(s%iso_color)
        o%vector_color = rgba_to_c(s%vector_color)
        o%outline_color = rgba_to_c(s%outline_color)
        o%cutaway_tint = rgba_to_c(s%cutaway_tint)
        if (allocated(s%cutaways)) then
            if (size(s%cutaways) > 0) then
                if (mod(size(s%cutaways), 6) /= 0 .or. size(s%cutaways) > 12) then
                    o%num_cutaways = -1  ! the C library refuses it by name
                else
                    cuts = levels(size(levels) - 11:size(levels))
                    cuts(1:size(s%cutaways)) = real(s%cutaways, c_double)
                    levels(size(levels) - 11:size(levels)) = cuts
                    o%cutaways = c_loc(levels(size(levels) - 11))
                    o%num_cutaways = int(size(s%cutaways) / 6, c_int32_t)
                end if
            end if
        end if
    end subroutine

    subroutine text_opts_to_c(t, o)
        type(mio_text_settings), intent(in) :: t
        type(mio_text_opts_t), intent(out) :: o
        o%encoding = int(t%encoding, c_int32_t)
        o%color_depth = int(t%color_depth, c_int32_t)
        o%format = int(t%format, c_int32_t)
        o%cols = int(t%cols, c_int32_t)
        o%rows = int(t%rows, c_int32_t)
        o%cell_pixel_width = int(t%cell_pixel_width, c_int32_t)
        o%cell_pixel_height = int(t%cell_pixel_height, c_int32_t)
        o%tmux = merge(1_c_int32_t, 0_c_int32_t, t%tmux)
        o%notes = merge(1_c_int32_t, 0_c_int32_t, t%notes)
        o%cell_aspect = real(t%cell_aspect, c_double)
    end subroutine

    !> Render this mesh into an RGBA frame (volume cells are drawn through their
    !> boundary skin). `settings` defaults to `mio_render_settings()`. See
    !> doc/tui.md. Free the frame with `call frame%free()`.
    function mesh_render(self, settings, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        type(mio_render_settings), intent(in), optional :: settings
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_frame) :: out
        type(mio_render_settings) :: s
        type(mio_render_opts_t) :: o
        character(kind=c_char, len=STRBUF_LEN), target :: bufs(9)
        real(c_double), target :: levels(1036)

        if (present(settings)) s = settings
        if (allocated(s%iso_levels)) then
            if (size(s%iso_levels) > size(levels) - 12) then
                call handle_failure('render', 'at most 1024 iso_levels', stat, errmsg)
                return
            end if
        end if
        call render_opts_to_c(s, o, bufs, levels)
        out%handle = c_mio_render(self%handle, o)
        if (.not. c_associated(out%handle)) then
            call handle_failure('render', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Render this mesh as text sized to `text%cols` x `text%rows` cells: terminal
    !> cells with ANSI colour, plain glyphs or an HTML page, or a graphics-protocol
    !> image. Returns '' on failure.
    function mesh_render_text(self, settings, text, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        type(mio_render_settings), intent(in), optional :: settings
        type(mio_text_settings), intent(in), optional :: text
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: out
        type(mio_render_settings) :: s
        type(mio_text_settings) :: t
        type(mio_render_opts_t) :: o
        type(mio_text_opts_t) :: to
        character(kind=c_char, len=STRBUF_LEN), target :: bufs(9)
        real(c_double), target :: levels(1036)
        character(c_char), allocatable :: buf(:)
        integer(c_int64_t) :: n

        out = ''
        if (present(settings)) s = settings
        if (present(text)) t = text
        if (allocated(s%iso_levels)) then
            if (size(s%iso_levels) > size(levels) - 12) then
                call handle_failure('render_text', 'at most 1024 iso_levels', stat, errmsg)
                return
            end if
        end if
        call render_opts_to_c(s, o, bufs, levels)
        call text_opts_to_c(t, to)
        allocate (buf(1))
        n = c_mio_render_text(self%handle, o, to, buf, 0_c_int64_t)
        if (n < 0 .or. n >= int(huge(0) - 1, c_int64_t)) then
            call handle_failure('render_text', mio_error_message(), stat, errmsg)
            return
        end if
        deallocate (buf)
        allocate (buf(int(n) + 1))
        n = c_mio_render_text(self%handle, o, to, buf, n + 1)
        if (n < 0) then
            call handle_failure('render_text', mio_error_message(), stat, errmsg)
            return
        end if
        out = from_c_buf(buf, int(n))
        call clear_status(stat, errmsg)
    end function

    !> Render this mesh to a file chosen by extension: `.png`, `.txt`, `.ansi`,
    !> `.html` or `.cast`. Text forms are sized by `text%cols` x `text%rows`.
    subroutine mesh_write_snapshot(self, path, settings, text, png_compress, cast_frames, &
                                   cast_fps, cast_degrees, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: path
        type(mio_render_settings), intent(in), optional :: settings
        type(mio_text_settings), intent(in), optional :: text
        integer, intent(in), optional :: png_compress, cast_frames
        real(real64), intent(in), optional :: cast_fps, cast_degrees
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_render_settings) :: s
        type(mio_text_settings) :: t
        type(mio_render_opts_t) :: o
        type(mio_text_opts_t) :: to
        type(mio_snapshot_opts_t) :: so
        character(kind=c_char, len=STRBUF_LEN), target :: bufs(9)
        real(c_double), target :: levels(1036)
        integer(c_int) :: rc

        if (present(settings)) s = settings
        if (present(text)) t = text
        if (allocated(s%iso_levels)) then
            if (size(s%iso_levels) > size(levels) - 12) then
                call handle_failure('write_snapshot', 'at most 1024 iso_levels', stat, errmsg)
                return
            end if
        end if
        call render_opts_to_c(s, o, bufs, levels)
        call text_opts_to_c(t, to)
        if (present(png_compress)) so%png_compress = int(png_compress, c_int32_t)
        if (present(cast_frames)) so%cast_frames = int(cast_frames, c_int32_t)
        if (present(cast_fps)) so%cast_fps = real(cast_fps, c_double)
        if (present(cast_degrees)) so%cast_degrees = real(cast_degrees, c_double)
        rc = c_mio_write_snapshot(trim(path)//c_null_char, self%handle, o, to, so)
        if (rc /= 0) then
            call handle_failure('write_snapshot', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end subroutine

    subroutine frame_free(self)
        class(mio_frame), intent(inout) :: self
        if (c_associated(self%handle)) call c_mio_frame_free(self%handle)
        self%handle = c_null_ptr
    end subroutine

    function frame_is_valid(self) result(ok)
        class(mio_frame), intent(in) :: self
        logical :: ok
        ok = c_associated(self%handle)
    end function

    function frame_width(self) result(n)
        class(mio_frame), intent(in) :: self
        integer :: n
        n = int(c_mio_frame_width(self%handle))
    end function

    function frame_height(self) result(n)
        class(mio_frame), intent(in) :: self
        integer :: n
        n = int(c_mio_frame_height(self%handle))
    end function

    !> The pixels as `rgba(4, width, height)`, each channel 0..255, row 1 on top.
    subroutine frame_get_rgba(self, rgba)
        class(mio_frame), intent(in) :: self
        integer, allocatable, intent(out) :: rgba(:, :, :)
        integer(c_int8_t), pointer :: raw(:)
        integer :: w, h, i
        w = self%width()
        h = self%height()
        allocate (rgba(4, w, h))
        if (w * h == 0) return
        call c_f_pointer(c_mio_frame_rgba(self%handle), raw, [4 * w * h])
        do i = 0, w * h - 1
            rgba(1:4, mod(i, w) + 1, i / w + 1) = iand(int(raw(4 * i + 1:4 * i + 4)), 255)
        end do
    end subroutine

    !> The input cell drawn at each pixel as `ids(width, height)`, -1 for none.
    subroutine frame_get_cell_ids(self, ids)
        class(mio_frame), intent(in) :: self
        integer(int64), allocatable, intent(out) :: ids(:, :)
        integer(c_int64_t), pointer :: raw(:)
        integer :: w, h, i
        w = self%width()
        h = self%height()
        allocate (ids(w, h))
        if (w * h == 0) return
        call c_f_pointer(c_mio_frame_cell_ids(self%handle), raw, [w * h])
        do i = 0, w * h - 1
            ids(mod(i, w) + 1, i / w + 1) = int(raw(i + 1), int64)
        end do
    end subroutine

    !> `.true.` when a data field was mapped, with its range in `vmin`, `vmax`.
    function frame_range(self, vmin, vmax) result(colored)
        class(mio_frame), intent(in) :: self
        real(real64), intent(out), optional :: vmin, vmax
        logical :: colored
        real(c_double) :: lo, hi
        colored = c_mio_frame_range(self%handle, lo, hi) == 1
        if (present(vmin)) vmin = real(lo, real64)
        if (present(vmax)) vmax = real(hi, real64)
    end function

    function frame_num_notes(self) result(n)
        class(mio_frame), intent(in) :: self
        integer :: n
        n = int(max(c_mio_frame_num_notes(self%handle), 0_c_int64_t))
    end function

    !> Note `index` (1-based): the colour range, ticks and keys of the picture.
    function frame_note(self, index) result(s)
        class(mio_frame), intent(in) :: self
        integer, intent(in) :: index
        character(:), allocatable :: s
        character(c_char), allocatable :: buf(:)
        integer(c_int64_t) :: n
        s = ''
        allocate (buf(1))
        n = c_mio_frame_note(self%handle, int(index - 1, c_int64_t), buf, 0_c_int64_t)
        if (n < 0 .or. n >= int(huge(0) - 1, c_int64_t)) return
        deallocate (buf)
        allocate (buf(int(n) + 1))
        n = c_mio_frame_note(self%handle, int(index - 1, c_int64_t), buf, n + 1)
        if (n >= 0) s = from_c_buf(buf, int(n))
    end function

    !> Encode the frame as text. A cell encoding needs a frame whose size is a
    !> multiple of its cell grid; `mio_mesh%render_text` sizes it for you.
    function frame_text(self, text, stat, errmsg) result(out)
        class(mio_frame), intent(in) :: self
        type(mio_text_settings), intent(in), optional :: text
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: out
        type(mio_text_settings) :: t
        type(mio_text_opts_t) :: to
        character(c_char), allocatable :: buf(:)
        integer(c_int64_t) :: n
        out = ''
        if (present(text)) t = text
        call text_opts_to_c(t, to)
        allocate (buf(1))
        n = c_mio_frame_text(self%handle, to, buf, 0_c_int64_t)
        if (n < 0 .or. n >= int(huge(0) - 1, c_int64_t)) then
            call handle_failure('frame_text', mio_error_message(), stat, errmsg)
            return
        end if
        deallocate (buf)
        allocate (buf(int(n) + 1))
        n = c_mio_frame_text(self%handle, to, buf, n + 1)
        if (n < 0) then
            call handle_failure('frame_text', mio_error_message(), stat, errmsg)
            return
        end if
        out = from_c_buf(buf, int(n))
        call clear_status(stat, errmsg)
    end function

    !> The frame as an RGBA PNG: `compress` 0 (default) stores the pixels with no
    !> zlib, the same bytes everywhere; 1..9 compress through zlib where the build
    !> has it. The bytes are returned as 0..255 values.
    subroutine frame_png(self, bytes, compress, stat, errmsg)
        class(mio_frame), intent(in) :: self
        integer, allocatable, intent(out) :: bytes(:)
        integer, intent(in), optional :: compress
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int8_t), allocatable :: raw(:)
        integer(c_int64_t) :: n
        integer(c_int32_t) :: level
        level = 0
        if (present(compress)) level = int(compress, c_int32_t)
        allocate (bytes(0))
        allocate (raw(1))
        n = c_mio_frame_png(self%handle, level, raw, 0_c_int64_t)
        if (n < 0 .or. n >= int(huge(0) - 1, c_int64_t)) then
            call handle_failure('frame_png', mio_error_message(), stat, errmsg)
            return
        end if
        deallocate (raw)
        allocate (raw(int(n)))
        n = c_mio_frame_png(self%handle, level, raw, n)
        if (n < 0) then
            call handle_failure('frame_png', mio_error_message(), stat, errmsg)
            return
        end if
        deallocate (bytes)
        allocate (bytes(int(n)))
        bytes = iand(int(raw), 255)
        call clear_status(stat, errmsg)
    end subroutine

    !> The MIO_COLOR_* depth a terminal advertises, from the values of NO_COLOR,
    !> COLORTERM and TERM (omit one that is unset).
    function mio_detect_color_depth(no_color, color_term, term) result(depth)
        character(*), intent(in), optional :: no_color, color_term, term
        integer :: depth
        character(kind=c_char, len=STRBUF_LEN), target :: a, b, c
        type(c_ptr) :: pa, pb, pc
        pa = c_null_ptr
        pb = c_null_ptr
        pc = c_null_ptr
        if (present(no_color)) then
            a = trim(no_color)//c_null_char
            pa = c_loc(a(1:1))
        end if
        if (present(color_term)) then
            b = trim(color_term)//c_null_char
            pb = c_loc(b(1:1))
        end if
        if (present(term)) then
            c = trim(term)//c_null_char
            pc = c_loc(c(1:1))
        end if
        depth = int(c_mio_detect_color_depth(pa, pb, pc))
    end function

end module meshioplusplus
