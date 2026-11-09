/* Combining / comparing, statistics and the opaque-result operations.
 * A port of the middle part of the R binding's src/ops.c (R_mio_merge through
 * R_mio_partition_labels); see mio_oct.h for the layout rules.
 *
 * Ownership rule for the opaque C results (split, partition, reorder, refine,
 * decimate, convert_cells, ...): read the maps and counters first, then TAKE
 * the mesh(es) out of the result, then free it, and only then wrap the meshes
 * into handles. Every mesh handed to MATLAB therefore owns its handle and
 * cannot dangle when the result goes away, and no raise can leave a result
 * allocated. This mirrors the R and Julia bindings.
 */

#include "mio_oct.h"

#include <math.h>
#include <string.h>

/* --- local helpers ------------------------------------------------------- */

/* The cellstr argument `a` ({} / [] for none) as a name array. */
static const char **oc_names(const mxArray *a, int64_t *n, const char *what) {
    return (const char **)mio_o_names(a, n, what);
}

/* Read a string through the C getter protocol WITHOUT raising (so the caller
 * can free a native result first). NULL on failure. */
static char *oc_try_getcstring(mio_o_str_getter getter, void *ctx) {
    int64_t n = getter(ctx, NULL, 0);
    if (n < 0) return NULL;
    char *buf = (char *)mxCalloc((size_t)n + 1, 1);
    if (n > 0) {
        if (getter(ctx, buf, n + 1) < 0) {
            mxFree(buf);
            return NULL;
        }
        buf[n] = '\0';
    }
    return buf;
}

/* --- combining / comparing ----------------------------------------------- */

static mxArray *O_mio_merge(const mxArray **a) {
    const mxArray *meshes = a[0];
    if (!mxIsCell(meshes) || mxGetNumberOfElements(meshes) < 1) {
        mio_o_error("merge", "`meshes` must be a non-empty cell array of mio.MeshHandle objects");
    }
    int64_t n = (int64_t)mxGetNumberOfElements(meshes);
    const mio_mesh **handles = (const mio_mesh **)mxMalloc((size_t)n * sizeof(mio_mesh *));
    for (int64_t i = 0; i < n; ++i) handles[i] = mio_o_mesh(mxGetCell(meshes, (mwIndex)i));
    mio_mesh *out = mio_merge(handles, n, mio_o_bool(a[1], "weld"), mio_o_double(a[2], "atol"),
                              mio_o_bool(a[3], "source_tag"), mio_o_int(a[4], "data_policy"),
                              mio_o_bool(a[5], "drop_duplicate_cells"));
    mxFree(handles);
    if (out == NULL) mio_o_fail("merge");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_interpolate(const mxArray **a) {
    int64_t count = 0;
    const char *const *names = mio_o_names(a[3], &count, "arrays");
    mio_mesh *out = mio_interpolate(mio_o_mesh(a[0]), mio_o_mesh(a[1]), mio_o_opt_string(a[2]),
                                    names, count, mio_o_bool(a[4], "extrapolate"),
                                    mio_o_double(a[5], "default_value"),
                                    mio_o_opt_string(a[6]));
    if (out == NULL) mio_o_fail("interpolate");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_conservative_interpolate(const mxArray **a) {
    int64_t count = 0;
    const char *const *names = mio_o_names(a[2], &count, "arrays");
    mio_mesh *out = mio_conservative_interpolate(mio_o_mesh(a[0]), mio_o_mesh(a[1]), names, count,
                                                 mio_o_double(a[3], "default_value"),
                                                 mio_o_opt_string(a[4]));
    if (out == NULL) mio_o_fail("conservative_interpolate");
    return mio_o_wrap_mesh(out);
}

/* error_size (roadmap 6.10.2): a per-point size field from a per-cell error
 * indicator. Returns struct(mesh, num_cells, num_skipped, num_clamped,
 * min_size, max_size). */
static mxArray *O_mio_error_size(const mxArray **a) {
    mio_error_size_opts opts;
    mio_error_size_opts_init(&opts);
    opts.target_error = mio_o_double(a[1], "target_error");
    opts.error_array = mio_o_string(a[2], "error_array");
    opts.order = mio_o_double(a[3], "order");
    opts.min_size = mio_o_double(a[4], "min_size");
    opts.max_size = mio_o_double(a[5], "max_size");
    opts.output = mio_o_string(a[6], "output");
    opts.overwrite = mio_o_bool(a[7], "overwrite");
    mio_error_size_report report;
    mio_mesh *out = mio_error_size(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("error_size");
    const char *names[] = {"mesh",        "num_cells", "num_skipped",
                           "num_clamped", "min_size",  "max_size"};
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    mio_o_scalar((double)report.num_cells),
                    mio_o_scalar((double)report.num_skipped),
                    mio_o_scalar((double)report.num_clamped),
                    mio_o_scalar(report.min_size),
                    mio_o_scalar(report.max_size)};
    return mio_o_struct(6, names, v);
}

/* check_conformity (roadmap 3.2.1): hanging nodes, T-junctions, partly
 * overlapping and unwelded facets. Returns struct(mesh, dimension, ...,
 * interface_mode, warnings); the .m wrapper splits the warnings. */
static mxArray *O_mio_check_conformity(const mxArray **a) {
    mio_conformity_opts opts;
    mio_conformity_opts_init(&opts);
    opts.tolerance = mio_o_double(a[1], "tolerance");
    opts.attach_flags = mio_o_bool(a[2], "attach_flags");
    mio_conformity_report report;
    mio_mesh *out = mio_check_conformity(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("check_conformity");
    static const char *modes[] = {"conforming", "weld", "proximity"};
    const char *mode = (report.interface_mode >= 0 && report.interface_mode <= 2)
                           ? modes[report.interface_mode]
                           : "unknown";
    /* Read the warnings before wrapping: nothing native is live while it raises. */
    mxArray *w = mio_o_getbuffer(mio_check_conformity_warnings, "check_conformity");
    const char *names[] = {"mesh",
                           "dimension",
                           "num_hanging_nodes",
                           "num_t_junctions",
                           "num_partial_overlaps",
                           "num_unwelded_facets",
                           "num_boundary_facets",
                           "num_non_manifold_facets",
                           "conforming",
                           "interface_mode",
                           "warnings"};
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    mio_o_scalar((double)report.dimension),
                    mio_o_scalar((double)report.num_hanging_nodes),
                    mio_o_scalar((double)report.num_t_junctions),
                    mio_o_scalar((double)report.num_partial_overlaps),
                    mio_o_scalar((double)report.num_unwelded_facets),
                    mio_o_scalar((double)report.num_boundary_facets),
                    mio_o_scalar((double)report.num_non_manifold_facets),
                    mio_o_logical(report.conforming != 0),
                    mio_o_str(mode),
                    w};
    return mio_o_struct(11, names, v);
}

/* check_tetrahedra: inverted, degenerate and overlapping tetrahedra. The mesh is
 * not changed. Returns struct(num_tetrahedra, ..., valid, warnings); `zero_volume_factor`
 * [] keeps the C default and the .m wrapper splits the warnings. */
static mxArray *O_mio_check_tetrahedra(const mxArray **a) {
    mio_check_tetrahedra_opts opts;
    mio_check_tetrahedra_opts_init(&opts);
    if (!mio_o_is_null(a[1])) opts.zero_volume_factor = mio_o_double(a[1], "zero_volume_factor");
    mio_check_tetrahedra_report report;
    memset(&report, 0, sizeof(report));
    mio_o_check(mio_check_tetrahedra(mio_o_mesh(a[0]), &opts, &report), "check_tetrahedra");
    mxArray *w = mio_o_getbuffer(mio_check_tetrahedra_warnings, "check_tetrahedra");
    const char *names[] = {"num_tetrahedra", "num_inverted", "num_degenerate", "num_overlaps",
                           "min_volume",     "max_volume",   "valid",          "warnings"};
    mxArray *v[] = {mio_o_scalar((double)report.num_tetrahedra),
                    mio_o_scalar((double)report.num_inverted),
                    mio_o_scalar((double)report.num_degenerate),
                    mio_o_scalar((double)report.num_overlaps),
                    mio_o_scalar(report.min_volume),
                    mio_o_scalar(report.max_volume),
                    mio_o_logical(report.valid != 0),
                    w};
    return mio_o_struct(8, names, v);
}

/* boundary_layer (roadmap 6.9): prisms grown from a surface. Returns
 * struct(mesh, num_cells, ..., num_front_crossings). */
static mxArray *O_mio_boundary_layer(const mxArray **a) {
    mio_boundary_layer_opts opts;
    mio_boundary_layer_opts_init(&opts);
    opts.first_height = mio_o_double(a[1], "first_height");
    opts.num_layers = (int32_t)mio_o_int(a[2], "num_layers");
    opts.growth = mio_o_double(a[3], "growth");
    opts.reverse = mio_o_bool(a[4], "reverse");
    opts.collision_fraction = mio_o_double(a[5], "collision_fraction");
    opts.smooth_iterations = (int32_t)mio_o_int(a[6], "smooth_iterations");
    opts.max_shrink = (int32_t)mio_o_int(a[7], "max_shrink");
    {
        const char *c = mio_o_string(a[8], "collapse");
        if (strcmp(c, "points") == 0)
            opts.collapse = MIO_BOUNDARY_LAYER_COLLAPSE_POINTS;
        else if (strcmp(c, "layers") == 0)
            opts.collapse = MIO_BOUNDARY_LAYER_COLLAPSE_LAYERS;
        else
            mio_o_error("boundary_layer",
                        "meshio++: boundary_layer: collapse must be 'points' or 'layers'");
    }
    opts.normal_smooth_iterations = (int32_t)mio_o_int(a[9], "normal_smooth_iterations");
    opts.front_check = mio_o_bool(a[10], "front_check");
    mio_boundary_layer_report report;
    mio_mesh *out = mio_boundary_layer(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("boundary_layer");
    const char *names[] = {"mesh",          "num_cells",       "num_wedges",
                           "num_hexahedra", "num_pyramids",    "num_tetrahedra",
                           "num_shrunk",    "num_collapsed",   "num_not_visible",
                           "min_scaled_jacobian", "min_thickness", "max_thickness",
                           "cap_closed",    "num_fewer_layers", "num_front_crossings"};
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    mio_o_scalar((double)report.num_cells),
                    mio_o_scalar((double)report.num_wedges),
                    mio_o_scalar((double)report.num_hexahedra),
                    mio_o_scalar((double)report.num_pyramids),
                    mio_o_scalar((double)report.num_tetrahedra),
                    mio_o_scalar((double)report.num_shrunk),
                    mio_o_scalar((double)report.num_collapsed),
                    mio_o_scalar((double)report.num_not_visible),
                    mio_o_scalar(report.min_scaled_jacobian),
                    mio_o_scalar(report.min_thickness),
                    mio_o_scalar(report.max_thickness),
                    mio_o_logical(report.cap_closed != 0),
                    mio_o_scalar((double)report.num_fewer_layers),
                    mio_o_scalar((double)report.num_front_crossings)};
    return mio_o_struct(15, names, v);
}

/* mio_o_getstring adapters: the warnings of this thread's last extrude,
 * revolve and mirror calls. */
/* The struct an extrude or revolve call returns. */
static mxArray *oc_sweep_struct(mio_mesh *out, const mio_sweep_report *report, mxArray *warnings) {
    const char *names[] = {"mesh",          "num_cells",      "num_hexahedra", "num_wedges",
                           "num_pyramids",  "num_tetrahedra", "num_quads",     "num_triangles",
                           "num_flipped",   "num_on_axis",    "min_scaled_jacobian",
                           "full",          "warnings"};
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    mio_o_scalar((double)report->num_cells),
                    mio_o_scalar((double)report->num_hexahedra),
                    mio_o_scalar((double)report->num_wedges),
                    mio_o_scalar((double)report->num_pyramids),
                    mio_o_scalar((double)report->num_tetrahedra),
                    mio_o_scalar((double)report->num_quads),
                    mio_o_scalar((double)report->num_triangles),
                    mio_o_scalar((double)report->num_flipped),
                    mio_o_scalar((double)report->num_on_axis),
                    mio_o_scalar(report->min_scaled_jacobian),
                    mio_o_logical(report->full != 0),
                    warnings};
    return mio_o_struct(13, names, v);
}

/* extrude (roadmap 3.1.2): a surface swept into a volume, or a curve into a
 * surface, along a vector. Returns struct(mesh, num_cells, ..., warnings). */
static mxArray *O_mio_extrude(const mxArray **a) {
    mio_extrude_opts opts;
    mio_extrude_opts_init(&opts);
    mio_o_vec3(a[1], opts.direction, "direction");
    opts.num_layers = (int32_t)mio_o_int(a[2], "num_layers");
    opts.ratio = mio_o_double(a[3], "ratio");
    size_t noffs = 0;
    double *offs = mio_o_doubles(a[4], &noffs, "offsets");
    opts.num_offsets = (int64_t)noffs;
    opts.offsets = noffs > 0 ? offs : NULL;
    opts.record_layers = mio_o_bool(a[5], "record_layers");
    opts.max_cells = (int64_t)mio_o_double(a[6], "max_cells");
    mio_sweep_report report;
    mio_mesh *out = mio_extrude(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("extrude");
    mxArray *w = mio_o_getbuffer(mio_extrude_warnings, "extrude");
    return oc_sweep_struct(out, &report, w);
}

/* revolve (roadmap 3.1.3): the same around an axis. */
static mxArray *O_mio_revolve(const mxArray **a) {
    mio_revolve_opts opts;
    mio_revolve_opts_init(&opts);
    mio_o_vec3(a[1], opts.axis, "axis");
    mio_o_vec3(a[2], opts.axis_point, "axis_point");
    opts.angle = mio_o_double(a[3], "angle");
    opts.num_layers = (int32_t)mio_o_int(a[4], "num_layers");
    opts.tolerance = mio_o_double(a[5], "tolerance");
    opts.record_layers = mio_o_bool(a[6], "record_layers");
    opts.max_cells = (int64_t)mio_o_double(a[7], "max_cells");
    mio_sweep_report report;
    mio_mesh *out = mio_revolve(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("revolve");
    mxArray *w = mio_o_getbuffer(mio_revolve_warnings, "revolve");
    return oc_sweep_struct(out, &report, w);
}

/* mirror (roadmap 3.2.5): reflect across a plane, every cell right side out.
 * Returns struct(mesh, num_flipped_cells, ..., warnings). */
static mxArray *O_mio_mirror(const mxArray **a) {
    mio_mirror_opts opts;
    mio_mirror_opts_init(&opts);
    mio_o_vec3(a[1], opts.origin, "origin");
    mio_o_vec3(a[2], opts.normal, "normal");
    opts.keep_original = mio_o_bool(a[3], "keep_original");
    opts.weld = mio_o_bool(a[4], "weld");
    opts.tolerance = mio_o_double(a[5], "tolerance");
    opts.reflect_vector_data = mio_o_bool(a[6], "reflect_vector_data");
    mio_mirror_report report;
    mio_mesh *out = mio_mirror(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("mirror");
    mxArray *w = mio_o_getbuffer(mio_mirror_warnings, "mirror");
    const char *names[] = {"mesh",              "num_flipped_cells", "num_points_on_plane",
                           "num_dropped_cells", "num_dropped_sides", "warnings"};
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    mio_o_scalar((double)report.num_flipped_cells),
                    mio_o_scalar((double)report.num_points_on_plane),
                    mio_o_scalar((double)report.num_dropped_cells),
                    mio_o_scalar((double)report.num_dropped_sides),
                    w};
    return mio_o_struct(6, names, v);
}

/* find_intersections(mesh, other, region, region_other, segments, max_pairs).
 * `pairs` is (num_pairs x 2), 1-based global cell ids, one pair per row --
 * the same shape as in R. R refuses more than INT_MAX pairs (its matrix
 * limit); a MATLAB matrix has no such limit, so there is no check here. */
static mxArray *O_mio_find_intersections(const mxArray **a) {
    mio_find_intersections_opts opts;
    mio_find_intersections_opts_init(&opts);
    opts.region = mio_o_opt_string(a[2]);
    opts.region_b = mio_o_opt_string(a[3]);
    opts.segments = mio_o_bool(a[4], "segments");
    opts.max_pairs = (int64_t)mio_o_double(a[5], "max_pairs");
    const mio_mesh *b = mio_o_is_null(a[1]) ? NULL : mio_o_mesh(a[1]);
    mio_find_intersections_report report;
    mio_find_intersections_result *result =
        mio_find_intersections(mio_o_mesh(a[0]), b, &opts, &report);
    if (result == NULL) mio_o_fail("find_intersections");
    int64_t count = 0;
    const int64_t *pairs = mio_find_intersections_result_pairs(result, &count);
    if (count < 0 || (count > 0 && pairs == NULL)) count = 0;
    mxArray *pm = mio_o_rows_from_c(pairs, MIO_INT64, count, 2, 1.0);
    mio_mesh *flagged = mio_find_intersections_result_take_mesh(result);
    mio_mesh *flagged_b = mio_find_intersections_result_take_mesh_b(result);
    mio_mesh *segs = mio_find_intersections_result_take_segments(result);
    mio_find_intersections_result_free(result);
    if (flagged == NULL) {
        if (flagged_b != NULL) mio_mesh_free(flagged_b);
        if (segs != NULL) mio_mesh_free(segs);
        mio_o_fail("find_intersections take_mesh");
    }
    mxArray *w = mio_o_getbuffer(mio_find_intersections_warnings, "find_intersections");
    const char *names[] = {"mesh",          "other",       "segments",     "pairs",
                           "num_pairs",     "num_triangle_pairs",          "num_coplanar",
                           "num_segments",  "num_cells_a", "num_cells_b",  "truncated",
                           "clean",         "warnings"};
    mxArray *v[13];
    v[0] = mio_o_wrap_mesh(flagged);
    v[1] = flagged_b ? mio_o_wrap_mesh(flagged_b) : mio_o_empty();
    v[2] = segs ? mio_o_wrap_mesh(segs) : mio_o_empty();
    v[3] = pm;
    v[4] = mio_o_scalar((double)report.num_pairs);
    v[5] = mio_o_scalar((double)report.num_triangle_pairs);
    v[6] = mio_o_scalar((double)report.num_coplanar);
    v[7] = mio_o_scalar((double)report.num_segments);
    v[8] = mio_o_scalar((double)report.num_cells_a);
    v[9] = mio_o_scalar((double)report.num_cells_b);
    v[10] = mio_o_logical(report.truncated != 0);
    v[11] = mio_o_logical(report.clean != 0);
    v[12] = w;
    return mio_o_struct(13, names, v);
}

/* thickness(mesh, other, direction, max_distance, threshold, region,
 * region_other, method, cone_angle, cone_rays, two_sided,
 * max_normal_deviation). `min_point` is a 1 x 3 row (one point). */
static mxArray *O_mio_thickness(const mxArray **a) {
    mio_thickness_opts opts;
    mio_thickness_opts_init(&opts);
    const char *dir = mio_o_string(a[2], "direction");
    if (strcmp(dir, "auto") == 0)
        opts.direction = MIO_THICKNESS_AUTO;
    else if (strcmp(dir, "inward") == 0)
        opts.direction = MIO_THICKNESS_INWARD;
    else if (strcmp(dir, "outward") == 0)
        opts.direction = MIO_THICKNESS_OUTWARD;
    else
        mio_o_error("thickness", "`direction` must be 'auto', 'inward' or 'outward'");
    opts.max_distance = mio_o_double(a[3], "max_distance");
    opts.threshold = mio_o_double(a[4], "threshold");
    opts.region = mio_o_opt_string(a[5]);
    opts.region_target = mio_o_opt_string(a[6]);
    const char *meth = mio_o_string(a[7], "method");
    if (strcmp(meth, "ray") == 0)
        opts.method = MIO_THICKNESS_RAY;
    else if (strcmp(meth, "sphere") == 0)
        opts.method = MIO_THICKNESS_SPHERE;
    else if (strcmp(meth, "cone") == 0)
        opts.method = MIO_THICKNESS_CONE;
    else
        mio_o_error("thickness", "`method` must be 'ray', 'sphere' or 'cone'");
    opts.cone_angle = mio_o_double(a[8], "cone_angle");
    opts.cone_rays = (int64_t)mio_o_double(a[9], "cone_rays");
    opts.two_sided = mio_o_bool(a[10], "two_sided");
    opts.max_normal_deviation = mio_o_double(a[11], "max_normal_deviation");
    const mio_mesh *target = mio_o_is_null(a[1]) ? NULL : mio_o_mesh(a[1]);
    mio_thickness_report r;
    mio_mesh *out = mio_thickness(mio_o_mesh(a[0]), target, &opts, &r);
    if (out == NULL) mio_o_fail("thickness");
    mxArray *w = mio_o_getbuffer(mio_thickness_warnings, "thickness");
    const char *names[] = {"mesh",      "direction", "num_cells",    "num_hit",  "num_missed",
                           "num_below", "min",       "max",          "mean",     "min_cell",
                           "min_point", "max_distance", "warnings", "method", "p2", "p95"};
    static const char *const kMethods[] = {"ray", "sphere", "cone"};
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    mio_o_str(r.direction == MIO_THICKNESS_OUTWARD ? "outward" : "inward"),
                    mio_o_scalar((double)r.num_cells),
                    mio_o_scalar((double)r.num_hit),
                    mio_o_scalar((double)r.num_missed),
                    mio_o_scalar((double)r.num_below),
                    mio_o_scalar(r.min),
                    mio_o_scalar(r.max),
                    mio_o_scalar(r.mean),
                    mio_o_scalar(r.min_cell < 0 ? 0.0 : (double)(r.min_cell + 1)),
                    mio_o_row3(r.min_point),
                    mio_o_scalar(r.max_distance),
                    w,
                    mio_o_str(kMethods[r.method >= 0 && r.method <= 2 ? r.method : 0]),
                    mio_o_scalar(r.p2),
                    mio_o_scalar(r.p95)};
    return mio_o_struct(16, names, v);
}

/* deviation(scan, reference, type, params, tolerance, bins, range_min,
 * range_max, max_normal_angle, reference_region, register); reference []
 * measures against the primitive (type 1-5, ten params). */
static mxArray *O_mio_deviation(const mxArray **a) {
    mio_deviation_opts opts;
    mio_deviation_opts_init(&opts);
    opts.tolerance = mio_o_double(a[4], "tolerance");
    opts.bins = (int64_t)mio_o_double(a[5], "bins");
    opts.range_min = mio_o_double(a[6], "range_min");
    opts.range_max = mio_o_double(a[7], "range_max");
    opts.max_normal_angle = mio_o_double(a[8], "max_normal_angle");
    opts.reference_region = mio_o_opt_string(a[9]);
    opts.register_first = mio_o_bool(a[10], "register");
    mio_deviation_report r;
    mio_mesh *out;
    if (mio_o_is_null(a[1])) {
        size_t n = 0;
        double *p = mio_o_doubles(a[3], &n, "params");
        if (n != 10) mio_o_error("deviation", "`params` must be ten numbers");
        out = mio_deviation_primitive(mio_o_mesh(a[0]), (int32_t)mio_o_double(a[2], "type"), p,
                                      &opts, &r);
    } else {
        out = mio_deviation(mio_o_mesh(a[0]), mio_o_mesh(a[1]), &opts, &r);
    }
    if (out == NULL) mio_o_fail("deviation");
    mxArray *w = mio_o_getbuffer(mio_deviation_warnings, "deviation");
    const int64_t nb = r.num_bins;
    mxArray *counts = mio_o_column(nb);
    mxArray *edges = mio_o_column(nb + 1);
    int64_t *tmp = (int64_t *)mxMalloc((size_t)(nb > 0 ? nb : 1) * sizeof(int64_t));
    mio_deviation_histogram(tmp, mxGetPr(edges), nb);
    for (int64_t k = 0; k < nb; ++k)
        mxGetPr(counts)[k] = (double)tmp[k];
    mxFree(tmp);
    mxArray *pct = mio_o_matrix(1, 5);
    for (int k = 0; k < 5; ++k)
        mxGetPr(pct)[k] = r.percentiles[k];
    const char *names[] = {"mesh",          "num_points",        "num_measured",
                           "min",           "max",               "mean",
                           "rms",           "sigma",             "percentiles",
                           "histogram_edges", "histogram_counts", "num_outside_range",
                           "share_within",  "suggested_range",   "num_without_normal",
                           "warnings",      "registered",        "matrix",
                           "registration_rms"};
    mxArray *mat = mio_o_matrix(4, 4);
    double dm[16];
    mio_deviation_matrix(dm);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) mxGetPr(mat)[j * 4 + i] = dm[4 * i + j];
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    mio_o_scalar((double)r.num_points),
                    mio_o_scalar((double)r.num_measured),
                    mio_o_scalar(r.min),
                    mio_o_scalar(r.max),
                    mio_o_scalar(r.mean),
                    mio_o_scalar(r.rms),
                    mio_o_scalar(r.sigma),
                    pct,
                    edges,
                    counts,
                    mio_o_scalar((double)r.num_outside_range),
                    mio_o_scalar(r.share_within),
                    mio_o_scalar(r.suggested_range),
                    mio_o_scalar((double)r.num_without_normal),
                    w,
                    mio_o_logical(r.registered != 0),
                    mat,
                    mio_o_scalar(r.registration_rms)};
    return mio_o_struct(19, names, v);
}

/* deviation_regions(mesh, tolerance, field, prefix, min_area); `regions` is an
 * n x 1 struct array. */
static mxArray *O_mio_deviation_regions(const mxArray **a) {
    mio_deviation_regions_opts opts;
    mio_deviation_regions_opts_init(&opts);
    opts.tolerance = mio_o_double(a[1], "tolerance");
    opts.field = mio_o_opt_string(a[2]);
    opts.prefix = mio_o_opt_string(a[3]);
    opts.min_area = mio_o_double(a[4], "min_area");
    mio_deviation_regions_report r;
    mio_mesh *out = mio_deviation_regions(mio_o_mesh(a[0]), &opts, &r);
    if (out == NULL) mio_o_fail("deviation_regions");
    const char *fields[] = {"name",   "sign",   "num_cells", "area",  "extreme",
                            "mean",   "volume", "centroid",  "bounds"};
    mxArray *regions = mio_o_struct_array(r.num_regions, 9, fields);
    for (int64_t k = 0; k < r.num_regions; ++k) {
        mio_deviation_region_info g;
        if (mio_deviation_region_get(k, &g) != MIO_OK) mio_o_fail("deviation_regions");
        char buf[512];
        mio_deviation_region_name(k, buf, (int64_t)sizeof(buf));
        mxArray *bounds = mio_o_matrix(1, 6);
        for (int c = 0; c < 6; ++c)
            mxGetPr(bounds)[c] = g.bounds[c];
        mxSetField(regions, (mwIndex)k, "name", mio_o_str(buf));
        mxSetField(regions, (mwIndex)k, "sign", mio_o_scalar((double)g.sign));
        mxSetField(regions, (mwIndex)k, "num_cells", mio_o_scalar((double)g.num_cells));
        mxSetField(regions, (mwIndex)k, "area", mio_o_scalar(g.area));
        mxSetField(regions, (mwIndex)k, "extreme", mio_o_scalar(g.extreme));
        mxSetField(regions, (mwIndex)k, "mean", mio_o_scalar(g.mean));
        mxSetField(regions, (mwIndex)k, "volume", mio_o_scalar(g.volume));
        mxSetField(regions, (mwIndex)k, "centroid", mio_o_row3(g.centroid));
        mxSetField(regions, (mwIndex)k, "bounds", bounds);
    }
    mxArray *w = mio_o_getbuffer(mio_deviation_regions_warnings, "deviation_regions");
    const char *names[] = {"mesh",        "regions",     "num_excess",
                           "num_missing", "num_dropped", "warnings"};
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    regions,
                    mio_o_scalar((double)r.num_excess),
                    mio_o_scalar((double)r.num_missing),
                    mio_o_scalar((double)r.num_dropped),
                    w};
    return mio_o_struct(6, names, v);
}

/* exaggerate(mesh, field, scale, clamp, smooth_fraction). */
static mxArray *O_mio_exaggerate(const mxArray **a) {
    mio_exaggerate_opts opts;
    mio_exaggerate_opts_init(&opts);
    opts.field = mio_o_opt_string(a[1]);
    opts.scale = mio_o_double(a[2], "scale");
    opts.clamp = mio_o_double(a[3], "clamp");
    opts.smooth_fraction = mio_o_double(a[4], "smooth_fraction");
    mio_exaggerate_report r;
    mio_mesh *out = mio_exaggerate(mio_o_mesh(a[0]), &opts, &r);
    if (out == NULL) mio_o_fail("exaggerate");
    mxArray *w = mio_o_getbuffer(mio_exaggerate_warnings, "exaggerate");
    const char *names[] = {"mesh", "clamp", "radius", "num_moved", "warnings"};
    mxArray *v[] = {mio_o_wrap_mesh(out), mio_o_scalar(r.clamp), mio_o_scalar(r.radius),
                    mio_o_scalar((double)r.num_moved), w};
    return mio_o_struct(5, names, v);
}

/* colorize(mesh, field, cmap, vmin, vmax, symmetric, steps, nan_color,
 * component); component is 0-based here (-1: the magnitude). */
static mxArray *O_mio_colorize(const mxArray **a) {
    mio_colorize_opts opts;
    mio_colorize_opts_init(&opts);
    opts.field = mio_o_string(a[1], "field");
    opts.cmap = mio_o_opt_string(a[2]);
    opts.vmin = mio_o_double(a[3], "vmin");
    opts.vmax = mio_o_double(a[4], "vmax");
    opts.symmetric = mio_o_double(a[5], "symmetric") != 0.0 ? 1 : 0;
    opts.steps = (int64_t)mio_o_double(a[6], "steps");
    opts.nan_color = mio_o_opt_string(a[7]);
    opts.component = (int64_t)mio_o_double(a[8], "component");
    mio_colorize_report r;
    mio_mesh *out = mio_colorize(mio_o_mesh(a[0]), &opts, &r);
    if (out == NULL) mio_o_fail("colorize");
    const char *names[] = {"mesh", "vmin", "vmax", "location"};
    mxArray *v[] = {mio_o_wrap_mesh(out), mio_o_scalar(r.vmin), mio_o_scalar(r.vmax),
                    mio_o_str(r.location == 1 ? "cell" : "point")};
    return mio_o_struct(4, names, v);
}

static mxArray *O_mio_has_netgen(const mxArray **a) {
    (void)a;
    return mio_o_logical(mio_has_netgen() != 0);
}

/* The backend names tetrahedralize takes besides "auto" (roadmap 7.7.2). */
static mxArray *O_mio_tetrahedralize_backends(const mxArray **a) {
    (void)a;
    const int64_t n = mio_tet_mesher_count();
    if (n < 0) mio_o_fail("tetrahedralize_backends");
    mxArray *out = mio_o_cellstr(n);
    for (int64_t i = 0; i < n; ++i) {
        const char *name = mio_tet_mesher_name(i);
        mxSetCell(out, (mwIndex)i, mio_o_str(name != NULL ? name : ""));
    }
    return out;
}

/* remesh_quads (roadmap 7.7.3): field-aligned quads through an installed
 * QuadriFlow or Instant Meshes. */
static mxArray *O_mio_remesh_quads(const mxArray **a) {
    mio_remesh_quads_opts opts;
    mio_remesh_quads_opts_init(&opts);
    opts.target_faces = (int64_t)mio_o_double(a[1], "target_faces");
    opts.backend = mio_o_string(a[2], "backend");
    opts.sharp = mio_o_bool(a[3], "sharp");
    opts.boundary = mio_o_bool(a[4], "boundary");
    opts.adaptive = mio_o_bool(a[5], "adaptive");
    opts.dominant = mio_o_bool(a[6], "dominant");
    opts.seed = (int32_t)mio_o_int(a[7], "seed");
    opts.max_cells = (int64_t)mio_o_double(a[8], "max_cells");
    mio_remesh_quads_report report;
    mio_mesh *out = mio_remesh_quads(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("remesh_quads");
    const char *names[] = {"mesh",
                           "num_quads",
                           "num_triangles",
                           "num_points",
                           "min_scaled_jacobian",
                           "mean_scaled_jacobian",
                           "hausdorff_to_surface",
                           "backend"};
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    mio_o_scalar((double)report.num_quads),
                    mio_o_scalar((double)report.num_triangles),
                    mio_o_scalar((double)report.num_points),
                    mio_o_scalar(report.min_scaled_jacobian),
                    mio_o_scalar(report.mean_scaled_jacobian),
                    mio_o_scalar(report.hausdorff_to_surface),
                    mio_o_str(report.backend_name != NULL ? report.backend_name : "")};
    return mio_o_struct(8, names, v);
}

/* The quad remeshers remesh_quads takes by name. */
static mxArray *O_mio_remesh_quads_backends(const mxArray **a) {
    (void)a;
    const int64_t n = mio_quad_remesher_count();
    if (n < 0) mio_o_fail("remesh_quads_backends");
    mxArray *out = mio_o_cellstr(n);
    for (int64_t i = 0; i < n; ++i) {
        const char *name = mio_quad_remesher_name(i);
        mxSetCell(out, (mwIndex)i, mio_o_str(name != NULL ? name : ""));
    }
    return out;
}

/* cut_level_set and cut_quadrature (doc/cut.md): the part level sets select,
 * through CutCells. `planes` is flat, six numbers per plane (the .m wrapper
 * flattens an (n x 6) matrix or a cell of 6-vectors). */
static void oc_cut_opts(mio_cut_opts *opts, const mxArray *level_sets, const mxArray *planes,
                        const mxArray *selection, const mxArray *include_uncut,
                        const mxArray *backend, const mxArray *triangulate,
                        const mxArray *record_parent_ids) {
    mio_cut_opts_init(opts);
    int64_t n = 0;
    const char *const *names = mio_o_names(level_sets, &n, "level_sets");
    opts->level_sets = (const char **)names;
    opts->num_level_sets = n;
    size_t np = 0;
    double *p = mio_o_doubles(planes, &np, "planes");
    if (np % 6 != 0) {
        mio_o_error("cut",
                    "meshio++: cut: planes must hold six numbers per plane (origin, then normal)");
    }
    opts->planes = np > 0 ? p : NULL;
    opts->num_planes = (int64_t)(np / 6);
    opts->selection = mio_o_string(selection, "selection");
    opts->backend = mio_o_string(backend, "backend");
    opts->include_uncut = mio_o_bool(include_uncut, "include_uncut");
    opts->triangulate = mio_o_bool(triangulate, "triangulate");
    opts->record_parent_ids = mio_o_bool(record_parent_ids, "record_parent_ids");
}

/* cut_level_set(mesh, level_sets, planes, selection, include_uncut, backend,
 *               degree, triangulate, record_parent_ids) */
static mxArray *O_mio_cut_level_set(const mxArray **a) {
    mio_cut_opts opts;
    oc_cut_opts(&opts, a[1], a[2], a[3], a[4], a[5], a[7], a[8]);
    opts.degree = (int32_t)mio_o_int(a[6], "degree");
    mio_cut_report report;
    mio_mesh *out = mio_cut_level_set(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("cut_level_set");
    const char *names[] = {"mesh",       "num_cells", "num_cut_cells", "num_uncut_cells",
                           "num_pieces", "dimension", "measure",       "backend"};
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    mio_o_scalar((double)report.num_cells),
                    mio_o_scalar((double)report.num_cut_cells),
                    mio_o_scalar((double)report.num_uncut_cells),
                    mio_o_scalar((double)report.num_pieces),
                    mio_o_scalar((double)report.dimension),
                    mio_o_scalar(report.measure),
                    mio_o_str(report.backend_name != NULL ? report.backend_name : "")};
    return mio_o_struct(8, names, v);
}

/* cut_quadrature(mesh, level_sets, planes, selection, order, include_uncut,
 *                backend, triangulate, record_parent_ids) */
static mxArray *O_mio_cut_quadrature(const mxArray **a) {
    mio_cut_opts opts;
    oc_cut_opts(&opts, a[1], a[2], a[3], a[5], a[6], a[7], a[8]);
    opts.order = (int32_t)mio_o_int(a[4], "order");
    mio_cut_report report;
    mio_mesh *out = mio_cut_quadrature(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("cut_quadrature");
    const char *names[] = {"mesh", "num_rules", "num_points", "sum_weights", "dimension",
                           "backend"};
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    mio_o_scalar((double)report.num_rules),
                    mio_o_scalar((double)report.num_points),
                    mio_o_scalar(report.measure),
                    mio_o_scalar((double)report.dimension),
                    mio_o_str(report.backend_name != NULL ? report.backend_name : "")};
    return mio_o_struct(6, names, v);
}

static mxArray *O_mio_has_cutcells(const mxArray **a) {
    (void)a;
    return mio_o_logical(mio_has_cutcells() != 0);
}

static mxArray *O_mio_has_mmg(const mxArray **a) {
    (void)a;
    return mio_o_logical(mio_has_mmg() != 0);
}

static mxArray *O_mio_has_mmg_elas(const mxArray **a) {
    (void)a;
    return mio_o_logical(mio_has_mmg_elas() != 0);
}

/* The report of an MMG operation, with its mesh, as a struct; `extra` more
 * fields (mesh_domain's counts) are appended. `quality` is the 3 x 1 column
 * (min, mean, max). */
static mxArray *oc_mmg_result(mio_mesh *out, const mio_mmg_report *report, int extra,
                              const char **extra_names, mxArray **extra_values) {
    const char *names[11] = {"mesh",          "status",           "library",
                             "num_points_in", "num_points_out",   "num_elements_in",
                             "num_elements_out", "quality",       "num_warnings"};
    mxArray *v[11];
    mxArray *q = mio_o_column(3);
    mxGetPr(q)[0] = report->quality_min;
    mxGetPr(q)[1] = report->quality_mean;
    mxGetPr(q)[2] = report->quality_max;
    v[0] = mio_o_wrap_mesh(out);
    v[1] = mio_o_str(report->status == 0 ? "success" : "low_failure");
    v[2] = mio_o_str(report->library);
    v[3] = mio_o_scalar((double)report->num_points_in);
    v[4] = mio_o_scalar((double)report->num_points_out);
    v[5] = mio_o_scalar((double)report->num_elements_in);
    v[6] = mio_o_scalar((double)report->num_elements_out);
    v[7] = q;
    v[8] = mio_o_scalar((double)report->num_warnings);
    for (int i = 0; i < extra && i < 2; ++i) {
        names[9 + i] = extra_names[i];
        v[9 + i] = extra_values[i];
    }
    return mio_o_struct(9 + (extra < 2 ? extra : 2), names, v);
}

/* The shared MMG parameters: sizes, the ridge angle, 1-based frozen points and
   the required region names. */
static void oc_mmg_params(mio_mmg_params *p, const mxArray *hmin, const mxArray *hmax,
                          const mxArray *hausd, const mxArray *hgrad,
                          const mxArray *angle_detection, const mxArray *frozen,
                          const mxArray *required) {
    mio_mmg_params_init(p);
    p->hmin = mio_o_double(hmin, "hmin");
    p->hmax = mio_o_double(hmax, "hmax");
    p->hausd = mio_o_double(hausd, "hausd");
    p->hgrad = mio_o_double(hgrad, "hgrad");
    p->angle_detection = mio_o_double(angle_detection, "angle_detection");
    int64_t n_frozen = 0;
    p->frozen = mio_o_ids_0based(frozen, &n_frozen, "frozen");
    p->num_frozen = n_frozen;
    int64_t n_required = 0;
    const char **names = oc_names(required, &n_required, "required");
    if (n_required > 0) {
        p->required = names;
        p->num_required = n_required;
    }
}

/* adapt(mesh, hsiz, size, metric, hmin, hmax, hausd, hgrad, angle_detection,
 *       frozen, required) */
static mxArray *O_mio_adapt(const mxArray **a) {
    mio_adapt_opts opts;
    mio_adapt_opts_init(&opts);
    oc_mmg_params(&opts.params, a[4], a[5], a[6], a[7], a[8], a[9], a[10]);
    opts.hsiz = mio_o_double(a[1], "hsiz");
    const char *s = mio_o_string(a[2], "size");
    const char *m = mio_o_string(a[3], "metric");
    opts.size = s[0] ? s : NULL;
    opts.metric = m[0] ? m : NULL;
    mio_mmg_report report;
    mio_mesh *out = mio_adapt(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("adapt");
    return oc_mmg_result(out, &report, 0, NULL, NULL);
}

/* discretize_levelset (roadmap 3.2.2). `materials` is [] or a numeric k x 4
   matrix of (ref, split, ref_inside, ref_outside) rows -- the same column-major
   k x 4 matrix as in R, so it is read the same way. */
static mxArray *O_mio_discretize_levelset(const mxArray **a) {
    mio_discretize_levelset_opts opts;
    mio_discretize_levelset_opts_init(&opts);
    oc_mmg_params(&opts.params, a[5], a[6], a[7], a[8], a[9], a[10], a[11]);
    opts.level_set = mio_o_string(a[1], "level_set");
    opts.isovalue = mio_o_double(a[2], "isovalue");
    opts.remove_small_components = mio_o_double(a[3], "remove_small_components");
    const mxArray *materials = a[4];
    if (!mio_o_is_null(materials)) {
        if (!mxIsDouble(materials) || mxIsComplex(materials) ||
            mxGetNumberOfElements(materials) % 4 != 0)
            mio_o_error("discretize_levelset", "materials must be a numeric matrix with 4 columns");
        size_t k = mxGetNumberOfElements(materials) / 4;
        const double *m = mxGetPr(materials);
        mio_mmg_material *mats =
            (mio_mmg_material *)mxCalloc(k > 0 ? k : 1, sizeof(mio_mmg_material));
        for (size_t i = 0; i < k; ++i) {
            mats[i].ref = (int64_t)m[i];
            mats[i].split = m[i + k] != 0.0;
            mats[i].ref_inside = (int64_t)m[i + 2 * k];
            mats[i].ref_outside = (int64_t)m[i + 3 * k];
        }
        opts.materials = mats;
        opts.num_materials = (int64_t)k;
    }
    mio_mmg_report report;
    mio_mesh *out = mio_discretize_levelset(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("discretize_levelset");
    return oc_mmg_result(out, &report, 0, NULL, NULL);
}

/* mesh_domain (roadmap 3.2.3). `holes` is [] or a numeric k x 2 matrix of
   points, one inside each hole. Returns the MMG report plus num_subdomains and
   num_holes. */
static mxArray *O_mio_mesh_domain(const mxArray **a) {
    mio_mesh_domain_opts opts;
    mio_mesh_domain_opts_init(&opts);
    oc_mmg_params(&opts.params, a[4], a[5], a[6], a[7], a[8], a[9], a[10]);
    opts.hsiz = mio_o_double(a[1], "hsiz");
    const char *s = mio_o_string(a[2], "size");
    opts.size = s[0] ? s : NULL;
    const mxArray *holes = a[3];
    if (!mio_o_is_null(holes)) {
        if (!mxIsDouble(holes) || mxIsComplex(holes) || mxGetNumberOfElements(holes) % 2 != 0)
            mio_o_error("mesh_domain", "holes must be a numeric matrix with 2 columns");
        size_t k = mxGetNumberOfElements(holes) / 2;
        const double *h = mxGetPr(holes);
        double *xy = (double *)mxMalloc((k > 0 ? 2 * k : 1) * sizeof(double));
        for (size_t i = 0; i < k; ++i) {
            xy[2 * i] = h[i];
            xy[2 * i + 1] = h[i + k];
        }
        opts.holes = xy;
        opts.num_holes = (int64_t)k;
    }
    mio_mmg_report report;
    int64_t num_subdomains = 0, num_holes = 0;
    mio_mesh *out =
        mio_mesh_domain(mio_o_mesh(a[0]), &opts, &report, &num_subdomains, &num_holes);
    if (out == NULL) mio_o_fail("mesh_domain");
    const char *extra_names[] = {"num_subdomains", "num_holes"};
    mxArray *extra[] = {mio_o_scalar((double)num_subdomains), mio_o_scalar((double)num_holes)};
    return oc_mmg_result(out, &report, 2, extra_names, extra);
}

/* move (roadmap 3.2.4): MMG's Lagrangian mode, with ELAS. */
static mxArray *O_mio_move(const mxArray **a) {
    mio_move_opts opts;
    mio_move_opts_init(&opts);
    oc_mmg_params(&opts.params, a[4], a[5], a[6], a[7], a[8], a[9], a[10]);
    opts.displacement = mio_o_string(a[1], "displacement");
    int64_t n_moving = 0;
    const char **names = oc_names(a[2], &n_moving, "moving");
    if (n_moving > 0) {
        opts.moving = names;
        opts.num_moving = n_moving;
    }
    opts.lagrangian = (int32_t)mio_o_double(a[3], "lagrangian");
    mio_mmg_report report;
    mio_mesh *out = mio_move(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("move");
    return oc_mmg_result(out, &report, 0, NULL, NULL);
}

static mxArray *O_mio_metric_from_hessian(const mxArray **a) {
    mio_mesh *out = mio_metric_from_hessian(
        mio_o_mesh(a[0]), mio_o_string(a[1], "hessian"), mio_o_double(a[2], "error"),
        mio_o_double(a[3], "hmin"), mio_o_double(a[4], "hmax"), mio_o_string(a[5], "output"));
    if (out == NULL) mio_o_fail("metric_from_hessian");
    return mio_o_wrap_mesh(out);
}

/* tetrahedralize (roadmap 6.7): tetrahedra filling a closed triangle surface.
 * Returns struct(mesh, num_tetrahedra, ..., backend). */
static mxArray *O_mio_tetrahedralize(const mxArray **a) {
    mio_tetrahedralize_opts opts;
    mio_tetrahedralize_opts_init(&opts);
    {
        const char *b = mio_o_string(a[1], "backend");
        if (strcmp(b, "auto") == 0)
            opts.backend = MIO_TETRAHEDRALIZE_BACKEND_AUTO;
        else if (strcmp(b, "netgen") == 0)
            opts.backend = MIO_TETRAHEDRALIZE_BACKEND_NETGEN;
        else if (strcmp(b, "octree") == 0)
            opts.backend = MIO_TETRAHEDRALIZE_BACKEND_OCTREE;
        else if (b[0] != '\0') {
            /* A named mesher (the installed tetgen or gmsh); an unknown name
             * fails in the library, listing the available ones. */
            opts.backend = MIO_TETRAHEDRALIZE_BACKEND_NAMED;
            opts.provider = b;
        } else
            mio_o_error("tetrahedralize",
                        "meshio++: tetrahedralize: backend must be 'auto', 'netgen', 'octree' "
                        "or a mesher's name");
    }
    opts.max_size = mio_o_double(a[2], "max_size");
    opts.max_volume = mio_o_double(a[3], "max_volume");
    opts.grading = mio_o_double(a[4], "grading");
    opts.optimize_steps = (int32_t)mio_o_int(a[5], "optimize_steps");
    opts.max_cells = (int64_t)mio_o_double(a[6], "max_cells");
    opts.face_array = mio_o_string(a[7], "face_array");
    mio_tetrahedralize_report report;
    mio_mesh *out = mio_tetrahedralize(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("tetrahedralize");
    const char *names[] = {"mesh",
                           "num_tetrahedra",
                           "num_points",
                           "num_new_points",
                           "num_regions",
                           "num_boundary_triangles",
                           "boundary_preserved",
                           "hausdorff_to_surface",
                           "min_dihedral_deg",
                           "max_dihedral_deg",
                           "volume",
                           "backend"};
    mxArray *v[] = {mio_o_wrap_mesh(out),
                    mio_o_scalar((double)report.num_tetrahedra),
                    mio_o_scalar((double)report.num_points),
                    mio_o_scalar((double)report.num_new_points),
                    mio_o_scalar((double)report.num_regions),
                    mio_o_scalar((double)report.num_boundary_triangles),
                    mio_o_logical(report.boundary_preserved != 0),
                    mio_o_scalar(report.hausdorff_to_surface),
                    mio_o_scalar(report.min_dihedral_deg),
                    mio_o_scalar(report.max_dihedral_deg),
                    mio_o_scalar(report.volume),
                    mio_o_str(report.backend_name != NULL ? report.backend_name
                              : report.backend == MIO_TETRAHEDRALIZE_BACKEND_NETGEN
                                  ? "netgen"
                                  : "octree")};
    return mio_o_struct(12, names, v);
}

static mxArray *O_mio_undo_green(const mxArray **a) {
    int64_t ngroups = 0, nremoved = 0;
    mio_mesh *out = mio_undo_green(mio_o_mesh(a[0]), mio_o_mesh(a[1]), &ngroups, &nremoved);
    if (out == NULL) mio_o_fail("undo_green");
    const char *names[] = {"mesh", "num_groups_undone", "num_cells_removed"};
    mxArray *v[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)ngroups),
                    mio_o_scalar((double)nremoved)};
    return mio_o_struct(3, names, v);
}

static mxArray *O_mio_meshes_equal(const mxArray **a) {
    int equal = 0;
    mio_o_check(mio_meshes_equal(mio_o_mesh(a[0]), mio_o_mesh(a[1]), mio_o_double(a[2], "atol"),
                                 mio_o_double(a[3], "rtol"), mio_o_bool(a[4], "unordered"),
                                 &equal),
                "meshes_equal");
    return mio_o_logical(equal != 0);
}

/* diff(a, b, atol, rtol, unordered) -> struct(verdict, points, blocks);
 * `blocks` holds one column entry per compared block. */
static mxArray *O_mio_diff(const mxArray **a) {
    mio_diff_result *result = NULL;
    mio_o_check(mio_diff(mio_o_mesh(a[0]), mio_o_mesh(a[1]), mio_o_double(a[2], "atol"),
                         mio_o_double(a[3], "rtol"), mio_o_bool(a[4], "unordered"), &result),
                "diff");

    static const char *const kVerdicts[] = {"identical", "equal_within_tolerance",
                                            "different"};
    int v = (int)mio_diff_result_verdict(result);
    if (v < 0 || v > 2) v = 2;

    double mabs = 0.0, mrel = 0.0;
    int64_t worst = -1, nex = 0;
    if (mio_diff_result_point_summary(result, &mabs, &mrel, &worst, &nex) != MIO_OK) {
        mio_diff_result_free(result);
        mio_o_fail("diff point summary");
    }

    int64_t nb = mio_diff_result_num_block_diffs(result);
    if (nb < 0) nb = 0;
    mxArray *tm = mio_o_logical_column(nb);
    mxArray *cm = mio_o_logical_column(nb);
    mxArray *cc = mio_o_column(nb);
    for (int64_t i = 0; i < nb; ++i) {
        int t = 0, c = 0;
        int64_t k = 0;
        if (mio_diff_result_block(result, i, &t, &c, &k) != MIO_OK) {
            mio_diff_result_free(result);
            mio_o_fail("diff block");
        }
        mxGetLogicals(tm)[i] = t != 0;
        mxGetLogicals(cm)[i] = c != 0;
        mxGetPr(cc)[i] = (double)k;
    }
    mio_diff_result_free(result);

    /* worst_index is a flat index into the coordinate buffer: 1-based here,
     * with the C API's -1 "no offender" becoming 0. */
    const char *pnames[] = {"max_abs_error", "max_rel_error", "worst_index", "num_exceeding"};
    mxArray *pvalues[] = {mio_o_scalar(mabs), mio_o_scalar(mrel),
                          mio_o_scalar(worst < 0 ? 0.0 : (double)(worst + 1)),
                          mio_o_scalar((double)nex)};
    mxArray *pts = mio_o_struct(4, pnames, pvalues);

    const char *bnames[] = {"type_mismatch", "count_mismatch", "conn_mismatch_count"};
    mxArray *bvalues[] = {tm, cm, cc};
    mxArray *blocks = mio_o_struct(3, bnames, bvalues);

    const char *names[] = {"verdict", "points", "blocks"};
    mxArray *values[] = {mio_o_str(kVerdicts[v]), pts, blocks};
    return mio_o_struct(3, names, values);
}

/* --- statistics ---------------------------------------------------------- */

/* stats(mesh). bbox_min, bbox_max, extent and centroid are 1 x 3 rows (each a
 * point, or a displacement, in space). */
static mxArray *O_mio_stats(const mxArray **a) {
    mio_stats_report r;
    memset(&r, 0, sizeof(r));
    mio_o_check(mio_stats(mio_o_mesh(a[0]), &r), "stats");
    const char *names[] = {"num_points", "num_cells",       "bbox_min",        "bbox_max",
                           "extent",     "centroid",        "total_area",      "signed_volume",
                           "unsigned_volume", "num_inverted"};
    mxArray *values[] = {mio_o_scalar((double)r.num_points),
                         mio_o_scalar((double)r.num_cells),
                         mio_o_row3(r.bbox_min),
                         mio_o_row3(r.bbox_max),
                         mio_o_row3(r.extent),
                         mio_o_row3(r.centroid),
                         mio_o_scalar(r.total_area),
                         mio_o_scalar(r.signed_volume),
                         mio_o_scalar(r.unsigned_volume),
                         mio_o_scalar((double)r.num_inverted)};
    return mio_o_struct(10, names, values);
}

static mxArray *O_mio_compute_bandwidth(const mxArray **a) {
    int64_t bw = mio_compute_bandwidth(mio_o_mesh(a[0]));
    if (bw < 0) mio_o_fail("compute_bandwidth");
    return mio_o_scalar((double)bw);
}

/* --- opaque-result operations -------------------------------------------- */

static mxArray *O_mio_reorder(const mxArray **a) {
    mio_reorder_result *r = mio_reorder(mio_o_mesh(a[0]), mio_o_string(a[1], "method"));
    if (r == NULL) mio_o_fail("reorder");

    const void *data = NULL;
    mio_dtype dt;
    int64_t n = 0;
    if (mio_reorder_result_node_perm(r, &data, &dt, &n) != MIO_OK) {
        mio_reorder_result_free(r);
        mio_o_fail("reorder node_perm");
    }
    mxArray *node_perm = mio_o_shift_map((const int64_t *)data, n);

    int64_t nperms = mio_reorder_result_num_cell_perms(r);
    if (nperms < 0) nperms = 0;
    mxArray *cell_perms = mio_o_cell(nperms);
    for (int64_t b = 0; b < nperms; ++b) {
        const void *d = NULL;
        int64_t len = 0;
        if (mio_reorder_result_cell_perm(r, b, &d, &dt, &len) != MIO_OK) {
            mio_reorder_result_free(r);
            mio_o_fail("reorder cell_perm");
        }
        mxSetCell(cell_perms, (mwIndex)b, mio_o_shift_map((const int64_t *)d, len));
    }

    mio_mesh *out = mio_reorder_result_take_mesh(r);
    mio_reorder_result_free(r);
    if (out == NULL) mio_o_fail("reorder take_mesh");
    const char *names[] = {"mesh", "node_perm", "cell_perms"};
    mxArray *values[] = {mio_o_wrap_mesh(out), node_perm, cell_perms};
    return mio_o_struct(3, names, values);
}

typedef struct {
    const mio_split_result *result;
    int64_t index;
} oc_split_ctx;

static int64_t oc_split_key_getter(void *c, char *buf, int64_t buflen) {
    oc_split_ctx *x = (oc_split_ctx *)c;
    return mio_split_result_key(x->result, x->index, buf, buflen);
}

/* split(mesh, by, tag_name) -> struct(meshes, keys): an n x 1 cell of mesh
 * handles and the n x 1 cellstr of their keys -- R's named list of pieces.
 * The .m wrapper returns them as two outputs. */
static mxArray *O_mio_split(const mxArray **a) {
    mio_split_result *r =
        mio_split(mio_o_mesh(a[0]), mio_o_string(a[1], "by"), mio_o_opt_string(a[2]));
    if (r == NULL) mio_o_fail("split");
    int64_t n = mio_split_result_count(r);
    if (n < 0) {
        mio_split_result_free(r);
        mio_o_fail("split count");
    }
    mxArray *keys = mio_o_cellstr(n);
    mio_mesh **taken = (mio_mesh **)mxCalloc((size_t)(n > 0 ? n : 1), sizeof(mio_mesh *));
    for (int64_t i = 0; i < n; ++i) {
        oc_split_ctx ctx = {r, i};
        char *k = oc_try_getcstring(oc_split_key_getter, &ctx);
        /* TAKE, never borrow: the piece must outlive the result handle. */
        mio_mesh *piece = k != NULL ? mio_split_result_take_mesh(r, i) : NULL;
        if (piece == NULL) {
            for (int64_t j = 0; j < i; ++j) mio_mesh_free(taken[j]);
            mio_split_result_free(r);
            mio_o_fail(k == NULL ? "split key" : "split take_mesh");
        }
        mxSetCell(keys, (mwIndex)i, mxCreateString(k));
        mxFree(k);
        taken[i] = piece;
    }
    mio_split_result_free(r);
    mxArray *pieces = mio_o_cell(n);
    for (int64_t i = 0; i < n; ++i) mxSetCell(pieces, (mwIndex)i, mio_o_wrap_mesh(taken[i]));
    mxFree(taken);
    const char *names[] = {"meshes", "keys"};
    mxArray *values[] = {pieces, keys};
    return mio_o_struct(2, names, values);
}

/* convert_cells / refine / decimate / decimate_volume share a result shape:
 * one mesh, a point map, and one cell map per input block. `kind` selects
 * which family of C getters to use (0 convert_cells, 1 refine, 2 decimate,
 * 3 decimate_volume). Returns 0 on failure (nothing is raised, so the caller can
 * free the result first); on success *point_map and *cell_maps are set. */
static int oc_result_maps(void *r, int kind, mxArray **point_map, mxArray **cell_maps) {
    const void *data = NULL;
    mio_dtype dt;
    int64_t n = 0;
    mio_status st;
    if (kind == 0) {
        st = mio_convert_cells_result_point_map((mio_convert_cells_result *)r, &data, &dt, &n);
    } else if (kind == 1) {
        st = mio_refine_result_point_map((mio_refine_result *)r, &data, &dt, &n);
    } else if (kind == 2) {
        st = mio_decimate_result_point_map((mio_decimate_result *)r, &data, &dt, &n);
    } else {
        st = mio_decimate_volume_result_point_map((mio_decimate_volume_result *)r, &data, &dt,
                                                  &n);
    }
    if (st != MIO_OK) return 0;
    mxArray *pm = mio_o_shift_map((const int64_t *)data, n);

    int64_t nmaps;
    if (kind == 0) {
        nmaps = mio_convert_cells_result_num_cell_maps((mio_convert_cells_result *)r);
    } else if (kind == 1) {
        nmaps = mio_refine_result_num_cell_maps((mio_refine_result *)r);
    } else if (kind == 2) {
        nmaps = mio_decimate_result_num_cell_maps((mio_decimate_result *)r);
    } else {
        nmaps = mio_decimate_volume_result_num_cell_maps((mio_decimate_volume_result *)r);
    }
    if (nmaps < 0) nmaps = 0;
    mxArray *cm = mio_o_cell(nmaps);
    for (int64_t b = 0; b < nmaps; ++b) {
        const void *d = NULL;
        int64_t len = 0;
        if (kind == 0) {
            st = mio_convert_cells_result_cell_map((mio_convert_cells_result *)r, b, &d, &dt,
                                                   &len);
        } else if (kind == 1) {
            st = mio_refine_result_cell_map((mio_refine_result *)r, b, &d, &dt, &len);
        } else if (kind == 2) {
            st = mio_decimate_result_cell_map((mio_decimate_result *)r, b, &d, &dt, &len);
        } else {
            st = mio_decimate_volume_result_cell_map((mio_decimate_volume_result *)r, b, &d,
                                                     &dt, &len);
        }
        if (st != MIO_OK) {
            mxDestroyArray(pm);
            mxDestroyArray(cm);
            return 0;
        }
        mxSetCell(cm, (mwIndex)b, mio_o_shift_map((const int64_t *)d, len));
    }
    *point_map = pm;
    *cell_maps = cm;
    return 1;
}

/* subdivide's cell maps only -- unlike convert_cells/refine/decimate, there
 * is no point map to read first (subdivide never prunes or renumbers a
 * point), so this does not share oc_result_maps() above. NULL on failure. */
static mxArray *oc_subdivide_cell_maps(mio_subdivide_result *r) {
    int64_t nmaps = mio_subdivide_result_num_cell_maps(r);
    if (nmaps < 0) nmaps = 0;
    mxArray *cell_maps = mio_o_cell(nmaps);
    for (int64_t b = 0; b < nmaps; ++b) {
        const void *d = NULL;
        mio_dtype dt;
        int64_t len = 0;
        if (mio_subdivide_result_cell_map(r, b, &d, &dt, &len) != MIO_OK) {
            mxDestroyArray(cell_maps);
            return NULL;
        }
        mxSetCell(cell_maps, (mwIndex)b, mio_o_shift_map((const int64_t *)d, len));
    }
    return cell_maps;
}

static mxArray *O_mio_subdivide(const mxArray **a) {
    mio_subdivide_result *r =
        mio_subdivide(mio_o_mesh(a[0]), mio_o_bool(a[1], "record_parent_ids"));
    if (r == NULL) mio_o_fail("subdivide");
    mxArray *cm = oc_subdivide_cell_maps(r);
    if (cm == NULL) {
        mio_subdivide_result_free(r);
        mio_o_fail("subdivide maps");
    }
    mio_mesh *out = mio_subdivide_result_take_mesh(r);
    mio_subdivide_result_free(r);
    if (out == NULL) mio_o_fail("subdivide take_mesh");
    const char *names[] = {"mesh", "cell_maps"};
    mxArray *values[] = {mio_o_wrap_mesh(out), cm};
    return mio_o_struct(2, names, values);
}

/* agglomerate's cell map is a single FLAT array (unlike subdivide's
 * per-block cell_maps), so it needs no per-block loop at all -- just one
 * mio_o_shift_map() call. */
static mxArray *O_mio_agglomerate(const mxArray **a) {
    mio_agglomerate_opts opts;
    mio_agglomerate_opts_init(&opts);
    opts.target_group_size = mio_o_int64(a[1], "target_group_size");
    opts.merge_coplanar_faces = mio_o_bool(a[2], "merge_coplanar_faces");
    opts.coplanar_angle = mio_o_double(a[3], "coplanar_angle");
    opts.min_sphericity = mio_o_double(a[4], "min_sphericity");
    int64_t merged = 0, rejected = 0;
    mio_agglomerate_result *r = mio_agglomerate_ex(mio_o_mesh(a[0]), &opts, &merged, &rejected);
    if (r == NULL) mio_o_fail("agglomerate");
    const void *d = NULL;
    mio_dtype dt;
    int64_t len = 0;
    if (mio_agglomerate_result_cell_map(r, &d, &dt, &len) != MIO_OK) {
        mio_agglomerate_result_free(r);
        mio_o_fail("agglomerate cell_map");
    }
    mxArray *cm = mio_o_shift_map((const int64_t *)d, len);
    mio_mesh *out = mio_agglomerate_result_take_mesh(r);
    mio_agglomerate_result_free(r);
    if (out == NULL) mio_o_fail("agglomerate take_mesh");
    const char *names[] = {"mesh", "cell_map", "num_faces_merged", "num_rejected"};
    mxArray *values[] = {mio_o_wrap_mesh(out), cm, mio_o_scalar((double)merged),
                         mio_o_scalar((double)rejected)};
    return mio_o_struct(4, names, values);
}

static mxArray *O_mio_blend_steps(const mxArray **a) {
    mio_mesh *out = mio_blend_steps(mio_o_mesh(a[0]), mio_o_mesh(a[1]), mio_o_double(a[2], "w"),
                                    mio_o_bool(a[3], "blend_points"));
    if (out == NULL) mio_o_fail("blend_steps");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_convert_cells(const mxArray **a) {
    mio_convert_cells_result *r = mio_convert_cells(
        mio_o_mesh(a[0]), mio_o_string(a[1], "mode"), mio_o_bool(a[2], "record_parent_ids"));
    if (r == NULL) mio_o_fail("convert_cells");
    mxArray *pm = NULL, *cm = NULL;
    if (!oc_result_maps(r, 0, &pm, &cm)) {
        mio_convert_cells_result_free(r);
        mio_o_fail("convert_cells maps");
    }
    mio_mesh *out = mio_convert_cells_result_take_mesh(r);
    mio_convert_cells_result_free(r);
    if (out == NULL) mio_o_fail("convert_cells take_mesh");
    const char *names[] = {"mesh", "point_map", "cell_maps"};
    mxArray *values[] = {mio_o_wrap_mesh(out), pm, cm};
    return mio_o_struct(3, names, values);
}

/* tet_to_poly's vertex map is a single FLAT array, like agglomerate's cell
 * map, returned 1-based. */
static mxArray *O_mio_tet_to_poly(const mxArray **a) {
    mio_tet_to_poly_opts opts;
    mio_tet_to_poly_opts_init(&opts);
    opts.feature_angle = mio_o_double(a[1], "feature_angle");
    opts.point_data_to_cell_data = mio_o_bool(a[2], "point_data_to_cell_data");
    opts.record_vertex_ids = mio_o_bool(a[3], "record_vertex_ids");
    mio_tet_to_poly_result *r = mio_tet_to_poly(mio_o_mesh(a[0]), &opts);
    if (r == NULL) mio_o_fail("tet_to_poly");
    const void *d = NULL;
    mio_dtype dt;
    int64_t len = 0, nnew = 0, nbnd = 0, nsplit = 0;
    if (mio_tet_to_poly_result_vertex_map(r, &d, &dt, &len) != MIO_OK ||
        mio_tet_to_poly_result_counts(r, &nnew, &nbnd, &nsplit) != MIO_OK) {
        mio_tet_to_poly_result_free(r);
        mio_o_fail("tet_to_poly vertex_map");
    }
    mxArray *vm = mio_o_shift_map((const int64_t *)d, len);
    mio_mesh *out = mio_tet_to_poly_result_take_mesh(r);
    mio_tet_to_poly_result_free(r);
    if (out == NULL) mio_o_fail("tet_to_poly take_mesh");
    const char *names[] = {"mesh", "vertex_map", "num_new_points", "num_boundary_faces",
                           "num_split_vertices"};
    mxArray *values[] = {mio_o_wrap_mesh(out), vm, mio_o_scalar((double)nnew),
                         mio_o_scalar((double)nbnd), mio_o_scalar((double)nsplit)};
    return mio_o_struct(5, names, values);
}

/* poly_to_tet returns a convert_cells result, so its maps read the same way. */
static mxArray *O_mio_poly_to_tet(const mxArray **a) {
    int64_t nexact = 0, nfanned = 0;
    mio_convert_cells_result *r = mio_poly_to_tet(
        mio_o_mesh(a[0]), mio_o_bool(a[1], "record_parent_ids"), &nexact, &nfanned);
    if (r == NULL) mio_o_fail("poly_to_tet");
    mxArray *pm = NULL, *cm = NULL;
    if (!oc_result_maps(r, 0, &pm, &cm)) {
        mio_convert_cells_result_free(r);
        mio_o_fail("poly_to_tet maps");
    }
    mio_mesh *out = mio_convert_cells_result_take_mesh(r);
    mio_convert_cells_result_free(r);
    if (out == NULL) mio_o_fail("poly_to_tet take_mesh");
    const char *names[] = {"mesh", "point_map", "cell_maps", "num_exact_tets",
                           "num_fanned_cells"};
    mxArray *values[] = {mio_o_wrap_mesh(out), pm, cm, mio_o_scalar((double)nexact),
                         mio_o_scalar((double)nfanned)};
    return mio_o_struct(5, names, values);
}

/* The mio_refine_closure / mio_refine_compare code for a name, or -1. */
static int oc_refine_closure_code(const char *name) {
    if (name == NULL || name[0] == '\0') return MIO_REFINE_CLOSURE_REDGREEN;
    if (strcmp(name, "redgreen") == 0 || strcmp(name, "red-green") == 0 ||
        strcmp(name, "green") == 0)
        return MIO_REFINE_CLOSURE_REDGREEN;
    if (strcmp(name, "propagate") == 0 || strcmp(name, "red") == 0)
        return MIO_REFINE_CLOSURE_PROPAGATE;
    if (strcmp(name, "balanced") == 0 || strcmp(name, "2:1") == 0)
        return MIO_REFINE_CLOSURE_BALANCED;
    return -1;
}

static int oc_refine_compare_code(const char *name) {
    if (name == NULL || name[0] == '\0') return MIO_REFINE_LT;
    if (strcmp(name, "<") == 0 || strcmp(name, "lt") == 0) return MIO_REFINE_LT;
    if (strcmp(name, "<=") == 0 || strcmp(name, "le") == 0) return MIO_REFINE_LE;
    if (strcmp(name, ">") == 0 || strcmp(name, "gt") == 0) return MIO_REFINE_GT;
    if (strcmp(name, ">=") == 0 || strcmp(name, "ge") == 0) return MIO_REFINE_GE;
    if (strcmp(name, "==") == 0 || strcmp(name, "=") == 0 || strcmp(name, "eq") == 0)
        return MIO_REFINE_EQ;
    if (strcmp(name, "!=") == 0 || strcmp(name, "ne") == 0) return MIO_REFINE_NE;
    return -1;
}

/* refine(mesh, levels, record_parent_ids, cells, region, where_array,
 *        where_op, where_value, closure, record_levels, record_hierarchy) */
static mxArray *O_mio_refine(const mxArray **a) {
    mio_refine_opts opts;
    mio_refine_opts_init(&opts);
    opts.levels = (int32_t)mio_o_int(a[1], "levels");
    opts.record_parent_ids = mio_o_bool(a[2], "record_parent_ids");
    opts.record_levels = mio_o_bool(a[9], "record_levels");
    opts.record_hierarchy = mio_o_bool(a[10], "record_hierarchy");
    opts.region = mio_o_opt_string(a[4]);
    opts.predicate_array = mio_o_opt_string(a[5]);
    if (!mio_o_is_null(a[7]) && mxGetNumberOfElements(a[7]) > 0)
        opts.predicate_value = mio_o_double(a[7], "where_value");
    opts.closure = oc_refine_closure_code(mio_o_opt_string(a[8]));
    if (opts.closure < 0) mio_o_error("refine", "meshio++: refine: unknown closure");
    opts.predicate_op = oc_refine_compare_code(mio_o_opt_string(a[6]));
    if (opts.predicate_op < 0) mio_o_error("refine", "meshio++: refine: unknown comparison");

    /* Cell indices arrive as double (exact to 2^53) or any integer class;
     * either way they are 1-based on this side. */
    int64_t n_ids = 0;
    int64_t *ids = mio_o_ids_0based(a[3], &n_ids, "cells");
    if (n_ids > 0) {
        opts.cells = ids;
        opts.num_cells = n_ids;
    }

    mio_refine_result *r = mio_refine_ex(mio_o_mesh(a[0]), &opts);
    if (r == NULL) mio_o_fail("refine");
    mxArray *pm = NULL, *cm = NULL;
    if (!oc_result_maps(r, 1, &pm, &cm)) {
        mio_refine_result_free(r);
        mio_o_fail("refine maps");
    }
    mio_mesh *out = mio_refine_result_take_mesh(r);
    mio_refine_result_free(r);
    if (out == NULL) mio_o_fail("refine take_mesh");
    const char *names[] = {"mesh", "point_map", "cell_maps"};
    mxArray *values[] = {mio_o_wrap_mesh(out), pm, cm};
    return mio_o_struct(3, names, values);
}

static mxArray *O_mio_tessellate(const mxArray **a) {
    mio_tessellate_opts opts;
    mio_tessellate_report report;
    mio_tessellate_opts_init(&opts);
    opts.levels = mio_o_int64(a[1], "levels");
    opts.curved = mio_o_bool(a[2], "curved");
    opts.fields = mio_o_bool(a[3], "fields");
    opts.record_stencil = mio_o_bool(a[4], "record_stencil");

    mio_mesh *out = mio_tessellate(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("tessellate");
    const char *names[] = {"mesh",
                           "num_source_points",
                           "num_source_cells",
                           "num_curved_source_cells",
                           "num_pass_through_source_cells",
                           "num_points",
                           "num_cells",
                           "num_boundary_facets",
                           "num_interior_facets",
                           "num_facets_with_bad_count"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)report.num_source_points),
                         mio_o_scalar((double)report.num_source_cells),
                         mio_o_scalar((double)report.num_curved_source_cells),
                         mio_o_scalar((double)report.num_pass_through_source_cells),
                         mio_o_scalar((double)report.num_points),
                         mio_o_scalar((double)report.num_cells),
                         mio_o_scalar((double)report.num_boundary_facets),
                         mio_o_scalar((double)report.num_interior_facets),
                         mio_o_scalar((double)report.num_facets_with_bad_count)};
    return mio_o_struct(10, names, values);
}

/* decimate(mesh, target_ratio, target_faces, max_error, placement,
 *          preserve_boundary, preserve_features, feature_angle, frozen) */
static mxArray *O_mio_decimate(const mxArray **a) {
    mio_decimate_opts opts;
    mio_decimate_opts_init(&opts);
    opts.target_ratio = mio_o_double(a[1], "target_ratio");
    opts.target_faces = mio_o_int64(a[2], "target_faces");
    opts.max_error = mio_o_double(a[3], "max_error");
    opts.placement = mio_o_opt_string(a[4]);
    opts.preserve_boundary = mio_o_bool(a[5], "preserve_boundary");
    opts.preserve_features = mio_o_bool(a[6], "preserve_features");
    opts.feature_angle = mio_o_double(a[7], "feature_angle");
    int64_t n_frozen = 0;
    opts.frozen = mio_o_ids_0based(a[8], &n_frozen, "frozen");
    opts.num_frozen = n_frozen;

    mio_decimate_result *r = mio_decimate_ex(mio_o_mesh(a[0]), &opts);
    if (r == NULL) mio_o_fail("decimate");
    mxArray *pm = NULL, *cm = NULL;
    if (!oc_result_maps(r, 2, &pm, &cm)) {
        mio_decimate_result_free(r);
        mio_o_fail("decimate maps");
    }
    mxArray *fr = mio_o_scalar((double)mio_decimate_result_faces_removed(r));
    mxArray *pr = mio_o_scalar((double)mio_decimate_result_points_removed(r));
    mxArray *cr = mio_o_scalar((double)mio_decimate_result_collapses_rejected(r));
    mxArray *me = mio_o_scalar(mio_decimate_result_max_error_applied(r));
    mio_mesh *out = mio_decimate_result_take_mesh(r);
    mio_decimate_result_free(r);
    if (out == NULL) mio_o_fail("decimate take_mesh");
    const char *names[] = {"mesh",           "point_map",          "cell_maps",
                           "faces_removed",  "points_removed",     "collapses_rejected",
                           "max_error_applied"};
    mxArray *values[] = {mio_o_wrap_mesh(out), pm, cm, fr, pr, cr, me};
    return mio_o_struct(7, names, values);
}

/* decimate_volume(mesh, ratio, target_cells, max_error, placement,
 *                 preserve_boundary, preserve_features, feature_angle, frozen) */
static mxArray *O_mio_decimate_volume(const mxArray **a) {
    mio_decimate_volume_opts opts;
    mio_decimate_volume_opts_init(&opts);
    opts.target_ratio = mio_o_double(a[1], "target_ratio");
    opts.target_cells = mio_o_int64(a[2], "target_cells");
    opts.max_error = mio_o_double(a[3], "max_error");
    opts.placement = mio_o_opt_string(a[4]);
    opts.preserve_boundary = mio_o_bool(a[5], "preserve_boundary");
    opts.preserve_features = mio_o_bool(a[6], "preserve_features");
    opts.feature_angle = mio_o_double(a[7], "feature_angle");
    int64_t n_frozen = 0;
    opts.frozen = mio_o_ids_0based(a[8], &n_frozen, "frozen");
    opts.num_frozen = n_frozen;

    mio_decimate_volume_result *r = mio_decimate_volume_ex(mio_o_mesh(a[0]), &opts);
    if (r == NULL) mio_o_fail("decimate_volume");
    mxArray *pm = NULL, *cm = NULL;
    if (!oc_result_maps(r, 3, &pm, &cm)) {
        mio_decimate_volume_result_free(r);
        mio_o_fail("decimate_volume maps");
    }
    mxArray *tr = mio_o_scalar((double)mio_decimate_volume_result_tets_removed(r));
    mxArray *pr = mio_o_scalar((double)mio_decimate_volume_result_points_removed(r));
    mxArray *cr = mio_o_scalar((double)mio_decimate_volume_result_collapses_rejected(r));
    mxArray *me = mio_o_scalar(mio_decimate_volume_result_max_error_applied(r));
    mio_mesh *out = mio_decimate_volume_result_take_mesh(r);
    mio_decimate_volume_result_free(r);
    if (out == NULL) mio_o_fail("decimate_volume take_mesh");
    const char *names[] = {"mesh",          "point_map",          "cell_maps",
                           "tets_removed",  "points_removed",     "collapses_rejected",
                           "max_error_applied"};
    mxArray *values[] = {mio_o_wrap_mesh(out), pm, cm, tr, pr, cr, me};
    return mio_o_struct(7, names, values);
}

static const char *kPartitionFields[] = {"part_id", "mesh", "point_map", "cell_maps"};

/* partition(mesh, nparts, method, imbalance, mode, seed, record_ids,
 *           ghost_layers, weights_key) -> an n x 1 struct array with fields
 * part_id (0-based, as in R), mesh, point_map and cell_maps. Every map and
 * mesh is read out of the result first; the meshes are wrapped after it is
 * freed. */
static mxArray *O_mio_partition(const mxArray **a) {
    mio_partition_result *r = mio_partition(
        mio_o_mesh(a[0]), mio_o_int(a[1], "nparts"), mio_o_opt_string(a[2]),
        mio_o_double(a[3], "imbalance"), mio_o_opt_string(a[4]), mio_o_int(a[5], "seed"),
        mio_o_bool(a[6], "record_ids"), mio_o_int(a[7], "ghost_layers"), mio_o_opt_string(a[8]));
    if (r == NULL) mio_o_fail("partition");

    int64_t n = mio_partition_result_num_pieces(r);
    int64_t nmaps = mio_partition_result_num_cell_maps(r);
    if (n < 0) n = 0;
    if (nmaps < 0) nmaps = 0;

    mxArray *out = mio_o_struct_array(n, 4, kPartitionFields);
    mio_mesh **taken = (mio_mesh **)mxCalloc((size_t)(n > 0 ? n : 1), sizeof(mio_mesh *));
    for (int64_t i = 0; i < n; ++i) {
        const void *d = NULL;
        mio_dtype dt;
        int64_t len = 0;
        const char *failed = NULL;
        if (mio_partition_result_point_map(r, i, &d, &dt, &len) != MIO_OK) {
            failed = "partition point_map";
        } else {
            mxSetField(out, (mwIndex)i, "point_map", mio_o_shift_map((const int64_t *)d, len));
            mxArray *cms = mio_o_cell(nmaps);
            mxSetField(out, (mwIndex)i, "cell_maps", cms);
            for (int64_t b = 0; b < nmaps && failed == NULL; ++b) {
                const void *d2 = NULL;
                int64_t l2 = 0;
                if (mio_partition_result_cell_map(r, i, b, &d2, &dt, &l2) != MIO_OK) {
                    failed = "partition cell_map";
                } else {
                    mxSetCell(cms, (mwIndex)b, mio_o_shift_map((const int64_t *)d2, l2));
                }
            }
        }
        if (failed == NULL) {
            mxSetField(out, (mwIndex)i, "part_id",
                       mio_o_scalar((double)mio_partition_result_part_id(r, i)));
            taken[i] = mio_partition_result_take_mesh(r, i);
            if (taken[i] == NULL) failed = "partition take_mesh";
        }
        if (failed != NULL) {
            for (int64_t j = 0; j < i; ++j) mio_mesh_free(taken[j]);
            mio_partition_result_free(r);
            mxDestroyArray(out);
            mio_o_fail(failed);
        }
    }
    mio_partition_result_free(r);
    for (int64_t i = 0; i < n; ++i) mxSetField(out, (mwIndex)i, "mesh", mio_o_wrap_mesh(taken[i]));
    mxFree(taken);
    return out;
}

/* partition_labels(mesh, nparts, method, imbalance, mode, seed, weights_key,
 *                  num_cells) -> a num_cells x 1 column of part ids. */
static mxArray *O_mio_partition_labels(const mxArray **a) {
    int64_t n = mio_o_int64(a[7], "num_cells");
    int64_t *labels = (int64_t *)mxMalloc((size_t)(n > 0 ? n : 1) * sizeof(int64_t));
    mio_status st = mio_partition_labels(
        mio_o_mesh(a[0]), mio_o_int(a[1], "nparts"), mio_o_opt_string(a[2]),
        mio_o_double(a[3], "imbalance"), mio_o_opt_string(a[4]), mio_o_int(a[5], "seed"),
        mio_o_opt_string(a[6]), labels, n);
    if (st != MIO_OK) mio_o_fail("partition_labels");
    /* Part IDs, not indices: deliberately NOT shifted to 1-based. */
    mxArray *out = mio_o_int64_column(labels, n);
    mxFree(labels);
    return out;
}

/* --- table --------------------------------------------------------------- */

static const mio_o_entry kEntries[] = {
    MIO_O_ENTRY(mio_merge, 6),
    MIO_O_ENTRY(mio_interpolate, 7),
    MIO_O_ENTRY(mio_conservative_interpolate, 5),
    MIO_O_ENTRY(mio_error_size, 8),
    MIO_O_ENTRY(mio_check_conformity, 3),
    MIO_O_ENTRY(mio_check_tetrahedra, 2),
    MIO_O_ENTRY(mio_boundary_layer, 11),
    MIO_O_ENTRY(mio_extrude, 7),
    MIO_O_ENTRY(mio_revolve, 8),
    MIO_O_ENTRY(mio_mirror, 7),
    MIO_O_ENTRY(mio_find_intersections, 6),
    MIO_O_ENTRY(mio_thickness, 12),
    MIO_O_ENTRY(mio_deviation, 11),
    MIO_O_ENTRY(mio_deviation_regions, 5),
    MIO_O_ENTRY(mio_exaggerate, 5),
    MIO_O_ENTRY(mio_colorize, 9),
    MIO_O_ENTRY(mio_has_netgen, 0),
    MIO_O_ENTRY(mio_tetrahedralize_backends, 0),
    MIO_O_ENTRY(mio_remesh_quads, 9),
    MIO_O_ENTRY(mio_remesh_quads_backends, 0),
    MIO_O_ENTRY(mio_cut_level_set, 9),
    MIO_O_ENTRY(mio_cut_quadrature, 9),
    MIO_O_ENTRY(mio_has_cutcells, 0),
    MIO_O_ENTRY(mio_has_mmg, 0),
    MIO_O_ENTRY(mio_has_mmg_elas, 0),
    MIO_O_ENTRY(mio_adapt, 11),
    MIO_O_ENTRY(mio_discretize_levelset, 12),
    MIO_O_ENTRY(mio_mesh_domain, 11),
    MIO_O_ENTRY(mio_move, 11),
    MIO_O_ENTRY(mio_metric_from_hessian, 6),
    MIO_O_ENTRY(mio_tetrahedralize, 8),
    MIO_O_ENTRY(mio_undo_green, 2),
    MIO_O_ENTRY(mio_meshes_equal, 5),
    MIO_O_ENTRY(mio_diff, 5),
    MIO_O_ENTRY(mio_stats, 1),
    MIO_O_ENTRY(mio_compute_bandwidth, 1),
    MIO_O_ENTRY(mio_reorder, 2),
    MIO_O_ENTRY(mio_split, 3),
    MIO_O_ENTRY(mio_subdivide, 2),
    MIO_O_ENTRY(mio_agglomerate, 5),
    MIO_O_ENTRY(mio_blend_steps, 4),
    MIO_O_ENTRY(mio_convert_cells, 3),
    MIO_O_ENTRY(mio_tet_to_poly, 4),
    MIO_O_ENTRY(mio_poly_to_tet, 2),
    MIO_O_ENTRY(mio_refine, 11),
    MIO_O_ENTRY(mio_tessellate, 5),
    MIO_O_ENTRY(mio_decimate, 9),
    MIO_O_ENTRY(mio_decimate_volume, 9),
    MIO_O_ENTRY(mio_partition, 9),
    MIO_O_ENTRY(mio_partition_labels, 8),
};

const mio_o_entry *mio_o_table_ops_combine(size_t *n) {
    *n = sizeof(kEntries) / sizeof(kEntries[0]);
    return kEntries;
}
