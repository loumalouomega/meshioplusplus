/* The operations layer: computations *on* a mesh, not file formats.
 *
 * Ownership rule for the opaque C results (split, partition, reorder, refine,
 * decimate, convert_cells): read the maps and counters first, then TAKE the
 * mesh(es) out of the result, then free it. Every mesh handed to R therefore
 * owns its handle and cannot dangle when the result goes away. This mirrors
 * the Julia binding exactly.
 */

#include "mio_r.h"

static SEXP quality_list(const mio_surface_quality *q); /* defined below; used earlier */

#include <limits.h>
#include <math.h>
#include <string.h>

/* --- surface / skin / quality ------------------------------------------- */

SEXP R_mio_extract_surface(SEXP mesh, SEXP record_parent_ids) {
    mio_mesh *out = mio_extract_surface(mio_r_mesh(mesh),
                                        mio_r_bool(record_parent_ids, "record_parent_ids"));
    if (out == NULL) mio_r_fail("extract_surface");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_extract_skin(SEXP mesh, SEXP linearize) {
    mio_mesh *out = mio_extract_skin(mio_r_mesh(mesh), mio_r_bool(linearize, "linearize"));
    if (out == NULL) mio_r_fail("extract_skin");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_attach_quality(SEXP mesh) {
    mio_mesh *out = mio_attach_quality(mio_r_mesh(mesh));
    if (out == NULL) mio_r_fail("attach_quality");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_quality_counts(SEXP mesh) {
    int64_t nc = 0, ni = 0, nd = 0;
    mio_r_check(mio_quality_counts(mio_r_mesh(mesh), &nc, &ni, &nd), "quality_counts");
    const char *names[] = {"num_cells", "num_inverted", "num_degenerate"};
    const double numbers[] = {(double)nc, (double)ni, (double)nd};
    return mio_r_report_list(3, names, NULL, numbers);
}

/* --- geometry-preserving transforms ------------------------------------- */

SEXP R_mio_transform(SEXP mesh, SEXP matrix, SEXP rotate_vector_data) {
    SEXP m = PROTECT(Rf_coerceVector(matrix, REALSXP));
    if (Rf_xlength(m) != 16) {
        UNPROTECT(1);
        Rf_error("`matrix` must be a 4x4 affine matrix (16 values)");
    }
    /* The C API takes the matrix ROW-major, while an R matrix is column-major,
     * so this one is transposed on the way out. That is correct here and only
     * here: a 4x4 affine matrix is a mathematical object, not a mesh array. */
    double rowmajor[16];
    const double *src = REAL(m);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            rowmajor[i * 4 + j] = src[j * 4 + i];
        }
    }
    mio_mesh *out = mio_transform(mio_r_mesh(mesh), rowmajor,
                                  mio_r_bool(rotate_vector_data, "rotate_vector_data"));
    UNPROTECT(1);
    if (out == NULL) mio_r_fail("transform");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_clean(SEXP mesh, SEXP weld, SEXP atol, SEXP remove_orphans, SEXP drop_degenerate,
                 SEXP drop_duplicate_cells) {
    int64_t welded = 0, orphan = 0, degen = 0, dup = 0;
    mio_mesh *out = mio_clean(mio_r_mesh(mesh), mio_r_bool(weld, "weld"),
                              mio_r_double(atol, "atol"),
                              mio_r_bool(remove_orphans, "remove_orphans"),
                              mio_r_bool(drop_degenerate, "drop_degenerate"),
                              mio_r_bool(drop_duplicate_cells, "drop_duplicate_cells"),
                              &welded, &orphan, &degen, &dup);
    if (out == NULL) mio_r_fail("clean");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh", "points_welded", "points_removed_orphan",
                           "cells_dropped_degenerate", "cells_dropped_duplicate"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)welded, (double)orphan, (double)degen, (double)dup};
    SEXP res = mio_r_report_list(5, names, values, numbers);
    UNPROTECT(1);
    return res;
}

/* An R id vector (integer or double, 1-based, NULL for none) -> a freshly
 * R_alloc'd 0-based int64_t array; *n_ids receives its length (0 if `ids`
 * is R_NilValue or empty). R_alloc's lifetime is tied to this .Call, so no
 * explicit free is needed -- the mio_refine_opts precedent this mirrors. */
static int64_t *r_ids_0based(SEXP ids, R_xlen_t *n_ids) {
    *n_ids = 0;
    if (ids == R_NilValue || Rf_length(ids) == 0) return NULL;
    R_xlen_t n = Rf_length(ids);
    int64_t *out = (int64_t *)R_alloc((size_t)n, sizeof(int64_t));
    for (R_xlen_t i = 0; i < n; ++i) {
        double v = (TYPEOF(ids) == REALSXP) ? REAL(ids)[i] : (double)INTEGER(ids)[i];
        out[i] = (int64_t)v - 1;
    }
    *n_ids = n;
    return out;
}

SEXP R_mio_smooth(SEXP mesh, SEXP method, SEXP iterations, SEXP lambda, SEXP mu,
                  SEXP fix_boundary, SEXP preserve_features, SEXP feature_angle,
                  SEXP guard_inversion, SEXP frozen) {
    mio_smooth_opts opts;
    int64_t moved = 0, skipped = 0;
    double disp = 0.0;
    R_xlen_t n_frozen = 0;

    mio_smooth_opts_init(&opts);
    opts.method = mio_r_opt_string(method);
    opts.iterations = (int32_t)mio_r_int(iterations, "iterations");
    opts.lambda = mio_r_double(lambda, "lambda");
    opts.mu = mio_r_double(mu, "mu");
    opts.fix_boundary = mio_r_bool(fix_boundary, "fix_boundary");
    opts.preserve_features = mio_r_bool(preserve_features, "preserve_features");
    opts.feature_angle = mio_r_double(feature_angle, "feature_angle");
    opts.guard_inversion = mio_r_bool(guard_inversion, "guard_inversion");
    opts.frozen = r_ids_0based(frozen, &n_frozen);
    opts.num_frozen = (int64_t)n_frozen;

    mio_mesh *out = mio_smooth_ex(mio_r_mesh(mesh), &opts, &moved, &disp, &skipped);
    if (out == NULL) mio_r_fail("smooth");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh", "num_nodes_moved", "max_displacement",
                           "num_skipped_inversion"};
    const SEXP values[] = {mo, NULL, NULL, NULL};
    const double numbers[] = {(double)moved, disp, (double)skipped};
    SEXP res = mio_r_report_list(4, names, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_optimize_volume(SEXP mesh, SEXP max_iterations, SEXP relocate, SEXP flip,
                           SEXP preserve_boundary, SEXP min_improvement, SEXP objective,
                           SEXP max_ring_size) {
    mio_optimize_volume_opts opts;
    mio_optimize_volume_report rep;
    const char *obj;
    mio_optimize_volume_opts_init(&opts);
    opts.max_iterations = mio_r_int(max_iterations, "max_iterations");
    opts.relocate = mio_r_bool(relocate, "relocate");
    opts.flip = mio_r_bool(flip, "flip");
    opts.preserve_boundary = mio_r_bool(preserve_boundary, "preserve_boundary");
    opts.min_improvement = mio_r_double(min_improvement, "min_improvement");
    obj = mio_r_string(objective, "objective");
    if (strcmp(obj, "worst") == 0)
        opts.objective = MIO_OPTIMIZE_VOLUME_WORST;
    else if (strcmp(obj, "energy") == 0)
        opts.objective = MIO_OPTIMIZE_VOLUME_ENERGY;
    else
        Rf_error("meshio++: optimize_volume: objective must be \"worst\" or \"energy\"");
    opts.max_ring_size = mio_r_int(max_ring_size, "max_ring_size");
    mio_mesh *out = mio_optimize_volume_ex(mio_r_mesh(mesh), &opts, &rep);
    if (out == NULL) mio_r_fail("optimize_volume");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh",
                           "num_flips",
                           "num_23_flips",
                           "num_32_flips",
                           "num_vertices_moved",
                           "num_tets",
                           "min_quality_before",
                           "min_quality_after",
                           "num_edge_removals",
                           "num_multi_face_removals",
                           "num_worst_smoothed",
                           "min_dihedral_before",
                           "min_dihedral_after"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                           NULL};
    const double numbers[] = {(double)rep.num_flips, (double)rep.num_23_flips,
                              (double)rep.num_32_flips, (double)rep.num_vertices_moved,
                              (double)rep.num_tets, rep.min_quality_before, rep.min_quality_after,
                              (double)rep.num_edge_removals, (double)rep.num_multi_face_removals,
                              (double)rep.num_worst_smoothed, rep.min_dihedral_before,
                              rep.min_dihedral_after};
    SEXP res = mio_r_report_list(13, names, values, numbers);
    UNPROTECT(1);
    return res;
}

/* --- subsetting --------------------------------------------------------- */

SEXP R_mio_crop_bbox(SEXP mesh, SEXP lo, SEXP hi, SEXP mode, SEXP record_ids) {
    double l[3], h[3];
    mio_r_vec3(lo, l, "lo");
    mio_r_vec3(hi, h, "hi");
    mio_mesh *out = mio_crop_bbox(mio_r_mesh(mesh), l, h, mio_r_int(mode, "mode"),
                                  mio_r_bool(record_ids, "record_ids"));
    if (out == NULL) mio_r_fail("crop_bbox");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_crop_predicate(SEXP mesh, SEXP array, SEXP compare, SEXP value,
                          SEXP record_ids) {
    mio_mesh *out = mio_crop_predicate(mio_r_mesh(mesh), mio_r_string(array, "array"),
                                       mio_r_int(compare, "compare"),
                                       mio_r_double(value, "value"),
                                       mio_r_bool(record_ids, "record_ids"));
    if (out == NULL) mio_r_fail("crop_predicate");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_crop_plane(SEXP mesh, SEXP point, SEXP normal, SEXP mode, SEXP record_ids) {
    double p[3], n[3];
    mio_r_vec3(point, p, "point");
    mio_r_vec3(normal, n, "normal");
    mio_mesh *out = mio_crop_plane(mio_r_mesh(mesh), p, n, mio_r_int(mode, "mode"),
                                   mio_r_bool(record_ids, "record_ids"));
    if (out == NULL) mio_r_fail("crop_plane");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_slice(SEXP mesh, SEXP origin, SEXP normal, SEXP record_parent_ids) {
    double o[3], n[3];
    mio_r_vec3(origin, o, "origin");
    mio_r_vec3(normal, n, "normal");
    mio_mesh *out = mio_slice(mio_r_mesh(mesh), o, n,
                              mio_r_bool(record_parent_ids, "record_parent_ids"));
    if (out == NULL) mio_r_fail("slice");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_isosurface(SEXP mesh, SEXP array, SEXP isovalues, SEXP component,
                      SEXP record_parent_ids) {
    SEXP v = PROTECT(Rf_coerceVector(isovalues, REALSXP));
    R_xlen_t n = Rf_xlength(v);
    if (n < 1) {
        UNPROTECT(1);
        Rf_error("at least one isovalue is required");
    }
    mio_mesh *out = mio_isosurface(mio_r_mesh(mesh), mio_r_string(array, "array"), REAL(v),
                                   (int)n, mio_r_int(component, "component"),
                                   mio_r_bool(record_parent_ids, "record_parent_ids"));
    UNPROTECT(1);
    if (out == NULL) mio_r_fail("isosurface");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_gradient(SEXP mesh, SEXP array, SEXP op, SEXP method, SEXP location, SEXP output,
                    SEXP component, SEXP overwrite) {
    int64_t skipped = 0, fallback = 0;
    mio_mesh *out = mio_gradient(mio_r_mesh(mesh), mio_r_string(array, "array"),
                                 mio_r_opt_string(op), mio_r_opt_string(method),
                                 mio_r_opt_string(location), mio_r_opt_string(output),
                                 mio_r_int(component, "component"),
                                 mio_r_bool(overwrite, "overwrite"), &skipped, &fallback);
    if (out == NULL) mio_r_fail("gradient");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    /* R has no native int64, so the counters come back as doubles (exact well
       past any plausible cell count). */
    const char *names[] = {"mesh", "num_skipped", "num_fallback"};
    const SEXP values[] = {mo, NULL, NULL};
    const double numbers[] = {(double)skipped, (double)fallback};
    SEXP res = mio_r_report_list(3, names, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_hessian(SEXP mesh, SEXP array, SEXP method, SEXP location, SEXP output,
                   SEXP overwrite) {
    int64_t skipped = 0, fallback = 0;
    mio_mesh *out = mio_hessian(mio_r_mesh(mesh), mio_r_string(array, "array"),
                                mio_r_opt_string(method), mio_r_opt_string(location),
                                mio_r_opt_string(output), mio_r_bool(overwrite, "overwrite"),
                                &skipped, &fallback);
    if (out == NULL) mio_r_fail("hessian");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    /* R has no native int64, so the counters come back as doubles -- the
       mio_gradient wrapper's rule. */
    const char *names[] = {"mesh", "num_skipped", "num_fallback"};
    const SEXP values[] = {mo, NULL, NULL};
    const double numbers[] = {(double)skipped, (double)fallback};
    SEXP res = mio_r_report_list(3, names, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_estimate_error(SEXP mesh, SEXP array, SEXP method, SEXP marking, SEXP marking_value,
                          SEXP output, SEXP marked, SEXP overwrite) {
    double global_error = 0.0;
    int64_t skipped = 0, nmarked = 0;
    mio_mesh *out = mio_estimate_error(
        mio_r_mesh(mesh), mio_r_string(array, "array"), mio_r_opt_string(method),
        mio_r_opt_string(marking), mio_r_double(marking_value, "marking_value"),
        mio_r_opt_string(output), mio_r_opt_string(marked),
        mio_r_bool(overwrite, "overwrite"), &global_error, &skipped, &nmarked);
    if (out == NULL) mio_r_fail("estimate_error");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    /* R has no native int64, so the counters come back as doubles (exact well
       past any plausible cell count) -- the mio_gradient wrapper's rule. */
    const char *names[] = {"mesh", "global_error", "num_skipped", "num_marked"};
    const SEXP values[] = {mo, NULL, NULL, NULL};
    const double numbers[] = {global_error, (double)skipped, (double)nmarked};
    SEXP res = mio_r_report_list(4, names, values, numbers);
    UNPROTECT(1);
    return res;
}

/* Goes through mio_remesh_ex/mio_remesh_opts rather than the flat mio_remesh
   -- the R_mio_refine/mio_refine_ex precedent above, needed because
   mio_remesh is a flat C function with no room to grow (it already changed
   once). Unlike Julia/Fortran, R links the real C header, so mio_remesh_opts
   is used directly -- no hand-mirrored struct to keep in sync. */
SEXP R_mio_remesh(SEXP mesh, SEXP num_clusters, SEXP subdivide, SEXP subsample_ratio,
                  SEXP max_subdivide, SEXP max_iterations, SEXP max_repair_passes, SEXP metric,
                  SEXP gradation, SEXP preserve_boundary, SEXP max_anisotropy) {
    mio_remesh_opts opts;
    mio_remesh_report report;
    mio_mesh *out;

    mio_remesh_opts_init(&opts);
    opts.num_clusters = (int64_t)mio_r_int(num_clusters, "num_clusters");
    opts.subdivide = (int32_t)mio_r_int(subdivide, "subdivide");
    opts.max_subdivide = (int32_t)mio_r_int(max_subdivide, "max_subdivide");
    opts.subsample_ratio = mio_r_double(subsample_ratio, "subsample_ratio");
    opts.max_iterations = (int32_t)mio_r_int(max_iterations, "max_iterations");
    opts.max_repair_passes = (int32_t)mio_r_int(max_repair_passes, "max_repair_passes");
    opts.metric = mio_r_opt_string(metric);
    opts.gradation = mio_r_double(gradation, "gradation");
    opts.preserve_boundary = mio_r_bool(preserve_boundary, "preserve_boundary") ? 1 : 0;
    opts.max_anisotropy = mio_r_double(max_anisotropy, "max_anisotropy");

    out = mio_remesh_ex(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("remesh");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    /* R has no native int64, so the counters come back as doubles -- the
       mio_gradient wrapper's rule. */
    const char *names[] = {"mesh", "num_clusters", "num_iterations", "subdivide_applied",
                           "num_isolated_clusters", "num_non_manifold_vertices"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)report.num_clusters, (double)report.num_iterations,
                              (double)report.subdivide_applied,
                              (double)report.num_isolated_clusters,
                              (double)report.num_non_manifold_vertices};
    SEXP res = mio_r_report_list(6, names, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_compute_curvature(SEXP mesh, SEXP mean, SEXP gaussian, SEXP dual_area,
                             SEXP include_boundary, SEXP record_area, SEXP record_principal,
                             SEXP region) {
    mio_curvature_opts opts;
    mio_curvature_report report;
    mio_mesh *out;

    mio_curvature_opts_init(&opts);
    opts.mean = mio_r_bool(mean, "mean") ? 1 : 0;
    opts.gaussian = mio_r_bool(gaussian, "gaussian") ? 1 : 0;
    const char *da = mio_r_opt_string(dual_area);
    if (da != NULL && strcmp(da, "barycentric") == 0) {
        opts.dual_area = MIO_CURVATURE_BARYCENTRIC;
    } else if (da != NULL && strcmp(da, "mixed-voronoi") != 0) {
        Rf_error("meshio++: curvature: unknown dual area '%s' "
                 "(expected 'mixed-voronoi' or 'barycentric')", da);
    }
    opts.include_boundary = mio_r_bool(include_boundary, "include_boundary") ? 1 : 0;
    opts.record_area = mio_r_bool(record_area, "record_area") ? 1 : 0;
    opts.record_principal = mio_r_bool(record_principal, "record_principal") ? 1 : 0;
    opts.region = mio_r_opt_string(region);

    out = mio_compute_curvature(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("curvature");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP q = PROTECT(quality_list(&report.quality));
    /* R has no native int64, so the counters come back as doubles -- the
       mio_gradient wrapper's rule. */
    const char *names[] = {"mesh", "quality", "num_boundary", "num_isolated",
                           "num_degenerate", "total_angle_defect"};
    const SEXP values[] = {mo, q, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)report.num_boundary, (double)report.num_isolated,
                              (double)report.num_degenerate, report.total_angle_defect};
    SEXP res = mio_r_report_list(6, names, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_compute_normals(SEXP mesh, SEXP point_normals, SEXP cell_normals, SEXP weight,
                           SEXP split_angle, SEXP record_parent_ids, SEXP region) {
    mio_normals_opts opts;
    mio_normals_report report;
    mio_mesh *out;

    mio_normals_opts_init(&opts);
    opts.point_normals = mio_r_bool(point_normals, "point_normals") ? 1 : 0;
    opts.cell_normals = mio_r_bool(cell_normals, "cell_normals") ? 1 : 0;
    const char *w = mio_r_opt_string(weight);
    if (w != NULL && strcmp(w, "area") == 0) {
        opts.weight = MIO_SDF_WEIGHT_AREA;
    } else if (w != NULL && strcmp(w, "angle") != 0) {
        Rf_error("meshio++: normals: unknown weight '%s' (expected 'angle' or 'area')", w);
    }
    /* NULL means "no split"; a number is the crease angle in degrees. */
    if (split_angle != R_NilValue) {
        opts.split = 1;
        opts.split_angle = mio_r_double(split_angle, "split_angle");
    }
    opts.record_parent_ids = mio_r_bool(record_parent_ids, "record_parent_ids") ? 1 : 0;
    opts.region = mio_r_opt_string(region);

    out = mio_compute_normals(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("normals");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP q = PROTECT(quality_list(&report.quality));
    /* R has no native int64, so the counters come back as doubles. */
    const char *names[] = {"mesh", "quality", "num_isolated", "num_undefined",
                           "num_degenerate", "num_split_points", "num_added_points"};
    const SEXP values[] = {mo, q, NULL, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)report.num_isolated, (double)report.num_undefined,
                              (double)report.num_degenerate, (double)report.num_split_points,
                              (double)report.num_added_points};
    SEXP res = mio_r_report_list(7, names, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_feature_edges(SEXP mesh, SEXP feature_angle, SEXP feature, SEXP boundary,
                         SEXP non_manifold, SEXP inconsistent, SEXP region) {
    mio_feature_edges_opts opts;
    mio_feature_edges_report report;
    mio_mesh *out;

    mio_feature_edges_opts_init(&opts);
    opts.feature_angle = mio_r_double(feature_angle, "feature_angle");
    opts.feature = mio_r_bool(feature, "feature") ? 1 : 0;
    opts.boundary = mio_r_bool(boundary, "boundary") ? 1 : 0;
    opts.non_manifold = mio_r_bool(non_manifold, "non_manifold") ? 1 : 0;
    opts.inconsistent = mio_r_bool(inconsistent, "inconsistent") ? 1 : 0;
    opts.region = mio_r_opt_string(region);

    out = mio_feature_edges(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("feature_edges");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh", "num_feature", "num_boundary", "num_non_manifold",
                           "num_inconsistent"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)report.num_feature, (double)report.num_boundary,
                              (double)report.num_non_manifold, (double)report.num_inconsistent};
    SEXP res = mio_r_report_list(5, names, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_segment(SEXP mesh, SEXP feature_angle, SEXP normal_tolerance,
                   SEXP curvature_tolerance, SEXP min_patch_cells, SEXP region, SEXP prefix,
                   SEXP method, SEXP tolerance) {
    mio_segment_opts opts;
    mio_segment_report report;
    mio_mesh *out;

    mio_segment_opts_init(&opts);
    opts.feature_angle = mio_r_double(feature_angle, "feature_angle");
    opts.normal_tolerance = mio_r_double(normal_tolerance, "normal_tolerance");
    opts.curvature_tolerance = mio_r_double(curvature_tolerance, "curvature_tolerance");
    opts.min_patch_cells = (int64_t)mio_r_int(min_patch_cells, "min_patch_cells");
    opts.region = mio_r_opt_string(region);
    opts.prefix = mio_r_opt_string(prefix);
    {
        const char *m = mio_r_string(method, "method");
        if (strcmp(m, "angle") == 0)
            opts.method = MIO_SEGMENT_ANGLE;
        else if (strcmp(m, "primitive") == 0)
            opts.method = MIO_SEGMENT_PRIMITIVE;
        else
            Rf_error("meshio++: segment: method must be \"angle\" or \"primitive\"");
    }
    opts.tolerance = mio_r_double(tolerance, "tolerance");

    out = mio_segment(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("segment");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh", "num_patches", "num_merged", "num_feature_edges",
                           "num_boundary_edges"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)report.num_patches, (double)report.num_merged,
                              (double)report.num_feature_edges, (double)report.num_boundary_edges};
    SEXP res = mio_r_report_list(5, names, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_fit_primitives(SEXP mesh, SEXP patch_data, SEXP tolerance, SEXP relative,
                          SEXP min_inlier_fraction, SEXP types, SEXP max_iterations,
                          SEXP region, SEXP seeds, SEXP feature_angle, SEXP rim_angle,
                          SEXP output_region, SEXP used_points) {
    mio_fit_primitives_ex_opts opts;
    mio_fit_primitives_report report;
    mio_mesh *out;

    mio_fit_primitives_ex_opts_init(&opts);
    opts.base.patch_data = mio_r_string(patch_data, "patch_data");
    opts.base.tolerance = mio_r_double(tolerance, "tolerance");
    opts.base.relative = mio_r_bool(relative, "relative") ? 1 : 0;
    opts.base.min_inlier_fraction = mio_r_double(min_inlier_fraction, "min_inlier_fraction");
    opts.base.types = (uint32_t)mio_r_int(types, "types");
    opts.base.max_iterations = (int64_t)mio_r_int(max_iterations, "max_iterations");
    opts.region = mio_r_opt_string(region);
    if (seeds != R_NilValue) {
        if (!Rf_isReal(seeds) || XLENGTH(seeds) % 3 != 0)
            Rf_error("`seeds` must hold x, y, z triples");
        opts.seeds = REAL(seeds);
        opts.num_seeds = (int64_t)(XLENGTH(seeds) / 3);
    }
    opts.feature_angle = mio_r_double(feature_angle, "feature_angle");
    opts.rim_angle = mio_r_double(rim_angle, "rim_angle");
    opts.output_region = mio_r_opt_string(output_region);
    opts.used_points = (int32_t)mio_r_int(used_points, "used_points");

    out = mio_fit_primitives_ex(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("fit_primitives");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP warnings = PROTECT(mio_r_getbuffer(mio_fit_primitives_warnings, "fit_primitives"));
    const char *names[] = {"mesh",      "num_patches", "num_plane",    "num_sphere",
                           "num_cylinder", "num_cone", "num_torus",    "num_unfitted",
                           "num_feature_cells", "warnings"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, warnings};
    const double numbers[] = {(double)report.num_patches, (double)report.num_plane,
                              (double)report.num_sphere, (double)report.num_cylinder,
                              (double)report.num_cone, (double)report.num_torus,
                              (double)report.num_unfitted, (double)report.num_feature_cells};
    SEXP res = mio_r_report_list(10, names, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_has_occt(void) {
    return Rf_ScalarLogical(mio_has_occt() != 0);
}

/* A RANSAC seed: NULL takes the default; a number must be a whole number in
 * [0, 2^64), held exactly below 2^53 (R numbers are doubles). */
static uint64_t ops_refine_seed(SEXP s) {
    double v;
    if (s == R_NilValue)
        return MIO_REFINE_DEFAULT_SEED;
    if (Rf_length(s) != 1 || (TYPEOF(s) != REALSXP && TYPEOF(s) != INTSXP))
        Rf_error("`ransac_seed` must be NULL or a single number");
    v = Rf_asReal(s);
    if (ISNAN(v) || v < 0.0 || v >= 18446744073709551616.0 || v != floor(v))
        Rf_error("`ransac_seed` must be a whole number in [0, 2^64)");
    return (uint64_t)v;
}

/* The mio_refine_stop code as a string. */
static const char* ops_refine_stop_name(int32_t stop) {
    switch (stop) {
        case MIO_REFINE_UNCHANGED:
            return "unchanged";
        case MIO_REFINE_SMALL_GAIN:
            return "small_gain";
        case MIO_REFINE_ITERATION_CAP:
            return "iteration_cap";
        case MIO_REFINE_ENERGY_ROSE:
            return "energy_rose";
        default:
            return "unknown";
    }
}

SEXP R_mio_refine_patches(SEXP mesh, SEXP tolerance, SEXP relative, SEXP normal_weight,
                          SEXP epsilon, SEXP feature_angle, SEXP max_iterations,
                          SEXP min_patch_cells, SEXP ransac_seed, SEXP types, SEXP allow_split,
                          SEXP allow_merge, SEXP prefix) {
    mio_refine_patches_opts opts;
    mio_refine_patches_report report;
    mio_mesh* out;

    mio_refine_patches_opts_init(&opts);
    opts.tolerance = mio_r_double(tolerance, "tolerance");
    opts.relative = mio_r_bool(relative, "relative") ? 1 : 0;
    opts.normal_weight = mio_r_double(normal_weight, "normal_weight");
    opts.epsilon = mio_r_double(epsilon, "epsilon");
    opts.feature_angle = mio_r_double(feature_angle, "feature_angle");
    opts.max_iterations = mio_r_int64(max_iterations, "max_iterations");
    opts.min_patch_cells = mio_r_int64(min_patch_cells, "min_patch_cells");
    opts.ransac_seed = ops_refine_seed(ransac_seed);
    opts.types = (uint32_t)mio_r_int(types, "types");
    opts.allow_split = mio_r_bool(allow_split, "allow_split") ? 1 : 0;
    opts.allow_merge = mio_r_bool(allow_merge, "allow_merge") ? 1 : 0;
    opts.prefix = mio_r_opt_string(prefix);

    out = mio_refine_patches(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL)
        mio_r_fail("refine_patches");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP st = PROTECT(Rf_mkString(ops_refine_stop_name(report.stop)));
    const char* names[] = {"mesh",
                           "num_patches",
                           "num_iterations",
                           "num_splits",
                           "num_merges",
                           "num_absorbed",
                           "num_primitive",
                           "num_freeform",
                           "num_facets",
                           "tolerance",
                           "explained_area_fraction",
                           "initial_explained_area_fraction",
                           "hausdorff_max",
                           "hausdorff_rms",
                           "initial_energy",
                           "energy",
                           "stop"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, st};
    const double numbers[] = {(double)report.num_patches, (double)report.num_iterations,
                              (double)report.num_splits, (double)report.num_merges,
                              (double)report.num_absorbed, (double)report.num_primitive,
                              (double)report.num_freeform, (double)report.num_facets,
                              report.tolerance, report.explained_area_fraction,
                              report.initial_explained_area_fraction, report.hausdorff_max,
                              report.hausdorff_rms, report.initial_energy, report.energy};
    SEXP res = mio_r_report_list(17, names, values, numbers);
    UNPROTECT(2);
    return res;
}

/* The 0-based position of a keyword option in a NULL-terminated list. */
static int32_t ops_fit_cad_keyword(SEXP value, const char* name, const char* const* words) {
    const char* v = mio_r_string(value, name);
    int32_t i;
    for (i = 0; words[i] != NULL; ++i)
        if (strcmp(v, words[i]) == 0)
            return i;
    Rf_error("meshio++: fit_cad: %s: unknown value \"%s\"", name, v);
    return 0;
}

SEXP R_mio_fit_cad(SEXP mesh, SEXP step_path, SEXP tolerance, SEXP exact_edges, SEXP types,
                   SEXP mesh_deflection, SEXP max_demotion_rounds, SEXP freeform,
                   SEXP freeform_degree, SEXP freeform_max_control, SEXP freeform_split,
                   SEXP fillets, SEXP threads, SEXP merge_faces, SEXP backend,
                   SEXP refine_segmentation, SEXP max_refine_iterations, SEXP refine_epsilon,
                   SEXP normal_weight, SEXP allow_merge, SEXP allow_split, SEXP ransac_seed,
                   SEXP refine_tolerance, SEXP relative_tolerance, SEXP freeform_general,
                   SEXP freeform_tolerance, SEXP freeform_smoothing, SEXP max_freeform_poles,
                   SEXP max_freeform_iterations, SEXP max_freeform_splits, SEXP sew_tolerance,
                   SEXP extend_surfaces, SEXP snap_junctions, SEXP face_deviation,
                   SEXP freeform_mode, SEXP freeform_parameterization, SEXP extraordinary_vertex,
                   SEXP local_refine, SEXP freeform_g1_tolerance_deg,
                   SEXP freeform_max_control_points, SEXP layout_target_patch_size,
                   SEXP layout_max_valence_deviation, SEXP resegment_rounds,
                   SEXP layout_snap_to_features, SEXP validate_by_surface_sampling,
                   SEXP constrain_primitives) {
    mio_fit_cad_opts opts;
    mio_fit_cad_opts_ex opts_ex;
    mio_fit_cad_opts_ex2 opts_ex2;
    mio_fit_cad_report report;
    mio_fit_cad_report_ex report_ex;
    mio_fit_cad_report_ex2 report_ex2;
    mio_fit_cad_report_ex3 report_ex3;
    mio_mesh *out;
    static const char* const modes[] = {"legacy", "adaptive", "network", "auto", NULL};
    static const char* const params[] = {"auto", "projection", "lscm", "harmonic", NULL};
    static const char* const evs[] = {"g0", "g1", "g1x", "g1xx", NULL};
    static const char* const refines[] = {"knots", "thb", NULL};

    mio_fit_cad_opts_init(&opts);
    mio_fit_cad_opts_ex_init(&opts_ex);
    mio_fit_cad_opts_ex2_init(&opts_ex2);
    opts_ex2.freeform_mode = ops_fit_cad_keyword(freeform_mode, "freeform_mode", modes);
    opts_ex2.freeform_parameterization =
        ops_fit_cad_keyword(freeform_parameterization, "freeform_parameterization", params);
    opts_ex2.extraordinary_vertex =
        ops_fit_cad_keyword(extraordinary_vertex, "extraordinary_vertex", evs);
    opts_ex2.local_refine = ops_fit_cad_keyword(local_refine, "local_refine", refines);
    opts_ex2.freeform_g1_tolerance_deg =
        mio_r_double(freeform_g1_tolerance_deg, "freeform_g1_tolerance_deg");
    opts_ex2.freeform_max_control_points =
        mio_r_int64(freeform_max_control_points, "freeform_max_control_points");
    opts_ex2.layout_target_patch_size =
        mio_r_double(layout_target_patch_size, "layout_target_patch_size");
    opts_ex2.layout_max_valence_deviation =
        mio_r_int64(layout_max_valence_deviation, "layout_max_valence_deviation");
    opts_ex2.resegment_rounds = mio_r_int64(resegment_rounds, "resegment_rounds");
    opts_ex2.layout_snap_to_features =
        mio_r_bool(layout_snap_to_features, "layout_snap_to_features") ? 1 : 0;
    opts_ex2.validate_by_surface_sampling =
        mio_r_bool(validate_by_surface_sampling, "validate_by_surface_sampling") ? 1 : 0;
    opts_ex2.constrain_primitives =
        mio_r_bool(constrain_primitives, "constrain_primitives") ? 1 : 0;
    opts_ex.refine_segmentation = mio_r_bool(refine_segmentation, "refine_segmentation") ? 1 : 0;
    opts_ex.refine.max_iterations = mio_r_int64(max_refine_iterations, "max_refine_iterations");
    opts_ex.refine.epsilon = mio_r_double(refine_epsilon, "refine_epsilon");
    opts_ex.refine.normal_weight = mio_r_double(normal_weight, "normal_weight");
    opts_ex.refine.allow_merge = mio_r_bool(allow_merge, "allow_merge") ? 1 : 0;
    opts_ex.refine.allow_split = mio_r_bool(allow_split, "allow_split") ? 1 : 0;
    opts_ex.refine.ransac_seed = ops_refine_seed(ransac_seed);
    opts_ex.refine.tolerance = mio_r_double(refine_tolerance, "refine_tolerance");
    opts_ex.refine.relative = mio_r_bool(relative_tolerance, "relative_tolerance") ? 1 : 0;
    opts_ex.freeform_general = mio_r_bool(freeform_general, "freeform_general") ? 1 : 0;
    opts_ex.freeform_tolerance = mio_r_double(freeform_tolerance, "freeform_tolerance");
    opts_ex.freeform_smoothing = mio_r_double(freeform_smoothing, "freeform_smoothing");
    opts_ex.max_freeform_poles = mio_r_int64(max_freeform_poles, "max_freeform_poles");
    opts_ex.max_freeform_iterations =
        mio_r_int64(max_freeform_iterations, "max_freeform_iterations");
    opts_ex.max_freeform_splits = mio_r_int64(max_freeform_splits, "max_freeform_splits");
    opts_ex.sew_tolerance = mio_r_double(sew_tolerance, "sew_tolerance");
    opts_ex.extend_surfaces = mio_r_bool(extend_surfaces, "extend_surfaces") ? 1 : 0;
    opts_ex.snap_junctions = mio_r_bool(snap_junctions, "snap_junctions") ? 1 : 0;
    opts_ex.face_deviation = mio_r_bool(face_deviation, "face_deviation") ? 1 : 0;
    opts.tolerance = mio_r_double(tolerance, "tolerance");
    opts.exact_edges = mio_r_bool(exact_edges, "exact_edges") ? 1 : 0;
    opts.types = (uint32_t)mio_r_int(types, "types");
    opts.mesh_deflection = mio_r_double(mesh_deflection, "mesh_deflection");
    opts.max_demotion_rounds = (int64_t)mio_r_int(max_demotion_rounds, "max_demotion_rounds");
    opts.freeform = mio_r_bool(freeform, "freeform") ? 1 : 0;
    opts.freeform_degree = (int32_t)mio_r_int(freeform_degree, "freeform_degree");
    opts.freeform_max_control = (int64_t)mio_r_int(freeform_max_control, "freeform_max_control");
    opts.freeform_split = mio_r_bool(freeform_split, "freeform_split") ? 1 : 0;
    opts.fillets = mio_r_bool(fillets, "fillets") ? 1 : 0;
    opts.threads = mio_r_bool(threads, "threads") ? 1 : 0;
    opts.merge_faces = mio_r_bool(merge_faces, "merge_faces") ? 1 : 0;
    {
        const char *b = mio_r_string(backend, "backend");
        if (strcmp(b, "auto") == 0)
            opts.backend = MIO_FIT_CAD_AUTO;
        else if (strcmp(b, "occt") == 0)
            opts.backend = MIO_FIT_CAD_OCCT;
        else if (strcmp(b, "native") == 0)
            opts.backend = MIO_FIT_CAD_NATIVE;
        else
            Rf_error("meshio++: fit_cad: backend must be \"auto\", \"occt\" or \"native\"");
    }

    out = mio_fit_cad3(mio_r_mesh(mesh), mio_r_string(step_path, "step_path"), &opts, &opts_ex,
                       &opts_ex2, &report, &report_ex, &report_ex2, &report_ex3);
    if (out == NULL) mio_r_fail("fit_cad");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP valid = PROTECT(Rf_ScalarLogical(report.valid != 0));
    SEXP engine = PROTECT(
        Rf_mkString(report_ex.backend == MIO_FIT_CAD_NATIVE ? "native" : "occt"));
    SEXP refined = PROTECT(Rf_ScalarLogical(report_ex2.refined != 0));
    SEXP closed = PROTECT(Rf_ScalarLogical(report_ex2.closed != 0));
    SEXP kept = PROTECT(Rf_ScalarLogical(report_ex2.refinement_kept != 0));
    SEXP checked = PROTECT(Rf_ScalarLogical(report_ex3.brep_checked != 0));
    SEXP mode = PROTECT(Rf_mkString(report_ex3.freeform_mode >= 0 && report_ex3.freeform_mode <= 3
                                        ? modes[report_ex3.freeform_mode]
                                        : "legacy"));
    SEXP g1 = PROTECT(Rf_ScalarReal(report_ex3.max_g1_angle_deg));
    const char* names[] = {"mesh",
                           "valid",
                           "tolerance",
                           "volume",
                           "skin_volume",
                           "volume_ratio",
                           "max_deviation",
                           "num_faces",
                           "num_plane",
                           "num_sphere",
                           "num_cylinder",
                           "num_cone",
                           "num_torus",
                           "num_faceted_faces",
                           "num_patches",
                           "num_demoted",
                           "num_exact_edges",
                           "num_approx_edges",
                           "num_solids",
                           "num_bspline",
                           "num_split_patches",
                           "num_fillet",
                           "num_helical",
                           "num_variable_fillet",
                           "num_blend",
                           "num_rounded_thread",
                           "num_vertex_blend",
                           "num_chamfer",
                           "num_conic_blend",
                           "num_tapered_thread",
                           "num_variable_pitch_thread",
                           "backend",
                           "refined",
                           "refine_iterations",
                           "num_refine_splits",
                           "num_refine_merges",
                           "initial_explained_area_fraction",
                           "explained_area_fraction",
                           "closed",
                           "num_freeform_general",
                           "num_freeform_refinements",
                           "avg_freeform_poles",
                           "num_extended_edges",
                           "num_snapped_junctions",
                           "hausdorff_max",
                           "hausdorff_rms",
                           "refinement_kept",
                           "brep_checked",
                           "freeform_mode",
                           "facet_area_fraction",
                           "num_freeform_adaptive",
                           "num_network_patches",
                           "num_network_edges",
                           "num_extraordinary_vertices",
                           "num_resegmented",
                           "num_escalations",
                           "max_g1_angle_deg",
                           "max_freeform_distance",
                           "rms_freeform_distance",
                           "num_face_records",
                           "num_constrained_primitives",
                           "num_constraint_groups"};
    const SEXP values[] = {mo, valid, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, engine, refined, NULL, NULL,
                           NULL, NULL, NULL, closed, NULL, NULL, NULL, NULL, NULL, NULL, NULL, kept,
                           checked, mode, NULL, NULL, NULL, NULL, NULL, NULL, NULL, g1, NULL, NULL,
                           NULL, NULL, NULL};
    const double numbers[] = {report.tolerance, report.volume, report.skin_volume,
                              report.volume_ratio, report.max_deviation, (double)report.num_faces,
                              (double)report.num_plane, (double)report.num_sphere,
                              (double)report.num_cylinder, (double)report.num_cone,
                              (double)report.num_torus, (double)report.num_faceted_faces,
                              (double)report.num_patches, (double)report.num_demoted,
                              (double)report.num_exact_edges, (double)report.num_approx_edges,
                              (double)report.num_solids, (double)report.num_bspline,
                              (double)report.num_split_patches, (double)report.num_fillet,
                              (double)report.num_helical, (double)report.num_variable_fillet,
                              (double)report.num_blend, (double)report.num_rounded_thread,
                              (double)report_ex.num_vertex_blend, (double)report_ex.num_chamfer,
                              (double)report_ex.num_conic_blend,
                              (double)report_ex.num_tapered_thread,
                              (double)report_ex.num_variable_pitch_thread,
                              (double)report_ex2.refine_iterations,
                              (double)report_ex2.num_refine_splits,
                              (double)report_ex2.num_refine_merges,
                              report_ex2.initial_explained_area_fraction,
                              report_ex2.explained_area_fraction,
                              (double)report_ex2.num_freeform_general,
                              (double)report_ex2.num_freeform_refinements,
                              report_ex2.avg_freeform_poles, (double)report_ex2.num_extended_edges,
                              (double)report_ex2.num_snapped_junctions, report_ex2.hausdorff_max,
                              report_ex2.hausdorff_rms, report_ex3.facet_area_fraction,
                              (double)report_ex3.num_freeform_adaptive,
                              (double)report_ex3.num_network_patches,
                              (double)report_ex3.num_network_edges,
                              (double)report_ex3.num_extraordinary_vertices,
                              (double)report_ex3.num_resegmented,
                              (double)report_ex3.num_escalations, report_ex3.max_freeform_distance,
                              report_ex3.rms_freeform_distance, (double)report_ex3.num_face_records,
                              (double)report_ex3.num_constrained_primitives,
                              (double)report_ex3.num_constraint_groups};
    SEXP res = mio_r_report_list(62, names, values, numbers);
    UNPROTECT(9);
    return res;
}

SEXP R_mio_region_adjacency(SEXP mesh, SEXP regions) {
    R_xlen_t n = 0;
    mio_region_selector *selectors = NULL;
    if (regions != R_NilValue) {
        if (TYPEOF(regions) != STRSXP)
            Rf_error("`regions` must be NULL or a character vector");
        n = Rf_xlength(regions);
        if (n > 0) selectors = (mio_region_selector *)R_alloc((size_t)n, sizeof(*selectors));
        for (R_xlen_t i = 0; i < n; ++i) {
            const char *name = CHAR(STRING_ELT(regions, i));
            mio_region_selector_init(&selectors[i], name);
            selectors[i].kind = MIO_REGION_CELL;
        }
    }
    mio_mesh *out = mio_region_adjacency(mio_r_mesh(mesh), selectors, (int64_t)n);
    if (out == NULL) mio_r_fail("region_adjacency");
    return mio_r_wrap_mesh(out);
}

/* These ids cross the flat C ABI in 0-based `(cell, local_facet)` pairs.
 * R's region convention shifts only the cell column; facet ordinals stay 0-based. */
static SEXP mio_r_side_matrix(const int64_t *entries, int64_t count) {
    if (count < 0 || count > INT_MAX / 2) Rf_error("interface result is too large for an R matrix");
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, 2, (int)count));
    double *dst = REAL(out);
    for (int64_t i = 0; i < 2 * count; ++i)
        dst[i] = (i % 2 == 0) ? (double)(entries[i] + 1) : (double)entries[i];
    UNPROTECT(1);
    return out;
}

SEXP R_mio_find_interface(SEXP mesh_a, SEXP region_a, SEXP region_b, SEXP mesh_b,
                          SEXP mode, SEXP master, SEXP gap_tolerance,
                          SEXP angle_tolerance, SEXP overlap_tolerance) {
    mio_find_interface_opts opts;
    mio_find_interface_opts_init(&opts);
    const char *mode_name = mio_r_opt_string(mode);
    if (mode_name != NULL && mode_name[0] != '\0') {
        if (strcmp(mode_name, "conforming") == 0) opts.mode = MIO_INTERFACE_CONFORMING;
        else if (strcmp(mode_name, "proximity") == 0) opts.mode = MIO_INTERFACE_PROXIMITY;
        else Rf_error("`mode` must be 'conforming' or 'proximity'");
    }
    const char *master_name = mio_r_opt_string(master);
    if (master_name != NULL && master_name[0] != '\0') {
        if (strcmp(master_name, "a") == 0) opts.master = MIO_INTERFACE_MASTER_A;
        else if (strcmp(master_name, "b") == 0) opts.master = MIO_INTERFACE_MASTER_B;
        else Rf_error("`master` must be 'a' or 'b'");
    }
    opts.gap_tolerance = mio_r_double(gap_tolerance, "gap_tolerance");
    opts.angle_tolerance = mio_r_double(angle_tolerance, "angle_tolerance");
    opts.overlap_tolerance = mio_r_double(overlap_tolerance, "overlap_tolerance");

    mio_region_selector selector_a, selector_b;
    mio_region_selector_init(&selector_a, mio_r_string(region_a, "region_a"));
    mio_region_selector_init(&selector_b, mio_r_string(region_b, "region_b"));
    selector_a.kind = selector_b.kind = MIO_REGION_CELL;
    const mio_mesh *b = mesh_b == R_NilValue ? NULL : mio_r_mesh(mesh_b);
    mio_find_interface_result *result = mio_find_interface(
        mio_r_mesh(mesh_a), &selector_a, b, &selector_b, &opts);
    if (result == NULL) mio_r_fail("find_interface");

    int64_t pairs = 0, unmatched_a = 0, unmatched_b = 0, count_a = 0, count_b = 0;
    double area = 0.0, max_gap = 0.0;
    if (mio_find_interface_result_report(result, &pairs, &area, &max_gap,
                                         &unmatched_a, &unmatched_b) != MIO_OK) {
        mio_find_interface_result_free(result);
        mio_r_fail("find_interface");
    }
    const int64_t *side_a = mio_find_interface_result_side_a(result, &count_a);
    const int64_t *side_b = mio_find_interface_result_side_b(result, &count_b);
    mio_mesh *out_mesh = mio_find_interface_result_take_mesh(result);
    if (out_mesh == NULL) {
        mio_find_interface_result_free(result);
        mio_r_fail("find_interface");
    }

    SEXP mesh_value = PROTECT(mio_r_wrap_mesh(out_mesh));
    SEXP side_a_value = PROTECT(mio_r_side_matrix(side_a, count_a));
    SEXP side_b_value = PROTECT(mio_r_side_matrix(side_b, count_b));
    mio_find_interface_result_free(result);
    const char *names[] = {"mesh", "num_pairs", "area", "max_gap", "unmatched_a",
                           "unmatched_b", "side_a", "side_b"};
    const SEXP values[] = {mesh_value, NULL, NULL, NULL, NULL, NULL, side_a_value, side_b_value};
    const double numbers[] = {(double)pairs, area, max_gap, (double)unmatched_a,
                              (double)unmatched_b};
    SEXP out = mio_r_report_list(8, names, values, numbers);
    UNPROTECT(3);
    return out;
}

static SEXP mio_r_copy_i64(const int64_t *data, int64_t count) {
    if (count == 0) return Rf_allocVector(REALSXP, 0);
    return mio_r_copy_as_real(data, MIO_INT64, (R_xlen_t)count);
}

static SEXP mio_r_copy_f64(const double *data, int64_t count) {
    if (count == 0) return Rf_allocVector(REALSXP, 0);
    return mio_r_copy_as_real(data, MIO_FLOAT64, (R_xlen_t)count);
}

static SEXP mio_r_copy_xyz(const double *data, int64_t count) {
    if (count < 0 || count > INT_MAX / 3) Rf_error("contact result is too large for an R matrix");
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, 3, (int)count));
    if (count > 0) memcpy(REAL(out), data, (size_t)(3 * count) * sizeof(double));
    UNPROTECT(1);
    return out;
}

SEXP R_mio_contact_pairs(SEXP slave_mesh, SEXP slave_points, SEXP master_mesh,
                         SEXP master_cells, SEXP tolerance, SEXP require_complete) {
    mio_contact_pairs_opts opts;
    mio_contact_pairs_opts_init(&opts);
    opts.tolerance = mio_r_double(tolerance, "tolerance");
    opts.require_complete = mio_r_bool(require_complete, "require_complete") ? 1 : 0;
    mio_region_selector point_selector, cell_selector;
    mio_region_selector_init(&point_selector, mio_r_string(slave_points, "slave_points"));
    mio_region_selector_init(&cell_selector, mio_r_string(master_cells, "master_cells"));
    point_selector.kind = MIO_REGION_POINT;
    cell_selector.kind = MIO_REGION_CELL;
    const mio_mesh *master = master_mesh == R_NilValue ? NULL : mio_r_mesh(master_mesh);
    mio_contact_pairs_result *result = mio_contact_pairs(
        mio_r_mesh(slave_mesh), &point_selector, master, &cell_selector, &opts);
    if (result == NULL) mio_r_fail("contact_pairs");
    int64_t count = 0, unmatched_count = 0;
    if (mio_contact_pairs_result_info(result, &count, &unmatched_count) != MIO_OK) {
        mio_contact_pairs_result_free(result);
        mio_r_fail("contact_pairs");
    }

    SEXP slave = PROTECT(mio_r_shift_map(mio_contact_pairs_slave_point(result), count));
    SEXP master_cell = PROTECT(mio_r_shift_map(mio_contact_pairs_master_cell(result), count));
    SEXP master_facet = PROTECT(mio_r_copy_i64(mio_contact_pairs_master_facet(result), count));
    SEXP master_subfacet = PROTECT(mio_r_copy_i64(mio_contact_pairs_master_subfacet(result), count));
    SEXP local_coordinates = PROTECT(mio_r_copy_xyz(
        mio_contact_pairs_local_coordinates(result), count));
    SEXP closest_point = PROTECT(mio_r_copy_xyz(mio_contact_pairs_closest_point(result), count));
    SEXP gap = PROTECT(mio_r_copy_f64(mio_contact_pairs_gap(result), count));
    SEXP normal = PROTECT(mio_r_copy_xyz(mio_contact_pairs_normal(result), count));
    SEXP unmatched = PROTECT(mio_r_shift_map(mio_contact_pairs_unmatched(result), unmatched_count));
    mio_contact_pairs_result_free(result);
    const char *names[] = {"slave_point", "master_cell", "master_facet", "master_subfacet",
                           "local_coordinates", "closest_point", "gap", "normal", "unmatched"};
    SEXP values[] = {slave, master_cell, master_facet, master_subfacet, local_coordinates,
                     closest_point, gap, normal, unmatched};
    SEXP out = PROTECT(mio_r_named_list(9, names, values));
    UNPROTECT(10);
    return out;
}

SEXP R_mio_split_interface(SEXP mesh, SEXP side, SEXP add_cohesive) {
    mio_split_interface_opts opts;
    mio_split_interface_opts_init(&opts);
    opts.add_cohesive = mio_r_bool(add_cohesive, "add_cohesive") ? 1 : 0;
    mio_region_selector selector;
    mio_region_selector_init(&selector, mio_r_string(side, "side"));
    selector.kind = MIO_REGION_SIDE;
    int64_t duplicated = 0, cohesive = 0;
    mio_mesh *out_mesh = mio_split_interface(mio_r_mesh(mesh), &selector, &opts,
                                             &duplicated, &cohesive);
    if (out_mesh == NULL) mio_r_fail("split_interface");
    SEXP mesh_value = PROTECT(mio_r_wrap_mesh(out_mesh));
    const char *names[] = {"mesh", "num_duplicated_points", "num_cohesive_cells"};
    const SEXP values[] = {mesh_value, NULL, NULL};
    const double numbers[] = {(double)duplicated, (double)cohesive};
    SEXP out = mio_r_report_list(3, names, values, numbers);
    UNPROTECT(1);
    return out;
}

SEXP R_mio_check_quality(SEXP mesh, SEXP require, SEXP max_inverted, SEXP max_degenerate) {
    mio_quality_gate_report r;
    char buf[8192];
    mio_r_check(mio_check_quality(mio_r_mesh(mesh), mio_r_opt_string(require),
                                  (int64_t)mio_r_double(max_inverted, "max_inverted"),
                                  (int64_t)mio_r_double(max_degenerate, "max_degenerate"), &r,
                                  buf, (int64_t)sizeof(buf)),
                "check_quality");
    SEXP v0 = PROTECT(Rf_ScalarLogical(r.passed != 0));
    SEXP v6 = PROTECT(Rf_mkString(buf));
    const char *names[] = {"passed", "num_checks", "num_failed", "num_cells",
                           "num_inverted", "num_degenerate", "summary"};
    const SEXP values[] = {v0, NULL, NULL, NULL, NULL, NULL, v6};
    const double numbers[] = {(double)r.num_checks, (double)r.num_failed, (double)r.num_cells,
                              (double)r.num_inverted, (double)r.num_degenerate};
    SEXP res = mio_r_report_list(7, names, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_hausdorff_distance(SEXP a, SEXP b, SEXP face_samples, SEXP region_a,
                              SEXP region_b) {
    mio_hausdorff_opts opts;
    mio_hausdorff_report r;

    mio_hausdorff_opts_init(&opts);
    opts.face_samples = (int64_t)mio_r_double(face_samples, "face_samples");
    opts.region_a = mio_r_opt_string(region_a);
    opts.region_b = mio_r_opt_string(region_b);
    mio_r_check(mio_hausdorff_distance(mio_r_mesh(a), mio_r_mesh(b), &opts, &r),
                "hausdorff_distance");
    SEXP wa = PROTECT(Rf_allocVector(REALSXP, 3));
    SEXP wb = PROTECT(Rf_allocVector(REALSXP, 3));
    for (int k = 0; k < 3; ++k) {
        REAL(wa)[k] = r.worst_point_a[k];
        REAL(wb)[k] = r.worst_point_b[k];
    }
    const char *names[] = {"distance",    "a_to_b",        "b_to_a",        "mean_a_to_b",
                           "rms_a_to_b",  "mean_b_to_a",   "rms_b_to_a",    "num_samples_a",
                           "num_samples_b", "worst_point_a", "worst_point_b"};
    const SEXP values[] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, wa, wb};
    const double numbers[] = {r.distance, r.a_to_b, r.b_to_a, r.mean_a_to_b, r.rms_a_to_b,
                              r.mean_b_to_a, r.rms_b_to_a, (double)r.num_samples_a,
                              (double)r.num_samples_b};
    SEXP res = mio_r_report_list(11, names, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_repair(SEXP mesh, SEXP fix_orientation, SEXP orient_outward, SEXP fill_holes,
                  SEXP split_non_manifold, SEXP max_hole_edges, SEXP weld_tolerance,
                  SEXP record_provenance) {
    mio_repair_opts opts;
    mio_repair_report report;
    mio_mesh *out;

    mio_repair_opts_init(&opts);
    opts.fix_orientation = mio_r_bool(fix_orientation, "fix_orientation") ? 1 : 0;
    opts.orient_outward = mio_r_bool(orient_outward, "orient_outward") ? 1 : 0;
    opts.fill_holes = mio_r_bool(fill_holes, "fill_holes") ? 1 : 0;
    opts.split_non_manifold = mio_r_bool(split_non_manifold, "split_non_manifold") ? 1 : 0;
    opts.max_hole_edges = (int64_t)mio_r_double(max_hole_edges, "max_hole_edges");
    opts.weld_tolerance = mio_r_double(weld_tolerance, "weld_tolerance");
    opts.record_provenance = mio_r_bool(record_provenance, "record_provenance") ? 1 : 0;

    out = mio_repair(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("repair");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP qb = PROTECT(quality_list(&report.quality_before));
    SEXP qa = PROTECT(quality_list(&report.quality_after));
    /* R has no native int64, so the counters come back as doubles. */
    const int64_t counters[] = {report.num_flipped,        report.num_components,
                                report.largest_component,  report.num_oriented_outward,
                                report.num_unorientable,   report.num_vertices_split,
                                report.num_holes_detected, report.num_holes_filled,
                                report.num_holes_skipped,  report.num_faces_added,
                                report.num_points_added,   report.points_welded};
    const char *names[] = {"mesh", "quality_before", "quality_after", "num_flipped",
                           "num_components", "largest_component", "num_oriented_outward",
                           "num_unorientable", "num_vertices_split", "num_holes_detected",
                           "num_holes_filled", "num_holes_skipped", "num_faces_added",
                           "num_points_added", "points_welded"};
    SEXP values[15];
    values[0] = mo;
    values[1] = qb;
    values[2] = qa;
    for (int k = 0; k < 12; ++k) values[3 + k] = PROTECT(Rf_ScalarReal((double)counters[k]));
    SEXP res = PROTECT(mio_r_named_list(15, names, values));
    UNPROTECT(16);
    return res;
}

SEXP R_mio_shrinkwrap(SEXP mesh, SEXP target, SEXP offset, SEXP max_distance, SEXP weights,
                      SEXP target_region, SEXP normal_weight, SEXP record_distance,
                      SEXP record_closest_cell) {
    mio_shrinkwrap_opts opts;
    mio_shrinkwrap_report report;
    mio_mesh *out;

    mio_shrinkwrap_opts_init(&opts);
    opts.offset = mio_r_double(offset, "offset");
    opts.max_distance = mio_r_double(max_distance, "max_distance");
    opts.weights = mio_r_opt_string(weights);
    opts.target_region = mio_r_opt_string(target_region);
    const char *nw = mio_r_opt_string(normal_weight);
    if (nw != NULL && strcmp(nw, "area") == 0) {
        opts.normal_weight = MIO_SDF_WEIGHT_AREA;
    } else if (nw != NULL && strcmp(nw, "angle") != 0) {
        Rf_error("meshio++: shrinkwrap: unknown normal weight '%s' (expected 'angle' or 'area')",
                 nw);
    }
    opts.record_distance = mio_r_bool(record_distance, "record_distance") ? 1 : 0;
    opts.record_closest_cell = mio_r_bool(record_closest_cell, "record_closest_cell") ? 1 : 0;

    out = mio_shrinkwrap(mio_r_mesh(mesh), mio_r_mesh(target), &opts, &report);
    if (out == NULL) mio_r_fail("shrinkwrap");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP q = PROTECT(quality_list(&report.quality));
    const char *names[] = {"mesh", "quality", "num_projected", "num_missed", "num_skipped",
                           "max_displacement"};
    const SEXP values[] = {mo, q, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)report.num_projected, (double)report.num_missed,
                              (double)report.num_skipped, report.max_displacement};
    SEXP res = mio_r_report_list(6, names, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_sobolev_deform(SEXP mesh, SEXP array, SEXP length_scale, SEXP fixed_points_array,
                          SEXP fix_boundary, SEXP record_filtered, SEXP max_iterations,
                          SEXP tolerance) {
    mio_sobolev_opts opts;
    mio_sobolev_report report;
    mio_mesh *out;

    mio_sobolev_opts_init(&opts);
    opts.array = mio_r_opt_string(array);
    opts.length_scale = mio_r_double(length_scale, "length_scale");
    opts.fixed_points_array = mio_r_opt_string(fixed_points_array);
    opts.fix_boundary = mio_r_bool(fix_boundary, "fix_boundary") ? 1 : 0;
    opts.record_filtered = mio_r_bool(record_filtered, "record_filtered") ? 1 : 0;
    opts.max_iterations = mio_r_int(max_iterations, "max_iterations");
    opts.tolerance = mio_r_double(tolerance, "tolerance");

    out = mio_sobolev_deform(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("sobolev_deform");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP cv = PROTECT(Rf_ScalarLogical(report.converged != 0));
    const char *names[] = {"mesh", "num_iterations", "residual", "converged", "num_fixed",
                           "num_isolated", "max_displacement"};
    const SEXP values[] = {mo, NULL, NULL, cv, NULL, NULL, NULL};
    const double numbers[] = {(double)report.num_iterations, report.residual,
                              (double)report.num_fixed, (double)report.num_isolated,
                              report.max_displacement};
    SEXP res = mio_r_report_list(7, names, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_prepare_metric(SEXP mesh, SEXP arrays, SEXP output, SEXP target_vertices,
                          SEXP norm, SEXP hmin, SEXP hmax, SEXP gradation) {
    mio_metric_opts opts;
    mio_metric_opts_init(&opts);
    R_xlen_t n = 0;
    const char **names = mio_r_cstrings(arrays, &n);
    opts.arrays = names;
    opts.num_arrays = (int64_t)n;
    opts.output = mio_r_string(output, "output");
    opts.target_vertices = Rf_asReal(target_vertices);
    opts.norm = Rf_asReal(norm);
    opts.hmin = Rf_asReal(hmin);
    opts.hmax = Rf_asReal(hmax);
    opts.gradation = Rf_asReal(gradation);
    mio_metric_report rep;
    mio_mesh *out = mio_prepare_metric(mio_r_mesh(mesh), &opts, &rep);
    if (out == NULL) mio_r_fail("prepare_metric");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *keys[] = {"mesh", "complexity_before", "complexity_after", "num_graded"};
    const SEXP values[] = {mo, NULL, NULL, NULL};
    const double numbers[] = {rep.complexity_before, rep.complexity_after, (double)rep.num_graded};
    SEXP res = mio_r_report_list(4, keys, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_untangle(SEXP mesh, SEXP fix_boundary, SEXP rings, SEXP smooth, SEXP max_sweeps) {
    mio_untangle_opts opts;
    mio_untangle_opts_init(&opts);
    opts.fix_boundary = mio_r_bool(fix_boundary, "fix_boundary") ? 1 : 0;
    opts.rings = (int64_t)Rf_asReal(rings);
    opts.smooth = mio_r_bool(smooth, "smooth") ? 1 : 0;
    opts.max_sweeps = (int64_t)Rf_asReal(max_sweeps);
    mio_untangle_report rep;
    mio_mesh *out = mio_untangle(mio_r_mesh(mesh), &opts, &rep);
    if (out == NULL) mio_r_fail("untangle");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *keys[] = {"mesh", "num_inverted_before", "num_inverted_after",
                          "min_quality_before", "min_quality_after", "sweeps"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)rep.num_inverted_before, (double)rep.num_inverted_after,
                              rep.min_quality_before, rep.min_quality_after, (double)rep.sweeps};
    SEXP res = mio_r_report_list(6, keys, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_morph(SEXP mesh, SEXP array, SEXP fix_boundary, SEXP constrained_array,
                 SEXP stiffening, SEXP steps, SEXP allow_inverted) {
    mio_morph_opts opts;
    mio_morph_opts_init(&opts);
    opts.array = mio_r_string(array, "array");
    if (!Rf_isNull(constrained_array))
        opts.constrained_array = mio_r_string(constrained_array, "constrained_array");
    opts.fix_boundary = mio_r_bool(fix_boundary, "fix_boundary") ? 1 : 0;
    opts.stiffening = Rf_asReal(stiffening);
    opts.steps = (int64_t)Rf_asReal(steps);
    opts.allow_inverted = mio_r_bool(allow_inverted, "allow_inverted") ? 1 : 0;
    mio_morph_report rep;
    mio_mesh *out = mio_morph(mio_r_mesh(mesh), &opts, &rep);
    if (out == NULL) mio_r_fail("morph");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP e = PROTECT(Rf_ScalarLogical(rep.converged != 0));
    const char *keys[] = {"mesh", "num_constrained", "num_free", "num_inverted",
                          "max_displacement", "converged"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, e};
    const double numbers[] = {(double)rep.num_constrained, (double)rep.num_free,
                              (double)rep.num_inverted, rep.max_displacement};
    SEXP res = mio_r_report_list(6, keys, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_laplacian_eigenmodes(SEXP mesh, SEXP num_modes, SEXP tolerance,
                                SEXP max_iterations) {
    mio_eigenmode_opts opts;
    mio_eigenmode_opts_init(&opts);
    opts.num_modes = (int64_t)Rf_asReal(num_modes);
    opts.tolerance = Rf_asReal(tolerance);
    opts.max_iterations = (int64_t)Rf_asReal(max_iterations);
    mio_eigenmode_report rep;
    mio_mesh *out = mio_laplacian_eigenmodes(mio_r_mesh(mesh), &opts, &rep);
    if (out == NULL) mio_r_fail("laplacian_eigenmodes");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP c = PROTECT(Rf_ScalarLogical(rep.converged != 0));
    const char *keys[] = {"mesh", "iterations", "max_residual", "converged"};
    const SEXP values[] = {mo, NULL, NULL, c};
    const double numbers[] = {(double)rep.iterations, rep.max_residual};
    SEXP res = mio_r_report_list(4, keys, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_register_rigid(SEXP source, SEXP target, SEXP method, SEXP max_iterations,
                          SEXP tolerance, SEXP max_distance, SEXP initial, SEXP initial_mode,
                          SEXP symmetric_check, SEXP source_region, SEXP target_region,
                          SEXP adaptive_trim) {
    mio_registration_ex_opts opts;
    mio_registration_ex_opts_init(&opts);
    const char *m = mio_r_string(method, "method");
    if (strcmp(m, "point_to_plane") == 0) {
        opts.base.method = MIO_REGISTER_POINT_TO_PLANE;
    } else if (strcmp(m, "point_to_point") != 0) {
        Rf_error("`method` must be \"point_to_point\" or \"point_to_plane\"");
    }
    opts.base.max_iterations = (int64_t)Rf_asReal(max_iterations);
    opts.base.tolerance = Rf_asReal(tolerance);
    opts.base.max_distance = Rf_asReal(max_distance);
    if (!Rf_isNull(initial)) {
        SEXP x = PROTECT(Rf_coerceVector(initial, REALSXP));
        if (Rf_xlength(x) != 16) {
            UNPROTECT(1);
            Rf_error("`initial` must be a 4 x 4 matrix");
        }
        /* R stores column-major: element (i, j) at j * 4 + i */
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) opts.base.initial[4 * i + j] = REAL(x)[j * 4 + i];
        UNPROTECT(1);
    }
    const char *im = mio_r_string(initial_mode, "initial");
    if (strcmp(im, "auto") == 0) {
        opts.initial_mode = MIO_REGISTER_INITIAL_AUTO;
    } else if (strcmp(im, "matrix") != 0) {
        Rf_error("`initial` must be a 4 x 4 matrix, NULL or \"auto\"");
    }
    opts.symmetric_check = mio_r_bool(symmetric_check, "symmetric_check");
    opts.source_region = mio_r_opt_string(source_region);
    opts.target_region = mio_r_opt_string(target_region);
    opts.adaptive_trim = mio_r_double(adaptive_trim, "adaptive_trim");
    mio_registration_report rep;
    if (mio_register_rigid_ex(mio_r_mesh(source), mio_r_mesh(target), &opts, &rep) != MIO_OK)
        mio_r_fail("register_rigid");
    SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, 4, 4));
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) REAL(mat)[j * 4 + i] = rep.matrix[4 * i + j];
    SEXP a = PROTECT(Rf_ScalarReal(rep.initial_rms));
    SEXP e = PROTECT(Rf_ScalarLogical(rep.converged != 0));
    SEXP cands = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)rep.num_candidates));
    if (rep.num_candidates > 0) mio_register_rigid_candidates(REAL(cands), rep.num_candidates);
    SEXP warnings = PROTECT(mio_r_getbuffer(mio_register_rigid_warnings, "register_rigid"));
    const char *keys[] = {"matrix",       "initial_rms",      "rms",           "iterations",
                          "num_pairs",    "converged",        "num_rejected",  "chosen_candidate",
                          "candidate_rms", "warnings"};
    const SEXP values[] = {mat, a, NULL, NULL, NULL, e, NULL, NULL, cands, warnings};
    const double numbers[] = {rep.rms, (double)rep.iterations, (double)rep.num_pairs,
                              (double)rep.num_rejected, (double)rep.chosen_candidate};
    SEXP res = mio_r_report_list(10, keys, values, numbers);
    UNPROTECT(5);
    return res;
}

/* n x 3 R matrices (column-major) of picked pairs; `weights` NULL or n. */
SEXP R_mio_register_landmarks(SEXP source, SEXP target, SEXP weights) {
    SEXP s = PROTECT(Rf_coerceVector(source, REALSXP));
    SEXP t = PROTECT(Rf_coerceVector(target, REALSXP));
    const R_xlen_t n = Rf_xlength(s) / 3;
    if (Rf_xlength(s) != 3 * n || Rf_xlength(t) != 3 * n) {
        UNPROTECT(2);
        Rf_error("`source` and `target` must be n x 3 matrices of the same size");
    }
    double *src = (double *)R_alloc((size_t)(n > 0 ? 3 * n : 1), sizeof(double));
    double *dst = (double *)R_alloc((size_t)(n > 0 ? 3 * n : 1), sizeof(double));
    for (R_xlen_t i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) {
            src[3 * i + k] = REAL(s)[k * n + i];
            dst[3 * i + k] = REAL(t)[k * n + i];
        }
    const double *w = NULL;
    if (!Rf_isNull(weights)) {
        SEXP wv = PROTECT(Rf_coerceVector(weights, REALSXP));
        if (Rf_xlength(wv) != n) {
            UNPROTECT(3);
            Rf_error("`weights` must have one entry per pair");
        }
        double *wc = (double *)R_alloc((size_t)(n > 0 ? n : 1), sizeof(double));
        for (R_xlen_t i = 0; i < n; ++i) wc[i] = REAL(wv)[i];
        w = wc;
        UNPROTECT(1);
    }
    SEXP res_v = PROTECT(Rf_allocVector(REALSXP, n));
    double m[16], rms = 0.0;
    if (mio_register_landmarks(src, dst, w, (int64_t)n, m, n > 0 ? REAL(res_v) : NULL, &rms) !=
        MIO_OK)
        mio_r_fail("register_landmarks");
    SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, 4, 4));
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) REAL(mat)[j * 4 + i] = m[4 * i + j];
    const char *keys[] = {"matrix", "residuals", "rms"};
    const SEXP values[] = {mat, res_v, NULL};
    const double numbers[] = {rms};
    SEXP out = mio_r_report_list(3, keys, values, numbers);
    UNPROTECT(4);
    return out;
}

SEXP R_mio_find_symmetry_plane(SEXP mesh, SEXP seed_origin, SEXP seed_normal, SEXP tolerance,
                               SEXP min_share, SEXP max_samples) {
    mio_symmetry_plane_opts opts;
    mio_symmetry_plane_opts_init(&opts);
    if (!Rf_isNull(seed_origin) || !Rf_isNull(seed_normal)) {
        if (Rf_xlength(seed_origin) != 3 || Rf_xlength(seed_normal) != 3)
            Rf_error("`seed` must be a list of a 3-vector origin and a 3-vector normal");
        SEXP so = PROTECT(Rf_coerceVector(seed_origin, REALSXP));
        SEXP sn = PROTECT(Rf_coerceVector(seed_normal, REALSXP));
        for (int k = 0; k < 3; ++k) {
            opts.seed_origin[k] = REAL(so)[k];
            opts.seed_normal[k] = REAL(sn)[k];
        }
        UNPROTECT(2);
    }
    opts.tolerance = mio_r_double(tolerance, "tolerance");
    opts.min_share = mio_r_double(min_share, "min_share");
    opts.max_samples = (int64_t)mio_r_double(max_samples, "max_samples");
    mio_symmetry_plane_report r;
    if (mio_find_symmetry_plane(mio_r_mesh(mesh), &opts, &r) != MIO_OK)
        mio_r_fail("find_symmetry_plane");
    SEXP o = PROTECT(Rf_allocVector(REALSXP, 3));
    SEXP nn = PROTECT(Rf_allocVector(REALSXP, 3));
    for (int k = 0; k < 3; ++k) {
        REAL(o)[k] = r.origin[k];
        REAL(nn)[k] = r.normal[k];
    }
    const char *keys[] = {"origin", "normal", "rms", "share", "tolerance", "chosen_candidate"};
    const SEXP values[] = {o, nn, NULL, NULL, NULL, NULL};
    const double numbers[] = {r.rms, r.share, r.tolerance, (double)r.chosen_candidate};
    SEXP res = mio_r_report_list(6, keys, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_auto_align(SEXP mesh, SEXP angle_tolerance, SEXP subdivisions, SEXP region,
                      SEXP rotate_vector_data) {
    mio_auto_align_opts opts;
    mio_auto_align_opts_init(&opts);
    opts.angle_tolerance = mio_r_double(angle_tolerance, "angle_tolerance");
    opts.subdivisions = (int64_t)mio_r_double(subdivisions, "subdivisions");
    opts.region = mio_r_opt_string(region);
    opts.rotate_vector_data = mio_r_bool(rotate_vector_data, "rotate_vector_data");
    mio_auto_align_report r;
    mio_mesh *out = mio_auto_align(mio_r_mesh(mesh), &opts, &r);
    if (out == NULL) mio_r_fail("auto_align");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, 4, 4));
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) REAL(mat)[j * 4 + i] = r.matrix[4 * i + j];
    SEXP mode = PROTECT(Rf_mkString(r.mode == MIO_AUTO_ALIGN_AXIS ? "axis" : "face"));
    SEXP pr = PROTECT(Rf_allocVector(REALSXP, 3));
    SEXP se = PROTECT(Rf_allocVector(REALSXP, 3));
    for (int k = 0; k < 3; ++k) {
        REAL(pr)[k] = r.primary[k];
        REAL(se)[k] = r.secondary[k];
    }
    const char *keys[] = {"mesh",  "matrix",     "mode",       "primary",     "secondary",
                          "share", "face_share", "axis_share", "opening_area"};
    const SEXP values[] = {mo, mat, mode, pr, se, NULL, NULL, NULL, NULL};
    const double numbers[] = {r.share, r.face_share, r.axis_share, r.opening_area};
    SEXP res = mio_r_report_list(9, keys, values, numbers);
    UNPROTECT(5);
    return res;
}

SEXP R_mio_estimate_normals(SEXP mesh, SEXP k, SEXP orientation, SEXP viewpoint,
                            SEXP record_variation) {
    mio_point_normals_opts opts;
    mio_point_normals_opts_init(&opts);
    opts.k = (int64_t)Rf_asReal(k);
    const char *o = mio_r_string(orientation, "orientation");
    if (strcmp(o, "consistent") == 0) {
        opts.orientation = MIO_NORMALS_CONSISTENT;
    } else if (strcmp(o, "viewpoint") == 0) {
        opts.orientation = MIO_NORMALS_VIEWPOINT;
    } else if (strcmp(o, "none") == 0) {
        opts.orientation = MIO_NORMALS_NONE;
    } else {
        Rf_error("`orientation` must be \"consistent\", \"viewpoint\" or \"none\"");
    }
    SEXP v = PROTECT(Rf_coerceVector(viewpoint, REALSXP));
    if (Rf_xlength(v) != 3) {
        UNPROTECT(1);
        Rf_error("`viewpoint` takes three numbers");
    }
    for (int a = 0; a < 3; ++a) opts.viewpoint[a] = REAL(v)[a];
    UNPROTECT(1);
    opts.record_variation = mio_r_bool(record_variation, "record_variation") ? 1 : 0;
    int64_t nd = 0, nc = 0;
    mio_mesh *out = mio_estimate_normals(mio_r_mesh(mesh), &opts, &nd, &nc);
    if (out == NULL) mio_r_fail("estimate_normals");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP a = PROTECT(Rf_ScalarReal((double)nd));
    const char *keys[] = {"mesh", "num_degenerate", "num_components"};
    const SEXP values[] = {mo, a, NULL};
    const double numbers[] = {(double)nc};
    SEXP res = mio_r_report_list(3, keys, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_reconstruct_surface(SEXP mesh, SEXP normals, SEXP k, SEXP cell_size,
                               SEXP resolution, SEXP max_distance, SEXP estimate_normals,
                               SEXP normals_k, SEXP max_cells) {
    mio_reconstruct_surface_opts opts;
    mio_reconstruct_surface_opts_init(&opts);
    opts.normals = mio_r_string(normals, "normals");
    opts.k = (int64_t)Rf_asReal(k);
    opts.cell_size = Rf_asReal(cell_size);
    opts.resolution = (int64_t)Rf_asReal(resolution);
    opts.max_distance = Rf_asReal(max_distance);
    opts.estimate_normals = mio_r_bool(estimate_normals, "estimate_normals") ? 1 : 0;
    opts.normals_k = (int64_t)Rf_asReal(normals_k);
    opts.max_cells = (int64_t)Rf_asReal(max_cells);
    mio_reconstruct_surface_report rep;
    mio_mesh *out = mio_reconstruct_surface(mio_r_mesh(mesh), &opts, &rep);
    if (out == NULL) mio_r_fail("reconstruct_surface");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP i = PROTECT(Rf_ScalarLogical(rep.estimated_normals != 0));
    const char *keys[] = {"mesh",          "num_samples",        "spacing",
                          "cell_size",     "max_distance",       "num_lattice_points",
                          "num_undefined", "num_triangles",      "num_components",
                          "estimated_normals"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, i};
    const double numbers[] = {(double)rep.num_samples, rep.spacing, rep.cell_size, rep.max_distance,
                              (double)rep.num_lattice_points, (double)rep.num_undefined,
                              (double)rep.num_triangles, (double)rep.num_components};
    SEXP res = mio_r_report_list(10, keys, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_convex_hull(SEXP mesh, SEXP dimension) {
    const int d = Rf_asInteger(dimension);
    if (d != 0 && d != 2 && d != 3) Rf_error("`dimension` must be 0, 2 or 3");
    mio_hull_report rep;
    mio_mesh *out = mio_convex_hull(mio_r_mesh(mesh), (int32_t)d, &rep);
    if (out == NULL) mio_r_fail("convex_hull");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP a = PROTECT(Rf_ScalarInteger(rep.dimension));
    const char *keys[] = {"mesh", "dimension", "num_vertices", "num_facets",
                          "area", "volume", "perimeter"};
    const SEXP values[] = {mo, a, NULL, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)rep.num_vertices, (double)rep.num_facets, rep.area,
                              rep.volume, rep.perimeter};
    SEXP res = mio_r_report_list(7, keys, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_delaunay_2d(SEXP mesh) {
    int64_t nt = 0, nd = 0, nb = 0;
    mio_mesh *out = mio_delaunay_2d(mio_r_mesh(mesh), &nt, &nd, &nb);
    if (out == NULL) mio_r_fail("delaunay_2d");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *keys[] = {"mesh", "num_triangles", "num_duplicates", "num_boundary"};
    const SEXP values[] = {mo, NULL, NULL, NULL};
    const double numbers[] = {(double)nt, (double)nd, (double)nb};
    SEXP res = mio_r_report_list(4, keys, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_convert_coordinates(SEXP mesh, SEXP system, SEXP axis, SEXP origin, SEXP reference,
                               SEXP arrays, SEXP inverse, SEXP record_coordinates) {
    mio_coordinate_opts opts;
    mio_coordinate_opts_init(&opts);
    const char *sys = mio_r_string(system, "system");
    if (strcmp(sys, "spherical") == 0) {
        opts.system = MIO_COORDS_SPHERICAL;
    } else if (strcmp(sys, "cylindrical") != 0) {
        Rf_error("`system` must be \"cylindrical\" or \"spherical\"");
    }
    SEXP vecs[3] = {axis, origin, reference};
    double *dst[3] = {opts.axis, opts.origin, opts.reference};
    for (int v = 0; v < 3; ++v) {
        SEXP x = PROTECT(Rf_coerceVector(vecs[v], REALSXP));
        if (Rf_xlength(x) != 3) {
            UNPROTECT(1);
            Rf_error("`axis`, `origin` and `reference` take three numbers");
        }
        for (int k = 0; k < 3; ++k) dst[v][k] = REAL(x)[k];
        UNPROTECT(1);
    }
    R_xlen_t n = 0;
    const char **names = mio_r_cstrings(arrays, &n);
    opts.arrays = names;
    opts.num_arrays = (int64_t)n;
    opts.inverse = mio_r_bool(inverse, "inverse") ? 1 : 0;
    opts.record_coordinates = mio_r_bool(record_coordinates, "record_coordinates") ? 1 : 0;
    int64_t na = 0, ns = 0;
    mio_mesh *out = mio_convert_coordinates(mio_r_mesh(mesh), &opts, &na, &ns);
    if (out == NULL) mio_r_fail("convert_coordinates");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *keys[] = {"mesh", "num_arrays", "num_singular"};
    const SEXP values[] = {mo, NULL, NULL};
    const double numbers[] = {(double)na, (double)ns};
    SEXP res = mio_r_report_list(3, keys, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_geodesic_distance(SEXP mesh, SEXP sources, SEXP source_region, SEXP time_factor,
                             SEXP tolerance, SEXP max_iterations) {
    mio_geodesic_opts opts;
    mio_geodesic_report report;
    mio_geodesic_opts_init(&opts);
    SEXP src = PROTECT(Rf_coerceVector(sources, REALSXP));
    const R_xlen_t n = Rf_xlength(src);
    int64_t *ids = n > 0 ? (int64_t *)R_alloc((size_t)n, sizeof(int64_t)) : NULL;
    for (R_xlen_t i = 0; i < n; ++i) ids[i] = (int64_t)REAL(src)[i] - 1; /* 1-based in R */
    opts.sources = ids;
    opts.num_sources = (int64_t)n;
    opts.source_region = mio_r_opt_string(source_region);
    opts.time_factor = mio_r_double(time_factor, "time_factor");
    opts.tolerance = mio_r_double(tolerance, "tolerance");
    opts.max_iterations = (int64_t)mio_r_double(max_iterations, "max_iterations");
    mio_mesh *out = mio_geodesic_distance(mio_r_mesh(mesh), &opts, &report);
    UNPROTECT(1);
    if (out == NULL) mio_r_fail("geodesic_distance");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP cv = PROTECT(Rf_ScalarLogical(report.converged != 0));
    const char *names[] = {"mesh", "num_sources", "num_unreached", "max_distance", "time_step",
                           "converged"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, cv};
    const double numbers[] = {(double)report.num_sources, (double)report.num_unreached,
                              report.max_distance, report.time_step};
    SEXP res = mio_r_report_list(6, names, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_parameterize(SEXP mesh, SEXP padding, SEXP max_depth, SEXP max_stretch) {
    mio_parameterize_opts opts;
    mio_parameterize_report report;
    mio_parameterize_opts_init(&opts);
    opts.padding = mio_r_double(padding, "padding");
    opts.max_depth = (int64_t)mio_r_double(max_depth, "max_depth");
    opts.max_stretch = mio_r_double(max_stretch, "max_stretch");
    mio_mesh* out = mio_parameterize(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL)
        mio_r_fail("parameterize");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char* names[] = {"mesh",    "num_charts",         "num_flipped",
                           "stretch", "packing_efficiency", "num_seam_points"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)report.num_charts, (double)report.num_flipped, report.stretch,
                              report.packing_efficiency, (double)report.num_seam_points};
    SEXP res = mio_r_report_list(6, names, values, numbers);
    UNPROTECT(1);
    return res;
}

/* --- combining / comparing ---------------------------------------------- */

SEXP R_mio_merge(SEXP meshes, SEXP weld, SEXP atol, SEXP source_tag, SEXP data_policy,
                 SEXP drop_duplicate_cells) {
    if (TYPEOF(meshes) != VECSXP || Rf_length(meshes) < 1) {
        Rf_error("`meshes` must be a non-empty list of mio_mesh objects");
    }
    R_xlen_t n = Rf_length(meshes);
    SEXP buf = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)(n * sizeof(mio_mesh *))));
    const mio_mesh **handles = (const mio_mesh **)RAW(buf);
    for (R_xlen_t i = 0; i < n; ++i) {
        handles[i] = mio_r_mesh(VECTOR_ELT(meshes, i));
    }
    mio_mesh *out = mio_merge(handles, (int64_t)n, mio_r_bool(weld, "weld"),
                              mio_r_double(atol, "atol"),
                              mio_r_bool(source_tag, "source_tag"),
                              mio_r_int(data_policy, "data_policy"),
                              mio_r_bool(drop_duplicate_cells, "drop_duplicate_cells"));
    UNPROTECT(1);
    if (out == NULL) mio_r_fail("merge");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_interpolate(SEXP source, SEXP target, SEXP method, SEXP arrays, SEXP extrapolate,
                       SEXP default_value, SEXP on_conflict) {
    SEXP shelter;
    int64_t count = 0;
    const char *const *names = mio_r_names(arrays, &count, &shelter);
    mio_mesh *out = mio_interpolate(mio_r_mesh(source), mio_r_mesh(target),
                                    mio_r_opt_string(method), names, count,
                                    mio_r_bool(extrapolate, "extrapolate"),
                                    mio_r_double(default_value, "default_value"),
                                    mio_r_opt_string(on_conflict));
    UNPROTECT(1); /* shelter */
    if (out == NULL) mio_r_fail("interpolate");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_conservative_interpolate(SEXP source, SEXP target, SEXP arrays, SEXP default_value,
                                    SEXP on_conflict) {
    SEXP shelter;
    int64_t count = 0;
    const char *const *names = mio_r_names(arrays, &count, &shelter);
    mio_mesh *out = mio_conservative_interpolate(mio_r_mesh(source), mio_r_mesh(target), names,
                                                 count, mio_r_double(default_value, "default_value"),
                                                 mio_r_opt_string(on_conflict));
    UNPROTECT(1); /* shelter */
    if (out == NULL) mio_r_fail("conservative_interpolate");
    return mio_r_wrap_mesh(out);
}

/* error_size (roadmap 6.10.2): a per-point size field from a per-cell error
 * indicator. Returns list(mesh, num_cells, num_skipped, num_clamped, min_size,
 * max_size). */
SEXP R_mio_error_size(SEXP mesh, SEXP target_error, SEXP error_array, SEXP order,
                      SEXP min_size, SEXP max_size, SEXP output, SEXP overwrite) {
    mio_error_size_opts opts;
    mio_error_size_opts_init(&opts);
    opts.target_error = mio_r_double(target_error, "target_error");
    opts.error_array = mio_r_string(error_array, "error_array");
    opts.order = mio_r_double(order, "order");
    opts.min_size = mio_r_double(min_size, "min_size");
    opts.max_size = mio_r_double(max_size, "max_size");
    opts.output = mio_r_string(output, "output");
    opts.overwrite = mio_r_bool(overwrite, "overwrite");
    mio_error_size_report report;
    mio_mesh *out = mio_error_size(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("error_size");
    SEXP v0 = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh",        "num_cells", "num_skipped",
                           "num_clamped", "min_size",  "max_size"};
    const SEXP values[] = {v0, NULL, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)report.num_cells, (double)report.num_skipped,
                              (double)report.num_clamped, report.min_size, report.max_size};
    SEXP res = mio_r_report_list(6, names, values, numbers);
    UNPROTECT(1);
    return res;
}

/* check_conformity (roadmap 3.2.1): hanging nodes, T-junctions, partly
 * overlapping and unwelded facets. Returns list(mesh, dimension, ...,
 * interface_mode, warnings). */
SEXP R_mio_check_conformity(SEXP mesh, SEXP tolerance, SEXP attach_flags) {
    mio_conformity_opts opts;
    mio_conformity_opts_init(&opts);
    opts.tolerance = mio_r_double(tolerance, "tolerance");
    opts.attach_flags = mio_r_bool(attach_flags, "attach_flags");
    mio_conformity_report report;
    mio_mesh *out = mio_check_conformity(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("check_conformity");
    static const char *modes[] = {"conforming", "weld", "proximity"};
    const char *mode = (report.interface_mode >= 0 && report.interface_mode <= 2)
                           ? modes[report.interface_mode]
                           : "unknown";
    SEXP v0 = PROTECT(mio_r_wrap_mesh(out));
    SEXP v8 = PROTECT(Rf_ScalarLogical(report.conforming != 0));
    SEXP v9 = PROTECT(Rf_mkString(mode));
    SEXP v10 = PROTECT(mio_r_getbuffer(mio_check_conformity_warnings, "check_conformity"));
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
    const SEXP values[] = {v0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, v8, v9, v10};
    const double numbers[] = {(double)report.dimension, (double)report.num_hanging_nodes,
                              (double)report.num_t_junctions, (double)report.num_partial_overlaps,
                              (double)report.num_unwelded_facets,
                              (double)report.num_boundary_facets,
                              (double)report.num_non_manifold_facets};
    SEXP res = mio_r_report_list(11, names, values, numbers);
    UNPROTECT(4);
    return res;
}

/* check_tetrahedra: inverted, degenerate and overlapping tetrahedra. The mesh is
 * not changed. Returns list(num_tetrahedra, ..., valid, warnings); `zero_volume_factor`
 * NULL keeps the C default. */
SEXP R_mio_check_tetrahedra(SEXP mesh, SEXP zero_volume_factor) {
    mio_check_tetrahedra_opts opts;
    mio_check_tetrahedra_opts_init(&opts);
    if (!Rf_isNull(zero_volume_factor))
        opts.zero_volume_factor = mio_r_double(zero_volume_factor, "zero_volume_factor");
    mio_check_tetrahedra_report report;
    memset(&report, 0, sizeof(report));
    mio_r_check(mio_check_tetrahedra(mio_r_mesh(mesh), &opts, &report), "check_tetrahedra");
    SEXP v6 = PROTECT(Rf_ScalarLogical(report.valid != 0));
    SEXP v7 = PROTECT(mio_r_getbuffer(mio_check_tetrahedra_warnings, "check_tetrahedra"));
    const char *names[] = {"num_tetrahedra", "num_inverted", "num_degenerate", "num_overlaps",
                           "min_volume",     "max_volume",   "valid",          "warnings"};
    const SEXP values[] = {NULL, NULL, NULL, NULL, NULL, NULL, v6, v7};
    const double numbers[] = {(double)report.num_tetrahedra, (double)report.num_inverted,
                              (double)report.num_degenerate, (double)report.num_overlaps,
                              report.min_volume, report.max_volume};
    SEXP res = mio_r_report_list(8, names, values, numbers);
    UNPROTECT(2);
    return res;
}

/* boundary_layer (roadmap 6.9): prisms grown from a surface. Returns
 * list(mesh, num_cells, ..., cap_closed). */
SEXP R_mio_boundary_layer(SEXP mesh, SEXP first_height, SEXP num_layers, SEXP growth,
                          SEXP reverse, SEXP collision_fraction, SEXP smooth_iterations,
                          SEXP max_shrink, SEXP collapse, SEXP normal_smooth_iterations,
                          SEXP front_check) {
    mio_boundary_layer_opts opts;
    mio_boundary_layer_opts_init(&opts);
    opts.first_height = mio_r_double(first_height, "first_height");
    opts.num_layers = (int32_t)mio_r_int(num_layers, "num_layers");
    opts.growth = mio_r_double(growth, "growth");
    opts.reverse = mio_r_bool(reverse, "reverse");
    opts.collision_fraction = mio_r_double(collision_fraction, "collision_fraction");
    opts.smooth_iterations = (int32_t)mio_r_int(smooth_iterations, "smooth_iterations");
    opts.max_shrink = (int32_t)mio_r_int(max_shrink, "max_shrink");
    {
        const char *c = mio_r_string(collapse, "collapse");
        if (strcmp(c, "points") == 0)
            opts.collapse = MIO_BOUNDARY_LAYER_COLLAPSE_POINTS;
        else if (strcmp(c, "layers") == 0)
            opts.collapse = MIO_BOUNDARY_LAYER_COLLAPSE_LAYERS;
        else
            Rf_error("meshio++: boundary_layer: collapse must be \"points\" or \"layers\"");
    }
    opts.normal_smooth_iterations =
        (int32_t)mio_r_int(normal_smooth_iterations, "normal_smooth_iterations");
    opts.front_check = mio_r_bool(front_check, "front_check");
    mio_boundary_layer_report report;
    mio_mesh *out = mio_boundary_layer(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("boundary_layer");
    SEXP v0 = PROTECT(mio_r_wrap_mesh(out));
    SEXP v12 = PROTECT(Rf_ScalarLogical(report.cap_closed != 0));
    const char *names[] = {"mesh",          "num_cells",       "num_wedges",
                           "num_hexahedra", "num_pyramids",    "num_tetrahedra",
                           "num_shrunk",    "num_collapsed",   "num_not_visible",
                           "min_scaled_jacobian", "min_thickness", "max_thickness",
                           "cap_closed",    "num_fewer_layers", "num_front_crossings"};
    const SEXP values[] = {v0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                           v12, NULL, NULL};
    const double numbers[] = {(double)report.num_cells, (double)report.num_wedges,
                              (double)report.num_hexahedra, (double)report.num_pyramids,
                              (double)report.num_tetrahedra, (double)report.num_shrunk,
                              (double)report.num_collapsed, (double)report.num_not_visible,
                              report.min_scaled_jacobian, report.min_thickness,
                              report.max_thickness, (double)report.num_fewer_layers,
                              (double)report.num_front_crossings};
    SEXP res = mio_r_report_list(15, names, values, numbers);
    UNPROTECT(2);
    return res;
}

/* The list an extrude or revolve call returns. */
static SEXP r_sweep_list(mio_mesh *out, const mio_sweep_report *report, SEXP warnings) {
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP full = PROTECT(Rf_ScalarLogical(report->full != 0));
    const char *names[] = {"mesh",          "num_cells",      "num_hexahedra", "num_wedges",
                           "num_pyramids",  "num_tetrahedra", "num_quads",     "num_triangles",
                           "num_flipped",   "num_on_axis",    "min_scaled_jacobian",
                           "full",          "warnings"};
    const SEXP values[] = {mo,   NULL, NULL, NULL, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, full, warnings};
    const double numbers[] = {
        (double)report->num_cells,     (double)report->num_hexahedra,
        (double)report->num_wedges,    (double)report->num_pyramids,
        (double)report->num_tetrahedra, (double)report->num_quads,
        (double)report->num_triangles, (double)report->num_flipped,
        (double)report->num_on_axis,   report->min_scaled_jacobian};
    SEXP res = mio_r_report_list(13, names, values, numbers);
    UNPROTECT(2);
    return res;
}

/* extrude (roadmap 3.1.2): a surface swept into a volume, or a curve into a
 * surface, along a vector. Returns list(mesh, num_cells, ..., warnings). */
SEXP R_mio_extrude(SEXP mesh, SEXP direction, SEXP num_layers, SEXP ratio, SEXP offsets,
                   SEXP record_layers, SEXP max_cells) {
    mio_extrude_opts opts;
    mio_extrude_opts_init(&opts);
    mio_r_vec3(direction, opts.direction, "direction");
    opts.num_layers = (int32_t)mio_r_int(num_layers, "num_layers");
    opts.ratio = mio_r_double(ratio, "ratio");
    SEXP offs = PROTECT(Rf_coerceVector(offsets, REALSXP));
    opts.num_offsets = (int64_t)Rf_xlength(offs);
    opts.offsets = opts.num_offsets > 0 ? REAL(offs) : NULL;
    opts.record_layers = mio_r_bool(record_layers, "record_layers");
    opts.max_cells = (int64_t)mio_r_double(max_cells, "max_cells");
    mio_sweep_report report;
    mio_mesh *out = mio_extrude(mio_r_mesh(mesh), &opts, &report);
    UNPROTECT(1);
    if (out == NULL) mio_r_fail("extrude");
    SEXP w = PROTECT(mio_r_getbuffer(mio_extrude_warnings, "extrude"));
    SEXP res = PROTECT(r_sweep_list(out, &report, w));
    UNPROTECT(2);
    return res;
}

/* revolve (roadmap 3.1.3): the same around an axis. */
SEXP R_mio_revolve(SEXP mesh, SEXP axis, SEXP axis_point, SEXP angle, SEXP num_layers,
                   SEXP tolerance, SEXP record_layers, SEXP max_cells) {
    mio_revolve_opts opts;
    mio_revolve_opts_init(&opts);
    mio_r_vec3(axis, opts.axis, "axis");
    mio_r_vec3(axis_point, opts.axis_point, "axis_point");
    opts.angle = mio_r_double(angle, "angle");
    opts.num_layers = (int32_t)mio_r_int(num_layers, "num_layers");
    opts.tolerance = mio_r_double(tolerance, "tolerance");
    opts.record_layers = mio_r_bool(record_layers, "record_layers");
    opts.max_cells = (int64_t)mio_r_double(max_cells, "max_cells");
    mio_sweep_report report;
    mio_mesh *out = mio_revolve(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("revolve");
    SEXP w = PROTECT(mio_r_getbuffer(mio_revolve_warnings, "revolve"));
    SEXP res = PROTECT(r_sweep_list(out, &report, w));
    UNPROTECT(2);
    return res;
}

/* mirror (roadmap 3.2.5): reflect across a plane, every cell right side out.
 * Returns list(mesh, num_flipped_cells, ..., warnings). */
SEXP R_mio_mirror(SEXP mesh, SEXP origin, SEXP normal, SEXP keep_original, SEXP weld,
                  SEXP tolerance, SEXP reflect_vector_data) {
    mio_mirror_opts opts;
    mio_mirror_opts_init(&opts);
    mio_r_vec3(origin, opts.origin, "origin");
    mio_r_vec3(normal, opts.normal, "normal");
    opts.keep_original = mio_r_bool(keep_original, "keep_original");
    opts.weld = mio_r_bool(weld, "weld");
    opts.tolerance = mio_r_double(tolerance, "tolerance");
    opts.reflect_vector_data = mio_r_bool(reflect_vector_data, "reflect_vector_data");
    mio_mirror_report report;
    mio_mesh *out = mio_mirror(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("mirror");
    SEXP v0 = PROTECT(mio_r_wrap_mesh(out));
    SEXP v5 = PROTECT(mio_r_getbuffer(mio_mirror_warnings, "mirror"));
    const char *names[] = {"mesh",              "num_flipped_cells", "num_points_on_plane",
                           "num_dropped_cells", "num_dropped_sides", "warnings"};
    const SEXP values[] = {v0, NULL, NULL, NULL, NULL, v5};
    const double numbers[] = {(double)report.num_flipped_cells, (double)report.num_points_on_plane,
                              (double)report.num_dropped_cells, (double)report.num_dropped_sides};
    SEXP res = mio_r_report_list(6, names, values, numbers);
    UNPROTECT(2);
    return res;
}


SEXP R_mio_find_intersections(SEXP mesh, SEXP other, SEXP region, SEXP region_other,
                              SEXP segments, SEXP max_pairs) {
    mio_find_intersections_opts opts;
    mio_find_intersections_opts_init(&opts);
    opts.region = mio_r_opt_string(region);
    opts.region_b = mio_r_opt_string(region_other);
    opts.segments = mio_r_bool(segments, "segments");
    opts.max_pairs = (int64_t)mio_r_double(max_pairs, "max_pairs");
    const mio_mesh *b = other == R_NilValue ? NULL : mio_r_mesh(other);
    mio_find_intersections_report report;
    mio_find_intersections_result *result =
        mio_find_intersections(mio_r_mesh(mesh), b, &opts, &report);
    if (result == NULL) mio_r_fail("find_intersections");
    int64_t count = 0;
    const int64_t *pairs = mio_find_intersections_result_pairs(result, &count);
    if (count > INT_MAX) {
        mio_find_intersections_result_free(result);
        Rf_error("find_intersections: too many pairs for an R matrix; set max_pairs");
    }
    mio_mesh *flagged = mio_find_intersections_result_take_mesh(result);
    mio_mesh *flagged_b = mio_find_intersections_result_take_mesh_b(result);
    mio_mesh *segs = mio_find_intersections_result_take_segments(result);
    SEXP v[13];
    v[0] = PROTECT(mio_r_wrap_mesh(flagged));
    v[1] = PROTECT(flagged_b ? mio_r_wrap_mesh(flagged_b) : R_NilValue);
    v[2] = PROTECT(segs ? mio_r_wrap_mesh(segs) : R_NilValue);
    v[3] = PROTECT(Rf_allocMatrix(REALSXP, (int)count, 2));
    for (int64_t i = 0; i < count; ++i) {
        REAL(v[3])[i] = (double)(pairs[2 * i] + 1);
        REAL(v[3])[count + i] = (double)(pairs[2 * i + 1] + 1);
    }
    mio_find_intersections_result_free(result);
    v[10] = PROTECT(Rf_ScalarLogical(report.truncated != 0));
    v[11] = PROTECT(Rf_ScalarLogical(report.clean != 0));
    v[12] = PROTECT(mio_r_getbuffer(mio_find_intersections_warnings, "find_intersections"));
    const char *names[] = {"mesh",          "other",       "segments",     "pairs",
                           "num_pairs",     "num_triangle_pairs",          "num_coplanar",
                           "num_segments",  "num_cells_a", "num_cells_b",  "truncated",
                           "clean",         "warnings"};
    const SEXP values[] = {v[0], v[1], v[2], v[3], NULL, NULL, NULL,
                           NULL, NULL, NULL, v[10], v[11], v[12]};
    const double numbers[] = {(double)report.num_pairs,    (double)report.num_triangle_pairs,
                              (double)report.num_coplanar, (double)report.num_segments,
                              (double)report.num_cells_a,  (double)report.num_cells_b};
    SEXP res = mio_r_report_list(13, names, values, numbers);
    UNPROTECT(7);
    return res;
}

SEXP R_mio_thickness(SEXP mesh, SEXP other, SEXP direction, SEXP max_distance, SEXP threshold,
                     SEXP region, SEXP region_other, SEXP method, SEXP cone_angle,
                     SEXP cone_rays, SEXP two_sided, SEXP max_normal_deviation) {
    mio_thickness_opts opts;
    mio_thickness_opts_init(&opts);
    const char *dir = mio_r_string(direction, "direction");
    if (strcmp(dir, "auto") == 0) opts.direction = MIO_THICKNESS_AUTO;
    else if (strcmp(dir, "inward") == 0) opts.direction = MIO_THICKNESS_INWARD;
    else if (strcmp(dir, "outward") == 0) opts.direction = MIO_THICKNESS_OUTWARD;
    else Rf_error("`direction` must be 'auto', 'inward' or 'outward'");
    opts.max_distance = mio_r_double(max_distance, "max_distance");
    opts.threshold = mio_r_double(threshold, "threshold");
    opts.region = mio_r_opt_string(region);
    opts.region_target = mio_r_opt_string(region_other);
    const char *meth = mio_r_string(method, "method");
    if (strcmp(meth, "ray") == 0) opts.method = MIO_THICKNESS_RAY;
    else if (strcmp(meth, "sphere") == 0) opts.method = MIO_THICKNESS_SPHERE;
    else if (strcmp(meth, "cone") == 0) opts.method = MIO_THICKNESS_CONE;
    else Rf_error("`method` must be 'ray', 'sphere' or 'cone'");
    opts.cone_angle = mio_r_double(cone_angle, "cone_angle");
    opts.cone_rays = (int64_t)mio_r_double(cone_rays, "cone_rays");
    opts.two_sided = mio_r_bool(two_sided, "two_sided");
    opts.max_normal_deviation = mio_r_double(max_normal_deviation, "max_normal_deviation");
    const mio_mesh *target = other == R_NilValue ? NULL : mio_r_mesh(other);
    mio_thickness_report r;
    mio_mesh *out = mio_thickness(mio_r_mesh(mesh), target, &opts, &r);
    if (out == NULL) mio_r_fail("thickness");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP dir_name =
        PROTECT(Rf_mkString(r.direction == MIO_THICKNESS_OUTWARD ? "outward" : "inward"));
    SEXP min_point = PROTECT(Rf_allocVector(REALSXP, 3));
    for (int k = 0; k < 3; ++k)
        REAL(min_point)[k] = r.min_point[k];
    SEXP warnings = PROTECT(mio_r_getbuffer(mio_thickness_warnings, "thickness"));
    static const char *const kMethods[] = {"ray", "sphere", "cone"};
    SEXP meth_name = PROTECT(Rf_mkString(kMethods[r.method >= 0 && r.method <= 2 ? r.method : 0]));
    const char *names[] = {"mesh",     "direction", "num_cells",    "num_hit",  "num_missed",
                           "num_below", "min",      "max",          "mean",     "min_cell",
                           "min_point", "max_distance", "warnings", "method", "p2", "p95"};
    const SEXP values[] = {mo,   dir_name, NULL,      NULL, NULL, NULL,    NULL,
                           NULL, NULL,      NULL,      min_point,  NULL,    warnings,
                           meth_name, NULL, NULL};
    const double numbers[] = {(double)r.num_cells,
                              (double)r.num_hit,
                              (double)r.num_missed,
                              (double)r.num_below,
                              r.min,
                              r.max,
                              r.mean,
                              r.min_cell < 0 ? 0.0 : (double)(r.min_cell + 1),
                              r.max_distance,
                              r.p2,
                              r.p95};
    SEXP res = mio_r_report_list(16, names, values, numbers);
    UNPROTECT(5);
    return res;
}

SEXP R_mio_has_netgen(void) {
    return Rf_ScalarLogical(mio_has_netgen() != 0);
}

/* The backend names tetrahedralize takes besides "auto" (roadmap 7.7.2). */
SEXP R_mio_tetrahedralize_backends(void) {
    const int64_t n = mio_tet_mesher_count();
    if (n < 0) mio_r_fail("tetrahedralize_backends");
    SEXP out = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)n));
    for (int64_t i = 0; i < n; ++i) {
        const char *name = mio_tet_mesher_name(i);
        SET_STRING_ELT(out, (R_xlen_t)i, Rf_mkChar(name != NULL ? name : ""));
    }
    UNPROTECT(1);
    return out;
}

/* remesh_quads (roadmap 7.7.3): field-aligned quads through an installed
 * QuadriFlow or Instant Meshes. */
SEXP R_mio_remesh_quads(SEXP mesh, SEXP target_faces, SEXP backend, SEXP sharp, SEXP boundary,
                        SEXP adaptive, SEXP dominant, SEXP seed, SEXP max_cells) {
    mio_remesh_quads_opts opts;
    mio_remesh_quads_opts_init(&opts);
    opts.target_faces = (int64_t)mio_r_double(target_faces, "target_faces");
    opts.backend = mio_r_string(backend, "backend");
    opts.sharp = mio_r_bool(sharp, "sharp");
    opts.boundary = mio_r_bool(boundary, "boundary");
    opts.adaptive = mio_r_bool(adaptive, "adaptive");
    opts.dominant = mio_r_bool(dominant, "dominant");
    opts.seed = (int32_t)mio_r_int(seed, "seed");
    opts.max_cells = (int64_t)mio_r_double(max_cells, "max_cells");
    mio_remesh_quads_report report;
    mio_mesh *out = mio_remesh_quads(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("remesh_quads");
    SEXP v0 = PROTECT(mio_r_wrap_mesh(out));
    SEXP v7 = PROTECT(Rf_mkString(report.backend_name != NULL ? report.backend_name : ""));
    const char *names[] = {"mesh",
                           "num_quads",
                           "num_triangles",
                           "num_points",
                           "min_scaled_jacobian",
                           "mean_scaled_jacobian",
                           "hausdorff_to_surface",
                           "backend"};
    const SEXP values[] = {v0, NULL, NULL, NULL, NULL, NULL, NULL, v7};
    const double numbers[] = {(double)report.num_quads, (double)report.num_triangles,
                              (double)report.num_points, report.min_scaled_jacobian,
                              report.mean_scaled_jacobian, report.hausdorff_to_surface};
    SEXP res = mio_r_report_list(8, names, values, numbers);
    UNPROTECT(2);
    return res;
}

/* The quad remeshers remesh_quads takes by name. */
SEXP R_mio_remesh_quads_backends(void) {
    const int64_t n = mio_quad_remesher_count();
    if (n < 0) mio_r_fail("remesh_quads_backends");
    SEXP out = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)n));
    for (int64_t i = 0; i < n; ++i) {
        const char *name = mio_quad_remesher_name(i);
        SET_STRING_ELT(out, (R_xlen_t)i, Rf_mkChar(name != NULL ? name : ""));
    }
    UNPROTECT(1);
    return out;
}

/* cut_level_set and cut_quadrature (doc/cut.md): the part level sets select,
 * through CutCells. `planes` is flat, six numbers per plane. */
static void mio_r_cut_opts(mio_cut_opts *opts, SEXP level_sets, SEXP planes, SEXP selection,
                           SEXP include_uncut, SEXP backend, SEXP triangulate,
                           SEXP record_parent_ids) {
    mio_cut_opts_init(opts);
    R_xlen_t n = 0;
    const char **names = mio_r_cstrings(level_sets, &n);
    opts->level_sets = names;
    opts->num_level_sets = (int64_t)n;
    const R_xlen_t np = Rf_xlength(planes);
    if (np % 6 != 0) Rf_error("meshio++: cut: planes must hold six numbers per plane (origin, then normal)");
    opts->planes = np > 0 ? REAL(planes) : NULL;
    opts->num_planes = (int64_t)(np / 6);
    opts->selection = mio_r_string(selection, "selection");
    opts->backend = mio_r_string(backend, "backend");
    opts->include_uncut = mio_r_bool(include_uncut, "include_uncut");
    opts->triangulate = mio_r_bool(triangulate, "triangulate");
    opts->record_parent_ids = mio_r_bool(record_parent_ids, "record_parent_ids");
}

SEXP R_mio_cut_level_set(SEXP mesh, SEXP level_sets, SEXP planes, SEXP selection,
                         SEXP include_uncut, SEXP backend, SEXP degree, SEXP triangulate,
                         SEXP record_parent_ids) {
    mio_cut_opts opts;
    mio_r_cut_opts(&opts, level_sets, planes, selection, include_uncut, backend, triangulate,
                   record_parent_ids);
    opts.degree = (int32_t)mio_r_int(degree, "degree");
    mio_cut_report report;
    mio_mesh *out = mio_cut_level_set(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("cut_level_set");
    SEXP v0 = PROTECT(mio_r_wrap_mesh(out));
    SEXP v5 = PROTECT(Rf_ScalarInteger((int)report.dimension));
    SEXP v7 = PROTECT(Rf_mkString(report.backend_name != NULL ? report.backend_name : ""));
    const char *names[] = {"mesh",       "num_cells", "num_cut_cells", "num_uncut_cells",
                           "num_pieces", "dimension", "measure",       "backend"};
    const SEXP values[] = {v0, NULL, NULL, NULL, NULL, v5, NULL, v7};
    const double numbers[] = {(double)report.num_cells, (double)report.num_cut_cells,
                              (double)report.num_uncut_cells, (double)report.num_pieces,
                              report.measure};
    SEXP res = mio_r_report_list(8, names, values, numbers);
    UNPROTECT(3);
    return res;
}

SEXP R_mio_cut_quadrature(SEXP mesh, SEXP level_sets, SEXP planes, SEXP selection, SEXP order,
                          SEXP include_uncut, SEXP backend, SEXP triangulate,
                          SEXP record_parent_ids) {
    mio_cut_opts opts;
    mio_r_cut_opts(&opts, level_sets, planes, selection, include_uncut, backend, triangulate,
                   record_parent_ids);
    opts.order = (int32_t)mio_r_int(order, "order");
    mio_cut_report report;
    mio_mesh *out = mio_cut_quadrature(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("cut_quadrature");
    SEXP v0 = PROTECT(mio_r_wrap_mesh(out));
    SEXP v4 = PROTECT(Rf_ScalarInteger((int)report.dimension));
    SEXP v5 = PROTECT(Rf_mkString(report.backend_name != NULL ? report.backend_name : ""));
    const char *names[] = {"mesh", "num_rules", "num_points", "sum_weights", "dimension", "backend"};
    const SEXP values[] = {v0, NULL, NULL, NULL, v4, v5};
    const double numbers[] = {(double)report.num_rules, (double)report.num_points, report.measure};
    SEXP res = mio_r_report_list(6, names, values, numbers);
    UNPROTECT(3);
    return res;
}

SEXP R_mio_has_cutcells(void) {
    return Rf_ScalarLogical(mio_has_cutcells() != 0);
}

SEXP R_mio_has_mmg(void) {
    return Rf_ScalarLogical(mio_has_mmg() != 0);
}

SEXP R_mio_has_mmg_elas(void) {
    return Rf_ScalarLogical(mio_has_mmg_elas() != 0);
}

/* The report of an MMG operation, with its mesh, as a named list. */
/* An MMG operation's result list; mesh_domain appends its two counts
 * (`n_extra` of them, named by `extra_names`) after the report's fields. */
static SEXP r_mmg_result_ex(mio_mesh *out, const mio_mmg_report *report,
                            const char *const *extra_names, const double *extra, int n_extra) {
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP status = PROTECT(Rf_mkString(report->status == 0 ? "success" : "low_failure"));
    SEXP library = PROTECT(Rf_mkString(report->library));
    SEXP quality = PROTECT(Rf_allocVector(REALSXP, 3));
    REAL(quality)[0] = report->quality_min;
    REAL(quality)[1] = report->quality_mean;
    REAL(quality)[2] = report->quality_max;
    const char *names[] = {"mesh",          "status",           "library",
                           "num_points_in", "num_points_out",   "num_elements_in",
                           "num_elements_out", "quality",       "num_warnings",
                           NULL,               NULL};
    const SEXP values[] = {mo, status, library, NULL, NULL, NULL, NULL, quality, NULL, NULL, NULL};
    double numbers[] = {(double)report->num_points_in,   (double)report->num_points_out,
                        (double)report->num_elements_in, (double)report->num_elements_out,
                        (double)report->num_warnings,    0.0,
                        0.0};
    for (int i = 0; i < n_extra; ++i) {
        names[9 + i] = extra_names[i];
        numbers[5 + i] = extra[i];
    }
    SEXP res = mio_r_report_list(9 + n_extra, names, values, numbers);
    UNPROTECT(4);
    return res;
}

static SEXP r_mmg_result(mio_mesh *out, const mio_mmg_report *report) {
    return r_mmg_result_ex(out, report, NULL, NULL, 0);
}

/* The shared MMG parameters: sizes, the ridge angle, 1-based frozen points and
   the required region names. */
static void r_mmg_params(mio_mmg_params *p, SEXP hmin, SEXP hmax, SEXP hausd, SEXP hgrad,
                         SEXP angle_detection, SEXP frozen, SEXP required) {
    mio_mmg_params_init(p);
    p->hmin = mio_r_double(hmin, "hmin");
    p->hmax = mio_r_double(hmax, "hmax");
    p->hausd = mio_r_double(hausd, "hausd");
    p->hgrad = mio_r_double(hgrad, "hgrad");
    p->angle_detection = mio_r_double(angle_detection, "angle_detection");
    R_xlen_t n_frozen = 0;
    p->frozen = r_ids_0based(frozen, &n_frozen);
    p->num_frozen = (int64_t)n_frozen;
    R_xlen_t n_required = Rf_length(required);
    if (n_required > 0) {
        const char **names = (const char **)R_alloc((size_t)n_required, sizeof(char *));
        for (R_xlen_t i = 0; i < n_required; ++i) names[i] = CHAR(STRING_ELT(required, i));
        p->required = names;
        p->num_required = (int64_t)n_required;
    }
}

SEXP R_mio_adapt(SEXP mesh, SEXP hsiz, SEXP size, SEXP metric, SEXP hmin, SEXP hmax, SEXP hausd,
                 SEXP hgrad, SEXP angle_detection, SEXP frozen, SEXP required) {
    mio_adapt_opts opts;
    mio_adapt_opts_init(&opts);
    r_mmg_params(&opts.params, hmin, hmax, hausd, hgrad, angle_detection, frozen, required);
    opts.hsiz = mio_r_double(hsiz, "hsiz");
    const char *s = mio_r_string(size, "size");
    const char *m = mio_r_string(metric, "metric");
    opts.size = s[0] ? s : NULL;
    opts.metric = m[0] ? m : NULL;
    mio_mmg_report report;
    mio_mesh *out = mio_adapt(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("adapt");
    return r_mmg_result(out, &report);
}

/* discretize_levelset (roadmap 3.2.2). `materials` is NULL or a numeric k x 4
   matrix of (ref, split, ref_inside, ref_outside) rows. */
SEXP R_mio_discretize_levelset(SEXP mesh, SEXP level_set, SEXP isovalue, SEXP rmc, SEXP materials,
                               SEXP hmin, SEXP hmax, SEXP hausd, SEXP hgrad, SEXP angle_detection,
                               SEXP frozen, SEXP required) {
    mio_discretize_levelset_opts opts;
    mio_discretize_levelset_opts_init(&opts);
    r_mmg_params(&opts.params, hmin, hmax, hausd, hgrad, angle_detection, frozen, required);
    opts.level_set = mio_r_string(level_set, "level_set");
    opts.isovalue = mio_r_double(isovalue, "isovalue");
    opts.remove_small_components = mio_r_double(rmc, "remove_small_components");
    if (!Rf_isNull(materials)) {
        if (!Rf_isReal(materials) || Rf_length(materials) % 4 != 0)
            Rf_error("materials must be a numeric matrix with 4 columns");
        R_xlen_t k = Rf_length(materials) / 4;
        const double* m = REAL(materials);
        mio_mmg_material* mats =
            (mio_mmg_material*)R_alloc((size_t)(k > 0 ? k : 1), sizeof(mio_mmg_material));
        for (R_xlen_t i = 0; i < k; ++i) {
            memset(&mats[i], 0, sizeof(mats[i]));
            mats[i].ref = (int64_t)m[i];
            mats[i].split = m[i + k] != 0.0;
            mats[i].ref_inside = (int64_t)m[i + 2 * k];
            mats[i].ref_outside = (int64_t)m[i + 3 * k];
        }
        opts.materials = mats;
        opts.num_materials = (int64_t)k;
    }
    mio_mmg_report report;
    mio_mesh* out = mio_discretize_levelset(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL)
        mio_r_fail("discretize_levelset");
    return r_mmg_result(out, &report);
}

/* mesh_domain (roadmap 3.2.3). `holes` is NULL or a numeric k x 2 matrix of
   points, one inside each hole. Returns the MMG report plus num_subdomains and
   num_holes. */
SEXP R_mio_mesh_domain(SEXP boundary, SEXP hsiz, SEXP size, SEXP holes, SEXP hmin, SEXP hmax,
                       SEXP hausd, SEXP hgrad, SEXP angle_detection, SEXP frozen, SEXP required) {
    mio_mesh_domain_opts opts;
    mio_mesh_domain_opts_init(&opts);
    r_mmg_params(&opts.params, hmin, hmax, hausd, hgrad, angle_detection, frozen, required);
    opts.hsiz = mio_r_double(hsiz, "hsiz");
    const char* s = mio_r_string(size, "size");
    opts.size = s[0] ? s : NULL;
    if (!Rf_isNull(holes)) {
        if (!Rf_isReal(holes) || Rf_length(holes) % 2 != 0)
            Rf_error("holes must be a numeric matrix with 2 columns");
        R_xlen_t k = Rf_length(holes) / 2;
        const double* h = REAL(holes);
        double* xy = (double*)R_alloc((size_t)(k > 0 ? 2 * k : 1), sizeof(double));
        for (R_xlen_t i = 0; i < k; ++i) {
            xy[2 * i] = h[i];
            xy[2 * i + 1] = h[i + k];
        }
        opts.holes = xy;
        opts.num_holes = (int64_t)k;
    }
    mio_mmg_report report;
    int64_t num_subdomains = 0, num_holes = 0;
    mio_mesh* out =
        mio_mesh_domain(mio_r_mesh(boundary), &opts, &report, &num_subdomains, &num_holes);
    if (out == NULL)
        mio_r_fail("mesh_domain");
    const char *extra_names[] = {"num_subdomains", "num_holes"};
    const double extra[] = {(double)num_subdomains, (double)num_holes};
    return r_mmg_result_ex(out, &report, extra_names, extra, 2);
}

/* move (roadmap 3.2.4): MMG's Lagrangian mode, with ELAS. */
SEXP R_mio_move(SEXP mesh, SEXP displacement, SEXP moving, SEXP lagrangian, SEXP hmin, SEXP hmax,
                SEXP hausd, SEXP hgrad, SEXP angle_detection, SEXP frozen, SEXP required) {
    mio_move_opts opts;
    mio_move_opts_init(&opts);
    r_mmg_params(&opts.params, hmin, hmax, hausd, hgrad, angle_detection, frozen, required);
    opts.displacement = mio_r_string(displacement, "displacement");
    R_xlen_t n_moving = Rf_length(moving);
    if (n_moving > 0) {
        const char** names = (const char**)R_alloc((size_t)n_moving, sizeof(char*));
        for (R_xlen_t i = 0; i < n_moving; ++i)
            names[i] = CHAR(STRING_ELT(moving, i));
        opts.moving = names;
        opts.num_moving = (int64_t)n_moving;
    }
    opts.lagrangian = (int32_t)mio_r_double(lagrangian, "lagrangian");
    mio_mmg_report report;
    mio_mesh* out = mio_move(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL)
        mio_r_fail("move");
    return r_mmg_result(out, &report);
}

SEXP R_mio_metric_from_hessian(SEXP mesh, SEXP hessian, SEXP error, SEXP hmin, SEXP hmax,
                               SEXP output) {
    mio_mesh *out = mio_metric_from_hessian(
        mio_r_mesh(mesh), mio_r_string(hessian, "hessian"), mio_r_double(error, "error"),
        mio_r_double(hmin, "hmin"), mio_r_double(hmax, "hmax"), mio_r_string(output, "output"));
    if (out == NULL) mio_r_fail("metric_from_hessian");
    return mio_r_wrap_mesh(out);
}

/* tetrahedralize (roadmap 6.7): tetrahedra filling a closed triangle surface.
 * Returns list(mesh, num_tetrahedra, ..., backend). */
SEXP R_mio_tetrahedralize(SEXP mesh, SEXP backend, SEXP max_size, SEXP max_volume, SEXP grading,
                          SEXP optimize_steps, SEXP max_cells, SEXP face_array) {
    mio_tetrahedralize_opts opts;
    mio_tetrahedralize_opts_init(&opts);
    {
        const char *b = mio_r_string(backend, "backend");
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
            Rf_error("meshio++: tetrahedralize: backend must be \"auto\", \"netgen\", \"octree\" "
                     "or a mesher's name");
    }
    opts.max_size = mio_r_double(max_size, "max_size");
    opts.max_volume = mio_r_double(max_volume, "max_volume");
    opts.grading = mio_r_double(grading, "grading");
    opts.optimize_steps = (int32_t)mio_r_int(optimize_steps, "optimize_steps");
    opts.max_cells = (int64_t)mio_r_double(max_cells, "max_cells");
    opts.face_array = mio_r_string(face_array, "face_array");
    mio_tetrahedralize_report report;
    mio_mesh *out = mio_tetrahedralize(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("tetrahedralize");
    SEXP v0 = PROTECT(mio_r_wrap_mesh(out));
    SEXP v6 = PROTECT(Rf_ScalarLogical(report.boundary_preserved != 0));
    SEXP v11 = PROTECT(Rf_mkString(report.backend_name != NULL ? report.backend_name
                                : report.backend == MIO_TETRAHEDRALIZE_BACKEND_NETGEN
                                    ? "netgen"
                                    : "octree"));
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
    const SEXP values[] = {v0, NULL, NULL, NULL, NULL, NULL, v6, NULL, NULL, NULL, NULL, v11};
    const double numbers[] = {(double)report.num_tetrahedra, (double)report.num_points, (double)report.num_new_points, (double)report.num_regions, (double)report.num_boundary_triangles, report.hausdorff_to_surface, report.min_dihedral_deg, report.max_dihedral_deg, report.volume};
    SEXP res = mio_r_report_list(12, names, values, numbers);
    UNPROTECT(3);
    return res;
}

SEXP R_mio_undo_green(SEXP coarse, SEXP fine) {
    int64_t ngroups = 0, nremoved = 0;
    mio_mesh *out =
        mio_undo_green(mio_r_mesh(coarse), mio_r_mesh(fine), &ngroups, &nremoved);
    if (out == NULL) mio_r_fail("undo_green");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh", "num_groups_undone", "num_cells_removed"};
    const SEXP values[] = {mo, NULL, NULL};
    const double numbers[] = {(double)ngroups, (double)nremoved};
    SEXP res = mio_r_report_list(3, names, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_meshes_equal(SEXP a, SEXP b, SEXP atol, SEXP rtol, SEXP unordered) {
    int equal = 0;
    mio_r_check(mio_meshes_equal(mio_r_mesh(a), mio_r_mesh(b), mio_r_double(atol, "atol"),
                                 mio_r_double(rtol, "rtol"),
                                 mio_r_bool(unordered, "unordered"), &equal),
                "meshes_equal");
    return Rf_ScalarLogical(equal != 0);
}

SEXP R_mio_diff(SEXP a, SEXP b, SEXP atol, SEXP rtol, SEXP unordered) {
    mio_diff_result *result = NULL;
    mio_r_check(mio_diff(mio_r_mesh(a), mio_r_mesh(b), mio_r_double(atol, "atol"),
                         mio_r_double(rtol, "rtol"), mio_r_bool(unordered, "unordered"),
                         &result),
                "diff");

    static const char *const kVerdicts[] = {"identical", "equal_within_tolerance",
                                            "different"};
    int v = (int)mio_diff_result_verdict(result);
    if (v < 0 || v > 2) v = 2;

    double mabs = 0.0, mrel = 0.0;
    int64_t worst = -1, nex = 0;
    if (mio_diff_result_point_summary(result, &mabs, &mrel, &worst, &nex) != MIO_OK) {
        mio_diff_result_free(result);
        mio_r_fail("diff point summary");
    }

    int64_t nb = mio_diff_result_num_block_diffs(result);
    if (nb < 0) nb = 0;
    SEXP tm = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t)nb));
    SEXP cm = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t)nb));
    SEXP cc = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)nb));
    for (int64_t i = 0; i < nb; ++i) {
        int t = 0, c = 0;
        int64_t k = 0;
        if (mio_diff_result_block(result, i, &t, &c, &k) != MIO_OK) {
            mio_diff_result_free(result);
            UNPROTECT(3);
            mio_r_fail("diff block");
        }
        LOGICAL(tm)[i] = t != 0;
        LOGICAL(cm)[i] = c != 0;
        REAL(cc)[i] = (double)k;
    }
    mio_diff_result_free(result);

    SEXP verdict = PROTECT(Rf_mkString(kVerdicts[v]));
    SEXP sa = PROTECT(Rf_ScalarReal(mabs));
    SEXP sr = PROTECT(Rf_ScalarReal(mrel));
    /* worst_index is a flat index into the coordinate buffer: 1-based here,
     * with the C API's -1 "no offender" becoming 0. */
    SEXP sw = PROTECT(Rf_ScalarReal(worst < 0 ? 0.0 : (double)(worst + 1)));
    SEXP sn = PROTECT(Rf_ScalarReal((double)nex));

    const char *pnames[] = {"max_abs_error", "max_rel_error", "worst_index",
                            "num_exceeding"};
    SEXP pvalues[] = {sa, sr, sw, sn};
    SEXP pts = PROTECT(mio_r_named_list(4, pnames, pvalues));

    const char *bnames[] = {"type_mismatch", "count_mismatch", "conn_mismatch_count"};
    SEXP bvalues[] = {tm, cm, cc};
    SEXP blocks = PROTECT(mio_r_named_list(3, bnames, bvalues));

    const char *names[] = {"verdict", "points", "blocks"};
    SEXP values[] = {verdict, pts, blocks};
    SEXP out = PROTECT(mio_r_named_list(3, names, values));
    UNPROTECT(11); /* tm cm cc verdict sa sr sw sn pts blocks out */
    return out;
}

/* --- statistics --------------------------------------------------------- */

SEXP R_mio_stats(SEXP mesh) {
    mio_stats_report r;
    memset(&r, 0, sizeof(r));
    mio_r_check(mio_stats(mio_r_mesh(mesh), &r), "stats");

    SEXP bmin = PROTECT(Rf_allocVector(REALSXP, 3));
    SEXP bmax = PROTECT(Rf_allocVector(REALSXP, 3));
    SEXP ext = PROTECT(Rf_allocVector(REALSXP, 3));
    SEXP cen = PROTECT(Rf_allocVector(REALSXP, 3));
    memcpy(REAL(bmin), r.bbox_min, 3 * sizeof(double));
    memcpy(REAL(bmax), r.bbox_max, 3 * sizeof(double));
    memcpy(REAL(ext), r.extent, 3 * sizeof(double));
    memcpy(REAL(cen), r.centroid, 3 * sizeof(double));

    const char *names[] = {"num_points", "num_cells",       "bbox_min",        "bbox_max",
                           "extent",     "centroid",        "total_area",      "signed_volume",
                           "unsigned_volume", "num_inverted"};
    const SEXP values[] = {NULL, NULL, bmin, bmax, ext, cen, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)r.num_points, (double)r.num_cells,    r.total_area,
                              r.signed_volume,      r.unsigned_volume,      (double)r.num_inverted};
    SEXP out = mio_r_report_list(10, names, values, numbers);
    UNPROTECT(4);
    return out;
}

SEXP R_mio_compute_bandwidth(SEXP mesh) {
    int64_t bw = mio_compute_bandwidth(mio_r_mesh(mesh));
    if (bw < 0) mio_r_fail("compute_bandwidth");
    return Rf_ScalarReal((double)bw);
}

/* --- opaque-result operations ------------------------------------------- */

SEXP R_mio_reorder(SEXP mesh, SEXP method) {
    mio_reorder_result *r = mio_reorder(mio_r_mesh(mesh), mio_r_string(method, "method"));
    if (r == NULL) mio_r_fail("reorder");

    const void *data = NULL;
    mio_dtype dt;
    int64_t n = 0;
    if (mio_reorder_result_node_perm(r, &data, &dt, &n) != MIO_OK) {
        mio_reorder_result_free(r);
        mio_r_fail("reorder node_perm");
    }
    SEXP node_perm = PROTECT(mio_r_shift_map((const int64_t *)data, n));

    int64_t nperms = mio_reorder_result_num_cell_perms(r);
    if (nperms < 0) nperms = 0;
    SEXP cell_perms = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nperms));
    for (int64_t b = 0; b < nperms; ++b) {
        const void *d = NULL;
        int64_t len = 0;
        if (mio_reorder_result_cell_perm(r, b, &d, &dt, &len) != MIO_OK) {
            mio_reorder_result_free(r);
            UNPROTECT(2);
            mio_r_fail("reorder cell_perm");
        }
        SEXP p = PROTECT(mio_r_shift_map((const int64_t *)d, len));
        SET_VECTOR_ELT(cell_perms, (R_xlen_t)b, p);
        UNPROTECT(1);
    }

    mio_mesh *out = mio_reorder_result_take_mesh(r);
    mio_reorder_result_free(r);
    if (out == NULL) {
        UNPROTECT(2);
        mio_r_fail("reorder take_mesh");
    }
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh", "node_perm", "cell_perms"};
    SEXP values[] = {mo, node_perm, cell_perms};
    SEXP res = PROTECT(mio_r_named_list(3, names, values));
    UNPROTECT(4);
    return res;
}

typedef struct {
    const mio_split_result *result;
    int64_t index;
} split_ctx;

static int64_t split_key_getter(void *c, char *buf, int64_t buflen) {
    split_ctx *x = (split_ctx *)c;
    return mio_split_result_key(x->result, x->index, buf, buflen);
}

SEXP R_mio_split(SEXP mesh, SEXP by, SEXP tag_name) {
    mio_split_result *r = mio_split(mio_r_mesh(mesh), mio_r_string(by, "by"),
                                    mio_r_opt_string(tag_name));
    if (r == NULL) mio_r_fail("split");
    int64_t n = mio_split_result_count(r);
    if (n < 0) {
        mio_split_result_free(r);
        mio_r_fail("split count");
    }
    SEXP pieces = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)n));
    SEXP keys = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)n));
    for (int64_t i = 0; i < n; ++i) {
        split_ctx ctx = {r, i};
        SEXP k = PROTECT(mio_r_getstring(split_key_getter, &ctx, "split key"));
        SET_STRING_ELT(keys, (R_xlen_t)i, STRING_ELT(k, 0));
        UNPROTECT(1);
        /* TAKE, never borrow: the piece must outlive the result handle. */
        mio_mesh *piece = mio_split_result_take_mesh(r, i);
        if (piece == NULL) {
            mio_split_result_free(r);
            UNPROTECT(2);
            mio_r_fail("split take_mesh");
        }
        SEXP p = PROTECT(mio_r_wrap_mesh(piece));
        SET_VECTOR_ELT(pieces, (R_xlen_t)i, p);
        UNPROTECT(1);
    }
    mio_split_result_free(r);
    Rf_setAttrib(pieces, R_NamesSymbol, keys);
    UNPROTECT(2);
    return pieces;
}

/* convert_cells / refine / decimate / decimate_volume share a result shape:
 * one mesh, a point map, and one cell map per input block. `kind` selects
 * which family of C getters to use (0 convert_cells, 1 refine, 2 decimate,
 * 3 decimate_volume). */
static SEXP result_maps(void *r, int kind, SEXP *point_map, SEXP *cell_maps) {
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
    if (st != MIO_OK) return NULL;
    *point_map = PROTECT(mio_r_shift_map((const int64_t *)data, n));

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
    *cell_maps = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nmaps));
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
            UNPROTECT(2);
            return NULL;
        }
        SEXP p = PROTECT(mio_r_shift_map((const int64_t *)d, len));
        SET_VECTOR_ELT(*cell_maps, (R_xlen_t)b, p);
        UNPROTECT(1);
    }
    return *cell_maps; /* non-NULL marker; two PROTECTs are outstanding */
}

/* subdivide's cell maps only -- unlike convert_cells/refine/decimate, there
 * is no point map to read first (subdivide never prunes or renumbers a
 * point), so this does not share result_maps() above. Leaves one PROTECT
 * outstanding on success (the returned list), matching result_maps()'s own
 * convention of leaving its caller to UNPROTECT. */
static SEXP subdivide_cell_maps(mio_subdivide_result *r) {
    int64_t nmaps = mio_subdivide_result_num_cell_maps(r);
    if (nmaps < 0) nmaps = 0;
    SEXP cell_maps = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nmaps));
    for (int64_t b = 0; b < nmaps; ++b) {
        const void *d = NULL;
        mio_dtype dt;
        int64_t len = 0;
        mio_status st = mio_subdivide_result_cell_map(r, b, &d, &dt, &len);
        if (st != MIO_OK) {
            UNPROTECT(1);
            return NULL;
        }
        SEXP p = PROTECT(mio_r_shift_map((const int64_t *)d, len));
        SET_VECTOR_ELT(cell_maps, (R_xlen_t)b, p);
        UNPROTECT(1);
    }
    return cell_maps; /* one PROTECT outstanding */
}

SEXP R_mio_subdivide(SEXP mesh, SEXP record_parent_ids) {
    mio_subdivide_result *r =
        mio_subdivide(mio_r_mesh(mesh), mio_r_bool(record_parent_ids, "record_parent_ids"));
    if (r == NULL) mio_r_fail("subdivide");
    SEXP cm = subdivide_cell_maps(r);
    if (cm == NULL) {
        mio_subdivide_result_free(r);
        mio_r_fail("subdivide maps");
    }
    mio_mesh *out = mio_subdivide_result_take_mesh(r);
    mio_subdivide_result_free(r);
    if (out == NULL) {
        UNPROTECT(1);
        mio_r_fail("subdivide take_mesh");
    }
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh", "cell_maps"};
    SEXP values[] = {mo, cm};
    SEXP res = PROTECT(mio_r_named_list(2, names, values));
    UNPROTECT(3);
    return res;
}

/* agglomerate's cell map is a single FLAT array (unlike subdivide's
 * per-block cell_maps), so it needs no per-block loop at all -- just one
 * mio_r_shift_map() call. */
SEXP R_mio_agglomerate(SEXP mesh, SEXP target_group_size, SEXP merge_coplanar_faces,
                       SEXP coplanar_angle, SEXP min_sphericity) {
    mio_agglomerate_opts opts;
    mio_agglomerate_opts_init(&opts);
    opts.target_group_size = mio_r_int64(target_group_size, "target_group_size");
    opts.merge_coplanar_faces = mio_r_bool(merge_coplanar_faces, "merge_coplanar_faces");
    opts.coplanar_angle = mio_r_double(coplanar_angle, "coplanar_angle");
    opts.min_sphericity = mio_r_double(min_sphericity, "min_sphericity");
    int64_t merged = 0, rejected = 0;
    mio_agglomerate_result *r = mio_agglomerate_ex(mio_r_mesh(mesh), &opts, &merged, &rejected);
    if (r == NULL) mio_r_fail("agglomerate");
    const void *d = NULL;
    mio_dtype dt;
    int64_t len = 0;
    mio_status st = mio_agglomerate_result_cell_map(r, &d, &dt, &len);
    if (st != MIO_OK) {
        mio_agglomerate_result_free(r);
        mio_r_fail("agglomerate cell_map");
    }
    SEXP cm = PROTECT(mio_r_shift_map((const int64_t *)d, len));
    mio_mesh *out = mio_agglomerate_result_take_mesh(r);
    mio_agglomerate_result_free(r);
    if (out == NULL) {
        UNPROTECT(1);
        mio_r_fail("agglomerate take_mesh");
    }
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh", "cell_map", "num_faces_merged", "num_rejected"};
    const SEXP values[] = {mo, cm, NULL, NULL};
    const double numbers[] = {(double)merged, (double)rejected};
    SEXP res = mio_r_report_list(4, names, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_blend_steps(SEXP a, SEXP b, SEXP w, SEXP blend_points) {
    mio_mesh *out = mio_blend_steps(mio_r_mesh(a), mio_r_mesh(b), mio_r_double(w, "w"),
                                    mio_r_bool(blend_points, "blend_points"));
    if (out == NULL) mio_r_fail("blend_steps");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_convert_cells(SEXP mesh, SEXP mode, SEXP record_parent_ids) {
    mio_convert_cells_result *r =
        mio_convert_cells(mio_r_mesh(mesh), mio_r_string(mode, "mode"),
                          mio_r_bool(record_parent_ids, "record_parent_ids"));
    if (r == NULL) mio_r_fail("convert_cells");
    SEXP pm = R_NilValue, cm = R_NilValue;
    if (result_maps(r, 0, &pm, &cm) == NULL) {
        mio_convert_cells_result_free(r);
        mio_r_fail("convert_cells maps");
    }
    mio_mesh *out = mio_convert_cells_result_take_mesh(r);
    mio_convert_cells_result_free(r);
    if (out == NULL) {
        UNPROTECT(2);
        mio_r_fail("convert_cells take_mesh");
    }
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh", "point_map", "cell_maps"};
    SEXP values[] = {mo, pm, cm};
    SEXP res = PROTECT(mio_r_named_list(3, names, values));
    UNPROTECT(4);
    return res;
}

/* tet_to_poly's vertex map is a single FLAT array, like agglomerate's cell
 * map, returned 1-based. */
SEXP R_mio_tet_to_poly(SEXP mesh, SEXP feature_angle, SEXP point_data_to_cell_data,
                       SEXP record_vertex_ids) {
    mio_tet_to_poly_opts opts;
    mio_tet_to_poly_opts_init(&opts);
    opts.feature_angle = mio_r_double(feature_angle, "feature_angle");
    opts.point_data_to_cell_data =
        mio_r_bool(point_data_to_cell_data, "point_data_to_cell_data");
    opts.record_vertex_ids = mio_r_bool(record_vertex_ids, "record_vertex_ids");
    mio_tet_to_poly_result *r = mio_tet_to_poly(mio_r_mesh(mesh), &opts);
    if (r == NULL) mio_r_fail("tet_to_poly");
    const void *d = NULL;
    mio_dtype dt;
    int64_t len = 0, nnew = 0, nbnd = 0, nsplit = 0;
    if (mio_tet_to_poly_result_vertex_map(r, &d, &dt, &len) != MIO_OK ||
        mio_tet_to_poly_result_counts(r, &nnew, &nbnd, &nsplit) != MIO_OK) {
        mio_tet_to_poly_result_free(r);
        mio_r_fail("tet_to_poly vertex_map");
    }
    SEXP vm = PROTECT(mio_r_shift_map((const int64_t *)d, len));
    mio_mesh *out = mio_tet_to_poly_result_take_mesh(r);
    mio_tet_to_poly_result_free(r);
    if (out == NULL) {
        UNPROTECT(1);
        mio_r_fail("tet_to_poly take_mesh");
    }
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP a = PROTECT(Rf_ScalarReal((double)nnew));
    const char *names[] = {"mesh", "vertex_map", "num_new_points", "num_boundary_faces",
                           "num_split_vertices"};
    const SEXP values[] = {mo, vm, a, NULL, NULL};
    const double numbers[] = {(double)nbnd, (double)nsplit};
    SEXP res = mio_r_report_list(5, names, values, numbers);
    UNPROTECT(3);
    return res;
}

/* poly_to_tet returns a convert_cells result, so its maps read the same way. */
SEXP R_mio_poly_to_tet(SEXP mesh, SEXP record_parent_ids) {
    int64_t nexact = 0, nfanned = 0;
    mio_convert_cells_result *r =
        mio_poly_to_tet(mio_r_mesh(mesh), mio_r_bool(record_parent_ids, "record_parent_ids"),
                        &nexact, &nfanned);
    if (r == NULL) mio_r_fail("poly_to_tet");
    SEXP pm = R_NilValue, cm = R_NilValue;
    if (result_maps(r, 0, &pm, &cm) == NULL) {
        mio_convert_cells_result_free(r);
        mio_r_fail("poly_to_tet maps");
    }
    mio_mesh *out = mio_convert_cells_result_take_mesh(r);
    mio_convert_cells_result_free(r);
    if (out == NULL) {
        UNPROTECT(2);
        mio_r_fail("poly_to_tet take_mesh");
    }
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP a = PROTECT(Rf_ScalarReal((double)nexact));
    const char *names[] = {"mesh", "point_map", "cell_maps", "num_exact_tets",
                           "num_fanned_cells"};
    const SEXP values[] = {mo, pm, cm, a, NULL};
    const double numbers[] = {(double)nfanned};
    SEXP res = mio_r_report_list(5, names, values, numbers);
    UNPROTECT(4);
    return res;
}

/* The mio_refine_closure / mio_refine_compare code for a name, or -1. */
static int refine_closure_code(const char *name) {
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

static int refine_compare_code(const char *name) {
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

SEXP R_mio_refine(SEXP mesh, SEXP levels, SEXP record_parent_ids, SEXP cells, SEXP region,
                  SEXP where_array, SEXP where_op, SEXP where_value, SEXP closure,
                  SEXP record_levels, SEXP record_hierarchy) {
    mio_refine_opts opts;
    int64_t *ids = NULL;
    R_xlen_t n_ids = 0;
    mio_refine_result *r = NULL;

    mio_refine_opts_init(&opts);
    opts.levels = (int32_t)mio_r_int(levels, "levels");
    opts.record_parent_ids = mio_r_bool(record_parent_ids, "record_parent_ids");
    opts.record_levels = mio_r_bool(record_levels, "record_levels");
    opts.record_hierarchy = mio_r_bool(record_hierarchy, "record_hierarchy");
    opts.region = mio_r_opt_string(region);
    opts.predicate_array = mio_r_opt_string(where_array);
    if (where_value != R_NilValue && Rf_length(where_value) > 0)
        opts.predicate_value = Rf_asReal(where_value);
    opts.closure = refine_closure_code(mio_r_opt_string(closure));
    if (opts.closure < 0) Rf_error("meshio++: refine: unknown closure");
    opts.predicate_op = refine_compare_code(mio_r_opt_string(where_op));
    if (opts.predicate_op < 0) Rf_error("meshio++: refine: unknown comparison");

    if (cells != R_NilValue && Rf_length(cells) > 0) {
        /* R has no native int64, so cell indices arrive as double (exact to
         * 2^53) or integer; either way they are 1-based on this side. */
        n_ids = Rf_length(cells);
        ids = (int64_t *)R_alloc((size_t)n_ids, sizeof(int64_t));
        for (R_xlen_t i = 0; i < n_ids; ++i) {
            double v = (TYPEOF(cells) == REALSXP) ? REAL(cells)[i] : (double)INTEGER(cells)[i];
            ids[i] = (int64_t)v - 1;
        }
        opts.cells = ids;
        opts.num_cells = (int64_t)n_ids;
    }

    r = mio_refine_ex(mio_r_mesh(mesh), &opts);
    if (r == NULL) mio_r_fail("refine");
    SEXP pm = R_NilValue, cm = R_NilValue;
    if (result_maps(r, 1, &pm, &cm) == NULL) {
        mio_refine_result_free(r);
        mio_r_fail("refine maps");
    }
    mio_mesh *out = mio_refine_result_take_mesh(r);
    mio_refine_result_free(r);
    if (out == NULL) {
        UNPROTECT(2);
        mio_r_fail("refine take_mesh");
    }
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh", "point_map", "cell_maps"};
    SEXP values[] = {mo, pm, cm};
    SEXP res = PROTECT(mio_r_named_list(3, names, values));
    UNPROTECT(4);
    return res;
}

SEXP R_mio_tessellate(SEXP mesh, SEXP levels, SEXP curved, SEXP fields, SEXP record_stencil) {
    mio_tessellate_opts opts;
    mio_tessellate_report report;
    mio_mesh *out;

    mio_tessellate_opts_init(&opts);
    opts.levels = mio_r_int64(levels, "levels");
    opts.curved = mio_r_bool(curved, "curved");
    opts.fields = mio_r_bool(fields, "fields");
    opts.record_stencil = mio_r_bool(record_stencil, "record_stencil");

    out = mio_tessellate(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("tessellate");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
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
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)report.num_source_points, (double)report.num_source_cells,
                              (double)report.num_curved_source_cells,
                              (double)report.num_pass_through_source_cells,
                              (double)report.num_points, (double)report.num_cells,
                              (double)report.num_boundary_facets,
                              (double)report.num_interior_facets,
                              (double)report.num_facets_with_bad_count};
    SEXP res = mio_r_report_list(10, names, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_decimate(SEXP mesh, SEXP target_ratio, SEXP target_faces, SEXP max_error,
                    SEXP placement, SEXP preserve_boundary, SEXP preserve_features,
                    SEXP feature_angle, SEXP frozen) {
    mio_decimate_opts opts;
    R_xlen_t n_frozen = 0;

    mio_decimate_opts_init(&opts);
    opts.target_ratio = mio_r_double(target_ratio, "target_ratio");
    opts.target_faces = mio_r_int64(target_faces, "target_faces");
    opts.max_error = mio_r_double(max_error, "max_error");
    opts.placement = mio_r_opt_string(placement);
    opts.preserve_boundary = mio_r_bool(preserve_boundary, "preserve_boundary");
    opts.preserve_features = mio_r_bool(preserve_features, "preserve_features");
    opts.feature_angle = mio_r_double(feature_angle, "feature_angle");
    opts.frozen = r_ids_0based(frozen, &n_frozen);
    opts.num_frozen = (int64_t)n_frozen;

    mio_decimate_result *r = mio_decimate_ex(mio_r_mesh(mesh), &opts);
    if (r == NULL) mio_r_fail("decimate");
    SEXP pm = R_NilValue, cm = R_NilValue;
    if (result_maps(r, 2, &pm, &cm) == NULL) {
        mio_decimate_result_free(r);
        mio_r_fail("decimate maps");
    }
    const double numbers[] = {(double)mio_decimate_result_faces_removed(r),
                              (double)mio_decimate_result_points_removed(r),
                              (double)mio_decimate_result_collapses_rejected(r),
                              mio_decimate_result_max_error_applied(r)};
    mio_mesh *out = mio_decimate_result_take_mesh(r);
    mio_decimate_result_free(r);
    if (out == NULL) {
        UNPROTECT(2);
        mio_r_fail("decimate take_mesh");
    }
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh",           "point_map",          "cell_maps",
                           "faces_removed",  "points_removed",     "collapses_rejected",
                           "max_error_applied"};
    const SEXP values[] = {mo, pm, cm, NULL, NULL, NULL, NULL};
    SEXP res = mio_r_report_list(7, names, values, numbers);
    UNPROTECT(3);
    return res;
}

/* decimate_volume(mesh, ratio, target_cells, max_error, placement,
 *                 preserve_boundary, preserve_features, feature_angle, frozen) --
 * the tet counterpart of R_mio_decimate, through mio_decimate_volume_ex so the
 * `frozen` pin list crosses too. */
SEXP R_mio_decimate_volume(SEXP mesh, SEXP target_ratio, SEXP target_cells, SEXP max_error,
                           SEXP placement, SEXP preserve_boundary, SEXP preserve_features,
                           SEXP feature_angle, SEXP frozen) {
    mio_decimate_volume_opts opts;
    R_xlen_t n_frozen = 0;

    mio_decimate_volume_opts_init(&opts);
    opts.target_ratio = mio_r_double(target_ratio, "target_ratio");
    opts.target_cells = mio_r_int64(target_cells, "target_cells");
    opts.max_error = mio_r_double(max_error, "max_error");
    opts.placement = mio_r_opt_string(placement);
    opts.preserve_boundary = mio_r_bool(preserve_boundary, "preserve_boundary");
    opts.preserve_features = mio_r_bool(preserve_features, "preserve_features");
    opts.feature_angle = mio_r_double(feature_angle, "feature_angle");
    opts.frozen = r_ids_0based(frozen, &n_frozen);
    opts.num_frozen = (int64_t)n_frozen;

    mio_decimate_volume_result *r = mio_decimate_volume_ex(mio_r_mesh(mesh), &opts);
    if (r == NULL) mio_r_fail("decimate_volume");
    SEXP pm = R_NilValue, cm = R_NilValue;
    if (result_maps(r, 3, &pm, &cm) == NULL) {
        mio_decimate_volume_result_free(r);
        mio_r_fail("decimate_volume maps");
    }
    const double numbers[] = {(double)mio_decimate_volume_result_tets_removed(r),
                              (double)mio_decimate_volume_result_points_removed(r),
                              (double)mio_decimate_volume_result_collapses_rejected(r),
                              mio_decimate_volume_result_max_error_applied(r)};
    mio_mesh *out = mio_decimate_volume_result_take_mesh(r);
    mio_decimate_volume_result_free(r);
    if (out == NULL) {
        UNPROTECT(2);
        mio_r_fail("decimate_volume take_mesh");
    }
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh",          "point_map",          "cell_maps",
                           "tets_removed",  "points_removed",     "collapses_rejected",
                           "max_error_applied"};
    const SEXP values[] = {mo, pm, cm, NULL, NULL, NULL, NULL};
    SEXP res = mio_r_report_list(7, names, values, numbers);
    UNPROTECT(3);
    return res;
}

SEXP R_mio_partition(SEXP mesh, SEXP nparts, SEXP method, SEXP imbalance, SEXP mode,
                     SEXP seed, SEXP record_ids, SEXP ghost_layers, SEXP weights_key) {
    mio_partition_result *r =
        mio_partition(mio_r_mesh(mesh), mio_r_int(nparts, "nparts"),
                      mio_r_opt_string(method), mio_r_double(imbalance, "imbalance"),
                      mio_r_opt_string(mode), mio_r_int(seed, "seed"),
                      mio_r_bool(record_ids, "record_ids"),
                      mio_r_int(ghost_layers, "ghost_layers"), mio_r_opt_string(weights_key));
    if (r == NULL) mio_r_fail("partition");

    int64_t n = mio_partition_result_num_pieces(r);
    int64_t nmaps = mio_partition_result_num_cell_maps(r);
    if (n < 0) n = 0;
    if (nmaps < 0) nmaps = 0;

    SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)n));
    for (int64_t i = 0; i < n; ++i) {
        const void *d = NULL;
        mio_dtype dt;
        int64_t len = 0;
        if (mio_partition_result_point_map(r, i, &d, &dt, &len) != MIO_OK) {
            mio_partition_result_free(r);
            UNPROTECT(1);
            mio_r_fail("partition point_map");
        }
        SEXP pm = PROTECT(mio_r_shift_map((const int64_t *)d, len));
        SEXP cms = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nmaps));
        for (int64_t b = 0; b < nmaps; ++b) {
            const void *d2 = NULL;
            int64_t l2 = 0;
            if (mio_partition_result_cell_map(r, i, b, &d2, &dt, &l2) != MIO_OK) {
                mio_partition_result_free(r);
                UNPROTECT(2); /* pm, cms */
                mio_r_fail("partition cell_map");
            }
            SEXP p = PROTECT(mio_r_shift_map((const int64_t *)d2, l2));
            SET_VECTOR_ELT(cms, (R_xlen_t)b, p);
            UNPROTECT(1);
        }
        SEXP pid = PROTECT(Rf_ScalarInteger(mio_partition_result_part_id(r, i)));
        mio_mesh *piece = mio_partition_result_take_mesh(r, i);
        if (piece == NULL) {
            mio_partition_result_free(r);
            UNPROTECT(4);
            mio_r_fail("partition take_mesh");
        }
        SEXP mo = PROTECT(mio_r_wrap_mesh(piece));
        const char *names[] = {"part_id", "mesh", "point_map", "cell_maps"};
        SEXP values[] = {pid, mo, pm, cms};
        SEXP item = PROTECT(mio_r_named_list(4, names, values));
        SET_VECTOR_ELT(out, (R_xlen_t)i, item);
        UNPROTECT(5);
    }
    mio_partition_result_free(r);
    UNPROTECT(1);
    return out;
}

SEXP R_mio_partition_labels(SEXP mesh, SEXP nparts, SEXP method, SEXP imbalance, SEXP mode,
                            SEXP seed, SEXP weights_key, SEXP num_cells) {
    int64_t n = mio_r_int64(num_cells, "num_cells");
    SEXP buf = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)((n > 0 ? n : 1) * sizeof(int64_t))));
    int64_t *labels = (int64_t *)RAW(buf);
    mio_status st = mio_partition_labels(
        mio_r_mesh(mesh), mio_r_int(nparts, "nparts"), mio_r_opt_string(method),
        mio_r_double(imbalance, "imbalance"), mio_r_opt_string(mode), mio_r_int(seed, "seed"),
        mio_r_opt_string(weights_key), labels, n);
    if (st != MIO_OK) {
        UNPROTECT(1);
        mio_r_fail("partition_labels");
    }
    SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)n));
    /* Part IDs, not indices: deliberately NOT shifted to 1-based. */
    for (int64_t i = 0; i < n; ++i) REAL(out)[i] = (double)labels[i];
    UNPROTECT(2);
    return out;
}

/* --- data operations ---------------------------------------------------- */

SEXP R_mio_sets_to_data(SEXP mesh, SEXP location, SEXP data_name, SEXP join_char, SEXP order) {
    SEXP shelter;
    int64_t count = 0;
    const char* const* nm = mio_r_names(order, &count, &shelter);
    mio_mesh* out =
        mio_sets_to_data(mio_r_mesh(mesh), (mio_data_location)mio_r_int(location, "location"),
                         mio_r_opt_string(data_name), mio_r_opt_string(join_char), nm, count);
    UNPROTECT(1);
    if (out == NULL)
        mio_r_fail("sets_to_data");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_data_to_sets(SEXP mesh, SEXP location, SEXP key) {
    mio_mesh* out =
        mio_data_to_sets(mio_r_mesh(mesh), (mio_data_location)mio_r_int(location, "location"),
                         mio_r_string(key, "key"));
    if (out == NULL)
        mio_r_fail("data_to_sets");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_data_drop(SEXP mesh, SEXP location, SEXP names, SEXP ignore_missing) {
    SEXP shelter;
    int64_t count = 0;
    const char *const *nm = mio_r_names(names, &count, &shelter);
    mio_mesh *out = mio_data_drop(mio_r_mesh(mesh),
                                  (mio_data_location)mio_r_int(location, "location"), nm,
                                  count, mio_r_bool(ignore_missing, "ignore_missing"));
    UNPROTECT(1);
    if (out == NULL) mio_r_fail("data_drop");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_data_keep(SEXP mesh, SEXP location, SEXP names, SEXP ignore_missing) {
    SEXP shelter;
    int64_t count = 0;
    const char *const *nm = mio_r_names(names, &count, &shelter);
    mio_mesh *out = mio_data_keep(mio_r_mesh(mesh),
                                  (mio_data_location)mio_r_int(location, "location"), nm,
                                  count, mio_r_bool(ignore_missing, "ignore_missing"));
    UNPROTECT(1);
    if (out == NULL) mio_r_fail("data_keep");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_data_rename(SEXP mesh, SEXP location, SEXP from, SEXP to) {
    mio_mesh *out = mio_data_rename(mio_r_mesh(mesh),
                                    (mio_data_location)mio_r_int(location, "location"),
                                    mio_r_string(from, "from"), mio_r_string(to, "to"));
    if (out == NULL) mio_r_fail("data_rename");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_data_point_to_cell(SEXP mesh, SEXP names, SEXP suffix) {
    SEXP shelter;
    int64_t count = 0;
    const char *const *nm = mio_r_names(names, &count, &shelter);
    mio_mesh *out = mio_data_point_to_cell(mio_r_mesh(mesh), nm, count,
                                           mio_r_opt_string(suffix));
    UNPROTECT(1);
    if (out == NULL) mio_r_fail("data_point_to_cell");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_data_cell_to_point(SEXP mesh, SEXP names, SEXP weight, SEXP suffix) {
    SEXP shelter;
    int64_t count = 0;
    const char *const *nm = mio_r_names(names, &count, &shelter);
    mio_mesh *out = mio_data_cell_to_point(
        mio_r_mesh(mesh), nm, count, (mio_cell_point_weight)mio_r_int(weight, "weight"),
        mio_r_opt_string(suffix));
    UNPROTECT(1);
    if (out == NULL) mio_r_fail("data_cell_to_point");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_data_calc(SEXP mesh, SEXP expression, SEXP location, SEXP output_name,
                     SEXP overwrite) {
    mio_mesh *out = mio_data_calc(mio_r_mesh(mesh), mio_r_string(expression, "expression"),
                                  (mio_data_location)mio_r_int(location, "location"),
                                  mio_r_string(output_name, "output_name"),
                                  mio_r_bool(overwrite, "overwrite"));
    if (out == NULL) mio_r_fail("data_calc");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_data_condition(SEXP mesh, SEXP location, SEXP names, SEXP mode, SEXP lo, SEXP hi,
                          SEXP scope, SEXP nan_policy, SEXP nan_replacement, SEXP suffix) {
    SEXP shelter;
    int64_t count = 0;
    const char *const *nm = mio_r_names(names, &count, &shelter);
    mio_mesh *out = mio_data_condition(
        mio_r_mesh(mesh), (mio_data_location)mio_r_int(location, "location"), nm, count,
        (mio_condition_mode)mio_r_int(mode, "mode"), mio_r_double(lo, "lo"),
        mio_r_double(hi, "hi"), (mio_condition_scope)mio_r_int(scope, "scope"),
        (mio_nan_policy)mio_r_int(nan_policy, "nan_policy"),
        mio_r_double(nan_replacement, "nan_replacement"), mio_r_opt_string(suffix));
    UNPROTECT(1);
    if (out == NULL) mio_r_fail("data_condition");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_tensor_invariants(SEXP mesh, SEXP location, SEXP names, SEXP outputs, SEXP prefix,
                             SEXP suffix, SEXP overwrite) {
    SEXP shelter;
    int64_t count = 0;
    const char *const *nm = mio_r_names(names, &count, &shelter);
    mio_mesh *out = mio_tensor_invariants(
        mio_r_mesh(mesh), (mio_data_location)mio_r_int(location, "location"), nm, count,
        mio_r_int(outputs, "outputs"), mio_r_opt_string(prefix), mio_r_opt_string(suffix),
        mio_r_bool(overwrite, "overwrite"));
    UNPROTECT(1);
    if (out == NULL) mio_r_fail("tensor_invariants");
    return mio_r_wrap_mesh(out);
}

typedef struct {
    const mio_data_info *info;
    int64_t index;
} info_ctx;

static int64_t info_name_getter(void *c, char *buf, int64_t buflen) {
    info_ctx *x = (info_ctx *)c;
    return mio_data_info_name(x->info, x->index, buf, buflen);
}

SEXP R_mio_data_info(SEXP mesh) {
    mio_data_info *info =
        (mio_data_info *)mio_r_check_ptr(mio_data_info_create(mio_r_mesh(mesh)), "data_info");
    int64_t n = mio_data_info_count(info);
    if (n < 0) {
        mio_data_info_free(info);
        mio_r_fail("data_info count");
    }
    static const char *const kLocations[] = {"point", "cell", "field"};

    SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)n));
    for (int64_t i = 0; i < n; ++i) {
        info_ctx ctx = {info, i};
        SEXP name = PROTECT(mio_r_getstring(info_name_getter, &ctx, "data_info name"));
        mio_data_array_info e;
        memset(&e, 0, sizeof(e));
        if (mio_data_info_entry(info, i, &e) != MIO_OK) {
            mio_data_info_free(info);
            UNPROTECT(1); /* name */
            mio_r_fail("data_info entry");
        }
        SEXP loc = PROTECT(Rf_mkString(
            (e.location >= 0 && e.location <= 2) ? kLocations[e.location] : "unknown"));
        SEXP dtype = PROTECT(Rf_mkString(mio_r_dtype_name((mio_dtype)e.dtype)));
        SEXP nb = PROTECT(Rf_ScalarReal((double)e.num_blocks));
        SEXP ne = PROTECT(Rf_ScalarReal((double)e.num_entries));
        SEXP ncomp = PROTECT(Rf_ScalarReal((double)e.num_components));
        SEXP nv = PROTECT(Rf_ScalarReal((double)e.num_values));
        SEXP mn = PROTECT(Rf_ScalarReal(e.min));
        SEXP mx = PROTECT(Rf_ScalarReal(e.max));
        SEXP mean = PROTECT(Rf_ScalarReal(e.mean));
        SEXP nnan = PROTECT(Rf_ScalarReal((double)e.num_nan));
        SEXP ninf = PROTECT(Rf_ScalarReal((double)e.num_inf));
        SEXP nfin = PROTECT(Rf_ScalarReal((double)e.num_finite));
        SEXP incons = PROTECT(Rf_ScalarLogical(e.inconsistent_blocks != 0));

        int64_t nc = e.num_components > 0 ? e.num_components : 0;
        SEXP comps = PROTECT(Rf_allocMatrix(REALSXP, (int)nc, 3));
        for (int64_t c = 0; c < nc; ++c) {
            double lo = 0.0, hi = 0.0, mu = 0.0;
            if (mio_data_info_component(info, i, c, &lo, &hi, &mu) != MIO_OK) {
                mio_data_info_free(info);
                UNPROTECT(15); /* this entry's 15 values */
                mio_r_fail("data_info component");
            }
            REAL(comps)[c] = lo;
            REAL(comps)[nc + c] = hi;
            REAL(comps)[2 * nc + c] = mu;
        }
        SEXP cdimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SEXP ccols = PROTECT(Rf_allocVector(STRSXP, 3));
        SET_STRING_ELT(ccols, 0, Rf_mkChar("min"));
        SET_STRING_ELT(ccols, 1, Rf_mkChar("max"));
        SET_STRING_ELT(ccols, 2, Rf_mkChar("mean"));
        SET_VECTOR_ELT(cdimnames, 0, R_NilValue);
        SET_VECTOR_ELT(cdimnames, 1, ccols);
        Rf_setAttrib(comps, R_DimNamesSymbol, cdimnames);
        UNPROTECT(2); /* cdimnames, ccols */

        const char *names[] = {"name",        "location",   "dtype",     "num_blocks",
                               "num_entries", "num_components", "num_values", "min",
                               "max",         "mean",       "num_nan",   "num_inf",
                               "num_finite",  "inconsistent_blocks", "components"};
        SEXP values[] = {name, loc,  dtype, nb,   ne,     ncomp, nv, mn,
                         mx,   mean, nnan,  ninf, nfin,   incons, comps};
        SEXP item = PROTECT(mio_r_named_list(15, names, values));
        SET_VECTOR_ELT(out, (R_xlen_t)i, item);
        UNPROTECT(16); /* the 15 values + item */
    }
    mio_data_info_free(info);
    UNPROTECT(1);
    return out;
}

typedef struct {
    const mio_data_integrate *result;
    int64_t index;
} integrate_name_ctx;

typedef struct {
    const mio_data_integrate *result;
    int64_t index;
    int64_t region;
} integrate_region_name_ctx;

static int64_t integrate_name_getter(void *c, char *buf, int64_t buflen) {
    integrate_name_ctx *x = (integrate_name_ctx *)c;
    return mio_data_integrate_name(x->result, x->index, buf, buflen);
}

static int64_t integrate_region_name_getter(void *c, char *buf, int64_t buflen) {
    integrate_region_name_ctx *x = (integrate_region_name_ctx *)c;
    return mio_data_integrate_region_name(x->result, x->index, x->region, buf, buflen);
}

/* Build the (num_components x 4) total/mean/domain_measure/num_nan matrix for
 * one array's whole-mesh entry (region < 0) or one of its named Cell
 * regions (region >= 0). Frees `result` and raises via mio_r_fail (which
 * longjmps, unwinding the protection stack) on any C API failure -- the
 * caller never sees a NULL back. Returns unprotected, like mio_r_named_list;
 * the caller must PROTECT it immediately. */
static SEXP integrate_components_matrix(mio_data_integrate *result, int64_t index, int64_t region,
                                        int64_t num_components) {
    int64_t nc = num_components > 0 ? num_components : 0;
    SEXP m = PROTECT(Rf_allocMatrix(REALSXP, (int)nc, 4));
    for (int64_t c = 0; c < nc; ++c) {
        double total = 0.0, mean = 0.0, dm = 0.0;
        int64_t nan = 0;
        mio_status s = (region < 0) ? mio_data_integrate_component(result, index, c, &total,
                                                                    &mean, &dm, &nan)
                                    : mio_data_integrate_region_component(
                                          result, index, region, c, &total, &mean, &dm, &nan);
        if (s != MIO_OK) {
            mio_data_integrate_free(result);
            UNPROTECT(1); /* m */
            mio_r_fail("data_integrate component");
        }
        REAL(m)[c] = total;
        REAL(m)[nc + c] = mean;
        REAL(m)[2 * nc + c] = dm;
        REAL(m)[3 * nc + c] = (double)nan;
    }
    SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP cols = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(cols, 0, Rf_mkChar("total"));
    SET_STRING_ELT(cols, 1, Rf_mkChar("mean"));
    SET_STRING_ELT(cols, 2, Rf_mkChar("domain_measure"));
    SET_STRING_ELT(cols, 3, Rf_mkChar("num_nan"));
    SET_VECTOR_ELT(dimnames, 0, R_NilValue);
    SET_VECTOR_ELT(dimnames, 1, cols);
    Rf_setAttrib(m, R_DimNamesSymbol, dimnames);
    UNPROTECT(3); /* dimnames, cols, m */
    return m;
}

/* Cell-measure-weighted total/mean of one or more cell_data arrays --
 * gradient's integration counterpart. `names` NULL means every cell_data
 * array. Reported for the whole mesh (`domain`) and independently for every
 * named Cell region (`regions`) -- a cell in two regions contributes fully
 * to both, one in none contributes to neither. A point_data-only name fails,
 * naming data_point_to_cell as the fix. See doc/field_integration.md. */
SEXP R_mio_data_integrate(SEXP mesh, SEXP names) {
    SEXP shelter;
    int64_t ncount = 0;
    const char *const *nm = mio_r_names(names, &ncount, &shelter);
    mio_data_integrate *result = (mio_data_integrate *)mio_r_check_ptr(
        mio_data_integrate_create(mio_r_mesh(mesh), nm, ncount), "data_integrate");
    UNPROTECT(1); /* shelter */

    int64_t n = mio_data_integrate_count(result);
    if (n < 0) {
        mio_data_integrate_free(result);
        mio_r_fail("data_integrate count");
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)n));
    for (int64_t i = 0; i < n; ++i) {
        integrate_name_ctx nctx = {result, i};
        SEXP name =
            PROTECT(mio_r_getstring(integrate_name_getter, &nctx, "data_integrate name"));

        mio_field_integral_info e;
        memset(&e, 0, sizeof(e));
        if (mio_data_integrate_entry(result, i, &e) != MIO_OK) {
            mio_data_integrate_free(result);
            UNPROTECT(2); /* name, out */
            mio_r_fail("data_integrate entry");
        }
        SEXP ncomp = PROTECT(Rf_ScalarReal((double)e.num_components));
        SEXP dcells = PROTECT(Rf_ScalarReal((double)e.num_cells));
        SEXP dskip = PROTECT(Rf_ScalarReal((double)e.num_skipped));
        SEXP dcomps = PROTECT(integrate_components_matrix(result, i, -1, e.num_components));

        const char *domain_names[] = {"num_cells", "num_skipped", "components"};
        SEXP domain_values[] = {dcells, dskip, dcomps};
        SEXP domain = PROTECT(mio_r_named_list(3, domain_names, domain_values));

        int64_t nregions = mio_data_integrate_region_count(result, i);
        if (nregions < 0) {
            mio_data_integrate_free(result);
            UNPROTECT(7); /* name, ncomp, dcells, dskip, dcomps, domain, out */
            mio_r_fail("data_integrate region count");
        }

        SEXP regions = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nregions));
        for (int64_t r = 0; r < nregions; ++r) {
            integrate_region_name_ctx rctx = {result, i, r};
            SEXP rname = PROTECT(mio_r_getstring(integrate_region_name_getter, &rctx,
                                                 "data_integrate region name"));

            mio_field_integral_info re;
            memset(&re, 0, sizeof(re));
            if (mio_data_integrate_region_entry(result, i, r, &re) != MIO_OK) {
                mio_data_integrate_free(result);
                UNPROTECT(9); /* rname + the 7 above + regions */
                mio_r_fail("data_integrate region entry");
            }
            SEXP rcells = PROTECT(Rf_ScalarReal((double)re.num_cells));
            SEXP rskip = PROTECT(Rf_ScalarReal((double)re.num_skipped));
            SEXP rcomps = PROTECT(integrate_components_matrix(result, i, r, re.num_components));

            const char *region_names[] = {"name", "num_cells", "num_skipped", "components"};
            SEXP region_values[] = {rname, rcells, rskip, rcomps};
            SEXP region_item = PROTECT(mio_r_named_list(4, region_names, region_values));
            SET_VECTOR_ELT(regions, (R_xlen_t)r, region_item);
            UNPROTECT(5); /* rname, rcells, rskip, rcomps, region_item */
        }

        const char *array_names[] = {"name", "num_components", "domain", "regions"};
        SEXP array_values[] = {name, ncomp, domain, regions};
        SEXP item = PROTECT(mio_r_named_list(4, array_names, array_values));
        SET_VECTOR_ELT(out, (R_xlen_t)i, item);
        UNPROTECT(8); /* name, ncomp, dcells, dskip, dcomps, domain, regions, item */
    }
    mio_data_integrate_free(result);
    UNPROTECT(1); /* out */
    return out;
}

/* --- settings pipeline (v9.11.0) ---------------------------------------- */

SEXP R_mio_pipeline_run_file(SEXP settings_path) {
    mio_r_check(mio_pipeline_run_file(mio_r_string(settings_path, "settings_path")),
                "pipeline_run_file");
    return R_NilValue;
}

static void r_pipeline_report_finalizer(SEXP handle) {
    mio_pipeline_report_free((mio_pipeline_report *)R_ExternalPtrAddr(handle));
    R_ClearExternalPtr(handle);
}

static SEXP r_pipeline_report(mio_pipeline_report *report) {
    if (report == NULL) mio_r_fail("pipeline_report");
    SEXP owner = PROTECT(R_MakeExternalPtr(report, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(owner, r_pipeline_report_finalizer, TRUE);
    int64_t n = mio_pipeline_report_json(report, NULL, 0);
    if (n < 0) mio_r_fail("pipeline_report");
    char *buf = (char *)R_alloc((size_t)n + 1, 1);
    if (mio_pipeline_report_json(report, buf, n + 1) < 0) mio_r_fail("pipeline_report");
    SEXP out = PROTECT(Rf_ScalarString(Rf_mkCharCE(buf, CE_UTF8)));
    r_pipeline_report_finalizer(owner);
    UNPROTECT(2);
    return out;
}

SEXP R_mio_pipeline_run_file_report(SEXP path) {
    return r_pipeline_report(mio_pipeline_run_file_report(mio_r_string(path, "path")));
}
SEXP R_mio_pipeline_run_json_report(SEXP text) {
    return r_pipeline_report(mio_pipeline_run_json_report(mio_r_string(text, "text")));
}
SEXP R_mio_sequence_pipeline_run_file_report(SEXP path) {
    return r_pipeline_report(mio_sequence_pipeline_run_file_report(mio_r_string(path, "path")));
}
SEXP R_mio_sequence_pipeline_run_json_report(SEXP text) {
    return r_pipeline_report(mio_sequence_pipeline_run_json_report(mio_r_string(text, "text")));
}

SEXP R_mio_pipeline_run_json(SEXP json_text) {
    mio_r_check(mio_pipeline_run_json(mio_r_string(json_text, "json_text")),
                "pipeline_run_json");
    return R_NilValue;
}

SEXP R_mio_pipeline_has_json(void) {
    return Rf_ScalarLogical(mio_pipeline_has_json() != 0);
}

/* --- regular grids and signed distance ---------------------------------- */

/* R has no int64, so the counts come back as doubles (exact well past any
 * plausible cell count) -- the convention every other counter here follows. */

static void fill_sdf_opts(mio_sdf_opts *o, SEXP sign, SEXP location, SEXP band,
                          SEXP record_inside, SEXP watertight_check) {
    mio_sdf_opts_init(o);
    o->sign = (int32_t)mio_r_int(sign, "sign");
    o->location = (int32_t)mio_r_int(location, "location");
    o->band = mio_r_double(band, "band");
    o->record_inside = mio_r_bool(record_inside, "record_inside");
    o->watertight_check = (int32_t)mio_r_int(watertight_check, "watertight_check");
}

static SEXP quality_list(const mio_surface_quality *q) {
    SEXP w = PROTECT(Rf_ScalarLogical(q->watertight != 0));
    const char *names[] = {"boundary_edges", "non_manifold_edges", "inconsistent_pairs",
                           "degenerate_triangles", "watertight"};
    const SEXP values[] = {NULL, NULL, NULL, NULL, w};
    const double numbers[] = {(double)q->boundary_edges, (double)q->non_manifold_edges,
                              (double)q->inconsistent_pairs, (double)q->degenerate_triangles};
    SEXP out = mio_r_report_list(5, names, values, numbers);
    UNPROTECT(1);
    return out;
}

SEXP R_mio_grid(SEXP dims, SEXP origin, SEXP spacing, SEXP max_cells) {
    SEXP d = PROTECT(Rf_coerceVector(dims, REALSXP));
    if (Rf_xlength(d) != 3) {
        UNPROTECT(1);
        Rf_error("`dims` must have exactly 3 elements");
    }
    int64_t cdims[3];
    for (int i = 0; i < 3; ++i) cdims[i] = (int64_t)REAL(d)[i];
    UNPROTECT(1);

    double o[3], s[3];
    mio_r_vec3(origin, o, "origin");
    mio_r_vec3(spacing, s, "spacing");
    mio_mesh *out = mio_grid(cdims, o, s, mio_r_int64(max_cells, "max_cells"));
    if (out == NULL) mio_r_fail("grid");
    return mio_r_wrap_mesh(out);
}

/* The parametric generators share one options struct; `shape` picks the entry
 * point. cell_type and split arrive as names and are mapped here so the R layer
 * and the other language bindings refuse the same spellings. */
static int32_t mio_r_gen_cell(const char *name) {
    if (name[0] == '\0') return MIO_GEN_DEFAULT;
    if (strcmp(name, "quad") == 0) return MIO_GEN_QUAD;
    if (strcmp(name, "triangle") == 0) return MIO_GEN_TRIANGLE;
    if (strcmp(name, "hexahedron") == 0) return MIO_GEN_HEXAHEDRON;
    if (strcmp(name, "tetra") == 0) return MIO_GEN_TETRA;
    if (strcmp(name, "wedge") == 0) return MIO_GEN_WEDGE;
    Rf_error("unknown cell_type '%s'", name);
    return 0;
}

static int32_t mio_r_gen_split(const char *name) {
    if (strcmp(name, "zigzag") == 0) return MIO_GEN_SPLIT_ZIGZAG;
    if (strcmp(name, "up") == 0) return MIO_GEN_SPLIT_UP;
    if (strcmp(name, "down") == 0) return MIO_GEN_SPLIT_DOWN;
    if (strcmp(name, "center") == 0) return MIO_GEN_SPLIT_CENTER;
    Rf_error("unknown split '%s'", name);
    return 0;
}

SEXP R_mio_generate(SEXP shape, SEXP counts, SEXP size, SEXP origin, SEXP radius, SEXP height,
                    SEXP cell_type, SEXP split, SEXP order, SEXP max_cells, SEXP radius2,
                    SEXP angle, SEXP num_twists, SEXP complete, SEXP ratio, SEXP seed) {
    mio_gen_opts opts;
    mio_gen_opts_init(&opts);
    SEXP c = PROTECT(Rf_coerceVector(counts, REALSXP));
    if (Rf_xlength(c) < 1 || Rf_xlength(c) > 3) {
        UNPROTECT(1);
        Rf_error("`counts` must have 1 to 3 elements");
    }
    for (R_xlen_t i = 0; i < Rf_xlength(c); ++i) opts.counts[i] = (int64_t)REAL(c)[i];
    for (R_xlen_t i = Rf_xlength(c); i < 3; ++i) opts.counts[i] = 1;
    UNPROTECT(1);
    mio_r_vec3(size, opts.size, "size");
    mio_r_vec3(origin, opts.origin, "origin");
    opts.radius = mio_r_double(radius, "radius");
    opts.height = mio_r_double(height, "height");
    opts.radius2 = mio_r_double(radius2, "radius2");
    opts.angle = mio_r_double(angle, "angle");
    opts.cell_type = mio_r_gen_cell(mio_r_string(cell_type, "cell_type"));
    opts.split = mio_r_gen_split(mio_r_string(split, "split"));
    opts.order = (int32_t)mio_r_int(order, "order");
    opts.max_cells = mio_r_int64(max_cells, "max_cells");
    opts.num_twists = (int32_t)mio_r_int(num_twists, "num_twists");
    opts.complete = mio_r_bool(complete, "complete") ? 1 : 0;
    opts.ratio = mio_r_double(ratio, "ratio");
    opts.seed = (uint64_t)mio_r_double(seed, "seed");

    const char *name = mio_r_string(shape, "shape");
    mio_mesh *out = NULL;
    if (strcmp(name, "rectangle") == 0) out = mio_rectangle(&opts);
    else if (strcmp(name, "box") == 0) out = mio_box(&opts);
    else if (strcmp(name, "disk") == 0) out = mio_disk(&opts);
    else if (strcmp(name, "uv_sphere") == 0) out = mio_uv_sphere(&opts);
    else if (strcmp(name, "icosa_sphere") == 0) out = mio_icosa_sphere(&opts);
    else if (strcmp(name, "cylinder") == 0) out = mio_cylinder(&opts);
    else if (strcmp(name, "reference_triangle") == 0) out = mio_reference_triangle(&opts);
    else if (strcmp(name, "ngon") == 0) out = mio_ngon(&opts);
    else if (strcmp(name, "annulus") == 0) out = mio_annulus(&opts);
    else if (strcmp(name, "lshape") == 0) out = mio_lshape(&opts);
    else if (strcmp(name, "plate_with_hole") == 0) out = mio_plate_with_hole(&opts);
    else if (strcmp(name, "honeycomb") == 0) out = mio_honeycomb(&opts);
    else if (strcmp(name, "line") == 0) out = mio_line(&opts);
    else if (strcmp(name, "circle") == 0) out = mio_circle(&opts);
    else if (strcmp(name, "helix") == 0) out = mio_helix(&opts);
    else if (strcmp(name, "sphere_points") == 0) out = mio_sphere_points(&opts);
    else if (strcmp(name, "box_points") == 0) out = mio_box_points(&opts);
    else if (strcmp(name, "tetra_sphere") == 0) out = mio_tetra_sphere(&opts);
    else if (strcmp(name, "octa_sphere") == 0) out = mio_octa_sphere(&opts);
    else if (strcmp(name, "cube_sphere") == 0) out = mio_cube_sphere(&opts);
    else if (strcmp(name, "geo_sphere") == 0) out = mio_geo_sphere(&opts);
    else if (strcmp(name, "torus") == 0) out = mio_torus(&opts);
    else if (strcmp(name, "moebius") == 0) out = mio_moebius(&opts);
    else if (strcmp(name, "tube") == 0) out = mio_tube(&opts);
    else if (strcmp(name, "cone") == 0) out = mio_cone(&opts);
    else if (strcmp(name, "capsule") == 0) out = mio_capsule(&opts);
    else if (strcmp(name, "box_surface") == 0) out = mio_box_surface(&opts);
    else if (strcmp(name, "platonic") == 0) out = mio_platonic(&opts);
    else if (strcmp(name, "ball") == 0) out = mio_ball(&opts);
    else if (strcmp(name, "spherical_shell") == 0) out = mio_spherical_shell(&opts);
    else if (strcmp(name, "pipe") == 0) out = mio_pipe(&opts);
    else if (strcmp(name, "fichera") == 0) out = mio_fichera(&opts);
    else if (strcmp(name, "solid_torus") == 0) out = mio_solid_torus(&opts);
    else Rf_error("unknown shape '%s'", name);
    if (out == NULL) mio_r_fail(name);
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_reference_cells(SEXP types, SEXP pairs) {
    mio_mesh *out = mio_reference_cells(mio_r_string(types, "types"), mio_r_bool(pairs, "pairs"));
    if (out == NULL) mio_r_fail("reference_cells");
    return mio_r_wrap_mesh(out);
}

SEXP R_mio_voxelize(SEXP mesh, SEXP resolution, SEXP cell_size, SEXP bounds, SEXP padding,
                    SEXP padding_relative, SEXP fill, SEXP sign, SEXP attach_occupancy,
                    SEXP max_cells, SEXP watertight_check, SEXP plane) {
    mio_voxel_opts opts;
    mio_voxel_opts_init(&opts);

    /* The buffers the option pointers reference must outlive the call. */
    int64_t res[3];
    double bnd[6];
    if (Rf_xlength(resolution) > 0) {
        SEXP r = PROTECT(Rf_coerceVector(resolution, REALSXP));
        if (Rf_xlength(r) != 3) {
            UNPROTECT(1);
            Rf_error("`resolution` must have exactly 3 elements");
        }
        for (int i = 0; i < 3; ++i) res[i] = (int64_t)REAL(r)[i];
        UNPROTECT(1);
        opts.resolution = res;
    }
    if (Rf_xlength(bounds) > 0) {
        SEXP b = PROTECT(Rf_coerceVector(bounds, REALSXP));
        if (Rf_xlength(b) != 6) {
            UNPROTECT(1);
            Rf_error("`bounds` must have exactly 6 elements");
        }
        memcpy(bnd, REAL(b), 6 * sizeof(double));
        UNPROTECT(1);
        opts.bounds = bnd;
    }
    opts.cell_size = mio_r_double(cell_size, "cell_size");
    opts.padding = mio_r_double(padding, "padding");
    opts.padding_relative = mio_r_double(padding_relative, "padding_relative");
    opts.fill = (int32_t)mio_r_int(fill, "fill");
    opts.sign = (int32_t)mio_r_int(sign, "sign");
    opts.attach_occupancy = mio_r_bool(attach_occupancy, "attach_occupancy");
    opts.max_cells = mio_r_int64(max_cells, "max_cells");
    opts.watertight_check = (int32_t)mio_r_int(watertight_check, "watertight_check");
    opts.plane = (int32_t)mio_r_int(plane, "plane");

    int64_t cdims[3] = {0, 0, 0}, occupied = 0;
    double origin[3] = {0, 0, 0}, spacing[3] = {0, 0, 0};
    mio_mesh *out = mio_voxelize(mio_r_mesh(mesh), &opts, cdims, origin, spacing, &occupied);
    if (out == NULL) mio_r_fail("voxelize");

    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP sd = PROTECT(Rf_allocVector(REALSXP, 3));
    for (int i = 0; i < 3; ++i) REAL(sd)[i] = (double)cdims[i];
    SEXP so = PROTECT(Rf_allocVector(REALSXP, 3));
    memcpy(REAL(so), origin, 3 * sizeof(double));
    SEXP ss = PROTECT(Rf_allocVector(REALSXP, 3));
    memcpy(REAL(ss), spacing, 3 * sizeof(double));
    const char *names[] = {"mesh", "dims", "origin", "spacing", "num_occupied"};
    const SEXP values[] = {mo, sd, so, ss, NULL};
    const double numbers[] = {(double)occupied};
    SEXP res_list = mio_r_report_list(5, names, values, numbers);
    UNPROTECT(4);
    return res_list;
}

SEXP R_mio_octree(SEXP mesh, SEXP root_cells, SEXP min_level, SEXP max_level, SEXP bounds,
                  SEXP padding, SEXP padding_relative, SEXP dimension, SEXP sample_distance,
                  SEXP balance, SEXP fill, SEXP sign, SEXP attach_distance, SEXP even_odd,
                  SEXP max_cells, SEXP watertight_check, SEXP transition, SEXP snap,
                  SEXP surface_level, SEXP curvature_threshold, SEXP feature_angle,
                  SEXP size_array, SEXP plane) {
    mio_octree_opts opts;
    mio_octree_opts_init(&opts);

    /* The buffers the option pointers reference must outlive the call. */
    int64_t roots[3];
    double bnd[6];
    if (Rf_xlength(root_cells) > 0) {
        SEXP r = PROTECT(Rf_coerceVector(root_cells, REALSXP));
        if (Rf_xlength(r) != 3) {
            UNPROTECT(1);
            Rf_error("`root_cells` must have exactly 3 elements");
        }
        for (int i = 0; i < 3; ++i) roots[i] = (int64_t)REAL(r)[i];
        UNPROTECT(1);
        opts.root_cells = roots;
    }
    if (Rf_xlength(bounds) > 0) {
        SEXP b = PROTECT(Rf_coerceVector(bounds, REALSXP));
        if (Rf_xlength(b) != 6) {
            UNPROTECT(1);
            Rf_error("`bounds` must have exactly 6 elements");
        }
        memcpy(bnd, REAL(b), 6 * sizeof(double));
        UNPROTECT(1);
        opts.bounds = bnd;
    }
    opts.min_level = mio_r_int64(min_level, "min_level");
    opts.max_level = mio_r_int64(max_level, "max_level");
    opts.padding = mio_r_double(padding, "padding");
    opts.padding_relative = mio_r_double(padding_relative, "padding_relative");
    opts.dimension = mio_r_int64(dimension, "dimension");
    opts.sample_distance = mio_r_double(sample_distance, "sample_distance");
    opts.balance = (int32_t)mio_r_int(balance, "balance");
    opts.fill = (int32_t)mio_r_int(fill, "fill");
    opts.sign = (int32_t)mio_r_int(sign, "sign");
    opts.attach_distance = mio_r_bool(attach_distance, "attach_distance");
    opts.even_odd = mio_r_bool(even_odd, "even_odd");
    opts.max_cells = mio_r_int64(max_cells, "max_cells");
    opts.watertight_check = (int32_t)mio_r_int(watertight_check, "watertight_check");
    opts.transition = (int32_t)mio_r_int(transition, "transition");
    opts.snap = mio_r_bool(snap, "snap");
    opts.surface_level = mio_r_int64(surface_level, "surface_level");
    opts.curvature_threshold = mio_r_double(curvature_threshold, "curvature_threshold");
    opts.feature_angle = mio_r_double(feature_angle, "feature_angle");
    opts.plane = (int32_t)mio_r_int(plane, "plane");
    if (TYPEOF(size_array) == STRSXP && Rf_xlength(size_array) == 1 &&
        CHAR(STRING_ELT(size_array, 0))[0] != '\0')
        opts.size_array = CHAR(STRING_ELT(size_array, 0));

    mio_octree_report report;
    mio_mesh *out = mio_octree(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("octree");

    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP so = PROTECT(Rf_allocVector(REALSXP, 3));
    memcpy(REAL(so), report.origin, 3 * sizeof(double));
    const char *names[] = {"mesh", "dimension", "num_leaves", "num_hanging",
                           "finest_level", "origin", "root_size", "num_snapped",
                           "num_refused"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, so, NULL, NULL, NULL};
    const double numbers[] = {(double)report.dimension,    (double)report.num_leaves,
                              (double)report.num_hanging,  (double)report.finest_level,
                              report.root_size,            (double)report.num_snapped,
                              (double)report.num_refused};
    SEXP res_list = mio_r_report_list(9, names, values, numbers);
    UNPROTECT(2);
    return res_list;
}


/* Every option is NULL for "keep mio_octree_mesh_opts_init's default", so the
 * defaults come from the C side; the R wrapper maps the enum strings. */
SEXP R_mio_octree_mesh(SEXP mesh, SEXP output, SEXP conformity, SEXP skin_grading, SEXP fill,
                       SEXP bounds, SEXP size_array, SEXP surface_region, SEXP min_size,
                       SEXP max_size, SEXP padding, SEXP padding_relative,
                       SEXP curvature_angle, SEXP feature_angle, SEXP corner_angle,
                       SEXP feature_size, SEXP gap_cells, SEXP surface_size,
                       SEXP winding_threshold, SEXP winding_beta, SEXP alpha_long,
                       SEXP alpha_short, SEXP max_level, SEXP grading_layers,
                       SEXP fit_iterations, SEXP untangle_iterations, SEXP max_cells,
                       SEXP max_tets, SEXP gap_outside, SEXP preserve_features,
                       SEXP dihedral_floor_guard, SEXP keep_outside, SEXP attach_distance,
                       SEXP attach_quality, SEXP pillow, SEXP hex_grading, SEXP volume_fix,
                       SEXP zero_volume_factor, SEXP cap_flatness, SEXP crossing_margin,
                       SEXP crossing_gap) {
    mio_octree_mesh_opts opts;
    mio_octree_mesh_opts_init(&opts);

#define R_OCTREE_MESH_ENUM(f) \
    if (!Rf_isNull(f)) opts.f = (int32_t)mio_r_int(f, #f)
#define R_OCTREE_MESH_DOUBLE(f) \
    if (!Rf_isNull(f)) opts.f = mio_r_double(f, #f)
#define R_OCTREE_MESH_INT64(f) \
    if (!Rf_isNull(f)) opts.f = mio_r_int64(f, #f)
#define R_OCTREE_MESH_BOOL(f) \
    if (!Rf_isNull(f)) opts.f = (int32_t)mio_r_bool(f, #f)
    R_OCTREE_MESH_ENUM(output);
    R_OCTREE_MESH_ENUM(conformity);
    R_OCTREE_MESH_ENUM(skin_grading);
    R_OCTREE_MESH_ENUM(fill);
    R_OCTREE_MESH_ENUM(hex_grading);
    R_OCTREE_MESH_DOUBLE(min_size);
    R_OCTREE_MESH_DOUBLE(max_size);
    R_OCTREE_MESH_DOUBLE(padding);
    R_OCTREE_MESH_DOUBLE(padding_relative);
    R_OCTREE_MESH_DOUBLE(curvature_angle);
    R_OCTREE_MESH_DOUBLE(feature_angle);
    R_OCTREE_MESH_DOUBLE(corner_angle);
    R_OCTREE_MESH_DOUBLE(feature_size);
    R_OCTREE_MESH_DOUBLE(gap_cells);
    R_OCTREE_MESH_DOUBLE(surface_size);
    R_OCTREE_MESH_DOUBLE(winding_threshold);
    R_OCTREE_MESH_DOUBLE(winding_beta);
    R_OCTREE_MESH_DOUBLE(alpha_long);
    R_OCTREE_MESH_DOUBLE(alpha_short);
    R_OCTREE_MESH_INT64(max_level);
    R_OCTREE_MESH_INT64(grading_layers);
    R_OCTREE_MESH_INT64(fit_iterations);
    R_OCTREE_MESH_INT64(untangle_iterations);
    R_OCTREE_MESH_INT64(max_cells);
    R_OCTREE_MESH_INT64(max_tets);
    R_OCTREE_MESH_BOOL(gap_outside);
    R_OCTREE_MESH_BOOL(preserve_features);
    R_OCTREE_MESH_BOOL(dihedral_floor_guard);
    R_OCTREE_MESH_BOOL(keep_outside);
    R_OCTREE_MESH_BOOL(attach_distance);
    R_OCTREE_MESH_BOOL(attach_quality);
    R_OCTREE_MESH_BOOL(pillow);
    R_OCTREE_MESH_BOOL(volume_fix);
    R_OCTREE_MESH_DOUBLE(zero_volume_factor);
    R_OCTREE_MESH_DOUBLE(cap_flatness);
    R_OCTREE_MESH_DOUBLE(crossing_margin);
    R_OCTREE_MESH_DOUBLE(crossing_gap);
#undef R_OCTREE_MESH_ENUM
#undef R_OCTREE_MESH_DOUBLE
#undef R_OCTREE_MESH_INT64
#undef R_OCTREE_MESH_BOOL

    /* The buffers the option pointers reference must outlive the call; the
     * names are owned by the caller's character vectors. */
    double bnd[6];
    if (Rf_xlength(bounds) > 0) {
        SEXP b = PROTECT(Rf_coerceVector(bounds, REALSXP));
        if (Rf_xlength(b) != 6) {
            UNPROTECT(1);
            Rf_error("`bounds` must have exactly 6 elements");
        }
        memcpy(bnd, REAL(b), 6 * sizeof(double));
        UNPROTECT(1);
        opts.bounds = bnd;
    }
    const char *size_name = mio_r_opt_string(size_array);
    if (size_name != NULL && size_name[0] != '\0') opts.size_array = size_name;
    const char *region_name = mio_r_opt_string(surface_region);
    if (region_name != NULL && region_name[0] != '\0') opts.surface_region = region_name;

    mio_octree_mesh_report report;
    mio_mesh *out = mio_octree_mesh(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("octree_mesh");

    /* Wrap the mesh first: once it has a finalizer nothing below can leak it. */
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP warnings = PROTECT(mio_r_getbuffer(mio_octree_mesh_warnings, "octree_mesh"));
    SEXP stages = PROTECT(mio_r_getbuffer(mio_octree_mesh_stage_names, "octree_mesh"));
    static const char *outputs[] = {"castellated", "tet", "hex"};
    const char *output_name =
        (report.output >= 0 && report.output <= 2) ? outputs[report.output] : "unknown";

    enum { R_OCTREE_MESH_N = 53 };
    SEXP res = PROTECT(Rf_allocVector(VECSXP, R_OCTREE_MESH_N));
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, R_OCTREE_MESH_N));
    int k = 0;
#define R_OCTREE_MESH_SET(name, value)                \
    do {                                              \
        SET_STRING_ELT(nm, k, Rf_mkChar(name));       \
        SET_VECTOR_ELT(res, k, (value));              \
        ++k;                                          \
    } while (0)
#define R_OCTREE_MESH_NUM(f) R_OCTREE_MESH_SET(#f, Rf_ScalarReal((double)report.f))
    R_OCTREE_MESH_SET("mesh", mo);
    R_OCTREE_MESH_SET("output", Rf_mkString(output_name));
    R_OCTREE_MESH_SET("guarantee_applies", Rf_ScalarLogical(report.guarantee_applies != 0));
    {
        SEXP so = Rf_allocVector(REALSXP, 3);
        SET_VECTOR_ELT(res, k, so); /* protected by `res` from here on */
        memcpy(REAL(so), report.origin, 3 * sizeof(double));
        SET_STRING_ELT(nm, k, Rf_mkChar("origin"));
        ++k;
    }
    R_OCTREE_MESH_NUM(root_size);
    {
        SEXP sr = Rf_allocVector(REALSXP, 3);
        SET_VECTOR_ELT(res, k, sr);
        for (int i = 0; i < 3; ++i) REAL(sr)[i] = (double)report.root_cells[i];
        SET_STRING_ELT(nm, k, Rf_mkChar("root_cells"));
        ++k;
    }
    R_OCTREE_MESH_NUM(max_level);
    R_OCTREE_MESH_NUM(finest_level);
    R_OCTREE_MESH_NUM(surface_level);
    R_OCTREE_MESH_NUM(min_size);
    R_OCTREE_MESH_NUM(max_size);
    R_OCTREE_MESH_NUM(num_leaves);
    R_OCTREE_MESH_NUM(num_points);
    R_OCTREE_MESH_NUM(num_cells);
    R_OCTREE_MESH_NUM(num_hexahedra);
    R_OCTREE_MESH_NUM(num_tetrahedra);
    R_OCTREE_MESH_NUM(num_inside);
    R_OCTREE_MESH_NUM(num_boundary_cells);
    R_OCTREE_MESH_NUM(num_surrogate_faces);
    R_OCTREE_MESH_NUM(num_hanging);
    R_OCTREE_MESH_NUM(num_cut_points);
    R_OCTREE_MESH_NUM(num_warped);
    R_OCTREE_MESH_NUM(num_warp_refused);
    R_OCTREE_MESH_NUM(min_dihedral);
    R_OCTREE_MESH_NUM(max_dihedral);
    R_OCTREE_MESH_NUM(min_scaled_jacobian);
    R_OCTREE_MESH_NUM(hausdorff_to_skin);
    R_OCTREE_MESH_NUM(num_inverted);
    R_OCTREE_MESH_NUM(num_boundary_edges_open);
    R_OCTREE_MESH_NUM(num_non_manifold_boundary_edges);
    R_OCTREE_MESH_NUM(num_corners);
    R_OCTREE_MESH_NUM(num_corners_snapped);
    R_OCTREE_MESH_NUM(num_feature_chains);
    R_OCTREE_MESH_NUM(num_feature_chains_missed);
    R_OCTREE_MESH_NUM(num_feature_moves);
    R_OCTREE_MESH_NUM(num_optimized_moves);
    R_OCTREE_MESH_NUM(num_flips);
    R_OCTREE_MESH_NUM(num_winding_rechecks);
    R_OCTREE_MESH_NUM(num_winding_skipped);
    R_OCTREE_MESH_NUM(num_gaps_filled);
    R_OCTREE_MESH_NUM(num_gaps_skipped);
    R_OCTREE_MESH_NUM(num_collapsed);
    R_OCTREE_MESH_NUM(num_reordered);
    R_OCTREE_MESH_NUM(num_overlaps);
    R_OCTREE_MESH_NUM(num_caps_removed);
    R_OCTREE_MESH_NUM(num_split);
    R_OCTREE_MESH_NUM(num_slivers);
    R_OCTREE_MESH_NUM(num_untangled);
    R_OCTREE_MESH_NUM(num_warnings);
    /* One string, newline-separated; the R wrapper splits it. */
    R_OCTREE_MESH_SET("warnings", warnings);
    /* The stage names, newline-separated (the R wrapper splits them), their
     * seconds, and the leaves at each level. */
    R_OCTREE_MESH_SET("stage_names", stages);
    {
        const int64_t count = mio_octree_mesh_stage_seconds(NULL, 0);
        if (count < 0) mio_r_fail("octree_mesh");
        SEXP ss = Rf_allocVector(REALSXP, (R_xlen_t)count);
        SET_VECTOR_ELT(res, k, ss);
        SET_STRING_ELT(nm, k, Rf_mkChar("stage_seconds"));
        ++k;
        if (count > 0 && mio_octree_mesh_stage_seconds(REAL(ss), count) < 0)
            mio_r_fail("octree_mesh");
    }
    {
        const int64_t count = mio_octree_mesh_leaves_per_level(NULL, 0);
        if (count < 0) mio_r_fail("octree_mesh");
        SEXP lv = Rf_allocVector(REALSXP, (R_xlen_t)count);
        SET_VECTOR_ELT(res, k, lv);
        SET_STRING_ELT(nm, k, Rf_mkChar("leaves_per_level"));
        ++k;
        if (count > 0) {
            int64_t *tmp = (int64_t *)R_alloc((size_t)count, sizeof(int64_t));
            if (mio_octree_mesh_leaves_per_level(tmp, count) < 0) mio_r_fail("octree_mesh");
            for (int64_t i = 0; i < count; ++i) REAL(lv)[i] = (double)tmp[i];
        }
    }
#undef R_OCTREE_MESH_NUM
#undef R_OCTREE_MESH_SET
    if (k != R_OCTREE_MESH_N) Rf_error("octree_mesh: report has %d fields, expected %d", k,
                                       (int)R_OCTREE_MESH_N);
    Rf_setAttrib(res, R_NamesSymbol, nm);
    UNPROTECT(5);
    return res;
}

SEXP R_mio_compute_sdf(SEXP mesh, SEXP structure, SEXP resolution, SEXP cell_size,
                       SEXP bounds, SEXP padding, SEXP padding_relative,
                       SEXP root_resolution, SEXP max_depth, SEXP band_cells,
                       SEXP record_levels, SEXP max_cells, SEXP sign, SEXP location,
                       SEXP band, SEXP watertight_check) {
    mio_compute_sdf_opts opts;
    mio_compute_sdf_opts_init(&opts);

    /* The buffers the option pointers reference must outlive the call. */
    int64_t res[3];
    double bnd[6];
    if (Rf_xlength(resolution) > 0) {
        SEXP r = PROTECT(Rf_coerceVector(resolution, REALSXP));
        if (Rf_xlength(r) != 3) {
            UNPROTECT(1);
            Rf_error("`resolution` must have exactly 3 elements");
        }
        for (int i = 0; i < 3; ++i) res[i] = (int64_t)REAL(r)[i];
        UNPROTECT(1);
        opts.resolution = res;
    }
    if (Rf_xlength(bounds) > 0) {
        SEXP b = PROTECT(Rf_coerceVector(bounds, REALSXP));
        if (Rf_xlength(b) != 6) {
            UNPROTECT(1);
            Rf_error("`bounds` must have exactly 6 elements");
        }
        memcpy(bnd, REAL(b), 6 * sizeof(double));
        UNPROTECT(1);
        opts.bounds = bnd;
    }
    opts.structure = (int32_t)mio_r_int(structure, "structure");
    opts.cell_size = mio_r_double(cell_size, "cell_size");
    opts.padding = mio_r_double(padding, "padding");
    opts.padding_relative = mio_r_double(padding_relative, "padding_relative");
    opts.root_resolution = mio_r_int64(root_resolution, "root_resolution");
    opts.max_depth = mio_r_int64(max_depth, "max_depth");
    opts.band_cells = mio_r_double(band_cells, "band_cells");
    opts.record_levels = mio_r_bool(record_levels, "record_levels");
    opts.max_cells = mio_r_int64(max_cells, "max_cells");
    fill_sdf_opts(&opts.distance, sign, location, band, Rf_ScalarLogical(0),
                  watertight_check);

    int64_t cdims[3] = {0, 0, 0}, depth = 0, banded = 0;
    double origin[3] = {0, 0, 0}, spacing[3] = {0, 0, 0};
    mio_surface_quality q;
    memset(&q, 0, sizeof(q));
    mio_mesh *out = mio_compute_sdf(mio_r_mesh(mesh), &opts, cdims, origin, spacing, &depth,
                                    &banded, &q);
    if (out == NULL) mio_r_fail("compute_sdf");

    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP sd = PROTECT(Rf_allocVector(REALSXP, 3));
    for (int i = 0; i < 3; ++i) REAL(sd)[i] = (double)cdims[i];
    SEXP so = PROTECT(Rf_allocVector(REALSXP, 3));
    memcpy(REAL(so), origin, 3 * sizeof(double));
    SEXP ss = PROTECT(Rf_allocVector(REALSXP, 3));
    memcpy(REAL(ss), spacing, 3 * sizeof(double));
    SEXP ql = PROTECT(quality_list(&q));
    const char *names[] = {"mesh", "dims", "origin", "spacing", "max_depth", "num_banded",
                           "quality"};
    const SEXP values[] = {mo, sd, so, ss, NULL, NULL, ql};
    const double numbers[] = {(double)depth, (double)banded};
    SEXP res_list = mio_r_report_list(7, names, values, numbers);
    UNPROTECT(5);
    return res_list;
}

/* Retetrahedralize a volume mesh (or closed surface) by isosurface stuffing
   -- remesh's volumetric sibling. Same "R links the real C header, no
   hand-mirrored struct" note as R_mio_remesh above; `distance` is set
   directly rather than through fill_sdf_opts, whose location/band/
   record_inside fields remesh_volume never reads. */
SEXP R_mio_remesh_volume(SEXP mesh, SEXP resolution, SEXP cell_size, SEXP bounds,
                         SEXP padding, SEXP padding_relative, SEXP max_cells, SEXP max_tets,
                         SEXP warp_fraction, SEXP sign, SEXP watertight_check,
                         SEXP surface_region, SEXP grid_cell_size, SEXP max_winding_work) {
    mio_remesh_volume_opts opts;
    mio_remesh_volume_report report;
    mio_mesh *out;

    mio_remesh_volume_opts_init(&opts);

    /* The buffers the option pointers reference must outlive the call, the
       R_mio_compute_sdf idiom above. */
    int64_t res[3];
    double bnd[6];
    if (Rf_xlength(resolution) > 0) {
        SEXP r = PROTECT(Rf_coerceVector(resolution, REALSXP));
        if (Rf_xlength(r) != 3) {
            UNPROTECT(1);
            Rf_error("`resolution` must have exactly 3 elements");
        }
        for (int i = 0; i < 3; ++i) res[i] = (int64_t)REAL(r)[i];
        UNPROTECT(1);
        opts.resolution = res;
    }
    if (Rf_xlength(bounds) > 0) {
        SEXP b = PROTECT(Rf_coerceVector(bounds, REALSXP));
        if (Rf_xlength(b) != 6) {
            UNPROTECT(1);
            Rf_error("`bounds` must have exactly 6 elements");
        }
        memcpy(bnd, REAL(b), 6 * sizeof(double));
        UNPROTECT(1);
        opts.bounds = bnd;
    }
    opts.cell_size = mio_r_double(cell_size, "cell_size");
    opts.padding = mio_r_double(padding, "padding");
    opts.padding_relative = mio_r_double(padding_relative, "padding_relative");
    opts.max_cells = mio_r_int64(max_cells, "max_cells");
    opts.max_tets = mio_r_int64(max_tets, "max_tets");
    opts.warp_fraction = mio_r_double(warp_fraction, "warp_fraction");
    opts.distance.sign = (int32_t)mio_r_int(sign, "sign");
    opts.distance.watertight_check = (int32_t)mio_r_int(watertight_check, "watertight_check");
    opts.distance.surface_region = mio_r_opt_string(surface_region);
    opts.distance.grid_cell_size = mio_r_double(grid_cell_size, "grid_cell_size");
    opts.distance.max_winding_work = mio_r_double(max_winding_work, "max_winding_work");

    out = mio_remesh_volume_ex(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("remesh_volume");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP ql = PROTECT(quality_list(&report.input_quality));
    const char *names[] = {"mesh", "num_tets", "num_vertices_warped", "num_tets_rejected",
                           "num_non_manifold_edges", "input_quality"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL, ql};
    const double numbers[] = {(double)report.num_tets, (double)report.num_vertices_warped,
                              (double)report.num_tets_rejected,
                              (double)report.num_non_manifold_edges};
    SEXP res_list = mio_r_report_list(6, names, values, numbers);
    UNPROTECT(2);
    return res_list;
}


SEXP R_mio_surface_watertight_check(SEXP mesh) {
    mio_surface_quality q;
    memset(&q, 0, sizeof(q));
    mio_r_check(mio_surface_watertight_check(mio_r_mesh(mesh), &q), "surface_watertight_check");
    return quality_list(&q);
}

SEXP R_mio_sample_distance(SEXP mesh, SEXP points, SEXP sign, SEXP band,
                           SEXP watertight_check) {
    if (!Rf_isMatrix(points)) Rf_error("`points` must be a (3 x n) matrix");
    SEXP p = PROTECT(Rf_coerceVector(points, REALSXP));
    SEXP dim = Rf_getAttrib(points, R_DimSymbol);
    if (INTEGER(dim)[0] != 3) {
        UNPROTECT(1);
        Rf_error("`points` must be a (3 x n) matrix");
    }
    int64_t n = INTEGER(dim)[1];

    mio_sdf_opts opts;
    fill_sdf_opts(&opts, sign, Rf_ScalarInteger(0), band, Rf_ScalarLogical(0),
                  watertight_check);

    SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)n));
    /* Column-major R == row-major C, so the buffer goes across as it stands. */
    mio_status st = mio_sample_distance(mio_r_mesh(mesh), REAL(p), n, &opts, REAL(out));
    if (st != MIO_OK) {
        UNPROTECT(2);
        mio_r_fail("sample_distance");
    }
    UNPROTECT(2);
    return out;
}

SEXP R_mio_distance_to_surface(SEXP mesh, SEXP surface, SEXP sign, SEXP location, SEXP band,
                               SEXP record_inside, SEXP watertight_check) {
    mio_sdf_opts opts;
    fill_sdf_opts(&opts, sign, location, band, record_inside, watertight_check);

    int64_t banded = 0;
    mio_surface_quality q;
    memset(&q, 0, sizeof(q));
    mio_mesh *out =
        mio_distance_to_surface(mio_r_mesh(mesh), mio_r_mesh(surface), &opts, &banded, &q);
    if (out == NULL) mio_r_fail("distance_to_surface");

    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP ql = PROTECT(quality_list(&q));
    const char *names[] = {"mesh", "num_banded", "quality"};
    const SEXP values[] = {mo, NULL, ql};
    const double numbers[] = {(double)banded};
    SEXP res = mio_r_report_list(3, names, values, numbers);
    UNPROTECT(2);
    return res;
}

/* --- Provenance (see doc/provenance.md) --------------------------------
 *
 * Rides the flat C ABI. R has no RAII, so the begin/end pair is exposed as
 * two entry points and the R-side `mio_with_provenance()` wrapper pairs them
 * with `on.exit()` -- the same shape `mio_xdmf_series()` already uses to give
 * a C handle an R-shaped lifetime. */

SEXP R_mio_provenance_begin(SEXP mode) {
    if (mio_provenance_scope_begin(Rf_asInteger(mode)) != MIO_OK)
        mio_r_fail("provenance_scope_begin");
    return R_NilValue;
}

SEXP R_mio_provenance_end(void) {
    mio_provenance_scope_end();
    return R_NilValue;
}

SEXP R_mio_provenance_note(SEXP category, SEXP detail) {
    mio_provenance_note(mio_r_string(category, "category"), mio_r_string(detail, "detail"));
    return R_NilValue;
}

SEXP R_mio_provenance_set_source(SEXP path, SEXP format) {
    mio_provenance_set_source(mio_r_string(path, "path"), mio_r_string(format, "format"));
    return R_NilValue;
}

SEXP R_mio_provenance_set_target(SEXP format, SEXP encoding, SEXP codec, SEXP float_format) {
    mio_provenance_set_target(mio_r_string(format, "format"), mio_r_string(encoding, "encoding"),
                              mio_r_string(codec, "codec"),
                              mio_r_string(float_format, "float_format"));
    return R_NilValue;
}

/* A numeric vector of exactly `n` values into `dst`; returns 0 when `s` is
 * NULL or empty (the criterion is off). */
static int select_vec(SEXP s, double *dst, R_xlen_t n, const char *what) {
    R_xlen_t i;
    SEXP r;
    if (Rf_isNull(s) || Rf_xlength(s) == 0) return 0;
    if (!Rf_isNumeric(s) || Rf_xlength(s) != n)
        Rf_error("meshio++: select: %s must be %d numbers", what, (int)n);
    r = PROTECT(Rf_coerceVector(s, REALSXP));
    for (i = 0; i < n; ++i) dst[i] = REAL(r)[i];
    UNPROTECT(1);
    return 1;
}

SEXP R_mio_select(SEXP mesh, SEXP name, SEXP kind, SEXP tag, SEXP mode, SEXP boundary_only,
                  SEXP box, SEXP plane, SEXP sphere, SEXP where_array, SEXP where_op,
                  SEXP where_value, SEXP quality_metric, SEXP quality_op, SEXP quality_value,
                  SEXP component, SEXP seed, SEXP feature_angle, SEXP normal_tolerance) {
    mio_select_opts opts;
    mio_select_report report;
    mio_mesh *out;
    double v[6];
    const char *k, *m;

    mio_select_opts_init(&opts);
    opts.name = mio_r_string(name, "name");
    k = mio_r_string(kind, "kind");
    if (strcmp(k, "point") == 0)
        opts.kind = MIO_REGION_POINT;
    else if (strcmp(k, "cell") == 0)
        opts.kind = MIO_REGION_CELL;
    else if (strcmp(k, "side") == 0)
        opts.kind = MIO_REGION_SIDE;
    else
        Rf_error("meshio++: select: kind must be \"point\", \"cell\" or \"side\"");
    opts.tag = (int64_t)mio_r_double(tag, "tag");
    m = mio_r_string(mode, "mode");
    if (strcmp(m, "all") == 0)
        opts.mode = MIO_SELECT_ALL;
    else if (strcmp(m, "any") == 0)
        opts.mode = MIO_SELECT_ANY;
    else if (strcmp(m, "centroid") == 0)
        opts.mode = MIO_SELECT_CENTROID;
    else
        Rf_error("meshio++: select: mode must be \"all\", \"any\" or \"centroid\"");
    opts.boundary_only = mio_r_bool(boundary_only, "boundary_only") ? 1 : 0;
    if (select_vec(box, v, 6, "box")) {
        opts.use_box = 1;
        memcpy(opts.box_lo, v, 3 * sizeof(double));
        memcpy(opts.box_hi, v + 3, 3 * sizeof(double));
    }
    if (select_vec(plane, v, 6, "plane")) {
        opts.use_plane = 1;
        memcpy(opts.plane_point, v, 3 * sizeof(double));
        memcpy(opts.plane_normal, v + 3, 3 * sizeof(double));
    }
    if (select_vec(sphere, v, 4, "sphere")) {
        opts.use_sphere = 1;
        memcpy(opts.sphere_center, v, 3 * sizeof(double));
        opts.sphere_radius = v[3];
    }
    opts.where_array = mio_r_opt_string(where_array);
    opts.where_compare = refine_compare_code(mio_r_opt_string(where_op));
    if (opts.where_compare < 0) Rf_error("meshio++: select: unknown comparison in where_op");
    opts.where_value = mio_r_double(where_value, "where_value");
    opts.quality_metric = mio_r_opt_string(quality_metric);
    opts.quality_compare = refine_compare_code(mio_r_opt_string(quality_op));
    if (opts.quality_compare < 0) Rf_error("meshio++: select: unknown comparison in quality_op");
    opts.quality_value = mio_r_double(quality_value, "quality_value");
    /* 1-based in R, like every other index here. */
    if (!Rf_isNull(component) && Rf_xlength(component) > 0) {
        opts.component = (int64_t)mio_r_double(component, "component") - 1;
        if (opts.component < 0) Rf_error("meshio++: select: component must be >= 1");
    }
    if (select_vec(seed, v, 3, "seed")) {
        opts.use_seed = 1;
        memcpy(opts.seed, v, 3 * sizeof(double));
    }
    opts.feature_angle = mio_r_double(feature_angle, "feature_angle");
    opts.normal_tolerance = mio_r_double(normal_tolerance, "normal_tolerance");

    out = mio_select(mio_r_mesh(mesh), &opts, &report);
    if (out == NULL) mio_r_fail("select");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const char *names[] = {"mesh", "num_candidates", "num_selected", "num_grown",
                           "num_components"};
    const SEXP values[] = {mo, NULL, NULL, NULL, NULL};
    const double numbers[] = {(double)report.num_candidates, (double)report.num_selected,
                              (double)report.num_grown, (double)report.num_components};
    SEXP res = mio_r_report_list(5, names, values, numbers);
    UNPROTECT(1);
    return res;
}

SEXP R_mio_deviation(SEXP scan, SEXP reference, SEXP ptype, SEXP params, SEXP tolerance,
                     SEXP bins, SEXP range_min, SEXP range_max, SEXP max_normal_angle,
                     SEXP reference_region, SEXP register_first) {
    mio_deviation_opts opts;
    mio_deviation_opts_init(&opts);
    opts.tolerance = mio_r_double(tolerance, "tolerance");
    opts.bins = (int64_t)mio_r_double(bins, "bins");
    opts.range_min = mio_r_double(range_min, "range_min");
    opts.range_max = mio_r_double(range_max, "range_max");
    opts.max_normal_angle = mio_r_double(max_normal_angle, "max_normal_angle");
    opts.reference_region = mio_r_opt_string(reference_region);
    opts.register_first = mio_r_bool(register_first, "register");
    mio_deviation_report r;
    mio_mesh *out;
    if (reference == R_NilValue) {
        if (!Rf_isReal(params) || XLENGTH(params) != 10)
            Rf_error("`params` must be ten numbers");
        out = mio_deviation_primitive(mio_r_mesh(scan), (int32_t)mio_r_double(ptype, "type"),
                                      REAL(params), &opts, &r);
    } else {
        out = mio_deviation(mio_r_mesh(scan), mio_r_mesh(reference), &opts, &r);
    }
    if (out == NULL) mio_r_fail("deviation");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP pct = PROTECT(Rf_allocVector(REALSXP, 5));
    for (int k = 0; k < 5; ++k)
        REAL(pct)[k] = r.percentiles[k];
    const int64_t nb = r.num_bins;
    SEXP counts = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)nb));
    SEXP edges = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)(nb + 1)));
    int64_t *tmp = (int64_t *)R_alloc((size_t)(nb > 0 ? nb : 1), sizeof(int64_t));
    mio_deviation_histogram(tmp, REAL(edges), nb);
    for (int64_t k = 0; k < nb; ++k)
        REAL(counts)[k] = (double)tmp[k];
    SEXP warnings = PROTECT(mio_r_getbuffer(mio_deviation_warnings, "deviation"));
    SEXP reg = PROTECT(Rf_ScalarLogical(r.registered != 0));
    SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, 4, 4));
    double m[16];
    mio_deviation_matrix(m);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) REAL(mat)[j * 4 + i] = m[4 * i + j];
    const char *names[] = {"mesh",          "num_points",        "num_measured",
                           "min",           "max",               "mean",
                           "rms",           "sigma",             "percentiles",
                           "histogram_edges", "histogram_counts", "num_outside_range",
                           "share_within",  "suggested_range",   "num_without_normal",
                           "warnings",      "registered",        "matrix",
                           "registration_rms"};
    const SEXP values[] = {mo,   NULL, NULL,  NULL,   NULL, NULL, NULL, NULL,
                           pct,  edges, counts, NULL, NULL, NULL, NULL, warnings,
                           reg,  mat,  NULL};
    const double numbers[] = {(double)r.num_points,
                              (double)r.num_measured,
                              r.min,
                              r.max,
                              r.mean,
                              r.rms,
                              r.sigma,
                              (double)r.num_outside_range,
                              r.share_within,
                              r.suggested_range,
                              (double)r.num_without_normal,
                              r.registration_rms};
    SEXP res = mio_r_report_list(19, names, values, numbers);
    UNPROTECT(7);
    return res;
}

SEXP R_mio_deviation_regions(SEXP mesh, SEXP tolerance, SEXP field, SEXP prefix,
                             SEXP min_area) {
    mio_deviation_regions_opts opts;
    mio_deviation_regions_opts_init(&opts);
    opts.tolerance = mio_r_double(tolerance, "tolerance");
    opts.min_area = mio_r_double(min_area, "min_area");
    opts.field = mio_r_opt_string(field);
    opts.prefix = mio_r_opt_string(prefix);
    mio_deviation_regions_report r;
    mio_mesh *out = mio_deviation_regions(mio_r_mesh(mesh), &opts, &r);
    if (out == NULL) mio_r_fail("deviation_regions");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const R_xlen_t n = (R_xlen_t)r.num_regions;
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n));
    SEXP sign = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP ncells = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP area = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP extreme = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP mean = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP volume = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP centroid = PROTECT(Rf_allocMatrix(REALSXP, (int)n, 3));
    SEXP bounds = PROTECT(Rf_allocMatrix(REALSXP, (int)n, 6));
    for (R_xlen_t k = 0; k < n; ++k) {
        mio_deviation_region_info g;
        if (mio_deviation_region_get((int64_t)k, &g) != MIO_OK) mio_r_fail("deviation_regions");
        char buf[512];
        mio_deviation_region_name((int64_t)k, buf, (int64_t)sizeof(buf));
        SET_STRING_ELT(names, k, Rf_mkChar(buf));
        REAL(sign)[k] = (double)g.sign;
        REAL(ncells)[k] = (double)g.num_cells;
        REAL(area)[k] = g.area;
        REAL(extreme)[k] = g.extreme;
        REAL(mean)[k] = g.mean;
        REAL(volume)[k] = g.volume;
        for (int c = 0; c < 3; ++c)
            REAL(centroid)[k + (R_xlen_t)c * n] = g.centroid[c];
        for (int c = 0; c < 6; ++c)
            REAL(bounds)[k + (R_xlen_t)c * n] = g.bounds[c];
    }
    SEXP warnings = PROTECT(mio_r_getbuffer(mio_deviation_regions_warnings, "deviation_regions"));
    const char *keys[] = {"mesh",    "name",   "sign",     "num_cells",  "area",
                          "extreme", "mean",   "volume",   "centroid",   "bounds",
                          "num_excess", "num_missing", "num_dropped", "warnings"};
    const SEXP values[] = {mo,   names, sign,     ncells, area, extreme, mean,
                           volume, centroid, bounds, NULL, NULL, NULL, warnings};
    const double numbers[] = {(double)r.num_excess, (double)r.num_missing,
                              (double)r.num_dropped};
    SEXP res = mio_r_report_list(14, keys, values, numbers);
    UNPROTECT(11);
    return res;
}

SEXP R_mio_exaggerate(SEXP mesh, SEXP field, SEXP scale, SEXP clamp, SEXP smooth_fraction) {
    mio_exaggerate_opts opts;
    mio_exaggerate_opts_init(&opts);
    opts.field = mio_r_opt_string(field);
    opts.scale = mio_r_double(scale, "scale");
    opts.clamp = mio_r_double(clamp, "clamp");
    opts.smooth_fraction = mio_r_double(smooth_fraction, "smooth_fraction");
    mio_exaggerate_report r;
    mio_mesh *out = mio_exaggerate(mio_r_mesh(mesh), &opts, &r);
    if (out == NULL) mio_r_fail("exaggerate");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP warnings = PROTECT(mio_r_getbuffer(mio_exaggerate_warnings, "exaggerate"));
    const char *names[] = {"mesh", "clamp", "radius", "num_moved", "warnings"};
    const SEXP values[] = {mo, NULL, NULL, NULL, warnings};
    const double numbers[] = {r.clamp, r.radius, (double)r.num_moved};
    SEXP res = mio_r_report_list(5, names, values, numbers);
    UNPROTECT(2);
    return res;
}

SEXP R_mio_colorize(SEXP mesh, SEXP field, SEXP cmap, SEXP vmin, SEXP vmax, SEXP symmetric,
                    SEXP steps, SEXP nan_color, SEXP component) {
    mio_colorize_opts opts;
    mio_colorize_opts_init(&opts);
    opts.field = mio_r_string(field, "field");
    opts.cmap = mio_r_opt_string(cmap);
    opts.vmin = mio_r_double(vmin, "vmin");
    opts.vmax = mio_r_double(vmax, "vmax");
    opts.symmetric = Rf_asLogical(symmetric) == TRUE ? 1 : 0;
    opts.steps = (int64_t)mio_r_double(steps, "steps");
    opts.nan_color = mio_r_opt_string(nan_color);
    opts.component = (int64_t)mio_r_double(component, "component");
    mio_colorize_report r;
    mio_mesh *out = mio_colorize(mio_r_mesh(mesh), &opts, &r);
    if (out == NULL) mio_r_fail("colorize");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP loc = PROTECT(Rf_mkString(r.location == 1 ? "cell" : "point"));
    const char *names[] = {"mesh", "vmin", "vmax", "location"};
    const SEXP values[] = {mo, NULL, NULL, loc};
    const double numbers[] = {r.vmin, r.vmax};
    SEXP res = mio_r_report_list(4, names, values, numbers);
    UNPROTECT(2);
    return res;
}
