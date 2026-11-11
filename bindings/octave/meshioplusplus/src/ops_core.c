/* The operations layer, first part: computations *on* a mesh, not file
 * formats -- surface extraction, transforms, cleaning and smoothing,
 * subsetting, differential operators, segmentation and CAD fitting, interface
 * search, repair and the geometry kernels through mio_parameterize.
 *
 * A port of lines 1-1596 of the R binding's src/ops.c; see mio_oct.h for the
 * layout rules. Ownership rule for the opaque C results (find_interface,
 * contact_pairs): read the maps and counters first, then TAKE the mesh out of
 * the result, then free it. Every mesh handed to Octave therefore owns its
 * handle and cannot dangle when the result goes away. This mirrors the Julia
 * and R bindings exactly.
 *
 * R has no native int64 (nor does this binding expose one for counters), so
 * every counter comes back as a double -- exact well past any plausible cell
 * count.
 */

#include "mio_oct.h"

#include <math.h>
#include <string.h>

/* --- local helpers ------------------------------------------------------- */

/* Exactly `n` numbers into `out`, raising `msg` (R's own wording) otherwise. */
static void oc_vec_exact(const mxArray *a, double *out, size_t n, const char *what,
                         const char *msg) {
    if (a == NULL || !(mxIsNumeric(a) || mxIsLogical(a)) || mxGetNumberOfElements(a) != n) {
        mio_o_error(what, "%s", msg);
    }
    mio_o_vecn(a, out, n, what);
}

/* --- surface / skin / quality -------------------------------------------- */

static mxArray *O_mio_extract_surface(const mxArray **a) {
    mio_mesh *out =
        mio_extract_surface(mio_o_mesh(a[0]), mio_o_bool(a[1], "record_parent_ids"));
    if (out == NULL) mio_o_fail("extract_surface");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_extract_skin(const mxArray **a) {
    mio_mesh *out = mio_extract_skin(mio_o_mesh(a[0]), mio_o_bool(a[1], "linearize"));
    if (out == NULL) mio_o_fail("extract_skin");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_attach_quality(const mxArray **a) {
    mio_mesh *out = mio_attach_quality(mio_o_mesh(a[0]));
    if (out == NULL) mio_o_fail("attach_quality");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_quality_counts(const mxArray **a) {
    int64_t nc = 0, ni = 0, nd = 0;
    mio_o_check(mio_quality_counts(mio_o_mesh(a[0]), &nc, &ni, &nd), "quality_counts");
    const char *names[] = {"num_cells", "num_inverted", "num_degenerate"};
    mxArray *values[] = {mio_o_scalar((double)nc), mio_o_scalar((double)ni),
                         mio_o_scalar((double)nd)};
    return mio_o_struct(3, names, values);
}

/* --- geometry-preserving transforms -------------------------------------- */

/* transform(mesh, matrix, rotate_vector_data) */
static mxArray *O_mio_transform(const mxArray **a) {
    double m[16];
    oc_vec_exact(a[1], m, 16, "transform", "`matrix` must be a 4x4 affine matrix (16 values)");
    /* The C API takes the matrix ROW-major, while an Octave matrix is
     * column-major, so this one is transposed on the way out. That is correct
     * here and only here: a 4x4 affine matrix is a mathematical object, not a
     * mesh array. */
    double rowmajor[16];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            rowmajor[i * 4 + j] = m[j * 4 + i];
        }
    }
    mio_mesh *out =
        mio_transform(mio_o_mesh(a[0]), rowmajor, mio_o_bool(a[2], "rotate_vector_data"));
    if (out == NULL) mio_o_fail("transform");
    return mio_o_wrap_mesh(out);
}

/* clean(mesh, weld, atol, remove_orphans, drop_degenerate, drop_duplicate_cells) */
static mxArray *O_mio_clean(const mxArray **a) {
    int64_t welded = 0, orphan = 0, degen = 0, dup = 0;
    mio_mesh *out = mio_clean(mio_o_mesh(a[0]), mio_o_bool(a[1], "weld"),
                              mio_o_double(a[2], "atol"), mio_o_bool(a[3], "remove_orphans"),
                              mio_o_bool(a[4], "drop_degenerate"),
                              mio_o_bool(a[5], "drop_duplicate_cells"), &welded, &orphan,
                              &degen, &dup);
    if (out == NULL) mio_o_fail("clean");
    const char *names[] = {"mesh", "points_welded", "points_removed_orphan",
                           "cells_dropped_degenerate", "cells_dropped_duplicate"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)welded),
                         mio_o_scalar((double)orphan), mio_o_scalar((double)degen),
                         mio_o_scalar((double)dup)};
    return mio_o_struct(5, names, values);
}

/* smooth(mesh, method, iterations, lambda, mu, fix_boundary, preserve_features,
 *        feature_angle, guard_inversion, frozen) */
static mxArray *O_mio_smooth(const mxArray **a) {
    mio_smooth_opts opts;
    int64_t moved = 0, skipped = 0, n_frozen = 0;
    double disp = 0.0;

    mio_smooth_opts_init(&opts);
    opts.method = mio_o_opt_string(a[1]);
    opts.iterations = (int32_t)mio_o_int(a[2], "iterations");
    opts.lambda = mio_o_double(a[3], "lambda");
    opts.mu = mio_o_double(a[4], "mu");
    opts.fix_boundary = mio_o_bool(a[5], "fix_boundary");
    opts.preserve_features = mio_o_bool(a[6], "preserve_features");
    opts.feature_angle = mio_o_double(a[7], "feature_angle");
    opts.guard_inversion = mio_o_bool(a[8], "guard_inversion");
    /* 1-based point ids, [] for none; shifted to 0-based across the C ABI, as
     * every index array in this binding. An out-of-range id raises by name in
     * the core. */
    opts.frozen = mio_o_ids_0based(a[9], &n_frozen, "frozen");
    opts.num_frozen = n_frozen;

    mio_mesh *out = mio_smooth_ex(mio_o_mesh(a[0]), &opts, &moved, &disp, &skipped);
    if (out == NULL) mio_o_fail("smooth");
    const char *names[] = {"mesh", "num_nodes_moved", "max_displacement",
                           "num_skipped_inversion"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)moved), mio_o_scalar(disp),
                         mio_o_scalar((double)skipped)};
    return mio_o_struct(4, names, values);
}

/* optimize_volume(mesh, max_iterations, relocate, flip, preserve_boundary,
 *                 min_improvement, objective, max_ring_size) */
static mxArray *O_mio_optimize_volume(const mxArray **a) {
    mio_optimize_volume_opts opts;
    mio_optimize_volume_report rep;
    mio_optimize_volume_opts_init(&opts);
    opts.max_iterations = mio_o_int(a[1], "max_iterations");
    opts.relocate = mio_o_bool(a[2], "relocate");
    opts.flip = mio_o_bool(a[3], "flip");
    opts.preserve_boundary = mio_o_bool(a[4], "preserve_boundary");
    opts.min_improvement = mio_o_double(a[5], "min_improvement");
    const char *obj = mio_o_string(a[6], "objective");
    if (strcmp(obj, "worst") == 0)
        opts.objective = MIO_OPTIMIZE_VOLUME_WORST;
    else if (strcmp(obj, "energy") == 0)
        opts.objective = MIO_OPTIMIZE_VOLUME_ENERGY;
    else
        mio_o_error("optimize_volume",
                    "meshio++: optimize_volume: objective must be \"worst\" or \"energy\"");
    opts.max_ring_size = mio_o_int(a[7], "max_ring_size");
    mio_mesh *out = mio_optimize_volume_ex(mio_o_mesh(a[0]), &opts, &rep);
    if (out == NULL) mio_o_fail("optimize_volume");
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
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)rep.num_flips),
                         mio_o_scalar((double)rep.num_23_flips),
                         mio_o_scalar((double)rep.num_32_flips),
                         mio_o_scalar((double)rep.num_vertices_moved),
                         mio_o_scalar((double)rep.num_tets),
                         mio_o_scalar(rep.min_quality_before),
                         mio_o_scalar(rep.min_quality_after),
                         mio_o_scalar((double)rep.num_edge_removals),
                         mio_o_scalar((double)rep.num_multi_face_removals),
                         mio_o_scalar((double)rep.num_worst_smoothed),
                         mio_o_scalar(rep.min_dihedral_before),
                         mio_o_scalar(rep.min_dihedral_after)};
    return mio_o_struct(13, names, values);
}

/* adapt_volume(mesh, hmin, hmax, hgrad, hausd, hausd_relative, size_array,
 * preserve_boundary, preserve_features, quality_metric, quality_min,
 * max_iterations, collapse, split, swap, smooth, boundary_sizing) */
static mxArray *O_mio_adapt_volume(const mxArray **a) {
    mio_adapt_volume_opts opts;
    mio_adapt_volume_report rep;
    mio_adapt_volume_opts_init(&opts);
    opts.hmin = mio_o_double(a[1], "hmin");
    opts.hmax = mio_o_double(a[2], "hmax");
    opts.hgrad = mio_o_double(a[3], "hgrad");
    opts.hausd = mio_o_double(a[4], "hausd");
    opts.hausd_relative = mio_o_bool(a[5], "hausd_relative");
    const char *size_name = mio_o_string(a[6], "size_array");
    if (size_name[0] != '\0') opts.size_array = size_name;
    opts.preserve_boundary = mio_o_bool(a[7], "preserve_boundary");
    opts.preserve_features = mio_o_bool(a[8], "preserve_features");
    const char *metric = mio_o_string(a[9], "quality_metric");
    if (strcmp(metric, "min_dihedral") == 0)
        opts.quality_metric = MIO_ADAPT_VOLUME_MIN_DIHEDRAL;
    else if (strcmp(metric, "mean_ratio") == 0)
        opts.quality_metric = MIO_ADAPT_VOLUME_MEAN_RATIO;
    else if (strcmp(metric, "amips") == 0)
        opts.quality_metric = MIO_ADAPT_VOLUME_AMIPS;
    else
        mio_o_error("adapt_volume", "meshio++: adapt_volume: quality_metric must be "
                                    "\"min_dihedral\", \"mean_ratio\" or \"amips\"");
    opts.quality_min = mio_o_double(a[10], "quality_min");
    opts.max_iterations = mio_o_int(a[11], "max_iterations");
    opts.collapse = mio_o_bool(a[12], "collapse");
    opts.split = mio_o_bool(a[13], "split");
    opts.swap = mio_o_bool(a[14], "swap");
    opts.smooth = mio_o_bool(a[15], "smooth");
    const char *sizing = mio_o_string(a[16], "boundary_sizing");
    if (strcmp(sizing, "input") == 0)
        opts.boundary_sizing = MIO_ADAPT_VOLUME_SIZING_INPUT;
    else if (strcmp(sizing, "geometry") == 0)
        opts.boundary_sizing = MIO_ADAPT_VOLUME_SIZING_GEOMETRY;
    else
        mio_o_error("adapt_volume",
                    "meshio++: adapt_volume: boundary_sizing must be \"input\" or \"geometry\"");
    mio_mesh *out = mio_adapt_volume(mio_o_mesh(a[0]), &opts, &rep);
    if (out == NULL) mio_o_fail("adapt_volume");
    const char *names[] = {"mesh",
                           "num_tets_before",
                           "num_tets_after",
                           "num_points_before",
                           "num_points_after",
                           "num_collapsed",
                           "num_split",
                           "num_swapped",
                           "num_smoothed",
                           "num_rejected_by_link",
                           "num_rejected_by_inversion",
                           "num_rejected_by_quality",
                           "num_rejected_by_length",
                           "num_rejected_by_class",
                           "num_rejected_by_hausdorff",
                           "num_rejected_by_normal",
                           "num_corners",
                           "num_ridge_vertices",
                           "num_chains",
                           "num_patches",
                           "iterations",
                           "min_dihedral_before",
                           "min_dihedral_after",
                           "max_dihedral_before",
                           "max_dihedral_after",
                           "mean_quality_before",
                           "mean_quality_after",
                           "max_boundary_deviation",
                           "hmin",
                           "hmax",
                           "hausd",
                           "volume_before",
                           "volume_after"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)rep.num_tets_before),
                         mio_o_scalar((double)rep.num_tets_after),
                         mio_o_scalar((double)rep.num_points_before),
                         mio_o_scalar((double)rep.num_points_after),
                         mio_o_scalar((double)rep.num_collapsed),
                         mio_o_scalar((double)rep.num_split),
                         mio_o_scalar((double)rep.num_swapped),
                         mio_o_scalar((double)rep.num_smoothed),
                         mio_o_scalar((double)rep.num_rejected_by_link),
                         mio_o_scalar((double)rep.num_rejected_by_inversion),
                         mio_o_scalar((double)rep.num_rejected_by_quality),
                         mio_o_scalar((double)rep.num_rejected_by_length),
                         mio_o_scalar((double)rep.num_rejected_by_class),
                         mio_o_scalar((double)rep.num_rejected_by_hausdorff),
                         mio_o_scalar((double)rep.num_rejected_by_normal),
                         mio_o_scalar((double)rep.num_corners),
                         mio_o_scalar((double)rep.num_ridge_vertices),
                         mio_o_scalar((double)rep.num_chains),
                         mio_o_scalar((double)rep.num_patches),
                         mio_o_scalar((double)rep.iterations),
                         mio_o_scalar(rep.min_dihedral_before),
                         mio_o_scalar(rep.min_dihedral_after),
                         mio_o_scalar(rep.max_dihedral_before),
                         mio_o_scalar(rep.max_dihedral_after),
                         mio_o_scalar(rep.mean_quality_before),
                         mio_o_scalar(rep.mean_quality_after),
                         mio_o_scalar(rep.max_boundary_deviation),
                         mio_o_scalar(rep.hmin),
                         mio_o_scalar(rep.hmax),
                         mio_o_scalar(rep.hausd),
                         mio_o_scalar(rep.volume_before),
                         mio_o_scalar(rep.volume_after)};
    return mio_o_struct(33, names, values);
}

/* --- subsetting ---------------------------------------------------------- */

/* crop_bbox(mesh, lo, hi, mode, record_ids) */
static mxArray *O_mio_crop_bbox(const mxArray **a) {
    double l[3], h[3];
    mio_o_vec3(a[1], l, "lo");
    mio_o_vec3(a[2], h, "hi");
    mio_mesh *out = mio_crop_bbox(mio_o_mesh(a[0]), l, h, mio_o_int(a[3], "mode"),
                                  mio_o_bool(a[4], "record_ids"));
    if (out == NULL) mio_o_fail("crop_bbox");
    return mio_o_wrap_mesh(out);
}

/* crop_predicate(mesh, array, compare, value, record_ids) */
static mxArray *O_mio_crop_predicate(const mxArray **a) {
    mio_mesh *out = mio_crop_predicate(mio_o_mesh(a[0]), mio_o_string(a[1], "array"),
                                       mio_o_int(a[2], "compare"), mio_o_double(a[3], "value"),
                                       mio_o_bool(a[4], "record_ids"));
    if (out == NULL) mio_o_fail("crop_predicate");
    return mio_o_wrap_mesh(out);
}

/* crop_plane(mesh, point, normal, mode, record_ids) */
static mxArray *O_mio_crop_plane(const mxArray **a) {
    double p[3], n[3];
    mio_o_vec3(a[1], p, "point");
    mio_o_vec3(a[2], n, "normal");
    mio_mesh *out = mio_crop_plane(mio_o_mesh(a[0]), p, n, mio_o_int(a[3], "mode"),
                                   mio_o_bool(a[4], "record_ids"));
    if (out == NULL) mio_o_fail("crop_plane");
    return mio_o_wrap_mesh(out);
}

/* slice(mesh, origin, normal, record_parent_ids, chain, min_length) */
static mxArray *O_mio_slice(const mxArray **a) {
    double o[3], n[3];
    mio_o_vec3(a[1], o, "origin");
    mio_o_vec3(a[2], n, "normal");
    mio_mesh *out = mio_slice_ex(mio_o_mesh(a[0]), o, n, mio_o_bool(a[3], "record_parent_ids"),
                                 mio_o_bool(a[4], "chain"), mio_o_double(a[5], "min_length"));
    if (out == NULL) mio_o_fail("slice");
    return mio_o_wrap_mesh(out);
}

/* profile_solid(sketch, features, step_path, scan, provider, deflection)
 * (roadmap 7.7.4); features in their text form, scan a mesh or []. */
static mxArray *O_mio_profile_solid(const mxArray **a) {
    const char *text = mio_o_string(a[1], "features");
    const int64_t n = mio_profile_features_from_text(text, NULL, 0);
    if (n < 0) mio_o_fail("profile_solid");
    mio_profile_feature *feats =
        (mio_profile_feature *)mxMalloc((size_t)(n > 0 ? n : 1) * sizeof(mio_profile_feature));
    if (mio_profile_features_from_text(text, feats, n) < 0) mio_o_fail("profile_solid");
    mio_profile_solid_report r;
    mio_mesh *out = mio_profile_solid(mio_o_mesh(a[0]), feats, n, mio_o_opt_string(a[2]),
                                      mio_o_opt_string(a[4]), mio_o_double(a[5], "deflection"),
                                      mio_o_is_null(a[3]) ? NULL : mio_o_mesh(a[3]), &r);
    if (out == NULL) mio_o_fail("profile_solid");
    mxArray *w = mio_o_getbuffer(mio_profile_solid_warnings, "profile_solid");
    const char *keys[] = {"mesh", "volume", "area", "num_faces", "valid", "warnings"};
    mxArray *values[] = {mio_o_wrap_mesh(out),           mio_o_scalar(r.volume),
                         mio_o_scalar(r.area),           mio_o_scalar((double)r.num_faces),
                         mio_o_logical(r.valid != 0),    w};
    return mio_o_struct(6, keys, values);
}

/* A flat double vector of (a, b, kind) triples as int64; NULL when empty. */
static const int64_t *oc_triples(const mxArray *a, int64_t *n) {
    *n = 0;
    if (mio_o_is_null(a)) return NULL;
    size_t len = 0;
    const double *x = mio_o_doubles(a, &len, "relations");
    if (len == 0) return NULL;
    if (len % 3 != 0) mio_o_error("constrain_primitives", "relations are (a, b, kind) triples");
    int64_t *out = (int64_t *)mxMalloc(len * sizeof(int64_t));
    for (size_t i = 0; i < len; ++i) out[i] = (int64_t)x[i];
    *n = (int64_t)(len / 3);
    return out;
}

/* constrain_primitives(mesh, face_tolerance, angle_tolerance, distance_tolerance,
 * min_fraction, force, forbid) (roadmap 7.7.3); force and forbid flat triples. */
static mxArray *O_mio_constrain_primitives(const mxArray **a) {
    mio_constrain_primitives_opts opts;
    mio_constrain_primitives_opts_init(&opts);
    opts.face_tolerance = mio_o_double(a[1], "face_tolerance");
    opts.angle_tolerance_deg = mio_o_double(a[2], "angle_tolerance");
    opts.distance_tolerance = mio_o_double(a[3], "distance_tolerance");
    opts.min_fraction = mio_o_double(a[4], "min_fraction");
    opts.force = oc_triples(a[5], &opts.num_force);
    opts.forbid = oc_triples(a[6], &opts.num_forbid);
    mio_constrain_primitives_report r;
    mio_mesh *out = mio_constrain_primitives(mio_o_mesh(a[0]), &opts, &r);
    if (out == NULL) mio_o_fail("constrain_primitives");
    const int64_t n = r.num_relations;
    mxArray *pa = mio_o_column(n), *pb = mio_o_column(n), *kind = mio_o_column(n);
    mxArray *ang = mio_o_column(n), *dst = mio_o_column(n), *fc = mio_o_column(n);
    for (int64_t k = 0; k < n; ++k) {
        int64_t ia = 0, ib = 0;
        int32_t kd = 0, forced = 0;
        double angle = 0.0, distance = 0.0;
        mio_o_check(mio_constrain_primitives_relation(k, &ia, &ib, &kd, &angle, &distance, &forced),
                    "constrain_primitives relation");
        mxGetPr(pa)[k] = (double)ia;
        mxGetPr(pb)[k] = (double)ib;
        mxGetPr(kind)[k] = (double)kd;
        mxGetPr(ang)[k] = angle;
        mxGetPr(dst)[k] = distance;
        mxGetPr(fc)[k] = (double)forced;
    }
    mxArray *unmet = mio_o_matrix(r.num_unmet, 3);
    for (int64_t k = 0; k < r.num_unmet; ++k) {
        int64_t ia = 0, ib = 0;
        int32_t kd = 0;
        mio_o_check(mio_constrain_primitives_unmet(k, &ia, &ib, &kd), "constrain_primitives unmet");
        mxGetPr(unmet)[k] = (double)ia;
        mxGetPr(unmet)[k + r.num_unmet] = (double)ib;
        mxGetPr(unmet)[k + 2 * r.num_unmet] = (double)kd;
    }
    const char *keys[] = {"mesh",     "num_constrained", "num_groups", "face_tolerance",
                          "patch_a",  "patch_b",         "relation",   "angle_deg",
                          "distance", "forced",          "unmet"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)r.num_constrained),
                         mio_o_scalar((double)r.num_groups),
                         mio_o_scalar(r.face_tolerance),
                         pa, pb, kind, ang, dst, fc, unmet};
    return mio_o_struct(11, keys, values);
}

/* face_support(solid, scan, tolerance, measured_share, inferred_share, hole_data)
 * (roadmap 7.7.1). */
static mxArray *O_mio_face_support(const mxArray **a) {
    mio_face_support_opts opts;
    mio_face_support_opts_init(&opts);
    opts.tolerance = mio_o_double(a[2], "tolerance");
    opts.measured_share = mio_o_double(a[3], "measured_share");
    opts.inferred_share = mio_o_double(a[4], "inferred_share");
    opts.hole_data = mio_o_opt_string(a[5]);
    mio_face_support_report r;
    mio_mesh *out = mio_face_support(mio_o_mesh(a[0]), mio_o_mesh(a[1]), &opts, &r);
    if (out == NULL) mio_o_fail("face_support");
    const char *keys[] = {"mesh",      "num_measured", "num_partial", "num_inferred",
                          "num_faces", "share_within", "tolerance"};
    mxArray *values[] = {mio_o_wrap_mesh(out),           mio_o_scalar((double)r.num_measured),
                         mio_o_scalar((double)r.num_partial), mio_o_scalar((double)r.num_inferred),
                         mio_o_scalar((double)r.num_faces),   mio_o_scalar(r.share_within),
                         mio_o_scalar(r.tolerance)};
    return mio_o_struct(7, keys, values);
}

/* fit_cad_iterate(scan, step_path, part_class, target_within, target_explained,
 * segment_tolerance, fit_tolerance, max_rounds, worst_faces, refine_factor,
 * backend) (roadmap 7.7.2); part_class and backend as their C enum values. */
static mxArray *O_mio_fit_cad_iterate(const mxArray **a) {
    mio_fit_cad_iterate_opts opts;
    mio_fit_cad_iterate_opts_init(&opts);
    opts.part_class = mio_o_int(a[2], "part_class");
    opts.target_within = mio_o_double(a[3], "target_within");
    opts.target_explained = mio_o_double(a[4], "target_explained");
    opts.segment_tolerance = mio_o_double(a[5], "segment_tolerance");
    opts.fit_tolerance = mio_o_double(a[6], "fit_tolerance");
    opts.max_rounds = mio_o_int64(a[7], "max_rounds");
    opts.worst_faces = mio_o_int64(a[8], "worst_faces");
    opts.refine_factor = mio_o_double(a[9], "refine_factor");
    opts.backend = mio_o_int(a[10], "backend");
    mio_fit_cad_iterate_report r;
    mio_mesh *out = mio_fit_cad_iterate(mio_o_mesh(a[0]), mio_o_opt_string(a[1]), &opts, &r);
    if (out == NULL) mio_o_fail("fit_cad_iterate");
    mxArray *share = mio_o_column(r.num_rounds), *expl = mio_o_column(r.num_rounds);
    for (int64_t k = 0; k < r.num_rounds; ++k)
        mio_o_check(mio_fit_cad_iterate_round(k, &mxGetPr(share)[k], &mxGetPr(expl)[k], NULL,
                                              NULL, NULL),
                    "fit_cad_iterate round");
    const char *keys[] = {"mesh",          "num_rounds",    "best_round",   "target_met",
                          "stop_reason",   "target_within", "share_within", "explained",
                          "hausdorff_max", "num_faces",     "round_share_within",
                          "round_explained"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)r.num_rounds),
                         mio_o_scalar((double)r.best_round),
                         mio_o_scalar((double)r.target_met),
                         mio_o_scalar((double)r.stop_reason),
                         mio_o_scalar(r.target_within),
                         mio_o_scalar(r.share_within),
                         mio_o_scalar(r.explained),
                         mio_o_scalar(r.hausdorff_max),
                         mio_o_scalar((double)r.num_faces),
                         share,
                         expl};
    return mio_o_struct(12, keys, values);
}

/* The struct offset_surface, thicken and hollow return (roadmap 7.8.4). */
static mxArray *mio_o_offset_struct(mio_mesh *out, const mio_offset_report *r, const char *what) {
    if (out == NULL) mio_o_fail(what);
    mxArray *w = mio_o_getbuffer(mio_offset_warnings, what);
    const char *keys[] = {"mesh",          "num_not_visible", "num_clamped", "num_folded",
                          "num_rim_quads", "max_move",        "warnings"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)r->num_not_visible),
                         mio_o_scalar((double)r->num_clamped),
                         mio_o_scalar((double)r->num_folded),
                         mio_o_scalar((double)r->num_rim_quads),
                         mio_o_scalar(r->max_move),
                         w};
    return mio_o_struct(7, keys, values);
}

/* offset_surface(mesh, distance, smooth_iterations) (roadmap 7.8.4). */
static mxArray *O_mio_offset_surface(const mxArray **a) {
    mio_offset_report r;
    mio_mesh *out = mio_offset_surface(mio_o_mesh(a[0]), mio_o_double(a[1], "distance"),
                                       mio_o_int(a[2], "smooth_iterations"), &r);
    return mio_o_offset_struct(out, &r, "offset_surface");
}

/* thicken(mesh, thickness, side, smooth_iterations); side 0 inward, 1
 * outward, 2 both (roadmap 7.8.4). */
static mxArray *O_mio_thicken(const mxArray **a) {
    mio_offset_report r;
    mio_mesh *out = mio_thicken(mio_o_mesh(a[0]), mio_o_double(a[1], "thickness"),
                                mio_o_int(a[2], "side"), mio_o_int(a[3], "smooth_iterations"), &r);
    return mio_o_offset_struct(out, &r, "thicken");
}

/* hollow(mesh, wall, smooth_iterations) (roadmap 7.8.4). */
static mxArray *O_mio_hollow(const mxArray **a) {
    mio_offset_report r;
    mio_mesh *out = mio_hollow(mio_o_mesh(a[0]), mio_o_double(a[1], "wall"),
                               mio_o_int(a[2], "smooth_iterations"), &r);
    return mio_o_offset_struct(out, &r, "hollow");
}

/* clip(mesh, origin, normal, keep_below, cap) (roadmap 7.8.3). */
static mxArray *O_mio_clip(const mxArray **a) {
    double o[3], n[3];
    mio_o_vec3(a[1], o, "origin");
    mio_o_vec3(a[2], n, "normal");
    mio_clip_report r;
    mio_mesh *out = mio_clip(mio_o_mesh(a[0]), o, n, mio_o_bool(a[3], "keep_below"),
                             mio_o_bool(a[4], "cap"), &r);
    if (out == NULL) mio_o_fail("clip");
    mxArray *w = mio_o_getbuffer(mio_clip_warnings, "clip");
    const char *keys[] = {"mesh",      "num_kept",        "num_split", "num_cap_triangles",
                          "num_loops", "num_open_chains", "warnings"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)r.num_kept),
                         mio_o_scalar((double)r.num_split),
                         mio_o_scalar((double)r.num_cap_triangles),
                         mio_o_scalar((double)r.num_loops),
                         mio_o_scalar((double)r.num_open_chains),
                         w};
    return mio_o_struct(7, keys, values);
}

/* flip_normals(mesh, region) (roadmap 7.8.2). */
static mxArray *O_mio_flip_normals(const mxArray **a) {
    int64_t n = 0;
    mio_mesh *out = mio_flip_normals(mio_o_mesh(a[0]), mio_o_opt_string(a[1]), &n);
    if (out == NULL) mio_o_fail("flip_normals");
    const char *keys[] = {"mesh", "num_flipped"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)n)};
    return mio_o_struct(2, keys, values);
}

/* remove_islands(mesh, connectivity, min_diagonal, min_area, min_cells,
 * relative, dry_run) (roadmap 7.8.1). */
static mxArray *O_mio_remove_islands(const mxArray **a) {
    mio_remove_islands_opts opts;
    mio_remove_islands_opts_init(&opts);
    opts.connectivity = mio_o_int(a[1], "connectivity");
    opts.min_diagonal = mio_o_double(a[2], "min_diagonal");
    opts.min_area = mio_o_double(a[3], "min_area");
    opts.min_cells = (int64_t)mio_o_double(a[4], "min_cells");
    opts.relative = mio_o_bool(a[5], "relative");
    opts.dry_run = mio_o_bool(a[6], "dry_run");
    mio_remove_islands_report r;
    mio_mesh *out = mio_remove_islands(mio_o_mesh(a[0]), &opts, &r);
    if (out == NULL) mio_o_fail("remove_islands");
    const char *keys[] = {"mesh", "num_components", "num_removed", "num_cells_removed",
                          "num_points_removed"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)r.num_components),
                         mio_o_scalar((double)r.num_removed),
                         mio_o_scalar((double)r.num_cells_removed),
                         mio_o_scalar((double)r.num_points_removed)};
    return mio_o_struct(5, keys, values);
}

/* measure(mesh, spec) (roadmap 7.2): the report as its JSON text. */
static mxArray *O_mio_measure(const mxArray **a) {
    mio_measure_report r;
    mio_mesh *out = mio_measure(mio_o_mesh(a[0]), mio_o_string(a[1], "spec"), &r);
    if (out == NULL) mio_o_fail("measure");
    mxArray *js = mio_o_getbuffer(mio_measure_report_json, "measure");
    mxArray *w = mio_o_getbuffer(mio_measure_warnings, "measure");
    const char *keys[] = {"mesh",     "num_elements", "num_dimensions", "num_tolerances",
                          "num_fail", "report_json",  "warnings"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)r.num_elements),
                         mio_o_scalar((double)r.num_dimensions),
                         mio_o_scalar((double)r.num_tolerances),
                         mio_o_scalar((double)r.num_fail),
                         js,
                         w};
    return mio_o_struct(7, keys, values);
}

/* align_datums(mesh, primary, secondary, origin) (roadmap 7.3.4). */
static mxArray *O_mio_align_datums(const mxArray **a) {
    double m[16];
    mio_mesh *out = mio_align_datums(mio_o_mesh(a[0]), mio_o_string(a[1], "primary"),
                                     mio_o_string(a[2], "secondary"), mio_o_opt_string(a[3]), m);
    if (out == NULL) mio_o_fail("align_datums");
    mxArray *mat = mio_o_matrix(4, 4);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) mxGetPr(mat)[i + 4 * j] = m[4 * i + j];
    const char *keys[] = {"mesh", "matrix"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mat};
    return mio_o_struct(2, keys, values);
}

/* measure_element_primitive(mesh, name) (roadmap 7.4.2). */
static mxArray *O_mio_measure_element_primitive(const mxArray **a) {
    int32_t type = 0;
    double p[10];
    mio_o_check(mio_measure_element_primitive(mio_o_mesh(a[0]), mio_o_string(a[1], "name"),
                                              &type, p),
                "measure_element_primitive");
    mxArray *params = mio_o_matrix(1, 10);
    for (int k = 0; k < 10; ++k) mxGetPr(params)[k] = p[k];
    const char *keys[] = {"type", "params"};
    mxArray *values[] = {mio_o_scalar((double)type), params};
    return mio_o_struct(2, keys, values);
}

/* write_step_elements(path, mesh, patches, segments) (roadmap 7.6.5): patches
 * a numeric vector of labels (empty: every fitted patch). */
static mxArray *O_mio_write_step_elements(const mxArray **a) {
    const mwSize n = mxIsEmpty(a[2]) ? 0 : mxGetNumberOfElements(a[2]);
    if (n > 0 && (!mxIsDouble(a[2]) || mxIsComplex(a[2])))
        mexErrMsgIdAndTxt("meshioplusplus:args", "write_step_elements: patches must be real numbers");
    int64_t *ids = (int64_t *)mxCalloc((size_t)(n > 0 ? n : 1), sizeof(int64_t));
    for (mwSize i = 0; i < n; ++i) ids[i] = (int64_t)mxGetPr(a[2])[i];
    mio_step_elements_report r;
    mio_o_check(mio_write_step_elements(mio_o_string(a[0], "path"), mio_o_mesh(a[1]),
                                        n > 0 ? ids : NULL, (int64_t)n,
                                        (int64_t)mio_o_double(a[3], "segments"), &r),
                "write_step_elements");
    mxArray *w = mio_o_getbuffer(mio_write_step_elements_warnings, "write_step_elements");
    const char *keys[] = {"num_planes", "num_spheres", "num_cylinders", "num_cones",
                          "num_tori",   "num_skipped", "warnings"};
    mxArray *values[] = {mio_o_scalar((double)r.num_planes),    mio_o_scalar((double)r.num_spheres),
                         mio_o_scalar((double)r.num_cylinders), mio_o_scalar((double)r.num_cones),
                         mio_o_scalar((double)r.num_tori),      mio_o_scalar((double)r.num_skipped),
                         w};
    return mio_o_struct(7, keys, values);
}

/* silhouette(mesh, normal, clearance, round_angle, arc_segments, extrude,
 * region, merge_tolerance) (roadmap 7.6.3). */
static mxArray *O_mio_silhouette(const mxArray **a) {
    mio_silhouette_opts opts;
    mio_silhouette_opts_init(&opts);
    mio_o_vec3(a[1], opts.normal, "normal");
    opts.clearance = mio_o_double(a[2], "clearance");
    opts.round_angle = mio_o_double(a[3], "round_angle");
    opts.arc_segments = (int64_t)mio_o_double(a[4], "arc_segments");
    opts.extrude = mio_o_bool(a[5], "extrude");
    opts.region = mio_o_opt_string(a[6]);
    opts.merge_tolerance = mio_o_double(a[7], "merge_tolerance");
    mio_silhouette_report r;
    mio_mesh *out = mio_silhouette(mio_o_mesh(a[0]), &opts, &r);
    if (out == NULL) mio_o_fail("silhouette");
    mxArray *w = mio_o_getbuffer(mio_silhouette_warnings, "silhouette");
    const char *keys[] = {"mesh",      "num_loops", "num_outer",  "num_holes",  "num_triangles",
                          "area",      "perimeter", "min_height", "max_height", "warnings"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)r.num_loops),
                         mio_o_scalar((double)r.num_outer),
                         mio_o_scalar((double)r.num_holes),
                         mio_o_scalar((double)r.num_triangles),
                         mio_o_scalar(r.area),
                         mio_o_scalar(r.perimeter),
                         mio_o_scalar(r.min_height),
                         mio_o_scalar(r.max_height),
                         w};
    return mio_o_struct(10, keys, values);
}

/* fit_curves(sections, tolerance, snap_angle, normal, min_arc_points,
 * min_arc_sweep): `pieces` is an n x 1 struct array. */
static mxArray *O_mio_fit_curves(const mxArray **a) {
    mio_fit_curves_opts opts;
    mio_fit_curves_opts_init(&opts);
    opts.tolerance = mio_o_double(a[1], "tolerance");
    opts.snap_angle = mio_o_double(a[2], "snap_angle");
    mio_o_vec3(a[3], opts.normal, "normal");
    opts.min_arc_points = (int64_t)mio_o_double(a[4], "min_arc_points");
    opts.min_arc_sweep = mio_o_double(a[5], "min_arc_sweep");
    mio_fit_curves_report r;
    mio_mesh *out = mio_fit_curves(mio_o_mesh(a[0]), &opts, &r);
    if (out == NULL) mio_o_fail("fit_curves");
    static const char *const kFields[] = {"type",     "start",    "end",        "centre",
                                          "radius",   "sweep",    "residual",   "num_points",
                                          "chain",    "frame"};
    static const char *const kTypes[] = {"line", "arc", "circle"};
    mxArray *pieces = mxCreateStructMatrix((mwSize)r.num_pieces, 1, 10, kFields);
    for (int64_t i = 0; i < r.num_pieces; ++i) {
        mio_sketch_piece p;
        if (mio_fit_curves_piece(i, &p) != MIO_OK) mio_o_fail("fit_curves");
        mxArray *st = mio_o_matrix(1, 2), *en = mio_o_matrix(1, 2), *ce = mio_o_matrix(1, 2);
        for (int k = 0; k < 2; ++k) {
            mxGetPr(st)[k] = p.start[k];
            mxGetPr(en)[k] = p.end[k];
            mxGetPr(ce)[k] = p.centre[k];
        }
        mxSetField(pieces, (mwIndex)i, "type",
                   mio_o_str(kTypes[(p.type >= 1 && p.type <= 3) ? p.type - 1 : 0]));
        mxSetField(pieces, (mwIndex)i, "start", st);
        mxSetField(pieces, (mwIndex)i, "end", en);
        mxSetField(pieces, (mwIndex)i, "centre", ce);
        mxSetField(pieces, (mwIndex)i, "radius", mio_o_scalar(p.radius));
        mxSetField(pieces, (mwIndex)i, "sweep", mio_o_scalar(p.sweep));
        mxSetField(pieces, (mwIndex)i, "residual", mio_o_scalar(p.residual));
        mxSetField(pieces, (mwIndex)i, "num_points", mio_o_scalar((double)p.num_points));
        mxSetField(pieces, (mwIndex)i, "chain", mio_o_scalar((double)p.chain));
        mxSetField(pieces, (mwIndex)i, "frame", mio_o_scalar((double)p.frame));
    }
    mxArray *w = mio_o_getbuffer(mio_fit_curves_warnings, "fit_curves");
    const char *keys[] = {"mesh",        "pieces",       "num_lines", "num_arcs",
                          "num_circles", "max_residual", "tolerance", "warnings"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         pieces,
                         mio_o_scalar((double)r.num_lines),
                         mio_o_scalar((double)r.num_arcs),
                         mio_o_scalar((double)r.num_circles),
                         mio_o_scalar(r.max_residual),
                         mio_o_scalar(r.tolerance),
                         w};
    return mio_o_struct(8, keys, values);
}

/* sections(mesh, direction, offsets, min_length, record_parent_ids) */
static mxArray *O_mio_sections(const mxArray **a) {
    double d[3];
    mio_o_vec3(a[1], d, "direction");
    size_t n = 0;
    double *v = mio_o_doubles(a[2], &n, "offsets");
    mio_mesh *out = mio_sections(mio_o_mesh(a[0]), d, v, (int64_t)n,
                                 mio_o_double(a[3], "min_length"),
                                 mio_o_bool(a[4], "record_parent_ids"));
    if (out == NULL) mio_o_fail("sections");
    return mio_o_wrap_mesh(out);
}

/* isosurface(mesh, array, isovalues, component, record_parent_ids) */
static mxArray *O_mio_isosurface(const mxArray **a) {
    size_t n = 0;
    double *v = mio_o_doubles(a[2], &n, "isovalues");
    if (n < 1) mio_o_error("isosurface", "at least one isovalue is required");
    mio_mesh *out = mio_isosurface(mio_o_mesh(a[0]), mio_o_string(a[1], "array"), v, (int)n,
                                   mio_o_int(a[3], "component"),
                                   mio_o_bool(a[4], "record_parent_ids"));
    if (out == NULL) mio_o_fail("isosurface");
    return mio_o_wrap_mesh(out);
}

/* gradient(mesh, array, op, method, location, output, component, overwrite) */
static mxArray *O_mio_gradient(const mxArray **a) {
    int64_t skipped = 0, fallback = 0;
    mio_mesh *out = mio_gradient(mio_o_mesh(a[0]), mio_o_string(a[1], "array"),
                                 mio_o_opt_string(a[2]), mio_o_opt_string(a[3]),
                                 mio_o_opt_string(a[4]), mio_o_opt_string(a[5]),
                                 mio_o_int(a[6], "component"), mio_o_bool(a[7], "overwrite"),
                                 &skipped, &fallback);
    if (out == NULL) mio_o_fail("gradient");
    const char *names[] = {"mesh", "num_skipped", "num_fallback"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)skipped),
                         mio_o_scalar((double)fallback)};
    return mio_o_struct(3, names, values);
}

/* hessian(mesh, array, method, location, output, overwrite) */
static mxArray *O_mio_hessian(const mxArray **a) {
    int64_t skipped = 0, fallback = 0;
    mio_mesh *out = mio_hessian(mio_o_mesh(a[0]), mio_o_string(a[1], "array"),
                                mio_o_opt_string(a[2]), mio_o_opt_string(a[3]),
                                mio_o_opt_string(a[4]), mio_o_bool(a[5], "overwrite"), &skipped,
                                &fallback);
    if (out == NULL) mio_o_fail("hessian");
    const char *names[] = {"mesh", "num_skipped", "num_fallback"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)skipped),
                         mio_o_scalar((double)fallback)};
    return mio_o_struct(3, names, values);
}

/* estimate_error(mesh, array, method, marking, marking_value, output, marked,
 *                overwrite) */
static mxArray *O_mio_estimate_error(const mxArray **a) {
    double global_error = 0.0;
    int64_t skipped = 0, nmarked = 0;
    mio_mesh *out = mio_estimate_error(
        mio_o_mesh(a[0]), mio_o_string(a[1], "array"), mio_o_opt_string(a[2]),
        mio_o_opt_string(a[3]), mio_o_double(a[4], "marking_value"), mio_o_opt_string(a[5]),
        mio_o_opt_string(a[6]), mio_o_bool(a[7], "overwrite"), &global_error, &skipped,
        &nmarked);
    if (out == NULL) mio_o_fail("estimate_error");
    const char *names[] = {"mesh", "global_error", "num_skipped", "num_marked"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar(global_error),
                         mio_o_scalar((double)skipped), mio_o_scalar((double)nmarked)};
    return mio_o_struct(4, names, values);
}

/* Goes through mio_remesh_ex/mio_remesh_opts rather than the flat mio_remesh
   -- the mio_refine_ex precedent, needed because mio_remesh is a flat C
   function with no room to grow (it already changed once). Like R, this links
   the real C header, so mio_remesh_opts is used directly -- no hand-mirrored
   struct to keep in sync.

   remesh(mesh, num_clusters, subdivide, subsample_ratio, max_subdivide,
          max_iterations, max_repair_passes, metric, gradation,
          preserve_boundary, max_anisotropy) */
static mxArray *O_mio_remesh(const mxArray **a) {
    mio_remesh_opts opts;
    mio_remesh_report report;

    mio_remesh_opts_init(&opts);
    opts.num_clusters = (int64_t)mio_o_int(a[1], "num_clusters");
    opts.subdivide = (int32_t)mio_o_int(a[2], "subdivide");
    opts.max_subdivide = (int32_t)mio_o_int(a[4], "max_subdivide");
    opts.subsample_ratio = mio_o_double(a[3], "subsample_ratio");
    opts.max_iterations = (int32_t)mio_o_int(a[5], "max_iterations");
    opts.max_repair_passes = (int32_t)mio_o_int(a[6], "max_repair_passes");
    opts.metric = mio_o_opt_string(a[7]);
    opts.gradation = mio_o_double(a[8], "gradation");
    opts.preserve_boundary = mio_o_bool(a[9], "preserve_boundary") ? 1 : 0;
    opts.max_anisotropy = mio_o_double(a[10], "max_anisotropy");

    mio_mesh *out = mio_remesh_ex(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("remesh");
    const char *names[] = {"mesh", "num_clusters", "num_iterations", "subdivide_applied",
                           "num_isolated_clusters", "num_non_manifold_vertices"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)report.num_clusters),
                         mio_o_scalar((double)report.num_iterations),
                         mio_o_scalar((double)report.subdivide_applied),
                         mio_o_scalar((double)report.num_isolated_clusters),
                         mio_o_scalar((double)report.num_non_manifold_vertices)};
    return mio_o_struct(6, names, values);
}

/* compute_curvature(mesh, mean, gaussian, dual_area, include_boundary,
 *                   record_area, record_principal, region) */
static mxArray *O_mio_compute_curvature(const mxArray **a) {
    mio_curvature_opts opts;
    mio_curvature_report report;

    mio_curvature_opts_init(&opts);
    opts.mean = mio_o_bool(a[1], "mean") ? 1 : 0;
    opts.gaussian = mio_o_bool(a[2], "gaussian") ? 1 : 0;
    const char *da = mio_o_opt_string(a[3]);
    if (da != NULL && strcmp(da, "barycentric") == 0) {
        opts.dual_area = MIO_CURVATURE_BARYCENTRIC;
    } else if (da != NULL && strcmp(da, "mixed-voronoi") != 0) {
        mio_o_error("curvature",
                    "meshio++: curvature: unknown dual area '%s' "
                    "(expected 'mixed-voronoi' or 'barycentric')",
                    da);
    }
    opts.include_boundary = mio_o_bool(a[4], "include_boundary") ? 1 : 0;
    opts.record_area = mio_o_bool(a[5], "record_area") ? 1 : 0;
    opts.record_principal = mio_o_bool(a[6], "record_principal") ? 1 : 0;
    opts.region = mio_o_opt_string(a[7]);

    mio_mesh *out = mio_compute_curvature(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("curvature");
    const char *names[] = {"mesh", "quality", "num_boundary", "num_isolated",
                           "num_degenerate", "total_angle_defect"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_surface_quality(&report.quality),
                         mio_o_scalar((double)report.num_boundary),
                         mio_o_scalar((double)report.num_isolated),
                         mio_o_scalar((double)report.num_degenerate),
                         mio_o_scalar(report.total_angle_defect)};
    return mio_o_struct(6, names, values);
}

/* compute_normals(mesh, point_normals, cell_normals, weight, split_angle,
 *                 record_parent_ids, region) */
static mxArray *O_mio_compute_normals(const mxArray **a) {
    mio_normals_opts opts;
    mio_normals_report report;

    mio_normals_opts_init(&opts);
    opts.point_normals = mio_o_bool(a[1], "point_normals") ? 1 : 0;
    opts.cell_normals = mio_o_bool(a[2], "cell_normals") ? 1 : 0;
    const char *w = mio_o_opt_string(a[3]);
    if (w != NULL && strcmp(w, "area") == 0) {
        opts.weight = MIO_SDF_WEIGHT_AREA;
    } else if (w != NULL && strcmp(w, "angle") != 0) {
        mio_o_error("normals", "meshio++: normals: unknown weight '%s' (expected 'angle' or 'area')",
                    w);
    }
    /* [] means "no split"; a number is the crease angle in degrees. */
    if (!mio_o_is_null(a[4])) {
        opts.split = 1;
        opts.split_angle = mio_o_double(a[4], "split_angle");
    }
    opts.record_parent_ids = mio_o_bool(a[5], "record_parent_ids") ? 1 : 0;
    opts.region = mio_o_opt_string(a[6]);

    mio_mesh *out = mio_compute_normals(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("normals");
    const char *names[] = {"mesh", "quality", "num_isolated", "num_undefined",
                           "num_degenerate", "num_split_points", "num_added_points"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_surface_quality(&report.quality),
                         mio_o_scalar((double)report.num_isolated),
                         mio_o_scalar((double)report.num_undefined),
                         mio_o_scalar((double)report.num_degenerate),
                         mio_o_scalar((double)report.num_split_points),
                         mio_o_scalar((double)report.num_added_points)};
    return mio_o_struct(7, names, values);
}

/* feature_edges(mesh, feature_angle, feature, boundary, non_manifold,
 *               inconsistent, region) */
static mxArray *O_mio_feature_edges(const mxArray **a) {
    mio_feature_edges_opts opts;
    mio_feature_edges_report report;

    mio_feature_edges_opts_init(&opts);
    opts.feature_angle = mio_o_double(a[1], "feature_angle");
    opts.feature = mio_o_bool(a[2], "feature") ? 1 : 0;
    opts.boundary = mio_o_bool(a[3], "boundary") ? 1 : 0;
    opts.non_manifold = mio_o_bool(a[4], "non_manifold") ? 1 : 0;
    opts.inconsistent = mio_o_bool(a[5], "inconsistent") ? 1 : 0;
    opts.region = mio_o_opt_string(a[6]);

    mio_mesh *out = mio_feature_edges(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("feature_edges");
    const char *names[] = {"mesh", "num_feature", "num_boundary", "num_non_manifold",
                           "num_inconsistent"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)report.num_feature),
                         mio_o_scalar((double)report.num_boundary),
                         mio_o_scalar((double)report.num_non_manifold),
                         mio_o_scalar((double)report.num_inconsistent)};
    return mio_o_struct(5, names, values);
}

/* segment(mesh, feature_angle, normal_tolerance, curvature_tolerance,
 *         min_patch_cells, region, prefix, method, tolerance) */
static mxArray *O_mio_segment(const mxArray **a) {
    mio_segment_opts opts;
    mio_segment_report report;

    mio_segment_opts_init(&opts);
    opts.feature_angle = mio_o_double(a[1], "feature_angle");
    opts.normal_tolerance = mio_o_double(a[2], "normal_tolerance");
    opts.curvature_tolerance = mio_o_double(a[3], "curvature_tolerance");
    opts.min_patch_cells = (int64_t)mio_o_int(a[4], "min_patch_cells");
    opts.region = mio_o_opt_string(a[5]);
    opts.prefix = mio_o_opt_string(a[6]);
    {
        const char *m = mio_o_string(a[7], "method");
        if (strcmp(m, "angle") == 0)
            opts.method = MIO_SEGMENT_ANGLE;
        else if (strcmp(m, "primitive") == 0)
            opts.method = MIO_SEGMENT_PRIMITIVE;
        else
            mio_o_error("segment", "meshio++: segment: method must be \"angle\" or \"primitive\"");
    }
    opts.tolerance = mio_o_double(a[8], "tolerance");

    mio_mesh *out = mio_segment(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("segment");
    const char *names[] = {"mesh", "num_patches", "num_merged", "num_feature_edges",
                           "num_boundary_edges"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)report.num_patches),
                         mio_o_scalar((double)report.num_merged),
                         mio_o_scalar((double)report.num_feature_edges),
                         mio_o_scalar((double)report.num_boundary_edges)};
    return mio_o_struct(5, names, values);
}

/* fit_primitives(mesh, patch_data, tolerance, relative, min_inlier_fraction,
 *                types, max_iterations, region, seeds, feature_angle,
 *                rim_angle, output_region, used_points, criterion) -- `types` the bit
 * mask and `seeds` the flat x, y, z list the .m wrapper builds. */
static mxArray *O_mio_fit_primitives(const mxArray **a) {
    mio_fit_primitives_ex_opts opts;
    mio_fit_primitives_report report;

    mio_fit_primitives_ex_opts_init(&opts);
    opts.base.patch_data = mio_o_string(a[1], "patch_data");
    opts.base.tolerance = mio_o_double(a[2], "tolerance");
    opts.base.relative = mio_o_bool(a[3], "relative") ? 1 : 0;
    opts.base.min_inlier_fraction = mio_o_double(a[4], "min_inlier_fraction");
    opts.base.types = (uint32_t)mio_o_int(a[5], "types");
    opts.base.max_iterations = (int64_t)mio_o_int(a[6], "max_iterations");
    opts.region = mio_o_opt_string(a[7]);
    if (!mio_o_is_null(a[8])) {
        size_t n = 0;
        double *seeds = mio_o_doubles(a[8], &n, "seeds");
        if (n % 3 != 0) mio_o_error("fit_primitives", "`seeds` must hold x, y, z triples");
        opts.seeds = seeds;
        opts.num_seeds = (int64_t)(n / 3);
    }
    opts.feature_angle = mio_o_double(a[9], "feature_angle");
    opts.rim_angle = mio_o_double(a[10], "rim_angle");
    opts.output_region = mio_o_opt_string(a[11]);
    opts.used_points = (int32_t)mio_o_int(a[12], "used_points");
    opts.criterion = (int32_t)mio_o_int(a[13], "criterion");

    mio_mesh *out = mio_fit_primitives_ex(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("fit_primitives");
    mxArray *w = mio_o_getbuffer(mio_fit_primitives_warnings, "fit_primitives");
    const char *names[] = {"mesh",         "num_patches", "num_plane",    "num_sphere",
                           "num_cylinder", "num_cone",    "num_torus",    "num_unfitted",
                           "num_feature_cells", "warnings"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)report.num_patches),
                         mio_o_scalar((double)report.num_plane),
                         mio_o_scalar((double)report.num_sphere),
                         mio_o_scalar((double)report.num_cylinder),
                         mio_o_scalar((double)report.num_cone),
                         mio_o_scalar((double)report.num_torus),
                         mio_o_scalar((double)report.num_unfitted),
                         mio_o_scalar((double)report.num_feature_cells),
                         w};
    return mio_o_struct(10, names, values);
}

static mxArray *O_mio_has_occt(const mxArray **a) {
    (void)a;
    return mio_o_logical(mio_has_occt() != 0);
}

/* A RANSAC seed: [] takes the default; a number must be a whole number in
 * [0, 2^64), held exactly below 2^53 (Octave numbers are doubles). */
static uint64_t oc_refine_seed(const mxArray *s) {
    if (mio_o_is_null(s)) return MIO_REFINE_DEFAULT_SEED;
    if (mxGetNumberOfElements(s) != 1 || !(mxIsNumeric(s) || mxIsLogical(s)) || mxIsComplex(s)) {
        mio_o_error("ransac_seed", "`ransac_seed` must be [] or a single number");
    }
    double v = mxGetScalar(s);
    if (isnan(v) || v < 0.0 || v >= 18446744073709551616.0 || v != floor(v)) {
        mio_o_error("ransac_seed", "`ransac_seed` must be a whole number in [0, 2^64)");
    }
    return (uint64_t)v;
}

/* The mio_refine_stop code as a string. */
static const char *oc_refine_stop_name(int32_t stop) {
    switch (stop) {
    case MIO_REFINE_UNCHANGED: return "unchanged";
    case MIO_REFINE_SMALL_GAIN: return "small_gain";
    case MIO_REFINE_ITERATION_CAP: return "iteration_cap";
    case MIO_REFINE_ENERGY_ROSE: return "energy_rose";
    default: return "unknown";
    }
}

/* refine_patches(mesh, tolerance, relative, normal_weight, epsilon,
 *                feature_angle, max_iterations, min_patch_cells, ransac_seed,
 *                types, allow_split, allow_merge, prefix) */
static mxArray *O_mio_refine_patches(const mxArray **a) {
    mio_refine_patches_opts opts;
    mio_refine_patches_report report;

    mio_refine_patches_opts_init(&opts);
    opts.tolerance = mio_o_double(a[1], "tolerance");
    opts.relative = mio_o_bool(a[2], "relative") ? 1 : 0;
    opts.normal_weight = mio_o_double(a[3], "normal_weight");
    opts.epsilon = mio_o_double(a[4], "epsilon");
    opts.feature_angle = mio_o_double(a[5], "feature_angle");
    opts.max_iterations = mio_o_int64(a[6], "max_iterations");
    opts.min_patch_cells = mio_o_int64(a[7], "min_patch_cells");
    opts.ransac_seed = oc_refine_seed(a[8]);
    opts.types = (uint32_t)mio_o_int(a[9], "types");
    opts.allow_split = mio_o_bool(a[10], "allow_split") ? 1 : 0;
    opts.allow_merge = mio_o_bool(a[11], "allow_merge") ? 1 : 0;
    opts.prefix = mio_o_opt_string(a[12]);

    mio_mesh *out = mio_refine_patches(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("refine_patches");
    const char *names[] = {"mesh",
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
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)report.num_patches),
                         mio_o_scalar((double)report.num_iterations),
                         mio_o_scalar((double)report.num_splits),
                         mio_o_scalar((double)report.num_merges),
                         mio_o_scalar((double)report.num_absorbed),
                         mio_o_scalar((double)report.num_primitive),
                         mio_o_scalar((double)report.num_freeform),
                         mio_o_scalar((double)report.num_facets),
                         mio_o_scalar(report.tolerance),
                         mio_o_scalar(report.explained_area_fraction),
                         mio_o_scalar(report.initial_explained_area_fraction),
                         mio_o_scalar(report.hausdorff_max),
                         mio_o_scalar(report.hausdorff_rms),
                         mio_o_scalar(report.initial_energy),
                         mio_o_scalar(report.energy),
                         mio_o_str(oc_refine_stop_name(report.stop))};
    return mio_o_struct(17, names, values);
}

/* The 0-based position of a keyword option in a NULL-terminated list. */
static int32_t oc_fit_cad_keyword(const mxArray *value, const char *name,
                                  const char *const *words) {
    const char *v = mio_o_string(value, name);
    for (int32_t i = 0; words[i] != NULL; ++i)
        if (strcmp(v, words[i]) == 0) return i;
    mio_o_error("fit_cad", "meshio++: fit_cad: %s: unknown value \"%s\"", name, v);
    return 0;
}

/* fit_cad(mesh, step_path, tolerance, exact_edges, types, mesh_deflection,
 *         max_demotion_rounds, freeform, freeform_degree, freeform_max_control,
 *         freeform_split, fillets, threads, merge_faces, backend,
 *         refine_segmentation, max_refine_iterations, refine_epsilon,
 *         normal_weight, allow_merge, allow_split, ransac_seed,
 *         refine_tolerance, relative_tolerance, freeform_general,
 *         freeform_tolerance, freeform_smoothing, max_freeform_poles,
 *         max_freeform_iterations, max_freeform_splits, sew_tolerance,
 *         extend_surfaces, snap_junctions, face_deviation, freeform_mode,
 *         freeform_parameterization, extraordinary_vertex, local_refine,
 *         freeform_g1_tolerance_deg, freeform_max_control_points,
 *         layout_target_patch_size, layout_max_valence_deviation,
 *         resegment_rounds, layout_snap_to_features,
 *         validate_by_surface_sampling, constrain_primitives, validation_report,
 *         validation_acceptance) -- 48 arguments, R's order exactly. */
static mxArray *O_mio_fit_cad(const mxArray **a) {
    mio_fit_cad_opts opts;
    mio_fit_cad_opts_ex opts_ex;
    mio_fit_cad_opts_ex2 opts_ex2;
    mio_fit_cad_report report;
    mio_fit_cad_report_ex report_ex;
    mio_fit_cad_report_ex2 report_ex2;
    mio_fit_cad_report_ex3 report_ex3;
    mio_fit_cad_report_ex4 report_ex4;
    static const char *const modes[] = {"legacy", "adaptive", "network", "auto", NULL};
    static const char *const params[] = {"auto", "projection", "lscm", "harmonic", NULL};
    static const char *const evs[] = {"g0", "g1", "g1x", "g1xx", NULL};
    static const char *const refines[] = {"knots", "thb", NULL};

    mio_fit_cad_opts_init(&opts);
    mio_fit_cad_opts_ex_init(&opts_ex);
    mio_fit_cad_opts_ex2_init(&opts_ex2);
    opts_ex2.freeform_mode = oc_fit_cad_keyword(a[34], "freeform_mode", modes);
    opts_ex2.freeform_parameterization =
        oc_fit_cad_keyword(a[35], "freeform_parameterization", params);
    opts_ex2.extraordinary_vertex = oc_fit_cad_keyword(a[36], "extraordinary_vertex", evs);
    opts_ex2.local_refine = oc_fit_cad_keyword(a[37], "local_refine", refines);
    opts_ex2.freeform_g1_tolerance_deg = mio_o_double(a[38], "freeform_g1_tolerance_deg");
    opts_ex2.freeform_max_control_points = mio_o_int64(a[39], "freeform_max_control_points");
    opts_ex2.layout_target_patch_size = mio_o_double(a[40], "layout_target_patch_size");
    opts_ex2.layout_max_valence_deviation = mio_o_int64(a[41], "layout_max_valence_deviation");
    opts_ex2.resegment_rounds = mio_o_int64(a[42], "resegment_rounds");
    opts_ex2.layout_snap_to_features = mio_o_bool(a[43], "layout_snap_to_features") ? 1 : 0;
    opts_ex2.validate_by_surface_sampling =
        mio_o_bool(a[44], "validate_by_surface_sampling") ? 1 : 0;
    opts_ex2.constrain_primitives = mio_o_bool(a[45], "constrain_primitives") ? 1 : 0;
    opts_ex2.validation_report = mio_o_bool(a[46], "validation_report") ? 1 : 0;
    opts_ex2.validation_acceptance = mio_o_bool(a[47], "validation_acceptance") ? 1 : 0;
    opts_ex.refine_segmentation = mio_o_bool(a[15], "refine_segmentation") ? 1 : 0;
    opts_ex.refine.max_iterations = mio_o_int64(a[16], "max_refine_iterations");
    opts_ex.refine.epsilon = mio_o_double(a[17], "refine_epsilon");
    opts_ex.refine.normal_weight = mio_o_double(a[18], "normal_weight");
    opts_ex.refine.allow_merge = mio_o_bool(a[19], "allow_merge") ? 1 : 0;
    opts_ex.refine.allow_split = mio_o_bool(a[20], "allow_split") ? 1 : 0;
    opts_ex.refine.ransac_seed = oc_refine_seed(a[21]);
    opts_ex.refine.tolerance = mio_o_double(a[22], "refine_tolerance");
    opts_ex.refine.relative = mio_o_bool(a[23], "relative_tolerance") ? 1 : 0;
    opts_ex.freeform_general = mio_o_bool(a[24], "freeform_general") ? 1 : 0;
    opts_ex.freeform_tolerance = mio_o_double(a[25], "freeform_tolerance");
    opts_ex.freeform_smoothing = mio_o_double(a[26], "freeform_smoothing");
    opts_ex.max_freeform_poles = mio_o_int64(a[27], "max_freeform_poles");
    opts_ex.max_freeform_iterations = mio_o_int64(a[28], "max_freeform_iterations");
    opts_ex.max_freeform_splits = mio_o_int64(a[29], "max_freeform_splits");
    opts_ex.sew_tolerance = mio_o_double(a[30], "sew_tolerance");
    opts_ex.extend_surfaces = mio_o_bool(a[31], "extend_surfaces") ? 1 : 0;
    opts_ex.snap_junctions = mio_o_bool(a[32], "snap_junctions") ? 1 : 0;
    opts_ex.face_deviation = mio_o_bool(a[33], "face_deviation") ? 1 : 0;
    opts.tolerance = mio_o_double(a[2], "tolerance");
    opts.exact_edges = mio_o_bool(a[3], "exact_edges") ? 1 : 0;
    opts.types = (uint32_t)mio_o_int(a[4], "types");
    opts.mesh_deflection = mio_o_double(a[5], "mesh_deflection");
    opts.max_demotion_rounds = (int64_t)mio_o_int(a[6], "max_demotion_rounds");
    opts.freeform = mio_o_bool(a[7], "freeform") ? 1 : 0;
    opts.freeform_degree = (int32_t)mio_o_int(a[8], "freeform_degree");
    opts.freeform_max_control = (int64_t)mio_o_int(a[9], "freeform_max_control");
    opts.freeform_split = mio_o_bool(a[10], "freeform_split") ? 1 : 0;
    opts.fillets = mio_o_bool(a[11], "fillets") ? 1 : 0;
    opts.threads = mio_o_bool(a[12], "threads") ? 1 : 0;
    opts.merge_faces = mio_o_bool(a[13], "merge_faces") ? 1 : 0;
    {
        const char *b = mio_o_string(a[14], "backend");
        if (strcmp(b, "auto") == 0)
            opts.backend = MIO_FIT_CAD_AUTO;
        else if (strcmp(b, "occt") == 0)
            opts.backend = MIO_FIT_CAD_OCCT;
        else if (strcmp(b, "native") == 0)
            opts.backend = MIO_FIT_CAD_NATIVE;
        else
            mio_o_error("fit_cad",
                        "meshio++: fit_cad: backend must be \"auto\", \"occt\" or \"native\"");
    }

    mio_mesh *out = mio_fit_cad4(mio_o_mesh(a[0]), mio_o_string(a[1], "step_path"), &opts,
                                 &opts_ex, &opts_ex2, &report, &report_ex, &report_ex2,
                                 &report_ex3, &report_ex4);
    if (out == NULL) mio_o_fail("fit_cad");
    const char *names[] = {"mesh",
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
                           "num_constraint_groups",
                           "validation_reported",
                           "cad_to_mesh_max",
                           "cad_to_mesh_p95",
                           "cad_to_mesh_p99",
                           "mesh_to_cad_max",
                           "mesh_to_cad_p95",
                           "mesh_to_cad_p99",
                           "normal_deviation_max_deg",
                           "normal_deviation_p95_deg",
                           "step_round_trip_checked",
                           "step_round_trip_valid",
                           "step_round_trip_volume_ratio",
                           "total_freeform_control_points"};
    mxArray *values[] = {
        mio_o_wrap_mesh(out),
        mio_o_logical(report.valid != 0),
        mio_o_scalar(report.tolerance),
        mio_o_scalar(report.volume),
        mio_o_scalar(report.skin_volume),
        mio_o_scalar(report.volume_ratio),
        mio_o_scalar(report.max_deviation),
        mio_o_scalar((double)report.num_faces),
        mio_o_scalar((double)report.num_plane),
        mio_o_scalar((double)report.num_sphere),
        mio_o_scalar((double)report.num_cylinder),
        mio_o_scalar((double)report.num_cone),
        mio_o_scalar((double)report.num_torus),
        mio_o_scalar((double)report.num_faceted_faces),
        mio_o_scalar((double)report.num_patches),
        mio_o_scalar((double)report.num_demoted),
        mio_o_scalar((double)report.num_exact_edges),
        mio_o_scalar((double)report.num_approx_edges),
        mio_o_scalar((double)report.num_solids),
        mio_o_scalar((double)report.num_bspline),
        mio_o_scalar((double)report.num_split_patches),
        mio_o_scalar((double)report.num_fillet),
        mio_o_scalar((double)report.num_helical),
        mio_o_scalar((double)report.num_variable_fillet),
        mio_o_scalar((double)report.num_blend),
        mio_o_scalar((double)report.num_rounded_thread),
        mio_o_scalar((double)report_ex.num_vertex_blend),
        mio_o_scalar((double)report_ex.num_chamfer),
        mio_o_scalar((double)report_ex.num_conic_blend),
        mio_o_scalar((double)report_ex.num_tapered_thread),
        mio_o_scalar((double)report_ex.num_variable_pitch_thread),
        mio_o_str(report_ex.backend == MIO_FIT_CAD_NATIVE ? "native" : "occt"),
        mio_o_logical(report_ex2.refined != 0),
        mio_o_scalar((double)report_ex2.refine_iterations),
        mio_o_scalar((double)report_ex2.num_refine_splits),
        mio_o_scalar((double)report_ex2.num_refine_merges),
        mio_o_scalar(report_ex2.initial_explained_area_fraction),
        mio_o_scalar(report_ex2.explained_area_fraction),
        mio_o_logical(report_ex2.closed != 0),
        mio_o_scalar((double)report_ex2.num_freeform_general),
        mio_o_scalar((double)report_ex2.num_freeform_refinements),
        mio_o_scalar(report_ex2.avg_freeform_poles),
        mio_o_scalar((double)report_ex2.num_extended_edges),
        mio_o_scalar((double)report_ex2.num_snapped_junctions),
        mio_o_scalar(report_ex2.hausdorff_max),
        mio_o_scalar(report_ex2.hausdorff_rms),
        mio_o_logical(report_ex2.refinement_kept != 0),
        mio_o_logical(report_ex3.brep_checked != 0),
        mio_o_str(report_ex3.freeform_mode >= 0 && report_ex3.freeform_mode <= 3
                      ? modes[report_ex3.freeform_mode]
                      : "legacy"),
        mio_o_scalar(report_ex3.facet_area_fraction),
        mio_o_scalar((double)report_ex3.num_freeform_adaptive),
        mio_o_scalar((double)report_ex3.num_network_patches),
        mio_o_scalar((double)report_ex3.num_network_edges),
        mio_o_scalar((double)report_ex3.num_extraordinary_vertices),
        mio_o_scalar((double)report_ex3.num_resegmented),
        mio_o_scalar((double)report_ex3.num_escalations),
        mio_o_scalar(report_ex3.max_g1_angle_deg),
        mio_o_scalar(report_ex3.max_freeform_distance),
        mio_o_scalar(report_ex3.rms_freeform_distance),
        mio_o_scalar((double)report_ex3.num_face_records),
        mio_o_scalar((double)report_ex3.num_constrained_primitives),
        mio_o_scalar((double)report_ex3.num_constraint_groups),
        mio_o_logical(report_ex4.validation_reported != 0),
        mio_o_scalar(report_ex4.cad_to_mesh_max),
        mio_o_scalar(report_ex4.cad_to_mesh_p95),
        mio_o_scalar(report_ex4.cad_to_mesh_p99),
        mio_o_scalar(report_ex4.mesh_to_cad_max),
        mio_o_scalar(report_ex4.mesh_to_cad_p95),
        mio_o_scalar(report_ex4.mesh_to_cad_p99),
        mio_o_scalar(report_ex4.normal_deviation_max_deg),
        mio_o_scalar(report_ex4.normal_deviation_p95_deg),
        mio_o_logical(report_ex4.step_round_trip_checked != 0),
        mio_o_logical(report_ex4.step_round_trip_valid != 0),
        mio_o_scalar(report_ex4.step_round_trip_volume_ratio),
        mio_o_scalar((double)report_ex4.total_freeform_control_points)};
    return mio_o_struct(75, names, values);
}

/* --- interfaces between regions ------------------------------------------ */

/* region_adjacency(mesh, regions) -- [] or {} for every Cell region. */
static mxArray *O_mio_region_adjacency(const mxArray **a) {
    int64_t n = 0;
    mio_region_selector *selectors = NULL;
    if (!mio_o_is_null(a[1])) {
        if (!mxIsCell(a[1]) && !mxIsChar(a[1])) {
            mio_o_error("region_adjacency", "`regions` must be [] or a cell array of strings");
        }
        const char *const *names = mio_o_names(a[1], &n, "regions");
        if (n > 0) {
            selectors = (mio_region_selector *)mxCalloc((size_t)n, sizeof(*selectors));
        }
        for (int64_t i = 0; i < n; ++i) {
            mio_region_selector_init(&selectors[i], names[i]);
            selectors[i].kind = MIO_REGION_CELL;
        }
    }
    mio_mesh *out = mio_region_adjacency(mio_o_mesh(a[0]), selectors, n);
    if (out == NULL) mio_o_fail("region_adjacency");
    return mio_o_wrap_mesh(out);
}

/* These ids cross the flat C ABI in 0-based `(cell, local_facet)` pairs. The
 * region convention shifts only the cell column; facet ordinals stay 0-based.
 * One pair per ROW: (count x 2). */
static mxArray *oc_side_matrix(const int64_t *entries, int64_t count) {
    if (count < 0) count = 0;
    mxArray *out = mio_o_matrix(count, 2);
    double *dst = mxGetPr(out);
    for (int64_t i = 0; i < count; ++i) {
        dst[i] = (double)(entries[2 * i] + 1);
        dst[count + i] = (double)entries[2 * i + 1];
    }
    return out;
}

/* find_interface(mesh_a, region_a, region_b, mesh_b, mode, master,
 *                gap_tolerance, angle_tolerance, overlap_tolerance) */
static mxArray *O_mio_find_interface(const mxArray **a) {
    mio_find_interface_opts opts;
    mio_find_interface_opts_init(&opts);
    const char *mode_name = mio_o_opt_string(a[4]);
    if (mode_name != NULL && mode_name[0] != '\0') {
        if (strcmp(mode_name, "conforming") == 0)
            opts.mode = MIO_INTERFACE_CONFORMING;
        else if (strcmp(mode_name, "proximity") == 0)
            opts.mode = MIO_INTERFACE_PROXIMITY;
        else
            mio_o_error("find_interface", "`mode` must be 'conforming' or 'proximity'");
    }
    const char *master_name = mio_o_opt_string(a[5]);
    if (master_name != NULL && master_name[0] != '\0') {
        if (strcmp(master_name, "a") == 0)
            opts.master = MIO_INTERFACE_MASTER_A;
        else if (strcmp(master_name, "b") == 0)
            opts.master = MIO_INTERFACE_MASTER_B;
        else
            mio_o_error("find_interface", "`master` must be 'a' or 'b'");
    }
    opts.gap_tolerance = mio_o_double(a[6], "gap_tolerance");
    opts.angle_tolerance = mio_o_double(a[7], "angle_tolerance");
    opts.overlap_tolerance = mio_o_double(a[8], "overlap_tolerance");

    mio_region_selector selector_a, selector_b;
    mio_region_selector_init(&selector_a, mio_o_string(a[1], "region_a"));
    mio_region_selector_init(&selector_b, mio_o_string(a[2], "region_b"));
    selector_a.kind = selector_b.kind = MIO_REGION_CELL;
    const mio_mesh *b = mio_o_is_null(a[3]) ? NULL : mio_o_mesh(a[3]);
    const mio_mesh *ma = mio_o_mesh(a[0]);
    mio_find_interface_result *result =
        mio_find_interface(ma, &selector_a, b, &selector_b, &opts);
    if (result == NULL) mio_o_fail("find_interface");

    int64_t pairs = 0, unmatched_a = 0, unmatched_b = 0, count_a = 0, count_b = 0;
    double area = 0.0, max_gap = 0.0;
    if (mio_find_interface_result_report(result, &pairs, &area, &max_gap, &unmatched_a,
                                         &unmatched_b) != MIO_OK) {
        mio_find_interface_result_free(result);
        mio_o_fail("find_interface");
    }
    const int64_t *side_a = mio_find_interface_result_side_a(result, &count_a);
    const int64_t *side_b = mio_find_interface_result_side_b(result, &count_b);
    mxArray *side_a_value = oc_side_matrix(side_a, count_a);
    mxArray *side_b_value = oc_side_matrix(side_b, count_b);
    mio_mesh *out_mesh = mio_find_interface_result_take_mesh(result);
    mio_find_interface_result_free(result);
    if (out_mesh == NULL) mio_o_fail("find_interface");

    const char *names[] = {"mesh",        "num_pairs",   "area",   "max_gap",
                           "unmatched_a", "unmatched_b", "side_a", "side_b"};
    mxArray *values[] = {mio_o_wrap_mesh(out_mesh),
                         mio_o_scalar((double)pairs),
                         mio_o_scalar(area),
                         mio_o_scalar(max_gap),
                         mio_o_scalar((double)unmatched_a),
                         mio_o_scalar((double)unmatched_b),
                         side_a_value,
                         side_b_value};
    return mio_o_struct(8, names, values);
}

/* contact_pairs(slave_mesh, slave_points, master_mesh, master_cells,
 *               tolerance, require_complete) */
static mxArray *O_mio_contact_pairs(const mxArray **a) {
    mio_contact_pairs_opts opts;
    mio_contact_pairs_opts_init(&opts);
    opts.tolerance = mio_o_double(a[4], "tolerance");
    opts.require_complete = mio_o_bool(a[5], "require_complete") ? 1 : 0;
    mio_region_selector point_selector, cell_selector;
    mio_region_selector_init(&point_selector, mio_o_string(a[1], "slave_points"));
    mio_region_selector_init(&cell_selector, mio_o_string(a[3], "master_cells"));
    point_selector.kind = MIO_REGION_POINT;
    cell_selector.kind = MIO_REGION_CELL;
    const mio_mesh *master = mio_o_is_null(a[2]) ? NULL : mio_o_mesh(a[2]);
    const mio_mesh *slave = mio_o_mesh(a[0]);
    mio_contact_pairs_result *result =
        mio_contact_pairs(slave, &point_selector, master, &cell_selector, &opts);
    if (result == NULL) mio_o_fail("contact_pairs");
    int64_t count = 0, unmatched_count = 0;
    if (mio_contact_pairs_result_info(result, &count, &unmatched_count) != MIO_OK) {
        mio_contact_pairs_result_free(result);
        mio_o_fail("contact_pairs");
    }

    /* Ids 1-based (a 0 sentinel for "none"), facet ordinals 0-based, and the
     * 3-vectors one pair per ROW: (count x 3). */
    mxArray *vs = mio_o_shift_map(mio_contact_pairs_slave_point(result), count);
    mxArray *vmc = mio_o_shift_map(mio_contact_pairs_master_cell(result), count);
    mxArray *vmf = mio_o_copy_as_real(mio_contact_pairs_master_facet(result), MIO_INT64, count);
    mxArray *vms =
        mio_o_copy_as_real(mio_contact_pairs_master_subfacet(result), MIO_INT64, count);
    mxArray *vlc =
        mio_o_rows_from_c(mio_contact_pairs_local_coordinates(result), MIO_FLOAT64, count, 3, 0.0);
    mxArray *vcp =
        mio_o_rows_from_c(mio_contact_pairs_closest_point(result), MIO_FLOAT64, count, 3, 0.0);
    mxArray *vg = mio_o_copy_as_real(mio_contact_pairs_gap(result), MIO_FLOAT64, count);
    mxArray *vn = mio_o_rows_from_c(mio_contact_pairs_normal(result), MIO_FLOAT64, count, 3, 0.0);
    mxArray *vu = mio_o_shift_map(mio_contact_pairs_unmatched(result), unmatched_count);
    mio_contact_pairs_result_free(result);
    const char *names[] = {"slave_point",   "master_cell", "master_facet",
                           "master_subfacet", "local_coordinates", "closest_point",
                           "gap",           "normal",      "unmatched"};
    mxArray *values[] = {vs, vmc, vmf, vms, vlc, vcp, vg, vn, vu};
    return mio_o_struct(9, names, values);
}

/* split_interface(mesh, side, add_cohesive) */
static mxArray *O_mio_split_interface(const mxArray **a) {
    mio_split_interface_opts opts;
    mio_split_interface_opts_init(&opts);
    opts.add_cohesive = mio_o_bool(a[2], "add_cohesive") ? 1 : 0;
    mio_region_selector selector;
    mio_region_selector_init(&selector, mio_o_string(a[1], "side"));
    selector.kind = MIO_REGION_SIDE;
    int64_t duplicated = 0, cohesive = 0;
    mio_mesh *out_mesh =
        mio_split_interface(mio_o_mesh(a[0]), &selector, &opts, &duplicated, &cohesive);
    if (out_mesh == NULL) mio_o_fail("split_interface");
    const char *names[] = {"mesh", "num_duplicated_points", "num_cohesive_cells"};
    mxArray *values[] = {mio_o_wrap_mesh(out_mesh), mio_o_scalar((double)duplicated),
                         mio_o_scalar((double)cohesive)};
    return mio_o_struct(3, names, values);
}

/* --- quality gate, distances and repair ---------------------------------- */

/* check_quality(mesh, require, max_inverted, max_degenerate) */
static mxArray *O_mio_check_quality(const mxArray **a) {
    mio_quality_gate_report r;
    char buf[8192];
    buf[0] = '\0';
    mio_o_check(mio_check_quality(mio_o_mesh(a[0]), mio_o_opt_string(a[1]),
                                  (int64_t)mio_o_double(a[2], "max_inverted"),
                                  (int64_t)mio_o_double(a[3], "max_degenerate"), &r, buf,
                                  (int64_t)sizeof(buf)),
                "check_quality");
    const char *names[] = {"passed",       "num_checks",     "num_failed", "num_cells",
                           "num_inverted", "num_degenerate", "summary"};
    mxArray *values[] = {mio_o_logical(r.passed != 0),     mio_o_scalar((double)r.num_checks),
                         mio_o_scalar((double)r.num_failed), mio_o_scalar((double)r.num_cells),
                         mio_o_scalar((double)r.num_inverted),
                         mio_o_scalar((double)r.num_degenerate), mio_o_str(buf)};
    return mio_o_struct(7, names, values);
}

/* hausdorff_distance(a, b, face_samples, region_a, region_b) */
static mxArray *O_mio_hausdorff_distance(const mxArray **a) {
    mio_hausdorff_opts opts;
    mio_hausdorff_report r;

    mio_hausdorff_opts_init(&opts);
    opts.face_samples = (int64_t)mio_o_double(a[2], "face_samples");
    opts.region_a = mio_o_opt_string(a[3]);
    opts.region_b = mio_o_opt_string(a[4]);
    mio_o_check(mio_hausdorff_distance(mio_o_mesh(a[0]), mio_o_mesh(a[1]), &opts, &r),
                "hausdorff_distance");
    const char *names[] = {"distance",      "a_to_b",        "b_to_a",       "mean_a_to_b",
                           "rms_a_to_b",    "mean_b_to_a",   "rms_b_to_a",   "num_samples_a",
                           "num_samples_b", "worst_point_a", "worst_point_b"};
    mxArray *values[] = {mio_o_scalar(r.distance),
                         mio_o_scalar(r.a_to_b),
                         mio_o_scalar(r.b_to_a),
                         mio_o_scalar(r.mean_a_to_b),
                         mio_o_scalar(r.rms_a_to_b),
                         mio_o_scalar(r.mean_b_to_a),
                         mio_o_scalar(r.rms_b_to_a),
                         mio_o_scalar((double)r.num_samples_a),
                         mio_o_scalar((double)r.num_samples_b),
                         mio_o_row3(r.worst_point_a),
                         mio_o_row3(r.worst_point_b)};
    return mio_o_struct(11, names, values);
}

/* repair(mesh, fix_orientation, orient_outward, fill_holes, split_non_manifold,
 *        max_hole_edges, weld_tolerance, record_provenance, fill_style) */
static mxArray *O_mio_repair(const mxArray **a) {
    mio_repair_opts opts;
    mio_repair_report report;

    mio_repair_opts_init(&opts);
    opts.fix_orientation = mio_o_bool(a[1], "fix_orientation") ? 1 : 0;
    opts.orient_outward = mio_o_bool(a[2], "orient_outward") ? 1 : 0;
    opts.fill_holes = mio_o_bool(a[3], "fill_holes") ? 1 : 0;
    opts.split_non_manifold = mio_o_bool(a[4], "split_non_manifold") ? 1 : 0;
    opts.max_hole_edges = (int64_t)mio_o_double(a[5], "max_hole_edges");
    opts.weld_tolerance = mio_o_double(a[6], "weld_tolerance");
    opts.record_provenance = mio_o_bool(a[7], "record_provenance") ? 1 : 0;
    opts.fill_style = mio_o_int(a[8], "fill_style"); /* roadmap 7.8.6 */

    mio_mesh *out = mio_repair(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("repair");
    const int64_t counters[] = {report.num_flipped,        report.num_components,
                                report.largest_component,  report.num_oriented_outward,
                                report.num_unorientable,   report.num_vertices_split,
                                report.num_holes_detected, report.num_holes_filled,
                                report.num_holes_skipped,  report.num_faces_added,
                                report.num_points_added,   report.points_welded};
    const char *names[] = {"mesh",
                           "quality_before",
                           "quality_after",
                           "num_flipped",
                           "num_components",
                           "largest_component",
                           "num_oriented_outward",
                           "num_unorientable",
                           "num_vertices_split",
                           "num_holes_detected",
                           "num_holes_filled",
                           "num_holes_skipped",
                           "num_faces_added",
                           "num_points_added",
                           "points_welded"};
    mxArray *values[15];
    values[0] = mio_o_wrap_mesh(out);
    values[1] = mio_o_surface_quality(&report.quality_before);
    values[2] = mio_o_surface_quality(&report.quality_after);
    for (int k = 0; k < 12; ++k) values[3 + k] = mio_o_scalar((double)counters[k]);
    return mio_o_struct(15, names, values);
}

/* shrinkwrap(mesh, target, offset, max_distance, weights, target_region,
 *            normal_weight, record_distance, record_closest_cell) */
static mxArray *O_mio_shrinkwrap(const mxArray **a) {
    mio_shrinkwrap_opts opts;
    mio_shrinkwrap_report report;

    mio_shrinkwrap_opts_init(&opts);
    opts.offset = mio_o_double(a[2], "offset");
    opts.max_distance = mio_o_double(a[3], "max_distance");
    opts.weights = mio_o_opt_string(a[4]);
    opts.target_region = mio_o_opt_string(a[5]);
    const char *nw = mio_o_opt_string(a[6]);
    if (nw != NULL && strcmp(nw, "area") == 0) {
        opts.normal_weight = MIO_SDF_WEIGHT_AREA;
    } else if (nw != NULL && strcmp(nw, "angle") != 0) {
        mio_o_error("shrinkwrap",
                    "meshio++: shrinkwrap: unknown normal weight '%s' (expected 'angle' or "
                    "'area')",
                    nw);
    }
    opts.record_distance = mio_o_bool(a[7], "record_distance") ? 1 : 0;
    opts.record_closest_cell = mio_o_bool(a[8], "record_closest_cell") ? 1 : 0;

    mio_mesh *out = mio_shrinkwrap(mio_o_mesh(a[0]), mio_o_mesh(a[1]), &opts, &report);
    if (out == NULL) mio_o_fail("shrinkwrap");
    const char *names[] = {"mesh",       "quality",     "num_projected",
                           "num_missed", "num_skipped", "max_displacement"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_surface_quality(&report.quality),
                         mio_o_scalar((double)report.num_projected),
                         mio_o_scalar((double)report.num_missed),
                         mio_o_scalar((double)report.num_skipped),
                         mio_o_scalar(report.max_displacement)};
    return mio_o_struct(6, names, values);
}

/* sobolev_deform(mesh, array, length_scale, fixed_points_array, fix_boundary,
 *                record_filtered, max_iterations, tolerance) */
static mxArray *O_mio_sobolev_deform(const mxArray **a) {
    mio_sobolev_opts opts;
    mio_sobolev_report report;

    mio_sobolev_opts_init(&opts);
    opts.array = mio_o_opt_string(a[1]);
    opts.length_scale = mio_o_double(a[2], "length_scale");
    opts.fixed_points_array = mio_o_opt_string(a[3]);
    opts.fix_boundary = mio_o_bool(a[4], "fix_boundary") ? 1 : 0;
    opts.record_filtered = mio_o_bool(a[5], "record_filtered") ? 1 : 0;
    opts.max_iterations = mio_o_int(a[6], "max_iterations");
    opts.tolerance = mio_o_double(a[7], "tolerance");

    mio_mesh *out = mio_sobolev_deform(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("sobolev_deform");
    const char *names[] = {"mesh",      "num_iterations", "residual",        "converged",
                           "num_fixed", "num_isolated",   "max_displacement"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)report.num_iterations),
                         mio_o_scalar(report.residual),
                         mio_o_logical(report.converged != 0),
                         mio_o_scalar((double)report.num_fixed),
                         mio_o_scalar((double)report.num_isolated),
                         mio_o_scalar(report.max_displacement)};
    return mio_o_struct(7, names, values);
}

/* --- metrics, untangling, morphing, spectra ------------------------------ */

/* prepare_metric(mesh, arrays, output, target_vertices, norm, hmin, hmax,
 *                gradation) */
static mxArray *O_mio_prepare_metric(const mxArray **a) {
    mio_metric_opts opts;
    mio_metric_opts_init(&opts);
    int64_t n = 0;
    const char *const *names = mio_o_names(a[1], &n, "arrays");
    opts.arrays = names;
    opts.num_arrays = n;
    opts.output = mio_o_string(a[2], "output");
    opts.target_vertices = mio_o_double(a[3], "target_vertices");
    opts.norm = mio_o_double(a[4], "norm");
    opts.hmin = mio_o_double(a[5], "hmin");
    opts.hmax = mio_o_double(a[6], "hmax");
    opts.gradation = mio_o_double(a[7], "gradation");
    mio_metric_report rep;
    mio_mesh *out = mio_prepare_metric(mio_o_mesh(a[0]), &opts, &rep);
    if (out == NULL) mio_o_fail("prepare_metric");
    const char *keys[] = {"mesh", "complexity_before", "complexity_after", "num_graded"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar(rep.complexity_before),
                         mio_o_scalar(rep.complexity_after), mio_o_scalar((double)rep.num_graded)};
    return mio_o_struct(4, keys, values);
}

/* untangle(mesh, fix_boundary, rings, smooth, max_sweeps) */
static mxArray *O_mio_untangle(const mxArray **a) {
    mio_untangle_opts opts;
    mio_untangle_opts_init(&opts);
    opts.fix_boundary = mio_o_bool(a[1], "fix_boundary") ? 1 : 0;
    opts.rings = (int64_t)mio_o_double(a[2], "rings");
    opts.smooth = mio_o_bool(a[3], "smooth") ? 1 : 0;
    opts.max_sweeps = (int64_t)mio_o_double(a[4], "max_sweeps");
    mio_untangle_report rep;
    mio_mesh *out = mio_untangle(mio_o_mesh(a[0]), &opts, &rep);
    if (out == NULL) mio_o_fail("untangle");
    const char *keys[] = {"mesh",          "num_inverted_before", "num_inverted_after",
                          "min_quality_before", "min_quality_after", "sweeps"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)rep.num_inverted_before),
                         mio_o_scalar((double)rep.num_inverted_after),
                         mio_o_scalar(rep.min_quality_before),
                         mio_o_scalar(rep.min_quality_after),
                         mio_o_scalar((double)rep.sweeps)};
    return mio_o_struct(6, keys, values);
}

/* morph(mesh, array, fix_boundary, constrained_array, stiffening, steps,
 *       allow_inverted) */
static mxArray *O_mio_morph(const mxArray **a) {
    mio_morph_opts opts;
    mio_morph_opts_init(&opts);
    opts.array = mio_o_string(a[1], "array");
    if (!mio_o_is_null(a[3])) opts.constrained_array = mio_o_string(a[3], "constrained_array");
    opts.fix_boundary = mio_o_bool(a[2], "fix_boundary") ? 1 : 0;
    opts.stiffening = mio_o_double(a[4], "stiffening");
    opts.steps = (int64_t)mio_o_double(a[5], "steps");
    opts.allow_inverted = mio_o_bool(a[6], "allow_inverted") ? 1 : 0;
    mio_morph_report rep;
    mio_mesh *out = mio_morph(mio_o_mesh(a[0]), &opts, &rep);
    if (out == NULL) mio_o_fail("morph");
    const char *keys[] = {"mesh",         "num_constrained",  "num_free",
                          "num_inverted", "max_displacement", "converged"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)rep.num_constrained),
                         mio_o_scalar((double)rep.num_free),
                         mio_o_scalar((double)rep.num_inverted),
                         mio_o_scalar(rep.max_displacement),
                         mio_o_logical(rep.converged != 0)};
    return mio_o_struct(6, keys, values);
}

/* laplacian_eigenmodes(mesh, num_modes, tolerance, max_iterations) */
static mxArray *O_mio_laplacian_eigenmodes(const mxArray **a) {
    mio_eigenmode_opts opts;
    mio_eigenmode_opts_init(&opts);
    opts.num_modes = (int64_t)mio_o_double(a[1], "num_modes");
    opts.tolerance = mio_o_double(a[2], "tolerance");
    opts.max_iterations = (int64_t)mio_o_double(a[3], "max_iterations");
    mio_eigenmode_report rep;
    mio_mesh *out = mio_laplacian_eigenmodes(mio_o_mesh(a[0]), &opts, &rep);
    if (out == NULL) mio_o_fail("laplacian_eigenmodes");
    const char *keys[] = {"mesh", "iterations", "max_residual", "converged"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)rep.iterations),
                         mio_o_scalar(rep.max_residual), mio_o_logical(rep.converged != 0)};
    return mio_o_struct(4, keys, values);
}

/* --- registration, point clouds, hulls ----------------------------------- */

/* register_rigid(source, target, method, max_iterations, tolerance,
 *                max_distance, initial, initial_mode, symmetric_check,
 *                source_region, target_region, adaptive_trim) */
static mxArray *O_mio_register_rigid(const mxArray **a) {
    mio_registration_ex_opts opts;
    mio_registration_ex_opts_init(&opts);
    const char *m = mio_o_string(a[2], "method");
    if (strcmp(m, "point_to_plane") == 0) {
        opts.base.method = MIO_REGISTER_POINT_TO_PLANE;
    } else if (strcmp(m, "point_to_point") != 0) {
        mio_o_error("register_rigid", "`method` must be \"point_to_point\" or \"point_to_plane\"");
    }
    opts.base.max_iterations = (int64_t)mio_o_double(a[3], "max_iterations");
    opts.base.tolerance = mio_o_double(a[4], "tolerance");
    opts.base.max_distance = mio_o_double(a[5], "max_distance");
    if (!mio_o_is_null(a[6])) {
        double x[16];
        oc_vec_exact(a[6], x, 16, "register_rigid", "`initial` must be a 4 x 4 matrix");
        /* Octave stores column-major: element (i, j) at j * 4 + i */
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) opts.base.initial[4 * i + j] = x[j * 4 + i];
    }
    const char *im = mio_o_string(a[7], "initial");
    if (strcmp(im, "auto") == 0) {
        opts.initial_mode = MIO_REGISTER_INITIAL_AUTO;
    } else if (strcmp(im, "matrix") != 0) {
        mio_o_error("register_rigid", "`initial` must be a 4 x 4 matrix, [] or \"auto\"");
    }
    opts.symmetric_check = mio_o_bool(a[8], "symmetric_check");
    opts.source_region = mio_o_opt_string(a[9]);
    opts.target_region = mio_o_opt_string(a[10]);
    opts.adaptive_trim = mio_o_double(a[11], "adaptive_trim");
    mio_registration_report rep;
    if (mio_register_rigid_ex(mio_o_mesh(a[0]), mio_o_mesh(a[1]), &opts, &rep) != MIO_OK)
        mio_o_fail("register_rigid");
    mxArray *mat = mio_o_matrix(4, 4);
    double *dst = mxGetPr(mat);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) dst[j * 4 + i] = rep.matrix[4 * i + j];
    mxArray *cands = mio_o_column(rep.num_candidates);
    if (rep.num_candidates > 0)
        mio_register_rigid_candidates(mxGetPr(cands), rep.num_candidates);
    mxArray *w = mio_o_getbuffer(mio_register_rigid_warnings, "register_rigid");
    const char *keys[] = {"matrix",       "initial_rms",      "rms",          "iterations",
                          "num_pairs",    "converged",        "num_rejected", "chosen_candidate",
                          "candidate_rms", "warnings"};
    mxArray *values[] = {mat,
                         mio_o_scalar(rep.initial_rms),
                         mio_o_scalar(rep.rms),
                         mio_o_scalar((double)rep.iterations),
                         mio_o_scalar((double)rep.num_pairs),
                         mio_o_logical(rep.converged != 0),
                         mio_o_scalar((double)rep.num_rejected),
                         mio_o_scalar((double)rep.chosen_candidate),
                         cands,
                         w};
    return mio_o_struct(10, keys, values);
}

/* register_landmarks(source, target, weights): n x 3 matrices, weights [] or n. */
static mxArray *O_mio_register_landmarks(const mxArray **a) {
    int64_t ns = 0, cs = 0, nt = 0, ct = 0;
    double *src = mio_o_rows_to_c(a[0], &ns, &cs, 0.0, "source");
    double *dst = mio_o_rows_to_c(a[1], &nt, &ct, 0.0, "target");
    if (cs != 3 || ct != 3 || ns != nt)
        mio_o_error("register_landmarks", "`source` and `target` must be n x 3 matrices");
    const double *w = NULL;
    if (!mio_o_is_null(a[2])) {
        size_t nw = 0;
        w = mio_o_doubles(a[2], &nw, "weights");
        if ((int64_t)nw != ns)
            mio_o_error("register_landmarks", "`weights` must have one entry per pair");
    }
    mxArray *res = mio_o_column(ns);
    double m[16], rms = 0.0;
    if (mio_register_landmarks(src, dst, w, ns, m, ns > 0 ? mxGetPr(res) : NULL, &rms) != MIO_OK)
        mio_o_fail("register_landmarks");
    mxArray *mat = mio_o_matrix(4, 4);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) mxGetPr(mat)[j * 4 + i] = m[4 * i + j];
    const char *keys[] = {"matrix", "residuals", "rms"};
    mxArray *values[] = {mat, res, mio_o_scalar(rms)};
    return mio_o_struct(3, keys, values);
}

/* find_symmetry_plane(mesh, seed_origin, seed_normal, tolerance, min_share,
 *                     max_samples) */
static mxArray *O_mio_find_symmetry_plane(const mxArray **a) {
    mio_symmetry_plane_opts opts;
    mio_symmetry_plane_opts_init(&opts);
    if (!mio_o_is_null(a[1]) || !mio_o_is_null(a[2])) {
        mio_o_vecn(a[1], opts.seed_origin, 3, "seed origin");
        mio_o_vecn(a[2], opts.seed_normal, 3, "seed normal");
    }
    opts.tolerance = mio_o_double(a[3], "tolerance");
    opts.min_share = mio_o_double(a[4], "min_share");
    opts.max_samples = (int64_t)mio_o_double(a[5], "max_samples");
    mio_symmetry_plane_report r;
    if (mio_find_symmetry_plane(mio_o_mesh(a[0]), &opts, &r) != MIO_OK)
        mio_o_fail("find_symmetry_plane");
    mxArray *o = mio_o_matrix(1, 3);
    mxArray *nn = mio_o_matrix(1, 3);
    for (int k = 0; k < 3; ++k) {
        mxGetPr(o)[k] = r.origin[k];
        mxGetPr(nn)[k] = r.normal[k];
    }
    const char *keys[] = {"origin", "normal", "rms", "share", "tolerance", "chosen_candidate"};
    mxArray *values[] = {o,
                         nn,
                         mio_o_scalar(r.rms),
                         mio_o_scalar(r.share),
                         mio_o_scalar(r.tolerance),
                         mio_o_scalar((double)r.chosen_candidate)};
    return mio_o_struct(6, keys, values);
}

/* auto_align(mesh, angle_tolerance, subdivisions, region, rotate_vector_data) */
static mxArray *O_mio_auto_align(const mxArray **a) {
    mio_auto_align_opts opts;
    mio_auto_align_opts_init(&opts);
    opts.angle_tolerance = mio_o_double(a[1], "angle_tolerance");
    opts.subdivisions = (int64_t)mio_o_double(a[2], "subdivisions");
    opts.region = mio_o_opt_string(a[3]);
    opts.rotate_vector_data = mio_o_bool(a[4], "rotate_vector_data");
    mio_auto_align_report r;
    mio_mesh *out = mio_auto_align(mio_o_mesh(a[0]), &opts, &r);
    if (out == NULL) mio_o_fail("auto_align");
    mxArray *mat = mio_o_matrix(4, 4);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) mxGetPr(mat)[j * 4 + i] = r.matrix[4 * i + j];
    mxArray *pr = mio_o_matrix(1, 3);
    mxArray *se = mio_o_matrix(1, 3);
    for (int k = 0; k < 3; ++k) {
        mxGetPr(pr)[k] = r.primary[k];
        mxGetPr(se)[k] = r.secondary[k];
    }
    const char *keys[] = {"mesh",  "matrix",     "mode",       "primary",     "secondary",
                          "share", "face_share", "axis_share", "opening_area"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mat,
                         mio_o_str(r.mode == MIO_AUTO_ALIGN_AXIS ? "axis" : "face"),
                         pr,
                         se,
                         mio_o_scalar(r.share),
                         mio_o_scalar(r.face_share),
                         mio_o_scalar(r.axis_share),
                         mio_o_scalar(r.opening_area)};
    return mio_o_struct(9, keys, values);
}

/* estimate_normals(mesh, k, orientation, viewpoint, record_variation) */
static mxArray *O_mio_estimate_normals(const mxArray **a) {
    mio_point_normals_opts opts;
    mio_point_normals_opts_init(&opts);
    opts.k = (int64_t)mio_o_double(a[1], "k");
    const char *o = mio_o_string(a[2], "orientation");
    if (strcmp(o, "consistent") == 0) {
        opts.orientation = MIO_NORMALS_CONSISTENT;
    } else if (strcmp(o, "viewpoint") == 0) {
        opts.orientation = MIO_NORMALS_VIEWPOINT;
    } else if (strcmp(o, "none") == 0) {
        opts.orientation = MIO_NORMALS_NONE;
    } else {
        mio_o_error("estimate_normals",
                    "`orientation` must be \"consistent\", \"viewpoint\" or \"none\"");
    }
    oc_vec_exact(a[3], opts.viewpoint, 3, "estimate_normals", "`viewpoint` takes three numbers");
    opts.record_variation = mio_o_bool(a[4], "record_variation") ? 1 : 0;
    int64_t nd = 0, nc = 0;
    mio_mesh *out = mio_estimate_normals(mio_o_mesh(a[0]), &opts, &nd, &nc);
    if (out == NULL) mio_o_fail("estimate_normals");
    const char *keys[] = {"mesh", "num_degenerate", "num_components"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)nd), mio_o_scalar((double)nc)};
    return mio_o_struct(3, keys, values);
}

/* reconstruct_surface(mesh, normals, k, cell_size, resolution, max_distance,
 *                     estimate_normals, normals_k, max_cells) */
static mxArray *O_mio_reconstruct_surface(const mxArray **a) {
    mio_reconstruct_surface_opts opts;
    mio_reconstruct_surface_opts_init(&opts);
    opts.normals = mio_o_string(a[1], "normals");
    opts.k = (int64_t)mio_o_double(a[2], "k");
    opts.cell_size = mio_o_double(a[3], "cell_size");
    opts.resolution = (int64_t)mio_o_double(a[4], "resolution");
    opts.max_distance = mio_o_double(a[5], "max_distance");
    opts.estimate_normals = mio_o_bool(a[6], "estimate_normals") ? 1 : 0;
    opts.normals_k = (int64_t)mio_o_double(a[7], "normals_k");
    opts.max_cells = (int64_t)mio_o_double(a[8], "max_cells");
    mio_reconstruct_surface_report rep;
    mio_mesh *out = mio_reconstruct_surface(mio_o_mesh(a[0]), &opts, &rep);
    if (out == NULL) mio_o_fail("reconstruct_surface");
    const char *keys[] = {"mesh",          "num_samples",   "spacing",
                          "cell_size",     "max_distance",  "num_lattice_points",
                          "num_undefined", "num_triangles", "num_components",
                          "estimated_normals"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)rep.num_samples),
                         mio_o_scalar(rep.spacing),
                         mio_o_scalar(rep.cell_size),
                         mio_o_scalar(rep.max_distance),
                         mio_o_scalar((double)rep.num_lattice_points),
                         mio_o_scalar((double)rep.num_undefined),
                         mio_o_scalar((double)rep.num_triangles),
                         mio_o_scalar((double)rep.num_components),
                         mio_o_logical(rep.estimated_normals != 0)};
    return mio_o_struct(10, keys, values);
}

/* convex_hull(mesh, dimension) */
static mxArray *O_mio_convex_hull(const mxArray **a) {
    const int d = mio_o_int(a[1], "dimension");
    if (d != 0 && d != 2 && d != 3) mio_o_error("convex_hull", "`dimension` must be 0, 2 or 3");
    mio_hull_report rep;
    mio_mesh *out = mio_convex_hull(mio_o_mesh(a[0]), (int32_t)d, &rep);
    if (out == NULL) mio_o_fail("convex_hull");
    const char *keys[] = {"mesh", "dimension", "num_vertices", "num_facets",
                          "area", "volume",    "perimeter"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)rep.dimension),
                         mio_o_scalar((double)rep.num_vertices),
                         mio_o_scalar((double)rep.num_facets),
                         mio_o_scalar(rep.area),
                         mio_o_scalar(rep.volume),
                         mio_o_scalar(rep.perimeter)};
    return mio_o_struct(7, keys, values);
}

static mxArray *O_mio_delaunay_2d(const mxArray **a) {
    int64_t nt = 0, nd = 0, nb = 0;
    mio_mesh *out = mio_delaunay_2d(mio_o_mesh(a[0]), &nt, &nd, &nb);
    if (out == NULL) mio_o_fail("delaunay_2d");
    const char *keys[] = {"mesh", "num_triangles", "num_duplicates", "num_boundary"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)nt), mio_o_scalar((double)nd),
                         mio_o_scalar((double)nb)};
    return mio_o_struct(4, keys, values);
}

/* convert_coordinates(mesh, system, axis, origin, reference, arrays, inverse,
 *                     record_coordinates) */
static mxArray *O_mio_convert_coordinates(const mxArray **a) {
    mio_coordinate_opts opts;
    mio_coordinate_opts_init(&opts);
    const char *sys = mio_o_string(a[1], "system");
    if (strcmp(sys, "spherical") == 0) {
        opts.system = MIO_COORDS_SPHERICAL;
    } else if (strcmp(sys, "cylindrical") != 0) {
        mio_o_error("convert_coordinates", "`system` must be \"cylindrical\" or \"spherical\"");
    }
    double *dst[3] = {opts.axis, opts.origin, opts.reference};
    for (int v = 0; v < 3; ++v) {
        oc_vec_exact(a[2 + v], dst[v], 3, "convert_coordinates",
                     "`axis`, `origin` and `reference` take three numbers");
    }
    int64_t n = 0;
    const char *const *names = mio_o_names(a[5], &n, "arrays");
    opts.arrays = names;
    opts.num_arrays = n;
    opts.inverse = mio_o_bool(a[6], "inverse") ? 1 : 0;
    opts.record_coordinates = mio_o_bool(a[7], "record_coordinates") ? 1 : 0;
    int64_t na = 0, ns = 0;
    mio_mesh *out = mio_convert_coordinates(mio_o_mesh(a[0]), &opts, &na, &ns);
    if (out == NULL) mio_o_fail("convert_coordinates");
    const char *keys[] = {"mesh", "num_arrays", "num_singular"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)na), mio_o_scalar((double)ns)};
    return mio_o_struct(3, keys, values);
}

/* geodesic_distance(mesh, sources, source_region, time_factor, tolerance,
 *                   max_iterations) */
static mxArray *O_mio_geodesic_distance(const mxArray **a) {
    mio_geodesic_opts opts;
    mio_geodesic_report report;
    mio_geodesic_opts_init(&opts);
    int64_t n = 0;
    opts.sources = mio_o_ids_0based(a[1], &n, "sources"); /* 1-based here */
    opts.num_sources = n;
    opts.source_region = mio_o_opt_string(a[2]);
    opts.time_factor = mio_o_double(a[3], "time_factor");
    opts.tolerance = mio_o_double(a[4], "tolerance");
    opts.max_iterations = (int64_t)mio_o_double(a[5], "max_iterations");
    mio_mesh *out = mio_geodesic_distance(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("geodesic_distance");
    const char *names[] = {"mesh",         "num_sources", "num_unreached",
                           "max_distance", "time_step",   "converged"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)report.num_sources),
                         mio_o_scalar((double)report.num_unreached),
                         mio_o_scalar(report.max_distance),
                         mio_o_scalar(report.time_step),
                         mio_o_logical(report.converged != 0)};
    return mio_o_struct(6, names, values);
}

/* parameterize(mesh, padding, max_depth, max_stretch) */
static mxArray *O_mio_parameterize(const mxArray **a) {
    mio_parameterize_opts opts;
    mio_parameterize_report report;
    mio_parameterize_opts_init(&opts);
    opts.padding = mio_o_double(a[1], "padding");
    opts.max_depth = (int64_t)mio_o_double(a[2], "max_depth");
    opts.max_stretch = mio_o_double(a[3], "max_stretch");
    mio_mesh *out = mio_parameterize(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("parameterize");
    const char *names[] = {"mesh",    "num_charts",         "num_flipped",
                           "stretch", "packing_efficiency", "num_seam_points"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)report.num_charts),
                         mio_o_scalar((double)report.num_flipped),
                         mio_o_scalar(report.stretch),
                         mio_o_scalar(report.packing_efficiency),
                         mio_o_scalar((double)report.num_seam_points)};
    return mio_o_struct(6, names, values);
}

/* --- table --------------------------------------------------------------- */

static const mio_o_entry kEntries[] = {
    MIO_O_ENTRY(mio_extract_surface, 2),
    MIO_O_ENTRY(mio_extract_skin, 2),
    MIO_O_ENTRY(mio_attach_quality, 1),
    MIO_O_ENTRY(mio_quality_counts, 1),
    MIO_O_ENTRY(mio_transform, 3),
    MIO_O_ENTRY(mio_clean, 6),
    MIO_O_ENTRY(mio_smooth, 10),
    MIO_O_ENTRY(mio_optimize_volume, 8),
    MIO_O_ENTRY(mio_adapt_volume, 17),
    MIO_O_ENTRY(mio_crop_bbox, 5),
    MIO_O_ENTRY(mio_crop_predicate, 5),
    MIO_O_ENTRY(mio_crop_plane, 5),
    MIO_O_ENTRY(mio_slice, 6),
    MIO_O_ENTRY(mio_sections, 5),
    MIO_O_ENTRY(mio_fit_curves, 6),
    MIO_O_ENTRY(mio_silhouette, 8),
    MIO_O_ENTRY(mio_write_step_elements, 4),
    MIO_O_ENTRY(mio_measure, 2),
    MIO_O_ENTRY(mio_remove_islands, 7),
    MIO_O_ENTRY(mio_flip_normals, 2),
    MIO_O_ENTRY(mio_clip, 5),
    MIO_O_ENTRY(mio_offset_surface, 3),
    MIO_O_ENTRY(mio_constrain_primitives, 7),
    MIO_O_ENTRY(mio_profile_solid, 6),
    MIO_O_ENTRY(mio_face_support, 6),
    MIO_O_ENTRY(mio_fit_cad_iterate, 11),
    MIO_O_ENTRY(mio_thicken, 4),
    MIO_O_ENTRY(mio_hollow, 3),
    MIO_O_ENTRY(mio_align_datums, 4),
    MIO_O_ENTRY(mio_measure_element_primitive, 2),
    MIO_O_ENTRY(mio_isosurface, 5),
    MIO_O_ENTRY(mio_gradient, 8),
    MIO_O_ENTRY(mio_hessian, 6),
    MIO_O_ENTRY(mio_estimate_error, 8),
    MIO_O_ENTRY(mio_remesh, 11),
    MIO_O_ENTRY(mio_compute_curvature, 8),
    MIO_O_ENTRY(mio_compute_normals, 7),
    MIO_O_ENTRY(mio_feature_edges, 7),
    MIO_O_ENTRY(mio_segment, 9),
    MIO_O_ENTRY(mio_fit_primitives, 14),
    MIO_O_ENTRY(mio_has_occt, 0),
    MIO_O_ENTRY(mio_refine_patches, 13),
    MIO_O_ENTRY(mio_fit_cad, 48),
    MIO_O_ENTRY(mio_region_adjacency, 2),
    MIO_O_ENTRY(mio_find_interface, 9),
    MIO_O_ENTRY(mio_contact_pairs, 6),
    MIO_O_ENTRY(mio_split_interface, 3),
    MIO_O_ENTRY(mio_check_quality, 4),
    MIO_O_ENTRY(mio_hausdorff_distance, 5),
    MIO_O_ENTRY(mio_repair, 9),
    MIO_O_ENTRY(mio_shrinkwrap, 9),
    MIO_O_ENTRY(mio_sobolev_deform, 8),
    MIO_O_ENTRY(mio_prepare_metric, 8),
    MIO_O_ENTRY(mio_untangle, 5),
    MIO_O_ENTRY(mio_morph, 7),
    MIO_O_ENTRY(mio_laplacian_eigenmodes, 4),
    MIO_O_ENTRY(mio_register_rigid, 12),
    MIO_O_ENTRY(mio_register_landmarks, 3),
    MIO_O_ENTRY(mio_find_symmetry_plane, 6),
    MIO_O_ENTRY(mio_auto_align, 5),
    MIO_O_ENTRY(mio_estimate_normals, 5),
    MIO_O_ENTRY(mio_reconstruct_surface, 9),
    MIO_O_ENTRY(mio_convex_hull, 2),
    MIO_O_ENTRY(mio_delaunay_2d, 1),
    MIO_O_ENTRY(mio_convert_coordinates, 8),
    MIO_O_ENTRY(mio_geodesic_distance, 6),
    MIO_O_ENTRY(mio_parameterize, 4),
};

const mio_o_entry *mio_o_table_ops_core(size_t *n) {
    *n = sizeof(kEntries) / sizeof(kEntries[0]);
    return kEntries;
}
