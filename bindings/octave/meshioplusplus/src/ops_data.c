/* Data operations, the settings pipeline, regular grids and the parametric
 * generators, voxelize / octree / signed distance, remesh_volume, provenance
 * and select. A port of the last part of the R binding's src/ops.c (from
 * R_mio_sets_to_data on); see mio_oct.h for the layout rules. */

#include "mio_oct.h"

#include <stdio.h>
#include <string.h>

/* --- small helpers ------------------------------------------------------- */

/* A numeric argument of exactly `n` values into `out`, or 0 when it is [] /
 * empty -- the R binding's `numeric(0)` for "not given". */
static int od_opt_vecn(const mxArray *a, double *out, size_t n, const char *what) {
    if (a == NULL || mxGetNumberOfElements(a) == 0) return 0;
    mio_o_vecn(a, out, n, what);
    return 1;
}

/* Raise with the C library's current message, copied first so freeing a native
 * object in between cannot disturb it. Frees `report` (may be NULL). */
static void od_fail_report(mio_pipeline_report *report, const char *what) {
    char msg[1024];
    const char *e = mio_last_error();
    snprintf(msg, sizeof(msg), "%s",
             (e != NULL && e[0] != '\0') ? e : "meshio++: pipeline report failed");
    if (report != NULL) mio_pipeline_report_free(report);
    mio_o_error(what, "%s", msg);
}

/* --- data operations ----------------------------------------------------- */

/* sets_to_data(mesh, location, data_name, join_char, order) */
static mxArray *O_mio_sets_to_data(const mxArray **a) {
    int64_t count = 0;
    const char *const *nm = mio_o_names(a[4], &count, "order");
    mio_mesh *out =
        mio_sets_to_data(mio_o_mesh(a[0]), (mio_data_location)mio_o_int(a[1], "location"),
                         mio_o_opt_string(a[2]), mio_o_opt_string(a[3]), nm, count);
    if (out == NULL) mio_o_fail("sets_to_data");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_data_to_sets(const mxArray **a) {
    mio_mesh *out =
        mio_data_to_sets(mio_o_mesh(a[0]), (mio_data_location)mio_o_int(a[1], "location"),
                         mio_o_string(a[2], "key"));
    if (out == NULL) mio_o_fail("data_to_sets");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_data_drop(const mxArray **a) {
    int64_t count = 0;
    const char *const *nm = mio_o_names(a[2], &count, "names");
    mio_mesh *out = mio_data_drop(mio_o_mesh(a[0]),
                                  (mio_data_location)mio_o_int(a[1], "location"), nm, count,
                                  mio_o_bool(a[3], "ignore_missing"));
    if (out == NULL) mio_o_fail("data_drop");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_data_keep(const mxArray **a) {
    int64_t count = 0;
    const char *const *nm = mio_o_names(a[2], &count, "names");
    mio_mesh *out = mio_data_keep(mio_o_mesh(a[0]),
                                  (mio_data_location)mio_o_int(a[1], "location"), nm, count,
                                  mio_o_bool(a[3], "ignore_missing"));
    if (out == NULL) mio_o_fail("data_keep");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_data_rename(const mxArray **a) {
    mio_mesh *out = mio_data_rename(mio_o_mesh(a[0]),
                                    (mio_data_location)mio_o_int(a[1], "location"),
                                    mio_o_string(a[2], "from"), mio_o_string(a[3], "to"));
    if (out == NULL) mio_o_fail("data_rename");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_data_point_to_cell(const mxArray **a) {
    int64_t count = 0;
    const char *const *nm = mio_o_names(a[1], &count, "names");
    mio_mesh *out = mio_data_point_to_cell(mio_o_mesh(a[0]), nm, count, mio_o_opt_string(a[2]));
    if (out == NULL) mio_o_fail("data_point_to_cell");
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_data_cell_to_point(const mxArray **a) {
    int64_t count = 0;
    const char *const *nm = mio_o_names(a[1], &count, "names");
    mio_mesh *out = mio_data_cell_to_point(
        mio_o_mesh(a[0]), nm, count, (mio_cell_point_weight)mio_o_int(a[2], "weight"),
        mio_o_opt_string(a[3]));
    if (out == NULL) mio_o_fail("data_cell_to_point");
    return mio_o_wrap_mesh(out);
}

/* data_calc(mesh, expression, location, output_name, overwrite) */
static mxArray *O_mio_data_calc(const mxArray **a) {
    mio_mesh *out = mio_data_calc(mio_o_mesh(a[0]), mio_o_string(a[1], "expression"),
                                  (mio_data_location)mio_o_int(a[2], "location"),
                                  mio_o_string(a[3], "output_name"),
                                  mio_o_bool(a[4], "overwrite"));
    if (out == NULL) mio_o_fail("data_calc");
    return mio_o_wrap_mesh(out);
}

/* data_condition(mesh, location, names, mode, lo, hi, scope, nan_policy,
 *                nan_replacement, suffix) */
static mxArray *O_mio_data_condition(const mxArray **a) {
    int64_t count = 0;
    const char *const *nm = mio_o_names(a[2], &count, "names");
    mio_mesh *out = mio_data_condition(
        mio_o_mesh(a[0]), (mio_data_location)mio_o_int(a[1], "location"), nm, count,
        (mio_condition_mode)mio_o_int(a[3], "mode"), mio_o_double(a[4], "lo"),
        mio_o_double(a[5], "hi"), (mio_condition_scope)mio_o_int(a[6], "scope"),
        (mio_nan_policy)mio_o_int(a[7], "nan_policy"),
        mio_o_double(a[8], "nan_replacement"), mio_o_opt_string(a[9]));
    if (out == NULL) mio_o_fail("data_condition");
    return mio_o_wrap_mesh(out);
}

/* tensor_invariants(mesh, location, names, outputs, prefix, suffix, overwrite) */
static mxArray *O_mio_tensor_invariants(const mxArray **a) {
    int64_t count = 0;
    const char *const *nm = mio_o_names(a[2], &count, "names");
    mio_mesh *out = mio_tensor_invariants(
        mio_o_mesh(a[0]), (mio_data_location)mio_o_int(a[1], "location"), nm, count,
        mio_o_int(a[3], "outputs"), mio_o_opt_string(a[4]), mio_o_opt_string(a[5]),
        mio_o_bool(a[6], "overwrite"));
    if (out == NULL) mio_o_fail("tensor_invariants");
    return mio_o_wrap_mesh(out);
}

typedef struct {
    const mio_data_info *info;
    int64_t index;
} od_info_ctx;

static int64_t od_info_name_getter(void *c, char *buf, int64_t buflen) {
    od_info_ctx *x = (od_info_ctx *)c;
    return mio_data_info_name(x->info, x->index, buf, buflen);
}

static const char *kDataInfoFields[] = {
    "name",        "location",       "dtype",      "num_blocks", "num_entries",
    "num_components", "num_values",  "min",        "max",        "mean",
    "num_nan",     "num_inf",        "num_finite", "inconsistent_blocks", "components"};

/* One summary per array, as an n x 1 struct array. `components` is
 * (num_components x 3): columns min, max, mean -- R's matrix with those
 * column names. */
static mxArray *O_mio_data_info(const mxArray **a) {
    mio_data_info *info =
        (mio_data_info *)mio_o_check_ptr(mio_data_info_create(mio_o_mesh(a[0])), "data_info");
    int64_t n = mio_data_info_count(info);
    if (n < 0) {
        mio_data_info_free(info);
        mio_o_fail("data_info count");
    }
    static const char *const kLocations[] = {"point", "cell", "field"};

    mxArray *out = mio_o_struct_array(n, 15, kDataInfoFields);
    for (int64_t i = 0; i < n; ++i) {
        /* Read the name through the C protocol directly: a getter failure
         * would otherwise leak `info`. */
        int64_t len = mio_data_info_name(info, i, NULL, 0);
        mio_data_array_info e;
        memset(&e, 0, sizeof(e));
        if (len < 0 || mio_data_info_entry(info, i, &e) != MIO_OK) {
            mio_data_info_free(info);
            mio_o_fail("data_info entry");
        }
        char *nm = (char *)mxCalloc((size_t)len + 1, 1);
        od_info_ctx ctx = {info, i};
        if (len > 0) od_info_name_getter(&ctx, nm, len + 1);

        int64_t nc = e.num_components > 0 ? e.num_components : 0;
        mxArray *comps = mio_o_matrix(nc, 3);
        double *cd = mxGetPr(comps);
        for (int64_t c = 0; c < nc; ++c) {
            double lo = 0.0, hi = 0.0, mu = 0.0;
            if (mio_data_info_component(info, i, c, &lo, &hi, &mu) != MIO_OK) {
                mio_data_info_free(info);
                mio_o_fail("data_info component");
            }
            cd[c] = lo;
            cd[nc + c] = hi;
            cd[2 * nc + c] = mu;
        }
        mwIndex k = (mwIndex)i;
        mxSetField(out, k, "name", mxCreateString(nm));
        mxFree(nm);
        mxSetField(out, k, "location",
                   mio_o_str((e.location >= 0 && e.location <= 2) ? kLocations[e.location]
                                                                  : "unknown"));
        mxSetField(out, k, "dtype", mio_o_str(mio_o_dtype_name((mio_dtype)e.dtype)));
        mxSetField(out, k, "num_blocks", mio_o_scalar((double)e.num_blocks));
        mxSetField(out, k, "num_entries", mio_o_scalar((double)e.num_entries));
        mxSetField(out, k, "num_components", mio_o_scalar((double)e.num_components));
        mxSetField(out, k, "num_values", mio_o_scalar((double)e.num_values));
        mxSetField(out, k, "min", mio_o_scalar(e.min));
        mxSetField(out, k, "max", mio_o_scalar(e.max));
        mxSetField(out, k, "mean", mio_o_scalar(e.mean));
        mxSetField(out, k, "num_nan", mio_o_scalar((double)e.num_nan));
        mxSetField(out, k, "num_inf", mio_o_scalar((double)e.num_inf));
        mxSetField(out, k, "num_finite", mio_o_scalar((double)e.num_finite));
        mxSetField(out, k, "inconsistent_blocks", mio_o_logical(e.inconsistent_blocks != 0));
        mxSetField(out, k, "components", comps);
    }
    mio_data_info_free(info);
    return out;
}

/* The C protocol for a name, read directly so a failure can free `result`
 * before raising. region < 0 reads the array's name, else that region's. */
static char *od_integrate_name(mio_data_integrate *result, int64_t index, int64_t region) {
    int64_t len = region < 0 ? mio_data_integrate_name(result, index, NULL, 0)
                             : mio_data_integrate_region_name(result, index, region, NULL, 0);
    if (len < 0) {
        mio_data_integrate_free(result);
        mio_o_fail("data_integrate name");
    }
    char *buf = (char *)mxCalloc((size_t)len + 1, 1);
    if (len > 0) {
        int64_t got = region < 0
                          ? mio_data_integrate_name(result, index, buf, len + 1)
                          : mio_data_integrate_region_name(result, index, region, buf, len + 1);
        if (got < 0) {
            mio_data_integrate_free(result);
            mio_o_fail("data_integrate name");
        }
    }
    return buf;
}

/* Build the (num_components x 4) total/mean/domain_measure/num_nan matrix for
 * one array's whole-mesh entry (region < 0) or one of its named Cell
 * regions (region >= 0). Frees `result` and raises on any C API failure --
 * the caller never sees a NULL back. */
static mxArray *od_integrate_components(mio_data_integrate *result, int64_t index,
                                        int64_t region, int64_t num_components) {
    int64_t nc = num_components > 0 ? num_components : 0;
    mxArray *m = mio_o_matrix(nc, 4);
    double *d = mxGetPr(m);
    for (int64_t c = 0; c < nc; ++c) {
        double total = 0.0, mean = 0.0, dm = 0.0;
        int64_t nan = 0;
        mio_status s = (region < 0) ? mio_data_integrate_component(result, index, c, &total,
                                                                    &mean, &dm, &nan)
                                    : mio_data_integrate_region_component(
                                          result, index, region, c, &total, &mean, &dm, &nan);
        if (s != MIO_OK) {
            mio_data_integrate_free(result);
            mio_o_fail("data_integrate component");
        }
        d[c] = total;
        d[nc + c] = mean;
        d[2 * nc + c] = dm;
        d[3 * nc + c] = (double)nan;
    }
    return m;
}

static const char *kIntegrateFields[] = {"name", "num_components", "domain", "regions"};
static const char *kIntegrateRegionFields[] = {"name", "num_cells", "num_skipped",
                                               "components"};

/* Cell-measure-weighted total/mean of one or more cell_data arrays --
 * gradient's integration counterpart. `names` [] means every cell_data
 * array. Reported for the whole mesh (`domain`) and independently for every
 * named Cell region (`regions`) -- a cell in two regions contributes fully
 * to both, one in none contributes to neither. A point_data-only name fails,
 * naming data_point_to_cell as the fix. See doc/field_integration.md. */
static mxArray *O_mio_data_integrate(const mxArray **a) {
    int64_t ncount = 0;
    const char *const *nm = mio_o_names(a[1], &ncount, "names");
    mio_data_integrate *result = (mio_data_integrate *)mio_o_check_ptr(
        mio_data_integrate_create(mio_o_mesh(a[0]), nm, ncount), "data_integrate");

    int64_t n = mio_data_integrate_count(result);
    if (n < 0) {
        mio_data_integrate_free(result);
        mio_o_fail("data_integrate count");
    }

    mxArray *out = mio_o_struct_array(n, 4, kIntegrateFields);
    for (int64_t i = 0; i < n; ++i) {
        char *name = od_integrate_name(result, i, -1);
        mio_field_integral_info e;
        memset(&e, 0, sizeof(e));
        if (mio_data_integrate_entry(result, i, &e) != MIO_OK) {
            mio_data_integrate_free(result);
            mio_o_fail("data_integrate entry");
        }
        const char *domain_names[] = {"num_cells", "num_skipped", "components"};
        mxArray *domain_values[] = {mio_o_scalar((double)e.num_cells),
                                    mio_o_scalar((double)e.num_skipped),
                                    od_integrate_components(result, i, -1, e.num_components)};
        mxArray *domain = mio_o_struct(3, domain_names, domain_values);

        int64_t nregions = mio_data_integrate_region_count(result, i);
        if (nregions < 0) {
            mio_data_integrate_free(result);
            mio_o_fail("data_integrate region count");
        }
        mxArray *regions = mio_o_struct_array(nregions, 4, kIntegrateRegionFields);
        for (int64_t r = 0; r < nregions; ++r) {
            char *rname = od_integrate_name(result, i, r);
            mio_field_integral_info re;
            memset(&re, 0, sizeof(re));
            if (mio_data_integrate_region_entry(result, i, r, &re) != MIO_OK) {
                mio_data_integrate_free(result);
                mio_o_fail("data_integrate region entry");
            }
            mwIndex rk = (mwIndex)r;
            mxSetField(regions, rk, "name", mxCreateString(rname));
            mxFree(rname);
            mxSetField(regions, rk, "num_cells", mio_o_scalar((double)re.num_cells));
            mxSetField(regions, rk, "num_skipped", mio_o_scalar((double)re.num_skipped));
            mxSetField(regions, rk, "components",
                       od_integrate_components(result, i, r, re.num_components));
        }

        mwIndex k = (mwIndex)i;
        mxSetField(out, k, "name", mxCreateString(name));
        mxFree(name);
        mxSetField(out, k, "num_components", mio_o_scalar((double)e.num_components));
        mxSetField(out, k, "domain", domain);
        mxSetField(out, k, "regions", regions);
    }
    mio_data_integrate_free(result);
    return out;
}

/* --- settings pipeline (v9.11.0) ----------------------------------------- */

static mxArray *O_mio_pipeline_run_file(const mxArray **a) {
    mio_o_check(mio_pipeline_run_file(mio_o_string(a[0], "settings_path")),
                "pipeline_run_file");
    return NULL;
}

/* The report's JSON as a char row vector; the report is freed on every path
 * (R hands it to a finalizer; here it is freed before any raise). */
static mxArray *od_pipeline_report(mio_pipeline_report *report) {
    if (report == NULL) mio_o_fail("pipeline_report");
    int64_t n = mio_pipeline_report_json(report, NULL, 0);
    if (n < 0) od_fail_report(report, "pipeline_report");
    char *buf = (char *)mxCalloc((size_t)n + 1, 1);
    if (mio_pipeline_report_json(report, buf, n + 1) < 0) od_fail_report(report, "pipeline_report");
    buf[n] = '\0';
    mio_pipeline_report_free(report);
    mxArray *out = mxCreateString(buf);
    mxFree(buf);
    return out;
}

static mxArray *O_mio_pipeline_run_file_report(const mxArray **a) {
    return od_pipeline_report(mio_pipeline_run_file_report(mio_o_string(a[0], "path")));
}
static mxArray *O_mio_pipeline_run_json_report(const mxArray **a) {
    return od_pipeline_report(mio_pipeline_run_json_report(mio_o_string(a[0], "text")));
}
static mxArray *O_mio_sequence_pipeline_run_file_report(const mxArray **a) {
    return od_pipeline_report(mio_sequence_pipeline_run_file_report(mio_o_string(a[0], "path")));
}
static mxArray *O_mio_sequence_pipeline_run_json_report(const mxArray **a) {
    return od_pipeline_report(mio_sequence_pipeline_run_json_report(mio_o_string(a[0], "text")));
}

static mxArray *O_mio_pipeline_run_json(const mxArray **a) {
    mio_o_check(mio_pipeline_run_json(mio_o_string(a[0], "json_text")), "pipeline_run_json");
    return NULL;
}

static mxArray *O_mio_pipeline_has_json(const mxArray **a) {
    (void)a;
    return mio_o_logical(mio_pipeline_has_json() != 0);
}

/* --- regular grids and signed distance ----------------------------------- */

/* MATLAB has int64, but counts come back as doubles (exact well past any
 * plausible cell count) -- the convention every other counter here follows.
 * Length-3 results (dims, origin, spacing) are 3 x 1 columns. */

static void od_fill_sdf_opts(mio_sdf_opts *o, const mxArray *sign, int location,
                             const mxArray *band, int record_inside,
                             const mxArray *watertight_check) {
    mio_sdf_opts_init(o);
    o->sign = (int32_t)mio_o_int(sign, "sign");
    o->location = (int32_t)location;
    o->band = mio_o_double(band, "band");
    o->record_inside = record_inside;
    o->watertight_check = (int32_t)mio_o_int(watertight_check, "watertight_check");
}

static mxArray *od_quality_struct(const mio_surface_quality *q) {
    const char *names[] = {"boundary_edges", "non_manifold_edges", "inconsistent_pairs",
                           "degenerate_triangles", "watertight"};
    mxArray *values[] = {mio_o_scalar((double)q->boundary_edges),
                         mio_o_scalar((double)q->non_manifold_edges),
                         mio_o_scalar((double)q->inconsistent_pairs),
                         mio_o_scalar((double)q->degenerate_triangles),
                         mio_o_logical(q->watertight != 0)};
    return mio_o_struct(5, names, values);
}

static mxArray *od_column3(const double *v) {
    mxArray *out = mio_o_column(3);
    memcpy(mxGetPr(out), v, 3 * sizeof(double));
    return out;
}

static mxArray *od_column3_i64(const int64_t *v) {
    mxArray *out = mio_o_column(3);
    for (int i = 0; i < 3; ++i) mxGetPr(out)[i] = (double)v[i];
    return out;
}

/* grid(dims, origin, spacing, max_cells) */
static mxArray *O_mio_grid(const mxArray **a) {
    double d[3];
    mio_o_vecn(a[0], d, 3, "dims");
    int64_t cdims[3];
    for (int i = 0; i < 3; ++i) cdims[i] = (int64_t)d[i];
    double o[3], s[3];
    mio_o_vec3(a[1], o, "origin");
    mio_o_vec3(a[2], s, "spacing");
    mio_mesh *out = mio_grid(cdims, o, s, mio_o_int64(a[3], "max_cells"));
    if (out == NULL) mio_o_fail("grid");
    return mio_o_wrap_mesh(out);
}

/* The parametric generators share one options struct; `shape` picks the entry
 * point. cell_type and split arrive as names and are mapped here so the
 * Octave layer and the other language bindings refuse the same spellings. */
static int32_t od_gen_cell(const char *name) {
    if (name[0] == '\0') return MIO_GEN_DEFAULT;
    if (strcmp(name, "quad") == 0) return MIO_GEN_QUAD;
    if (strcmp(name, "triangle") == 0) return MIO_GEN_TRIANGLE;
    if (strcmp(name, "hexahedron") == 0) return MIO_GEN_HEXAHEDRON;
    if (strcmp(name, "tetra") == 0) return MIO_GEN_TETRA;
    if (strcmp(name, "wedge") == 0) return MIO_GEN_WEDGE;
    mio_o_error("generate", "unknown cell_type '%s'", name);
    return 0;
}

static int32_t od_gen_split(const char *name) {
    if (strcmp(name, "zigzag") == 0) return MIO_GEN_SPLIT_ZIGZAG;
    if (strcmp(name, "up") == 0) return MIO_GEN_SPLIT_UP;
    if (strcmp(name, "down") == 0) return MIO_GEN_SPLIT_DOWN;
    if (strcmp(name, "center") == 0) return MIO_GEN_SPLIT_CENTER;
    mio_o_error("generate", "unknown split '%s'", name);
    return 0;
}

typedef mio_mesh *(*od_gen_fn)(const mio_gen_opts *);

static const struct {
    const char *name;
    od_gen_fn fn;
} kGenerators[] = {
    {"rectangle", mio_rectangle},
    {"box", mio_box},
    {"disk", mio_disk},
    {"uv_sphere", mio_uv_sphere},
    {"icosa_sphere", mio_icosa_sphere},
    {"cylinder", mio_cylinder},
    {"reference_triangle", mio_reference_triangle},
    {"ngon", mio_ngon},
    {"annulus", mio_annulus},
    {"lshape", mio_lshape},
    {"plate_with_hole", mio_plate_with_hole},
    {"honeycomb", mio_honeycomb},
    {"line", mio_line},
    {"circle", mio_circle},
    {"helix", mio_helix},
    {"sphere_points", mio_sphere_points},
    {"box_points", mio_box_points},
    {"tetra_sphere", mio_tetra_sphere},
    {"octa_sphere", mio_octa_sphere},
    {"cube_sphere", mio_cube_sphere},
    {"geo_sphere", mio_geo_sphere},
    {"torus", mio_torus},
    {"moebius", mio_moebius},
    {"tube", mio_tube},
    {"cone", mio_cone},
    {"capsule", mio_capsule},
    {"box_surface", mio_box_surface},
    {"platonic", mio_platonic},
    {"ball", mio_ball},
    {"spherical_shell", mio_spherical_shell},
    {"pipe", mio_pipe},
    {"fichera", mio_fichera},
    {"solid_torus", mio_solid_torus},
};

/* generate(shape, counts, size, origin, radius, height, cell_type, split,
 *          order, max_cells, radius2, angle, num_twists, complete, ratio, seed) */
static mxArray *O_mio_generate(const mxArray **a) {
    mio_gen_opts opts;
    mio_gen_opts_init(&opts);
    size_t nc = 0;
    double *c = mio_o_doubles(a[1], &nc, "counts");
    if (nc < 1 || nc > 3) mio_o_error("generate", "`counts` must have 1 to 3 elements");
    for (size_t i = 0; i < nc; ++i) opts.counts[i] = (int64_t)c[i];
    for (size_t i = nc; i < 3; ++i) opts.counts[i] = 1;
    mxFree(c);
    mio_o_vec3(a[2], opts.size, "size");
    mio_o_vec3(a[3], opts.origin, "origin");
    opts.radius = mio_o_double(a[4], "radius");
    opts.height = mio_o_double(a[5], "height");
    opts.radius2 = mio_o_double(a[10], "radius2");
    opts.angle = mio_o_double(a[11], "angle");
    opts.cell_type = od_gen_cell(mio_o_string(a[6], "cell_type"));
    opts.split = od_gen_split(mio_o_string(a[7], "split"));
    opts.order = (int32_t)mio_o_int(a[8], "order");
    opts.max_cells = mio_o_int64(a[9], "max_cells");
    opts.num_twists = (int32_t)mio_o_int(a[12], "num_twists");
    opts.complete = mio_o_bool(a[13], "complete") ? 1 : 0;
    opts.ratio = mio_o_double(a[14], "ratio");
    opts.seed = (uint64_t)mio_o_double(a[15], "seed");

    const char *name = mio_o_string(a[0], "shape");
    od_gen_fn fn = NULL;
    for (size_t i = 0; i < sizeof(kGenerators) / sizeof(kGenerators[0]); ++i) {
        if (strcmp(name, kGenerators[i].name) == 0) fn = kGenerators[i].fn;
    }
    if (fn == NULL) mio_o_error("generate", "unknown shape '%s'", name);
    mio_mesh *out = fn(&opts);
    if (out == NULL) mio_o_fail(name);
    return mio_o_wrap_mesh(out);
}

static mxArray *O_mio_reference_cells(const mxArray **a) {
    mio_mesh *out = mio_reference_cells(mio_o_string(a[0], "types"), mio_o_bool(a[1], "pairs"));
    if (out == NULL) mio_o_fail("reference_cells");
    return mio_o_wrap_mesh(out);
}

/* voxelize(mesh, resolution, cell_size, bounds, padding, padding_relative,
 *          fill, sign, attach_occupancy, max_cells, watertight_check, plane) */
static mxArray *O_mio_voxelize(const mxArray **a) {
    mio_voxel_opts opts;
    mio_voxel_opts_init(&opts);

    /* The buffers the option pointers reference must outlive the call. */
    int64_t res[3];
    double bnd[6], r[3];
    if (od_opt_vecn(a[1], r, 3, "resolution")) {
        for (int i = 0; i < 3; ++i) res[i] = (int64_t)r[i];
        opts.resolution = res;
    }
    if (od_opt_vecn(a[3], bnd, 6, "bounds")) opts.bounds = bnd;
    opts.cell_size = mio_o_double(a[2], "cell_size");
    opts.padding = mio_o_double(a[4], "padding");
    opts.padding_relative = mio_o_double(a[5], "padding_relative");
    opts.fill = (int32_t)mio_o_int(a[6], "fill");
    opts.sign = (int32_t)mio_o_int(a[7], "sign");
    opts.attach_occupancy = mio_o_bool(a[8], "attach_occupancy");
    opts.max_cells = mio_o_int64(a[9], "max_cells");
    opts.watertight_check = (int32_t)mio_o_int(a[10], "watertight_check");
    opts.plane = (int32_t)mio_o_int(a[11], "plane");

    int64_t cdims[3] = {0, 0, 0}, occupied = 0;
    double origin[3] = {0, 0, 0}, spacing[3] = {0, 0, 0};
    mio_mesh *out = mio_voxelize(mio_o_mesh(a[0]), &opts, cdims, origin, spacing, &occupied);
    if (out == NULL) mio_o_fail("voxelize");

    const char *names[] = {"mesh", "dims", "origin", "spacing", "num_occupied"};
    mxArray *values[] = {mio_o_wrap_mesh(out), od_column3_i64(cdims), od_column3(origin),
                         od_column3(spacing), mio_o_scalar((double)occupied)};
    return mio_o_struct(5, names, values);
}

/* octree(mesh, root_cells, min_level, max_level, bounds, padding,
 *        padding_relative, dimension, sample_distance, balance, fill, sign,
 *        attach_distance, even_odd, max_cells, watertight_check, transition,
 *        snap, surface_level, curvature_threshold, feature_angle, size_array,
 *        plane) */
static mxArray *O_mio_octree(const mxArray **a) {
    mio_octree_opts opts;
    mio_octree_opts_init(&opts);

    /* The buffers the option pointers reference must outlive the call. */
    int64_t roots[3];
    double bnd[6], r[3];
    if (od_opt_vecn(a[1], r, 3, "root_cells")) {
        for (int i = 0; i < 3; ++i) roots[i] = (int64_t)r[i];
        opts.root_cells = roots;
    }
    if (od_opt_vecn(a[4], bnd, 6, "bounds")) opts.bounds = bnd;
    opts.min_level = mio_o_int64(a[2], "min_level");
    opts.max_level = mio_o_int64(a[3], "max_level");
    opts.padding = mio_o_double(a[5], "padding");
    opts.padding_relative = mio_o_double(a[6], "padding_relative");
    opts.dimension = mio_o_int64(a[7], "dimension");
    opts.sample_distance = mio_o_double(a[8], "sample_distance");
    opts.balance = (int32_t)mio_o_int(a[9], "balance");
    opts.fill = (int32_t)mio_o_int(a[10], "fill");
    opts.sign = (int32_t)mio_o_int(a[11], "sign");
    opts.attach_distance = mio_o_bool(a[12], "attach_distance");
    opts.even_odd = mio_o_bool(a[13], "even_odd");
    opts.max_cells = mio_o_int64(a[14], "max_cells");
    opts.watertight_check = (int32_t)mio_o_int(a[15], "watertight_check");
    opts.transition = (int32_t)mio_o_int(a[16], "transition");
    opts.snap = mio_o_bool(a[17], "snap");
    opts.surface_level = mio_o_int64(a[18], "surface_level");
    opts.curvature_threshold = mio_o_double(a[19], "curvature_threshold");
    opts.feature_angle = mio_o_double(a[20], "feature_angle");
    opts.plane = (int32_t)mio_o_int(a[22], "plane");
    const char *size_array = mio_o_opt_string(a[21]);
    if (size_array != NULL && size_array[0] != '\0') opts.size_array = size_array;

    mio_octree_report report;
    mio_mesh *out = mio_octree(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("octree");

    const char *names[] = {"mesh",         "dimension", "num_leaves",  "num_hanging",
                           "finest_level", "origin",    "root_size",   "num_snapped",
                           "num_refused"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)report.dimension),
                         mio_o_scalar((double)report.num_leaves),
                         mio_o_scalar((double)report.num_hanging),
                         mio_o_scalar((double)report.finest_level),
                         od_column3(report.origin),
                         mio_o_scalar(report.root_size),
                         mio_o_scalar((double)report.num_snapped),
                         mio_o_scalar((double)report.num_refused)};
    return mio_o_struct(9, names, values);
}

/* mio_o_getstring adapter: the warnings of this thread's last octree_mesh call. */
static int64_t od_octree_mesh_warnings(void *ctx, char *buf, int64_t buflen) {
    (void)ctx;
    return mio_octree_mesh_warnings(buf, buflen);
}

/* Every option is [] for "keep mio_octree_mesh_opts_init's default", so the
 * defaults come from the C side; the .m wrapper maps the enum strings. The
 * argument order is R's .Call order:
 *   mesh, output, conformity, skin_grading, fill, bounds, size_array,
 *   surface_region, min_size, max_size, padding, padding_relative,
 *   curvature_angle, feature_angle, corner_angle, feature_size, gap_cells,
 *   surface_size, winding_threshold, winding_beta, alpha_long, alpha_short,
 *   max_level, grading_layers, fit_iterations, untangle_iterations,
 *   max_cells, max_tets, gap_outside, preserve_features,
 *   dihedral_floor_guard, keep_outside, attach_distance, attach_quality,
 *   pillow, hex_grading */
static mxArray *O_mio_octree_mesh(const mxArray **a) {
    mio_octree_mesh_opts opts;
    mio_octree_mesh_opts_init(&opts);

#define OD_ENUM(i, f) \
    if (!mio_o_is_null(a[i])) opts.f = (int32_t)mio_o_int(a[i], #f)
#define OD_DOUBLE(i, f) \
    if (!mio_o_is_null(a[i])) opts.f = mio_o_double(a[i], #f)
#define OD_INT64(i, f) \
    if (!mio_o_is_null(a[i])) opts.f = mio_o_int64(a[i], #f)
#define OD_BOOL(i, f) \
    if (!mio_o_is_null(a[i])) opts.f = (int32_t)mio_o_bool(a[i], #f)
    OD_ENUM(1, output);
    OD_ENUM(2, conformity);
    OD_ENUM(3, skin_grading);
    OD_ENUM(4, fill);
    OD_DOUBLE(8, min_size);
    OD_DOUBLE(9, max_size);
    OD_DOUBLE(10, padding);
    OD_DOUBLE(11, padding_relative);
    OD_DOUBLE(12, curvature_angle);
    OD_DOUBLE(13, feature_angle);
    OD_DOUBLE(14, corner_angle);
    OD_DOUBLE(15, feature_size);
    OD_DOUBLE(16, gap_cells);
    OD_DOUBLE(17, surface_size);
    OD_DOUBLE(18, winding_threshold);
    OD_DOUBLE(19, winding_beta);
    OD_DOUBLE(20, alpha_long);
    OD_DOUBLE(21, alpha_short);
    OD_INT64(22, max_level);
    OD_INT64(23, grading_layers);
    OD_INT64(24, fit_iterations);
    OD_INT64(25, untangle_iterations);
    OD_INT64(26, max_cells);
    OD_INT64(27, max_tets);
    OD_BOOL(28, gap_outside);
    OD_BOOL(29, preserve_features);
    OD_BOOL(30, dihedral_floor_guard);
    OD_BOOL(31, keep_outside);
    OD_BOOL(32, attach_distance);
    OD_BOOL(33, attach_quality);
    OD_BOOL(34, pillow);
    OD_ENUM(35, hex_grading);
#undef OD_ENUM
#undef OD_DOUBLE
#undef OD_INT64
#undef OD_BOOL

    /* The buffers the option pointers reference must outlive the call; the
     * names are mxMalloc'd copies, alive until the MEX call returns. */
    double bnd[6];
    if (od_opt_vecn(a[5], bnd, 6, "bounds")) opts.bounds = bnd;
    const char *size_name = mio_o_opt_string(a[6]);
    if (size_name != NULL && size_name[0] != '\0') opts.size_array = size_name;
    const char *region_name = mio_o_opt_string(a[7]);
    if (region_name != NULL && region_name[0] != '\0') opts.surface_region = region_name;

    mio_octree_mesh_report report;
    mio_mesh *out = mio_octree_mesh(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("octree_mesh");

    /* Wrap the mesh first: once it is a handle object nothing below can leak
     * it. */
    mxArray *mo = mio_o_wrap_mesh(out);
    mxArray *warnings = mio_o_getstring(od_octree_mesh_warnings, NULL, "octree_mesh");
    static const char *outputs[] = {"castellated", "tet", "hex"};
    const char *output_name =
        (report.output >= 0 && report.output <= 2) ? outputs[report.output] : "unknown";

    enum { OD_OCTREE_MESH_N = 41 };
    const char *names[OD_OCTREE_MESH_N];
    mxArray *values[OD_OCTREE_MESH_N];
    int k = 0;
#define OD_SET(name, value)  \
    do {                     \
        names[k] = (name);   \
        values[k] = (value); \
        ++k;                 \
    } while (0)
#define OD_NUM(f) OD_SET(#f, mio_o_scalar((double)report.f))
    OD_SET("mesh", mo);
    OD_SET("output", mio_o_str(output_name));
    OD_SET("guarantee_applies", mio_o_logical(report.guarantee_applies != 0));
    OD_SET("origin", od_column3(report.origin));
    OD_NUM(root_size);
    OD_SET("root_cells", od_column3_i64(report.root_cells));
    OD_NUM(max_level);
    OD_NUM(finest_level);
    OD_NUM(surface_level);
    OD_NUM(min_size);
    OD_NUM(max_size);
    OD_NUM(num_leaves);
    OD_NUM(num_points);
    OD_NUM(num_cells);
    OD_NUM(num_hexahedra);
    OD_NUM(num_tetrahedra);
    OD_NUM(num_inside);
    OD_NUM(num_boundary_cells);
    OD_NUM(num_surrogate_faces);
    OD_NUM(num_hanging);
    OD_NUM(num_cut_points);
    OD_NUM(num_warped);
    OD_NUM(num_warp_refused);
    OD_NUM(min_dihedral);
    OD_NUM(max_dihedral);
    OD_NUM(min_scaled_jacobian);
    OD_NUM(hausdorff_to_skin);
    OD_NUM(num_inverted);
    OD_NUM(num_boundary_edges_open);
    OD_NUM(num_non_manifold_boundary_edges);
    OD_NUM(num_corners);
    OD_NUM(num_corners_snapped);
    OD_NUM(num_feature_chains);
    OD_NUM(num_feature_chains_missed);
    OD_NUM(num_feature_moves);
    OD_NUM(num_optimized_moves);
    OD_NUM(num_flips);
    OD_NUM(num_winding_rechecks);
    OD_NUM(num_winding_skipped);
    OD_NUM(num_warnings);
    /* One string, newline-separated; the .m wrapper splits it. */
    OD_SET("warnings", warnings);
#undef OD_NUM
#undef OD_SET
    if (k != OD_OCTREE_MESH_N) {
        mio_o_error("octree_mesh", "octree_mesh: report has %d fields, expected %d", k,
                    (int)OD_OCTREE_MESH_N);
    }
    return mio_o_struct(OD_OCTREE_MESH_N, names, values);
}

/* compute_sdf(mesh, structure, resolution, cell_size, bounds, padding,
 *             padding_relative, root_resolution, max_depth, band_cells,
 *             record_levels, max_cells, sign, location, band, watertight_check) */
static mxArray *O_mio_compute_sdf(const mxArray **a) {
    mio_compute_sdf_opts opts;
    mio_compute_sdf_opts_init(&opts);

    /* The buffers the option pointers reference must outlive the call. */
    int64_t res[3];
    double bnd[6], r[3];
    if (od_opt_vecn(a[2], r, 3, "resolution")) {
        for (int i = 0; i < 3; ++i) res[i] = (int64_t)r[i];
        opts.resolution = res;
    }
    if (od_opt_vecn(a[4], bnd, 6, "bounds")) opts.bounds = bnd;
    opts.structure = (int32_t)mio_o_int(a[1], "structure");
    opts.cell_size = mio_o_double(a[3], "cell_size");
    opts.padding = mio_o_double(a[5], "padding");
    opts.padding_relative = mio_o_double(a[6], "padding_relative");
    opts.root_resolution = mio_o_int64(a[7], "root_resolution");
    opts.max_depth = mio_o_int64(a[8], "max_depth");
    opts.band_cells = mio_o_double(a[9], "band_cells");
    opts.record_levels = mio_o_bool(a[10], "record_levels");
    opts.max_cells = mio_o_int64(a[11], "max_cells");
    od_fill_sdf_opts(&opts.distance, a[12], mio_o_int(a[13], "location"), a[14], 0, a[15]);

    int64_t cdims[3] = {0, 0, 0}, depth = 0, banded = 0;
    double origin[3] = {0, 0, 0}, spacing[3] = {0, 0, 0};
    mio_surface_quality q;
    memset(&q, 0, sizeof(q));
    mio_mesh *out = mio_compute_sdf(mio_o_mesh(a[0]), &opts, cdims, origin, spacing, &depth,
                                    &banded, &q);
    if (out == NULL) mio_o_fail("compute_sdf");

    const char *names[] = {"mesh",      "dims",       "origin", "spacing",
                           "max_depth", "num_banded", "quality"};
    mxArray *values[] = {mio_o_wrap_mesh(out),         od_column3_i64(cdims),
                         od_column3(origin),           od_column3(spacing),
                         mio_o_scalar((double)depth),  mio_o_scalar((double)banded),
                         od_quality_struct(&q)};
    return mio_o_struct(7, names, values);
}

/* Retetrahedralize a volume mesh (or closed surface) by isosurface stuffing
 * -- remesh's volumetric sibling. `distance` is set directly rather than
 * through od_fill_sdf_opts, whose location/band/record_inside fields
 * remesh_volume never reads.
 * remesh_volume(mesh, resolution, cell_size, bounds, padding,
 *               padding_relative, max_cells, max_tets, warp_fraction, sign,
 *               watertight_check, surface_region, grid_cell_size,
 *               max_winding_work) */
static mxArray *O_mio_remesh_volume(const mxArray **a) {
    mio_remesh_volume_opts opts;
    mio_remesh_volume_report report;
    mio_remesh_volume_opts_init(&opts);

    /* The buffers the option pointers reference must outlive the call, the
     * compute_sdf idiom above. */
    int64_t res[3];
    double bnd[6], r[3];
    if (od_opt_vecn(a[1], r, 3, "resolution")) {
        for (int i = 0; i < 3; ++i) res[i] = (int64_t)r[i];
        opts.resolution = res;
    }
    if (od_opt_vecn(a[3], bnd, 6, "bounds")) opts.bounds = bnd;
    opts.cell_size = mio_o_double(a[2], "cell_size");
    opts.padding = mio_o_double(a[4], "padding");
    opts.padding_relative = mio_o_double(a[5], "padding_relative");
    opts.max_cells = mio_o_int64(a[6], "max_cells");
    opts.max_tets = mio_o_int64(a[7], "max_tets");
    opts.warp_fraction = mio_o_double(a[8], "warp_fraction");
    opts.distance.sign = (int32_t)mio_o_int(a[9], "sign");
    opts.distance.watertight_check = (int32_t)mio_o_int(a[10], "watertight_check");
    opts.distance.surface_region = mio_o_opt_string(a[11]);
    opts.distance.grid_cell_size = mio_o_double(a[12], "grid_cell_size");
    opts.distance.max_winding_work = mio_o_double(a[13], "max_winding_work");

    mio_mesh *out = mio_remesh_volume_ex(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("remesh_volume");
    const char *names[] = {"mesh", "num_tets", "num_vertices_warped", "num_tets_rejected",
                           "num_non_manifold_edges", "input_quality"};
    mxArray *values[] = {mio_o_wrap_mesh(out),
                         mio_o_scalar((double)report.num_tets),
                         mio_o_scalar((double)report.num_vertices_warped),
                         mio_o_scalar((double)report.num_tets_rejected),
                         mio_o_scalar((double)report.num_non_manifold_edges),
                         od_quality_struct(&report.input_quality)};
    return mio_o_struct(6, names, values);
}

static mxArray *O_mio_surface_watertight_check(const mxArray **a) {
    mio_surface_quality q;
    memset(&q, 0, sizeof(q));
    mio_o_check(mio_surface_watertight_check(mio_o_mesh(a[0]), &q), "surface_watertight_check");
    return od_quality_struct(&q);
}

/* sample_distance(mesh, points, sign, band, watertight_check) -- `points` is
 * (n x 3), one query point per row; the result an n x 1 column. */
static mxArray *O_mio_sample_distance(const mxArray **a) {
    const mxArray *pts = a[1];
    if (pts == NULL || mxGetNumberOfDimensions(pts) != 2 || mxGetN(pts) != 3) {
        mio_o_error("sample_distance", "`points` must be an (n x 3) matrix");
    }
    int64_t n = 0, k = 0;
    /* rows_to_c transposes to the row-major (n, 3) buffer the C API reads. */
    double *p = mio_o_rows_to_c(pts, &n, &k, 0.0, "points");

    mio_sdf_opts opts;
    od_fill_sdf_opts(&opts, a[2], 0, a[3], 0, a[4]);

    mxArray *out = mio_o_column(n);
    mio_o_check(mio_sample_distance(mio_o_mesh(a[0]), p, n, &opts, mxGetPr(out)),
                "sample_distance");
    mxFree(p);
    return out;
}

/* distance_to_surface(mesh, surface, sign, location, band, record_inside,
 *                     watertight_check) */
static mxArray *O_mio_distance_to_surface(const mxArray **a) {
    mio_sdf_opts opts;
    od_fill_sdf_opts(&opts, a[2], mio_o_int(a[3], "location"), a[4],
                     mio_o_bool(a[5], "record_inside"), a[6]);

    int64_t banded = 0;
    mio_surface_quality q;
    memset(&q, 0, sizeof(q));
    mio_mesh *out =
        mio_distance_to_surface(mio_o_mesh(a[0]), mio_o_mesh(a[1]), &opts, &banded, &q);
    if (out == NULL) mio_o_fail("distance_to_surface");

    const char *names[] = {"mesh", "num_banded", "quality"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)banded),
                         od_quality_struct(&q)};
    return mio_o_struct(3, names, values);
}

/* --- Provenance (see doc/provenance.md) ----------------------------------
 *
 * Rides the flat C ABI. The begin/end pair is exposed as two entry points and
 * the .m-side `mio.with_provenance()` wrapper pairs them with onCleanup -- so
 * the scope closes even when the wrapped function handle throws. */

static mxArray *O_mio_provenance_begin(const mxArray **a) {
    if (mio_provenance_scope_begin(mio_o_int(a[0], "mode")) != MIO_OK)
        mio_o_fail("provenance_scope_begin");
    return NULL;
}

static mxArray *O_mio_provenance_end(const mxArray **a) {
    (void)a;
    mio_provenance_scope_end();
    return NULL;
}

static mxArray *O_mio_provenance_note(const mxArray **a) {
    mio_provenance_note(mio_o_string(a[0], "category"), mio_o_string(a[1], "detail"));
    return NULL;
}

static mxArray *O_mio_provenance_set_source(const mxArray **a) {
    mio_provenance_set_source(mio_o_string(a[0], "path"), mio_o_string(a[1], "format"));
    return NULL;
}

static mxArray *O_mio_provenance_set_target(const mxArray **a) {
    mio_provenance_set_target(mio_o_string(a[0], "format"), mio_o_string(a[1], "encoding"),
                              mio_o_string(a[2], "codec"),
                              mio_o_string(a[3], "float_format"));
    return NULL;
}

/* --- select -------------------------------------------------------------- */

/* The comparison names refine and select share (the R binding's
 * refine_compare_code; "" is "<"). */
static int od_compare_code(const char *name) {
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

/* A numeric vector of exactly `n` values into `dst`; returns 0 when `s` is
 * [] or empty (the criterion is off). */
static int od_select_vec(const mxArray *s, double *dst, size_t n, const char *what) {
    if (s == NULL || mxGetNumberOfElements(s) == 0) return 0;
    if (!(mxIsNumeric(s) || mxIsLogical(s)) || mxGetNumberOfElements(s) != n)
        mio_o_error("select", "meshio++: select: %s must be %d numbers", what, (int)n);
    mio_o_vecn(s, dst, n, what);
    return 1;
}

/* select(mesh, name, kind, tag, mode, boundary_only, box, plane, sphere,
 *        where_array, where_op, where_value, quality_metric, quality_op,
 *        quality_value, component, seed, feature_angle, normal_tolerance) */
static mxArray *O_mio_select(const mxArray **a) {
    mio_select_opts opts;
    mio_select_report report;
    double v[6];

    mio_select_opts_init(&opts);
    opts.name = mio_o_string(a[1], "name");
    const char *k = mio_o_string(a[2], "kind");
    if (strcmp(k, "point") == 0)
        opts.kind = MIO_REGION_POINT;
    else if (strcmp(k, "cell") == 0)
        opts.kind = MIO_REGION_CELL;
    else if (strcmp(k, "side") == 0)
        opts.kind = MIO_REGION_SIDE;
    else
        mio_o_error("select", "meshio++: select: kind must be 'point', 'cell' or 'side'");
    opts.tag = (int64_t)mio_o_double(a[3], "tag");
    const char *m = mio_o_string(a[4], "mode");
    if (strcmp(m, "all") == 0)
        opts.mode = MIO_SELECT_ALL;
    else if (strcmp(m, "any") == 0)
        opts.mode = MIO_SELECT_ANY;
    else if (strcmp(m, "centroid") == 0)
        opts.mode = MIO_SELECT_CENTROID;
    else
        mio_o_error("select", "meshio++: select: mode must be 'all', 'any' or 'centroid'");
    opts.boundary_only = mio_o_bool(a[5], "boundary_only") ? 1 : 0;
    if (od_select_vec(a[6], v, 6, "box")) {
        opts.use_box = 1;
        memcpy(opts.box_lo, v, 3 * sizeof(double));
        memcpy(opts.box_hi, v + 3, 3 * sizeof(double));
    }
    if (od_select_vec(a[7], v, 6, "plane")) {
        opts.use_plane = 1;
        memcpy(opts.plane_point, v, 3 * sizeof(double));
        memcpy(opts.plane_normal, v + 3, 3 * sizeof(double));
    }
    if (od_select_vec(a[8], v, 4, "sphere")) {
        opts.use_sphere = 1;
        memcpy(opts.sphere_center, v, 3 * sizeof(double));
        opts.sphere_radius = v[3];
    }
    opts.where_array = mio_o_opt_string(a[9]);
    opts.where_compare = od_compare_code(mio_o_opt_string(a[10]));
    if (opts.where_compare < 0) mio_o_error("select", "meshio++: select: unknown comparison in where_op");
    opts.where_value = mio_o_double(a[11], "where_value");
    opts.quality_metric = mio_o_opt_string(a[12]);
    opts.quality_compare = od_compare_code(mio_o_opt_string(a[13]));
    if (opts.quality_compare < 0)
        mio_o_error("select", "meshio++: select: unknown comparison in quality_op");
    opts.quality_value = mio_o_double(a[14], "quality_value");
    /* 1-based here, like every other index. */
    if (!mio_o_is_null(a[15]) && mxGetNumberOfElements(a[15]) > 0) {
        opts.component = (int64_t)mio_o_double(a[15], "component") - 1;
        if (opts.component < 0) mio_o_error("select", "meshio++: select: component must be >= 1");
    }
    if (od_select_vec(a[16], v, 3, "seed")) {
        opts.use_seed = 1;
        memcpy(opts.seed, v, 3 * sizeof(double));
    }
    opts.feature_angle = mio_o_double(a[17], "feature_angle");
    opts.normal_tolerance = mio_o_double(a[18], "normal_tolerance");

    mio_mesh *out = mio_select(mio_o_mesh(a[0]), &opts, &report);
    if (out == NULL) mio_o_fail("select");
    const char *names[] = {"mesh", "num_candidates", "num_selected", "num_grown",
                           "num_components"};
    mxArray *values[] = {mio_o_wrap_mesh(out), mio_o_scalar((double)report.num_candidates),
                         mio_o_scalar((double)report.num_selected),
                         mio_o_scalar((double)report.num_grown),
                         mio_o_scalar((double)report.num_components)};
    return mio_o_struct(5, names, values);
}

/* --- table --------------------------------------------------------------- */

static const mio_o_entry kEntries[] = {
    MIO_O_ENTRY(mio_sets_to_data, 5),
    MIO_O_ENTRY(mio_data_to_sets, 3),
    MIO_O_ENTRY(mio_data_drop, 4),
    MIO_O_ENTRY(mio_data_keep, 4),
    MIO_O_ENTRY(mio_data_rename, 4),
    MIO_O_ENTRY(mio_data_point_to_cell, 3),
    MIO_O_ENTRY(mio_data_cell_to_point, 4),
    MIO_O_ENTRY(mio_data_calc, 5),
    MIO_O_ENTRY(mio_data_condition, 10),
    MIO_O_ENTRY(mio_tensor_invariants, 7),
    MIO_O_ENTRY(mio_data_info, 1),
    MIO_O_ENTRY(mio_data_integrate, 2),
    MIO_O_ENTRY(mio_pipeline_run_file, 1),
    MIO_O_ENTRY(mio_pipeline_run_file_report, 1),
    MIO_O_ENTRY(mio_pipeline_run_json_report, 1),
    MIO_O_ENTRY(mio_sequence_pipeline_run_file_report, 1),
    MIO_O_ENTRY(mio_sequence_pipeline_run_json_report, 1),
    MIO_O_ENTRY(mio_pipeline_run_json, 1),
    MIO_O_ENTRY(mio_pipeline_has_json, 0),
    MIO_O_ENTRY(mio_grid, 4),
    MIO_O_ENTRY(mio_generate, 16),
    MIO_O_ENTRY(mio_reference_cells, 2),
    MIO_O_ENTRY(mio_voxelize, 12),
    MIO_O_ENTRY(mio_octree, 23),
    MIO_O_ENTRY(mio_octree_mesh, 36),
    MIO_O_ENTRY(mio_compute_sdf, 16),
    MIO_O_ENTRY(mio_remesh_volume, 14),
    MIO_O_ENTRY(mio_surface_watertight_check, 1),
    MIO_O_ENTRY(mio_sample_distance, 5),
    MIO_O_ENTRY(mio_distance_to_surface, 7),
    MIO_O_ENTRY(mio_provenance_begin, 1),
    MIO_O_ENTRY(mio_provenance_end, 0),
    MIO_O_ENTRY(mio_provenance_note, 2),
    MIO_O_ENTRY(mio_provenance_set_source, 2),
    MIO_O_ENTRY(mio_provenance_set_target, 4),
    MIO_O_ENTRY(mio_select, 19),
};

const mio_o_entry *mio_o_table_ops_data(size_t *n) {
    *n = sizeof(kEntries) / sizeof(kEntries[0]);
    return kEntries;
}
