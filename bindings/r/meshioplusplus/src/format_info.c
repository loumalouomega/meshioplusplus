/* Format side channels: what a format carries that a mesh cannot hold (MDPA's
 * tables, geometries, Mesh blocks and constraints), kept by
 * `mio_read_with_info()` in its own external pointer -- own tag, own finalizer,
 * exactly as `mio_mesh` -- for `mio_write_with_info()` to put back.
 *
 * `mio_mdpa_info()` copies the whole MDPA side channel into plain R values in
 * one call, so nothing but the handle ever outlives a .Call. */

#include "mio_r.h"

#include <string.h>

SEXP mio_r_info_tag = NULL;

static void mio_r_info_finalizer(SEXP x) {
    mio_format_info *info = (mio_format_info *)R_ExternalPtrAddr(x);
    if (info != NULL) {
        mio_format_info_free(info);
        R_ClearExternalPtr(x);
    }
}

static SEXP wrap_info(mio_format_info *info) {
    SEXP ptr = PROTECT(R_MakeExternalPtr(info, mio_r_info_tag, R_NilValue));
    R_RegisterCFinalizerEx(ptr, mio_r_info_finalizer, TRUE);
    SEXP cls = PROTECT(Rf_mkString("mio_format_info"));
    Rf_setAttrib(ptr, R_ClassSymbol, cls);
    UNPROTECT(2); /* cls, ptr */
    return ptr;
}

/* NULL for an R NULL (no side channel); an error for anything else foreign. */
static const mio_format_info *info_of(SEXP x, int allow_null) {
    if (allow_null && x == R_NilValue) {
        return NULL;
    }
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != mio_r_info_tag) {
        Rf_error("expected a mio_format_info object");
    }
    const mio_format_info *info = (const mio_format_info *)R_ExternalPtrAddr(x);
    if (info == NULL) {
        Rf_error("this mio_format_info has already been released");
    }
    return info;
}

SEXP R_mio_read_with_info(SEXP path, SEXP format, SEXP lenient) {
    mio_read_opts opts;
    mio_read_opts_init(&opts);
    opts.lenient = mio_r_bool(lenient, "lenient");
    mio_format_info *info = NULL;
    mio_mesh *m = mio_read_with_info(mio_r_string(path, "path"), mio_r_opt_string(format), &opts,
                                     &info);
    if (m == NULL) mio_r_fail("read_with_info");
    SEXP values[2];
    values[0] = PROTECT(mio_r_wrap_mesh(m));
    values[1] = info != NULL ? wrap_info(info) : R_NilValue;
    PROTECT(values[1]);
    const char *names[2] = {"mesh", "info"};
    SEXP out = mio_r_named_list(2, names, values);
    UNPROTECT(2);
    return out;
}

SEXP R_mio_write_with_info(SEXP mesh, SEXP info, SEXP path, SEXP format) {
    mio_r_check(mio_write_with_info(mio_r_string(path, "path"), mio_r_mesh(mesh),
                                    mio_r_opt_string(format), info_of(info, 1)),
                "write_with_info");
    return R_NilValue;
}

SEXP R_mio_format_info_release(SEXP x) {
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != mio_r_info_tag) {
        Rf_error("expected a mio_format_info object");
    }
    mio_format_info *info = (mio_format_info *)R_ExternalPtrAddr(x);
    if (info != NULL) {
        mio_format_info_free(info);
        R_ClearExternalPtr(x);
    }
    return R_NilValue;
}

SEXP R_mio_format_info_is_open(SEXP x) {
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != mio_r_info_tag) {
        return Rf_ScalarLogical(FALSE);
    }
    return Rf_ScalarLogical(R_ExternalPtrAddr(x) != NULL);
}

/* --- string getters ------------------------------------------------------ */

typedef struct {
    const mio_format_info *info;
    int32_t section;
    int64_t index;
    int64_t entry;
    int32_t field;
} mdpa_ctx;

static int64_t get_format(void *ctx, char *buf, int64_t buflen) {
    return mio_format_info_format(((mdpa_ctx *)ctx)->info, buf, buflen);
}

static int64_t get_item_string(void *ctx, char *buf, int64_t buflen) {
    mdpa_ctx *c = (mdpa_ctx *)ctx;
    return mio_mdpa_info_string(c->info, c->section, c->index, c->field, buf, buflen);
}

static int64_t get_data_string(void *ctx, char *buf, int64_t buflen) {
    mdpa_ctx *c = (mdpa_ctx *)ctx;
    return mio_mdpa_info_data_string(c->info, c->section, c->index, c->entry, c->field, buf,
                                     buflen);
}

SEXP R_mio_format_info_format(SEXP x) {
    mdpa_ctx c = {info_of(x, 0), 0, 0, 0, 0};
    return mio_r_getstring(get_format, &c, "format_info_format");
}

static SEXP item_string(const mio_format_info *info, int32_t section, int64_t index,
                        int32_t field) {
    mdpa_ctx c = {info, section, index, 0, field};
    return mio_r_getstring(get_item_string, &c, "mdpa_info");
}

static int64_t item_int(const mio_format_info *info, int32_t section, int64_t index,
                        int32_t field) {
    int64_t v = 0;
    mio_r_check(mio_mdpa_info_int(info, section, index, field, &v), "mdpa_info");
    return v;
}

static int64_t count_of(const mio_format_info *info, int32_t section) {
    int64_t n = mio_mdpa_info_count(info, section);
    if (n < 0) mio_r_fail("mdpa_info");
    return n;
}

/* --- arrays -------------------------------------------------------------- */

/* A borrowed C array as an R value: a vector for 1-D, and for 2-D either the
 * same memory as a (cols, rows) matrix (`transpose` 0 -- the connectivity
 * convention) or a (rows, cols) matrix (`transpose` 1 -- a table). `shift`
 * is added to every value (1 for point references). */
static SEXP array_to_r(const void *data, mio_dtype dtype, int32_t ndim, const int64_t *shape,
                       int transpose, double shift) {
    R_xlen_t n = 1;
    for (int32_t k = 0; k < ndim; ++k) n *= (R_xlen_t)shape[k];
    if (ndim == 0) n = 0;
    SEXP flat = PROTECT(mio_r_copy_as_real(data, dtype, n));
    double *p = REAL(flat);
    for (R_xlen_t i = 0; i < n; ++i) p[i] += shift;
    if (ndim != 2) {
        UNPROTECT(1);
        return flat;
    }
    const int rows = (int)shape[0], cols = (int)shape[1];
    SEXP out;
    if (!transpose) {
        out = PROTECT(Rf_allocMatrix(REALSXP, cols, rows));
        if (n > 0) memcpy(REAL(out), p, (size_t)n * sizeof(double));
    } else {
        out = PROTECT(Rf_allocMatrix(REALSXP, rows, cols));
        double *q = REAL(out);
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c) q[r + (R_xlen_t)c * rows] = p[(R_xlen_t)r * cols + c];
    }
    UNPROTECT(2);
    return out;
}

static SEXP item_array(const mio_format_info *info, int32_t section, int64_t index,
                       int32_t field, int transpose, double shift) {
    const void *data = NULL;
    mio_dtype dtype = MIO_FLOAT64;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_r_check(mio_mdpa_info_array(info, section, index, field, &data, &dtype, &ndim, shape),
                "mdpa_info");
    return array_to_r(data, dtype, ndim, shape, transpose, shift);
}

/* The key/value entries of one item, as a named list (number, string, or a
 * (rows, cols) matrix for an inline table). */
static SEXP item_data(const mio_format_info *info, int32_t section, int64_t index) {
    int64_t n = mio_mdpa_info_data_count(info, section, index);
    if (n < 0) mio_r_fail("mdpa_info");
    SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)n));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)n));
    for (int64_t e = 0; e < n; ++e) {
        mdpa_ctx c = {info, section, index, e, 0};
        SEXP key = PROTECT(mio_r_getstring(get_data_string, &c, "mdpa_info"));
        SET_STRING_ELT(names, (R_xlen_t)e, STRING_ELT(key, 0));
        const int32_t kind = mio_mdpa_info_data_kind(info, section, index, e);
        if (kind < 0) mio_r_fail("mdpa_info");
        SEXP value;
        if (kind == MIO_MDPA_VALUE_TEXT) {
            c.field = 1;
            value = PROTECT(mio_r_getstring(get_data_string, &c, "mdpa_info"));
        } else {
            const void *data = NULL;
            mio_dtype dtype = MIO_FLOAT64;
            int32_t ndim = 0;
            int64_t shape[MIO_MAX_NDIM] = {0};
            mio_r_check(mio_mdpa_info_data_array(info, section, index, e, &data, &dtype, &ndim,
                                                 shape),
                        "mdpa_info");
            value = PROTECT(array_to_r(data, dtype, ndim, shape, 1, 0.0));
        }
        SET_VECTOR_ELT(out, (R_xlen_t)e, value);
        UNPROTECT(2); /* value, key */
    }
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}

/* A list of `n` items, each built by `make(info, i)`. */
typedef SEXP (*item_maker)(const mio_format_info *info, int64_t i);

static SEXP items(const mio_format_info *info, int32_t section, item_maker make) {
    const int64_t n = count_of(info, section);
    SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)n));
    for (int64_t i = 0; i < n; ++i) SET_VECTOR_ELT(out, (R_xlen_t)i, make(info, i));
    UNPROTECT(1);
    return out;
}

static SEXP make_property(const mio_format_info *info, int64_t i) {
    SEXP v[2];
    v[0] = PROTECT(Rf_ScalarReal((double)item_int(info, MIO_MDPA_PROPERTIES, i, 0)));
    v[1] = PROTECT(item_data(info, MIO_MDPA_PROPERTIES, i));
    const char *names[2] = {"id", "values"};
    SEXP out = mio_r_named_list(2, names, v);
    UNPROTECT(2);
    return out;
}

static SEXP make_entity_name(const mio_format_info *info, int64_t i) {
    SEXP v[2];
    v[0] = PROTECT(item_string(info, MIO_MDPA_ENTITY_NAMES, i, 0));
    v[1] = PROTECT(Rf_ScalarLogical(item_int(info, MIO_MDPA_ENTITY_NAMES, i, 0) == 1));
    const char *names[2] = {"name", "is_condition"};
    SEXP out = mio_r_named_list(2, names, v);
    UNPROTECT(2);
    return out;
}

static SEXP make_skipped(const mio_format_info *info, int64_t i) {
    return item_string(info, MIO_MDPA_SKIPPED, i, 0);
}

static SEXP make_table(const mio_format_info *info, int64_t i) {
    SEXP v[2];
    v[0] = PROTECT(item_string(info, MIO_MDPA_TABLES, i, 0));
    v[1] = PROTECT(item_array(info, MIO_MDPA_TABLES, i, 0, 1, 0.0));
    const char *names[2] = {"header", "values"};
    SEXP out = mio_r_named_list(2, names, v);
    UNPROTECT(2);
    return out;
}

static SEXP make_geometry(const mio_format_info *info, int64_t i) {
    SEXP v[4];
    v[0] = PROTECT(item_string(info, MIO_MDPA_GEOMETRIES, i, 0));
    v[1] = PROTECT(item_string(info, MIO_MDPA_GEOMETRIES, i, 1));
    v[2] = PROTECT(item_array(info, MIO_MDPA_GEOMETRIES, i, 0, 0, 1.0)); /* 1-based */
    v[3] = PROTECT(item_array(info, MIO_MDPA_GEOMETRIES, i, 1, 0, 0.0));
    const char *names[4] = {"name", "type", "connectivity", "ids"};
    SEXP out = mio_r_named_list(4, names, v);
    UNPROTECT(4);
    return out;
}

static SEXP make_mesh_block(const mio_format_info *info, int64_t i) {
    SEXP v[5];
    v[0] = PROTECT(Rf_ScalarReal((double)item_int(info, MIO_MDPA_MESH_BLOCKS, i, 0)));
    v[1] = PROTECT(item_data(info, MIO_MDPA_MESH_BLOCKS, i));
    v[2] = PROTECT(item_array(info, MIO_MDPA_MESH_BLOCKS, i, 0, 0, 1.0)); /* 1-based */
    v[3] = PROTECT(item_array(info, MIO_MDPA_MESH_BLOCKS, i, 1, 0, 0.0));
    v[4] = PROTECT(item_array(info, MIO_MDPA_MESH_BLOCKS, i, 2, 0, 0.0));
    const char *names[5] = {"id", "data", "nodes", "element_ids", "condition_ids"};
    SEXP out = mio_r_named_list(5, names, v);
    UNPROTECT(5);
    return out;
}

static SEXP make_submodelpart(const mio_format_info *info, int64_t i) {
    SEXP v[5];
    v[0] = PROTECT(item_string(info, MIO_MDPA_SUBMODELPARTS, i, 0));
    v[1] = PROTECT(item_data(info, MIO_MDPA_SUBMODELPARTS, i));
    v[2] = PROTECT(item_array(info, MIO_MDPA_SUBMODELPARTS, i, 0, 0, 0.0));
    v[3] = PROTECT(item_array(info, MIO_MDPA_SUBMODELPARTS, i, 1, 0, 0.0));
    v[4] = PROTECT(item_array(info, MIO_MDPA_SUBMODELPARTS, i, 2, 0, 0.0));
    const char *names[5] = {"name", "data", "tables", "geometry_ids", "constraint_ids"};
    SEXP out = mio_r_named_list(5, names, v);
    UNPROTECT(5);
    return out;
}

static SEXP make_raw_block(const mio_format_info *info, int64_t i) {
    SEXP v[3];
    v[0] = PROTECT(item_string(info, MIO_MDPA_RAW_BLOCKS, i, 0));
    v[1] = PROTECT(item_string(info, MIO_MDPA_RAW_BLOCKS, i, 1));
    v[2] = PROTECT(item_string(info, MIO_MDPA_RAW_BLOCKS, i, 2));
    const char *names[3] = {"header", "body", "terminator"};
    SEXP out = mio_r_named_list(3, names, v);
    UNPROTECT(3);
    return out;
}

SEXP R_mio_mdpa_info(SEXP x) {
    const mio_format_info *info = info_of(x, 0);
    SEXP v[9];
    v[0] = PROTECT(items(info, MIO_MDPA_PROPERTIES, make_property));
    v[1] = PROTECT(items(info, MIO_MDPA_ENTITY_NAMES, make_entity_name));
    v[2] = PROTECT(items(info, MIO_MDPA_SKIPPED, make_skipped));
    v[3] = PROTECT(item_data(info, MIO_MDPA_MODEL_PART_DATA, 0));
    v[4] = PROTECT(items(info, MIO_MDPA_TABLES, make_table));
    v[5] = PROTECT(items(info, MIO_MDPA_GEOMETRIES, make_geometry));
    v[6] = PROTECT(items(info, MIO_MDPA_MESH_BLOCKS, make_mesh_block));
    v[7] = PROTECT(items(info, MIO_MDPA_SUBMODELPARTS, make_submodelpart));
    v[8] = PROTECT(items(info, MIO_MDPA_RAW_BLOCKS, make_raw_block));
    const char *names[9] = {"properties",    "entity_names", "skipped",
                            "model_part_data", "tables",       "geometries",
                            "mesh_blocks",   "submodelparts", "raw_blocks"};
    SEXP out = mio_r_named_list(9, names, v);
    UNPROTECT(9);
    return out;
}

static SEXP gmsh_array(const mio_format_info *info, int32_t section, int64_t i, int32_t field) {
    const void *data = NULL;
    mio_dtype dtype = MIO_FLOAT64;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_r_check(mio_gmsh_info_array(info, section, i, field, &data, &dtype, &ndim, shape),
                "gmsh_info");
    return array_to_r(data, dtype, ndim, shape, 0, field == 2 ? 1.0 : 0.0);
}

SEXP R_mio_gmsh_info(SEXP x) {
    const mio_format_info *info = info_of(x, 0);
    int64_t nb = mio_gmsh_info_count(info, MIO_GMSH_BOUNDING_ENTITIES);
    int64_t np = mio_gmsh_info_count(info, MIO_GMSH_PERIODIC);
    if (nb < 0 || np < 0) mio_r_fail("gmsh_info");
    SEXP v[2];
    v[0] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nb));
    v[1] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)np));
    for (int64_t i = 0; i < nb; ++i)
        SET_VECTOR_ELT(v[0], i, gmsh_array(info, MIO_GMSH_BOUNDING_ENTITIES, i, 0));
    for (int64_t i = 0; i < np; ++i) {
        SEXP fields[3];
        for (int32_t f = 0; f < 3; ++f)
            fields[f] = PROTECT(gmsh_array(info, MIO_GMSH_PERIODIC, i, f));
        const char *names[3] = {"entity", "affine", "node_pairs"};
        SET_VECTOR_ELT(v[1], i, mio_r_named_list(3, names, fields));
        UNPROTECT(3);
    }
    const char *names[2] = {"bounding_entities", "periodic"};
    SEXP out = mio_r_named_list(2, names, v);
    UNPROTECT(2);
    return out;
}

/* The NURBS side channel of an MFEM NURBS mesh (v16.40.0), copied: per patch
 * its degrees, knot vectors, control net (space dimension x points), weights
 * and field control values; the field names; the patch topology; the cell map.
 * Vertices, patches and knot spans come back 1-based; attributes and the
 * edges' knot-vector ids as the file spelled them. */
static SEXP nurbs_array(const mio_format_info *info, int32_t section, int64_t i, int32_t field,
                        double shift) {
    const void *data = NULL;
    mio_dtype dtype = MIO_FLOAT64;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_r_check(mio_nurbs_info_array(info, section, i, field, &data, &dtype, &ndim, shape),
                "nurbs_info");
    return array_to_r(data, dtype, ndim, shape, 0, shift);
}

static SEXP nurbs_int(const mio_format_info *info, int32_t section, int64_t i, int32_t field) {
    int64_t v = 0;
    mio_r_check(mio_nurbs_info_int(info, section, i, field, &v), "nurbs_info");
    return Rf_ScalarReal((double)v);
}

SEXP R_mio_nurbs_info(SEXP x) {
    const mio_format_info *info = info_of(x, 0);
    const int64_t np = mio_nurbs_info_count(info, MIO_NURBS_PATCHES);
    const int64_t nf = mio_nurbs_info_count(info, MIO_NURBS_FIELDS);
    const int64_t nt = mio_nurbs_info_count(info, MIO_NURBS_TOPOLOGY);
    const int64_t nc = mio_nurbs_info_count(info, MIO_NURBS_CELL_MAP);
    if (np < 0 || nf < 0 || nt < 0 || nc < 0) mio_r_fail("nurbs_info");
    SEXP v[6];
    v[0] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)np));
    for (int64_t p = 0; p < np; ++p) {
        int64_t dim = 0;
        mio_r_check(mio_nurbs_info_int(info, MIO_NURBS_PATCHES, p, 0, &dim), "nurbs_info");
        SEXP f[5];
        f[0] = PROTECT(nurbs_array(info, MIO_NURBS_PATCHES, p, 0, 0.0));
        f[1] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)dim));
        for (int64_t d = 0; d < dim; ++d)
            SET_VECTOR_ELT(f[1], d, nurbs_array(info, MIO_NURBS_PATCHES, p, (int32_t)(3 + d), 0.0));
        f[2] = PROTECT(nurbs_array(info, MIO_NURBS_PATCHES, p, 1, 0.0));
        f[3] = PROTECT(nurbs_array(info, MIO_NURBS_PATCHES, p, 2, 0.0));
        f[4] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nf));
        for (int64_t k = 0; k < nf; ++k)
            SET_VECTOR_ELT(f[4], k, nurbs_array(info, MIO_NURBS_PATCHES, p, (int32_t)(6 + k), 0.0));
        const char *names[5] = {"degrees", "knots", "control_points", "weights", "fields"};
        SET_VECTOR_ELT(v[0], p, mio_r_named_list(5, names, f));
        UNPROTECT(5);
    }
    v[1] = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)nf));
    for (int64_t k = 0; k < nf; ++k) {
        const int64_t n = mio_nurbs_info_string(info, MIO_NURBS_FIELDS, k, 0, NULL, 0);
        if (n < 0) mio_r_fail("nurbs_info");
        char *buf = R_alloc((size_t)n + 1, 1);
        mio_nurbs_info_string(info, MIO_NURBS_FIELDS, k, 0, buf, n + 1);
        SET_STRING_ELT(v[1], k, Rf_mkCharLen(buf, (int)n));
    }
    if (nt > 0) {
        SEXP t[7];
        t[0] = PROTECT(nurbs_int(info, MIO_NURBS_TOPOLOGY, 0, 0));
        t[1] = PROTECT(nurbs_int(info, MIO_NURBS_TOPOLOGY, 0, 1));
        t[2] = PROTECT(nurbs_array(info, MIO_NURBS_TOPOLOGY, 0, 0, 1.0));
        t[3] = PROTECT(nurbs_array(info, MIO_NURBS_TOPOLOGY, 0, 1, 0.0));
        t[4] = PROTECT(nurbs_array(info, MIO_NURBS_TOPOLOGY, 0, 2, 1.0));
        t[5] = PROTECT(nurbs_array(info, MIO_NURBS_TOPOLOGY, 0, 3, 0.0));
        t[6] = PROTECT(nurbs_array(info, MIO_NURBS_TOPOLOGY, 0, 4, 0.0));
        if (Rf_isMatrix(t[6])) { /* (knot vector, vertex, vertex) per column */
            double *e = REAL(t[6]);
            for (R_xlen_t c = 0; c < Rf_ncols(t[6]); ++c) {
                e[3 * c + 1] += 1.0;
                e[3 * c + 2] += 1.0;
            }
        }
        const char *names[7] = {"dimension",         "num_vertices",        "patch_vertices",
                                "patch_attributes",  "boundary_vertices",   "boundary_attributes",
                                "edges"};
        v[2] = mio_r_named_list(7, names, t);
        UNPROTECT(7);
    } else {
        v[2] = R_NilValue;
    }
    PROTECT(v[2]);
    for (int32_t f = 0; f < 3; ++f)
        v[3 + f] = nc > 0 ? nurbs_array(info, MIO_NURBS_CELL_MAP, 0, f, 1.0)
                          : Rf_allocVector(REALSXP, 0);
    for (int32_t f = 0; f < 3; ++f) PROTECT(v[3 + f]);
    const char *names[6] = {"patches",    "field_names", "topology",
                            "cell_patch", "cell_span",   "boundary_cell_patch"};
    SEXP out = mio_r_named_list(6, names, v);
    UNPROTECT(6);
    return out;
}

/* The isogeometric layer (v16.43.0): refine, sample and evaluate the NURBS
 * model in a handle without going through its sampling. */
static void iga_int3(SEXP x, const char *what, int32_t out[3]) {
    SEXP v = PROTECT(Rf_coerceVector(x, INTSXP));
    const R_xlen_t n = Rf_xlength(v);
    if (n < 1 || n > 3) Rf_error("'%s' must be one integer or one per direction (at most 3)", what);
    for (int d = 0; d < 3; ++d) {
        const int e = INTEGER(v)[d < n ? d : n - 1];
        if (e == NA_INTEGER) Rf_error("'%s' must not contain NA", what);
        out[d] = (int32_t)e;
    }
    UNPROTECT(1);
}

static mio_iga_opts iga_opts(SEXP subdivisions, SEXP elevate, SEXP sample_subdivisions,
                             SEXP trim_samples, SEXP edges, SEXP trim) {
    mio_iga_opts opts;
    mio_iga_opts_init(&opts);
    iga_int3(subdivisions, "subdivisions", opts.subdivisions);
    iga_int3(elevate, "elevate", opts.elevate);
    opts.sample_subdivisions = (int32_t)mio_r_int(sample_subdivisions, "sample_subdivisions");
    opts.trim_samples = (int32_t)mio_r_int(trim_samples, "trim_samples");
    opts.edges = mio_r_bool(edges, "edges") ? 1 : 0;
    opts.trim = mio_r_bool(trim, "trim") ? 1 : 0;
    return opts;
}

static SEXP iga_mesh_and_info(mio_mesh *mesh, mio_format_info *info, const char *what) {
    if (mesh == NULL) mio_r_fail(what);
    SEXP values[2];
    values[0] = PROTECT(mio_r_wrap_mesh(mesh));
    values[1] = info != NULL ? wrap_info(info) : R_NilValue;
    PROTECT(values[1]);
    const char *names[2] = {"mesh", "info"};
    SEXP out = mio_r_named_list(2, names, values);
    UNPROTECT(2);
    return out;
}

SEXP R_mio_iga_refine(SEXP x, SEXP subdivisions, SEXP elevate, SEXP sample_subdivisions,
                      SEXP trim_samples, SEXP edges, SEXP trim) {
    const mio_format_info *info = info_of(x, 0);
    const mio_iga_opts opts =
        iga_opts(subdivisions, elevate, sample_subdivisions, trim_samples, edges, trim);
    mio_format_info *refined = NULL;
    mio_mesh *mesh = mio_iga_refine(info, &opts, &refined);
    return iga_mesh_and_info(mesh, refined, "iga_refine");
}

SEXP R_mio_iga_sample(SEXP x, SEXP sample_subdivisions, SEXP trim_samples, SEXP edges,
                      SEXP trim) {
    const mio_format_info *info = info_of(x, 0);
    SEXP one = PROTECT(Rf_ScalarInteger(1));
    SEXP zero = PROTECT(Rf_ScalarInteger(0));
    const mio_iga_opts opts = iga_opts(one, zero, sample_subdivisions, trim_samples, edges, trim);
    UNPROTECT(2);
    mio_format_info *mapped = NULL;
    mio_mesh *mesh = mio_iga_sample(info, &opts, &mapped);
    return iga_mesh_and_info(mesh, mapped, "iga_sample");
}

/* --- mesh generation from CAD (v17.13.0; doc/mesh_cad.md) ----------------- */

static int64_t r_mesh_cad_warnings(void *ctx, char *buf, int64_t buflen) {
    (void)ctx;
    return mio_mesh_cad_surface_warnings(buf, buflen);
}

static int64_t r_mesh_cad_reason(void *ctx, char *buf, int64_t buflen) {
    return mio_mesh_cad_surface_face_reason(*(const int64_t *)ctx, buf, buflen);
}

/* The corners (v17.15.0) of the last call as named column vectors: 1-based
 * face and point, the parameters, the point, the angles, the reach and
 * whether it is a pole. */
static SEXP r_mesh_cad_corners(void) {
    const int64_t n = mio_mesh_cad_surface_num_corners();
    static const char *names[] = {"face", "point", "u", "v", "x", "y", "z",
                                  "angle_deg", "min_forced_deg", "reach", "pole"};
    enum { R_CORNER_N = 11 };
    SEXP cols = PROTECT(Rf_allocVector(VECSXP, R_CORNER_N));
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, R_CORNER_N));
    for (int c = 0; c < R_CORNER_N; ++c) {
        SET_STRING_ELT(nm, c, Rf_mkChar(names[c]));
        SET_VECTOR_ELT(cols, c, Rf_allocVector(c == 10 ? LGLSXP : REALSXP, (R_xlen_t)n));
    }
    for (int64_t k = 0; k < n; ++k) {
        mio_mesh_cad_corner corner;
        if (mio_mesh_cad_surface_corner(k, &corner) != MIO_OK) mio_r_fail("mesh_cad_surface");
        const double v[10] = {(double)corner.face + 1.0, (double)corner.point + 1.0, corner.u,
                              corner.v, corner.x[0], corner.x[1], corner.x[2], corner.angle_deg,
                              corner.min_forced_deg, corner.reach};
        for (int c = 0; c < 10; ++c) REAL(VECTOR_ELT(cols, c))[k] = v[c];
        LOGICAL(VECTOR_ELT(cols, 10))[k] = corner.pole != 0;
    }
    Rf_setAttrib(cols, R_NamesSymbol, nm);
    UNPROTECT(2);
    return cols;
}

/* A named list of equal-length columns. */
static SEXP r_columns(int n, const char** names, SEXP* cols) {
    SEXP out = PROTECT(Rf_allocVector(VECSXP, n));
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, n));
    for (int i = 0; i < n; ++i) {
        SET_VECTOR_ELT(out, i, cols[i]);
        SET_STRING_ELT(nm, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(out, R_NamesSymbol, nm);
    UNPROTECT(2);
    return out;
}

/* A reason string of record `k` through `fn` (the C string rule). */
static SEXP r_reason(int64_t (*fn)(int64_t, char*, int64_t), int64_t k) {
    const int64_t n = fn(k, NULL, 0);
    if (n <= 0)
        return Rf_mkChar("");
    char* buf = (char*)R_alloc((size_t)n + 1, 1);
    fn(k, buf, n + 1);
    return Rf_mkCharLen(buf, (int)n);
}

/* v17.20.0: the last call's reparametrised faces (6.4.4). */
static SEXP r_mesh_cad_reparametrised(void) {
    const R_xlen_t n = (R_xlen_t)mio_mesh_cad_surface_num_reparametrised();
    SEXP face = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP cond = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP disc = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP method = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP reason = PROTECT(Rf_allocVector(STRSXP, n));
    for (R_xlen_t i = 0; i < n; ++i) {
        int64_t f = -1;
        double c = 0.0, d = 0.0;
        int32_t m = 0;
        mio_mesh_cad_surface_reparametrised((int64_t)i, &f, &c, &d, &m);
        REAL(face)[i] = (double)f + 1.0;
        REAL(cond)[i] = c;
        REAL(disc)[i] = d;
        REAL(method)[i] = (double)m;
        SET_STRING_ELT(reason, i, r_reason(mio_mesh_cad_surface_reparam_reason, (int64_t)i));
    }
    const char* names[] = {"face", "condition_number", "condition_number_disc", "method", "reason"};
    SEXP cols[] = {face, cond, disc, method, reason};
    SEXP out = r_columns(5, names, cols);
    UNPROTECT(5);
    return out;
}

/* v17.20.0: the slivers found (6.5.2); host 0 for none. */
static SEXP r_mesh_cad_merged_faces(void) {
    const R_xlen_t n = (R_xlen_t)mio_mesh_cad_surface_num_merged_faces();
    SEXP sliver = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP host = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP width = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP merged = PROTECT(Rf_allocVector(LGLSXP, n));
    SEXP reason = PROTECT(Rf_allocVector(STRSXP, n));
    for (R_xlen_t i = 0; i < n; ++i) {
        int64_t s = -1, h = -1;
        double w = 0.0;
        int32_t m = 0;
        mio_mesh_cad_surface_merged_face((int64_t)i, &s, &h, &w, &m);
        REAL(sliver)[i] = (double)s + 1.0;
        REAL(host)[i] = (double)h + 1.0;
        REAL(width)[i] = w;
        LOGICAL(merged)[i] = m != 0;
        SET_STRING_ELT(reason, i, r_reason(mio_mesh_cad_surface_merged_face_reason, (int64_t)i));
    }
    const char* names[] = {"sliver", "host", "width", "merged", "reason"};
    SEXP cols[] = {sliver, host, width, merged, reason};
    SEXP out = r_columns(5, names, cols);
    UNPROTECT(5);
    return out;
}

/* v17.20.0: where a declared tolerance accepted a gap (6.5.1). */
static SEXP r_mesh_cad_tolerance_uses(void) {
    const R_xlen_t n = (R_xlen_t)mio_mesh_cad_surface_num_tolerance_uses();
    SEXP edge = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP declared = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP gap = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP where = PROTECT(Rf_allocVector(STRSXP, n));
    static const char* kWhere[] = {"closed_start", "corner", "projection", "pairing"};
    for (R_xlen_t i = 0; i < n; ++i) {
        int64_t e = -1;
        double d = 0.0, g = 0.0;
        int32_t w = 0;
        mio_mesh_cad_surface_tolerance_use((int64_t)i, &e, &d, &g, &w);
        REAL(edge)[i] = (double)e + 1.0;
        REAL(declared)[i] = d;
        REAL(gap)[i] = g;
        SET_STRING_ELT(where, i, Rf_mkChar(kWhere[w >= 0 && w < 4 ? w : 2]));
    }
    const char* names[] = {"edge", "declared", "gap", "where"};
    SEXP cols[] = {edge, declared, gap, where};
    SEXP out = r_columns(4, names, cols);
    UNPROTECT(4);
    return out;
}

/* v17.20.0: the corners lopped (6.2.3), 1-based face and point. */
static SEXP r_mesh_cad_lopped_corners(void) {
    const R_xlen_t n = (R_xlen_t)mio_mesh_cad_surface_num_lopped_corners();
    SEXP face = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP point = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP angle = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP lfs0 = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP radius = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP turning = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP bound = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP forced = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP min_forced = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP within = PROTECT(Rf_allocVector(LGLSXP, n));
    SEXP guaranteed = PROTECT(Rf_allocVector(LGLSXP, n));
    SEXP reason = PROTECT(Rf_allocVector(STRSXP, n));
    for (R_xlen_t i = 0; i < n; ++i) {
        int64_t f = -1, p = -1, nf = 0;
        double a = 0.0, l = 0.0, r = 0.0, t = 0.0, b = 0.0, m = 0.0;
        int32_t w = 0, g = 0;
        mio_mesh_cad_surface_lopped_corner((int64_t)i, &f, &p, &a, &l, &r, &t, &b, &nf, &m, &w,
                                           &g);
        REAL(face)[i] = (double)f + 1.0;
        REAL(point)[i] = (double)p + 1.0;
        REAL(angle)[i] = a;
        REAL(lfs0)[i] = l;
        REAL(radius)[i] = r;
        REAL(turning)[i] = t;
        REAL(bound)[i] = b;
        REAL(forced)[i] = (double)nf;
        REAL(min_forced)[i] = m;
        LOGICAL(within)[i] = w != 0;
        LOGICAL(guaranteed)[i] = g != 0;
        SET_STRING_ELT(reason, i,
                       r_reason(mio_mesh_cad_surface_lopped_corner_reason, (int64_t)i));
    }
    const char* names[] = {"face",          "point",          "angle_deg",   "lfs0",
                           "lop_radius",    "turning_deg",    "stated_bound_deg",
                           "num_forced",    "min_forced_deg", "within_two_radii",
                           "guaranteed",    "reason"};
    SEXP cols[] = {face,   point,      angle,  lfs0,   radius,     turning,
                   bound,  forced,     min_forced, within, guaranteed, reason};
    SEXP out = r_columns(12, names, cols);
    UNPROTECT(12);
    return out;
}

/* project_to_cad (roadmap 6.10.1): every point of a mesh onto the CAD face or
 * edge its cells name. Returns list(mesh, num_projected, ..., outside_trim). */
SEXP R_mio_project_to_cad(SEXP mesh, SEXP x, SEXP face_array, SEXP edge_array, SEXP max_distance,
                          SEXP trim_samples) {
    const mio_format_info* info = info_of(x, 0);
    mio_project_to_cad_opts opts;
    mio_project_to_cad_opts_init(&opts);
    opts.face_array = mio_r_string(face_array, "face_array");
    opts.edge_array = mio_r_string(edge_array, "edge_array");
    opts.max_distance = mio_r_double(max_distance, "max_distance");
    opts.trim_samples = (int32_t)mio_r_int(trim_samples, "trim_samples");
    mio_project_to_cad_report report;
    mio_mesh* out = mio_project_to_cad(mio_r_mesh(mesh), info, &opts, &report);
    if (out == NULL)
        mio_r_fail("project_to_cad");
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    const int64_t n = mio_project_to_cad_outside_trim(NULL, 0);
    SEXP outside = PROTECT(Rf_allocVector(REALSXP, n > 0 ? (R_xlen_t)n : 0));
    if (n > 0) {
        int64_t* buf = (int64_t*)R_alloc((size_t)n, sizeof(int64_t));
        mio_project_to_cad_outside_trim(buf, n);
        for (int64_t i = 0; i < n; ++i)
            REAL(outside)[i] = (double)buf[i] + 1.0;
    }
    SEXP v1 = PROTECT(Rf_ScalarReal((double)report.num_projected));
    SEXP v2 = PROTECT(Rf_ScalarReal((double)report.num_on_edges));
    SEXP v3 = PROTECT(Rf_ScalarReal((double)report.num_corners));
    SEXP v4 = PROTECT(Rf_ScalarReal((double)report.num_ambiguous));
    SEXP v5 = PROTECT(Rf_ScalarReal((double)report.num_too_far));
    SEXP v6 = PROTECT(Rf_ScalarReal((double)report.num_failed));
    SEXP v7 = PROTECT(Rf_ScalarReal((double)report.num_outside_trim));
    SEXP v8 = PROTECT(Rf_ScalarReal(report.max_move));
    SEXP v9 = PROTECT(Rf_ScalarReal((double)report.num_on_curves));
    const char* names[] = {"mesh",          "num_projected", "num_on_edges", "num_on_curves",
                           "num_corners",   "num_ambiguous", "num_too_far",  "num_failed",
                           "num_outside_trim", "max_move",   "outside_trim"};
    SEXP values[] = {mo, v1, v2, v9, v3, v4, v5, v6, v7, v8, outside};
    SEXP res = PROTECT(mio_r_named_list(11, names, values));
    UNPROTECT(12);
    return res;
}

SEXP R_mio_mesh_cad_surface(
    SEXP x, SEXP target_size, SEXP max_chordal_deviation, SEXP min_points_per_edge, SEXP tolerance,
    SEXP record_face_ids, SEXP edges, SEXP match_unreferenced_trims, SEXP strict, SEXP max_points,
    SEXP min_angle, SEXP smooth_iterations, SEXP smooth_method, SEXP order, SEXP grading,
    SEXP curvature_deviation, SEXP cells_across_gaps, SEXP pole_rings, SEXP face_size_faces,
    SEXP face_sizes, SEXP edge_size_edges, SEXP edge_sizes, SEXP background, SEXP background_array,
    SEXP metric, SEXP virtual_edges, SEXP reparam_condition_number, SEXP model_tolerances,
    SEXP virtual_faces, SEXP lop_corners, SEXP metric_background, SEXP metric_array,
    SEXP metric_error, SEXP untangle, SEXP quads, SEXP quad_pairs) {
    const mio_format_info *info = info_of(x, 0);
    mio_mesh_cad_surface_opts opts;
    mio_mesh_cad_surface_opts_init(&opts);
    opts.target_size = mio_r_double(target_size, "target_size");
    opts.max_chordal_deviation = mio_r_double(max_chordal_deviation, "max_chordal_deviation");
    opts.min_points_per_edge = mio_r_int64(min_points_per_edge, "min_points_per_edge");
    opts.tolerance = mio_r_double(tolerance, "tolerance");
    opts.record_face_ids = mio_r_bool(record_face_ids, "record_face_ids") ? 1 : 0;
    opts.edges = mio_r_bool(edges, "edges") ? 1 : 0;
    opts.match_unreferenced_trims =
        mio_r_bool(match_unreferenced_trims, "match_unreferenced_trims") ? 1 : 0;
    opts.strict = mio_r_bool(strict, "strict") ? 1 : 0;
    opts.max_points = mio_r_int64(max_points, "max_points");
    opts.min_angle_deg = mio_r_double(min_angle, "min_angle");
    opts.smooth_iterations = (int32_t)mio_r_int64(smooth_iterations, "smooth_iterations");
    opts.smooth_method = (int32_t)mio_r_int64(smooth_method, "smooth_method");
    opts.order = (int32_t)mio_r_int64(order, "order");
    opts.grading = mio_r_double(grading, "grading");
    opts.curvature_deviation = mio_r_double(curvature_deviation, "curvature_deviation");
    opts.cells_across_gaps = mio_r_double(cells_across_gaps, "cells_across_gaps");
    opts.pole_rings = (int32_t)mio_r_int64(pole_rings, "pole_rings");
    /* v17.19.0: given sizes (1-based indices from R), the background, the
     * metric and virtual edges. */
    mio_mesh_cad_surface_opts_ex ex;
    mio_mesh_cad_surface_opts_ex_init(&ex);
    ex.metric = (int32_t)mio_r_int64(metric, "metric");
    ex.virtual_edges = mio_r_bool(virtual_edges, "virtual_edges") ? 1 : 0;
    /* v17.20.0: the reparametrisation, model tolerances, virtual faces and
     * corner lopping. */
    ex.reparam_condition_number =
        mio_r_double(reparam_condition_number, "reparam_condition_number");
    ex.model_tolerances = mio_r_bool(model_tolerances, "model_tolerances") ? 1 : 0;
    ex.virtual_faces = mio_r_bool(virtual_faces, "virtual_faces") ? 1 : 0;
    ex.lop_corners = mio_r_bool(lop_corners, "lop_corners") ? 1 : 0;
    const R_xlen_t nfs = Rf_xlength(face_size_faces), nes = Rf_xlength(edge_size_edges);
    if (Rf_xlength(face_sizes) != nfs || Rf_xlength(edge_sizes) != nes)
        Rf_error("mesh_cad_surface: the size indices and sizes must have the same length");
    int64_t* fidx = (int64_t*)R_alloc(nfs > 0 ? (size_t)nfs : 1, sizeof(int64_t));
    double* fsz = (double*)R_alloc(nfs > 0 ? (size_t)nfs : 1, sizeof(double));
    int64_t* eidx = (int64_t*)R_alloc(nes > 0 ? (size_t)nes : 1, sizeof(int64_t));
    double* esz = (double*)R_alloc(nes > 0 ? (size_t)nes : 1, sizeof(double));
    for (R_xlen_t i = 0; i < nfs; ++i) {
        /* Numeric vectors: the R wrapper coerces them. */
        fidx[i] = (int64_t)REAL(face_size_faces)[i] - 1;
        fsz[i] = REAL(face_sizes)[i];
    }
    for (R_xlen_t i = 0; i < nes; ++i) {
        eidx[i] = (int64_t)REAL(edge_size_edges)[i] - 1;
        esz[i] = REAL(edge_sizes)[i];
    }
    ex.face_size_faces = fidx;
    ex.face_sizes = fsz;
    ex.num_face_sizes = (int64_t)nfs;
    ex.edge_size_edges = eidx;
    ex.edge_sizes = esz;
    ex.num_edge_sizes = (int64_t)nes;
    if (!Rf_isNull(background))
        ex.background = mio_r_mesh(background);
    ex.background_array = mio_r_string(background_array, "background_array");
    /* Roadmap 6.10.2: a metric from a field's Hessian. */
    if (!Rf_isNull(metric_background))
        ex.metric_background = mio_r_mesh(metric_background);
    ex.metric_array = mio_r_string(metric_array, "metric_array");
    ex.metric_error = mio_r_double(metric_error, "metric_error");
    /* Roadmap 6.6.3 and 6.8: untangling and quads (a mio_mesh_cad_quads). */
    ex.untangle = mio_r_bool(untangle, "untangle") ? 1 : 0;
    ex.quads = (int32_t)mio_r_int64(quads, "quads");
    ex.quad_pairs = (int32_t)mio_r_int64(quad_pairs, "quad_pairs");
    mio_mesh_cad_surface_report report;
    mio_mesh* out = mio_mesh_cad_surface2(info, &opts, &ex, &report);
    if (out == NULL) mio_r_fail("mesh_cad_surface");

    /* Wrap the mesh first: once it has a finalizer nothing below can leak it. */
    SEXP mo = PROTECT(mio_r_wrap_mesh(out));
    SEXP warnings = PROTECT(mio_r_getstring(r_mesh_cad_warnings, NULL, "mesh_cad_surface"));
    static const char *statuses[] = {"ok",           "empty",           "open_loop",
                                     "self_intersecting", "degenerate", "projection_failed",
                                     "budget_exceeded", "unsupported",  "other"};
    const R_xlen_t nf = (R_xlen_t)report.num_faces;
    SEXP status = PROTECT(Rf_allocVector(STRSXP, nf));
    SEXP reason = PROTECT(Rf_allocVector(STRSXP, nf));
    for (R_xlen_t f = 0; f < nf; ++f) {
        const int32_t st = mio_mesh_cad_surface_face_status((int64_t)f);
        SET_STRING_ELT(status, f, Rf_mkChar(st >= 0 && st <= 8 ? statuses[st] : "other"));
        int64_t face = (int64_t)f;
        SEXP why = PROTECT(mio_r_getstring(r_mesh_cad_reason, &face, "mesh_cad_surface"));
        SET_STRING_ELT(reason, f, STRING_ELT(why, 0));
        UNPROTECT(1);
    }

    SEXP histogram = PROTECT(Rf_allocVector(REALSXP, MIO_MESH_CAD_ANGLE_BINS));
    {
        int64_t bins[MIO_MESH_CAD_ANGLE_BINS] = {0};
        mio_mesh_cad_surface_angle_histogram(bins, MIO_MESH_CAD_ANGLE_BINS);
        for (int b = 0; b < MIO_MESH_CAD_ANGLE_BINS; ++b) REAL(histogram)[b] = (double)bins[b];
    }
    SEXP corners = PROTECT(r_mesh_cad_corners());
    /* v17.19.0: per face, the 1-based edges that limited its given size; and
     * each virtual edge's chain of 1-based edges. */
    SEXP limits = PROTECT(Rf_allocVector(VECSXP, nf));
    for (R_xlen_t f = 0; f < nf; ++f) {
        const int64_t n = mio_mesh_cad_surface_face_size_limits((int64_t)f, NULL, 0);
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n > 0 ? (R_xlen_t)n : 0));
        if (n > 0) {
            int64_t* buf = (int64_t*)R_alloc((size_t)n, sizeof(int64_t));
            mio_mesh_cad_surface_face_size_limits((int64_t)f, buf, n);
            for (int64_t i = 0; i < n; ++i)
                REAL(v)[i] = (double)buf[i] + 1.0;
        }
        SET_VECTOR_ELT(limits, f, v);
        UNPROTECT(1);
    }
    const int64_t nv = mio_mesh_cad_surface_num_virtual_edges();
    SEXP virtual_edges_out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nv));
    for (int64_t k = 0; k < nv; ++k) {
        const int64_t n = mio_mesh_cad_surface_virtual_edge(k, NULL, 0);
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n > 0 ? (R_xlen_t)n : 0));
        if (n > 0) {
            int64_t* buf = (int64_t*)R_alloc((size_t)n, sizeof(int64_t));
            mio_mesh_cad_surface_virtual_edge(k, buf, n);
            for (int64_t i = 0; i < n; ++i)
                REAL(v)[i] = (double)buf[i] + 1.0;
        }
        SET_VECTOR_ELT(virtual_edges_out, (R_xlen_t)k, v);
        UNPROTECT(1);
    }

    /* v17.20.0: the faces over the reparametrisation threshold, the slivers
     * found and the declared tolerances used, as named columns (1-based). */
    SEXP reparam = PROTECT(r_mesh_cad_reparametrised());
    SEXP merged = PROTECT(r_mesh_cad_merged_faces());
    SEXP tol_uses = PROTECT(r_mesh_cad_tolerance_uses());
    SEXP lopped = PROTECT(r_mesh_cad_lopped_corners());

    int64_t num_quads = 0, num_leftover = 0;
    double quad_sj = 0.0;
    mio_mesh_cad_surface_quad_stats(&num_quads, &num_leftover, &quad_sj);

    enum { R_MESH_CAD_N = 38 };
    SEXP res = PROTECT(Rf_allocVector(VECSXP, R_MESH_CAD_N));
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, R_MESH_CAD_N));
    int k = 0;
#define R_MESH_CAD_SET(name, value)              \
    do {                                         \
        SET_STRING_ELT(nm, k, Rf_mkChar(name));  \
        SET_VECTOR_ELT(res, k, (value));         \
        ++k;                                     \
    } while (0)
#define R_MESH_CAD_NUM(f) R_MESH_CAD_SET(#f, Rf_ScalarReal((double)report.f))
    R_MESH_CAD_SET("mesh", mo);
    R_MESH_CAD_NUM(num_points);
    R_MESH_CAD_NUM(num_triangles);
    R_MESH_CAD_NUM(num_edge_segments);
    R_MESH_CAD_NUM(num_faces);
    R_MESH_CAD_NUM(num_failed);
    R_MESH_CAD_NUM(num_edges);
    R_MESH_CAD_NUM(num_open_edges);
    R_MESH_CAD_NUM(num_non_manifold_edges);
    R_MESH_CAD_NUM(num_misoriented_edges);
    R_MESH_CAD_NUM(num_degenerate_dropped);
    R_MESH_CAD_NUM(tolerance);
    R_MESH_CAD_NUM(min_angle_deg);
    R_MESH_CAD_NUM(max_angle_deg);
    R_MESH_CAD_NUM(min_edge);
    R_MESH_CAD_NUM(max_edge);
    R_MESH_CAD_NUM(max_edge_gap);
    R_MESH_CAD_SET("face_status", status);
    R_MESH_CAD_SET("face_reason", reason);
    /* One string, newline-separated; the R wrapper splits it. */
    R_MESH_CAD_SET("warnings", warnings);
    R_MESH_CAD_NUM(num_refinement_points);
    R_MESH_CAD_NUM(num_faces_below_bound);
    R_MESH_CAD_NUM(num_smoothed);
    R_MESH_CAD_NUM(min_input_angle_deg);
    R_MESH_CAD_SET("angle_histogram", histogram);
    /* Named columns; the R wrapper makes a data frame. */
    R_MESH_CAD_SET("corners", corners);
    R_MESH_CAD_NUM(num_curved_straightened);
    R_MESH_CAD_SET("face_size_limits", limits);
    R_MESH_CAD_SET("virtual_edges", virtual_edges_out);
    R_MESH_CAD_SET("reparametrised", reparam);
    R_MESH_CAD_SET("merged_faces", merged);
    R_MESH_CAD_SET("tolerance_uses", tol_uses);
    R_MESH_CAD_SET("lopped_corners", lopped);
    R_MESH_CAD_SET("num_untangled", Rf_ScalarReal((double)mio_mesh_cad_surface_num_untangled()));
    R_MESH_CAD_SET("num_quads", Rf_ScalarReal((double)num_quads));
    R_MESH_CAD_SET("num_leftover_triangles", Rf_ScalarReal((double)num_leftover));
    R_MESH_CAD_SET("min_quad_scaled_jacobian", Rf_ScalarReal(quad_sj));
    R_MESH_CAD_SET("min_metric_angle", Rf_ScalarReal(mio_mesh_cad_surface_min_metric_angle()));
#undef R_MESH_CAD_NUM
#undef R_MESH_CAD_SET
    if (k != R_MESH_CAD_N)
        Rf_error("mesh_cad_surface: report has %d fields, expected %d", k, (int)R_MESH_CAD_N);
    Rf_setAttrib(res, R_NamesSymbol, nm);
    UNPROTECT(14);
    return res;
}

/* One array of an evaluation as an R array with its dimensions reversed (the
 * C layout (N, d, sdim) is R's (sdim, d, N)); R_NilValue when it is absent. */
static SEXP iga_array(const mio_iga_evaluation *ev, int32_t which, int64_t index, int as_logical) {
    const void *data = NULL;
    mio_dtype dtype = MIO_FLOAT64;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    if (mio_iga_evaluation_array(ev, which, index, &data, &dtype, &ndim, shape) != MIO_OK)
        return R_NilValue;
    R_xlen_t n = 1;
    for (int32_t k = 0; k < ndim; ++k) n *= (R_xlen_t)shape[k];
    SEXP flat = PROTECT(mio_r_copy_as_real(data, dtype, n));
    SEXP out = flat;
    if (as_logical) {
        out = PROTECT(Rf_coerceVector(flat, LGLSXP));
        UNPROTECT(2);
        return out;
    }
    if (ndim >= 2) {
        SEXP dims = PROTECT(Rf_allocVector(INTSXP, ndim));
        for (int32_t k = 0; k < ndim; ++k) INTEGER(dims)[k] = (int)shape[ndim - 1 - k];
        Rf_setAttrib(flat, R_DimSymbol, dims);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return out;
}

SEXP R_mio_iga_evaluate(SEXP x, SEXP patch, SEXP params, SEXP derivatives, SEXP trim_samples) {
    const mio_format_info *info = info_of(x, 0);
    if (!Rf_isMatrix(params)) Rf_error("'params' must be a (d, N) matrix, one column per point");
    SEXP p = PROTECT(Rf_coerceVector(params, REALSXP));
    const int64_t n = (int64_t)Rf_ncols(params);
    mio_iga_evaluation *ev =
        mio_iga_evaluate(info, (int64_t)mio_r_int(patch, "patch") - 1, REAL(p), n,
                         (int32_t)mio_r_int(derivatives, "derivatives"),
                         (int32_t)mio_r_int(trim_samples, "trim_samples"));
    UNPROTECT(1);
    if (ev == NULL) mio_r_fail("iga_evaluate");
    const int64_t nf = mio_iga_evaluation_num_fields(ev);
    SEXP v[5];
    v[0] = PROTECT(iga_array(ev, MIO_IGA_POINTS, 0, 0));
    v[1] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nf));
    for (int64_t k = 0; k < nf; ++k) SET_VECTOR_ELT(v[1], k, iga_array(ev, MIO_IGA_FIELD, k, 0));
    v[2] = PROTECT(iga_array(ev, MIO_IGA_DERIVATIVES, 0, 0));
    v[3] = PROTECT(iga_array(ev, MIO_IGA_NORMALS, 0, 0));
    v[4] = PROTECT(iga_array(ev, MIO_IGA_INSIDE, 0, 1));
    mio_iga_evaluation_free(ev);
    const char *names[5] = {"points", "fields", "derivatives", "normals", "inside"};
    SEXP out = mio_r_named_list(5, names, v);
    UNPROTECT(5);
    return out;
}

/* The B-rep side channel of a Kratos CAD JSON (or OpenCASCADE STEP/IGES)
 * handle, copied: patches, faces with their loops and trims, edges and
 * vertices with their topology, the brep groups and the file's other members.
 * Nested lists come from the C API's flat sections (each owner names its first
 * item and count). Indices (a face's patch, a brep's lists) come back 1-based;
 * ids, trim indices and flags as the file spelled them. */
typedef struct {
    const mio_format_info *info;
    int32_t section;
    int64_t index;
    int32_t field;
} brep_ctx;

static int64_t brep_get_string(void *ctx, char *buf, int64_t buflen) {
    brep_ctx *c = (brep_ctx *)ctx;
    return mio_brep_info_string(c->info, c->section, c->index, c->field, buf, buflen);
}

static SEXP brep_str(const mio_format_info *info, int32_t section, int64_t i, int32_t field) {
    brep_ctx c = {info, section, i, field};
    return mio_r_getstring(brep_get_string, &c, "brep_info");
}

static int64_t brep_ival(const mio_format_info *info, int32_t section, int64_t i, int32_t field) {
    int64_t v = 0;
    mio_r_check(mio_brep_info_int(info, section, i, field, &v), "brep_info");
    return v;
}

static SEXP brep_num(const mio_format_info *info, int32_t section, int64_t i, int32_t field) {
    return Rf_ScalarReal((double)brep_ival(info, section, i, field));
}

static SEXP brep_flag(const mio_format_info *info, int32_t section, int64_t i, int32_t field) {
    return Rf_ScalarLogical(brep_ival(info, section, i, field) != 0);
}

static SEXP brep_array(const mio_format_info *info, int32_t section, int64_t i, int32_t field,
                       double shift) {
    const void *data = NULL;
    mio_dtype dtype = MIO_FLOAT64;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_r_check(mio_brep_info_array(info, section, i, field, &data, &dtype, &ndim, shape),
                "brep_info");
    return array_to_r(data, dtype, ndim, shape, 0, shift);
}

static SEXP brep_curve(const mio_format_info *info, int32_t section, int64_t i) {
    SEXP f[9];
    f[0] = PROTECT(brep_num(info, section, i, 10));
    f[1] = PROTECT(brep_array(info, section, i, 12, 0.0));
    f[2] = PROTECT(brep_flag(info, section, i, 12));
    f[3] = PROTECT(brep_flag(info, section, i, 11));
    f[4] = PROTECT(brep_array(info, section, i, 10, 0.0));
    f[5] = PROTECT(brep_array(info, section, i, 11, 0.0));
    f[6] = PROTECT(brep_array(info, section, i, 13, 0.0));
    f[7] = PROTECT(brep_array(info, section, i, 14, 0.0));
    f[8] = PROTECT(brep_str(info, section, i, 10));
    const char *names[9] = {"degree", "knots",   "reduced_knots", "rational",    "control_points",
                            "weights", "ids",    "active_range",  "extra"};
    SEXP out = mio_r_named_list(9, names, f);
    UNPROTECT(9);
    return out;
}

static SEXP brep_trim(const mio_format_info *info, int64_t t) {
    SEXP f[4];
    f[0] = PROTECT(brep_num(info, MIO_BREP_TRIMS, t, 2));
    f[1] = PROTECT(brep_flag(info, MIO_BREP_TRIMS, t, 3));
    f[2] = PROTECT(brep_curve(info, MIO_BREP_TRIMS, t));
    f[3] = PROTECT(brep_str(info, MIO_BREP_TRIMS, t, 0));
    const char *names[4] = {"trim_index", "direction", "curve", "extra"};
    SEXP out = mio_r_named_list(4, names, f);
    UNPROTECT(4);
    return out;
}

static SEXP brep_refs(const mio_format_info *info, int64_t first, int64_t count) {
    SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)count));
    for (int64_t k = 0; k < count; ++k) {
        const int64_t r = first + k;
        SEXP f[6];
        f[0] = PROTECT(brep_num(info, MIO_BREP_REFS, r, 0));
        f[1] = PROTECT(brep_str(info, MIO_BREP_REFS, r, 0));
        f[2] = PROTECT(brep_num(info, MIO_BREP_REFS, r, 1));
        f[3] = PROTECT(brep_num(info, MIO_BREP_REFS, r, 2));
        f[4] = PROTECT(brep_array(info, MIO_BREP_REFS, r, 0, 0.0));
        f[5] = PROTECT(brep_str(info, MIO_BREP_REFS, r, 1));
        const char *names[6] = {"id",                 "name", "trim_index",
                                "relative_direction", "local_coordinates", "extra"};
        SET_VECTOR_ELT(out, k, mio_r_named_list(6, names, f));
        UNPROTECT(6);
    }
    UNPROTECT(1);
    return out;
}

static int64_t brep_count(const mio_format_info *info, int32_t section) {
    const int64_t n = mio_brep_info_count(info, section);
    if (n < 0) mio_r_fail("brep_info");
    return n;
}

SEXP R_mio_brep_info(SEXP x) {
    const mio_format_info *info = info_of(x, 0);
    SEXP v[6];
    const int64_t np = brep_count(info, MIO_BREP_PATCHES);
    v[0] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)np));
    for (int64_t p = 0; p < np; ++p) {
        SEXP f[4];
        f[0] = PROTECT(Rf_allocVector(REALSXP, 2));
        REAL(f[0])[0] = (double)brep_ival(info, MIO_BREP_PATCHES, p, 0);
        REAL(f[0])[1] = (double)brep_ival(info, MIO_BREP_PATCHES, p, 1);
        f[1] = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(f[1], 0, brep_array(info, MIO_BREP_PATCHES, p, 2, 0.0));
        SET_VECTOR_ELT(f[1], 1, brep_array(info, MIO_BREP_PATCHES, p, 3, 0.0));
        f[2] = PROTECT(brep_array(info, MIO_BREP_PATCHES, p, 0, 0.0));
        f[3] = PROTECT(brep_array(info, MIO_BREP_PATCHES, p, 1, 0.0));
        const char *names[4] = {"degrees", "knots", "control_points", "weights"};
        SET_VECTOR_ELT(v[0], p, mio_r_named_list(4, names, f));
        UNPROTECT(4);
    }
    const int64_t nb = brep_count(info, MIO_BREP_BREPS);
    v[1] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nb));
    for (int64_t b = 0; b < nb; ++b) {
        SEXP f[6];
        f[0] = PROTECT(brep_num(info, MIO_BREP_BREPS, b, 0));
        f[1] = PROTECT(brep_str(info, MIO_BREP_BREPS, b, 0));
        f[2] = PROTECT(brep_array(info, MIO_BREP_BREPS, b, 0, 1.0));
        f[3] = PROTECT(brep_array(info, MIO_BREP_BREPS, b, 1, 1.0));
        f[4] = PROTECT(brep_array(info, MIO_BREP_BREPS, b, 2, 1.0));
        f[5] = PROTECT(brep_str(info, MIO_BREP_BREPS, b, 1));
        const char *names[6] = {"id", "name", "faces", "edges", "vertices", "extra"};
        SET_VECTOR_ELT(v[1], b, mio_r_named_list(6, names, f));
        UNPROTECT(6);
    }
    const int64_t nf = brep_count(info, MIO_BREP_FACES);
    v[2] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nf));
    for (int64_t fi = 0; fi < nf; ++fi) {
        const int64_t first_loop = brep_ival(info, MIO_BREP_FACES, fi, 6);
        const int64_t loops = brep_ival(info, MIO_BREP_FACES, fi, 7);
        const int64_t first_trim = brep_ival(info, MIO_BREP_FACES, fi, 8);
        const int64_t trims = brep_ival(info, MIO_BREP_FACES, fi, 9);
        SEXP f[12];
        f[0] = PROTECT(brep_num(info, MIO_BREP_FACES, fi, 0));
        f[1] = PROTECT(brep_str(info, MIO_BREP_FACES, fi, 0));
        f[2] = PROTECT(Rf_ScalarReal((double)brep_ival(info, MIO_BREP_FACES, fi, 1) + 1.0));
        f[3] = PROTECT(brep_flag(info, MIO_BREP_FACES, fi, 2));
        f[4] = PROTECT(brep_flag(info, MIO_BREP_FACES, fi, 3));
        f[5] = PROTECT(brep_flag(info, MIO_BREP_FACES, fi, 4));
        f[6] = PROTECT(brep_flag(info, MIO_BREP_FACES, fi, 5));
        f[7] = PROTECT(brep_array(info, MIO_BREP_FACES, fi, 0, 0.0));
        f[8] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)loops));
        for (int64_t l = 0; l < loops; ++l) {
            const int64_t li = first_loop + l;
            const int64_t t0 = brep_ival(info, MIO_BREP_LOOPS, li, 1);
            const int64_t tn = brep_ival(info, MIO_BREP_LOOPS, li, 2);
            SEXP lf[3];
            lf[0] = PROTECT(brep_str(info, MIO_BREP_LOOPS, li, 0));
            lf[1] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)tn));
            for (int64_t t = 0; t < tn; ++t) SET_VECTOR_ELT(lf[1], t, brep_trim(info, t0 + t));
            lf[2] = PROTECT(brep_str(info, MIO_BREP_LOOPS, li, 1));
            const char *lnames[3] = {"type", "trims", "extra"};
            SET_VECTOR_ELT(f[8], l, mio_r_named_list(3, lnames, lf));
            UNPROTECT(3);
        }
        int64_t embedded = 0;
        for (int64_t t = 0; t < trims; ++t)
            if (brep_ival(info, MIO_BREP_TRIMS, first_trim + t, 1) < 0) ++embedded;
        f[9] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)embedded));
        for (int64_t t = 0, k = 0; t < trims; ++t)
            if (brep_ival(info, MIO_BREP_TRIMS, first_trim + t, 1) < 0)
                SET_VECTOR_ELT(f[9], k++, brep_trim(info, first_trim + t));
        f[10] = PROTECT(brep_str(info, MIO_BREP_FACES, fi, 2));
        f[11] = PROTECT(brep_str(info, MIO_BREP_FACES, fi, 1));
        const char *names[12] = {"id",           "name",          "patch",
                                 "swapped_normal", "trimmed",     "rational",
                                 "reduced_knots", "control_point_ids", "loops",
                                 "embedded_edges", "surface_extra", "extra"};
        SET_VECTOR_ELT(v[2], fi, mio_r_named_list(12, names, f));
        UNPROTECT(12);
    }
    const int64_t ne = brep_count(info, MIO_BREP_EDGES);
    v[3] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)ne));
    for (int64_t e = 0; e < ne; ++e) {
        SEXP f[5];
        f[0] = PROTECT(brep_num(info, MIO_BREP_EDGES, e, 0));
        f[1] = PROTECT(brep_str(info, MIO_BREP_EDGES, e, 0));
        f[2] = PROTECT(brep_ival(info, MIO_BREP_EDGES, e, 1) ? brep_curve(info, MIO_BREP_EDGES, e)
                                                             : R_NilValue);
        f[3] = PROTECT(brep_refs(info, brep_ival(info, MIO_BREP_EDGES, e, 2),
                                 brep_ival(info, MIO_BREP_EDGES, e, 3)));
        f[4] = PROTECT(brep_str(info, MIO_BREP_EDGES, e, 1));
        const char *names[5] = {"id", "name", "curve", "topology", "extra"};
        SET_VECTOR_ELT(v[3], e, mio_r_named_list(5, names, f));
        UNPROTECT(5);
    }
    const int64_t nv = brep_count(info, MIO_BREP_VERTICES);
    v[4] = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)nv));
    for (int64_t k = 0; k < nv; ++k) {
        SEXP f[4];
        f[0] = PROTECT(brep_num(info, MIO_BREP_VERTICES, k, 0));
        f[1] = PROTECT(brep_str(info, MIO_BREP_VERTICES, k, 0));
        f[2] = PROTECT(brep_refs(info, brep_ival(info, MIO_BREP_VERTICES, k, 1),
                                 brep_ival(info, MIO_BREP_VERTICES, k, 2)));
        f[3] = PROTECT(brep_str(info, MIO_BREP_VERTICES, k, 1));
        const char *names[4] = {"id", "name", "topology", "extra"};
        SET_VECTOR_ELT(v[4], k, mio_r_named_list(4, names, f));
        UNPROTECT(4);
    }
    v[5] = PROTECT(brep_str(info, MIO_BREP_FILE, 0, 0));
    const char *names[6] = {"patches", "breps", "faces", "edges", "vertices", "extra"};
    SEXP out = mio_r_named_list(6, names, v);
    UNPROTECT(6);
    return out;
}
