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
