/* Lifecycle, file I/O, the setter/getter pairs, and named regions. */

#include "mio_r.h"

#include <string.h>

/* --- introspection ------------------------------------------------------ */

SEXP R_mio_version(void) { return Rf_mkString(mio_version()); }
SEXP R_mio_mesh_backend(void) { return Rf_mkString(mio_mesh_backend()); }
SEXP R_mio_last_error(void) { return Rf_mkString(mio_last_error()); }

SEXP R_mio_format_readable(SEXP format) {
    return Rf_ScalarLogical(mio_format_readable(mio_r_string(format, "format")) != 0);
}

SEXP R_mio_format_writable(SEXP format) {
    return Rf_ScalarLogical(mio_format_writable(mio_r_string(format, "format")) != 0);
}

typedef struct {
    const char *path;
} sniff_ctx;

static int64_t sniff_getter(void *ctx, char *buf, int64_t buflen) {
    return mio_sniff_format(((sniff_ctx *)ctx)->path, buf, buflen);
}

SEXP R_mio_sniff_format(SEXP path) {
    sniff_ctx ctx = {mio_r_string(path, "path")};
    return mio_r_getstring(sniff_getter, &ctx, "sniff_format");
}

SEXP R_mio_cell_type_num_nodes(SEXP name) {
    mio_cell_type t = mio_cell_type_from_name(mio_r_string(name, "name"));
    return Rf_ScalarInteger(mio_cell_type_num_nodes(t));
}

SEXP R_mio_cell_type_dimension(SEXP name) {
    mio_cell_type t = mio_cell_type_from_name(mio_r_string(name, "name"));
    return Rf_ScalarInteger(mio_cell_type_dimension(t));
}

SEXP R_mio_reader_supports_options(SEXP format) {
    int r = mio_reader_supports_options(mio_r_string(format, "format"));
    if (r < 0) mio_r_fail("reader_supports_options");
    return Rf_ScalarLogical(r != 0);
}

/* --- lifecycle and file I/O --------------------------------------------- */

SEXP R_mio_mesh_create(void) {
    mio_mesh *m = (mio_mesh *)mio_r_check_ptr(mio_mesh_create(), "mesh_create");
    return mio_r_wrap_mesh(m);
}

SEXP R_mio_mesh_release(SEXP x) {
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != mio_r_mesh_tag) {
        Rf_error("expected a mio_mesh object");
    }
    mio_mesh *m = (mio_mesh *)R_ExternalPtrAddr(x);
    if (m != NULL) {
        mio_mesh_free(m);
        R_ClearExternalPtr(x);
    }
    return R_NilValue;
}

SEXP R_mio_mesh_is_open(SEXP x) {
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != mio_r_mesh_tag) {
        return Rf_ScalarLogical(FALSE);
    }
    return Rf_ScalarLogical(R_ExternalPtrAddr(x) != NULL);
}

SEXP R_mio_read(SEXP path, SEXP format, SEXP points_only, SEXP metadata_only, SEXP arrays,
                SEXP mmap_mode, SEXP time_step, SEXP lenient) {
    const char *p = mio_r_string(path, "path");
    const char *f = mio_r_opt_string(format);

    mio_read_opts opts;
    mio_read_opts_init(&opts);
    opts.points_only = mio_r_bool(points_only, "points_only");
    opts.metadata_only = mio_r_bool(metadata_only, "metadata_only");
    opts.mmap_mode = mio_r_int(mmap_mode, "mmap_mode");
    opts.time_step = mio_r_int(time_step, "time_step");
    opts.lenient = mio_r_bool(lenient, "lenient");

    /* NULL means "every array"; a valid pointer with count 0 means "no arrays
     * at all". The distinction is load-bearing at the ABI, so an R NULL and an
     * R character(0) are deliberately different requests here. */
    SEXP shelter;
    int64_t count = 0;
    const char *const *names = mio_r_names(arrays, &count, &shelter);
    static const char *const kNoArrays[1] = {NULL};
    if (arrays != R_NilValue && names == NULL) {
        names = kNoArrays; /* character(0): a non-NULL pointer, zero entries */
    }
    opts.arrays = names;
    opts.num_arrays = count;

    mio_mesh *m = mio_read_ex(p, f, &opts);
    UNPROTECT(1); /* shelter */
    if (m == NULL) mio_r_fail("read");
    return mio_r_wrap_mesh(m);
}

SEXP R_mio_write(SEXP mesh, SEXP path, SEXP format) {
    mio_r_check(mio_write(mio_r_string(path, "path"), mio_r_mesh(mesh),
                          mio_r_opt_string(format)),
                "write");
    return R_NilValue;
}

SEXP R_mio_convert(SEXP in_path, SEXP in_format, SEXP out_path, SEXP out_format) {
    mio_r_check(mio_convert(mio_r_string(in_path, "in_path"), mio_r_opt_string(in_format),
                            mio_r_string(out_path, "out_path"),
                            mio_r_opt_string(out_format)),
                "convert");
    return R_NilValue;
}

/* --- file metadata ------------------------------------------------------ */

typedef struct {
    const mio_read_metadata *meta;
    int64_t index;
    int location;
} meta_ctx;

static int64_t meta_block_type_getter(void *c, char *buf, int64_t buflen) {
    meta_ctx *x = (meta_ctx *)c;
    return mio_read_metadata_cell_block_type(x->meta, x->index, buf, buflen);
}

static int64_t meta_name_getter(void *c, char *buf, int64_t buflen) {
    meta_ctx *x = (meta_ctx *)c;
    return mio_read_metadata_name(x->meta, x->location, x->index, buf, buflen);
}

static int64_t meta_region_name_getter(void *c, char *buf, int64_t buflen) {
    meta_ctx *x = (meta_ctx *)c;
    return mio_read_metadata_region_name(x->meta, x->index, buf, buflen);
}

/* The mesh.c-local region_kind_name() (defined further down, near R_mio_regions)
 * is reused here rather than duplicated -- forward-declared since this function
 * is defined earlier in the file. */
static const char *region_kind_name(int32_t kind);

/* One named region's shape, without its entries -- the read_metadata
 * counterpart of R_mio_regions(). */
static SEXP meta_regions(const mio_read_metadata *meta) {
    int64_t n = mio_read_metadata_num_regions(meta);
    if (n < 0) n = 0;
    SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)n));
    for (int64_t i = 0; i < n; ++i) {
        meta_ctx ctx = {meta, i, 0};
        SEXP name = PROTECT(mio_r_getstring(meta_region_name_getter, &ctx, "region name"));
        mio_region_info info;
        if (mio_read_metadata_region_info(meta, i, &info) != MIO_OK) {
            UNPROTECT(1); /* name */
            continue;     /* skip a region the core could not describe */
        }
        SEXP kind = PROTECT(Rf_mkString(region_kind_name(info.kind)));
        SEXP dim = PROTECT(Rf_ScalarInteger(info.dim));
        SEXP tag = PROTECT(Rf_ScalarReal((double)info.tag));
        SEXP num_entries = PROTECT(Rf_ScalarReal((double)info.num_entries));
        const char *item_names[] = {"name", "kind", "dim", "tag", "num_entries"};
        SEXP item_values[] = {name, kind, dim, tag, num_entries};
        SEXP item = PROTECT(mio_r_named_list(5, item_names, item_values));
        SET_VECTOR_ELT(out, (R_xlen_t)i, item);
        UNPROTECT(6);
    }
    UNPROTECT(1);
    return out;
}

static SEXP meta_names(const mio_read_metadata *meta, int location) {
    int64_t n = mio_read_metadata_num_names(meta, location);
    if (n < 0) n = 0;
    SEXP out = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)n));
    for (int64_t i = 0; i < n; ++i) {
        meta_ctx ctx = {meta, i, location};
        SEXP s = PROTECT(mio_r_getstring(meta_name_getter, &ctx, "metadata name"));
        SET_STRING_ELT(out, (R_xlen_t)i, STRING_ELT(s, 0));
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return out;
}

SEXP R_mio_read_metadata(SEXP path, SEXP format) {
    mio_read_metadata *meta = mio_read_metadata_create(mio_r_string(path, "path"),
                                                       mio_r_opt_string(format));
    if (meta == NULL) mio_r_fail("read_metadata");

    SEXP out = R_NilValue;
    int64_t nblocks = mio_read_metadata_num_cell_blocks(meta);
    if (nblocks < 0) nblocks = 0;

    SEXP types = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)nblocks));
    SEXP counts = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)nblocks));
    SEXP npcs = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)nblocks));
    SEXP ragged = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t)nblocks));
    for (int64_t i = 0; i < nblocks; ++i) {
        int64_t nc = 0, npc = 0;
        int is_ragged = 0;
        mio_status st = mio_read_metadata_cell_block(meta, i, &nc, &npc, &is_ragged);
        if (st != MIO_OK) {
            mio_read_metadata_free(meta);
            UNPROTECT(4);
            mio_r_fail("read_metadata cell block");
        }
        REAL(counts)[i] = (double)nc;
        REAL(npcs)[i] = (double)npc;
        LOGICAL(ragged)[i] = is_ragged != 0;
        meta_ctx ctx = {meta, i, 0};
        SEXP s = PROTECT(mio_r_getstring(meta_block_type_getter, &ctx, "block type"));
        SET_STRING_ELT(types, (R_xlen_t)i, STRING_ELT(s, 0));
        UNPROTECT(1);
    }

    double lo[3], hi[3];
    int have_bbox = mio_read_metadata_bbox(meta, lo, hi) == MIO_OK;
    SEXP bbox = R_NilValue;
    if (have_bbox) {
        bbox = PROTECT(Rf_allocMatrix(REALSXP, 3, 2));
        memcpy(REAL(bbox), lo, 3 * sizeof(double));
        memcpy(REAL(bbox) + 3, hi, 3 * sizeof(double));
    } else {
        bbox = PROTECT(R_NilValue);
    }

    SEXP pn = PROTECT(meta_names(meta, MIO_DATA_POINT));
    SEXP cn = PROTECT(meta_names(meta, MIO_DATA_CELL));
    SEXP fn = PROTECT(meta_names(meta, MIO_DATA_FIELD));

    SEXP npoints = PROTECT(Rf_ScalarReal((double)mio_read_metadata_num_points(meta)));
    SEXP pdim = PROTECT(Rf_ScalarReal((double)mio_read_metadata_point_dim(meta)));
    SEXP ncells = PROTECT(Rf_ScalarReal((double)mio_read_metadata_num_cells(meta)));
    SEXP nblk = PROTECT(Rf_ScalarReal((double)nblocks));
    SEXP fell = PROTECT(Rf_ScalarLogical(mio_read_metadata_fell_back(meta) == 1));

    /* The file's recorded time-series values; length 0 for a format with no
     * time concept. This is the count `time_step` may name. */
    int64_t nsteps = mio_read_metadata_num_time_values(meta);
    if (nsteps < 0) nsteps = 0;
    SEXP times = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)nsteps));
    if (nsteps > 0) mio_read_metadata_time_values(meta, REAL(times), nsteps);

    /* The file's named regions, without their entries. Empty on a native
     * metadata path (none of those formats currently map regions); cheap
     * whenever the summary came from an already-read mesh. */
    SEXP regions = PROTECT(meta_regions(meta));

    mio_read_metadata_free(meta);

    const char *names[] = {"num_points",      "point_dim",       "num_cells",
                           "num_cell_blocks", "cell_block_types", "cell_block_num_cells",
                           "cell_block_nodes_per_cell", "cell_block_is_ragged",
                           "point_data_names", "cell_data_names", "field_data_names",
                           "bbox",            "fell_back",       "time_values", "regions"};
    SEXP values[] = {npoints, pdim, ncells, nblk, types, counts, npcs,
                     ragged,  pn,   cn,     fn,   bbox,  fell,   times, regions};
    out = PROTECT(mio_r_named_list(15, names, values));
    UNPROTECT(16); /* the 15 values + out */
    return out;
}

/* --- counts ------------------------------------------------------------- */

SEXP R_mio_num_points(SEXP mesh) {
    return Rf_ScalarReal((double)mio_mesh_num_points(mio_r_mesh(mesh)));
}

SEXP R_mio_point_dim(SEXP mesh) {
    return Rf_ScalarReal((double)mio_mesh_point_dim(mio_r_mesh(mesh)));
}

SEXP R_mio_num_cell_blocks(SEXP mesh) {
    return Rf_ScalarReal((double)mio_mesh_num_cell_blocks(mio_r_mesh(mesh)));
}

SEXP R_mio_cell_block_info(SEXP mesh, SEXP block) {
    int64_t nc = 0, npc = 0;
    int32_t ragged = 0;
    mio_r_check(mio_mesh_cell_block_info(mio_r_mesh(mesh), mio_r_int64(block, "block"), &nc,
                                         &npc, &ragged),
                "cell_block_info");
    SEXP a = PROTECT(Rf_ScalarReal((double)nc));
    SEXP b = PROTECT(Rf_ScalarReal((double)npc));
    SEXP c = PROTECT(Rf_ScalarLogical(ragged != 0));
    const char *names[] = {"num_cells", "nodes_per_cell", "is_ragged"};
    SEXP values[] = {a, b, c};
    SEXP out = PROTECT(mio_r_named_list(3, names, values));
    UNPROTECT(4);
    return out;
}

typedef struct {
    const mio_mesh *mesh;
    int64_t index;
} mesh_ctx;

static int64_t block_type_getter(void *c, char *buf, int64_t buflen) {
    mesh_ctx *x = (mesh_ctx *)c;
    return mio_mesh_cell_block_type(x->mesh, x->index, buf, buflen);
}

SEXP R_mio_cell_block_type(SEXP mesh, SEXP block) {
    mesh_ctx ctx = {mio_r_mesh(mesh), mio_r_int64(block, "block")};
    return mio_r_getstring(block_type_getter, &ctx, "cell_block_type");
}

/* --- points ------------------------------------------------------------- */

SEXP R_mio_points(SEXP mesh) {
    const mio_mesh *m = mio_r_mesh(mesh);
    const void *data = NULL;
    mio_dtype dt;
    mio_r_check(mio_mesh_get_points(m, &data, &dt), "get_points");
    int64_t n = mio_mesh_num_points(m);
    int64_t dim = mio_mesh_point_dim(m);

    /* Rf_allocMatrix(REALSXP, dim, n) is column-major, i.e. byte-identical to
     * the C API's row-major (n, dim). One memcpy, no transpose -- the same
     * identity the Julia and Fortran bindings rely on. */
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int)dim, (int)n));
    SEXP flat = PROTECT(mio_r_copy_as_real(data, dt, (R_xlen_t)(n * dim)));
    memcpy(REAL(out), REAL(flat), (size_t)(n * dim) * sizeof(double));
    UNPROTECT(2);
    return out;
}

/* Shared by the 1-based and 0-based connectivity accessors. */
static SEXP conn_matrix(SEXP mesh, SEXP block, int shift) {
    const mio_mesh *m = mio_r_mesh(mesh);
    int64_t b = mio_r_int64(block, "block");
    int64_t nc = 0, npc = 0;
    int32_t ragged = 0;
    mio_r_check(mio_mesh_cell_block_info(m, b, &nc, &npc, &ragged), "cell_block_info");
    if (ragged) {
        Rf_error("cell block %d is ragged (polygons or polyhedra of varying size); "
                 "its connectivity is not reachable through the meshio++ C API",
                 (int)(b + 1));
    }
    const void *conn = NULL;
    mio_dtype dt;
    mio_r_check(mio_mesh_cell_block_conn(m, b, &conn, &dt), "cell_block_conn");

    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int)npc, (int)nc));
    SEXP flat = PROTECT(mio_r_copy_as_real(conn, dt, (R_xlen_t)(nc * npc)));
    double *dst = REAL(out);
    const double *src = REAL(flat);
    for (R_xlen_t i = 0; i < (R_xlen_t)(nc * npc); ++i) {
        dst[i] = src[i] + (shift ? 1.0 : 0.0);
    }
    UNPROTECT(2);
    return out;
}

SEXP R_mio_connectivity(SEXP mesh, SEXP block) { return conn_matrix(mesh, block, 1); }
SEXP R_mio_connectivity_raw(SEXP mesh, SEXP block) { return conn_matrix(mesh, block, 0); }

/* --- data arrays -------------------------------------------------------- */

static int64_t point_data_name_getter(void *c, char *buf, int64_t buflen) {
    mesh_ctx *x = (mesh_ctx *)c;
    return mio_mesh_point_data_name(x->mesh, x->index, buf, buflen);
}
static int64_t cell_data_name_getter(void *c, char *buf, int64_t buflen) {
    mesh_ctx *x = (mesh_ctx *)c;
    return mio_mesh_cell_data_name(x->mesh, x->index, buf, buflen);
}
static int64_t field_data_name_getter(void *c, char *buf, int64_t buflen) {
    mesh_ctx *x = (mesh_ctx *)c;
    return mio_mesh_field_data_name(x->mesh, x->index, buf, buflen);
}

static SEXP data_names(const mio_mesh *m, int64_t n, mio_r_str_getter getter) {
    if (n < 0) n = 0;
    SEXP out = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)n));
    for (int64_t i = 0; i < n; ++i) {
        mesh_ctx ctx = {m, i};
        SEXP s = PROTECT(mio_r_getstring(getter, &ctx, "data name"));
        SET_STRING_ELT(out, (R_xlen_t)i, STRING_ELT(s, 0));
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return out;
}

SEXP R_mio_point_data_names(SEXP mesh) {
    const mio_mesh *m = mio_r_mesh(mesh);
    return data_names(m, mio_mesh_num_point_data(m), point_data_name_getter);
}
SEXP R_mio_cell_data_names(SEXP mesh) {
    const mio_mesh *m = mio_r_mesh(mesh);
    return data_names(m, mio_mesh_num_cell_data(m), cell_data_name_getter);
}
SEXP R_mio_field_data_names(SEXP mesh) {
    const mio_mesh *m = mio_r_mesh(mesh);
    return data_names(m, mio_mesh_num_field_data(m), field_data_name_getter);
}

SEXP R_mio_cell_data_num_blocks(SEXP mesh, SEXP name) {
    int64_t n = mio_mesh_cell_data_num_blocks(mio_r_mesh(mesh), mio_r_string(name, "name"));
    if (n < 0) mio_r_fail("cell_data_num_blocks");
    return Rf_ScalarReal((double)n);
}

/* Shape a returned data array. The C API reports it row-major as
 * (rows, components...); the column-major R array over the same memory has
 * the reversed dim vector, so a scalar field is a plain vector and an
 * (n, 3) vector field is a 3 x n matrix. */
static SEXP shape_data(const void *data, mio_dtype dt, int32_t ndim, const int64_t *shape) {
    R_xlen_t total = 1;
    for (int32_t i = 0; i < ndim; ++i) total *= (R_xlen_t)shape[i];
    SEXP out = PROTECT(mio_r_copy_as_real(data, dt, total));
    if (ndim > 1) {
        SEXP dim = PROTECT(Rf_allocVector(INTSXP, ndim));
        for (int32_t i = 0; i < ndim; ++i) {
            INTEGER(dim)[i] = (int)shape[ndim - 1 - i]; /* reversed */
        }
        Rf_setAttrib(out, R_DimSymbol, dim);
        UNPROTECT(1);
    }
    SEXP dtn = PROTECT(Rf_mkString(mio_r_dtype_name(dt)));
    Rf_setAttrib(out, Rf_install("dtype"), dtn);
    UNPROTECT(2);
    return out;
}

SEXP R_mio_point_data(SEXP mesh, SEXP name) {
    const void *data = NULL;
    mio_dtype dt;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_r_check(mio_mesh_get_point_data(mio_r_mesh(mesh), mio_r_string(name, "name"), &data,
                                        &dt, &ndim, shape),
                "get_point_data");
    return shape_data(data, dt, ndim, shape);
}

SEXP R_mio_cell_data(SEXP mesh, SEXP name, SEXP block) {
    const void *data = NULL;
    mio_dtype dt;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_r_check(mio_mesh_get_cell_data(mio_r_mesh(mesh), mio_r_string(name, "name"),
                                       mio_r_int64(block, "block"), &data, &dt, &ndim, shape),
                "get_cell_data");
    return shape_data(data, dt, ndim, shape);
}

SEXP R_mio_field_data(SEXP mesh, SEXP name) {
    const void *data = NULL;
    mio_dtype dt;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_r_check(mio_mesh_get_field_data(mio_r_mesh(mesh), mio_r_string(name, "name"), &data,
                                        &dt, &ndim, shape),
                "get_field_data");
    return shape_data(data, dt, ndim, shape);
}

/* --- setters (all COPY caller memory) ----------------------------------- */

SEXP R_mio_set_points(SEXP mesh, SEXP points) {
    if (!Rf_isMatrix(points)) Rf_error("`points` must be a (dim x num_points) matrix");
    SEXP p = PROTECT(Rf_coerceVector(points, REALSXP));
    SEXP dim = Rf_getAttrib(points, R_DimSymbol);
    int64_t d = INTEGER(dim)[0];
    int64_t n = INTEGER(dim)[1];
    /* Column-major R == row-major C, so the buffer goes across as it stands. */
    mio_r_check(mio_mesh_set_points(mio_r_mesh(mesh), MIO_FLOAT64, n, d, REAL(p)),
                "set_points");
    UNPROTECT(1);
    return R_NilValue;
}

SEXP R_mio_add_cell_block(SEXP mesh, SEXP cell_type, SEXP conn) {
    if (!Rf_isMatrix(conn)) {
        Rf_error("`conn` must be a (nodes_per_cell x num_cells) matrix");
    }
    SEXP c = PROTECT(Rf_coerceVector(conn, REALSXP));
    SEXP dim = Rf_getAttrib(conn, R_DimSymbol);
    int64_t npc = INTEGER(dim)[0];
    int64_t nc = INTEGER(dim)[1];
    R_xlen_t total = (R_xlen_t)(npc * nc);

    /* 1-based in R, 0-based at the ABI: the shift happens here, in the copy. */
    SEXP buf = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)(total * sizeof(int64_t))));
    int64_t *dst = (int64_t *)RAW(buf);
    const double *src = REAL(c);
    for (R_xlen_t i = 0; i < total; ++i) {
        if (ISNA(src[i]) || src[i] < 1.0) {
            UNPROTECT(2);
            Rf_error("connectivity is 1-based here; got %g", src[i]);
        }
        dst[i] = (int64_t)src[i] - 1;
    }
    mio_r_check(mio_mesh_add_cell_block(mio_r_mesh(mesh),
                                        mio_r_string(cell_type, "cell_type"), nc, npc,
                                        MIO_INT64, dst),
                "add_cell_block");
    UNPROTECT(2);
    return R_NilValue;
}

/* Shared by the three data setters. An R array shaped (components..., rows) is
 * the C row-major (rows, components...), so the C shape is the reversed dim
 * vector. */
static SEXP add_data(SEXP mesh, SEXP name, SEXP data, int which) {
    SEXP d = PROTECT(Rf_coerceVector(data, REALSXP));
    SEXP dim = Rf_getAttrib(data, R_DimSymbol);
    int32_t ndim;
    int64_t shape[MIO_MAX_NDIM] = {0};
    if (dim == R_NilValue) {
        ndim = 1;
        shape[0] = (int64_t)Rf_xlength(d);
    } else {
        ndim = (int32_t)Rf_length(dim);
        if (ndim > MIO_MAX_NDIM) {
            UNPROTECT(1);
            Rf_error("arrays of rank > %d are not supported", MIO_MAX_NDIM);
        }
        for (int32_t i = 0; i < ndim; ++i) {
            shape[i] = (int64_t)INTEGER(dim)[ndim - 1 - i]; /* reversed */
        }
    }
    const char *nm = mio_r_string(name, "name");
    mio_mesh *m = mio_r_mesh(mesh);
    mio_status st;
    if (which == 0) {
        st = mio_mesh_add_point_data(m, nm, MIO_FLOAT64, ndim, shape, REAL(d));
    } else if (which == 1) {
        st = mio_mesh_append_cell_data(m, nm, MIO_FLOAT64, ndim, shape, REAL(d));
    } else {
        st = mio_mesh_add_field_data(m, nm, MIO_FLOAT64, ndim, shape, REAL(d));
    }
    UNPROTECT(1);
    mio_r_check(st, "add data");
    return R_NilValue;
}

SEXP R_mio_add_point_data(SEXP mesh, SEXP name, SEXP data) {
    return add_data(mesh, name, data, 0);
}
SEXP R_mio_append_cell_data(SEXP mesh, SEXP name, SEXP data) {
    return add_data(mesh, name, data, 1);
}
SEXP R_mio_add_field_data(SEXP mesh, SEXP name, SEXP data) {
    return add_data(mesh, name, data, 2);
}

/* --- named regions ------------------------------------------------------ */

typedef struct {
    const mio_regions *regions;
    int64_t index;
} regions_ctx;

static int64_t region_name_getter(void *c, char *buf, int64_t buflen) {
    regions_ctx *x = (regions_ctx *)c;
    return mio_regions_name(x->regions, x->index, buf, buflen);
}

static const char *region_kind_name(int32_t kind) {
    switch (kind) {
    case MIO_REGION_POINT: return "point";
    case MIO_REGION_CELL: return "cell";
    case MIO_REGION_SIDE: return "side";
    default: return "unknown";
    }
}

SEXP R_mio_regions(SEXP mesh) {
    mio_regions *regions =
        (mio_regions *)mio_r_check_ptr(mio_regions_create(mio_r_mesh(mesh)), "regions");
    int64_t n = mio_regions_count(regions);
    if (n < 0) {
        mio_regions_free(regions);
        mio_r_fail("regions_count");
    }
    SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)n));
    for (int64_t i = 0; i < n; ++i) {
        regions_ctx ctx = {regions, i};
        SEXP name = PROTECT(mio_r_getstring(region_name_getter, &ctx, "region name"));
        mio_region_info info;
        if (mio_regions_info(regions, i, &info) != MIO_OK) {
            mio_regions_free(regions);
            UNPROTECT(1); /* name */
            mio_r_fail("regions_info");
        }
        SEXP kind = PROTECT(Rf_mkString(region_kind_name(info.kind)));
        SEXP dim = PROTECT(Rf_ScalarInteger(info.dim));
        SEXP tag = PROTECT(Rf_ScalarReal((double)info.tag));

        SEXP entries = PROTECT(Rf_allocMatrix(REALSXP, (int)info.stride,
                                              (int)info.num_entries));
        if (info.num_entries > 0) {
            int64_t count = 0;
            const int64_t *src = mio_regions_entries(regions, i, &count);
            if (src == NULL) {
                mio_regions_free(regions);
                UNPROTECT(5); /* name, kind, dim, tag, entries */
                mio_r_fail("regions_entries");
            }
            double *dst = REAL(entries);
            /* +1 on the index rows. For a SIDE region row 2 is a FACET ordinal
             * within the cell type, not a mesh index, so it is NOT shifted --
             * the same rule the Fortran and Julia bindings follow. */
            for (int64_t e = 0; e < info.num_entries; ++e) {
                for (int64_t r = 0; r < info.stride; ++r) {
                    int64_t v = src[e * info.stride + r];
                    int is_facet = (info.kind == MIO_REGION_SIDE && r == 1);
                    dst[e * info.stride + r] = is_facet ? (double)v : (double)(v + 1);
                }
            }
        }
        const char *names[] = {"name", "kind", "dim", "tag", "entries"};
        SEXP values[] = {name, kind, dim, tag, entries};
        SEXP item = PROTECT(mio_r_named_list(5, names, values));
        SET_VECTOR_ELT(out, (R_xlen_t)i, item);
        UNPROTECT(6);
    }
    mio_regions_free(regions);
    UNPROTECT(1);
    return out;
}

SEXP R_mio_add_region(SEXP mesh, SEXP name, SEXP kind, SEXP entries, SEXP dim, SEXP tag) {
    const char *k = mio_r_string(kind, "kind");
    mio_region_kind kv;
    if (strcmp(k, "point") == 0) {
        kv = MIO_REGION_POINT;
    } else if (strcmp(k, "cell") == 0) {
        kv = MIO_REGION_CELL;
    } else if (strcmp(k, "side") == 0) {
        kv = MIO_REGION_SIDE;
    } else {
        Rf_error("`kind` must be \"point\", \"cell\" or \"side\", got \"%s\"", k);
    }

    SEXP e = PROTECT(Rf_coerceVector(entries, REALSXP));
    SEXP edim = Rf_getAttrib(entries, R_DimSymbol);
    int64_t stride = 1, nent = Rf_xlength(e);
    if (edim != R_NilValue) {
        if (Rf_length(edim) != 2) {
            UNPROTECT(1);
            Rf_error("`entries` must be a vector or a (stride x n) matrix");
        }
        stride = INTEGER(edim)[0];
        nent = INTEGER(edim)[1];
    }
    int64_t expected = (kv == MIO_REGION_SIDE) ? 2 : 1;
    if (stride != expected) {
        UNPROTECT(1);
        Rf_error("a \"%s\" region needs %d value(s) per entry, got %d", k, (int)expected,
                 (int)stride);
    }

    R_xlen_t total = (R_xlen_t)(stride * nent);
    SEXP buf = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)(total * sizeof(int64_t))));
    int64_t *dst = (int64_t *)RAW(buf);
    const double *src = REAL(e);
    for (int64_t i = 0; i < nent; ++i) {
        for (int64_t r = 0; r < stride; ++r) {
            R_xlen_t k2 = (R_xlen_t)(i * stride + r);
            int is_facet = (kv == MIO_REGION_SIDE && r == 1);
            if (is_facet) {
                dst[k2] = (int64_t)src[k2];
            } else {
                if (ISNA(src[k2]) || src[k2] < 1.0) {
                    UNPROTECT(2);
                    Rf_error("region entries are 1-based here; got %g", src[k2]);
                }
                dst[k2] = (int64_t)src[k2] - 1;
            }
        }
    }
    mio_status st = mio_mesh_add_region(mio_r_mesh(mesh), mio_r_string(name, "name"), kv,
                                        (int32_t)mio_r_int(dim, "dim"),
                                        mio_r_int64(tag, "tag"), dst, total);
    UNPROTECT(2);
    mio_r_check(st, "add_region");
    return R_NilValue;
}
