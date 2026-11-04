/* Lifecycle, file I/O, the setter/getter pairs, and named regions.
 * A port of the R binding's src/mesh.c; see mio_oct.h for the layout rules. */

#include "mio_oct.h"

#include <math.h>
#include <string.h>

/* --- introspection ------------------------------------------------------- */

static mxArray *O_mio_version(const mxArray **a) {
    (void)a;
    return mio_o_str(mio_version());
}

static mxArray *O_mio_mesh_backend(const mxArray **a) {
    (void)a;
    return mio_o_str(mio_mesh_backend());
}

static mxArray *O_mio_last_error(const mxArray **a) {
    (void)a;
    return mio_o_str(mio_last_error());
}

static mxArray *O_mio_load_plugins(const mxArray **a) {
    return mio_o_scalar((double)mio_load_plugins(mio_o_string(a[0], "paths")));
}

static mxArray *O_mio_format_readable(const mxArray **a) {
    return mio_o_logical(mio_format_readable(mio_o_string(a[0], "format")) != 0);
}

static mxArray *O_mio_format_writable(const mxArray **a) {
    return mio_o_logical(mio_format_writable(mio_o_string(a[0], "format")) != 0);
}

typedef struct {
    const char *path;
} sniff_ctx;

static int64_t sniff_getter(void *ctx, char *buf, int64_t buflen) {
    return mio_sniff_format(((sniff_ctx *)ctx)->path, buf, buflen);
}

static mxArray *O_mio_sniff_format(const mxArray **a) {
    sniff_ctx ctx = {mio_o_string(a[0], "path")};
    return mio_o_getstring(sniff_getter, &ctx, "sniff_format");
}

static mxArray *O_mio_cell_type_num_nodes(const mxArray **a) {
    mio_cell_type t = mio_cell_type_from_name(mio_o_string(a[0], "name"));
    return mio_o_scalar((double)mio_cell_type_num_nodes(t));
}

static mxArray *O_mio_cell_type_dimension(const mxArray **a) {
    mio_cell_type t = mio_cell_type_from_name(mio_o_string(a[0], "name"));
    return mio_o_scalar((double)mio_cell_type_dimension(t));
}

static mxArray *O_mio_reader_supports_options(const mxArray **a) {
    int r = mio_reader_supports_options(mio_o_string(a[0], "format"));
    if (r < 0) mio_o_fail("reader_supports_options");
    return mio_o_logical(r != 0);
}

/* --- lifecycle and file I/O ---------------------------------------------- */

static mxArray *O_mio_mesh_create(const mxArray **a) {
    (void)a;
    return mio_o_wrap_mesh((mio_mesh *)mio_o_check_ptr(mio_mesh_create(), "mesh_create"));
}

static mxArray *O_mio_mesh_release(const mxArray **a) {
    mio_o_release(a[0], MIO_O_MESH);
    return NULL;
}

static mxArray *O_mio_mesh_is_open(const mxArray **a) {
    return mio_o_logical(mio_o_is_open(a[0], MIO_O_MESH));
}

/* read(path, format, points_only, metadata_only, arrays, mmap_mode, time_step,
 *      lenient, piece, drop_ghosts) */
static mxArray *O_mio_read(const mxArray **a) {
    const char *p = mio_o_string(a[0], "path");
    const char *f = mio_o_opt_string(a[1]);

    mio_read_opts opts;
    mio_read_opts_init(&opts);
    opts.points_only = mio_o_bool(a[2], "points_only");
    opts.metadata_only = mio_o_bool(a[3], "metadata_only");
    opts.mmap_mode = mio_o_int(a[5], "mmap_mode");
    opts.time_step = mio_o_int(a[6], "time_step");
    opts.lenient = mio_o_bool(a[7], "lenient");
    /* [] merges every piece of a partitioned file; an integer keeps that piece.
     * piece_set is what selects: piece alone (0) must never mean "merge". */
    if (!mio_o_is_null(a[8])) {
        opts.piece = mio_o_int(a[8], "piece");
        opts.piece_set = 1;
    }
    opts.drop_ghosts = mio_o_bool(a[9], "drop_ghosts");

    /* [] means "every array"; {} means "no arrays at all". The distinction is
     * load-bearing at the ABI, exactly like R's NULL vs character(0). */
    int64_t count = 0;
    const char *const *names = mio_o_names(a[4], &count, "arrays");
    static const char *const kNoArrays[1] = {NULL};
    if (!mio_o_is_null(a[4]) && names == NULL) names = kNoArrays;
    opts.arrays = names;
    opts.num_arrays = count;

    mio_mesh *m = mio_read_ex(p, f, &opts);
    if (m == NULL) mio_o_fail("read");
    return mio_o_wrap_mesh(m);
}

/* write(mesh, path, format, encoding, codec, float_format) */
static mxArray *O_mio_write(const mxArray **a) {
    mio_write_opts opts;
    mio_write_opts_init(&opts);
    opts.encoding = mio_o_int(a[3], "encoding");
    opts.codec = mio_o_int(a[4], "codec");
    opts.float_format = mio_o_opt_string(a[5]);
    mio_o_check(mio_write_ex(mio_o_string(a[1], "path"), mio_o_mesh(a[0]),
                             mio_o_opt_string(a[2]), &opts),
                "write");
    return NULL;
}

static mxArray *O_mio_convert(const mxArray **a) {
    mio_o_check(mio_convert(mio_o_string(a[0], "in_path"), mio_o_opt_string(a[1]),
                            mio_o_string(a[2], "out_path"), mio_o_opt_string(a[3])),
                "convert");
    return NULL;
}

/* write_gltf(mesh, path, options) -- `options` is a 1x19 cell in the R
 * binding's list order. */
static mxArray *O_mio_write_gltf(const mxArray **a) {
    const mxArray *o = a[2];
    if (!mxIsCell(o) || mxGetNumberOfElements(o) != 19) {
        mio_o_error("write_gltf", "invalid glTF option list");
    }
#define OPT(i) mxGetCell(o, (mwIndex)(i))
    mio_gltf_opts opts;
    mio_gltf_opts_init(&opts);
    opts.container = mio_o_int(OPT(0), "container");
    opts.up_axis = mio_o_int(OPT(1), "up_axis");
    opts.normal_weight = mio_o_int(OPT(2), "normal_weight");
    opts.normals = mio_o_bool(OPT(3), "normals");
    opts.fields = mio_o_bool(OPT(4), "fields");
    opts.recenter = mio_o_bool(OPT(5), "recenter");
    opts.by_region = mio_o_bool(OPT(6), "by_region");
    opts.unlit = mio_o_bool(OPT(7), "unlit");
    opts.component = mio_o_int(OPT(8), "component") - 1;
    opts.component_set = mio_o_bool(OPT(9), "component_set");
    opts.vmin_set = mio_o_bool(OPT(10), "vmin_set");
    opts.vmax_set = mio_o_bool(OPT(11), "vmax_set");
    opts.split_angle = mio_o_double(OPT(12), "split_angle");
    opts.scale = mio_o_double(OPT(13), "scale");
    opts.vmin = mio_o_double(OPT(14), "vmin");
    opts.vmax = mio_o_double(OPT(15), "vmax");
    opts.color_by = mio_o_opt_string(OPT(16));
    opts.cmap = mio_o_opt_string(OPT(17));
    opts.nan_color = mio_o_opt_string(OPT(18));
#undef OPT
    mio_o_check(mio_write_gltf(mio_o_string(a[1], "path"), mio_o_mesh(a[0]), &opts),
                "write_gltf");
    return NULL;
}

static mxArray *O_mio_med_mesh_names(const mxArray **a) {
    const char *p = mio_o_string(a[0], "path");
    int64_t n = mio_med_mesh_count(p);
    if (n < 0) mio_o_fail("med_mesh_names");
    mxArray *out = mio_o_cellstr(n);
    for (int64_t i = 0; i < n; ++i) {
        int64_t len = mio_med_mesh_name(p, i, NULL, 0);
        if (len < 0) mio_o_fail("med_mesh_name");
        char *buf = (char *)mxCalloc((size_t)len + 1, 1);
        if (mio_med_mesh_name(p, i, buf, len + 1) < 0) mio_o_fail("med_mesh_name");
        mxSetCell(out, (mwIndex)i, mxCreateString(buf));
        mxFree(buf);
    }
    return out;
}

static mxArray *O_mio_med_read_named(const mxArray **a) {
    mio_read_opts opts;
    mio_read_opts_init(&opts);
    opts.time_step = mio_o_int64(a[2], "time_step");
    opts.lenient = mio_o_bool(a[3], "lenient");
    mio_mesh *mesh =
        mio_med_read_named(mio_o_string(a[0], "path"), mio_o_string(a[1], "name"), &opts);
    if (mesh == NULL) mio_o_fail("med_read_named");
    return mio_o_wrap_mesh(mesh);
}

static mxArray *O_mio_med_write_multi(const mxArray **a) {
    const mxArray *meshes = a[1];
    int64_t count = 0;
    const char *const *names = mio_o_names(a[2], &count, "names");
    if (!mxIsCell(meshes) || (int64_t)mxGetNumberOfElements(meshes) != count || count == 0) {
        mio_o_error("med_write_multi", "MED: provide a nonempty mesh list and one name per mesh");
    }
    const mio_mesh **handles = (const mio_mesh **)mxMalloc((size_t)count * sizeof(mio_mesh *));
    for (int64_t i = 0; i < count; ++i) handles[i] = mio_o_mesh(mxGetCell(meshes, (mwIndex)i));
    mio_o_check(mio_med_write_multi(mio_o_string(a[0], "path"), handles, names, count,
                                    mio_o_string(a[3], "version")),
                "med_write_multi");
    return NULL;
}

/* --- file metadata ------------------------------------------------------- */

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

static const char *region_kind_name(int32_t kind) {
    switch (kind) {
    case MIO_REGION_POINT: return "point";
    case MIO_REGION_CELL: return "cell";
    case MIO_REGION_SIDE: return "side";
    default: return "unknown";
    }
}

static const char *kRegionMetaFields[] = {"name", "kind", "dim", "tag", "num_entries"};

/* The file's named regions without their entries, as an n x 1 struct array. */
static mxArray *meta_regions(const mio_read_metadata *meta) {
    int64_t n = mio_read_metadata_num_regions(meta);
    if (n < 0) n = 0;
    /* Count the describable ones first: a region the core cannot describe is
     * skipped, as in the R binding. */
    mio_region_info *infos = (mio_region_info *)mxCalloc((size_t)(n > 0 ? n : 1),
                                                         sizeof(mio_region_info));
    int64_t *keep = (int64_t *)mxCalloc((size_t)(n > 0 ? n : 1), sizeof(int64_t));
    int64_t kept = 0;
    for (int64_t i = 0; i < n; ++i) {
        if (mio_read_metadata_region_info(meta, i, &infos[kept]) == MIO_OK) keep[kept++] = i;
    }
    mxArray *out = mio_o_struct_array(kept, 5, kRegionMetaFields);
    for (int64_t k = 0; k < kept; ++k) {
        meta_ctx ctx = {meta, keep[k], 0};
        mxSetField(out, (mwIndex)k, "name",
                   mio_o_getstring(meta_region_name_getter, &ctx, "region name"));
        mxSetField(out, (mwIndex)k, "kind", mio_o_str(region_kind_name(infos[k].kind)));
        mxSetField(out, (mwIndex)k, "dim", mio_o_scalar((double)infos[k].dim));
        mxSetField(out, (mwIndex)k, "tag", mio_o_scalar((double)infos[k].tag));
        mxSetField(out, (mwIndex)k, "num_entries", mio_o_scalar((double)infos[k].num_entries));
    }
    mxFree(infos);
    mxFree(keep);
    return out;
}

static mxArray *meta_names(const mio_read_metadata *meta, int location) {
    int64_t n = mio_read_metadata_num_names(meta, location);
    if (n < 0) n = 0;
    mxArray *out = mio_o_cellstr(n);
    for (int64_t i = 0; i < n; ++i) {
        meta_ctx ctx = {meta, i, location};
        mxSetCell(out, (mwIndex)i, mio_o_getstring(meta_name_getter, &ctx, "metadata name"));
    }
    return out;
}

static mxArray *O_mio_read_metadata(const mxArray **a) {
    mio_read_metadata *meta =
        mio_read_metadata_create(mio_o_string(a[0], "path"), mio_o_opt_string(a[1]));
    if (meta == NULL) mio_o_fail("read_metadata");

    int64_t nblocks = mio_read_metadata_num_cell_blocks(meta);
    if (nblocks < 0) nblocks = 0;

    mxArray *types = mio_o_cellstr(nblocks);
    mxArray *counts = mio_o_column(nblocks);
    mxArray *npcs = mio_o_column(nblocks);
    mxArray *ragged = mio_o_logical_column(nblocks);
    for (int64_t i = 0; i < nblocks; ++i) {
        int64_t nc = 0, npc = 0;
        int is_ragged = 0;
        if (mio_read_metadata_cell_block(meta, i, &nc, &npc, &is_ragged) != MIO_OK) {
            mio_read_metadata_free(meta);
            mio_o_fail("read_metadata cell block");
        }
        mxGetPr(counts)[i] = (double)nc;
        mxGetPr(npcs)[i] = (double)npc;
        mxGetLogicals(ragged)[i] = is_ragged != 0;
        meta_ctx ctx = {meta, i, 0};
        /* A getter failure here would leak `meta`; read the name through the
         * C protocol directly so the free can happen first. */
        int64_t len = mio_read_metadata_cell_block_type(meta, i, NULL, 0);
        if (len < 0) {
            mio_read_metadata_free(meta);
            mio_o_fail("block type");
        }
        (void)ctx;
        char *buf = (char *)mxCalloc((size_t)len + 1, 1);
        if (len > 0) meta_block_type_getter(&ctx, buf, len + 1);
        mxSetCell(types, (mwIndex)i, mxCreateString(buf));
        mxFree(buf);
    }

    double lo[3], hi[3];
    mxArray *bbox;
    if (mio_read_metadata_bbox(meta, lo, hi) == MIO_OK) {
        /* 2 x 3: row 1 the minimum corner, row 2 the maximum. */
        bbox = mio_o_matrix(2, 3);
        double *b = mxGetPr(bbox);
        for (int j = 0; j < 3; ++j) {
            b[j * 2] = lo[j];
            b[j * 2 + 1] = hi[j];
        }
    } else {
        bbox = mio_o_empty();
    }

    mxArray *pn = meta_names(meta, MIO_DATA_POINT);
    mxArray *cn = meta_names(meta, MIO_DATA_CELL);
    mxArray *fn = meta_names(meta, MIO_DATA_FIELD);

    mxArray *npoints = mio_o_scalar((double)mio_read_metadata_num_points(meta));
    mxArray *pdim = mio_o_scalar((double)mio_read_metadata_point_dim(meta));
    mxArray *ncells = mio_o_scalar((double)mio_read_metadata_num_cells(meta));
    mxArray *nblk = mio_o_scalar((double)nblocks);
    mxArray *fell = mio_o_logical(mio_read_metadata_fell_back(meta) == 1);

    int64_t nsteps = mio_read_metadata_num_time_values(meta);
    if (nsteps < 0) nsteps = 0;
    mxArray *times = mio_o_column(nsteps);
    if (nsteps > 0) mio_read_metadata_time_values(meta, mxGetPr(times), nsteps);

    mxArray *regions = meta_regions(meta);
    mio_read_metadata_free(meta);

    const char *names[] = {"num_points",       "point_dim",
                           "num_cells",        "num_cell_blocks",
                           "cell_block_types", "cell_block_num_cells",
                           "cell_block_nodes_per_cell", "cell_block_is_ragged",
                           "point_data_names", "cell_data_names",
                           "field_data_names", "bbox",
                           "fell_back",        "time_values",
                           "regions"};
    mxArray *values[] = {npoints, pdim, ncells, nblk, types, counts, npcs, ragged,
                         pn,      cn,   fn,     bbox, fell,  times,  regions};
    return mio_o_struct(15, names, values);
}

/* --- counts -------------------------------------------------------------- */

static mxArray *O_mio_num_points(const mxArray **a) {
    return mio_o_scalar((double)mio_mesh_num_points(mio_o_mesh(a[0])));
}

static mxArray *O_mio_point_dim(const mxArray **a) {
    return mio_o_scalar((double)mio_mesh_point_dim(mio_o_mesh(a[0])));
}

static mxArray *O_mio_num_cell_blocks(const mxArray **a) {
    return mio_o_scalar((double)mio_mesh_num_cell_blocks(mio_o_mesh(a[0])));
}

/* cell_block_info(mesh, block0) */
static mxArray *O_mio_cell_block_info(const mxArray **a) {
    mio_cell_block_info info;
    memset(&info, 0, sizeof(info));
    mio_o_check(mio_mesh_cell_block_info_ex(mio_o_mesh(a[0]), mio_o_int64(a[1], "block"), &info),
                "cell_block_info");
    const char *names[] = {"num_cells",     "nodes_per_cell", "is_ragged",
                           "is_polyhedron", "num_faces",      "num_nodes"};
    mxArray *values[] = {mio_o_scalar((double)info.num_cells),
                         mio_o_scalar((double)info.nodes_per_cell),
                         mio_o_logical(info.is_ragged != 0),
                         mio_o_logical(info.is_polyhedron != 0),
                         mio_o_scalar((double)info.num_faces),
                         mio_o_scalar((double)info.num_nodes)};
    return mio_o_struct(6, names, values);
}

typedef struct {
    const mio_mesh *mesh;
    int64_t index;
} mesh_ctx;

static int64_t block_type_getter(void *c, char *buf, int64_t buflen) {
    mesh_ctx *x = (mesh_ctx *)c;
    return mio_mesh_cell_block_type(x->mesh, x->index, buf, buflen);
}

static mxArray *O_mio_cell_block_type(const mxArray **a) {
    mesh_ctx ctx = {mio_o_mesh(a[0]), mio_o_int64(a[1], "block")};
    return mio_o_getstring(block_type_getter, &ctx, "cell_block_type");
}

/* --- points -------------------------------------------------------------- */

static mxArray *O_mio_points(const mxArray **a) {
    const mio_mesh *m = mio_o_mesh(a[0]);
    const void *data = NULL;
    mio_dtype dt;
    mio_o_check(mio_mesh_get_points(m, &data, &dt), "get_points");
    int64_t n = mio_mesh_num_points(m);
    int64_t dim = mio_mesh_point_dim(m);
    /* (num_points x dim) double: one transposing copy. */
    return mio_o_rows_from_c(data, dt, n, dim, 0.0);
}

/* Shared by the 1-based and 0-based connectivity accessors. */
static mxArray *conn_matrix(const mxArray **a, int shift) {
    const mio_mesh *m = mio_o_mesh(a[0]);
    int64_t b = mio_o_int64(a[1], "block");
    int64_t nc = 0, npc = 0;
    int32_t ragged = 0;
    mio_o_check(mio_mesh_cell_block_info(m, b, &nc, &npc, &ragged), "cell_block_info");
    if (ragged) {
        mio_o_error("ragged",
                    "cell block %d is ragged (polygons or polyhedra of varying size), so it "
                    "has no connectivity matrix; read it with mio.polygon_block() or "
                    "mio.polyhedron_block()",
                    (int)(b + 1));
    }
    const void *conn = NULL;
    mio_dtype dt;
    mio_o_check(mio_mesh_cell_block_conn(m, b, &conn, &dt), "cell_block_conn");
    /* (num_cells x nodes_per_cell), shifted while it is copied. */
    return mio_o_rows_from_c(conn, dt, nc, npc, shift ? 1.0 : 0.0);
}

static mxArray *O_mio_connectivity(const mxArray **a) { return conn_matrix(a, 1); }
static mxArray *O_mio_connectivity_raw(const mxArray **a) { return conn_matrix(a, 0); }

/* --- ragged (polygon / polyhedron) connectivity -------------------------- *
 *
 * The C ABI carries these as flat CSR; MATLAB has cell arrays, so that is what
 * comes back: one 1-based ROW vector per cell (one row of a connectivity
 * matrix), and for polyhedra one cell of face rows per cell. */

static mxArray *poly_span(const int64_t *nodes, int64_t from, int64_t to) {
    mxArray *out = mio_o_matrix(1, to - from);
    double *dst = mxGetPr(out);
    for (int64_t i = from; i < to; ++i) dst[i - from] = (double)(nodes[i] + 1);
    return out;
}

static mio_poly_conn *open_poly(const mio_mesh *m, int64_t b, mio_poly_conn_shape *shape) {
    mio_poly_conn *pc = (mio_poly_conn *)mio_o_check_ptr(mio_poly_conn_create(m, b), "poly_conn");
    memset(shape, 0, sizeof(*shape));
    if (mio_poly_conn_get_shape(pc, shape) != MIO_OK) {
        mio_poly_conn_free(pc);
        mio_o_fail("poly_conn_get_shape");
    }
    return pc;
}

static mxArray *O_mio_polygon_block(const mxArray **a) {
    const mio_mesh *m = mio_o_mesh(a[0]);
    int64_t b = mio_o_int64(a[1], "block");
    mio_poly_conn_shape shape;
    mio_poly_conn *pc = open_poly(m, b, &shape);
    if (shape.is_polyhedron) {
        mio_poly_conn_free(pc);
        mio_o_error("polygon_block", "cell block %d is a polyhedron block; use mio.polyhedron_block()",
                    (int)(b + 1));
    }
    const int64_t *nodes = mio_poly_conn_nodes(pc, NULL);
    const int64_t *rows = mio_poly_conn_face_offsets(pc, NULL);
    if (shape.num_cells > 0 && (nodes == NULL || rows == NULL)) {
        mio_poly_conn_free(pc);
        mio_o_fail("poly_conn arrays");
    }
    mxArray *out = mio_o_cell(shape.num_cells);
    for (int64_t c = 0; c < shape.num_cells; ++c) {
        mxSetCell(out, (mwIndex)c, poly_span(nodes, rows[c], rows[c + 1]));
    }
    mio_poly_conn_free(pc);
    return out;
}

static mxArray *O_mio_polyhedron_block(const mxArray **a) {
    const mio_mesh *m = mio_o_mesh(a[0]);
    int64_t b = mio_o_int64(a[1], "block");
    mio_poly_conn_shape shape;
    mio_poly_conn *pc = open_poly(m, b, &shape);
    if (!shape.is_polyhedron) {
        mio_poly_conn_free(pc);
        mio_o_error("polyhedron_block",
                    "cell block %d is a 1-level polygon block; use mio.polygon_block()",
                    (int)(b + 1));
    }
    const int64_t *nodes = mio_poly_conn_nodes(pc, NULL);
    const int64_t *faces = mio_poly_conn_face_offsets(pc, NULL);
    const int64_t *cells = mio_poly_conn_cell_offsets(pc, NULL);
    if (shape.num_cells > 0 && (nodes == NULL || faces == NULL || cells == NULL)) {
        mio_poly_conn_free(pc);
        mio_o_fail("poly_conn arrays");
    }
    mxArray *out = mio_o_cell(shape.num_cells);
    for (int64_t c = 0; c < shape.num_cells; ++c) {
        mxArray *cell = mio_o_cell(cells[c + 1] - cells[c]);
        for (int64_t f = cells[c]; f < cells[c + 1]; ++f) {
            mxSetCell(cell, (mwIndex)(f - cells[c]), poly_span(nodes, faces[f], faces[f + 1]));
        }
        mxSetCell(out, (mwIndex)c, cell);
    }
    mio_poly_conn_free(pc);
    return out;
}

/* Append the 0-based nodes of one numeric vector at `*at`. */
static void flatten_row(const mxArray *row, int64_t *nodes, int64_t *at, const char *what) {
    int64_t k = 0;
    int64_t *idx = mio_o_index_vector(row, &k, what);
    memcpy(nodes + *at, idx, (size_t)k * sizeof(int64_t));
    *at += k;
    mxFree(idx);
}

static mxArray *O_mio_add_polygon_block(const mxArray **a) {
    const mxArray *rows = a[2];
    if (!mxIsCell(rows)) {
        mio_o_error("add_polygon_block", "`rows` must be a cell array of numeric vectors, one per cell");
    }
    int64_t n = (int64_t)mxGetNumberOfElements(rows);
    int64_t total = 0;
    for (int64_t i = 0; i < n; ++i) {
        const mxArray *r = mxGetCell(rows, (mwIndex)i);
        total += r != NULL ? (int64_t)mxGetNumberOfElements(r) : 0;
    }
    int64_t *offsets = (int64_t *)mxMalloc((size_t)(n + 1) * sizeof(int64_t));
    int64_t *nodes = (int64_t *)mxMalloc((size_t)(total > 0 ? total : 1) * sizeof(int64_t));
    offsets[0] = 0;
    int64_t at = 0;
    for (int64_t i = 0; i < n; ++i) {
        const mxArray *r = mxGetCell(rows, (mwIndex)i);
        if (r == NULL) mio_o_error("add_polygon_block", "cell %d is empty", (int)(i + 1));
        flatten_row(r, nodes, &at, "connectivity");
        offsets[i + 1] = at;
    }
    mio_o_check(mio_mesh_add_polygon_block(mio_o_mesh(a[0]), mio_o_string(a[1], "cell_type"), n,
                                           offsets, nodes, total),
                "add_polygon_block");
    return NULL;
}

static mxArray *O_mio_add_polyhedron_block(const mxArray **a) {
    const mxArray *cells = a[2];
    if (!mxIsCell(cells)) {
        mio_o_error("add_polyhedron_block",
                    "`cells` must be a cell array of cells, each a cell array of face vectors");
    }
    int64_t ncells = (int64_t)mxGetNumberOfElements(cells);
    int64_t nfaces = 0, total = 0;
    for (int64_t c = 0; c < ncells; ++c) {
        const mxArray *cell = mxGetCell(cells, (mwIndex)c);
        if (cell == NULL || !mxIsCell(cell)) {
            mio_o_error("add_polyhedron_block", "cell %d must be a cell array of face vectors",
                        (int)(c + 1));
        }
        int64_t nf = (int64_t)mxGetNumberOfElements(cell);
        nfaces += nf;
        for (int64_t f = 0; f < nf; ++f) {
            const mxArray *face = mxGetCell(cell, (mwIndex)f);
            total += face != NULL ? (int64_t)mxGetNumberOfElements(face) : 0;
        }
    }
    int64_t *cell_offsets = (int64_t *)mxMalloc((size_t)(ncells + 1) * sizeof(int64_t));
    int64_t *face_offsets = (int64_t *)mxMalloc((size_t)(nfaces + 1) * sizeof(int64_t));
    int64_t *nodes = (int64_t *)mxMalloc((size_t)(total > 0 ? total : 1) * sizeof(int64_t));
    cell_offsets[0] = 0;
    face_offsets[0] = 0;
    int64_t face_at = 0, at = 0;
    for (int64_t c = 0; c < ncells; ++c) {
        const mxArray *cell = mxGetCell(cells, (mwIndex)c);
        int64_t nf = (int64_t)mxGetNumberOfElements(cell);
        for (int64_t f = 0; f < nf; ++f) {
            const mxArray *face = mxGetCell(cell, (mwIndex)f);
            if (face == NULL) {
                mio_o_error("add_polyhedron_block", "face %d of cell %d is empty", (int)(f + 1),
                            (int)(c + 1));
            }
            flatten_row(face, nodes, &at, "connectivity");
            face_offsets[++face_at] = at;
        }
        cell_offsets[c + 1] = face_at;
    }
    mio_o_check(mio_mesh_add_polyhedron_block(mio_o_mesh(a[0]), mio_o_string(a[1], "cell_type"),
                                              ncells, cell_offsets, nfaces, face_offsets, nodes,
                                              total),
                "add_polyhedron_block");
    return NULL;
}

/* --- data arrays --------------------------------------------------------- */

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

static mxArray *data_names(const mio_mesh *m, int64_t n, mio_o_str_getter getter) {
    if (n < 0) n = 0;
    mxArray *out = mio_o_cellstr(n);
    for (int64_t i = 0; i < n; ++i) {
        mesh_ctx ctx = {m, i};
        mxSetCell(out, (mwIndex)i, mio_o_getstring(getter, &ctx, "data name"));
    }
    return out;
}

static mxArray *O_mio_point_data_names(const mxArray **a) {
    const mio_mesh *m = mio_o_mesh(a[0]);
    return data_names(m, mio_mesh_num_point_data(m), point_data_name_getter);
}
static mxArray *O_mio_cell_data_names(const mxArray **a) {
    const mio_mesh *m = mio_o_mesh(a[0]);
    return data_names(m, mio_mesh_num_cell_data(m), cell_data_name_getter);
}
static mxArray *O_mio_field_data_names(const mxArray **a) {
    const mio_mesh *m = mio_o_mesh(a[0]);
    return data_names(m, mio_mesh_num_field_data(m), field_data_name_getter);
}

static mxArray *O_mio_cell_data_num_blocks(const mxArray **a) {
    int64_t n = mio_mesh_cell_data_num_blocks(mio_o_mesh(a[0]), mio_o_string(a[1], "name"));
    if (n < 0) mio_o_fail("cell_data_num_blocks");
    return mio_o_scalar((double)n);
}

/* A data array keeps its stored class and its C shape: (n) is n x 1, (n, k)
 * is n x k. */
static mxArray *O_mio_point_data(const mxArray **a) {
    const void *data = NULL;
    mio_dtype dt;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_o_check(mio_mesh_get_point_data(mio_o_mesh(a[0]), mio_o_string(a[1], "name"), &data, &dt,
                                        &ndim, shape),
                "get_point_data");
    return mio_o_typed_from_c(data, dt, ndim, shape);
}

static mxArray *O_mio_cell_data(const mxArray **a) {
    const void *data = NULL;
    mio_dtype dt;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_o_check(mio_mesh_get_cell_data(mio_o_mesh(a[0]), mio_o_string(a[1], "name"),
                                       mio_o_int64(a[2], "block"), &data, &dt, &ndim, shape),
                "get_cell_data");
    return mio_o_typed_from_c(data, dt, ndim, shape);
}

static mxArray *O_mio_field_data(const mxArray **a) {
    const void *data = NULL;
    mio_dtype dt;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_o_check(mio_mesh_get_field_data(mio_o_mesh(a[0]), mio_o_string(a[1], "name"), &data, &dt,
                                        &ndim, shape),
                "get_field_data");
    return mio_o_typed_from_c(data, dt, ndim, shape);
}

/* --- setters (all COPY caller memory) ------------------------------------ */

static mxArray *O_mio_set_points(const mxArray **a) {
    int64_t n = 0, d = 0;
    double *p = mio_o_rows_to_c(a[1], &n, &d, 0.0, "points");
    mio_o_check(mio_mesh_set_points(mio_o_mesh(a[0]), MIO_FLOAT64, n, d, p), "set_points");
    return NULL;
}

static mxArray *O_mio_add_cell_block(const mxArray **a) {
    int64_t nc = 0, npc = 0;
    int64_t *conn = mio_o_index_rows(a[2], &nc, &npc, "connectivity");
    mio_o_check(mio_mesh_add_cell_block(mio_o_mesh(a[0]), mio_o_string(a[1], "cell_type"), nc,
                                        npc, MIO_INT64, conn),
                "add_cell_block");
    return NULL;
}

static mxArray *add_data(const mxArray **a, int which) {
    mio_dtype dt;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    void *buf = mio_o_data_to_c(a[2], &dt, &ndim, shape, "data");
    const char *nm = mio_o_string(a[1], "name");
    mio_mesh *m = mio_o_mesh(a[0]);
    mio_status st;
    if (which == 0) {
        st = mio_mesh_add_point_data(m, nm, dt, ndim, shape, buf);
    } else if (which == 1) {
        st = mio_mesh_append_cell_data(m, nm, dt, ndim, shape, buf);
    } else {
        st = mio_mesh_add_field_data(m, nm, dt, ndim, shape, buf);
    }
    mio_o_check(st, "add_data");
    return NULL;
}

static mxArray *O_mio_add_point_data(const mxArray **a) { return add_data(a, 0); }
static mxArray *O_mio_append_cell_data(const mxArray **a) { return add_data(a, 1); }
static mxArray *O_mio_add_field_data(const mxArray **a) { return add_data(a, 2); }

/* --- named regions ------------------------------------------------------- */

typedef struct {
    const mio_regions *regions;
    int64_t index;
} regions_ctx;

static int64_t region_name_getter(void *c, char *buf, int64_t buflen) {
    regions_ctx *x = (regions_ctx *)c;
    return mio_regions_name(x->regions, x->index, buf, buflen);
}

static const char *kRegionFields[] = {"name", "kind", "dim", "tag", "entries"};

/* The mesh's regions as an n x 1 struct array. `entries` is
 * (num_entries x stride): 1-based indices, except a side region's second
 * column, which is a FACET ordinal within the cell type and is NOT shifted --
 * the rule every binding follows. */
static mxArray *O_mio_regions(const mxArray **a) {
    mio_regions *regions =
        (mio_regions *)mio_o_check_ptr(mio_regions_create(mio_o_mesh(a[0])), "regions");
    int64_t n = mio_regions_count(regions);
    if (n < 0) {
        mio_regions_free(regions);
        mio_o_fail("regions_count");
    }
    mxArray *out = mio_o_struct_array(n, 5, kRegionFields);
    for (int64_t i = 0; i < n; ++i) {
        regions_ctx ctx = {regions, i};
        int64_t len = mio_regions_name(regions, i, NULL, 0);
        mio_region_info info;
        if (len < 0 || mio_regions_info(regions, i, &info) != MIO_OK) {
            mio_regions_free(regions);
            mio_o_fail("regions_info");
        }
        char *nm = (char *)mxCalloc((size_t)len + 1, 1);
        if (len > 0) region_name_getter(&ctx, nm, len + 1);
        mxArray *entries = mio_o_matrix(info.num_entries, info.stride);
        if (info.num_entries > 0) {
            int64_t count = 0;
            const int64_t *src = mio_regions_entries(regions, i, &count);
            if (src == NULL) {
                mio_regions_free(regions);
                mio_o_fail("regions_entries");
            }
            double *dst = mxGetPr(entries);
            for (int64_t e = 0; e < info.num_entries; ++e) {
                for (int64_t r = 0; r < info.stride; ++r) {
                    int64_t v = src[e * info.stride + r];
                    int is_facet = (info.kind == MIO_REGION_SIDE && r == 1);
                    dst[r * info.num_entries + e] = is_facet ? (double)v : (double)(v + 1);
                }
            }
        }
        mxSetField(out, (mwIndex)i, "name", mxCreateString(nm));
        mxFree(nm);
        mxSetField(out, (mwIndex)i, "kind", mio_o_str(region_kind_name(info.kind)));
        mxSetField(out, (mwIndex)i, "dim", mio_o_scalar((double)info.dim));
        mxSetField(out, (mwIndex)i, "tag", mio_o_scalar((double)info.tag));
        mxSetField(out, (mwIndex)i, "entries", entries);
    }
    mio_regions_free(regions);
    return out;
}

static mio_region_kind parse_kind(const char *k) {
    if (strcmp(k, "point") == 0) return MIO_REGION_POINT;
    if (strcmp(k, "cell") == 0) return MIO_REGION_CELL;
    if (strcmp(k, "side") == 0) return MIO_REGION_SIDE;
    mio_o_error("region_kind", "`kind` must be 'point', 'cell' or 'side', got '%s'", k);
    return MIO_REGION_POINT;
}

/* add_region(mesh, name, kind, entries, dim, tag) -- entries a vector, or an
 * (n x stride) matrix. */
static mxArray *O_mio_add_region(const mxArray **a) {
    const char *k = mio_o_string(a[2], "kind");
    mio_region_kind kv = parse_kind(k);
    const mxArray *e = a[3];
    if (!(mxIsNumeric(e) || mxIsLogical(e)) || mxGetNumberOfDimensions(e) > 2) {
        mio_o_error("add_region", "`entries` must be a vector or an (n x stride) matrix");
    }
    int64_t expected = (kv == MIO_REGION_SIDE) ? 2 : 1;
    int64_t rows = (int64_t)mxGetM(e), cols = (int64_t)mxGetN(e);
    int64_t nent, stride;
    if (expected == 1 && (rows == 1 || cols == 1)) {
        nent = rows * cols;
        stride = 1;
    } else {
        nent = rows;
        stride = cols;
    }
    if (mxIsEmpty(e)) {
        nent = 0;
        stride = expected;
    }
    if (stride != expected) {
        mio_o_error("add_region", "a '%s' region needs %d value(s) per entry, got %d", k,
                    (int)expected, (int)stride);
    }
    size_t total_n = 0;
    double *src = mio_o_doubles(e, &total_n, "entries");
    int64_t total = stride * nent;
    int64_t *dst = (int64_t *)mxMalloc((size_t)(total > 0 ? total : 1) * sizeof(int64_t));
    for (int64_t i = 0; i < nent; ++i) {
        for (int64_t r = 0; r < stride; ++r) {
            double v = src[r * nent + i]; /* column-major (nent x stride) */
            int is_facet = (kv == MIO_REGION_SIDE && r == 1);
            if (is_facet) {
                dst[i * stride + r] = (int64_t)v;
            } else {
                if (!isfinite(v) || v < 1.0) {
                    mio_o_error("add_region", "region entries are 1-based here; got %g", v);
                }
                dst[i * stride + r] = (int64_t)v - 1;
            }
        }
    }
    mio_o_check(mio_mesh_add_region(mio_o_mesh(a[0]), mio_o_string(a[1], "name"), kv,
                                    (int32_t)mio_o_int(a[4], "dim"), mio_o_int64(a[5], "tag"), dst,
                                    total),
                "add_region");
    return NULL;
}

/* --- region editing and periodic pairs ----------------------------------- */

static int32_t region_kind_or_any(const mxArray *kind) {
    const char *k = mio_o_opt_string(kind);
    if (k == NULL || k[0] == '\0') return -1;
    return (int32_t)parse_kind(k);
}

/* edit_regions(mesh, op, inputs, output, kind, dim, tag, keep_inputs) */
static mxArray *O_mio_edit_regions(const mxArray **a) {
    static const char *const ops[] = {"union",  "intersection", "difference",
                                      "rename", "retag",        "delete"};
    const char *o = mio_o_string(a[1], "op");
    int32_t code = -1;
    for (int i = 0; i < 6; ++i)
        if (strcmp(o, ops[i]) == 0) code = i;
    if (strcmp(o, "intersect") == 0) code = 1;
    if (code < 0) {
        mio_o_error("edit_regions",
                    "meshio++: edit_regions: unknown operation '%s' (expected union, "
                    "intersection, difference, rename, retag or delete)",
                    o);
    }
    const int32_t k = region_kind_or_any(a[4]);
    int64_t n = 0;
    const char *const *names = mio_o_names(a[2], &n, "inputs");
    mio_region_selector *sel =
        (mio_region_selector *)mxCalloc(n > 0 ? (size_t)n : 1, sizeof(mio_region_selector));
    for (int64_t i = 0; i < n; ++i) {
        mio_region_selector_init(&sel[i], names[i]);
        sel[i].kind = k;
    }
    const char *out_name = mio_o_opt_string(a[3]);
    if (out_name != NULL && out_name[0] == '\0') out_name = NULL;
    const int64_t d = mio_o_is_null(a[5]) ? MIO_REGION_ANY : (int64_t)mio_o_double(a[5], "dim");
    const int64_t t = mio_o_is_null(a[6]) ? MIO_REGION_ANY : (int64_t)mio_o_double(a[6], "tag");
    mio_mesh *out = mio_edit_regions(mio_o_mesh(a[0]), code, sel, n, out_name, d, t,
                                     mio_o_bool(a[7], "keep_inputs") ? 1 : 0);
    if (out == NULL) mio_o_fail("edit_regions");
    return mio_o_wrap_mesh(out);
}

/* remove_region(mesh, index1) */
static mxArray *O_mio_remove_region(const mxArray **a) {
    const int64_t i = (int64_t)mio_o_double(a[1], "index");
    mio_o_check(mio_mesh_remove_region(mio_o_mesh(a[0]), i - 1), "remove_region");
    return NULL;
}

/* match_periodic_nodes(mesh, slave, master, matrix4x4, atol, require_complete) */
static mxArray *O_mio_match_periodic_nodes(const mxArray **a) {
    mio_periodic_opts opts;
    mio_periodic_opts_init(&opts);
    double m[16];
    mio_o_vecn(a[3], m, 16, "matrix");
    /* Row-major on the C side, column-major here: transpose (see transform). */
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) opts.matrix[i * 4 + j] = m[j * 4 + i];
    opts.atol = mio_o_double(a[4], "atol");
    opts.require_complete = mio_o_bool(a[5], "require_complete") ? 1 : 0;
    mio_region_selector s, t;
    mio_region_selector_init(&s, mio_o_string(a[1], "slave"));
    mio_region_selector_init(&t, mio_o_string(a[2], "master"));
    mio_periodic_pairs *pairs = mio_match_periodic_nodes(mio_o_mesh(a[0]), &s, &t, &opts);
    if (pairs == NULL) mio_o_fail("match_periodic_nodes");
    int64_t np = 0, nu = 0, nf = 0, cnt = 0;
    double res = 0.0;
    if (mio_periodic_pairs_info(pairs, &np, &nu, &nf, &res) != MIO_OK) {
        mio_periodic_pairs_free(pairs);
        mio_o_fail("match_periodic_nodes");
    }
    const int64_t *sp = mio_periodic_pairs_slave(pairs, &cnt);
    mxArray *vs = mio_o_shift_map(sp, cnt);
    const int64_t *mp = mio_periodic_pairs_master(pairs, &cnt);
    mxArray *vm = mio_o_shift_map(mp, cnt);
    const int64_t *up = mio_periodic_pairs_unmatched(pairs, &cnt);
    mxArray *vu = mio_o_shift_map(up, cnt);
    mio_periodic_pairs_free(pairs);
    const char *names[] = {"slave", "master", "unmatched", "num_fixed", "max_residual"};
    mxArray *values[] = {vs, vm, vu, mio_o_scalar((double)nf), mio_o_scalar(res)};
    return mio_o_struct(5, names, values);
}

static mxArray *O_mio_adjacency(const mxArray **a) {
    mio_adjacency *g = mio_adjacency_create(mio_o_mesh(a[0]), mio_o_string(a[1], "kind"));
    if (g == NULL) mio_o_fail("adjacency");
    int64_t cnt = 0;
    const int64_t *pp = mio_adjacency_indptr(g, &cnt);
    mxArray *vp = mio_o_shift_map(pp, cnt);
    const int64_t *pi = mio_adjacency_indices(g, &cnt);
    mxArray *vi = mio_o_shift_map(pi, cnt);
    mio_adjacency_free(g);
    const char *names[] = {"indptr", "indices"};
    mxArray *values[] = {vp, vi};
    return mio_o_struct(2, names, values);
}

/* write_graph(mesh, path, kind, format) */
static mxArray *O_mio_write_graph(const mxArray **a) {
    int64_t nv = 0, ne = 0, md = 0;
    mio_o_check(mio_write_graph(mio_o_mesh(a[0]), mio_o_string(a[2], "kind"),
                                mio_o_string(a[1], "path"), mio_o_opt_string(a[3]), &nv, &ne, &md),
                "write_graph");
    const char *names[] = {"num_vertices", "num_edges", "max_degree"};
    mxArray *values[] = {mio_o_scalar((double)nv), mio_o_scalar((double)ne),
                         mio_o_scalar((double)md)};
    return mio_o_struct(3, names, values);
}

/* --- table --------------------------------------------------------------- */

static const mio_o_entry kEntries[] = {
    MIO_O_ENTRY(mio_version, 0),
    MIO_O_ENTRY(mio_mesh_backend, 0),
    MIO_O_ENTRY(mio_last_error, 0),
    MIO_O_ENTRY(mio_load_plugins, 1),
    MIO_O_ENTRY(mio_format_readable, 1),
    MIO_O_ENTRY(mio_format_writable, 1),
    MIO_O_ENTRY(mio_sniff_format, 1),
    MIO_O_ENTRY(mio_cell_type_num_nodes, 1),
    MIO_O_ENTRY(mio_cell_type_dimension, 1),
    MIO_O_ENTRY(mio_reader_supports_options, 1),
    MIO_O_ENTRY(mio_mesh_create, 0),
    MIO_O_ENTRY(mio_mesh_release, 1),
    MIO_O_ENTRY(mio_mesh_is_open, 1),
    MIO_O_ENTRY(mio_read, 10),
    MIO_O_ENTRY(mio_write, 6),
    MIO_O_ENTRY(mio_convert, 4),
    MIO_O_ENTRY(mio_write_gltf, 3),
    MIO_O_ENTRY(mio_med_mesh_names, 1),
    MIO_O_ENTRY(mio_med_read_named, 4),
    MIO_O_ENTRY(mio_med_write_multi, 4),
    MIO_O_ENTRY(mio_read_metadata, 2),
    MIO_O_ENTRY(mio_num_points, 1),
    MIO_O_ENTRY(mio_point_dim, 1),
    MIO_O_ENTRY(mio_num_cell_blocks, 1),
    MIO_O_ENTRY(mio_cell_block_info, 2),
    MIO_O_ENTRY(mio_cell_block_type, 2),
    MIO_O_ENTRY(mio_points, 1),
    MIO_O_ENTRY(mio_connectivity, 2),
    MIO_O_ENTRY(mio_connectivity_raw, 2),
    MIO_O_ENTRY(mio_polygon_block, 2),
    MIO_O_ENTRY(mio_polyhedron_block, 2),
    MIO_O_ENTRY(mio_add_polygon_block, 3),
    MIO_O_ENTRY(mio_add_polyhedron_block, 3),
    MIO_O_ENTRY(mio_point_data_names, 1),
    MIO_O_ENTRY(mio_cell_data_names, 1),
    MIO_O_ENTRY(mio_field_data_names, 1),
    MIO_O_ENTRY(mio_cell_data_num_blocks, 2),
    MIO_O_ENTRY(mio_point_data, 2),
    MIO_O_ENTRY(mio_cell_data, 3),
    MIO_O_ENTRY(mio_field_data, 2),
    MIO_O_ENTRY(mio_set_points, 2),
    MIO_O_ENTRY(mio_add_cell_block, 3),
    MIO_O_ENTRY(mio_add_point_data, 3),
    MIO_O_ENTRY(mio_append_cell_data, 3),
    MIO_O_ENTRY(mio_add_field_data, 3),
    MIO_O_ENTRY(mio_regions, 1),
    MIO_O_ENTRY(mio_add_region, 6),
    MIO_O_ENTRY(mio_edit_regions, 8),
    MIO_O_ENTRY(mio_remove_region, 2),
    MIO_O_ENTRY(mio_match_periodic_nodes, 6),
    MIO_O_ENTRY(mio_adjacency, 2),
    MIO_O_ENTRY(mio_write_graph, 4),
};

const mio_o_entry *mio_o_table_mesh(size_t *n) {
    *n = sizeof(kEntries) / sizeof(kEntries[0]);
    return kEntries;
}
