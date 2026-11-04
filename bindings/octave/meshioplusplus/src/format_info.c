/* Format side channels: what a format carries that a mesh cannot hold (MDPA's
 * tables, geometries, Mesh blocks and constraints), kept by
 * `mio.read_with_info()` in its own handle kind (mio.FormatInfo) -- own
 * registry kind, own free function, exactly as mio_mesh -- for
 * `mio.write_with_info()` to put back.
 *
 * `mio.mdpa_info()` copies the whole MDPA side channel into plain MATLAB
 * values in one call, so nothing but the handle ever outlives a MEX call.
 *
 * A port of the R binding's src/format_info.c; see mio_oct.h for the layout
 * rules. Every per-entity array is (entities x components): where R hands a
 * C (rows, cols) array back as a (cols, rows) matrix (the connectivity
 * convention) or as a (rows, cols) table, both are (rows x cols) here. */

#include "mio_oct.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* NULL for [] (no side channel); an error for anything else foreign. */
static const mio_format_info *ofi_info_of(const mxArray *x, int allow_null) {
    if (allow_null && mio_o_is_null(x)) return NULL;
    return (const mio_format_info *)mio_o_unwrap(x, MIO_O_FORMAT_INFO, "info");
}

static mxArray *ofi_wrap_info(mio_format_info *info) {
    return info != NULL ? mio_o_wrap(MIO_O_FORMAT_INFO, info) : mio_o_empty();
}

/* struct('mesh', handle, 'info', handle or []) */
static mxArray *ofi_mesh_and_info(mio_mesh *mesh, mio_format_info *info, const char *what) {
    if (mesh == NULL) {
        if (info != NULL) mio_format_info_free(info);
        mio_o_fail(what);
    }
    mxArray *values[2];
    values[0] = mio_o_wrap_mesh(mesh);
    values[1] = ofi_wrap_info(info);
    const char *names[2] = {"mesh", "info"};
    return mio_o_struct(2, names, values);
}

/* read_with_info(path, format, lenient, deflection, angular_deflection, heal,
 *                sew_tolerance) */
static mxArray *O_mio_read_with_info(const mxArray **a) {
    mio_read_opts opts;
    mio_read_opts_init(&opts);
    opts.lenient = mio_o_bool(a[2], "lenient");
    mio_cad_read_opts cad;
    mio_cad_read_opts_init(&cad);
    cad.deflection = mio_o_double(a[3], "deflection");
    cad.angular_deflection = mio_o_double(a[4], "angular_deflection");
    cad.heal = mio_o_bool(a[5], "heal");
    cad.sew_tolerance = mio_o_double(a[6], "sew_tolerance");
    mio_format_info *info = NULL;
    mio_mesh *m = mio_read_with_info2(mio_o_string(a[0], "path"), mio_o_opt_string(a[1]), &opts,
                                      &cad, &info);
    return ofi_mesh_and_info(m, info, "read_with_info");
}

/* write_with_info(mesh, info, path, format) */
static mxArray *O_mio_write_with_info(const mxArray **a) {
    mio_o_check(mio_write_with_info(mio_o_string(a[2], "path"), mio_o_mesh(a[0]),
                                    mio_o_opt_string(a[3]), ofi_info_of(a[1], 1)),
                "write_with_info");
    return NULL;
}

/* Idempotent: an already-released handle is a no-op; a foreign object raises,
 * as in R. */
static mxArray *O_mio_format_info_release(const mxArray **a) {
    mio_o_release(a[0], MIO_O_FORMAT_INFO);
    return NULL;
}

static mxArray *O_mio_format_info_is_open(const mxArray **a) {
    return mio_o_logical(mio_o_is_open(a[0], MIO_O_FORMAT_INFO));
}

/* --- string getters ------------------------------------------------------ */

typedef struct {
    const mio_format_info *info;
    int32_t section;
    int64_t index;
    int64_t entry;
    int32_t field;
} ofi_mdpa_ctx;

static int64_t ofi_get_format(void *ctx, char *buf, int64_t buflen) {
    return mio_format_info_format(((ofi_mdpa_ctx *)ctx)->info, buf, buflen);
}

static int64_t ofi_get_item_string(void *ctx, char *buf, int64_t buflen) {
    ofi_mdpa_ctx *c = (ofi_mdpa_ctx *)ctx;
    return mio_mdpa_info_string(c->info, c->section, c->index, c->field, buf, buflen);
}

static int64_t ofi_get_data_string(void *ctx, char *buf, int64_t buflen) {
    ofi_mdpa_ctx *c = (ofi_mdpa_ctx *)ctx;
    return mio_mdpa_info_data_string(c->info, c->section, c->index, c->entry, c->field, buf,
                                     buflen);
}

static mxArray *O_mio_format_info_format(const mxArray **a) {
    ofi_mdpa_ctx c = {ofi_info_of(a[0], 0), 0, 0, 0, 0};
    return mio_o_getstring(ofi_get_format, &c, "format_info_format");
}

static mxArray *ofi_item_string(const mio_format_info *info, int32_t section, int64_t index,
                                int32_t field) {
    ofi_mdpa_ctx c = {info, section, index, 0, field};
    return mio_o_getstring(ofi_get_item_string, &c, "mdpa_info");
}

static int64_t ofi_item_int(const mio_format_info *info, int32_t section, int64_t index,
                            int32_t field) {
    int64_t v = 0;
    mio_o_check(mio_mdpa_info_int(info, section, index, field, &v), "mdpa_info");
    return v;
}

static int64_t ofi_count_of(const mio_format_info *info, int32_t section) {
    int64_t n = mio_mdpa_info_count(info, section);
    if (n < 0) mio_o_fail("mdpa_info");
    return n;
}

/* --- arrays -------------------------------------------------------------- */

/* A borrowed C array as a MATLAB double: a column for 1-D, (rows x cols) for
 * 2-D -- rows are entities, so R's two conventions (the (cols, rows)
 * connectivity matrix and the (rows, cols) table) are one here. Higher ranks
 * come back flat, as in R. `shift` is added to every value (1 for point
 * references). */
static mxArray *ofi_array(const void *data, mio_dtype dtype, int32_t ndim, const int64_t *shape,
                          double shift) {
    if (ndim == 2) return mio_o_rows_from_c(data, dtype, shape[0], shape[1], shift);
    int64_t n = 1;
    for (int32_t k = 0; k < ndim; ++k) n *= shape[k];
    if (ndim == 0) n = 0;
    mxArray *out = mio_o_column(n);
    if (n > 0) mio_o_fill_real(mxGetPr(out), data, dtype, (size_t)n, shift);
    return out;
}

static mxArray *ofi_item_array(const mio_format_info *info, int32_t section, int64_t index,
                               int32_t field, double shift) {
    const void *data = NULL;
    mio_dtype dtype = MIO_FLOAT64;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_o_check(mio_mdpa_info_array(info, section, index, field, &data, &dtype, &ndim, shape),
                "mdpa_info");
    return ofi_array(data, dtype, ndim, shape, shift);
}

/* A struct field name for an MDPA key. Keys are Kratos variable names, which
 * are valid identifiers; anything else is made one (invalid characters become
 * '_', a leading non-letter gets an 'x' prefix, at most 63 characters) and a
 * repeated key gets a _2, _3, ... suffix, since a struct cannot hold two fields
 * of one name the way an R named list can. */
static char *ofi_field_name(const char *key, char **taken, int64_t ntaken) {
    size_t len = strlen(key);
    char *base = (char *)mxCalloc(len + 2, 1);
    size_t at = 0;
    int lead_ok = (key[0] >= 'a' && key[0] <= 'z') || (key[0] >= 'A' && key[0] <= 'Z');
    if (!lead_ok) base[at++] = 'x';
    for (size_t i = 0; i < len; ++i) {
        char ch = key[i];
        int ok = (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') || (ch >= '0' && ch <= '9') ||
                 ch == '_';
        base[at++] = ok ? ch : '_';
    }
    base[at] = '\0';
    if (at > 63) base[63] = '\0';
    char *name = (char *)mxCalloc(strlen(base) + 24, 1);
    snprintf(name, strlen(base) + 24, "%s", base);
    for (int suffix = 2;; ++suffix) {
        int clash = 0;
        for (int64_t k = 0; k < ntaken; ++k) {
            if (strcmp(taken[k], name) == 0) {
                clash = 1;
                break;
            }
        }
        if (!clash) break;
        char tail[24];
        snprintf(tail, sizeof(tail), "_%d", suffix);
        size_t room = 63 - strlen(tail);
        snprintf(name, strlen(base) + 24, "%.*s%s", (int)room, base, tail);
    }
    mxFree(base);
    return name;
}

/* The key/value entries of one item, as a 1x1 struct (number, string, or a
 * (rows x cols) matrix for an inline table). */
static mxArray *ofi_item_data(const mio_format_info *info, int32_t section, int64_t index) {
    int64_t n = mio_mdpa_info_data_count(info, section, index);
    if (n < 0) mio_o_fail("mdpa_info");
    char **names = (char **)mxCalloc((size_t)(n > 0 ? n : 1), sizeof(char *));
    mxArray **values = (mxArray **)mxCalloc((size_t)(n > 0 ? n : 1), sizeof(mxArray *));
    for (int64_t e = 0; e < n; ++e) {
        ofi_mdpa_ctx c = {info, section, index, e, 0};
        char *key = mio_o_getcstring(ofi_get_data_string, &c, "mdpa_info");
        names[e] = ofi_field_name(key, names, e);
        mxFree(key);
        const int32_t kind = mio_mdpa_info_data_kind(info, section, index, e);
        if (kind < 0) mio_o_fail("mdpa_info");
        if (kind == MIO_MDPA_VALUE_TEXT) {
            c.field = 1;
            values[e] = mio_o_getstring(ofi_get_data_string, &c, "mdpa_info");
        } else {
            const void *data = NULL;
            mio_dtype dtype = MIO_FLOAT64;
            int32_t ndim = 0;
            int64_t shape[MIO_MAX_NDIM] = {0};
            mio_o_check(mio_mdpa_info_data_array(info, section, index, e, &data, &dtype, &ndim,
                                                 shape),
                        "mdpa_info");
            values[e] = ofi_array(data, dtype, ndim, shape, 0.0);
        }
    }
    return mio_o_struct((int)n, (const char **)names, values);
}

/* --- MDPA ---------------------------------------------------------------- */

static mxArray *O_mio_mdpa_info(const mxArray **a) {
    const mio_format_info *info = ofi_info_of(a[0], 0);
    int64_t n;

    static const char *kProp[] = {"id", "values"};
    n = ofi_count_of(info, MIO_MDPA_PROPERTIES);
    mxArray *props = mio_o_struct_array(n, 2, kProp);
    for (int64_t i = 0; i < n; ++i) {
        mxSetField(props, (mwIndex)i, "id",
                   mio_o_scalar((double)ofi_item_int(info, MIO_MDPA_PROPERTIES, i, 0)));
        mxSetField(props, (mwIndex)i, "values", ofi_item_data(info, MIO_MDPA_PROPERTIES, i));
    }

    static const char *kEntity[] = {"name", "is_condition"};
    n = ofi_count_of(info, MIO_MDPA_ENTITY_NAMES);
    mxArray *entities = mio_o_struct_array(n, 2, kEntity);
    for (int64_t i = 0; i < n; ++i) {
        mxSetField(entities, (mwIndex)i, "name",
                   ofi_item_string(info, MIO_MDPA_ENTITY_NAMES, i, 0));
        mxSetField(entities, (mwIndex)i, "is_condition",
                   mio_o_logical(ofi_item_int(info, MIO_MDPA_ENTITY_NAMES, i, 0) == 1));
    }

    n = ofi_count_of(info, MIO_MDPA_SKIPPED);
    mxArray *skipped = mio_o_cellstr(n);
    for (int64_t i = 0; i < n; ++i) {
        mxSetCell(skipped, (mwIndex)i, ofi_item_string(info, MIO_MDPA_SKIPPED, i, 0));
    }

    mxArray *model_part_data = ofi_item_data(info, MIO_MDPA_MODEL_PART_DATA, 0);

    static const char *kTable[] = {"header", "values"};
    n = ofi_count_of(info, MIO_MDPA_TABLES);
    mxArray *tables = mio_o_struct_array(n, 2, kTable);
    for (int64_t i = 0; i < n; ++i) {
        mxSetField(tables, (mwIndex)i, "header", ofi_item_string(info, MIO_MDPA_TABLES, i, 0));
        mxSetField(tables, (mwIndex)i, "values",
                   ofi_item_array(info, MIO_MDPA_TABLES, i, 0, 0.0));
    }

    static const char *kGeom[] = {"name", "type", "connectivity", "ids"};
    n = ofi_count_of(info, MIO_MDPA_GEOMETRIES);
    mxArray *geoms = mio_o_struct_array(n, 4, kGeom);
    for (int64_t i = 0; i < n; ++i) {
        mxSetField(geoms, (mwIndex)i, "name", ofi_item_string(info, MIO_MDPA_GEOMETRIES, i, 0));
        mxSetField(geoms, (mwIndex)i, "type", ofi_item_string(info, MIO_MDPA_GEOMETRIES, i, 1));
        mxSetField(geoms, (mwIndex)i, "connectivity", /* 1-based */
                   ofi_item_array(info, MIO_MDPA_GEOMETRIES, i, 0, 1.0));
        mxSetField(geoms, (mwIndex)i, "ids",
                   ofi_item_array(info, MIO_MDPA_GEOMETRIES, i, 1, 0.0));
    }

    static const char *kMesh[] = {"id", "data", "nodes", "element_ids", "condition_ids"};
    n = ofi_count_of(info, MIO_MDPA_MESH_BLOCKS);
    mxArray *meshes = mio_o_struct_array(n, 5, kMesh);
    for (int64_t i = 0; i < n; ++i) {
        mxSetField(meshes, (mwIndex)i, "id",
                   mio_o_scalar((double)ofi_item_int(info, MIO_MDPA_MESH_BLOCKS, i, 0)));
        mxSetField(meshes, (mwIndex)i, "data", ofi_item_data(info, MIO_MDPA_MESH_BLOCKS, i));
        mxSetField(meshes, (mwIndex)i, "nodes", /* 1-based */
                   ofi_item_array(info, MIO_MDPA_MESH_BLOCKS, i, 0, 1.0));
        mxSetField(meshes, (mwIndex)i, "element_ids",
                   ofi_item_array(info, MIO_MDPA_MESH_BLOCKS, i, 1, 0.0));
        mxSetField(meshes, (mwIndex)i, "condition_ids",
                   ofi_item_array(info, MIO_MDPA_MESH_BLOCKS, i, 2, 0.0));
    }

    static const char *kSmp[] = {"name", "data", "tables", "geometry_ids", "constraint_ids"};
    n = ofi_count_of(info, MIO_MDPA_SUBMODELPARTS);
    mxArray *smps = mio_o_struct_array(n, 5, kSmp);
    for (int64_t i = 0; i < n; ++i) {
        mxSetField(smps, (mwIndex)i, "name", ofi_item_string(info, MIO_MDPA_SUBMODELPARTS, i, 0));
        mxSetField(smps, (mwIndex)i, "data", ofi_item_data(info, MIO_MDPA_SUBMODELPARTS, i));
        mxSetField(smps, (mwIndex)i, "tables",
                   ofi_item_array(info, MIO_MDPA_SUBMODELPARTS, i, 0, 0.0));
        mxSetField(smps, (mwIndex)i, "geometry_ids",
                   ofi_item_array(info, MIO_MDPA_SUBMODELPARTS, i, 1, 0.0));
        mxSetField(smps, (mwIndex)i, "constraint_ids",
                   ofi_item_array(info, MIO_MDPA_SUBMODELPARTS, i, 2, 0.0));
    }

    static const char *kRaw[] = {"header", "body", "terminator"};
    n = ofi_count_of(info, MIO_MDPA_RAW_BLOCKS);
    mxArray *raws = mio_o_struct_array(n, 3, kRaw);
    for (int64_t i = 0; i < n; ++i) {
        mxSetField(raws, (mwIndex)i, "header", ofi_item_string(info, MIO_MDPA_RAW_BLOCKS, i, 0));
        mxSetField(raws, (mwIndex)i, "body", ofi_item_string(info, MIO_MDPA_RAW_BLOCKS, i, 1));
        mxSetField(raws, (mwIndex)i, "terminator",
                   ofi_item_string(info, MIO_MDPA_RAW_BLOCKS, i, 2));
    }

    const char *names[9] = {"properties",      "entity_names", "skipped",
                            "model_part_data", "tables",       "geometries",
                            "mesh_blocks",     "submodelparts", "raw_blocks"};
    mxArray *values[9] = {props, entities, skipped, model_part_data, tables,
                          geoms, meshes,   smps,    raws};
    return mio_o_struct(9, names, values);
}

/* --- Gmsh ---------------------------------------------------------------- */

static mxArray *ofi_gmsh_array(const mio_format_info *info, int32_t section, int64_t i,
                               int32_t field) {
    const void *data = NULL;
    mio_dtype dtype = MIO_FLOAT64;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_o_check(mio_gmsh_info_array(info, section, i, field, &data, &dtype, &ndim, shape),
                "gmsh_info");
    /* node_pairs (field 2) are point rows: 1-based. */
    return ofi_array(data, dtype, ndim, shape, field == 2 ? 1.0 : 0.0);
}

static mxArray *O_mio_gmsh_info(const mxArray **a) {
    const mio_format_info *info = ofi_info_of(a[0], 0);
    int64_t nb = mio_gmsh_info_count(info, MIO_GMSH_BOUNDING_ENTITIES);
    int64_t np = mio_gmsh_info_count(info, MIO_GMSH_PERIODIC);
    if (nb < 0 || np < 0) mio_o_fail("gmsh_info");
    mxArray *bounding = mio_o_cell(nb);
    for (int64_t i = 0; i < nb; ++i) {
        mxSetCell(bounding, (mwIndex)i, ofi_gmsh_array(info, MIO_GMSH_BOUNDING_ENTITIES, i, 0));
    }
    static const char *kPeriodic[] = {"entity", "affine", "node_pairs"};
    mxArray *periodic = mio_o_struct_array(np, 3, kPeriodic);
    for (int64_t i = 0; i < np; ++i) {
        for (int32_t f = 0; f < 3; ++f) {
            mxSetField(periodic, (mwIndex)i, kPeriodic[f],
                       ofi_gmsh_array(info, MIO_GMSH_PERIODIC, i, f));
        }
    }
    const char *names[2] = {"bounding_entities", "periodic"};
    mxArray *values[2] = {bounding, periodic};
    return mio_o_struct(2, names, values);
}

/* --- MFEM NURBS ---------------------------------------------------------- *
 *
 * The NURBS side channel of an MFEM NURBS mesh (v16.40.0), copied: per patch
 * its degrees, knot vectors, control net (points x space dimension), weights
 * and field control values; the field names; the patch topology; the cell map.
 * Vertices, patches and knot spans come back 1-based; attributes and the
 * edges' knot-vector ids as the file spelled them. */

static mxArray *ofi_nurbs_array(const mio_format_info *info, int32_t section, int64_t i,
                                int32_t field, double shift, int32_t *ndim_out) {
    const void *data = NULL;
    mio_dtype dtype = MIO_FLOAT64;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_o_check(mio_nurbs_info_array(info, section, i, field, &data, &dtype, &ndim, shape),
                "nurbs_info");
    if (ndim_out != NULL) *ndim_out = ndim;
    return ofi_array(data, dtype, ndim, shape, shift);
}

static mxArray *ofi_nurbs_int(const mio_format_info *info, int32_t section, int64_t i,
                              int32_t field) {
    int64_t v = 0;
    mio_o_check(mio_nurbs_info_int(info, section, i, field, &v), "nurbs_info");
    return mio_o_scalar((double)v);
}

static mxArray *O_mio_nurbs_info(const mxArray **a) {
    const mio_format_info *info = ofi_info_of(a[0], 0);
    const int64_t np = mio_nurbs_info_count(info, MIO_NURBS_PATCHES);
    const int64_t nf = mio_nurbs_info_count(info, MIO_NURBS_FIELDS);
    const int64_t nt = mio_nurbs_info_count(info, MIO_NURBS_TOPOLOGY);
    const int64_t nc = mio_nurbs_info_count(info, MIO_NURBS_CELL_MAP);
    if (np < 0 || nf < 0 || nt < 0 || nc < 0) mio_o_fail("nurbs_info");

    static const char *kPatch[] = {"degrees", "knots", "control_points", "weights", "fields"};
    mxArray *patches = mio_o_struct_array(np, 5, kPatch);
    for (int64_t p = 0; p < np; ++p) {
        int64_t dim = 0;
        mio_o_check(mio_nurbs_info_int(info, MIO_NURBS_PATCHES, p, 0, &dim), "nurbs_info");
        mxSetField(patches, (mwIndex)p, "degrees",
                   ofi_nurbs_array(info, MIO_NURBS_PATCHES, p, 0, 0.0, NULL));
        mxArray *knots = mio_o_cell(dim);
        for (int64_t d = 0; d < dim; ++d) {
            mxSetCell(knots, (mwIndex)d,
                      ofi_nurbs_array(info, MIO_NURBS_PATCHES, p, (int32_t)(3 + d), 0.0, NULL));
        }
        mxSetField(patches, (mwIndex)p, "knots", knots);
        mxSetField(patches, (mwIndex)p, "control_points",
                   ofi_nurbs_array(info, MIO_NURBS_PATCHES, p, 1, 0.0, NULL));
        mxSetField(patches, (mwIndex)p, "weights",
                   ofi_nurbs_array(info, MIO_NURBS_PATCHES, p, 2, 0.0, NULL));
        mxArray *fields = mio_o_cell(nf);
        for (int64_t k = 0; k < nf; ++k) {
            mxSetCell(fields, (mwIndex)k,
                      ofi_nurbs_array(info, MIO_NURBS_PATCHES, p, (int32_t)(6 + k), 0.0, NULL));
        }
        mxSetField(patches, (mwIndex)p, "fields", fields);
    }

    mxArray *field_names = mio_o_cellstr(nf);
    for (int64_t k = 0; k < nf; ++k) {
        const int64_t n = mio_nurbs_info_string(info, MIO_NURBS_FIELDS, k, 0, NULL, 0);
        if (n < 0) mio_o_fail("nurbs_info");
        char *buf = (char *)mxCalloc((size_t)n + 1, 1);
        mio_nurbs_info_string(info, MIO_NURBS_FIELDS, k, 0, buf, n + 1);
        mxSetCell(field_names, (mwIndex)k, mxCreateString(buf));
        mxFree(buf);
    }

    mxArray *topology;
    if (nt > 0) {
        mxArray *t[7];
        int32_t edges_ndim = 0;
        t[0] = ofi_nurbs_int(info, MIO_NURBS_TOPOLOGY, 0, 0);
        t[1] = ofi_nurbs_int(info, MIO_NURBS_TOPOLOGY, 0, 1);
        t[2] = ofi_nurbs_array(info, MIO_NURBS_TOPOLOGY, 0, 0, 1.0, NULL);
        t[3] = ofi_nurbs_array(info, MIO_NURBS_TOPOLOGY, 0, 1, 0.0, NULL);
        t[4] = ofi_nurbs_array(info, MIO_NURBS_TOPOLOGY, 0, 2, 1.0, NULL);
        t[5] = ofi_nurbs_array(info, MIO_NURBS_TOPOLOGY, 0, 3, 0.0, NULL);
        t[6] = ofi_nurbs_array(info, MIO_NURBS_TOPOLOGY, 0, 4, 0.0, &edges_ndim);
        if (edges_ndim == 2 && mxGetN(t[6]) == 3) {
            /* (knot vector, vertex, vertex) per row: the vertices 1-based, the
             * knot-vector id as the file spelled it. */
            double *e = mxGetPr(t[6]);
            const size_t rows = mxGetM(t[6]);
            for (size_t r = 0; r < rows; ++r) {
                e[rows + r] += 1.0;
                e[2 * rows + r] += 1.0;
            }
        }
        const char *names[7] = {"dimension",        "num_vertices",      "patch_vertices",
                                "patch_attributes", "boundary_vertices", "boundary_attributes",
                                "edges"};
        topology = mio_o_struct(7, names, t);
    } else {
        topology = mio_o_empty();
    }

    mxArray *cellmap[3];
    for (int32_t f = 0; f < 3; ++f) {
        cellmap[f] = nc > 0 ? ofi_nurbs_array(info, MIO_NURBS_CELL_MAP, 0, f, 1.0, NULL)
                            : mio_o_column(0);
    }
    const char *names[6] = {"patches",    "field_names", "topology",
                            "cell_patch", "cell_span",   "boundary_cell_patch"};
    mxArray *values[6] = {patches, field_names, topology, cellmap[0], cellmap[1], cellmap[2]};
    return mio_o_struct(6, names, values);
}

/* --- the isogeometric layer (v16.43.0) ----------------------------------- *
 *
 * Refine, sample and evaluate the NURBS model in a handle without going
 * through its sampling. */

static void ofi_iga_int3(const mxArray *x, const char *what, int32_t out[3]) {
    size_t n = 0;
    double *v = mio_o_doubles(x, &n, what);
    if (n < 1 || n > 3) {
        mio_o_error("iga", "'%s' must be one integer or one per direction (at most 3)", what);
    }
    for (size_t d = 0; d < 3; ++d) {
        const double e = v[d < n ? d : n - 1];
        if (isnan(e)) mio_o_error("iga", "'%s' must not contain NaN", what);
        out[d] = (int32_t)e;
    }
    mxFree(v);
}

static mxArray *O_mio_iga_refine(const mxArray **a) {
    const mio_format_info *info = ofi_info_of(a[0], 0);
    mio_iga_opts opts;
    mio_iga_opts_init(&opts);
    ofi_iga_int3(a[1], "subdivisions", opts.subdivisions);
    ofi_iga_int3(a[2], "elevate", opts.elevate);
    opts.sample_subdivisions = (int32_t)mio_o_int(a[3], "sample_subdivisions");
    opts.trim_samples = (int32_t)mio_o_int(a[4], "trim_samples");
    opts.edges = mio_o_bool(a[5], "edges") ? 1 : 0;
    opts.trim = mio_o_bool(a[6], "trim") ? 1 : 0;
    mio_format_info *refined = NULL;
    mio_mesh *mesh = mio_iga_refine(info, &opts, &refined);
    return ofi_mesh_and_info(mesh, refined, "iga_refine");
}

static mxArray *O_mio_iga_sample(const mxArray **a) {
    const mio_format_info *info = ofi_info_of(a[0], 0);
    mio_iga_opts opts;
    mio_iga_opts_init(&opts);
    for (int d = 0; d < 3; ++d) {
        opts.subdivisions[d] = 1;
        opts.elevate[d] = 0;
    }
    opts.sample_subdivisions = (int32_t)mio_o_int(a[1], "sample_subdivisions");
    opts.trim_samples = (int32_t)mio_o_int(a[2], "trim_samples");
    opts.edges = mio_o_bool(a[3], "edges") ? 1 : 0;
    opts.trim = mio_o_bool(a[4], "trim") ? 1 : 0;
    mio_format_info *mapped = NULL;
    mio_mesh *mesh = mio_iga_sample(info, &opts, &mapped);
    return ofi_mesh_and_info(mesh, mapped, "iga_sample");
}

/* One array of an evaluation as a double array with the C dims kept -- rows
 * are points, so (N, d, sdim) stays N x d x sdim (R reverses them to
 * (sdim, d, N)); [] when it is absent. */
static mxArray *ofi_iga_array(const mio_iga_evaluation *ev, int32_t which, int64_t index,
                              int as_logical) {
    const void *data = NULL;
    mio_dtype dtype = MIO_FLOAT64;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    if (mio_iga_evaluation_array(ev, which, index, &data, &dtype, &ndim, shape) != MIO_OK) {
        return mio_o_empty();
    }
    int64_t n = 1;
    for (int32_t k = 0; k < ndim; ++k) n *= shape[k];
    if (as_logical) {
        mxArray *out = mio_o_logical_column(n);
        if (n > 0) {
            double *tmp = (double *)mxMalloc((size_t)n * sizeof(double));
            mio_o_fill_real(tmp, data, dtype, (size_t)n, 0.0);
            mxLogical *dst = mxGetLogicals(out);
            for (int64_t i = 0; i < n; ++i) dst[i] = tmp[i] != 0.0;
            mxFree(tmp);
        }
        return out;
    }
    if (dtype == MIO_FLOAT64) return mio_o_typed_from_c(data, dtype, ndim, shape);
    double *tmp = (double *)mxMalloc((size_t)(n > 0 ? n : 1) * sizeof(double));
    if (n > 0) mio_o_fill_real(tmp, data, dtype, (size_t)n, 0.0);
    mxArray *out = mio_o_typed_from_c(tmp, MIO_FLOAT64, ndim, shape);
    mxFree(tmp);
    return out;
}

/* iga_evaluate(info, patch, params, derivatives, trim_samples) -- params is
 * (N x d), one row per point: exactly the row-major (N, d) the C API takes. */
static mxArray *O_mio_iga_evaluate(const mxArray **a) {
    const mio_format_info *info = ofi_info_of(a[0], 0);
    const int64_t patch = (int64_t)mio_o_int(a[1], "patch") - 1;
    if (mxGetNumberOfDimensions(a[2]) != 2) {
        mio_o_error("iga_evaluate", "'params' must be an (N x d) matrix, one row per point");
    }
    int64_t n = 0, d = 0;
    double *p = mio_o_rows_to_c(a[2], &n, &d, 0.0, "params");
    /* The C API reads N * (the patch's parametric dimension) doubles: check
     * the column count so a short matrix is an error, not an overread. A
     * B-rep face is 2-D; an MFEM patch states its dimension. */
    int64_t want = 2;
    int64_t pd = 0;
    if (patch >= 0 && mio_nurbs_info_int(info, MIO_NURBS_PATCHES, patch, 0, &pd) == MIO_OK &&
        pd > 0) {
        want = pd;
    }
    if (n > 0 && d != want) {
        mio_o_error("iga_evaluate",
                    "'params' must be an (N x %d) matrix, one row per point; got %d columns",
                    (int)want, (int)d);
    }
    mio_iga_evaluation *ev = mio_iga_evaluate(info, patch, p, n,
                                              (int32_t)mio_o_int(a[3], "derivatives"),
                                              (int32_t)mio_o_int(a[4], "trim_samples"));
    if (ev == NULL) mio_o_fail("iga_evaluate");
    const int64_t nf = mio_iga_evaluation_num_fields(ev);
    mxArray *v[5];
    v[0] = ofi_iga_array(ev, MIO_IGA_POINTS, 0, 0);
    v[1] = mio_o_cell(nf > 0 ? nf : 0);
    for (int64_t k = 0; k < nf; ++k) {
        mxSetCell(v[1], (mwIndex)k, ofi_iga_array(ev, MIO_IGA_FIELD, k, 0));
    }
    v[2] = ofi_iga_array(ev, MIO_IGA_DERIVATIVES, 0, 0);
    v[3] = ofi_iga_array(ev, MIO_IGA_NORMALS, 0, 0);
    v[4] = ofi_iga_array(ev, MIO_IGA_INSIDE, 0, 1);
    mio_iga_evaluation_free(ev);
    const char *names[5] = {"points", "fields", "derivatives", "normals", "inside"};
    return mio_o_struct(5, names, v);
}

/* --- mesh generation from CAD (v17.13.0; doc/mesh_cad.md) ----------------- */

static int64_t ofi_mesh_cad_warnings(void *ctx, char *buf, int64_t buflen) {
    (void)ctx;
    return mio_mesh_cad_surface_warnings(buf, buflen);
}

static int64_t ofi_mesh_cad_reason(void *ctx, char *buf, int64_t buflen) {
    return mio_mesh_cad_surface_face_reason(*(const int64_t *)ctx, buf, buflen);
}

/* A reason string of record `k` through `fn` (the C string rule). */
static mxArray *ofi_reason(int64_t (*fn)(int64_t, char *, int64_t), int64_t k) {
    const int64_t n = fn(k, NULL, 0);
    if (n <= 0) return mio_o_str("");
    char *buf = (char *)mxCalloc((size_t)n + 1, 1);
    fn(k, buf, n + 1);
    mxArray *out = mxCreateString(buf);
    mxFree(buf);
    return out;
}

/* The R binding's data frames are 1x1 structs of equal-length n x 1 columns
 * here (numbers double, flags logical, strings an n x 1 cellstr):
 * `r.corners.angle_deg` is the column, `numel(r.corners.face)` the row
 * count. */

/* The corners (v17.15.0) of the last call: 1-based face and point, the
 * parameters, the point, the angles, the reach and whether it is a pole. */
static mxArray *ofi_mesh_cad_corners(void) {
    const int64_t n = mio_mesh_cad_surface_num_corners();
    const int64_t rows = n > 0 ? n : 0;
    static const char *names[] = {"face", "point", "u", "v", "x", "y", "z",
                                  "angle_deg", "min_forced_deg", "reach", "pole"};
    enum { OFI_CORNER_N = 11 };
    mxArray *cols[OFI_CORNER_N];
    for (int c = 0; c < OFI_CORNER_N; ++c) {
        cols[c] = c == 10 ? mio_o_logical_column(rows) : mio_o_column(rows);
    }
    for (int64_t k = 0; k < rows; ++k) {
        mio_mesh_cad_corner corner;
        if (mio_mesh_cad_surface_corner(k, &corner) != MIO_OK) mio_o_fail("mesh_cad_surface");
        const double v[10] = {(double)corner.face + 1.0, (double)corner.point + 1.0, corner.u,
                              corner.v, corner.x[0], corner.x[1], corner.x[2], corner.angle_deg,
                              corner.min_forced_deg, corner.reach};
        for (int c = 0; c < 10; ++c) mxGetPr(cols[c])[k] = v[c];
        mxGetLogicals(cols[10])[k] = corner.pole != 0;
    }
    return mio_o_struct(OFI_CORNER_N, names, cols);
}

/* v17.20.0: the last call's reparametrised faces (6.4.4). */
static mxArray *ofi_mesh_cad_reparametrised(void) {
    const int64_t n = mio_mesh_cad_surface_num_reparametrised();
    const int64_t rows = n > 0 ? n : 0;
    mxArray *face = mio_o_column(rows), *cond = mio_o_column(rows), *disc = mio_o_column(rows);
    mxArray *method = mio_o_column(rows), *reason = mio_o_cellstr(rows);
    for (int64_t i = 0; i < rows; ++i) {
        int64_t f = -1;
        double c = 0.0, d = 0.0;
        int32_t m = 0;
        mio_mesh_cad_surface_reparametrised(i, &f, &c, &d, &m);
        mxGetPr(face)[i] = (double)f + 1.0;
        mxGetPr(cond)[i] = c;
        mxGetPr(disc)[i] = d;
        mxGetPr(method)[i] = (double)m;
        mxSetCell(reason, (mwIndex)i, ofi_reason(mio_mesh_cad_surface_reparam_reason, i));
    }
    const char *names[] = {"face", "condition_number", "condition_number_disc", "method",
                           "reason"};
    mxArray *cols[] = {face, cond, disc, method, reason};
    return mio_o_struct(5, names, cols);
}

/* v17.20.0: the slivers found (6.5.2); host 0 for none. */
static mxArray *ofi_mesh_cad_merged_faces(void) {
    const int64_t n = mio_mesh_cad_surface_num_merged_faces();
    const int64_t rows = n > 0 ? n : 0;
    mxArray *sliver = mio_o_column(rows), *host = mio_o_column(rows), *width = mio_o_column(rows);
    mxArray *merged = mio_o_logical_column(rows), *reason = mio_o_cellstr(rows);
    for (int64_t i = 0; i < rows; ++i) {
        int64_t s = -1, h = -1;
        double w = 0.0;
        int32_t m = 0;
        mio_mesh_cad_surface_merged_face(i, &s, &h, &w, &m);
        mxGetPr(sliver)[i] = (double)s + 1.0;
        mxGetPr(host)[i] = (double)h + 1.0;
        mxGetPr(width)[i] = w;
        mxGetLogicals(merged)[i] = m != 0;
        mxSetCell(reason, (mwIndex)i, ofi_reason(mio_mesh_cad_surface_merged_face_reason, i));
    }
    const char *names[] = {"sliver", "host", "width", "merged", "reason"};
    mxArray *cols[] = {sliver, host, width, merged, reason};
    return mio_o_struct(5, names, cols);
}

/* v17.20.0: where a declared tolerance accepted a gap (6.5.1). */
static mxArray *ofi_mesh_cad_tolerance_uses(void) {
    const int64_t n = mio_mesh_cad_surface_num_tolerance_uses();
    const int64_t rows = n > 0 ? n : 0;
    mxArray *edge = mio_o_column(rows), *declared = mio_o_column(rows), *gap = mio_o_column(rows);
    mxArray *where = mio_o_cellstr(rows);
    static const char *kWhere[] = {"closed_start", "corner", "projection", "pairing"};
    for (int64_t i = 0; i < rows; ++i) {
        int64_t e = -1;
        double d = 0.0, g = 0.0;
        int32_t w = 0;
        mio_mesh_cad_surface_tolerance_use(i, &e, &d, &g, &w);
        mxGetPr(edge)[i] = (double)e + 1.0;
        mxGetPr(declared)[i] = d;
        mxGetPr(gap)[i] = g;
        mxSetCell(where, (mwIndex)i, mio_o_str(kWhere[w >= 0 && w < 4 ? w : 2]));
    }
    const char *names[] = {"edge", "declared", "gap", "where"};
    mxArray *cols[] = {edge, declared, gap, where};
    return mio_o_struct(4, names, cols);
}

/* v17.20.0: the corners lopped (6.2.3), 1-based face and point. */
static mxArray *ofi_mesh_cad_lopped_corners(void) {
    const int64_t n = mio_mesh_cad_surface_num_lopped_corners();
    const int64_t rows = n > 0 ? n : 0;
    mxArray *face = mio_o_column(rows), *point = mio_o_column(rows), *angle = mio_o_column(rows);
    mxArray *lfs0 = mio_o_column(rows), *radius = mio_o_column(rows);
    mxArray *turning = mio_o_column(rows), *bound = mio_o_column(rows);
    mxArray *forced = mio_o_column(rows), *min_forced = mio_o_column(rows);
    mxArray *within = mio_o_logical_column(rows), *guaranteed = mio_o_logical_column(rows);
    mxArray *reason = mio_o_cellstr(rows);
    for (int64_t i = 0; i < rows; ++i) {
        int64_t f = -1, p = -1, nf = 0;
        double an = 0.0, l = 0.0, r = 0.0, t = 0.0, b = 0.0, m = 0.0;
        int32_t w = 0, g = 0;
        mio_mesh_cad_surface_lopped_corner(i, &f, &p, &an, &l, &r, &t, &b, &nf, &m, &w, &g);
        mxGetPr(face)[i] = (double)f + 1.0;
        mxGetPr(point)[i] = (double)p + 1.0;
        mxGetPr(angle)[i] = an;
        mxGetPr(lfs0)[i] = l;
        mxGetPr(radius)[i] = r;
        mxGetPr(turning)[i] = t;
        mxGetPr(bound)[i] = b;
        mxGetPr(forced)[i] = (double)nf;
        mxGetPr(min_forced)[i] = m;
        mxGetLogicals(within)[i] = w != 0;
        mxGetLogicals(guaranteed)[i] = g != 0;
        mxSetCell(reason, (mwIndex)i, ofi_reason(mio_mesh_cad_surface_lopped_corner_reason, i));
    }
    const char *names[] = {"face",       "point",          "angle_deg",        "lfs0",
                           "lop_radius", "turning_deg",    "stated_bound_deg", "num_forced",
                           "min_forced_deg", "within_two_radii", "guaranteed", "reason"};
    mxArray *cols[] = {face,  point,  angle,      lfs0,   radius,     turning,
                       bound, forced, min_forced, within, guaranteed, reason};
    return mio_o_struct(12, names, cols);
}

/* A list of 1-based index vectors from a (k, buf, n) getter. */
static mxArray *ofi_index_lists(int64_t count, int64_t (*fn)(int64_t, int64_t *, int64_t)) {
    mxArray *out = mio_o_cell(count > 0 ? count : 0);
    for (int64_t k = 0; k < count; ++k) {
        const int64_t n = fn(k, NULL, 0);
        mxArray *v = mio_o_column(n > 0 ? n : 0);
        if (n > 0) {
            int64_t *buf = (int64_t *)mxMalloc((size_t)n * sizeof(int64_t));
            fn(k, buf, n);
            for (int64_t i = 0; i < n; ++i) mxGetPr(v)[i] = (double)buf[i] + 1.0;
            mxFree(buf);
        }
        mxSetCell(out, (mwIndex)k, v);
    }
    return out;
}

/* project_to_cad (roadmap 6.10.1): every point of a mesh onto the CAD face or
 * edge its cells name. Returns struct(mesh, num_projected, ..., outside_trim).
 * project_to_cad(mesh, info, face_array, edge_array, max_distance,
 *                trim_samples) */
static mxArray *O_mio_project_to_cad(const mxArray **a) {
    const mio_format_info *info = ofi_info_of(a[1], 0);
    mio_project_to_cad_opts opts;
    mio_project_to_cad_opts_init(&opts);
    opts.face_array = mio_o_string(a[2], "face_array");
    opts.edge_array = mio_o_string(a[3], "edge_array");
    opts.max_distance = mio_o_double(a[4], "max_distance");
    opts.trim_samples = (int32_t)mio_o_int(a[5], "trim_samples");
    mio_project_to_cad_report report;
    memset(&report, 0, sizeof(report));
    mio_mesh *out = mio_project_to_cad(mio_o_mesh(a[0]), info, &opts, &report);
    if (out == NULL) mio_o_fail("project_to_cad");
    mxArray *mo = mio_o_wrap_mesh(out);
    const int64_t n = mio_project_to_cad_outside_trim(NULL, 0);
    mxArray *outside = mio_o_column(n > 0 ? n : 0);
    if (n > 0) {
        int64_t *buf = (int64_t *)mxMalloc((size_t)n * sizeof(int64_t));
        mio_project_to_cad_outside_trim(buf, n);
        for (int64_t i = 0; i < n; ++i) mxGetPr(outside)[i] = (double)buf[i] + 1.0;
        mxFree(buf);
    }
    const char *names[] = {"mesh",          "num_projected", "num_on_edges", "num_on_curves",
                           "num_corners",   "num_ambiguous", "num_too_far",  "num_failed",
                           "num_outside_trim", "max_move",   "outside_trim"};
    mxArray *values[] = {mo,
                         mio_o_scalar((double)report.num_projected),
                         mio_o_scalar((double)report.num_on_edges),
                         mio_o_scalar((double)report.num_on_curves),
                         mio_o_scalar((double)report.num_corners),
                         mio_o_scalar((double)report.num_ambiguous),
                         mio_o_scalar((double)report.num_too_far),
                         mio_o_scalar((double)report.num_failed),
                         mio_o_scalar((double)report.num_outside_trim),
                         mio_o_scalar(report.max_move),
                         outside};
    return mio_o_struct(11, names, values);
}

/* mesh_cad_surface(info, target_size, max_chordal_deviation, min_points_per_edge,
 *   tolerance, record_face_ids, edges, match_unreferenced_trims, strict,
 *   max_points, min_angle, smooth_iterations, smooth_method, order, grading,
 *   curvature_deviation, cells_across_gaps, pole_rings, face_size_faces,
 *   face_sizes, edge_size_edges, edge_sizes, background, background_array,
 *   metric, virtual_edges, reparam_condition_number, model_tolerances,
 *   virtual_faces, lop_corners, metric_background, metric_array, metric_error,
 *   untangle, quads, quad_pairs) -- 36 arguments, the enums already 0-based. */
static mxArray *O_mio_mesh_cad_surface(const mxArray **a) {
    const mio_format_info *info = ofi_info_of(a[0], 0);
    mio_mesh_cad_surface_opts opts;
    mio_mesh_cad_surface_opts_init(&opts);
    opts.target_size = mio_o_double(a[1], "target_size");
    opts.max_chordal_deviation = mio_o_double(a[2], "max_chordal_deviation");
    opts.min_points_per_edge = mio_o_int64(a[3], "min_points_per_edge");
    opts.tolerance = mio_o_double(a[4], "tolerance");
    opts.record_face_ids = mio_o_bool(a[5], "record_face_ids") ? 1 : 0;
    opts.edges = mio_o_bool(a[6], "edges") ? 1 : 0;
    opts.match_unreferenced_trims = mio_o_bool(a[7], "match_unreferenced_trims") ? 1 : 0;
    opts.strict = mio_o_bool(a[8], "strict") ? 1 : 0;
    opts.max_points = mio_o_int64(a[9], "max_points");
    opts.min_angle_deg = mio_o_double(a[10], "min_angle");
    opts.smooth_iterations = (int32_t)mio_o_int64(a[11], "smooth_iterations");
    opts.smooth_method = (int32_t)mio_o_int64(a[12], "smooth_method");
    opts.order = (int32_t)mio_o_int64(a[13], "order");
    opts.grading = mio_o_double(a[14], "grading");
    opts.curvature_deviation = mio_o_double(a[15], "curvature_deviation");
    opts.cells_across_gaps = mio_o_double(a[16], "cells_across_gaps");
    opts.pole_rings = (int32_t)mio_o_int64(a[17], "pole_rings");
    /* v17.19.0: given sizes (1-based indices from Octave), the background, the
     * metric and virtual edges. */
    mio_mesh_cad_surface_opts_ex ex;
    mio_mesh_cad_surface_opts_ex_init(&ex);
    ex.metric = (int32_t)mio_o_int64(a[24], "metric");
    ex.virtual_edges = mio_o_bool(a[25], "virtual_edges") ? 1 : 0;
    /* v17.20.0: the reparametrisation, model tolerances, virtual faces and
     * corner lopping. */
    ex.reparam_condition_number = mio_o_double(a[26], "reparam_condition_number");
    ex.model_tolerances = mio_o_bool(a[27], "model_tolerances") ? 1 : 0;
    ex.virtual_faces = mio_o_bool(a[28], "virtual_faces") ? 1 : 0;
    ex.lop_corners = mio_o_bool(a[29], "lop_corners") ? 1 : 0;
    size_t nfs = 0, nfsz = 0, nes = 0, nesz = 0;
    /* Numeric vectors: the .m wrapper splits the (index, size) pairs. */
    double *fsf = mio_o_doubles(a[18], &nfs, "face_size_faces");
    double *fsz = mio_o_doubles(a[19], &nfsz, "face_sizes");
    double *esf = mio_o_doubles(a[20], &nes, "edge_size_edges");
    double *esz = mio_o_doubles(a[21], &nesz, "edge_sizes");
    if (nfsz != nfs || nesz != nes) {
        mio_o_error("mesh_cad_surface",
                    "mesh_cad_surface: the size indices and sizes must have the same length");
    }
    int64_t *fidx = (int64_t *)mxMalloc((nfs > 0 ? nfs : 1) * sizeof(int64_t));
    int64_t *eidx = (int64_t *)mxMalloc((nes > 0 ? nes : 1) * sizeof(int64_t));
    for (size_t i = 0; i < nfs; ++i) fidx[i] = (int64_t)fsf[i] - 1;
    for (size_t i = 0; i < nes; ++i) eidx[i] = (int64_t)esf[i] - 1;
    ex.face_size_faces = fidx;
    ex.face_sizes = fsz;
    ex.num_face_sizes = (int64_t)nfs;
    ex.edge_size_edges = eidx;
    ex.edge_sizes = esz;
    ex.num_edge_sizes = (int64_t)nes;
    if (!mio_o_is_null(a[22])) ex.background = mio_o_mesh(a[22]);
    ex.background_array = mio_o_string(a[23], "background_array");
    /* Roadmap 6.10.2: a metric from a field's Hessian. */
    if (!mio_o_is_null(a[30])) ex.metric_background = mio_o_mesh(a[30]);
    ex.metric_array = mio_o_string(a[31], "metric_array");
    ex.metric_error = mio_o_double(a[32], "metric_error");
    /* Roadmap 6.6.3 and 6.8: untangling and quads (a mio_mesh_cad_quads). */
    ex.untangle = mio_o_bool(a[33], "untangle") ? 1 : 0;
    ex.quads = (int32_t)mio_o_int64(a[34], "quads");
    ex.quad_pairs = (int32_t)mio_o_int64(a[35], "quad_pairs");
    mio_mesh_cad_surface_report report;
    memset(&report, 0, sizeof(report));
    mio_mesh *out = mio_mesh_cad_surface2(info, &opts, &ex, &report);
    if (out == NULL) mio_o_fail("mesh_cad_surface");

    /* Wrap the mesh first: once it is registered nothing below can leak it. */
    mxArray *mo = mio_o_wrap_mesh(out);
    mxArray *warnings = mio_o_getstring(ofi_mesh_cad_warnings, NULL, "mesh_cad_surface");
    static const char *statuses[] = {"ok",           "empty",           "open_loop",
                                     "self_intersecting", "degenerate", "projection_failed",
                                     "budget_exceeded", "unsupported",  "other"};
    const int64_t nf = report.num_faces > 0 ? report.num_faces : 0;
    mxArray *status = mio_o_cellstr(nf);
    mxArray *reason = mio_o_cellstr(nf);
    for (int64_t f = 0; f < nf; ++f) {
        const int32_t st = mio_mesh_cad_surface_face_status(f);
        mxSetCell(status, (mwIndex)f, mio_o_str(st >= 0 && st <= 8 ? statuses[st] : "other"));
        int64_t face = f;
        mxSetCell(reason, (mwIndex)f, mio_o_getstring(ofi_mesh_cad_reason, &face, "mesh_cad_surface"));
    }

    mxArray *histogram = mio_o_column(MIO_MESH_CAD_ANGLE_BINS);
    {
        int64_t bins[MIO_MESH_CAD_ANGLE_BINS] = {0};
        mio_mesh_cad_surface_angle_histogram(bins, MIO_MESH_CAD_ANGLE_BINS);
        for (int b = 0; b < MIO_MESH_CAD_ANGLE_BINS; ++b) mxGetPr(histogram)[b] = (double)bins[b];
    }
    mxArray *corners = ofi_mesh_cad_corners();
    /* v17.19.0: per face, the 1-based edges that limited its given size; and
     * each virtual edge's chain of 1-based edges. */
    mxArray *limits = ofi_index_lists(nf, mio_mesh_cad_surface_face_size_limits);
    mxArray *virtual_edges_out =
        ofi_index_lists(mio_mesh_cad_surface_num_virtual_edges(), mio_mesh_cad_surface_virtual_edge);

    /* v17.20.0: the faces over the reparametrisation threshold, the slivers
     * found and the declared tolerances used, as named columns (1-based). */
    mxArray *reparam = ofi_mesh_cad_reparametrised();
    mxArray *merged = ofi_mesh_cad_merged_faces();
    mxArray *tol_uses = ofi_mesh_cad_tolerance_uses();
    mxArray *lopped = ofi_mesh_cad_lopped_corners();

    int64_t num_quads = 0, num_leftover = 0;
    double quad_sj = 0.0;
    mio_mesh_cad_surface_quad_stats(&num_quads, &num_leftover, &quad_sj);

    enum { OFI_MESH_CAD_N = 38 };
    const char *names[OFI_MESH_CAD_N];
    mxArray *values[OFI_MESH_CAD_N];
    int k = 0;
#define OFI_SET(name, value)                                                                       \
    do {                                                                                           \
        names[k] = (name);                                                                         \
        values[k] = (value);                                                                       \
        ++k;                                                                                       \
    } while (0)
#define OFI_NUM(f) OFI_SET(#f, mio_o_scalar((double)report.f))
    OFI_SET("mesh", mo);
    OFI_NUM(num_points);
    OFI_NUM(num_triangles);
    OFI_NUM(num_edge_segments);
    OFI_NUM(num_faces);
    OFI_NUM(num_failed);
    OFI_NUM(num_edges);
    OFI_NUM(num_open_edges);
    OFI_NUM(num_non_manifold_edges);
    OFI_NUM(num_misoriented_edges);
    OFI_NUM(num_degenerate_dropped);
    OFI_NUM(tolerance);
    OFI_NUM(min_angle_deg);
    OFI_NUM(max_angle_deg);
    OFI_NUM(min_edge);
    OFI_NUM(max_edge);
    OFI_NUM(max_edge_gap);
    OFI_SET("face_status", status);
    OFI_SET("face_reason", reason);
    /* One string, newline-separated; the .m wrapper splits it. */
    OFI_SET("warnings", warnings);
    OFI_NUM(num_refinement_points);
    OFI_NUM(num_faces_below_bound);
    OFI_NUM(num_smoothed);
    OFI_NUM(min_input_angle_deg);
    OFI_SET("angle_histogram", histogram);
    /* Named columns (R makes a data frame). */
    OFI_SET("corners", corners);
    OFI_NUM(num_curved_straightened);
    OFI_SET("face_size_limits", limits);
    OFI_SET("virtual_edges", virtual_edges_out);
    OFI_SET("reparametrised", reparam);
    OFI_SET("merged_faces", merged);
    OFI_SET("tolerance_uses", tol_uses);
    OFI_SET("lopped_corners", lopped);
    OFI_SET("num_untangled", mio_o_scalar((double)mio_mesh_cad_surface_num_untangled()));
    OFI_SET("num_quads", mio_o_scalar((double)num_quads));
    OFI_SET("num_leftover_triangles", mio_o_scalar((double)num_leftover));
    OFI_SET("min_quad_scaled_jacobian", mio_o_scalar(quad_sj));
    OFI_SET("min_metric_angle", mio_o_scalar(mio_mesh_cad_surface_min_metric_angle()));
#undef OFI_NUM
#undef OFI_SET
    if (k != OFI_MESH_CAD_N) {
        mio_o_error("mesh_cad_surface", "mesh_cad_surface: report has %d fields, expected %d", k,
                    (int)OFI_MESH_CAD_N);
    }
    return mio_o_struct(OFI_MESH_CAD_N, names, values);
}

/* --- B-rep --------------------------------------------------------------- *
 *
 * The B-rep side channel of a Kratos CAD JSON (or OpenCASCADE STEP/IGES)
 * handle, copied: patches, faces with their loops and trims, edges and
 * vertices with their topology, the brep groups and the file's other members.
 * Nested lists come from the C API's flat sections (each owner names its first
 * item and count) and become struct arrays here. Indices (a face's patch, a
 * brep's lists) come back 1-based; ids, trim indices and flags as the file
 * spelled them. */

typedef struct {
    const mio_format_info *info;
    int32_t section;
    int64_t index;
    int32_t field;
} ofi_brep_ctx;

static int64_t ofi_brep_get_string(void *ctx, char *buf, int64_t buflen) {
    ofi_brep_ctx *c = (ofi_brep_ctx *)ctx;
    return mio_brep_info_string(c->info, c->section, c->index, c->field, buf, buflen);
}

static mxArray *ofi_brep_str(const mio_format_info *info, int32_t section, int64_t i,
                             int32_t field) {
    ofi_brep_ctx c = {info, section, i, field};
    return mio_o_getstring(ofi_brep_get_string, &c, "brep_info");
}

static int64_t ofi_brep_ival(const mio_format_info *info, int32_t section, int64_t i,
                             int32_t field) {
    int64_t v = 0;
    mio_o_check(mio_brep_info_int(info, section, i, field, &v), "brep_info");
    return v;
}

static mxArray *ofi_brep_num(const mio_format_info *info, int32_t section, int64_t i,
                             int32_t field) {
    return mio_o_scalar((double)ofi_brep_ival(info, section, i, field));
}

static mxArray *ofi_brep_flag(const mio_format_info *info, int32_t section, int64_t i,
                              int32_t field) {
    return mio_o_logical(ofi_brep_ival(info, section, i, field) != 0);
}

static mxArray *ofi_brep_array(const mio_format_info *info, int32_t section, int64_t i,
                               int32_t field, double shift) {
    const void *data = NULL;
    mio_dtype dtype = MIO_FLOAT64;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {0};
    mio_o_check(mio_brep_info_array(info, section, i, field, &data, &dtype, &ndim, shape),
                "brep_info");
    return ofi_array(data, dtype, ndim, shape, shift);
}

static mxArray *ofi_brep_curve(const mio_format_info *info, int32_t section, int64_t i) {
    mxArray *f[9];
    f[0] = ofi_brep_num(info, section, i, 10);
    f[1] = ofi_brep_array(info, section, i, 12, 0.0);
    f[2] = ofi_brep_flag(info, section, i, 12);
    f[3] = ofi_brep_flag(info, section, i, 11);
    f[4] = ofi_brep_array(info, section, i, 10, 0.0);
    f[5] = ofi_brep_array(info, section, i, 11, 0.0);
    f[6] = ofi_brep_array(info, section, i, 13, 0.0);
    f[7] = ofi_brep_array(info, section, i, 14, 0.0);
    f[8] = ofi_brep_str(info, section, i, 10);
    const char *names[9] = {"degree",  "knots", "reduced_knots", "rational", "control_points",
                            "weights", "ids",   "active_range",  "extra"};
    return mio_o_struct(9, names, f);
}

static const char *kOfiTrimFields[] = {"trim_index", "direction", "curve", "extra"};

static void ofi_brep_trim_into(mxArray *arr, mwIndex at, const mio_format_info *info, int64_t t) {
    mxSetField(arr, at, "trim_index", ofi_brep_num(info, MIO_BREP_TRIMS, t, 2));
    mxSetField(arr, at, "direction", ofi_brep_flag(info, MIO_BREP_TRIMS, t, 3));
    mxSetField(arr, at, "curve", ofi_brep_curve(info, MIO_BREP_TRIMS, t));
    mxSetField(arr, at, "extra", ofi_brep_str(info, MIO_BREP_TRIMS, t, 0));
}

static mxArray *ofi_brep_refs(const mio_format_info *info, int64_t first, int64_t count) {
    static const char *names[6] = {"id",                 "name",              "trim_index",
                                   "relative_direction", "local_coordinates", "extra"};
    mxArray *out = mio_o_struct_array(count, 6, names);
    for (int64_t k = 0; k < count; ++k) {
        const int64_t r = first + k;
        mxSetField(out, (mwIndex)k, "id", ofi_brep_num(info, MIO_BREP_REFS, r, 0));
        mxSetField(out, (mwIndex)k, "name", ofi_brep_str(info, MIO_BREP_REFS, r, 0));
        mxSetField(out, (mwIndex)k, "trim_index", ofi_brep_num(info, MIO_BREP_REFS, r, 1));
        mxSetField(out, (mwIndex)k, "relative_direction", ofi_brep_num(info, MIO_BREP_REFS, r, 2));
        mxSetField(out, (mwIndex)k, "local_coordinates",
                   ofi_brep_array(info, MIO_BREP_REFS, r, 0, 0.0));
        mxSetField(out, (mwIndex)k, "extra", ofi_brep_str(info, MIO_BREP_REFS, r, 1));
    }
    return out;
}

static int64_t ofi_brep_count(const mio_format_info *info, int32_t section) {
    const int64_t n = mio_brep_info_count(info, section);
    if (n < 0) mio_o_fail("brep_info");
    return n;
}

static mxArray *O_mio_brep_info(const mxArray **a) {
    const mio_format_info *info = ofi_info_of(a[0], 0);

    static const char *kPatch[] = {"degrees", "knots", "control_points", "weights"};
    const int64_t np = ofi_brep_count(info, MIO_BREP_PATCHES);
    mxArray *patches = mio_o_struct_array(np, 4, kPatch);
    for (int64_t p = 0; p < np; ++p) {
        mxArray *deg = mio_o_column(2);
        mxGetPr(deg)[0] = (double)ofi_brep_ival(info, MIO_BREP_PATCHES, p, 0);
        mxGetPr(deg)[1] = (double)ofi_brep_ival(info, MIO_BREP_PATCHES, p, 1);
        mxArray *knots = mio_o_cell(2);
        mxSetCell(knots, 0, ofi_brep_array(info, MIO_BREP_PATCHES, p, 2, 0.0));
        mxSetCell(knots, 1, ofi_brep_array(info, MIO_BREP_PATCHES, p, 3, 0.0));
        mxSetField(patches, (mwIndex)p, "degrees", deg);
        mxSetField(patches, (mwIndex)p, "knots", knots);
        mxSetField(patches, (mwIndex)p, "control_points",
                   ofi_brep_array(info, MIO_BREP_PATCHES, p, 0, 0.0));
        mxSetField(patches, (mwIndex)p, "weights",
                   ofi_brep_array(info, MIO_BREP_PATCHES, p, 1, 0.0));
    }

    static const char *kBrep[] = {"id", "name", "faces", "edges", "vertices", "extra"};
    const int64_t nb = ofi_brep_count(info, MIO_BREP_BREPS);
    mxArray *breps = mio_o_struct_array(nb, 6, kBrep);
    for (int64_t b = 0; b < nb; ++b) {
        mxSetField(breps, (mwIndex)b, "id", ofi_brep_num(info, MIO_BREP_BREPS, b, 0));
        mxSetField(breps, (mwIndex)b, "name", ofi_brep_str(info, MIO_BREP_BREPS, b, 0));
        mxSetField(breps, (mwIndex)b, "faces", ofi_brep_array(info, MIO_BREP_BREPS, b, 0, 1.0));
        mxSetField(breps, (mwIndex)b, "edges", ofi_brep_array(info, MIO_BREP_BREPS, b, 1, 1.0));
        mxSetField(breps, (mwIndex)b, "vertices",
                   ofi_brep_array(info, MIO_BREP_BREPS, b, 2, 1.0));
        mxSetField(breps, (mwIndex)b, "extra", ofi_brep_str(info, MIO_BREP_BREPS, b, 1));
    }

    static const char *kFace[] = {"id",             "name",          "patch",
                                  "swapped_normal", "trimmed",       "rational",
                                  "reduced_knots",  "control_point_ids", "loops",
                                  "embedded_edges", "surface_extra", "extra"};
    static const char *kLoop[] = {"type", "trims", "extra"};
    const int64_t nf = ofi_brep_count(info, MIO_BREP_FACES);
    mxArray *faces = mio_o_struct_array(nf, 12, kFace);
    for (int64_t fi = 0; fi < nf; ++fi) {
        const mwIndex at = (mwIndex)fi;
        const int64_t first_loop = ofi_brep_ival(info, MIO_BREP_FACES, fi, 6);
        const int64_t loops = ofi_brep_ival(info, MIO_BREP_FACES, fi, 7);
        const int64_t first_trim = ofi_brep_ival(info, MIO_BREP_FACES, fi, 8);
        const int64_t trims = ofi_brep_ival(info, MIO_BREP_FACES, fi, 9);
        mxSetField(faces, at, "id", ofi_brep_num(info, MIO_BREP_FACES, fi, 0));
        mxSetField(faces, at, "name", ofi_brep_str(info, MIO_BREP_FACES, fi, 0));
        mxSetField(faces, at, "patch",
                   mio_o_scalar((double)ofi_brep_ival(info, MIO_BREP_FACES, fi, 1) + 1.0));
        mxSetField(faces, at, "swapped_normal", ofi_brep_flag(info, MIO_BREP_FACES, fi, 2));
        mxSetField(faces, at, "trimmed", ofi_brep_flag(info, MIO_BREP_FACES, fi, 3));
        mxSetField(faces, at, "rational", ofi_brep_flag(info, MIO_BREP_FACES, fi, 4));
        mxSetField(faces, at, "reduced_knots", ofi_brep_flag(info, MIO_BREP_FACES, fi, 5));
        mxSetField(faces, at, "control_point_ids",
                   ofi_brep_array(info, MIO_BREP_FACES, fi, 0, 0.0));
        mxArray *loop_arr = mio_o_struct_array(loops, 3, kLoop);
        for (int64_t l = 0; l < loops; ++l) {
            const int64_t li = first_loop + l;
            const int64_t t0 = ofi_brep_ival(info, MIO_BREP_LOOPS, li, 1);
            const int64_t tn = ofi_brep_ival(info, MIO_BREP_LOOPS, li, 2);
            mxArray *trim_arr = mio_o_struct_array(tn, 4, kOfiTrimFields);
            for (int64_t t = 0; t < tn; ++t) ofi_brep_trim_into(trim_arr, (mwIndex)t, info, t0 + t);
            mxSetField(loop_arr, (mwIndex)l, "type", ofi_brep_str(info, MIO_BREP_LOOPS, li, 0));
            mxSetField(loop_arr, (mwIndex)l, "trims", trim_arr);
            mxSetField(loop_arr, (mwIndex)l, "extra", ofi_brep_str(info, MIO_BREP_LOOPS, li, 1));
        }
        mxSetField(faces, at, "loops", loop_arr);
        int64_t embedded = 0;
        for (int64_t t = 0; t < trims; ++t) {
            if (ofi_brep_ival(info, MIO_BREP_TRIMS, first_trim + t, 1) < 0) ++embedded;
        }
        mxArray *emb = mio_o_struct_array(embedded, 4, kOfiTrimFields);
        for (int64_t t = 0, k = 0; t < trims; ++t) {
            if (ofi_brep_ival(info, MIO_BREP_TRIMS, first_trim + t, 1) < 0) {
                ofi_brep_trim_into(emb, (mwIndex)k++, info, first_trim + t);
            }
        }
        mxSetField(faces, at, "embedded_edges", emb);
        mxSetField(faces, at, "surface_extra", ofi_brep_str(info, MIO_BREP_FACES, fi, 2));
        mxSetField(faces, at, "extra", ofi_brep_str(info, MIO_BREP_FACES, fi, 1));
    }

    static const char *kEdge[] = {"id", "name", "curve", "topology", "extra"};
    const int64_t ne = ofi_brep_count(info, MIO_BREP_EDGES);
    mxArray *edges = mio_o_struct_array(ne, 5, kEdge);
    for (int64_t e = 0; e < ne; ++e) {
        mxSetField(edges, (mwIndex)e, "id", ofi_brep_num(info, MIO_BREP_EDGES, e, 0));
        mxSetField(edges, (mwIndex)e, "name", ofi_brep_str(info, MIO_BREP_EDGES, e, 0));
        mxSetField(edges, (mwIndex)e, "curve",
                   ofi_brep_ival(info, MIO_BREP_EDGES, e, 1) ? ofi_brep_curve(info, MIO_BREP_EDGES, e)
                                                             : mio_o_empty());
        mxSetField(edges, (mwIndex)e, "topology",
                   ofi_brep_refs(info, ofi_brep_ival(info, MIO_BREP_EDGES, e, 2),
                                 ofi_brep_ival(info, MIO_BREP_EDGES, e, 3)));
        mxSetField(edges, (mwIndex)e, "extra", ofi_brep_str(info, MIO_BREP_EDGES, e, 1));
    }

    static const char *kVertex[] = {"id", "name", "topology", "extra"};
    const int64_t nv = ofi_brep_count(info, MIO_BREP_VERTICES);
    mxArray *vertices = mio_o_struct_array(nv, 4, kVertex);
    for (int64_t k = 0; k < nv; ++k) {
        mxSetField(vertices, (mwIndex)k, "id", ofi_brep_num(info, MIO_BREP_VERTICES, k, 0));
        mxSetField(vertices, (mwIndex)k, "name", ofi_brep_str(info, MIO_BREP_VERTICES, k, 0));
        mxSetField(vertices, (mwIndex)k, "topology",
                   ofi_brep_refs(info, ofi_brep_ival(info, MIO_BREP_VERTICES, k, 1),
                                 ofi_brep_ival(info, MIO_BREP_VERTICES, k, 2)));
        mxSetField(vertices, (mwIndex)k, "extra", ofi_brep_str(info, MIO_BREP_VERTICES, k, 1));
    }

    const char *names[6] = {"patches", "breps", "faces", "edges", "vertices", "extra"};
    mxArray *values[6] = {patches, breps, faces, edges, vertices,
                          ofi_brep_str(info, MIO_BREP_FILE, 0, 0)};
    return mio_o_struct(6, names, values);
}

/* --- table --------------------------------------------------------------- */

static const mio_o_entry kEntries[] = {
    MIO_O_ENTRY(mio_read_with_info, 7),
    MIO_O_ENTRY(mio_write_with_info, 4),
    MIO_O_ENTRY(mio_format_info_release, 1),
    MIO_O_ENTRY(mio_format_info_is_open, 1),
    MIO_O_ENTRY(mio_format_info_format, 1),
    MIO_O_ENTRY(mio_gmsh_info, 1),
    MIO_O_ENTRY(mio_nurbs_info, 1),
    MIO_O_ENTRY(mio_iga_refine, 7),
    MIO_O_ENTRY(mio_iga_sample, 5),
    MIO_O_ENTRY(mio_project_to_cad, 6),
    MIO_O_ENTRY(mio_mesh_cad_surface, 36),
    MIO_O_ENTRY(mio_iga_evaluate, 5),
    MIO_O_ENTRY(mio_brep_info, 1),
    MIO_O_ENTRY(mio_mdpa_info, 1),
};

const mio_o_entry *mio_o_table_format_info(size_t *n) {
    *n = sizeof(kEntries) / sizeof(kEntries[0]);
    return kEntries;
}
