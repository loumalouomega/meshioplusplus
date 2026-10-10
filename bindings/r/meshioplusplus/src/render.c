/* Software rendering (v16.33.0, field rendering v16.34.0): a deterministic
 * software rasterizer and its terminal, HTML and PNG encodings. See
 * doc/tui.md.
 *
 * Every entry point takes the options as one named R list, so adding a field
 * to mio_render_opts is a row in a table here, not a change of signature. An
 * unknown name is an error naming the ones that exist: a typo never silently
 * does nothing. */

#include "mio_r.h"

#include <stddef.h>
#include <string.h>

typedef enum { F_STR, F_DBL, F_INT, F_BOOL, F_COLOR } field_kind;

typedef struct {
    const char *name;
    field_kind kind;
    size_t offset;
} field_def;

#define FIELD(n, k) {#n, k, offsetof(mio_render_opts, n)}

static const field_def render_fields[] = {
    FIELD(view, F_STR),           FIELD(color_by, F_STR),      FIELD(cmap, F_STR),
    FIELD(reduce, F_STR),         FIELD(expr, F_STR),          FIELD(vectors, F_STR),
    FIELD(warp, F_STR),           FIELD(quality_metric, F_STR), FIELD(streamlines, F_STR),
    FIELD(pixel_aspect, F_DBL),   FIELD(azimuth, F_DBL),       FIELD(elevation, F_DBL),
    FIELD(roll, F_DBL),           FIELD(fov_deg, F_DBL),       FIELD(zoom, F_DBL),
    FIELD(ambient, F_DBL),        FIELD(split_angle, F_DBL),   FIELD(feature_angle, F_DBL),
    FIELD(point_radius, F_DBL),   FIELD(scale_threshold, F_DBL), FIELD(vector_length, F_DBL),
    FIELD(warp_scale, F_DBL),     FIELD(stream_length, F_DBL),
    FIELD(width, F_INT),          FIELD(height, F_INT),        FIELD(supersample, F_INT),
    FIELD(component, F_INT),      FIELD(isolines, F_INT),      FIELD(vector_count, F_INT),
    FIELD(stream_seeds, F_INT),
    FIELD(perspective, F_BOOL),   FIELD(two_sided, F_BOOL),    FIELD(colorbar, F_BOOL),
    FIELD(axes, F_BOOL),          FIELD(scale_bar, F_BOOL),    FIELD(symmetric, F_BOOL),
    FIELD(categorical, F_BOOL),   FIELD(color_regions, F_BOOL), FIELD(category_edges, F_BOOL),
    FIELD(warp_outline, F_BOOL),
    FIELD(edge_color, F_COLOR),   FIELD(fill_color, F_COLOR), FIELD(line_color, F_COLOR),
    FIELD(background, F_COLOR),   FIELD(nan_color, F_COLOR),  FIELD(iso_color, F_COLOR),
    FIELD(vector_color, F_COLOR), FIELD(outline_color, F_COLOR), FIELD(cutaway_tint, F_COLOR),
    FIELD(stream_color, F_COLOR),
};
#define NUM_RENDER_FIELDS (sizeof(render_fields) / sizeof(render_fields[0]))

static int lookup(const char *what, const char *const *table, int n, const char *value) {
    int i;
    for (i = 0; i < n; ++i)
        if (strcmp(table[i], value) == 0) return i;
    Rf_error("meshio++: render: %s must be one of %s, %s, ..., not '%s'", what, table[0],
             table[1], value);
    return 0;
}

/* A colour: a number 0..2^32-1 as 0xRRGGBBAA (see mio_rgba), or c(r, g, b[, a]). */
static uint32_t color_of(SEXP s, const char *name) {
    if (TYPEOF(s) != REALSXP && TYPEOF(s) != INTSXP)
        Rf_error("meshio++: render: '%s' must be a number or c(r, g, b[, a])", name);
    if (XLENGTH(s) == 1) {
        double v = Rf_asReal(s);
        if (!(v >= 0.0 && v <= 4294967295.0))
            Rf_error("meshio++: render: '%s' must lie in [0, 2^32 - 1]", name);
        return (uint32_t)v;
    }
    if (XLENGTH(s) == 3 || XLENGTH(s) == 4) {
        SEXP d = PROTECT(Rf_coerceVector(s, REALSXP));
        uint32_t out = 0;
        R_xlen_t i;
        for (i = 0; i < 4; ++i) {
            double c = i < XLENGTH(d) ? REAL(d)[i] : 255.0;
            if (!(c >= 0.0 && c <= 255.0))
                Rf_error("meshio++: render: '%s' channels must lie in 0..255", name);
            out = (out << 8) | (uint32_t)c;
        }
        UNPROTECT(1);
        return out;
    }
    Rf_error("meshio++: render: '%s' must be a number or c(r, g, b[, a])", name);
    return 0;
}

/* Fill `o` from the named list `opts`. iso_levels is borrowed from a PROTECT-ed
 * REALSXP, counted in `*nprotect`: the caller UNPROTECTs that many objects once
 * the C API call has returned. */
static void fill_render_opts(SEXP opts, mio_render_opts *o, int *nprotect) {
    static const char *const shading_names[] = {"none", "flat", "smooth"};
    static const char *const edges_names[] = {"none", "all", "feature"};
    static const char *const scale_names[] = {"linear", "log", "symlog"};
    static const char *const diagnostic_names[] = {"none", "quality", "inverted", "degenerate",
                                                   "orientation", "free_edges", "edge_length"};
    SEXP names;
    R_xlen_t i;
    *nprotect = 0;
    mio_render_opts_init(o);
    if (opts == R_NilValue) return;
    if (TYPEOF(opts) != VECSXP) Rf_error("meshio++: render: options must be a named list");
    names = Rf_getAttrib(opts, R_NamesSymbol);
    if (XLENGTH(opts) > 0 && TYPEOF(names) != STRSXP)
        Rf_error("meshio++: render: every option must be named");
    for (i = 0; i < XLENGTH(opts); ++i) {
        const char *name = CHAR(STRING_ELT(names, i));
        SEXP v = VECTOR_ELT(opts, i);
        size_t k;
        int known = 0;
        if (v == R_NilValue) continue;
        for (k = 0; k < NUM_RENDER_FIELDS; ++k) {
            const field_def *f = &render_fields[k];
            char *base = (char *)o;
            if (strcmp(f->name, name) != 0) continue;
            known = 1;
            switch (f->kind) {
                case F_STR: *(const char **)(base + f->offset) = mio_r_string(v, name); break;
                case F_DBL: *(double *)(base + f->offset) = mio_r_double(v, name); break;
                case F_INT: *(int32_t *)(base + f->offset) = (int32_t)mio_r_int(v, name); break;
                case F_BOOL: *(int32_t *)(base + f->offset) = mio_r_bool(v, name) ? 1 : 0; break;
                case F_COLOR: *(uint32_t *)(base + f->offset) = color_of(v, name); break;
            }
        }
        if (known) continue;
        if (strcmp(name, "vmin") == 0) { o->vmin = mio_r_double(v, name); o->has_vmin = 1; }
        else if (strcmp(name, "vmax") == 0) { o->vmax = mio_r_double(v, name); o->has_vmax = 1; }
        else if (strcmp(name, "clip") == 0) {
            SEXP d;
            if (XLENGTH(v) != 2) Rf_error("meshio++: render: clip must be c(low, high) percentiles");
            d = PROTECT(Rf_coerceVector(v, REALSXP));
            if (!ISNAN(REAL(d)[0])) { o->clip_low = REAL(d)[0]; o->has_clip_low = 1; }
            if (!ISNAN(REAL(d)[1])) { o->clip_high = REAL(d)[1]; o->has_clip_high = 1; }
            UNPROTECT(1);
        } else if (strcmp(name, "pan") == 0) {
            SEXP d;
            if (XLENGTH(v) != 2) Rf_error("meshio++: render: pan must be c(x, y)");
            d = PROTECT(Rf_coerceVector(v, REALSXP));
            o->pan_x = REAL(d)[0];
            o->pan_y = REAL(d)[1];
            UNPROTECT(1);
        } else if (strcmp(name, "light_dir") == 0) {
            SEXP d;
            if (XLENGTH(v) != 3) Rf_error("meshio++: render: light_dir must be c(x, y, z)");
            d = PROTECT(Rf_coerceVector(v, REALSXP));
            o->light_x = REAL(d)[0];
            o->light_y = REAL(d)[1];
            o->light_z = REAL(d)[2];
            UNPROTECT(1);
        } else if (strcmp(name, "projection") == 0) {
            const char *p = mio_r_string(v, name);
            if (strcmp(p, "perspective") == 0) o->perspective = 1;
            else if (strcmp(p, "orthographic") == 0) o->perspective = 0;
            else Rf_error("meshio++: render: projection must be 'orthographic' or 'perspective'");
        } else if (strcmp(name, "shading") == 0) {
            o->shading = lookup("shading", shading_names, 3, mio_r_string(v, name));
        } else if (strcmp(name, "edges") == 0) {
            o->edges = lookup("edges", edges_names, 3, mio_r_string(v, name));
        } else if (strcmp(name, "scale") == 0) {
            o->scale = lookup("scale", scale_names, 3, mio_r_string(v, name));
        } else if (strcmp(name, "diagnostic") == 0) {
            o->diagnostic = lookup("diagnostic", diagnostic_names, 7, mio_r_string(v, name));
        } else if (strcmp(name, "cutaway") == 0) {
            /* One or two planes of six numbers: a point, then the normal of the side
             * kept. A plain vector of 6 or 12 numbers, or a list of such vectors. */
            SEXP d;
            if (TYPEOF(v) == VECSXP) {
                R_xlen_t j, total = 0;
                int ok = XLENGTH(v) >= 1 && XLENGTH(v) <= 2;
                for (j = 0; ok && j < XLENGTH(v); ++j) ok = XLENGTH(VECTOR_ELT(v, j)) == 6;
                if (!ok) Rf_error("meshio++: render: cutaway is one or two planes of six numbers each");
                d = PROTECT(Rf_allocVector(REALSXP, 6 * XLENGTH(v)));
                ++*nprotect;
                for (j = 0; j < XLENGTH(v); ++j) {
                    SEXP plane = PROTECT(Rf_coerceVector(VECTOR_ELT(v, j), REALSXP));
                    R_xlen_t t;
                    for (t = 0; t < 6; ++t) REAL(d)[total++] = REAL(plane)[t];
                    UNPROTECT(1);
                }
            } else {
                if (XLENGTH(v) != 6 && XLENGTH(v) != 12)
                    Rf_error("meshio++: render: cutaway is one or two planes of six numbers each");
                d = PROTECT(Rf_coerceVector(v, REALSXP));
                ++*nprotect;
            }
            o->cutaways = REAL(d);
            o->num_cutaways = (int32_t)(XLENGTH(d) / 6);
        } else if (strcmp(name, "iso_levels") == 0) {
            SEXP d = PROTECT(Rf_coerceVector(v, REALSXP));
            ++*nprotect;
            o->iso_levels = REAL(d);
            o->num_iso_levels = (int64_t)XLENGTH(d);
        } else {
            Rf_error("meshio++: render: unknown option '%s'", name);
        }
    }
}

static void fill_text_opts(SEXP t, mio_text_opts *o) {
    static const char *const encodings[] = {"halfblock", "quadrant", "sextant", "braille",
                                            "ascii", "kitty", "iterm2", "sixel"};
    static const char *const depths[] = {"truecolor", "256", "16", "mono"};
    static const char *const formats[] = {"ansi", "plain", "html"};
    SEXP names;
    R_xlen_t i;
    mio_text_opts_init(o);
    if (t == R_NilValue) return;
    if (TYPEOF(t) != VECSXP) Rf_error("meshio++: render: text options must be a named list");
    names = Rf_getAttrib(t, R_NamesSymbol);
    for (i = 0; i < XLENGTH(t); ++i) {
        const char *name = CHAR(STRING_ELT(names, i));
        SEXP v = VECTOR_ELT(t, i);
        if (v == R_NilValue) continue;
        if (strcmp(name, "encoding") == 0) o->encoding = lookup("encoding", encodings, 8, mio_r_string(v, name));
        else if (strcmp(name, "color_depth") == 0) o->color_depth = lookup("color_depth", depths, 4, mio_r_string(v, name));
        else if (strcmp(name, "format") == 0) o->format = lookup("format", formats, 3, mio_r_string(v, name));
        else if (strcmp(name, "cols") == 0) o->cols = mio_r_int(v, name);
        else if (strcmp(name, "rows") == 0) o->rows = mio_r_int(v, name);
        else if (strcmp(name, "cell_aspect") == 0) o->cell_aspect = mio_r_double(v, name);
        else if (strcmp(name, "cell_pixel_width") == 0) o->cell_pixel_width = mio_r_int(v, name);
        else if (strcmp(name, "cell_pixel_height") == 0) o->cell_pixel_height = mio_r_int(v, name);
        else if (strcmp(name, "tmux") == 0) o->tmux = mio_r_bool(v, name) ? 1 : 0;
        else if (strcmp(name, "notes") == 0) o->notes = mio_r_bool(v, name) ? 1 : 0;
        else Rf_error("meshio++: render: unknown text option '%s'", name);
    }
}

typedef struct {
    const mio_mesh *mesh;
    const mio_render_opts *ropts;
    const mio_text_opts *topts;
} text_ctx;

static int64_t text_getter(void *ctx, char *buf, int64_t buflen) {
    text_ctx *c = (text_ctx *)ctx;
    return mio_render_text(c->mesh, c->ropts, c->topts, buf, buflen);
}

SEXP R_mio_render(SEXP mesh, SEXP opts) {
    mio_render_opts o;
    mio_frame *f;
    int np;
    fill_render_opts(opts, &o, &np);
    f = mio_render(mio_r_mesh(mesh), &o);
    UNPROTECT(np);
    if (f == NULL) mio_r_fail("render");
    {
        const int w = mio_frame_width(f), h = mio_frame_height(f);
        const R_xlen_t n = (R_xlen_t)w * h;
        SEXP rgba = PROTECT(Rf_allocVector(INTSXP, 4 * n));
        SEXP ids = PROTECT(Rf_allocVector(REALSXP, n));
        SEXP dim, range, notes, sw, sh;
        const uint8_t *px = mio_frame_rgba(f);
        const int64_t *id = mio_frame_cell_ids(f);
        double vmin = 0.0, vmax = 0.0;
        int colored;
        R_xlen_t i, nn;
        for (i = 0; i < 4 * n; ++i) INTEGER(rgba)[i] = px[i];
        for (i = 0; i < n; ++i) REAL(ids)[i] = (double)id[i];
        dim = PROTECT(Rf_allocVector(INTSXP, 3));
        INTEGER(dim)[0] = 4; INTEGER(dim)[1] = w; INTEGER(dim)[2] = h;
        Rf_setAttrib(rgba, R_DimSymbol, dim);
        UNPROTECT(1);
        dim = PROTECT(Rf_allocVector(INTSXP, 2));
        INTEGER(dim)[0] = w; INTEGER(dim)[1] = h;
        Rf_setAttrib(ids, R_DimSymbol, dim);
        UNPROTECT(1);
        colored = mio_frame_range(f, &vmin, &vmax);
        if (colored == 1) {
            range = PROTECT(Rf_allocVector(REALSXP, 2));
            REAL(range)[0] = vmin; REAL(range)[1] = vmax;
        } else {
            range = PROTECT(R_NilValue);
        }
        nn = (R_xlen_t)mio_frame_num_notes(f);
        notes = PROTECT(Rf_allocVector(STRSXP, nn));
        for (i = 0; i < nn; ++i) {
            int64_t len = mio_frame_note(f, (int64_t)i, NULL, 0);
            char *buf = R_alloc((size_t)len + 1, 1);
            if (len < 0 || mio_frame_note(f, (int64_t)i, buf, len + 1) < 0) {
                mio_frame_free(f);
                mio_r_fail("render");
            }
            SET_STRING_ELT(notes, i, Rf_mkCharCE(buf, CE_UTF8));
        }
        sw = PROTECT(Rf_ScalarInteger(w));
        sh = PROTECT(Rf_ScalarInteger(h));
        mio_frame_free(f);
        {
            const char *names[] = {"width", "height", "rgba", "cell_ids", "range", "notes"};
            SEXP values[] = {sw, sh, rgba, ids, range, notes};
            SEXP res = PROTECT(mio_r_named_list(6, names, values));
            UNPROTECT(7);
            return res;
        }
    }
}

SEXP R_mio_render_text(SEXP mesh, SEXP opts, SEXP text) {
    mio_render_opts o;
    mio_text_opts t;
    text_ctx ctx;
    SEXP out;
    int np;
    fill_render_opts(opts, &o, &np);
    fill_text_opts(text, &t);
    ctx.mesh = mio_r_mesh(mesh);
    ctx.ropts = &o;
    ctx.topts = &t;
    {
        /* The text is UTF-8 whatever the session's locale (block characters). */
        int64_t n = text_getter(&ctx, NULL, 0);
        char *buf;
        if (n < 0) mio_r_fail("render_text");
        buf = R_alloc((size_t)n + 1, 1);
        if (n > 0 && text_getter(&ctx, buf, n + 1) < 0) mio_r_fail("render_text");
        buf[n] = '\0';
        out = PROTECT(Rf_ScalarString(Rf_mkCharLenCE(buf, (int)n, CE_UTF8)));
    }
    UNPROTECT(np + 1);
    return out;
}

SEXP R_mio_render_png(SEXP mesh, SEXP opts, SEXP compress) {
    mio_render_opts o;
    mio_frame *f;
    int np;
    int64_t n;
    SEXP out;
    fill_render_opts(opts, &o, &np);
    f = mio_render(mio_r_mesh(mesh), &o);
    UNPROTECT(np);
    if (f == NULL) mio_r_fail("render_png");
    n = mio_frame_png(f, (int32_t)mio_r_int(compress, "compress"), NULL, 0);
    if (n < 0) { mio_frame_free(f); mio_r_fail("render_png"); }
    out = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)n));
    if (mio_frame_png(f, (int32_t)mio_r_int(compress, "compress"), RAW(out), n) < 0) {
        mio_frame_free(f);
        mio_r_fail("render_png");
    }
    mio_frame_free(f);
    UNPROTECT(1);
    return out;
}

SEXP R_mio_write_snapshot(SEXP path, SEXP mesh, SEXP opts, SEXP text, SEXP png_compress,
                          SEXP cast_frames, SEXP cast_fps, SEXP cast_degrees) {
    mio_render_opts o;
    mio_text_opts t;
    mio_snapshot_opts s;
    int np;
    mio_status st;
    fill_render_opts(opts, &o, &np);
    fill_text_opts(text, &t);
    mio_snapshot_opts_init(&s);
    s.png_compress = mio_r_int(png_compress, "png_compress");
    s.cast_frames = mio_r_int(cast_frames, "cast_frames");
    s.cast_fps = mio_r_double(cast_fps, "cast_fps");
    s.cast_degrees = mio_r_double(cast_degrees, "cast_degrees");
    st = mio_write_snapshot(mio_r_string(path, "path"), mio_r_mesh(mesh), &o, &t, &s);
    UNPROTECT(np);
    mio_r_check(st, "write_snapshot");
    return R_NilValue;
}

SEXP R_mio_detect_color_depth(SEXP no_color, SEXP color_term, SEXP term) {
    static const char *const depths[] = {"truecolor", "256", "16", "mono"};
    int32_t d = mio_detect_color_depth(mio_r_opt_string(no_color), mio_r_opt_string(color_term),
                                       mio_r_opt_string(term));
    return Rf_mkString(depths[d < 0 || d > 3 ? 3 : d]);
}
