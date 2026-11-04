/* Software rendering (v16.33.0, field rendering v16.34.0): a deterministic
 * software rasterizer and its terminal, HTML and PNG encodings. See
 * doc/tui.md. A port of the R binding's src/render.c; see mio_oct.h for the
 * layout rules.
 *
 * Every entry point takes the options as one 1 x 2k cell of name/value pairs
 * (the counterpart of R's named list; [] for none), so adding a field to
 * mio_render_opts is a row in a table here, not a change of signature. An
 * unknown name is an error naming the ones that exist: a typo never silently
 * does nothing. A [] value is skipped, like R's NULL list element.
 *
 * Text (the terminal text and the frame's notes) is UTF-8 and crosses as
 * uint8 byte rows; the .m wrappers decode it, so block characters survive in
 * MATLAB whatever its locale. */

#include "mio_oct.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

typedef enum { ORN_STR, ORN_DBL, ORN_INT, ORN_BOOL, ORN_COLOR } orn_kind;

typedef struct {
    const char *name;
    orn_kind kind;
    size_t offset;
} orn_field;

#define ORN_FIELD(n, k) {#n, k, offsetof(mio_render_opts, n)}

static const orn_field orn_fields[] = {
    ORN_FIELD(view, ORN_STR),           ORN_FIELD(color_by, ORN_STR),
    ORN_FIELD(cmap, ORN_STR),           ORN_FIELD(reduce, ORN_STR),
    ORN_FIELD(expr, ORN_STR),           ORN_FIELD(vectors, ORN_STR),
    ORN_FIELD(warp, ORN_STR),           ORN_FIELD(quality_metric, ORN_STR),
    ORN_FIELD(streamlines, ORN_STR),    ORN_FIELD(pixel_aspect, ORN_DBL),
    ORN_FIELD(azimuth, ORN_DBL),        ORN_FIELD(elevation, ORN_DBL),
    ORN_FIELD(roll, ORN_DBL),           ORN_FIELD(fov_deg, ORN_DBL),
    ORN_FIELD(zoom, ORN_DBL),           ORN_FIELD(ambient, ORN_DBL),
    ORN_FIELD(split_angle, ORN_DBL),    ORN_FIELD(feature_angle, ORN_DBL),
    ORN_FIELD(point_radius, ORN_DBL),   ORN_FIELD(scale_threshold, ORN_DBL),
    ORN_FIELD(vector_length, ORN_DBL),  ORN_FIELD(warp_scale, ORN_DBL),
    ORN_FIELD(stream_length, ORN_DBL),  ORN_FIELD(width, ORN_INT),
    ORN_FIELD(height, ORN_INT),         ORN_FIELD(supersample, ORN_INT),
    ORN_FIELD(component, ORN_INT),      ORN_FIELD(isolines, ORN_INT),
    ORN_FIELD(vector_count, ORN_INT),   ORN_FIELD(stream_seeds, ORN_INT),
    ORN_FIELD(perspective, ORN_BOOL),   ORN_FIELD(two_sided, ORN_BOOL),
    ORN_FIELD(colorbar, ORN_BOOL),      ORN_FIELD(axes, ORN_BOOL),
    ORN_FIELD(scale_bar, ORN_BOOL),     ORN_FIELD(symmetric, ORN_BOOL),
    ORN_FIELD(categorical, ORN_BOOL),   ORN_FIELD(color_regions, ORN_BOOL),
    ORN_FIELD(category_edges, ORN_BOOL), ORN_FIELD(warp_outline, ORN_BOOL),
    ORN_FIELD(edge_color, ORN_COLOR),   ORN_FIELD(fill_color, ORN_COLOR),
    ORN_FIELD(line_color, ORN_COLOR),   ORN_FIELD(background, ORN_COLOR),
    ORN_FIELD(nan_color, ORN_COLOR),    ORN_FIELD(iso_color, ORN_COLOR),
    ORN_FIELD(vector_color, ORN_COLOR), ORN_FIELD(outline_color, ORN_COLOR),
    ORN_FIELD(cutaway_tint, ORN_COLOR), ORN_FIELD(stream_color, ORN_COLOR),
};
#define ORN_NUM_FIELDS (sizeof(orn_fields) / sizeof(orn_fields[0]))

static int orn_lookup(const char *what, const char *const *table, int n, const char *value) {
    for (int i = 0; i < n; ++i)
        if (strcmp(table[i], value) == 0) return i;
    mio_o_error("render", "meshio++: render: %s must be one of %s, %s, ..., not '%s'", what,
                table[0], table[1], value);
    return 0;
}

/* A colour: a number 0..2^32-1 as 0xRRGGBBAA (see mio.rgba), or [r g b (a)]. */
static uint32_t orn_color_of(const mxArray *s, const char *name) {
    if (!(mxIsNumeric(s) || mxIsLogical(s)) || mxIsComplex(s))
        mio_o_error("render", "meshio++: render: '%s' must be a number or [r g b (a)]", name);
    size_t n = 0;
    double *d = mio_o_doubles(s, &n, name);
    if (n == 1) {
        double v = d[0];
        if (!(v >= 0.0 && v <= 4294967295.0))
            mio_o_error("render", "meshio++: render: '%s' must lie in [0, 2^32 - 1]", name);
        return (uint32_t)v;
    }
    if (n == 3 || n == 4) {
        uint32_t out = 0;
        for (size_t i = 0; i < 4; ++i) {
            double c = i < n ? d[i] : 255.0;
            if (!(c >= 0.0 && c <= 255.0))
                mio_o_error("render", "meshio++: render: '%s' channels must lie in 0..255", name);
            out = (out << 8) | (uint32_t)c;
        }
        return out;
    }
    mio_o_error("render", "meshio++: render: '%s' must be a number or [r g b (a)]", name);
    return 0;
}

/* Exactly `n` doubles of a numeric option, else the error `msg`. */
static double *orn_fixed(const mxArray *v, size_t n, const char *name, const char *msg) {
    if (!(mxIsNumeric(v) || mxIsLogical(v)) || mxGetNumberOfElements(v) != n)
        mio_o_error("render", "meshio++: render: %s", msg);
    size_t got = 0;
    return mio_o_doubles(v, &got, name);
}

/* The name/value pairs of an options argument: [] for none, or a cell with an
 * even number of elements whose odd ones are names. */
static size_t orn_pairs(const mxArray *opts, const char *what) {
    if (mio_o_is_null(opts)) return 0;
    if (!mxIsCell(opts) || mxGetNumberOfElements(opts) % 2 != 0)
        mio_o_error("render", "meshio++: render: %s must be name/value pairs", what);
    size_t n = mxGetNumberOfElements(opts) / 2;
    for (size_t i = 0; i < n; ++i) {
        const mxArray *nm = mxGetCell(opts, (mwIndex)(2 * i));
        if (nm == NULL || !mxIsChar(nm) || mxIsEmpty(nm))
            mio_o_error("render", "meshio++: render: every option must be named");
    }
    return n;
}

/* Fill `o` from the pairs in `opts`. Every buffer it points into (strings,
 * cutaways, iso levels) is mxMalloc'd and lives until the MEX call returns,
 * so no bookkeeping is needed -- the counterpart of R's PROTECT count. */
static void orn_fill_render_opts(const mxArray *opts, mio_render_opts *o) {
    static const char *const shading_names[] = {"none", "flat", "smooth"};
    static const char *const edges_names[] = {"none", "all", "feature"};
    static const char *const scale_names[] = {"linear", "log", "symlog"};
    static const char *const diagnostic_names[] = {"none", "quality", "inverted", "degenerate",
                                                   "orientation", "free_edges", "edge_length"};
    mio_render_opts_init(o);
    size_t npairs = orn_pairs(opts, "options");
    for (size_t i = 0; i < npairs; ++i) {
        const char *name = mio_o_string(mxGetCell(opts, (mwIndex)(2 * i)), "option name");
        const mxArray *v = mxGetCell(opts, (mwIndex)(2 * i + 1));
        int known = 0;
        if (v == NULL || mio_o_is_null(v)) continue;
        for (size_t k = 0; k < ORN_NUM_FIELDS; ++k) {
            const orn_field *f = &orn_fields[k];
            char *base = (char *)o;
            if (strcmp(f->name, name) != 0) continue;
            known = 1;
            switch (f->kind) {
            case ORN_STR: *(const char **)(base + f->offset) = mio_o_string(v, name); break;
            case ORN_DBL: *(double *)(base + f->offset) = mio_o_double(v, name); break;
            case ORN_INT: *(int32_t *)(base + f->offset) = (int32_t)mio_o_int(v, name); break;
            case ORN_BOOL: *(int32_t *)(base + f->offset) = mio_o_bool(v, name) ? 1 : 0; break;
            case ORN_COLOR: *(uint32_t *)(base + f->offset) = orn_color_of(v, name); break;
            }
        }
        if (known) continue;
        if (strcmp(name, "vmin") == 0) {
            o->vmin = mio_o_double(v, name);
            o->has_vmin = 1;
        } else if (strcmp(name, "vmax") == 0) {
            o->vmax = mio_o_double(v, name);
            o->has_vmax = 1;
        } else if (strcmp(name, "clip") == 0) {
            /* [low high] percentiles; NaN leaves that end unclipped. */
            double *d = orn_fixed(v, 2, name, "clip must be [low high] percentiles");
            if (!isnan(d[0])) {
                o->clip_low = d[0];
                o->has_clip_low = 1;
            }
            if (!isnan(d[1])) {
                o->clip_high = d[1];
                o->has_clip_high = 1;
            }
        } else if (strcmp(name, "pan") == 0) {
            double *d = orn_fixed(v, 2, name, "pan must be [x y]");
            o->pan_x = d[0];
            o->pan_y = d[1];
        } else if (strcmp(name, "light_dir") == 0) {
            double *d = orn_fixed(v, 3, name, "light_dir must be [x y z]");
            o->light_x = d[0];
            o->light_y = d[1];
            o->light_z = d[2];
        } else if (strcmp(name, "projection") == 0) {
            const char *p = mio_o_string(v, name);
            if (strcmp(p, "perspective") == 0)
                o->perspective = 1;
            else if (strcmp(p, "orthographic") == 0)
                o->perspective = 0;
            else
                mio_o_error("render",
                            "meshio++: render: projection must be 'orthographic' or 'perspective'");
        } else if (strcmp(name, "shading") == 0) {
            o->shading = orn_lookup("shading", shading_names, 3, mio_o_string(v, name));
        } else if (strcmp(name, "edges") == 0) {
            o->edges = orn_lookup("edges", edges_names, 3, mio_o_string(v, name));
        } else if (strcmp(name, "scale") == 0) {
            o->scale = orn_lookup("scale", scale_names, 3, mio_o_string(v, name));
        } else if (strcmp(name, "diagnostic") == 0) {
            o->diagnostic = orn_lookup("diagnostic", diagnostic_names, 7, mio_o_string(v, name));
        } else if (strcmp(name, "cutaway") == 0) {
            /* One or two planes of six numbers: a point, then the normal of the side
             * kept. A plain vector of 6 or 12 numbers, or a cell of such vectors. */
            static const char *const bad =
                "cutaway is one or two planes of six numbers each";
            double *d;
            size_t total = 0;
            if (mxIsCell(v)) {
                size_t np = mxGetNumberOfElements(v);
                int ok = np >= 1 && np <= 2;
                for (size_t j = 0; ok && j < np; ++j) {
                    const mxArray *plane = mxGetCell(v, (mwIndex)j);
                    ok = plane != NULL && (mxIsNumeric(plane) || mxIsLogical(plane)) &&
                         mxGetNumberOfElements(plane) == 6;
                }
                if (!ok) mio_o_error("render", "meshio++: render: %s", bad);
                d = (double *)mxMalloc(6 * np * sizeof(double));
                for (size_t j = 0; j < np; ++j) {
                    size_t got = 0;
                    double *plane = mio_o_doubles(mxGetCell(v, (mwIndex)j), &got, name);
                    for (size_t t = 0; t < 6; ++t) d[total++] = plane[t];
                }
            } else {
                size_t nv = mxGetNumberOfElements(v);
                if (!(mxIsNumeric(v) || mxIsLogical(v)) || (nv != 6 && nv != 12))
                    mio_o_error("render", "meshio++: render: %s", bad);
                d = mio_o_doubles(v, &total, name);
            }
            o->cutaways = d;
            o->num_cutaways = (int32_t)(total / 6);
        } else if (strcmp(name, "iso_levels") == 0) {
            size_t n = 0;
            o->iso_levels = mio_o_doubles(v, &n, name);
            o->num_iso_levels = (int64_t)n;
        } else {
            mio_o_error("render", "meshio++: render: unknown option '%s'", name);
        }
    }
}

static void orn_fill_text_opts(const mxArray *t, mio_text_opts *o) {
    static const char *const encodings[] = {"halfblock", "quadrant", "sextant", "braille",
                                            "ascii",     "kitty",    "iterm2",  "sixel"};
    static const char *const depths[] = {"truecolor", "256", "16", "mono"};
    static const char *const formats[] = {"ansi", "plain", "html"};
    mio_text_opts_init(o);
    size_t npairs = orn_pairs(t, "text options");
    for (size_t i = 0; i < npairs; ++i) {
        const char *name = mio_o_string(mxGetCell(t, (mwIndex)(2 * i)), "option name");
        const mxArray *v = mxGetCell(t, (mwIndex)(2 * i + 1));
        if (v == NULL || mio_o_is_null(v)) continue;
        if (strcmp(name, "encoding") == 0)
            o->encoding = orn_lookup("encoding", encodings, 8, mio_o_string(v, name));
        else if (strcmp(name, "color_depth") == 0)
            o->color_depth = orn_lookup("color_depth", depths, 4, mio_o_string(v, name));
        else if (strcmp(name, "format") == 0)
            o->format = orn_lookup("format", formats, 3, mio_o_string(v, name));
        else if (strcmp(name, "cols") == 0)
            o->cols = mio_o_int(v, name);
        else if (strcmp(name, "rows") == 0)
            o->rows = mio_o_int(v, name);
        else if (strcmp(name, "cell_aspect") == 0)
            o->cell_aspect = mio_o_double(v, name);
        else if (strcmp(name, "cell_pixel_width") == 0)
            o->cell_pixel_width = mio_o_int(v, name);
        else if (strcmp(name, "cell_pixel_height") == 0)
            o->cell_pixel_height = mio_o_int(v, name);
        else if (strcmp(name, "tmux") == 0)
            o->tmux = mio_o_bool(v, name) ? 1 : 0;
        else if (strcmp(name, "notes") == 0)
            o->notes = mio_o_bool(v, name) ? 1 : 0;
        else
            mio_o_error("render", "meshio++: render: unknown text option '%s'", name);
    }
}

/* A 1 x n uint8 row holding `n` bytes. */
static mxArray *orn_bytes(const void *data, int64_t n) {
    mxArray *out = mxCreateNumericMatrix(1, (mwSize)(n > 0 ? n : 0), mxUINT8_CLASS, mxREAL);
    if (n > 0) memcpy(mxGetData(out), data, (size_t)n);
    return out;
}

/* render(mesh, opts) -> struct width, height, rgba, cell_ids, range, notes.
 *
 * rgba is height x width x 4 uint8, straight alpha, row 1 the TOP of the
 * picture (the frame's row order), so image(f.rgba(:,:,1:3)) and imwrite
 * show it upright. cell_ids is height x width in the same orientation: the
 * 1-based global (block-major) cell drawn at each pixel, 0 for none -- the
 * index-map convention of this binding. */
static mxArray *O_mio_render(const mxArray **a) {
    mio_render_opts o;
    orn_fill_render_opts(a[1], &o);
    mio_frame *f = mio_render(mio_o_mesh(a[0]), &o);
    if (f == NULL) mio_o_fail("render");

    const int64_t w = mio_frame_width(f), h = mio_frame_height(f);
    const uint8_t *px = mio_frame_rgba(f);
    const int64_t *id = mio_frame_cell_ids(f);

    mwSize dims[3] = {(mwSize)h, (mwSize)w, 4};
    mxArray *rgba = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
    mxArray *ids = mio_o_matrix(h, w);
    uint8_t *dst = (uint8_t *)mxGetData(rgba);
    double *idd = mxGetPr(ids);
    for (int64_t y = 0; y < h; ++y) {
        for (int64_t x = 0; x < w; ++x) {
            const int64_t src = y * w + x;   /* row-major, row 0 on top */
            const int64_t cm = x * h + y;    /* column-major (y, x) */
            for (int64_t c = 0; c < 4; ++c) dst[c * w * h + cm] = px[src * 4 + c];
            idd[cm] = id[src] < 0 ? 0.0 : (double)(id[src] + 1);
        }
    }

    double vmin = 0.0, vmax = 0.0;
    mxArray *range;
    if (mio_frame_range(f, &vmin, &vmax) == 1) {
        range = mio_o_matrix(1, 2);
        mxGetPr(range)[0] = vmin;
        mxGetPr(range)[1] = vmax;
    } else {
        range = mio_o_empty();
    }

    int64_t nn = mio_frame_num_notes(f);
    if (nn < 0) nn = 0;
    mxArray *notes = mio_o_cell(nn);
    for (int64_t i = 0; i < nn; ++i) {
        int64_t len = mio_frame_note(f, i, NULL, 0);
        if (len < 0) {
            mio_frame_free(f);
            mio_o_fail("render");
        }
        char *buf = (char *)mxCalloc((size_t)len + 1, 1);
        if (mio_frame_note(f, i, buf, len + 1) < 0) {
            mio_frame_free(f);
            mio_o_fail("render");
        }
        mxSetCell(notes, (mwIndex)i, orn_bytes(buf, len));
        mxFree(buf);
    }
    mio_frame_free(f);

    const char *names[] = {"width", "height", "rgba", "cell_ids", "range", "notes"};
    mxArray *values[] = {mio_o_scalar((double)w), mio_o_scalar((double)h), rgba, ids, range,
                         notes};
    return mio_o_struct(6, names, values);
}

typedef struct {
    const mio_mesh *mesh;
    const mio_render_opts *ropts;
    const mio_text_opts *topts;
} orn_text_ctx;

static int64_t orn_text_getter(void *ctx, char *buf, int64_t buflen) {
    orn_text_ctx *c = (orn_text_ctx *)ctx;
    return mio_render_text(c->mesh, c->ropts, c->topts, buf, buflen);
}

/* render_text(mesh, opts, text) -> the UTF-8 text as a 1 x n uint8 row. */
static mxArray *O_mio_render_text(const mxArray **a) {
    mio_render_opts o;
    mio_text_opts t;
    orn_fill_render_opts(a[1], &o);
    orn_fill_text_opts(a[2], &t);
    orn_text_ctx ctx = {mio_o_mesh(a[0]), &o, &t};
    int64_t n = orn_text_getter(&ctx, NULL, 0);
    if (n < 0) mio_o_fail("render_text");
    char *buf = (char *)mxMalloc((size_t)n + 1);
    if (n > 0 && orn_text_getter(&ctx, buf, n + 1) < 0) mio_o_fail("render_text");
    mxArray *out = orn_bytes(buf, n);
    mxFree(buf);
    return out;
}

/* render_png(mesh, opts, compress) -> the PNG file bytes as an n x 1 uint8. */
static mxArray *O_mio_render_png(const mxArray **a) {
    mio_render_opts o;
    orn_fill_render_opts(a[1], &o);
    const int32_t compress = (int32_t)mio_o_int(a[2], "compress");
    mio_frame *f = mio_render(mio_o_mesh(a[0]), &o);
    if (f == NULL) mio_o_fail("render_png");
    int64_t n = mio_frame_png(f, compress, NULL, 0);
    if (n < 0) {
        mio_frame_free(f);
        mio_o_fail("render_png");
    }
    mxArray *out = mxCreateNumericMatrix((mwSize)n, 1, mxUINT8_CLASS, mxREAL);
    if (mio_frame_png(f, compress, (uint8_t *)mxGetData(out), n) < 0) {
        mio_frame_free(f);
        mxDestroyArray(out);
        mio_o_fail("render_png");
    }
    mio_frame_free(f);
    return out;
}

/* write_snapshot(path, mesh, opts, text, png_compress, cast_frames, cast_fps,
 *                cast_degrees) */
static mxArray *O_mio_write_snapshot(const mxArray **a) {
    mio_render_opts o;
    mio_text_opts t;
    mio_snapshot_opts s;
    orn_fill_render_opts(a[2], &o);
    orn_fill_text_opts(a[3], &t);
    mio_snapshot_opts_init(&s);
    s.png_compress = mio_o_int(a[4], "png_compress");
    s.cast_frames = mio_o_int(a[5], "cast_frames");
    s.cast_fps = mio_o_double(a[6], "cast_fps");
    s.cast_degrees = mio_o_double(a[7], "cast_degrees");
    mio_o_check(mio_write_snapshot(mio_o_string(a[0], "path"), mio_o_mesh(a[1]), &o, &t, &s),
                "write_snapshot");
    return NULL;
}

static mxArray *O_mio_detect_color_depth(const mxArray **a) {
    static const char *const depths[] = {"truecolor", "256", "16", "mono"};
    int32_t d = mio_detect_color_depth(mio_o_opt_string(a[0]), mio_o_opt_string(a[1]),
                                       mio_o_opt_string(a[2]));
    return mio_o_str(depths[d < 0 || d > 3 ? 3 : d]);
}

static const mio_o_entry kEntries[] = {
    MIO_O_ENTRY(mio_render, 2),          MIO_O_ENTRY(mio_render_text, 3),
    MIO_O_ENTRY(mio_render_png, 3),      MIO_O_ENTRY(mio_write_snapshot, 8),
    MIO_O_ENTRY(mio_detect_color_depth, 3),
};

const mio_o_entry *mio_o_table_render(size_t *n) {
    *n = sizeof(kEntries) / sizeof(kEntries[0]);
    return kEntries;
}
