/* The transient (time-series) XDMF writer.
 * A port of the R binding's src/xdmf_series.c; see mio_oct.h for the layout
 * rules.
 *
 * The one writer `mio.write()` cannot express: a series is a stateful
 * multi-call object -- the grid goes out once and each solve appends a cheap
 * step -- so it gets its own handle kind (MIO_O_XDMF_SERIES, class
 * mio.XdmfSeries), exactly as a mesh does. A foreign or already-released
 * handle is an error, never a dereference. The registry frees a series with
 * mio_xdmf_series_free(), which finalizes the .xdmf if the caller never did; a
 * failure there cannot be reported from a destructor, which is why
 * mio.xdmf_series_finalize() is exposed separately. */

#include "mio_oct.h"

#include <string.h>

static mio_xdmf_series *oxs_series_of(const mxArray *x) {
    return (mio_xdmf_series *)mio_o_unwrap(x, MIO_O_XDMF_SERIES, "series");
}

/* xdmf_series_create(path, data_format, gzip_level, mode, auto_flush) */
static mxArray *O_mio_xdmf_series_create(const mxArray **a) {
    const char *fmt = mio_o_opt_string(a[1]); /* NULL means "HDF" */
    const char *md = mio_o_string(a[3], "mode");
    mio_xdmf_series_opts opts;
    mio_xdmf_series_opts_init(&opts);
    opts.data_format = fmt;
    opts.gzip_level = (int32_t)mio_o_int(a[2], "gzip_level");
    if (strcmp(md, "append") == 0) {
        opts.mode = MIO_XDMF_SERIES_APPEND;
    } else if (strcmp(md, "truncate") != 0) {
        mio_o_error("xdmf_series", "mode must be 'truncate' or 'append'");
    }
    opts.auto_flush = mio_o_bool(a[4], "auto_flush");
    mio_xdmf_series *series = mio_xdmf_series_create_ex(mio_o_string(a[0], "path"), &opts);
    if (series == NULL) mio_o_fail("xdmf_series");
    return mio_o_wrap(MIO_O_XDMF_SERIES, series);
}

static mxArray *O_mio_xdmf_series_write_points_cells(const mxArray **a) {
    mio_o_check(mio_xdmf_series_write_points_cells(oxs_series_of(a[0]), mio_o_mesh(a[1])),
                "xdmf_series_write_points_cells");
    return NULL;
}

static mxArray *O_mio_xdmf_series_write_data(const mxArray **a) {
    mio_o_check(mio_xdmf_series_write_data(oxs_series_of(a[0]), mio_o_double(a[1], "time"),
                                           mio_o_mesh(a[2])),
                "xdmf_series_write_data");
    return NULL;
}

static mxArray *O_mio_xdmf_series_finalize(const mxArray **a) {
    mio_o_check(mio_xdmf_series_finalize(oxs_series_of(a[0])), "xdmf_series_finalize");
    return NULL;
}

static mxArray *O_mio_xdmf_series_flush(const mxArray **a) {
    mio_o_check(mio_xdmf_series_flush(oxs_series_of(a[0])), "xdmf_series_flush");
    return NULL;
}

static mxArray *O_mio_xdmf_series_finalized(const mxArray **a) {
    int32_t f = mio_xdmf_series_finalized(oxs_series_of(a[0]));
    if (f < 0) mio_o_fail("xdmf_series_finalized");
    return mio_o_logical(f == 1);
}

static mxArray *O_mio_xdmf_series_num_steps(const mxArray **a) {
    int64_t n = mio_xdmf_series_num_steps(oxs_series_of(a[0]));
    if (n < 0) mio_o_fail("xdmf_series_num_steps");
    /* Counts cross as double; a step count is far inside its exact range. */
    return mio_o_scalar((double)n);
}

/* Idempotent: an already-released series is a no-op; a non-series raises. */
static mxArray *O_mio_xdmf_series_release(const mxArray **a) {
    mio_o_release(a[0], MIO_O_XDMF_SERIES);
    return NULL;
}

static mxArray *O_mio_xdmf_series_is_open(const mxArray **a) {
    return mio_o_logical(mio_o_is_open(a[0], MIO_O_XDMF_SERIES));
}

static const mio_o_entry kEntries[] = {
    MIO_O_ENTRY(mio_xdmf_series_create, 5),
    MIO_O_ENTRY(mio_xdmf_series_flush, 1),
    MIO_O_ENTRY(mio_xdmf_series_finalized, 1),
    MIO_O_ENTRY(mio_xdmf_series_write_points_cells, 2),
    MIO_O_ENTRY(mio_xdmf_series_write_data, 3),
    MIO_O_ENTRY(mio_xdmf_series_finalize, 1),
    MIO_O_ENTRY(mio_xdmf_series_num_steps, 1),
    MIO_O_ENTRY(mio_xdmf_series_release, 1),
    MIO_O_ENTRY(mio_xdmf_series_is_open, 1),
};

const mio_o_entry *mio_o_table_xdmf_series(size_t *n) {
    *n = sizeof(kEntries) / sizeof(kEntries[0]);
    return kEntries;
}
