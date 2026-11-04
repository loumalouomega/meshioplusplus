/* Owning Exodus series; a distinct handle kind (MIO_O_EXODUS_SERIES, class
 * mio.ExodusSeries) rejects mesh/XDMF/foreign handles. A port of the R
 * binding's src/exodus_series.c; see mio_oct.h for the layout rules.
 *
 * Without netCDF the C library refuses to create a series, and its error
 * names the MESHIOPLUSPLUS_WITH_NETCDF dependency; that message is raised
 * unchanged. */

#include "mio_oct.h"

static mio_exodus_series *oex_series_of(const mxArray *x) {
    return (mio_exodus_series *)mio_o_unwrap(x, MIO_O_EXODUS_SERIES, "series");
}

static mxArray *O_mio_exodus_series_create(const mxArray **a) {
    mio_exodus_series *s = mio_exodus_series_create(mio_o_string(a[0], "path"));
    if (s == NULL) mio_o_fail("exodus_series");
    return mio_o_wrap(MIO_O_EXODUS_SERIES, s);
}

static mxArray *O_mio_exodus_series_write_points_cells(const mxArray **a) {
    mio_o_check(mio_exodus_series_write_points_cells(oex_series_of(a[0]), mio_o_mesh(a[1])),
                "exodus_series grid");
    return NULL;
}

static mxArray *O_mio_exodus_series_write_data(const mxArray **a) {
    mio_o_check(mio_exodus_series_write_data(oex_series_of(a[0]), mio_o_double(a[1], "time"),
                                             mio_o_mesh(a[2])),
                "exodus_series data");
    return NULL;
}

static mxArray *O_mio_exodus_series_flush(const mxArray **a) {
    mio_o_check(mio_exodus_series_flush(oex_series_of(a[0])), "exodus_series flush");
    return NULL;
}

static mxArray *O_mio_exodus_series_finalize(const mxArray **a) {
    mio_o_check(mio_exodus_series_finalize(oex_series_of(a[0])), "exodus_series finalize");
    return NULL;
}

static mxArray *O_mio_exodus_series_num_steps(const mxArray **a) {
    int64_t n = mio_exodus_series_num_steps(oex_series_of(a[0]));
    if (n < 0) mio_o_fail("exodus_series steps");
    return mio_o_scalar((double)n);
}

static mxArray *O_mio_exodus_series_finalized(const mxArray **a) {
    int32_t f = mio_exodus_series_finalized(oex_series_of(a[0]));
    if (f < 0) mio_o_fail("exodus_series finalized");
    return mio_o_logical(f == 1);
}

/* Idempotent: an already-released series is a no-op; a non-series raises. */
static mxArray *O_mio_exodus_series_release(const mxArray **a) {
    mio_o_release(a[0], MIO_O_EXODUS_SERIES);
    return NULL;
}

static mxArray *O_mio_exodus_series_is_open(const mxArray **a) {
    return mio_o_logical(mio_o_is_open(a[0], MIO_O_EXODUS_SERIES));
}

static const mio_o_entry kEntries[] = {
    MIO_O_ENTRY(mio_exodus_series_create, 1),
    MIO_O_ENTRY(mio_exodus_series_write_points_cells, 2),
    MIO_O_ENTRY(mio_exodus_series_write_data, 3),
    MIO_O_ENTRY(mio_exodus_series_flush, 1),
    MIO_O_ENTRY(mio_exodus_series_finalize, 1),
    MIO_O_ENTRY(mio_exodus_series_num_steps, 1),
    MIO_O_ENTRY(mio_exodus_series_finalized, 1),
    MIO_O_ENTRY(mio_exodus_series_release, 1),
    MIO_O_ENTRY(mio_exodus_series_is_open, 1),
};

const mio_o_entry *mio_o_table_exodus_series(size_t *n) {
    *n = sizeof(kEntries) / sizeof(kEntries[0]);
    return kEntries;
}
