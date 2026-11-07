/* Owning Exodus series; a distinct tag rejects mesh/XDMF/foreign pointers. */
#include "mio_r.h"

static SEXP exodus_tag(void) { return Rf_install("mio_exodus_series_handle"); }

static mio_exodus_series *exodus_of(SEXP x) {
    return (mio_exodus_series *)mio_r_handle(x, exodus_tag(), "mio_exodus_series");
}

static void exodus_finalizer(SEXP x) {
    mio_exodus_series_free((mio_exodus_series *)R_ExternalPtrAddr(x));
    R_ClearExternalPtr(x);
}

SEXP R_mio_exodus_series_create(SEXP path) {
    mio_exodus_series *s = mio_exodus_series_create(mio_r_string(path, "path"));
    if (!s) mio_r_fail("exodus_series");
    return mio_r_wrap_handle(s, exodus_tag(), exodus_finalizer, "mio_exodus_series");
}

SEXP R_mio_exodus_series_write_points_cells(SEXP series, SEXP mesh) {
    mio_r_check(mio_exodus_series_write_points_cells(exodus_of(series), mio_r_mesh(mesh)), "exodus_series grid");
    return R_NilValue;
}
SEXP R_mio_exodus_series_write_data(SEXP series, SEXP time, SEXP mesh) {
    mio_r_check(mio_exodus_series_write_data(exodus_of(series), mio_r_double(time, "time"), mio_r_mesh(mesh)), "exodus_series data");
    return R_NilValue;
}
SEXP R_mio_exodus_series_flush(SEXP series) {
    mio_r_check(mio_exodus_series_flush(exodus_of(series)), "exodus_series flush");
    return R_NilValue;
}
SEXP R_mio_exodus_series_finalize(SEXP series) {
    mio_r_check(mio_exodus_series_finalize(exodus_of(series)), "exodus_series finalize");
    return R_NilValue;
}
SEXP R_mio_exodus_series_num_steps(SEXP series) {
    int64_t n = mio_exodus_series_num_steps(exodus_of(series));
    if (n < 0) mio_r_fail("exodus_series steps");
    return Rf_ScalarReal((double)n);
}
SEXP R_mio_exodus_series_finalized(SEXP series) {
    int32_t f = mio_exodus_series_finalized(exodus_of(series));
    if (f < 0) mio_r_fail("exodus_series finalized");
    return Rf_ScalarLogical(f == 1);
}
SEXP R_mio_exodus_series_release(SEXP series) {
    if (TYPEOF(series) != EXTPTRSXP || R_ExternalPtrTag(series) != exodus_tag())
        Rf_error("expected a mio_exodus_series object");
    exodus_finalizer(series);
    return R_NilValue;
}
SEXP R_mio_exodus_series_is_open(SEXP series) {
    return Rf_ScalarLogical(TYPEOF(series) == EXTPTRSXP && R_ExternalPtrTag(series) == exodus_tag() && R_ExternalPtrAddr(series) != NULL);
}
