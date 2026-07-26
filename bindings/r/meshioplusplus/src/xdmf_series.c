/* The transient (time-series) XDMF writer.
 *
 * The one writer `mio_write()` cannot express: a series is a stateful
 * multi-call object -- the grid goes out once and each solve appends a cheap
 * step -- so it gets its own external pointer here, with its own tag and
 * finalizer, exactly as `mio_mesh` does. A foreign or already-released pointer
 * is an R error, never a dereference. */

#include "mio_r.h"

SEXP mio_r_series_tag = NULL;

static void mio_r_series_finalizer(SEXP x) {
    mio_xdmf_series *series = (mio_xdmf_series *)R_ExternalPtrAddr(x);
    if (series != NULL) {
        /* Finalizes the .xdmf if the caller never did; a failure there cannot
         * be reported from a GC finalizer, which is why
         * mio_xdmf_series_finalize() is exposed separately. */
        mio_xdmf_series_free(series);
        R_ClearExternalPtr(x);
    }
}

static SEXP wrap_series(mio_xdmf_series *series) {
    SEXP ptr = PROTECT(R_MakeExternalPtr(series, mio_r_series_tag, R_NilValue));
    R_RegisterCFinalizerEx(ptr, mio_r_series_finalizer, TRUE);
    SEXP cls = PROTECT(Rf_mkString("mio_xdmf_series"));
    Rf_setAttrib(ptr, R_ClassSymbol, cls);
    UNPROTECT(2); /* cls, ptr */
    return ptr;
}

static mio_xdmf_series *series_of(SEXP x) {
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != mio_r_series_tag) {
        Rf_error("expected a mio_xdmf_series object");
    }
    mio_xdmf_series *series = (mio_xdmf_series *)R_ExternalPtrAddr(x);
    if (series == NULL) {
        Rf_error("this mio_xdmf_series has already been released");
    }
    return series;
}

SEXP R_mio_xdmf_series_create(SEXP path, SEXP data_format, SEXP gzip_level, SEXP mode,
                              SEXP auto_flush) {
    const char *fmt = mio_r_opt_string(data_format); /* NULL means "HDF" */
    const char *md = mio_r_string(mode, "mode");
    mio_xdmf_series_opts opts;
    mio_xdmf_series_opts_init(&opts);
    opts.data_format = fmt;
    opts.gzip_level = (int32_t)mio_r_int(gzip_level, "gzip_level");
    if (strcmp(md, "append") == 0) {
        opts.mode = MIO_XDMF_SERIES_APPEND;
    } else if (strcmp(md, "truncate") != 0) {
        Rf_error("mode must be \"truncate\" or \"append\"");
    }
    opts.auto_flush = mio_r_bool(auto_flush, "auto_flush");
    mio_xdmf_series *series = mio_xdmf_series_create_ex(mio_r_string(path, "path"), &opts);
    if (series == NULL) mio_r_fail("xdmf_series");
    return wrap_series(series);
}

SEXP R_mio_xdmf_series_write_points_cells(SEXP series, SEXP mesh) {
    mio_r_check(mio_xdmf_series_write_points_cells(series_of(series), mio_r_mesh(mesh)),
                "xdmf_series_write_points_cells");
    return R_NilValue;
}

SEXP R_mio_xdmf_series_write_data(SEXP series, SEXP time, SEXP mesh) {
    mio_r_check(mio_xdmf_series_write_data(series_of(series), mio_r_double(time, "time"),
                                           mio_r_mesh(mesh)),
                "xdmf_series_write_data");
    return R_NilValue;
}

SEXP R_mio_xdmf_series_finalize(SEXP series) {
    mio_r_check(mio_xdmf_series_finalize(series_of(series)), "xdmf_series_finalize");
    return R_NilValue;
}

SEXP R_mio_xdmf_series_flush(SEXP series) {
    mio_r_check(mio_xdmf_series_flush(series_of(series)), "xdmf_series_flush");
    return R_NilValue;
}

SEXP R_mio_xdmf_series_finalized(SEXP series) {
    int32_t f = mio_xdmf_series_finalized(series_of(series));
    if (f < 0) mio_r_fail("xdmf_series_finalized");
    return Rf_ScalarLogical(f == 1);
}

SEXP R_mio_xdmf_series_num_steps(SEXP series) {
    int64_t n = mio_xdmf_series_num_steps(series_of(series));
    if (n < 0) mio_r_fail("xdmf_series_num_steps");
    /* R has no native int64; a step count is far inside double's exact range. */
    return Rf_ScalarReal((double)n);
}

SEXP R_mio_xdmf_series_release(SEXP x) {
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != mio_r_series_tag) {
        Rf_error("expected a mio_xdmf_series object");
    }
    mio_xdmf_series *series = (mio_xdmf_series *)R_ExternalPtrAddr(x);
    if (series != NULL) {
        mio_xdmf_series_free(series);
        R_ClearExternalPtr(x);
    }
    return R_NilValue;
}

SEXP R_mio_xdmf_series_is_open(SEXP x) {
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != mio_r_series_tag) {
        return Rf_ScalarLogical(FALSE);
    }
    return Rf_ScalarLogical(R_ExternalPtrAddr(x) != NULL);
}
