/* Multi-file / transient datasets.
 *
 * A sequence is an ordered PLAN over a set of files (or the steps inside one
 * multi-step file), read one step at a time. It owns no mesh, deliberately:
 * caching what it handed out is exactly the accumulation the streaming
 * guarantee forbids, which is why `mio_sequence_read()` returns an OWNED mesh.
 *
 * Its own tag and finalizer, exactly as `mio_mesh` and `mio_xdmf_series` have,
 * so the three can never be passed for one another. A foreign or
 * already-released pointer is an R error, never a dereference. */

#include "mio_r.h"

SEXP mio_r_sequence_tag = NULL;

static void mio_r_sequence_finalizer(SEXP x) {
    mio_sequence *seq = (mio_sequence *)R_ExternalPtrAddr(x);
    if (seq != NULL) {
        mio_sequence_free(seq);
        R_ClearExternalPtr(x);
    }
}

static SEXP wrap_sequence(mio_sequence *seq) {
    SEXP ptr = PROTECT(R_MakeExternalPtr(seq, mio_r_sequence_tag, R_NilValue));
    R_RegisterCFinalizerEx(ptr, mio_r_sequence_finalizer, TRUE);
    SEXP cls = PROTECT(Rf_mkString("mio_sequence"));
    Rf_setAttrib(ptr, R_ClassSymbol, cls);
    UNPROTECT(2); /* cls, ptr */
    return ptr;
}

static mio_sequence *sequence_of(SEXP x) {
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != mio_r_sequence_tag) {
        Rf_error("expected a mio_sequence object");
    }
    mio_sequence *seq = (mio_sequence *)R_ExternalPtrAddr(x);
    if (seq == NULL) {
        Rf_error("this mio_sequence has already been released");
    }
    return seq;
}

/* "auto"/"file"/"filename"/"index" -> the C enum. */
static int32_t sequence_time_from(SEXP time_from) {
    const char *s = mio_r_string(time_from, "time_from");
    if (strcmp(s, "auto") == 0)
        return 0;
    if (strcmp(s, "file") == 0)
        return 1;
    if (strcmp(s, "filename") == 0)
        return 2;
    if (strcmp(s, "index") == 0)
        return 3;
    Rf_error("time_from must be \"auto\", \"file\", \"filename\" or \"index\"");
    return 0; /* unreachable */
}

SEXP R_mio_sequence_open(SEXP pattern, SEXP format, SEXP times, SEXP time_from, SEXP sort) {
    mio_sequence_opts opts;
    mio_sequence_opts_init(&opts);
    opts.format = mio_r_opt_string(format);
    opts.time_from = sequence_time_from(time_from);
    opts.sort = (int32_t)mio_r_bool(sort, "sort");

    double *tbuf = NULL;
    R_xlen_t ntimes = 0;
    if (times != R_NilValue) {
        SEXP nums = PROTECT(Rf_coerceVector(times, REALSXP));
        ntimes = Rf_xlength(nums);
        tbuf = REAL(nums);
        opts.times = tbuf;
        opts.num_times = (int64_t)ntimes;
        mio_sequence *seq = mio_sequence_open_ex(mio_r_string(pattern, "pattern"), &opts);
        UNPROTECT(1); /* nums */
        if (seq == NULL)
            mio_r_fail("sequence_open");
        return wrap_sequence(seq);
    }
    mio_sequence *seq = mio_sequence_open_ex(mio_r_string(pattern, "pattern"), &opts);
    if (seq == NULL)
        mio_r_fail("sequence_open");
    return wrap_sequence(seq);
}

SEXP R_mio_sequence_open_list(SEXP paths, SEXP format, SEXP times, SEXP time_from, SEXP sort) {
    if (TYPEOF(paths) != STRSXP || Rf_xlength(paths) == 0) {
        Rf_error("paths must be a non-empty character vector");
    }
    R_xlen_t n = Rf_xlength(paths);
    const char **cpaths = (const char **)R_alloc((size_t)n, sizeof(const char *));
    for (R_xlen_t i = 0; i < n; ++i) {
        cpaths[i] = CHAR(STRING_ELT(paths, i));
    }

    mio_sequence_opts opts;
    mio_sequence_opts_init(&opts);
    opts.format = mio_r_opt_string(format);
    opts.time_from = sequence_time_from(time_from);
    opts.sort = (int32_t)mio_r_bool(sort, "sort");

    if (times != R_NilValue) {
        SEXP nums = PROTECT(Rf_coerceVector(times, REALSXP));
        opts.times = REAL(nums);
        opts.num_times = (int64_t)Rf_xlength(nums);
        mio_sequence *seq = mio_sequence_open_list(cpaths, (int64_t)n, &opts);
        UNPROTECT(1); /* nums */
        if (seq == NULL)
            mio_r_fail("sequence_open_list");
        return wrap_sequence(seq);
    }
    mio_sequence *seq = mio_sequence_open_list(cpaths, (int64_t)n, &opts);
    if (seq == NULL)
        mio_r_fail("sequence_open_list");
    return wrap_sequence(seq);
}

SEXP R_mio_sequence_count(SEXP seq) {
    int64_t n = mio_sequence_count(sequence_of(seq));
    if (n < 0)
        mio_r_fail("sequence_count");
    return Rf_ScalarReal((double)n);
}

SEXP R_mio_sequence_path(SEXP seq, SEXP index) {
    mio_sequence *s = sequence_of(seq);
    int64_t i = (int64_t)mio_r_int(index, "index") - 1; /* R is 1-based */
    int64_t need = mio_sequence_path(s, i, NULL, 0);
    if (need < 0)
        mio_r_fail("sequence_path");
    char *buf = (char *)R_alloc((size_t)need + 1, 1);
    if (mio_sequence_path(s, i, buf, need + 1) < 0)
        mio_r_fail("sequence_path");
    return Rf_mkString(buf);
}

SEXP R_mio_sequence_step(SEXP seq, SEXP index) {
    int64_t k = mio_sequence_step(sequence_of(seq), (int64_t)mio_r_int(index, "index") - 1);
    if (k < 0)
        mio_r_fail("sequence_step");
    return Rf_ScalarReal((double)k);
}

SEXP R_mio_sequence_time(SEXP seq, SEXP index) {
    double t = 0.0;
    mio_r_check(mio_sequence_time(sequence_of(seq), (int64_t)mio_r_int(index, "index") - 1, &t),
                "sequence_time");
    return Rf_ScalarReal(t);
}

SEXP R_mio_sequence_time_source(SEXP seq, SEXP index) {
    int32_t code =
        mio_sequence_time_source(sequence_of(seq), (int64_t)mio_r_int(index, "index") - 1);
    if (code < 0)
        mio_r_fail("sequence_time_source");
    static const char *names[] = {"explicit", "file", "filename", "index"};
    return Rf_mkString(names[code]);
}

SEXP R_mio_sequence_read(SEXP seq, SEXP index) {
    mio_mesh *m = mio_sequence_read(sequence_of(seq), (int64_t)mio_r_int(index, "index") - 1);
    if (m == NULL)
        mio_r_fail("sequence_read");
    return mio_r_wrap_mesh(m);
}

SEXP R_mio_sequence_free(SEXP seq) {
    if (TYPEOF(seq) == EXTPTRSXP && R_ExternalPtrTag(seq) == mio_r_sequence_tag) {
        mio_sequence *s = (mio_sequence *)R_ExternalPtrAddr(seq);
        if (s != NULL) {
            mio_sequence_free(s);
            R_ClearExternalPtr(seq);
        }
    }
    return R_NilValue;
}

SEXP R_mio_sequence_to_timeseries(SEXP seq, SEXP out_path, SEXP out_format, SEXP ascii) {
    /* The transient writer drives XdmfTimeSeriesWriter directly rather than
     * the registry, so `ascii` is the one write option with anywhere to go:
     * it selects XDMF's "XML" data format (no HDF5 needed) over the default
     * "HDF" -- the option a build without HDF5 support (this package's own
     * notebook environment among them, see doc/r.md) needs. */
    if (!mio_r_bool(ascii, "ascii")) {
        mio_r_check(mio_sequence_to_timeseries(sequence_of(seq), mio_r_string(out_path, "out_path"),
                                               mio_r_opt_string(out_format)),
                    "sequence_to_timeseries");
        return R_NilValue;
    }
    mio_write_opts opts;
    mio_write_opts_init(&opts);
    opts.encoding = MIO_ENCODING_ASCII;
    mio_r_check(mio_sequence_to_timeseries_ex(sequence_of(seq), mio_r_string(out_path, "out_path"),
                                              mio_r_opt_string(out_format), &opts),
                "sequence_to_timeseries");
    return R_NilValue;
}

SEXP R_mio_timeseries_to_sequence(SEXP in_path, SEXP in_format, SEXP out_pattern, SEXP out_format) {
    mio_r_check(mio_timeseries_to_sequence(
                    mio_r_string(in_path, "in_path"), mio_r_opt_string(in_format),
                    mio_r_string(out_pattern, "out_pattern"), mio_r_opt_string(out_format)),
                "timeseries_to_sequence");
    return R_NilValue;
}

SEXP R_mio_sequence_pipeline_run_file(SEXP settings_path) {
    mio_r_check(mio_sequence_pipeline_run_file(mio_r_string(settings_path, "settings_path")),
                "sequence_pipeline_run_file");
    return R_NilValue;
}

SEXP R_mio_sequence_pipeline_run_json(SEXP json_text) {
    mio_r_check(mio_sequence_pipeline_run_json(mio_r_string(json_text, "json_text")),
                "sequence_pipeline_run_json");
    return R_NilValue;
}
