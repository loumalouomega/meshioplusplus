/* Multi-file / transient datasets.
 *
 * A sequence is an ordered PLAN over a set of files (or the steps inside one
 * multi-step file), read one step at a time. It owns no mesh, deliberately:
 * caching what it handed out is exactly the accumulation the streaming
 * guarantee forbids, which is why `mio.sequence_read()` returns an OWNED mesh.
 *
 * Its own handle kind (mio.SequenceHandle), exactly as mio_mesh and
 * mio_xdmf_series have, so the three can never be passed for one another. A
 * foreign or already-released handle is an error, never a dereference.
 *
 * A port of the R binding's src/sequence.c; see mio_oct.h. */

#include "mio_oct.h"

#include <string.h>

static mio_sequence *oseq_of(const mxArray *x) {
    return (mio_sequence *)mio_o_unwrap(x, MIO_O_SEQUENCE, "seq");
}

/* "auto"/"file"/"filename"/"index" -> the C enum. */
static int32_t oseq_time_from(const mxArray *time_from) {
    const char *s = mio_o_string(time_from, "time_from");
    if (strcmp(s, "auto") == 0) return 0;
    if (strcmp(s, "file") == 0) return 1;
    if (strcmp(s, "filename") == 0) return 2;
    if (strcmp(s, "index") == 0) return 3;
    mio_o_error("time_from", "time_from must be \"auto\", \"file\", \"filename\" or \"index\"");
    return 0; /* unreachable */
}

/* Shared option parsing of the two openers: (format, times, time_from, sort). */
static void oseq_opts(const mxArray **a, mio_sequence_opts *opts) {
    mio_sequence_opts_init(opts);
    opts->format = mio_o_opt_string(a[0]);
    opts->time_from = oseq_time_from(a[2]);
    opts->sort = (int32_t)mio_o_bool(a[3], "sort");
    if (!mio_o_is_null(a[1])) {
        size_t n = 0;
        opts->times = mio_o_doubles(a[1], &n, "times");
        opts->num_times = (int64_t)n;
    }
}

/* sequence_open(pattern, format, times, time_from, sort) */
static mxArray *O_mio_sequence_open(const mxArray **a) {
    mio_sequence_opts opts;
    oseq_opts(a + 1, &opts);
    mio_sequence *seq = mio_sequence_open_ex(mio_o_string(a[0], "pattern"), &opts);
    if (seq == NULL) mio_o_fail("sequence_open");
    return mio_o_wrap(MIO_O_SEQUENCE, seq);
}

/* sequence_open_list(paths, format, times, time_from, sort) */
static mxArray *O_mio_sequence_open_list(const mxArray **a) {
    int64_t n = 0;
    const char *const *paths = mio_o_names(a[0], &n, "paths");
    if (paths == NULL || n == 0) {
        mio_o_error("sequence_open_list", "paths must be a non-empty cell array of strings");
    }
    mio_sequence_opts opts;
    oseq_opts(a + 1, &opts);
    mio_sequence *seq = mio_sequence_open_list(paths, n, &opts);
    if (seq == NULL) mio_o_fail("sequence_open_list");
    return mio_o_wrap(MIO_O_SEQUENCE, seq);
}

static mxArray *O_mio_sequence_count(const mxArray **a) {
    int64_t n = mio_sequence_count(oseq_of(a[0]));
    if (n < 0) mio_o_fail("sequence_count");
    return mio_o_scalar((double)n);
}

typedef struct {
    const mio_sequence *seq;
    int64_t index;
} oseq_path_ctx;

static int64_t oseq_path_getter(void *ctx, char *buf, int64_t buflen) {
    oseq_path_ctx *c = (oseq_path_ctx *)ctx;
    return mio_sequence_path(c->seq, c->index, buf, buflen);
}

static mxArray *O_mio_sequence_path(const mxArray **a) {
    oseq_path_ctx c = {oseq_of(a[0]), (int64_t)mio_o_int(a[1], "index") - 1}; /* 1-based */
    return mio_o_getstring(oseq_path_getter, &c, "sequence_path");
}

static mxArray *O_mio_sequence_step(const mxArray **a) {
    int64_t k = mio_sequence_step(oseq_of(a[0]), (int64_t)mio_o_int(a[1], "index") - 1);
    if (k < 0) mio_o_fail("sequence_step");
    return mio_o_scalar((double)k);
}

static mxArray *O_mio_sequence_time(const mxArray **a) {
    double t = 0.0;
    mio_o_check(mio_sequence_time(oseq_of(a[0]), (int64_t)mio_o_int(a[1], "index") - 1, &t),
                "sequence_time");
    return mio_o_scalar(t);
}

static mxArray *O_mio_sequence_time_source(const mxArray **a) {
    int32_t code = mio_sequence_time_source(oseq_of(a[0]), (int64_t)mio_o_int(a[1], "index") - 1);
    if (code < 0 || code > 3) mio_o_fail("sequence_time_source");
    static const char *names[] = {"explicit", "file", "filename", "index"};
    return mio_o_str(names[code]);
}

static mxArray *O_mio_sequence_read(const mxArray **a) {
    mio_mesh *m = mio_sequence_read(oseq_of(a[0]), (int64_t)mio_o_int(a[1], "index") - 1);
    if (m == NULL) mio_o_fail("sequence_read");
    return mio_o_wrap_mesh(m);
}

/* Idempotent, and -- as in R -- a no-op for anything that is not a sequence
 * handle. */
static mxArray *O_mio_sequence_free(const mxArray **a) {
    if (mio_o_kind_of(a[0]) == (int)MIO_O_SEQUENCE) mio_o_release(a[0], MIO_O_SEQUENCE);
    return NULL;
}

/* sequence_to_timeseries(seq, out_path, out_format, ascii) */
static mxArray *O_mio_sequence_to_timeseries(const mxArray **a) {
    /* The transient writer drives XdmfTimeSeriesWriter directly rather than
     * the registry, so `ascii` is the one write option with anywhere to go:
     * it selects XDMF's "XML" data format (no HDF5 needed) over the default
     * "HDF" -- the option a build without HDF5 support needs. */
    mio_sequence *seq = oseq_of(a[0]);
    const char *out_path = mio_o_string(a[1], "out_path");
    const char *out_format = mio_o_opt_string(a[2]);
    if (!mio_o_bool(a[3], "ascii")) {
        mio_o_check(mio_sequence_to_timeseries(seq, out_path, out_format),
                    "sequence_to_timeseries");
        return NULL;
    }
    mio_write_opts opts;
    mio_write_opts_init(&opts);
    /* The token, as the R binding and c_api.cpp spell it: the header's later
     * `#define MIO_ENCODING_ASCII` (a mio_text_opts value) shadows the
     * mio_write_encoding enumerator, and the C API compares against the same
     * expansion. */
    opts.encoding = MIO_ENCODING_ASCII;
    mio_o_check(mio_sequence_to_timeseries_ex(seq, out_path, out_format, &opts),
                "sequence_to_timeseries");
    return NULL;
}

/* sequence_resample(seq, out_path, out_format, times, method, clamp,
 *                   blend_points) */
static mxArray *O_mio_sequence_resample(const mxArray **a) {
    if (!mxIsNumeric(a[3])) mio_o_error("sequence_resample", "meshio++: times must be numeric");
    size_t n = 0;
    double *t = mio_o_doubles(a[3], &n, "times");
    mio_resample_opts opts;
    mio_resample_opts_init(&opts);
    opts.times = t;
    opts.num_times = (int64_t)n;
    opts.method = (int32_t)mio_o_int64(a[4], "method");
    opts.extrapolate = mio_o_bool(a[5], "clamp");
    opts.blend_points = mio_o_bool(a[6], "blend_points");
    mio_o_check(mio_sequence_resample(oseq_of(a[0]), mio_o_string(a[1], "out_path"),
                                      mio_o_opt_string(a[2]), &opts),
                "sequence_resample");
    return NULL;
}

/* timeseries_to_sequence(in_path, in_format, out_pattern, out_format) */
static mxArray *O_mio_timeseries_to_sequence(const mxArray **a) {
    mio_o_check(mio_timeseries_to_sequence(mio_o_string(a[0], "in_path"), mio_o_opt_string(a[1]),
                                           mio_o_string(a[2], "out_pattern"),
                                           mio_o_opt_string(a[3])),
                "timeseries_to_sequence");
    return NULL;
}

static mxArray *O_mio_sequence_pipeline_run_file(const mxArray **a) {
    mio_o_check(mio_sequence_pipeline_run_file(mio_o_string(a[0], "settings_path")),
                "sequence_pipeline_run_file");
    return NULL;
}

static mxArray *O_mio_sequence_pipeline_run_json(const mxArray **a) {
    mio_o_check(mio_sequence_pipeline_run_json(mio_o_string(a[0], "json_text")),
                "sequence_pipeline_run_json");
    return NULL;
}

/* --- table --------------------------------------------------------------- */

static const mio_o_entry kEntries[] = {
    MIO_O_ENTRY(mio_sequence_open, 5),
    MIO_O_ENTRY(mio_sequence_open_list, 5),
    MIO_O_ENTRY(mio_sequence_count, 1),
    MIO_O_ENTRY(mio_sequence_path, 2),
    MIO_O_ENTRY(mio_sequence_step, 2),
    MIO_O_ENTRY(mio_sequence_time, 2),
    MIO_O_ENTRY(mio_sequence_time_source, 2),
    MIO_O_ENTRY(mio_sequence_read, 2),
    MIO_O_ENTRY(mio_sequence_free, 1),
    MIO_O_ENTRY(mio_sequence_to_timeseries, 4),
    MIO_O_ENTRY(mio_sequence_resample, 7),
    MIO_O_ENTRY(mio_timeseries_to_sequence, 4),
    MIO_O_ENTRY(mio_sequence_pipeline_run_file, 1),
    MIO_O_ENTRY(mio_sequence_pipeline_run_json, 1),
};

const mio_o_entry *mio_o_table_sequence(size_t *n) {
    *n = sizeof(kEntries) / sizeof(kEntries[0]);
    return kEntries;
}
