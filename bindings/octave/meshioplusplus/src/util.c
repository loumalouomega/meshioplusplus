/* Handle registry, error translation and the small marshalling helpers. */

#include "mio_oct.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

/* --- handles ------------------------------------------------------------- *
 *
 * A handle object holds a uint64 id = (generation << 32) | slot. The slot
 * indexes a table the MEX owns; the generation is bumped every time a slot is
 * freed, so a stale id (a released handle, or an object loaded back from a
 * .mat file into a fresh session) never matches a live entry -- it is rejected
 * rather than dereferenced, the counterpart of the R binding's
 * external-pointer tags. A forged id that happens to match a live slot of the
 * right kind refers to a live, valid object, so it cannot crash either. */

typedef struct {
    void *ptr;
    uint32_t generation;
    int kind; /* -1 when free */
} mio_o_slot;

static mio_o_slot *g_slots = NULL;
static uint32_t g_num_slots = 0;
static uint32_t g_live = 0;

static void free_exodus(void *p) { mio_exodus_series_free((mio_exodus_series *)p); }
static void free_mesh(void *p) { mio_mesh_free((mio_mesh *)p); }
static void free_info(void *p) { mio_format_info_free((mio_format_info *)p); }
static void free_series(void *p) { mio_xdmf_series_free((mio_xdmf_series *)p); }
static void free_sequence(void *p) { mio_sequence_free((mio_sequence *)p); }

typedef void (*mio_o_free_fn)(void *);

static const struct {
    const char *class_name;
    mio_o_free_fn free_fn;
} g_kinds[MIO_O_NUM_KINDS] = {
    {"mio.MeshHandle", free_mesh},
    {"mio.FormatInfo", free_info},
    {"mio.XdmfSeries", free_series},
    {"mio.SequenceHandle", free_sequence},
    {"mio.ExodusSeries", free_exodus},
};

static uint32_t registry_add(int kind, void *ptr) {
    for (uint32_t i = 0; i < g_num_slots; ++i) {
        if (g_slots[i].kind < 0) {
            g_slots[i].kind = kind;
            g_slots[i].ptr = ptr;
            ++g_live;
            return i;
        }
    }
    uint32_t n = g_num_slots == 0 ? 64 : g_num_slots * 2;
    /* Plain realloc, NOT mxRealloc: the table outlives every MEX call. */
    mio_o_slot *grown = (mio_o_slot *)realloc(g_slots, (size_t)n * sizeof(mio_o_slot));
    if (grown == NULL) {
        g_kinds[kind].free_fn(ptr);
        mexErrMsgIdAndTxt("meshioplusplus:handle", "meshio++: out of memory for handles");
    }
    for (uint32_t i = g_num_slots; i < n; ++i) {
        grown[i].ptr = NULL;
        grown[i].generation = 1;
        grown[i].kind = -1;
    }
    g_slots = grown;
    uint32_t slot = g_num_slots;
    g_num_slots = n;
    g_slots[slot].kind = kind;
    g_slots[slot].ptr = ptr;
    ++g_live;
    return slot;
}

static void registry_free_slot(uint32_t slot) {
    void *ptr = g_slots[slot].ptr;
    int kind = g_slots[slot].kind;
    g_slots[slot].ptr = NULL;
    g_slots[slot].kind = -1;
    ++g_slots[slot].generation;
    --g_live;
    if (ptr != NULL) g_kinds[kind].free_fn(ptr);
    if (g_live == 0 && mexIsLocked()) mexUnlock();
}

/* The id stored in a handle object, or 0 when `obj` is not one of `kind`. */
static uint64_t handle_id(const mxArray *obj, mio_o_kind kind) {
    if (obj == NULL || !mxIsClass(obj, g_kinds[kind].class_name)) return 0;
    mxArray *p = mxGetProperty(obj, 0, "id");
    if (p == NULL) return 0;
    uint64_t id = 0;
    if (mxGetClassID(p) == mxUINT64_CLASS && mxGetNumberOfElements(p) == 1) {
        id = *(const uint64_t *)mxGetData(p);
    }
    mxDestroyArray(p);
    return id;
}

/* The live slot for `id` of `kind`, or -1. */
static int64_t live_slot(uint64_t id, mio_o_kind kind) {
    if (id == 0) return -1;
    uint32_t slot = (uint32_t)(id & 0xffffffffu);
    uint32_t gen = (uint32_t)(id >> 32);
    if (slot >= g_num_slots) return -1;
    if (g_slots[slot].kind != (int)kind || g_slots[slot].generation != gen) return -1;
    return (int64_t)slot;
}

mxArray *mio_o_wrap(mio_o_kind kind, void *ptr) {
    if (ptr == NULL) mio_o_fail("wrap");
    uint32_t slot = registry_add((int)kind, ptr);
    if (!mexIsLocked()) mexLock();
    uint64_t id = ((uint64_t)g_slots[slot].generation << 32) | slot;
    mxArray *arg = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *(uint64_t *)mxGetData(arg) = id;
    mxArray *out = NULL;
    /* mexCallMATLABWithTrap, so a failing constructor frees the native object
     * here instead of leaking it behind a raised error. */
    mxArray *err = mexCallMATLABWithTrap(1, &out, 1, &arg, g_kinds[kind].class_name);
    mxDestroyArray(arg);
    if (err != NULL) {
        registry_free_slot(slot);
        mexErrMsgIdAndTxt("meshioplusplus:handle", "meshio++: could not construct a %s object",
                          g_kinds[kind].class_name);
    }
    return out;
}

void *mio_o_unwrap(const mxArray *obj, mio_o_kind kind, const char *what) {
    const char *cls = g_kinds[kind].class_name;
    if (obj == NULL || !mxIsClass(obj, cls)) {
        mexErrMsgIdAndTxt("meshioplusplus:handle", "`%s` must be a %s object", what, cls);
    }
    int64_t slot = live_slot(handle_id(obj, kind), kind);
    if (slot < 0) {
        mexErrMsgIdAndTxt("meshioplusplus:released",
                          "this %s has already been released (or belongs to another session)",
                          cls);
    }
    return g_slots[slot].ptr;
}

int mio_o_release(const mxArray *obj, mio_o_kind kind) {
    const char *cls = g_kinds[kind].class_name;
    if (obj == NULL || !mxIsClass(obj, cls)) {
        mexErrMsgIdAndTxt("meshioplusplus:handle", "expected a %s object", cls);
    }
    int64_t slot = live_slot(handle_id(obj, kind), kind);
    if (slot < 0) return 0;
    registry_free_slot((uint32_t)slot);
    return 1;
}

int mio_o_is_open(const mxArray *obj, mio_o_kind kind) {
    return live_slot(handle_id(obj, kind), kind) >= 0;
}

void mio_o_release_all(void) {
    for (uint32_t i = 0; i < g_num_slots; ++i) {
        if (g_slots[i].kind >= 0) registry_free_slot(i);
    }
    free(g_slots);
    g_slots = NULL;
    g_num_slots = 0;
}

int mio_o_kind_of(const mxArray *obj) {
    if (obj == NULL) return -1;
    for (int k = 0; k < MIO_O_NUM_KINDS; ++k) {
        if (mxIsClass(obj, g_kinds[k].class_name)) return k;
    }
    return -1;
}

mxArray *mio_o_wrap_mesh(mio_mesh *mesh) { return mio_o_wrap(MIO_O_MESH, mesh); }

mio_mesh *mio_o_mesh(const mxArray *obj) {
    return (mio_mesh *)mio_o_unwrap(obj, MIO_O_MESH, "mesh");
}

/* --- errors -------------------------------------------------------------- */

static void raise_with(const char *what, const char *msg) {
    char id[128];
    /* A MATLAB message identifier is component:component, each a valid
     * identifier; anything else in `what` becomes an underscore. */
    size_t at = (size_t)snprintf(id, sizeof(id), "meshioplusplus:");
    for (const char *c = what; *c != '\0' && at + 1 < sizeof(id); ++c) {
        char ch = *c;
        int ok = (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') ||
                 (ch >= '0' && ch <= '9' && at > 15) || ch == '_';
        id[at++] = ok ? ch : '_';
    }
    id[at] = '\0';
    if (at == 15) snprintf(id + at, sizeof(id) - at, "error");
    mexErrMsgIdAndTxt(id, "%s", msg);
}

void mio_o_fail(const char *what) {
    const char *msg = mio_last_error();
    if (msg != NULL && msg[0] != '\0') {
        raise_with(what, msg);
    }
    char buf[256];
    snprintf(buf, sizeof(buf), "meshio++: %s failed", what);
    raise_with(what, buf);
}

void mio_o_error(const char *what, const char *fmt, ...) {
    char buf[1024];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    raise_with(what, buf);
}

void mio_o_check(mio_status status, const char *what) {
    if (status != MIO_OK) mio_o_fail(what);
}

void *mio_o_check_ptr(void *ptr, const char *what) {
    if (ptr == NULL) mio_o_fail(what);
    return ptr;
}

/* --- argument helpers ---------------------------------------------------- */

int mio_o_is_null(const mxArray *a) {
    return a == NULL || (mxIsDouble(a) && mxIsEmpty(a));
}

static char *to_cstring(const mxArray *a) {
#ifndef HAVE_OCTAVE
    /* MATLAB: mxArrayToString uses the user's locale; paths and names cross as
     * UTF-8 like every other binding. Octave's char arrays already are. */
    return mxArrayToUTF8String(a);
#else
    return mxArrayToString(a);
#endif
}

const char *mio_o_opt_string(const mxArray *a) {
    if (a == NULL || !mxIsChar(a)) return NULL;
    if (mxIsEmpty(a)) {
        char *s = (char *)mxCalloc(1, 1);
        return s;
    }
    char *s = to_cstring(a);
    if (s == NULL) mio_o_error("string", "meshio++: could not convert a string argument");
    return s;
}

const char *mio_o_string(const mxArray *a, const char *what) {
    const char *v = mio_o_opt_string(a);
    if (v == NULL) mio_o_error("argument", "`%s` must be a string", what);
    return v;
}

static double scalar_value(const mxArray *a, const char *what, const char *kind) {
    if (a == NULL || mxGetNumberOfElements(a) < 1 || !(mxIsNumeric(a) || mxIsLogical(a)) ||
        mxIsComplex(a)) {
        mio_o_error("argument", "`%s` must be a single %s", what, kind);
    }
    return mxGetScalar(a);
}

int mio_o_int(const mxArray *a, const char *what) {
    double v = scalar_value(a, what, "number");
    if (isnan(v)) mio_o_error("argument", "`%s` must not be NaN", what);
    return (int)v;
}

double mio_o_double(const mxArray *a, const char *what) {
    return scalar_value(a, what, "number");
}

int mio_o_bool(const mxArray *a, const char *what) {
    double v = scalar_value(a, what, "logical (true or false)");
    if (isnan(v)) mio_o_error("argument", "`%s` must not be NaN", what);
    return v != 0.0 ? 1 : 0;
}

int64_t mio_o_int64(const mxArray *a, const char *what) {
    if (a != NULL && mxGetNumberOfElements(a) == 1 && mxGetClassID(a) == mxINT64_CLASS) {
        return *(const int64_t *)mxGetData(a);
    }
    double v = scalar_value(a, what, "number");
    if (isnan(v)) mio_o_error("argument", "`%s` must not be NaN", what);
    return (int64_t)v;
}

/* Element `i` of any real numeric / logical array, as double. */
static double element(const mxArray *a, size_t i) {
    const void *d = mxGetData(a);
    switch (mxGetClassID(a)) {
    case mxDOUBLE_CLASS: return ((const double *)d)[i];
    case mxSINGLE_CLASS: return (double)((const float *)d)[i];
    case mxINT8_CLASS: return (double)((const int8_t *)d)[i];
    case mxUINT8_CLASS: return (double)((const uint8_t *)d)[i];
    case mxINT16_CLASS: return (double)((const int16_t *)d)[i];
    case mxUINT16_CLASS: return (double)((const uint16_t *)d)[i];
    case mxINT32_CLASS: return (double)((const int32_t *)d)[i];
    case mxUINT32_CLASS: return (double)((const uint32_t *)d)[i];
    case mxINT64_CLASS: return (double)((const int64_t *)d)[i];
    case mxUINT64_CLASS: return (double)((const uint64_t *)d)[i];
    case mxLOGICAL_CLASS: return ((const mxLogical *)d)[i] ? 1.0 : 0.0;
    default: return NAN;
    }
}

static void require_numeric(const mxArray *a, const char *what) {
    if (a == NULL || !(mxIsNumeric(a) || mxIsLogical(a)) || mxIsComplex(a) || mxIsSparse(a)) {
        mio_o_error("argument", "`%s` must be a real numeric array", what);
    }
}

void mio_o_vec3(const mxArray *a, double *out, const char *what) {
    mio_o_vecn(a, out, 3, what);
}

void mio_o_vecn(const mxArray *a, double *out, size_t n, const char *what) {
    require_numeric(a, what);
    size_t got = mxGetNumberOfElements(a);
    if (got != n) {
        mio_o_error("argument", "`%s` must have exactly %d elements, got %d", what, (int)n,
                    (int)got);
    }
    for (size_t i = 0; i < n; ++i) out[i] = element(a, i);
}

double *mio_o_doubles(const mxArray *a, size_t *n, const char *what) {
    require_numeric(a, what);
    size_t k = mxGetNumberOfElements(a);
    double *out = (double *)mxMalloc((k > 0 ? k : 1) * sizeof(double));
    for (size_t i = 0; i < k; ++i) out[i] = element(a, i);
    *n = k;
    return out;
}

const char *const *mio_o_names(const mxArray *a, int64_t *count, const char *what) {
    *count = 0;
    if (a == NULL || mio_o_is_null(a)) return NULL;
    if (mxIsChar(a)) {
        const char **ptrs = (const char **)mxMalloc(sizeof(char *));
        ptrs[0] = mio_o_string(a, what);
        *count = 1;
        return ptrs;
    }
    if (!mxIsCell(a)) mio_o_error("argument", "`%s` must be a cell array of strings", what);
    size_t n = mxGetNumberOfElements(a);
    if (n == 0) return NULL;
    const char **ptrs = (const char **)mxMalloc(n * sizeof(char *));
    for (size_t i = 0; i < n; ++i) {
        const mxArray *e = mxGetCell(a, (mwIndex)i);
        if (e == NULL || !mxIsChar(e)) {
            mio_o_error("argument", "`%s` must be a cell array of strings (element %d is not)",
                        what, (int)(i + 1));
        }
        ptrs[i] = mio_o_string(e, what);
    }
    *count = (int64_t)n;
    return ptrs;
}

/* (rows, cols) of a 2-D argument; a vector is (n, 1). */
static void matrix_dims(const mxArray *a, int64_t *rows, int64_t *cols, const char *what) {
    if (mxGetNumberOfDimensions(a) > 2) {
        mio_o_error("argument", "`%s` must be a 2-D matrix", what);
    }
    *rows = (int64_t)mxGetM(a);
    *cols = (int64_t)mxGetN(a);
}

double *mio_o_rows_to_c(const mxArray *a, int64_t *rows, int64_t *cols, double offset,
                        const char *what) {
    require_numeric(a, what);
    int64_t n, k;
    matrix_dims(a, &n, &k, what);
    double *out = (double *)mxMalloc((size_t)(n * k > 0 ? n * k : 1) * sizeof(double));
    for (int64_t i = 0; i < n; ++i)
        for (int64_t j = 0; j < k; ++j) out[i * k + j] = element(a, (size_t)(j * n + i)) - offset;
    *rows = n;
    *cols = k;
    return out;
}

static int64_t checked_index(double v, const char *what) {
    if (!isfinite(v) || v < 1.0 || v != floor(v)) {
        mio_o_error("argument", "`%s` is 1-based here; got %g", what, v);
    }
    return (int64_t)v - 1;
}

int64_t *mio_o_index_rows(const mxArray *a, int64_t *rows, int64_t *cols, const char *what) {
    require_numeric(a, what);
    int64_t n, k;
    matrix_dims(a, &n, &k, what);
    int64_t *out = (int64_t *)mxMalloc((size_t)(n * k > 0 ? n * k : 1) * sizeof(int64_t));
    for (int64_t i = 0; i < n; ++i)
        for (int64_t j = 0; j < k; ++j)
            out[i * k + j] = checked_index(element(a, (size_t)(j * n + i)), what);
    *rows = n;
    *cols = k;
    return out;
}

int64_t *mio_o_index_vector(const mxArray *a, int64_t *n, const char *what) {
    require_numeric(a, what);
    size_t k = mxGetNumberOfElements(a);
    int64_t *out = (int64_t *)mxMalloc((k > 0 ? k : 1) * sizeof(int64_t));
    for (size_t i = 0; i < k; ++i) out[i] = checked_index(element(a, i), what);
    *n = (int64_t)k;
    return out;
}

static int dtype_of_class(mxClassID c, mio_dtype *out) {
    switch (c) {
    case mxDOUBLE_CLASS: *out = MIO_FLOAT64; return 1;
    case mxSINGLE_CLASS: *out = MIO_FLOAT32; return 1;
    case mxINT8_CLASS: *out = MIO_INT8; return 1;
    case mxUINT8_CLASS: *out = MIO_UINT8; return 1;
    case mxINT16_CLASS: *out = MIO_INT16; return 1;
    case mxUINT16_CLASS: *out = MIO_UINT16; return 1;
    case mxINT32_CLASS: *out = MIO_INT32; return 1;
    case mxUINT32_CLASS: *out = MIO_UINT32; return 1;
    case mxINT64_CLASS: *out = MIO_INT64; return 1;
    case mxUINT64_CLASS: *out = MIO_UINT64; return 1;
    case mxLOGICAL_CLASS: *out = MIO_UINT8; return 1;
    default: return 0;
    }
}

static size_t dtype_size(mio_dtype dt) {
    switch (dt) {
    case MIO_INT8:
    case MIO_UINT8: return 1;
    case MIO_INT16:
    case MIO_UINT16: return 2;
    case MIO_FLOAT32:
    case MIO_INT32:
    case MIO_UINT32: return 4;
    default: return 8;
    }
}

/* Copy a column-major array with MATLAB dims `dims` (ndim of them) into a
 * row-major buffer of the same dims, element size `es`. Iterates the
 * destination in row-major order and computes the column-major source index. */
static void transpose_copy(void *dst, const void *src, size_t es, int ndim, const int64_t *dims,
                           int to_row_major) {
    size_t total = 1;
    for (int i = 0; i < ndim; ++i) total *= (size_t)dims[i];
    if (total == 0) return;
    if (ndim <= 1) {
        memcpy(dst, src, total * es);
        return;
    }
    int64_t idx[MIO_MAX_NDIM] = {0};
    size_t cm_stride[MIO_MAX_NDIM];
    cm_stride[0] = 1;
    for (int i = 1; i < ndim; ++i) cm_stride[i] = cm_stride[i - 1] * (size_t)dims[i - 1];
    for (size_t r = 0; r < total; ++r) {
        size_t c = 0;
        for (int i = 0; i < ndim; ++i) c += (size_t)idx[i] * cm_stride[i];
        if (to_row_major) {
            memcpy((char *)dst + r * es, (const char *)src + c * es, es);
        } else {
            memcpy((char *)dst + c * es, (const char *)src + r * es, es);
        }
        /* advance the row-major multi-index (last dim fastest) */
        for (int i = ndim - 1; i >= 0; --i) {
            if (++idx[i] < dims[i]) break;
            idx[i] = 0;
        }
    }
}

void *mio_o_data_to_c(const mxArray *a, mio_dtype *dtype, int32_t *ndim, int64_t *shape,
                      const char *what) {
    require_numeric(a, what);
    if (!dtype_of_class(mxGetClassID(a), dtype)) {
        mio_o_error("argument", "`%s` has an unsupported class", what);
    }
    mwSize nd = mxGetNumberOfDimensions(a);
    const mwSize *dims = mxGetDimensions(a);
    if (nd > MIO_MAX_NDIM) {
        mio_o_error("argument", "arrays of rank > %d are not supported", MIO_MAX_NDIM);
    }
    int64_t d[MIO_MAX_NDIM] = {0};
    int32_t n = (int32_t)nd;
    for (int32_t i = 0; i < n; ++i) d[i] = (int64_t)dims[i];
    /* A vector -- row or column -- is a 1-D array of its length. Trailing
     * singleton dims beyond 2 are already dropped by MATLAB. */
    if (n == 2 && (d[1] == 1 || d[0] == 1)) {
        d[0] = d[0] * d[1];
        n = 1;
    }
    for (int32_t i = 0; i < n; ++i) shape[i] = d[i];
    *ndim = n;
    size_t es = dtype_size(*dtype);
    size_t total = mxGetNumberOfElements(a);
    void *out = mxMalloc((total > 0 ? total : 1) * es);
    if (mxIsLogical(a)) {
        const mxLogical *src = (const mxLogical *)mxGetData(a);
        uint8_t *tmp = (uint8_t *)mxMalloc(total > 0 ? total : 1);
        for (size_t i = 0; i < total; ++i) tmp[i] = src[i] ? 1 : 0;
        transpose_copy(out, tmp, 1, n, d, 1);
    } else {
        transpose_copy(out, mxGetData(a), es, n, d, 1);
    }
    return out;
}

/* --- result helpers ------------------------------------------------------ */

mxArray *mio_o_str(const char *s) { return mxCreateString(s != NULL ? s : ""); }
mxArray *mio_o_scalar(double v) { return mxCreateDoubleScalar(v); }
mxArray *mio_o_logical(int v) { return mxCreateLogicalScalar(v ? true : false); }
mxArray *mio_o_empty(void) { return mxCreateDoubleMatrix(0, 0, mxREAL); }
mxArray *mio_o_cellstr(int64_t n) { return mxCreateCellMatrix((mwSize)n, 1); }
mxArray *mio_o_cell(int64_t n) { return mxCreateCellMatrix((mwSize)n, 1); }
mxArray *mio_o_column(int64_t n) { return mxCreateDoubleMatrix((mwSize)(n > 0 ? n : 0), 1, mxREAL); }

mxArray *mio_o_matrix(int64_t rows, int64_t cols) {
    return mxCreateDoubleMatrix((mwSize)(rows > 0 ? rows : 0), (mwSize)(cols > 0 ? cols : 0),
                                mxREAL);
}

mxArray *mio_o_logical_column(int64_t n) {
    return mxCreateLogicalMatrix((mwSize)(n > 0 ? n : 0), 1);
}

char *mio_o_getcstring(mio_o_str_getter getter, void *ctx, const char *what) {
    int64_t n = getter(ctx, NULL, 0);
    if (n < 0) mio_o_fail(what);
    char *buf = (char *)mxCalloc((size_t)n + 1, 1);
    if (n > 0) {
        if (getter(ctx, buf, n + 1) < 0) mio_o_fail(what);
        buf[n] = '\0';
    }
    return buf;
}

mxArray *mio_o_getstring(mio_o_str_getter getter, void *ctx, const char *what) {
    char *buf = mio_o_getcstring(getter, ctx, what);
    mxArray *out = mxCreateString(buf);
    mxFree(buf);
    return out;
}

void mio_o_fill_real(double *dst, const void *data, mio_dtype dtype, size_t n, double offset) {
    if (n == 0) return;
    if (dtype == MIO_FLOAT64 && offset == 0.0) {
        memcpy(dst, data, n * sizeof(double));
        return;
    }
#define MIO_O_FILL(T)                                                                              \
    do {                                                                                           \
        const T *p = (const T *)data;                                                              \
        for (size_t i = 0; i < n; ++i) dst[i] = (double)p[i] + offset;                             \
    } while (0)
    switch (dtype) {
    case MIO_FLOAT32: MIO_O_FILL(float); break;
    case MIO_FLOAT64: MIO_O_FILL(double); break;
    case MIO_INT8: MIO_O_FILL(int8_t); break;
    case MIO_INT16: MIO_O_FILL(int16_t); break;
    case MIO_INT32: MIO_O_FILL(int32_t); break;
    case MIO_INT64: MIO_O_FILL(int64_t); break;
    case MIO_UINT8: MIO_O_FILL(uint8_t); break;
    case MIO_UINT16: MIO_O_FILL(uint16_t); break;
    case MIO_UINT32: MIO_O_FILL(uint32_t); break;
    case MIO_UINT64: MIO_O_FILL(uint64_t); break;
    default: mio_o_error("dtype", "meshio++: unknown dtype %d", (int)dtype);
    }
#undef MIO_O_FILL
}

mxArray *mio_o_copy_as_real(const void *data, mio_dtype dtype, int64_t n) {
    mxArray *out = mio_o_column(n);
    if (n > 0) mio_o_fill_real(mxGetPr(out), data, dtype, (size_t)n, 0.0);
    return out;
}

mxArray *mio_o_rows_from_c(const void *data, mio_dtype dtype, int64_t rows, int64_t cols,
                           double offset) {
    mxArray *out = mio_o_matrix(rows, cols);
    if (rows <= 0 || cols <= 0) return out;
    size_t total = (size_t)(rows * cols);
    double *tmp = (double *)mxMalloc(total * sizeof(double));
    mio_o_fill_real(tmp, data, dtype, total, offset);
    double *dst = mxGetPr(out);
    for (int64_t i = 0; i < rows; ++i)
        for (int64_t j = 0; j < cols; ++j) dst[j * rows + i] = tmp[i * cols + j];
    mxFree(tmp);
    return out;
}

static mxClassID class_of_dtype(mio_dtype dt) {
    switch (dt) {
    case MIO_FLOAT32: return mxSINGLE_CLASS;
    case MIO_FLOAT64: return mxDOUBLE_CLASS;
    case MIO_INT8: return mxINT8_CLASS;
    case MIO_INT16: return mxINT16_CLASS;
    case MIO_INT32: return mxINT32_CLASS;
    case MIO_INT64: return mxINT64_CLASS;
    case MIO_UINT8: return mxUINT8_CLASS;
    case MIO_UINT16: return mxUINT16_CLASS;
    case MIO_UINT32: return mxUINT32_CLASS;
    case MIO_UINT64: return mxUINT64_CLASS;
    default: mio_o_error("dtype", "meshio++: unknown dtype %d", (int)dt); return mxUNKNOWN_CLASS;
    }
}

mxArray *mio_o_typed_from_c(const void *data, mio_dtype dtype, int32_t ndim,
                            const int64_t *shape) {
    mwSize dims[MIO_MAX_NDIM > 2 ? MIO_MAX_NDIM : 2];
    int64_t d[MIO_MAX_NDIM > 2 ? MIO_MAX_NDIM : 2];
    int nd = ndim;
    if (nd <= 0) {
        /* rank-0 scalar: one element */
        nd = 2;
        d[0] = 1;
        d[1] = 1;
    } else if (nd == 1) {
        nd = 2;
        d[0] = shape[0];
        d[1] = 1;
    } else {
        for (int i = 0; i < nd; ++i) d[i] = shape[i];
    }
    for (int i = 0; i < nd; ++i) dims[i] = (mwSize)d[i];
    mxArray *out = mxCreateNumericArray((mwSize)nd, dims, class_of_dtype(dtype), mxREAL);
    size_t total = mxGetNumberOfElements(out);
    if (total == 0 || data == NULL) return out;
    transpose_copy(mxGetData(out), data, dtype_size(dtype), nd, d, 0);
    return out;
}

mxArray *mio_o_shift_map(const int64_t *data, int64_t n) {
    if (n < 0) n = 0;
    mxArray *out = mio_o_column(n);
    double *dst = mxGetPr(out);
    for (int64_t i = 0; i < n; ++i) dst[i] = data[i] < 0 ? 0.0 : (double)(data[i] + 1);
    return out;
}

mxArray *mio_o_int64_column(const int64_t *data, int64_t n) {
    if (n < 0) n = 0;
    mxArray *out = mio_o_column(n);
    double *dst = mxGetPr(out);
    for (int64_t i = 0; i < n; ++i) dst[i] = (double)data[i];
    return out;
}

mxArray *mio_o_struct(int n, const char **names, mxArray **values) {
    mxArray *out = mxCreateStructMatrix(1, 1, n, names);
    for (int i = 0; i < n; ++i) {
        mxSetFieldByNumber(out, 0, i, values[i] != NULL ? values[i] : mio_o_empty());
    }
    return out;
}

mxArray *mio_o_struct_array(int64_t n, int nfields, const char **names) {
    return mxCreateStructMatrix((mwSize)(n > 0 ? n : 0), 1, nfields, names);
}

const char *mio_o_dtype_name(mio_dtype dtype) {
    switch (dtype) {
    case MIO_FLOAT32: return "float32";
    case MIO_FLOAT64: return "float64";
    case MIO_INT8: return "int8";
    case MIO_INT16: return "int16";
    case MIO_INT32: return "int32";
    case MIO_INT64: return "int64";
    case MIO_UINT8: return "uint8";
    case MIO_UINT16: return "uint16";
    case MIO_UINT32: return "uint32";
    case MIO_UINT64: return "uint64";
    default: return "unknown";
    }
}
