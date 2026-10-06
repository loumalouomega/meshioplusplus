/* Handle management, error translation and the small marshalling helpers. */

#include "mio_r.h"

#include <string.h>

SEXP mio_r_mesh_tag = NULL;

/* --- handles ------------------------------------------------------------ */

static void mio_r_mesh_finalizer(SEXP x) {
    mio_mesh *mesh = (mio_mesh *)R_ExternalPtrAddr(x);
    if (mesh != NULL) {
        mio_mesh_free(mesh);
        R_ClearExternalPtr(x);
    }
}

SEXP mio_r_wrap_mesh(mio_mesh *mesh) {
    SEXP ptr = PROTECT(R_MakeExternalPtr(mesh, mio_r_mesh_tag, R_NilValue));
    /* onexit = TRUE: run the finalizer at the end of the session too, not only
     * when the GC happens to collect. */
    R_RegisterCFinalizerEx(ptr, mio_r_mesh_finalizer, TRUE);
    SEXP cls = PROTECT(Rf_mkString("mio_mesh"));
    Rf_setAttrib(ptr, R_ClassSymbol, cls);
    UNPROTECT(2); /* cls, ptr */
    return ptr;
}

mio_mesh *mio_r_mesh(SEXP x) {
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != mio_r_mesh_tag) {
        Rf_error("expected a mio_mesh object");
    }
    mio_mesh *mesh = (mio_mesh *)R_ExternalPtrAddr(x);
    if (mesh == NULL) {
        Rf_error("this mio_mesh has already been released");
    }
    return mesh;
}

/* --- errors ------------------------------------------------------------- */

void mio_r_fail(const char *what) {
    const char *msg = mio_last_error();
    if (msg != NULL && msg[0] != '\0') {
        Rf_error("%s", msg);
    } else {
        Rf_error("meshio++: %s failed", what);
    }
}

void mio_r_check(mio_status status, const char *what) {
    if (status != MIO_OK) {
        mio_r_fail(what);
    }
}

void *mio_r_check_ptr(void *ptr, const char *what) {
    if (ptr == NULL) {
        mio_r_fail(what);
    }
    return ptr;
}

/* --- argument helpers --------------------------------------------------- */

const char *mio_r_opt_string(SEXP s) {
    if (s == R_NilValue || TYPEOF(s) != STRSXP || Rf_length(s) < 1) {
        return NULL;
    }
    SEXP e = STRING_ELT(s, 0);
    if (e == NA_STRING) {
        return NULL;
    }
    return CHAR(e);
}

const char *mio_r_string(SEXP s, const char *what) {
    const char *v = mio_r_opt_string(s);
    if (v == NULL) {
        Rf_error("`%s` must be a single non-NA string", what);
    }
    return v;
}

int mio_r_int(SEXP s, const char *what) {
    if (Rf_length(s) < 1) {
        Rf_error("`%s` must be a single number", what);
    }
    double v = Rf_asReal(s);
    if (ISNA(v) || ISNAN(v)) {
        Rf_error("`%s` must not be NA", what);
    }
    return (int)v;
}

double mio_r_double(SEXP s, const char *what) {
    if (Rf_length(s) < 1) {
        Rf_error("`%s` must be a single number", what);
    }
    return Rf_asReal(s);
}

int mio_r_bool(SEXP s, const char *what) {
    if (Rf_length(s) < 1) {
        Rf_error("`%s` must be TRUE or FALSE", what);
    }
    int v = Rf_asLogical(s);
    if (v == NA_LOGICAL) {
        Rf_error("`%s` must not be NA", what);
    }
    return v ? 1 : 0;
}

int64_t mio_r_int64(SEXP s, const char *what) {
    if (Rf_length(s) < 1) {
        Rf_error("`%s` must be a single number", what);
    }
    double v = Rf_asReal(s);
    if (ISNA(v) || ISNAN(v)) {
        Rf_error("`%s` must not be NA", what);
    }
    return (int64_t)v;
}

void mio_r_vec3(SEXP s, double *out, const char *what) {
    if (TYPEOF(s) != REALSXP && TYPEOF(s) != INTSXP) {
        Rf_error("`%s` must be numeric", what);
    }
    if (Rf_length(s) != 3) {
        Rf_error("`%s` must have exactly 3 elements, got %d", what, (int)Rf_length(s));
    }
    SEXP r = PROTECT(Rf_coerceVector(s, REALSXP));
    memcpy(out, REAL(r), 3 * sizeof(double));
    UNPROTECT(1);
}

const char *const *mio_r_names(SEXP s, int64_t *count, SEXP *shelter) {
    if (s == R_NilValue || TYPEOF(s) != STRSXP || Rf_length(s) == 0) {
        /* Always leave exactly one PROTECT on the stack, so every caller can
         * UNPROTECT(1) unconditionally. */
        *shelter = PROTECT(R_NilValue);
        *count = 0;
        return NULL;
    }
    R_xlen_t n = Rf_length(s);
    /* A raw vector holds the pointer array. It is left PROTECT-ed and handed
     * back through *shelter -- unprotecting it here and re-protecting in the
     * caller would leave a GC window in between. The CHARSXP bodies it points
     * at are owned by `s`, which the caller keeps alive for the same span. */
    SEXP buf = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)(n * sizeof(char *))));
    const char **ptrs = (const char **)RAW(buf);
    for (R_xlen_t i = 0; i < n; ++i) {
        SEXP e = STRING_ELT(s, i);
        ptrs[i] = (e == NA_STRING) ? "" : CHAR(e);
    }
    *shelter = buf; /* still protected; the caller UNPROTECTs */
    *count = (int64_t)n;
    return ptrs;
}

/* --- result helpers ----------------------------------------------------- */

SEXP mio_r_getstring(mio_r_str_getter getter, void *ctx, const char *what) {
    int64_t n = getter(ctx, NULL, 0);
    if (n < 0) {
        mio_r_fail(what);
    }
    if (n == 0) {
        return Rf_mkString("");
    }
    char *buf = (char *)R_alloc((size_t)n + 1, 1);
    int64_t m = getter(ctx, buf, n + 1);
    if (m < 0) {
        mio_r_fail(what);
    }
    buf[n] = '\0';
    return Rf_mkString(buf);
}

/* R has no native 64-bit integer type, so every numeric array the C API hands
 * back -- including the int64 connectivity, region entries and index maps --
 * comes into R as `double`. That is exact to 2^53, far beyond any mesh this
 * library will meet, and it avoids a hard dependency on the bit64 package for
 * a theoretical case. Documented as a named limitation. */
void mio_r_fill_real(double *dst, const void *data, mio_dtype dtype, R_xlen_t n, double offset) {
    if (dtype == MIO_FLOAT64 && offset == 0.0) {
        memcpy(dst, data, (size_t)n * sizeof(double));
        return;
    }
#define MIO_R_FILL(T)                                                                              \
    do {                                                                                           \
        const T *p = (const T *)data;                                                              \
        for (R_xlen_t i = 0; i < n; ++i) dst[i] = (double)p[i] + offset;                           \
    } while (0)
    switch (dtype) {
    case MIO_FLOAT32: MIO_R_FILL(float); break;
    case MIO_FLOAT64: MIO_R_FILL(double); break;
    case MIO_INT8: MIO_R_FILL(int8_t); break;
    case MIO_INT16: MIO_R_FILL(int16_t); break;
    case MIO_INT32: MIO_R_FILL(int32_t); break;
    case MIO_INT64: MIO_R_FILL(int64_t); break;
    case MIO_UINT8: MIO_R_FILL(uint8_t); break;
    case MIO_UINT16: MIO_R_FILL(uint16_t); break;
    case MIO_UINT32: MIO_R_FILL(uint32_t); break;
    case MIO_UINT64: MIO_R_FILL(uint64_t); break;
    default:
        Rf_error("meshio++: unknown dtype %d", (int)dtype);
    }
#undef MIO_R_FILL
}

SEXP mio_r_copy_as_real(const void *data, mio_dtype dtype, R_xlen_t n) {
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    mio_r_fill_real(REAL(out), data, dtype, n, 0.0);
    UNPROTECT(1);
    return out;
}

SEXP mio_r_shift_map(const int64_t *data, int64_t n) {
    if (n < 0) n = 0;
    SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)n));
    double *dst = REAL(out);
    for (int64_t i = 0; i < n; ++i) {
        dst[i] = data[i] < 0 ? 0.0 : (double)(data[i] + 1);
    }
    UNPROTECT(1);
    return out;
}

SEXP mio_r_named_list(int n, const char **names, SEXP *values) {
    SEXP out = PROTECT(Rf_allocVector(VECSXP, n));
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, n));
    for (int i = 0; i < n; ++i) {
        SET_VECTOR_ELT(out, i, values[i]);
        SET_STRING_ELT(nm, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(out, R_NamesSymbol, nm);
    UNPROTECT(2);
    return out;
}

/* The dtype name is attached to returned data arrays so the caller can see how
 * the value was stored, even though it arrives as a double. */
const char *mio_r_dtype_name(mio_dtype dtype) {
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
