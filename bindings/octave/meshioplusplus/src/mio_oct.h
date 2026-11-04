/* Shared declarations for the meshio++ Octave / MATLAB binding.
 *
 * The binding is plain C over the MEX API (mex.h / matrix.h) -- the one native
 * interface GNU Octave and MATLAB both implement -- layered on the installed
 * flat C API, exactly like the R binding sits on .Call. Nothing here includes
 * Octave's C++ headers (oct.h), so the same sources build with `mkoctfile
 * --mex` and with MATLAB's `mex`.
 *
 * One MEX file, `meshioplusplus_mex`, exports every entry point. Its first
 * argument names the command; mio_mex.c looks it up in the per-file tables
 * (`mio_o_entries_*`) and checks the argument count, the counterpart of R's
 * registered CallMethodDef table. Each handler has the shape
 *
 *     mxArray *O_mio_x(const mxArray **a);   // a[0] .. a[nargs-1]
 *
 * and returns its result (NULL means "no value": the gateway returns []).
 *
 * The facts that govern every handler:
 *
 *  - ROWS ARE ENTITIES. MATLAB and Octave arrays are column-major like R's,
 *    but the MATLAB idiom puts one point / cell / entry per ROW
 *    (`trisurf(conn, p(:,1), p(:,2), p(:,3))`). Every per-entity array crosses
 *    as (num_entities x components): the helpers below transpose while they
 *    copy, so the C side and the user both see their natural layout. Plain
 *    vectors are columns (n x 1). Mathematical matrices (an affine 4x4, a 3x3
 *    tensor) keep their meaning, so the R binding's transpose for those is
 *    kept verbatim.
 *  - COPY-ONLY. Every accessor copies; MATLAB arrays are runtime-managed.
 *  - 1-BASED. Connectivity, block, region and index-map entries are 1-based
 *    doubles; index maps turn the C API's -1 sentinel into 0.
 *  - NATIVE TYPES for data arrays. Unlike R, MATLAB has every integer width,
 *    so a point/cell/field array keeps its stored class (int32, uint8,
 *    single, ...). Indices and counts are double, exact to 2^53.
 *  - ERRORS raise through mexErrMsgIdAndTxt (identifier
 *    "meshioplusplus:<what>") and never return. Temporary memory comes from
 *    mxMalloc / mxCalloc, which both runtimes release automatically when a MEX
 *    call ends or errors -- the R_alloc counterpart. Native objects the C API
 *    allocated (a poly_conn, a result) must still be freed by the handler
 *    before it raises, exactly as in the R shim.
 *  - HANDLES are classdef objects (mio.MeshHandle, mio.FormatInfo, ...) built from
 *    inside the MEX. They carry an opaque uint64 id into a registry that the
 *    MEX owns; an id that is unknown, stale (released, or loaded back from a
 *    .mat file) or of the wrong kind raises an error instead of being
 *    dereferenced. The MEX stays locked while any handle is live, so `clear
 *    mex` cannot strand one.
 */

#ifndef MESHIOPLUSPLUS_OCTAVE_H
#define MESHIOPLUSPLUS_OCTAVE_H

#include <stddef.h>
#include <stdint.h>

#include "mex.h"
#include <meshioplusplus/meshioplusplus.h>

/* --- dispatch table ------------------------------------------------------ */

typedef mxArray *(*mio_o_fn)(const mxArray **a);

typedef struct {
    const char *name; /* command name, e.g. "mio_read" */
    mio_o_fn fn;
    int nargs; /* exact argument count after the command name */
} mio_o_entry;

/* Every source file `<file>.c` ends with its table:
 *
 *     static const mio_o_entry kEntries[] = { MIO_O_ENTRY(mio_x, 2), ... };
 *     const mio_o_entry *mio_o_table_<file>(size_t *n) {
 *         *n = sizeof(kEntries) / sizeof(kEntries[0]);
 *         return kEntries;
 *     }
 *
 * and mio_mex.c lists it. MIO_O_ENTRY(mio_x, n) binds the command "mio_x" to
 * the handler O_mio_x taking n arguments. */
#define MIO_O_ENTRY(name, n) {#name, &O_##name, n}

/* --- handles ------------------------------------------------------------- */

typedef enum {
    MIO_O_MESH = 0,
    MIO_O_FORMAT_INFO,
    MIO_O_XDMF_SERIES,
    MIO_O_SEQUENCE,
    MIO_O_EXODUS_SERIES,
    MIO_O_NUM_KINDS
} mio_o_kind;

/* Register `ptr` (owned; freed with the kind's free function on release) and
 * return a new classdef object of the kind's class. On failure to build the
 * object the pointer is freed before raising. */
mxArray *mio_o_wrap(mio_o_kind kind, void *ptr);

/* The live pointer behind a handle object of `kind`, or an error naming
 * `what`. Never returns NULL. */
void *mio_o_unwrap(const mxArray *obj, mio_o_kind kind, const char *what);

/* Free and unregister. Idempotent: an already-released handle is a no-op.
 * Returns 1 when something was freed. A non-handle or wrong-kind argument
 * raises. */
int mio_o_release(const mxArray *obj, mio_o_kind kind);

/* 1 while the handle is live, 0 once released or for any non-handle. */
int mio_o_is_open(const mxArray *obj, mio_o_kind kind);

/* The kind of a handle object, or -1 for anything else. */
int mio_o_kind_of(const mxArray *obj);

/* Release every live handle (mexAtExit). */
void mio_o_release_all(void);

/* Convenience wrappers for the mesh kind. */
mxArray *mio_o_wrap_mesh(mio_mesh *mesh);
mio_mesh *mio_o_mesh(const mxArray *obj);

/* --- errors -------------------------------------------------------------- */

/* Raise carrying the C API's thread-local message. Never returns. */
void mio_o_fail(const char *what);

/* Raise a formatted error with identifier meshioplusplus:<what>. */
void mio_o_error(const char *what, const char *fmt, ...);

/* Raise unless `status` is MIO_OK. */
void mio_o_check(mio_status status, const char *what);

/* Raise when `ptr` is NULL. */
void *mio_o_check_ptr(void *ptr, const char *what);

/* --- argument helpers ---------------------------------------------------- */

/* The R binding's NULL is MATLAB's [] (an empty double). `{}` and '' are NOT
 * null: they are an empty list and an empty string. */
int mio_o_is_null(const mxArray *a);

/* A char row vector, or NULL for [] / any non-char argument. The returned
 * string is mxMalloc'd (auto-freed). '' gives "". */
const char *mio_o_opt_string(const mxArray *a);

/* A char argument that must be present. */
const char *mio_o_string(const mxArray *a, const char *what);

int mio_o_int(const mxArray *a, const char *what);
double mio_o_double(const mxArray *a, const char *what);
int mio_o_bool(const mxArray *a, const char *what);
int64_t mio_o_int64(const mxArray *a, const char *what);

/* Copy a length-3 numeric argument into `out`. */
void mio_o_vec3(const mxArray *a, double *out, const char *what);

/* Copy exactly `n` numeric values into `out` (any numeric class). */
void mio_o_vecn(const mxArray *a, double *out, size_t n, const char *what);

/* Every element of a numeric array as a fresh mxMalloc'd double buffer (column
 * order); `*n` receives the element count. */
double *mio_o_doubles(const mxArray *a, size_t *n, const char *what);

/* A cellstr (or a single char row) as a `const char* const*` array, mxMalloc'd.
 * [] gives NULL with count 0; {} gives NULL with count 0 too -- use
 * mio_o_is_null() to tell "absent" from "empty" where the ABI cares. */
const char *const *mio_o_names(const mxArray *a, int64_t *count, const char *what);

/* An (n x k) numeric matrix as a row-major (n, k) buffer, converted to double
 * and with `offset` subtracted from each element (1.0 turns 1-based indices
 * into 0-based). `*rows`/`*cols` receive n and k. */
double *mio_o_rows_to_c(const mxArray *a, int64_t *rows, int64_t *cols, double offset,
                        const char *what);

/* An (n x k) array of 1-based indices as a row-major int64 (n, k) buffer,
 * shifted to 0-based. Every value must be a finite integer >= 1. */
int64_t *mio_o_index_rows(const mxArray *a, int64_t *rows, int64_t *cols, const char *what);

/* A vector of 1-based indices as a 0-based int64 buffer. */
int64_t *mio_o_index_vector(const mxArray *a, int64_t *n, const char *what);

/* A typed (n, comps...) data array for a setter: keeps a numeric argument's
 * class (double, single, int*, uint*, logical as uint8), returning a row-major
 * buffer laid out as the C API expects plus its dtype and shape. A vector
 * becomes shape (n). */
void *mio_o_data_to_c(const mxArray *a, mio_dtype *dtype, int32_t *ndim, int64_t *shape,
                      const char *what);

/* --- result helpers ------------------------------------------------------ */

mxArray *mio_o_str(const char *s);          /* char row vector ("" -> '') */
mxArray *mio_o_scalar(double v);            /* double scalar */
mxArray *mio_o_logical(int v);              /* logical scalar */
mxArray *mio_o_empty(void);                 /* [] -- the R NULL */
mxArray *mio_o_cellstr(int64_t n);          /* n x 1 cell, filled by the caller */
mxArray *mio_o_cell(int64_t n);             /* n x 1 cell */
mxArray *mio_o_column(int64_t n);           /* n x 1 double, zeroed */
mxArray *mio_o_matrix(int64_t rows, int64_t cols); /* rows x cols double, zeroed */
mxArray *mio_o_logical_column(int64_t n);   /* n x 1 logical */

/* Run the C API's string-getter protocol into a fresh char row vector. */
typedef int64_t (*mio_o_str_getter)(void *ctx, char *buf, int64_t buflen);
mxArray *mio_o_getstring(mio_o_str_getter getter, void *ctx, const char *what);

/* The same, into a mxMalloc'd C string. */
char *mio_o_getcstring(mio_o_str_getter getter, void *ctx, const char *what);

/* Convert `n` elements of a typed C buffer into doubles, adding `offset`. */
void mio_o_fill_real(double *dst, const void *data, mio_dtype dtype, size_t n, double offset);

/* A typed C buffer as an n x 1 double column. */
mxArray *mio_o_copy_as_real(const void *data, mio_dtype dtype, int64_t n);

/* A row-major (rows, cols) C buffer as a rows x cols double matrix, adding
 * `offset` to each element (1.0 for the 1-based shift). */
mxArray *mio_o_rows_from_c(const void *data, mio_dtype dtype, int64_t rows, int64_t cols,
                           double offset);

/* A row-major data array of the given dtype and shape as a MATLAB array of the
 * matching class with the SAME dims: (n) -> n x 1, (n, k) -> n x k,
 * (n, a, b) -> n x a x b. */
mxArray *mio_o_typed_from_c(const void *data, mio_dtype dtype, int32_t ndim,
                            const int64_t *shape);

/* An int64 index map as an n x 1 double column, shifted to 1-based with the
 * -1 "absent" sentinel as 0. */
mxArray *mio_o_shift_map(const int64_t *data, int64_t n);

/* A plain int64 buffer as an n x 1 double column (no shift). */
mxArray *mio_o_int64_column(const int64_t *data, int64_t n);

/* A 1 x 1 struct with the given fields. The values are adopted (they must be
 * freshly created arrays, or NULL for []). */
mxArray *mio_o_struct(int n, const char **names, mxArray **values);

/* An n x 1 struct array with the given field names, values set by the caller
 * with mxSetField. */
mxArray *mio_o_struct_array(int64_t n, int nfields, const char **names);

/* Name of a `mio_dtype`. */
const char *mio_o_dtype_name(mio_dtype dtype);

#endif /* MESHIOPLUSPLUS_OCTAVE_H */
