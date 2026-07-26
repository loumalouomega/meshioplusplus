/* Shared declarations for the meshio++ R binding.
 *
 * The binding is plain C over R's own C API (.Call / SEXP / PROTECT), not
 * Rcpp: it keeps the dependency footprint at exactly zero and matches the flat
 * C ABI the rest of meshio++'s bindings sit on. The whole surface is scalars,
 * vectors and matrices, where Rf_allocVector / Rf_allocMatrix plus a memcpy is
 * all that is needed.
 *
 * Two facts govern everything here:
 *
 *  - COLUMN-MAJOR. R matrices are column-major and the C core row-major, so an
 *    R `(dim, num_points)` matrix is the SAME memory as the C API's row-major
 *    `(num_points, dim)`. Nothing is transposed anywhere.
 *  - COPY-ONLY. R vectors are R-managed, so the C API's zero-copy borrow
 *    cannot survive into R without ALTREP machinery that is out of scope.
 *    Every accessor copies, and the docs say so plainly rather than implying
 *    parity with the Julia binding.
 */

#ifndef MESHIOPLUSPLUS_R_H
#define MESHIOPLUSPLUS_R_H

#include <R.h>
#include <Rinternals.h>
#include <meshioplusplus/meshioplusplus.h>

/* --- handles ------------------------------------------------------------ */

/* Tag identifying our external pointers, so a foreign pointer is rejected
 * rather than dereferenced. */
extern SEXP mio_r_mesh_tag;

/* The same, for the transient XDMF writer's own opaque handle (xdmf_series.c).
 * A separate tag so a mesh and a series cannot be passed for one another. */
extern SEXP mio_r_series_tag;

/* Wrap an owning `mio_mesh*` in an external pointer with a registered
 * finalizer calling mio_mesh_free. The returned SEXP is NOT protected: a
 * caller that allocates again before returning it must PROTECT it first. */
SEXP mio_r_wrap_mesh(mio_mesh *mesh);

/* Extract the handle from an external pointer, erroring if it is not one of
 * ours or has already been released. */
mio_mesh *mio_r_mesh(SEXP x);

/* --- errors ------------------------------------------------------------- */

/* Raise an R condition carrying the C API's thread-local message. Never
 * returns (Rf_error does not). A `mio_status` code never reaches R. */
void mio_r_fail(const char *what);

/* Raise unless `status` is MIO_OK. */
void mio_r_check(mio_status status, const char *what);

/* Raise when `ptr` is NULL. */
void *mio_r_check_ptr(void *ptr, const char *what);

/* --- argument helpers --------------------------------------------------- */

/* A character(1) argument, or NULL for R_NilValue / NA / an empty vector. */
const char *mio_r_opt_string(SEXP s);

/* A character(1) argument that must be present. */
const char *mio_r_string(SEXP s, const char *what);

int mio_r_int(SEXP s, const char *what);
double mio_r_double(SEXP s, const char *what);
int mio_r_bool(SEXP s, const char *what);
int64_t mio_r_int64(SEXP s, const char *what);

/* Copy a length-3 numeric argument into `out`. */
void mio_r_vec3(SEXP s, double *out, const char *what);

/* Build a `const char* const*` array from an R character vector. `*count`
 * receives the length; `s` may be R_NilValue, giving a NULL pointer and a zero
 * count.
 *
 * The pointer array lives in a RAWSXP returned through `*shelter`, which is
 * left PROTECT-ed -- the caller MUST UNPROTECT(1) it once the C API call has
 * returned. (The CHARSXP bodies it points at are owned by `s`, which the
 * caller keeps alive for the same span.) */
const char *const *mio_r_names(SEXP s, int64_t *count, SEXP *shelter);

/* Name of a `mio_dtype`, attached to returned data arrays so the caller can
 * see how a value was stored even though it arrives as a double. */
const char *mio_r_dtype_name(mio_dtype dtype);

/* --- result helpers ----------------------------------------------------- */

/* Run the C API's string-getter protocol into a fresh R character(1). */
typedef int64_t (*mio_r_str_getter)(void *ctx, char *buf, int64_t buflen);
SEXP mio_r_getstring(mio_r_str_getter getter, void *ctx, const char *what);

/* Copy a typed C buffer of `n` elements into a fresh REALSXP. */
SEXP mio_r_copy_as_real(const void *data, mio_dtype dtype, R_xlen_t n);

/* Copy an int64 index map into a REALSXP, shifting to 1-based and turning the
 * C API's -1 "pruned / absent" sentinel into 0 -- never a valid 1-based index,
 * so it stays unambiguous. This is the Fortran and Julia bindings' rule; the
 * three must agree. */
SEXP mio_r_shift_map(const int64_t *data, int64_t n);

/* A named list of `n` elements; `names` must have `n` entries. */
SEXP mio_r_named_list(int n, const char **names, SEXP *values);

#endif /* MESHIOPLUSPLUS_R_H */
