/* The single MEX entry point: meshioplusplus_mex(command, args...).
 *
 * The command tables live in the source file that implements them (mesh.c,
 * ops_*.c, ...), one per file like the R binding's src/ layout. They are
 * merged and sorted once, on the first call, and looked up by binary search;
 * the argument count is checked before a handler runs, so a handler can index
 * a[0..nargs-1] without further checks -- the role R's registered
 * CallMethodDef table plays. */

#include "mio_oct.h"

#include <stdlib.h>
#include <string.h>

typedef const mio_o_entry *(*mio_o_table_fn)(size_t *n);

extern const mio_o_entry *mio_o_table_mesh(size_t *n);
extern const mio_o_entry *mio_o_table_ops_core(size_t *n);
extern const mio_o_entry *mio_o_table_ops_combine(size_t *n);
extern const mio_o_entry *mio_o_table_ops_data(size_t *n);
extern const mio_o_entry *mio_o_table_format_info(size_t *n);
extern const mio_o_entry *mio_o_table_sequence(size_t *n);
extern const mio_o_entry *mio_o_table_xdmf_series(size_t *n);
extern const mio_o_entry *mio_o_table_exodus_series(size_t *n);
extern const mio_o_entry *mio_o_table_render(size_t *n);

static const mio_o_table_fn kTables[] = {
    mio_o_table_mesh,        mio_o_table_ops_core,      mio_o_table_ops_combine,
    mio_o_table_ops_data,    mio_o_table_format_info,   mio_o_table_sequence,
    mio_o_table_xdmf_series, mio_o_table_exodus_series, mio_o_table_render,
};

static mio_o_entry *g_entries = NULL;
static size_t g_num_entries = 0;

static int compare_entries(const void *a, const void *b) {
    return strcmp(((const mio_o_entry *)a)->name, ((const mio_o_entry *)b)->name);
}

static void at_exit(void) {
    mio_o_release_all();
    free(g_entries);
    g_entries = NULL;
    g_num_entries = 0;
}

static void build_table(void) {
    size_t total = 0;
    for (size_t t = 0; t < sizeof(kTables) / sizeof(kTables[0]); ++t) {
        size_t n = 0;
        (void)kTables[t](&n);
        total += n;
    }
    /* Plain malloc: the table lives as long as the MEX. */
    g_entries = (mio_o_entry *)malloc((total > 0 ? total : 1) * sizeof(mio_o_entry));
    if (g_entries == NULL) mexErrMsgIdAndTxt("meshioplusplus:init", "meshio++: out of memory");
    size_t at = 0;
    for (size_t t = 0; t < sizeof(kTables) / sizeof(kTables[0]); ++t) {
        size_t n = 0;
        const mio_o_entry *e = kTables[t](&n);
        for (size_t i = 0; i < n; ++i) g_entries[at++] = e[i];
    }
    g_num_entries = total;
    qsort(g_entries, g_num_entries, sizeof(mio_o_entry), compare_entries);
    for (size_t i = 1; i < g_num_entries; ++i) {
        if (strcmp(g_entries[i - 1].name, g_entries[i].name) == 0) {
            mexErrMsgIdAndTxt("meshioplusplus:init", "meshio++: duplicate command '%s'",
                              g_entries[i].name);
        }
    }
    mexAtExit(at_exit);
}

static const mio_o_entry *find(const char *name) {
    size_t lo = 0, hi = g_num_entries;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int c = strcmp(name, g_entries[mid].name);
        if (c == 0) return &g_entries[mid];
        if (c < 0) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    return NULL;
}

/* Built-in commands that are not ports of an R entry point. */
static int builtin(const char *cmd, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (strcmp(cmd, "__commands__") == 0) {
        /* Every registered command, sorted: the parity test reads this. */
        mxArray *out = mxCreateCellMatrix((mwSize)g_num_entries, 1);
        for (size_t i = 0; i < g_num_entries; ++i) {
            mxSetCell(out, (mwIndex)i, mxCreateString(g_entries[i].name));
        }
        plhs[0] = out;
        return 1;
    }
    if (strcmp(cmd, "__release__") == 0) {
        /* Any handle kind: what mio.Handle/delete and mio.release call. */
        int released = 0;
        if (nrhs == 2) {
            int kind = mio_o_kind_of(prhs[1]);
            if (kind >= 0) released = mio_o_release(prhs[1], (mio_o_kind)kind);
        }
        if (nlhs > 0) plhs[0] = mxCreateLogicalScalar(released ? true : false);
        return 1;
    }
    if (strcmp(cmd, "__is_open__") == 0) {
        int open = 0;
        if (nrhs == 2) {
            int kind = mio_o_kind_of(prhs[1]);
            if (kind >= 0) open = mio_o_is_open(prhs[1], (mio_o_kind)kind);
        }
        plhs[0] = mxCreateLogicalScalar(open ? true : false);
        return 1;
    }
    return 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (g_entries == NULL) build_table();
    if (nrhs < 1 || !mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("meshioplusplus:usage",
                          "usage: meshioplusplus_mex(command, args...) -- call the mio.* "
                          "functions instead");
    }
    char cmd[128];
    if (mxGetString(prhs[0], cmd, sizeof(cmd)) != 0) {
        mexErrMsgIdAndTxt("meshioplusplus:usage", "meshio++: command name too long");
    }
    if (builtin(cmd, nlhs, plhs, nrhs, prhs)) return;
    const mio_o_entry *e = find(cmd);
    if (e == NULL) {
        mexErrMsgIdAndTxt("meshioplusplus:usage", "meshio++: unknown command '%s'", cmd);
    }
    if (nrhs - 1 != e->nargs) {
        mexErrMsgIdAndTxt("meshioplusplus:usage", "meshio++: '%s' takes %d argument(s), got %d",
                          cmd, e->nargs, nrhs - 1);
    }
    mxArray *out = e->fn(prhs + 1);
    if (out == NULL) {
        if (nlhs > 0) plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
    } else if (nlhs > 0 || out != NULL) {
        plhs[0] = out; /* ans */
    }
}
