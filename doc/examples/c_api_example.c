/* Minimal meshio++ C API consumer: build a tet mesh, write it, read it back. */
#include <meshioplusplus/meshioplusplus.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    const double points[15] = {0, 0, 0, 1.1, 0.2, 0.3, 0.4, 1.2, 0.5, 0.6, 0.7, 1.3, 1.4, 1.5, 1.6};
    const int64_t conn[8] = {0, 1, 2, 3, 1, 2, 3, 4};

    printf("meshio++ %s (backend: %s)\n", mio_version(), mio_mesh_backend());

    mio_mesh* m = mio_mesh_create();
    if (mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, points) != MIO_OK ||
        mio_mesh_add_cell_block(m, "tetra", 2, 4, MIO_INT64, conn) != MIO_OK ||
        mio_write("/tmp/mio_example.vtu", m, NULL) != MIO_OK) {
        fprintf(stderr, "build/write failed: %s\n", mio_last_error());
        return 1;
    }
    mio_mesh_free(m);

    mio_mesh* r = mio_read("/tmp/mio_example.vtu", NULL);
    if (!r) {
        fprintf(stderr, "read failed: %s\n", mio_last_error());
        return 1;
    }
    const void* pts = NULL;
    mio_dtype dt;
    if (mio_mesh_get_points(r, &pts, &dt) != MIO_OK || dt != MIO_FLOAT64 ||
        mio_mesh_num_points(r) != 5 || ((const double*)pts)[4] != 0.2) {
        fprintf(stderr, "verification failed\n");
        return 1;
    }
    char type[32];
    mio_mesh_cell_block_type(r, 0, type, sizeof type);
    if (strcmp(type, "tetra") != 0) {
        fprintf(stderr, "unexpected cell type '%s'\n", type);
        return 1;
    }
    mio_mesh_free(r);
    remove("/tmp/mio_example.vtu");
    printf("example.c: OK\n");
    return 0;
}
