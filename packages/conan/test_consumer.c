/* Smallest possible meshio++ C API consumer: link the packaged library and
 * exercise a couple of symbols. Success = the config-package resolved. */
#include <meshioplusplus/meshioplusplus.h>
#include <stdio.h>

int main(void) {
    printf("meshio++ %s (backend: %s)\n", mio_version(), mio_mesh_backend());
    mio_mesh* m = mio_mesh_create();
    if (!m) {
        fprintf(stderr, "mio_mesh_create failed\n");
        return 1;
    }
    mio_mesh_free(m);
    printf("test_package: OK\n");
    return 0;
}
