// SPDX-License-Identifier: MIT
// Fault injection for the real HDF5 custom mutator. This test target replaces
// object traversal so every structure-aware mutation fails in its child.
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <hdf5.h>

namespace {

[[noreturn]] void mutator_probe_fail() {
    std::fprintf(stderr, "mutator probe: injected failure\n");
    const char* mode = std::getenv("MIO_MUTATOR_PROBE_FAILURE");
    if (mode && std::strcmp(mode, "sanitizer") == 0) {
        volatile int* p_value = new int(0);
        delete p_value;
        *p_value = 1;  // intentional ASan report in the mutator child
    }
    std::raise(SIGABRT);
    std::abort();
}

}  // namespace

#if H5_VERSION_GE(1, 12, 0)
extern "C" herr_t H5Ovisit3(hid_t, H5_index_t, H5_iter_order_t, H5O_iterate2_t, void*, unsigned) {
    mutator_probe_fail();
}
#elif H5_VERSION_GE(1, 10, 3)
extern "C" herr_t H5Ovisit2(hid_t, H5_index_t, H5_iter_order_t, H5O_iterate_t, void*, unsigned) {
    mutator_probe_fail();
}
#else
extern "C" herr_t H5Ovisit(hid_t, H5_index_t, H5_iter_order_t, H5O_iterate_t, void*) {
    mutator_probe_fail();
}
#endif

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t*, std::size_t) {
    return 0;
}
