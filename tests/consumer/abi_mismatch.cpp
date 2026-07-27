//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
/**
 * @file abi_mismatch.cpp
 * @brief The ABI-version link guard's negative test: a TU compiled against one
 * meshio++'s headers and linked against another meshio++'s library.
 *
 * The twin of `no_backend_macro.cpp`, for the second axis. Deliberately NOT part
 * of the `consumer_smoketest` target and not referenced by `CMakeLists.txt`: it
 * is compiled by hand, twice, against two different header trees.
 *
 * ### How the mismatch is staged, and why this way
 *
 * `MESHIOPLUSPLUS_ABI_VERSION` lives in `abi_version.hpp` and is deliberately not
 * overridable with `-D` -- a knob to change it would be a knob to silence the
 * guard while still miscompiling. So the mismatch is staged the way it actually
 * happens in the wild: **copy the installed headers and edit the copy**, which is
 * exactly a consumer whose include path points at an older or newer prefix than
 * its `-L` does.
 *
 *   cp -r <prefix>/include /tmp/h
 *   sed -i 's/MESHIOPLUSPLUS_ABI_VERSION 3/MESHIOPLUSPLUS_ABI_VERSION 99/' \
 *          /tmp/h/meshioplusplus/abi_version.hpp
 *
 *   # must FAIL, naming abi_version_is_99
 *   g++ -std=c++20 -DMESHIOPLUSPLUS_MESH_BACKEND_KRATOS -I /tmp/h abi_mismatch.cpp \
 *       -L <prefix>/lib -lmeshioplusplus_core_kratos -o /tmp/abi_bad
 *
 *   # positive control: must SUCCEED and run
 *   g++ -std=c++20 -DMESHIOPLUSPLUS_MESH_BACKEND_KRATOS -I <prefix>/include \
 *       abi_mismatch.cpp -L <prefix>/lib -lmeshioplusplus_core_kratos -o /tmp/abi_ok
 *
 * Both halves matter. Without the positive control, a typo'd include path would
 * satisfy the negative one for entirely the wrong reason -- the failure mode that
 * let the backend guard ship inert for a full release.
 *
 * The backend macro IS passed here, on purpose: it isolates the axis under test.
 * Without it this TU would fail on `mesh_backend_is_kratos` instead and prove
 * nothing about the ABI guard.
 */

// System includes
#include <cstdio>

// Project includes
#include "meshioplusplus/mesh.hpp"

int main() {
    // The guard fires from including mesh.hpp at all, not from what is called.
    // Touching the object keeps the compiler from dropping the whole thing
    // before the sentinel's relocation is emitted.
    meshioplusplus::Mesh mesh;
    std::printf("abi %d, points: %zu\n", MESHIOPLUSPLUS_ABI_VERSION,
                static_cast<std::size_t>(mesh.NumPoints()));
    return 0;
}
