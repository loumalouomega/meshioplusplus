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
 * @file no_backend_macro.cpp
 * @brief The mesh-backend link guard's negative test: a TU that gets the include
 * path but **not** the backend definitions.
 *
 * Deliberately NOT part of the `consumer_smoketest` target and not referenced by
 * `CMakeLists.txt`. It is compiled by hand, with `-I<prefix>/include` only, to
 * simulate the consumer `detail/mesh_backend_check.hpp` exists for: pkg-config,
 * a hand-written makefile, or an IDE that picked up the include path without the
 * imported target's `INTERFACE_COMPILE_DEFINITIONS`. A CMake consumer cannot
 * reach this state, which is exactly why `tests/consumer/main.cpp` cannot host
 * the check -- it `#error`s on a missing macro before the linker is ever run.
 *
 * With no macro defined, `mesh.hpp` silently falls through to its MESHIO
 * default, so this TU's `Mesh` is a different type from the one a NATIVE/KRATOS
 * library was built with. Linking it against such a library **must fail**, and
 * the failure must name `mesh_backend_is_<backend>` rather than emitting raw
 * mangled undefined references to whatever the TU happened to call.
 *
 * Expected uses (both halves matter -- see the CI step; without the positive
 * control a typo'd include path would satisfy the negative one for the wrong
 * reason):
 *
 *   # must FAIL, naming mesh_backend_is_kratos
 *   g++ -std=c++20 -I <prefix>/include no_backend_macro.cpp \
 *       -L <prefix>/lib -lmeshioplusplus_core_kratos -o /tmp/nb
 *
 *   # must SUCCEED and run
 *   g++ -std=c++20 -DMESHIOPLUSPLUS_MESH_BACKEND_KRATOS -I <prefix>/include \
 *       no_backend_macro.cpp -L <prefix>/lib -lmeshioplusplus_core_kratos -o /tmp/nb_ok
 */

// System includes
#include <cstdio>

// Project includes
#include "meshioplusplus/mesh.hpp"

int main() {
    // Any use of Mesh will do: the guard fires from including mesh.hpp at all,
    // not from what is called. Touching the object keeps a compiler from
    // deciding the whole thing is dead before the reference is emitted.
    meshioplusplus::Mesh mesh;
    std::printf("points: %zu\n", static_cast<std::size_t>(mesh.NumPoints()));
    return 0;
}
