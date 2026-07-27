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
 * @file abi_version_check.cpp
 * @brief The one definition of the ABI sentinel declared in
 * `detail/abi_version_check.hpp`.
 *
 * Which symbol this is depends entirely on the `MESHIOPLUSPLUS_ABI_VERSION` this
 * translation unit was compiled with, so a library built from one set of headers
 * defines a different name from the one a consumer compiled against another set
 * references, and the two fail to link. See the header for the full rationale,
 * and `src/cpp/src/mesh_backend_check.cpp` for the same pattern on the backend
 * axis.
 */

// Project includes
#include "meshioplusplus/detail/abi_version_check.hpp"

namespace meshioplusplus::detail {

void MESHIOPLUSPLUS_ABI_SYM(MESHIOPLUSPLUS_ABI_VERSION)() {}

}  // namespace meshioplusplus::detail
