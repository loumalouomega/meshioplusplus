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

// The meshioplusplus-mmg extension: importing it registers the MMG remesher
// with the base wheel's shared core. The meshes never cross this boundary; the
// base `_core` converts the arguments and the core reaches the provider
// in-process.

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "../mmg_provider.hpp"
#include "meshioplusplus/version.hpp"

#ifndef MESHIOPLUSPLUS_MMG_VERSION
#define MESHIOPLUSPLUS_MMG_VERSION "unknown"
#endif

namespace py = pybind11;

PYBIND11_MODULE(_mmg, m) {
    m.doc() = "MMG remesher for meshio++";
    meshioplusplus::providers::register_mmg();
    m.attr("__version__") = MESHIOPLUSPLUS_VERSION_STRING;
    m.attr("__mmg_version__") = MESHIOPLUSPLUS_MMG_VERSION;
}
