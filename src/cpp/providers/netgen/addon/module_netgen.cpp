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

// The meshioplusplus-netgen extension: importing it registers the Netgen
// tetrahedral mesher with the base wheel's shared core. The meshes never cross
// this boundary; the base `_core` converts the arguments and the core reaches
// the provider in-process.

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "../netgen_provider.hpp"
#include "meshioplusplus/version.hpp"

#ifndef MESHIOPLUSPLUS_NETGEN_VERSION
#define MESHIOPLUSPLUS_NETGEN_VERSION "unknown"
#endif

namespace py = pybind11;

PYBIND11_MODULE(_netgen, m) {
    m.doc() = "Netgen tetrahedral mesher for meshio++";
    meshioplusplus::providers::register_netgen();
    m.attr("__version__") = MESHIOPLUSPLUS_VERSION_STRING;
    m.attr("__netgen_version__") = MESHIOPLUSPLUS_NETGEN_VERSION;
}
