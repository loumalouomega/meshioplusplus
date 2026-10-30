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

// The meshioplusplus-occt extension: importing it registers the OpenCASCADE
// CAD kernel with the base wheel's shared core. The meshes never cross this
// boundary; the base `_core` converts the arguments and the core reaches the
// provider in-process.

// External includes
#include <pybind11/pybind11.h>

#include <Standard_Version.hxx>

// Project includes
#include "meshioplusplus/version.hpp"
#include "../occt_provider.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_occt, m) {
    m.doc() = "OpenCASCADE CAD kernel for meshio++";
    meshioplusplus::providers::register_occt();
    m.attr("__version__") = MESHIOPLUSPLUS_VERSION_STRING;
    m.attr("__occt_version__") = OCC_VERSION_COMPLETE;
}
