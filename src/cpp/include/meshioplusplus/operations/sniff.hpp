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
#pragma once

/**
 * @file operations/sniff.hpp
 * @brief Content-based mesh-format detection, part of the dependency-free
 * mesh-operations layer.
 *
 * `sniff_format` reads the leading bytes / first lines of a file and returns
 * the meshio++ format name on a *confident* signature match, or `""` when it
 * cannot tell. It is used as a fallback by the read paths (C API, WASM, and
 * the Python reader) when the format cannot be inferred from the extension —
 * never on a write path. It is deliberately conservative: signatures shared by
 * several formats (e.g. the generic HDF5 magic used by med/h5m/cgns/hmf, or a
 * headerless binary STL) yield `""` rather than a guess.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {

/**
 * @brief Guess a mesh file's format from its contents.
 * @param rPath path to an existing, readable file
 * @return a meshio++ format name (e.g. `"vtu"`, `"gmsh"`) on a confident
 *         signature match, or `""` if the format cannot be determined (or the
 *         file cannot be opened)
 */
MESHIOPLUSPLUS_API std::string sniff_format(const std::string& rPath);

}  // namespace meshioplusplus
