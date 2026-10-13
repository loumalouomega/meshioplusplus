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
 * @file formats/iges.hpp
 * @brief IGES (`.igs`, `.iges`) reader, through the optional OpenCASCADE
 * backend: the twin of `read_step` (formats/step.hpp) for IGES files.
 *
 * Read-only. The file's shape is triangulated by `BRepMesh` at
 * `CadReadOptions::mDeflection` into one `triangle` block with the `cad:*`
 * cell data of `BrepInfo`; the overload taking a `BrepInfo` also returns the
 * shape as NURBS faces, trimming loops and edges, which `write_cad_json`
 * writes as Kratos CAD JSON. Without `-DMESHIOPLUSPLUS_WITH_OCCT=ON` both
 * throw `Unsupported` naming the flag. See doc/formats/iges.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Read an IGES file as its triangulation.
/// @throws Unsupported in a build without OpenCASCADE; ReadError on a file
///         OpenCASCADE cannot read or transfer.
MESHIOPLUSPLUS_API Mesh read_iges(const std::string& rPath);

/// Read an IGES file as its triangulation and its B-rep.
MESHIOPLUSPLUS_API Mesh read_iges(const std::string& rPath, BrepInfo& rInfo,
                                  const CadReadOptions& rOptions = {});

}  // namespace meshioplusplus
