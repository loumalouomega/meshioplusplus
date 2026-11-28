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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file operations/mesh_boolean.hpp
 * @brief Watertight surface mesh boolean operations (`Union`, `Difference`,
 * `Intersection`) dispatched through a registered `BooleanProvider`
 * (e.g. Manifold, Apache-2.0; roadmap 9.1).
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Boolean set operation on two watertight triangle surfaces.
enum class BooleanOp : std::int32_t {
    Union = 0,
    Difference = 1,
    Intersection = 2,
};

/// Parse `"union"`, `"difference"` / `"subtract"`, or `"intersection"` / `"intersect"`.
MESHIOPLUSPLUS_API BooleanOp boolean_op_from_name(const std::string& rName);

/// Options for `mesh_boolean`.
struct MeshBooleanOptions {
    BooleanOp mOp = BooleanOp::Union;
    /// Optional provider name; empty selects the default registered `BooleanProvider`.
    std::string mBackend;
};

/// Result of `mesh_boolean`.
struct MeshBooleanResult {
    Mesh mMesh;
    std::string mBackend;
    std::int64_t mNumTriangles = 0;
};

/// Perform a surface boolean operation on two closed triangle meshes.
MESHIOPLUSPLUS_API MeshBooleanResult mesh_boolean(const Mesh& rA, const Mesh& rB,
                                                  const MeshBooleanOptions& rOptions = {});

/// Whether a `BooleanProvider` is currently registered.
MESHIOPLUSPLUS_API bool has_boolean_provider() noexcept;

}  // namespace meshioplusplus
