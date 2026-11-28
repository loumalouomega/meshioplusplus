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
 * @file operations/units.hpp
 * @brief Per-field physical unit metadata (`field_data["meshio:unit:<array>"]`)
 * and consistent unit-system conversion (`convert_units`; roadmap 9.10.1).
 *
 * Each declared field unit is stored as 8 float64 values:
 * `[scale_to_SI, m, kg, s, K, A, mol, cd]`.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Prefix for per-array unit entries in `field_data`: `meshio:unit:<array>`.
inline constexpr const char* kFieldUnitPrefix = "meshio:unit:";

/// A physical unit represented as an SI scale factor + 7 SI base exponents
/// `[m, kg, s, K, A, mol, cd]`.
struct FieldUnit {
    double mScale = 1.0;
    std::array<double, 7> mExponents{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
};

/// Parse a common unit name (`"m"`, `"mm"`, `"Pa"`, `"MPa"`, `"GPa"`, `"N"`,
/// `"m/s"`, `"kg/m^3"`, `"t/mm^3"`, `"1"`, ...) into its SI scale and exponents.
MESHIOPLUSPLUS_API FieldUnit field_unit_from_name(const std::string& rName);

/// Attach `field_data["meshio:unit:" + rArrayName]` as Float64(8).
MESHIOPLUSPLUS_API void set_field_unit(Mesh& rMesh, const std::string& rArrayName,
                                       const FieldUnit& rUnit);

/// Read `field_data["meshio:unit:" + rArrayName]`; returns false if absent or malformed.
MESHIOPLUSPLUS_API bool get_field_unit(const Mesh& rMesh, const std::string& rArrayName,
                                       FieldUnit& rOut);

/// Options for `convert_units`.
struct ConvertUnitsOptions {
    /// Target consistent unit system: `"SI"`, `"mm-t-s"`, `"mm-kg-ms"`, `"cm-g-s"`, `"mm-g-ms"`.
    std::string mSystem = "SI";
    /// When true (default), throw if any point_data or cell_data array lacks a declared unit.
    /// When false, leave undeclared arrays untouched and record a warning.
    bool mStrict = true;
};

/// Result of `convert_units`.
struct ConvertUnitsResult {
    Mesh mMesh;
    std::int64_t mNumArraysConverted = 0;
    std::vector<std::string> mWarnings;
};

/// Convert mesh coordinates (`meshio:length_unit`) and all declared fields to `rOptions.mSystem`.
MESHIOPLUSPLUS_API ConvertUnitsResult convert_units(const Mesh& rMesh,
                                                    const ConvertUnitsOptions& rOptions = {});

}  // namespace meshioplusplus
