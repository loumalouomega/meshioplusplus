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
 * @file detail/length_unit.hpp
 * @brief Length units as metadata (roadmap 7.1.4): the `field_data` key a mesh
 * states its length unit under, and the accessors the readers, writers and the
 * transform share.
 *
 * `field_data["meshio:length_unit"]` holds one float64: the number of **metres
 * per mesh unit** (`0.001` for millimetres). The UNV and Radioss readers fill
 * it from the unit the file states, their writers state it, and a transform
 * with a unit scale (`ScaleUnits`, `--scale-units`) updates it. A mesh without
 * the key states no unit. The Python twin is `meshioplusplus._units`.
 */

// System includes
#include <cmath>
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace detail {

/// The `field_data` key of a mesh's length unit (metres per unit, one float64).
inline constexpr const char* kLengthUnitName = "meshio:length_unit";

/**
 * @brief The metres per unit @p rMesh states, or `-1` when it states none (no
 * key, or a value that is not one positive finite number).
 */
inline double mesh_length_unit(const Mesh& rMesh) {
    if (!rMesh.HasFieldData(kLengthUnitName))
        return -1.0;
    const NDArray& value = rMesh.FieldData(kLengthUnitName);
    if (value.Size() != 1)
        return -1.0;
    double v = -1.0;
    switch (value.Dtype()) {
        case DType::Float64:
            v = value.As<double>()[0];
            break;
        case DType::Float32:
            v = static_cast<double>(value.As<float>()[0]);
            break;
        default:
            return -1.0;
    }
    return std::isfinite(v) && v > 0.0 ? v : -1.0;
}

/** @brief Record @p MetresPerUnit (positive) as @p rMesh's length unit. */
inline void set_mesh_length_unit(Mesh& rMesh, double MetresPerUnit) {
    NDArray value(DType::Float64, {std::size_t{1}});
    value.As<double>()[0] = MetresPerUnit;
    rMesh.AddFieldData(kLengthUnitName, std::move(value));
}

/**
 * @brief After coordinates were multiplied by @p UnitScale (a unit conversion,
 * `0.001` for millimetres to metres), the unit @p rMesh is in: its stated unit
 * divided by the scale. A mesh that states no unit is left alone.
 */
inline void rescale_mesh_length_unit(Mesh& rMesh, double UnitScale) {
    const double unit = mesh_length_unit(rMesh);
    if (unit > 0.0 && std::isfinite(UnitScale) && UnitScale > 0.0)
        set_mesh_length_unit(rMesh, unit / UnitScale);
}

}  // namespace detail
}  // namespace meshioplusplus
