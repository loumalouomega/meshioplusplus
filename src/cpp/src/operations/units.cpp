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

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/units.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/mesh_api.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kUnitsPrefix = "meshio++: convert_units: ";

struct UnitSystemDef {
    const char* mpName;
    double mLength;  // metres per system length unit
    double mMass;    // kg per system mass unit
    double mTime;    // seconds per system time unit
};

const UnitSystemDef* units_find_system(const std::string& rName) {
    static const UnitSystemDef kSystems[] = {
        {"SI", 1.0, 1.0, 1.0},      {"si", 1.0, 1.0, 1.0},         {"m-kg-s", 1.0, 1.0, 1.0},
        {"mm-t-s", 1e-3, 1e3, 1.0}, {"mm-kg-ms", 1e-3, 1.0, 1e-3}, {"cm-g-s", 1e-2, 1e-3, 1.0},
        {"cgs", 1e-2, 1e-3, 1.0},   {"mm-g-ms", 1e-3, 1e-3, 1e-3},
    };
    for (const auto& s : kSystems) {
        if (rName == s.mpName)
            return &s;
    }
    return nullptr;
}

double units_target_scale(const UnitSystemDef& rSys, const FieldUnit& rUnit) {
    return std::pow(rSys.mLength, rUnit.mExponents[0]) * std::pow(rSys.mMass, rUnit.mExponents[1]) *
           std::pow(rSys.mTime, rUnit.mExponents[2]);
}

NDArray units_scale_array(const NDArray& rArr, double Factor) {
    NDArray out = detail::data_owned_copy(rArr);
    if (Factor == 1.0)
        return out;
    if (out.Dtype() == DType::Float64) {
        auto* p = out.As<double>();
        for (std::size_t i = 0; i < out.Size(); ++i)
            p[i] *= Factor;
    } else if (out.Dtype() == DType::Float32) {
        auto* p = out.As<float>();
        for (std::size_t i = 0; i < out.Size(); ++i)
            p[i] = static_cast<float>(static_cast<double>(p[i]) * Factor);
    }
    return out;
}

}  // namespace

FieldUnit field_unit_from_name(const std::string& rName) {
    // Exponents: {m, kg, s, K, A, mol, cd}
    static const std::pair<const char*, FieldUnit> kTable[] = {
        {"1", {1.0, {0, 0, 0, 0, 0, 0, 0}}},
        {"dimensionless", {1.0, {0, 0, 0, 0, 0, 0, 0}}},
        {"m", {1.0, {1, 0, 0, 0, 0, 0, 0}}},
        {"mm", {1e-3, {1, 0, 0, 0, 0, 0, 0}}},
        {"cm", {1e-2, {1, 0, 0, 0, 0, 0, 0}}},
        {"km", {1e3, {1, 0, 0, 0, 0, 0, 0}}},
        {"in", {0.0254, {1, 0, 0, 0, 0, 0, 0}}},
        {"ft", {0.3048, {1, 0, 0, 0, 0, 0, 0}}},
        {"kg", {1.0, {0, 1, 0, 0, 0, 0, 0}}},
        {"g", {1e-3, {0, 1, 0, 0, 0, 0, 0}}},
        {"t", {1e3, {0, 1, 0, 0, 0, 0, 0}}},
        {"s", {1.0, {0, 0, 1, 0, 0, 0, 0}}},
        {"ms", {1e-3, {0, 0, 1, 0, 0, 0, 0}}},
        {"K", {1.0, {0, 0, 0, 1, 0, 0, 0}}},
        {"A", {1.0, {0, 0, 0, 0, 1, 0, 0}}},
        {"Pa", {1.0, {-1, 1, -2, 0, 0, 0, 0}}},
        {"kPa", {1e3, {-1, 1, -2, 0, 0, 0, 0}}},
        {"MPa", {1e6, {-1, 1, -2, 0, 0, 0, 0}}},
        {"GPa", {1e9, {-1, 1, -2, 0, 0, 0, 0}}},
        {"bar", {1e5, {-1, 1, -2, 0, 0, 0, 0}}},
        {"psi", {6894.757293168361, {-1, 1, -2, 0, 0, 0, 0}}},
        {"N", {1.0, {1, 1, -2, 0, 0, 0, 0}}},
        {"kN", {1e3, {1, 1, -2, 0, 0, 0, 0}}},
        {"m/s", {1.0, {1, 0, -1, 0, 0, 0, 0}}},
        {"mm/s", {1e-3, {1, 0, -1, 0, 0, 0, 0}}},
        {"m/s^2", {1.0, {1, 0, -2, 0, 0, 0, 0}}},
        {"mm/s^2", {1e-3, {1, 0, -2, 0, 0, 0, 0}}},
        {"kg/m^3", {1.0, {-3, 1, 0, 0, 0, 0, 0}}},
        {"t/mm^3", {1e12, {-3, 1, 0, 0, 0, 0, 0}}},
        {"J", {1.0, {2, 1, -2, 0, 0, 0, 0}}},
        {"W", {1.0, {2, 1, -3, 0, 0, 0, 0}}},
    };
    for (const auto& [name, unit] : kTable) {
        if (rName == name)
            return unit;
    }
    throw std::invalid_argument(std::string(kUnitsPrefix) + "unknown unit name '" + rName + "'");
}

void set_field_unit(Mesh& rMesh, const std::string& rArrayName, const FieldUnit& rUnit) {
    if (!std::isfinite(rUnit.mScale) || rUnit.mScale <= 0.0)
        throw std::invalid_argument(std::string(kUnitsPrefix) +
                                    "unit scale must be positive and finite");
    NDArray arr = NDArray::Uninit(DType::Float64, {std::size_t{8}});
    auto* p = arr.As<double>();
    p[0] = rUnit.mScale;
    for (std::size_t k = 0; k < 7; ++k)
        p[k + 1] = rUnit.mExponents[k];
    rMesh.AddFieldData(std::string(kFieldUnitPrefix) + rArrayName, std::move(arr));
}

bool get_field_unit(const Mesh& rMesh, const std::string& rArrayName, FieldUnit& rOut) {
    const std::string key = std::string(kFieldUnitPrefix) + rArrayName;
    if (!rMesh.HasFieldData(key))
        return false;
    const NDArray& arr = rMesh.FieldData(key);
    if (arr.Size() != 8 || arr.Dtype() != DType::Float64)
        return false;
    const double* p = arr.As<double>();
    if (!std::isfinite(p[0]) || p[0] <= 0.0)
        return false;
    rOut.mScale = p[0];
    for (std::size_t k = 0; k < 7; ++k)
        rOut.mExponents[k] = p[k + 1];
    return true;
}

ConvertUnitsResult convert_units(const Mesh& rMesh, const ConvertUnitsOptions& rOptions) {
    const UnitSystemDef* p_sys = units_find_system(rOptions.mSystem);
    if (!p_sys)
        throw std::invalid_argument(
            std::string(kUnitsPrefix) + "unknown unit system '" + rOptions.mSystem +
            "' (expected 'SI', 'mm-t-s', 'mm-kg-ms', 'cm-g-s', or 'mm-g-ms')");

    ConvertUnitsResult res;
    Mesh out = detail::clone_mesh(rMesh);

    // 1. Scale coordinates if meshio:length_unit is declared.
    const double src_len = detail::mesh_length_unit(rMesh);
    if (src_len > 0.0) {
        const double factor = src_len / p_sys->mLength;
        out.AssignPoints(units_scale_array(rMesh.Points(), factor));
        detail::set_mesh_length_unit(out, p_sys->mLength);
    }

    // 2. Scale point_data arrays with declared units.
    for (const std::string& name : rMesh.PointDataNames()) {
        FieldUnit u;
        if (!get_field_unit(rMesh, name, u)) {
            if (rOptions.mStrict)
                throw std::invalid_argument(std::string(kUnitsPrefix) + "point_data array '" +
                                            name + "' has no declared unit (set '" +
                                            kFieldUnitPrefix + name + "' or pass strict=False)");
            res.mWarnings.push_back("point_data '" + name +
                                    "' has no declared unit; left unchanged");
            continue;
        }
        const double tgt_scale = units_target_scale(*p_sys, u);
        const double factor = u.mScale / tgt_scale;
        out.AddPointData(name, units_scale_array(rMesh.PointData(name), factor));
        FieldUnit next = u;
        next.mScale = tgt_scale;
        set_field_unit(out, name, next);
        ++res.mNumArraysConverted;
    }

    // 3. Scale cell_data arrays with declared units.
    for (const std::string& name : rMesh.CellDataNames()) {
        FieldUnit u;
        if (!get_field_unit(rMesh, name, u)) {
            if (rOptions.mStrict)
                throw std::invalid_argument(std::string(kUnitsPrefix) + "cell_data array '" + name +
                                            "' has no declared unit (set '" + kFieldUnitPrefix +
                                            name + "' or pass strict=False)");
            res.mWarnings.push_back("cell_data '" + name +
                                    "' has no declared unit; left unchanged");
            continue;
        }
        const double tgt_scale = units_target_scale(*p_sys, u);
        const double factor = u.mScale / tgt_scale;
        const std::size_t nb = rMesh.CellDataNumBlocks(name);
        std::vector<NDArray> blocks;
        blocks.reserve(nb);
        for (std::size_t b = 0; b < nb; ++b)
            blocks.push_back(units_scale_array(rMesh.CellData(name, b), factor));
        out.AddCellData(name, std::move(blocks));
        FieldUnit next = u;
        next.mScale = tgt_scale;
        set_field_unit(out, name, next);
        ++res.mNumArraysConverted;
    }

    res.mMesh = std::move(out);
    return res;
}

}  // namespace meshioplusplus
