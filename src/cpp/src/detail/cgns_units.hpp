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
 * @file detail/cgns_units.hpp
 * @brief SIDS physical units of a CGNS `DataArray_t` and `meshio:unit:<array>`
 * (roadmap 9.10.1), shared by the raw HDF5 path (formats/cgns.cpp) and the
 * cgnslib one (formats/cgns_mll.cpp).
 *
 * A field unit is written as the array's own `DimensionalUnits_t` (always SI:
 * Kilogram, Meter, Second, Kelvin, Radian), `DimensionalExponents_t` (Mass,
 * Length, Time, Temperature, Angle), `AdditionalExponents_t` (ElectricCurrent,
 * Substance, LuminousIntensity) when any is nonzero, and `DataConversion_t`
 * (ConversionScale, ConversionOffset = 0) when the scale is not 1: SIDS reads
 * `raw * ConversionScale + ConversionOffset` in the DimensionalUnits.
 */

// System includes
#include <array>
#include <cmath>
#include <string>
#include <utility>

// Project includes
#include "meshioplusplus/operations/units.hpp"

namespace meshioplusplus {
namespace detail {

/// The SI unit names a written `DimensionalUnits_t` carries, in SIDS order.
inline const std::array<const char*, 5>& cgns_si_unit_names() {
    static const std::array<const char*, 5> kNames = {"Kilogram", "Meter", "Second", "Kelvin",
                                                      "Radian"};
    return kNames;
}

/// SI per unit of a SIDS unit name of slot @p Slot (0 Mass .. 4 Angle); 0 for
/// `Null`/`UserDefined`/an offset temperature scale, which have no factor.
inline double cgns_unit_factor(std::size_t Slot, const std::string& rName) {
    static const std::pair<const char*, double> kMass[] = {
        {"Kilogram", 1.0}, {"Gram", 1e-3}, {"Slug", 14.593902937206364}, {"PoundMass", 0.45359237}};
    static const std::pair<const char*, double> kLength[] = {{"Meter", 1.0},
                                                             {"Centimeter", 0.01},
                                                             {"Millimeter", 0.001},
                                                             {"Foot", 0.3048},
                                                             {"Inch", 0.0254}};
    static const std::pair<const char*, double> kTime[] = {{"Second", 1.0}};
    static const std::pair<const char*, double> kTemp[] = {{"Kelvin", 1.0}, {"Rankine", 5.0 / 9.0}};
    static const std::pair<const char*, double> kAngle[] = {
        {"Radian", 1.0}, {"Degree", 3.14159265358979323846 / 180.0}};
    const auto find = [&](const auto& rTable) {
        for (const auto& [n, v] : rTable)
            if (rName == n)
                return v;
        return 0.0;
    };
    switch (Slot) {
        case 0:
            return find(kMass);
        case 1:
            return find(kLength);
        case 2:
            return find(kTime);
        case 3:
            return find(kTemp);
        default:
            return find(kAngle);
    }
}

/// The field unit SIDS nodes state: @p rExponents (Mass, Length, Time,
/// Temperature, Angle), @p rAdditional (ElectricCurrent, Substance,
/// LuminousIntensity), the conversion scale, and the array's DimensionalUnits
/// names (empty strings for SI). False when a slot with a nonzero exponent
/// names a unit without a factor.
inline bool cgns_unit_from_sids(const std::array<double, 5>& rExponents,
                                const std::array<double, 3>& rAdditional, double Scale,
                                const std::array<std::string, 5>& rUnits, FieldUnit& rOut) {
    if (!(Scale > 0.0) || !std::isfinite(Scale))
        return false;
    double factor = Scale;
    for (std::size_t s = 0; s < 5; ++s) {
        if (rExponents[s] == 0.0)
            continue;
        const double f = rUnits[s].empty() ? 1.0 : cgns_unit_factor(s, rUnits[s]);
        if (!(f > 0.0))
            return false;
        factor *= std::pow(f, rExponents[s]);
    }
    rOut.mScale = factor;
    rOut.mExponents = {rExponents[1],  rExponents[0],  rExponents[2], rExponents[3],
                       rAdditional[0], rAdditional[1], rAdditional[2]};
    return true;
}

/// The SIDS exponents (Mass, Length, Time, Temperature, Angle) of @p rUnit.
inline std::array<double, 5> cgns_sids_exponents(const FieldUnit& rUnit) {
    return {rUnit.mExponents[1], rUnit.mExponents[0], rUnit.mExponents[2], rUnit.mExponents[3],
            0.0};
}

/// The SIDS additional exponents (ElectricCurrent, Substance, LuminousIntensity) of @p rUnit.
inline std::array<double, 3> cgns_sids_additional(const FieldUnit& rUnit) {
    return {rUnit.mExponents[4], rUnit.mExponents[5], rUnit.mExponents[6]};
}

}  // namespace detail
}  // namespace meshioplusplus
