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
 * @file detail/sfc_quantize.hpp
 * @brief The grid quantization ahead of a space-filling-curve key, shared by
 * `partition` and `reorder`: per-axis scales from a bounding box, and the
 * clamped rounding of a scaled coordinate.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <cstddef>
#include <cstdint>

namespace meshioplusplus {
namespace detail {

/// Per axis `d < Dim`, the scale mapping `[pLo[d], pHi[d]]` onto
/// `[0, 2^Bits - 1]`; 0 on a flat axis. `pScale` holds three values, the
/// axes from `Dim` on left as they are.
inline void sfc_axis_scales(const double* pLo, const double* pHi, std::size_t Dim, int Bits,
                            double* pScale) {
    const double qmax = static_cast<double>((std::int64_t(1) << Bits) - 1);
    for (std::size_t d = 0; d < Dim; ++d) {
        const double span = pHi[d] - pLo[d];
        pScale[d] = span > 0.0 ? qmax / span : 0.0;
    }
}

/// The grid index of a scaled coordinate `T`: rounded, clamped to
/// `[0, 2^Bits - 1]`.
inline std::uint32_t sfc_quantize(double T, int Bits) {
    std::int64_t qi = static_cast<std::int64_t>(T + 0.5);
    if (qi < 0)
        qi = 0;
    const std::int64_t maxq = (std::int64_t(1) << Bits) - 1;
    if (qi > maxq)
        qi = maxq;
    return static_cast<std::uint32_t>(qi);
}

}  // namespace detail
}  // namespace meshioplusplus
