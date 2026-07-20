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
 * @file detail/space_filling.hpp
 * @brief Space-filling-curve key functions (Morton, Hilbert) shared by the
 * `reorder` and `partition` operations.
 *
 * These are the pure integer transforms only: a quantized 3D point (21 bits
 * per axis) to a 63-bit scalar curve distance. The Mesh-coupled parts —
 * bounding-box computation and coordinate quantization — stay with each
 * operation (reorder quantizes node coordinates, partition cell centroids)
 * so that reorder's output remains byte-identical to what it was before the
 * hoist. Both key functions are locality-preserving bijections; the absolute
 * distance values need not match any external convention, they only need to
 * be stable, which the tests pin.
 */

// System includes
#include <cstdint>

namespace meshioplusplus {
namespace detail {

/// Quantization width per axis: 3 * 21 = 63 bits fit in a uint64 key.
inline constexpr int sfc_bits = 21;

/// Spread the low 21 bits of `x` so bit i lands at position 3*i (Morton).
inline std::uint64_t sfc_part1by2(std::uint64_t x) {
    x &= 0x1fffffULL;
    x = (x | (x << 32)) & 0x1f00000000ffffULL;
    x = (x | (x << 16)) & 0x1f0000ff0000ffULL;
    x = (x | (x << 8)) & 0x100f00f00f00f00fULL;
    x = (x | (x << 4)) & 0x10c30c30c30c30c3ULL;
    x = (x | (x << 2)) & 0x1249249249249249ULL;
    return x;
}

/// Morton (Z-order) key of a quantized 3D point.
inline std::uint64_t sfc_morton_key(const std::uint32_t q[3]) {
    return sfc_part1by2(q[0]) | (sfc_part1by2(q[1]) << 1) | (sfc_part1by2(q[2]) << 2);
}

/// Hilbert distance of a 3D quantized point via Skilling's AxesToTranspose
/// transform (an in-place Gray-code + rotation), then interleaving the
/// transpose columns MSB-first into a scalar distance.
inline std::uint64_t sfc_hilbert_key(const std::uint32_t q[3], int bits) {
    std::uint32_t X[3] = {q[0], q[1], q[2]};
    const int n = 3;
    std::uint32_t M = 1u << (bits - 1);
    std::uint32_t P, Q, t;
    int i;
    // Inverse undo excess work.
    for (Q = M; Q > 1; Q >>= 1) {
        P = Q - 1;
        for (i = 0; i < n; i++) {
            if (X[i] & Q) {
                X[0] ^= P;  // invert
            } else {
                t = (X[0] ^ X[i]) & P;  // exchange
                X[0] ^= t;
                X[i] ^= t;
            }
        }
    }
    // Gray encode.
    for (i = 1; i < n; i++)
        X[i] ^= X[i - 1];
    t = 0;
    for (Q = M; Q > 1; Q >>= 1)
        if (X[n - 1] & Q)
            t ^= Q - 1;
    for (i = 0; i < n; i++)
        X[i] ^= t;
    // Interleave transpose columns, most-significant bit first.
    std::uint64_t d = 0;
    for (int b = bits - 1; b >= 0; --b)
        for (i = 0; i < n; ++i)
            d = (d << 1) | static_cast<std::uint64_t>((X[i] >> b) & 1u);
    return d;
}

}  // namespace detail
}  // namespace meshioplusplus
