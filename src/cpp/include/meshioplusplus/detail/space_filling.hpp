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
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

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

/**
 * @brief Stable argsort of curve keys: the indices `0 .. n-1` ordered by
 * `rKeys[i]`, ties in ascending index -- exactly what `std::stable_sort` with
 * an indirect `rKeys[a] < rKeys[b]` comparator returns, as a least-significant-
 * digit radix sort (11-bit digits; a digit every key shares is skipped). It
 * replaces that comparison sort in `reorder` and `partition` (roadmap §4):
 * linear passes over two index buffers instead of O(n log n) indirect loads.
 */
inline std::vector<std::int64_t> sfc_stable_argsort(const std::vector<std::uint64_t>& rKeys) {
    const std::size_t n = rKeys.size();
    std::vector<std::int64_t> order(n);
    std::iota(order.begin(), order.end(), std::int64_t{0});
    if (n < 2)
        return order;
    constexpr int kDigitBits = 11;
    constexpr std::size_t kBuckets = std::size_t{1} << kDigitBits;
    std::vector<std::int64_t> next(n);
    std::vector<std::size_t> count(kBuckets);
    for (int shift = 0; shift < 64; shift += kDigitBits) {
        std::fill(count.begin(), count.end(), 0);
        for (std::uint64_t k : rKeys)
            ++count[(k >> shift) & (kBuckets - 1)];
        if (count[(rKeys[0] >> shift) & (kBuckets - 1)] == n)
            continue;  // every key has this digit: the pass would not move anything
        std::size_t sum = 0;
        for (std::size_t& c : count) {
            const std::size_t here = c;
            c = sum;
            sum += here;
        }
        for (std::int64_t idx : order)
            next[count[(rKeys[static_cast<std::size_t>(idx)] >> shift) & (kBuckets - 1)]++] = idx;
        order.swap(next);
    }
    return order;
}

}  // namespace detail
}  // namespace meshioplusplus
