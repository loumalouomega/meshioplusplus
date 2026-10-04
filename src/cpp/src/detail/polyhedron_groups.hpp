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

#pragma once

/**
 * @file detail/polyhedron_groups.hpp
 * @brief Split a staged polyhedron CSR into `polyhedron<N>` groups, N being a
 * cell's unique node count, in first-seen order.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed
 * header names it. It is the convention the EnSight, CGNS, OpenFOAM and VTU
 * readers share. The input is the CSR triple a reader staged for a whole
 * section -- no vector per cell or per face -- and each group comes back as a
 * CSR triple ready for `AddPolyhedronBlock` (doc/benchmarks.md, "CSR ragged
 * readers and `reorder`").
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// One `polyhedron<N>` group: its cells (indices into the staged section, in
/// section order) and their CSR triple.
struct PolyhedronGroup {
    std::size_t mNodeCount = 0;
    std::vector<std::size_t> mCells;
    std::vector<std::int64_t> mFlat;
    std::vector<std::int64_t> mRowOffsets{0};
    std::vector<std::int64_t> mFaceOffsets{0};
};

/**
 * @brief Group the polyhedra `(Flat, RowOffsets, FaceOffsets)` by unique node
 * count, groups in first-seen order. The arrays are consumed: when one group
 * holds every cell it takes them whole, with no copy.
 * @pre The offsets were checked (`check_csr`) or built by the caller.
 */
inline std::vector<PolyhedronGroup> group_polyhedra_by_node_count(
    std::vector<std::int64_t> Flat, std::vector<std::int64_t> RowOffsets,
    std::vector<std::int64_t> FaceOffsets) {
    const std::size_t ncells = FaceOffsets.empty() ? 0 : FaceOffsets.size() - 1;
    std::vector<std::size_t> node_counts(ncells);
    std::vector<std::int64_t> uniq;
    for (std::size_t c = 0; c < ncells; ++c) {
        const auto first = Flat.begin() + RowOffsets[static_cast<std::size_t>(FaceOffsets[c])];
        const auto last = Flat.begin() + RowOffsets[static_cast<std::size_t>(FaceOffsets[c + 1])];
        uniq.assign(first, last);
        std::sort(uniq.begin(), uniq.end());
        node_counts[c] =
            static_cast<std::size_t>(std::unique(uniq.begin(), uniq.end()) - uniq.begin());
    }
    std::vector<std::size_t> order;
    std::map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t c = 0; c < ncells; ++c) {
        if (groups.find(node_counts[c]) == groups.end())
            order.push_back(node_counts[c]);
        groups[node_counts[c]].push_back(c);
    }
    std::vector<PolyhedronGroup> out;
    out.reserve(order.size());
    for (const std::size_t n : order) {
        PolyhedronGroup g;
        g.mNodeCount = n;
        g.mCells = std::move(groups[n]);
        if (g.mCells.size() == ncells) {
            g.mFlat = std::move(Flat);
            g.mRowOffsets = std::move(RowOffsets);
            g.mFaceOffsets = std::move(FaceOffsets);
        } else {
            for (const std::size_t c : g.mCells) {
                const auto f0 = static_cast<std::size_t>(FaceOffsets[c]);
                const auto f1 = static_cast<std::size_t>(FaceOffsets[c + 1]);
                const auto base = static_cast<std::int64_t>(g.mFlat.size());
                g.mFlat.insert(g.mFlat.end(), Flat.begin() + RowOffsets[f0],
                               Flat.begin() + RowOffsets[f1]);
                for (std::size_t f = f0; f < f1; ++f)
                    g.mRowOffsets.push_back(base + RowOffsets[f + 1] - RowOffsets[f0]);
                g.mFaceOffsets.push_back(static_cast<std::int64_t>(g.mRowOffsets.size() - 1));
            }
        }
        out.push_back(std::move(g));
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
