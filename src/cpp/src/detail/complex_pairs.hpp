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
 * @file detail/complex_pairs.hpp
 * @brief Complex-field naming convention (`<base>_real` / `<base>_imag`).
 *
 * Complex fields are stored as matched pairs of real-valued arrays named
 * `<base>_real` and `<base>_imag` (the convention shared by UNV, Nastran OP2,
 * Nastran HDF5 and EnSight). Core-private: no installed header includes this.
 */

// System includes
#include <algorithm>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus::detail {

inline constexpr std::string_view kComplexRealSuffix = "_real";
inline constexpr std::string_view kComplexImagSuffix = "_imag";

/// A matched `<base>_real` / `<base>_imag` pair.
struct ComplexPair {
    std::string mBase;
    std::string mReal;
    std::string mImag;
};

/// Returns every matched `<base>_real` / `<base>_imag` pair in `rNames`,
/// sorted by `mBase`. `<base>` must be non-empty.
inline std::vector<ComplexPair> find_complex_pairs(const std::vector<std::string>& rNames) {
    std::unordered_set<std::string_view> present(rNames.begin(), rNames.end());
    std::vector<ComplexPair> pairs;
    for (const std::string& name : rNames) {
        if (name.size() <= kComplexRealSuffix.size())
            continue;
        if (std::string_view(name).ends_with(kComplexRealSuffix)) {
            const std::string base = name.substr(0, name.size() - kComplexRealSuffix.size());
            const std::string imag = base + std::string(kComplexImagSuffix);
            if (present.count(imag) != 0)
                pairs.push_back({base, name, imag});
        }
    }
    std::sort(pairs.begin(), pairs.end(),
              [](const ComplexPair& rA, const ComplexPair& rB) { return rA.mBase < rB.mBase; });
    return pairs;
}

/// Maps each `<base>_real` name to its matching `ComplexPair` index, and marks
/// every `<base>_imag` partner so writers can emit each pair once in sorted
/// order and skip the companion `_imag` key.
struct ComplexPairIndex {
    std::vector<ComplexPair> mPairs;
    std::unordered_map<std::string, std::size_t> mByReal;
    std::unordered_set<std::string> mImagNames;

    static ComplexPairIndex Build(const std::vector<std::string>& rNames) {
        ComplexPairIndex idx;
        idx.mPairs = find_complex_pairs(rNames);
        for (std::size_t i = 0; i < idx.mPairs.size(); ++i) {
            idx.mByReal.emplace(idx.mPairs[i].mReal, i);
            idx.mImagNames.insert(idx.mPairs[i].mImag);
        }
        return idx;
    }

    const ComplexPair* PairForReal(const std::string& rName) const {
        auto it = mByReal.find(rName);
        return it == mByReal.end() ? nullptr : &mPairs[it->second];
    }

    bool IsImagMember(const std::string& rName) const { return mImagNames.count(rName) != 0; }
};

/// Returns true if `rMesh` carries both `<rBase>_real` and `<rBase>_imag` at `Loc`,
/// writing the two array names to `rReal` and `rImag`.
inline bool find_complex_pair_for_base(const Mesh& rMesh, DataLocation Loc,
                                       const std::string& rBase, std::string& rReal,
                                       std::string& rImag) {
    if (rBase.empty())
        return false;
    const std::string re = rBase + std::string(kComplexRealSuffix);
    const std::string im = rBase + std::string(kComplexImagSuffix);
    bool has_re = false;
    bool has_im = false;
    switch (Loc) {
        case DataLocation::Point:
            has_re = rMesh.HasPointData(re);
            has_im = rMesh.HasPointData(im);
            break;
        case DataLocation::Cell:
            has_re = rMesh.HasCellData(re);
            has_im = rMesh.HasCellData(im);
            break;
        case DataLocation::Field:
            has_re = rMesh.HasFieldData(re);
            has_im = rMesh.HasFieldData(im);
            break;
    }
    if (has_re && has_im) {
        rReal = re;
        rImag = im;
        return true;
    }
    return false;
}

}  // namespace meshioplusplus::detail
