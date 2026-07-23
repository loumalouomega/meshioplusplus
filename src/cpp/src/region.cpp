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
// The Region kind vocabulary and the canonicalization that makes region
// equality exact across backends and thread counts. See region.hpp.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {

const char* region_kind_name(RegionKind Kind) {
    switch (Kind) {
        case RegionKind::Point:
            return "point";
        case RegionKind::Cell:
            return "cell";
        case RegionKind::Side:
            return "side";
    }
    return "point";
}

RegionKind region_kind_from_name(const std::string& rName) {
    if (rName == "point")
        return RegionKind::Point;
    if (rName == "cell")
        return RegionKind::Cell;
    if (rName == "side")
        return RegionKind::Side;
    throw std::invalid_argument("meshio++: unknown region kind '" + rName +
                                "' (expected 'point', 'cell' or 'side')");
}

void Region::Canonicalize() {
    const std::size_t stride = Stride();

    // Widen whatever integer dtype arrived (and tolerate a view over
    // Python-owned memory) into an owned Int64 buffer we may sort in place.
    const std::size_t n_values = mEntries.Size();
    std::vector<std::int64_t> flat(n_values);
    for (std::size_t i = 0; i < n_values; ++i)
        flat[i] = detail::read_int(mEntries, i);

    // A trailing partial row cannot be interpreted; drop it rather than read
    // past the end. (Only reachable through a hand-built malformed region.)
    const std::size_t n_rows = stride == 0 ? 0 : flat.size() / stride;

    if (stride == 1) {
        flat.resize(n_rows);
        std::sort(flat.begin(), flat.end());
        flat.erase(std::unique(flat.begin(), flat.end()), flat.end());
    } else {
        std::vector<std::pair<std::int64_t, std::int64_t>> rows(n_rows);
        for (std::size_t r = 0; r < n_rows; ++r)
            rows[r] = {flat[r * 2], flat[r * 2 + 1]};
        std::sort(rows.begin(), rows.end());
        rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
        flat.resize(rows.size() * 2);
        for (std::size_t r = 0; r < rows.size(); ++r) {
            flat[r * 2] = rows[r].first;
            flat[r * 2 + 1] = rows[r].second;
        }
    }

    const std::size_t out_rows = stride == 0 ? 0 : flat.size() / stride;
    std::vector<std::size_t> shape;
    if (stride == 1)
        shape = {out_rows};
    else
        shape = {out_rows, stride};

    NDArray out = NDArray::Uninit(DType::Int64, std::move(shape));
    if (!flat.empty())
        std::copy(flat.begin(), flat.end(), out.As<std::int64_t>());
    mEntries = std::move(out);
}

bool regions_equal(const Region& rA, const Region& rB) {
    if (rA.mName != rB.mName || rA.mKind != rB.mKind || rA.mDim != rB.mDim || rA.mTag != rB.mTag)
        return false;
    if (rA.mEntries.Size() != rB.mEntries.Size())
        return false;
    const std::size_t n = rA.mEntries.Size();
    for (std::size_t i = 0; i < n; ++i)
        if (detail::read_int(rA.mEntries, i) != detail::read_int(rB.mEntries, i))
            return false;
    return true;
}

}  // namespace meshioplusplus
