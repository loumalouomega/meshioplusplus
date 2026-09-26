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
// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/detail/degenerate_solid.hpp"

namespace meshioplusplus {
namespace detail {

CollapsedBrick collapse_brick(const std::array<std::int64_t, 8>& rNodes) {
    const std::array<std::int64_t, 8>& n = rNodes;
    if (n[3] == n[4] && n[4] == n[5] && n[5] == n[6] && n[6] == n[7])
        return {"tetra", {n[0], n[1], n[2], n[3]}};
    if (n[2] == n[3] && n[4] == n[5] && n[5] == n[6] && n[6] == n[7])
        return {"tetra", {n[0], n[1], n[2], n[4]}};
    if (n[4] == n[5] && n[5] == n[6] && n[6] == n[7])
        return {"pyramid", {n[0], n[1], n[2], n[3], n[4]}};
    if (n[2] == n[3] && n[6] == n[7])
        return {"wedge", {n[0], n[1], n[2], n[4], n[5], n[6]}};
    if (n[4] == n[5] && n[6] == n[7])
        return {"wedge", {n[0], n[4], n[1], n[3], n[6], n[2]}};
    // One side edge collapsed in both the bottom and the top face (Radioss
    // writes `1 2 3 1 5 6 7 5`): the wedge keeps the faces' cyclic order.
    for (std::size_t i = 0; i < 4; ++i) {
        const std::size_t j = (i + 1) % 4, k = (i + 2) % 4, l = (i + 3) % 4;
        if (n[i] == n[j] && n[i + 4] == n[j + 4] && n[j] != n[k] && n[k] != n[l] && n[l] != n[j])
            return {"wedge", {n[j], n[k], n[l], n[j + 4], n[k + 4], n[l + 4]}};
    }
    return {"hexahedron", std::vector<std::int64_t>(n.begin(), n.end())};
}

std::array<std::int64_t, 8> expand_brick(std::string_view Type, const std::int64_t* pRow) {
    const std::int64_t* r = pRow;
    if (Type == "tetra")
        return {r[0], r[1], r[2], r[3], r[3], r[3], r[3], r[3]};
    if (Type == "pyramid")
        return {r[0], r[1], r[2], r[3], r[4], r[4], r[4], r[4]};
    if (Type == "wedge")
        return {r[0], r[2], r[5], r[3], r[1], r[1], r[4], r[4]};
    return {r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]};
}

}  // namespace detail
}  // namespace meshioplusplus
