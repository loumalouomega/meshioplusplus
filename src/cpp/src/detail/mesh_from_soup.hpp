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
 * @file detail/mesh_from_soup.hpp
 * @brief A one-block mesh from point and fixed-corner cell arrays, as the
 * octree meshers hand their working state to the mesh operations.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace detail {

/// A mesh of the 3-D points `rPoints` (any indexable `[p][k]` rows) and one
/// `rType` block of the `K`-corner cells `rCells`, through the uniform mesh API.
template <class TPoints, std::size_t K>
inline Mesh mesh_from_soup(const TPoints& rPoints,
                           const std::vector<std::array<std::int64_t, K>>& rCells,
                           const std::string& rType) {
    Mesh m;
    NDArray xyzs = NDArray::Uninit(DType::Float64, {rPoints.size(), std::size_t{3}});
    for (std::size_t p = 0; p < rPoints.size(); ++p)
        for (std::size_t k = 0; k < 3; ++k)
            xyzs.As<double>()[p * 3 + k] = rPoints[p][k];
    NDArray conn = NDArray::Uninit(DType::Int64, {rCells.size(), K});
    for (std::size_t n = 0; n < rCells.size(); ++n)
        for (std::size_t c = 0; c < K; ++c)
            conn.As<std::int64_t>()[n * K + c] = rCells[n][c];
    m.AssignPoints(std::move(xyzs));
    m.AddCellBlock(rType, std::move(conn));
    return m;
}

}  // namespace detail
}  // namespace meshioplusplus
