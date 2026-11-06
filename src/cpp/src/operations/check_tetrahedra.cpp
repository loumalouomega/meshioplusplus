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
// check_tetrahedra. See operations/check_tetrahedra.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/check_tetrahedra.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/ndarray.hpp"

// Project includes (private, not installed)
#include "../detail/tet_checks.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

TetrahedraCheckResult check_tetrahedra(const Mesh& rMesh, const TetrahedraCheckOptions& rOptions) {
    if (!(rOptions.mZeroVolumeFactor >= 0.0) || !std::isfinite(rOptions.mZeroVolumeFactor))
        throw std::invalid_argument(
            "meshio++: check_tetrahedra: zero_volume_factor must be finite and non-negative");
    TetrahedraCheckResult result;

    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    std::vector<detail::Vec3> xyz(rMesh.NumPoints());
    for (std::size_t i = 0; i < xyz.size(); ++i)
        xyz[i] = detail::read_point(points, dim, static_cast<std::int64_t>(i));

    std::vector<detail::TetNodes> tets;
    std::size_t block = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t index = block++;
        const std::string type(cb.Type());
        if (cb.NumCells() == 0 || (type != "tetra" && type != "tetra10"))
            continue;
        if (cb.IsRagged()) {
            result.mWarnings.push_back("block " + std::to_string(index) + " ('" + type +
                                       "') is ragged and was not examined");
            continue;
        }
        const detail::Int64View conn(cb.Conn());
        const std::size_t width = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            detail::TetNodes tet;
            for (std::size_t k = 0; k < 4; ++k) {
                const std::int64_t node = conn[c * width + k];
                if (node < 0 || static_cast<std::size_t>(node) >= xyz.size())
                    throw std::invalid_argument(
                        "meshio++: check_tetrahedra: a cell refers to a point out of range");
                tet[k] = node;
            }
            tets.push_back(tet);
        }
    }

    result.mNumTetrahedra = static_cast<std::int64_t>(tets.size());
    if (tets.empty()) {
        result.mWarnings.push_back("the mesh has no tetrahedra; nothing was examined");
        return result;
    }
    result.mMinVolume = std::numeric_limits<double>::infinity();
    result.mMaxVolume = -std::numeric_limits<double>::infinity();
    for (const detail::TetNodes& tet : tets) {
        const auto [volume, tolerance] =
            detail::tet_volume_and_tolerance(xyz, tet, rOptions.mZeroVolumeFactor);
        result.mMinVolume = std::min(result.mMinVolume, volume);
        result.mMaxVolume = std::max(result.mMaxVolume, volume);
        if (volume < -tolerance)
            ++result.mNumInverted;
        else if (volume <= tolerance)
            ++result.mNumDegenerate;
    }
    result.mNumOverlaps = detail::count_tet_overlaps(xyz, tets);
    result.mValid =
        result.mNumInverted == 0 && result.mNumDegenerate == 0 && result.mNumOverlaps == 0;
    return result;
}

}  // namespace meshioplusplus
