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
 * @file detail/node_table.hpp
 * @brief The `.node` table shared by the TetGen and Triangle writers.
 *
 * A **core-private** header: no installed header names it, and it adds
 * nothing to the API or the ABI.
 */

// System includes
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/mesh.hpp"
#include "typed_view.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Write a marker/ref value: integral values as integers, else `%.16e`.
 */
inline void node_table_write_marker(std::ostream& rOs, double v) {
    const double r = std::nearbyint(v);
    if (v == r && std::fabs(v) < 9.2e18) {
        rOs << static_cast<std::int64_t>(r);
    } else {
        char buf[40];
        snprintf_c(buf, sizeof(buf), "%.16e", v);
        rOs << buf;
    }
}

/**
 * @brief Split point_data keys into attribute columns and at most one marker
 * column (the first ":ref" key, or else the first key).
 */
inline void node_table_split_keys(const Mesh& rMesh, std::vector<std::string>& rAttrKeys,
                                  std::vector<std::string>& rRefKeys) {
    rAttrKeys = rMesh.PointDataNames();  // sorted: deterministic column order
    if (rAttrKeys.empty())
        return;
    for (const auto& k : rAttrKeys)
        if (k.find(":ref") != std::string::npos) {
            rRefKeys.push_back(k);
            break;
        }
    if (!rRefKeys.empty()) {
        rAttrKeys.erase(std::remove(rAttrKeys.begin(), rAttrKeys.end(), rRefKeys[0]),
                        rAttrKeys.end());
    } else {
        rRefKeys.push_back(rAttrKeys.front());
        rAttrKeys.erase(rAttrKeys.begin());
    }
}

/**
 * @brief Write the `.node` rows: index, @p Dim coordinates and the attribute
 * columns as `%.16e`, then the marker columns through
 * `node_table_write_marker`.
 */
template <int Dim>
void node_table_write_rows(std::ostream& rOs, const Mesh& rMesh,
                           const std::vector<std::string>& rAttrKeys,
                           const std::vector<std::string>& rRefKeys) {
    const NDArray& points = rMesh.Points();
    const std::int64_t np = static_cast<std::int64_t>(rMesh.NumPoints());
    char fbuf[40];
    const DoubleView point_values(points);
    // One view per attribute / reference array, made before the row loop.
    std::vector<std::optional<DoubleView>> attrs(rAttrKeys.size());
    for (std::size_t a = 0; a < rAttrKeys.size(); ++a)
        attrs[a].emplace(rMesh.PointData(rAttrKeys[a]));
    std::vector<std::optional<DoubleView>> refs(rRefKeys.size());
    for (std::size_t a = 0; a < rRefKeys.size(); ++a)
        refs[a].emplace(rMesh.PointData(rRefKeys[a]));
    for (std::int64_t i = 0; i < np; ++i) {
        rOs << i;
        for (int c = 0; c < Dim; ++c) {
            snprintf_c(fbuf, sizeof(fbuf), "%.16e", point_values[i * Dim + c]);
            rOs << " " << fbuf;
        }
        for (const auto& attr : attrs) {
            snprintf_c(fbuf, sizeof(fbuf), "%.16e", (*attr)[i]);
            rOs << " " << fbuf;
        }
        for (const auto& ref : refs) {
            rOs << " ";
            node_table_write_marker(rOs, (*ref)[i]);
        }
        rOs << "\n";
    }
}

}  // namespace detail
}  // namespace meshioplusplus
