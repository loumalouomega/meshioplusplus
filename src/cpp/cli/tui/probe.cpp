// SPDX-License-Identifier: MIT
/// @file probe.cpp
/// @brief Implementation of probe.hpp.

#include "probe.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus::cli::tui {

namespace {

// Four significant digits, in the C locale.
std::string number(double Value) {
    char buffer[40];
    detail::snprintf_c(buffer, sizeof(buffer), "%.4g", Value);
    return buffer;
}

// Values per row of a data array: 1 for a scalar, the product of the trailing
// dimensions otherwise.
std::size_t components(const NDArray& rArray) {
    std::size_t k = 1;
    for (std::size_t d = 1; d < rArray.Ndim(); ++d)
        k *= rArray.Shape()[d];
    return std::max<std::size_t>(k, 1);
}

double magnitude(const NDArray& rArray, std::size_t Row, std::size_t K) {
    double sum = 0.0;
    for (std::size_t c = 0; c < K; ++c) {
        const double v = detail::read_double(rArray, Row * K + c);
        sum += v * v;
    }
    return std::sqrt(sum);
}

}  // namespace

Probe probe_cell(const Mesh& rMesh, std::int64_t Cell) {
    Probe probe;
    if (Cell < 0)
        return probe;
    std::int64_t base = 0;
    std::size_t block = 0;
    bool found = false;
    std::string type;
    std::vector<std::int64_t> nodes;
    bool polyhedron = false;
    std::size_t local = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t n = static_cast<std::int64_t>(cb.NumCells());
        if (Cell < base + n) {
            found = true;
            local = static_cast<std::size_t>(Cell - base);
            type = cb.Type();
            if (cb.IsPolyhedron()) {
                polyhedron = true;
            } else if (cb.IsRagged()) {
                const std::int64_t* row = cb.Row(local);
                nodes.assign(row, row + cb.RowSize(local));
            } else {
                const std::size_t k = cb.NodesPerCell();
                for (std::size_t j = 0; j < k; ++j)
                    nodes.push_back(detail::read_int(cb.Conn(), local * k + j));
            }
            break;
        }
        base += n;
        ++block;
    }
    if (!found)
        return probe;
    probe.mCell = Cell;

    std::string head = "cell " + std::to_string(Cell) + "  block " + std::to_string(block) + " " +
                       type + " #" + std::to_string(local);
    if (polyhedron) {
        head += "  (a polyhedron)";
    } else {
        head += "  nodes";
        for (std::size_t i = 0; i < nodes.size() && i < 8; ++i)
            head += " " + std::to_string(nodes[i]);
        if (nodes.size() > 8)
            head += " ... (" + std::to_string(nodes.size()) + ")";
    }
    probe.mLines.push_back(head);

    // Cell data: this cell's value in every array.
    std::string cell_line;
    for (const std::string& name : rMesh.CellDataNames()) {
        const NDArray& array = rMesh.CellData(name, block);
        const std::size_t k = components(array);
        if ((local + 1) * k > array.Size())
            continue;
        cell_line += (cell_line.empty() ? "" : "  ") + name + "=";
        if (k == 1) {
            const double v = detail::read_double(array, local);
            cell_line += number(v);
            probe.mValues[name] = v;
        } else {
            cell_line += "(";
            for (std::size_t c = 0; c < k && c < 6; ++c)
                cell_line += (c ? ", " : "") + number(detail::read_double(array, local * k + c));
            cell_line += k > 6 ? ", ...)" : ")";
            probe.mValues[name] = magnitude(array, local, k);
        }
    }
    if (!cell_line.empty())
        probe.mLines.push_back("cell data: " + cell_line);

    // Point data across the nodes: the mean (and the spread of a scalar).
    std::string point_line;
    if (!nodes.empty()) {
        for (const std::string& name : rMesh.PointDataNames()) {
            const NDArray& array = rMesh.PointData(name);
            const std::size_t k = components(array);
            double lo = 0.0, hi = 0.0, sum = 0.0;
            std::size_t count = 0;
            for (std::int64_t node : nodes) {
                if (node < 0 || (static_cast<std::size_t>(node) + 1) * k > array.Size())
                    continue;
                const double v = k == 1 ? detail::read_double(array, static_cast<std::size_t>(node))
                                        : magnitude(array, static_cast<std::size_t>(node), k);
                lo = count == 0 ? v : std::min(lo, v);
                hi = count == 0 ? v : std::max(hi, v);
                sum += v;
                ++count;
            }
            if (count == 0)
                continue;
            const double mean = sum / static_cast<double>(count);
            probe.mValues[name] = mean;
            point_line += (point_line.empty() ? "" : "  ") + name + (k == 1 ? " " : " |.| ") +
                          number(lo) + ".." + number(hi) + " (mean " + number(mean) + ")";
        }
    }
    if (!point_line.empty())
        probe.mLines.push_back("at the nodes: " + point_line);

    // The named cell regions that hold the cell.
    std::string regions;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const meshioplusplus::Region& region = rMesh.Region(r);
        if (region.mKind != RegionKind::Cell)
            continue;
        for (std::size_t e = 0; e < region.mEntries.Size(); ++e)
            if (detail::read_int(region.mEntries, e) == Cell) {
                regions += (regions.empty() ? "" : ", ") + region.mName;
                break;
            }
    }
    if (!regions.empty())
        probe.mLines.push_back("regions: " + regions);
    return probe;
}

std::string probe_difference(const Probe& rA, const Probe& rB) {
    std::string out;
    for (const auto& [name, a] : rA.mValues) {
        const auto it = rB.mValues.find(name);
        if (it == rB.mValues.end())
            continue;
        const double d = it->second - a;
        out += (out.empty() ? "" : ", ") + name + " " + (d >= 0.0 ? "+" : "") + number(d);
    }
    return out;
}

}  // namespace meshioplusplus::cli::tui
