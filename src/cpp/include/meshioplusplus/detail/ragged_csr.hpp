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
 * @file detail/ragged_csr.hpp
 * @brief The CSR form every mesh backend stores ragged (polygon, polyhedron)
 * cell blocks in, and the conversions and checks the backends share.
 *
 * A polygon block is `(flat, rowOffsets)`: row `i`'s node ids are
 * `flat[rowOffsets[i] .. rowOffsets[i+1])`. A polyhedron block adds
 * `faceOffsets`: cell `c`'s faces are rows `faceOffsets[c] ..
 * faceOffsets[c+1])`. Readers that build these arrays directly hand them to
 * `AddPolygonBlock`/`AddPolyhedronBlock` with no per-cell allocation; the
 * nested-vector overloads convert through `csr_from_rows`/`csr_from_cells`.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace meshioplusplus {
namespace detail {

/**
 * @brief Check that @p rOffsets describes a buffer of @p NumValues values:
 * it starts at 0, never decreases and ends at @p NumValues.
 * @throws std::invalid_argument naming @p pWhat otherwise.
 */
inline void check_csr(std::size_t NumValues, const std::vector<std::int64_t>& rOffsets,
                      const char* pWhat) {
    const auto fail = [&](const std::string& rWhy) {
        throw std::invalid_argument(std::string("meshio++: ragged block: ") + pWhat + " offsets " +
                                    rWhy);
    };
    if (rOffsets.empty() || rOffsets.front() != 0)
        fail("must start at 0");
    for (std::size_t i = 1; i < rOffsets.size(); ++i)
        if (rOffsets[i] < rOffsets[i - 1])
            fail("must not decrease");
    if (static_cast<std::uint64_t>(rOffsets.back()) != NumValues)
        fail("must end at the " + std::to_string(NumValues) + " values they index, not at " +
             std::to_string(rOffsets.back()));
}

/// Flatten polygon rows into `(rFlat, rRowOffsets)`.
inline void csr_from_rows(const std::vector<std::vector<std::int64_t>>& rRows,
                          std::vector<std::int64_t>& rFlat,
                          std::vector<std::int64_t>& rRowOffsets) {
    std::size_t total = 0;
    for (const auto& r_row : rRows)
        total += r_row.size();
    rFlat.clear();
    rFlat.reserve(total);
    rRowOffsets.assign(1, 0);
    rRowOffsets.reserve(rRows.size() + 1);
    for (const auto& r_row : rRows) {
        rFlat.insert(rFlat.end(), r_row.begin(), r_row.end());
        rRowOffsets.push_back(static_cast<std::int64_t>(rFlat.size()));
    }
}

/// Flatten polyhedron cells (lists of faces) into `(rFlat, rRowOffsets, rFaceOffsets)`.
inline void csr_from_cells(const std::vector<std::vector<std::vector<std::int64_t>>>& rCells,
                           std::vector<std::int64_t>& rFlat, std::vector<std::int64_t>& rRowOffsets,
                           std::vector<std::int64_t>& rFaceOffsets) {
    std::size_t nfaces = 0, total = 0;
    for (const auto& r_cell : rCells) {
        nfaces += r_cell.size();
        for (const auto& r_face : r_cell)
            total += r_face.size();
    }
    rFlat.clear();
    rFlat.reserve(total);
    rRowOffsets.assign(1, 0);
    rRowOffsets.reserve(nfaces + 1);
    rFaceOffsets.assign(1, 0);
    rFaceOffsets.reserve(rCells.size() + 1);
    for (const auto& r_cell : rCells) {
        for (const auto& r_face : r_cell) {
            rFlat.insert(rFlat.end(), r_face.begin(), r_face.end());
            rRowOffsets.push_back(static_cast<std::int64_t>(rFlat.size()));
        }
        rFaceOffsets.push_back(static_cast<std::int64_t>(rRowOffsets.size() - 1));
    }
}

/**
 * @brief Append a verbatim copy of the ragged block viewed by @p rCells (a mesh
 * backend's `CellView`) to @p rOut: CSR to CSR, with no per-cell vector.
 */
template <class CellViewT, class MeshT>
void append_ragged_copy(const CellViewT& rCells, MeshT& rOut) {
    const std::size_t nc = rCells.NumCells();
    std::vector<std::int64_t> flat, rows{0};
    if (rCells.IsPolyhedron()) {
        std::vector<std::int64_t> faces{0};
        faces.reserve(nc + 1);
        for (std::size_t c = 0; c < nc; ++c) {
            for (std::size_t f = 0; f < rCells.NumFaces(c); ++f) {
                const auto face = rCells.Face(c, f);
                flat.insert(flat.end(), face.first, face.first + face.second);
                rows.push_back(static_cast<std::int64_t>(flat.size()));
            }
            faces.push_back(static_cast<std::int64_t>(rows.size() - 1));
        }
        rOut.AddPolyhedronBlock(std::string(rCells.Type()), std::move(flat), std::move(rows),
                                std::move(faces));
        return;
    }
    rows.reserve(nc + 1);
    for (std::size_t c = 0; c < nc; ++c) {
        flat.insert(flat.end(), rCells.Row(c), rCells.Row(c) + rCells.RowSize(c));
        rows.push_back(static_cast<std::int64_t>(flat.size()));
    }
    rOut.AddPolygonBlock(std::string(rCells.Type()), std::move(flat), std::move(rows));
}

}  // namespace detail
}  // namespace meshioplusplus
