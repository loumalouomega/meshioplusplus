//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace meshioplusplus::detail {

/**
 * @brief The cells of an `Ni x Nj x Nk` lattice of points numbered with `i` fastest: hexahedra,
 *        quads or lines by how many of the dimensions are above one (a single point is a
 *        `vertex`). The twin of `_grid.lattice_cells` in the Python package.
 *
 * Connectivity is 0-based over the lattice's own points, `rNodes` per cell.
 */
inline void lattice_cells(std::int64_t Ni, std::int64_t Nj, std::int64_t Nk, std::string& rType,
                          std::vector<std::int64_t>& rConn, std::size_t& rNodes) {
    const auto at = [&](std::int64_t I, std::int64_t J, std::int64_t K) {
        return I + Ni * (J + Nj * K);
    };
    const int big = (Ni > 1) + (Nj > 1) + (Nk > 1);
    rConn.clear();
    if (big == 3) {
        rType = "hexahedron";
        rNodes = 8;
        for (std::int64_t k = 0; k + 1 < Nk; ++k)
            for (std::int64_t j = 0; j + 1 < Nj; ++j)
                for (std::int64_t i = 0; i + 1 < Ni; ++i)
                    for (const auto& c : {std::array<std::int64_t, 3>{i, j, k},
                                          {i + 1, j, k},
                                          {i + 1, j + 1, k},
                                          {i, j + 1, k},
                                          {i, j, k + 1},
                                          {i + 1, j, k + 1},
                                          {i + 1, j + 1, k + 1},
                                          {i, j + 1, k + 1}})
                        rConn.push_back(at(c[0], c[1], c[2]));
        return;
    }
    // Flatten the lattice to the (up to) two axes above one, slowest first.
    std::vector<std::int64_t> dims;
    std::vector<std::int64_t> strides;
    const std::int64_t d[3] = {Ni, Nj, Nk};
    const std::int64_t st[3] = {1, Ni, Ni * Nj};
    for (int a = 2; a >= 0; --a)
        if (d[a] > 1) {
            dims.push_back(d[a]);
            strides.push_back(st[a]);
        }
    if (big == 2) {
        rType = "quad";
        rNodes = 4;
        for (std::int64_t a = 0; a + 1 < dims[0]; ++a)
            for (std::int64_t b = 0; b + 1 < dims[1]; ++b) {
                const auto p = [&](std::int64_t A, std::int64_t B) {
                    return A * strides[0] + B * strides[1];
                };
                rConn.insert(rConn.end(), {p(a, b), p(a, b + 1), p(a + 1, b + 1), p(a + 1, b)});
            }
    } else if (big == 1) {
        rType = "line";
        rNodes = 2;
        for (std::int64_t a = 0; a + 1 < dims[0]; ++a)
            rConn.insert(rConn.end(), {a * strides[0], (a + 1) * strides[0]});
    } else {
        rType = "vertex";
        rNodes = 1;
        rConn.push_back(0);
    }
}

}  // namespace meshioplusplus::detail
