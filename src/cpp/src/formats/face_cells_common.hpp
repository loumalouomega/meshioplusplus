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
 * @file formats/face_cells_common.hpp
 * @brief Volume cells from their faces: the kernel shared by the readers of
 *        face-based formats (OpenFOAM `polyMesh`, ANSYS Fluent `.msh`).
 *
 * A cell is given as its faces, each a list of point ids wound so that its
 * right-hand normal points out of the cell. Tetrahedra, pyramids, wedges and
 * hexahedra come back in meshio++'s (VTK) node order with a positive volume --
 * the orientation is checked geometrically, so a flipped input face does not
 * invert the cell; anything else is reported as a polyhedron and the caller
 * keeps its outward faces. `polygon_from_edges` does the 2-D counterpart.
 *
 * A **format-private** header (the `gid_common.hpp` precedent): it sits beside
 * the `.cpp` files, is never installed, and lives in a *named* namespace
 * because the amalgamation concatenates every source into one translation
 * unit. Twin of `src/python/meshioplusplus/_face_cells.py`.
 */

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace meshioplusplus {
namespace face_cells {

using Face = std::vector<std::int64_t>;
using P3 = std::vector<std::array<double, 3>>;

inline double triple(const std::array<double, 3>& rA, const std::array<double, 3>& rB,
                     const std::array<double, 3>& rC) {
    // a . (b x c)
    double cx = rB[1] * rC[2] - rB[2] * rC[1];
    double cy = rB[2] * rC[0] - rB[0] * rC[2];
    double cz = rB[0] * rC[1] - rB[1] * rC[0];
    return rA[0] * cx + rA[1] * cy + rA[2] * cz;
}

inline std::array<double, 3> sub(const std::array<double, 3>& rA, const std::array<double, 3>& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}

inline std::size_t unique_node_count(const std::vector<Face>& rFaces) {
    std::unordered_set<std::int64_t> s;
    for (const auto& f : rFaces)
        for (std::int64_t v : f)
            s.insert(v);
    return s.size();
}

inline std::unordered_map<std::int64_t, std::unordered_set<std::int64_t>> node_adjacency(
    const std::vector<Face>& rFaces) {
    std::unordered_map<std::int64_t, std::unordered_set<std::int64_t>> adj;
    for (const auto& f : rFaces) {
        std::size_t m = f.size();
        for (std::size_t i = 0; i < m; ++i) {
            std::int64_t a = f[i], b = f[(i + 1) % m];
            adj[a].insert(b);
            adj[b].insert(a);
        }
    }
    return adj;
}

// Returns the ordered top ring, or empty if ambiguous.
inline std::vector<std::int64_t> match_top(const Face& rBottom,
                                           const std::vector<Face>& rOriented) {
    auto adj = node_adjacency(rOriented);
    std::unordered_set<std::int64_t> base(rBottom.begin(), rBottom.end());
    std::vector<std::int64_t> top;
    for (std::int64_t b : rBottom) {
        std::vector<std::int64_t> cand;
        for (std::int64_t x : adj[b])
            if (!base.count(x))
                cand.push_back(x);
        if (cand.size() != 1)
            return {};
        top.push_back(cand[0]);
    }
    return top;
}

inline Face build_tetra(const std::vector<Face>& rOriented, const P3& rP) {
    const Face& base = rOriented[0];
    std::unordered_set<std::int64_t> all;
    for (const auto& f : rOriented)
        for (std::int64_t v : f)
            all.insert(v);
    for (std::int64_t v : base)
        all.erase(v);
    if (all.empty())
        return {};
    std::int64_t apex = *all.begin();
    Face n = {base[0], base[1], base[2], apex};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[2]], rP[n[0]]), sub(rP[n[3]], rP[n[0]])) < 0)
        n = {base[0], base[2], base[1], apex};
    return n;
}

inline Face build_pyramid(const std::vector<Face>& rOriented, const P3& rP) {
    Face quad;
    for (const auto& f : rOriented)
        if (f.size() == 4) {
            quad = f;
            break;
        }
    std::unordered_set<std::int64_t> all;
    for (const auto& f : rOriented)
        for (std::int64_t v : f)
            all.insert(v);
    for (std::int64_t v : quad)
        all.erase(v);
    std::int64_t apex = *all.begin();
    Face n = {quad[0], quad[1], quad[2], quad[3], apex};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[3]], rP[n[0]]), sub(rP[n[4]], rP[n[0]])) < 0)
        n = {quad[0], quad[3], quad[2], quad[1], apex};
    return n;
}

inline Face build_wedge(const std::vector<Face>& rOriented, const P3& rP) {
    Face bottom;
    for (const auto& f : rOriented)
        if (f.size() == 3) {
            bottom = f;
            break;
        }
    std::vector<std::int64_t> top = match_top(bottom, rOriented);
    if (top.empty())
        return {};
    Face n = {bottom[0], bottom[1], bottom[2], top[0], top[1], top[2]};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[2]], rP[n[0]]), sub(rP[n[3]], rP[n[0]])) < 0)
        n = {bottom[0], bottom[2], bottom[1], top[0], top[2], top[1]};
    return n;
}

inline Face build_hexahedron(const std::vector<Face>& rOriented, const P3& rP) {
    Face bottom;
    for (const auto& f : rOriented)
        if (f.size() == 4) {
            bottom = f;
            break;
        }
    std::vector<std::int64_t> top = match_top(bottom, rOriented);
    if (top.empty())
        return {};
    Face n = {bottom[0], bottom[1], bottom[2], bottom[3], top[0], top[1], top[2], top[3]};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[3]], rP[n[0]]), sub(rP[n[4]], rP[n[0]])) < 0)
        n = {bottom[0], bottom[3], bottom[2], bottom[1], top[0], top[3], top[2], top[1]};
    return n;
}

// Classify a cell. Returns {meshio type, connectivity}. For "polyhedron" the
// connectivity is empty (the caller keeps the oriented faces).
inline std::pair<std::string, Face> reconstruct_cell(const std::vector<Face>& rOriented,
                                                     const P3& rP) {
    // A face needs three corners; a malformed list is skipped by the caller
    // (an empty connectivity), never indexed.
    for (const auto& f : rOriented)
        if (f.size() < 3)
            return {"invalid", {}};
    std::size_t nf = rOriented.size();
    std::size_t np = unique_node_count(rOriented);
    if (nf == 4 && np == 4)
        return {"tetra", build_tetra(rOriented, rP)};
    if (nf == 5 && np == 5)
        return {"pyramid", build_pyramid(rOriented, rP)};
    if (nf == 5 && np == 6)
        return {"wedge", build_wedge(rOriented, rP)};
    if (nf == 6 && np == 8)
        return {"hexahedron", build_hexahedron(rOriented, rP)};
    return {"polyhedron", {}};
}

// Chain a 2-D cell's boundary edges into one ring of point ids, counter-
// clockwise in the xy plane. Empty when the edges do not close a single loop.
inline Face polygon_from_edges(const std::vector<std::array<std::int64_t, 2>>& rEdges,
                               const P3& rP) {
    if (rEdges.size() < 3)
        return {};
    std::unordered_map<std::int64_t, std::vector<std::int64_t>> next;
    for (const auto& e : rEdges) {
        // Point ids index rP below: one outside it is a malformed face list.
        for (const std::int64_t v : e)
            if (v < 0 || static_cast<std::size_t>(v) >= rP.size())
                return {};
        next[e[0]].push_back(e[1]);
        next[e[1]].push_back(e[0]);
    }
    for (const auto& [node, nbrs] : next)
        if (nbrs.size() != 2)
            return {};
    const std::int64_t start = rEdges[0][0];
    Face ring = {start, rEdges[0][1]};
    while (ring.size() < next.size()) {
        const auto& nbrs = next[ring.back()];
        const std::int64_t c = nbrs[0] != ring[ring.size() - 2] ? nbrs[0] : nbrs[1];
        if (c == start)
            return {};
        ring.push_back(c);
    }
    const auto& last = next[ring.back()];
    if (ring.size() != rEdges.size() || (last[0] != start && last[1] != start))
        return {};
    double area = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto& a = rP[static_cast<std::size_t>(ring[i])];
        const auto& b = rP[static_cast<std::size_t>(ring[(i + 1) % ring.size()])];
        area += a[0] * b[1] - b[0] * a[1];
    }
    if (area < 0)
        std::reverse(ring.begin() + 1, ring.end());
    return ring;
}

}  // namespace face_cells
}  // namespace meshioplusplus
