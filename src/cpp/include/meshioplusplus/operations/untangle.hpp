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
 * @file operations/untangle.hpp
 * @brief Untangling inverted linear volume cells (roadmap 7.4.2): move nodes
 * until every tetrahedron and every hexahedron corner has a positive signed
 * volume, keeping the connectivity. Meshes come out tangled from large
 * morphs, from explicit dynamics result files, and from hand edits;
 * `smooth`'s inversion guard only keeps a valid mesh valid.
 *
 * `untangle(mesh, options)` minimises, node by node, the untangling and
 * smoothing functional of J. M. Escobar, E. Rodriguez, R. Montenegro, G.
 * Montero and J. M. Gonzalez-Yuste ("Simultaneous untangling and smoothing of
 * tetrahedral meshes", Comput. Methods Appl. Mech. Engrg. 192, 2003): every
 * corner simplex `c` (a tetrahedron itself; each of a hexahedron's eight
 * corners with its three edges, P. M. Knupp's corner Jacobians) contributes
 * `|S|^2 / (3 h(det S)^(2/3))`, `S` the corner's edge matrix against the ideal
 * one (a regular tetrahedron, a cube corner) in units of the mesh's mean edge
 * length, and `h(s) = (s + sqrt(s^2 + 4 delta^2)) / 2` a regularised signed
 * volume that stays positive while `s` is negative (`delta` from the most
 * negative corner of the sweep, 0 once none is), so the functional is finite
 * across inversions and pushes every corner positive. A sweep visits the free
 * nodes in ascending order (Gauss-Seidel), each taking up to four steps of
 * steepest descent (a central-difference gradient, a backtracking line search
 * that accepts only a decrease).
 *
 * By default only the nodes within `mRings` edge rings of an inverted cell
 * move (recomputed each sweep), and it stops once no cell is inverted;
 * `mSmooth` moves every free node until the functional settles, smoothing as
 * it untangles. The boundary nodes (of faces used by one cell) are fixed by
 * default, `mFrozen` pins more, and so are the nodes of any other volume
 * cell (`wedge`, `pyramid`, higher order), which it does not measure. Points
 * move; connectivity, data and regions pass through. Serial, so the result is
 * the same on every backend.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// How `untangle` moves the nodes.
struct UntangleOptions {
    /// Keep the boundary nodes fixed.
    bool mFixBoundary = true;
    /// More nodes to keep fixed: empty, or one entry per point (nonzero pins).
    std::vector<std::uint8_t> mFrozen;
    /// Edge rings around the inverted cells whose nodes may move.
    std::int64_t mRings = 2;
    /// Move every free node until the functional settles (untangle and smooth).
    bool mSmooth = false;
    /// Sweeps at most; at least 1.
    std::int64_t mMaxSweeps = 100;
};

/// What `untangle` did.
struct UntangleResult {
    /// The mesh with its nodes moved.
    Mesh mMesh;
    /// Tetrahedra and hexahedra with a corner of non-positive volume, before
    /// and after.
    std::int64_t mNumInvertedBefore = 0;
    std::int64_t mNumInvertedAfter = 0;
    /// The smallest corner scaled Jacobian (determinant over the edge
    /// lengths' product), before and after.
    double mMinQualityBefore = 0.0;
    double mMinQualityAfter = 0.0;
    /// Sweeps run, and node moves kept.
    std::int64_t mSweeps = 0;
    std::int64_t mMoves = 0;
};

/**
 * @brief Untangle a mesh's tetrahedra and hexahedra (see the file comment).
 * @throws std::invalid_argument on a mesh with no tetrahedron or hexahedron,
 *         a `mFrozen` mask of the wrong length, `mMaxSweeps` below 1 or a
 *         negative `mRings`.
 */
MESHIOPLUSPLUS_API UntangleResult untangle(const Mesh& rMesh, const UntangleOptions& rOptions = {});

}  // namespace meshioplusplus
