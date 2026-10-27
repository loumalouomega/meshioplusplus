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
 * @file detail/layer_sweep.hpp
 * @brief Layers of cells swept from a surface along a direction per point
 * (roadmap 6.9.2; the sweep `boundary_layer` uses and `extrude` will).
 *
 * Every point `p` of a triangle and quad surface is copied at `L` levels,
 * `x_p + o_{p,k} d_p` for `k = 1..L` with offsets increasing in `k`; between
 * levels `k - 1` and `k` (level 0 the surface itself) each triangle gives a
 * `wedge` and each quad a `hexahedron`, ordered so that its Jacobian is
 * positive when the face's normal (right-hand rule) points along the
 * directions. A point whose offset at level `k` equals the one below (0 at
 * the surface) does not grow in layer `k`: its copy there is the one below,
 * and its corner of every cell of that layer is *collapsed*, so a point
 * whose offsets stop growing after `m` levels grows `m` layers (Garimella and
 * Shephard's fewer layers at a point), and one whose offsets are all zero
 * stays the surface point. A face with a corner collapsed in any layer is
 * split into triangles (a quad along its first diagonal) in every layer, and
 * a wedge with one collapsed corner becomes a `pyramid` on its opposite side
 * face, with two a `tetra`, with three nothing; the faces the cells share
 * stay conforming, since every column through a collapsed point is split
 * the same way in every layer. The surface (`wall`) and the last level
 * (`cap`) come out as faces, split like their columns.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

namespace meshioplusplus::detail {

struct LayerSweepInput {
    /// The surface's points, three coordinates each.
    const double* pPoints = nullptr;
    std::int64_t mNumPoints = 0;
    /// Its faces: a triangle has -1 last.
    std::vector<std::array<std::int64_t, 4>> mFaces;
    /// Per point, the direction (three coordinates each).
    std::vector<double> mDirections;
    /// Per point, its `mLayers` offsets along its direction, level by level,
    /// not decreasing; an offset equal to the one below (0 at the surface)
    /// does not grow that layer, and all zero is a collapsed point.
    std::vector<double> mOffsets;
    int mLayers = 1;
};

/// One cell or face of the sweep: its node ids (unused slots -1), the
/// surface face it comes from and its layer (1..L; 0 the wall, L + 1 the cap).
struct LayerSweepCell {
    std::array<std::int64_t, 8> mNodes{};
    std::int64_t mFace = 0;
    std::int32_t mLayer = 0;
};

struct LayerSweepResult {
    /// The points: the surface's, then the levels' own (three coordinates).
    std::vector<double> mPoints;
    /// Per point, the surface point it is a copy of.
    std::vector<std::int64_t> mOrigin;
    std::vector<LayerSweepCell> mHexahedra, mWedges, mPyramids, mTetrahedra;
    /// The wall and cap faces, triangles and quads apart.
    std::vector<LayerSweepCell> mWallTriangles, mWallQuads, mCapTriangles, mCapQuads;
};

LayerSweepResult layer_sweep(const LayerSweepInput& rIn);

/// The least corner Jacobian of a sweep cell over its corner edges' lengths
/// (1 for a right prism or a cube; at most 0 when it folds): the measure the
/// sweep's caller judges a cell by. `Kind` is 8, 6, 5 or 4 nodes.
double layer_cell_scaled_jacobian(const double* pX, const LayerSweepCell& rCell, int Kind);

}  // namespace meshioplusplus::detail
