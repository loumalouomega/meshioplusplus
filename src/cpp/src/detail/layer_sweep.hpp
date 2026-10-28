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
 * `extrude` and `revolve` (roadmap 3.1.2, 3.1.3) ask for three more things,
 * all off by default so `boundary_layer`'s sweep is unchanged:
 *
 *  - **edges**: each directed edge `a -> b` sweeps to a side face per layer,
 *    `[a, b, b', a']` (a triangle when one end is collapsed, nothing when
 *    both are). An edge that runs along its face's winding gives the face
 *    of the column beside it that faces outward, so the border edges of an
 *    open surface give its side walls, and a curve's own edges its surface.
 *  - **periodic**: the last level is the surface itself again, so a full
 *    turn closes on itself; there is then no wall and no cap.
 *  - **collapsed-edge wedges**: a quad with two adjacent corners collapsed in
 *    every layer, and the other two in none (an edge on a revolution's
 *    axis), gives one wedge per layer instead of a split column; its wall
 *    and cap stay quads.
 *
 * The id of every point at every level is returned too, so a caller can put
 * the levels' points somewhere other than along a straight direction (a
 * rotation) and tell which corners collapsed.
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
    /// Directed edges to sweep into side faces (`mSideQuads`, `mSideTriangles`).
    std::vector<std::array<std::int64_t, 2>> mEdges;
    /// The last level is the surface again: no points of its own, no wall, no cap.
    bool mPeriodic = false;
    /// A quad with an edge collapsed in every layer gives wedges (see above).
    bool mCollapsedEdgeWedges = false;
};

/// One cell or face of the sweep: its node ids (unused slots -1), the
/// surface face (or edge) it comes from and its layer (1..L; 0 the wall, L + 1
/// the cap).
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
    /// The edges' side faces, by edge then layer; `mFace` is the edge index.
    std::vector<LayerSweepCell> mSideQuads, mSideTriangles;
    /// The id of point P at level K (1..L) is `mLevelIds[(K - 1) * n + P]`;
    /// level 0 is P itself. Equal ids at K - 1 and K: P is collapsed in K.
    std::vector<std::int64_t> mLevelIds;
};

LayerSweepResult layer_sweep(const LayerSweepInput& rIn);

/// The least corner Jacobian of a sweep cell over its corner edges' lengths
/// (1 for a right prism or a cube; at most 0 when it folds): the measure the
/// sweep's caller judges a cell by. `Kind` is 8, 6, 5 or 4 nodes.
double layer_cell_scaled_jacobian(const double* pX, const LayerSweepCell& rCell, int Kind);

}  // namespace meshioplusplus::detail
