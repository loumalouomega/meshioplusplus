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
// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

// Project includes
#include "layer_sweep.hpp"

namespace meshioplusplus::detail {

namespace {

using LsV3 = std::array<double, 3>;

LsV3 ls_at(const double* pX, std::int64_t I) {
    const auto i = static_cast<std::size_t>(I);
    return {pX[3 * i], pX[3 * i + 1], pX[3 * i + 2]};
}

// det(A - O, B - O, C - O) over the three lengths; 0 for a zero edge.
double ls_corner(const LsV3& rO, const LsV3& rA, const LsV3& rB, const LsV3& rC) {
    const LsV3 a{rA[0] - rO[0], rA[1] - rO[1], rA[2] - rO[2]};
    const LsV3 b{rB[0] - rO[0], rB[1] - rO[1], rB[2] - rO[2]};
    const LsV3 c{rC[0] - rO[0], rC[1] - rO[1], rC[2] - rO[2]};
    const double la = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    const double lb = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    const double lc = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    const double l = la * lb * lc;
    if (!(l > 0.0))
        return 0.0;
    const double det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
                       a[2] * (b[0] * c[1] - b[1] * c[0]);
    return det / l;
}

}  // namespace

double layer_cell_scaled_jacobian(const double* pX, const LayerSweepCell& rCell, int Kind) {
    std::array<LsV3, 8> c{};
    for (int i = 0; i < Kind; ++i)
        c[static_cast<std::size_t>(i)] = ls_at(pX, rCell.mNodes[static_cast<std::size_t>(i)]);
    double sj = std::numeric_limits<double>::infinity();
    if (Kind == 8) {
        for (int i = 0; i < 4; ++i) {
            const auto b = static_cast<std::size_t>(i), n = static_cast<std::size_t>((i + 1) % 4),
                       p = static_cast<std::size_t>((i + 3) % 4);
            sj = std::min(sj, ls_corner(c[b], c[n], c[p], c[b + 4]));
            sj = std::min(sj, ls_corner(c[b + 4], c[p + 4], c[n + 4], c[b]));
        }
    } else if (Kind == 6) {
        static constexpr int kT[6][3] = {{1, 2, 3}, {2, 0, 4}, {0, 1, 5},
                                         {5, 4, 0}, {3, 5, 1}, {4, 3, 2}};
        for (std::size_t i = 0; i < 6; ++i)
            sj = std::min(sj, ls_corner(c[i], c[static_cast<std::size_t>(kT[i][0])],
                                        c[static_cast<std::size_t>(kT[i][1])],
                                        c[static_cast<std::size_t>(kT[i][2])]));
    } else if (Kind == 5) {
        for (std::size_t i = 0; i < 4; ++i)
            sj = std::min(sj, ls_corner(c[i], c[(i + 1) % 4], c[(i + 3) % 4], c[4]));
    } else if (Kind == 4) {
        sj = ls_corner(c[0], c[1], c[2], c[3]);
    } else {
        throw std::invalid_argument("meshio++: layer_sweep: a cell has 4, 5, 6 or 8 nodes");
    }
    return sj;
}

LayerSweepResult layer_sweep(const LayerSweepInput& rIn) {
    const std::int64_t n = rIn.mNumPoints;
    const int layers = rIn.mLayers;
    if (layers < 1)
        throw std::invalid_argument("meshio++: layer_sweep: at least one layer");
    const auto nn = static_cast<std::size_t>(n);
    const auto nl = static_cast<std::size_t>(layers);
    if (rIn.mDirections.size() != 3 * nn || rIn.mOffsets.size() != nl * nn)
        throw std::invalid_argument("meshio++: layer_sweep: a direction and offsets per point");
    LayerSweepResult out;
    out.mPoints.assign(rIn.pPoints, rIn.pPoints + 3 * nn);
    out.mOrigin.resize(nn);
    for (std::size_t p = 0; p < nn; ++p)
        out.mOrigin[p] = static_cast<std::int64_t>(p);
    // A point grows in layer K when its offset there differs from the one
    // below (0 at the surface); otherwise its copy at K is the one below, and
    // its corner of every cell of that layer is collapsed. A point that grows
    // in no layer is the surface point at every level.
    std::vector<char> grows(nn * nl, 0);
    for (std::size_t p = 0; p < nn; ++p)
        for (std::size_t k = 0; k < nl; ++k) {
            const double below = k == 0 ? 0.0 : rIn.mOffsets[p * nl + k - 1];
            grows[k * nn + p] = rIn.mOffsets[p * nl + k] != below ? 1 : 0;
        }
    // The id of point P at level K (0 the surface).
    std::vector<std::int64_t> level_id(nn * nl, -1);
    for (std::size_t k = 0; k < nl; ++k)
        for (std::size_t p = 0; p < nn; ++p) {
            if (rIn.mPeriodic && k + 1 == nl) {
                // The last level of a full turn is the surface again.
                level_id[k * nn + p] = static_cast<std::int64_t>(p);
                continue;
            }
            if (!grows[k * nn + p]) {
                level_id[k * nn + p] =
                    k == 0 ? static_cast<std::int64_t>(p) : level_id[(k - 1) * nn + p];
                continue;
            }
            level_id[k * nn + p] = static_cast<std::int64_t>(out.mOrigin.size());
            out.mOrigin.push_back(static_cast<std::int64_t>(p));
            for (std::size_t c = 0; c < 3; ++c)
                out.mPoints.push_back(rIn.pPoints[3 * p + c] +
                                      rIn.mOffsets[p * nl + k] * rIn.mDirections[3 * p + c]);
        }
    const auto at = [&](std::int64_t P, int K) {
        return K == 0
                   ? P
                   : level_id[static_cast<std::size_t>(K - 1) * nn + static_cast<std::size_t>(P)];
    };
    // Collapsed in layer K (between levels K - 1 and K).
    const auto is_collapsed = [&](std::int64_t P, int K) {
        return grows[static_cast<std::size_t>(K - 1) * nn + static_cast<std::size_t>(P)] == 0;
    };
    // A triangle (A, B, C) between levels K - 1 and K: a wedge, or what its
    // collapsed corners leave of it.
    const auto triangle = [&](std::int64_t A, std::int64_t B, std::int64_t C, std::int64_t F,
                              int K) {
        const std::array<std::int64_t, 3> v{A, B, C};
        int ncol = 0;
        for (const std::int64_t x : v)
            ncol += is_collapsed(x, K) ? 1 : 0;
        LayerSweepCell cell;
        cell.mNodes.fill(-1);
        cell.mFace = F;
        cell.mLayer = K;
        if (ncol == 0) {
            for (std::size_t i = 0; i < 3; ++i) {
                cell.mNodes[i] = at(v[i], K - 1);
                cell.mNodes[i + 3] = at(v[i], K);
            }
            out.mWedges.push_back(cell);
            return;
        }
        if (ncol == 3)
            return;
        // Rotate so the collapsed ones come first, in cyclic order.
        std::size_t s = 0;
        if (ncol == 1) {
            while (!is_collapsed(v[s], K))
                ++s;
            const std::int64_t x = v[s], y = v[(s + 1) % 3], z = v[(s + 2) % 3];
            cell.mNodes = {at(y, K - 1), at(y, K), at(z, K), at(z, K - 1),
                           at(x, K - 1), -1,       -1,       -1};
            out.mPyramids.push_back(cell);
        } else {
            while (is_collapsed(v[s], K))
                ++s;
            const std::int64_t z = v[s], x = v[(s + 1) % 3], y = v[(s + 2) % 3];
            cell.mNodes = {at(x, K - 1), at(y, K - 1), at(z, K - 1), at(z, K), -1, -1, -1, -1};
            out.mTetrahedra.push_back(cell);
        }
    };
    for (std::size_t f = 0; f < rIn.mFaces.size(); ++f) {
        const auto& r_f = rIn.mFaces[f];
        const bool quad = r_f[3] >= 0;
        const int corners = quad ? 4 : 3;
        // A quad with a corner collapsed in any layer is split in every
        // layer, so the columns beside it stay conforming.
        bool any = false;
        for (int i = 0; i < corners; ++i)
            for (int k = 1; k <= layers && !any; ++k)
                any = is_collapsed(r_f[static_cast<std::size_t>(i)], k);
        const auto face = static_cast<std::int64_t>(f);
        // An edge collapsed in every layer and the rest in none: corner S
        // and S + 1 on the edge.
        int edge_at = -1;
        if (quad && any && rIn.mCollapsedEdgeWedges)
            for (int s = 0; s < 4 && edge_at < 0; ++s) {
                bool ok = true;
                for (int i = 0; i < 4 && ok; ++i) {
                    const bool on_edge = i == s || i == (s + 1) % 4;
                    for (int k = 1; k <= layers && ok; ++k)
                        ok = is_collapsed(r_f[static_cast<std::size_t>(i)], k) == on_edge;
                }
                if (ok)
                    edge_at = s;
            }
        const bool split = quad && any && edge_at < 0;
        for (int k = 1; k <= layers; ++k) {
            if (edge_at >= 0) {
                const auto s = static_cast<std::size_t>(edge_at);
                const std::int64_t x = r_f[s], y = r_f[(s + 1) % 4], z = r_f[(s + 2) % 4],
                                   w = r_f[(s + 3) % 4];
                LayerSweepCell cell;
                cell.mFace = face;
                cell.mLayer = k;
                cell.mNodes = {x, at(w, k - 1), at(w, k), y, at(z, k - 1), at(z, k), -1, -1};
                out.mWedges.push_back(cell);
            } else if (quad && !split) {
                LayerSweepCell cell;
                cell.mFace = face;
                cell.mLayer = k;
                for (std::size_t i = 0; i < 4; ++i) {
                    cell.mNodes[i] = at(r_f[i], k - 1);
                    cell.mNodes[i + 4] = at(r_f[i], k);
                }
                out.mHexahedra.push_back(cell);
            } else if (split) {
                triangle(r_f[0], r_f[1], r_f[2], face, k);
                triangle(r_f[0], r_f[2], r_f[3], face, k);
            } else {
                triangle(r_f[0], r_f[1], r_f[2], face, k);
            }
        }
        // The wall and the cap, split as their column.
        const auto put = [&](int K, std::vector<LayerSweepCell>& rTris,
                             std::vector<LayerSweepCell>& rQuads, std::int32_t Layer) {
            LayerSweepCell c;
            c.mNodes.fill(-1);
            c.mFace = face;
            c.mLayer = Layer;
            if (quad && !split) {
                for (std::size_t i = 0; i < 4; ++i)
                    c.mNodes[i] = at(r_f[i], K);
                rQuads.push_back(c);
            } else if (split) {
                c.mNodes[0] = at(r_f[0], K);
                c.mNodes[1] = at(r_f[1], K);
                c.mNodes[2] = at(r_f[2], K);
                rTris.push_back(c);
                c.mNodes[1] = at(r_f[2], K);
                c.mNodes[2] = at(r_f[3], K);
                rTris.push_back(c);
            } else {
                for (std::size_t i = 0; i < 3; ++i)
                    c.mNodes[i] = at(r_f[i], K);
                rTris.push_back(c);
            }
        };
        if (!rIn.mPeriodic) {
            put(0, out.mWallTriangles, out.mWallQuads, 0);
            put(layers, out.mCapTriangles, out.mCapQuads, layers + 1);
        }
    }
    // The edges' side faces.
    for (std::size_t e = 0; e < rIn.mEdges.size(); ++e) {
        const std::int64_t a = rIn.mEdges[e][0], b = rIn.mEdges[e][1];
        for (int k = 1; k <= layers; ++k) {
            const bool ca = is_collapsed(a, k), cb = is_collapsed(b, k);
            if (ca && cb)
                continue;
            LayerSweepCell c;
            c.mNodes.fill(-1);
            c.mFace = static_cast<std::int64_t>(e);
            c.mLayer = k;
            if (!ca && !cb) {
                c.mNodes[0] = at(a, k - 1);
                c.mNodes[1] = at(b, k - 1);
                c.mNodes[2] = at(b, k);
                c.mNodes[3] = at(a, k);
                out.mSideQuads.push_back(c);
            } else {
                c.mNodes[0] = at(a, k - 1);
                c.mNodes[1] = at(b, k - 1);
                c.mNodes[2] = ca ? at(b, k) : at(a, k);
                out.mSideTriangles.push_back(c);
            }
        }
    }
    out.mLevelIds = std::move(level_id);
    return out;
}

}  // namespace meshioplusplus::detail
