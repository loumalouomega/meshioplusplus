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
// A UV atlas of a triangle surface: charts cut and flattened by fit_cad's
// kernels, shelf-packed into the unit square. See operations/parameterize.hpp
// for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/parameterize.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "../detail/cut_graph.hpp"
#include "../detail/param_disc.hpp"
#include "../detail/patch_flatten.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kPrmOp = "meshio++: parameterize";
// Pieces larger than this are halved without trying to flatten them first.
constexpr std::size_t kPrmMaxFlattenTriangles = 8192;

// One flattened chart: its triangles (global triangle indices), its points
// (input point ids) and their plane coordinates.
struct PrmChart {
    std::vector<std::size_t> mTris;
    std::vector<std::int64_t> mPoints;
    std::vector<double> mU;
    std::vector<double> mV;
    std::size_t mFlips = 0;
    double mW = 0.0;  // bounding box once aligned and moved to the origin
    double mH = 0.0;
    double mX = 0.0;  // placement in the atlas
    double mY = 0.0;
};

// A piece of the surface while it is cut: triangle indices into the global list.
struct PrmPiece {
    std::vector<std::size_t> mTris;
    std::int64_t mDepth = 0;
};

using PrmTri = std::array<std::int64_t, 3>;

double prm_tri_area3(const std::vector<detail::PatchVec>& rPts, const detail::PatchTri& rT) {
    const detail::PatchVec& a = rPts[rT[0]];
    const detail::PatchVec& b = rPts[rT[1]];
    const detail::PatchVec& c = rPts[rT[2]];
    const double e1[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double e2[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    const double n[3] = {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                         e1[0] * e2[1] - e1[1] * e2[0]};
    return 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
}

double prm_signed_area2(const std::vector<double>& rU, const std::vector<double>& rV,
                        const detail::PatchTri& rT) {
    return 0.5 * ((rU[rT[1]] - rU[rT[0]]) * (rV[rT[2]] - rV[rT[0]]) -
                  (rU[rT[2]] - rU[rT[0]]) * (rV[rT[1]] - rV[rT[0]]));
}

// The local (compact) points and triangles of a piece.
void prm_local(const std::vector<detail::PatchVec>& rXyz, const std::vector<PrmTri>& rTris,
               const std::vector<std::size_t>& rPiece, std::vector<std::int64_t>& rGlobal,
               std::vector<detail::PatchVec>& rPts, std::vector<detail::PatchTri>& rLocal) {
    rGlobal.clear();
    rPts.clear();
    rLocal.clear();
    for (std::size_t t : rPiece)
        rGlobal.insert(rGlobal.end(), rTris[t].begin(), rTris[t].end());
    std::sort(rGlobal.begin(), rGlobal.end());
    rGlobal.erase(std::unique(rGlobal.begin(), rGlobal.end()), rGlobal.end());
    for (std::int64_t p : rGlobal)
        rPts.push_back(rXyz[static_cast<std::size_t>(p)]);
    auto local_of = [&](std::int64_t P) {
        return static_cast<std::size_t>(std::lower_bound(rGlobal.begin(), rGlobal.end(), P) -
                                        rGlobal.begin());
    };
    for (std::size_t t : rPiece)
        rLocal.push_back({local_of(rTris[t][0]), local_of(rTris[t][1]), local_of(rTris[t][2])});
}

// One chart per triangle, laid out exactly: a fallback for a piece that has no
// map at all.
void prm_single_triangles(const std::vector<detail::PatchVec>& rXyz,
                          const std::vector<PrmTri>& rTris, const std::vector<std::size_t>& rPiece,
                          std::vector<PrmChart>& rCharts) {
    for (std::size_t t : rPiece) {
        PrmChart chart;
        chart.mTris = {t};
        const detail::PatchVec& a = rXyz[static_cast<std::size_t>(rTris[t][0])];
        const detail::PatchVec& b = rXyz[static_cast<std::size_t>(rTris[t][1])];
        const detail::PatchVec& c = rXyz[static_cast<std::size_t>(rTris[t][2])];
        double e1[3], e2[3];
        for (std::size_t k = 0; k < 3; ++k) {
            e1[k] = b[k] - a[k];
            e2[k] = c[k] - a[k];
        }
        const double l1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
        const double x = l1 > 0.0 ? (e2[0] * e1[0] + e2[1] * e1[1] + e2[2] * e1[2]) / l1 : 0.0;
        const double l2sq = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
        const double y = std::sqrt(std::max(0.0, l2sq - x * x));
        // A triangle used twice by a corner (degenerate) keeps distinct points.
        chart.mPoints = {rTris[t][0], rTris[t][1], rTris[t][2]};
        chart.mU = {0.0, l1, x};
        chart.mV = {0.0, 0.0, y};
        chart.mFlips = y > 0.0 && l1 > 0.0 ? 0 : 1;
        rCharts.push_back(std::move(chart));
    }
}

// Sander et al.'s area-weighted L2 stretch of one triangle's map, as its
// squared value times its area (`rArea3` out); 0 for a flipped one.
double prm_stretch_sq_area(const detail::PatchVec& rQ1, const detail::PatchVec& rQ2,
                           const detail::PatchVec& rQ3, double S1, double S2, double S3, double T1,
                           double T2, double T3, double& rArea3) {
    const double e1[3] = {rQ2[0] - rQ1[0], rQ2[1] - rQ1[1], rQ2[2] - rQ1[2]};
    const double e2[3] = {rQ3[0] - rQ1[0], rQ3[1] - rQ1[1], rQ3[2] - rQ1[2]};
    const double n[3] = {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                         e1[0] * e2[1] - e1[1] * e2[0]};
    rArea3 = 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    const double a2 = 0.5 * ((S2 - S1) * (T3 - T1) - (S3 - S1) * (T2 - T1));
    if (!(a2 > 0.0) || !(rArea3 > 0.0))
        return 0.0;
    double ss[3], st[3];
    for (std::size_t k = 0; k < 3; ++k) {
        ss[k] = (rQ1[k] * (T2 - T3) + rQ2[k] * (T3 - T1) + rQ3[k] * (T1 - T2)) / (2.0 * a2);
        st[k] = (rQ1[k] * (S3 - S2) + rQ2[k] * (S1 - S3) + rQ3[k] * (S2 - S1)) / (2.0 * a2);
    }
    const double a = ss[0] * ss[0] + ss[1] * ss[1] + ss[2] * ss[2];
    const double c = st[0] * st[0] + st[1] * st[1] + st[2] * st[2];
    return rArea3 * 0.5 * (a + c);
}

// A chart's L2 stretch, its coordinates indexed like the piece's local points.
double prm_chart_stretch(const std::vector<double>& rU, const std::vector<double>& rV,
                         const std::vector<detail::PatchVec>& rPts,
                         const std::vector<detail::PatchTri>& rLocal) {
    double num = 0.0, den = 0.0;
    for (const detail::PatchTri& t : rLocal) {
        double a3 = 0.0;
        num += prm_stretch_sq_area(rPts[t[0]], rPts[t[1]], rPts[t[2]], rU[t[0]], rU[t[1]], rU[t[2]],
                                   rV[t[0]], rV[t[1]], rV[t[2]], a3);
        den += a3;
    }
    return den > 0.0 ? std::sqrt(num / den) : 0.0;
}

// The map of a piece made a chart: global points, the coordinates, and the
// flipped triangles counted against the majority orientation.
PrmChart prm_make_chart(const std::vector<std::size_t>& rPiece,
                        const std::vector<std::int64_t>& rGlobal,
                        const std::vector<detail::PatchVec>& rPts,
                        const std::vector<detail::PatchTri>& rLocal,
                        const detail::DiscParam& rMap) {
    PrmChart chart;
    chart.mTris = rPiece;
    chart.mPoints = rGlobal;
    chart.mU = rMap.mU;
    chart.mV = rMap.mV;
    // Wind like the surface: the majority orientation counter-clockwise.
    double signed_sum = 0.0;
    for (const detail::PatchTri& t : rLocal)
        signed_sum += prm_signed_area2(chart.mU, chart.mV, t);
    if (signed_sum < 0.0)
        for (double& u : chart.mU)
            u = -u;
    // One texel density everywhere: the area in the plane is the area on the surface.
    double area2 = 0.0, area3 = 0.0;
    for (const detail::PatchTri& t : rLocal) {
        const double a = prm_signed_area2(chart.mU, chart.mV, t);
        if (a <= 0.0)
            ++chart.mFlips;
        area2 += std::fabs(a);
        area3 += prm_tri_area3(rPts, t);
    }
    if (area2 > 0.0 && area3 > 0.0) {
        const double s = std::sqrt(area3 / area2);
        for (std::size_t i = 0; i < chart.mU.size(); ++i) {
            chart.mU[i] *= s;
            chart.mV[i] *= s;
        }
    }
    return chart;
}

// Where input point @p P sits in a chart's point list: sorted, except a
// single-triangle chart's three.
std::size_t prm_index(const PrmChart& rChart, std::int64_t P) {
    const auto& pts = rChart.mPoints;
    if (pts.size() <= 3)
        return static_cast<std::size_t>(std::find(pts.begin(), pts.end(), P) - pts.begin());
    return static_cast<std::size_t>(std::lower_bound(pts.begin(), pts.end(), P) - pts.begin());
}

// Turn a chart to its principal axes (the longer one along u) and move its
// bounding box to the origin.
void prm_align(PrmChart& rChart) {
    const std::size_t n = rChart.mU.size();
    if (n == 0)
        return;
    double mu = 0.0, mv = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        mu += rChart.mU[i];
        mv += rChart.mV[i];
    }
    mu /= static_cast<double>(n);
    mv /= static_cast<double>(n);
    double cuu = 0.0, cvv = 0.0, cuv = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double du = rChart.mU[i] - mu, dv = rChart.mV[i] - mv;
        cuu += du * du;
        cvv += dv * dv;
        cuv += du * dv;
    }
    const double theta = 0.5 * std::atan2(2.0 * cuv, cuu - cvv);
    const double c = std::cos(theta), s = std::sin(theta);
    double lo_u = std::numeric_limits<double>::infinity(), lo_v = lo_u;
    double hi_u = -lo_u, hi_v = -lo_u;
    for (std::size_t i = 0; i < n; ++i) {
        const double du = rChart.mU[i] - mu, dv = rChart.mV[i] - mv;
        rChart.mU[i] = c * du + s * dv;
        rChart.mV[i] = -s * du + c * dv;
        lo_u = std::min(lo_u, rChart.mU[i]);
        hi_u = std::max(hi_u, rChart.mU[i]);
        lo_v = std::min(lo_v, rChart.mV[i]);
        hi_v = std::max(hi_v, rChart.mV[i]);
    }
    for (std::size_t i = 0; i < n; ++i) {
        rChart.mU[i] -= lo_u;
        rChart.mV[i] -= lo_v;
    }
    rChart.mW = hi_u - lo_u;
    rChart.mH = hi_v - lo_v;
}

}  // namespace

ParameterizeResult parameterize(const Mesh& rMesh, const ParameterizeOptions& rOptions) {
    const std::string op = kPrmOp;
    if (!(rOptions.mPadding >= 0.0) || !std::isfinite(rOptions.mPadding))
        throw std::invalid_argument(op + ": padding must be a finite number >= 0");
    if (rOptions.mMaxDepth < 0)
        throw std::invalid_argument(op + ": max_depth must be >= 0");
    if (!(rOptions.mMaxStretch >= 1.0))
        throw std::invalid_argument(op + ": max_stretch must be >= 1");

    // --- the surface: every triangle block; lines and vertices ride along ---
    const std::size_t np = rMesh.NumPoints();
    std::vector<PrmTri> tris;
    std::vector<std::size_t> tri_block;  // per triangle, its block
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string& type = cb.Type();
        if (cb.IsPolyhedron())
            throw std::invalid_argument(
                op + ": a volume mesh; take its surface with extract_surface first");
        const int d = cb.IsRagged() ? 2 : cell_type_dimension(cell_type_from_name(type));
        if (type == "triangle" && !cb.IsRagged()) {
            const NDArray& conn = cb.Conn();
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                PrmTri t;
                for (std::size_t k = 0; k < 3; ++k) {
                    t[k] = detail::read_int(conn, c * 3 + k);
                    if (t[k] < 0 || static_cast<std::size_t>(t[k]) >= np)
                        throw std::invalid_argument(op +
                                                    ": a triangle names a point outside "
                                                    "the mesh");
                }
                tris.push_back(t);
                tri_block.push_back(b);
            }
        } else if (d == 3) {
            throw std::invalid_argument(op + ": '" + type +
                                        "' cells make a volume; take its surface with "
                                        "extract_surface first");
        } else if (d == 2 || type.rfind("polygon", 0) == 0) {
            throw std::invalid_argument(
                op + ": '" + type +
                "' cells are not linear triangles; triangulate them with "
                "convert_cells(mode=\"simplexify\") (or \"linearize\" a higher-order one) first");
        }
    }
    if (tris.empty())
        throw std::invalid_argument(op + ": the mesh has no triangle");

    std::vector<detail::PatchVec> xyz(np, {0.0, 0.0, 0.0});
    {
        const NDArray& points = rMesh.Points();
        const std::size_t dim = rMesh.PointDim();
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t k = 0; k < dim && k < 3; ++k)
                xyz[i][k] = detail::read_double(points, i * dim + k);
    }

    // --- the edge-connected pieces, in order of their lowest triangle ---
    std::vector<std::size_t> parent(tris.size());
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    auto find = [&](std::size_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    {
        std::vector<std::array<std::int64_t, 3>> edges;  // lo, hi, triangle
        edges.reserve(3 * tris.size());
        for (std::size_t t = 0; t < tris.size(); ++t)
            for (std::size_t k = 0; k < 3; ++k) {
                const std::int64_t a = tris[t][k], b = tris[t][(k + 1) % 3];
                edges.push_back({std::min(a, b), std::max(a, b), static_cast<std::int64_t>(t)});
            }
        std::sort(edges.begin(), edges.end());
        for (std::size_t i = 1; i < edges.size(); ++i)
            if (edges[i][0] == edges[i - 1][0] && edges[i][1] == edges[i - 1][1]) {
                const std::size_t ra = find(static_cast<std::size_t>(edges[i][2]));
                const std::size_t rb = find(static_cast<std::size_t>(edges[i - 1][2]));
                if (ra != rb)
                    parent[std::max(ra, rb)] = std::min(ra, rb);
            }
    }
    std::vector<PrmPiece> queue;
    {
        std::vector<std::int64_t> slot(tris.size(), -1);
        for (std::size_t t = 0; t < tris.size(); ++t) {
            const std::size_t r = find(t);
            if (slot[r] < 0) {
                slot[r] = static_cast<std::int64_t>(queue.size());
                queue.push_back({});
            }
            queue[static_cast<std::size_t>(slot[r])].mTris.push_back(t);
        }
    }

    // --- cut and flatten, depth first in a fixed order ---
    std::vector<PrmChart> charts;
    std::vector<std::int64_t> global;
    std::vector<detail::PatchVec> pts;
    std::vector<detail::PatchTri> local;
    std::reverse(queue.begin(), queue.end());  // a stack: the first piece on top
    while (!queue.empty()) {
        PrmPiece piece = std::move(queue.back());
        queue.pop_back();
        prm_local(xyz, tris, piece.mTris, global, pts, local);
        const bool disc = detail::param_is_disc(local);
        detail::DiscParam map;
        const bool at_cap = piece.mDepth >= rOptions.mMaxDepth;
        bool have_chart = false;
        PrmChart best;
        // A very large piece is halved before any solve: its conformal map is
        // slow, and rarely stretches little enough to keep.
        const bool too_big = piece.mTris.size() > kPrmMaxFlattenTriangles && !at_cap;
        if (disc && !too_big) {
            map = detail::param_disc(pts, local);
            if (map.mOk) {
                // A map that folds nothing but stretches too much is halved
                // again while the depth allows.
                best = prm_make_chart(piece.mTris, global, pts, local, map);
                have_chart = true;
                if (at_cap || piece.mTris.size() < 2 ||
                    prm_chart_stretch(best.mU, best.mV, pts, local) <= rOptions.mMaxStretch) {
                    charts.push_back(std::move(best));
                    continue;
                }
            }
        }
        // Split: a genus-0 piece that is not a disc along its cut graph,
        // anything else in halves.
        std::vector<std::int64_t> labels;
        std::size_t parts = 0;
        if (!disc)
            labels = detail::patch_split_to_discs(pts, local, 1, parts);
        if (parts < 2 && !at_cap)
            labels = detail::patch_bisect(pts, local, 1, parts);
        if (parts >= 2 && labels.size() == piece.mTris.size()) {
            std::vector<PrmPiece> split(parts);
            for (std::size_t i = 0; i < labels.size(); ++i)
                if (labels[i] >= 0)
                    split[static_cast<std::size_t>(labels[i])].mTris.push_back(piece.mTris[i]);
            // A triangle the split left out (-1) becomes its own chart.
            std::vector<std::size_t> loose;
            for (std::size_t i = 0; i < labels.size(); ++i)
                if (labels[i] < 0)
                    loose.push_back(piece.mTris[i]);
            prm_single_triangles(xyz, tris, loose, charts);
            for (std::size_t k = parts; k-- > 0;)
                if (!split[k].mTris.empty()) {
                    split[k].mDepth = piece.mDepth + 1;
                    queue.push_back(std::move(split[k]));
                }
            continue;
        }
        // Nothing splits it further: keep the best map, or one chart per triangle.
        if (have_chart) {
            charts.push_back(std::move(best));
            continue;
        }
        bool usable = disc && map.mU.size() == pts.size() && map.mV.size() == pts.size();
        for (std::size_t i = 0; usable && i < map.mU.size(); ++i)
            usable = std::isfinite(map.mU[i]) && std::isfinite(map.mV[i]);
        if (usable)
            charts.push_back(prm_make_chart(piece.mTris, global, pts, local, map));
        else
            prm_single_triangles(xyz, tris, piece.mTris, charts);
    }

    ParameterizeResult result;
    result.mNumCharts = static_cast<std::int64_t>(charts.size());

    // --- stretch, at each chart's own (area-preserving) scale ---
    {
        double num = 0.0, den = 0.0;
        for (const PrmChart& chart : charts) {
            for (std::size_t t : chart.mTris) {
                std::array<std::size_t, 3> loc{};
                for (std::size_t k = 0; k < 3; ++k)
                    loc[k] = prm_index(chart, tris[t][k]);
                double a3 = 0.0;
                num += prm_stretch_sq_area(
                    xyz[static_cast<std::size_t>(tris[t][0])],
                    xyz[static_cast<std::size_t>(tris[t][1])],
                    xyz[static_cast<std::size_t>(tris[t][2])], chart.mU[loc[0]], chart.mU[loc[1]],
                    chart.mU[loc[2]], chart.mV[loc[0]], chart.mV[loc[1]], chart.mV[loc[2]], a3);
                den += a3;
            }
            result.mNumFlipped += static_cast<std::int64_t>(chart.mFlips);
        }
        result.mStretch = den > 0.0 ? std::sqrt(num / den) : 0.0;
    }

    // --- pack: shelves, tallest first ---
    double total_area3 = 0.0;
    for (const PrmTri& t : tris)
        total_area3 +=
            prm_tri_area3(xyz, {static_cast<std::size_t>(t[0]), static_cast<std::size_t>(t[1]),
                                static_cast<std::size_t>(t[2])});
    const double gap = rOptions.mPadding * std::sqrt(total_area3);
    for (PrmChart& chart : charts)
        prm_align(chart);
    std::vector<std::size_t> order(charts.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return charts[a].mH > charts[b].mH; });
    double box_area = 0.0, widest = 0.0;
    for (const PrmChart& chart : charts) {
        box_area += (chart.mW + gap) * (chart.mH + gap);
        widest = std::max(widest, chart.mW + gap);
    }
    // Shelves of a given width; the side of the square the result needs. The
    // width is chosen among a fixed ladder for the smallest square, so the
    // atlas fills the square instead of growing tall.
    auto shelve = [&](double RowWidth, bool Place) {
        double x = 0.0, y = 0.0, shelf = 0.0, used_w = 0.0;
        for (std::size_t i : order) {
            PrmChart& chart = charts[i];
            if (x > 0.0 && x + chart.mW > RowWidth) {
                y += shelf + gap;
                x = 0.0;
                shelf = 0.0;
            }
            if (Place) {
                chart.mX = x;
                chart.mY = y;
            }
            used_w = std::max(used_w, x + chart.mW);
            x += chart.mW + gap;
            shelf = std::max(shelf, chart.mH);
        }
        return std::max(used_w, y + shelf);
    };
    const double base_width = std::max(widest, std::sqrt(box_area));
    double best_width = base_width;
    double side = shelve(base_width, false);
    for (int k = -20; k <= 40; ++k) {
        const double width = std::max(widest, base_width * (1.0 + 0.025 * k));
        const double s = shelve(width, false);
        if (s < side) {
            side = s;
            best_width = width;
        }
    }
    side = shelve(best_width, true);
    const double scale = side > 0.0 ? 1.0 / side : 1.0;
    {
        double covered = 0.0;
        for (const PrmChart& chart : charts)
            for (std::size_t t : chart.mTris) {
                std::array<std::size_t, 3> loc{};
                for (std::size_t k = 0; k < 3; ++k)
                    loc[k] = prm_index(chart, tris[t][k]);
                covered += std::fabs(prm_signed_area2(chart.mU, chart.mV, loc));
            }
        result.mPackingEfficiency = covered * scale * scale;
    }

    // --- split the points at the seams ---
    // Per input point, the (chart, output point) pairs, in the order met.
    std::vector<std::vector<std::pair<std::size_t, std::int64_t>>> copies(np);
    std::vector<std::size_t> chart_of_tri(tris.size(), 0);
    for (std::size_t c = 0; c < charts.size(); ++c)
        for (std::size_t t : charts[c].mTris)
            chart_of_tri[t] = c;
    std::vector<std::int64_t> parent_point(np);
    std::iota(parent_point.begin(), parent_point.end(), std::int64_t{0});
    std::vector<double> uv(2 * np, std::nan(""));
    std::vector<PrmTri> out_tris(tris.size());
    auto chart_uv = [&](std::size_t C, std::int64_t P, double& rU, double& rV) {
        const PrmChart& chart = charts[C];
        const std::size_t i = prm_index(chart, P);
        rU = (chart.mX + chart.mU[i]) * scale;
        rV = (chart.mY + chart.mV[i]) * scale;
    };
    for (std::size_t t = 0; t < tris.size(); ++t) {
        const std::size_t c = chart_of_tri[t];
        for (std::size_t k = 0; k < 3; ++k) {
            const std::int64_t p = tris[t][k];
            auto& list = copies[static_cast<std::size_t>(p)];
            std::int64_t id = -1;
            for (const auto& [chart, out] : list)
                if (chart == c)
                    id = out;
            if (id < 0) {
                id = list.empty() ? p : static_cast<std::int64_t>(parent_point.size());
                if (!list.empty()) {
                    parent_point.push_back(p);
                    uv.push_back(0.0);
                    uv.push_back(0.0);
                }
                list.emplace_back(c, id);
                chart_uv(c, p, uv[2 * static_cast<std::size_t>(id)],
                         uv[2 * static_cast<std::size_t>(id) + 1]);
            }
            out_tris[t][k] = id;
        }
    }
    const std::size_t n_out = parent_point.size();
    result.mNumSeamPoints = static_cast<std::int64_t>(n_out - np);

    // --- the output mesh ---
    Mesh& out = result.mMesh;
    out.AssignPoints(detail::subset_gather_rows(rMesh.Points(), parent_point));
    std::size_t tri_at = 0;
    std::vector<NDArray> chart_blocks;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t n = cb.NumCells();
        NDArray chart_ids = NDArray::Uninit(DType::Int64, {n});
        std::int64_t* p_chart = chart_ids.As<std::int64_t>();
        if (cb.Type() == "triangle" && !cb.IsRagged()) {
            NDArray conn = NDArray::Uninit(DType::Int64, {n, 3});
            std::int64_t* p_conn = conn.As<std::int64_t>();
            for (std::size_t c = 0; c < n; ++c, ++tri_at) {
                for (std::size_t k = 0; k < 3; ++k)
                    p_conn[3 * c + k] = out_tris[tri_at][k];
                p_chart[c] = static_cast<std::int64_t>(chart_of_tri[tri_at]);
            }
            out.AddCellBlock(cb.Type(), std::move(conn));
        } else {
            out.AddCellBlock(cb.Type(), detail::data_owned_copy(cb.Conn()));
            std::fill(p_chart, p_chart + n, std::int64_t{-1});
        }
        chart_blocks.push_back(std::move(chart_ids));
    }
    for (const std::string& name : rMesh.PointDataNames())
        if (name != kUvName && name != kParameterizeParentPointName)
            out.AddPointData(name, detail::subset_gather_rows(rMesh.PointData(name), parent_point));
    for (const std::string& name : rMesh.CellDataNames()) {
        if (name == kParameterizeChartName)
            continue;
        std::vector<NDArray> blocks;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            blocks.push_back(detail::data_owned_copy(rMesh.CellData(name, b)));
        out.AddCellData(name, std::move(blocks));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));
    {
        NDArray uv_array = NDArray::Uninit(DType::Float64, {n_out, 2});
        std::copy(uv.begin(), uv.end(), uv_array.As<double>());
        out.AddPointData(kUvName, std::move(uv_array));
        NDArray parents = NDArray::Uninit(DType::Int64, {n_out});
        std::copy(parent_point.begin(), parent_point.end(), parents.As<std::int64_t>());
        out.AddPointData(kParameterizeParentPointName, std::move(parents));
        out.AddCellData(kParameterizeChartName, std::move(chart_blocks));
    }
    // Every input point keeps its index and every cell its place, so regions
    // carry as they are.
    detail::RegionRemap maps;
    maps.mOpName = "parameterize";
    detail::remap_regions(rMesh, out, maps);
    return result;
}

}  // namespace meshioplusplus
