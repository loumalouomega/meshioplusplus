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

// The kernel-free STEP B-rep writer. See detail/step_brep.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/ndarray.hpp"

// Project includes (private, not installed)
#include "primitive_geom.hpp"
#include "step_brep.hpp"
#include "step_stream.hpp"

namespace meshioplusplus {
namespace detail {
namespace {

constexpr const char* kSbPrefix = "meshio++: step: ";
constexpr double kSbPi = 3.14159265358979323846;
constexpr int kSbParams = kPrimitiveParams;

double sb_dot(const double* pA, const double* pB) {
    return pA[0] * pB[0] + pA[1] * pB[1] + pA[2] * pB[2];
}

double sb_dist(const double* pA, const double* pB) {
    const double d[3] = {pA[0] - pB[0], pA[1] - pB[1], pA[2] - pB[2]};
    return std::sqrt(sb_dot(d, d));
}

// A unit direction perpendicular to unit pA: the axis least aligned with it,
// made orthogonal.
void sb_perpendicular(const double* pA, double* pOut) {
    int k = 0;
    for (int j = 1; j < 3; ++j)
        if (std::fabs(pA[j]) < std::fabs(pA[k]))
            k = j;
    double e[3] = {0.0, 0.0, 0.0};
    e[k] = 1.0;
    const double h = sb_dot(e, pA);
    for (int j = 0; j < 3; ++j)
        pOut[j] = e[j] - h * pA[j];
    const double len = std::sqrt(sb_dot(pOut, pOut));
    for (int j = 0; j < 3; ++j)
        pOut[j] /= len;
}

// Distance from pP to the segment pA-pB.
double sb_segment_distance(const double* pP, const double* pA, const double* pB) {
    const double ab[3] = {pB[0] - pA[0], pB[1] - pA[1], pB[2] - pA[2]};
    const double ap[3] = {pP[0] - pA[0], pP[1] - pA[1], pP[2] - pA[2]};
    const double len2 = sb_dot(ab, ab);
    double t = len2 > 0.0 ? sb_dot(ap, ab) / len2 : 0.0;
    t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
    const double q[3] = {pA[0] + t * ab[0], pA[1] + t * ab[1], pA[2] + t * ab[2]};
    return sb_dist(pP, q);
}

std::uint64_t sb_key(std::int64_t U, std::int64_t V) {
    const std::uint64_t lo = static_cast<std::uint64_t>(std::min(U, V));
    const std::uint64_t hi = static_cast<std::uint64_t>(std::max(U, V));
    return (lo << 32) ^ hi ^ (hi >> 32 << 48);
}

// One row of the fit table.
struct SbRow {
    std::int64_t mLabel = -1;
    std::int64_t mType = 0;
    double mParams[kSbParams] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double mTolerance = 0.0;
    bool mAccepted = false;  // a finite primitive of type 1..5
    bool mAnalytic = false;  // written on its surface (this round)
    bool mSense = true;      // its surface's normal agrees with the facets'
    // A closed sphere or torus patch with no boundary has no STEP face: it is
    // split in two by a plane through its centre, the far half going to its
    // twin row.
    std::int64_t mTwin = -1;  // this row's twin (on an input row)
    bool mSplit = false;
    double mSplitNormal[3] = {1.0, 0.0, 0.0};
};

// A chain of boundary mesh edges between two faces.
struct SbChain {
    std::vector<std::int64_t> mVerts;
    std::int64_t mFwd = -1;  // the group whose facets walk mVerts forwards
    std::int64_t mBwd = -1;
};

// An EDGE_CURVE: verts [mFrom, mTo] of a chain, its kind and geometry.
struct SbSegment {
    std::size_t mChain = 0;
    std::size_t mFrom = 0;
    std::size_t mTo = 0;
    int mKind = 0;  // 0 line, 1 circle, 2 polyline B-spline
    double mCentre[3] = {0, 0, 0};
    double mAxis[3] = {0, 0, 1};
    double mRadius = 0.0;
    std::int64_t mId = 0;
};

// Douglas-Peucker breakpoints of verts [I0, I1] (both kept).
void sb_split(const std::vector<std::int64_t>& rV, const std::vector<double>& rPos, std::size_t I0,
              std::size_t I1, double Tol, std::vector<std::size_t>& rOut) {
    double worst = -1.0;
    std::size_t at = I0;
    for (std::size_t i = I0 + 1; i < I1; ++i) {
        const double d = sb_segment_distance(&rPos[static_cast<std::size_t>(rV[i]) * 3],
                                             &rPos[static_cast<std::size_t>(rV[I0]) * 3],
                                             &rPos[static_cast<std::size_t>(rV[I1]) * 3]);
        if (d > worst) {
            worst = d;
            at = i;
        }
    }
    if (worst > Tol) {
        sb_split(rV, rPos, I0, at, Tol, rOut);
        rOut.push_back(at);
        sb_split(rV, rPos, at, I1, Tol, rOut);
    }
}

}  // namespace

StepBrepResult step_brep_write(const std::string& rPath, const Mesh& rSurface,
                               const ClosedFacets& rShell, const StepBrepOptions& rOptions) {
    StepBrepResult res;
    const std::size_t nf = rShell.NumFacets();
    const std::size_t npts = rShell.mXyz.size() / 3;
    const std::vector<double>& xyz = rShell.mXyz;
    std::vector<std::int64_t> start = rShell.mStart;
    std::vector<std::int64_t> nodes = rShell.mNodes;
    std::vector<double> normals = rShell.mNormals;

    // Pieces, and the ones wound inwards turned outwards.
    std::vector<std::size_t> parent(nf);
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    auto find = [&](std::size_t X) {
        while (parent[X] != X) {
            parent[X] = parent[parent[X]];
            X = parent[X];
        }
        return X;
    };
    {
        std::unordered_map<std::uint64_t, std::size_t> owner;
        owner.reserve(nodes.size());
        for (std::size_t f = 0; f < nf; ++f) {
            const std::int64_t b = start[f], s = start[f + 1] - b;
            for (std::int64_t k = 0; k < s; ++k) {
                const std::uint64_t key = sb_key(nodes[static_cast<std::size_t>(b + k)],
                                                 nodes[static_cast<std::size_t>(b + (k + 1) % s)]);
                auto ins = owner.emplace(key, f);
                if (!ins.second) {
                    const std::size_t ra = find(f), rb = find(ins.first->second);
                    if (ra != rb)
                        parent[std::max(ra, rb)] = std::min(ra, rb);
                }
            }
        }
    }
    std::vector<std::size_t> piece(nf);
    std::vector<std::size_t> roots;
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t r = find(f);
        if (r == f)
            roots.push_back(f);
        piece[f] = static_cast<std::size_t>(std::lower_bound(roots.begin(), roots.end(), r) -
                                            roots.begin());
    }
    const std::size_t npieces = roots.size();
    std::vector<double> volume(npieces, 0.0);
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t b = start[f], s = start[f + 1] - b;
        const double* p0 = &xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(b)]) * 3];
        for (std::int64_t k = 1; k + 1 < s; ++k) {
            const double* p1 =
                &xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(b + k)]) * 3];
            const double* p2 =
                &xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(b + k + 1)]) * 3];
            volume[piece[f]] +=
                (p0[0] * (p1[1] * p2[2] - p1[2] * p2[1]) - p0[1] * (p1[0] * p2[2] - p1[2] * p2[0]) +
                 p0[2] * (p1[0] * p2[1] - p1[1] * p2[0])) /
                6.0;
        }
    }
    for (std::size_t f = 0; f < nf; ++f)
        if (volume[piece[f]] < 0.0) {
            std::reverse(nodes.begin() + start[f], nodes.begin() + start[f + 1]);
            for (int k = 0; k < 3; ++k)
                normals[f * 3 + static_cast<std::size_t>(k)] =
                    -normals[f * 3 + static_cast<std::size_t>(k)];
        }
    for (double v : volume)
        res.mVolume += std::fabs(v);

    // Per undirected edge, the facets walking it lo -> hi and hi -> lo.
    std::unordered_map<std::uint64_t, std::array<std::int64_t, 2>> walkers;
    walkers.reserve(nodes.size());
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t b = start[f], s = start[f + 1] - b;
        for (std::int64_t k = 0; k < s; ++k) {
            const std::int64_t u = nodes[static_cast<std::size_t>(b + k)];
            const std::int64_t v = nodes[static_cast<std::size_t>(b + (k + 1) % s)];
            auto ins = walkers.emplace(sb_key(u, v), std::array<std::int64_t, 2>{-1, -1});
            ins.first->second[u < v ? 0 : 1] = static_cast<std::int64_t>(f);
        }
    }
    auto walker = [&](std::int64_t U, std::int64_t V) -> std::int64_t {
        auto it = walkers.find(sb_key(U, V));
        if (it == walkers.end())
            return -1;
        return U < V ? it->second[0] : it->second[1];
    };
    auto ring_next = [&](std::size_t F, std::int64_t V) -> std::int64_t {
        const std::int64_t b = start[F], n = start[F + 1] - b;
        for (std::int64_t k = 0; k < n; ++k)
            if (nodes[static_cast<std::size_t>(b + k)] == V)
                return nodes[static_cast<std::size_t>(b + (k + 1) % n)];
        return -1;
    };

    double lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
        lo[k] = std::numeric_limits<double>::infinity();
        hi[k] = -lo[k];
    }
    for (std::int64_t v : nodes)
        for (int k = 0; k < 3; ++k) {
            lo[k] =
                std::min(lo[k], xyz[static_cast<std::size_t>(v) * 3 + static_cast<std::size_t>(k)]);
            hi[k] =
                std::max(hi[k], xyz[static_cast<std::size_t>(v) * 3 + static_cast<std::size_t>(k)]);
        }
    const double ext[3] = {hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
    const double diag = std::sqrt(sb_dot(ext, ext));
    const double eps = 1e-6 * diag;

    // The fit table and each facet's row.
    std::vector<SbRow> rows;
    std::vector<std::int64_t> row_of(nf, -1);
    std::vector<std::int64_t> label_of(nf, -1);
    const bool has_patches = rSurface.HasCellData("segment:patch");
    if (has_patches) {
        const std::vector<std::int64_t> bases = block_bases(rSurface);
        std::vector<std::int64_t> label(static_cast<std::size_t>(total_cells(bases)), -1);
        for (std::size_t b = 0; b < rSurface.NumCellBlocks(); ++b) {
            const NDArray& a = rSurface.CellData("segment:patch", b);
            const std::size_t nc = rSurface.Cells(b).NumCells();
            for (std::size_t c = 0; c < nc && c < a.Size(); ++c)
                label[static_cast<std::size_t>(bases[b]) + c] = read_int(a, c);
        }
        for (std::size_t f = 0; f < nf; ++f)
            label_of[f] = label[static_cast<std::size_t>(rShell.mCell[f])];
    }
    const bool has_table =
        rSurface.HasFieldData("fit:patch_id") && rSurface.HasFieldData("fit:patch_type") &&
        rSurface.HasFieldData("fit:patch_params") && rSurface.HasFieldData("fit:patch_tolerance");
    if (rOptions.mFitted && has_patches && has_table) {
        const NDArray& ids = rSurface.FieldData("fit:patch_id");
        const NDArray& types = rSurface.FieldData("fit:patch_type");
        const NDArray& params = rSurface.FieldData("fit:patch_params");
        const NDArray& tols = rSurface.FieldData("fit:patch_tolerance");
        const std::size_t m = ids.Size();
        if (types.Size() != m || tols.Size() != m ||
            params.Size() != m * static_cast<std::size_t>(kSbParams))
            throw WriteError(std::string(kSbPrefix) +
                             "the fit:patch_* arrays disagree on the number of patches");
        rows.resize(m);
        std::unordered_map<std::int64_t, std::int64_t> by_label;
        for (std::size_t r = 0; r < m; ++r) {
            SbRow& row = rows[r];
            row.mLabel = read_int(ids, r);
            row.mType = read_int(types, r);
            bool finite = true;
            for (int k = 0; k < kSbParams; ++k) {
                row.mParams[k] = read_double(params, r * kSbParams + static_cast<std::size_t>(k));
                finite = finite && std::isfinite(row.mParams[k]);
            }
            row.mTolerance = read_double(tols, r);
            row.mAccepted = row.mType >= 1 && row.mType <= 5 && finite;
            by_label.emplace(row.mLabel, static_cast<std::int64_t>(r));
        }
        for (std::size_t f = 0; f < nf; ++f) {
            auto it = by_label.find(label_of[f]);
            if (it != by_label.end())
                row_of[f] = it->second;
        }
    }
    res.mNumPatches = static_cast<std::int64_t>(rows.size());
    const std::size_t num_input_rows = rows.size();
    for (std::size_t r = 0; r < num_input_rows; ++r) {
        if (!rows[r].mAccepted || (rows[r].mType != 2 && rows[r].mType != 5))
            continue;
        SbRow twin = rows[r];
        twin.mAccepted = false;
        rows[r].mTwin = static_cast<std::int64_t>(rows.size());
        if (rows[r].mType == 5)
            sb_perpendicular(rows[r].mParams + 3, rows[r].mSplitNormal);
        rows.push_back(twin);
    }
    double tol = rOptions.mTolerance;
    if (!(tol > 0.0)) {
        double largest = 0.0;
        for (const SbRow& r : rows)
            if (r.mAccepted && std::isfinite(r.mTolerance))
                largest = std::max(largest, r.mTolerance);
        tol = largest > 0.0 ? 3.0 * largest : 1e-3 * diag;
    }
    res.mTolerance = tol;
    const std::int64_t m = static_cast<std::int64_t>(rows.size());

    // A row is analytic when every point of its facets lies within the
    // tolerance; its sense from the facets' normals.
    {
        std::vector<double> worst(rows.size(), 0.0), agree(rows.size(), 0.0);
        std::vector<std::uint8_t> present(rows.size(), 0);
        for (std::size_t f = 0; f < nf; ++f) {
            const std::int64_t r = row_of[f];
            if (r < 0 || !rows[static_cast<std::size_t>(r)].mAccepted)
                continue;
            const SbRow& row = rows[static_cast<std::size_t>(r)];
            present[static_cast<std::size_t>(r)] = 1;
            double c[3] = {0.0, 0.0, 0.0};
            const std::int64_t b = start[f], s = start[f + 1] - b;
            for (std::int64_t k = 0; k < s; ++k) {
                const double* p =
                    &xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(b + k)]) * 3];
                worst[static_cast<std::size_t>(r)] =
                    std::max(worst[static_cast<std::size_t>(r)],
                             primitive_point_distance(row.mType, row.mParams, p));
                for (int j = 0; j < 3; ++j)
                    c[j] += p[j] / static_cast<double>(s);
            }
            double n[3];
            if (primitive_normal(row.mType, row.mParams, c, n))
                agree[static_cast<std::size_t>(r)] += sb_dot(n, &normals[f * 3]);
        }
        for (std::size_t r = 0; r < rows.size(); ++r) {
            rows[r].mAnalytic = rows[r].mAccepted && present[r] && worst[r] <= tol;
            rows[r].mSense = agree[r] >= 0.0;
            if (rows[r].mAccepted && present[r] && !rows[r].mAnalytic)
                ++res.mNumDemoted;
        }
    }

    // Build rounds: a row whose loops fail is demoted and the faces rebuilt.
    for (;;) {
        // Groups: analytic rows, then coplanar groups `m + seed`.
        std::vector<std::int64_t> group(nf, -1);
        for (std::size_t f = 0; f < nf; ++f) {
            if (row_of[f] < 0 || !rows[static_cast<std::size_t>(row_of[f])].mAnalytic)
                continue;
            const SbRow& row = rows[static_cast<std::size_t>(row_of[f])];
            group[f] = row_of[f];
            if (row.mSplit) {
                double side = 0.0;
                for (std::int64_t k = start[f]; k < start[f + 1]; ++k)
                    for (int j = 0; j < 3; ++j)
                        side +=
                            (xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(k)]) * 3 +
                                 static_cast<std::size_t>(j)] -
                             row.mParams[j]) *
                            row.mSplitNormal[j];
                if (side < 0.0)
                    group[f] = row.mTwin;
            }
        }
        {
            std::vector<double> area(nf, 0.0);
            for (std::size_t f = 0; f < nf; ++f) {
                const std::int64_t b = start[f], s = start[f + 1] - b;
                double twice[3] = {0.0, 0.0, 0.0};
                for (std::int64_t k = 0; k < s; ++k) {
                    const double* a =
                        &xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(b + k)]) * 3];
                    const double* e = &xyz[static_cast<std::size_t>(
                                               nodes[static_cast<std::size_t>(b + (k + 1) % s)]) *
                                           3];
                    twice[0] += a[1] * e[2] - a[2] * e[1];
                    twice[1] += a[2] * e[0] - a[0] * e[2];
                    twice[2] += a[0] * e[1] - a[1] * e[0];
                }
                area[f] = 0.5 * std::sqrt(sb_dot(twice, twice));
            }
            std::vector<std::size_t> order(nf);
            std::iota(order.begin(), order.end(), std::size_t{0});
            std::stable_sort(order.begin(), order.end(),
                             [&](std::size_t A, std::size_t B) { return area[A] > area[B]; });
            std::vector<std::size_t> queue;
            for (std::size_t seed : order) {
                if (group[seed] != -1)
                    continue;
                const std::int64_t id = m + static_cast<std::int64_t>(seed);
                const double* p0 =
                    &xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(start[seed])]) *
                         3];
                const double* n0 = &normals[seed * 3];
                group[seed] = id;
                queue.assign(1, seed);
                for (std::size_t q = 0; q < queue.size(); ++q) {
                    const std::size_t f = queue[q];
                    const std::int64_t b = start[f], s = start[f + 1] - b;
                    for (std::int64_t k = 0; k < s; ++k) {
                        const std::int64_t u = nodes[static_cast<std::size_t>(b + k)];
                        const std::int64_t v = nodes[static_cast<std::size_t>(b + (k + 1) % s)];
                        const std::int64_t other = walker(v, u);
                        if (other < 0)
                            continue;
                        const std::size_t o = static_cast<std::size_t>(other);
                        if (group[o] != -1 || sb_dot(&normals[o * 3], n0) <= 0.0)
                            continue;
                        bool flat = true;
                        for (std::int64_t j = start[o]; j < start[o + 1] && flat; ++j) {
                            const double* p =
                                &xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(j)]) *
                                     3];
                            const double d[3] = {p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]};
                            flat = std::fabs(sb_dot(d, n0)) <= eps;
                        }
                        if (!flat)
                            continue;
                        group[o] = id;
                        queue.push_back(o);
                    }
                }
            }
        }
        auto analytic_group = [&](std::int64_t G) { return G >= 0 && G < m; };
        // A coplanar group every point of which lies on an analytic
        // neighbour's surface joins that face (the lowest such row).
        {
            std::map<std::int64_t, std::vector<std::size_t>> members;
            for (std::size_t f = 0; f < nf; ++f)
                if (!analytic_group(group[f]))
                    members[group[f]].push_back(f);
            for (const auto& kv : members) {
                std::int64_t best = -1;
                for (std::size_t f : kv.second) {
                    const std::int64_t b = start[f], s = start[f + 1] - b;
                    for (std::int64_t k = 0; k < s; ++k) {
                        const std::int64_t u = nodes[static_cast<std::size_t>(b + k)];
                        const std::int64_t v = nodes[static_cast<std::size_t>(b + (k + 1) % s)];
                        const std::int64_t o = walker(v, u);
                        if (o < 0 || !analytic_group(group[static_cast<std::size_t>(o)]))
                            continue;
                        const std::int64_t g = group[static_cast<std::size_t>(o)];
                        if (best >= 0 && g >= best)
                            continue;
                        const SbRow& row = rows[static_cast<std::size_t>(g)];
                        bool on = true;
                        for (std::size_t h : kv.second)
                            for (std::int64_t j = start[h]; j < start[h + 1] && on; ++j)
                                on = primitive_point_distance(
                                         row.mType, row.mParams,
                                         &xyz[static_cast<std::size_t>(
                                                  nodes[static_cast<std::size_t>(j)]) *
                                              3]) <= tol;
                        if (on)
                            best = g;
                    }
                }
                if (best >= 0)
                    for (std::size_t f : kv.second)
                        group[f] = best;
            }
        }

        // The groups around each point, and the point positions.
        std::vector<std::vector<std::int64_t>> groups_at(npts);
        for (std::size_t f = 0; f < nf; ++f)
            for (std::int64_t k = start[f]; k < start[f + 1]; ++k) {
                std::vector<std::int64_t>& g =
                    groups_at[static_cast<std::size_t>(nodes[static_cast<std::size_t>(k)])];
                if (std::find(g.begin(), g.end(), group[f]) == g.end())
                    g.push_back(group[f]);
            }
        for (auto& g : groups_at)
            std::sort(g.begin(), g.end());
        std::vector<double> pos = xyz;
        for (std::size_t v = 0; v < npts; ++v) {
            const std::vector<std::int64_t>& g = groups_at[v];
            if (g.empty() || !std::all_of(g.begin(), g.end(), analytic_group))
                continue;
            double q[3] = {xyz[v * 3], xyz[v * 3 + 1], xyz[v * 3 + 2]};
            for (int it = 0; it < (g.size() > 1 ? 8 : 1); ++it)
                for (std::int64_t r : g) {
                    double out[3];
                    if (primitive_project(rows[static_cast<std::size_t>(r)].mType,
                                          rows[static_cast<std::size_t>(r)].mParams, q, out))
                        for (int j = 0; j < 3; ++j)
                            q[j] = out[j];
                }
            if (sb_dist(q, &xyz[v * 3]) <= tol)
                for (int j = 0; j < 3; ++j)
                    pos[v * 3 + static_cast<std::size_t>(j)] = q[j];
        }

        // Boundary edges and their chains.
        std::vector<std::pair<std::int64_t, std::int64_t>> bedges;
        for (const auto& kv : walkers) {
            const std::int64_t fa = kv.second[0], fb = kv.second[1];
            if (fa < 0 || fb < 0)
                throw WriteError(std::string(kSbPrefix) + "the surface is not closed");
            if (group[static_cast<std::size_t>(fa)] == group[static_cast<std::size_t>(fb)])
                continue;
            std::int64_t u = -1, v = -1;
            // Recover the endpoints from the facet that walks lo -> hi.
            const std::int64_t b = start[static_cast<std::size_t>(fa)],
                               s = start[static_cast<std::size_t>(fa) + 1] - b;
            for (std::int64_t k = 0; k < s; ++k) {
                const std::int64_t x = nodes[static_cast<std::size_t>(b + k)];
                const std::int64_t y = nodes[static_cast<std::size_t>(b + (k + 1) % s)];
                if (x < y && sb_key(x, y) == kv.first) {
                    u = x;
                    v = y;
                    break;
                }
            }
            bedges.emplace_back(u, v);
        }
        std::sort(bedges.begin(), bedges.end());
        std::vector<std::vector<std::int64_t>> nbr(npts);
        for (const auto& e : bedges) {
            nbr[static_cast<std::size_t>(e.first)].push_back(e.second);
            nbr[static_cast<std::size_t>(e.second)].push_back(e.first);
        }
        for (auto& n : nbr)
            std::sort(n.begin(), n.end());
        auto corner = [&](std::int64_t V) {
            const std::size_t v = static_cast<std::size_t>(V);
            return nbr[v].size() != 2 || groups_at[v].size() >= 3;
        };
        std::unordered_map<std::uint64_t, std::uint8_t> taken;
        taken.reserve(bedges.size());
        for (const auto& e : bedges)
            taken.emplace(sb_key(e.first, e.second), 0);
        auto take = [&](std::int64_t U, std::int64_t V) {
            auto it = taken.find(sb_key(U, V));
            if (it == taken.end() || it->second)
                return false;
            it->second = 1;
            return true;
        };
        std::vector<SbChain> chains;
        auto walk = [&](std::int64_t Start, std::int64_t First, bool Loop) {
            SbChain c;
            c.mVerts = {Start, First};
            std::int64_t prev = Start, cur = First;
            while (!(Loop ? cur == Start : corner(cur))) {
                std::int64_t next = -1;
                for (std::int64_t n : nbr[static_cast<std::size_t>(cur)])
                    if (n != prev && take(cur, n)) {
                        next = n;
                        break;
                    }
                if (next < 0) {
                    if (Loop && take(cur, Start))
                        c.mVerts.push_back(Start);
                    break;
                }
                c.mVerts.push_back(next);
                prev = cur;
                cur = next;
            }
            const std::int64_t f = walker(c.mVerts[0], c.mVerts[1]);
            const std::int64_t g = walker(c.mVerts[1], c.mVerts[0]);
            c.mFwd = group[static_cast<std::size_t>(f)];
            c.mBwd = group[static_cast<std::size_t>(g)];
            chains.push_back(std::move(c));
        };
        for (std::int64_t v = 0; v < static_cast<std::int64_t>(npts); ++v) {
            if (nbr[static_cast<std::size_t>(v)].empty() || !corner(v))
                continue;
            for (std::int64_t n : nbr[static_cast<std::size_t>(v)])
                if (take(v, n))
                    walk(v, n, false);
        }
        for (const auto& e : bedges)
            if (take(e.first, e.second))
                walk(e.first, e.second, true);

        // Segments: the EDGE_CURVEs each chain is cut into.
        auto planar_group = [&](std::int64_t G) {
            return !analytic_group(G) || rows[static_cast<std::size_t>(G)].mType == 1;
        };
        std::vector<SbSegment> segs;
        for (std::size_t ci = 0; ci < chains.size(); ++ci) {
            const SbChain& c = chains[ci];
            const std::vector<std::int64_t>& vs = c.mVerts;
            const bool closed = vs.front() == vs.back();
            const std::size_t last = vs.size() - 1;
            auto add = [&](std::size_t A, std::size_t B, int Kind) {
                SbSegment s;
                s.mChain = ci;
                s.mFrom = A;
                s.mTo = B;
                s.mKind = Kind;
                segs.push_back(s);
                return segs.size() - 1;
            };
            auto lines = [&](double LineTol) {
                std::vector<std::size_t> cut{0};
                if (closed) {
                    // Split a closed chain at its point farthest from the start.
                    std::size_t far_index = 1;
                    double best = -1.0;
                    for (std::size_t i = 1; i < last; ++i) {
                        const double d = sb_dist(&pos[static_cast<std::size_t>(vs[i]) * 3],
                                                 &pos[static_cast<std::size_t>(vs[0]) * 3]);
                        if (d > best) {
                            best = d;
                            far_index = i;
                        }
                    }
                    sb_split(vs, pos, 0, far_index, LineTol, cut);
                    cut.push_back(far_index);
                    sb_split(vs, pos, far_index, last, LineTol, cut);
                } else {
                    sb_split(vs, pos, 0, last, LineTol, cut);
                }
                cut.push_back(last);
                for (std::size_t i = 0; i + 1 < cut.size(); ++i)
                    add(cut[i], cut[i + 1], 0);
            };
            const bool flat_side = !analytic_group(c.mFwd) || !analytic_group(c.mBwd);
            if (planar_group(c.mFwd) && planar_group(c.mBwd)) {
                lines(flat_side ? eps : tol);
                continue;
            }
            if (flat_side) {
                for (std::size_t i = 0; i < last; ++i)
                    add(i, i + 1, 0);
                continue;
            }
            // Two analytic faces, one of them curved: a circle about a shared
            // axis, a straight line, or a polyline B-spline.
            bool done = false;
            std::vector<std::array<double, 6>> axes;
            for (std::int64_t g : {c.mFwd, c.mBwd}) {
                const SbRow& r = rows[static_cast<std::size_t>(g)];
                if (r.mType >= 3)
                    axes.push_back({r.mParams[0], r.mParams[1], r.mParams[2], r.mParams[3],
                                    r.mParams[4], r.mParams[5]});
            }
            {
                const SbRow& a = rows[static_cast<std::size_t>(c.mFwd)];
                const SbRow& b = rows[static_cast<std::size_t>(c.mBwd)];
                const SbRow* sphere = a.mType == 2 ? &a : (b.mType == 2 ? &b : nullptr);
                const SbRow* other = sphere == &a ? &b : &a;
                if (sphere != nullptr && other->mType == 1)
                    axes.push_back({sphere->mParams[0], sphere->mParams[1], sphere->mParams[2],
                                    other->mParams[3], other->mParams[4], other->mParams[5]});
                if (a.mType == 2 && b.mType == 2) {
                    double d[3] = {b.mParams[0] - a.mParams[0], b.mParams[1] - a.mParams[1],
                                   b.mParams[2] - a.mParams[2]};
                    const double len = std::sqrt(sb_dot(d, d));
                    if (len > 0.0)
                        axes.push_back({a.mParams[0], a.mParams[1], a.mParams[2], d[0] / len,
                                        d[1] / len, d[2] / len});
                }
            }
            for (const auto& ax : axes) {
                const double* x = ax.data();
                const double* d = ax.data() + 3;
                std::vector<double> h(vs.size()), rho(vs.size());
                double hm = 0.0, rm = 0.0;
                for (std::size_t i = 0; i < vs.size(); ++i) {
                    const double* p = &pos[static_cast<std::size_t>(vs[i]) * 3];
                    const double v[3] = {p[0] - x[0], p[1] - x[1], p[2] - x[2]};
                    h[i] = sb_dot(v, d);
                    const double r[3] = {v[0] - h[i] * d[0], v[1] - h[i] * d[1],
                                         v[2] - h[i] * d[2]};
                    rho[i] = std::sqrt(sb_dot(r, r));
                    hm += h[i];
                    rm += rho[i];
                }
                hm /= static_cast<double>(vs.size());
                rm /= static_cast<double>(vs.size());
                bool fits = rm > tol;
                for (std::size_t i = 0; i < vs.size() && fits; ++i)
                    fits = std::fabs(h[i] - hm) <= tol && std::fabs(rho[i] - rm) <= tol;
                if (!fits)
                    continue;
                const std::size_t k = add(0, last, 1);
                SbSegment& s = segs[k];
                for (int j = 0; j < 3; ++j) {
                    s.mCentre[j] = x[j] + hm * d[j];
                    s.mAxis[j] = d[j];
                }
                s.mRadius = rm;
                done = true;
                break;
            }
            if (done)
                continue;
            if (!closed) {
                bool straight = true;
                for (std::size_t i = 1; i < last && straight; ++i)
                    straight =
                        sb_segment_distance(&pos[static_cast<std::size_t>(vs[i]) * 3],
                                            &pos[static_cast<std::size_t>(vs[0]) * 3],
                                            &pos[static_cast<std::size_t>(vs[last]) * 3]) <= tol;
                if (straight) {
                    add(0, last, 0);
                    continue;
                }
            }
            add(0, last, 2);
        }

        // Each boundary edge's segment, and whether the segment runs lo -> hi.
        std::unordered_map<std::uint64_t, std::pair<std::size_t, bool>> seg_of;
        seg_of.reserve(bedges.size());
        for (std::size_t si = 0; si < segs.size(); ++si) {
            const SbSegment& s = segs[si];
            const std::vector<std::int64_t>& vs = chains[s.mChain].mVerts;
            for (std::size_t i = s.mFrom; i < s.mTo; ++i)
                seg_of[sb_key(vs[i], vs[i + 1])] = {si, vs[i] < vs[i + 1]};
        }

        // Loops: each group's boundary half-edges walked through the fans.
        std::map<std::int64_t, std::vector<std::pair<std::int64_t, std::int64_t>>> half;
        for (const auto& e : bedges) {
            half[group[static_cast<std::size_t>(walker(e.first, e.second))]].push_back(e);
            half[group[static_cast<std::size_t>(walker(e.second, e.first))]].emplace_back(e.second,
                                                                                          e.first);
        }
        auto fan_next = [&](std::int64_t G, std::int64_t X, std::int64_t V) -> std::int64_t {
            std::int64_t f = walker(X, V);
            for (std::size_t guard = 0; f >= 0 && guard <= nf; ++guard) {
                const std::int64_t w = ring_next(static_cast<std::size_t>(f), V);
                if (w < 0)
                    return -1;
                const std::int64_t across = walker(w, V);
                if (across < 0 || group[static_cast<std::size_t>(across)] != G)
                    return w;
                f = across;
            }
            return -1;
        };
        // Per group, its loops as (segment, along) lists.
        std::map<std::int64_t, std::vector<std::vector<std::pair<std::size_t, bool>>>> loops;
        std::vector<std::int64_t> broken;
        for (auto& kv : half) {
            const std::int64_t g = kv.first;
            std::vector<std::pair<std::int64_t, std::int64_t>>& hs = kv.second;
            std::sort(hs.begin(), hs.end());
            std::map<std::pair<std::int64_t, std::int64_t>, std::uint8_t> used;
            for (const auto& h : hs)
                used.emplace(h, 0);
            bool bad = false;
            for (const auto& h0 : hs) {
                if (used[h0])
                    continue;
                std::vector<std::pair<std::int64_t, std::int64_t>> loop;
                std::pair<std::int64_t, std::int64_t> h = h0;
                for (std::size_t guard = 0; guard <= hs.size(); ++guard) {
                    auto it = used.find(h);
                    if (it == used.end() || it->second) {
                        bad = true;
                        break;
                    }
                    it->second = 1;
                    loop.push_back(h);
                    const std::int64_t w = fan_next(g, h.first, h.second);
                    if (w < 0) {
                        bad = true;
                        break;
                    }
                    h = {h.second, w};
                    if (h == h0)
                        break;
                }
                if (bad || h != h0) {
                    bad = true;
                    break;
                }
                // Rotate to a segment start, then emit one entry per segment.
                auto starts = [&](const std::pair<std::int64_t, std::int64_t>& rH) {
                    const auto& so = seg_of.at(sb_key(rH.first, rH.second));
                    const SbSegment& s = segs[so.first];
                    const std::vector<std::int64_t>& vs = chains[s.mChain].mVerts;
                    const bool along = (rH.first < rH.second) == so.second;
                    return rH.first == (along ? vs[s.mFrom] : vs[s.mTo]);
                };
                std::size_t at = loop.size();
                for (std::size_t i = 0; i < loop.size(); ++i)
                    if (starts(loop[i])) {
                        at = i;
                        break;
                    }
                if (at == loop.size()) {
                    bad = true;
                    break;
                }
                std::rotate(loop.begin(), loop.begin() + static_cast<std::ptrdiff_t>(at),
                            loop.end());
                std::vector<std::pair<std::size_t, bool>> entry;
                for (std::size_t i = 0; i < loop.size(); ++i) {
                    if (!starts(loop[i]))
                        continue;
                    const auto& so = seg_of.at(sb_key(loop[i].first, loop[i].second));
                    entry.emplace_back(so.first, (loop[i].first < loop[i].second) == so.second);
                }
                loops[g].push_back(std::move(entry));
            }
            if (bad)
                broken.push_back(g);
        }
        bool retry = false;
        for (std::int64_t g : broken) {
            if (!analytic_group(g))
                throw WriteError(std::string(kSbPrefix) +
                                 "a coplanar face's boundary does not close into loops");
            for (std::size_t r = 0; r < rows.size(); ++r)
                if (static_cast<std::int64_t>(r) == g || rows[r].mTwin == g)
                    rows[r].mAnalytic = false;
            if (rows[static_cast<std::size_t>(g)].mTwin >= 0)
                rows[static_cast<std::size_t>(rows[static_cast<std::size_t>(g)].mTwin)].mAnalytic =
                    false;
            ++res.mNumDemoted;
            retry = true;
        }
        // A face with no boundary at all (a whole sphere or torus) has no STEP
        // face: split it in two, else write it as coplanar groups.
        std::vector<std::int64_t> unbounded;
        for (std::size_t f = 0; f < nf; ++f)
            if (analytic_group(group[f]) && loops.find(group[f]) == loops.end())
                unbounded.push_back(group[f]);
        std::sort(unbounded.begin(), unbounded.end());
        unbounded.erase(std::unique(unbounded.begin(), unbounded.end()), unbounded.end());
        for (std::int64_t g : unbounded) {
            if (!rows[static_cast<std::size_t>(g)].mAnalytic)
                continue;
            SbRow& row = rows[static_cast<std::size_t>(g)];
            if (row.mTwin >= 0 && !row.mSplit) {
                row.mSplit = true;
                rows[static_cast<std::size_t>(row.mTwin)].mAnalytic = true;
                rows[static_cast<std::size_t>(row.mTwin)].mSense = row.mSense;
            } else {
                row.mAnalytic = false;
                if (row.mTwin >= 0)
                    rows[static_cast<std::size_t>(row.mTwin)].mAnalytic = false;
                ++res.mNumDemoted;
            }
            retry = true;
        }
        if (retry)
            continue;

        // Write.
        StepStream st(rPath, "meshio++ B-rep");
        std::map<std::int64_t, std::int64_t> vertex_id;
        auto vertex = [&](std::int64_t V) {
            auto it = vertex_id.find(V);
            if (it != vertex_id.end())
                return it->second;
            const std::int64_t p =
                st.Ent("CARTESIAN_POINT(''," +
                       StepStream::Triple(&pos[static_cast<std::size_t>(V) * 3]) + ")");
            const std::int64_t id = st.Ent("VERTEX_POINT(''," + StepStream::Ref(p) + ")");
            vertex_id.emplace(V, id);
            return id;
        };
        auto direction = [&](const double* pD) {
            return st.Ent("DIRECTION(''," + StepStream::Triple(pD) + ")");
        };
        auto placement = [&](const double* pAt, const double* pZ) {
            double x[3];
            sb_perpendicular(pZ, x);
            const std::int64_t p = st.Ent("CARTESIAN_POINT(''," + StepStream::Triple(pAt) + ")");
            const std::int64_t z = direction(pZ);
            const std::int64_t xd = direction(x);
            return st.Ent("AXIS2_PLACEMENT_3D(''," + StepStream::Ref(p) + "," + StepStream::Ref(z) +
                          "," + StepStream::Ref(xd) + ")");
        };
        for (SbSegment& s : segs) {
            const std::vector<std::int64_t>& vs = chains[s.mChain].mVerts;
            const std::int64_t a = vs[s.mFrom], b = vs[s.mTo];
            const std::int64_t va = vertex(a), vb = vertex(b);
            std::int64_t curve = 0;
            bool sense = true;
            if (s.mKind == 0) {
                const double* pa = &pos[static_cast<std::size_t>(a) * 3];
                const double* pb = &pos[static_cast<std::size_t>(b) * 3];
                double d[3] = {pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]};
                const double len = std::sqrt(sb_dot(d, d));
                if (!(len > 0.0))
                    throw WriteError(std::string(kSbPrefix) +
                                     "an edge has zero length (coincident "
                                     "points; run clean first)");
                for (double& x : d)
                    x /= len;
                const std::int64_t p = st.Ent("CARTESIAN_POINT(''," + StepStream::Triple(pa) + ")");
                const std::int64_t vec = st.Ent("VECTOR(''," + StepStream::Ref(direction(d)) + "," +
                                                step_real(len) + ")");
                curve = st.Ent("LINE(''," + StepStream::Ref(p) + "," + StepStream::Ref(vec) + ")");
                ++res.mNumLineEdges;
            } else if (s.mKind == 1) {
                curve = st.Ent("CIRCLE(''," + StepStream::Ref(placement(s.mCentre, s.mAxis)) + "," +
                               step_real(s.mRadius) + ")");
                double turn = 0.0;
                for (std::size_t i = s.mFrom; i < s.mTo; ++i) {
                    const double* p = &pos[static_cast<std::size_t>(vs[i]) * 3];
                    const double* q = &pos[static_cast<std::size_t>(vs[i + 1]) * 3];
                    const double u[3] = {p[0] - s.mCentre[0], p[1] - s.mCentre[1],
                                         p[2] - s.mCentre[2]};
                    const double w[3] = {q[0] - s.mCentre[0], q[1] - s.mCentre[1],
                                         q[2] - s.mCentre[2]};
                    const double x[3] = {u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
                                         u[0] * w[1] - u[1] * w[0]};
                    turn += sb_dot(x, s.mAxis);
                }
                sense = turn >= 0.0;
                ++res.mNumCircleEdges;
            } else {
                std::string pts, mults, knots;
                const std::size_t n = s.mTo - s.mFrom + 1;
                for (std::size_t i = s.mFrom; i <= s.mTo; ++i) {
                    const std::int64_t p =
                        st.Ent("CARTESIAN_POINT(''," +
                               StepStream::Triple(&pos[static_cast<std::size_t>(vs[i]) * 3]) + ")");
                    pts += (i == s.mFrom ? "" : ",") + StepStream::Ref(p);
                }
                for (std::size_t i = 0; i < n; ++i) {
                    mults += (i == 0 ? "" : ",") + std::string(i == 0 || i + 1 == n ? "2" : "1");
                    knots += (i == 0 ? "" : ",") + step_real(static_cast<double>(i));
                }
                curve = st.Ent("B_SPLINE_CURVE_WITH_KNOTS('',1,(" + pts +
                               "),.POLYLINE_FORM.,.F.,.F.,(" + mults + "),(" + knots +
                               "),.UNSPECIFIED.)");
                ++res.mNumSplineEdges;
            }
            s.mId = st.Ent("EDGE_CURVE(''," + StepStream::Ref(va) + "," + StepStream::Ref(vb) +
                           "," + StepStream::Ref(curve) + (sense ? ",.T.)" : ",.F.)"));
        }

        // Faces, numbered in group order.
        std::map<std::int64_t, std::int64_t> face_index;
        std::vector<std::vector<std::int64_t>> shell_faces(npieces);
        std::map<std::int64_t, std::size_t> seed_of;  // group -> a facet of it
        for (std::size_t f = 0; f < nf; ++f)
            seed_of.emplace(group[f], f);
        for (const auto& kv : loops) {
            const std::int64_t g = kv.first;
            const std::size_t f0 = seed_of.at(g);
            const bool analytic = analytic_group(g);
            const SbRow* row = analytic ? &rows[static_cast<std::size_t>(g)] : nullptr;
            // The surface.
            std::int64_t surface = 0;
            bool same_sense = true;
            double plane_n[3] = {0.0, 0.0, 0.0};
            const bool planar = !analytic || row->mType == 1;
            if (!analytic) {
                const std::size_t seed = static_cast<std::size_t>(g - m);
                for (int j = 0; j < 3; ++j)
                    plane_n[j] = normals[seed * 3 + static_cast<std::size_t>(j)];
                surface = st.Ent("PLANE(''," +
                                 StepStream::Ref(placement(
                                     &xyz[static_cast<std::size_t>(
                                              nodes[static_cast<std::size_t>(start[seed])]) *
                                          3],
                                     plane_n)) +
                                 ")");
            } else {
                const double* p = row->mParams;
                same_sense = row->mSense;
                switch (row->mType) {
                    case 1:
                        for (int j = 0; j < 3; ++j)
                            plane_n[j] = (same_sense ? 1.0 : -1.0) * p[3 + j];
                        surface = st.Ent("PLANE(''," + StepStream::Ref(placement(p, p + 3)) + ")");
                        break;
                    case 2: {
                        const double z[3] = {0.0, 0.0, 1.0};
                        surface =
                            st.Ent("SPHERICAL_SURFACE(''," + StepStream::Ref(placement(p, z)) +
                                   "," + step_real(p[6]) + ")");
                        break;
                    }
                    case 3:
                        surface = st.Ent("CYLINDRICAL_SURFACE(''," +
                                         StepStream::Ref(placement(p, p + 3)) + "," +
                                         step_real(p[6]) + ")");
                        break;
                    case 4: {
                        // Placed at the mean height of the patch's points, where
                        // its radius is positive.
                        double hm = 0.0;
                        std::size_t count = 0;
                        for (std::size_t f = 0; f < nf; ++f) {
                            if (group[f] != g)
                                continue;
                            for (std::int64_t k = start[f]; k < start[f + 1]; ++k) {
                                const double* q = &xyz[static_cast<std::size_t>(
                                                           nodes[static_cast<std::size_t>(k)]) *
                                                       3];
                                const double v[3] = {q[0] - p[0], q[1] - p[1], q[2] - p[2]};
                                hm += sb_dot(v, p + 3);
                                ++count;
                            }
                        }
                        hm /= static_cast<double>(std::max<std::size_t>(count, 1));
                        const double a = p[8] * kSbPi / 180.0;
                        const double at[3] = {p[0] + hm * p[3], p[1] + hm * p[4], p[2] + hm * p[5]};
                        surface =
                            st.Ent("CONICAL_SURFACE(''," + StepStream::Ref(placement(at, p + 3)) +
                                   "," + step_real(hm * std::tan(a)) + "," + step_real(a) + ")");
                        break;
                    }
                    default:
                        surface =
                            st.Ent("TOROIDAL_SURFACE(''," + StepStream::Ref(placement(p, p + 3)) +
                                   "," + step_real(p[6]) + "," + step_real(p[7]) + ")");
                        break;
                }
            }
            // The bounds: on a plane, the loop of largest area is the outer one.
            std::size_t outer = kv.second.size();
            if (planar) {
                double best = -std::numeric_limits<double>::infinity();
                for (std::size_t li = 0; li < kv.second.size(); ++li) {
                    double twice[3] = {0.0, 0.0, 0.0};
                    std::vector<std::int64_t> ring;
                    for (const auto& e : kv.second[li]) {
                        const SbSegment& s = segs[e.first];
                        const std::vector<std::int64_t>& vs = chains[s.mChain].mVerts;
                        if (e.second)
                            for (std::size_t i = s.mFrom; i < s.mTo; ++i)
                                ring.push_back(vs[i]);
                        else
                            for (std::size_t i = s.mTo; i > s.mFrom; --i)
                                ring.push_back(vs[i]);
                    }
                    for (std::size_t i = 0; i < ring.size(); ++i) {
                        const double* a = &pos[static_cast<std::size_t>(ring[i]) * 3];
                        const double* b =
                            &pos[static_cast<std::size_t>(ring[(i + 1) % ring.size()]) * 3];
                        twice[0] += a[1] * b[2] - a[2] * b[1];
                        twice[1] += a[2] * b[0] - a[0] * b[2];
                        twice[2] += a[0] * b[1] - a[1] * b[0];
                    }
                    const double area = sb_dot(twice, plane_n);
                    if (area > best) {
                        best = area;
                        outer = li;
                    }
                }
            }
            std::string bounds;
            for (std::size_t li = 0; li < kv.second.size(); ++li) {
                std::string list;
                for (const auto& e : kv.second[li]) {
                    const std::int64_t oe =
                        st.Ent("ORIENTED_EDGE('',*,*," + StepStream::Ref(segs[e.first].mId) +
                               (e.second ? ",.T.)" : ",.F.)"));
                    list += (list.empty() ? "" : ",") + StepStream::Ref(oe);
                }
                const std::int64_t el = st.Ent("EDGE_LOOP('',(" + list + "))");
                const std::int64_t bound =
                    st.Ent(std::string(li == outer ? "FACE_OUTER_BOUND" : "FACE_BOUND") + "(''," +
                           StepStream::Ref(el) + ",.T.)");
                bounds += (bounds.empty() ? "" : ",") + StepStream::Ref(bound);
            }
            const std::int64_t face =
                st.Ent("ADVANCED_FACE('',(" + bounds + ")," + StepStream::Ref(surface) +
                       (same_sense ? ",.T.)" : ",.F.)"));
            face_index.emplace(g, res.mNumFaces++);
            shell_faces[piece[f0]].push_back(face);
            if (!analytic) {
                ++res.mNumFacetedFaces;
            } else {
                switch (row->mType) {
                    case 1:
                        ++res.mNumPlane;
                        break;
                    case 2:
                        ++res.mNumSphere;
                        break;
                    case 3:
                        ++res.mNumCylinder;
                        break;
                    case 4:
                        ++res.mNumCone;
                        break;
                    default:
                        ++res.mNumTorus;
                        break;
                }
            }
        }
        std::vector<std::int64_t> solids;
        for (std::size_t c = 0; c < npieces; ++c) {
            if (shell_faces[c].empty())
                continue;
            std::string list;
            for (std::int64_t f : shell_faces[c])
                list += (list.empty() ? "" : ",") + StepStream::Ref(f);
            const std::int64_t shell = st.Ent("CLOSED_SHELL('',(" + list + "))");
            solids.push_back(st.Ent("MANIFOLD_SOLID_BREP('solid" + std::to_string(c + 1) + "'," +
                                    StepStream::Ref(shell) + ")"));
        }
        res.mNumSolids = static_cast<std::int64_t>(solids.size());
        st.Finish(solids);

        // Per facet: its face, type and patch; the largest deviation.
        res.mFaceOfFacet.assign(nf, -1);
        res.mTypeOfFacet.assign(nf, 0);
        res.mPatchOfFacet.assign(nf, -1);
        for (std::size_t f = 0; f < nf; ++f) {
            const std::int64_t g = group[f];
            auto it = face_index.find(g);
            res.mFaceOfFacet[f] = it == face_index.end() ? -1 : it->second;
            res.mPatchOfFacet[f] = label_of[f];
            if (analytic_group(g)) {
                const SbRow& row = rows[static_cast<std::size_t>(g)];
                res.mTypeOfFacet[f] = row.mType;
                for (std::int64_t k = start[f]; k < start[f + 1]; ++k)
                    res.mMaxDeviation = std::max(
                        res.mMaxDeviation,
                        primitive_point_distance(
                            row.mType, row.mParams,
                            &xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(k)]) *
                                 3]));
            }
        }
        res.mValid = true;
        return res;
    }
}

}  // namespace detail
}  // namespace meshioplusplus
