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

// The C++ twin of `_tessellation.py`. Every step mirrors the numpy reference in
// the same order -- the lattice, the shape functions, the point keys and the
// order points are first created in -- because two choices depend on that
// order: which corner a sub-quad's diagonal starts from (the first-created
// one) and which cell's evaluation a shared point keeps (the first cell's).
// Coordinates, weights and interpolated values then agree with the numpy twin
// within rounding (numpy's matrix products sum in their own order), and
// everything integer -- point order, connectivity, stencils, keys -- exactly.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/tessellate.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"
// Project includes (private, not installed)
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

constexpr std::size_t kTsWidth = kTessellateStencilWidth;

// --- the five curved types ----------------------------------------------------

enum class TsBase { Triangle, Quad, Tetra, Hexahedron };

struct TsType {
    TsBase mBase = TsBase::Triangle;
    CellType mBaseType = CellType::Triangle;
    int mNumCorners = 0;
    int mNumNodes = 0;
    int mDim = 0;
    bool mFaces = false;  // quad9 / hexahedron27: face-centre nodes
    bool mBody = false;   // hexahedron27: a body-centre node
    bool mQuad8 = false;  // the serendipity basis
};

bool ts_curved_type(std::string_view Name, TsType& rOut) {
    if (Name == "triangle6") {
        rOut = {TsBase::Triangle, CellType::Triangle, 3, 6, 2, false, false, false};
    } else if (Name == "quad8") {
        rOut = {TsBase::Quad, CellType::Quad, 4, 8, 2, false, false, true};
    } else if (Name == "quad9") {
        rOut = {TsBase::Quad, CellType::Quad, 4, 9, 2, true, false, false};
    } else if (Name == "tetra10") {
        rOut = {TsBase::Tetra, CellType::Tetra, 4, 10, 3, false, false, false};
    } else if (Name == "hexahedron27") {
        rOut = {TsBase::Hexahedron, CellType::Hexahedron, 8, 27, 3, true, true, false};
    } else {
        return false;
    }
    return true;
}

bool ts_is_simplex(TsBase Base) {
    return Base == TsBase::Triangle || Base == TsBase::Tetra;
}

// The reference elements of `_REFERENCE_CORNERS`: the unit simplex, and
// [-1, 1]^d for the tensor bases.
std::vector<std::array<double, 3>> ts_reference_corners(TsBase Base) {
    switch (Base) {
        case TsBase::Triangle:
            return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
        case TsBase::Quad:
            return {{-1.0, -1.0, 0.0}, {1.0, -1.0, 0.0}, {1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}};
        case TsBase::Tetra:
            return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        case TsBase::Hexahedron:
        default:
            return {{-1.0, -1.0, -1.0}, {1.0, -1.0, -1.0}, {1.0, 1.0, -1.0}, {-1.0, 1.0, -1.0},
                    {-1.0, -1.0, 1.0},  {1.0, -1.0, 1.0},  {1.0, 1.0, 1.0},  {-1.0, 1.0, 1.0}};
    }
}

// The quad faces of a hexahedron in cyclic corner order; a face's lattice
// points are keyed by walking its corners around this cycle.
constexpr std::array<std::array<int, 4>, 6> kTsHexFaces = {{
    {0, 1, 2, 3},
    {4, 5, 6, 7},
    {0, 1, 5, 4},
    {1, 2, 6, 5},
    {2, 3, 7, 6},
    {3, 0, 4, 7},
}};

// --- one curved type's reference lattice ------------------------------------

struct TsLattice {
    TsType mType;
    int mDim = 0;
    std::int64_t mN = 1;           // the key quantum: divisions per edge
    std::vector<double> mXi;       // (M, dim)
    std::vector<int> mCells;       // (K, corners)
    std::vector<double> mCornerW;  // (M, corners)
    std::vector<double> mNodeW;    // (M, nodes)
    std::size_t NumPoints() const { return mXi.size() / static_cast<std::size_t>(mDim); }
};

// Node parametric coordinates: corners, edge midpoints (in the mid-node slot
// order of `cell_refine_edges`), face centres (`cell_refine_quad_faces`) and
// the body centre.
std::vector<std::array<double, 3>> ts_node_coords(const TsType& rType) {
    const auto corners = ts_reference_corners(rType.mBase);
    std::vector<std::array<double, 3>> out(corners.begin(), corners.end());
    for (const auto& e : detail::cell_refine_edges(rType.mBaseType)) {
        std::array<double, 3> m{};
        for (int a = 0; a < 3; ++a)
            m[a] = 0.5 * (corners[e[0]][a] + corners[e[1]][a]);
        out.push_back(m);
    }
    if (rType.mFaces) {
        for (const auto& f : detail::cell_refine_quad_faces(rType.mBaseType)) {
            std::array<double, 3> m{};
            for (int a = 0; a < 3; ++a)
                m[a] = (((corners[f[0]][a] + corners[f[1]][a]) + corners[f[2]][a]) +
                        corners[f[3]][a]) /
                       4.0;
            out.push_back(m);
        }
    }
    if (rType.mBody) {
        std::array<double, 3> m{};
        for (int a = 0; a < 3; ++a) {
            double s = 0.0;
            for (const auto& c : corners)
                s += c[a];
            m[a] = s / static_cast<double>(corners.size());
        }
        out.push_back(m);
    }
    if (static_cast<int>(out.size()) != rType.mNumNodes)
        throw std::logic_error("tessellate: node table does not match the cell type");
    return out;
}

double ts_lagrange1d(double T, double XiVal) {
    if (XiVal < -0.5)
        return 0.5 * T * (T - 1.0);
    if (XiVal > 0.5)
        return 0.5 * T * (T + 1.0);
    return 1.0 - T * T;
}

// Corner weights at one reference point: barycentric for a simplex, the
// multilinear tensor weight otherwise.
void ts_corner_weights(const TsType& rType, const double* pXi, double* pOut) {
    const int dim = rType.mDim;
    if (ts_is_simplex(rType.mBase)) {
        double sum = 0.0;
        for (int a = 0; a < dim; ++a)
            sum += pXi[a];
        pOut[0] = 1.0 - sum;
        for (int a = 0; a < dim; ++a)
            pOut[a + 1] = pXi[a];
        return;
    }
    const auto corners = ts_reference_corners(rType.mBase);
    for (int c = 0; c < rType.mNumCorners; ++c) {
        double w = 1.0;
        for (int a = 0; a < dim; ++a)
            w *= corners[c][a] > 0 ? (1.0 + pXi[a]) * 0.5 : (1.0 - pXi[a]) * 0.5;
        pOut[c] = w;
    }
}

// Quadratic shape functions at one reference point.
void ts_shape_functions(const TsType& rType, const std::vector<std::array<double, 3>>& rNodes,
                        const double* pXi, double* pOut) {
    const int dim = rType.mDim;
    if (rType.mQuad8) {
        const double x = pXi[0], y = pXi[1];
        for (int k = 0; k < 8; ++k) {
            const double xk = rNodes[k][0], yk = rNodes[k][1];
            if (std::fabs(xk) < 0.5)
                pOut[k] = 0.5 * (1.0 - x * x) * (1.0 + y * yk);
            else if (std::fabs(yk) < 0.5)
                pOut[k] = 0.5 * (1.0 + x * xk) * (1.0 - y * y);
            else
                pOut[k] = 0.25 * (1.0 + x * xk) * (1.0 + y * yk) * (x * xk + y * yk - 1.0);
        }
        return;
    }
    if (!ts_is_simplex(rType.mBase)) {
        for (int k = 0; k < rType.mNumNodes; ++k) {
            double w = 1.0;
            for (int a = 0; a < dim; ++a)
                w *= ts_lagrange1d(pXi[a], rNodes[k][a]);
            pOut[k] = w;
        }
        return;
    }
    std::array<double, 4> lam{};
    ts_corner_weights(rType, pXi, lam.data());
    const int nc = rType.mNumCorners;
    for (int i = 0; i < nc; ++i)
        pOut[i] = lam[i] * (2.0 * lam[i] - 1.0);
    const auto& edges = detail::cell_refine_edges(rType.mBaseType);
    for (std::size_t k = 0; k < edges.size(); ++k)
        pOut[nc + static_cast<int>(k)] = 4.0 * lam[edges[k][0]] * lam[edges[k][1]];
}

// A simplex base's `levels` recursive same-type splits, with points numbered in
// the order they are first created and deduplicated by their exact dyadic
// coordinates.
void ts_recursive_lattice(const TsType& rType, std::int64_t Levels, TsLattice& rOut) {
    const int dim = rType.mDim;
    const double scale = std::ldexp(1.0, static_cast<int>(std::max<std::int64_t>(Levels, 0)));
    std::vector<std::array<double, 3>> points;
    std::unordered_map<std::uint64_t, int> index_of;
    // Coordinates are integers in [0, 2^levels] once scaled; levels is capped
    // well below 21, so three of them pack into one 64-bit key.
    auto get_index = [&](const std::array<double, 3>& rP) {
        std::uint64_t key = 0;
        for (int a = 0; a < dim; ++a)
            key = (key << 21) | static_cast<std::uint64_t>(std::nearbyint(rP[a] * scale));
        auto it = index_of.find(key);
        if (it != index_of.end())
            return it->second;
        const int idx = static_cast<int>(points.size());
        index_of.emplace(key, idx);
        points.push_back(rP);
        return idx;
    };
    auto mid = [](const std::array<double, 3>& rA, const std::array<double, 3>& rB) {
        return std::array<double, 3>{0.5 * (rA[0] + rB[0]), 0.5 * (rA[1] + rB[1]),
                                     0.5 * (rA[2] + rB[2])};
    };
    const auto corners = ts_reference_corners(rType.mBase);
    const int nc = rType.mNumCorners;
    std::vector<int> cells;
    for (int i = 0; i < nc; ++i)
        cells.push_back(get_index(corners[i]));
    for (std::int64_t level = 0; level < Levels; ++level) {
        std::vector<int> next;
        for (std::size_t c = 0; c < cells.size(); c += nc) {
            std::vector<std::array<double, 3>> cp(nc);
            for (int i = 0; i < nc; ++i)
                cp[i] = points[cells[c + i]];
            std::vector<std::array<std::array<double, 3>, 4>> children;
            if (rType.mBase == TsBase::Triangle) {
                const auto m01 = mid(cp[0], cp[1]), m12 = mid(cp[1], cp[2]),
                           m20 = mid(cp[2], cp[0]);
                children = {{cp[0], m01, m20, {}},
                            {cp[1], m12, m01, {}},
                            {cp[2], m20, m12, {}},
                            {m01, m12, m20, {}}};
            } else {
                const auto m01 = mid(cp[0], cp[1]), m12 = mid(cp[1], cp[2]);
                const auto m02 = mid(cp[0], cp[2]), m03 = mid(cp[0], cp[3]);
                const auto m13 = mid(cp[1], cp[3]), m23 = mid(cp[2], cp[3]);
                children = {{cp[0], m01, m02, m03}, {cp[1], m01, m12, m13}, {cp[2], m02, m12, m23},
                            {cp[3], m03, m13, m23}, {m01, m23, m02, m03},   {m01, m23, m03, m13},
                            {m01, m23, m13, m12},   {m01, m23, m12, m02}};
            }
            for (const auto& child : children)
                for (int i = 0; i < nc; ++i)
                    next.push_back(get_index(child[i]));
        }
        cells = std::move(next);
    }
    rOut.mXi.resize(points.size() * dim);
    for (std::size_t p = 0; p < points.size(); ++p)
        for (int a = 0; a < dim; ++a)
            rOut.mXi[p * dim + a] = points[p][a];
    rOut.mCells = std::move(cells);
}

// A tensor base's (n + 1)-per-axis grid, numbered and connected exactly as the
// numpy twin's meshgrid ravel and index lambdas do.
void ts_tensor_lattice(const TsType& rType, std::int64_t Levels, TsLattice& rOut) {
    const std::int64_t n = std::max<std::int64_t>(Levels, 1);
    const std::size_t np1 = static_cast<std::size_t>(n + 1);
    std::vector<double> t(np1);
    const double step = 2.0 / static_cast<double>(n);
    for (std::size_t i = 0; i < np1; ++i)
        t[i] = static_cast<double>(i) * step + -1.0;
    t[np1 - 1] = 1.0;
    std::vector<double>& xi = rOut.mXi;
    std::vector<int>& cells = rOut.mCells;
    if (rType.mBase == TsBase::Quad) {
        for (std::size_t i = 0; i < np1; ++i)
            for (std::size_t j = 0; j < np1; ++j) {
                xi.push_back(t[i]);
                xi.push_back(t[j]);
            }
        auto idx = [np1](std::size_t i, std::size_t j) { return static_cast<int>(i + np1 * j); };
        for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i)
            for (std::size_t j = 0; j < static_cast<std::size_t>(n); ++j)
                for (int v : {idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)})
                    cells.push_back(v);
        return;
    }
    for (std::size_t i = 0; i < np1; ++i)
        for (std::size_t j = 0; j < np1; ++j)
            for (std::size_t k = 0; k < np1; ++k) {
                xi.push_back(t[i]);
                xi.push_back(t[j]);
                xi.push_back(t[k]);
            }
    auto idx3 = [np1](std::size_t i, std::size_t j, std::size_t k) {
        return static_cast<int>(i + np1 * (j + np1 * k));
    };
    const std::size_t nn = static_cast<std::size_t>(n);
    for (std::size_t i = 0; i < nn; ++i)
        for (std::size_t j = 0; j < nn; ++j)
            for (std::size_t k = 0; k < nn; ++k)
                for (int v : {idx3(i, j, k), idx3(i + 1, j, k), idx3(i + 1, j + 1, k),
                              idx3(i, j + 1, k), idx3(i, j, k + 1), idx3(i + 1, j, k + 1),
                              idx3(i + 1, j + 1, k + 1), idx3(i, j + 1, k + 1)})
                    cells.push_back(v);
}

TsLattice ts_build_lattice(const TsType& rType, std::int64_t Levels) {
    TsLattice lat;
    lat.mType = rType;
    lat.mDim = rType.mDim;
    if (ts_is_simplex(rType.mBase)) {
        ts_recursive_lattice(rType, Levels, lat);
        lat.mN = std::int64_t{1} << std::max<std::int64_t>(Levels, 0);
    } else {
        ts_tensor_lattice(rType, Levels, lat);
        lat.mN = std::max<std::int64_t>(Levels, 1);
    }
    const auto nodes = ts_node_coords(rType);
    const std::size_t m = lat.NumPoints();
    lat.mCornerW.resize(m * rType.mNumCorners);
    lat.mNodeW.resize(m * rType.mNumNodes);
    for (std::size_t p = 0; p < m; ++p) {
        const double* xi = lat.mXi.data() + p * lat.mDim;
        ts_corner_weights(rType, xi, lat.mCornerW.data() + p * rType.mNumCorners);
        ts_shape_functions(rType, nodes, xi, lat.mNodeW.data() + p * rType.mNumNodes);
    }
    return lat;
}

// --- point keys ---------------------------------------------------------------

// A point's identity: a short tuple of integers, compared lexicographically
// with a proper prefix first -- Python's tuple order, which the output point
// table is sorted by.
struct TsKey {
    std::array<std::int64_t, 8> mV{};
    std::uint8_t mLen = 0;
    void Push(std::int64_t V) { mV[mLen++] = V; }
    bool operator==(const TsKey& rO) const {
        return mLen == rO.mLen && std::equal(mV.begin(), mV.begin() + mLen, rO.mV.begin());
    }
    bool operator<(const TsKey& rO) const {
        return std::lexicographical_compare(mV.begin(), mV.begin() + mLen, rO.mV.begin(),
                                            rO.mV.begin() + rO.mLen);
    }
};

struct TsKeyHash {
    std::size_t operator()(const TsKey& rK) const {
        std::uint64_t h = 1469598103934665603ULL ^ rK.mLen;
        for (std::uint8_t i = 0; i < rK.mLen; ++i) {
            h ^= static_cast<std::uint64_t>(rK.mV[i]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            h *= 1099511628211ULL;
        }
        return static_cast<std::size_t>(h);
    }
};

std::int64_t ts_quantize(double V) {
    return static_cast<std::int64_t>(std::nearbyint(V));
}

// The shared point table, in creation order.
struct TsPoints {
    int mDim = 3;
    std::unordered_map<TsKey, std::int64_t, TsKeyHash> mIndex;
    std::vector<TsKey> mKeys;
    std::vector<double> mXyz;            // (P, dim)
    std::vector<std::int64_t> mStencil;  // (P, 27)
    std::vector<double> mWeights;        // (P, 27)

    // Returns the point's index, creating it (by calling `rMake`, which writes
    // its coordinates, stencil and weights) the first time the key is seen.
    template <class TMake>
    std::int64_t Resolve(const TsKey& rKey, TMake&& rMake) {
        auto it = mIndex.find(rKey);
        if (it != mIndex.end())
            return it->second;
        const std::int64_t idx = static_cast<std::int64_t>(mKeys.size());
        mIndex.emplace(rKey, idx);
        mKeys.push_back(rKey);
        mXyz.resize(mXyz.size() + mDim, 0.0);
        mStencil.resize(mStencil.size() + kTsWidth, -1);
        mWeights.resize(mWeights.size() + kTsWidth, 0.0);
        const std::size_t u = static_cast<std::size_t>(idx);
        rMake(mXyz.data() + u * mDim, mStencil.data() + u * kTsWidth,
              mWeights.data() + u * kTsWidth);
        return idx;
    }
};

// One curved cell's nodes, read from the source mesh.
struct TsCell {
    const std::int64_t* pNodes = nullptr;  // every node id, corners first
    std::vector<double> mXyz;              // (nodes, dim), Float64
    std::int64_t mGlobal = 0;              // global block-major cell index
};

// Write a stencil (nodes whose weight exceeds 1e-12) and the weighted
// position for a combination `pW` of the cell's nodes.
void ts_fill_point(const TsLattice& rLat, const TsCell& rCell, const double* pW, int Dim,
                   double* pXyz, std::int64_t* pStencil, double* pWeights) {
    const int nn = rLat.mType.mNumNodes;
    for (int a = 0; a < Dim; ++a) {
        double s = 0.0;
        for (int k = 0; k < nn; ++k)
            s += pW[k] * rCell.mXyz[static_cast<std::size_t>(k) * Dim + a];
        pXyz[a] = s;
    }
    std::size_t w = 0;
    for (int k = 0; k < nn && w < kTsWidth; ++k) {
        if (std::fabs(pW[k]) > 1e-12) {
            pStencil[w] = rCell.pNodes[k];
            pWeights[w] = pW[k];
            ++w;
        }
    }
}

// The canonical traversal of a corner cycle: start at the smallest id and walk
// toward its smaller neighbour.
void ts_canonical_cycle(const std::int64_t* pIds, int N, std::int64_t* pOut) {
    int start = 0;
    for (int i = 1; i < N; ++i)
        if (pIds[i] < pIds[start])
            start = i;
    const std::int64_t nxt = pIds[(start + 1) % N];
    const std::int64_t prv = pIds[(start - 1 + N) % N];
    const int step = nxt < prv ? 1 : -1;
    for (int k = 0; k < N; ++k)
        pOut[k] = pIds[((start + step * k) % N + N) % N];
}

// The key of reference lattice point `Local` of one cell.
TsKey ts_point_key(const TsLattice& rLat, std::size_t Local, const TsCell& rCell) {
    const TsType& type = rLat.mType;
    const int nc = type.mNumCorners;
    const double* cw = rLat.mCornerW.data() + Local * nc;
    std::array<int, 8> nz{};
    int nnz = 0;
    for (int c = 0; c < nc; ++c)
        if (std::fabs(cw[c]) > 1e-9)
            nz[nnz++] = c;
    const double n = static_cast<double>(rLat.mN);
    const std::int64_t* gid = rCell.pNodes;
    TsKey key;
    if (nnz == 1) {
        key.Push(0);
        key.Push(gid[nz[0]]);
    } else if (nnz == 2) {
        std::int64_t ga = gid[nz[0]], gb = gid[nz[1]];
        double wb = cw[nz[1]];
        if (ga > gb) {
            std::swap(ga, gb);
            wb = cw[nz[0]];
        }
        key.Push(1);
        key.Push(ga);
        key.Push(gb);
        key.Push(ts_quantize(wb * n));
    } else if ((type.mBase == TsBase::Hexahedron && nnz == 4) ||
               (type.mBase == TsBase::Tetra && nnz == 3)) {
        // A face point. A hexahedron's face corners are taken in their cyclic
        // order around the face, not in ascending local order (which, for four
        // of the six faces, is not a cycle at all).
        std::array<int, 4> order{};
        if (nnz == 4) {
            for (const auto& f : kTsHexFaces) {
                std::array<int, 4> sorted = f;
                std::sort(sorted.begin(), sorted.end());
                if (std::equal(sorted.begin(), sorted.end(), nz.begin())) {
                    order = f;
                    break;
                }
            }
        } else {
            for (int i = 0; i < 3; ++i)
                order[i] = nz[i];
        }
        std::array<std::int64_t, 4> ids{}, canon{};
        for (int i = 0; i < nnz; ++i)
            ids[i] = gid[order[i]];
        ts_canonical_cycle(ids.data(), nnz, canon.data());
        auto w_of = [&](std::int64_t g) {
            // The weight of the corner whose global id is `g` (the last such
            // corner, as the numpy twin's dict keeps).
            double w = 0.0;
            for (int i = 0; i < nnz; ++i)
                if (gid[order[i]] == g)
                    w = cw[order[i]];
            return w;
        };
        double u, v;
        if (nnz == 3) {
            u = w_of(canon[1]);
            v = w_of(canon[2]);
        } else {
            u = w_of(canon[1]) + w_of(canon[2]);
            v = w_of(canon[2]) + w_of(canon[3]);
        }
        key.Push(2);
        key.Push(nnz);
        for (int i = 0; i < nnz; ++i)
            key.Push(canon[i]);
        key.Push(ts_quantize(u * n));
        key.Push(ts_quantize(v * n));
    } else {
        key.Push(3);
        key.Push(0);
        key.Push(rCell.mGlobal);
        const double* xi = rLat.mXi.data() + Local * rLat.mDim;
        for (int a = 0; a < rLat.mDim; ++a)
            key.Push(ts_quantize(xi[a] * n));
    }
    return key;
}

// Split a quad (point indices in creation order) along the diagonal from its
// smallest index to the corner two steps away.
void ts_quad_split(const std::array<std::int64_t, 4>& rIds, std::vector<std::int64_t>& rOut,
                   std::int64_t Apex) {
    int p = 0;
    for (int i = 1; i < 4; ++i)
        if (rIds[i] < rIds[p])
            p = i;
    const std::int64_t a = rIds[p], b = rIds[(p + 1) % 4], c = rIds[(p + 2) % 4],
                       d = rIds[(p + 3) % 4];
    for (std::int64_t v : {a, b, c})
        rOut.push_back(v);
    if (Apex >= 0)
        rOut.push_back(Apex);
    for (std::int64_t v : {a, c, d})
        rOut.push_back(v);
    if (Apex >= 0)
        rOut.push_back(Apex);
}

// Tessellate one curved cell; append its children to `rOut` (3 or 4 ids each).
void ts_tessellate_cell(const TsLattice& rLat, const TsCell& rCell, TsPoints& rPoints,
                        std::vector<std::int64_t>& rOut, std::int64_t& rNumChildren) {
    const int dim = rPoints.mDim;
    const int nn = rLat.mType.mNumNodes;
    const std::size_t m = rLat.NumPoints();
    std::vector<std::int64_t> resolved(m);
    for (std::size_t i = 0; i < m; ++i) {
        const TsKey key = ts_point_key(rLat, i, rCell);
        resolved[i] = rPoints.Resolve(key, [&](double* pXyz, std::int64_t* pS, double* pW) {
            ts_fill_point(rLat, rCell, rLat.mNodeW.data() + i * nn, dim, pXyz, pS, pW);
        });
    }
    const int nc = rLat.mType.mNumCorners;
    const std::size_t ncells = rLat.mCells.size() / static_cast<std::size_t>(nc);
    const TsBase base = rLat.mType.mBase;
    rNumChildren = 0;
    if (ts_is_simplex(base)) {
        for (std::size_t c = 0; c < ncells; ++c)
            for (int i = 0; i < nc; ++i)
                rOut.push_back(resolved[rLat.mCells[c * nc + i]]);
        rNumChildren = static_cast<std::int64_t>(ncells);
        return;
    }
    if (base == TsBase::Quad) {
        for (std::size_t c = 0; c < ncells; ++c) {
            std::array<std::int64_t, 4> ids{};
            for (int i = 0; i < 4; ++i)
                ids[i] = resolved[rLat.mCells[c * 4 + i]];
            ts_quad_split(ids, rOut, -1);
        }
        rNumChildren = static_cast<std::int64_t>(2 * ncells);
        return;
    }
    // Hexahedron: 12 tetrahedra per sub-hexahedron around a centroid apex.
    std::vector<double> apex_w(nn);
    for (std::size_t c = 0; c < ncells; ++c) {
        const int* local = rLat.mCells.data() + c * 8;
        TsKey apex_key;
        apex_key.Push(3);
        apex_key.Push(1);
        apex_key.Push(rCell.mGlobal);
        apex_key.Push(static_cast<std::int64_t>(c));
        for (int a = 0; a < 3; ++a) {
            double s = 0.0;
            for (int i = 0; i < 8; ++i)
                s += rLat.mXi[static_cast<std::size_t>(local[i]) * 3 + a];
            apex_key.Push(ts_quantize(s / 8.0 * static_cast<double>(rLat.mN)));
        }
        for (int k = 0; k < nn; ++k) {
            double s = 0.0;
            for (int i = 0; i < 8; ++i)
                s += rLat.mNodeW[static_cast<std::size_t>(local[i]) * nn + k];
            apex_w[k] = s / 8.0;
        }
        const std::int64_t apex =
            rPoints.Resolve(apex_key, [&](double* pXyz, std::int64_t* pS, double* pW) {
                ts_fill_point(rLat, rCell, apex_w.data(), dim, pXyz, pS, pW);
            });
        for (const auto& f : kTsHexFaces) {
            std::array<std::int64_t, 4> ids{};
            for (int i = 0; i < 4; ++i)
                ids[i] = resolved[local[f[i]]];
            ts_quad_split(ids, rOut, apex);
        }
    }
    rNumChildren = static_cast<std::int64_t>(12 * ncells);
}

// --- a pass-through block, staged with source point ids ----------------------

struct TsBlock {
    std::string mType;
    std::size_t mNodesPerCell = 0;
    std::vector<std::int64_t> mConn;  // rectangular, or the ragged CSR flat list
    std::vector<std::int64_t> mRowOffsets;
    std::vector<std::int64_t> mFaceOffsets;
    bool mIsRagged = false;
    bool mIsPolyhedron = false;
    std::vector<std::int64_t> mSourceCell;
    std::vector<std::int64_t> mSubIndex;
};

// Stage a pass-through block, creating its points (in the numpy twin's order:
// row by row) and storing their creation indices.
TsBlock ts_stage_passthrough(const Mesh::CellView& rBlock, const detail::DoubleView& rPoints,
                             TsPoints& rTable, std::int64_t Base) {
    TsBlock out;
    out.mType = std::string(rBlock.Type());
    const int dim = rTable.mDim;
    auto kept = [&](std::int64_t gid) {
        TsKey key;
        key.Push(0);
        key.Push(gid);
        return rTable.Resolve(key, [&](double* pXyz, std::int64_t* pS, double* pW) {
            for (int a = 0; a < dim; ++a)
                pXyz[a] = rPoints[static_cast<std::size_t>(gid) * dim + a];
            pS[0] = gid;
            pW[0] = 1.0;
        });
    };
    const std::size_t nc = rBlock.NumCells();
    if (rBlock.IsPolyhedron()) {
        out.mIsRagged = out.mIsPolyhedron = true;
        out.mRowOffsets.push_back(0);
        out.mFaceOffsets.push_back(0);
        for (std::size_t c = 0; c < nc; ++c) {
            for (std::size_t f = 0; f < rBlock.NumFaces(c); ++f) {
                const auto face = rBlock.Face(c, f);
                for (std::size_t i = 0; i < face.second; ++i)
                    out.mConn.push_back(kept(face.first[i]));
                out.mRowOffsets.push_back(static_cast<std::int64_t>(out.mConn.size()));
            }
            out.mFaceOffsets.push_back(static_cast<std::int64_t>(out.mRowOffsets.size() - 1));
        }
    } else if (rBlock.IsRagged()) {
        out.mIsRagged = true;
        out.mRowOffsets.push_back(0);
        for (std::size_t c = 0; c < nc; ++c) {
            const std::int64_t* row = rBlock.Row(c);
            for (std::size_t i = 0; i < rBlock.RowSize(c); ++i)
                out.mConn.push_back(kept(row[i]));
            out.mRowOffsets.push_back(static_cast<std::int64_t>(out.mConn.size()));
        }
    } else {
        const detail::Int64View conn(rBlock.Conn());
        const std::size_t npc = rBlock.NodesPerCell();
        out.mNodesPerCell = npc;
        out.mConn.reserve(nc * npc);
        for (std::size_t i = 0; i < nc * npc; ++i)
            out.mConn.push_back(kept(conn[i]));
    }
    out.mSourceCell.resize(nc);
    for (std::size_t c = 0; c < nc; ++c)
        out.mSourceCell[c] = Base + static_cast<std::int64_t>(c);
    out.mSubIndex.assign(nc, 0);
    return out;
}

void ts_emit_block(Mesh& rOut, TsBlock& rBlock, const std::vector<std::int64_t>& rRemap) {
    for (std::int64_t& id : rBlock.mConn)
        id = rRemap[static_cast<std::size_t>(id)];
    if (rBlock.mIsPolyhedron) {
        rOut.AddPolyhedronBlock(rBlock.mType, std::move(rBlock.mConn),
                                std::move(rBlock.mRowOffsets), std::move(rBlock.mFaceOffsets));
    } else if (rBlock.mIsRagged) {
        rOut.AddPolygonBlock(rBlock.mType, std::move(rBlock.mConn), std::move(rBlock.mRowOffsets));
    } else {
        const std::size_t npc = rBlock.mNodesPerCell;
        const std::size_t ncells = npc == 0 ? rBlock.mSourceCell.size() : rBlock.mConn.size() / npc;
        NDArray conn = NDArray::Uninit(DType::Int64, {ncells, npc});
        if (!rBlock.mConn.empty())
            std::memcpy(conn.Data(), rBlock.mConn.data(),
                        rBlock.mConn.size() * sizeof(std::int64_t));
        rOut.AddCellBlock(rBlock.mType, std::move(conn));
    }
}

NDArray ts_int64_array(const std::vector<std::int64_t>& rValues) {
    NDArray out = NDArray::Uninit(DType::Int64, {rValues.size()});
    if (!rValues.empty())
        std::memcpy(out.Data(), rValues.data(), rValues.size() * sizeof(std::int64_t));
    return out;
}

bool ts_reserved_name(const std::string& rName) {
    return rName.rfind("tessellate:", 0) == 0;
}

// Facet-use counts over every tetra block of the output.
void ts_watertight(const Mesh& rOut, TessellateResult& rRes) {
    std::vector<std::array<std::int64_t, 3>> facets;
    for (const auto cb : rOut.CellRange()) {
        if (cb.IsRagged() || std::string_view(cb.Type()) != "tetra")
            continue;
        const detail::Int64View conn(cb.Conn());
        static constexpr int kFaces[4][3] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            for (const auto& f : kFaces) {
                std::array<std::int64_t, 3> key{conn[c * 4 + f[0]], conn[c * 4 + f[1]],
                                                conn[c * 4 + f[2]]};
                std::sort(key.begin(), key.end());
                facets.push_back(key);
            }
    }
    std::sort(facets.begin(), facets.end());
    for (std::size_t i = 0; i < facets.size();) {
        std::size_t j = i;
        while (j < facets.size() && facets[j] == facets[i])
            ++j;
        const std::size_t count = j - i;
        if (count == 1)
            ++rRes.mNumBoundaryFacets;
        else if (count == 2)
            ++rRes.mNumInteriorFacets;
        else
            ++rRes.mNumFacetsWithBadCount;
        i = j;
    }
}

}  // namespace

TessellateResult tessellate(const Mesh& rMesh, const TessellateOptions& rOptions) {
    if (rOptions.mLevels < 0)
        throw std::invalid_argument("tessellate: levels must be >= 0, got " +
                                    std::to_string(rOptions.mLevels));
    // A simplex lattice has 2^levels divisions per edge and a tetrahedron
    // 8^levels children: 2^20 per edge is far past any useful tessellation and
    // keeps the lattice keys inside 21 bits.
    if (rOptions.mLevels > 20)
        throw std::invalid_argument("tessellate: levels must be <= 20, got " +
                                    std::to_string(rOptions.mLevels));

    const std::size_t num_points = rMesh.NumPoints();
    const int dim = num_points ? static_cast<int>(rMesh.PointDim()) : 3;
    const std::size_t nblocks = rMesh.NumCellBlocks();
    for (const std::string& name : rMesh.PointDataNames()) {
        if (rOptions.mFields && !ts_reserved_name(name) &&
            detail::rows(rMesh.PointData(name)) != num_points)
            throw std::invalid_argument("tessellate: point data '" + name + "' has " +
                                        std::to_string(detail::rows(rMesh.PointData(name))) +
                                        " rows, expected " + std::to_string(num_points));
    }

    std::vector<std::int64_t> bases(nblocks + 1, 0);
    for (std::size_t b = 0; b < nblocks; ++b)
        bases[b + 1] = bases[b] + static_cast<std::int64_t>(rMesh.Cells(b).NumCells());

    const NDArray& points_arr = rMesh.Points();
    const detail::DoubleView points(points_arr);

    TsPoints table;
    table.mDim = dim;
    std::unordered_map<std::string, TsLattice> lattices;

    // Block slots in output order: pass-through blocks keep their input order;
    // curved children go to the trailing triangle and tetra blocks.
    std::vector<TsBlock> passthrough;
    std::vector<std::int64_t> tri_conn, tet_conn;
    std::vector<std::int64_t> tri_src, tri_sub, tet_src, tet_sub;

    TessellateResult res;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        TsType type;
        const bool curved = rOptions.mCurved && !cb.IsRagged() &&
                            ts_curved_type(std::string_view(cb.Type()), type) &&
                            cb.NodesPerCell() == static_cast<std::size_t>(type.mNumNodes);
        if (!curved) {
            res.mNumPassThroughSourceCells += static_cast<std::int64_t>(cb.NumCells());
            passthrough.push_back(ts_stage_passthrough(cb, points, table, bases[b]));
            continue;
        }
        res.mNumCurvedSourceCells += static_cast<std::int64_t>(cb.NumCells());
        auto it = lattices.find(std::string(cb.Type()));
        if (it == lattices.end())
            it = lattices.emplace(std::string(cb.Type()), ts_build_lattice(type, rOptions.mLevels))
                     .first;
        const TsLattice& lat = it->second;
        const detail::Int64View conn(cb.Conn());
        const int nn = type.mNumNodes;
        const bool to_tet = type.mDim == 3;
        std::vector<std::int64_t>& dest = to_tet ? tet_conn : tri_conn;
        std::vector<std::int64_t>& dest_src = to_tet ? tet_src : tri_src;
        std::vector<std::int64_t>& dest_sub = to_tet ? tet_sub : tri_sub;
        TsCell cell;
        cell.mXyz.resize(static_cast<std::size_t>(nn) * dim);
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            cell.pNodes = conn.Data() + c * nn;
            cell.mGlobal = bases[b] + static_cast<std::int64_t>(c);
            for (int k = 0; k < nn; ++k) {
                const std::int64_t g = cell.pNodes[k];
                if (g < 0 || static_cast<std::size_t>(g) >= num_points)
                    throw std::invalid_argument("tessellate: cell node " + std::to_string(g) +
                                                " is outside the point table");
                for (int a = 0; a < dim; ++a)
                    cell.mXyz[static_cast<std::size_t>(k) * dim + a] =
                        points[static_cast<std::size_t>(g) * dim + a];
            }
            std::int64_t nchildren = 0;
            ts_tessellate_cell(lat, cell, table, dest, nchildren);
            for (std::int64_t k = 0; k < nchildren; ++k) {
                dest_src.push_back(cell.mGlobal);
                dest_sub.push_back(k);
            }
        }
    }

    // The output point table, sorted by key.
    const std::size_t np = table.mKeys.size();
    std::vector<std::int64_t> order(np);
    for (std::size_t i = 0; i < np; ++i)
        order[i] = static_cast<std::int64_t>(i);
    std::sort(order.begin(), order.end(), [&](std::int64_t a, std::int64_t b) {
        return table.mKeys[static_cast<std::size_t>(a)] < table.mKeys[static_cast<std::size_t>(b)];
    });
    std::vector<std::int64_t> remap(np, -1);
    for (std::size_t i = 0; i < np; ++i)
        remap[static_cast<std::size_t>(order[i])] = static_cast<std::int64_t>(i);

    Mesh out;
    std::vector<std::int64_t> source_point(np, -1);
    {
        NDArray pts = NDArray::Uninit(DType::Float64, {np, static_cast<std::size_t>(dim)});
        double* dst = pts.As<double>();
        res.mStencil.assign(np * kTsWidth, -1);
        res.mWeights.assign(np * kTsWidth, 0.0);
        for (std::size_t i = 0; i < np; ++i) {
            const std::size_t old = static_cast<std::size_t>(order[i]);
            std::memcpy(dst + i * dim, table.mXyz.data() + old * dim, sizeof(double) * dim);
            std::memcpy(res.mStencil.data() + i * kTsWidth, table.mStencil.data() + old * kTsWidth,
                        sizeof(std::int64_t) * kTsWidth);
            std::memcpy(res.mWeights.data() + i * kTsWidth, table.mWeights.data() + old * kTsWidth,
                        sizeof(double) * kTsWidth);
            if (table.mKeys[old].mV[0] == 0)
                source_point[i] = table.mKeys[old].mV[1];
        }
        out.AssignPoints(std::move(pts));
    }

    // Cells: pass-through blocks, then the triangles, then the tetrahedra.
    std::vector<std::vector<std::int64_t>> block_src, block_sub;
    for (TsBlock& blk : passthrough) {
        block_src.push_back(blk.mSourceCell);
        block_sub.push_back(blk.mSubIndex);
        ts_emit_block(out, blk, remap);
    }
    auto emit_simplices = [&](const char* pType, std::size_t Arity,
                              std::vector<std::int64_t>& rConn, std::vector<std::int64_t>& rSrc,
                              std::vector<std::int64_t>& rSub) {
        if (rSrc.empty())
            return;
        NDArray conn = NDArray::Uninit(DType::Int64, {rSrc.size(), Arity});
        std::int64_t* dst = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < rConn.size(); ++i)
            dst[i] = remap[static_cast<std::size_t>(rConn[i])];
        out.AddCellBlock(pType, std::move(conn));
        block_src.push_back(std::move(rSrc));
        block_sub.push_back(std::move(rSub));
    };
    emit_simplices("triangle", 3, tri_conn, tri_src, tri_sub);
    emit_simplices("tetra", 4, tet_conn, tet_src, tet_sub);

    // Point data.
    out.AddPointData(kTessellateSourcePointName, ts_int64_array(source_point));
    if (rOptions.mFields) {
        for (const std::string& name : rMesh.PointDataNames()) {
            if (ts_reserved_name(name))
                continue;
            const NDArray& a = rMesh.PointData(name);
            std::vector<std::size_t> shape = a.Shape();
            if (shape.empty())
                shape.push_back(num_points);
            const std::size_t width = num_points == 0 ? 1 : a.Size() / num_points;
            shape[0] = np;
            NDArray g = NDArray::Uninit(DType::Float64, std::move(shape));
            double* dst = g.As<double>();
            const detail::DoubleView src(a);
            parallel_for(np, [&](std::size_t p) {
                const std::int64_t* st = res.mStencil.data() + p * kTsWidth;
                const double* w = res.mWeights.data() + p * kTsWidth;
                for (std::size_t c = 0; c < width; ++c) {
                    double s = 0.0;
                    for (std::size_t k = 0; k < kTsWidth; ++k)
                        if (w[k] != 0.0)
                            s += src[static_cast<std::size_t>(st[k]) * width + c] * w[k];
                    dst[p * width + c] = s;
                }
            });
            out.AddPointData(name, std::move(g));
        }
    }
    if (rOptions.mRecordStencil) {
        NDArray st = NDArray::Uninit(DType::Int64, {np, kTsWidth});
        NDArray wt = NDArray::Uninit(DType::Float64, {np, kTsWidth});
        if (np) {
            std::memcpy(st.Data(), res.mStencil.data(), res.mStencil.size() * sizeof(std::int64_t));
            std::memcpy(wt.Data(), res.mWeights.data(), res.mWeights.size() * sizeof(double));
        }
        out.AddPointData(kTessellateStencilName, std::move(st));
        out.AddPointData(kTessellateWeightsName, std::move(wt));
    }

    // Cell data.
    {
        std::vector<NDArray> src_blocks, sub_blocks;
        for (std::size_t b = 0; b < block_src.size(); ++b) {
            src_blocks.push_back(ts_int64_array(block_src[b]));
            sub_blocks.push_back(ts_int64_array(block_sub[b]));
        }
        if (!src_blocks.empty()) {
            out.AddCellData(kTessellateSourceCellName, std::move(src_blocks));
            out.AddCellData(kTessellateSubIndexName, std::move(sub_blocks));
        }
    }
    if (rOptions.mFields && !block_src.empty()) {
        for (const std::string& name : rMesh.CellDataNames()) {
            if (ts_reserved_name(name))
                continue;
            // Rows are copied byte for byte, so every block must share one dtype
            // and one row shape (the numpy twin, which can upcast, takes the rest).
            bool complete = nblocks > 0 && rMesh.CellDataNumBlocks(name) == nblocks;
            const NDArray* first = complete ? &rMesh.CellData(name, 0) : nullptr;
            for (std::size_t b = 0; complete && b < nblocks; ++b) {
                const NDArray& a = rMesh.CellData(name, b);
                std::vector<std::size_t> tail = a.Shape(), tail0 = first->Shape();
                complete = !tail.empty() && !tail0.empty() &&
                           detail::rows(a) == rMesh.Cells(b).NumCells() &&
                           a.Dtype() == first->Dtype() &&
                           std::equal(tail.begin() + 1, tail.end(), tail0.begin() + 1, tail0.end());
            }
            if (!complete) {
                log::warn(
                    "tessellate: dropping cell data '{}': its blocks do not share one dtype "
                    "and row shape with one row per cell",
                    name);
                continue;
            }
            std::size_t row_bytes = dtype_size(first->Dtype());
            for (std::size_t k = 1; k < first->Shape().size(); ++k)
                row_bytes *= first->Shape()[k];
            std::vector<NDArray> blocks;
            for (std::size_t ob = 0; ob < block_src.size(); ++ob) {
                const std::vector<std::int64_t>& srcs = block_src[ob];
                std::vector<std::size_t> shape = first->Shape();
                shape[0] = srcs.size();
                NDArray dst = NDArray::Uninit(first->Dtype(), std::move(shape));
                std::byte* d = dst.Data();
                for (std::size_t i = 0; i < srcs.size(); ++i) {
                    const std::int64_t g = srcs[i];
                    const std::size_t b = static_cast<std::size_t>(
                        std::upper_bound(bases.begin(), bases.end(), g) - bases.begin() - 1);
                    const NDArray& a = rMesh.CellData(name, b);
                    std::memcpy(d + i * row_bytes,
                                a.Data() + static_cast<std::size_t>(g - bases[b]) * row_bytes,
                                row_bytes);
                }
                blocks.push_back(std::move(dst));
            }
            out.AddCellData(name, std::move(blocks));
        }
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    // Regions: point regions follow the kept points; cell regions keep their
    // pass-through cells; side regions are dropped.
    if (rMesh.NumRegions()) {
        std::vector<std::int64_t> point_map(num_points, -1);
        for (std::size_t i = 0; i < np; ++i)
            if (source_point[i] >= 0)
                point_map[static_cast<std::size_t>(source_point[i])] = static_cast<std::int64_t>(i);
        // A source cell maps to an output cell only when it passed through.
        std::vector<std::int64_t> cell_map(static_cast<std::size_t>(bases[nblocks]), -1);
        std::int64_t out_cell = 0;
        for (const TsBlock& blk : passthrough)
            for (std::int64_t g : blk.mSourceCell)
                cell_map[static_cast<std::size_t>(g)] = out_cell++;
        std::vector<std::string> dropped;
        for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
            const Region& reg = rMesh.Region(r);
            if (reg.mKind == RegionKind::Side) {
                dropped.push_back(reg.mName + " (side)");
                continue;
            }
            const detail::Int64View ent(reg.mEntries);
            const std::vector<std::int64_t>& map =
                reg.mKind == RegionKind::Point ? point_map : cell_map;
            std::vector<std::int64_t> kept;
            for (std::size_t i = 0; i < reg.NumEntries(); ++i) {
                const std::int64_t e = ent[i];
                if (e >= 0 && static_cast<std::size_t>(e) < map.size() &&
                    map[static_cast<std::size_t>(e)] >= 0)
                    kept.push_back(map[static_cast<std::size_t>(e)]);
            }
            if (reg.mKind == RegionKind::Cell && kept.size() != reg.NumEntries())
                dropped.push_back(reg.mName + " (cell, entries on a tessellated cell)");
            out.AddRegion(meshioplusplus::Region(reg.mName, reg.mKind, reg.mDim, reg.mTag,
                                                 ts_int64_array(kept)));
        }
        if (!dropped.empty()) {
            std::string list;
            for (const std::string& d : dropped)
                list += (list.empty() ? "" : ", ") + d;
            log::warn("tessellate: dropped region entries: {}", list);
        }
    }

    res.mNumSourcePoints = static_cast<std::int64_t>(num_points);
    res.mNumSourceCells = bases[nblocks];
    res.mNumPoints = static_cast<std::int64_t>(np);
    std::int64_t ncells = 0;
    for (const auto& s : block_src)
        ncells += static_cast<std::int64_t>(s.size());
    res.mNumCells = ncells;
    ts_watertight(out, res);
    res.mMesh = std::move(out);
    return res;
}

}  // namespace meshioplusplus
