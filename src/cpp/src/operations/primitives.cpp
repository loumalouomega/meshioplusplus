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
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

constexpr double kPrimPi = 3.14159265358979323846;

// Marks a node whose azimuth is meaningless (a pole, a disk centre): a mid-edge
// node on an edge to it takes the other end's azimuth.
constexpr std::uint16_t kPrimDegenerate = 0x8000;

// How a node's three parameters become a position. Every generator builds its
// nodes in a parameter space where the mid-edge average is exact (a lattice
// coordinate, an angle, a unit vector); the map is applied once, at the end, to
// corners and mid-edge nodes alike, so a quadratic cell lands on the true
// geometry rather than on a chord.
enum class PrimMap { Identity, Squircle, Polar, Sphere, Unit };

struct PrimBlock {
    CellType mType = CellType::Vertex;
    std::size_t mWidth = 0;
    std::vector<std::int64_t> mConn;
};

struct PrimSoup {
    PrimMap mMap = PrimMap::Identity;
    double mRadius = 1.0;
    std::array<double, 3> mCenter = {{0.0, 0.0, 0.0}};
    std::vector<double> mParam;
    std::vector<std::uint16_t> mMask;
    std::vector<PrimBlock> mBlocks;

    std::int64_t AddNode(double A, double B, double C, std::uint16_t Mask) {
        mParam.push_back(A);
        mParam.push_back(B);
        mParam.push_back(C);
        mMask.push_back(Mask);
        return static_cast<std::int64_t>(mMask.size() - 1);
    }

    PrimBlock& Block(CellType Type, std::size_t Width) {
        mBlocks.emplace_back();
        mBlocks.back().mType = Type;
        mBlocks.back().mWidth = Width;
        return mBlocks.back();
    }
};

struct PrimRegionSpec {
    const char* mName;
    std::uint16_t mBit;
};

[[noreturn]] void prim_fail(const char* pShape, const std::string& rMessage) {
    throw std::invalid_argument(std::string("meshio++: ") + pShape + ": " + rMessage);
}

void prim_check_count(const char* pShape, const char* pWhat, std::int64_t Value) {
    if (Value < 0)
        prim_fail(pShape, std::string(pWhat) + " must not be negative, got " +
                              std::to_string(Value));
}

void prim_check_positive(const char* pShape, const char* pWhat, double Value) {
    if (!(Value > 0.0))
        prim_fail(pShape, std::string(pWhat) + " must be positive, got " + std::to_string(Value));
}

void prim_check_order(const char* pShape, int Order) {
    if (Order != 1 && Order != 2)
        prim_fail(pShape, "order must be 1 or 2, got " + std::to_string(Order));
}

[[noreturn]] void prim_bad_choice(const char* pShape, const char* pWhat, const std::string& rGot,
                                  const char* pChoices) {
    prim_fail(pShape, std::string("unknown ") + pWhat + " '" + rGot +
                          "' (expected one of: " + pChoices + ")");
}

// Saturating product: the budget check must not wrap on an absurd request.
std::int64_t prim_mul(std::int64_t A, std::int64_t B) {
    if (A == 0 || B == 0)
        return 0;
    if (A > std::numeric_limits<std::int64_t>::max() / B)
        return std::numeric_limits<std::int64_t>::max();
    return A * B;
}

void prim_check_budget(const char* pShape, std::int64_t Cells, std::int64_t MaxCells) {
    if (MaxCells <= 0 || Cells <= MaxCells)
        return;
    prim_fail(pShape, "the requested mesh has " + std::to_string(Cells) +
                          " cells, above the limit of " + std::to_string(MaxCells) +
                          " (raise max_cells or coarsen the counts)");
}

Mesh prim_empty() {
    Mesh out;
    out.AssignPoints(NDArray(DType::Float64, {std::size_t{0}, std::size_t{3}}));
    return out;
}

// ---------------------------------------------------------------------------
// Parameter space -> position, and the mid-edge average in parameter space.
// ---------------------------------------------------------------------------

void prim_position(const PrimSoup& rSoup, std::size_t Node, double* pOut) {
    const double* p = rSoup.mParam.data() + 3 * Node;
    const double r = rSoup.mRadius;
    const double cx = rSoup.mCenter[0], cy = rSoup.mCenter[1], cz = rSoup.mCenter[2];
    switch (rSoup.mMap) {
        case PrimMap::Identity:
            pOut[0] = p[0];
            pOut[1] = p[1];
            pOut[2] = p[2];
            break;
        case PrimMap::Squircle:
            pOut[0] = cx + r * p[0] * std::sqrt(1.0 - p[1] * p[1] / 2.0);
            pOut[1] = cy + r * p[1] * std::sqrt(1.0 - p[0] * p[0] / 2.0);
            pOut[2] = cz + p[2];
            break;
        case PrimMap::Polar:
            pOut[0] = cx + r * p[0] * std::cos(p[1]);
            pOut[1] = cy + r * p[0] * std::sin(p[1]);
            pOut[2] = cz + p[2];
            break;
        case PrimMap::Sphere:
            if (rSoup.mMask[Node] & kPrimDegenerate) {
                // A pole is exactly on the axis; sin(pi) would leave 1e-16.
                pOut[0] = cx;
                pOut[1] = cy;
                pOut[2] = p[0] == 0.0 ? cz + r : cz - r;
            } else {
                const double st = std::sin(p[0]);
                pOut[0] = cx + r * st * std::cos(p[1]);
                pOut[1] = cy + r * st * std::sin(p[1]);
                pOut[2] = cz + r * std::cos(p[0]);
            }
            break;
        case PrimMap::Unit:
            pOut[0] = cx + r * p[0];
            pOut[1] = cy + r * p[1];
            pOut[2] = cz + r * p[2];
            break;
    }
}

// The azimuth halfway between two, along the shorter arc; the end without a
// meaningful azimuth defers to the other.
double prim_mid_angle(double A, bool DegenA, double B, bool DegenB) {
    if (DegenA)
        return B;
    if (DegenB)
        return A;
    const double two_pi = 2.0 * kPrimPi;
    double d = B - A;
    d -= two_pi * std::floor(d / two_pi + 0.5);
    return A + 0.5 * d;
}

// Parameters of the node midway along edge (A, B), A < B.
void prim_mid(const PrimSoup& rSoup, std::size_t A, std::size_t B, double* pOut) {
    const double* a = rSoup.mParam.data() + 3 * A;
    const double* b = rSoup.mParam.data() + 3 * B;
    const bool da = (rSoup.mMask[A] & kPrimDegenerate) != 0;
    const bool db = (rSoup.mMask[B] & kPrimDegenerate) != 0;
    switch (rSoup.mMap) {
        case PrimMap::Identity:
        case PrimMap::Squircle:
            pOut[0] = (a[0] + b[0]) * 0.5;
            pOut[1] = (a[1] + b[1]) * 0.5;
            pOut[2] = (a[2] + b[2]) * 0.5;
            break;
        case PrimMap::Polar:
            pOut[0] = (a[0] + b[0]) * 0.5;
            pOut[1] = prim_mid_angle(a[1], da, b[1], db);
            pOut[2] = (a[2] + b[2]) * 0.5;
            break;
        case PrimMap::Sphere:
            pOut[0] = (a[0] + b[0]) * 0.5;
            pOut[1] = prim_mid_angle(a[1], da, b[1], db);
            pOut[2] = 0.0;
            break;
        case PrimMap::Unit: {
            const double x = a[0] + b[0], y = a[1] + b[1], z = a[2] + b[2];
            const double len = std::sqrt(x * x + y * y + z * z);
            pOut[0] = x / len;
            pOut[1] = y / len;
            pOut[2] = z / len;
            break;
        }
    }
}

std::uint64_t prim_edge_key(std::int64_t A, std::int64_t B) {
    return (static_cast<std::uint64_t>(A) << 32) | static_cast<std::uint64_t>(B);
}

CellType prim_elevated(CellType Type) {
    switch (Type) {
        case CellType::Triangle:
            return CellType::Triangle6;
        case CellType::Quad:
            return CellType::Quad8;
        case CellType::Tetra:
            return CellType::Tetra10;
        case CellType::Hexahedron:
            return CellType::Hexahedron20;
        case CellType::Wedge:
            return CellType::Wedge15;
        default:
            return Type;
    }
}

// Quadratic cells: mid-edge node `corners + e` of a cell is appended on first
// appearance, in cell-major order over `cell_refine_edges`, and shared by every
// cell on the edge.
void prim_elevate(const char* pShape, PrimSoup& rSoup) {
    if (rSoup.mMask.size() >= (std::size_t{1} << 31))
        prim_fail(pShape, "too many nodes for a quadratic mesh");
    std::unordered_map<std::uint64_t, std::int64_t> ids;
    for (PrimBlock& blk : rSoup.mBlocks) {
        const auto& edges = detail::cell_refine_edges(blk.mType);
        if (edges.empty())
            continue;
        const std::size_t ncells = blk.mConn.size() / blk.mWidth;
        const std::size_t width = blk.mWidth + edges.size();
        std::vector<std::int64_t> conn(ncells * width);
        for (std::size_t c = 0; c < ncells; ++c) {
            const std::int64_t* row = blk.mConn.data() + c * blk.mWidth;
            std::int64_t* dst = conn.data() + c * width;
            for (std::size_t k = 0; k < blk.mWidth; ++k)
                dst[k] = row[k];
            for (std::size_t e = 0; e < edges.size(); ++e) {
                std::int64_t a = row[edges[e][0]], b = row[edges[e][1]];
                if (a > b)
                    std::swap(a, b);
                const auto key = prim_edge_key(a, b);
                auto found = ids.find(key);
                if (found == ids.end()) {
                    double mid[3];
                    prim_mid(rSoup, static_cast<std::size_t>(a), static_cast<std::size_t>(b), mid);
                    const std::int64_t id = rSoup.AddNode(mid[0], mid[1], mid[2], 0);
                    found = ids.emplace(key, id).first;
                }
                dst[blk.mWidth + e] = found->second;
            }
        }
        blk.mConn = std::move(conn);
        blk.mWidth = width;
        blk.mType = prim_elevated(blk.mType);
    }
}

// ---------------------------------------------------------------------------
// Boundary regions: a boundary facet joins the first region whose bit every one
// of its corners carries.
// ---------------------------------------------------------------------------

void prim_add_regions(Mesh& rMesh, const PrimSoup& rSoup,
                      const std::vector<PrimRegionSpec>& rSpecs, int FacetDim) {
    if (rSpecs.empty())
        return;
    const detail::FacetIndex index(rMesh);
    std::vector<std::vector<std::pair<std::int64_t, std::int64_t>>> entries(rSpecs.size());
    std::int64_t cell = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto view = rMesh.Cells(b);
        const CellType type = cell_type_from_name(view.Type());
        const std::size_t ncells = view.NumCells();
        const std::size_t width = view.NodesPerCell();
        const std::int64_t* conn = view.Conn().As<std::int64_t>();
        for (std::size_t c = 0; c < ncells; ++c, ++cell) {
            const std::int64_t* row = conn + c * width;
            const std::size_t nfacets =
                FacetDim == 2 ? detail::cell_faces(type).size() : detail::cell_edges(type).size();
            for (std::size_t f = 0; f < nfacets; ++f) {
                std::int64_t corners[4];
                std::size_t ncorners = 0;
                if (FacetDim == 2) {
                    const auto& def = detail::cell_faces(type)[f];
                    ncorners = def.mNumCorners;
                    for (std::size_t k = 0; k < ncorners; ++k)
                        corners[k] = row[def.mNodes[k]];
                } else {
                    const auto& def = detail::cell_edges(type)[f];
                    ncorners = def.mNumCorners;
                    for (std::size_t k = 0; k < ncorners; ++k)
                        corners[k] = row[def.mNodes[k]];
                }
                const detail::FacetHit* hit = index.Find(corners, ncorners);
                if (hit == nullptr || hit->mCount != 1)
                    continue;
                std::uint16_t common = 0xFFFF;
                for (std::size_t k = 0; k < ncorners; ++k)
                    common &= rSoup.mMask[static_cast<std::size_t>(corners[k])];
                // The lowest bit decides, so the order the specs are listed in is
                // the order of preference between regions sharing nodes.
                std::size_t best = rSpecs.size();
                for (std::size_t q = 0; q < rSpecs.size(); ++q)
                    if ((common & rSpecs[q].mBit) &&
                        (best == rSpecs.size() || rSpecs[q].mBit < rSpecs[best].mBit))
                        best = q;
                if (best != rSpecs.size())
                    entries[best].emplace_back(cell, static_cast<std::int64_t>(f));
            }
        }
    }
    for (std::size_t r = 0; r < rSpecs.size(); ++r) {
        NDArray ids = NDArray::Uninit(DType::Int64, {entries[r].size(), std::size_t{2}});
        std::int64_t* dst = ids.As<std::int64_t>();
        for (std::size_t i = 0; i < entries[r].size(); ++i) {
            dst[2 * i] = entries[r][i].first;
            dst[2 * i + 1] = entries[r][i].second;
        }
        rMesh.AddRegion(Region(rSpecs[r].mName, RegionKind::Side, FacetDim, -1, std::move(ids)));
    }
}

Mesh prim_finish(const char* pShape, PrimSoup& rSoup, int Order,
                 const std::vector<PrimRegionSpec>& rSpecs, int FacetDim) {
    if (Order == 2)
        prim_elevate(pShape, rSoup);
    Mesh out;
    const std::size_t n = rSoup.mMask.size();
    NDArray points = NDArray::Uninit(DType::Float64, {n, std::size_t{3}});
    double* dst = points.As<double>();
    for (std::size_t i = 0; i < n; ++i)
        prim_position(rSoup, i, dst + 3 * i);
    out.AssignPoints(std::move(points));
    for (PrimBlock& blk : rSoup.mBlocks) {
        if (blk.mConn.empty())
            continue;
        NDArray conn = NDArray::Uninit(DType::Int64, {blk.mConn.size() / blk.mWidth, blk.mWidth});
        std::copy(blk.mConn.begin(), blk.mConn.end(), conn.As<std::int64_t>());
        out.AddCellBlock(cell_type_name(blk.mType), std::move(conn));
    }
    prim_add_regions(out, rSoup, rSpecs, FacetDim);
    return out;
}

// ---------------------------------------------------------------------------
// Shared lattices.
// ---------------------------------------------------------------------------

// Hexahedron `hex[8]` of lattice cell (i, j, k) of a lattice `px` x `py` nodes.
void prim_hex_nodes(std::int64_t I, std::int64_t J, std::int64_t K, std::int64_t Px,
                    std::int64_t Py, std::int64_t* pHex) {
    const std::int64_t b = (K * Py + J) * Px + I;
    const std::int64_t t = b + Px * Py;
    pHex[0] = b;
    pHex[1] = b + 1;
    pHex[2] = b + Px + 1;
    pHex[3] = b + Px;
    pHex[4] = t;
    pHex[5] = t + 1;
    pHex[6] = t + Px + 1;
    pHex[7] = t + Px;
}

void prim_push(PrimBlock& rBlock, std::initializer_list<std::int64_t> Row) {
    rBlock.mConn.insert(rBlock.mConn.end(), Row.begin(), Row.end());
}

// The squircle disk lattice: (n + 1)^2 nodes, x fastest, at layer height z.
// Returns the id of the first node it added.
std::int64_t prim_add_squircle_layer(PrimSoup& rSoup, std::int64_t N, double Z,
                                     std::uint16_t LayerMask) {
    const std::int64_t first = static_cast<std::int64_t>(rSoup.mMask.size());
    for (std::int64_t j = 0; j <= N; ++j)
        for (std::int64_t i = 0; i <= N; ++i) {
            const double s = static_cast<double>(2 * i - N) / static_cast<double>(N);
            const double t = static_cast<double>(2 * j - N) / static_cast<double>(N);
            std::uint16_t mask = LayerMask;
            if (i == 0 || i == N || j == 0 || j == N)
                mask |= 4;
            rSoup.AddNode(s, t, Z, mask);
        }
    return first;
}

// The triangle disk: a centre node and N rings of 6k nodes at radius k / N,
// every ring node at the angle 2 pi m / 6k. Returns the id of the first node.
std::int64_t prim_add_polar_layer(PrimSoup& rSoup, std::int64_t N, double Z,
                                  std::uint16_t LayerMask) {
    const std::int64_t first = static_cast<std::int64_t>(rSoup.mMask.size());
    rSoup.AddNode(0.0, 0.0, Z, static_cast<std::uint16_t>(LayerMask | kPrimDegenerate));
    for (std::int64_t k = 1; k <= N; ++k)
        for (std::int64_t m = 0; m < 6 * k; ++m) {
            const double rho = static_cast<double>(k) / static_cast<double>(N);
            const double ang = (2.0 * kPrimPi * static_cast<double>(m)) / static_cast<double>(6 * k);
            rSoup.AddNode(rho, ang, Z, static_cast<std::uint16_t>(LayerMask | (k == N ? 4 : 0)));
        }
    return first;
}

std::int64_t prim_polar_nodes(std::int64_t N) {
    return 1 + 3 * N * (N + 1);
}

// Triangles of a polar layer whose first node id is Base, in the cell order
// documented in doc/primitives.md: ring by ring, sector by sector.
template <typename Emit>
void prim_polar_triangles(std::int64_t N, std::int64_t Base, Emit&& rEmit) {
    auto ring = [&](std::int64_t k, std::int64_t m) {
        return Base + 1 + 3 * k * (k - 1) + (m % (6 * k));
    };
    for (std::int64_t k = 1; k <= N; ++k)
        for (std::int64_t s = 0; s < 6; ++s) {
            auto outer = [&](std::int64_t t) { return ring(k, s * k + t); };
            auto inner = [&](std::int64_t t) { return k == 1 ? Base : ring(k - 1, s * (k - 1) + t); };
            for (std::int64_t t = 0; t < k; ++t) {
                rEmit(inner(t), outer(t), outer(t + 1));
                if (t < k - 1)
                    rEmit(inner(t), outer(t + 1), inner(t + 1));
            }
        }
}

}  // namespace

// ===========================================================================
// rectangle
// ===========================================================================

Mesh rectangle(const std::array<std::int64_t, 2>& rDims, const std::array<double, 2>& rSize,
               const std::array<double, 3>& rOrigin, const std::string& rCellType,
               const std::string& rSplit, int Order, std::int64_t MaxCells) {
    constexpr const char* shape = "rectangle";
    for (std::size_t k = 0; k < 2; ++k) {
        prim_check_count(shape, "cell counts", rDims[k]);
        if (rDims[k] > 0)
            prim_check_positive(shape, "size", rSize[k]);
    }
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (rSplit != "zigzag" && rSplit != "up" && rSplit != "down" && rSplit != "center")
        prim_bad_choice(shape, "split", rSplit, "zigzag, up, down, center");
    const std::int64_t nx = rDims[0], ny = rDims[1];
    if (nx == 0 || ny == 0)
        return prim_empty();
    const bool center = tri && rSplit == "center";
    const std::int64_t per_cell = !tri ? 1 : (center ? 4 : 2);
    prim_check_budget(shape, prim_mul(prim_mul(nx, ny), per_cell), MaxCells);

    PrimSoup soup;
    const std::int64_t px = nx + 1, py = ny + 1;
    const double lx = rSize[0], ly = rSize[1];
    for (std::int64_t j = 0; j < py; ++j)
        for (std::int64_t i = 0; i < px; ++i) {
            std::uint16_t mask = 0;
            if (i == 0)
                mask |= 1;
            if (i == nx)
                mask |= 2;
            if (j == 0)
                mask |= 4;
            if (j == ny)
                mask |= 8;
            soup.AddNode(rOrigin[0] + (lx * static_cast<double>(i)) / static_cast<double>(nx),
                         rOrigin[1] + (ly * static_cast<double>(j)) / static_cast<double>(ny),
                         rOrigin[2], mask);
        }
    if (center)
        for (std::int64_t j = 0; j < ny; ++j)
            for (std::int64_t i = 0; i < nx; ++i)
                soup.AddNode(
                    rOrigin[0] + (lx * static_cast<double>(2 * i + 1)) / static_cast<double>(2 * nx),
                    rOrigin[1] + (ly * static_cast<double>(2 * j + 1)) / static_cast<double>(2 * ny),
                    rOrigin[2], 0);

    PrimBlock& blk = soup.Block(tri ? CellType::Triangle : CellType::Quad, tri ? 3 : 4);
    blk.mConn.reserve(static_cast<std::size_t>(nx * ny * per_cell) * blk.mWidth);
    for (std::int64_t j = 0; j < ny; ++j)
        for (std::int64_t i = 0; i < nx; ++i) {
            const std::int64_t a = j * px + i, b = a + 1, c = a + px + 1, d = a + px;
            if (!tri) {
                prim_push(blk, {a, b, c, d});
            } else if (center) {
                const std::int64_t m = px * py + j * nx + i;
                prim_push(blk, {a, b, m});
                prim_push(blk, {b, c, m});
                prim_push(blk, {c, d, m});
                prim_push(blk, {d, a, m});
            } else {
                const bool up = rSplit == "up" || (rSplit == "zigzag" && ((i + j) % 2 == 0));
                if (up) {
                    prim_push(blk, {a, b, c});
                    prim_push(blk, {a, c, d});
                } else {
                    prim_push(blk, {a, b, d});
                    prim_push(blk, {b, c, d});
                }
            }
        }
    return prim_finish(shape, soup, Order,
                       {{"xmin", 1}, {"xmax", 2}, {"ymin", 4}, {"ymax", 8}}, 1);
}

// ===========================================================================
// box
// ===========================================================================

Mesh box(const std::array<std::int64_t, 3>& rDims, const std::array<double, 3>& rSize,
         const std::array<double, 3>& rOrigin, const std::string& rCellType, int Order,
         std::int64_t MaxCells) {
    constexpr const char* shape = "box";
    for (std::size_t k = 0; k < 3; ++k) {
        prim_check_count(shape, "cell counts", rDims[k]);
        if (rDims[k] > 0)
            prim_check_positive(shape, "size", rSize[k]);
    }
    prim_check_order(shape, Order);
    if (rCellType != "hexahedron" && rCellType != "tetra" && rCellType != "wedge")
        prim_bad_choice(shape, "cell_type", rCellType, "hexahedron, tetra, wedge");
    const std::int64_t nx = rDims[0], ny = rDims[1], nz = rDims[2];
    if (nx == 0 || ny == 0 || nz == 0)
        return prim_empty();
    const std::int64_t per_cell = rCellType == "hexahedron" ? 1 : (rCellType == "tetra" ? 6 : 2);
    prim_check_budget(shape, prim_mul(prim_mul(prim_mul(nx, ny), nz), per_cell), MaxCells);

    PrimSoup soup;
    const std::int64_t px = nx + 1, py = ny + 1, pz = nz + 1;
    for (std::int64_t k = 0; k < pz; ++k)
        for (std::int64_t j = 0; j < py; ++j)
            for (std::int64_t i = 0; i < px; ++i) {
                std::uint16_t mask = 0;
                if (i == 0)
                    mask |= 1;
                if (i == nx)
                    mask |= 2;
                if (j == 0)
                    mask |= 4;
                if (j == ny)
                    mask |= 8;
                if (k == 0)
                    mask |= 16;
                if (k == nz)
                    mask |= 32;
                soup.AddNode(
                    rOrigin[0] + (rSize[0] * static_cast<double>(i)) / static_cast<double>(nx),
                    rOrigin[1] + (rSize[1] * static_cast<double>(j)) / static_cast<double>(ny),
                    rOrigin[2] + (rSize[2] * static_cast<double>(k)) / static_cast<double>(nz),
                    mask);
            }

    // The six Kuhn tetrahedra around the main diagonal 0-6, all positively
    // oriented; every one contains the edge (0, 6), so neighbours conform.
    static const std::int64_t kKuhn[6][4] = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6},
                                             {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    const CellType type = rCellType == "hexahedron" ? CellType::Hexahedron
                          : rCellType == "tetra"    ? CellType::Tetra
                                                    : CellType::Wedge;
    PrimBlock& blk = soup.Block(type, type == CellType::Hexahedron ? 8 : (type == CellType::Tetra ? 4 : 6));
    for (std::int64_t k = 0; k < nz; ++k)
        for (std::int64_t j = 0; j < ny; ++j)
            for (std::int64_t i = 0; i < nx; ++i) {
                std::int64_t h[8];
                prim_hex_nodes(i, j, k, px, py, h);
                if (type == CellType::Hexahedron) {
                    prim_push(blk, {h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]});
                } else if (type == CellType::Tetra) {
                    for (const auto& t : kKuhn)
                        prim_push(blk, {h[t[0]], h[t[1]], h[t[2]], h[t[3]]});
                } else {
                    prim_push(blk, {h[0], h[1], h[2], h[4], h[5], h[6]});
                    prim_push(blk, {h[0], h[2], h[3], h[4], h[6], h[7]});
                }
            }
    return prim_finish(shape, soup, Order,
                       {{"xmin", 1}, {"xmax", 2}, {"ymin", 4}, {"ymax", 8}, {"zmin", 16},
                        {"zmax", 32}},
                       2);
}

// ===========================================================================
// disk
// ===========================================================================

Mesh disk(double Radius, std::int64_t N, const std::array<double, 3>& rCenter,
          const std::string& rCellType, int Order, std::int64_t MaxCells) {
    constexpr const char* shape = "disk";
    prim_check_positive(shape, "radius", Radius);
    prim_check_count(shape, "n", N);
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (N == 0)
        return prim_empty();
    prim_check_budget(shape, prim_mul(prim_mul(N, N), tri ? 6 : 1), MaxCells);

    PrimSoup soup;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    if (tri) {
        soup.mMap = PrimMap::Polar;
        const std::int64_t base = prim_add_polar_layer(soup, N, 0.0, 0);
        PrimBlock& blk = soup.Block(CellType::Triangle, 3);
        prim_polar_triangles(N, base, [&](std::int64_t A, std::int64_t B, std::int64_t C) {
            prim_push(blk, {A, B, C});
        });
    } else {
        soup.mMap = PrimMap::Squircle;
        prim_add_squircle_layer(soup, N, 0.0, 0);
        PrimBlock& blk = soup.Block(CellType::Quad, 4);
        const std::int64_t p = N + 1;
        for (std::int64_t j = 0; j < N; ++j)
            for (std::int64_t i = 0; i < N; ++i) {
                const std::int64_t a = j * p + i;
                prim_push(blk, {a, a + 1, a + p + 1, a + p});
            }
    }
    return prim_finish(shape, soup, Order, {{"rim", 4}}, 1);
}

// ===========================================================================
// uv_sphere
// ===========================================================================

Mesh uv_sphere(double Radius, std::int64_t Nu, std::int64_t Nv, const std::array<double, 3>& rCenter,
               const std::string& rCellType, int Order, std::int64_t MaxCells) {
    constexpr const char* shape = "uv_sphere";
    prim_check_positive(shape, "radius", Radius);
    prim_check_count(shape, "nu", Nu);
    prim_check_count(shape, "nv", Nv);
    prim_check_order(shape, Order);
    const bool quad = rCellType == "quad";
    if (!quad && rCellType != "triangle")
        prim_bad_choice(shape, "cell_type", rCellType, "triangle, quad");
    if (Nu == 0 || Nv == 0)
        return prim_empty();
    if (Nu < 3)
        prim_fail(shape, "nu must be at least 3, got " + std::to_string(Nu));
    if (Nv < 2)
        prim_fail(shape, "nv must be at least 2, got " + std::to_string(Nv));
    prim_check_budget(shape,
                      quad ? prim_mul(Nu, Nv) : prim_mul(prim_mul(2, Nu), Nv - 1), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Sphere;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    soup.AddNode(0.0, 0.0, 0.0, kPrimDegenerate);
    for (std::int64_t j = 1; j < Nv; ++j)
        for (std::int64_t i = 0; i < Nu; ++i)
            soup.AddNode((kPrimPi * static_cast<double>(j)) / static_cast<double>(Nv),
                         (2.0 * kPrimPi * static_cast<double>(i)) / static_cast<double>(Nu), 0.0, 0);
    const std::int64_t south = soup.AddNode(kPrimPi, 0.0, 0.0, kPrimDegenerate);
    auto ring = [&](std::int64_t j, std::int64_t i) { return 1 + (j - 1) * Nu + (i % Nu); };

    if (quad) {
        PrimBlock& body = soup.Block(CellType::Quad, 4);
        for (std::int64_t j = 1; j + 1 < Nv; ++j)
            for (std::int64_t i = 0; i < Nu; ++i)
                prim_push(body, {ring(j, i), ring(j + 1, i), ring(j + 1, i + 1), ring(j, i + 1)});
    }
    PrimBlock& tris = soup.Block(CellType::Triangle, 3);
    for (std::int64_t i = 0; i < Nu; ++i)
        prim_push(tris, {0, ring(1, i), ring(1, i + 1)});
    if (!quad)
        for (std::int64_t j = 1; j + 1 < Nv; ++j)
            for (std::int64_t i = 0; i < Nu; ++i) {
                const std::int64_t a = ring(j, i), b = ring(j, i + 1);
                const std::int64_t c = ring(j + 1, i + 1), d = ring(j + 1, i);
                prim_push(tris, {a, d, c});
                prim_push(tris, {a, c, b});
            }
    for (std::int64_t i = 0; i < Nu; ++i)
        prim_push(tris, {south, ring(Nv - 1, i + 1), ring(Nv - 1, i)});
    return prim_finish(shape, soup, Order, {}, 2);
}

// ===========================================================================
// icosa_sphere
// ===========================================================================

Mesh icosa_sphere(double Radius, std::int64_t Subdivisions, const std::array<double, 3>& rCenter,
                  int Order, std::int64_t MaxCells) {
    constexpr const char* shape = "icosa_sphere";
    prim_check_positive(shape, "radius", Radius);
    prim_check_count(shape, "subdivisions", Subdivisions);
    prim_check_order(shape, Order);
    std::int64_t cells = 20;
    for (std::int64_t s = 0; s < Subdivisions; ++s)
        cells = prim_mul(cells, 4);
    prim_check_budget(shape, cells, MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Unit;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    const double t = (1.0 + std::sqrt(5.0)) / 2.0;
    const double len = std::sqrt(1.0 + t * t);
    const double one = 1.0 / len, tau = t / len;
    const double verts[12][3] = {{-one, tau, 0}, {one, tau, 0},   {-one, -tau, 0}, {one, -tau, 0},
                                 {0, -one, tau}, {0, one, tau},   {0, -one, -tau}, {0, one, -tau},
                                 {tau, 0, -one}, {tau, 0, one},   {-tau, 0, -one}, {-tau, 0, one}};
    for (const auto& v : verts)
        soup.AddNode(v[0], v[1], v[2], 0);
    static const std::int64_t kFaces[20][3] = {
        {0, 11, 5}, {0, 5, 1},  {0, 1, 7},  {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4},
        {11, 10, 2}, {10, 7, 6}, {7, 1, 8},  {3, 9, 4},  {3, 4, 2},   {3, 2, 6}, {3, 6, 8},
        {3, 8, 9},  {4, 9, 5},  {2, 4, 11}, {6, 2, 10}, {8, 6, 7},   {9, 8, 1}};
    std::vector<std::int64_t> faces;
    for (const auto& f : kFaces)
        faces.insert(faces.end(), {f[0], f[1], f[2]});
    for (std::int64_t level = 0; level < Subdivisions; ++level) {
        std::unordered_map<std::uint64_t, std::int64_t> mids;
        auto midpoint = [&](std::int64_t A, std::int64_t B) {
            if (A > B)
                std::swap(A, B);
            const auto key = prim_edge_key(A, B);
            auto found = mids.find(key);
            if (found == mids.end()) {
                double mid[3];
                prim_mid(soup, static_cast<std::size_t>(A), static_cast<std::size_t>(B), mid);
                found = mids.emplace(key, soup.AddNode(mid[0], mid[1], mid[2], 0)).first;
            }
            return found->second;
        };
        std::vector<std::int64_t> next;
        next.reserve(faces.size() * 4);
        for (std::size_t f = 0; f < faces.size(); f += 3) {
            const std::int64_t v0 = faces[f], v1 = faces[f + 1], v2 = faces[f + 2];
            const std::int64_t m01 = midpoint(v0, v1);
            const std::int64_t m12 = midpoint(v1, v2);
            const std::int64_t m20 = midpoint(v2, v0);
            next.insert(next.end(), {v0, m01, m20, v1, m12, m01, v2, m20, m12, m01, m12, m20});
        }
        faces = std::move(next);
    }
    PrimBlock& blk = soup.Block(CellType::Triangle, 3);
    blk.mConn = std::move(faces);
    return prim_finish(shape, soup, Order, {}, 2);
}

// ===========================================================================
// cylinder
// ===========================================================================

Mesh cylinder(double Radius, double Height, std::int64_t Nxy, std::int64_t Nz,
              const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
              std::int64_t MaxCells) {
    constexpr const char* shape = "cylinder";
    prim_check_positive(shape, "radius", Radius);
    prim_check_positive(shape, "height", Height);
    prim_check_count(shape, "nxy", Nxy);
    prim_check_count(shape, "nz", Nz);
    prim_check_order(shape, Order);
    const bool wedge = rCellType == "wedge";
    if (!wedge && rCellType != "hexahedron")
        prim_bad_choice(shape, "cell_type", rCellType, "hexahedron, wedge");
    if (Nxy == 0 || Nz == 0)
        return prim_empty();
    prim_check_budget(shape, prim_mul(prim_mul(prim_mul(Nxy, Nxy), Nz), wedge ? 6 : 1), MaxCells);

    PrimSoup soup;
    soup.mMap = wedge ? PrimMap::Polar : PrimMap::Squircle;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    // Bits: 1 top, 2 bottom, 4 side (the rim of the disk).
    std::int64_t layer_nodes = 0;
    for (std::int64_t k = 0; k <= Nz; ++k) {
        const double z = (Height * static_cast<double>(2 * k - Nz)) / static_cast<double>(2 * Nz);
        std::uint16_t mask = 0;
        if (k == 0)
            mask |= 2;
        if (k == Nz)
            mask |= 1;
        if (wedge) {
            prim_add_polar_layer(soup, Nxy, z, mask);
            layer_nodes = prim_polar_nodes(Nxy);
        } else {
            prim_add_squircle_layer(soup, Nxy, z, mask);
            layer_nodes = (Nxy + 1) * (Nxy + 1);
        }
    }
    PrimBlock& blk = soup.Block(wedge ? CellType::Wedge : CellType::Hexahedron, wedge ? 6 : 8);
    for (std::int64_t k = 0; k < Nz; ++k) {
        const std::int64_t base = k * layer_nodes;
        if (wedge) {
            prim_polar_triangles(Nxy, base, [&](std::int64_t A, std::int64_t B, std::int64_t C) {
                prim_push(blk, {A, B, C, A + layer_nodes, B + layer_nodes, C + layer_nodes});
            });
        } else {
            const std::int64_t p = Nxy + 1;
            for (std::int64_t j = 0; j < Nxy; ++j)
                for (std::int64_t i = 0; i < Nxy; ++i) {
                    const std::int64_t a = base + j * p + i;
                    const std::int64_t t = a + layer_nodes;
                    prim_push(blk, {a, a + 1, a + p + 1, a + p, t, t + 1, t + p + 1, t + p});
                }
        }
    }
    return prim_finish(shape, soup, Order, {{"side", 4}, {"top", 1}, {"bottom", 2}}, 2);
}

// ===========================================================================
// reference_cells
// ===========================================================================

namespace {

struct PrimRefCell {
    std::vector<std::array<double, 3>> mPoints;
    std::vector<std::vector<std::int64_t>> mCells;
};

const std::vector<std::string>& prim_all_reference_types() {
    static const std::vector<std::string> kAll = {
        "vertex", "line",  "triangle", "quad",         "tetra",    "hexahedron",
        "wedge",  "pyramid", "line3",   "triangle6",   "quad8",    "tetra10",
        "hexahedron20", "wedge15", "pyramid13", "polyhedron"};
    return kAll;
}

const std::vector<std::string>& prim_pair_reference_types() {
    static const std::vector<std::string> kPairs = {"line", "triangle", "quad", "tetra",
                                                    "hexahedron", "wedge", "pyramid"};
    return kPairs;
}

// One cell of a linear type in its reference coordinates.
PrimRefCell prim_reference_linear(const std::string& rType) {
    PrimRefCell c;
    if (rType == "vertex") {
        c.mPoints = {{0, 0, 0}};
        c.mCells = {{0}};
    } else if (rType == "line") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}};
        c.mCells = {{0, 1}};
    } else if (rType == "triangle") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
        c.mCells = {{0, 1, 2}};
    } else if (rType == "quad") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
        c.mCells = {{0, 1, 2, 3}};
    } else if (rType == "tetra") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        c.mCells = {{0, 1, 2, 3}};
    } else if (rType == "hexahedron") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                     {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
        c.mCells = {{0, 1, 2, 3, 4, 5, 6, 7}};
    } else if (rType == "wedge") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
        c.mCells = {{0, 1, 2, 3, 4, 5}};
    } else if (rType == "pyramid") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1}};
        c.mCells = {{0, 1, 2, 3, 4}};
    }
    return c;
}

// Two cells of a linear type sharing a facet.
PrimRefCell prim_reference_pair(const std::string& rType) {
    PrimRefCell c;
    if (rType == "line") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}};
        c.mCells = {{0, 1}, {1, 2}};
    } else if (rType == "triangle") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
        c.mCells = {{0, 1, 2}, {1, 3, 2}};
    } else if (rType == "quad") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}};
        c.mCells = {{0, 1, 4, 3}, {1, 2, 5, 4}};
    } else if (rType == "tetra") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}};
        c.mCells = {{0, 1, 2, 3}, {2, 3, 1, 4}};
    } else if (rType == "hexahedron") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0},
                     {0, 0, 1}, {1, 0, 1}, {2, 0, 1}, {0, 1, 1}, {1, 1, 1}, {2, 1, 1}};
        c.mCells = {{0, 1, 4, 3, 6, 7, 10, 9}, {1, 2, 5, 4, 7, 8, 11, 10}};
    } else if (rType == "wedge") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1},
                     {0, 1, 1}, {0, 0, 2}, {1, 0, 2}, {0, 1, 2}};
        c.mCells = {{0, 1, 2, 3, 4, 5}, {3, 4, 5, 6, 7, 8}};
    } else if (rType == "pyramid") {
        c.mPoints = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1}, {0.5, 0.5, -1}};
        c.mCells = {{0, 1, 2, 3, 4}, {0, 3, 2, 1, 5}};
    }
    return c;
}

// A quadratic type's cell: the linear cell plus an exact midpoint per edge.
PrimRefCell prim_reference_quadratic(CellType Linear) {
    PrimRefCell c = prim_reference_linear(cell_type_name(Linear));
    const auto& edges = detail::cell_refine_edges(Linear);
    std::vector<std::int64_t> row = c.mCells[0];
    for (const auto& e : edges) {
        const auto& a = c.mPoints[static_cast<std::size_t>(c.mCells[0][e[0]])];
        const auto& b = c.mPoints[static_cast<std::size_t>(c.mCells[0][e[1]])];
        c.mPoints.push_back({(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5, (a[2] + b[2]) * 0.5});
        row.push_back(static_cast<std::int64_t>(c.mPoints.size() - 1));
    }
    c.mCells[0] = std::move(row);
    return c;
}

}  // namespace

Mesh reference_cells(const std::vector<std::string>& rTypes, bool Pairs) {
    constexpr const char* shape = "reference_cells";
    const auto& all = prim_all_reference_types();
    std::vector<std::string> types = rTypes;
    if (types.empty())
        types = Pairs ? prim_pair_reference_types() : all;
    for (std::size_t i = 0; i < types.size(); ++i) {
        if (std::find(all.begin(), all.end(), types[i]) == all.end()) {
            std::string choices;
            for (const auto& name : all)
                choices += (choices.empty() ? "" : ", ") + name;
            prim_bad_choice(shape, "cell type", types[i], choices.c_str());
        }
        if (std::find(types.begin(), types.begin() + static_cast<std::ptrdiff_t>(i), types[i]) !=
            types.begin() + static_cast<std::ptrdiff_t>(i))
            prim_fail(shape, "cell type '" + types[i] + "' is listed twice");
        if (Pairs) {
            const auto& ok = prim_pair_reference_types();
            if (std::find(ok.begin(), ok.end(), types[i]) == ok.end())
                prim_fail(shape, "pairs=True supports line, triangle, quad, tetra, hexahedron, "
                                 "wedge and pyramid, not '" + types[i] + "'");
        }
    }

    // Cells fit in [0, 1]^3 (a pair in [0, 2] x [0, 1] x [-1, 2]), so a stride
    // of 2 (3 for pairs) keeps neighbours apart.
    const double stride = Pairs ? 3.0 : 2.0;
    std::vector<std::array<double, 3>> points;
    struct Block {
        std::string mName;
        std::size_t mWidth;
        std::vector<std::int64_t> mConn;
        bool mPoly = false;
    };
    std::vector<Block> blocks;
    for (std::size_t m = 0; m < types.size(); ++m) {
        const std::string& name = types[m];
        const double dx = stride * static_cast<double>(m);
        PrimRefCell cell;
        bool poly = false;
        if (name == "polyhedron") {
            cell = prim_reference_linear("hexahedron");
            poly = true;
        } else if (Pairs) {
            cell = prim_reference_pair(name);
        } else {
            // A quadratic name maps back to its linear base.
            static const std::pair<const char*, const char*> kBase[] = {
                {"line3", "line"},     {"triangle6", "triangle"}, {"quad8", "quad"},
                {"tetra10", "tetra"},  {"hexahedron20", "hexahedron"}, {"wedge15", "wedge"},
                {"pyramid13", "pyramid"}};
            const char* base = nullptr;
            for (const auto& kv : kBase)
                if (name == kv.first)
                    base = kv.second;
            cell = base != nullptr ? prim_reference_quadratic(cell_type_from_name(base))
                                   : prim_reference_linear(name);
        }
        const std::int64_t offset = static_cast<std::int64_t>(points.size());
        for (const auto& p : cell.mPoints)
            points.push_back({p[0] + dx, p[1], p[2]});
        Block blk;
        blk.mName = name;
        blk.mWidth = cell.mCells[0].size();
        blk.mPoly = poly;
        for (const auto& row : cell.mCells)
            for (std::int64_t v : row)
                blk.mConn.push_back(v + offset);
        blocks.push_back(std::move(blk));
    }

    Mesh out;
    NDArray pts = NDArray::Uninit(DType::Float64, {points.size(), std::size_t{3}});
    double* dst = pts.As<double>();
    for (std::size_t i = 0; i < points.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            dst[3 * i + k] = points[i][k];
    out.AssignPoints(std::move(pts));
    for (Block& blk : blocks) {
        if (blk.mPoly) {
            // The hexahedron as one polyhedron: its six outward faces.
            const auto& faces = detail::cell_faces(CellType::Hexahedron);
            std::vector<std::int64_t> flat, row_offsets{0}, face_offsets{0};
            for (const auto& f : faces) {
                for (std::size_t k = 0; k < f.mNumCorners; ++k)
                    flat.push_back(blk.mConn[f.mNodes[k]]);
                row_offsets.push_back(static_cast<std::int64_t>(flat.size()));
            }
            face_offsets.push_back(static_cast<std::int64_t>(faces.size()));
            out.AddPolyhedronBlock("polyhedron8", std::move(flat), std::move(row_offsets),
                                   std::move(face_offsets));
            continue;
        }
        NDArray conn = NDArray::Uninit(DType::Int64, {blk.mConn.size() / blk.mWidth, blk.mWidth});
        std::copy(blk.mConn.begin(), blk.mConn.end(), conn.As<std::int64_t>());
        out.AddCellBlock(blk.mName, std::move(conn));
    }
    return out;
}

}  // namespace meshioplusplus
