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
#include <map>
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

#include "../formats/lagrange_common.hpp"

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
//
// The maps after Plate came with the second tier: Project pushes a point of a
// flat lattice (a cube face, an icosahedron face) radially onto the sphere;
// Torus, Moebius, Capsule, Graded and SolidTorus are surfaces and solids of
// revolution in angles; Ball blends the butterfly's inner cube into the sphere
// and Shell a cube-surface direction times a radial parameter.
enum class PrimMap {
    Identity,
    Squircle,
    Polar,
    Sphere,
    Unit,
    Plate,
    Project,
    Torus,
    Moebius,
    Capsule,
    Ball,
    Shell,
    Graded,
    SolidTorus
};

struct PrimBlock {
    CellType mType = CellType::Vertex;
    std::size_t mWidth = 0;
    std::vector<std::int64_t> mConn;
};

struct PrimSoup {
    PrimMap mMap = PrimMap::Identity;
    double mRadius = 1.0;
    std::array<double, 3> mCenter = {{0.0, 0.0, 0.0}};
    // PrimMap::Plate: the hole radius is mRadius; the half side, the number of
    // angular positions around the hole and the cells per 90 degree block.
    double mHalf = 1.0;
    double mPeriod = 1.0;
    double mNt = 1.0;
    // Torus / SolidTorus: the minor radius. Shell / Graded: the outer radius
    // (mRadius is the inner one), the grading ratio and the layer count.
    // Moebius: mHalf is the half width, mTwists the half-twists. Capsule: mHalf
    // is half the cylinder length.
    double mRadius2 = 0.5;
    double mRatio = 1.0;
    double mLayers = 1.0;
    std::int64_t mTwists = 1;
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

constexpr int kPrimMaxOrder = 10;

void prim_check_order(const char* pShape, int Order) {
    if (Order < 1 || Order > kPrimMaxOrder)
        prim_fail(pShape, "order must be between 1 and " + std::to_string(kPrimMaxOrder) +
                              ", got " + std::to_string(Order));
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

// The radial parameter t in [0, 1] of a graded wall: uniform for a ratio of 1,
// else layer thicknesses in geometric progression (q^(n t) - 1) / (q^n - 1).
double prim_grade(const PrimSoup& rSoup, double T) {
    if (rSoup.mRatio == 1.0)
        return T;
    if (T >= 1.0)
        return 1.0;
    return (std::pow(rSoup.mRatio, rSoup.mLayers * T) - 1.0) /
           (std::pow(rSoup.mRatio, rSoup.mLayers) - 1.0);
}

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
        case PrimMap::Plate: {
            // p[0] blends the hole circle (0) into the square (1) at the angular
            // position p[1] in [0, period]; p[2] is the height.
            double q = p[1];
            if (q >= rSoup.mPeriod)
                q -= rSoup.mPeriod;
            const double nt = rSoup.mNt;
            double b = std::floor(q / nt);
            b = b < 0.0 ? 0.0 : (b > 3.0 ? 3.0 : b);
            const double m = q - b * nt;
            const double ang = -kPrimPi / 4.0 + ((kPrimPi / 2.0) * q) / nt;
            const double hx = r * std::cos(ang), hy = r * std::sin(ang);
            // The east block's outer point, turned by b quarter turns (exact).
            const double ex = rSoup.mHalf;
            const double ey = rSoup.mHalf * ((2.0 * m) / nt - 1.0);
            double ox = ex, oy = ey;
            switch (static_cast<int>(b)) {
                case 1:
                    ox = -ey;
                    oy = ex;
                    break;
                case 2:
                    ox = -ex;
                    oy = -ey;
                    break;
                case 3:
                    ox = ey;
                    oy = -ex;
                    break;
                default:
                    break;
            }
            pOut[0] = cx + ((1.0 - p[0]) * hx + p[0] * ox);
            pOut[1] = cy + ((1.0 - p[0]) * hy + p[0] * oy);
            pOut[2] = cz + p[2];
            break;
        }
        case PrimMap::Project: {
            const double len = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            pOut[0] = cx + r * (p[0] / len);
            pOut[1] = cy + r * (p[1] / len);
            pOut[2] = cz + r * (p[2] / len);
            break;
        }
        case PrimMap::Torus: {
            // p[0] the angle around the axis, p[1] the angle around the tube.
            const double rho = r + rSoup.mRadius2 * std::cos(p[1]);
            pOut[0] = cx + rho * std::cos(p[0]);
            pOut[1] = cy + rho * std::sin(p[0]);
            pOut[2] = cz + rSoup.mRadius2 * std::sin(p[1]);
            break;
        }
        case PrimMap::Moebius: {
            // p[0] the angle along the band, p[1] in [-1, 1] across it.
            const double tw = (static_cast<double>(rSoup.mTwists) * p[0]) / 2.0;
            const double off = rSoup.mHalf * p[1];
            const double rho = r + off * std::cos(tw);
            pOut[0] = cx + rho * std::cos(p[0]);
            pOut[1] = cy + rho * std::sin(p[0]);
            pOut[2] = cz + off * std::sin(tw);
            break;
        }
        case PrimMap::Capsule: {
            // p[0] the arc length along the profile over the radius, from the
            // north pole: a quarter circle, the straight side, a quarter circle.
            const double h = rSoup.mHalf;
            const double side = (2.0 * h) / r;
            if (rSoup.mMask[Node] & kPrimDegenerate) {
                pOut[0] = cx;
                pOut[1] = cy;
                pOut[2] = p[0] == 0.0 ? cz + h + r : cz - h - r;
                break;
            }
            double rho, z;
            if (p[0] <= kPrimPi / 2.0) {
                rho = r * std::sin(p[0]);
                z = h + r * std::cos(p[0]);
            } else if (p[0] <= kPrimPi / 2.0 + side) {
                rho = r;
                z = h - r * (p[0] - kPrimPi / 2.0);
            } else {
                const double t = p[0] - side;
                rho = r * std::sin(t);
                z = -h + r * std::cos(t);
            }
            pOut[0] = cx + rho * std::cos(p[1]);
            pOut[1] = cy + rho * std::sin(p[1]);
            pOut[2] = cz + z;
            break;
        }
        case PrimMap::Ball: {
            // A lattice point of [-1, 1]^3: the inner cube (max norm up to 1/2)
            // scaled by the radius, blended outside it into the sphere along the
            // ray from the centre.
            const double s = std::max(std::max(std::fabs(p[0]), std::fabs(p[1])), std::fabs(p[2]));
            if (s <= 0.5) {
                pOut[0] = cx + r * p[0];
                pOut[1] = cy + r * p[1];
                pOut[2] = cz + r * p[2];
                break;
            }
            const double d0 = p[0] / s, d1 = p[1] / s, d2 = p[2] / s;
            const double len = std::sqrt(d0 * d0 + d1 * d1 + d2 * d2);
            const double lam = 2.0 * s - 1.0;
            pOut[0] = cx + r * ((1.0 - lam) * 0.5 * d0 + lam * (d0 / len));
            pOut[1] = cy + r * ((1.0 - lam) * 0.5 * d1 + lam * (d1 / len));
            pOut[2] = cz + r * ((1.0 - lam) * 0.5 * d2 + lam * (d2 / len));
            break;
        }
        case PrimMap::Shell: {
            // A cube-surface direction times t in [1, 2]: t - 1 is the radial
            // parameter, graded between the inner and the outer radius.
            const double s = std::max(std::max(std::fabs(p[0]), std::fabs(p[1])), std::fabs(p[2]));
            const double d0 = p[0] / s, d1 = p[1] / s, d2 = p[2] / s;
            const double len = std::sqrt(d0 * d0 + d1 * d1 + d2 * d2);
            const double rad = r + (rSoup.mRadius2 - r) * prim_grade(rSoup, s - 1.0);
            pOut[0] = cx + rad * (d0 / len);
            pOut[1] = cy + rad * (d1 / len);
            pOut[2] = cz + rad * (d2 / len);
            break;
        }
        case PrimMap::Graded: {
            const double rad = r + (rSoup.mRadius2 - r) * prim_grade(rSoup, p[0]);
            pOut[0] = cx + rad * std::cos(p[1]);
            pOut[1] = cy + rad * std::sin(p[1]);
            pOut[2] = cz + p[2];
            break;
        }
        case PrimMap::SolidTorus: {
            // The squircle disk (p[0], p[1]) as the cross-section, p[2] the angle
            // around the axis.
            const double a = rSoup.mRadius2 * p[0] * std::sqrt(1.0 - p[1] * p[1] / 2.0);
            const double b = rSoup.mRadius2 * p[1] * std::sqrt(1.0 - p[0] * p[0] / 2.0);
            pOut[0] = cx + (r + a) * std::cos(p[2]);
            pOut[1] = cy + (r + a) * std::sin(p[2]);
            pOut[2] = cz + b;
            break;
        }
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

// A weighted combination of K nodes in parameter space: the face and body
// centres of the full-Lagrange types, the lattice nodes of the VTK Lagrange
// family and the mid-edge nodes of the maps after Plate. Lattice coordinates
// average; an angle averages along the short way from the first node with a
// meaningful one (a pole or an apex has none); a unit vector is the normalised
// sum; and a Moebius node seen across the seam flips across the band.
void prim_combine(const PrimSoup& rSoup, const std::int64_t* pNodes, const double* pW,
                  std::size_t K, double* pOut) {
    const double two_pi = 2.0 * kPrimPi;
    auto param = [&](std::size_t I) {
        return rSoup.mParam.data() + 3 * static_cast<std::size_t>(pNodes[I]);
    };
    auto degenerate = [&](std::size_t I) {
        return (rSoup.mMask[static_cast<std::size_t>(pNodes[I])] & kPrimDegenerate) != 0;
    };
    double wsum = 0.0;
    for (std::size_t i = 0; i < K; ++i)
        wsum += pW[i];
    auto linear = [&](std::size_t C) {
        double acc = 0.0;
        for (std::size_t i = 0; i < K; ++i)
            acc += pW[i] * param(i)[C];
        return acc / wsum;
    };
    auto angular = [&](std::size_t C, double Period) {
        std::size_t ref = 0;
        while (ref < K && degenerate(ref))
            ++ref;
        if (ref == K)
            return param(0)[C];
        const double base = param(ref)[C];
        double acc = 0.0, ws = 0.0;
        for (std::size_t i = 0; i < K; ++i) {
            if (degenerate(i))
                continue;
            double d = param(i)[C] - base;
            d -= Period * std::floor(d / Period + 0.5);
            acc += pW[i] * (base + d);
            ws += pW[i];
        }
        return acc / ws;
    };
    switch (rSoup.mMap) {
        case PrimMap::Unit: {
            const double x = linear(0), y = linear(1), z = linear(2);
            const double len = std::sqrt(x * x + y * y + z * z);
            pOut[0] = x / len;
            pOut[1] = y / len;
            pOut[2] = z / len;
            return;
        }
        case PrimMap::Polar:
        case PrimMap::Graded:
        case PrimMap::Sphere:
        case PrimMap::Capsule:
            pOut[0] = linear(0);
            pOut[1] = angular(1, two_pi);
            pOut[2] = linear(2);
            return;
        case PrimMap::Torus:
            pOut[0] = angular(0, two_pi);
            pOut[1] = angular(1, two_pi);
            pOut[2] = linear(2);
            return;
        case PrimMap::SolidTorus:
            pOut[0] = linear(0);
            pOut[1] = linear(1);
            pOut[2] = angular(2, two_pi);
            return;
        case PrimMap::Plate: {
            pOut[0] = linear(0);
            double q = angular(1, rSoup.mPeriod);
            if (q < 0.0)
                q += rSoup.mPeriod;
            pOut[1] = q;
            pOut[2] = linear(2);
            return;
        }
        case PrimMap::Moebius: {
            const double base = param(0)[0];
            double au = 0.0, aw = 0.0;
            for (std::size_t i = 0; i < K; ++i) {
                double d = param(i)[0] - base;
                const double turns = std::floor(d / two_pi + 0.5);
                d -= two_pi * turns;
                const bool flip =
                    (static_cast<std::int64_t>(std::fabs(turns)) * rSoup.mTwists) % 2 != 0;
                au += pW[i] * (base + d);
                aw += pW[i] * (flip ? -param(i)[1] : param(i)[1]);
            }
            pOut[0] = au / wsum;
            pOut[1] = aw / wsum;
            pOut[2] = linear(2);
            return;
        }
        case PrimMap::Identity:
        case PrimMap::Squircle:
        case PrimMap::Project:
        case PrimMap::Ball:
        case PrimMap::Shell:
            break;
    }
    pOut[0] = linear(0);
    pOut[1] = linear(1);
    pOut[2] = linear(2);
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
        case PrimMap::Plate: {
            // The angular position is periodic: average along the shorter way.
            const double per = rSoup.mPeriod;
            double d = b[1] - a[1];
            if (d > 0.5 * per)
                d -= per;
            else if (d < -0.5 * per)
                d += per;
            double q = a[1] + 0.5 * d;
            if (q < 0.0)
                q += per;
            pOut[0] = (a[0] + b[0]) * 0.5;
            pOut[1] = q;
            pOut[2] = (a[2] + b[2]) * 0.5;
            break;
        }
        case PrimMap::Project:
        case PrimMap::Torus:
        case PrimMap::Moebius:
        case PrimMap::Capsule:
        case PrimMap::Ball:
        case PrimMap::Shell:
        case PrimMap::Graded:
        case PrimMap::SolidTorus: {
            const std::int64_t nodes[2] = {static_cast<std::int64_t>(A),
                                           static_cast<std::int64_t>(B)};
            const double w[2] = {1.0, 1.0};
            prim_combine(rSoup, nodes, w, 2, pOut);
            break;
        }
    }
}

std::uint64_t prim_edge_key(std::int64_t A, std::int64_t B) {
    return (static_cast<std::uint64_t>(A) << 32) | static_cast<std::uint64_t>(B);
}

CellType prim_elevated(CellType Type) {
    switch (Type) {
        case CellType::Line:
            return CellType::Line3;
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

// The centre node a full-Lagrange type adds over its serendipity twin, and the
// corners it averages, by local node index: quad9's node 8, hexahedron27's face
// centres 20-25 (cell_faces' own order) and body centre 26, wedge18's three
// quad-face centres 15-17.
std::vector<std::pair<std::size_t, std::vector<std::size_t>>> prim_centres(CellType Full) {
    std::vector<std::pair<std::size_t, std::vector<std::size_t>>> out;
    if (Full == CellType::Quad9) {
        out.push_back({8, {0, 1, 2, 3}});
        return out;
    }
    for (const auto& face : detail::cell_faces(Full))
        if (face.mNumNodes == 9)
            out.push_back(
                {face.mNodes[8], {face.mNodes[0], face.mNodes[1], face.mNodes[2], face.mNodes[3]}});
    if (Full == CellType::Hexahedron27)
        out.push_back({26, {0, 1, 2, 3, 4, 5, 6, 7}});
    std::sort(out.begin(), out.end());
    return out;
}

// Order 2 with Complete: the serendipity cells gain their face and body
// centres, appended after every mid-edge node in first-appearance order over
// the cells, each cell's in local node order, and shared by every cell on the
// face.
void prim_complete(PrimSoup& rSoup) {
    std::map<std::vector<std::int64_t>, std::int64_t> ids;
    for (PrimBlock& blk : rSoup.mBlocks) {
        CellType full;
        std::size_t width;
        switch (blk.mType) {
            case CellType::Quad8:
                full = CellType::Quad9;
                width = 9;
                break;
            case CellType::Hexahedron20:
                full = CellType::Hexahedron27;
                width = 27;
                break;
            case CellType::Wedge15:
                full = CellType::Wedge18;
                width = 18;
                break;
            default:
                continue;
        }
        const auto centres = prim_centres(full);
        const std::size_t ncells = blk.mConn.size() / blk.mWidth;
        std::vector<std::int64_t> conn(ncells * width);
        for (std::size_t c = 0; c < ncells; ++c) {
            const std::int64_t* row = blk.mConn.data() + c * blk.mWidth;
            std::int64_t* dst = conn.data() + c * width;
            std::copy(row, row + blk.mWidth, dst);
            for (const auto& [local, corners] : centres) {
                std::vector<std::int64_t> nodes;
                for (std::size_t k : corners)
                    nodes.push_back(row[k]);
                std::vector<std::int64_t> key = nodes;
                std::sort(key.begin(), key.end());
                auto found = ids.find(key);
                if (found == ids.end()) {
                    const std::vector<double> w(nodes.size(), 1.0);
                    double centre[3];
                    prim_combine(rSoup, nodes.data(), w.data(), nodes.size(), centre);
                    found =
                        ids.emplace(key, rSoup.AddNode(centre[0], centre[1], centre[2], 0)).first;
                }
                dst[local] = found->second;
            }
        }
        blk.mConn = std::move(conn);
        blk.mWidth = width;
        blk.mType = full;
    }
}

// Order 3 and up: every linear cell becomes its VTK Lagrange twin. A lattice
// node is named by the corners it combines (lagrange::lattice_weights), so two
// cells on one edge or face share it; it is appended on first appearance over
// the cells, each cell's in VTK order, its parameters the same weighted
// combination of the corners'.
void prim_lagrange(const char* pShape, PrimSoup& rSoup, int Order) {
    std::map<std::vector<std::int64_t>, std::int64_t> ids;
    for (PrimBlock& blk : rSoup.mBlocks) {
        lagrange::Shape shape;
        CellType out;
        switch (blk.mType) {
            case CellType::Line:
                shape = lagrange::Shape::Line;
                out = CellType::VtkLagrangeCurve;
                break;
            case CellType::Triangle:
                shape = lagrange::Shape::Triangle;
                out = CellType::VtkLagrangeTriangle;
                break;
            case CellType::Quad:
                shape = lagrange::Shape::Quad;
                out = CellType::VtkLagrangeQuadrilateral;
                break;
            case CellType::Tetra:
                shape = lagrange::Shape::Tetra;
                out = CellType::VtkLagrangeTetrahedron;
                break;
            case CellType::Hexahedron:
                shape = lagrange::Shape::Hexahedron;
                out = CellType::VtkLagrangeHexahedron;
                break;
            case CellType::Wedge:
                shape = lagrange::Shape::Wedge;
                out = CellType::VtkLagrangeWedge;
                break;
            default:
                prim_fail(pShape, std::string("order ") + std::to_string(Order) +
                                      " is not available for " + cell_type_name(blk.mType) +
                                      " cells");
        }
        const auto lattice = lagrange::vtk_lattice(shape, Order);
        std::vector<std::vector<std::pair<int, std::int64_t>>> weights;
        weights.reserve(lattice.size());
        for (const auto& ijk : lattice)
            weights.push_back(lagrange::lattice_weights(shape, Order, ijk));
        const std::size_t width = lattice.size();
        const std::size_t ncells = blk.mConn.size() / blk.mWidth;
        std::vector<std::int64_t> conn(ncells * width);
        std::vector<std::int64_t> nodes;
        std::vector<double> w;
        for (std::size_t c = 0; c < ncells; ++c) {
            const std::int64_t* row = blk.mConn.data() + c * blk.mWidth;
            std::int64_t* dst = conn.data() + c * width;
            for (std::size_t l = 0; l < width; ++l) {
                const auto& wl = weights[l];
                if (wl.size() == 1) {
                    dst[l] = row[wl[0].first];
                    continue;
                }
                std::vector<std::pair<std::int64_t, std::int64_t>> pairs;
                for (const auto& [corner, weight] : wl)
                    pairs.emplace_back(row[corner], weight);
                std::sort(pairs.begin(), pairs.end());
                std::vector<std::int64_t> key;
                for (const auto& [node, weight] : pairs) {
                    key.push_back(node);
                    key.push_back(weight);
                }
                auto found = ids.find(key);
                if (found == ids.end()) {
                    nodes.clear();
                    w.clear();
                    for (const auto& [corner, weight] : wl) {
                        nodes.push_back(row[corner]);
                        w.push_back(static_cast<double>(weight));
                    }
                    double param[3];
                    prim_combine(rSoup, nodes.data(), w.data(), nodes.size(), param);
                    found =
                        ids.emplace(std::move(key), rSoup.AddNode(param[0], param[1], param[2], 0))
                            .first;
                }
                dst[l] = found->second;
            }
        }
        blk.mConn = std::move(conn);
        blk.mWidth = width;
        blk.mType = out;
    }
}

// ---------------------------------------------------------------------------
// Boundary regions: a boundary facet joins the first region whose bit every one
// of its corners carries.
// ---------------------------------------------------------------------------

std::vector<Region> prim_side_regions(const Mesh& rMesh, const PrimSoup& rSoup,
                                      const std::vector<PrimRegionSpec>& rSpecs, int FacetDim) {
    std::vector<Region> out;
    if (rSpecs.empty())
        return out;
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
        out.emplace_back(rSpecs[r].mName, RegionKind::Side, FacetDim, -1, std::move(ids));
    }
    return out;
}

Mesh prim_build(const PrimSoup& rSoup) {
    Mesh out;
    const std::size_t n = rSoup.mMask.size();
    NDArray points = NDArray::Uninit(DType::Float64, {n, std::size_t{3}});
    double* dst = points.As<double>();
    for (std::size_t i = 0; i < n; ++i)
        prim_position(rSoup, i, dst + 3 * i);
    out.AssignPoints(std::move(points));
    for (const PrimBlock& blk : rSoup.mBlocks) {
        if (blk.mConn.empty())
            continue;
        NDArray conn = NDArray::Uninit(DType::Int64, {blk.mConn.size() / blk.mWidth, blk.mWidth});
        std::copy(blk.mConn.begin(), blk.mConn.end(), conn.As<std::int64_t>());
        out.AddCellBlock(cell_type_name(blk.mType), std::move(conn));
    }
    return out;
}

// Order 2 elevates (and with Complete completes) before the side regions are
// found on the result; order 3 and up finds them on the linear cells first,
// since the VTK Lagrange types have no facet table: a region names a facet by
// the linear type's local number, which is the Lagrange cell's too.
Mesh prim_finish(const char* pShape, PrimSoup& rSoup, int Order, bool Complete,
                 const std::vector<PrimRegionSpec>& rSpecs, int FacetDim) {
    std::vector<Region> regions;
    if (Order >= 3 && !rSpecs.empty())
        regions = prim_side_regions(prim_build(rSoup), rSoup, rSpecs, FacetDim);
    if (Order == 2) {
        prim_elevate(pShape, rSoup);
        if (Complete)
            prim_complete(rSoup);
    } else if (Order >= 3) {
        prim_lagrange(pShape, rSoup, Order);
    }
    Mesh out = prim_build(rSoup);
    if (Order < 3)
        regions = prim_side_regions(out, rSoup, rSpecs, FacetDim);
    for (Region& region : regions)
        out.AddRegion(std::move(region));
    return out;
}

// Cell regions from one label per cell (block-major, in soup order; empty
// blocks add no cells), label k naming rNames[k].
void prim_add_cell_regions(Mesh& rMesh, const std::vector<std::uint8_t>& rLabels,
                           const std::vector<const char*>& rNames, int Dim) {
    for (std::size_t k = 0; k < rNames.size(); ++k) {
        std::vector<std::int64_t> ids;
        for (std::size_t c = 0; c < rLabels.size(); ++c)
            if (rLabels[c] == k)
                ids.push_back(static_cast<std::int64_t>(c));
        NDArray arr = NDArray::Uninit(DType::Int64, {ids.size()});
        std::copy(ids.begin(), ids.end(), arr.As<std::int64_t>());
        rMesh.AddRegion(Region(rNames[k], RegionKind::Cell, Dim, -1, std::move(arr)));
    }
}

void prim_add_point_region(Mesh& rMesh, const char* pName, std::int64_t Point) {
    NDArray arr = NDArray::Uninit(DType::Int64, {std::size_t{1}});
    arr.As<std::int64_t>()[0] = Point;
    rMesh.AddRegion(Region(pName, RegionKind::Point, 0, -1, std::move(arr)));
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

// A quad as itself, or as the triangles (a b c) and (a c d).
void prim_push_quad(PrimBlock& rBlock, bool Tri, std::int64_t A, std::int64_t B, std::int64_t C,
                    std::int64_t D) {
    if (Tri) {
        prim_push(rBlock, {A, B, C});
        prim_push(rBlock, {A, C, D});
    } else {
        prim_push(rBlock, {A, B, C, D});
    }
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
// documented in doc/primitives.md: ring by ring, sector by sector. Ring k holds
// Sectors * k nodes (6 for the disk, the polygon's side count for an ngon).
template <typename Emit>
void prim_polar_triangles(std::int64_t N, std::int64_t Base, Emit&& rEmit,
                          std::int64_t Sectors = 6) {
    auto ring = [&](std::int64_t k, std::int64_t m) {
        return Base + 1 + (Sectors * k * (k - 1)) / 2 + (m % (Sectors * k));
    };
    for (std::int64_t k = 1; k <= N; ++k)
        for (std::int64_t s = 0; s < Sectors; ++s) {
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
               const std::string& rSplit, int Order, std::int64_t MaxCells, bool Complete) {
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
    return prim_finish(shape, soup, Order, Complete,
                       {{"xmin", 1}, {"xmax", 2}, {"ymin", 4}, {"ymax", 8}}, 1);
}

// ===========================================================================
// box
// ===========================================================================

Mesh box(const std::array<std::int64_t, 3>& rDims, const std::array<double, 3>& rSize,
         const std::array<double, 3>& rOrigin, const std::string& rCellType, int Order,
         std::int64_t MaxCells, bool Complete) {
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
    return prim_finish(
        shape, soup, Order, Complete,
        {{"xmin", 1}, {"xmax", 2}, {"ymin", 4}, {"ymax", 8}, {"zmin", 16}, {"zmax", 32}}, 2);
}

// ===========================================================================
// disk
// ===========================================================================

Mesh disk(double Radius, std::int64_t N, const std::array<double, 3>& rCenter,
          const std::string& rCellType, int Order, std::int64_t MaxCells, bool Complete) {
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
    return prim_finish(shape, soup, Order, Complete, {{"rim", 4}}, 1);
}

// ===========================================================================
// uv_sphere
// ===========================================================================

Mesh uv_sphere(double Radius, std::int64_t Nu, std::int64_t Nv,
               const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
               std::int64_t MaxCells, bool Complete) {
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
    return prim_finish(shape, soup, Order, Complete, {}, 2);
}

// ===========================================================================
// icosa_sphere
// ===========================================================================

namespace {
void prim_subdivide_sphere(PrimSoup& rSoup, std::vector<std::int64_t>& rFaces,
                           std::int64_t Subdivisions);
}  // namespace

Mesh icosa_sphere(double Radius, std::int64_t Subdivisions, const std::array<double, 3>& rCenter,
                  int Order, std::int64_t MaxCells, bool Complete) {
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
    prim_subdivide_sphere(soup, faces, Subdivisions);
    PrimBlock& blk = soup.Block(CellType::Triangle, 3);
    blk.mConn = std::move(faces);
    return prim_finish(shape, soup, Order, Complete, {}, 2);
}

namespace {

// Each level splits every triangle in four through its edge midpoints, pushed
// onto the sphere (PrimMap::Unit).
void prim_subdivide_sphere(PrimSoup& rSoup, std::vector<std::int64_t>& rFaces,
                           std::int64_t Subdivisions) {
    PrimSoup& soup = rSoup;
    std::vector<std::int64_t>& faces = rFaces;
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
}

}  // namespace

// ===========================================================================
// cylinder
// ===========================================================================

Mesh cylinder(double Radius, double Height, std::int64_t Nxy, std::int64_t Nz,
              const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
              std::int64_t MaxCells, bool Complete) {
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
    return prim_finish(shape, soup, Order, Complete, {{"side", 4}, {"top", 1}, {"bottom", 2}}, 2);
}

// ===========================================================================
// triangle
// ===========================================================================

Mesh reference_triangle(std::int64_t N, const std::array<double, 2>& rSize,
                        const std::array<double, 3>& rOrigin, int Order, std::int64_t MaxCells,
                        bool Complete) {
    constexpr const char* shape = "reference_triangle";
    prim_check_count(shape, "n", N);
    prim_check_order(shape, Order);
    if (N == 0)
        return prim_empty();
    for (std::size_t k = 0; k < 2; ++k)
        prim_check_positive(shape, "size", rSize[k]);
    prim_check_budget(shape, prim_mul(N, N), MaxCells);

    PrimSoup soup;
    // Row j holds N + 1 - j nodes; id(i, j) = j (N + 1) - j (j - 1) / 2 + i.
    auto id = [&](std::int64_t I, std::int64_t J) { return J * (N + 1) - (J * (J - 1)) / 2 + I; };
    for (std::int64_t j = 0; j <= N; ++j)
        for (std::int64_t i = 0; i <= N - j; ++i) {
            std::uint16_t mask = 0;
            if (j == 0)
                mask |= 1;
            if (i == 0)
                mask |= 2;
            if (i + j == N)
                mask |= 4;
            soup.AddNode(rOrigin[0] + (rSize[0] * static_cast<double>(i)) / static_cast<double>(N),
                         rOrigin[1] + (rSize[1] * static_cast<double>(j)) / static_cast<double>(N),
                         rOrigin[2], mask);
        }
    PrimBlock& blk = soup.Block(CellType::Triangle, 3);
    blk.mConn.reserve(static_cast<std::size_t>(N * N) * 3);
    for (std::int64_t j = 0; j < N; ++j)
        for (std::int64_t i = 0; i < N - j; ++i) {
            prim_push(blk, {id(i, j), id(i + 1, j), id(i, j + 1)});
            if (i < N - 1 - j)
                prim_push(blk, {id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)});
        }
    return prim_finish(shape, soup, Order, Complete, {{"ymin", 1}, {"xmin", 2}, {"hypotenuse", 4}},
                       1);
}

// ===========================================================================
// ngon
// ===========================================================================

Mesh ngon(std::int64_t Sides, double Radius, std::int64_t Rings,
          const std::array<double, 3>& rCenter, int Order, std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "ngon";
    prim_check_count(shape, "sides", Sides);
    prim_check_count(shape, "rings", Rings);
    prim_check_order(shape, Order);
    if (Sides == 0 || Rings == 0)
        return prim_empty();
    if (Sides < 3)
        prim_fail(shape, "sides must be at least 3, got " + std::to_string(Sides));
    prim_check_positive(shape, "radius", Radius);
    prim_check_budget(shape, prim_mul(Sides, prim_mul(Rings, Rings)), MaxCells);

    PrimSoup soup;
    std::vector<double> vx(static_cast<std::size_t>(Sides)), vy(static_cast<std::size_t>(Sides));
    for (std::int64_t s = 0; s < Sides; ++s) {
        const double ang = (2.0 * kPrimPi * static_cast<double>(s)) / static_cast<double>(Sides);
        vx[static_cast<std::size_t>(s)] = Radius * std::cos(ang);
        vy[static_cast<std::size_t>(s)] = Radius * std::sin(ang);
    }
    const double rings = static_cast<double>(Rings);
    soup.AddNode(rCenter[0], rCenter[1], rCenter[2], 0);
    for (std::int64_t k = 1; k <= Rings; ++k)
        for (std::int64_t s = 0; s < Sides; ++s)
            for (std::int64_t t = 0; t < k; ++t) {
                const std::size_t a = static_cast<std::size_t>(s);
                const std::size_t b = static_cast<std::size_t>((s + 1) % Sides);
                const double w0 = static_cast<double>(k - t), w1 = static_cast<double>(t);
                soup.AddNode(rCenter[0] + (w0 * vx[a] + w1 * vx[b]) / rings,
                             rCenter[1] + (w0 * vy[a] + w1 * vy[b]) / rings, rCenter[2],
                             k == Rings ? 1 : 0);
            }
    PrimBlock& blk = soup.Block(CellType::Triangle, 3);
    prim_polar_triangles(
        Rings, 0,
        [&](std::int64_t A, std::int64_t B, std::int64_t C) { prim_push(blk, {A, B, C}); }, Sides);
    return prim_finish(shape, soup, Order, Complete, {{"rim", 1}}, 1);
}

// ===========================================================================
// annulus
// ===========================================================================

Mesh annulus(double InnerRadius, double OuterRadius, std::int64_t Ntheta, std::int64_t Nr,
             double Angle, const std::array<double, 3>& rCenter, const std::string& rCellType,
             int Order, std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "annulus";
    prim_check_count(shape, "ntheta", Ntheta);
    prim_check_count(shape, "nr", Nr);
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (Ntheta == 0 || Nr == 0)
        return prim_empty();
    prim_check_positive(shape, "inner_radius", InnerRadius);
    if (!(OuterRadius > InnerRadius))
        prim_fail(shape, "outer_radius must exceed inner_radius, got " +
                             std::to_string(OuterRadius) + " and " + std::to_string(InnerRadius));
    const double two_pi = 2.0 * kPrimPi;
    if (!(Angle > 0.0) || Angle > two_pi + 1e-12)
        prim_fail(shape, "angle must be in (0, 2*pi], got " + std::to_string(Angle));
    const bool full = Angle >= two_pi - 1e-12;
    if (full && Ntheta < 3)
        prim_fail(shape,
                  "ntheta must be at least 3 for a full ring, got " + std::to_string(Ntheta));
    if (!full && Angle / static_cast<double>(Ntheta) >= kPrimPi)
        prim_fail(shape, "ntheta is too small: each cell must span less than pi, got angle " +
                             std::to_string(Angle) + " over " + std::to_string(Ntheta) + " cells");
    prim_check_budget(shape, prim_mul(prim_mul(Ntheta, Nr), tri ? 2 : 1), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Polar;
    soup.mRadius = 1.0;
    soup.mCenter = rCenter;
    const std::int64_t np = full ? Ntheta : Ntheta + 1;
    for (std::int64_t k = 0; k <= Nr; ++k) {
        const double rho =
            k == Nr ? OuterRadius
                    : InnerRadius + ((OuterRadius - InnerRadius) * static_cast<double>(k)) /
                                        static_cast<double>(Nr);
        for (std::int64_t i = 0; i < np; ++i) {
            std::uint16_t mask = 0;
            if (k == 0)
                mask |= 1;
            if (k == Nr)
                mask |= 2;
            if (!full && i == 0)
                mask |= 4;
            if (!full && i == Ntheta)
                mask |= 8;
            soup.AddNode(rho, (Angle * static_cast<double>(i)) / static_cast<double>(Ntheta), 0.0,
                         mask);
        }
    }
    PrimBlock& blk = soup.Block(tri ? CellType::Triangle : CellType::Quad, tri ? 3 : 4);
    for (std::int64_t k = 0; k < Nr; ++k)
        for (std::int64_t i = 0; i < Ntheta; ++i) {
            const std::int64_t i1 = (i + 1) % np;
            prim_push_quad(blk, tri, k * np + i, (k + 1) * np + i, (k + 1) * np + i1, k * np + i1);
        }
    std::vector<PrimRegionSpec> specs = {{"inner", 1}, {"outer", 2}};
    if (!full) {
        specs.push_back({"start", 4});
        specs.push_back({"end", 8});
    }
    return prim_finish(shape, soup, Order, Complete, specs, 1);
}

// ===========================================================================
// lshape
// ===========================================================================

Mesh lshape(std::int64_t N, double Size, const std::array<double, 3>& rCenter,
            const std::string& rCellType, int Order, std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "lshape";
    prim_check_count(shape, "n", N);
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (N == 0)
        return prim_empty();
    prim_check_positive(shape, "size", Size);
    prim_check_budget(shape, prim_mul(prim_mul(3, prim_mul(N, N)), tri ? 2 : 1), MaxCells);

    PrimSoup soup;
    const std::int64_t w = 2 * N + 1;
    // The corner of the square in the removed quadrant (x > 0, y < 0) has no nodes
    // strictly inside it; the lattice keeps its boundary.
    std::vector<std::int64_t> ids(static_cast<std::size_t>(w * w), -1);
    for (std::int64_t j = 0; j < w; ++j)
        for (std::int64_t i = 0; i < w; ++i) {
            if (i > N && j < N)
                continue;
            std::uint16_t mask = 0;
            if (i == 0)
                mask |= 1;
            if (i == 2 * N && j >= N)
                mask |= 2;
            if (j == 0 && i <= N)
                mask |= 4;
            if (j == 2 * N)
                mask |= 8;
            if (i == N && j <= N)
                mask |= 16;
            if (j == N && i >= N)
                mask |= 32;
            ids[static_cast<std::size_t>(j * w + i)] = soup.AddNode(
                rCenter[0] + (Size * static_cast<double>(i - N)) / static_cast<double>(2 * N),
                rCenter[1] + (Size * static_cast<double>(j - N)) / static_cast<double>(2 * N),
                rCenter[2], mask);
        }
    PrimBlock& blk = soup.Block(tri ? CellType::Triangle : CellType::Quad, tri ? 3 : 4);
    auto at = [&](std::int64_t I, std::int64_t J) {
        return ids[static_cast<std::size_t>(J * w + I)];
    };
    for (std::int64_t j = 0; j < 2 * N; ++j)
        for (std::int64_t i = 0; i < 2 * N; ++i) {
            if (i >= N && j < N)
                continue;
            prim_push_quad(blk, tri, at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
        }
    return prim_finish(
        shape, soup, Order, Complete,
        {{"xmin", 1}, {"xmax", 2}, {"ymin", 4}, {"ymax", 8}, {"notch_x", 16}, {"notch_y", 32}}, 1);
}

// ===========================================================================
// plate_with_hole
// ===========================================================================

Mesh plate_with_hole(double Radius, double Size, std::int64_t Nt, std::int64_t Nr,
                     const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
                     std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "plate_with_hole";
    prim_check_count(shape, "nt", Nt);
    prim_check_count(shape, "nr", Nr);
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (Nt == 0 || Nr == 0)
        return prim_empty();
    prim_check_positive(shape, "radius", Radius);
    prim_check_positive(shape, "size", Size);
    if (!(Radius < Size / 2.0))
        prim_fail(shape, "radius must be less than half the size, got radius " +
                             std::to_string(Radius) + " and size " + std::to_string(Size));
    const std::int64_t m = prim_mul(4, Nt);
    prim_check_budget(shape, prim_mul(prim_mul(m, Nr), tri ? 2 : 1), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Plate;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    soup.mHalf = Size / 2.0;
    soup.mPeriod = static_cast<double>(m);
    soup.mNt = static_cast<double>(Nt);
    // The block a position belongs to decides its outer side: east, north, west,
    // south. A block corner carries the side before it too.
    const std::uint16_t side_bit[4] = {4, 16, 2, 8};  // xmax, ymax, xmin, ymin
    for (std::int64_t k = 0; k <= Nr; ++k)
        for (std::int64_t p = 0; p < m; ++p) {
            std::uint16_t mask = 0;
            if (k == 0)
                mask |= 1;
            if (k == Nr) {
                const std::int64_t b = p / Nt;
                mask |= side_bit[b];
                if (p % Nt == 0)
                    mask |= side_bit[(b + 3) % 4];
            }
            soup.AddNode(static_cast<double>(k) / static_cast<double>(Nr), static_cast<double>(p),
                         0.0, mask);
        }
    PrimBlock& blk = soup.Block(tri ? CellType::Triangle : CellType::Quad, tri ? 3 : 4);
    for (std::int64_t k = 0; k < Nr; ++k)
        for (std::int64_t p = 0; p < m; ++p) {
            const std::int64_t p1 = (p + 1) % m;
            prim_push_quad(blk, tri, k * m + p, (k + 1) * m + p, (k + 1) * m + p1, k * m + p1);
        }
    return prim_finish(shape, soup, Order, Complete,
                       {{"hole", 1}, {"xmin", 2}, {"xmax", 4}, {"ymin", 8}, {"ymax", 16}}, 1);
}

// ===========================================================================
// honeycomb
// ===========================================================================

Mesh honeycomb(const std::array<std::int64_t, 2>& rDims, double Edge,
               const std::array<double, 3>& rOrigin, std::int64_t MaxCells) {
    constexpr const char* shape = "honeycomb";
    prim_check_count(shape, "cell counts", rDims[0]);
    prim_check_count(shape, "cell counts", rDims[1]);
    const std::int64_t nx = rDims[0], ny = rDims[1];
    if (nx == 0 || ny == 0)
        return prim_empty();
    prim_check_positive(shape, "edge", Edge);
    prim_check_budget(shape, prim_mul(nx, ny), MaxCells);

    // Pointy-top hexagons on the vertex lattice x = (sqrt(3) / 2) edge u,
    // y = (edge / 2) v: the centre of hexagon (i, j) is (2 i + 1 + (j & 1), 3 j + 2).
    const std::int64_t wu = 2 * nx + 2, wv = 3 * ny + 2;
    std::vector<std::int64_t> ids(static_cast<std::size_t>(wu * wv), -1);
    std::vector<double> coords;
    std::vector<std::int64_t> flat;
    flat.reserve(static_cast<std::size_t>(nx * ny) * 6);
    std::vector<std::int64_t> offsets{0};
    const double dx = (std::sqrt(3.0) * Edge) / 2.0, dy = Edge / 2.0;
    static const std::int64_t kCorner[6][2] = {{0, -2}, {1, -1}, {1, 1}, {0, 2}, {-1, 1}, {-1, -1}};
    for (std::int64_t j = 0; j < ny; ++j)
        for (std::int64_t i = 0; i < nx; ++i) {
            const std::int64_t cu = 2 * i + 1 + (j & 1), cv = 3 * j + 2;
            for (const auto& corner : kCorner) {
                const std::int64_t u = cu + corner[0], v = cv + corner[1];
                std::int64_t& slot = ids[static_cast<std::size_t>(v * wu + u)];
                if (slot < 0) {
                    slot = static_cast<std::int64_t>(coords.size() / 3);
                    coords.push_back(rOrigin[0] + dx * static_cast<double>(u));
                    coords.push_back(rOrigin[1] + dy * static_cast<double>(v));
                    coords.push_back(rOrigin[2]);
                }
                flat.push_back(slot);
            }
            offsets.push_back(static_cast<std::int64_t>(flat.size()));
        }

    Mesh out;
    NDArray points = NDArray::Uninit(DType::Float64, {coords.size() / 3, std::size_t{3}});
    std::copy(coords.begin(), coords.end(), points.As<double>());
    out.AssignPoints(std::move(points));

    // Facet f of a polygon is its edge from node f to node f + 1; an edge seen
    // once is on the boundary.
    std::unordered_map<std::uint64_t, int> seen;
    const std::size_t ncells = static_cast<std::size_t>(nx * ny);
    auto key = [&](std::size_t C, std::size_t F) {
        std::int64_t a = flat[C * 6 + F], b = flat[C * 6 + (F + 1) % 6];
        if (a > b)
            std::swap(a, b);
        return prim_edge_key(a, b);
    };
    for (std::size_t c = 0; c < ncells; ++c)
        for (std::size_t f = 0; f < 6; ++f)
            ++seen[key(c, f)];
    std::vector<std::int64_t> boundary;
    for (std::size_t c = 0; c < ncells; ++c)
        for (std::size_t f = 0; f < 6; ++f)
            if (seen[key(c, f)] == 1) {
                boundary.push_back(static_cast<std::int64_t>(c));
                boundary.push_back(static_cast<std::int64_t>(f));
            }
    out.AddPolygonBlock("polygon", std::move(flat), std::move(offsets));
    NDArray ids2 = NDArray::Uninit(DType::Int64, {boundary.size() / 2, std::size_t{2}});
    std::copy(boundary.begin(), boundary.end(), ids2.As<std::int64_t>());
    out.AddRegion(Region("boundary", RegionKind::Side, 1, -1, std::move(ids2)));
    return out;
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

// ===========================================================================
// line, circle, helix
// ===========================================================================

Mesh line(std::int64_t N, const std::array<double, 3>& rStart, const std::array<double, 3>& rEnd,
          int Order, std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "line";
    prim_check_count(shape, "n", N);
    prim_check_order(shape, Order);
    if (N == 0)
        return prim_empty();
    if (rStart == rEnd)
        prim_fail(shape, "start and end must differ");
    prim_check_budget(shape, N, MaxCells);

    PrimSoup soup;
    for (std::int64_t i = 0; i <= N; ++i) {
        double q[3];
        for (std::size_t k = 0; k < 3; ++k)
            q[k] = i == N ? rEnd[k]
                          : rStart[k] + ((rEnd[k] - rStart[k]) * static_cast<double>(i)) /
                                            static_cast<double>(N);
        soup.AddNode(q[0], q[1], q[2], 0);
    }
    PrimBlock& blk = soup.Block(CellType::Line, 2);
    for (std::int64_t i = 0; i < N; ++i)
        prim_push(blk, {i, i + 1});
    Mesh out = prim_finish(shape, soup, Order, Complete, {}, 0);
    prim_add_point_region(out, "start", 0);
    prim_add_point_region(out, "end", N);
    return out;
}

Mesh circle(double Radius, std::int64_t N, const std::array<double, 3>& rCenter, int Order,
            std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "circle";
    prim_check_positive(shape, "radius", Radius);
    prim_check_count(shape, "n", N);
    prim_check_order(shape, Order);
    if (N == 0)
        return prim_empty();
    if (N < 3)
        prim_fail(shape, "n must be at least 3, got " + std::to_string(N));
    prim_check_budget(shape, N, MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Polar;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    for (std::int64_t i = 0; i < N; ++i)
        soup.AddNode(1.0, (2.0 * kPrimPi * static_cast<double>(i)) / static_cast<double>(N), 0.0,
                     0);
    PrimBlock& blk = soup.Block(CellType::Line, 2);
    for (std::int64_t i = 0; i < N; ++i)
        prim_push(blk, {i, (i + 1) % N});
    return prim_finish(shape, soup, Order, Complete, {}, 0);
}

Mesh helix(double Radius, double Pitch, double Turns, std::int64_t N,
           const std::array<double, 3>& rCenter, int Order, std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "helix";
    prim_check_positive(shape, "radius", Radius);
    prim_check_positive(shape, "turns", Turns);
    prim_check_count(shape, "n", N);
    prim_check_order(shape, Order);
    if (N == 0)
        return prim_empty();
    if (!(Pitch != 0.0) || !std::isfinite(Pitch))
        prim_fail(shape, "pitch must be finite and nonzero, got " + std::to_string(Pitch));
    if (!(static_cast<double>(N) > 2.0 * Turns))
        prim_fail(shape, "n must exceed 2 * turns so no cell spans half a turn, got n " +
                             std::to_string(N) + " and turns " + std::to_string(Turns));
    prim_check_budget(shape, N, MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Polar;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    for (std::int64_t i = 0; i <= N; ++i)
        soup.AddNode(1.0, (2.0 * kPrimPi * Turns * static_cast<double>(i)) / static_cast<double>(N),
                     (Pitch * Turns * static_cast<double>(i)) / static_cast<double>(N), 0);
    PrimBlock& blk = soup.Block(CellType::Line, 2);
    for (std::int64_t i = 0; i < N; ++i)
        prim_push(blk, {i, i + 1});
    Mesh out = prim_finish(shape, soup, Order, Complete, {}, 0);
    prim_add_point_region(out, "start", 0);
    prim_add_point_region(out, "end", N);
    return out;
}

// ===========================================================================
// sphere_points, box_points
// ===========================================================================

namespace {

std::uint64_t prim_splitmix(std::uint64_t& rState) {
    rState += 0x9E3779B97F4A7C15ULL;
    std::uint64_t z = rState;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// A uniform draw in [0, 1): the top 53 bits, exact in a double.
double prim_draw(std::uint64_t& rState) {
    return static_cast<double>(prim_splitmix(rState) >> 11) * (1.0 / 9007199254740992.0);
}

Mesh prim_vertices(PrimSoup& rSoup) {
    const std::int64_t n = static_cast<std::int64_t>(rSoup.mMask.size());
    PrimBlock& blk = rSoup.Block(CellType::Vertex, 1);
    blk.mConn.resize(static_cast<std::size_t>(n));
    for (std::int64_t i = 0; i < n; ++i)
        blk.mConn[static_cast<std::size_t>(i)] = i;
    return prim_build(rSoup);
}

}  // namespace

Mesh sphere_points(std::int64_t N, double Radius, const std::array<double, 3>& rCenter,
                   std::int64_t MaxCells) {
    constexpr const char* shape = "sphere_points";
    prim_check_count(shape, "n", N);
    if (N == 0)
        return prim_empty();
    prim_check_positive(shape, "radius", Radius);
    prim_check_budget(shape, N, MaxCells);
    const double golden = kPrimPi * (3.0 - std::sqrt(5.0));
    PrimSoup soup;
    for (std::int64_t i = 0; i < N; ++i) {
        const double z = 1.0 - static_cast<double>(2 * i + 1) / static_cast<double>(N);
        const double rho = std::sqrt(1.0 - z * z);
        const double phi = golden * static_cast<double>(i);
        soup.AddNode(rCenter[0] + Radius * rho * std::cos(phi),
                     rCenter[1] + Radius * rho * std::sin(phi), rCenter[2] + Radius * z, 0);
    }
    return prim_vertices(soup);
}

Mesh box_points(std::int64_t N, const std::array<double, 3>& rSize,
                const std::array<double, 3>& rOrigin, std::uint64_t Seed, std::int64_t MaxCells) {
    constexpr const char* shape = "box_points";
    prim_check_count(shape, "n", N);
    if (N == 0)
        return prim_empty();
    for (std::size_t k = 0; k < 3; ++k)
        prim_check_positive(shape, "size", rSize[k]);
    prim_check_budget(shape, N, MaxCells);
    PrimSoup soup;
    std::uint64_t state = Seed;
    for (std::int64_t i = 0; i < N; ++i) {
        const double u = prim_draw(state);
        const double v = prim_draw(state);
        const double w = prim_draw(state);
        soup.AddNode(rOrigin[0] + rSize[0] * u, rOrigin[1] + rSize[1] * v,
                     rOrigin[2] + rSize[2] * w, 0);
    }
    return prim_vertices(soup);
}

// ===========================================================================
// tetra_sphere, octa_sphere
// ===========================================================================

namespace {

template <std::size_t NV, std::size_t NF>
Mesh prim_seeded_sphere(const char* pShape, const double (&rVerts)[NV][3],
                        const std::int64_t (&rFaces)[NF][3], double Radius,
                        std::int64_t Subdivisions, const std::array<double, 3>& rCenter, int Order,
                        std::int64_t MaxCells, bool Complete) {
    prim_check_positive(pShape, "radius", Radius);
    prim_check_count(pShape, "subdivisions", Subdivisions);
    prim_check_order(pShape, Order);
    std::int64_t cells = static_cast<std::int64_t>(NF);
    for (std::int64_t s = 0; s < Subdivisions; ++s)
        cells = prim_mul(cells, 4);
    prim_check_budget(pShape, cells, MaxCells);
    PrimSoup soup;
    soup.mMap = PrimMap::Unit;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    for (const auto& v : rVerts)
        soup.AddNode(v[0], v[1], v[2], 0);
    std::vector<std::int64_t> faces;
    for (const auto& f : rFaces)
        faces.insert(faces.end(), {f[0], f[1], f[2]});
    prim_subdivide_sphere(soup, faces, Subdivisions);
    PrimBlock& blk = soup.Block(CellType::Triangle, 3);
    blk.mConn = std::move(faces);
    return prim_finish(pShape, soup, Order, Complete, {}, 2);
}

}  // namespace

Mesh tetra_sphere(double Radius, std::int64_t Subdivisions, const std::array<double, 3>& rCenter,
                  int Order, std::int64_t MaxCells, bool Complete) {
    const double a = 1.0 / std::sqrt(3.0);
    const double verts[4][3] = {{a, a, a}, {a, -a, -a}, {-a, a, -a}, {-a, -a, a}};
    static const std::int64_t kFaces[4][3] = {{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2}};
    return prim_seeded_sphere("tetra_sphere", verts, kFaces, Radius, Subdivisions, rCenter, Order,
                              MaxCells, Complete);
}

Mesh octa_sphere(double Radius, std::int64_t Subdivisions, const std::array<double, 3>& rCenter,
                 int Order, std::int64_t MaxCells, bool Complete) {
    static const double kVerts[6][3] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                                        {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
    static const std::int64_t kFaces[8][3] = {{0, 2, 4}, {0, 3, 5}, {0, 4, 3}, {0, 5, 2},
                                              {1, 2, 5}, {1, 3, 4}, {1, 4, 2}, {1, 5, 3}};
    return prim_seeded_sphere("octa_sphere", kVerts, kFaces, Radius, Subdivisions, rCenter, Order,
                              MaxCells, Complete);
}

// ===========================================================================
// The boundary lattice of a box: cube_sphere, box_surface, spherical_shell
// ===========================================================================

namespace {

// The nodes of an nx x ny x nz lattice on the box boundary, numbered in lattice
// order (x fastest, then y, then z) with the interior skipped: the full bottom
// and top layers, and a ring on each layer between.
struct PrimCubeLattice {
    std::int64_t mNx, mNy, mNz;

    std::int64_t Full() const { return (mNx + 1) * (mNy + 1); }
    std::int64_t Ring() const { return 2 * (mNx + 1) + 2 * (mNy - 1); }
    std::int64_t Count() const { return 2 * Full() + (mNz - 1) * Ring(); }

    std::int64_t Id(std::int64_t I, std::int64_t J, std::int64_t K) const {
        if (K == 0)
            return J * (mNx + 1) + I;
        const std::int64_t base = Full() + (K - 1) * Ring();
        if (K == mNz)
            return base + J * (mNx + 1) + I;
        if (J == 0)
            return base + I;
        if (J == mNy)
            return base + (mNx + 1) + 2 * (mNy - 1) + I;
        return base + (mNx + 1) + 2 * (J - 1) + (I == 0 ? 0 : 1);
    }

    bool OnBoundary(std::int64_t I, std::int64_t J, std::int64_t K) const {
        return I == 0 || I == mNx || J == 0 || J == mNy || K == 0 || K == mNz;
    }

    // Every face quad, oriented outward, face by face: xmin, xmax, ymin, ymax,
    // zmin, zmax; on each, the second in-plane index outer.
    template <typename Emit>
    void Faces(Emit&& rEmit) const {
        const std::int64_t nx = mNx, ny = mNy, nz = mNz;
        for (std::int64_t k = 0; k < nz; ++k)
            for (std::int64_t j = 0; j < ny; ++j)
                rEmit(0, Id(0, j, k), Id(0, j, k + 1), Id(0, j + 1, k + 1), Id(0, j + 1, k));
        for (std::int64_t k = 0; k < nz; ++k)
            for (std::int64_t j = 0; j < ny; ++j)
                rEmit(1, Id(nx, j, k), Id(nx, j + 1, k), Id(nx, j + 1, k + 1), Id(nx, j, k + 1));
        for (std::int64_t k = 0; k < nz; ++k)
            for (std::int64_t i = 0; i < nx; ++i)
                rEmit(2, Id(i, 0, k), Id(i + 1, 0, k), Id(i + 1, 0, k + 1), Id(i, 0, k + 1));
        for (std::int64_t k = 0; k < nz; ++k)
            for (std::int64_t i = 0; i < nx; ++i)
                rEmit(3, Id(i, ny, k), Id(i, ny, k + 1), Id(i + 1, ny, k + 1), Id(i + 1, ny, k));
        for (std::int64_t j = 0; j < ny; ++j)
            for (std::int64_t i = 0; i < nx; ++i)
                rEmit(4, Id(i, j, 0), Id(i, j + 1, 0), Id(i + 1, j + 1, 0), Id(i + 1, j, 0));
        for (std::int64_t j = 0; j < ny; ++j)
            for (std::int64_t i = 0; i < nx; ++i)
                rEmit(5, Id(i, j, nz), Id(i + 1, j, nz), Id(i + 1, j + 1, nz), Id(i, j + 1, nz));
    }

    std::int64_t NumFaces() const {
        return 2 * (prim_mul(mNx, mNy) + prim_mul(mNy, mNz) + prim_mul(mNx, mNz));
    }
};

// The boundary nodes of the lattice in order, each through rAdd(i, j, k).
template <typename Add>
void prim_cube_nodes(const PrimCubeLattice& rLat, Add&& rAdd) {
    for (std::int64_t k = 0; k <= rLat.mNz; ++k)
        for (std::int64_t j = 0; j <= rLat.mNy; ++j)
            for (std::int64_t i = 0; i <= rLat.mNx; ++i)
                if (rLat.OnBoundary(i, j, k))
                    rAdd(i, j, k);
}

}  // namespace

Mesh cube_sphere(double Radius, std::int64_t N, const std::array<double, 3>& rCenter,
                 const std::string& rCellType, int Order, std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "cube_sphere";
    prim_check_positive(shape, "radius", Radius);
    prim_check_count(shape, "n", N);
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (N == 0)
        return prim_empty();
    const PrimCubeLattice lat{N, N, N};
    prim_check_budget(shape, prim_mul(lat.NumFaces(), tri ? 2 : 1), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Project;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    const double n = static_cast<double>(N);
    prim_cube_nodes(lat, [&](std::int64_t I, std::int64_t J, std::int64_t K) {
        soup.AddNode(static_cast<double>(2 * I - N) / n, static_cast<double>(2 * J - N) / n,
                     static_cast<double>(2 * K - N) / n, 0);
    });
    PrimBlock& blk = soup.Block(tri ? CellType::Triangle : CellType::Quad, tri ? 3 : 4);
    lat.Faces([&](int, std::int64_t A, std::int64_t B, std::int64_t C, std::int64_t D) {
        prim_push_quad(blk, tri, A, B, C, D);
    });
    return prim_finish(shape, soup, Order, Complete, {}, 2);
}

Mesh box_surface(const std::array<std::int64_t, 3>& rDims, const std::array<double, 3>& rSize,
                 const std::array<double, 3>& rOrigin, const std::string& rCellType, int Order,
                 std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "box_surface";
    for (std::size_t k = 0; k < 3; ++k) {
        prim_check_count(shape, "cell counts", rDims[k]);
        if (rDims[k] > 0)
            prim_check_positive(shape, "size", rSize[k]);
    }
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (rDims[0] == 0 || rDims[1] == 0 || rDims[2] == 0)
        return prim_empty();
    const PrimCubeLattice lat{rDims[0], rDims[1], rDims[2]};
    prim_check_budget(shape, prim_mul(lat.NumFaces(), tri ? 2 : 1), MaxCells);

    PrimSoup soup;
    prim_cube_nodes(lat, [&](std::int64_t I, std::int64_t J, std::int64_t K) {
        soup.AddNode(
            rOrigin[0] + (rSize[0] * static_cast<double>(I)) / static_cast<double>(rDims[0]),
            rOrigin[1] + (rSize[1] * static_cast<double>(J)) / static_cast<double>(rDims[1]),
            rOrigin[2] + (rSize[2] * static_cast<double>(K)) / static_cast<double>(rDims[2]), 0);
    });
    PrimBlock& blk = soup.Block(tri ? CellType::Triangle : CellType::Quad, tri ? 3 : 4);
    std::vector<std::uint8_t> labels;
    lat.Faces([&](int Face, std::int64_t A, std::int64_t B, std::int64_t C, std::int64_t D) {
        prim_push_quad(blk, tri, A, B, C, D);
        labels.insert(labels.end(), tri ? 2 : 1, static_cast<std::uint8_t>(Face));
    });
    Mesh out = prim_finish(shape, soup, Order, Complete, {}, 2);
    prim_add_cell_regions(out, labels, {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}, 2);
    return out;
}

// ===========================================================================
// geo_sphere
// ===========================================================================

Mesh geo_sphere(double Radius, std::int64_t N, const std::array<double, 3>& rCenter, int Order,
                std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "geo_sphere";
    prim_check_positive(shape, "radius", Radius);
    prim_check_count(shape, "n", N);
    prim_check_order(shape, Order);
    if (N == 0)
        return prim_empty();
    prim_check_budget(shape, prim_mul(20, prim_mul(N, N)), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Project;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    const double t = (1.0 + std::sqrt(5.0)) / 2.0;
    const double len = std::sqrt(1.0 + t * t);
    const double one = 1.0 / len, tau = t / len;
    const double verts[12][3] = {{-one, tau, 0}, {one, tau, 0}, {-one, -tau, 0}, {one, -tau, 0},
                                 {0, -one, tau}, {0, one, tau}, {0, -one, -tau}, {0, one, -tau},
                                 {tau, 0, -one}, {tau, 0, one}, {-tau, 0, -one}, {-tau, 0, one}};
    static const std::int64_t kFaces[20][3] = {
        {0, 11, 5},  {0, 5, 1},  {0, 1, 7},  {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4},
        {11, 10, 2}, {10, 7, 6}, {7, 1, 8},  {3, 9, 4},  {3, 4, 2},   {3, 2, 6}, {3, 6, 8},
        {3, 8, 9},   {4, 9, 5},  {2, 4, 11}, {6, 2, 10}, {8, 6, 7},   {9, 8, 1}};
    for (const auto& v : verts)
        soup.AddNode(v[0], v[1], v[2], 0);
    const double n = static_cast<double>(N);
    // The 30 edges in first-appearance order over the faces' (v0 v1, v1 v2,
    // v2 v0), then N - 1 nodes on each from its lower vertex.
    std::unordered_map<std::uint64_t, std::int64_t> edge_index;
    std::vector<std::pair<std::int64_t, std::int64_t>> edges;
    for (const auto& f : kFaces)
        for (int e = 0; e < 3; ++e) {
            std::int64_t a = f[e], b = f[(e + 1) % 3];
            if (a > b)
                std::swap(a, b);
            if (edge_index.emplace(prim_edge_key(a, b), static_cast<std::int64_t>(edges.size()))
                    .second)
                edges.emplace_back(a, b);
        }
    const std::int64_t edge_base = 12;
    for (const auto& [a, b] : edges)
        for (std::int64_t s = 1; s < N; ++s) {
            const double wa = static_cast<double>(N - s), wb = static_cast<double>(s);
            soup.AddNode((verts[a][0] * wa + verts[b][0] * wb) / n,
                         (verts[a][1] * wa + verts[b][1] * wb) / n,
                         (verts[a][2] * wa + verts[b][2] * wb) / n, 0);
        }
    const std::int64_t face_base = edge_base + 30 * (N - 1);
    const std::int64_t per_face = ((N - 1) * (N - 2)) / 2;
    for (const auto& f : kFaces)
        for (std::int64_t j = 1; j + 1 < N; ++j)
            for (std::int64_t i = 1; i + j < N; ++i) {
                const double w0 = static_cast<double>(N - i - j);
                const double w1 = static_cast<double>(i), w2 = static_cast<double>(j);
                soup.AddNode((verts[f[0]][0] * w0 + verts[f[1]][0] * w1 + verts[f[2]][0] * w2) / n,
                             (verts[f[0]][1] * w0 + verts[f[1]][1] * w1 + verts[f[2]][1] * w2) / n,
                             (verts[f[0]][2] * w0 + verts[f[1]][2] * w1 + verts[f[2]][2] * w2) / n,
                             0);
            }
    // The node at (N - i - j, i, j) of face F over its corners (v0, v1, v2).
    auto on_edge = [&](std::int64_t From, std::int64_t To, std::int64_t Steps) {
        const std::int64_t lo = std::min(From, To), hi = std::max(From, To);
        const std::int64_t e = edge_index.at(prim_edge_key(lo, hi));
        const std::int64_t s = From == lo ? Steps : N - Steps;
        return edge_base + e * (N - 1) + (s - 1);
    };
    auto id = [&](std::size_t F, std::int64_t I, std::int64_t J) -> std::int64_t {
        const auto& f = kFaces[F];
        if (I == 0 && J == 0)
            return f[0];
        if (I == N)
            return f[1];
        if (J == N)
            return f[2];
        if (J == 0)
            return on_edge(f[0], f[1], I);
        if (I == 0)
            return on_edge(f[0], f[2], J);
        if (I + J == N)
            return on_edge(f[1], f[2], J);
        return face_base + static_cast<std::int64_t>(F) * per_face + (J - 1) * (N - 1) -
               ((J - 1) * J) / 2 + (I - 1);
    };
    PrimBlock& blk = soup.Block(CellType::Triangle, 3);
    for (std::size_t f = 0; f < 20; ++f)
        for (std::int64_t j = 0; j < N; ++j)
            for (std::int64_t i = 0; i < N - j; ++i) {
                prim_push(blk, {id(f, i, j), id(f, i + 1, j), id(f, i, j + 1)});
                if (i < N - 1 - j)
                    prim_push(blk, {id(f, i + 1, j), id(f, i + 1, j + 1), id(f, i, j + 1)});
            }
    return prim_finish(shape, soup, Order, Complete, {}, 2);
}

// ===========================================================================
// torus, moebius, tube
// ===========================================================================

Mesh torus(double MajorRadius, double MinorRadius, std::int64_t Nu, std::int64_t Nv,
           const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
           std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "torus";
    prim_check_positive(shape, "major_radius", MajorRadius);
    prim_check_positive(shape, "minor_radius", MinorRadius);
    if (!(MinorRadius < MajorRadius))
        prim_fail(shape, "minor_radius must be less than major_radius, got " +
                             std::to_string(MinorRadius) + " and " + std::to_string(MajorRadius));
    prim_check_count(shape, "nu", Nu);
    prim_check_count(shape, "nv", Nv);
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (Nu == 0 || Nv == 0)
        return prim_empty();
    if (Nu < 3 || Nv < 3)
        prim_fail(shape, "nu and nv must be at least 3, got " + std::to_string(Nu) + " and " +
                             std::to_string(Nv));
    prim_check_budget(shape, prim_mul(prim_mul(Nu, Nv), tri ? 2 : 1), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Torus;
    soup.mRadius = MajorRadius;
    soup.mRadius2 = MinorRadius;
    soup.mCenter = rCenter;
    for (std::int64_t i = 0; i < Nu; ++i)
        for (std::int64_t j = 0; j < Nv; ++j)
            soup.AddNode((2.0 * kPrimPi * static_cast<double>(i)) / static_cast<double>(Nu),
                         (2.0 * kPrimPi * static_cast<double>(j)) / static_cast<double>(Nv), 0.0,
                         0);
    PrimBlock& blk = soup.Block(tri ? CellType::Triangle : CellType::Quad, tri ? 3 : 4);
    for (std::int64_t i = 0; i < Nu; ++i)
        for (std::int64_t j = 0; j < Nv; ++j) {
            const std::int64_t i1 = (i + 1) % Nu, j1 = (j + 1) % Nv;
            prim_push_quad(blk, tri, i * Nv + j, i1 * Nv + j, i1 * Nv + j1, i * Nv + j1);
        }
    return prim_finish(shape, soup, Order, Complete, {}, 2);
}

Mesh moebius(double Radius, double Width, std::int64_t Nu, std::int64_t Nv, std::int64_t NumTwists,
             const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
             std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "moebius";
    prim_check_positive(shape, "radius", Radius);
    prim_check_positive(shape, "width", Width);
    if (!(Width < 2.0 * Radius))
        prim_fail(shape, "width must be less than twice the radius, got " + std::to_string(Width) +
                             " and " + std::to_string(Radius));
    prim_check_count(shape, "nu", Nu);
    prim_check_count(shape, "nv", Nv);
    prim_check_count(shape, "num_twists", NumTwists);
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (Nu == 0 || Nv == 0)
        return prim_empty();
    if (Nu < 3)
        prim_fail(shape, "nu must be at least 3, got " + std::to_string(Nu));
    prim_check_budget(shape, prim_mul(prim_mul(Nu, Nv), tri ? 2 : 1), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Moebius;
    soup.mRadius = Radius;
    soup.mHalf = Width / 2.0;
    soup.mTwists = NumTwists;
    soup.mCenter = rCenter;
    const std::int64_t p = Nv + 1;
    for (std::int64_t i = 0; i < Nu; ++i)
        for (std::int64_t j = 0; j <= Nv; ++j)
            soup.AddNode((2.0 * kPrimPi * static_cast<double>(i)) / static_cast<double>(Nu),
                         static_cast<double>(2 * j - Nv) / static_cast<double>(Nv), 0.0,
                         (j == 0 || j == Nv) ? 1 : 0);
    // Across the seam, an odd count of half-twists meets the first ring upside
    // down: node j there is node Nv - j.
    const bool flip = NumTwists % 2 != 0;
    PrimBlock& blk = soup.Block(tri ? CellType::Triangle : CellType::Quad, tri ? 3 : 4);
    for (std::int64_t i = 0; i < Nu; ++i)
        for (std::int64_t j = 0; j < Nv; ++j) {
            const std::int64_t a = i * p + j, d = a + 1;
            std::int64_t b, c;
            if (i + 1 < Nu) {
                b = (i + 1) * p + j;
                c = b + 1;
            } else {
                b = flip ? Nv - j : j;
                c = flip ? Nv - j - 1 : j + 1;
            }
            prim_push_quad(blk, tri, a, b, c, d);
        }
    return prim_finish(shape, soup, Order, Complete, {{"boundary", 1}}, 1);
}

Mesh tube(double Radius, double Height, std::int64_t Ntheta, std::int64_t Nz,
          const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
          std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "tube";
    prim_check_positive(shape, "radius", Radius);
    prim_check_positive(shape, "height", Height);
    prim_check_count(shape, "ntheta", Ntheta);
    prim_check_count(shape, "nz", Nz);
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (Ntheta == 0 || Nz == 0)
        return prim_empty();
    if (Ntheta < 3)
        prim_fail(shape, "ntheta must be at least 3, got " + std::to_string(Ntheta));
    prim_check_budget(shape, prim_mul(prim_mul(Ntheta, Nz), tri ? 2 : 1), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Polar;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    for (std::int64_t k = 0; k <= Nz; ++k) {
        const double z = (Height * static_cast<double>(2 * k - Nz)) / static_cast<double>(2 * Nz);
        const std::uint16_t mask = (k == 0 ? 1 : 0) | (k == Nz ? 2 : 0);
        for (std::int64_t i = 0; i < Ntheta; ++i)
            soup.AddNode(1.0,
                         (2.0 * kPrimPi * static_cast<double>(i)) / static_cast<double>(Ntheta), z,
                         mask);
    }
    PrimBlock& blk = soup.Block(tri ? CellType::Triangle : CellType::Quad, tri ? 3 : 4);
    for (std::int64_t k = 0; k < Nz; ++k)
        for (std::int64_t i = 0; i < Ntheta; ++i) {
            const std::int64_t i1 = (i + 1) % Ntheta;
            prim_push_quad(blk, tri, k * Ntheta + i, k * Ntheta + i1, (k + 1) * Ntheta + i1,
                           (k + 1) * Ntheta + i);
        }
    return prim_finish(shape, soup, Order, Complete, {{"bottom", 1}, {"top", 2}}, 1);
}

// ===========================================================================
// cone, capsule: closed surfaces of revolution in rings between two poles
// ===========================================================================

namespace {

// Cells between a pole (node 0), Rings rings of Ntheta nodes (ring r from 1)
// and a second pole (node Last), as uv_sphere lays them out: quads between the
// rings and triangle fans at the poles (or triangles throughout), oriented
// outward when the rings run from the first pole to the second with the angle
// increasing counter-clockwise seen from the first. rLabel(r) labels the band
// between rings r and r + 1 (r = 0 the first fan, r = Rings the second).
template <typename Label>
void prim_ring_cells(PrimSoup& rSoup, std::int64_t Rings, std::int64_t Ntheta, std::int64_t Last,
                     bool Tri, Label&& rLabel, std::vector<std::uint8_t>& rLabels) {
    auto ring = [&](std::int64_t R, std::int64_t I) { return 1 + (R - 1) * Ntheta + (I % Ntheta); };
    if (!Tri) {
        PrimBlock& body = rSoup.Block(CellType::Quad, 4);
        for (std::int64_t r = 1; r < Rings; ++r)
            for (std::int64_t i = 0; i < Ntheta; ++i) {
                prim_push(body, {ring(r, i), ring(r + 1, i), ring(r + 1, i + 1), ring(r, i + 1)});
                rLabels.push_back(rLabel(r));
            }
    }
    PrimBlock& tris = rSoup.Block(CellType::Triangle, 3);
    for (std::int64_t i = 0; i < Ntheta; ++i) {
        prim_push(tris, {0, ring(1, i), ring(1, i + 1)});
        rLabels.push_back(rLabel(0));
    }
    if (Tri)
        for (std::int64_t r = 1; r < Rings; ++r)
            for (std::int64_t i = 0; i < Ntheta; ++i) {
                const std::int64_t a = ring(r, i), b = ring(r, i + 1);
                const std::int64_t c = ring(r + 1, i + 1), d = ring(r + 1, i);
                prim_push(tris, {a, d, c});
                prim_push(tris, {a, c, b});
                rLabels.insert(rLabels.end(), 2, rLabel(r));
            }
    for (std::int64_t i = 0; i < Ntheta; ++i) {
        prim_push(tris, {Last, ring(Rings, i + 1), ring(Rings, i)});
        rLabels.push_back(rLabel(Rings));
    }
}

std::int64_t prim_ring_cell_count(std::int64_t Rings, std::int64_t Ntheta, bool Tri) {
    return prim_mul(Ntheta, Tri ? 2 * Rings : Rings + 1);
}

}  // namespace

Mesh cone(double Radius, double Height, std::int64_t Ntheta, std::int64_t Nz,
          const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
          std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "cone";
    prim_check_positive(shape, "radius", Radius);
    prim_check_positive(shape, "height", Height);
    prim_check_count(shape, "ntheta", Ntheta);
    prim_check_count(shape, "nz", Nz);
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (Ntheta == 0 || Nz == 0)
        return prim_empty();
    if (Ntheta < 3)
        prim_fail(shape, "ntheta must be at least 3, got " + std::to_string(Ntheta));
    const std::int64_t rings = 2 * Nz - 1;
    prim_check_budget(shape, prim_ring_cell_count(rings, Ntheta, tri), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Polar;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    soup.AddNode(0.0, 0.0, Height, kPrimDegenerate);
    for (std::int64_t r = 1; r <= rings; ++r) {
        // Down the slant to the rim (r = Nz), then in across the base.
        const std::int64_t m = r <= Nz ? r : 2 * Nz - r;
        const double rho = static_cast<double>(m) / static_cast<double>(Nz);
        const double z =
            r <= Nz ? (Height * static_cast<double>(Nz - r)) / static_cast<double>(Nz) : 0.0;
        for (std::int64_t i = 0; i < Ntheta; ++i)
            soup.AddNode(
                rho, (2.0 * kPrimPi * static_cast<double>(i)) / static_cast<double>(Ntheta), z, 0);
    }
    const std::int64_t last = soup.AddNode(0.0, 0.0, 0.0, kPrimDegenerate);
    std::vector<std::uint8_t> labels;
    prim_ring_cells(
        soup, rings, Ntheta, last, tri,
        [&](std::int64_t R) { return static_cast<std::uint8_t>(R < Nz ? 0 : 1); }, labels);
    Mesh out = prim_finish(shape, soup, Order, Complete, {}, 2);
    prim_add_cell_regions(out, labels, {"side", "base"}, 2);
    return out;
}

Mesh capsule(double Radius, double Height, std::int64_t Ntheta, std::int64_t Nz, std::int64_t Ncap,
             const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
             std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "capsule";
    prim_check_positive(shape, "radius", Radius);
    prim_check_positive(shape, "height", Height);
    prim_check_count(shape, "ntheta", Ntheta);
    prim_check_count(shape, "nz", Nz);
    prim_check_count(shape, "ncap", Ncap);
    prim_check_order(shape, Order);
    const bool tri = rCellType == "triangle";
    if (!tri && rCellType != "quad")
        prim_bad_choice(shape, "cell_type", rCellType, "quad, triangle");
    if (Ntheta == 0 || Nz == 0 || Ncap == 0)
        return prim_empty();
    if (Ntheta < 3)
        prim_fail(shape, "ntheta must be at least 3, got " + std::to_string(Ntheta));
    const std::int64_t rings = 2 * Ncap + Nz - 1;
    prim_check_budget(shape, prim_ring_cell_count(rings, Ntheta, tri), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Capsule;
    soup.mRadius = Radius;
    soup.mHalf = Height / 2.0;
    soup.mCenter = rCenter;
    const double side = (2.0 * soup.mHalf) / Radius;
    soup.AddNode(0.0, 0.0, 0.0, kPrimDegenerate);
    for (std::int64_t r = 1; r <= rings; ++r) {
        double t;
        if (r <= Ncap)
            t = (kPrimPi * static_cast<double>(r)) / static_cast<double>(2 * Ncap);
        else if (r < Ncap + Nz)
            t = kPrimPi / 2.0 + (side * static_cast<double>(r - Ncap)) / static_cast<double>(Nz);
        else
            t = kPrimPi / 2.0 + side +
                (kPrimPi * static_cast<double>(r - Ncap - Nz)) / static_cast<double>(2 * Ncap);
        for (std::int64_t i = 0; i < Ntheta; ++i)
            soup.AddNode(t, (2.0 * kPrimPi * static_cast<double>(i)) / static_cast<double>(Ntheta),
                         0.0, 0);
    }
    const std::int64_t last = soup.AddNode(kPrimPi + side, 0.0, 0.0, kPrimDegenerate);
    std::vector<std::uint8_t> labels;
    prim_ring_cells(
        soup, rings, Ntheta, last, tri, [](std::int64_t) { return std::uint8_t{0}; }, labels);
    return prim_finish(shape, soup, Order, Complete, {}, 2);
}

// ===========================================================================
// platonic
// ===========================================================================

Mesh platonic(const std::string& rName, double Radius, const std::array<double, 3>& rCenter,
              std::int64_t MaxCells) {
    constexpr const char* shape = "platonic";
    prim_check_positive(shape, "radius", Radius);
    const double a = 1.0 / std::sqrt(3.0);
    std::vector<std::array<double, 3>> verts;
    std::vector<std::vector<std::int64_t>> faces;
    if (rName == "tetrahedron") {
        verts = {{a, a, a}, {a, -a, -a}, {-a, a, -a}, {-a, -a, a}};
        faces = {{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2}};
    } else if (rName == "cube") {
        for (int z = -1; z <= 1; z += 2)
            for (int y = -1; y <= 1; y += 2)
                for (int x = -1; x <= 1; x += 2)
                    verts.push_back({x * a, y * a, z * a});
        faces = {{0, 1, 5, 4}, {0, 2, 3, 1}, {0, 4, 6, 2},
                 {1, 3, 7, 5}, {2, 6, 7, 3}, {4, 5, 7, 6}};
    } else if (rName == "octahedron") {
        verts = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
        faces = {{0, 2, 4}, {0, 3, 5}, {0, 4, 3}, {0, 5, 2},
                 {1, 2, 5}, {1, 3, 4}, {1, 4, 2}, {1, 5, 3}};
    } else if (rName == "dodecahedron") {
        const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
        const double ip = 1.0 / phi;
        for (int z = -1; z <= 1; z += 2)
            for (int y = -1; y <= 1; y += 2)
                for (int x = -1; x <= 1; x += 2)
                    verts.push_back({x * a, y * a, z * a});
        for (int s1 = -1; s1 <= 1; s1 += 2)
            for (int s2 = -1; s2 <= 1; s2 += 2)
                verts.push_back({0.0, (s1 * ip) * a, (s2 * phi) * a});
        for (int s1 = -1; s1 <= 1; s1 += 2)
            for (int s2 = -1; s2 <= 1; s2 += 2)
                verts.push_back({(s1 * ip) * a, (s2 * phi) * a, 0.0});
        for (int s1 = -1; s1 <= 1; s1 += 2)
            for (int s2 = -1; s2 <= 1; s2 += 2)
                verts.push_back({(s1 * phi) * a, 0.0, (s2 * ip) * a});
        faces = {{0, 8, 1, 14, 12},  {0, 12, 4, 17, 16}, {0, 16, 2, 10, 8},  {1, 8, 10, 3, 18},
                 {1, 18, 19, 5, 14}, {2, 13, 15, 3, 10}, {2, 16, 17, 6, 13}, {3, 15, 7, 19, 18},
                 {4, 9, 11, 6, 17},  {4, 12, 14, 5, 9},  {5, 19, 7, 11, 9},  {6, 11, 7, 15, 13}};
    } else if (rName == "icosahedron") {
        const double t = (1.0 + std::sqrt(5.0)) / 2.0;
        const double len = std::sqrt(1.0 + t * t);
        const double one = 1.0 / len, tau = t / len;
        verts = {{-one, tau, 0}, {one, tau, 0}, {-one, -tau, 0}, {one, -tau, 0},
                 {0, -one, tau}, {0, one, tau}, {0, -one, -tau}, {0, one, -tau},
                 {tau, 0, -one}, {tau, 0, one}, {-tau, 0, -one}, {-tau, 0, one}};
        faces = {{0, 11, 5}, {0, 5, 1},  {0, 1, 7},   {0, 7, 10}, {0, 10, 11},
                 {1, 5, 9},  {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                 {3, 9, 4},  {3, 4, 2},  {3, 2, 6},   {3, 6, 8},  {3, 8, 9},
                 {4, 9, 5},  {2, 4, 11}, {6, 2, 10},  {8, 6, 7},  {9, 8, 1}};
    } else {
        prim_bad_choice(shape, "name", rName,
                        "tetrahedron, cube, octahedron, dodecahedron, icosahedron");
    }
    prim_check_budget(shape, static_cast<std::int64_t>(faces.size()), MaxCells);

    Mesh out;
    NDArray points = NDArray::Uninit(DType::Float64, {verts.size(), std::size_t{3}});
    double* dst = points.As<double>();
    for (std::size_t i = 0; i < verts.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            dst[3 * i + k] = rCenter[k] + Radius * verts[i][k];
    out.AssignPoints(std::move(points));
    const std::size_t width = faces[0].size();
    if (width == 5) {
        std::vector<std::int64_t> flat, offsets{0};
        for (const auto& f : faces) {
            flat.insert(flat.end(), f.begin(), f.end());
            offsets.push_back(static_cast<std::int64_t>(flat.size()));
        }
        out.AddPolygonBlock("polygon", std::move(flat), std::move(offsets));
    } else {
        NDArray conn = NDArray::Uninit(DType::Int64, {faces.size(), width});
        std::int64_t* c = conn.As<std::int64_t>();
        for (std::size_t f = 0; f < faces.size(); ++f)
            std::copy(faces[f].begin(), faces[f].end(), c + f * width);
        out.AddCellBlock(width == 3 ? "triangle" : "quad", std::move(conn));
    }
    return out;
}

// ===========================================================================
// ball, spherical_shell, pipe, fichera, solid_torus
// ===========================================================================

namespace {

void prim_push_hex(PrimBlock& rBlock, bool Tetra, const std::int64_t* pHex) {
    static const std::int64_t kKuhn[6][4] = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6},
                                             {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    if (Tetra) {
        for (const auto& t : kKuhn)
            prim_push(rBlock, {pHex[t[0]], pHex[t[1]], pHex[t[2]], pHex[t[3]]});
    } else {
        prim_push(rBlock, {pHex[0], pHex[1], pHex[2], pHex[3], pHex[4], pHex[5], pHex[6], pHex[7]});
    }
}

}  // namespace

Mesh ball(double Radius, std::int64_t N, const std::array<double, 3>& rCenter,
          const std::string& rCellType, int Order, std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "ball";
    prim_check_positive(shape, "radius", Radius);
    prim_check_count(shape, "n", N);
    prim_check_order(shape, Order);
    const bool tet = rCellType == "tetra";
    if (!tet && rCellType != "hexahedron")
        prim_bad_choice(shape, "cell_type", rCellType, "hexahedron, tetra");
    if (N == 0)
        return prim_empty();
    const std::int64_t m = 2 * N;
    prim_check_budget(shape, prim_mul(prim_mul(prim_mul(m, m), m), tet ? 6 : 1), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Ball;
    soup.mRadius = Radius;
    soup.mCenter = rCenter;
    const double n = static_cast<double>(N);
    for (std::int64_t k = 0; k <= m; ++k)
        for (std::int64_t j = 0; j <= m; ++j)
            for (std::int64_t i = 0; i <= m; ++i) {
                const bool rim = i == 0 || i == m || j == 0 || j == m || k == 0 || k == m;
                soup.AddNode(static_cast<double>(i - N) / n, static_cast<double>(j - N) / n,
                             static_cast<double>(k - N) / n, rim ? 1 : 0);
            }
    PrimBlock& blk = soup.Block(tet ? CellType::Tetra : CellType::Hexahedron, tet ? 4 : 8);
    for (std::int64_t k = 0; k < m; ++k)
        for (std::int64_t j = 0; j < m; ++j)
            for (std::int64_t i = 0; i < m; ++i) {
                std::int64_t h[8];
                prim_hex_nodes(i, j, k, m + 1, m + 1, h);
                prim_push_hex(blk, tet, h);
            }
    return prim_finish(shape, soup, Order, Complete, {{"surface", 1}}, 2);
}

Mesh spherical_shell(double InnerRadius, double OuterRadius, std::int64_t N, std::int64_t Nr,
                     double Ratio, const std::array<double, 3>& rCenter, int Order,
                     std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "spherical_shell";
    prim_check_count(shape, "n", N);
    prim_check_count(shape, "nr", Nr);
    prim_check_order(shape, Order);
    if (N == 0 || Nr == 0)
        return prim_empty();
    prim_check_positive(shape, "inner_radius", InnerRadius);
    if (!(OuterRadius > InnerRadius))
        prim_fail(shape, "outer_radius must exceed inner_radius, got " +
                             std::to_string(OuterRadius) + " and " + std::to_string(InnerRadius));
    prim_check_positive(shape, "ratio", Ratio);
    const PrimCubeLattice lat{N, N, N};
    prim_check_budget(shape, prim_mul(lat.NumFaces(), Nr), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Shell;
    soup.mRadius = InnerRadius;
    soup.mRadius2 = OuterRadius;
    soup.mRatio = Ratio;
    soup.mLayers = static_cast<double>(Nr);
    soup.mCenter = rCenter;
    const double n = static_cast<double>(N);
    for (std::int64_t k = 0; k <= Nr; ++k) {
        const double t = static_cast<double>(Nr + k) / static_cast<double>(Nr);
        const std::uint16_t mask = (k == 0 ? 1 : 0) | (k == Nr ? 2 : 0);
        prim_cube_nodes(lat, [&](std::int64_t I, std::int64_t J, std::int64_t K) {
            soup.AddNode((static_cast<double>(2 * I - N) / n) * t,
                         (static_cast<double>(2 * J - N) / n) * t,
                         (static_cast<double>(2 * K - N) / n) * t, mask);
        });
    }
    const std::int64_t layer = lat.Count();
    PrimBlock& blk = soup.Block(CellType::Hexahedron, 8);
    for (std::int64_t k = 0; k < Nr; ++k) {
        const std::int64_t lo = k * layer, hi = lo + layer;
        lat.Faces([&](int, std::int64_t A, std::int64_t B, std::int64_t C, std::int64_t D) {
            prim_push(blk, {A + lo, B + lo, C + lo, D + lo, A + hi, B + hi, C + hi, D + hi});
        });
    }
    return prim_finish(shape, soup, Order, Complete, {{"inner", 1}, {"outer", 2}}, 2);
}

Mesh pipe(double InnerRadius, double OuterRadius, double Height, std::int64_t Ntheta,
          std::int64_t Nr, std::int64_t Nz, double Ratio, const std::array<double, 3>& rCenter,
          int Order, std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "pipe";
    prim_check_count(shape, "ntheta", Ntheta);
    prim_check_count(shape, "nr", Nr);
    prim_check_count(shape, "nz", Nz);
    prim_check_order(shape, Order);
    if (Ntheta == 0 || Nr == 0 || Nz == 0)
        return prim_empty();
    prim_check_positive(shape, "inner_radius", InnerRadius);
    if (!(OuterRadius > InnerRadius))
        prim_fail(shape, "outer_radius must exceed inner_radius, got " +
                             std::to_string(OuterRadius) + " and " + std::to_string(InnerRadius));
    prim_check_positive(shape, "height", Height);
    prim_check_positive(shape, "ratio", Ratio);
    if (Ntheta < 3)
        prim_fail(shape, "ntheta must be at least 3, got " + std::to_string(Ntheta));
    prim_check_budget(shape, prim_mul(prim_mul(Ntheta, Nr), Nz), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::Graded;
    soup.mRadius = InnerRadius;
    soup.mRadius2 = OuterRadius;
    soup.mRatio = Ratio;
    soup.mLayers = static_cast<double>(Nr);
    soup.mCenter = rCenter;
    for (std::int64_t k = 0; k <= Nz; ++k) {
        const double z = (Height * static_cast<double>(2 * k - Nz)) / static_cast<double>(2 * Nz);
        for (std::int64_t m = 0; m <= Nr; ++m) {
            const std::uint16_t mask =
                (m == 0 ? 1 : 0) | (m == Nr ? 2 : 0) | (k == 0 ? 4 : 0) | (k == Nz ? 8 : 0);
            for (std::int64_t i = 0; i < Ntheta; ++i)
                soup.AddNode(static_cast<double>(m) / static_cast<double>(Nr),
                             (2.0 * kPrimPi * static_cast<double>(i)) / static_cast<double>(Ntheta),
                             z, mask);
        }
    }
    auto id = [&](std::int64_t K, std::int64_t M, std::int64_t I) {
        return (K * (Nr + 1) + M) * Ntheta + (I % Ntheta);
    };
    PrimBlock& blk = soup.Block(CellType::Hexahedron, 8);
    for (std::int64_t k = 0; k < Nz; ++k)
        for (std::int64_t m = 0; m < Nr; ++m)
            for (std::int64_t i = 0; i < Ntheta; ++i)
                prim_push(blk, {id(k, m, i), id(k, m + 1, i), id(k, m + 1, i + 1), id(k, m, i + 1),
                                id(k + 1, m, i), id(k + 1, m + 1, i), id(k + 1, m + 1, i + 1),
                                id(k + 1, m, i + 1)});
    return prim_finish(shape, soup, Order, Complete,
                       {{"inner", 1}, {"outer", 2}, {"bottom", 4}, {"top", 8}}, 2);
}

Mesh fichera(std::int64_t N, double Size, const std::array<double, 3>& rCenter,
             const std::string& rCellType, int Order, std::int64_t MaxCells, bool Complete) {
    constexpr const char* shape = "fichera";
    prim_check_count(shape, "n", N);
    prim_check_order(shape, Order);
    const bool tet = rCellType == "tetra";
    if (!tet && rCellType != "hexahedron")
        prim_bad_choice(shape, "cell_type", rCellType, "hexahedron, tetra");
    if (N == 0)
        return prim_empty();
    prim_check_positive(shape, "size", Size);
    prim_check_budget(shape, prim_mul(prim_mul(7, prim_mul(prim_mul(N, N), N)), tet ? 6 : 1),
                      MaxCells);

    PrimSoup soup;
    const std::int64_t w = 2 * N + 1, top = 2 * N;
    // Nodes strictly inside the removed octant (x, y, z > 0) are skipped; the
    // lattice keeps the notch walls.
    std::vector<std::int64_t> ids(static_cast<std::size_t>(w * w * w), -1);
    for (std::int64_t k = 0; k < w; ++k)
        for (std::int64_t j = 0; j < w; ++j)
            for (std::int64_t i = 0; i < w; ++i) {
                if (i > N && j > N && k > N)
                    continue;
                std::uint16_t mask = 0;
                if (i == 0)
                    mask |= 1;
                if (i == top)
                    mask |= 2;
                if (j == 0)
                    mask |= 4;
                if (j == top)
                    mask |= 8;
                if (k == 0)
                    mask |= 16;
                if (k == top)
                    mask |= 32;
                if (i == N && j >= N && k >= N)
                    mask |= 64;
                if (j == N && i >= N && k >= N)
                    mask |= 128;
                if (k == N && i >= N && j >= N)
                    mask |= 256;
                ids[static_cast<std::size_t>((k * w + j) * w + i)] = soup.AddNode(
                    rCenter[0] + (Size * static_cast<double>(i - N)) / static_cast<double>(2 * N),
                    rCenter[1] + (Size * static_cast<double>(j - N)) / static_cast<double>(2 * N),
                    rCenter[2] + (Size * static_cast<double>(k - N)) / static_cast<double>(2 * N),
                    mask);
            }
    auto at = [&](std::int64_t I, std::int64_t J, std::int64_t K) {
        return ids[static_cast<std::size_t>((K * w + J) * w + I)];
    };
    PrimBlock& blk = soup.Block(tet ? CellType::Tetra : CellType::Hexahedron, tet ? 4 : 8);
    for (std::int64_t k = 0; k < top; ++k)
        for (std::int64_t j = 0; j < top; ++j)
            for (std::int64_t i = 0; i < top; ++i) {
                if (i >= N && j >= N && k >= N)
                    continue;
                const std::int64_t h[8] = {at(i, j, k),
                                           at(i + 1, j, k),
                                           at(i + 1, j + 1, k),
                                           at(i, j + 1, k),
                                           at(i, j, k + 1),
                                           at(i + 1, j, k + 1),
                                           at(i + 1, j + 1, k + 1),
                                           at(i, j + 1, k + 1)};
                prim_push_hex(blk, tet, h);
            }
    return prim_finish(shape, soup, Order, Complete,
                       {{"xmin", 1},
                        {"xmax", 2},
                        {"ymin", 4},
                        {"ymax", 8},
                        {"zmin", 16},
                        {"zmax", 32},
                        {"notch_x", 64},
                        {"notch_y", 128},
                        {"notch_z", 256}},
                       2);
}

Mesh solid_torus(double MajorRadius, double MinorRadius, std::int64_t N, std::int64_t Nu,
                 const std::array<double, 3>& rCenter, int Order, std::int64_t MaxCells,
                 bool Complete) {
    constexpr const char* shape = "solid_torus";
    prim_check_positive(shape, "major_radius", MajorRadius);
    prim_check_positive(shape, "minor_radius", MinorRadius);
    if (!(MinorRadius < MajorRadius))
        prim_fail(shape, "minor_radius must be less than major_radius, got " +
                             std::to_string(MinorRadius) + " and " + std::to_string(MajorRadius));
    prim_check_count(shape, "n", N);
    prim_check_count(shape, "nu", Nu);
    prim_check_order(shape, Order);
    if (N == 0 || Nu == 0)
        return prim_empty();
    if (Nu < 3)
        prim_fail(shape, "nu must be at least 3, got " + std::to_string(Nu));
    prim_check_budget(shape, prim_mul(prim_mul(N, N), Nu), MaxCells);

    PrimSoup soup;
    soup.mMap = PrimMap::SolidTorus;
    soup.mRadius = MajorRadius;
    soup.mRadius2 = MinorRadius;
    soup.mCenter = rCenter;
    const std::int64_t p = N + 1, layer = p * p;
    for (std::int64_t i = 0; i < Nu; ++i) {
        const double u = (2.0 * kPrimPi * static_cast<double>(i)) / static_cast<double>(Nu);
        for (std::int64_t jj = 0; jj <= N; ++jj)
            for (std::int64_t ii = 0; ii <= N; ++ii) {
                const bool rim = ii == 0 || ii == N || jj == 0 || jj == N;
                soup.AddNode(static_cast<double>(2 * ii - N) / static_cast<double>(N),
                             static_cast<double>(2 * jj - N) / static_cast<double>(N), u,
                             rim ? 4 : 0);
            }
    }
    PrimBlock& blk = soup.Block(CellType::Hexahedron, 8);
    for (std::int64_t i = 0; i < Nu; ++i) {
        const std::int64_t i1 = (i + 1) % Nu;
        for (std::int64_t jj = 0; jj < N; ++jj)
            for (std::int64_t ii = 0; ii < N; ++ii) {
                // The cross-section quad runs clockwise in (s, t), so the hexahedron
                // swept counter-clockwise around the axis has a positive volume.
                const std::int64_t a = i * layer + jj * p + ii;
                const std::int64_t b = i1 * layer + jj * p + ii;
                prim_push(blk, {a, a + p, a + p + 1, a + 1, b, b + p, b + p + 1, b + 1});
            }
    }
    return prim_finish(shape, soup, Order, Complete, {{"surface", 4}}, 2);
}

// ===========================================================================
// The shorter signatures: Complete = false.
// ===========================================================================

Mesh rectangle(const std::array<std::int64_t, 2>& rDims, const std::array<double, 2>& rSize,
               const std::array<double, 3>& rOrigin, const std::string& rCellType,
               const std::string& rSplit, int Order, std::int64_t MaxCells) {
    return rectangle(rDims, rSize, rOrigin, rCellType, rSplit, Order, MaxCells, false);
}

Mesh box(const std::array<std::int64_t, 3>& rDims, const std::array<double, 3>& rSize,
         const std::array<double, 3>& rOrigin, const std::string& rCellType, int Order,
         std::int64_t MaxCells) {
    return box(rDims, rSize, rOrigin, rCellType, Order, MaxCells, false);
}

Mesh disk(double Radius, std::int64_t N, const std::array<double, 3>& rCenter,
          const std::string& rCellType, int Order, std::int64_t MaxCells) {
    return disk(Radius, N, rCenter, rCellType, Order, MaxCells, false);
}

Mesh uv_sphere(double Radius, std::int64_t Nu, std::int64_t Nv,
               const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
               std::int64_t MaxCells) {
    return uv_sphere(Radius, Nu, Nv, rCenter, rCellType, Order, MaxCells, false);
}

Mesh icosa_sphere(double Radius, std::int64_t Subdivisions, const std::array<double, 3>& rCenter,
                  int Order, std::int64_t MaxCells) {
    return icosa_sphere(Radius, Subdivisions, rCenter, Order, MaxCells, false);
}

Mesh cylinder(double Radius, double Height, std::int64_t Nxy, std::int64_t Nz,
              const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
              std::int64_t MaxCells) {
    return cylinder(Radius, Height, Nxy, Nz, rCenter, rCellType, Order, MaxCells, false);
}

Mesh reference_triangle(std::int64_t N, const std::array<double, 2>& rSize,
                        const std::array<double, 3>& rOrigin, int Order, std::int64_t MaxCells) {
    return reference_triangle(N, rSize, rOrigin, Order, MaxCells, false);
}

Mesh ngon(std::int64_t Sides, double Radius, std::int64_t Rings,
          const std::array<double, 3>& rCenter, int Order, std::int64_t MaxCells) {
    return ngon(Sides, Radius, Rings, rCenter, Order, MaxCells, false);
}

Mesh annulus(double InnerRadius, double OuterRadius, std::int64_t Ntheta, std::int64_t Nr,
             double Angle, const std::array<double, 3>& rCenter, const std::string& rCellType,
             int Order, std::int64_t MaxCells) {
    return annulus(InnerRadius, OuterRadius, Ntheta, Nr, Angle, rCenter, rCellType, Order, MaxCells,
                   false);
}

Mesh lshape(std::int64_t N, double Size, const std::array<double, 3>& rCenter,
            const std::string& rCellType, int Order, std::int64_t MaxCells) {
    return lshape(N, Size, rCenter, rCellType, Order, MaxCells, false);
}

Mesh plate_with_hole(double Radius, double Size, std::int64_t Nt, std::int64_t Nr,
                     const std::array<double, 3>& rCenter, const std::string& rCellType, int Order,
                     std::int64_t MaxCells) {
    return plate_with_hole(Radius, Size, Nt, Nr, rCenter, rCellType, Order, MaxCells, false);
}

}  // namespace meshioplusplus
