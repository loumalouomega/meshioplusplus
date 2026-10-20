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
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

using meshioplusplus::annulus;
using meshioplusplus::box;
using meshioplusplus::cylinder;
using meshioplusplus::disk;
using meshioplusplus::honeycomb;
using meshioplusplus::icosa_sphere;
using meshioplusplus::lshape;
using meshioplusplus::Mesh;
using meshioplusplus::ngon;
using meshioplusplus::plate_with_hole;
using meshioplusplus::rectangle;
using meshioplusplus::reference_cells;
using meshioplusplus::reference_triangle;
using meshioplusplus::uv_sphere;

namespace {

using Vec = std::array<double, 3>;

Vec point(const Mesh& rMesh, std::int64_t I) {
    const double* p = rMesh.Points().As<double>();
    return {p[3 * I], p[3 * I + 1], p[3 * I + 2]};
}

Vec sub(const Vec& A, const Vec& B) {
    return {A[0] - B[0], A[1] - B[1], A[2] - B[2]};
}

Vec cross(const Vec& A, const Vec& B) {
    return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
}

double dot(const Vec& A, const Vec& B) {
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}

double norm(const Vec& A) {
    return std::sqrt(dot(A, A));
}

double tet_volume(const Mesh& rMesh, const std::int64_t* pRow, const int (&rIdx)[4]) {
    const Vec a = point(rMesh, pRow[rIdx[0]]);
    return dot(cross(sub(point(rMesh, pRow[rIdx[1]]), a), sub(point(rMesh, pRow[rIdx[2]]), a)),
               sub(point(rMesh, pRow[rIdx[3]]), a)) /
           6.0;
}

// Signed volume of every cell of block B, by the Kuhn / prism decompositions of
// the corner nodes; the same oracle for the linear and the quadratic types.
std::vector<double> volumes(const Mesh& rMesh, std::size_t B) {
    const auto view = rMesh.Cells(B);
    const std::string type = view.Type();
    const std::int64_t* conn = view.Conn().As<std::int64_t>();
    const std::size_t w = view.NodesPerCell();
    static const int kHex[6][4] = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6},
                                   {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    static const int kWedge[3][4] = {{0, 1, 2, 3}, {1, 2, 3, 4}, {2, 3, 4, 5}};
    static const int kTet[1][4] = {{0, 1, 2, 3}};
    std::vector<double> out;
    for (std::size_t c = 0; c < view.NumCells(); ++c) {
        const std::int64_t* row = conn + c * w;
        double v = 0.0;
        if (type.rfind("hexahedron", 0) == 0)
            for (const auto& t : kHex)
                v += tet_volume(rMesh, row, t);
        else if (type.rfind("wedge", 0) == 0)
            for (const auto& t : kWedge)
                v += tet_volume(rMesh, row, t);
        else
            v = tet_volume(rMesh, row, kTet[0]);
        out.push_back(v);
    }
    return out;
}

// Twice the signed area, seen from +z, of every triangle or quad of block B.
std::vector<double> areas(const Mesh& rMesh, std::size_t B) {
    const auto view = rMesh.Cells(B);
    const std::int64_t* conn = view.Conn().As<std::int64_t>();
    const std::size_t w = view.NodesPerCell();
    const bool tri = view.Type().rfind("triangle", 0) == 0;
    std::vector<double> out;
    for (std::size_t c = 0; c < view.NumCells(); ++c) {
        const std::int64_t* row = conn + c * w;
        if (tri)
            out.push_back(cross(sub(point(rMesh, row[1]), point(rMesh, row[0])),
                                sub(point(rMesh, row[2]), point(rMesh, row[0])))[2]);
        else
            out.push_back(cross(sub(point(rMesh, row[2]), point(rMesh, row[0])),
                                sub(point(rMesh, row[3]), point(rMesh, row[1])))[2]);
    }
    return out;
}

const meshioplusplus::Region& region(const Mesh& rMesh, const std::string& rName) {
    const std::size_t i = rMesh.FindRegion(rName, meshioplusplus::RegionKind::Side);
    EXPECT_NE(i, Mesh::npos) << rName;
    return rMesh.Region(i);
}

// Surface cells whose normal points toward the origin.
std::size_t inward_cells(const Mesh& rMesh) {
    std::size_t bad = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto view = rMesh.Cells(b);
        const std::int64_t* conn = view.Conn().As<std::int64_t>();
        const std::size_t w = view.NodesPerCell();
        const bool tri = view.Type().rfind("triangle", 0) == 0;
        for (std::size_t c = 0; c < view.NumCells(); ++c) {
            const std::int64_t* row = conn + c * w;
            const Vec a = point(rMesh, row[0]), b1 = point(rMesh, row[1]),
                      c1 = point(rMesh, row[2]);
            Vec n, cen;
            if (tri) {
                n = cross(sub(b1, a), sub(c1, a));
                cen = {(a[0] + b1[0] + c1[0]) / 3, (a[1] + b1[1] + c1[1]) / 3,
                       (a[2] + b1[2] + c1[2]) / 3};
            } else {
                const Vec d = point(rMesh, row[3]);
                n = cross(sub(c1, a), sub(d, b1));
                cen = {(a[0] + b1[0] + c1[0] + d[0]) / 4, (a[1] + b1[1] + c1[1] + d[1]) / 4,
                       (a[2] + b1[2] + c1[2] + d[2]) / 4};
            }
            if (dot(n, cen) <= 0.0)
                ++bad;
        }
    }
    return bad;
}

std::size_t distinct_points(const Mesh& rMesh) {
    std::set<std::array<long long, 3>> seen;
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i) {
        const Vec p = point(rMesh, static_cast<std::int64_t>(i));
        seen.insert({std::llround(p[0] * 1e9), std::llround(p[1] * 1e9), std::llround(p[2] * 1e9)});
    }
    return seen.size();
}

bool same_bytes(const Mesh& rA, const Mesh& rB) {
    if (rA.NumPoints() != rB.NumPoints() || rA.NumCellBlocks() != rB.NumCellBlocks())
        return false;
    const double* a = rA.Points().As<double>();
    const double* b = rB.Points().As<double>();
    for (std::size_t i = 0; i < rA.NumPoints() * 3; ++i)
        if (a[i] != b[i])
            return false;
    for (std::size_t k = 0; k < rA.NumCellBlocks(); ++k) {
        if (rA.Cells(k).Type() != rB.Cells(k).Type() || rA.Cells(k).NumCells() != rB.Cells(k).NumCells())
            return false;
        const std::size_t n = rA.Cells(k).NumCells() * rA.Cells(k).NodesPerCell();
        const std::int64_t* ca = rA.Cells(k).Conn().As<std::int64_t>();
        const std::int64_t* cb = rB.Cells(k).Conn().As<std::int64_t>();
        for (std::size_t i = 0; i < n; ++i)
            if (ca[i] != cb[i])
                return false;
    }
    return true;
}

}  // namespace

// --- the shared contract ---------------------------------------------------

TEST(Primitives, ZeroCountsGiveAnEmptyMesh) {
    for (const Mesh& m :
         {rectangle({{0, 3}}), box({{2, 0, 2}}), disk(1.0, 0), uv_sphere(1.0, 0, 4),
          cylinder(1.0, 1.0, 4, 0), reference_triangle(0), ngon(0), ngon(5, 1.0, 0),
          annulus(0.5, 1.0, 0), annulus(0.5, 1.0, 8, 0), lshape(0), plate_with_hole(0.25, 1.0, 0),
          plate_with_hole(0.25, 1.0, 4, 0), honeycomb({{0, 3}}), honeycomb({{3, 0}})}) {
        EXPECT_EQ(m.NumPoints(), 0u);
        EXPECT_EQ(m.NumCellBlocks(), 0u);
    }
}

TEST(Primitives, BadArgumentsRaiseByName) {
    EXPECT_THROW(rectangle({{-1, 2}}), std::invalid_argument);
    EXPECT_THROW(rectangle({{2, 2}}, {{0.0, 1.0}}), std::invalid_argument);
    EXPECT_THROW(rectangle({{2, 2}}, {{1.0, 1.0}}, {{0, 0, 0}}, "hexahedron"), std::invalid_argument);
    EXPECT_THROW(rectangle({{2, 2}}, {{1.0, 1.0}}, {{0, 0, 0}}, "quad", "diag"), std::invalid_argument);
    EXPECT_THROW(rectangle({{2, 2}}, {{1.0, 1.0}}, {{0, 0, 0}}, "quad", "zigzag", 3),
                 std::invalid_argument);
    EXPECT_THROW(box({{2, 2, 2}}, {{1, 1, 1}}, {{0, 0, 0}}, "quad"), std::invalid_argument);
    EXPECT_THROW(disk(-1.0), std::invalid_argument);
    EXPECT_THROW(uv_sphere(1.0, 2, 4), std::invalid_argument);
    EXPECT_THROW(uv_sphere(1.0, 8, 1), std::invalid_argument);
    EXPECT_THROW(icosa_sphere(1.0, -1), std::invalid_argument);
    EXPECT_THROW(cylinder(1.0, 0.0), std::invalid_argument);
    EXPECT_THROW(reference_triangle(-1), std::invalid_argument);
    EXPECT_THROW(reference_triangle(3, {{1.0, 0.0}}), std::invalid_argument);
    EXPECT_THROW(ngon(2), std::invalid_argument);
    EXPECT_THROW(ngon(5, -1.0), std::invalid_argument);
    EXPECT_THROW(annulus(1.0, 0.5), std::invalid_argument);
    EXPECT_THROW(annulus(0.0, 1.0), std::invalid_argument);
    EXPECT_THROW(annulus(0.5, 1.0, 2), std::invalid_argument);
    EXPECT_THROW(annulus(0.5, 1.0, 8, 2, 7.0), std::invalid_argument);
    EXPECT_THROW(annulus(0.5, 1.0, 1, 2, 4.0), std::invalid_argument);
    EXPECT_THROW(annulus(0.5, 1.0, 8, 2, 1.0, {{0, 0, 0}}, "hexahedron"), std::invalid_argument);
    EXPECT_THROW(lshape(2, 0.0), std::invalid_argument);
    EXPECT_THROW(lshape(2, 2.0, {{0, 0, 0}}, "wedge"), std::invalid_argument);
    EXPECT_THROW(plate_with_hole(0.6, 1.0), std::invalid_argument);
    EXPECT_THROW(plate_with_hole(0.25, 1.0, 4, 2, {{0, 0, 0}}, "quad", 3), std::invalid_argument);
    EXPECT_THROW(honeycomb({{2, -1}}), std::invalid_argument);
    EXPECT_THROW(honeycomb({{2, 2}}, 0.0), std::invalid_argument);
    EXPECT_THROW(reference_cells({"blob"}), std::invalid_argument);
    EXPECT_THROW(reference_cells({"quad", "quad"}), std::invalid_argument);
    EXPECT_THROW(reference_cells({"tetra10"}, true), std::invalid_argument);
    try {
        rectangle({{2, 2}}, {{1.0, 1.0}}, {{0, 0, 0}}, "hex");
        FAIL();
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("meshio++: rectangle: unknown cell_type 'hex'"),
                  std::string::npos);
    }
}

TEST(Primitives, MaxCellsRefusesByName) {
    EXPECT_THROW(rectangle({{100, 100}}, {{1, 1}}, {{0, 0, 0}}, "quad", "zigzag", 1, 50),
                 std::invalid_argument);
    EXPECT_THROW(box({{10, 10, 10}}, {{1, 1, 1}}, {{0, 0, 0}}, "hexahedron", 1, 50),
                 std::invalid_argument);
    EXPECT_THROW(disk(1.0, 20, {{0, 0, 0}}, "quad", 1, 50), std::invalid_argument);
    EXPECT_THROW(uv_sphere(1.0, 64, 64, {{0, 0, 0}}, "triangle", 1, 50), std::invalid_argument);
    EXPECT_THROW(icosa_sphere(1.0, 6, {{0, 0, 0}}, 1, 50), std::invalid_argument);
    EXPECT_THROW(cylinder(1.0, 1.0, 10, 10, {{0, 0, 0}}, "hexahedron", 1, 50),
                 std::invalid_argument);
    EXPECT_THROW(reference_triangle(20, {{1, 1}}, {{0, 0, 0}}, 1, 50), std::invalid_argument);
    EXPECT_THROW(ngon(6, 1.0, 10, {{0, 0, 0}}, 1, 50), std::invalid_argument);
    EXPECT_THROW(annulus(0.5, 1.0, 20, 20, 6.283185307179586, {{0, 0, 0}}, "quad", 1, 50),
                 std::invalid_argument);
    EXPECT_THROW(lshape(10, 2.0, {{0, 0, 0}}, "quad", 1, 50), std::invalid_argument);
    EXPECT_THROW(plate_with_hole(0.25, 1.0, 8, 8, {{0, 0, 0}}, "quad", 1, 50),
                 std::invalid_argument);
    EXPECT_THROW(honeycomb({{10, 10}}, 1.0, {{0, 0, 0}}, 50), std::invalid_argument);
    EXPECT_NO_THROW(rectangle({{8, 8}}, {{1, 1}}, {{0, 0, 0}}, "quad", "zigzag", 1, 0));
}

TEST(Primitives, OutputIsDeterministic) {
    EXPECT_TRUE(same_bytes(uv_sphere(1.0, 9, 5, {{0, 0, 0}}, "quad", 2),
                           uv_sphere(1.0, 9, 5, {{0, 0, 0}}, "quad", 2)));
    EXPECT_TRUE(same_bytes(icosa_sphere(1.0, 3, {{0, 0, 0}}, 2), icosa_sphere(1.0, 3, {{0, 0, 0}}, 2)));
}

// --- rectangle and box -----------------------------------------------------

TEST(Primitives, BoxHexahedraMatchGridNumbering) {
    const Mesh b = box({{3, 2, 4}}, {{3.0, 2.0, 4.0}}, {{-1.0, 0.0, 2.0}});
    const Mesh g = meshioplusplus::grid({{3, 2, 4}}, {{-1.0, 0.0, 2.0}}, {{1.0, 1.0, 1.0}});
    ASSERT_EQ(b.NumPoints(), g.NumPoints());
    const double* pb = b.Points().As<double>();
    const double* pg = g.Points().As<double>();
    for (std::size_t i = 0; i < b.NumPoints() * 3; ++i)
        ASSERT_DOUBLE_EQ(pb[i], pg[i]);
    ASSERT_EQ(b.Cells(0).NumCells(), g.Cells(0).NumCells());
    const std::int64_t* cb = b.Cells(0).Conn().As<std::int64_t>();
    const std::int64_t* cg = g.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t i = 0; i < g.Cells(0).NumCells() * 8; ++i)
        ASSERT_EQ(cb[i], cg[i]);
}

TEST(Primitives, BoxCellTypesFillTheVolumePositively) {
    for (const auto& [type, per] : std::map<std::string, std::size_t>{
             {"hexahedron", 1}, {"tetra", 6}, {"wedge", 2}}) {
        const Mesh m = box({{2, 3, 2}}, {{1.0, 1.0, 1.0}}, {{0, 0, 0}}, type);
        ASSERT_EQ(m.NumCellBlocks(), 1u);
        EXPECT_EQ(m.Cells(0).NumCells(), 12 * per) << type;
        double total = 0.0;
        for (double v : volumes(m, 0)) {
            EXPECT_GT(v, 0.0) << type;
            total += v;
        }
        EXPECT_NEAR(total, 1.0, 1e-12) << type;
    }
}

TEST(Primitives, BoxRegionsPartitionTheBoundary) {
    for (const std::string type : {"hexahedron", "tetra", "wedge"}) {
        const Mesh m = box({{2, 2, 2}}, {{1, 1, 1}}, {{0, 0, 0}}, type);
        std::size_t total = 0;
        for (const char* name : {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}) {
            const auto& r = region(m, name);
            EXPECT_EQ(r.mDim, 2);
            total += r.NumEntries();
        }
        // 6 faces x 4 quads of the hexahedron box; two triangles per quad for
        // tetrahedra; a wedge keeps one quad per side face and two triangles
        // per top and bottom cell.
        const std::size_t want = type == "hexahedron" ? 24 : (type == "tetra" ? 48 : 32);
        EXPECT_EQ(total, want) << type;
    }
    const Mesh m = box({{2, 2, 2}});
    const auto& xmin = region(m, "xmin");
    ASSERT_EQ(xmin.NumEntries(), 4u);
    for (std::size_t i = 0; i < 4; ++i)
        EXPECT_EQ(xmin.mEntries.As<std::int64_t>()[2 * i + 1], 0);  // hexahedron facet 0
}

TEST(Primitives, RectangleSplitsAreAllPositiveAndFillTheArea) {
    for (const std::string split : {"zigzag", "up", "down", "center"}) {
        const Mesh m = rectangle({{4, 3}}, {{2.0, 1.5}}, {{0, 0, 0}}, "triangle", split);
        EXPECT_EQ(m.Cells(0).NumCells(), split == "center" ? 48u : 24u) << split;
        double total = 0.0;
        for (double a : areas(m, 0)) {
            EXPECT_GT(a, 0.0) << split;
            total += a / 2;
        }
        EXPECT_NEAR(total, 3.0, 1e-12) << split;
        EXPECT_EQ(region(m, "xmin").NumEntries(), 3u) << split;
        EXPECT_EQ(region(m, "ymax").NumEntries(), 4u) << split;
    }
}

TEST(Primitives, RectangleCoordinatesComeFromIndices) {
    const Mesh m = rectangle({{7, 3}}, {{0.7, 0.3}}, {{0.1, 0.2, 5.0}});
    const double* p = m.Points().As<double>();
    for (std::size_t j = 0; j <= 3; ++j)
        for (std::size_t i = 0; i <= 7; ++i) {
            const std::size_t v = j * 8 + i;
            EXPECT_EQ(p[3 * v], 0.1 + (0.7 * static_cast<double>(i)) / 7.0);
            EXPECT_EQ(p[3 * v + 1], 0.2 + (0.3 * static_cast<double>(j)) / 3.0);
            EXPECT_EQ(p[3 * v + 2], 5.0);
        }
    // The far edge is exactly origin + size.
    EXPECT_EQ(p[3 * 7], 0.1 + 0.7);
}

// --- disk and cylinder -----------------------------------------------------

TEST(Primitives, DiskRimLiesOnTheCircle) {
    for (const std::string type : {"quad", "triangle"}) {
        const Mesh m = disk(2.0, 4, {{1.0, -1.0, 0.5}}, type);
        for (double a : areas(m, 0))
            EXPECT_GT(a, 0.0) << type;
        const auto& rim = region(m, "rim");
        EXPECT_EQ(rim.mDim, 1);
        EXPECT_EQ(rim.NumEntries(), type == "quad" ? 16u : 24u);
        const std::int64_t* conn = m.Cells(0).Conn().As<std::int64_t>();
        const std::size_t w = m.Cells(0).NodesPerCell();
        for (std::size_t i = 0; i < rim.NumEntries(); ++i) {
            const std::int64_t c = rim.mEntries.As<std::int64_t>()[2 * i];
            const std::int64_t f = rim.mEntries.As<std::int64_t>()[2 * i + 1];
            for (std::int64_t k : {f, (f + 1) % static_cast<std::int64_t>(w)}) {
                const Vec p = point(m, conn[static_cast<std::size_t>(c) * w + static_cast<std::size_t>(k)]);
                EXPECT_NEAR(std::hypot(p[0] - 1.0, p[1] + 1.0), 2.0, 1e-12) << type;
                EXPECT_EQ(p[2], 0.5);
            }
        }
    }
}

TEST(Primitives, DiskTriangleCountIsSixNSquared) {
    for (std::int64_t n = 1; n <= 5; ++n) {
        const Mesh m = disk(1.0, n, {{0, 0, 0}}, "triangle");
        EXPECT_EQ(m.Cells(0).NumCells(), static_cast<std::size_t>(6 * n * n));
        EXPECT_EQ(m.NumPoints(), static_cast<std::size_t>(1 + 3 * n * (n + 1)));
    }
}

TEST(Primitives, CylinderCellsAndRegions) {
    for (const std::string type : {"hexahedron", "wedge"}) {
        const Mesh m = cylinder(1.0, 4.0, 3, 2, {{0, 0, 10.0}}, type);
        for (double v : volumes(m, 0))
            EXPECT_GT(v, 0.0) << type;
        const std::size_t side = region(m, "side").NumEntries();
        const std::size_t top = region(m, "top").NumEntries();
        EXPECT_EQ(top, region(m, "bottom").NumEntries());
        EXPECT_EQ(side, type == "hexahedron" ? 24u : 36u);
        EXPECT_EQ(top, type == "hexahedron" ? 9u : 54u);
        double zmin = 1e300, zmax = -1e300;
        for (std::size_t i = 0; i < m.NumPoints(); ++i) {
            const double z = point(m, static_cast<std::int64_t>(i))[2];
            zmin = std::min(zmin, z);
            zmax = std::max(zmax, z);
        }
        EXPECT_EQ(zmin, 8.0);
        EXPECT_EQ(zmax, 12.0);
    }
    // One cell: the caps and the side stay separate regions.
    const Mesh one = cylinder(1.0, 1.0, 1, 1);
    EXPECT_EQ(region(one, "side").NumEntries(), 4u);
    EXPECT_EQ(region(one, "top").NumEntries(), 1u);
    EXPECT_EQ(region(one, "bottom").NumEntries(), 1u);
}

// --- spheres ---------------------------------------------------------------

TEST(Primitives, UvSphereIsOutwardAndWelded) {
    for (const std::string type : {"triangle", "quad"}) {
        const Mesh m = uv_sphere(2.0, 12, 6, {{0, 0, 0}}, type);
        EXPECT_EQ(m.NumPoints(), static_cast<std::size_t>(2 + 12 * 5));
        EXPECT_EQ(distinct_points(m), m.NumPoints());
        EXPECT_EQ(inward_cells(m), 0u) << type;
        for (std::size_t i = 0; i < m.NumPoints(); ++i)
            EXPECT_NEAR(norm(point(m, static_cast<std::int64_t>(i))), 2.0, 1e-14);
        EXPECT_EQ(point(m, 0), (Vec{0.0, 0.0, 2.0}));
        EXPECT_EQ(point(m, static_cast<std::int64_t>(m.NumPoints()) - 1), (Vec{0.0, 0.0, -2.0}));
    }
}

TEST(Primitives, IcosaSphereCountsOrientationAndClosure) {
    for (std::int64_t s = 0; s <= 3; ++s) {
        const Mesh m = icosa_sphere(3.0, s);
        std::int64_t faces = 20;
        for (std::int64_t i = 0; i < s; ++i)
            faces *= 4;
        EXPECT_EQ(m.Cells(0).NumCells(), static_cast<std::size_t>(faces));
        EXPECT_EQ(m.NumPoints(), static_cast<std::size_t>(faces / 2 + 2));
        EXPECT_EQ(inward_cells(m), 0u);
        for (std::size_t i = 0; i < m.NumPoints(); ++i)
            EXPECT_NEAR(norm(point(m, static_cast<std::int64_t>(i))), 3.0, 1e-14);
        // Closed: every undirected edge is shared by exactly two triangles.
        std::map<std::pair<std::int64_t, std::int64_t>, int> edges;
        const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
        for (std::int64_t f = 0; f < faces; ++f)
            for (int e = 0; e < 3; ++e) {
                std::int64_t a = c[3 * f + e], b = c[3 * f + (e + 1) % 3];
                if (a > b)
                    std::swap(a, b);
                ++edges[{a, b}];
            }
        for (const auto& kv : edges)
            EXPECT_EQ(kv.second, 2);
    }
}

// --- order 2 ---------------------------------------------------------------

TEST(Primitives, QuadraticTypesAreNamedAndSharedBetweenNeighbours) {
    struct Case {
        Mesh mMesh;
        const char* mType;
    };
    const Case cases[] = {
        {rectangle({{3, 2}}, {{1, 1}}, {{0, 0, 0}}, "quad", "zigzag", 2), "quad8"},
        {rectangle({{3, 2}}, {{1, 1}}, {{0, 0, 0}}, "triangle", "center", 2), "triangle6"},
        {box({{2, 2, 2}}, {{1, 1, 1}}, {{0, 0, 0}}, "hexahedron", 2), "hexahedron20"},
        {box({{2, 2, 2}}, {{1, 1, 1}}, {{0, 0, 0}}, "tetra", 2), "tetra10"},
        {box({{2, 2, 2}}, {{1, 1, 1}}, {{0, 0, 0}}, "wedge", 2), "wedge15"},
    };
    for (const auto& c : cases) {
        ASSERT_EQ(c.mMesh.NumCellBlocks(), 1u) << c.mType;
        EXPECT_EQ(c.mMesh.Cells(0).Type(), c.mType);
        EXPECT_EQ(distinct_points(c.mMesh), c.mMesh.NumPoints()) << c.mType;
    }
}

TEST(Primitives, CurvedShapesPlaceMidNodesOnTheExactGeometry) {
    for (const Mesh& m : {uv_sphere(3.0, 7, 5, {{0, 0, 0}}, "triangle", 2),
                          uv_sphere(3.0, 7, 5, {{0, 0, 0}}, "quad", 2),
                          icosa_sphere(3.0, 2, {{0, 0, 0}}, 2)}) {
        EXPECT_EQ(distinct_points(m), m.NumPoints());
        for (std::size_t i = 0; i < m.NumPoints(); ++i)
            EXPECT_NEAR(norm(point(m, static_cast<std::int64_t>(i))), 3.0, 1e-13);
    }
    // A pole mid-node stays on its meridian: node 5 of a north-cap triangle6
    // (pole, ring i, ring i + 1) is the midpoint of ring i + 1 -> pole.
    const Mesh s = uv_sphere(1.0, 6, 4, {{0, 0, 0}}, "triangle", 2);
    const std::int64_t* conn = s.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t c = 0; c < 6; ++c) {
        const Vec mid = point(s, conn[c * 6 + 5]);
        const Vec ring = point(s, conn[c * 6 + 2]);
        EXPECT_NEAR(std::atan2(mid[1], mid[0]), std::atan2(ring[1], ring[0]), 1e-13);
    }
    // The curved rim of a quadratic disk is on the circle.
    const Mesh d = disk(2.0, 3, {{0, 0, 0}}, "quad", 2);
    const auto& rim = region(d, "rim");
    const std::int64_t* dc = d.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t i = 0; i < rim.NumEntries(); ++i) {
        const std::int64_t cell = rim.mEntries.As<std::int64_t>()[2 * i];
        const std::int64_t f = rim.mEntries.As<std::int64_t>()[2 * i + 1];
        const Vec mid = point(d, dc[cell * 8 + 4 + f]);
        EXPECT_NEAR(std::hypot(mid[0], mid[1]), 2.0, 1e-13);
    }
}

TEST(Primitives, QuadraticMidNodesAreAppendedInFirstAppearanceOrder) {
    const Mesh m = box({{2, 1, 1}}, {{2, 1, 1}}, {{0, 0, 0}}, "hexahedron", 2);
    const std::int64_t* conn = m.Cells(0).Conn().As<std::int64_t>();
    // 12 corners, then the first cell numbers its 12 edges 12..23, and the
    // second cell reuses the shared face's four.
    for (std::size_t k = 0; k < 12; ++k)
        EXPECT_EQ(conn[8 + k], static_cast<std::int64_t>(12 + k));
    std::int64_t top = 23;
    for (std::size_t k = 0; k < 12; ++k)
        top = std::max(top, conn[20 + 8 + k]);
    EXPECT_EQ(top, 31);
    EXPECT_EQ(m.NumPoints(), 12u + 20u);
}

// --- reference cells -------------------------------------------------------

TEST(Primitives, ReferenceCellsDefaultCoversEveryType) {
    const Mesh m = reference_cells();
    std::vector<std::string> types;
    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
        types.push_back(m.Cells(b).Type());
        EXPECT_EQ(m.Cells(b).NumCells(), 1u);
    }
    const std::vector<std::string> want = {
        "vertex",   "line",         "triangle", "quad",      "tetra",     "hexahedron",
        "wedge",    "pyramid",      "line3",    "triangle6", "quad8",     "tetra10",
        "hexahedron20", "wedge15",  "pyramid13", "polyhedron8"};
    EXPECT_EQ(types, want);
}

TEST(Primitives, ReferenceCellsArePositiveAndApart) {
    const Mesh m = reference_cells({"tetra", "hexahedron", "wedge"});
    EXPECT_NEAR(volumes(m, 0)[0], 1.0 / 6.0, 1e-15);
    EXPECT_NEAR(volumes(m, 1)[0], 1.0, 1e-15);
    EXPECT_NEAR(volumes(m, 2)[0], 0.5, 1e-15);
    double prev_hi = -1e300;
    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
        const std::int64_t* c = m.Cells(b).Conn().As<std::int64_t>();
        double lo = 1e300, hi = -1e300;
        for (std::size_t k = 0; k < m.Cells(b).NodesPerCell(); ++k) {
            lo = std::min(lo, point(m, c[k])[0]);
            hi = std::max(hi, point(m, c[k])[0]);
        }
        EXPECT_GT(lo, prev_hi);
        prev_hi = hi;
    }
}

TEST(Primitives, ReferencePairsShareAFacetAndStayPositive) {
    for (const std::string type : {"line", "triangle", "quad", "tetra", "hexahedron", "wedge", "pyramid"}) {
        const Mesh m = reference_cells({type}, true);
        ASSERT_EQ(m.Cells(0).NumCells(), 2u) << type;
        const std::size_t w = m.Cells(0).NodesPerCell();
        const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
        std::set<std::int64_t> a(c, c + w), b(c + w, c + 2 * w);
        std::size_t shared = 0;
        for (auto v : a)
            shared += b.count(v);
        const std::map<std::string, std::size_t> want = {{"line", 1},  {"triangle", 2}, {"quad", 2},
                                                         {"tetra", 3}, {"hexahedron", 4}, {"wedge", 3},
                                                         {"pyramid", 4}};
        EXPECT_EQ(shared, want.at(type)) << type;
        if (type == "tetra" || type == "hexahedron" || type == "wedge")
            for (double v : volumes(m, 0))
                EXPECT_GT(v, 0.0) << type;
    }
}

TEST(Primitives, ReferenceQuadraticMidpointsAreExact) {
    const Mesh m = reference_cells({"tetra10", "hexahedron20"});
    const std::int64_t* t = m.Cells(0).Conn().As<std::int64_t>();
    // tetra10 node 4 is the midpoint of edge (0, 1), node 6 of edge (0, 2).
    for (const auto& [node, a, b] : {std::array<int, 3>{4, 0, 1}, std::array<int, 3>{6, 0, 2},
                                     std::array<int, 3>{9, 2, 3}}) {
        const Vec mid = point(m, t[node]), pa = point(m, t[a]), pb = point(m, t[b]);
        for (int k = 0; k < 3; ++k)
            EXPECT_EQ(mid[k], (pa[k] + pb[k]) * 0.5);
    }
}

// --- planar domains --------------------------------------------------------

namespace {

// Distinct undirected edges of every cell of the mesh (corner nodes only).
std::size_t distinct_edges(const Mesh& rMesh) {
    std::set<std::pair<std::int64_t, std::int64_t>> edges;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto view = rMesh.Cells(b);
        for (std::size_t c = 0; c < view.NumCells(); ++c) {
            std::vector<std::int64_t> nodes;
            if (view.IsRagged()) {
                const std::int64_t* row = view.Row(c);
                nodes.assign(row, row + view.RowSize(c));
            } else {
                const std::int64_t* conn = view.Conn().As<std::int64_t>();
                const std::size_t w = view.NodesPerCell();
                const std::size_t corners = view.Type().rfind("triangle", 0) == 0 ? 3 : 4;
                nodes.assign(conn + c * w, conn + c * w + corners);
            }
            for (std::size_t k = 0; k < nodes.size(); ++k) {
                std::int64_t a = nodes[k], e = nodes[(k + 1) % nodes.size()];
                if (a > e)
                    std::swap(a, e);
                edges.insert({a, e});
            }
        }
    }
    return edges.size();
}

std::size_t total_cells(const Mesh& rMesh) {
    std::size_t n = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        n += rMesh.Cells(b).NumCells();
    return n;
}

// Euler characteristic V - E + F of the corner graph.
long long euler(const Mesh& rMesh) {
    return static_cast<long long>(rMesh.NumPoints()) -
           static_cast<long long>(distinct_edges(rMesh)) +
           static_cast<long long>(total_cells(rMesh));
}

double total_area(const Mesh& rMesh) {
    double total = 0.0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        for (double a : areas(rMesh, b))
            total += a / 2;  // areas() returns twice the signed area
    return total;
}

}  // namespace

TEST(PrimitivesPlanar, TriangleCountsAreaAndRegions) {
    for (std::int64_t n = 1; n <= 5; ++n) {
        const Mesh m = reference_triangle(n, {{2.0, 3.0}}, {{1.0, -1.0, 4.0}});
        EXPECT_EQ(m.Cells(0).NumCells(), static_cast<std::size_t>(n * n));
        EXPECT_EQ(m.NumPoints(), static_cast<std::size_t>((n + 1) * (n + 2) / 2));
        for (double a : areas(m, 0))
            EXPECT_GT(a, 0.0);
        EXPECT_NEAR(total_area(m), 3.0, 1e-12);
        for (const char* name : {"ymin", "xmin", "hypotenuse"})
            EXPECT_EQ(region(m, name).NumEntries(), static_cast<std::size_t>(n)) << name;
        EXPECT_EQ(euler(m), 1);
        EXPECT_EQ(point(m, 0), (Vec{1.0, -1.0, 4.0}));
    }
    // Barycentric order: row j holds the nodes with i = 0 .. n - j, x fastest.
    const Mesh m = reference_triangle(2);
    EXPECT_EQ(point(m, 1), (Vec{0.5, 0.0, 0.0}));
    EXPECT_EQ(point(m, 3), (Vec{0.0, 0.5, 0.0}));
    EXPECT_EQ(point(m, 5), (Vec{0.0, 1.0, 0.0}));
}

TEST(PrimitivesPlanar, TriangleQuadraticMidpointsAreChordMidpoints) {
    const Mesh m = reference_triangle(3, {{1.0, 1.0}}, {{0, 0, 0}}, 2);
    ASSERT_EQ(m.Cells(0).Type(), "triangle6");
    const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t i = 0; i < m.Cells(0).NumCells(); ++i)
        for (int e = 0; e < 3; ++e) {
            const Vec a = point(m, c[6 * i + e]), b = point(m, c[6 * i + (e + 1) % 3]);
            const Vec mid = point(m, c[6 * i + 3 + e]);
            for (int k = 0; k < 3; ++k)
                EXPECT_EQ(mid[k], (a[k] + b[k]) * 0.5);
        }
}

TEST(PrimitivesPlanar, NgonFansFromTheCentre) {
    for (const std::int64_t sides : {3, 5, 8}) {
        for (const std::int64_t rings : {1, 3}) {
            const Mesh m = ngon(sides, 2.0, rings, {{1.0, 1.0, 0.0}});
            EXPECT_EQ(m.Cells(0).NumCells(), static_cast<std::size_t>(sides * rings * rings));
            EXPECT_EQ(m.NumPoints(), static_cast<std::size_t>(1 + sides * rings * (rings + 1) / 2));
            for (double a : areas(m, 0))
                EXPECT_GT(a, 0.0);
            const double want = 0.5 * static_cast<double>(sides) * 4.0 *
                                std::sin(2.0 * M_PI / static_cast<double>(sides));
            EXPECT_NEAR(total_area(m), want, 1e-12);
            EXPECT_EQ(region(m, "rim").NumEntries(), static_cast<std::size_t>(sides * rings));
            EXPECT_EQ(distinct_points(m), m.NumPoints());
            EXPECT_EQ(euler(m), 1);
            EXPECT_EQ(point(m, 0), (Vec{1.0, 1.0, 0.0}));
        }
    }
    // The plain fan: the centre, then the vertices.
    const Mesh m = ngon(4, 1.0, 1);
    EXPECT_EQ(m.NumPoints(), 5u);
    EXPECT_NEAR(point(m, 1)[0], 1.0, 1e-15);
    EXPECT_NEAR(point(m, 2)[1], 1.0, 1e-15);
    // order 2: a side's midpoint is on the straight side.
    const Mesh q = ngon(6, 1.0, 2, {{0, 0, 0}}, 2);
    ASSERT_EQ(q.Cells(0).Type(), "triangle6");
}

TEST(PrimitivesPlanar, AnnulusFullRingIsWeldedAndOnItsCircles) {
    for (const std::string type : {"quad", "triangle"}) {
        const Mesh m = annulus(0.5, 1.5, 12, 3, 6.283185307179586, {{1.0, 2.0, 3.0}}, type);
        const std::size_t cells = type == "quad" ? 36u : 72u;
        EXPECT_EQ(m.Cells(0).NumCells(), cells);
        EXPECT_EQ(m.NumPoints(), 12u * 4u);
        EXPECT_EQ(distinct_points(m), m.NumPoints());
        for (double a : areas(m, 0))
            EXPECT_GT(a, 0.0);
        EXPECT_EQ(euler(m), 0);
        EXPECT_EQ(region(m, "inner").NumEntries(), 12u);
        EXPECT_EQ(region(m, "outer").NumEntries(), 12u);
        EXPECT_EQ(m.FindRegion("start", meshioplusplus::RegionKind::Side), Mesh::npos);
        EXPECT_EQ(m.FindRegion("end", meshioplusplus::RegionKind::Side), Mesh::npos);
        // 12-gon areas: 12 * 0.5 * r0 * r1 * sin(2 pi / 12) summed over the rings.
        double want = 0.0;
        for (int k = 0; k < 3; ++k) {
            const double r0 = 0.5 + k / 3.0, r1 = 0.5 + (k + 1) / 3.0;
            want += 6.0 * (r1 * r1 - r0 * r0) * std::sin(M_PI / 6.0) * 2.0 / 2.0;
        }
        EXPECT_NEAR(total_area(m), want, 1e-12) << type;
        for (std::size_t i = 0; i < m.NumPoints(); ++i) {
            const Vec p = point(m, static_cast<std::int64_t>(i));
            const double r = std::hypot(p[0] - 1.0, p[1] - 2.0);
            EXPECT_GE(r, 0.5 - 1e-12);
            EXPECT_LE(r, 1.5 + 1e-12);
            EXPECT_EQ(p[2], 3.0);
        }
    }
}

TEST(PrimitivesPlanar, AnnulusSectorKeepsItsSeamApart) {
    const Mesh m = annulus(1.0, 2.0, 6, 2, M_PI / 2.0);
    EXPECT_EQ(m.NumPoints(), 7u * 3u);
    EXPECT_EQ(m.Cells(0).NumCells(), 12u);
    EXPECT_EQ(region(m, "start").NumEntries(), 2u);
    EXPECT_EQ(region(m, "end").NumEntries(), 2u);
    EXPECT_EQ(region(m, "inner").NumEntries(), 6u);
    EXPECT_EQ(region(m, "outer").NumEntries(), 6u);
    EXPECT_EQ(euler(m), 1);
    for (double a : areas(m, 0))
        EXPECT_GT(a, 0.0);
    // The start ray is on +x, the end ray on +y.
    EXPECT_NEAR(point(m, 0)[0], 1.0, 1e-15);
    EXPECT_NEAR(point(m, 0)[1], 0.0, 1e-15);
    EXPECT_NEAR(point(m, 6)[0], 0.0, 1e-15);
    EXPECT_NEAR(point(m, 6)[1], 1.0, 1e-15);
}

TEST(PrimitivesPlanar, AnnulusQuadraticRimsAreOnTheCircles) {
    for (const double angle : {6.283185307179586, 4.0}) {
        const Mesh m = annulus(1.0, 3.0, 7, 2, angle, {{0.5, -0.5, 0}}, "quad", 2);
        ASSERT_EQ(m.Cells(0).Type(), "quad8");
        const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
        for (const char* name : {"inner", "outer"}) {
            const auto& r = region(m, name);
            const double want = std::string(name) == "inner" ? 1.0 : 3.0;
            for (std::size_t i = 0; i < r.NumEntries(); ++i) {
                const std::int64_t cell = r.mEntries.As<std::int64_t>()[2 * i];
                const std::int64_t f = r.mEntries.As<std::int64_t>()[2 * i + 1];
                for (std::int64_t node : {c[8 * cell + f], c[8 * cell + 4 + f]}) {
                    const Vec p = point(m, node);
                    EXPECT_NEAR(std::hypot(p[0] - 0.5, p[1] + 0.5), want, 1e-12) << name;
                }
            }
        }
    }
}

TEST(PrimitivesPlanar, LshapeLosesAQuadrantAndKeepsTheNotch) {
    for (const std::string type : {"quad", "triangle"}) {
        const Mesh m = lshape(3, 2.0, {{0, 0, 1.0}}, type);
        EXPECT_EQ(m.Cells(0).NumCells(), type == "quad" ? 27u : 54u);
        EXPECT_EQ(m.NumPoints(), 7u * 7u - 9u);
        for (double a : areas(m, 0))
            EXPECT_GT(a, 0.0);
        EXPECT_NEAR(total_area(m), 3.0, 1e-12);
        for (const char* name : {"notch_x", "notch_y"})
            EXPECT_EQ(region(m, name).NumEntries(), 3u) << name;
        EXPECT_EQ(region(m, "xmin").NumEntries(), 6u);
        EXPECT_EQ(region(m, "ymax").NumEntries(), 6u);
        EXPECT_EQ(region(m, "xmax").NumEntries(), 3u);
        EXPECT_EQ(region(m, "ymin").NumEntries(), 3u);
        EXPECT_EQ(euler(m), 1);
        for (std::size_t i = 0; i < m.NumPoints(); ++i) {
            const Vec p = point(m, static_cast<std::int64_t>(i));
            EXPECT_FALSE(p[0] > 1e-12 && p[1] < -1e-12);
            EXPECT_EQ(p[2], 1.0);
        }
    }
    // The re-entrant corner is a node at the centre.
    const Mesh m = lshape(2);
    bool found = false;
    for (std::size_t i = 0; i < m.NumPoints(); ++i) {
        const Vec p = point(m, static_cast<std::int64_t>(i));
        found = found || (p[0] == 0.0 && p[1] == 0.0);
    }
    EXPECT_TRUE(found);
}

TEST(PrimitivesPlanar, PlateWithHoleSurroundsTheHole) {
    for (const std::string type : {"quad", "triangle"}) {
        const Mesh m = plate_with_hole(0.25, 2.0, 3, 2, {{1.0, 1.0, 0.0}}, type);
        EXPECT_EQ(m.Cells(0).NumCells(), type == "quad" ? 24u : 48u);
        EXPECT_EQ(m.NumPoints(), 3u * 12u);
        EXPECT_EQ(distinct_points(m), m.NumPoints());
        for (double a : areas(m, 0))
            EXPECT_GT(a, 0.0);
        EXPECT_EQ(euler(m), 0);
        EXPECT_EQ(region(m, "hole").NumEntries(), 12u);
        for (const char* name : {"xmin", "xmax", "ymin", "ymax"})
            EXPECT_EQ(region(m, name).NumEntries(), 3u) << name;
        // The area is the square minus a 12-gon hole of circumradius 0.25 (the
        // hole boundary is a polygon in the linear mesh).
        const double hole = 0.5 * 12.0 * 0.0625 * std::sin(2.0 * M_PI / 12.0);
        EXPECT_NEAR(total_area(m), 4.0 - hole, 1e-12) << type;
        for (std::size_t i = 0; i < m.NumPoints(); ++i) {
            const Vec p = point(m, static_cast<std::int64_t>(i));
            EXPECT_LE(std::fabs(p[0] - 1.0), 1.0 + 1e-12);
            EXPECT_LE(std::fabs(p[1] - 1.0), 1.0 + 1e-12);
        }
    }
}

TEST(PrimitivesPlanar, PlateWithHoleQuadraticHoleIsCircular) {
    const Mesh m = plate_with_hole(0.5, 2.0, 2, 3, {{0, 0, 0}}, "quad", 2);
    ASSERT_EQ(m.Cells(0).Type(), "quad8");
    const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
    const auto& hole = region(m, "hole");
    ASSERT_EQ(hole.NumEntries(), 8u);
    for (std::size_t i = 0; i < hole.NumEntries(); ++i) {
        const std::int64_t cell = hole.mEntries.As<std::int64_t>()[2 * i];
        const std::int64_t f = hole.mEntries.As<std::int64_t>()[2 * i + 1];
        for (std::int64_t node : {c[8 * cell + f], c[8 * cell + 4 + f]}) {
            const Vec p = point(m, node);
            EXPECT_NEAR(std::hypot(p[0], p[1]), 0.5, 1e-12);
        }
    }
    // The outer side's midpoints stay on the straight side.
    const auto& east = region(m, "xmax");
    for (std::size_t i = 0; i < east.NumEntries(); ++i) {
        const std::int64_t cell = east.mEntries.As<std::int64_t>()[2 * i];
        const std::int64_t f = east.mEntries.As<std::int64_t>()[2 * i + 1];
        EXPECT_EQ(point(m, c[8 * cell + 4 + f])[0], 1.0);
    }
}

TEST(PrimitivesPlanar, HoneycombIsARaggedPolygonBlock) {
    for (const auto& dims :
         std::vector<std::array<std::int64_t, 2>>{{{1, 1}}, {{3, 2}}, {{4, 5}}}) {
        const Mesh m = honeycomb(dims, 2.0, {{1.0, 2.0, 0.0}});
        ASSERT_EQ(m.NumCellBlocks(), 1u);
        const auto view = m.Cells(0);
        EXPECT_TRUE(view.IsRagged());
        EXPECT_EQ(view.Type().rfind("polygon", 0), 0u);
        const std::size_t cells = static_cast<std::size_t>(dims[0] * dims[1]);
        EXPECT_EQ(view.NumCells(), cells);
        double area = 0.0;
        for (std::size_t c = 0; c < cells; ++c) {
            const std::int64_t* row = view.Row(c);
            ASSERT_EQ(view.RowSize(c), 6u);
            double a2 = 0.0;
            for (std::size_t k = 0; k < 6; ++k) {
                const Vec p = point(m, row[k]), q = point(m, row[(k + 1) % 6]);
                a2 += p[0] * q[1] - q[0] * p[1];
                EXPECT_NEAR(std::hypot(q[0] - p[0], q[1] - p[1]), 2.0, 1e-12);
            }
            EXPECT_GT(a2, 0.0);
            area += a2 / 2;
        }
        EXPECT_NEAR(area, static_cast<double>(cells) * 1.5 * std::sqrt(3.0) * 4.0, 1e-9);
        EXPECT_EQ(euler(m), 1);
        EXPECT_EQ(distinct_points(m), m.NumPoints());
        // The boundary is the edges counted once: 6 F - 2 (E - boundary) edge uses.
        const auto& b = region(m, "boundary");
        EXPECT_EQ(b.mDim, 1);
        EXPECT_EQ(b.NumEntries(), 2 * distinct_edges(m) - 6 * cells);
        double xmin = 1e300;
        for (std::size_t i = 0; i < m.NumPoints(); ++i)
            xmin = std::min(xmin, point(m, static_cast<std::int64_t>(i))[0]);
        EXPECT_EQ(xmin, 1.0);
    }
}

TEST(PrimitivesPlanar, HoneycombNodesAreNumberedByFirstAppearance) {
    const Mesh m = honeycomb({{2, 2}});
    const std::int64_t* row = m.Cells(0).Row(0);
    for (std::size_t k = 0; k < 6; ++k)
        EXPECT_EQ(row[k], static_cast<std::int64_t>(k));
    // 2 x 2 hexagons share 7 vertices per row of pairs: 16 nodes in all.
    EXPECT_EQ(m.NumPoints(), 16u);
}
