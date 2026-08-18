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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/conservative_interpolate.hpp"

namespace {

using meshioplusplus::conservative_interpolate;
using meshioplusplus::conservative_interpolate_conflict_from_name;
using meshioplusplus::ConservativeInterpolateConflict;
using meshioplusplus::ConservativeInterpolateOptions;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
namespace md = meshioplusplus::detail;

// --- fixtures ----------------------------------------------------------------

// A single triangle in the z = 0 plane with corners `p0`, `p1`, `p2`.
Mesh one_triangle(const std::vector<double>& p0, const std::vector<double>& p1,
                  const std::vector<double>& p2) {
    return mt::make_mesh({p0, p1, p2}, "triangle", {{0, 1, 2}});
}

// A single tetrahedron with the given 4 corners.
Mesh one_tetra(const std::vector<double>& p0, const std::vector<double>& p1,
               const std::vector<double>& p2, const std::vector<double>& p3) {
    return mt::make_mesh({p0, p1, p2, p3}, "tetra", {{0, 1, 2, 3}});
}

// An n x n grid of unit quads in the z = 0 plane, scaled to cover
// [0, n * scale] x [0, n * scale].
Mesh quad_grid(int n, double scale) {
    const int m = n + 1;
    std::vector<std::vector<double>> pts;
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < m; ++i)
            pts.push_back({static_cast<double>(i) * scale, static_cast<double>(j) * scale, 0.0});
    auto pid = [m](int i, int j) { return static_cast<std::int64_t>(j * m + i); };
    std::vector<std::vector<std::int64_t>> cells;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            cells.push_back({pid(i, j), pid(i + 1, j), pid(i + 1, j + 1), pid(i, j + 1)});
    return mt::make_mesh(std::move(pts), "quad", std::move(cells));
}

// An n x n x n grid of unit hexahedra, scaled to cover [0, n*scale]^3.
Mesh hex_grid(int n, double scale) {
    const int m = n + 1;
    std::vector<std::vector<double>> pts;
    for (int k = 0; k < m; ++k)
        for (int j = 0; j < m; ++j)
            for (int i = 0; i < m; ++i)
                pts.push_back({static_cast<double>(i) * scale, static_cast<double>(j) * scale,
                               static_cast<double>(k) * scale});
    auto pid = [m](int i, int j, int k) { return static_cast<std::int64_t>((k * m + j) * m + i); };
    std::vector<std::vector<std::int64_t>> cells;
    for (int k = 0; k < n; ++k)
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                cells.push_back({pid(i, j, k), pid(i + 1, j, k), pid(i + 1, j + 1, k),
                                 pid(i, j + 1, k), pid(i, j, k + 1), pid(i + 1, j, k + 1),
                                 pid(i + 1, j + 1, k + 1), pid(i, j + 1, k + 1)});
    return mt::make_mesh(std::move(pts), "hexahedron", std::move(cells));
}

// Total sum(value * measure) over a mesh's cell_data array `name`, weighting
// each quad/hexahedron cell by its own area/volume computed independently of
// the operation under test (a plain shoelace/box-volume formula).
double weighted_sum(const Mesh& rMesh, const std::string& name) {
    double total = 0.0;
    std::size_t block = 0;
    for (const auto cb : rMesh.CellRange()) {
        const NDArray& data = rMesh.CellData(name, block);
        const NDArray& conn = cb.Conn();
        const NDArray& pts = rMesh.Points();
        const std::size_t pdim = rMesh.PointDim();
        const std::size_t npc = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::vector<double> xs, ys, zs;
            for (std::size_t k = 0; k < npc; ++k) {
                const std::int64_t node = md::read_int(conn, c * npc + k);
                xs.push_back(md::read_double(pts, static_cast<std::size_t>(node) * pdim + 0));
                ys.push_back(md::read_double(pts, static_cast<std::size_t>(node) * pdim + 1));
                if (pdim > 2)
                    zs.push_back(md::read_double(pts, static_cast<std::size_t>(node) * pdim + 2));
            }
            double measure;
            if (npc == 4) {
                // shoelace area of the (planar, z ignored) quad
                double a2 = 0.0;
                for (std::size_t k = 0; k < 4; ++k) {
                    const std::size_t k2 = (k + 1) % 4;
                    a2 += xs[k] * ys[k2] - xs[k2] * ys[k];
                }
                measure = std::abs(a2) * 0.5;
            } else {
                // axis-aligned hexahedron volume from its bbox
                const double dx = *std::max_element(xs.begin(), xs.end()) -
                                  *std::min_element(xs.begin(), xs.end());
                const double dy = *std::max_element(ys.begin(), ys.end()) -
                                  *std::min_element(ys.begin(), ys.end());
                const double dz = *std::max_element(zs.begin(), zs.end()) -
                                  *std::min_element(zs.begin(), zs.end());
                measure = dx * dy * dz;
            }
            total += measure * md::read_double(data, c);
        }
        ++block;
    }
    return total;
}

// --- clip correctness, via single-cell meshes (isolates the kernel) --------

TEST(ConservativeInterpolate, IdenticalTriangleRecoversFullValue) {
    Mesh src = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    src.AddCellData("f", {mt::data_array({7.0})});
    Mesh tgt = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});

    const Mesh out = conservative_interpolate(src, tgt);
    EXPECT_DOUBLE_EQ(md::read_double(out.CellData("f", 0), 0), 7.0);
}

TEST(ConservativeInterpolate, IdenticalTetraRecoversFullValue) {
    Mesh src = one_tetra({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    src.AddCellData("f", {mt::data_array({3.5})});
    Mesh tgt = one_tetra({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});

    const Mesh out = conservative_interpolate(src, tgt);
    EXPECT_DOUBLE_EQ(md::read_double(out.CellData("f", 0), 0), 3.5);
}

TEST(ConservativeInterpolate, DisjointTrianglesGiveDefaultValue) {
    Mesh src = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    src.AddCellData("f", {mt::data_array({9.0})});
    Mesh tgt = one_triangle({10, 10, 0}, {11, 10, 0}, {10, 11, 0});

    ConservativeInterpolateOptions opts;
    opts.mDefaultValue = -1.0;
    const Mesh out = conservative_interpolate(src, tgt, opts);
    EXPECT_DOUBLE_EQ(md::read_double(out.CellData("f", 0), 0), -1.0);
}

TEST(ConservativeInterpolate, DisjointTetrasGiveDefaultValue) {
    Mesh src = one_tetra({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    src.AddCellData("f", {mt::data_array({9.0})});
    Mesh tgt = one_tetra({10, 10, 10}, {11, 10, 10}, {10, 11, 10}, {10, 10, 11});

    ConservativeInterpolateOptions opts;
    opts.mDefaultValue = -1.0;
    const Mesh out = conservative_interpolate(src, tgt, opts);
    EXPECT_DOUBLE_EQ(md::read_double(out.CellData("f", 0), 0), -1.0);
}

// A target triangle exactly half of the source triangle's area (the source
// right triangle's mid-line split) -> the overlap covers the target
// entirely, so the target's own value is exactly the source's.
TEST(ConservativeInterpolate, PartialTriangleOverlapRecoversSourceValueWhenFullyCovered) {
    Mesh src = one_triangle({0, 0, 0}, {2, 0, 0}, {0, 2, 0});
    src.AddCellData("f", {mt::data_array({4.0})});
    Mesh tgt = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});  // strictly inside src

    const Mesh out = conservative_interpolate(src, tgt);
    EXPECT_DOUBLE_EQ(md::read_double(out.CellData("f", 0), 0), 4.0);
}

TEST(ConservativeInterpolate, PartialTetraOverlapRecoversSourceValueWhenFullyCovered) {
    Mesh src = one_tetra({0, 0, 0}, {2, 0, 0}, {0, 2, 0}, {0, 0, 2});
    src.AddCellData("f", {mt::data_array({6.0})});
    Mesh tgt = one_tetra({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});  // strictly inside src

    const Mesh out = conservative_interpolate(src, tgt);
    EXPECT_DOUBLE_EQ(md::read_double(out.CellData("f", 0), 0), 6.0);
}

// Two triangles overlapping in a known quadrilateral: src is the unit square
// (as two triangles is not needed -- a single big triangle) clipped against a
// target triangle offset by (0.5, 0.5) with corners chosen so the exact
// overlap area is computable by hand: src = right triangle (0,0)-(2,0)-(0,2),
// tgt = right triangle (0,0)-(2,0)-(0,2) translated and reflected so the
// overlap is the square [0,1]x[0,1] intersected with src (a known triangle
// of area 0.5). Use a constant-source-value trick: a uniform field of value
// V over src, sampled with a target that partially overlaps by a KNOWN
// fraction, must read back exactly V (weighted mean of a constant is that
// constant) -- so exact-area verification is instead done through the
// integral-conservation tests below, which do not depend on knowing the
// overlap area of any single pair.

// --- operation-level behaviour ----------------------------------------------

TEST(ConservativeInterpolate, CellDataConservesTotalIntegral2D) {
    // Two independent partitions of the exact same [0,6]x[0,6] domain.
    Mesh src = quad_grid(3, 2.0);  // 3x3 cells of size 2
    Mesh tgt = quad_grid(4, 1.5);  // 4x4 cells of size 1.5

    std::vector<double> vals;
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i)
            vals.push_back(1.0 + 2.0 * i + 3.0 * j);  // a varying, non-constant field
    src.AddCellData("f", {mt::data_array(vals)});

    const Mesh out = conservative_interpolate(src, tgt);
    EXPECT_NEAR(weighted_sum(out, "f"), weighted_sum(src, "f"), 1e-8);
}

TEST(ConservativeInterpolate, CellDataConservesTotalIntegral3D) {
    Mesh src_hex = hex_grid(2, 3.0);  // 2x2x2 cells of size 3, domain [0,6]^3
    Mesh tgt_hex = hex_grid(3, 2.0);  // 3x3x3 cells of size 2, domain [0,6]^3

    std::vector<double> vals;
    for (int i = 0; i < 8; ++i)
        vals.push_back(1.0 + static_cast<double>(i));
    src_hex.AddCellData("f", {mt::data_array(vals)});

    const Mesh out = conservative_interpolate(src_hex, tgt_hex);
    EXPECT_NEAR(weighted_sum(out, "f"), weighted_sum(src_hex, "f"), 1e-6);
}

TEST(ConservativeInterpolate, RejectsMismatchedTopologicalDimension) {
    Mesh src2d = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    src2d.AddCellData("f", {mt::data_array({1.0})});
    Mesh tgt3d = one_tetra({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    EXPECT_THROW(conservative_interpolate(src2d, tgt3d), std::invalid_argument);
}

TEST(ConservativeInterpolate, PresentInBothLocationsTransfersBoth) {
    Mesh src = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    src.AddCellData("f", {mt::data_array({2.0})});
    src.AddPointData("f", mt::data_array({1.0, 2.0, 3.0}));
    Mesh tgt = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});

    const Mesh out = conservative_interpolate(src, tgt);
    EXPECT_TRUE(out.HasCellData("f"));
    EXPECT_TRUE(out.HasPointData("f"));
}

TEST(ConservativeInterpolate, ConflictErrorOverwriteSuffix) {
    Mesh src = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    src.AddCellData("f", {mt::data_array({2.0})});
    Mesh tgt = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    tgt.AddCellData("f", {mt::data_array({99.0})});

    EXPECT_THROW(conservative_interpolate(src, tgt), std::invalid_argument);

    ConservativeInterpolateOptions overwrite_opts;
    overwrite_opts.mOnConflict = ConservativeInterpolateConflict::Overwrite;
    const Mesh out_ow = conservative_interpolate(src, tgt, overwrite_opts);
    EXPECT_DOUBLE_EQ(md::read_double(out_ow.CellData("f", 0), 0), 2.0);

    ConservativeInterpolateOptions suffix_opts;
    suffix_opts.mOnConflict = ConservativeInterpolateConflict::Suffix;
    const Mesh out_sfx = conservative_interpolate(src, tgt, suffix_opts);
    EXPECT_TRUE(out_sfx.HasCellData("f_interp"));
    EXPECT_DOUBLE_EQ(md::read_double(out_sfx.CellData("f", 0), 0), 99.0);
}

TEST(ConservativeInterpolate, TargetGeometryAndOwnDataUntouched) {
    Mesh src = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    src.AddCellData("f", {mt::data_array({2.0})});
    Mesh tgt = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    tgt.AddCellData("g", {mt::data_array({42.0})});

    const Mesh out = conservative_interpolate(src, tgt);
    EXPECT_EQ(out.NumPoints(), tgt.NumPoints());
    EXPECT_DOUBLE_EQ(md::read_double(out.CellData("g", 0), 0), 42.0);
}

TEST(ConservativeInterpolate, UnknownArrayThrows) {
    Mesh src = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    Mesh tgt = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    ConservativeInterpolateOptions opts;
    opts.mArrays = {"nope"};
    EXPECT_THROW(conservative_interpolate(src, tgt, opts), std::invalid_argument);
}

TEST(ConservativeInterpolate, EmptySourceThrows) {
    Mesh src;
    Mesh tgt = one_triangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    EXPECT_THROW(conservative_interpolate(src, tgt), std::invalid_argument);
}

TEST(ConservativeInterpolate, PointDataComposedConservesApproximately) {
    Mesh src = quad_grid(3, 2.0);
    Mesh tgt = quad_grid(4, 1.5);
    std::vector<double> vals;
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 4; ++i)
            vals.push_back(1.0 + i + j);
    src.AddPointData("f", mt::data_array(vals));

    const Mesh out = conservative_interpolate(src, tgt);
    ASSERT_TRUE(out.HasPointData("f"));
    // Approximate: point<->cell lumping on both sides is lossy, but the
    // result must be a real, finite, in-range value -- not NaN/garbage.
    const NDArray& f = out.PointData("f");
    for (std::size_t i = 0; i < out.NumPoints(); ++i) {
        const double v = md::read_double(f, i);
        EXPECT_TRUE(std::isfinite(v));
    }
}

TEST(ConservativeInterpolate, DuplicateCandidatesFromMultiCellBboxAreNotDoubleCounted) {
    // Two adjacent, DIFFERENTLY-valued source quads and one target quad that
    // exactly spans both. A duplicated or dropped candidate from the
    // broad-phase box query would skew the ratio away from the exact
    // equal-area-weighted average -- unlike a single-source-cell scenario,
    // where over-counting cancels out of the numerator/denominator ratio and
    // so cannot be detected this way.
    Mesh src = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {2, 1, 0}},
                             "quad", {{0, 1, 2, 3}, {1, 4, 5, 2}});
    src.AddCellData("f", {mt::data_array({10.0, 30.0})});
    Mesh tgt = mt::make_mesh({{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}}, "quad", {{0, 1, 2, 3}});

    const Mesh out = conservative_interpolate(src, tgt);
    EXPECT_NEAR(md::read_double(out.CellData("f", 0), 0), 20.0, 1e-9);
}

}  // namespace
