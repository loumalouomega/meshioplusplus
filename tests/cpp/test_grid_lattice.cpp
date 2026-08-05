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
// The regular lattice: numbering, geometry and the public grid() constructor.
//
// The point of testing the lattice on its own -- rather than only through
// voxelize -- is that it has exactly one job and a closed-form answer, so a
// failure here localises immediately instead of looking like a distance bug.

// System includes
#include <array>
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
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

using meshioplusplus::attach_quality;
using meshioplusplus::compute_stats;
using meshioplusplus::grid;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
namespace d = meshioplusplus::detail;

namespace {

// The canonical hex lattice, written out here exactly as tests/cpp/test_refine.cpp
// writes it. It is duplicated rather than shared on purpose: it is the
// independent statement of the numbering that grid() must reproduce, and sharing
// one implementation between the claim and the check would make the comparison
// vacuous.
Mesh reference_hex_grid(std::size_t n) {
    const std::size_t s = n + 1;
    std::vector<std::vector<double>> pts;
    for (std::size_t k = 0; k < s; ++k)
        for (std::size_t j = 0; j < s; ++j)
            for (std::size_t i = 0; i < s; ++i)
                pts.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    const auto id = [s](std::size_t i, std::size_t j, std::size_t k) {
        return static_cast<std::int64_t>((k * s + j) * s + i);
    };
    std::vector<std::vector<std::int64_t>> cells;
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i)
                cells.push_back({id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k),
                                 id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1),
                                 id(i, j + 1, k + 1)});
    return mt::make_mesh(std::move(pts), "hexahedron", std::move(cells));
}

}  // namespace

TEST(GridLattice, CountsFollowTheCellCounts) {
    const Mesh m = grid({{3, 4, 5}});
    EXPECT_EQ(m.NumPoints(), 4u * 5u * 6u);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "hexahedron");
    EXPECT_EQ(m.Cells(0).NumCells(), 3u * 4u * 5u);
}

TEST(GridLattice, ReproducesTheRepositorysExistingNumbering) {
    // grid() must agree cell-for-cell and point-for-point with the fixture the
    // refine/partition/gradient/interpolate suites already pin their expected
    // values against. A different linearization would be a silent fork.
    for (std::size_t n : {1u, 2u, 3u}) {
        const Mesh want = reference_hex_grid(n);
        const Mesh got = grid({{static_cast<std::int64_t>(n), static_cast<std::int64_t>(n),
                                static_cast<std::int64_t>(n)}});
        ASSERT_EQ(got.NumPoints(), want.NumPoints()) << "n=" << n;
        const double* pg = got.Points().As<double>();
        const double* pw = want.Points().As<double>();
        for (std::size_t i = 0; i < want.NumPoints() * 3; ++i)
            ASSERT_DOUBLE_EQ(pg[i], pw[i]) << "n=" << n << " point component " << i;

        ASSERT_EQ(got.Cells(0).NumCells(), want.Cells(0).NumCells()) << "n=" << n;
        const std::int64_t* cg = got.Cells(0).Conn().As<std::int64_t>();
        const std::int64_t* cw = want.Cells(0).Conn().As<std::int64_t>();
        for (std::size_t i = 0; i < want.Cells(0).NumCells() * 8; ++i)
            ASSERT_EQ(cg[i], cw[i]) << "n=" << n << " conn entry " << i;
    }
}

TEST(GridLattice, EveryCellIsAPerfectHexahedron) {
    // The free structural oracle: a right parallelepiped has scaled Jacobian
    // exactly 1, so any transposed axis or mis-wound cell shows up as an
    // inverted cell rather than as a picture that looks slightly wrong.
    const Mesh m = grid({{3, 3, 3}}, {{0.0, 0.0, 0.0}}, {{2.0, 0.5, 1.5}});
    const Mesh q = attach_quality(m);
    const NDArray& sj = q.CellData("quality:scaled_jacobian", 0);
    ASSERT_EQ(sj.Size(), 27u);
    const double* v = sj.As<double>();
    for (std::size_t c = 0; c < sj.Size(); ++c)
        EXPECT_DOUBLE_EQ(v[c], 1.0) << "cell " << c;

    const meshioplusplus::StatsReport st = compute_stats(m);
    EXPECT_EQ(st.mNumInverted, 0);
    // 3x3x3 cells of 2 x 0.5 x 1.5 each.
    EXPECT_DOUBLE_EQ(st.mSignedVolume, 27.0 * 2.0 * 0.5 * 1.5);
}

TEST(GridLattice, TheHexahedronOracleActuallyFires) {
    // Exchanging the bottom and top rings mirrors every cell through its
    // mid-plane, which is what a swapped axis in the connectivity formula would
    // do. Both halves of the oracle above must notice: the signed volume flips
    // and the scaled Jacobian stops being 1.
    //
    // Note a *weaker* sabotage is not enough -- transposing two corners of one
    // ring twists the cell without necessarily negating its signed volume, and
    // an earlier version of this test used exactly that and passed. The oracle
    // detects mis-orientation, not every possible wrong ordering, and the test
    // has to say which.
    const Mesh good = reference_hex_grid(2);
    std::vector<std::vector<double>> pts;
    const double* p = good.Points().As<double>();
    for (std::size_t i = 0; i < good.NumPoints(); ++i)
        pts.push_back({p[i * 3], p[i * 3 + 1], p[i * 3 + 2]});
    std::vector<std::vector<std::int64_t>> cells;
    const std::int64_t* src = good.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t c = 0; c < good.Cells(0).NumCells(); ++c) {
        std::vector<std::int64_t> row(8);
        for (std::size_t i = 0; i < 8; ++i)
            row[i] = src[c * 8 + ((i + 4) % 8)];
        cells.push_back(std::move(row));
    }
    const Mesh bad = mt::make_mesh(std::move(pts), "hexahedron", std::move(cells));
    EXPECT_GT(compute_stats(bad).mNumInverted, 0) << "the inversion oracle does not fire";

    const Mesh q = attach_quality(bad);
    const NDArray& sj = q.CellData("quality:scaled_jacobian", 0);
    const double* v = sj.As<double>();
    bool all_one = true;
    for (std::size_t c = 0; c < sj.Size(); ++c)
        all_one = all_one && v[c] == 1.0;
    EXPECT_FALSE(all_one) << "the scaled-Jacobian oracle does not fire";
}

TEST(GridLattice, CoordinatesAreOriginPlusIndexTimesSpacing) {
    const Mesh m = grid({{4, 1, 1}}, {{-1.0, 2.0, 0.5}}, {{0.25, 3.0, 7.0}});
    const double* p = m.Points().As<double>();
    // Point (i, 0, 0) is at origin.x + i * 0.25 exactly -- an accumulating loop
    // would drift here and make the numpy twin's implementation observable.
    for (std::int64_t i = 0; i <= 4; ++i)
        EXPECT_DOUBLE_EQ(p[static_cast<std::size_t>(i) * 3 + 0],
                         -1.0 + static_cast<double>(i) * 0.25);
    EXPECT_DOUBLE_EQ(p[1], 2.0);
    EXPECT_DOUBLE_EQ(p[2], 0.5);
}

TEST(GridLattice, AnEmptyLatticeIsAnEmptyMeshNotAThrow) {
    const Mesh m = grid({{0, 0, 0}});
    EXPECT_EQ(m.NumPoints(), 0u);
    EXPECT_EQ(m.NumCellBlocks(), 0u);
}

TEST(GridLattice, FromBoundsDividesTheBoxExactly) {
    const d::LatticeSpec spec =
        d::lattice_from_bounds({{0.0, 0.0, 0.0}}, {{1.0, 2.0, 4.0}}, {{2, 4, 8}});
    EXPECT_DOUBLE_EQ(spec.mSpacing[0], 0.5);
    EXPECT_DOUBLE_EQ(spec.mSpacing[1], 0.5);
    EXPECT_DOUBLE_EQ(spec.mSpacing[2], 0.5);
    EXPECT_EQ(d::lattice_num_cells(spec), 64);
    EXPECT_EQ(d::lattice_num_points(spec), 3 * 5 * 9);
}

TEST(GridLattice, FromCellSizeCoversTheBoxRatherThanClippingIt) {
    // 1.0 does not divide 2.5, so the lattice must grow to 3 cells rather than
    // stop at 2 and leave a quarter of the box uncovered.
    const d::LatticeSpec spec =
        d::lattice_from_cell_size({{0.0, 0.0, 0.0}}, {{2.5, 1.0, 1.0}}, {{1.0, 1.0, 1.0}});
    EXPECT_EQ(spec.mDims[0], 3);
    EXPECT_DOUBLE_EQ(spec.mSpacing[0], 1.0);
    EXPECT_GE(spec.mOrigin[0] + static_cast<double>(spec.mDims[0]) * spec.mSpacing[0], 2.5);
}

TEST(GridLattice, PointBboxMatchesComputeStats) {
    const Mesh m = mt::tet_mesh();
    std::array<double, 3> lo{}, hi{};
    ASSERT_TRUE(d::point_bbox(m, lo, hi));
    const meshioplusplus::StatsReport st = compute_stats(m);
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_DOUBLE_EQ(lo[k], st.mBBoxMin[k]) << "axis " << k;
        EXPECT_DOUBLE_EQ(hi[k], st.mBBoxMax[k]) << "axis " << k;
    }
}

TEST(GridLattice, PointBboxReportsAnEmptyMeshRatherThanGuessing) {
    Mesh m;
    m.AssignPoints(NDArray(meshioplusplus::DType::Float64, {0u, 3u}));
    std::array<double, 3> lo{{7.0, 7.0, 7.0}}, hi{{9.0, 9.0, 9.0}};
    EXPECT_FALSE(d::point_bbox(m, lo, hi));
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_DOUBLE_EQ(lo[k], 0.0);
        EXPECT_DOUBLE_EQ(hi[k], 0.0);
    }
}

TEST(GridLattice, BadArgumentsRaiseByName) {
    EXPECT_THROW(grid({{-1, 1, 1}}), std::invalid_argument);
    EXPECT_THROW(grid({{1, 1, 1}}, {{0, 0, 0}}, {{0.0, 1.0, 1.0}}), std::invalid_argument);
    // The budget refusal, which exists because 512^3 is ~11.8 GB and is asked
    // for by typo more often than on purpose.
    try {
        grid({{100, 100, 100}}, {{0, 0, 0}}, {{1, 1, 1}}, 1000);
        FAIL() << "the cell budget did not refuse";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("max_cells"), std::string::npos) << e.what();
    }
}
