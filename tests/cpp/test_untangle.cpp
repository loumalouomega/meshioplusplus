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
// Untangling inverted tetrahedra and hexahedra (roadmap 7.4.2).

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/untangle.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

Mesh ut_hex_box(int N) {
    std::vector<std::vector<double>> pts;
    for (int k = 0; k <= N; ++k)
        for (int j = 0; j <= N; ++j)
            for (int i = 0; i <= N; ++i)
                pts.push_back({double(i) / N, double(j) / N, double(k) / N});
    const auto id = [&](int i, int j, int k) {
        return std::int64_t((k * (N + 1) + j) * (N + 1) + i);
    };
    std::vector<std::vector<std::int64_t>> hexes;
    for (int k = 0; k < N; ++k)
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i)
                hexes.push_back({id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k),
                                 id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1),
                                 id(i, j + 1, k + 1)});
    return mt::make_mesh(pts, "hexahedron", hexes);
}

Mesh ut_move(const Mesh& rMesh, std::size_t Node, double X, double Y, double Z) {
    Mesh m = rMesh;
    NDArray p(DType::Float64, {m.NumPoints(), std::size_t{3}});
    std::copy(rMesh.Points().As<double>(), rMesh.Points().As<double>() + 3 * m.NumPoints(),
              p.As<double>());
    p.As<double>()[3 * Node] = X;
    p.As<double>()[3 * Node + 1] = Y;
    p.As<double>()[3 * Node + 2] = Z;
    m.AssignPoints(std::move(p));
    return m;
}

}  // namespace

TEST(Untangle, TetCubeWithANodePushedThroughAFace) {
    // The roadmap probe: a tetrahedral cube, an interior node pushed past its
    // neighbours, untangled to every signed volume positive.
    ConvertCellsOptions co;
    co.mMode = ConvertCellsMode::Simplexify;
    const Mesh tets = convert_cells(ut_hex_box(4), co).mMesh;
    const std::size_t node = (1 * 5 + 1) * 5 + 1;  // (0.25, 0.25, 0.25)
    const Mesh tangled = ut_move(tets, node, 0.62, 0.58, 0.6);
    ASSERT_GT(compute_quality(tangled).mNumInverted, 0);
    const UntangleResult r = untangle(tangled);
    EXPECT_GT(r.mNumInvertedBefore, 0);
    EXPECT_EQ(r.mNumInvertedAfter, 0);
    EXPECT_GT(r.mMinQualityAfter, 0.0);
    EXPECT_EQ(compute_quality(r.mMesh).mNumInverted, 0);
    // The boundary did not move.
    const double* a = tangled.Points().As<double>();
    const double* b = r.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < tangled.NumPoints(); ++i) {
        bool boundary = false;
        for (int k = 0; k < 3; ++k)
            boundary = boundary || a[3 * i + k] == 0.0 || a[3 * i + k] == 1.0;
        if (boundary)
            for (int k = 0; k < 3; ++k)
                EXPECT_EQ(a[3 * i + k], b[3 * i + k]);
    }
}

TEST(Untangle, HexBoxAndValidMeshesAreLeftAlone) {
    const Mesh box = ut_hex_box(4);
    const UntangleResult valid = untangle(box);
    EXPECT_EQ(valid.mNumInvertedBefore, 0);
    EXPECT_EQ(valid.mSweeps, 0);
    EXPECT_EQ(valid.mMoves, 0);
    const std::size_t node = (2 * 5 + 2) * 5 + 2;  // the centre
    const Mesh tangled = ut_move(box, node, 0.5, 0.5, 0.8);
    const UntangleResult r = untangle(tangled);
    EXPECT_GT(r.mNumInvertedBefore, 0);
    EXPECT_EQ(r.mNumInvertedAfter, 0);
}

TEST(Untangle, Refusals) {
    EXPECT_THROW(untangle(mt::tri_mesh()), std::invalid_argument);
    UntangleOptions o;
    o.mMaxSweeps = 0;
    EXPECT_THROW(untangle(ut_hex_box(2), o), std::invalid_argument);
    o = UntangleOptions{};
    o.mFrozen = {1, 0};
    EXPECT_THROW(untangle(ut_hex_box(2), o), std::invalid_argument);
}
