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
// check_conformity: hand-built meshes whose defects are known by construction.

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/check_conformity.hpp"

using meshioplusplus::check_conformity;
using meshioplusplus::ConformityOptions;
using meshioplusplus::ConformityResult;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

namespace {

// The flags as a vector, or empty when the array is missing.
std::vector<std::int64_t> cc_flags(const ConformityResult& rR) {
    std::vector<std::int64_t> out;
    if (!rR.mMesh.HasPointData(meshioplusplus::kConformityHangingName))
        return out;
    const auto& a = rR.mMesh.PointData(meshioplusplus::kConformityHangingName);
    const std::int64_t* p = a.As<std::int64_t>();
    out.assign(p, p + rR.mMesh.NumPoints());
    return out;
}

// A hexahedron's corners in VTK order, from its lo corner and edge length.
std::vector<std::int64_t> cc_hex(const std::vector<std::vector<double>>& rPts,
                                 const std::array<double, 3>& rLo, double H) {
    std::vector<std::int64_t> row;
    const std::array<std::array<double, 3>, 8> corner = {
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
    for (const auto& c : corner) {
        for (std::size_t p = 0; p < rPts.size(); ++p)
            if (rPts[p][0] == rLo[0] + c[0] * H && rPts[p][1] == rLo[1] + c[1] * H &&
                rPts[p][2] == rLo[2] + c[2] * H) {
                row.push_back(static_cast<std::int64_t>(p));
                break;
            }
    }
    return row;
}

// A unit hexahedron at x in [0, 1] beside eight half-size ones filling
// [1, 2] x [0, 1]^2: the coarse face x = 1 has five hanging nodes on it.
Mesh cc_coarse_beside_fine() {
    std::vector<std::vector<double>> pts;
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j)
            for (int k = 0; k <= 2; ++k)
                pts.push_back({1.0 + 0.5 * i, 0.5 * j, 0.5 * k});
    for (double y : {0.0, 1.0})
        for (double z : {0.0, 1.0})
            pts.push_back({0.0, y, z});
    std::vector<std::vector<std::int64_t>> hexes;
    hexes.push_back(cc_hex(pts, {0.0, 0.0, 0.0}, 1.0));
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
                hexes.push_back(cc_hex(pts, {1.0 + 0.5 * i, 0.5 * j, 0.5 * k}, 0.5));
    return mt::make_mesh(std::move(pts), "hexahedron", std::move(hexes));
}

}  // namespace

TEST(CheckConformity, ASingleHexahedronConforms) {
    const ConformityResult r = check_conformity(mt::hex_mesh());
    EXPECT_EQ(r.mDimension, 3);
    EXPECT_TRUE(r.mConforming);
    EXPECT_EQ(r.mInterfaceMode, "conforming");
    EXPECT_EQ(r.mNumHangingNodes, 0);
    EXPECT_EQ(r.mNumPartialOverlaps, 0);
    EXPECT_EQ(r.mNumUnweldedFacets, 0);
    EXPECT_EQ(r.mNumBoundaryFacets, 6);
    for (std::int64_t f : cc_flags(r))
        EXPECT_EQ(f, 0);
}

TEST(CheckConformity, ACoarseFaceBesideFineOnesHasFiveHangingNodes) {
    const Mesh m = cc_coarse_beside_fine();
    const ConformityResult r = check_conformity(m);
    EXPECT_FALSE(r.mConforming);
    EXPECT_EQ(r.mInterfaceMode, "proximity");
    EXPECT_EQ(r.mNumHangingNodes, 5);
    EXPECT_EQ(r.mNumTJunctions, 0);
    // The coarse face against each of the four fine faces it holds.
    EXPECT_EQ(r.mNumPartialOverlaps, 4);
    const std::vector<std::int64_t> flags = cc_flags(r);
    ASSERT_EQ(flags.size(), m.NumPoints());
    const double* p = m.Points().As<double>();
    for (std::size_t i = 0; i < m.NumPoints(); ++i) {
        const bool on_face = p[i * 3] == 1.0;
        const bool corner = (p[i * 3 + 1] == 0.0 || p[i * 3 + 1] == 1.0) &&
                            (p[i * 3 + 2] == 0.0 || p[i * 3 + 2] == 1.0);
        EXPECT_EQ(flags[i], on_face && !corner ? 1 : 0) << i;
    }
}

TEST(CheckConformity, GeometryAndDataPassThroughBitIdentical) {
    Mesh m = cc_coarse_beside_fine();
    NDArray marker = NDArray::Uninit(meshioplusplus::DType::Float64, {m.NumPoints()});
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        marker.As<double>()[i] = static_cast<double>(i) * 0.5;
    m.AddPointData("marker", std::move(marker));
    const ConformityResult r = check_conformity(m);
    ASSERT_EQ(r.mMesh.NumPoints(), m.NumPoints());
    for (std::size_t i = 0; i < m.NumPoints() * 3; ++i)
        EXPECT_EQ(r.mMesh.Points().As<double>()[i], m.Points().As<double>()[i]);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), m.NumCellBlocks());
    const auto& a = r.mMesh.Cells(0).Conn();
    const auto& b = m.Cells(0).Conn();
    ASSERT_EQ(a.Size(), b.Size());
    for (std::size_t i = 0; i < a.Size(); ++i)
        EXPECT_EQ(a.As<std::int64_t>()[i], b.As<std::int64_t>()[i]);
    ASSERT_TRUE(r.mMesh.HasPointData("marker"));
    EXPECT_EQ(r.mMesh.PointData("marker").As<double>()[3], 1.5);

    ConformityOptions o;
    o.mAttachFlags = false;
    EXPECT_FALSE(check_conformity(m, o).mMesh.HasPointData(meshioplusplus::kConformityHangingName));
}

TEST(CheckConformity, ATJunctionOnAQuadEdge) {
    const Mesh m = mt::make_mesh({{0, 0, 0},
                                  {1, 0, 0},
                                  {1, 1, 0},
                                  {0, 1, 0},
                                  {2, 0, 0},
                                  {2, 0.5, 0},
                                  {2, 1, 0},
                                  {1, 0.5, 0}},
                                 "quad", {{0, 1, 2, 3}, {1, 4, 5, 7}, {7, 5, 6, 2}});
    const ConformityResult r = check_conformity(m);
    EXPECT_EQ(r.mDimension, 2);
    EXPECT_EQ(r.mNumTJunctions, 1);
    EXPECT_EQ(r.mNumHangingNodes, 0);
    EXPECT_EQ(r.mNumPartialOverlaps, 2);
    const std::vector<std::int64_t> flags = cc_flags(r);
    for (std::size_t i = 0; i < flags.size(); ++i)
        EXPECT_EQ(flags[i], i == 7 ? 2 : 0) << i;
}

TEST(CheckConformity, CoincidentFacesWithOtherIdsAreUnwelded) {
    // Two unit cubes side by side, each with its own eight points.
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<std::int64_t>> hexes(2);
    for (int c = 0; c < 2; ++c) {
        const std::vector<std::vector<double>> corners = {
            {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
        for (const auto& q : corners) {
            hexes[static_cast<std::size_t>(c)].push_back(static_cast<std::int64_t>(pts.size()));
            pts.push_back({q[0] + c, q[1], q[2]});
        }
    }
    const ConformityResult r =
        check_conformity(mt::make_mesh(std::move(pts), "hexahedron", std::move(hexes)));
    EXPECT_EQ(r.mNumUnweldedFacets, 1);
    EXPECT_EQ(r.mNumHangingNodes, 0);
    EXPECT_EQ(r.mNumPartialOverlaps, 0);
    EXPECT_FALSE(r.mConforming);
    EXPECT_EQ(r.mInterfaceMode, "weld");
}

TEST(CheckConformity, ShiftedNeighboursOverlapPartly) {
    // A unit cube and one shifted by half along y, each with its own points:
    // the faces at x = 1 overlap over half their area, and two corners of each
    // lie on the other's face edges.
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<std::int64_t>> hexes(2);
    for (int c = 0; c < 2; ++c) {
        const std::vector<std::vector<double>> corners = {
            {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
        for (const auto& q : corners) {
            hexes[static_cast<std::size_t>(c)].push_back(static_cast<std::int64_t>(pts.size()));
            pts.push_back({q[0] + c, q[1] + 0.5 * c, q[2]});
        }
    }
    const ConformityResult r =
        check_conformity(mt::make_mesh(std::move(pts), "hexahedron", std::move(hexes)));
    EXPECT_EQ(r.mNumPartialOverlaps, 1);
    EXPECT_EQ(r.mNumHangingNodes, 4);
    EXPECT_EQ(r.mInterfaceMode, "proximity");
}

TEST(CheckConformity, QuadraticMidEdgeNodesAreNotHanging) {
    // Two quad8 sharing an edge and its mid-edge node: every mid-edge node
    // lies on an edge of its own cell, which never flags it.
    const Mesh m = mt::make_mesh({{0, 0, 0},
                                  {1, 0, 0},
                                  {1, 1, 0},
                                  {0, 1, 0},
                                  {2, 0, 0},
                                  {2, 1, 0},
                                  {0.5, 0, 0},
                                  {1, 0.5, 0},
                                  {0.5, 1, 0},
                                  {0, 0.5, 0},
                                  {1.5, 0, 0},
                                  {2, 0.5, 0},
                                  {1.5, 1, 0}},
                                 "quad8", {{0, 1, 2, 3, 6, 7, 8, 9}, {1, 4, 5, 2, 10, 11, 12, 7}});
    const ConformityResult r = check_conformity(m);
    EXPECT_TRUE(r.mConforming);
    EXPECT_EQ(r.mNumTJunctions, 0);
    EXPECT_EQ(r.mNumBoundaryFacets, 6);
}

TEST(CheckConformity, LinesAloneAreNotExamined) {
    const ConformityResult r = check_conformity(mt::line_mesh());
    EXPECT_EQ(r.mDimension, 0);
    EXPECT_TRUE(r.mConforming);
    ASSERT_EQ(r.mWarnings.size(), 1u);
}

TEST(CheckConformity, RefusesABadTolerance) {
    ConformityOptions o;
    o.mTolerance = -1.0;
    EXPECT_THROW(check_conformity(mt::hex_mesh(), o), std::invalid_argument);
    o.mTolerance = std::numeric_limits<double>::infinity();
    EXPECT_THROW(check_conformity(mt::hex_mesh(), o), std::invalid_argument);
}
