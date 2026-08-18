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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/operations/decimate_volume.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::compute_stats;
using meshioplusplus::decimate_volume;
using meshioplusplus::DecimatePlacement;
using meshioplusplus::DecimateVolumeOptions;
using meshioplusplus::DecimateVolumeResult;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::RegionKind;

// A unit cube split into 6 positively-oriented tetrahedra sharing the main
// diagonal 0-6 (the standard FEM decomposition). Every vertex is a cube
// corner, so every vertex touches the boundary skin -- exercises the
// boundary-quadric/link-condition machinery on every candidate edge.
Mesh cube6_mesh() {
    return mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "tetra",
        {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}});
}

// Two unit cubes stacked along z, each split the same way, sharing the z=1
// face -- gives genuinely interior structure (the shared face's diagonal
// edges are not cube corners of the outer 2x1x1 box) alongside the boundary.
Mesh cube6x2_mesh() {
    Mesh lo = cube6_mesh();
    Mesh hi = mt::make_mesh(
        {{0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {0, 0, 2}, {1, 0, 2}, {1, 1, 2}, {0, 1, 2}},
        "tetra",
        {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}});
    // Weld the shared z=1 face by hand (4 point pairs), avoiding a dependency
    // on `merge`/`clean` from this test.
    Mesh out;
    std::vector<std::vector<double>> pts;
    for (std::size_t i = 0; i < 8; ++i)
        pts.push_back({lo.Points().As<double>()[i * 3], lo.Points().As<double>()[i * 3 + 1],
                       lo.Points().As<double>()[i * 3 + 2]});
    // hi's points 0..3 (z=1) map onto lo's 4..7; hi's 4..7 (z=2) are new ids 8..11.
    std::vector<std::int64_t> hi_map = {4, 5, 6, 7, 8, 9, 10, 11};
    for (std::size_t i = 4; i < 8; ++i)
        pts.push_back({hi.Points().As<double>()[i * 3], hi.Points().As<double>()[i * 3 + 1],
                       hi.Points().As<double>()[i * 3 + 2]});
    out.AssignPoints(mt::points_from(pts));
    std::vector<std::vector<std::int64_t>> cells;
    const std::int64_t* lc = lo.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t c = 0; c < 6; ++c)
        cells.push_back({lc[c * 4], lc[c * 4 + 1], lc[c * 4 + 2], lc[c * 4 + 3]});
    const std::int64_t* hc = hi.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t c = 0; c < 6; ++c)
        cells.push_back({hi_map[static_cast<std::size_t>(hc[c * 4])],
                         hi_map[static_cast<std::size_t>(hc[c * 4 + 1])],
                         hi_map[static_cast<std::size_t>(hc[c * 4 + 2])],
                         hi_map[static_cast<std::size_t>(hc[c * 4 + 3])]});
    out.AddCellBlock("tetra", mt::conn_from(cells));
    return out;
}

std::size_t alive_tets(const Mesh& rMesh) {
    std::size_t n = 0;
    for (const auto cb : rMesh.CellRange())
        n += cb.NumCells();
    return n;
}

void assert_no_repeated_corners(const Mesh& rMesh) {
    for (const auto cb : rMesh.CellRange()) {
        const NDArray& conn = cb.Conn();
        const std::int64_t* c = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < cb.NumCells(); ++i) {
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b)
                    ASSERT_NE(c[i * 4 + a], c[i * 4 + b])
                        << "tet " << i << " has a repeated corner";
        }
    }
}

// --- the topology oracle, written before any of the guard/greedy-loop code
// above was trusted: collapse one edge at a time (via successive
// target_cells = current-1 calls, each of which the greedy loop can only
// satisfy via a single commit) and check invariants no quality objective can
// accidentally satisfy by luck.
//
// Total volume is deliberately NOT one of them: unlike refine/subdivide,
// which partition the same geometry, decimation moves vertices, so surviving
// tets legitimately change volume as their corner is repositioned -- only
// TOPOLOGY (manifoldness, no repeated corner) and the inversion guard's own
// promise (no positively-oriented tet flips sign) are true invariants here.
TEST(DecimateVolume, TopologyOracleSurvivesEveryCandidateCollapse) {
    for (Mesh mesh : {cube6_mesh(), cube6x2_mesh()}) {
        std::size_t guard = 0;
        while (alive_tets(mesh) > 1 && guard++ < 64) {
            DecimateVolumeOptions o;
            o.mTargetCells = static_cast<std::int64_t>(alive_tets(mesh)) - 1;
            o.mPreserveBoundary = false;
            // Every vertex of these fixtures is a cube corner (a sharp
            // geometric feature); the feature-preservation default would pin
            // all of them and make the loop vacuous.
            o.mPreserveFeatures = false;
            DecimateVolumeResult r = decimate_volume(mesh, o);
            if (alive_tets(r.mMesh) == alive_tets(mesh))
                break;  // queue exhausted: no valid collapse left, not a bug

            assert_no_repeated_corners(r.mMesh);
            const meshioplusplus::detail::GlobalFaces gf =
                meshioplusplus::detail::build_global_faces(r.mMesh);
            EXPECT_EQ(gf.mNumNonManifold, 0);
            // Every input tet was positively oriented; the inversion guard's
            // entire promise is that no surviving tet flips sign.
            EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);
            EXPECT_EQ(alive_tets(r.mMesh), alive_tets(mesh) - r.mTetsRemoved);

            mesh = std::move(r.mMesh);
        }
    }
}

// --- validation --------------------------------------------------------------

TEST(DecimateVolume, RejectsOutOfScopeInputsByName) {
    EXPECT_THROW(decimate_volume(mt::hex_mesh()), std::invalid_argument);
    EXPECT_THROW(decimate_volume(mt::wedge_mesh()), std::invalid_argument);
    EXPECT_THROW(decimate_volume(mt::tri_mesh()), std::invalid_argument);
}

TEST(DecimateVolume, ValidatesTheStoppingCriteria) {
    Mesh mesh = cube6_mesh();
    EXPECT_THROW(decimate_volume(mesh, {}), std::invalid_argument);
    DecimateVolumeOptions both;
    both.mTargetRatio = 0.5;
    both.mTargetCells = 3;
    EXPECT_THROW(decimate_volume(mesh, both), std::invalid_argument);
    DecimateVolumeOptions bad_ratio;
    bad_ratio.mTargetRatio = 1.5;
    EXPECT_THROW(decimate_volume(mesh, bad_ratio), std::invalid_argument);
}

TEST(DecimateVolume, ValidatesTheFrozenMaskSize) {
    Mesh mesh = cube6_mesh();
    DecimateVolumeOptions o;
    o.mTargetCells = 1;
    o.mFrozen = std::vector<std::uint8_t>(3, 0);
    EXPECT_THROW(decimate_volume(mesh, o), std::invalid_argument);
}

// --- core behaviour ------------------------------------------------------------

TEST(DecimateVolume, TargetCellsIsHonoured) {
    Mesh mesh = cube6_mesh();
    DecimateVolumeOptions o;
    o.mTargetCells = 1;
    o.mPreserveFeatures = false;  // every vertex here is a sharp cube corner
    DecimateVolumeResult r = decimate_volume(mesh, o);
    EXPECT_LE(alive_tets(r.mMesh), 1u);
    EXPECT_GT(r.mTetsRemoved, 0);
}

TEST(DecimateVolume, NoSurvivingTetInvertsUnderDefaultOptions) {
    // Every input tet is positively oriented; the inversion guard's promise
    // is that decimation cannot flip one, even though (unlike refine's exact
    // partition) the surviving tets' individual volumes legitimately shift as
    // their collapsed corner is repositioned.
    Mesh mesh = cube6_mesh();
    DecimateVolumeOptions o;
    o.mTargetRatio = 0.4;
    o.mPreserveFeatures = false;
    DecimateVolumeResult r = decimate_volume(mesh, o);
    EXPECT_GT(alive_tets(r.mMesh), 0u);
    EXPECT_EQ(compute_stats(r.mMesh).mNumInverted, 0);
}

TEST(DecimateVolume, PreserveBoundaryPinsEveryVertexOfAnAllBoundaryMesh) {
    // Every vertex of cube6_mesh() touches the boundary, so pinning the
    // boundary leaves nothing collapsible: the result equals the input.
    Mesh mesh = cube6_mesh();
    DecimateVolumeOptions o;
    o.mTargetCells = 1;
    o.mPreserveBoundary = true;
    DecimateVolumeResult r = decimate_volume(mesh, o);
    EXPECT_EQ(alive_tets(r.mMesh), alive_tets(mesh));
    EXPECT_EQ(r.mTetsRemoved, 0);
}

TEST(DecimateVolume, FrozenMaskPreventsCollapse) {
    Mesh mesh = cube6_mesh();
    DecimateVolumeOptions o;
    o.mTargetCells = 1;
    o.mPreserveFeatures = false;
    o.mFrozen = std::vector<std::uint8_t>(mesh.NumPoints(), 1);
    DecimateVolumeResult r = decimate_volume(mesh, o);
    EXPECT_EQ(alive_tets(r.mMesh), alive_tets(mesh));
}

TEST(DecimateVolume, ResultIsStableAcrossRepeatedRuns) {
    Mesh mesh = cube6x2_mesh();
    DecimateVolumeOptions o;
    o.mTargetRatio = 0.3;
    o.mPreserveFeatures = false;
    DecimateVolumeResult a = decimate_volume(mesh, o);
    DecimateVolumeResult b = decimate_volume(mesh, o);
    EXPECT_EQ(alive_tets(a.mMesh), alive_tets(b.mMesh));
    ASSERT_EQ(a.mMesh.Points().Size(), b.mMesh.Points().Size());
    const double* pa = a.mMesh.Points().As<double>();
    const double* pb = b.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < a.mMesh.Points().Size(); ++i)
        EXPECT_EQ(pa[i], pb[i]);
}

TEST(DecimateVolume, PlacementModesAllProduceValidMeshes) {
    for (DecimatePlacement p :
         {DecimatePlacement::Optimal, DecimatePlacement::Midpoint, DecimatePlacement::Endpoint}) {
        Mesh mesh = cube6x2_mesh();
        DecimateVolumeOptions o;
        o.mTargetRatio = 0.5;
        o.mPreserveFeatures = false;
        o.mPlacement = p;
        DecimateVolumeResult r = decimate_volume(mesh, o);
        assert_no_repeated_corners(r.mMesh);
    }
}

TEST(DecimateVolume, DataIsCarriedAndCellMapsAreValid) {
    Mesh mesh = cube6_mesh();
    std::vector<double> pd(mesh.NumPoints());
    for (std::size_t i = 0; i < pd.size(); ++i)
        pd[i] = static_cast<double>(i);
    NDArray pa = NDArray::Uninit(meshioplusplus::DType::Float64, {mesh.NumPoints()});
    std::copy(pd.begin(), pd.end(), pa.As<double>());
    mesh.AddPointData("scalar", std::move(pa));

    NDArray cd = NDArray::Uninit(meshioplusplus::DType::Int64, {mesh.Cells(0).NumCells()});
    for (std::size_t i = 0; i < mesh.Cells(0).NumCells(); ++i)
        cd.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    mesh.AddCellData("tag", {std::move(cd)});

    DecimateVolumeOptions o;
    o.mTargetCells = 3;
    o.mPreserveFeatures = false;
    DecimateVolumeResult r = decimate_volume(mesh, o);
    EXPECT_TRUE(r.mMesh.HasPointData("scalar"));
    EXPECT_TRUE(r.mMesh.HasCellData("tag"));
    ASSERT_EQ(r.mCellMaps.size(), 1u);
    const std::int64_t* map = r.mCellMaps[0].As<std::int64_t>();
    for (std::size_t i = 0; i < r.mCellMaps[0].Size(); ++i)
        EXPECT_LT(map[i], static_cast<std::int64_t>(alive_tets(r.mMesh)));
}

TEST(DecimateVolume, RegionsSurviveTheNonInjectiveDirectMap) {
    Mesh mesh = cube6_mesh();
    meshioplusplus::Region region;
    region.mName = "all";
    region.mKind = RegionKind::Cell;
    NDArray entries = NDArray::Uninit(meshioplusplus::DType::Int64, {6});
    for (std::int64_t i = 0; i < 6; ++i)
        entries.As<std::int64_t>()[i] = i;
    region.mEntries = std::move(entries);
    mesh.AddRegion(region);

    DecimateVolumeOptions o;
    o.mTargetCells = 1;
    o.mPreserveFeatures = false;
    DecimateVolumeResult r = decimate_volume(mesh, o);
    ASSERT_TRUE(r.mMesh.HasRegion("all", RegionKind::Cell));
    const meshioplusplus::Region& out = r.mMesh.Region(r.mMesh.FindRegion("all", RegionKind::Cell));
    // Every original cell (whether it survived or was welded into another)
    // maps into the same region: the count must equal the survivor count,
    // deduplicated by Region::Canonicalize.
    EXPECT_EQ(static_cast<std::size_t>(out.mEntries.Size()), alive_tets(r.mMesh));
}

}  // namespace
