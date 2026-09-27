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

// Sequence time resampling (resample_plan, blend_steps, the driver) and the
// agglomerate follow-ups (sphericity gate, coplanar face merging) -- roadmap
// §5, v16.25.0.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/operations/blend.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/registry.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

NDArray f64(const std::vector<double>& rV) {
    NDArray a = NDArray::Uninit(DType::Float64, {rV.size()});
    std::copy(rV.begin(), rV.end(), a.As<double>());
    return a;
}

Mesh tri_with(double U) {
    Mesh m = mt::tri_mesh();
    std::vector<double> u(m.NumPoints(), U);
    m.AddPointData("u", f64(u));
    NDArray id = NDArray::Uninit(DType::Int32, {m.NumPoints()});
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        id.As<std::int32_t>()[i] = static_cast<std::int32_t>(U);
    m.AddPointData("id", std::move(id));
    return m;
}

// An nx x ny x nz grid of unit hexahedra.
Mesh hex_grid(int Nx, int Ny, int Nz) {
    auto id = [&](int i, int j, int k) {
        return static_cast<std::int64_t>((k * (Ny + 1) + j) * (Nx + 1) + i);
    };
    std::vector<std::vector<double>> pts;
    for (int k = 0; k <= Nz; ++k)
        for (int j = 0; j <= Ny; ++j)
            for (int i = 0; i <= Nx; ++i)
                pts.push_back({double(i), double(j), double(k)});
    std::vector<std::vector<std::int64_t>> hexes;
    for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                hexes.push_back({id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k),
                                 id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1),
                                 id(i, j + 1, k + 1)});
    return mt::make_mesh(pts, "hexahedron", hexes);
}

// Every face (as a sorted node set) of every polyhedron; a conforming mesh
// has each interior face exactly twice.
std::map<std::vector<std::int64_t>, int> face_uses(const Mesh& rMesh) {
    std::map<std::vector<std::int64_t>, int> uses;
    for (const auto cb : rMesh.CellRange()) {
        if (!cb.IsPolyhedron())
            continue;
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                const auto face = cb.Face(c, f);
                std::vector<std::int64_t> key(face.first, face.first + face.second);
                std::sort(key.begin(), key.end());
                ++uses[key];
            }
    }
    return uses;
}

class RsTempDir {
public:
    RsTempDir() {
        static std::atomic<unsigned> counter{0};
        mPath = std::filesystem::temp_directory_path() /
                ("meshio_rs_" + std::to_string(counter++) + "_" +
                 std::to_string(::testing::UnitTest::GetInstance()->random_seed()));
        std::filesystem::create_directories(mPath);
    }
    ~RsTempDir() {
        std::error_code ec;
        std::filesystem::remove_all(mPath, ec);
    }
    std::string operator/(const std::string& rName) const { return (mPath / rName).string(); }

private:
    std::filesystem::path mPath;
};

}  // namespace

// ------------------------------------------------------------ resample_plan

TEST(ResamplePlan, LinearNearestPrevious) {
    const std::vector<double> t = {0.0, 1.0, 3.0};
    auto lin =
        resample_plan(t, {0.25, 1.0, 2.5}, ResampleMethod::Linear, ResampleExtrapolate::Error);
    ASSERT_EQ(lin.size(), 3u);
    EXPECT_EQ(lin[0].mLo, 0u);
    EXPECT_EQ(lin[0].mHi, 1u);
    EXPECT_DOUBLE_EQ(lin[0].mWeight, 0.25);
    EXPECT_EQ(lin[1].mLo, 1u);  // an exact hit takes the step alone
    EXPECT_EQ(lin[1].mHi, 1u);
    EXPECT_DOUBLE_EQ(lin[2].mWeight, 0.75);
    auto near =
        resample_plan(t, {0.5, 1.9, 2.1}, ResampleMethod::Nearest, ResampleExtrapolate::Error);
    EXPECT_EQ(near[0].mLo, 0u);  // a tie goes to the earlier step
    EXPECT_EQ(near[1].mLo, 1u);
    EXPECT_EQ(near[2].mLo, 2u);
    auto prev = resample_plan(t, {2.9}, ResampleMethod::Previous, ResampleExtrapolate::Error);
    EXPECT_EQ(prev[0].mLo, 1u);
}

TEST(ResamplePlan, ExtrapolationAndValidation) {
    const std::vector<double> t = {0.0, 1.0};
    EXPECT_THROW(resample_plan(t, {1.5}, ResampleMethod::Linear, ResampleExtrapolate::Error),
                 std::invalid_argument);
    auto c = resample_plan(t, {-1.0, 2.0}, ResampleMethod::Linear, ResampleExtrapolate::Clamp);
    EXPECT_EQ(c[0].mLo, 0u);
    EXPECT_EQ(c[1].mLo, 1u);
    EXPECT_THROW(
        resample_plan({1.0, 1.0}, {1.0}, ResampleMethod::Linear, ResampleExtrapolate::Error),
        std::invalid_argument);
    EXPECT_THROW(resample_plan({}, {1.0}, ResampleMethod::Linear, ResampleExtrapolate::Error),
                 std::invalid_argument);
    const auto r = resample_times_range(0.0, 1.0, 0.1);
    ASSERT_EQ(r.size(), 11u);
    EXPECT_NEAR(r.back(), 1.0, 1e-12);
    EXPECT_THROW(resample_times_range(0.0, 1.0, 0.0), std::invalid_argument);
    EXPECT_EQ(resample_method_from_name("previous"), ResampleMethod::Previous);
    EXPECT_THROW(resample_method_from_name("cubic"), std::invalid_argument);
}

// -------------------------------------------------------------- blend_steps

TEST(BlendSteps, BlendsFloatsAndTakesIntegersFromTheNearerStep) {
    const Mesh out = blend_steps(tri_with(0.0), tri_with(8.0), 0.25);
    for (std::size_t i = 0; i < out.NumPoints(); ++i) {
        EXPECT_DOUBLE_EQ(out.PointData("u").As<double>()[i], 2.0);
        EXPECT_EQ(out.PointData("id").As<std::int32_t>()[i], 0);
    }
    const Mesh late = blend_steps(tri_with(0.0), tri_with(8.0), 0.75);
    EXPECT_EQ(late.PointData("id").As<std::int32_t>()[0], 8);
}

TEST(BlendSteps, RefusesDifferentTopologies) {
    EXPECT_THROW(blend_steps(tri_with(0.0), mt::quad_mesh(), 0.5), std::invalid_argument);
    Mesh other = tri_with(1.0);
    other.AddPointData("extra", f64(std::vector<double>(other.NumPoints(), 0.0)));
    EXPECT_THROW(blend_steps(tri_with(0.0), other, 0.5), std::invalid_argument);
}

// ---------------------------------------------------- the resampling driver

TEST(SequenceResample, WritesOneFilePerTargetTime) {
    RsTempDir dir;
    SequencePipeline p;
    for (int i = 0; i < 3; ++i) {
        const std::string path = dir / ("step_" + std::to_string(i) + ".vtu");
        registry_writers().at("vtu")(path, tri_with(10.0 * i));
        p.mInput.mPaths.push_back(path);
    }
    p.mInput.mTimes = {0.0, 1.0, 2.0};
    p.mOutput.mPath = dir / "out_{index}.vtu";
    SequenceResample rs;
    rs.mTimes = {0.5, 1.25, 2.0};
    p.mResample = rs;
    run_sequence_pipeline(p);
    const double want[] = {5.0, 12.5, 20.0};
    for (int k = 0; k < 3; ++k) {
        const Mesh m = registry_readers().at("vtu")(dir / ("out_" + std::to_string(k) + ".vtu"));
        EXPECT_NEAR(detail::read_double(m.PointData("u"), 0), want[k], 1e-12) << k;
    }
    // A target outside the range is an error unless clamped.
    p.mResample->mTimes = {3.0};
    EXPECT_THROW(run_sequence_pipeline(p), std::invalid_argument);
    p.mResample->mExtrapolate = ResampleExtrapolate::Clamp;
    EXPECT_NO_THROW(run_sequence_pipeline(p));
}

// ------------------------------------------------ agglomerate: shape gate

TEST(AgglomerateGate, SphericityKeepsGroupsCompact) {
    const Mesh bar = hex_grid(8, 1, 1);
    AgglomerateOptions plain;
    plain.mTargetGroupSize = 8;
    const AgglomerateResult a = agglomerate(bar, plain);
    EXPECT_EQ(a.mMesh.Cells(0).NumCells(), 1u);  // the whole bar is one group
    AgglomerateOptions gated = plain;
    gated.mMinSphericity = 0.7;
    const AgglomerateResult b = agglomerate(bar, gated);
    EXPECT_GT(b.mMesh.Cells(0).NumCells(), 1u);
    EXPECT_GT(b.mNumRejected, 0);
    EXPECT_NEAR(compute_stats(b.mMesh).mSignedVolume, 8.0, 1e-12);
    gated.mMinSphericity = 1.5;
    EXPECT_THROW(agglomerate(bar, gated), std::invalid_argument);
}

// ---------------------------------------------- agglomerate: coplanar merge

TEST(AgglomerateCoplanar, FusesTheFlatFacesOfOneGroup) {
    AgglomerateOptions o;
    o.mTargetGroupSize = 4;
    o.mMergeCoplanarFaces = true;
    const AgglomerateResult r = agglomerate(hex_grid(2, 2, 1), o);
    ASSERT_EQ(r.mMesh.Cells(0).NumCells(), 1u);
    // A 2x2x1 box: top and bottom 4 -> 1, the four sides 2 -> 1 each.
    EXPECT_EQ(r.mMesh.Cells(0).NumFaces(0), 6u);
    EXPECT_EQ(r.mNumFacesMerged, 3 + 3 + 4);
    EXPECT_NEAR(compute_stats(r.mMesh).mSignedVolume, 4.0, 1e-12);
    for (const auto& kv : face_uses(r.mMesh))
        EXPECT_EQ(kv.second, 1);
}

TEST(AgglomerateCoplanar, SharedFacesStayConforming) {
    AgglomerateOptions o;
    o.mTargetGroupSize = 4;
    o.mMergeCoplanarFaces = true;
    const Mesh m = hex_grid(4, 2, 2);
    const AgglomerateResult r = agglomerate(m, o);
    EXPECT_GT(r.mMesh.Cells(0).NumCells(), 1u);
    EXPECT_GT(r.mNumFacesMerged, 0);
    EXPECT_NEAR(compute_stats(r.mMesh).mSignedVolume, 16.0, 1e-12);
    // Every face is on the boundary (used once) or shared by exactly two
    // groups with the identical fused polygon (used twice).
    for (const auto& kv : face_uses(r.mMesh))
        EXPECT_TRUE(kv.second == 1 || kv.second == 2);
    AgglomerateOptions off = o;
    off.mMergeCoplanarFaces = false;
    EXPECT_EQ(agglomerate(m, off).mNumFacesMerged, 0);
    o.mCoplanarAngleDeg = 95.0;
    EXPECT_THROW(agglomerate(m, o), std::invalid_argument);
}
