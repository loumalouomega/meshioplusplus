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
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "analytic_surfaces.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

struct FpRow {
    std::int64_t mType;
    std::vector<double> mParams;
    double mRms;
};

std::vector<FpRow> rows(const Mesh& rMesh) {
    const NDArray& t = rMesh.FieldData("fit:patch_type");
    const NDArray& p = rMesh.FieldData("fit:patch_params");
    const NDArray& r = rMesh.FieldData("fit:patch_rms");
    std::vector<FpRow> out;
    for (std::size_t i = 0; i < t.Size(); ++i)
        out.push_back({t.As<std::int64_t>()[i],
                       std::vector<double>(p.As<double>() + i * 10, p.As<double>() + i * 10 + 10),
                       r.As<double>()[i]});
    return out;
}

FpRow fit_whole(const Mesh& rMesh, double Tol = 1e-3) {
    FitPrimitivesOptions o;
    o.mPatchData = "";
    o.mTolerance = Tol;
    const std::vector<FpRow> r = rows(fit_primitives(rMesh, o).mMesh);
    EXPECT_EQ(r.size(), 1u);
    return r.empty() ? FpRow{} : r[0];
}

// splitmix64: a portable seeded generator, so the noise is the same everywhere.
struct SplitMix {
    std::uint64_t mState;
    double Next() {
        std::uint64_t z = (mState += 0x9e3779b97f4a7c15ull);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
        z ^= z >> 31;
        return static_cast<double>(z >> 11) * (1.0 / 9007199254740992.0) * 2.0 - 1.0;
    }
};

Mesh jitter(const Mesh& rMesh, double Amp, std::uint64_t Seed) {
    SplitMix g{Seed};
    NDArray pts = NDArray::Uninit(DType::Float64, rMesh.Points().Shape());
    for (std::size_t i = 0; i < pts.Size(); ++i)
        pts.As<double>()[i] = rMesh.Points().As<double>()[i] + Amp * g.Next();
    Mesh out;
    out.AssignPoints(std::move(pts));
    NDArray conn = NDArray::Uninit(DType::Int64, rMesh.Cells(0).Conn().Shape());
    std::copy(rMesh.Cells(0).Conn().As<std::int64_t>(),
              rMesh.Cells(0).Conn().As<std::int64_t>() + conn.Size(), conn.As<std::int64_t>());
    out.AddCellBlock("quad", std::move(conn));
    return out;
}

}  // namespace

TEST(FitPrimitives, BoxGivesSixPlanes) {
    const Mesh seg = segment(analytic::box_skin(3)).mMesh;
    const FitPrimitivesResult r = fit_primitives(seg);
    EXPECT_EQ(r.mNumPatches, 6);
    EXPECT_EQ(r.mNumPlane, 6);
    for (const FpRow& row : rows(r.mMesh)) {
        EXPECT_EQ(row.mType, 1);
        EXPECT_LT(row.mRms, 1e-12);
        // A unit axis normal pointing out of the cube.
        double len = 0.0, dot = 0.0;
        for (int k = 0; k < 3; ++k) {
            len += row.mParams[3 + k] * row.mParams[3 + k];
            dot += row.mParams[3 + k] * (row.mParams[k] - 0.5);
        }
        EXPECT_NEAR(len, 1.0, 1e-12);
        EXPECT_NEAR(dot, 0.5, 1e-12);
    }
    const NDArray& t = r.mMesh.CellData(kFitTypeName, 0);
    for (std::size_t c = 0; c < t.Size(); ++c)
        EXPECT_EQ(t.As<std::int64_t>()[c], 1);
}

TEST(FitPrimitives, RecoversCylinderConeSphereTorus) {
    const FpRow cyl = fit_whole(analytic::cylinder(1.25, 2.0));
    EXPECT_EQ(cyl.mType, 3);
    EXPECT_NEAR(cyl.mParams[6], 1.25, 1e-9);
    EXPECT_NEAR(std::fabs(cyl.mParams[5]), 1.0, 1e-9);
    EXPECT_NEAR(cyl.mParams[0], 0.0, 1e-9);
    EXPECT_NEAR(cyl.mParams[1], 0.0, 1e-9);

    const FpRow cone = fit_whole(analytic::cone(20.0));
    EXPECT_EQ(cone.mType, 4);
    EXPECT_NEAR(cone.mParams[8], 20.0, 1e-7);
    for (int k = 0; k < 3; ++k)
        EXPECT_NEAR(cone.mParams[k], 0.0, 1e-7);  // the apex
    EXPECT_NEAR(cone.mParams[5], 1.0, 1e-9);      // opening along +z

    const FpRow sph = fit_whole(analytic::sphere_band(1.5));
    EXPECT_EQ(sph.mType, 2);
    EXPECT_NEAR(sph.mParams[6], 1.5, 1e-9);

    const FpRow tor = fit_whole(analytic::torus(3.0, 1.0));
    EXPECT_EQ(tor.mType, 5);
    EXPECT_NEAR(tor.mParams[6], 3.0, 1e-9);
    EXPECT_NEAR(tor.mParams[7], 1.0, 1e-9);
    EXPECT_NEAR(tor.mParams[5], 1.0, 1e-9);
}

TEST(FitPrimitives, NoiseRaisesTheResidualNotTheType) {
    const Mesh noisy = jitter(analytic::cylinder(1.0, 2.0), 1e-3, 7);
    const FpRow clean = fit_whole(analytic::cylinder(1.0, 2.0));
    const FpRow row = fit_whole(noisy, 5e-3);
    EXPECT_EQ(row.mType, 3);
    EXPECT_GT(row.mRms, clean.mRms);
    EXPECT_LT(row.mRms, 2e-3);
    EXPECT_NEAR(row.mParams[6], 1.0, 2e-3);
}

TEST(FitPrimitives, SaddleIsReportedUnfitted) {
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<std::int64_t>> quads;
    const int n = 8;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) {
            const double x = -1.0 + 2.0 * i / n, y = -1.0 + 2.0 * j / n;
            pts.push_back({x, y, 0.5 * (x * x - y * y)});
        }
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            const std::int64_t p = j * (n + 1) + i;
            quads.push_back({p, p + 1, p + n + 2, p + n + 1});
        }
    FitPrimitivesOptions o;
    o.mPatchData = "";
    o.mTolerance = 1e-4;
    const FitPrimitivesResult r = fit_primitives(mt::make_mesh(pts, "quad", quads), o);
    EXPECT_EQ(r.mNumUnfitted, 1);
    const FpRow row = rows(r.mMesh)[0];
    EXPECT_EQ(row.mType, 0);
    EXPECT_GT(row.mRms, 1e-4);  // the best attempt's residual, reported
    EXPECT_TRUE(std::isnan(row.mParams[0]));
}

TEST(FitPrimitives, TypesRestrictTheCandidates) {
    FitPrimitivesOptions o;
    o.mPatchData = "";
    o.mTypes = 1u << static_cast<int>(FitPrimitiveType::Plane);
    const FitPrimitivesResult r = fit_primitives(analytic::cylinder(), o);
    EXPECT_EQ(r.mNumUnfitted, 1);
    o.mTypes |= 1u << static_cast<int>(FitPrimitiveType::Cylinder);
    EXPECT_EQ(fit_primitives(analytic::cylinder(), o).mNumCylinder, 1);
}

TEST(FitPrimitives, RefusesBadOptionsAndMissingLabels) {
    EXPECT_THROW(fit_primitives(analytic::cylinder()), std::invalid_argument);  // no labels
    FitPrimitivesOptions o;
    o.mPatchData = "";
    o.mTolerance = 0.0;
    EXPECT_THROW(fit_primitives(analytic::cylinder(), o), std::invalid_argument);
    o = {};
    o.mPatchData = "";
    o.mTypes = 0;
    EXPECT_THROW(fit_primitives(analytic::cylinder(), o), std::invalid_argument);
    o = {};
    o.mPatchData = "";
    o.mMinInlierFraction = 1.5;
    EXPECT_THROW(fit_primitives(analytic::cylinder(), o), std::invalid_argument);
    o = {};
    o.mPatchData = "";
    EXPECT_THROW(fit_primitives(mt::tet_mesh(), o), std::invalid_argument);
}

TEST(FitPrimitives, ResultIsStableAcrossRuns) {
    const Mesh seg = segment(analytic::box_skin(4)).mMesh;
    const Mesh a = fit_primitives(seg).mMesh;
    const Mesh b = fit_primitives(seg).mMesh;
    const NDArray& pa = a.FieldData("fit:patch_params");
    const NDArray& pb = b.FieldData("fit:patch_params");
    ASSERT_EQ(pa.Size(), pb.Size());
    EXPECT_EQ(std::memcmp(pa.Data(), pb.Data(), pa.Nbytes()), 0);
}
