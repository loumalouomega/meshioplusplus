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

TEST(FitPrimitives, ATwoRimBoreIsACylinderNotASphere) {
    // One ring of quads: every point lies on the two rims, and two coaxial
    // circles lie on a sphere as well as on the cylinder.
    const FpRow band = fit_whole(analytic::cylinder(1.6, 0.5, 1));
    EXPECT_EQ(band.mType, 3);
    EXPECT_NEAR(band.mParams[6], 1.6, 1e-9);
    // A band of two rims of unequal radius is the cone, not a sphere either.
    const FpRow cone_band = fit_whole(analytic::cone(30.0, 1));
    EXPECT_EQ(cone_band.mType, 4);
    EXPECT_NEAR(cone_band.mParams[8], 30.0, 1e-6);
    // A real sphere band keeps its type.
    EXPECT_EQ(fit_whole(analytic::sphere_band(1.5)).mType, 2);
}

namespace {

double fd_value(const Mesh& rMesh, const char* pName) {
    return rMesh.FieldData(pName).As<double>()[0];
}

Mesh fp_fit_whole_mesh(const Mesh& rMesh, double Tol = 1e-3) {
    FitPrimitivesOptions o;
    o.mPatchData = "";
    o.mTolerance = Tol;
    return fit_primitives(rMesh, o).mMesh;
}

// Quads on the half cylinder x = R cos a, y = R sin a, a in [0, 180] degrees.
Mesh half_cylinder(double R, double H, int N, int M) {
    std::vector<std::vector<double>> pts;
    for (int i = 0; i <= M; ++i)
        for (int k = 0; k <= N; ++k) {
            const double a = analytic::kPi * k / N;
            pts.push_back({R * std::cos(a), R * std::sin(a), H * i / M});
        }
    std::vector<std::vector<std::int64_t>> quads;
    for (int i = 0; i < M; ++i)
        for (int k = 0; k < N; ++k) {
            const std::int64_t a = i * (N + 1) + k, b = (i + 1) * (N + 1) + k;
            quads.push_back({a, a + 1, b + 1, b});
        }
    return mt::make_mesh(pts, "quad", quads);
}

}  // namespace

TEST(FitPrimitivesForm, HalfCylinderArcAndZeroFormError) {
    const Mesh out = fp_fit_whole_mesh(half_cylinder(1.5, 2.0, 24, 4));
    EXPECT_EQ(out.FieldData("fit:patch_type").As<std::int64_t>()[0], 3);
    EXPECT_NEAR(fd_value(out, kFitPatchArcName), 180.0, 1e-9);
    EXPECT_NEAR(fd_value(out, kFitPatchFormErrorName), 0.0, 1e-9);
    EXPECT_NEAR(fd_value(out, kFitPatchLengthName), 2.0, 1e-9);
    EXPECT_NEAR(fd_value(out, kFitPatchDiameterMinName), 3.0, 1e-9);
    EXPECT_NEAR(fd_value(out, kFitPatchDiameterMaxName), 3.0, 1e-9);
    EXPECT_TRUE(std::isnan(fd_value(out, kFitPatchExtentName)));
}

TEST(FitPrimitivesForm, ClosedCylinderSpansTheFullCircle) {
    const Mesh out = fp_fit_whole_mesh(analytic::cylinder(1.0, 2.0, 6));
    EXPECT_EQ(out.FieldData("fit:patch_type").As<std::int64_t>()[0], 3);
    EXPECT_DOUBLE_EQ(fd_value(out, kFitPatchArcName), 360.0);
}

TEST(FitPrimitivesForm, SteppedCylinderReportsBothDiameters) {
    std::vector<std::pair<double, double>> prof;
    for (int i = 0; i <= 4; ++i)
        prof.push_back({1.0, 0.25 * i});
    for (int i = 4; i <= 8; ++i)
        prof.push_back({1.002, 0.25 * i});
    const Mesh out = fp_fit_whole_mesh(analytic::revolve(prof), 3e-3);
    ASSERT_EQ(out.FieldData("fit:patch_type").As<std::int64_t>()[0], 3);
    EXPECT_NEAR(fd_value(out, kFitPatchDiameterMinName), 2.0, 1e-6);
    EXPECT_NEAR(fd_value(out, kFitPatchDiameterMaxName), 2.004, 1e-6);
    EXPECT_GT(fd_value(out, kFitPatchFormErrorName), 0.0015);
}

TEST(FitPrimitivesForm, PlaneExtentsAreTheInPlaneSides) {
    std::vector<std::vector<double>> pts;
    for (int j = 0; j <= 2; ++j)
        for (int i = 0; i <= 4; ++i)
            pts.push_back({0.5 * i, 0.5 * j, 0.0});
    std::vector<std::vector<std::int64_t>> quads;
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 4; ++i)
            quads.push_back({j * 5 + i, j * 5 + i + 1, (j + 1) * 5 + i + 1, (j + 1) * 5 + i});
    const Mesh out = fp_fit_whole_mesh(mt::make_mesh(pts, "quad", quads));
    ASSERT_EQ(out.FieldData("fit:patch_type").As<std::int64_t>()[0], 1);
    const double* e = out.FieldData(kFitPatchExtentName).As<double>();
    EXPECT_NEAR(e[0], 2.0, 1e-12);
    EXPECT_NEAR(e[1], 1.0, 1e-12);
    EXPECT_NEAR(fd_value(out, kFitPatchFormErrorName), 0.0, 1e-12);
    EXPECT_TRUE(std::isnan(fd_value(out, kFitPatchArcName)));
}

// ---- roadmap 7.1.1 and 7.1.2 --------------------------------------------------------------

namespace {

std::int64_t fp_region_size(const Mesh& rMesh, const std::string& rName) {
    const std::size_t i = rMesh.FindRegion(rName, RegionKind::Cell);
    return i == Mesh::npos ? -1 : static_cast<std::int64_t>(rMesh.Region(i).NumEntries());
}

// A cylinder of radius R along z, two bands z in [0, 1] and [1.5, 2.5]: a bore
// a cross hole interrupts.
Mesh fp_split_bore(double R) {
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<std::int64_t>> quads;
    const int n = 48;
    for (double z0 : {0.0, 1.5}) {
        const std::int64_t base = static_cast<std::int64_t>(pts.size());
        for (int i = 0; i <= 4; ++i)
            for (int k = 0; k < n; ++k) {
                const double a = 2.0 * analytic::kPi * k / n;
                pts.push_back({R * std::cos(a), R * std::sin(a), z0 + 0.25 * i});
            }
        for (int i = 0; i < 4; ++i)
            for (int k = 0; k < n; ++k) {
                const std::int64_t k1 = (k + 1) % n;
                quads.push_back({base + i * n + k, base + i * n + k1, base + (i + 1) * n + k1,
                                 base + (i + 1) * n + k});
            }
    }
    return mt::make_mesh(pts, "quad", quads);
}

}  // namespace

TEST(FitPrimitivesSeeds, ASeedOnTheHoleFitsTheHoleOnly) {
    const Mesh plate = analytic::filleted_plate();
    FitPrimitivesOptions o;
    o.mSeeds = {{0.6, 0.0, 0.25}};
    o.mTypes = 1u << 3;  // a cylinder
    const FitPrimitivesResult r = fit_primitives(plate, o);
    EXPECT_EQ(r.mNumCylinder, 1);
    // The hole's cells: those with every corner at radius 0.6.
    std::int64_t hole = 0;
    const NDArray& conn = plate.Cells(0).Conn();
    const std::size_t npc = plate.Cells(0).NodesPerCell();
    for (std::size_t c = 0; c < plate.Cells(0).NumCells(); ++c) {
        bool on = true;
        for (std::size_t k = 0; k < npc; ++k) {
            const double* p =
                plate.Points().As<double>() + 3 * conn.As<std::int64_t>()[c * npc + k];
            on = on && std::fabs(std::hypot(p[0], p[1]) - 0.6) < 1e-12;
        }
        hole += on ? 1 : 0;
    }
    EXPECT_EQ(r.mNumFeatureCells, hole);
    EXPECT_EQ(fp_region_size(r.mMesh, "fit:feature"), r.mNumFeatureCells);
    const std::vector<FpRow> rr = rows(r.mMesh);
    ASSERT_EQ(rr.size(), 1u);
    EXPECT_NEAR(rr[0].mParams[6], 0.6, 1e-9);
    EXPECT_NEAR(std::fabs(rr[0].mParams[5]), 1.0, 1e-9);  // axis along z
    // The rim angle peels nothing off a hole bounded by sharp edges.
    EXPECT_EQ(r.mMesh.FieldData(kFitPatchNumRejectedName).As<std::int64_t>()[0], 0);
}

TEST(FitPrimitivesSeeds, ASeedOnTheFilletFitsTheFilletOnly) {
    const Mesh plate = analytic::filleted_plate();
    const double c = 2.0 - 0.8, s = std::sqrt(0.5);
    FitPrimitivesOptions o;
    o.mSeeds = {{c + 0.8 * s, c + 0.8 * s, 0.25}};
    o.mTypes = 1u << 3;
    const FitPrimitivesResult r = fit_primitives(plate, o);
    ASSERT_EQ(r.mNumCylinder, 1) << (r.mWarnings.empty() ? "" : r.mWarnings[0]);
    // The fillet's cells: those with every corner on its cylinder.
    std::int64_t fillet = 0;
    const NDArray& conn = plate.Cells(0).Conn();
    const std::size_t npc = plate.Cells(0).NodesPerCell();
    for (std::size_t k = 0; k < plate.Cells(0).NumCells(); ++k) {
        bool on = true;
        for (std::size_t j = 0; j < npc; ++j) {
            const double* p =
                plate.Points().As<double>() + 3 * conn.As<std::int64_t>()[k * npc + j];
            on = on && std::fabs(std::hypot(p[0] - c, p[1] - c) - 0.8) < 1e-12 &&
                 p[0] > c - 1e-12 && p[1] > c - 1e-12;
        }
        fillet += on ? 1 : 0;
    }
    EXPECT_EQ(r.mNumFeatureCells, fillet);
    const std::vector<FpRow> rr = rows(r.mMesh);
    EXPECT_NEAR(rr[0].mParams[6], 0.8, 1e-9);
    // The axis passes through (c, c).
    EXPECT_NEAR(rr[0].mParams[0], c, 1e-9);
    EXPECT_NEAR(rr[0].mParams[1], c, 1e-9);
}

TEST(FitPrimitivesSeeds, TwoSeedsAcrossACrossHoleFitOneCylinder) {
    FitPrimitivesOptions o;
    o.mSeeds = {{1.0, 0.0, 0.5}, {0.0, 1.0, 2.0}};
    o.mTypes = 1u << 3;
    const FitPrimitivesResult r = fit_primitives(fp_split_bore(1.0), o);
    EXPECT_EQ(r.mNumPatches, 1);
    EXPECT_EQ(r.mNumCylinder, 1);
    EXPECT_EQ(r.mNumFeatureCells, 2 * 4 * 48);
    EXPECT_NEAR(rows(r.mMesh)[0].mParams[6], 1.0, 1e-9);
    // One seed alone fits its band.
    o.mSeeds.pop_back();
    EXPECT_EQ(fit_primitives(fp_split_bore(1.0), o).mNumFeatureCells, 4 * 48);
}

TEST(FitPrimitivesSeeds, RegionModeAndBadOptions) {
    Mesh plate = analytic::filleted_plate();
    FitPrimitivesOptions o;
    o.mSeeds = {{0.6, 0.0, 0.25}};
    o.mTypes = 1u << 3;
    const FitPrimitivesResult seeded = fit_primitives(plate, o);
    FitPrimitivesOptions ro;
    ro.mRegion = "fit:feature";
    ro.mTypes = 1u << 3;
    const FitPrimitivesResult region = fit_primitives(seeded.mMesh, ro);
    EXPECT_EQ(region.mNumCylinder, 1);
    EXPECT_NEAR(rows(region.mMesh)[0].mParams[6], 0.6, 1e-9);
    ro.mRegion = "nope";
    EXPECT_THROW(fit_primitives(plate, ro), std::invalid_argument);
    o.mRimAngleDeg = 200.0;
    EXPECT_THROW(fit_primitives(plate, o), std::invalid_argument);
}

TEST(FitPrimitivesTrimming, OutliersDoNotDragTheRadius) {
    // A noisy cylinder (sigma 1e-3) with every 20th point 10 sigma outward.
    const double sigma = 1e-3;
    Mesh cyl = analytic::cylinder(1.0, 2.0, 20);
    NDArray pts = NDArray::Uninit(DType::Float64, cyl.Points().Shape());
    SplitMix g{42};
    for (std::size_t i = 0; i < cyl.NumPoints(); ++i) {
        const double* p = cyl.Points().As<double>() + 3 * i;
        const double r = std::hypot(p[0], p[1]);
        // A sum of uniforms: near-normal noise of standard deviation sigma.
        double noise = 0.0;
        for (int k = 0; k < 12; ++k)
            noise += 0.5 * g.Next();
        noise *= sigma / std::sqrt(1.0);
        const double dr = (i % 20 == 0) ? 10.0 * sigma : noise;
        double* q = pts.As<double>() + 3 * i;
        q[0] = p[0] * (r + dr) / r;
        q[1] = p[1] * (r + dr) / r;
        q[2] = p[2];
    }
    cyl.AssignPoints(std::move(pts));
    FitPrimitivesOptions o;
    o.mPatchData = "";
    o.mTypes = 1u << 3;
    o.mTolerance = 0.02;
    o.mMinInlierFraction = 0.5;
    o.mUsedPoints = FitUsedPoints::All;
    const double plain = rows(fit_primitives(cyl, o).mMesh)[0].mParams[6];
    o.mUsedPoints = FitUsedPoints::Sigma3;
    const FitPrimitivesResult trimmed = fit_primitives(cyl, o);
    const double robust = rows(trimmed.mMesh)[0].mParams[6];
    EXPECT_LT(std::fabs(robust - 1.0), sigma);
    EXPECT_GT(std::fabs(plain - 1.0), std::fabs(robust - 1.0));
    const std::int64_t rejected =
        trimmed.mMesh.FieldData(kFitPatchNumRejectedName).As<std::int64_t>()[0];
    EXPECT_GE(rejected, static_cast<std::int64_t>(cyl.NumPoints() / 20));
    EXPECT_LT(trimmed.mMesh.FieldData(kFitPatchSigmaName).As<double>()[0], 2.0 * sigma);
}

TEST(FitPrimitivesTrimming, AConsensusStartFitsAPlaneThroughGrossOutliers) {
    // A plane with a third of its points lifted far off it: least squares tilts,
    // the least median of squares start does not.
    std::vector<std::vector<double>> p;
    std::vector<std::vector<std::int64_t>> quads;
    for (int j = 0; j <= 10; ++j)
        for (int i = 0; i <= 10; ++i)
            p.push_back({0.1 * i, 0.1 * j, (i + j) % 3 == 0 && i < 6 ? 0.2 : 0.0});
    for (int j = 0; j < 10; ++j)
        for (int i = 0; i < 10; ++i)
            quads.push_back({j * 11 + i, j * 11 + i + 1, (j + 1) * 11 + i + 1, (j + 1) * 11 + i});
    const Mesh plate = mt::make_mesh(p, "quad", quads);
    FitPrimitivesOptions o;
    o.mPatchData = "";
    o.mTypes = 1u << 1;
    o.mTolerance = 0.5;
    o.mMinInlierFraction = 0.0;
    o.mUsedPoints = FitUsedPoints::Sigma3;
    const FitPrimitivesResult r = fit_primitives(plate, o);
    const std::vector<FpRow> rr = rows(r.mMesh);
    ASSERT_EQ(rr[0].mType, 1);
    EXPECT_NEAR(std::fabs(rr[0].mParams[5]), 1.0, 1e-12);  // the normal is z
    EXPECT_NEAR(rr[0].mParams[2], 0.0, 1e-12);
}
