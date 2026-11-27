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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/critical_points.hpp"
#include "meshioplusplus/operations/data_integrate.hpp"
#include "meshioplusplus/operations/deform.hpp"
#include "meshioplusplus/operations/denoise_surface.hpp"
#include "meshioplusplus/operations/eigenmodes.hpp"
#include "meshioplusplus/operations/evaluate_cells.hpp"
#include "meshioplusplus/operations/fvm_quality.hpp"
#include "meshioplusplus/operations/gradient.hpp"
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/operations/mesh_operator.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/quality_gate.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/shortest_path.hpp"
#include "meshioplusplus/operations/skeleton.hpp"
#include "meshioplusplus/operations/smooth_field.hpp"
#include "meshioplusplus/operations/spectral_descriptors.hpp"
#include "meshioplusplus/operations/transfer_operator.hpp"
#include "meshioplusplus/operations/transform.hpp"

#include "../../src/cpp/src/detail/reference_element.hpp"
#include "../../src/cpp/src/detail/sparse_solve.hpp"
#include "../../src/cpp/src/detail/surface_fem.hpp"

using namespace meshioplusplus;

namespace {

std::vector<double> csr_matvec(const MeshOperatorResult& rOp, const std::vector<double>& rX) {
    std::vector<double> y(static_cast<std::size_t>(rOp.mRows), 0.0);
    for (std::int64_t r = 0; r < rOp.mRows; ++r) {
        double acc = 0.0;
        for (std::int64_t k = rOp.mIndptr[static_cast<std::size_t>(r)];
             k < rOp.mIndptr[static_cast<std::size_t>(r) + 1]; ++k) {
            acc += rOp.mValues[static_cast<std::size_t>(k)] *
                   rX[static_cast<std::size_t>(rOp.mIndices[static_cast<std::size_t>(k)])];
        }
        y[static_cast<std::size_t>(r)] = acc;
    }
    return y;
}

}  // namespace

// --- 15.1 Integration by quadrature -----------------------------------------

TEST(Roadmap15, QuadratureIntegrateQuadraticCells) {
    // Unit square as a single Quad8 with quadratic scalar f(x, y) = x^2 + 3 y^2.
    // Exact integral over [0, 1]^2 = 1/3 + 1 = 4/3.
    Mesh m;
    NDArray pts(DType::Float64, {8, 3});
    const double coords[8][3] = {{0, 0, 0},   {1, 0, 0},   {1, 1, 0},   {0, 1, 0},
                                 {0.5, 0, 0}, {1, 0.5, 0}, {0.5, 1, 0}, {0, 0.5, 0}};
    double* pp = pts.As<double>();
    NDArray f(DType::Float64, {8});
    double* pf = f.As<double>();
    for (int i = 0; i < 8; ++i) {
        pp[i * 3 + 0] = coords[i][0];
        pp[i * 3 + 1] = coords[i][1];
        pp[i * 3 + 2] = coords[i][2];
        pf[i] = coords[i][0] * coords[i][0] + 3.0 * coords[i][1] * coords[i][1];
    }
    m.AssignPoints(std::move(pts));
    NDArray conn(DType::Int64, {1, 8});
    for (int i = 0; i < 8; ++i)
        conn.As<std::int64_t>()[i] = i;
    m.AddCellBlock("quad8", std::move(conn));
    m.AddPointData("f", std::move(f));

    DataIntegrateOptions opt;
    opt.mArrayNames = {"quad:f", "point:f"};
    const DataIntegrateReport rep = data_integrate(m, opt);
    ASSERT_EQ(rep.mArrays.size(), 2u);
    EXPECT_NEAR(rep.mArrays[0].mDomain.mTotalPerComponent[0], 4.0 / 3.0, 1e-12);
    // First-order point: mean is NOT exact for quadratic f.
    EXPECT_GT(std::abs(rep.mArrays[1].mDomain.mTotalPerComponent[0] - 4.0 / 3.0), 1e-2);
}

// --- 15.2.1 Tetrahedral operators & smooth_field ----------------------------

TEST(Roadmap15, TetraMeshOperatorAndSmoothField) {
    Mesh tet = box({2, 2, 2}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "tetra");
    const std::size_t n = tet.NumPoints();

    const MeshOperatorResult K = mesh_operator(tet, {MeshOperatorKind::Laplacian});
    const MeshOperatorResult ML = mesh_operator(tet, {MeshOperatorKind::MassLumped});
    const MeshOperatorResult MC = mesh_operator(tet, {MeshOperatorKind::MassConsistent});
    const MeshOperatorResult G = mesh_operator(tet, {MeshOperatorKind::Gradient});

    EXPECT_EQ(K.mRows, static_cast<std::int64_t>(n));
    EXPECT_EQ(K.mCols, static_cast<std::int64_t>(n));

    // Constant vector is in nullspace of K; 1^T M 1 = Volume = 1.0.
    std::vector<double> ones(n, 1.0);
    const auto k1 = csr_matvec(K, ones);
    for (double v : k1)
        EXPECT_NEAR(v, 0.0, 1e-12);

    const auto ml1 = csr_matvec(ML, ones);
    const auto mc1 = csr_matvec(MC, ones);
    double vl = 0.0, vc = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        vl += ml1[i];
        vc += mc1[i];
    }
    EXPECT_NEAR(vl, 1.0, 1e-12);
    EXPECT_NEAR(vc, 1.0, 1e-12);

    // Linear field f = 2x - 3y + 5z has exact gradient (2, -3, 5) on every tet.
    std::vector<double> lin(n);
    const NDArray& pts = tet.Points();
    for (std::size_t i = 0; i < n; ++i) {
        lin[i] = 2.0 * detail::read_double(pts, i * 3 + 0) -
                 3.0 * detail::read_double(pts, i * 3 + 1) +
                 5.0 * detail::read_double(pts, i * 3 + 2);
    }
    const auto gf = csr_matvec(G, lin);
    ASSERT_EQ(gf.size(), static_cast<std::size_t>(G.mRows));
    for (std::size_t c = 0; c < gf.size() / 3; ++c) {
        EXPECT_NEAR(gf[c * 3 + 0], 2.0, 1e-11);
        EXPECT_NEAR(gf[c * 3 + 1], -3.0, 1e-11);
        EXPECT_NEAR(gf[c * 3 + 2], 5.0, 1e-11);
    }

    // smooth_field preserves M-weighted mean on tetrahedra.
    NDArray f_arr(DType::Float64, {n});
    for (std::size_t i = 0; i < n; ++i)
        f_arr.As<double>()[i] = lin[i] * lin[i];
    tet.AddPointData("u", std::move(f_arr));
    SmoothFieldOptions sopt;
    sopt.mArray = "u";
    sopt.mLength = 0.25;
    const SmoothFieldResult sres = smooth_field(tet, sopt);
    EXPECT_TRUE(sres.mConverged);
    const double* pu = tet.PointData("u").As<double>();
    const double* ps = sres.mMesh.PointData("u:smooth").As<double>();
    double mean_in = 0.0, mean_out = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        mean_in += ml1[i] * pu[i];
        mean_out += ml1[i] * ps[i];
    }
    EXPECT_NEAR(mean_in, mean_out, 1e-10);
}

// --- 15.2.3 DEC spike verification (d1 * d0 == 0, d0^T *1 d0 == K_cot) ------

TEST(Roadmap15, DecSpikeExactness) {
    Mesh surf = icosa_sphere(1.0, 2);
    const detail::SurfaceFem fem = detail::build_surface_fem(surf, "dec_spike");
    const std::size_t nv = fem.mPoints.size();
    const std::size_t nf = fem.mTriangles.size();

    // Index oriented undirected edges (lo < hi).
    std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> edge_id;
    std::vector<std::pair<std::int64_t, std::int64_t>> edges;
    auto get_edge = [&](std::int64_t a, std::int64_t b, double& rSign) -> std::size_t {
        if (a < b) {
            rSign = 1.0;
            auto it = edge_id.find({a, b});
            if (it == edge_id.end()) {
                const std::size_t id = edges.size();
                edge_id[{a, b}] = id;
                edges.emplace_back(a, b);
                return id;
            }
            return it->second;
        }
        rSign = -1.0;
        auto it = edge_id.find({b, a});
        if (it == edge_id.end()) {
            const std::size_t id = edges.size();
            edge_id[{b, a}] = id;
            edges.emplace_back(b, a);
            return id;
        }
        return it->second;
    };

    std::vector<detail::SparseTriplet> d1_trip;
    std::vector<double> star1;
    for (std::size_t f = 0; f < nf; ++f) {
        const auto& t = fem.mTriangles[f];
        for (int k = 0; k < 3; ++k) {
            double sgn = 1.0;
            const std::size_t e = get_edge(t.mNodes[k], t.mNodes[(k + 1) % 3], sgn);
            d1_trip.push_back({f, e, sgn});
            if (star1.size() <= e)
                star1.resize(e + 1, 0.0);
            // Opposite vertex is (k + 2) % 3; cotangent weight is 0.5 * cot.
            star1[e] += 0.5 * t.mCot[(k + 2) % 3];
        }
    }
    const std::size_t ne = edges.size();
    std::vector<detail::SparseTriplet> d0_trip;
    d0_trip.reserve(2 * ne);
    for (std::size_t e = 0; e < ne; ++e) {
        d0_trip.push_back({e, static_cast<std::size_t>(edges[e].first), -1.0});
        d0_trip.push_back({e, static_cast<std::size_t>(edges[e].second), 1.0});
    }

    const detail::SparseCsr d0 = detail::csr_from_triplets(ne, nv, std::move(d0_trip));
    const detail::SparseCsr d1 = detail::csr_from_triplets(nf, ne, std::move(d1_trip));

    // Verify d1 * d0 == 0 identically on arbitrary 0-form.
    std::vector<double> f0(nv);
    for (std::size_t i = 0; i < nv; ++i)
        f0[i] = std::sin(static_cast<double>(i) * 0.37) + 1.0;
    std::vector<double> df(ne, 0.0), ddf(nf, 0.0);
    detail::csr_multiply(d0, f0.data(), df.data());
    detail::csr_multiply(d1, df.data(), ddf.data());
    for (double v : ddf)
        EXPECT_NEAR(v, 0.0, 1e-14);

    // Verify d0^T * star1 * d0 matches mesh_operator(Laplacian).
    const MeshOperatorResult K = mesh_operator(surf, {MeshOperatorKind::Laplacian});
    const auto kf = csr_matvec(K, f0);
    std::vector<double> dec_kf(nv, 0.0);
    for (std::size_t e = 0; e < ne; ++e) {
        const double wdf = star1[e] * df[e];
        dec_kf[static_cast<std::size_t>(edges[e].first)] -= wdf;
        dec_kf[static_cast<std::size_t>(edges[e].second)] += wdf;
    }
    for (std::size_t i = 0; i < nv; ++i)
        EXPECT_NEAR(dec_kf[i], kf[i], 1e-11);
}

// --- 15.5 Inter-mesh transfer operator --------------------------------------

TEST(Roadmap15, TransferOperatorNestedAndBarycentric) {
    Mesh coarse = box({2, 2, 2}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "tetra");
    RefineOptions ropt;
    ropt.mRecordHierarchy = true;
    Mesh fine = refine(coarse, ropt).mMesh;

    const MeshOperatorResult P = transfer_operator(coarse, fine);
    EXPECT_EQ(P.mRows, static_cast<std::int64_t>(fine.NumPoints()));
    EXPECT_EQ(P.mCols, static_cast<std::int64_t>(coarse.NumPoints()));
    EXPECT_EQ(P.mNumDegenerate, 0);

    // Partition of unity: every row sums to 1.
    std::vector<double> ones(coarse.NumPoints(), 1.0);
    const auto p1 = csr_matvec(P, ones);
    for (double v : p1)
        EXPECT_NEAR(v, 1.0, 1e-13);

    // Exact linear field reproduction on fine.
    std::vector<double> f_coarse(coarse.NumPoints());
    const NDArray& cp = coarse.Points();
    for (std::size_t i = 0; i < coarse.NumPoints(); ++i)
        f_coarse[i] = 1.0 + 2.0 * detail::read_double(cp, i * 3) -
                      3.0 * detail::read_double(cp, i * 3 + 1) +
                      4.0 * detail::read_double(cp, i * 3 + 2);
    const auto f_fine = csr_matvec(P, f_coarse);
    const NDArray& fp = fine.Points();
    for (std::size_t i = 0; i < fine.NumPoints(); ++i) {
        const double expected = 1.0 + 2.0 * detail::read_double(fp, i * 3) -
                                3.0 * detail::read_double(fp, i * 3 + 1) +
                                4.0 * detail::read_double(fp, i * 3 + 2);
        EXPECT_NEAR(f_fine[i], expected, 1e-12);
    }

    // Transpose restriction R = P^T has exact transposed dimensions.
    TransferOperatorOptions topt;
    topt.mMode = TransferOperatorMode::Barycentric;
    topt.mTranspose = true;
    const MeshOperatorResult R = transfer_operator(coarse, fine, topt);
    EXPECT_EQ(R.mRows, static_cast<std::int64_t>(coarse.NumPoints()));
    EXPECT_EQ(R.mCols, static_cast<std::int64_t>(fine.NumPoints()));
}

// --- 15.6 Spectral shape descriptors (HKS & WKS) ----------------------------

TEST(Roadmap15, SpectralDescriptorsRigidInvariance) {
    Mesh s1 = icosa_sphere(1.0, 2);
    Mesh s2 = icosa_sphere(1.0, 2);
    {
        NDArray pts = s2.Points();
        pts.MakeOwned();
        double* pp = pts.As<double>();
        const double ang = 47.0 * 3.14159265358979323846 / 180.0;
        const double c = std::cos(ang), s = std::sin(ang);
        for (std::size_t i = 0; i < s2.NumPoints(); ++i) {
            const double x = pp[i * 3 + 0], y = pp[i * 3 + 1], z = pp[i * 3 + 2];
            pp[i * 3 + 0] = c * x - s * y + 3.0;
            pp[i * 3 + 1] = s * x + c * y - 2.0;
            pp[i * 3 + 2] = z + 5.0;
        }
        s2.AssignPoints(std::move(pts));
    }

    SpectralDescriptorsOptions opt;
    opt.mKind = SpectralDescriptorKind::Hks;
    opt.mNumModes = 12;
    opt.mNumScales = 8;
    const SpectralDescriptorsResult r1 = spectral_descriptors(s1, opt);
    const SpectralDescriptorsResult r2 = spectral_descriptors(s2, opt);
    EXPECT_FALSE(r1.mReusedEigenmodes);
    EXPECT_EQ(r1.mScales.size(), 8u);

    const double* a1 = r1.mMesh.PointData(kSpectralHksName).As<double>();
    const double* a2 = r2.mMesh.PointData(kSpectralHksName).As<double>();
    for (std::size_t i = 0; i < s1.NumPoints() * 8; ++i)
        EXPECT_NEAR(a1[i], a2[i], 1e-5);

    // Reuse eigenmodes for WKS without re-running eigensolver.
    opt.mKind = SpectralDescriptorKind::Wks;
    const SpectralDescriptorsResult w1 = spectral_descriptors(r1.mMesh, opt);
    EXPECT_TRUE(w1.mReusedEigenmodes);
    EXPECT_TRUE(w1.mMesh.HasPointData(kSpectralWksName));
}

// --- 9.11 FVM face-pair mesh quality & quality_gate -------------------------

TEST(Roadmap9, FvmQualityAndGate) {
    // Orthogonal hex grid: non-orthogonality = 0 deg, skewness = 0, volume ratio = 1.
    Mesh hex = box({2, 2, 2}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "hexahedron");
    const FvmQualityReport rep = compute_fvm_quality(hex);
    EXPECT_EQ(rep.mNumCells, 8);
    EXPECT_EQ(rep.mNumInteriorFaces, 12);
    EXPECT_EQ(rep.mNumNegativePyramids, 0);
    EXPECT_NEAR(rep.mMetrics[0].second.mMax, 0.0, 1e-12);  // non_orthogonality
    EXPECT_NEAR(rep.mMetrics[1].second.mMax, 0.0, 1e-12);  // face_skewness
    EXPECT_NEAR(rep.mMetrics[2].second.mMax, 1.0, 1e-12);  // volume_ratio

    QualityGateOptions gopt;
    gopt.mThresholds =
        parse_quality_thresholds("fvm:non_orthogonality <= 70; fvm:volume_ratio <= 2");
    const QualityGateResult gres = check_quality(hex, gopt);
    EXPECT_TRUE(gres.mPassed);

    // Shear the hex mesh by x' = x + 1.5*z to create > 45 deg non-orthogonality.
    Mesh skewed = box({2, 2, 2}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "hexahedron");
    NDArray pts = skewed.Points();
    pts.MakeOwned();
    double* pp = pts.As<double>();
    for (std::size_t i = 0; i < skewed.NumPoints(); ++i)
        pp[i * 3 + 0] += 1.5 * pp[i * 3 + 2];
    skewed.AssignPoints(std::move(pts));
    const FvmQualityReport srep = compute_fvm_quality(skewed);
    EXPECT_GT(srep.mMetrics[0].second.mMax, 35.0);
}

// --- 9.12 Feature-preserving surface denoising ------------------------------

TEST(Roadmap9, DenoiseSurfaceReducesNoise) {
    Mesh clean = icosa_sphere(1.0, 2);
    Mesh noisy = icosa_sphere(1.0, 2);
    NDArray pts = noisy.Points();
    pts.MakeOwned();
    double* pp = pts.As<double>();
    std::mt19937_64 rng(42);
    std::normal_distribution<double> dist(0.0, 0.03);
    for (std::size_t i = 0; i < noisy.NumPoints(); ++i) {
        const double x = pp[i * 3 + 0], y = pp[i * 3 + 1], z = pp[i * 3 + 2];
        const double r = std::sqrt(x * x + y * y + z * z);
        const double dr = dist(rng);
        pp[i * 3 + 0] += dr * x / r;
        pp[i * 3 + 1] += dr * y / r;
        pp[i * 3 + 2] += dr * z / r;
    }
    noisy.AssignPoints(std::move(pts));

    auto rms_err = [](const Mesh& m) {
        const NDArray& p = m.Points();
        double s = 0.0;
        for (std::size_t i = 0; i < m.NumPoints(); ++i) {
            const double x = detail::read_double(p, i * 3 + 0);
            const double y = detail::read_double(p, i * 3 + 1);
            const double z = detail::read_double(p, i * 3 + 2);
            const double e = std::sqrt(x * x + y * y + z * z) - 1.0;
            s += e * e;
        }
        return std::sqrt(s / static_cast<double>(m.NumPoints()));
    };

    const double err_before = rms_err(noisy);
    DenoiseSurfaceOptions opt;
    opt.mMethod = DenoiseMethod::Bilateral;
    opt.mNormalIterations = 10;
    opt.mVertexIterations = 10;
    const DenoiseSurfaceResult dres = denoise_surface(noisy, opt);
    const double err_after = rms_err(dres.mMesh);
    EXPECT_LT(err_after, 0.5 * err_before);
}

// --- 9.13 Curve-skeleton extraction -----------------------------------------

TEST(Roadmap9, ExtractSkeletonTorus) {
    Mesh tor = torus(1.0, 0.25, 20, 12, {0.0, 0.0, 0.0}, "triangle");
    SkeletonOptions opt;
    opt.mMaxIterations = 10;
    const SkeletonResult sres = extract_skeleton(tor, opt);
    EXPECT_GT(sres.mNumNodes, 4);
    EXPECT_GT(sres.mNumSegments, 4);
    EXPECT_TRUE(sres.mMesh.HasPointData(kSkeletonRadiusName));
    EXPECT_TRUE(sres.mMesh.HasPointData(kSkeletonDegreeName));
    EXPECT_TRUE(sres.mMesh.HasFieldData(kSkeletonSurfaceMapName));
    const double* rad = sres.mMesh.PointData(kSkeletonRadiusName).As<double>();
    double mean_r = 0.0;
    for (std::int64_t i = 0; i < sres.mNumNodes; ++i)
        mean_r += rad[i];
    mean_r /= static_cast<double>(sres.mNumNodes);
    EXPECT_NEAR(mean_r, 0.25, 0.05);
}

// --- 9.14 Shortest geodesic paths & surface Voronoi -------------------------

TEST(Roadmap9, ShortestPathAndGeodesicVoronoi) {
    Mesh sph = icosa_sphere(1.0, 2);
    // Find vertices closest to north pole (0,0,1) and south pole (0,0,-1).
    ShortestPathOptions opt;
    opt.mHasSourcePoint = true;
    opt.mSourcePoint = {0.0, 0.0, 1.0};
    opt.mHasTargetPoint = true;
    opt.mTargetPoint = {0.0, 0.0, -1.0};
    opt.mMethod = ShortestPathMethod::Dijkstra;
    const ShortestPathResult dijk = shortest_path(sph, opt);
    EXPECT_TRUE(dijk.mReached);
    EXPECT_NEAR(dijk.mTotalLength, 3.14159265358979323846, 0.25);

    opt.mMethod = ShortestPathMethod::HeatTrace;
    const ShortestPathResult ht = shortest_path(sph, opt);
    EXPECT_TRUE(ht.mReached);
    EXPECT_NEAR(ht.mTotalLength, 3.14159265358979323846, 0.25);

    GeodesicVoronoiOptions vopt;
    vopt.mSeedPoints = {0.0, 0.0, 1.0, 0.0, 0.0, -1.0};
    const GeodesicVoronoiResult vres = geodesic_voronoi(sph, vopt);
    EXPECT_EQ(vres.mNumSources, 2);
    EXPECT_EQ(vres.mNumUnreached, 0);
    EXPECT_TRUE(vres.mMesh.HasPointData(kGeodesicNearestSourceName));
}

// --- 9.15 FFD & RBF morphing ------------------------------------------------

TEST(Roadmap9, DeformFfdAndRbf) {
    Mesh cube = box({2, 2, 2}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "hexahedron");
    FfdOptions fopt;
    fopt.mGridDims = {2, 2, 2};
    fopt.mDisplacements.assign(8 * 3, 0.0);
    // Stretch top (k=1) vertices by +0.2 in z.
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i)
            fopt.mDisplacements[((1 * 2 + j) * 2 + i) * 3 + 2] = 0.2;
    const DeformResult fres = deform_ffd(cube, fopt);
    EXPECT_EQ(fres.mNumInverted, 0);
    EXPECT_NEAR(fres.mMaxDisplacement, 0.2, 1e-12);

    RbfOptions ropt;
    ropt.mKernel = RbfKernel::Tps;
    ropt.mControlPoints = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1};
    ropt.mDisplacements = {0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0.2, 0.1, 0, 0.2};
    const DeformResult rres = deform_rbf(cube, ropt);
    EXPECT_EQ(rres.mNumInverted, 0);
    EXPECT_GT(rres.mMaxDisplacement, 0.1);
}

// --- 15.3 High-order evaluation, change_basis, critical_points, gradient ----

TEST(Roadmap15, HighOrderEvaluationChangeBasisAndCriticalPoints) {
    // Build a single quadratic Hexahedron27 on [0, 1]^3 with saddle field
    // f(x, y, z) = (x - 0.5)^2 + (y - 0.5)^2 - (z - 0.5)^2.
    const auto lat = lagrange::vtk_lattice(lagrange::Shape::Hexahedron, 2);
    Mesh m;
    NDArray pts(DType::Float64, {27, 3});
    NDArray f(DType::Float64, {27});
    double* pp = pts.As<double>();
    double* pf = f.As<double>();
    for (std::size_t i = 0; i < 27; ++i) {
        const double x = 0.5 * lat[i][0];
        const double y = 0.5 * lat[i][1];
        const double z = 0.5 * lat[i][2];
        pp[i * 3 + 0] = x;
        pp[i * 3 + 1] = y;
        pp[i * 3 + 2] = z;
        pf[i] = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) - (z - 0.5) * (z - 0.5);
    }
    m.AssignPoints(std::move(pts));
    NDArray conn(DType::Int64, {1, 27});
    for (std::size_t i = 0; i < 27; ++i)
        conn.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    m.AddCellBlock("VTK_LAGRANGE_HEXAHEDRON", std::move(conn));
    m.AddPointData("f", std::move(f));

    // 15.3.1 evaluate_cells
    EvaluateCellsOptions eopt;
    eopt.mCells = {0};
    eopt.mReferencePoints = {0.25, 0.75, 0.5};
    eopt.mFields = {"f"};
    eopt.mDerivatives = true;
    const EvaluateCellsResult eres = evaluate_cells(m, eopt);
    EXPECT_EQ(eres.mNumEvaluated, 1);
    EXPECT_NEAR(eres.mValues[0].second.As<double>()[0], 0.0625 + 0.0625, 1e-12);
    const double* g = eres.mGradients[0].second.As<double>();
    EXPECT_NEAR(g[0], -0.5, 1e-12);
    EXPECT_NEAR(g[1], 0.5, 1e-12);
    EXPECT_NEAR(g[2], 0.0, 1e-12);

    // 15.3.1 gradient on VtkLagrangeHexahedron (centroid (0.5, 0.5, 0.5) => grad = 0, not NaN!)
    GradientOptions gropt;
    gropt.mArrayName = "f";
    gropt.mLocation = DataLocation::Cell;
    const GradientResult grres = gradient(m, gropt);
    ASSERT_FALSE(grres.mMesh.CellDataNames().empty());
    const double* pg = grres.mMesh.CellData(grres.mMesh.CellDataNames()[0], 0).As<double>();
    EXPECT_NEAR(pg[0], 0.0, 1e-12);
    EXPECT_NEAR(pg[1], 0.0, 1e-12);
    EXPECT_NEAR(pg[2], 0.0, 1e-12);

    // 15.3.2 change_basis: Lagrange -> Bezier -> Lagrange is identity
    const ChangeBasisResult bez = change_basis(m, {BasisTarget::Bezier, true});
    EXPECT_EQ(bez.mNumConverted, 1);
    EXPECT_TRUE(bez.mMesh.HasCellData("HigherOrderDegrees"));
    const ChangeBasisResult back = change_basis(bez.mMesh, {BasisTarget::Lagrange, true});
    const double* pf_in = m.PointData("f").As<double>();
    const double* pf_back = back.mMesh.PointData("f").As<double>();
    for (std::size_t i = 0; i < 27; ++i)
        EXPECT_NEAR(pf_back[i], pf_in[i], 1e-11);

    // 15.3.3 critical_points: finds the saddle at (0.5, 0.5, 0.5)
    CriticalPointsOptions copt;
    copt.mArray = "f";
    const CriticalPointsResult cpres = critical_points(m, copt);
    EXPECT_EQ(cpres.mNumSaddles, 1);
    ASSERT_EQ(cpres.mMesh.NumPoints(), 1u);
    const double* cp_pts = cpres.mMesh.Points().As<double>();
    EXPECT_NEAR(cp_pts[0], 0.5, 1e-6);
    EXPECT_NEAR(cp_pts[1], 0.5, 1e-6);
    EXPECT_NEAR(cp_pts[2], 0.5, 1e-6);

    // 15.3.4 high-order isosurface on quadratic field f = 0.05 produces non-empty contour
    IsosurfaceOptions iopt;
    iopt.mArrayName = "f";
    iopt.mIsovalues = {0.05};
    const Mesh iso = isosurface(m, iopt);
    EXPECT_GT(iso.NumPoints(), 10u);
}
