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
/**
 * @file test_iga.cpp
 * @brief The isogeometric layer: knot insertion and degree elevation that
 *        leave the geometry, the fields and a trimmed face's trimming loops
 *        exactly where they were, evaluation with derivatives and normals,
 *        the generic NURBS sampler and its MFEM twin, and the error paths,
 *        on the native path directly (the Python tests cover parity with the
 *        numpy twin on every NURBS fixture).
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/formats/mfem.hpp"
#include "meshioplusplus/operations/iga.hpp"

namespace {

using meshioplusplus::BrepInfo;
using meshioplusplus::IgaRefineOptions;
using meshioplusplus::IgaSampleOptions;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::NurbsInfo;
using meshioplusplus::NurbsPatch;

NDArray doubles(const std::vector<std::size_t>& rShape, const std::vector<double>& rValues) {
    NDArray a(meshioplusplus::DType::Float64, rShape);
    std::copy(rValues.begin(), rValues.end(), a.As<double>());
    return a;
}

// A patch from its parts: control points are rows of `sdim`, weights optional.
NurbsPatch make_patch(std::vector<std::int32_t> Degrees, std::vector<std::vector<double>> Knots,
                      std::size_t Sdim, const std::vector<double>& rPoints,
                      const std::vector<double>& rWeights = {}) {
    NurbsPatch p;
    p.mDegrees = std::move(Degrees);
    p.mKnots = std::move(Knots);
    const std::size_t n = rPoints.size() / Sdim;
    p.mControlPoints = doubles({n, Sdim}, rPoints);
    p.mWeights = doubles({n}, rWeights.empty() ? std::vector<double>(n, 1.0) : rWeights);
    return p;
}

// The rational quadratic quarter circle of radius 1.
NurbsPatch arc() {
    const double s = std::sqrt(0.5);
    return make_patch({2}, {{0, 0, 0, 1, 1, 1}}, 2, {1, 0, 1, 1, 0, 1}, {1.0, s, 1.0});
}

// The quarter annulus, radii 1 and 2, order 2 x 1 and rational: a surface in 3-D.
NurbsPatch annulus() {
    const double s = std::sqrt(0.5);
    // rows along v (radius 1 then 2), three control points along u each
    std::vector<double> pts, w;
    for (const double r : {1.0, 2.0}) {
        pts.insert(pts.end(), {r, 0, 0, r, r, 0, 0, r, 0});
        w.insert(w.end(), {1.0, s, 1.0});
    }
    return make_patch({2, 1}, {{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}, 3, pts, w);
}

// A trilinear hexahedron with its corners pushed about, polynomial.
NurbsPatch hexahedron() {
    return make_patch(
        {1, 1, 1}, {{0, 0, 1, 1}, {0, 0, 1, 1}, {0, 0, 1, 1}}, 3,
        {0, 0, 0, 1, 0, 0, 0, 1, 0, 1.2, 1.1, 0, 0, 0, 1, 1, 0, 1.3, 0, 1, 1, 1.1, 1.2, 1.4});
}

NDArray grid(const NurbsPatch& rP, std::size_t PerDirection) {
    const std::size_t dim = rP.mDegrees.size();
    std::size_t total = 1;
    for (std::size_t d = 0; d < dim; ++d)
        total *= PerDirection;
    NDArray u(meshioplusplus::DType::Float64, {total, dim});
    double* pu = u.As<double>();
    for (std::size_t i = 0; i < total; ++i) {
        std::size_t rest = i;
        for (std::size_t d = 0; d < dim; ++d) {
            const std::size_t k = rest % PerDirection;
            rest /= PerDirection;
            const double lo = rP.mKnots[d][static_cast<std::size_t>(rP.mDegrees[d])];
            const double hi =
                rP.mKnots[d][rP.mKnots[d].size() - static_cast<std::size_t>(rP.mDegrees[d]) - 1];
            pu[i * dim + d] =
                lo + (hi - lo) * static_cast<double>(k) / static_cast<double>(PerDirection - 1);
        }
    }
    return u;
}

NurbsInfo one_patch(NurbsPatch P) {
    NurbsInfo info;
    info.mPatches.push_back(std::move(P));
    return info;
}

double max_diff(const NDArray& rA, const NDArray& rB) {
    EXPECT_EQ(rA.Size(), rB.Size());
    double d = 0.0;
    for (std::size_t i = 0; i < std::min(rA.Size(), rB.Size()); ++i)
        d = std::max(d, std::abs(rA.As<double>()[i] - rB.As<double>()[i]));
    return d;
}

TEST(Iga, KnotInsertionKeepsTheArcExact) {
    const NurbsPatch p = arc();
    const NurbsPatch q = meshioplusplus::iga_insert_knots(p, 0, {0.25, 0.5, 0.5, 0.75});
    EXPECT_EQ(q.mKnots[0], (std::vector<double>{0, 0, 0, 0.25, 0.5, 0.5, 0.75, 1, 1, 1}));
    EXPECT_EQ(q.mControlPoints.Shape()[0], 7u);
    EXPECT_EQ(q.mWeights.Size(), 7u);
    const NDArray u = grid(p, 41);
    const auto a = meshioplusplus::iga_evaluate(one_patch(p), 0, u);
    const auto b = meshioplusplus::iga_evaluate(one_patch(q), 0, u);
    EXPECT_LT(max_diff(a.mPoints, b.mPoints), 1e-14);
    for (std::size_t i = 0; i < b.mPoints.Shape()[0]; ++i)
        EXPECT_NEAR(std::hypot(b.mPoints.As<double>()[i * 2], b.mPoints.As<double>()[i * 2 + 1]),
                    1.0, 1e-14);
}

TEST(Iga, DegreeElevationKeepsTheArcExactAndRaisesContinuity) {
    const NurbsPatch p = arc();
    const NurbsPatch q = meshioplusplus::iga_elevate_degree(p, 0, 2);
    EXPECT_EQ(q.mDegrees[0], 4);
    EXPECT_EQ(q.mKnots[0], (std::vector<double>{0, 0, 0, 0, 0, 1, 1, 1, 1, 1}));
    EXPECT_EQ(q.mControlPoints.Shape()[0], 5u);  // a quartic Bezier segment
    const NDArray u = grid(p, 31);
    EXPECT_LT(max_diff(meshioplusplus::iga_evaluate(one_patch(p), 0, u).mPoints,
                       meshioplusplus::iga_evaluate(one_patch(q), 0, u).mPoints),
              1e-14);
}

TEST(Iga, ElevationOfAMultiSpanCurveKeepsEveryKnotsContinuity) {
    // C1 at 1: quadratic, one interior knot of multiplicity 1
    const NurbsPatch p = make_patch({2}, {{0, 0, 0, 1, 2, 2, 2}}, 2, {0, 0, 1, 2, 2, 0, 3, 2});
    const NurbsPatch q = meshioplusplus::iga_elevate_degree(p, 0, 2);
    EXPECT_EQ(q.mKnots[0], (std::vector<double>{0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2}));
    const NDArray u = grid(p, 51);
    EXPECT_LT(max_diff(meshioplusplus::iga_evaluate(one_patch(p), 0, u).mPoints,
                       meshioplusplus::iga_evaluate(one_patch(q), 0, u).mPoints),
              1e-14);
}

TEST(Iga, RefinesASurfaceAndItsDerivativesStay) {
    IgaRefineOptions o;
    o.mSubdivisions = {3, 2, 1};
    o.mElevate = {1, 2, 0};
    const NurbsInfo in = one_patch(annulus());
    const NurbsInfo out = meshioplusplus::iga_refine(in, o);
    const NurbsPatch& q = out.mPatches[0];
    EXPECT_EQ(q.mDegrees, (std::vector<std::int32_t>{3, 3}));
    const NDArray u = grid(in.mPatches[0], 9);
    const auto a = meshioplusplus::iga_evaluate(in, 0, u, 1);
    const auto b = meshioplusplus::iga_evaluate(out, 0, u, 1);
    EXPECT_LT(max_diff(a.mPoints, b.mPoints), 1e-13);
    EXPECT_LT(max_diff(a.mDerivatives, b.mDerivatives), 1e-12);
    EXPECT_LT(max_diff(a.mNormals, b.mNormals), 1e-12);
    // the annulus is a 2-D surface in the plane: unit normals along +-z, radius 1..2
    for (std::size_t i = 0; i < a.mNormals.Shape()[0]; ++i) {
        const double* n = a.mNormals.As<double>() + i * 3;
        EXPECT_NEAR(std::abs(n[2]), 1.0, 1e-13);
        EXPECT_NEAR(n[0], 0.0, 1e-13);
    }
    // cleared, because the sampling they described is gone
    EXPECT_EQ(out.mCellPatch.Size(), 0u);
    EXPECT_EQ(out.mBoundaryCellPatch.Size(), 0u);
}

TEST(Iga, DerivativesMatchCentralDifferencesInAVolume) {
    const NurbsInfo in = one_patch(hexahedron());
    IgaRefineOptions o;
    o.mSubdivisions = {2, 2, 2};
    o.mElevate = {1, 1, 1};
    const NurbsInfo out = meshioplusplus::iga_refine(in, o);
    const NDArray u = grid(in.mPatches[0], 5);
    const auto a = meshioplusplus::iga_evaluate(in, 0, u);
    const auto b = meshioplusplus::iga_evaluate(out, 0, u, 1);
    EXPECT_LT(max_diff(a.mPoints, b.mPoints), 1e-13);
    EXPECT_EQ(b.mNormals.Size(), 0u);  // a volume has no normal
    const std::size_t n = u.Shape()[0];
    const double h = 1e-6;
    for (std::size_t d = 0; d < 3; ++d) {
        NDArray fwd = u, bwd = u;
        for (std::size_t i = 0; i < n; ++i) {
            // stay inside the domain: central differences from inner points only
            double* f = fwd.As<double>() + i * 3 + d;
            double* r = bwd.As<double>() + i * 3 + d;
            const double base = *f;
            *f = std::min(base + h, 1.0);
            *r = std::max(base - h, 0.0);
        }
        const auto pf = meshioplusplus::iga_evaluate(out, 0, fwd);
        const auto pb = meshioplusplus::iga_evaluate(out, 0, bwd);
        for (std::size_t i = 0; i < n; ++i) {
            const double span = fwd.As<double>()[i * 3 + d] - bwd.As<double>()[i * 3 + d];
            for (std::size_t c = 0; c < 3; ++c) {
                const double fd =
                    (pf.mPoints.As<double>()[i * 3 + c] - pb.mPoints.As<double>()[i * 3 + c]) /
                    span;
                EXPECT_NEAR(b.mDerivatives.As<double>()[(i * 3 + d) * 3 + c], fd, 1e-6);
            }
        }
    }
}

TEST(Iga, FieldsRefineWithTheGeometry) {
    NurbsInfo in = one_patch(arc());
    in.mFieldNames = {"u"};
    // a field with two components, one control value per control point
    in.mPatches[0].mFields.push_back(doubles({3, 2}, {1, 10, 2, 20, 3, 30}));
    IgaRefineOptions o;
    o.mSubdivisions = {4, 1, 1};
    o.mElevate = {1, 0, 0};
    const NurbsInfo out = meshioplusplus::iga_refine(in, o);
    ASSERT_EQ(out.mPatches[0].mFields.size(), 1u);
    EXPECT_EQ(out.mPatches[0].mFields[0].Shape()[0], out.mPatches[0].mWeights.Size());
    EXPECT_EQ(out.mPatches[0].mFields[0].Shape()[1], 2u);
    const NDArray u = grid(in.mPatches[0], 17);
    const auto a = meshioplusplus::iga_evaluate(in, 0, u);
    const auto b = meshioplusplus::iga_evaluate(out, 0, u);
    ASSERT_EQ(b.mFields.size(), 1u);
    EXPECT_LT(max_diff(a.mFields[0], b.mFields[0]), 1e-13);
}

TEST(Iga, APolynomialPatchKeepsUnitWeights) {
    const NurbsPatch q = meshioplusplus::iga_elevate_degree(hexahedron(), 1, 1);
    for (std::size_t i = 0; i < q.mWeights.Size(); ++i)
        EXPECT_EQ(q.mWeights.As<double>()[i], 1.0);
}

TEST(Iga, SamplesAHandBuiltModelAtItsLattice) {
    NurbsInfo in = one_patch(annulus());
    in.mFieldNames = {"f"};
    std::vector<double> f;
    for (std::size_t i = 0; i < 6; ++i)
        f.push_back(static_cast<double>(i));
    in.mPatches[0].mFields.push_back(doubles({6, 1}, f));
    IgaSampleOptions so;
    so.mSubdivisions = 2;
    NurbsInfo mapped;
    const Mesh mesh = meshioplusplus::iga_sample(in, so, &mapped);
    // one span in u, one in v, split in 2 each: 4 quad9 cells, 25 nodes
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "quad9");
    EXPECT_EQ(mesh.Cells(0).NumCells(), 4u);
    EXPECT_EQ(mesh.Points().Shape()[0], 25u);
    EXPECT_TRUE(mesh.HasPointData("f"));
    EXPECT_EQ(mapped.mCellPatch.Size(), 4u);
    EXPECT_EQ(mapped.mCellSpan.Shape()[1], 3u);
    // every node of the annulus lies between radius 1 and 2 on z = 0
    for (std::size_t i = 0; i < 25; ++i) {
        const double* x = mesh.Points().As<double>() + i * 3;
        const double r = std::hypot(x[0], x[1]);
        EXPECT_GE(r, 1.0 - 1e-13);
        EXPECT_LE(r, 2.0 + 1e-13);
        EXPECT_EQ(x[2], 0.0);
    }
    so.mSubdivisions = 65;
    EXPECT_THROW(meshioplusplus::iga_sample(in, so), std::invalid_argument);
}

TEST(Iga, ErrorsAreInvalidArguments) {
    const NurbsInfo in = one_patch(arc());
    IgaRefineOptions o;
    o.mSubdivisions = {0, 1, 1};
    EXPECT_THROW(meshioplusplus::iga_refine(in, o), std::invalid_argument);
    o.mSubdivisions = {1, 1, 1};
    o.mElevate = {-1, 0, 0};
    EXPECT_THROW(meshioplusplus::iga_refine(in, o), std::invalid_argument);
    o.mElevate = {31, 0, 0};
    EXPECT_THROW(meshioplusplus::iga_refine(in, o), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_insert_knots(arc(), 1, {0.5}), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_insert_knots(arc(), 0, {0.0}), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_insert_knots(arc(), 0, {0.6, 0.4}), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_elevate_degree(arc(), 0, -1), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_evaluate(in, 1, grid(arc(), 3)), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_evaluate(in, 0, grid(arc(), 3), 2), std::invalid_argument);
    NDArray out_of_range = doubles({1, 1}, {1.5});
    EXPECT_THROW(meshioplusplus::iga_evaluate(in, 0, out_of_range), std::invalid_argument);
    NDArray wrong_width = doubles({1, 2}, {0.5, 0.5});
    EXPECT_THROW(meshioplusplus::iga_evaluate(in, 0, wrong_width), std::invalid_argument);
    NurbsInfo bad = in;
    bad.mPatches[0].mWeights.As<double>()[1] = 0.0;
    EXPECT_THROW(meshioplusplus::iga_refine(bad, IgaRefineOptions{}), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_refine(NurbsInfo{}, IgaRefineOptions{}),
                 std::invalid_argument);
}

// --- MFEM's patch topology ---------------------------------------------------

// The quarter annulus in MFEM's global form (see test_mfem.cpp).
std::string annulus_mesh() {
    const std::string w = "0.70710678118654757";
    const std::string path = mt::temp_path(".mesh");
    std::ofstream(path)
        << "MFEM NURBS mesh v1.0\n"
           "dimension\n2\n"
           "elements\n1\n5 3 0 1 2 3\n"
           "boundary\n0\n"
           "edges\n4\n0 0 1\n0 3 2\n1 0 3\n1 1 2\n"
           "vertices\n4\n"
           "knotvectors\n2\n2 3 0 0 0 1 1 1\n1 2 0 0 1 1\n"
           "weights\n1\n1\n1\n1\n"
        << w << "\n"
        << w
        << "\n"
           "FiniteElementSpace\nFiniteElementCollection: NURBS2\nVDim: 2\nOrdering: 1\n"
           "1 0\n0 1\n0 2\n2 0\n1 1\n2 2\n";
    return path;
}

TEST(Iga, ADefaultSampleOfAReadModelIsTheReadersMesh) {
    NurbsInfo info;
    const Mesh read =
        meshioplusplus::read_mfem(annulus_mesh(), {}, meshioplusplus::ReadOptions{}, info);
    ASSERT_EQ(info.mPatches.size(), 1u);
    NurbsInfo mapped;
    const Mesh sampled = meshioplusplus::iga_sample(info, {}, &mapped);
    EXPECT_EQ(max_diff(read.Points(), sampled.Points()), 0.0);
    EXPECT_EQ(read.NumCellBlocks(), sampled.NumCellBlocks());
    EXPECT_EQ(mapped.mCellPatch.Size(), info.mCellPatch.Size());
    EXPECT_EQ(mapped.mTopology.mNumVertices, info.mTopology.mNumVertices);
    // an explicit subdivision count is the generic sampler: no boundary cells
    IgaSampleOptions so;
    so.mSubdivisions = 1;
    EXPECT_EQ(meshioplusplus::iga_sample(info, so).NumCellBlocks(), 1u);
}

TEST(Iga, ARefinedReadModelWritesBackAsAnMfemNurbsMesh) {
    NurbsInfo info;
    meshioplusplus::read_mfem(annulus_mesh(), {}, meshioplusplus::ReadOptions{}, info);
    IgaRefineOptions o;
    o.mSubdivisions = {3, 3, 1};
    o.mElevate = {1, 1, 0};
    const NurbsInfo refined = meshioplusplus::iga_refine(info, o);
    NurbsInfo mapped;
    const Mesh mesh = meshioplusplus::iga_sample(refined, {}, &mapped);
    const std::string out = mt::temp_path(".mesh");
    meshioplusplus::write_mfem(out, mesh, false, mapped);
    NurbsInfo back;
    const Mesh again = meshioplusplus::read_mfem(out, {}, meshioplusplus::ReadOptions{}, back);
    EXPECT_EQ(max_diff(again.Points(), mesh.Points()), 0.0);
    ASSERT_EQ(back.mPatches.size(), 1u);
    EXPECT_EQ(back.mPatches[0].mDegrees, (std::vector<std::int32_t>{3, 2}));
    EXPECT_EQ(max_diff(back.mPatches[0].mControlPoints, mapped.mPatches[0].mControlPoints), 0.0);
    // refining directions differently would break the shared knot vectors
    o.mSubdivisions = {2, 3, 1};
    EXPECT_THROW(meshioplusplus::iga_refine(info, o), std::invalid_argument);
}

// --- trimmed B-reps ----------------------------------------------------------

// A unit-spaced 4 x 4 degree-1 plate, trimmed by a square loop and a diamond hole.
const char* kPlate = R"({
    "breps": [{
        "brep_id": 1,
        "faces": [{
            "brep_id": 2,
            "swapped_surface_normal": false,
            "surface": {
                "is_trimmed": true,
                "is_rational": false,
                "degrees": [1, 1],
                "knot_vectors": [[0, 4], [0, 4]],
                "control_points": [[1, [0, 0, 0, 1]], [2, [4, 0, 0, 1]],
                                   [3, [0, 4, 0, 1]], [4, [4, 4, 0, 1]]]
            },
            "boundary_loops": [
                {"loop_type": "outer", "trimming_curves": [
                    {"trim_index": 0, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[0, 0, 0, 1]], [[4, 0, 0, 1]]]}},
                    {"trim_index": 1, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[4, 0, 0, 1]], [[4, 4, 0, 1]]]}},
                    {"trim_index": 2, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[4, 4, 0, 1]], [[0, 4, 0, 1]]]}},
                    {"trim_index": 3, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[0, 4, 0, 1]], [[0, 0, 0, 1]]]}}]},
                {"loop_type": "inner", "trimming_curves": [
                    {"trim_index": 4, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 0.25, 0.5, 0.75, 1, 1], "active_range": [0, 1],
                        "control_points": [[[2, 0.8, 0, 1]], [[3.2, 2, 0, 1]], [[2, 3.2, 0, 1]],
                                           [[0.8, 2, 0, 1]], [[2, 0.8, 0, 1]]]}}]}],
            "embedded_edges": []
        }],
        "edges": [],
        "vertices": []
    }]
})";

BrepInfo plate() {
    const std::string path = mt::temp_path(".cad.json");
    std::ofstream(path) << kPlate;
    return meshioplusplus::read_cad_json_info(path);
}

TEST(Iga, ARefinedTrimmedFaceKeepsItsLoopsAndItsInsideTest) {
    const BrepInfo in = plate();
    IgaRefineOptions o;
    o.mSubdivisions = {3, 2, 1};
    o.mElevate = {1, 2, 0};
    const BrepInfo out = meshioplusplus::iga_refine(in, o);
    ASSERT_EQ(out.mFaces.size(), 1u);
    EXPECT_EQ(out.mFaces[0].mControlPointIds.Size(), 0u);
    EXPECT_EQ(out.mFaces[0].mLoops.size(), 2u);
    EXPECT_EQ(out.mPatches[0].mDegrees, (std::vector<std::int32_t>{2, 3}));
    // a grid of parameters: the plate is planar and unrefined geometry is the identity map
    NDArray u(meshioplusplus::DType::Float64, {25, 2});
    for (std::size_t j = 0; j < 5; ++j)
        for (std::size_t i = 0; i < 5; ++i) {
            u.As<double>()[(j * 5 + i) * 2] = 0.1 + 0.8 * static_cast<double>(i);
            u.As<double>()[(j * 5 + i) * 2 + 1] = 0.1 + 0.8 * static_cast<double>(j);
        }
    const auto a = meshioplusplus::iga_evaluate(in, 0, u, 1);
    const auto b = meshioplusplus::iga_evaluate(out, 0, u, 1);
    EXPECT_LT(max_diff(a.mPoints, b.mPoints), 1e-13);
    EXPECT_LT(max_diff(a.mNormals, b.mNormals), 1e-13);
    ASSERT_EQ(a.mInside.Size(), 25u);
    std::size_t inside = 0;
    for (std::size_t i = 0; i < 25; ++i) {
        EXPECT_EQ(a.mInside.As<std::uint8_t>()[i], b.mInside.As<std::uint8_t>()[i]);
        inside += a.mInside.As<std::uint8_t>()[i];
    }
    // the diamond |x-2| + |y-2| < 1.2 holds the centre of the grid and no corner
    EXPECT_EQ(a.mInside.As<std::uint8_t>()[12], 0);
    EXPECT_EQ(a.mInside.As<std::uint8_t>()[0], 1);
    EXPECT_LT(inside, 25u);
    // sampling the refined face is sampling CAD JSON as it does
    IgaSampleOptions so;
    so.mSubdivisions = 8;
    const Mesh mesh = meshioplusplus::iga_sample(out, so);
    EXPECT_GT(mesh.NumCellBlocks(), 0u);
}

TEST(Iga, FacesSharingAPatchMustAgreeOnRational) {
    BrepInfo in = plate();
    in.mFaces.push_back(in.mFaces[0]);
    in.mFaces[1].mRational = true;
    EXPECT_THROW(meshioplusplus::iga_refine(in, IgaRefineOptions{}), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::iga_evaluate(in, 5, grid(in.mPatches[0], 3)),
                 std::invalid_argument);
}

}  // namespace
