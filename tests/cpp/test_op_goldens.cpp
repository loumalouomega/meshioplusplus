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
 * @file test_op_goldens.cpp
 * @brief Pins the output of the operations roadmap §4's third part rewrites.
 *
 * Welding (`clean`, `merge`), the distance kernel's construction (`compute_sdf`,
 * `shrinkwrap`, `voxelize`, `compute_curvature`, `compute_normals`) and the
 * operations that had no parallel phase (`agglomerate`, `split`, `undo_green`,
 * `hessian`), and the operations its remaining items rewrite (`remesh`,
 * `remesh_volume`, `sample_distance`, `distance_to_surface`, the marching
 * kernel behind `isosurface` and `slice`, `convert_cells`, `interpolate`,
 * `repair`, `diff`, `sobolev_deform`, `gradient`, and the duplicate-cell
 * passes over mixed, polyhedral and ragged blocks). The digests were taken
 * from the implementation before each rewrite; the repeated-run checks hold
 * everywhere, the golden digests on x86-64 GCC/Clang only (no FMA contraction
 * there).
 */

#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "golden_fixtures.hpp"
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/curvature.hpp"
#include "meshioplusplus/operations/diff.hpp"
#include "meshioplusplus/operations/gradient.hpp"
#include "meshioplusplus/operations/hessian.hpp"
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/remesh.hpp"
#include "meshioplusplus/operations/remesh_volume.hpp"
#include "meshioplusplus/operations/repair.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/shrinkwrap.hpp"
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/operations/sobolev_deform.hpp"
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/undo_green.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

namespace {

using namespace golden;
namespace mio = meshioplusplus;

Mesh og_copy_points(const Mesh& rM, const std::vector<std::array<double, 3>>& rExtra) {
    const std::size_t n = rM.NumPoints();
    NDArray pts = NDArray::Uninit(DType::Float64, {n + rExtra.size(), 3});
    std::memcpy(pts.Data(), rM.Points().Data(), n * 3 * sizeof(double));
    for (std::size_t i = 0; i < rExtra.size(); ++i)
        for (int d = 0; d < 3; ++d)
            pts.As<double>()[(n + i) * 3 + static_cast<std::size_t>(d)] = rExtra[i][d];
    Mesh m;
    m.AssignPoints(std::move(pts));
    return m;
}

// Two unwelded copies of a jittered tet cube, the second one's points moved
// by a fraction of the weld tolerance, plus chains of near points (A ~ B,
// B ~ C, A !~ C), duplicate cells and degenerate cells.
Mesh og_weld_input() {
    const Mesh cube = ft_tet_cube(4, 0.2);
    const std::size_t n = cube.NumPoints();
    const double atol = 1.0 / 1024.0;
    std::vector<std::array<double, 3>> extra;
    const double* p = cube.Points().As<double>();
    for (std::size_t i = 0; i < n; ++i)
        extra.push_back({p[i * 3] + 0.25 * atol, p[i * 3 + 1], p[i * 3 + 2] - 0.125 * atol});
    for (std::size_t i = 0; i < n; i += 7) {  // chains along x
        extra.push_back({p[i * 3] + 0.75 * atol, p[i * 3 + 1], p[i * 3 + 2]});
        extra.push_back({p[i * 3] + 1.5 * atol, p[i * 3 + 1], p[i * 3 + 2]});
    }
    Mesh m = og_copy_points(cube, extra);
    const auto cb = cube.Cells(0);
    const std::int64_t* c = cb.Conn().As<std::int64_t>();
    const std::size_t nc = cb.NumCells();
    std::vector<std::int64_t> rows;
    for (std::size_t t = 0; t < nc; ++t)  // original cells
        rows.insert(rows.end(), c + t * 4, c + t * 4 + 4);
    for (std::size_t t = 0; t < nc; ++t)  // the copy, on the shifted points
        for (int k = 0; k < 4; ++k)
            rows.push_back(c[t * 4 + static_cast<std::size_t>(k)] + static_cast<std::int64_t>(n));
    for (std::size_t t = 0; t < nc; t += 5)  // exact duplicates, permuted
        rows.insert(rows.end(), {c[t * 4 + 2], c[t * 4], c[t * 4 + 3], c[t * 4 + 1]});
    for (std::size_t t = 0; t < nc; t += 11)  // repeated node
        rows.insert(rows.end(), {c[t * 4], c[t * 4], c[t * 4 + 2], c[t * 4 + 3]});
    NDArray conn = NDArray::Uninit(DType::Int64, {rows.size() / 4, 4});
    std::memcpy(conn.Data(), rows.data(), rows.size() * sizeof(std::int64_t));
    m.AddCellBlock("tetra", std::move(conn));
    return m;
}

std::uint64_t og_clean(bool weld) {
    mio::CleanOptions o;
    o.weld = weld;
    o.atol = 1.0 / 1024.0;
    auto r = mio::clean(og_weld_input(), o);
    MeshDigest d;
    d.Of(r.mMesh);
    d.Array(r.mPointMap);
    d.Arrays(r.mCellMaps);
    d.U64(static_cast<std::uint64_t>(r.mPointsWelded));
    d.U64(static_cast<std::uint64_t>(r.mPointsRemovedOrphan));
    d.U64(static_cast<std::uint64_t>(r.mCellsDroppedDegenerate));
    d.U64(static_cast<std::uint64_t>(r.mCellsDroppedDuplicate));
    return d.Value();
}

std::uint64_t og_merge(bool dedupe) {
    const Mesh a = ft_tet_cube(4, 0.2);
    const Mesh b = og_weld_input();
    mio::MergeOptions o;
    o.weld = true;
    o.atol = 1.0 / 1024.0;
    o.drop_duplicate_cells = dedupe;
    auto r = mio::merge({&a, &b}, o);
    MeshDigest d;
    d.Of(r.mMesh);
    d.Arrays(r.mPointMaps);
    d.Arrays(r.mCellMaps);
    return d.Value();
}

Mesh og_surface() {
    return mio::extract_surface(ft_tet_cube(6, 0.2));
}

std::uint64_t og_sdf() {
    mio::SdfOptions o;
    o.mResolution = std::array<std::int64_t, 3>{20, 20, 20};
    return ft_mesh(mio::compute_sdf(og_surface(), o).mMesh);
}

std::uint64_t og_shrinkwrap() {
    const Mesh surface = og_surface();
    NDArray pts = NDArray::Uninit(DType::Float64, surface.Points().Shape());
    const double* p = surface.Points().As<double>();
    for (std::size_t i = 0; i < pts.Size(); ++i)
        pts.As<double>()[i] = 0.5 + 1.25 * (p[i] - 0.5);
    Mesh inflated;
    inflated.AssignPoints(std::move(pts));
    NDArray conn =
        NDArray::Uninit(surface.Cells(0).Conn().Dtype(), surface.Cells(0).Conn().Shape());
    std::memcpy(conn.Data(), surface.Cells(0).Conn().Data(), conn.Nbytes());
    inflated.AddCellBlock(surface.Cells(0).Type(), std::move(conn));
    return ft_mesh(mio::shrinkwrap(inflated, surface).mMesh);
}

std::uint64_t og_voxelize() {
    mio::VoxelOptions o;
    o.mResolution = std::array<std::int64_t, 3>{12, 12, 12};
    return ft_mesh(mio::voxelize(og_surface(), o).mMesh);
}

std::uint64_t og_curvature() {
    return ft_mesh(mio::compute_curvature(og_surface()).mMesh);
}

std::uint64_t og_normals() {
    return ft_mesh(mio::compute_normals(og_surface()).mMesh);
}

std::uint64_t og_agglomerate() {
    auto r = mio::agglomerate(ft_tet_cube(5, 0.2));
    MeshDigest d;
    d.Of(r.mMesh);
    d.Array(r.mCellMap);
    return d.Value();
}

std::uint64_t og_split() {
    mio::MergeOptions o;
    o.source_tag = false;
    const Mesh a = ft_tet_cube(3, 0.2), b = ft_hex_grid(3, 0.1, false);
    auto r = mio::split(mio::merge({&a, &b}, o).mMesh, mio::SplitBy::Component);
    MeshDigest d;
    for (const auto& piece : r.mPieces) {
        d.Str(piece.mKey);
        d.Of(piece.mMesh);
        d.Array(piece.mPointMap);
        d.Arrays(piece.mCellMaps);
    }
    return d.Value();
}

std::uint64_t og_undo_green() {
    const Mesh coarse = ft_tet_cube(4, 0.2);
    mio::RefineOptions o;
    for (std::int64_t c = 0; c < 384; c += 9)
        o.mCells.push_back(c);
    o.mRecordHierarchy = true;
    o.mRecordLevels = true;
    const Mesh fine = mio::refine(coarse, o).mMesh;
    auto r = mio::undo_green(coarse, fine);
    MeshDigest d;
    d.Of(r.mMesh);
    d.Arrays(r.mCellMaps);
    d.U64(static_cast<std::uint64_t>(r.mNumGroupsUndone));
    return d.Value();
}

std::uint64_t og_hessian() {
    Mesh m = ft_tet_cube(5, 0.2);
    NDArray u = NDArray::Uninit(DType::Float64, {m.NumPoints()});
    const double* p = m.Points().As<double>();
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        u.As<double>()[i] = p[i * 3] * p[i * 3] + p[i * 3 + 1] * p[i * 3 + 2];
    m.AddPointData("u", std::move(u));
    mio::HessianOptions o;
    o.mArrayName = "u";
    o.mLocation = mio::DataLocation::Point;
    return ft_mesh(mio::hessian(m, o).mMesh);
}

// --- Roadmap §4 part 3b: the operations its remaining items rewrite -------

// A tet cube carrying the smooth point field u = x^2 + y z.
Mesh og_with_field(std::size_t n) {
    Mesh m = ft_tet_cube(n, 0.2);
    NDArray u = NDArray::Uninit(DType::Float64, {m.NumPoints()});
    const double* p = m.Points().As<double>();
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        u.As<double>()[i] = p[i * 3] * p[i * 3] + p[i * 3 + 1] * p[i * 3 + 2];
    m.AddPointData("u", std::move(u));
    return m;
}

void og_f64(MeshDigest& rD, double v) {
    rD.Bytes(&v, sizeof v);
}

// Query points around and inside the unit cube.
NDArray og_query_points() {
    NDArray q = ft_grid_points(5, 0.15);
    for (std::size_t i = 0; i < q.Size(); ++i)
        q.As<double>()[i] = 1.5 * q.As<double>()[i] - 0.25;
    return q;
}

std::uint64_t og_remesh() {
    mio::RemeshOptions o;
    o.mNumClusters = 120;
    auto r = mio::remesh(og_surface(), o);
    MeshDigest d;
    d.Of(r.mMesh);
    d.U64(static_cast<std::uint64_t>(r.mNumIterations));
    d.U64(static_cast<std::uint64_t>(r.mSubdivideApplied));
    return d.Value();
}

// The metric and gradation paths: vertex normals, curvature fits, tensors
// and face quadrics all feed the clustering.
std::uint64_t og_remesh_metrics() {
    MeshDigest d;
    for (auto metric : {mio::RemeshMetric::Quadric, mio::RemeshMetric::Anisotropic}) {
        mio::RemeshOptions o;
        o.mNumClusters = 120;
        o.mMetric = metric;
        o.mGradation = 0.5;
        d.Of(mio::remesh(og_surface(), o).mMesh);
    }
    return d.Value();
}

std::uint64_t og_remesh_volume() {
    mio::RemeshVolumeOptions o;
    o.mCellSize = 0.15;
    auto r = mio::remesh_volume(og_surface(), o);
    MeshDigest d;
    d.Of(r.mMesh);
    d.U64(static_cast<std::uint64_t>(r.mNumVerticesWarped));
    d.U64(static_cast<std::uint64_t>(r.mNumTetsRejected));
    return d.Value();
}

std::uint64_t og_sample_distance() {
    MeshDigest d;
    d.Array(mio::sample_distance(og_surface(), og_query_points()));
    return d.Value();
}

std::uint64_t og_distance_to_surface() {
    Mesh query;
    query.AssignPoints(og_query_points());
    mio::SurfaceDistanceOptions o;
    o.mRecordClosestCell = true;
    return ft_mesh(mio::distance_to_surface(query, og_surface(), o).mMesh);
}

std::uint64_t og_isosurface() {
    mio::IsosurfaceOptions o;
    o.mArrayName = "u";
    o.mIsovalues = {0.3, 0.6};
    o.mRecordParentIds = true;
    return ft_mesh(mio::isosurface(og_with_field(5), o));
}

std::uint64_t og_slice() {
    mio::SliceOptions o;
    o.mOrigin = {0.5, 0.45, 0.47};
    o.mNormal = {0.3, 0.2, 1.0};
    o.mRecordParentIds = true;
    MeshDigest d;
    d.Of(mio::slice(og_with_field(5), o));
    d.Of(mio::slice(ft_hex_grid(4, 0.1, true), o));
    return d.Value();
}

std::uint64_t og_convert_cells() {
    MeshDigest d;
    const auto run = [&](const Mesh& rM, mio::ConvertCellsMode mode) {
        mio::ConvertCellsOptions o;
        o.mMode = mode;
        o.mRecordParentIds = true;
        auto r = mio::convert_cells(rM, o);
        d.Of(r.mMesh);
        d.Array(r.mPointMap);
        d.Arrays(r.mCellMaps);
        return r.mMesh;
    };
    const Mesh hex = ft_hex_grid(3, 0.1, false);
    const Mesh tet = ft_tet_cube(3, 0.2);
    const Mesh hex2 = run(hex, mio::ConvertCellsMode::Elevate);
    run(tet, mio::ConvertCellsMode::Elevate);
    run(hex, mio::ConvertCellsMode::Simplexify);
    run(hex2, mio::ConvertCellsMode::Linearize);
    return d.Value();
}

std::uint64_t og_interpolate() {
    const Mesh src = og_with_field(5);
    Mesh tgt;
    tgt.AssignPoints(og_query_points());
    MeshDigest d;
    for (auto method : {mio::InterpolateMethod::Nearest, mio::InterpolateMethod::Barycentric}) {
        mio::InterpolateOptions o;
        o.mMethod = method;
        o.mExtrapolate = true;
        d.Of(mio::interpolate(src, tgt, o));
    }
    return d.Value();
}

// og_surface with every 7th triangle flipped and every 50th removed (holes).
Mesh og_broken_surface() {
    const Mesh s = og_surface();
    const auto cb = s.Cells(0);
    const std::int64_t* c = cb.Conn().As<std::int64_t>();
    std::vector<std::int64_t> rows;
    for (std::size_t t = 0; t < cb.NumCells(); ++t) {
        if (t % 50 == 3)
            continue;
        if (t % 7 == 0)
            rows.insert(rows.end(), {c[t * 3], c[t * 3 + 2], c[t * 3 + 1]});
        else
            rows.insert(rows.end(), c + t * 3, c + t * 3 + 3);
    }
    NDArray conn = NDArray::Uninit(DType::Int64, {rows.size() / 3, 3});
    std::memcpy(conn.Data(), rows.data(), rows.size() * sizeof(std::int64_t));
    Mesh m = og_copy_points(s, {});
    m.AddCellBlock("triangle", std::move(conn));
    return m;
}

std::uint64_t og_repair() {
    auto r = mio::repair(og_broken_surface());
    MeshDigest d;
    d.Of(r.mMesh);
    d.Array(r.mPointMap);
    d.Arrays(r.mCellMaps);
    for (std::int64_t v : {r.mNumFlipped, r.mNumComponents, r.mNumOrientedOutward,
                           r.mNumHolesFilled, r.mNumFacesAdded, r.mNumPointsAdded})
        d.U64(static_cast<std::uint64_t>(v));
    return d.Value();
}

std::uint64_t og_diff() {
    const Mesh a = og_with_field(4);
    Mesh b = og_with_field(4);
    NDArray moved = NDArray::Uninit(DType::Float64, b.Points().Shape());
    std::memcpy(moved.Data(), b.Points().Data(), moved.Nbytes());
    for (std::size_t i = 0; i < moved.Size(); i += 5)
        moved.As<double>()[i] += 1e-7 * static_cast<double>(i % 11);
    b.AssignPoints(std::move(moved));
    MeshDigest d;
    for (bool unordered : {false, true}) {
        mio::DiffOptions o;
        o.unordered = unordered;
        o.atol = 1e-9;
        const mio::DiffReport r = mio::diff(a, b, o);
        d.U64(static_cast<std::uint64_t>(r.mVerdict));
        d.U64(static_cast<std::uint64_t>(r.mCorrespondenceFailed));
        og_f64(d, r.mPoints.mMaxAbsError);
        og_f64(d, r.mPoints.mMaxRelError);
        d.U64(static_cast<std::uint64_t>(r.mPoints.mWorstIndex));
        d.U64(static_cast<std::uint64_t>(r.mPoints.mNumExceeding));
        for (const std::string& rMsg : r.mMessages)
            d.Str(rMsg);
    }
    return d.Value();
}

std::uint64_t og_sobolev() {
    Mesh m = og_surface();
    NDArray disp = NDArray::Uninit(DType::Float64, {m.NumPoints(), 3});
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            disp.As<double>()[i * 3 + k] = 0.01 * static_cast<double>((i * 7 + k * 3) % 13) - 0.06;
    m.AddPointData("d", std::move(disp));
    mio::SobolevOptions o;
    o.mArrayName = "d";
    o.mLengthScale = 0.3;
    auto r = mio::sobolev_deform(m, o);
    MeshDigest d;
    d.Of(r.mMesh);
    d.U64(static_cast<std::uint64_t>(r.mNumIterations));
    og_f64(d, r.mResidual);
    return d.Value();
}

std::uint64_t og_gradient() {
    MeshDigest d;
    const Mesh m = og_with_field(5);
    for (auto method : {mio::GradientMethod::GreenGauss, mio::GradientMethod::LeastSquares})
        for (auto loc : {mio::DataLocation::Cell, mio::DataLocation::Point}) {
            mio::GradientOptions o;
            o.mArrayName = "u";
            o.mMethod = method;
            o.mLocation = loc;
            d.Of(mio::gradient(m, o).mMesh);
        }
    return d.Value();
}

// Two unwelded copies of a mixed hexahedron/polyhedron grid plus a ragged
// polygon block, so welding makes every cell a duplicate on each path.
Mesh og_mixed_doubled() {
    Mesh g = ft_hex_grid(3, 0.1, true);
    g.AddPolygonBlock("polygon",
                      {{0, 1, 5, 4}, {1, 2, 6}, {2, 3, 7, 6, 10}, {4, 5, 9, 8}, {0, 4, 1}});
    mio::MergeOptions o;
    o.source_tag = false;
    return mio::merge({&g, &g}, o).mMesh;
}

std::uint64_t og_clean_mixed() {
    mio::CleanOptions o;
    o.weld = true;
    o.atol = 1.0 / 1024.0;
    auto r = mio::clean(og_mixed_doubled(), o);
    MeshDigest d;
    d.Of(r.mMesh);
    d.Array(r.mPointMap);
    d.Arrays(r.mCellMaps);
    d.U64(static_cast<std::uint64_t>(r.mCellsDroppedDuplicate));
    return d.Value();
}

std::uint64_t og_merge_mixed() {
    const Mesh a = og_mixed_doubled();
    Mesh b = ft_hex_grid(3, 0.1, true);
    b.AddPolygonBlock("polygon", {{4, 5, 1, 0}, {6, 2, 1}, {2, 3, 7, 6, 10}});
    mio::MergeOptions o;
    o.weld = true;
    o.atol = 1.0 / 1024.0;
    o.drop_duplicate_cells = true;
    auto r = mio::merge({&a, &b}, o);
    MeshDigest d;
    d.Of(r.mMesh);
    d.Arrays(r.mPointMaps);
    d.Arrays(r.mCellMaps);
    return d.Value();
}

struct OgCase {
    const char* mName;
    std::uint64_t (*mRun)();
    std::uint64_t mGolden;
};

const OgCase kOgCases[] = {
    {"clean_weld", [] { return og_clean(true); }, 0x5f6d105fde1e0997ull},
    {"clean_dedupe", [] { return og_clean(false); }, 0xfda9903f09cf8bc8ull},
    {"merge_weld", [] { return og_merge(false); }, 0x8d6b12cea13b9286ull},
    {"merge_dedupe", [] { return og_merge(true); }, 0x12ade067e19545bfull},
    {"sdf", og_sdf, 0x175cf3b488dec5f1ull},
    {"shrinkwrap", og_shrinkwrap, 0xb0848f2913e6c5c1ull},
    {"voxelize", og_voxelize, 0x6ba9809e780fc0adull},
    {"curvature", og_curvature, 0x57e791c665c3c096ull},
    {"normals", og_normals, 0xdd92696ea89a5ad9ull},
    {"agglomerate", og_agglomerate, 0x7c40650413522165ull},
    {"split", og_split, 0x258b588246f2ea38ull},
    {"undo_green", og_undo_green, 0x994d7afcbc6f56d9ull},
    {"hessian", og_hessian, 0x68a9b378a927d1b3ull},
    {"remesh", og_remesh, 0x79378985b87aac45ull},
    {"remesh_metrics", og_remesh_metrics, 0x8e4ec97dbb94abbfull},
    {"remesh_volume", og_remesh_volume, 0x7494ff7f3f992c87ull},
    {"sample_distance", og_sample_distance, 0x95506a1b94a08374ull},
    {"distance_to_surface", og_distance_to_surface, 0xd6ade37c6fb97ec8ull},
    {"isosurface", og_isosurface, 0x43e8aaa22e1a29f8ull},
    {"slice", og_slice, 0x125f56f49a888578ull},
    {"convert_cells", og_convert_cells, 0x2128f4d74274f7ddull},
    {"interpolate", og_interpolate, 0x7a243a2a64708882ull},
    {"repair", og_repair, 0x9d08b9e5c5d7c69ull},
    {"diff", og_diff, 0x5c8d33b762693b9full},
    {"sobolev_deform", og_sobolev, 0x489219385fb3f63aull},
    {"gradient", og_gradient, 0xc3a947740b3221bcull},
    {"clean_mixed", og_clean_mixed, 0x7a70d94cd8d6e8feull},
    {"merge_mixed", og_merge_mixed, 0x275bb867891ac42bull},
};

TEST(OpGoldens, ResultsAreStableAcrossRepeatedRuns) {
    for (const OgCase& c : kOgCases) {
        const std::uint64_t first = c.mRun();
        EXPECT_EQ(c.mRun(), first) << c.mName;
    }
}

// The fixtures reach the paths they pin: duplicates on every block kind,
// holes and flipped triangles for `repair`, a non-empty cut for marching.
TEST(OpGoldens, FixturesExerciseTheirPaths) {
    mio::CleanOptions co;
    co.weld = true;
    co.atol = 1.0 / 1024.0;
    const auto cleaned = mio::clean(og_mixed_doubled(), co);
    const Mesh g = ft_hex_grid(3, 0.1, true);
    EXPECT_EQ(cleaned.mCellsDroppedDuplicate,
              static_cast<std::int64_t>(g.Cells(0).NumCells() + g.Cells(1).NumCells() + 5));
    ASSERT_EQ(cleaned.mMesh.NumCellBlocks(), 3u);
    EXPECT_TRUE(cleaned.mMesh.Cells(1).IsPolyhedron());
    EXPECT_TRUE(cleaned.mMesh.Cells(2).IsRagged());
    const auto repaired = mio::repair(og_broken_surface());
    EXPECT_GT(repaired.mNumFlipped, 0);
    EXPECT_GT(repaired.mNumHolesFilled, 0);
    mio::IsosurfaceOptions io;
    io.mArrayName = "u";
    io.mIsovalues = {0.3};
    EXPECT_GT(mio::isosurface(og_with_field(5), io).NumPoints(), 0u);
}

#ifdef GOLDEN_PINNED
TEST(OpGoldens, ResultsMatchThePreviousImplementation) {
    for (const OgCase& c : kOgCases) {
        const std::uint64_t got = c.mRun();
        EXPECT_EQ(got, c.mGolden) << c.mName << ": 0x" << std::hex << got;
    }
}
#endif

}  // namespace
