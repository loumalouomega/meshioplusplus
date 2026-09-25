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
 * `hessian`). The digests were taken from the implementation before the
 * rewrite; the repeated-run checks hold everywhere, the golden digests on
 * x86-64 GCC/Clang only (no FMA contraction there).
 */

#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "golden_fixtures.hpp"
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/curvature.hpp"
#include "meshioplusplus/operations/hessian.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/shrinkwrap.hpp"
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
};

TEST(OpGoldens, ResultsAreStableAcrossRepeatedRuns) {
    for (const OgCase& c : kOgCases) {
        const std::uint64_t first = c.mRun();
        EXPECT_EQ(c.mRun(), first) << c.mName;
    }
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
