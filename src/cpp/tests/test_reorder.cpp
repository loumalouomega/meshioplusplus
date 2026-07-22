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
// GoogleTest suite for the mesh renumbering operation (operations/reorder.hpp).
// Runs under every mesh backend via the cpp-tests matrix.

// System includes
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "mesh_fixtures.hpp"

using meshioplusplus::compute_bandwidth;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::reorder;
using meshioplusplus::ReorderMethod;
using meshioplusplus::ReorderResult;

namespace {

// A path of `n` nodes with an even-first node labelling that gives a large
// bandwidth (~n/2): spatial position p sits at x=p but is labelled so that
// spatially-adjacent nodes are ~n/2 apart in index. RCM should collapse it.
Mesh bad_path_mesh(int n) {
    std::vector<std::int64_t> label(n);
    int half = (n + 1) / 2;
    for (int p = 0; p < n; ++p)
        label[p] = (p % 2 == 0) ? (p / 2) : (half + p / 2);
    std::vector<std::vector<double>> pts(n, std::vector<double>(3, 0.0));
    for (int p = 0; p < n; ++p)
        pts[static_cast<std::size_t>(label[p])] = {static_cast<double>(p), 0.0, 0.0};
    std::vector<std::vector<std::int64_t>> edges;
    for (int p = 0; p + 1 < n; ++p)
        edges.push_back({label[p], label[p + 1]});
    return mt::make_mesh(pts, "line", edges);
}

// An nx*ny grid of quads, row-major node numbering. `stride` relabels the
// nodes by (k*stride)%n (a coprime stride scrambles the numbering -> large
// bandwidth); stride==1 leaves the natural row-major numbering.
Mesh grid_mesh(int nx, int ny, std::int64_t stride = 1) {
    const std::int64_t n = static_cast<std::int64_t>((nx + 1) * (ny + 1));
    auto relabel = [n, stride](std::int64_t k) { return (k * stride) % n; };
    std::vector<std::vector<double>> pts(static_cast<std::size_t>(n), std::vector<double>(3, 0.0));
    std::int64_t k = 0;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i, ++k)
            pts[static_cast<std::size_t>(relabel(k))] = {static_cast<double>(i),
                                                         static_cast<double>(j), 0.0};
    auto idx = [nx, &relabel](int i, int j) {
        return relabel(static_cast<std::int64_t>(j * (nx + 1) + i));
    };
    std::vector<std::vector<std::int64_t>> quads;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            quads.push_back({idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)});
    return mt::make_mesh(pts, "quad", quads);
}

// Assert the result is a pure permutation of the input: node_perm is a
// bijection, points move by it, and connectivity un-permutes back to the input.
void check_pure_permutation(const Mesh& rIn, const ReorderResult& rRes) {
    const std::size_t n = rIn.NumPoints();
    ASSERT_EQ(rRes.mNodePermutation.Size(), n);
    const std::int64_t* np = rRes.mNodePermutation.As<std::int64_t>();

    std::vector<char> seen(n, 0);
    std::vector<std::int64_t> inv(n, -1);
    for (std::size_t i = 0; i < n; ++i) {
        const std::int64_t v = np[i];
        ASSERT_GE(v, 0);
        ASSERT_LT(static_cast<std::size_t>(v), n);
        ASSERT_EQ(seen[static_cast<std::size_t>(v)], 0) << "node_perm is not a bijection";
        seen[static_cast<std::size_t>(v)] = 1;
        inv[static_cast<std::size_t>(v)] = static_cast<std::int64_t>(i);
    }

    // points: out.points[np[i]] == in.points[i]
    const NDArray& pin = rIn.Points();
    const NDArray& pout = rRes.mMesh.Points();
    const std::size_t d = meshioplusplus::detail::cols(pin);
    ASSERT_EQ(rRes.mMesh.NumPoints(), n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < d; ++j) {
            const double a = meshioplusplus::detail::read_double(pin, i * d + j);
            const double b =
                meshioplusplus::detail::read_double(pout, static_cast<std::size_t>(np[i]) * d + j);
            EXPECT_DOUBLE_EQ(a, b) << "point " << i << " comp " << j;
        }

    // connectivity: map out's node labels back through inv, compare multisets.
    std::map<std::string, std::multiset<std::vector<std::int64_t>>> back;
    for (const auto cb : rRes.mMesh.CellRange()) {
        const NDArray& conn = cb.Conn();
        const std::size_t nc = cb.NumCells();
        const std::size_t k = meshioplusplus::detail::cols(conn);
        for (std::size_t r = 0; r < nc; ++r) {
            std::vector<std::int64_t> row(k);
            for (std::size_t j = 0; j < k; ++j)
                row[j] = inv[static_cast<std::size_t>(
                    meshioplusplus::detail::read_int(conn, r * k + j))];
            back[cb.Type()].insert(std::move(row));
        }
    }
    EXPECT_EQ(mt::cell_rows(rIn), back);
}

}  // namespace

TEST(Reorder, MethodFromName) {
    EXPECT_EQ(meshioplusplus::reorder_method_from_name("rcm"), ReorderMethod::RCM);
    EXPECT_EQ(meshioplusplus::reorder_method_from_name("morton"), ReorderMethod::Morton);
    EXPECT_EQ(meshioplusplus::reorder_method_from_name("hilbert"), ReorderMethod::Hilbert);
    EXPECT_THROW(meshioplusplus::reorder_method_from_name("nope"), std::invalid_argument);
}

TEST(Reorder, PurePermutationAllMethods) {
    const Mesh fixtures[] = {mt::tri_mesh(), mt::tet_mesh(), grid_mesh(4, 3), bad_path_mesh(16)};
    for (const Mesh& m : fixtures)
        for (ReorderMethod method :
             {ReorderMethod::RCM, ReorderMethod::Morton, ReorderMethod::Hilbert}) {
            ReorderResult res = reorder(m, method);
            check_pure_permutation(m, res);
        }
}

TEST(Reorder, RcmReducesBandwidth) {
    Mesh m = bad_path_mesh(64);
    const std::int64_t before = compute_bandwidth(m);
    ReorderResult res = reorder(m, ReorderMethod::RCM);
    const std::int64_t after = compute_bandwidth(res.mMesh);
    EXPECT_LT(after, before) << "RCM should strictly reduce this bad numbering";
    EXPECT_LE(after, 2);  // a path collapses to bandwidth ~1
}

TEST(Reorder, RcmReducesBandwidthOnShuffledGrid) {
    Mesh m = grid_mesh(6, 5, /*stride=*/13);  // scrambled numbering
    const std::int64_t before = compute_bandwidth(m);
    ReorderResult res = reorder(m, ReorderMethod::RCM);
    const std::int64_t after = compute_bandwidth(res.mMesh);
    EXPECT_LT(after, before);
    EXPECT_LE(after, 16);
}

TEST(Reorder, SfcLocalityOnGrid) {
    // A space-filling curve keeps spatially-close nodes close in index: the
    // mean spatial step between consecutive new indices stays far below the
    // ~4-5 a random ordering would give (grid spacing 1).
    Mesh m = grid_mesh(8, 8);
    const std::size_t n = m.NumPoints();
    for (ReorderMethod method : {ReorderMethod::Morton, ReorderMethod::Hilbert}) {
        ReorderResult res = reorder(m, method);
        check_pure_permutation(m, res);
        const std::int64_t* np = res.mNodePermutation.As<std::int64_t>();
        std::vector<std::int64_t> inv(n);
        for (std::size_t i = 0; i < n; ++i)
            inv[static_cast<std::size_t>(np[i])] = static_cast<std::int64_t>(i);
        const NDArray& pts = m.Points();
        const std::size_t d = meshioplusplus::detail::cols(pts);
        double total = 0.0;
        for (std::size_t k = 1; k < n; ++k) {
            double sq = 0.0;
            for (std::size_t j = 0; j < d; ++j) {
                const double a = meshioplusplus::detail::read_double(
                    pts, static_cast<std::size_t>(inv[k]) * d + j);
                const double b = meshioplusplus::detail::read_double(
                    pts, static_cast<std::size_t>(inv[k - 1]) * d + j);
                sq += (a - b) * (a - b);
            }
            total += std::sqrt(sq);
        }
        EXPECT_LT(total / static_cast<double>(n - 1), 3.0);
    }
}

TEST(Reorder, PreservesPointData) {
    Mesh m = grid_mesh(3, 3);
    const std::size_t n = m.NumPoints();
    NDArray pd = NDArray::Uninit(meshioplusplus::DType::Float64, {n, 1});
    for (std::size_t i = 0; i < n; ++i)
        pd.As<double>()[i] = static_cast<double>(i) + 0.5;
    m.AddPointData("tag", std::move(pd));

    ReorderResult res = reorder(m, ReorderMethod::RCM);
    ASSERT_TRUE(res.mMesh.HasPointData("tag"));
    const std::int64_t* np = res.mNodePermutation.As<std::int64_t>();
    const NDArray& out_pd = res.mMesh.PointData("tag");
    for (std::size_t i = 0; i < n; ++i)
        EXPECT_DOUBLE_EQ(out_pd.As<double>()[static_cast<std::size_t>(np[i])],
                         static_cast<double>(i) + 0.5);
}

TEST(Reorder, CellPermutationsPresent) {
    Mesh m = mt::tri_quad_mesh();
    ReorderResult res = reorder(m, ReorderMethod::RCM);
    EXPECT_EQ(res.mCellPermutations.size(), m.NumCellBlocks());
    std::size_t b = 0;
    for (const auto cb : m.CellRange())
        EXPECT_EQ(res.mCellPermutations[b++].Size(), cb.NumCells());
}
