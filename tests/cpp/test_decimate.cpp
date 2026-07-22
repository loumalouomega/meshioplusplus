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

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/skin.hpp"

namespace {

using meshioplusplus::compute_stats;
using meshioplusplus::decimate;
using meshioplusplus::DecimateOptions;
using meshioplusplus::DecimatePlacement;
using meshioplusplus::DecimateResult;
using meshioplusplus::extract_skin;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::refine;
using meshioplusplus::RefineOptions;

DecimateOptions target_opts(std::int64_t faces) {
    DecimateOptions o;
    o.mTargetFaces = faces;
    return o;
}

DecimateOptions ratio_opts(double ratio) {
    DecimateOptions o;
    o.mTargetRatio = ratio;
    return o;
}

// A regular (n+1)x(n+1) triangulated grid over the unit square (z = 0),
// 2*n*n triangles with a fixed diagonal.
Mesh grid_mesh(int n) {
    std::vector<std::vector<double>> pts;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i)
            pts.push_back({static_cast<double>(i) / n, static_cast<double>(j) / n, 0.0});
    std::vector<std::vector<std::int64_t>> cells;
    auto v = [n](int i, int j) { return static_cast<std::int64_t>(j * (n + 1) + i); };
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            cells.push_back({v(i, j), v(i + 1, j), v(i + 1, j + 1)});
            cells.push_back({v(i, j), v(i + 1, j + 1), v(i, j + 1)});
        }
    return mt::make_mesh(std::move(pts), "triangle", std::move(cells));
}

// The grid lifted onto the strictly convex bump z = x^2 + y^2. Unlike a ridged
// (ruled) sheet, which still admits zero-error collapses ALONG each crease,
// strict convexity gives every possible collapse a genuinely positive error.
Mesh bump_grid_mesh(int n) {
    Mesh flat = grid_mesh(n);
    std::vector<std::vector<double>> pts;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) {
            const double x = static_cast<double>(i) / n;
            const double y = static_cast<double>(j) / n;
            pts.push_back({x, y, x * x + y * y});
        }
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    NDArray conn = flat.Cells(0).Conn();
    conn.MakeOwned();
    m.AddCellBlock("triangle", std::move(conn));
    return m;
}

// An octahedron refined `levels` times with every point pushed onto the unit
// sphere. Same recipe as tests/python/test_decimate.py's `_sphere`.
Mesh sphere_mesh(int levels) {
    Mesh octa = mt::make_mesh(
        {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}, "triangle",
        {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4}, {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}});
    RefineOptions ro;
    ro.mLevels = levels;
    Mesh fine = refine(octa, ro).mMesh;
    const std::size_t n = fine.NumPoints();
    NDArray pts = NDArray::Uninit(meshioplusplus::DType::Float64, {n, 3});
    double* dst = pts.As<double>();
    for (std::size_t i = 0; i < n; ++i) {
        double x = meshioplusplus::detail::read_double(fine.Points(), i * 3);
        double y = meshioplusplus::detail::read_double(fine.Points(), i * 3 + 1);
        double z = meshioplusplus::detail::read_double(fine.Points(), i * 3 + 2);
        const double len = std::sqrt(x * x + y * y + z * z);
        dst[i * 3] = x / len;
        dst[i * 3 + 1] = y / len;
        dst[i * 3 + 2] = z / len;
    }
    Mesh m;
    m.AssignPoints(std::move(pts));
    NDArray conn = fine.Cells(0).Conn();
    conn.MakeOwned();
    m.AddCellBlock("triangle", std::move(conn));
    return m;
}

std::size_t total_cells(const Mesh& rMesh) {
    std::size_t total = 0;
    for (const auto cb : rMesh.CellRange())
        total += cb.NumCells();
    return total;
}

std::array<double, 3> point_at(const Mesh& rMesh, std::size_t i) {
    const std::size_t dim = rMesh.PointDim();
    std::array<double, 3> p = {0.0, 0.0, 0.0};
    for (std::size_t d = 0; d < dim && d < 3; ++d)
        p[d] = meshioplusplus::detail::read_double(rMesh.Points(), i * dim + d);
    return p;
}

bool has_point(const Mesh& rMesh, double x, double y, double z) {
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i) {
        const std::array<double, 3> p = point_at(rMesh, i);
        if (p[0] == x && p[1] == y && p[2] == z)
            return true;
    }
    return false;
}

// Every edge of a decimated surface must be used by at most two faces, no face
// may repeat a corner, and no two faces may share the same corner set.
void expect_manifold(const Mesh& rMesh) {
    std::map<std::pair<std::int64_t, std::int64_t>, int> edge_use;
    std::set<std::array<std::int64_t, 3>> face_keys;
    for (const auto cb : rMesh.CellRange()) {
        const NDArray& conn = cb.Conn();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::array<std::int64_t, 3> corners;
            for (std::size_t k = 0; k < 3; ++k)
                corners[k] = meshioplusplus::detail::read_int(conn, c * 3 + k);
            EXPECT_NE(corners[0], corners[1]);
            EXPECT_NE(corners[1], corners[2]);
            EXPECT_NE(corners[0], corners[2]);
            std::array<std::int64_t, 3> key = corners;
            std::sort(key.begin(), key.end());
            EXPECT_TRUE(face_keys.insert(key).second) << "duplicate face";
            for (std::size_t k = 0; k < 3; ++k) {
                std::int64_t a = corners[k];
                std::int64_t b = corners[(k + 1) % 3];
                if (b < a)
                    std::swap(a, b);
                ++edge_use[{a, b}];
            }
        }
    }
    for (const auto& [key, count] : edge_use)
        EXPECT_LE(count, 2) << "non-manifold edge (" << key.first << ", " << key.second << ")";
}

TEST(Decimate, PlanarGridCollapsesToBoundaryOnly) {
    const Mesh grid = grid_mesh(4);  // 32 faces, 16 boundary + 9 interior points
    const DecimateResult res = decimate(grid, target_opts(1));

    // The boundary is fully pinned, so the reachable minimum keeps all 16
    // outline vertices; the interior should be gone (or nearly so).
    EXPECT_LE(total_cells(res.mMesh), 16u);
    EXPECT_LT(total_cells(res.mMesh), 32u);
    EXPECT_EQ(res.mFacesRemoved, 32 - static_cast<std::int64_t>(total_cells(res.mMesh)));

    // The outline is preserved exactly: every input boundary point survives at
    // its exact coordinate, and the area of the sheet is unchanged.
    for (int k = 0; k <= 4; ++k) {
        const double t = static_cast<double>(k) / 4;
        EXPECT_TRUE(has_point(res.mMesh, t, 0.0, 0.0));
        EXPECT_TRUE(has_point(res.mMesh, t, 1.0, 0.0));
        EXPECT_TRUE(has_point(res.mMesh, 0.0, t, 0.0));
        EXPECT_TRUE(has_point(res.mMesh, 1.0, t, 0.0));
    }
    const double area = compute_stats(res.mMesh).mTotalArea;
    EXPECT_NEAR(area, 1.0, 1e-9);
    expect_manifold(res.mMesh);
}

TEST(Decimate, TargetLandsWithinOneCollapse) {
    const Mesh grid = grid_mesh(4);
    for (const std::int64_t target : {30, 24, 21, 20, 16}) {
        const DecimateResult res = decimate(grid, target_opts(target));
        const std::int64_t faces = static_cast<std::int64_t>(total_cells(res.mMesh));
        EXPECT_LE(faces, target) << "target " << target;
        EXPECT_GE(faces, target - 1) << "target " << target;
    }
}

TEST(Decimate, RatioMatchesFaceTarget) {
    const Mesh grid = grid_mesh(4);
    const DecimateResult by_ratio = decimate(grid, ratio_opts(0.625));  // 20 of 32
    const DecimateResult by_faces = decimate(grid, target_opts(20));
    EXPECT_EQ(total_cells(by_ratio.mMesh), total_cells(by_faces.mMesh));
}

TEST(Decimate, SphereKeepsShapeAndOrientation) {
    const Mesh sphere = sphere_mesh(3);  // 512 faces
    const double area_in = compute_stats(sphere).mTotalArea;

    const DecimateResult res = decimate(sphere, ratio_opts(0.25));
    const std::int64_t faces = static_cast<std::int64_t>(total_cells(res.mMesh));
    EXPECT_LE(faces, 128);
    EXPECT_GE(faces, 127);

    const auto stats = compute_stats(res.mMesh);
    EXPECT_NEAR(stats.mTotalArea, area_in, 0.25 * area_in);
    for (int d = 0; d < 3; ++d) {
        EXPECT_NEAR(stats.mBBoxMin[d], -1.0, 0.25);
        EXPECT_NEAR(stats.mBBoxMax[d], 1.0, 0.25);
    }

    // No face may point inward: the decimated sphere stays star-shaped around
    // the origin, so a flipped normal would show as a negative dot product.
    for (const auto cb : res.mMesh.CellRange()) {
        const NDArray& conn = cb.Conn();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::array<std::array<double, 3>, 3> p;
            for (std::size_t k = 0; k < 3; ++k)
                p[k] = point_at(res.mMesh, static_cast<std::size_t>(
                                               meshioplusplus::detail::read_int(conn, c * 3 + k)));
            const double e1x = p[1][0] - p[0][0], e1y = p[1][1] - p[0][1], e1z = p[1][2] - p[0][2];
            const double e2x = p[2][0] - p[0][0], e2y = p[2][1] - p[0][1], e2z = p[2][2] - p[0][2];
            const double nx = e1y * e2z - e1z * e2y;
            const double ny = e1z * e2x - e1x * e2z;
            const double nz = e1x * e2y - e1y * e2x;
            const double cx = (p[0][0] + p[1][0] + p[2][0]) / 3.0;
            const double cy = (p[0][1] + p[1][1] + p[2][1]) / 3.0;
            const double cz = (p[0][2] + p[1][2] + p[2][2]) / 3.0;
            EXPECT_GT(nx * cx + ny * cy + nz * cz, 0.0);
        }
    }
    expect_manifold(res.mMesh);
}

TEST(Decimate, MaxErrorCollapsesFlatButNotRidged) {
    // On the flat grid every collapse is (numerically) error-free, so the
    // max_error criterion decimates it all the way to the pinned frame.
    const Mesh flat = grid_mesh(4);
    DecimateOptions opts;
    opts.mMaxError = 1e-9;
    const DecimateResult res_flat = decimate(flat, opts);
    EXPECT_LE(total_cells(res_flat.mMesh), 16u);
    EXPECT_LE(res_flat.mMaxErrorApplied, 1e-9);

    // On the strictly convex bump every collapse has a genuinely positive
    // error, so a tiny budget stops decimation entirely.
    DecimateOptions bump_opts;
    bump_opts.mMaxError = 1e-12;
    bump_opts.mPreserveFeatures = false;
    const DecimateResult res_bump = decimate(bump_grid_mesh(4), bump_opts);
    EXPECT_EQ(res_bump.mFacesRemoved, 0);
    EXPECT_EQ(total_cells(res_bump.mMesh), 32u);
}

TEST(Decimate, CubeIsFullyPinnedByFeatures) {
    // Every vertex of a plain cube skin is a corner, so with features on the
    // queue is exhausted with nothing collapsible -- and that is a normal
    // return, not an error.
    const Mesh cube = extract_skin(mt::hex_mesh(), /*linearize=*/true);
    const DecimateResult res = decimate(cube, target_opts(1));
    EXPECT_EQ(res.mFacesRemoved, 0);
    EXPECT_EQ(total_cells(res.mMesh), 12u);  // 6 quads triangulated
    for (int x = 0; x <= 1; ++x)
        for (int y = 0; y <= 1; ++y)
            for (int z = 0; z <= 1; ++z)
                EXPECT_TRUE(has_point(res.mMesh, x, y, z));
}

TEST(Decimate, RefinedCubeKeepsCornersAndCreases) {
    RefineOptions two_levels;
    two_levels.mLevels = 2;
    const Mesh cube = refine(extract_skin(mt::hex_mesh(), /*linearize=*/true), two_levels).mMesh;
    const std::size_t faces_in = 2 * total_cells(cube);  // quads triangulate to 2

    const DecimateResult res = decimate(cube, ratio_opts(0.25));
    EXPECT_LT(total_cells(res.mMesh), faces_in);

    // Every input vertex on a cube edge (two coordinates in {0, 1}) is a
    // feature vertex and must survive at its exact position.
    for (std::size_t i = 0; i < cube.NumPoints(); ++i) {
        const std::array<double, 3> p = point_at(cube, i);
        int extreme = 0;
        for (int d = 0; d < 3; ++d)
            if (p[d] == 0.0 || p[d] == 1.0)
                ++extreme;
        if (extreme >= 2)
            EXPECT_TRUE(has_point(res.mMesh, p[0], p[1], p[2]))
                << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")";
    }
    expect_manifold(res.mMesh);
}

TEST(Decimate, FrozenVertexSurvives) {
    const Mesh grid = grid_mesh(4);
    DecimateOptions opts = target_opts(1);
    opts.mFrozen.assign(grid.NumPoints(), 0);
    opts.mFrozen[2 * 5 + 2] = 1;  // the centre vertex (0.5, 0.5)
    const DecimateResult res = decimate(grid, opts);
    EXPECT_TRUE(has_point(res.mMesh, 0.5, 0.5, 0.0));

    DecimateOptions bad = target_opts(1);
    bad.mFrozen.assign(3, 1);
    EXPECT_THROW(decimate(grid, bad), std::invalid_argument);
}

TEST(Decimate, MapsAreConsistentAcrossBlocks) {
    // All seven points of the strip are on its boundary, so nothing collapses;
    // the quad block still comes back triangulated, blocks stay 1:1.
    const Mesh strip = mt::tri_quad_mesh();
    const DecimateResult res = decimate(strip, target_opts(1));
    ASSERT_EQ(res.mMesh.NumCellBlocks(), 3u);
    EXPECT_EQ(std::string(res.mMesh.Cells(0).Type()), "triangle");
    EXPECT_EQ(std::string(res.mMesh.Cells(1).Type()), "triangle");
    EXPECT_EQ(res.mMesh.Cells(0).NumCells(), 2u);
    EXPECT_EQ(res.mMesh.Cells(1).NumCells(), 2u);  // 1 quad -> 2 triangles
    EXPECT_EQ(res.mMesh.Cells(2).NumCells(), 1u);

    ASSERT_EQ(res.mPointMap.Size(), strip.NumPoints());
    const std::int64_t* pm = res.mPointMap.As<std::int64_t>();
    for (std::size_t i = 0; i < strip.NumPoints(); ++i)
        EXPECT_EQ(pm[i], static_cast<std::int64_t>(i));

    ASSERT_EQ(res.mCellMaps.size(), 3u);
    EXPECT_EQ(res.mCellMaps[0].Size(), 2u);  // input triangle cells
    EXPECT_EQ(res.mCellMaps[1].Size(), 1u);  // ONE input quad cell
    EXPECT_EQ(res.mCellMaps[1].As<std::int64_t>()[0], 0);
    EXPECT_EQ(res.mCellMaps[2].Size(), 1u);
}

TEST(Decimate, MapsPointAtSurvivors) {
    const Mesh grid = grid_mesh(4);
    const DecimateResult res = decimate(grid, target_opts(16));
    ASSERT_EQ(res.mPointMap.Size(), grid.NumPoints());
    const std::int64_t* pm = res.mPointMap.As<std::int64_t>();
    const std::int64_t nout = static_cast<std::int64_t>(res.mMesh.NumPoints());
    for (std::size_t i = 0; i < grid.NumPoints(); ++i) {
        EXPECT_GE(pm[i], 0);  // every point has a live survivor
        EXPECT_LT(pm[i], nout);
    }
    ASSERT_EQ(res.mCellMaps.size(), 1u);
    const std::int64_t* cm = res.mCellMaps[0].As<std::int64_t>();
    const std::int64_t ncells = static_cast<std::int64_t>(res.mMesh.Cells(0).NumCells());
    std::int64_t mapped = 0;
    for (std::size_t c = 0; c < res.mCellMaps[0].Size(); ++c) {
        EXPECT_GE(cm[c], -1);
        EXPECT_LT(cm[c], ncells);
        if (cm[c] >= 0)
            ++mapped;
    }
    EXPECT_GT(mapped, 0);
    EXPECT_EQ(res.mPointsRemoved,
              static_cast<std::int64_t>(grid.NumPoints() - res.mMesh.NumPoints()));
}

TEST(Decimate, DataIsCarriedAndBlended) {
    Mesh sphere = sphere_mesh(2);
    const std::size_t n = sphere.NumPoints();
    std::vector<double> t_vals(n);
    std::vector<std::int32_t> id_vals(n);
    for (std::size_t i = 0; i < n; ++i) {
        t_vals[i] = point_at(sphere, i)[0];  // T = x, in [-1, 1]
        id_vals[i] = static_cast<std::int32_t>(i % 7);
    }
    sphere.AddPointData("T", mt::data_array(t_vals));
    {
        NDArray ids = mt::int_data_array(id_vals);
        sphere.AddPointData("id", std::move(ids));
    }
    std::vector<NDArray> tags;
    std::vector<std::int32_t> tag_vals(sphere.Cells(0).NumCells());
    for (std::size_t c = 0; c < tag_vals.size(); ++c)
        tag_vals[c] = static_cast<std::int32_t>(c % 3);
    tags.push_back(mt::int_data_array(tag_vals));
    sphere.AddCellData("tag", std::move(tags));
    sphere.AddFieldData("meta", mt::data_array({42.0}));

    DecimateOptions opts = ratio_opts(0.5);
    opts.mPreserveFeatures = false;  // the coarse sphere pins itself otherwise
    const DecimateResult res = decimate(sphere, opts);
    const Mesh& out = res.mMesh;

    // Blended float data stays inside the input hull; integer point data keeps
    // survivor rows (values from the input set); cell tags keep their rows.
    ASSERT_TRUE(out.HasPointData("T"));
    const NDArray& t_out = out.PointData("T");
    EXPECT_TRUE(meshioplusplus::detail::is_float_dtype(t_out.Dtype()));
    for (std::size_t i = 0; i < t_out.Size(); ++i) {
        const double v = meshioplusplus::detail::read_double(t_out, i);
        EXPECT_GE(v, -1.0 - 1e-12);
        EXPECT_LE(v, 1.0 + 1e-12);
    }
    ASSERT_TRUE(out.HasPointData("id"));
    const NDArray& id_out = out.PointData("id");
    EXPECT_FALSE(meshioplusplus::detail::is_float_dtype(id_out.Dtype()));
    for (std::size_t i = 0; i < id_out.Size(); ++i) {
        const std::int64_t v = meshioplusplus::detail::read_int(id_out, i);
        EXPECT_GE(v, 0);
        EXPECT_LT(v, 7);
    }
    ASSERT_TRUE(out.HasCellData("tag"));
    EXPECT_EQ(meshioplusplus::detail::rows(out.CellData("tag", 0)), out.Cells(0).NumCells());
    EXPECT_FALSE(out.HasCellData("convert:parent_cell"));
    ASSERT_TRUE(out.HasFieldData("meta"));
}

TEST(Decimate, PlacementModesAllProduceValidMeshes) {
    const Mesh sphere = sphere_mesh(2);
    for (const DecimatePlacement placement :
         {DecimatePlacement::Optimal, DecimatePlacement::Midpoint, DecimatePlacement::Endpoint}) {
        DecimateOptions opts = ratio_opts(0.5);
        opts.mPlacement = placement;
        opts.mPreserveFeatures = false;
        const DecimateResult res = decimate(sphere, opts);
        EXPECT_LT(total_cells(res.mMesh), total_cells(sphere));
        expect_manifold(res.mMesh);
    }
}

TEST(Decimate, EndpointPlacementKeepsInputCoordinates) {
    Mesh sphere = sphere_mesh(2);
    DecimateOptions opts = ratio_opts(0.5);
    opts.mPlacement = DecimatePlacement::Endpoint;
    opts.mPreserveFeatures = false;
    const DecimateResult res = decimate(sphere, opts);
    for (std::size_t i = 0; i < res.mMesh.NumPoints(); ++i) {
        const std::array<double, 3> p = point_at(res.mMesh, i);
        EXPECT_TRUE(has_point(sphere, p[0], p[1], p[2]));
    }
}

TEST(Decimate, ResultIsStableAcrossRepeatedRuns) {
    const Mesh sphere = sphere_mesh(2);
    DecimateOptions opts = ratio_opts(0.4);
    opts.mPreserveFeatures = false;
    const DecimateResult first = decimate(sphere, opts);
    for (int run = 0; run < 3; ++run) {
        const DecimateResult again = decimate(sphere, opts);
        ASSERT_EQ(again.mMesh.NumPoints(), first.mMesh.NumPoints());
        for (std::size_t i = 0; i < first.mMesh.Points().Size(); ++i)
            ASSERT_EQ(meshioplusplus::detail::read_double(again.mMesh.Points(), i),
                      meshioplusplus::detail::read_double(first.mMesh.Points(), i));
        ASSERT_EQ(mt::cell_rows(again.mMesh), mt::cell_rows(first.mMesh));
        ASSERT_EQ(again.mFacesRemoved, first.mFacesRemoved);
        ASSERT_EQ(again.mCollapsesRejected, first.mCollapsesRejected);
        ASSERT_EQ(again.mMaxErrorApplied, first.mMaxErrorApplied);
    }
}

TEST(Decimate, RejectsOutOfScopeInputsByName) {
    const auto message_of = [](const Mesh& rMesh) {
        try {
            decimate(rMesh, target_opts(1));
        } catch (const std::invalid_argument& e) {
            return std::string(e.what());
        }
        return std::string();
    };
    EXPECT_NE(message_of(mt::tet_mesh()).find("extract_surface"), std::string::npos);
    EXPECT_NE(message_of(mt::hex_mesh()).find("extract_surface"), std::string::npos);
    EXPECT_NE(message_of(mt::triangle6_mesh()).find("linearize"), std::string::npos);
    EXPECT_NE(message_of(mt::line_mesh()).find("surface"), std::string::npos);

    Mesh mixed = mt::tri_mesh();
    mixed.AddCellBlock("line", mt::conn_from({{0, 1}}));
    EXPECT_NE(message_of(mixed).find("line"), std::string::npos);
}

TEST(Decimate, ValidatesTheStoppingCriteria) {
    const Mesh grid = grid_mesh(2);
    EXPECT_THROW(decimate(grid, DecimateOptions{}), std::invalid_argument);
    {
        DecimateOptions two;
        two.mTargetRatio = 0.5;
        two.mTargetFaces = 4;
        EXPECT_THROW(decimate(grid, two), std::invalid_argument);
    }
    {
        DecimateOptions bad_ratio;
        bad_ratio.mTargetRatio = 1.5;
        EXPECT_THROW(decimate(grid, bad_ratio), std::invalid_argument);
    }
    EXPECT_THROW(meshioplusplus::decimate_placement_from_name("nearest"), std::invalid_argument);
    EXPECT_EQ(meshioplusplus::decimate_placement_from_name("midpoint"),
              DecimatePlacement::Midpoint);
}

}  // namespace
