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
#include <array>
#include <cmath>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/tessellate.hpp"
#include "meshioplusplus/region.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

using Vec3 = std::array<double, 3>;

// Reference node coordinates of a curved type, in the layout tessellate reads:
// corners, edge midpoints in mid-node slot order, face centres, body centre.
std::vector<Vec3> reference_nodes(const std::string& rType) {
    std::vector<Vec3> corners;
    CellType base;
    bool faces = false, body = false;
    if (rType == "triangle6") {
        corners = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
        base = CellType::Triangle;
    } else if (rType == "quad8" || rType == "quad9") {
        corners = {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}};
        base = CellType::Quad;
        faces = rType == "quad9";
    } else if (rType == "tetra10") {
        corners = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        base = CellType::Tetra;
    } else {
        corners = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                   {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
        base = CellType::Hexahedron;
        faces = body = true;
    }
    std::vector<Vec3> out = corners;
    for (const auto& e : detail::cell_refine_edges(base)) {
        Vec3 m{};
        for (int a = 0; a < 3; ++a)
            m[a] = 0.5 * (corners[e[0]][a] + corners[e[1]][a]);
        out.push_back(m);
    }
    if (faces)
        for (const auto& f : detail::cell_refine_quad_faces(base)) {
            Vec3 m{};
            for (int a = 0; a < 3; ++a)
                m[a] = 0.25 *
                       (corners[f[0]][a] + corners[f[1]][a] + corners[f[2]][a] + corners[f[3]][a]);
            out.push_back(m);
        }
    if (body)
        out.push_back({0, 0, 0});
    return out;
}

// x -> A x + b with a generic (asymmetric) A.
Vec3 affine(const Vec3& rX) {
    return {1.7 * rX[0] + 0.3 * rX[1] - 0.2 * rX[2] + 0.5,
            0.1 * rX[0] + 2.1 * rX[1] + 0.4 * rX[2] - 1.0,
            -0.3 * rX[0] + 0.2 * rX[1] + 1.3 * rX[2] + 2.0};
}

double affine_det() {
    const double a[3][3] = {{1.7, 0.3, -0.2}, {0.1, 2.1, 0.4}, {-0.3, 0.2, 1.3}};
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
           a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
           a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

// One cell of `rType` under the affine map, plus a linear field f = x + 2y + 3z.
Mesh affine_cell(const std::string& rType) {
    std::vector<std::vector<double>> pts;
    std::vector<std::int64_t> conn;
    std::vector<double> f;
    for (const Vec3& x : reference_nodes(rType)) {
        const Vec3 y = affine(x);
        pts.push_back({y[0], y[1], y[2]});
        f.push_back(y[0] + 2.0 * y[1] + 3.0 * y[2]);
        conn.push_back(static_cast<std::int64_t>(conn.size()));
    }
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddCellBlock(rType, mt::conn_from({conn}));
    NDArray fa = NDArray::Uninit(DType::Float64, {f.size()});
    std::copy(f.begin(), f.end(), fa.As<double>());
    m.AddPointData("f", std::move(fa));
    return m;
}

const double* pts_of(const Mesh& rMesh) {
    return rMesh.Points().As<double>();
}

double tet_volume(const double* pP, const std::int64_t* pC) {
    const double* a = pP + 3 * pC[0];
    const double* b = pP + 3 * pC[1];
    const double* c = pP + 3 * pC[2];
    const double* d = pP + 3 * pC[3];
    const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    const double w[3] = {d[0] - a[0], d[1] - a[1], d[2] - a[2]};
    return (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) +
            u[2] * (v[0] * w[1] - v[1] * w[0])) /
           6.0;
}

}  // namespace

TEST(Tessellate, EveryCurvedTypeReproducesAnAffineMap) {
    for (const std::string type : {"triangle6", "quad8", "quad9", "tetra10", "hexahedron27"}) {
        for (std::int64_t levels = 0; levels <= 3; ++levels) {
            const Mesh m = affine_cell(type);
            TessellateOptions o;
            o.mLevels = levels;
            const TessellateResult r = tessellate(m, o);
            const std::size_t np = r.mMesh.NumPoints();
            ASSERT_EQ(r.mNumPoints, static_cast<std::int64_t>(np));
            const double* p = pts_of(r.mMesh);
            const double* src = pts_of(m);
            const NDArray& f = r.mMesh.PointData("f");
            for (std::size_t i = 0; i < np; ++i) {
                // Partition of unity, and the stencil reproduces the point.
                double wsum = 0.0, x[3] = {0, 0, 0};
                for (std::size_t k = 0; k < kTessellateStencilWidth; ++k) {
                    const std::int64_t s = r.mStencil[i * kTessellateStencilWidth + k];
                    const double w = r.mWeights[i * kTessellateStencilWidth + k];
                    if (s < 0) {
                        EXPECT_EQ(w, 0.0);
                        continue;
                    }
                    wsum += w;
                    for (int a = 0; a < 3; ++a)
                        x[a] += w * src[3 * s + a];
                }
                EXPECT_NEAR(wsum, 1.0, 1e-12) << type << " level " << levels;
                for (int a = 0; a < 3; ++a)
                    EXPECT_NEAR(x[a], p[3 * i + a], 1e-12);
                // A linear field is interpolated exactly.
                EXPECT_NEAR(f.As<double>()[i], p[3 * i] + 2.0 * p[3 * i + 1] + 3.0 * p[3 * i + 2],
                            1e-11);
            }
            EXPECT_EQ(r.mNumFacetsWithBadCount, 0);
            EXPECT_EQ(r.mNumCurvedSourceCells, 1);
        }
    }
}

TEST(Tessellate, ChildCountsFollowTheLattice) {
    auto cells = [](const std::string& rType, std::int64_t Levels) {
        TessellateOptions o;
        o.mLevels = Levels;
        return tessellate(affine_cell(rType), o).mNumCells;
    };
    EXPECT_EQ(cells("triangle6", 0), 1);
    EXPECT_EQ(cells("triangle6", 2), 16);
    EXPECT_EQ(cells("tetra10", 2), 64);
    EXPECT_EQ(cells("quad9", 0), 2);  // max(levels, 1) divisions
    EXPECT_EQ(cells("quad8", 3), 18);
    EXPECT_EQ(cells("hexahedron27", 2), 96);
}

TEST(Tessellate, HexahedronVolumeIsPreserved) {
    TessellateOptions o;
    o.mLevels = 3;
    const TessellateResult r = tessellate(affine_cell("hexahedron27"), o);
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    const auto cb = r.mMesh.Cells(0);
    ASSERT_EQ(std::string(cb.Type()), "tetra");
    double vol = 0.0;
    for (std::size_t c = 0; c < cb.NumCells(); ++c)
        vol += std::fabs(tet_volume(pts_of(r.mMesh), cb.Conn().As<std::int64_t>() + 4 * c));
    EXPECT_NEAR(vol, 8.0 * std::fabs(affine_det()), 1e-10);
}

TEST(Tessellate, ASharedHexahedronFaceConformsInEveryOrientation) {
    const std::vector<Vec3> ref = reference_nodes("hexahedron27");
    std::vector<std::array<int, 8>> rotations;
    const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (const auto& perm : perms)
        for (int s = 0; s < 8; ++s) {
            double r[3][3] = {};
            for (int i = 0; i < 3; ++i)
                r[i][perm[i]] = (s >> i) & 1 ? -1.0 : 1.0;
            const double det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1]) -
                               r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0]) +
                               r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
            if (det < 0)
                continue;
            std::array<int, 8> order{};
            for (int c = 0; c < 8; ++c) {
                Vec3 moved{};
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        moved[i] += r[i][j] * ref[c][j];
                for (int k = 0; k < 8; ++k)
                    if (ref[k] == moved)
                        order[c] = k;
            }
            rotations.push_back(order);
        }
    ASSERT_EQ(rotations.size(), 24u);
    for (const auto& order : rotations) {
        // Cell A on [-1, 1]^3, cell B its copy shifted by 2 in x with its
        // corners renumbered by the rotation; shared nodes welded by position.
        std::vector<Vec3> nodes;
        std::vector<std::vector<std::int64_t>> conn(2);
        auto node_id = [&](const Vec3& rP) {
            for (std::size_t i = 0; i < nodes.size(); ++i)
                if (nodes[i] == rP)
                    return static_cast<std::int64_t>(i);
            nodes.push_back(rP);
            return static_cast<std::int64_t>(nodes.size() - 1);
        };
        for (const Vec3& x : ref)
            conn[0].push_back(node_id(x));
        // B's corner c sits where A's corner order[c] sits, shifted; its other
        // nodes follow from its corners exactly as the reference layout does.
        std::vector<Vec3> bc(8);
        for (int c = 0; c < 8; ++c)
            bc[c] = {ref[order[c]][0] + 2.0, ref[order[c]][1], ref[order[c]][2]};
        for (const Vec3& x : ref) {
            // Trilinear image of reference node x through B's corners.
            Vec3 y{};
            for (int c = 0; c < 8; ++c) {
                double w = 1.0;
                for (int a = 0; a < 3; ++a)
                    w *= ref[c][a] > 0 ? (1.0 + x[a]) * 0.5 : (1.0 - x[a]) * 0.5;
                for (int a = 0; a < 3; ++a)
                    y[a] += w * bc[c][a];
            }
            for (double& v : y)
                v = std::round(v * 1e9) / 1e9;
            conn[1].push_back(node_id(y));
        }
        std::vector<std::vector<double>> pts;
        for (const Vec3& p : nodes)
            pts.push_back({p[0], p[1], p[2]});
        Mesh m;
        m.AssignPoints(mt::points_from(pts));
        m.AddCellBlock("hexahedron27", mt::conn_from(conn));
        const TessellateResult r = tessellate(m);
        EXPECT_EQ(r.mNumFacetsWithBadCount, 0);
        std::set<std::array<double, 3>> unique;
        const double* p = pts_of(r.mMesh);
        for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i)
            unique.insert({std::round(p[3 * i] * 1e9), std::round(p[3 * i + 1] * 1e9),
                           std::round(p[3 * i + 2] * 1e9)});
        EXPECT_EQ(unique.size(), r.mMesh.NumPoints());
    }
}

TEST(Tessellate, ALinearMeshPassesThroughUnchanged) {
    const Mesh m = mt::tri_quad_mesh();
    const TessellateResult r = tessellate(m);
    ASSERT_EQ(r.mMesh.NumPoints(), m.NumPoints());
    ASSERT_EQ(r.mMesh.NumCellBlocks(), m.NumCellBlocks());
    const NDArray& sp = r.mMesh.PointData(kTessellateSourcePointName);
    for (std::size_t i = 0; i < m.NumPoints(); ++i) {
        EXPECT_EQ(detail::read_int(sp, i), static_cast<std::int64_t>(i));
        for (std::size_t a = 0; a < m.PointDim(); ++a)
            EXPECT_EQ(pts_of(r.mMesh)[i * m.PointDim() + a],
                      detail::read_double(m.Points(), i * m.PointDim() + a));
    }
    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
        EXPECT_EQ(std::string(r.mMesh.Cells(b).Type()), std::string(m.Cells(b).Type()));
        const NDArray& x = r.mMesh.Cells(b).Conn();
        const NDArray& y = m.Cells(b).Conn();
        for (std::size_t i = 0; i < y.Size(); ++i)
            EXPECT_EQ(detail::read_int(x, i), detail::read_int(y, i));
    }
    EXPECT_EQ(r.mNumCurvedSourceCells, 0);
}

TEST(Tessellate, CurvedFalseIsANoOp) {
    TessellateOptions o;
    o.mCurved = false;
    const TessellateResult r = tessellate(affine_cell("tetra10"), o);
    EXPECT_EQ(std::string(r.mMesh.Cells(0).Type()), "tetra10");
    EXPECT_EQ(r.mNumPassThroughSourceCells, 1);
}

TEST(Tessellate, RegionsAndOptionalArrays) {
    Mesh m = affine_cell("quad9");
    m.AddCellBlock("vertex", mt::conn_from({{4}}));
    NDArray pe = NDArray::Uninit(DType::Int64, {2});
    pe.As<std::int64_t>()[0] = 0;
    pe.As<std::int64_t>()[1] = 4;
    m.AddRegion(Region("pts", RegionKind::Point, 0, 1, std::move(pe)));
    NDArray ce = NDArray::Uninit(DType::Int64, {2});
    ce.As<std::int64_t>()[0] = 0;  // the quad9: tessellated, so dropped
    ce.As<std::int64_t>()[1] = 1;  // the vertex: passes through
    m.AddRegion(Region("cells", RegionKind::Cell, -1, 2, std::move(ce)));
    NDArray se = NDArray::Uninit(DType::Int64, {1, 2});
    se.As<std::int64_t>()[0] = 0;
    se.As<std::int64_t>()[1] = 0;
    m.AddRegion(Region("side", RegionKind::Side, -1, 3, std::move(se)));

    TessellateOptions o;
    o.mRecordStencil = true;
    const TessellateResult r = tessellate(m, o);
    EXPECT_TRUE(r.mMesh.HasPointData(kTessellateStencilName));
    EXPECT_TRUE(r.mMesh.HasPointData(kTessellateWeightsName));
    EXPECT_FALSE(r.mMesh.HasRegion("side"));
    const std::size_t pi = r.mMesh.FindRegion("pts", RegionKind::Point);
    ASSERT_NE(pi, Mesh::npos);
    const Region& pr = r.mMesh.Region(pi);
    ASSERT_EQ(pr.NumEntries(), 2u);
    const NDArray& sp = r.mMesh.PointData(kTessellateSourcePointName);
    std::set<std::int64_t> src;
    for (std::size_t i = 0; i < 2; ++i)
        src.insert(
            detail::read_int(sp, static_cast<std::size_t>(detail::read_int(pr.mEntries, i))));
    EXPECT_EQ(src, (std::set<std::int64_t>{0, 4}));
    const std::size_t ci = r.mMesh.FindRegion("cells", RegionKind::Cell);
    ASSERT_NE(ci, Mesh::npos);
    ASSERT_EQ(r.mMesh.Region(ci).NumEntries(), 1u);
    // The vertex block comes first in the output, so it is output cell 0.
    EXPECT_EQ(detail::read_int(r.mMesh.Region(ci).mEntries, 0), 0);

    TessellateOptions bare;
    bare.mFields = false;
    const TessellateResult b = tessellate(m, bare);
    EXPECT_FALSE(b.mMesh.HasPointData("f"));
    EXPECT_TRUE(b.mMesh.HasPointData(kTessellateSourcePointName));
}

TEST(Tessellate, AnUnusedPointIsDropped) {
    Mesh m = affine_cell("triangle6");
    std::vector<std::vector<double>> pts;
    const double* p = pts_of(m);
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        pts.push_back({p[3 * i], p[3 * i + 1], p[3 * i + 2]});
    pts.push_back({9.0, 9.0, 9.0});
    Mesh n;
    n.AssignPoints(mt::points_from(pts));
    n.AddCellBlock("triangle6", mt::conn_from({{0, 1, 2, 3, 4, 5}}));
    TessellateOptions o;
    o.mLevels = 0;
    const TessellateResult r = tessellate(n, o);
    EXPECT_EQ(r.mMesh.NumPoints(), 3u);
    EXPECT_EQ(r.mNumSourcePoints, 7);
}

TEST(Tessellate, BadArgumentsAreRefused) {
    TessellateOptions o;
    o.mLevels = -1;
    EXPECT_THROW(tessellate(affine_cell("quad9"), o), std::invalid_argument);
    o.mLevels = 21;
    EXPECT_THROW(tessellate(affine_cell("quad9"), o), std::invalid_argument);
    Mesh m = affine_cell("quad9");
    m.AddPointData("short", NDArray::Uninit(DType::Float64, {2}));
    EXPECT_THROW(tessellate(m), std::invalid_argument);
}
