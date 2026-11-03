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
// cut_level_set and cut_quadrature (doc/cut.md): exact measures of planar cuts
// on every cell family, several level sets, convergence on a sphere, quadratic
// level sets, quadrature, the data carried to the pieces, determinism and the
// refusals. Without CutCells only the refusal runs.

// System includes
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/cut.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/region.hpp"

namespace mio = meshioplusplus;

namespace {

int tcut_dim(const std::string& rType) {
    if (rType == "vertex")
        return 0;
    if (rType == "line" || rType == "VTK_LAGRANGE_CURVE")
        return 1;
    if (rType == "triangle" || rType == "quad" ||
        rType.find("QUADRILATERAL") != std::string::npos ||
        rType.find("TRIANGLE") != std::string::npos)
        return 2;
    return 3;
}

/// The summed measure of the cells of dimension `Dim`.
double tcut_measure(const mio::Mesh& rMesh, int Dim) {
    double s = 0.0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cells = rMesh.Cells(b);
        if (tcut_dim(cells.Type()) != Dim)
            continue;
        for (std::size_t c = 0; c < cells.NumCells(); ++c)
            s += mio::detail::cell_measure(rMesh.Points(), rMesh.PointDim(), cells, c);
    }
    return s;
}

/// The unit cube as six pyramids around its centre.
mio::Mesh tcut_pyramids() {
    mio::Mesh m;
    const double xyz[9][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},      {0, 0, 1},
                              {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {0.5, 0.5, 0.5}};
    mio::NDArray p = mio::NDArray::Uninit(mio::DType::Float64, {std::size_t{9}, std::size_t{3}});
    for (int i = 0; i < 9; ++i)
        for (int k = 0; k < 3; ++k)
            p.As<double>()[3 * i + k] = xyz[i][k];
    m.AssignPoints(std::move(p));
    const std::int64_t faces[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                      {2, 3, 7, 6}, {1, 2, 6, 5}, {0, 4, 7, 3}};
    mio::NDArray c = mio::NDArray::Uninit(mio::DType::Int64, {std::size_t{6}, std::size_t{5}});
    for (int f = 0; f < 6; ++f) {
        for (int k = 0; k < 4; ++k)
            c.As<std::int64_t>()[5 * f + k] = faces[f][k];
        c.As<std::int64_t>()[5 * f + 4] = 8;
    }
    m.AddCellBlock("pyramid", std::move(c));
    return m;
}

/// Adds point data `rName` = f(x, y, z) at every point.
template <typename F>
void tcut_field(mio::Mesh& rMesh, const std::string& rName, F Fn) {
    const std::size_t n = rMesh.NumPoints(), d = rMesh.PointDim();
    mio::NDArray a = mio::NDArray::Uninit(mio::DType::Float64, {n});
    for (std::size_t i = 0; i < n; ++i) {
        double x[3] = {0, 0, 0};
        for (std::size_t k = 0; k < d; ++k)
            x[k] = mio::detail::read_double(rMesh.Points(), i * d + k);
        a.As<double>()[i] = Fn(x[0], x[1], x[2]);
    }
    rMesh.AddPointData(rName, std::move(a));
}

mio::CutOptions tcut_plane(double X, const std::string& rSelection) {
    mio::CutOptions o;
    mio::CutPlane pl;
    pl.mOrigin = {X, 0.0, 0.0};
    pl.mNormal = {1.0, 0.0, 0.0};
    o.mPlanes.push_back(pl);
    o.mSelection = rSelection;
    return o;
}

#define TCUT_NEED_CUTCELLS()      \
    if (!mio::cut_has_cutcells()) \
    GTEST_SKIP() << "this build has no CutCells"

}  // namespace

TEST(Cut, WithoutCutCellsBothOperationsRefuseByName) {
    if (mio::cut_has_cutcells())
        GTEST_SKIP() << "CutCells is compiled in";
    const mio::Mesh box = mio::box({2, 2, 2});
    try {
        (void)mio::cut_level_set(box, tcut_plane(0.3, ""));
        FAIL() << "expected Unsupported";
    } catch (const mio::Unsupported& rE) {
        EXPECT_NE(std::string(rE.what()).find("MESHIOPLUSPLUS_WITH_CUTCELLS"), std::string::npos);
    }
    EXPECT_THROW((void)mio::cut_quadrature(box, tcut_plane(0.3, "")), mio::Unsupported);
}

TEST(Cut, APlaneCutGivesExactMeasuresOnEveryFamily) {
    TCUT_NEED_CUTCELLS();
    std::map<std::string, mio::Mesh> meshes;
    meshes["hexahedron"] = mio::box({3, 3, 3});
    meshes["tetra"] = mio::box({3, 3, 3}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "tetra");
    meshes["wedge"] = mio::box({3, 3, 3}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "wedge");
    meshes["pyramid"] = tcut_pyramids();
    meshes["quad"] = mio::rectangle({4, 4});
    meshes["triangle"] = mio::rectangle({4, 4}, {1.0, 1.0}, {0.0, 0.0, 0.0}, "triangle");
    for (const auto& [name, mesh] : meshes) {
        SCOPED_TRACE(name);
        const int dim = name == "quad" || name == "triangle" ? 2 : 3;
        for (const bool include_uncut : {true, false}) {
            auto o = tcut_plane(0.3, "plane0 < 0");
            o.mIncludeUncut = include_uncut;
            const auto in = mio::cut_level_set(mesh, o);
            if (include_uncut) {
                EXPECT_NEAR(in.mMeasure, 0.3, 1e-12);
                EXPECT_NEAR(tcut_measure(in.mMesh, dim), 0.3, 1e-12);
            }
            EXPECT_EQ(in.mDimension, dim);
            EXPECT_EQ(in.mBackend, "lut");
            EXPECT_GT(in.mNumCutCells, 0);
        }
        const auto out = mio::cut_level_set(mesh, tcut_plane(0.3, "plane0 > 0"));
        EXPECT_NEAR(tcut_measure(out.mMesh, dim), 0.7, 1e-12);
        const auto face = mio::cut_level_set(mesh, tcut_plane(0.3, "plane0 = 0"));
        EXPECT_EQ(face.mDimension, dim - 1);
        EXPECT_NEAR(face.mMeasure, 1.0, 1e-12);
        EXPECT_NEAR(tcut_measure(face.mMesh, dim - 1), 1.0, 1e-12);
        // triangulate splits every piece and every whole cell into simplices.
        auto t = tcut_plane(0.3, "plane0 > 0");
        t.mTriangulate = true;
        const auto simplices = mio::cut_level_set(mesh, t);
        for (std::size_t b = 0; b < simplices.mMesh.NumCellBlocks(); ++b)
            EXPECT_EQ(simplices.mMesh.Cells(b).Type(), dim == 3 ? "tetra" : "triangle");
        EXPECT_NEAR(tcut_measure(simplices.mMesh, dim), 0.7, 1e-12);
        // The pieces are welded: no point appears twice.
        std::size_t refs = 0;
        for (std::size_t b = 0; b < out.mMesh.NumCellBlocks(); ++b)
            refs += out.mMesh.Cells(b).NumCells() * out.mMesh.Cells(b).NodesPerCell();
        EXPECT_LT(out.mMesh.NumPoints(), refs);
    }
}

TEST(Cut, SeveralLevelSetsCombine) {
    TCUT_NEED_CUTCELLS();
    const mio::Mesh box = mio::box({2, 2, 2});
    mio::CutOptions o;
    mio::CutPlane a, b;
    a.mOrigin = {0.5, 0.0, 0.0};
    a.mNormal = {1.0, 0.0, 0.0};
    b.mOrigin = {0.0, 0.5, 0.0};
    b.mNormal = {0.0, 2.0, 0.0};
    o.mPlanes = {a, b};
    o.mSelection = "plane0 < 0 and plane1 < 0";
    EXPECT_NEAR(mio::cut_level_set(box, o).mMeasure, 0.25, 1e-12);
    o.mSelection = "plane0 < 0 or plane1 < 0";
    EXPECT_NEAR(mio::cut_level_set(box, o).mMeasure, 0.75, 1e-12);
    o.mSelection = "plane0 = 0 and plane1 < 0";
    const auto face = mio::cut_level_set(box, o);
    EXPECT_EQ(face.mDimension, 2);
    EXPECT_NEAR(face.mMeasure, 0.5, 1e-12);
}

TEST(Cut, ASphereConvergesAndAQuadraticLevelSetIsMoreAccurate) {
    TCUT_NEED_CUTCELLS();
    const double r = 0.3, exact = 4.0 / 3.0 * M_PI * r * r * r;
    const auto ball = [&](double x, double y, double z) {
        return std::sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) + (z - 0.5) * (z - 0.5)) - r;
    };
    double previous = 1.0;
    for (const std::int64_t n : {4, 8, 16}) {
        mio::Mesh m = mio::box({n, n, n}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "tetra");
        tcut_field(m, "ball", ball);
        mio::CutOptions o;
        o.mLevelSets = {"ball"};
        const auto res = mio::cut_level_set(m, o);
        const double err = std::fabs(res.mMeasure - exact);
        EXPECT_NEAR(tcut_measure(res.mMesh, 3), res.mMeasure, 1e-12);
        EXPECT_LT(err, previous / 2.5) << n;
        previous = err;
    }
    // The same mesh size with quadratic tetrahedra: the level set is quadratic.
    mio::Mesh lin = mio::box({4, 4, 4}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "tetra");
    mio::Mesh quad = mio::box({4, 4, 4}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "tetra", 2);
    ASSERT_EQ(quad.Cells(0).Type(), "tetra10");
    tcut_field(lin, "ball", ball);
    tcut_field(quad, "ball", ball);
    mio::CutOptions o;
    o.mLevelSets = {"ball"};
    const double e1 = std::fabs(mio::cut_quadrature(lin, o).mSumWeights - exact);
    o.mOrder = 4;
    const auto q2 = mio::cut_quadrature(quad, o);
    EXPECT_EQ(q2.mBackend, "quadrays");
    const double e2 = std::fabs(q2.mSumWeights - exact);
    EXPECT_LT(e2, e1 / 5.0) << e1 << " " << e2;
    // Curved pieces: Lagrange cells from the quadrays backend.
    o.mDegree = 2;
    const auto curved = mio::cut_level_set(quad, o);
    bool lagrange = false;
    for (std::size_t b = 0; b < curved.mMesh.NumCellBlocks(); ++b)
        lagrange = lagrange || curved.mMesh.Cells(b).Type().rfind("VTK_LAGRANGE_", 0) == 0;
    EXPECT_TRUE(lagrange);
    EXPECT_NEAR(curved.mMeasure, q2.mSumWeights, 1e-3 * exact);
}

TEST(Cut, QuadratureIntegratesOnThePartAndItsInterface) {
    TCUT_NEED_CUTCELLS();
    const mio::Mesh box = mio::box({2, 2, 2});
    auto o = tcut_plane(0.3, "plane0 < 0");
    o.mOrder = 2;
    const auto q = mio::cut_quadrature(box, o);
    EXPECT_NEAR(q.mSumWeights, 0.3, 1e-12);
    const auto& m = q.mMesh;
    ASSERT_TRUE(m.HasPointData("quadrature:weight"));
    ASSERT_TRUE(m.HasPointData("quadrature:parent_cell"));
    ASSERT_EQ(m.Cells(0).Type(), "vertex");
    double ix = 0.0, ixx = 0.0;
    for (std::size_t i = 0; i < m.NumPoints(); ++i) {
        const double w = mio::detail::read_double(m.PointData("quadrature:weight"), i);
        const double x = mio::detail::read_double(m.Points(), 3 * i);
        EXPECT_LE(x, 0.3 + 1e-12);
        ix += w * x;
        ixx += w * x * x;
    }
    EXPECT_NEAR(ix, 0.3 * 0.3 / 2.0, 1e-12);
    EXPECT_NEAR(ixx, 0.3 * 0.3 * 0.3 / 3.0, 1e-12);
    const auto f = mio::cut_quadrature(box, tcut_plane(0.3, "plane0 = 0"));
    EXPECT_EQ(f.mDimension, 2);
    EXPECT_NEAR(f.mSumWeights, 1.0, 1e-12);
    for (std::size_t i = 0; i < f.mMesh.NumPoints(); ++i)
        EXPECT_NEAR(mio::detail::read_double(f.mMesh.Points(), 3 * i), 0.3, 1e-12);
}

TEST(Cut, DataAndCellRegionsFollowTheParents) {
    TCUT_NEED_CUTCELLS();
    for (const std::string type : {"hexahedron", "wedge", "tetra"}) {
        SCOPED_TRACE(type);
        mio::Mesh m = mio::box({2, 2, 2}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, type);
        tcut_field(m, "f", [](double x, double y, double z) { return 2 * x + 3 * y - z + 1; });
        const std::size_t nc = m.Cells(0).NumCells();
        mio::NDArray id = mio::NDArray::Uninit(mio::DType::Int64, {nc});
        std::vector<std::int64_t> left;
        for (std::size_t c = 0; c < nc; ++c) {
            id.As<std::int64_t>()[c] = static_cast<std::int64_t>(100 + c);
            if (c % 3 == 0)
                left.push_back(static_cast<std::int64_t>(c));
        }
        m.AddCellData("id", {std::move(id)});
        mio::NDArray entries = mio::NDArray::Uninit(mio::DType::Int64, {left.size()});
        std::copy(left.begin(), left.end(), entries.As<std::int64_t>());
        m.AddRegion(mio::Region("left", mio::RegionKind::Cell, std::move(entries)));
        const auto res = mio::cut_level_set(m, tcut_plane(0.3, "plane0 < 0"));
        const mio::Mesh& out = res.mMesh;
        ASSERT_TRUE(out.HasPointData("f"));
        for (std::size_t i = 0; i < out.NumPoints(); ++i) {
            const double x = mio::detail::read_double(out.Points(), 3 * i);
            const double y = mio::detail::read_double(out.Points(), 3 * i + 1);
            const double z = mio::detail::read_double(out.Points(), 3 * i + 2);
            EXPECT_NEAR(mio::detail::read_double(out.PointData("f"), i), 2 * x + 3 * y - z + 1,
                        1e-10);
        }
        std::int64_t global = 0;
        ASSERT_TRUE(out.HasRegion("left"));
        const mio::Region& reg = out.Region(out.FindRegion("left", mio::RegionKind::Cell));
        std::vector<std::int64_t> in_left;
        for (std::size_t b = 0; b < out.NumCellBlocks(); ++b) {
            const std::size_t n = out.Cells(b).NumCells();
            for (std::size_t c = 0; c < n; ++c, ++global) {
                const std::int64_t parent =
                    mio::detail::read_int(out.CellData("cut:parent_cell", b), c);
                EXPECT_EQ(mio::detail::read_int(out.CellData("id", b), c), 100 + parent);
                if (parent % 3 == 0)
                    in_left.push_back(global);
            }
        }
        ASSERT_EQ(reg.NumEntries(), in_left.size());
        for (std::size_t i = 0; i < in_left.size(); ++i)
            EXPECT_EQ(reg.Entries()[i], in_left[i]);
    }
}

TEST(Cut, IsDeterministicAndLeavesLowerDimensionsOut) {
    TCUT_NEED_CUTCELLS();
    mio::Mesh m = mio::box({3, 3, 3}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, "tetra");
    tcut_field(m, "ball", [](double x, double y, double z) { return x * x + y * y + z * z - 0.5; });
    mio::NDArray tri = mio::NDArray::Uninit(mio::DType::Int64, {std::size_t{1}, std::size_t{3}});
    tri.As<std::int64_t>()[0] = 0;
    tri.As<std::int64_t>()[1] = 1;
    tri.As<std::int64_t>()[2] = 4;
    m.AddCellBlock("triangle", std::move(tri));
    mio::CutOptions o;
    o.mLevelSets = {"ball"};
    const auto a = mio::cut_level_set(m, o);
    const auto b = mio::cut_level_set(m, o);
    ASSERT_EQ(a.mMesh.NumPoints(), b.mMesh.NumPoints());
    for (std::size_t i = 0; i < a.mMesh.Points().Size(); ++i)
        ASSERT_EQ(a.mMesh.Points().As<double>()[i], b.mMesh.Points().As<double>()[i]);
    bool warned = false;
    for (const auto& r_w : a.mWarnings)
        warned = warned || r_w.find("'triangle'") != std::string::npos;
    EXPECT_TRUE(warned);
    for (std::size_t k = 0; k < a.mMesh.NumCellBlocks(); ++k)
        EXPECT_NE(a.mMesh.Cells(k).Type(), "triangle");
}

TEST(Cut, RefusesByName) {
    TCUT_NEED_CUTCELLS();
    mio::Mesh m = mio::box({2, 2, 2});
    tcut_field(m, "phi", [](double x, double, double) { return x - 0.4; });
    mio::CutOptions o;
    EXPECT_THROW((void)mio::cut_level_set(m, o), std::invalid_argument);  // no level set
    o.mLevelSets = {"nope"};
    EXPECT_THROW((void)mio::cut_level_set(m, o), std::invalid_argument);
    o.mLevelSets = {"phi"};
    o.mSelection = "psi < 0";
    try {
        (void)mio::cut_level_set(m, o);
        FAIL() << "expected invalid_argument";
    } catch (const std::invalid_argument& rE) {
        EXPECT_NE(std::string(rE.what()).find("psi"), std::string::npos) << rE.what();
    }
    o.mSelection.clear();
    o.mBackend = "marching";
    EXPECT_THROW((void)mio::cut_level_set(m, o), std::invalid_argument);
    o.mBackend = "auto";
    o.mOrder = -1;
    EXPECT_THROW((void)mio::cut_quadrature(m, o), std::invalid_argument);
    o.mOrder = 2;
    EXPECT_NO_THROW((void)mio::cut_quadrature(m, o));
    mio::Mesh poly;
    mio::NDArray p = mio::NDArray::Uninit(mio::DType::Float64, {std::size_t{4}, std::size_t{3}});
    const double xyz[12] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    std::copy(xyz, xyz + 12, p.As<double>());
    poly.AssignPoints(std::move(p));
    poly.AddPolyhedronBlock("polyhedron4", {{{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}}});
    tcut_field(poly, "phi", [](double x, double, double) { return x - 0.4; });
    EXPECT_THROW((void)mio::cut_level_set(poly, o), mio::Unsupported);
}
