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
// mirror: reflections whose cells come out right side out (roadmap 3.2.5),
// and the transform probe that showed a reflecting matrix does not.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/cell_flip.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/check_conformity.hpp"
#include "meshioplusplus/operations/mirror.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/transform.hpp"

using meshioplusplus::CellType;
using meshioplusplus::Mesh;
using meshioplusplus::MirrorOptions;
using meshioplusplus::MirrorResult;
using meshioplusplus::NDArray;
using meshioplusplus::RegionKind;

namespace {

using Vec3 = meshioplusplus::detail::Vec3;

// The signed volume of every 3-D cell, block-major.
std::vector<double> mir_volumes(const Mesh& rM) {
    std::vector<double> out;
    meshioplusplus::detail::CellRings rings;
    std::vector<Vec3> xyz;
    for (const auto cb : rM.CellRange())
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            if (meshioplusplus::detail::cell_rings(cb, c, rM.Points(), rM.PointDim(), rings, xyz))
                out.push_back(meshioplusplus::detail::poly_measure(rings, xyz.data()).mVolume);
    return out;
}

double mir_sum(const std::vector<double>& rV) {
    double s = 0.0;
    for (const double v : rV)
        s += v;
    return s;
}

Vec3 mir_point(const Mesh& rM, std::int64_t P) {
    Vec3 x{0.0, 0.0, 0.0};
    for (std::size_t d = 0; d < rM.PointDim() && d < 3; ++d)
        x[d] = meshioplusplus::detail::read_double(rM.Points(),
                                                   static_cast<std::size_t>(P) * rM.PointDim() + d);
    return x;
}

// Per 2-D cell, the z component of its corner ring's area vector.
std::vector<double> mir_areas_z(const Mesh& rM) {
    std::vector<double> out;
    for (const auto cb : rM.CellRange()) {
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::vector<std::int64_t> ring;
            if (cb.IsRagged()) {
                ring.assign(cb.Row(c), cb.Row(c) + cb.RowSize(c));
            } else {
                const CellType ct = meshioplusplus::cell_type_from_name(std::string(cb.Type()));
                const auto& r_edges = meshioplusplus::detail::cell_edges(ct);
                for (const auto& r_e : r_edges)
                    ring.push_back(meshioplusplus::detail::read_int(
                        cb.Conn(), c * cb.NodesPerCell() + r_e.mNodes[0]));
            }
            double a = 0.0;
            for (std::size_t k = 0; k < ring.size(); ++k) {
                const Vec3 p = mir_point(rM, ring[k]);
                const Vec3 q = mir_point(rM, ring[(k + 1) % ring.size()]);
                a += 0.5 * (p[0] * q[1] - p[1] * q[0]);
            }
            out.push_back(a);
        }
    }
    return out;
}

std::vector<std::int64_t> mir_entries(const meshioplusplus::Region& rR) {
    const std::int64_t* e = rR.Entries();
    return e ? std::vector<std::int64_t>(e, e + rR.mEntries.Size()) : std::vector<std::int64_t>{};
}

// The sorted node set of each Side entry of a region.
std::vector<std::vector<std::int64_t>> mir_side_nodes(const Mesh& rM,
                                                      const meshioplusplus::Region& rR) {
    std::vector<std::vector<std::int64_t>> out;
    const std::vector<std::int64_t> e = mir_entries(rR);
    for (std::size_t k = 0; k + 1 < e.size(); k += 2) {
        CellType t = CellType::Custom;
        std::vector<std::int64_t> nodes;
        EXPECT_TRUE(meshioplusplus::detail::facet_nodes(rM, e[k], e[k + 1], t, nodes));
        std::sort(nodes.begin(), nodes.end());
        out.push_back(nodes);
    }
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace

// The roadmap's verify-first probe: a reflecting transform inverts every cell,
// which is why `mirror` exists.
TEST(Mirror, ReflectingTransformInvertsEveryCell) {
    const Mesh box = meshioplusplus::box({{2, 2, 2}}, {{1.0, 1.0, 1.0}}, {{0.0, 0.0, 0.0}}, "tetra");
    const Mesh flipped =
        meshioplusplus::transform(box, meshioplusplus::transform_scale(-1.0, 1.0, 1.0));
    const auto q = meshioplusplus::compute_quality(flipped);
    EXPECT_EQ(q.mNumInverted, q.mNumCells);
    EXPECT_GT(q.mNumCells, 0);
}

TEST(Mirror, FlipTablesAreInvolutionsWithBijectiveFacetMaps) {
    using CT = CellType;
    for (const CT t : {CT::Vertex, CT::Line, CT::Line3, CT::Line6, CT::Triangle, CT::Triangle6,
                       CT::Triangle7, CT::Quad, CT::Quad8, CT::Quad9, CT::Tetra, CT::Tetra10,
                       CT::Pyramid, CT::Pyramid13, CT::Pyramid14, CT::Wedge, CT::Wedge15,
                       CT::Wedge18, CT::Hexahedron, CT::Hexahedron20, CT::Hexahedron27}) {
        const auto* p = meshioplusplus::detail::cell_flip_permutation(t);
        const auto* f = meshioplusplus::detail::cell_flip_facet_map(t);
        ASSERT_NE(p, nullptr) << meshioplusplus::cell_type_name(t);
        ASSERT_NE(f, nullptr);
        ASSERT_EQ(static_cast<int>(p->size()), meshioplusplus::cell_type_num_nodes(t));
        for (std::size_t i = 0; i < p->size(); ++i)
            EXPECT_EQ((*p)[static_cast<std::size_t>((*p)[i])], static_cast<int>(i))
                << meshioplusplus::cell_type_name(t);
        std::set<int> image(f->begin(), f->end());
        EXPECT_EQ(image.size(), f->size());
    }
    EXPECT_EQ(meshioplusplus::detail::cell_flip_permutation(CT::Triangle10), nullptr);
    EXPECT_EQ(meshioplusplus::detail::cell_flip_permutation(CT::VtkLagrangeTetrahedron), nullptr);
    // The facet maps of the linear solids, as derived by hand from cell_faces.
    EXPECT_EQ(*meshioplusplus::detail::cell_flip_facet_map(CT::Tetra),
              (std::vector<int>{2, 1, 0, 3}));
    EXPECT_EQ(*meshioplusplus::detail::cell_flip_facet_map(CT::Hexahedron),
              (std::vector<int>{2, 3, 0, 1, 4, 5}));
    EXPECT_EQ(*meshioplusplus::detail::cell_flip_facet_map(CT::Wedge),
              (std::vector<int>{0, 1, 4, 3, 2}));
    EXPECT_EQ(*meshioplusplus::detail::cell_flip_facet_map(CT::Pyramid),
              (std::vector<int>{0, 4, 3, 2, 1}));
    EXPECT_EQ(*meshioplusplus::detail::cell_flip_facet_map(CT::Quad),
              (std::vector<int>{3, 2, 1, 0}));
}

TEST(Mirror, BoxesComeOutRightSideOut) {
    for (const char* type : {"tetra", "hexahedron", "wedge"}) {
        for (int order : {1, 2}) {
            const Mesh box = meshioplusplus::box({{2, 3, 2}}, {{1.0, 2.0, 0.5}},
                                                 {{0.25, 0.0, 0.0}}, type, order);
            MirrorOptions o;
            o.mOrigin = {{0.1, 0.2, 0.3}};
            o.mNormal = {{1.0, 2.0, -0.5}};
            const MirrorResult r = meshioplusplus::mirror(box, o);
            const std::vector<double> before = mir_volumes(box), after = mir_volumes(r.mMesh);
            ASSERT_EQ(before.size(), after.size());
            for (const double v : after)
                EXPECT_GT(v, 0.0) << type << order;
            EXPECT_NEAR(mir_sum(after), mir_sum(before), 1e-12) << type;
            EXPECT_EQ(meshioplusplus::compute_quality(r.mMesh).mNumInverted, 0) << type << order;
            EXPECT_EQ(r.mNumFlippedCells, static_cast<std::int64_t>(before.size()));
            EXPECT_EQ(r.mMesh.NumPoints(), box.NumPoints());
            // Every Side region names the same nodes as before.
            for (const char* side : {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}) {
                const auto& r_in = box.Region(box.FindRegion(side, RegionKind::Side));
                const auto& r_out = r.mMesh.Region(r.mMesh.FindRegion(side, RegionKind::Side));
                EXPECT_EQ(mir_side_nodes(box, r_in), mir_side_nodes(r.mMesh, r_out))
                    << type << order << side;
            }
        }
    }
}

// Mid-edge and centre nodes stay where they belong: every facet's k-th mid
// node is the midpoint of its corners k and k + 1 after the flip, as before.
TEST(Mirror, QuadraticNodesStayInTheirSlots) {
    std::vector<Mesh> meshes;
    meshes.push_back(meshioplusplus::reference_cells(
        {"line3", "triangle6", "quad8", "tetra10", "hexahedron20", "wedge15", "pyramid13"}));
    meshes.push_back(mt::wedge18_mesh());
    meshes.push_back(mt::pyramid14_mesh());
    for (const Mesh& m : meshes) {
        MirrorOptions o;
        o.mNormal = {{0.3, -1.0, 0.2}};
        const MirrorResult r = meshioplusplus::mirror(m, o);
        const Mesh& out = r.mMesh;
        for (std::size_t b = 0; b < out.NumCellBlocks(); ++b) {
            const auto cb = out.Cells(b);
            const CellType ct = meshioplusplus::cell_type_from_name(std::string(cb.Type()));
            const std::size_t npc = cb.NodesPerCell();
            const auto node = [&](std::size_t C, std::size_t K) {
                return mir_point(out, meshioplusplus::detail::read_int(cb.Conn(), C * npc + K));
            };
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                const auto check = [&](const std::uint8_t* pNodes, std::size_t NC, std::size_t NN) {
                    for (std::size_t k = 0; k < NC && NC + k < NN; ++k) {
                        const Vec3 a = node(c, pNodes[k]), z = node(c, pNodes[(k + 1) % NC]);
                        const Vec3 mid = node(c, pNodes[NC + k]);
                        for (std::size_t d = 0; d < 3; ++d)
                            EXPECT_NEAR(mid[d], 0.5 * (a[d] + z[d]), 1e-12) << cb.Type();
                    }
                };
                if (meshioplusplus::cell_type_dimension(ct) == 3)
                    for (const auto& r_f : meshioplusplus::detail::cell_faces(ct))
                        check(r_f.mNodes.data(), r_f.mNumCorners,
                              std::min<std::size_t>(r_f.mNumNodes, 2u * r_f.mNumCorners));
                else if (meshioplusplus::cell_type_dimension(ct) == 2)
                    for (const auto& r_e : meshioplusplus::detail::cell_edges(ct)) {
                        if (r_e.mNumNodes < 3)
                            continue;
                        const Vec3 a = node(c, r_e.mNodes[0]), z = node(c, r_e.mNodes[1]);
                        const Vec3 mid = node(c, r_e.mNodes[2]);
                        for (std::size_t d = 0; d < 3; ++d)
                            EXPECT_NEAR(mid[d], 0.5 * (a[d] + z[d]), 1e-12) << cb.Type();
                    }
            }
            if (meshioplusplus::cell_type_dimension(ct) == 3)
                for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                    meshioplusplus::detail::CellRings rings;
                    std::vector<Vec3> xyz;
                    ASSERT_TRUE(meshioplusplus::detail::cell_rings(cb, c, out.Points(),
                                                                   out.PointDim(), rings, xyz));
                    EXPECT_GT(meshioplusplus::detail::poly_measure(rings, xyz.data()).mVolume,
                              0.0)
                        << cb.Type();
                }
        }
    }
}

TEST(Mirror, PolygonsAndPolyhedraAreRewound) {
    const Mesh hc = meshioplusplus::honeycomb({{3, 2}});
    for (const double a : mir_areas_z(hc))
        ASSERT_GT(a, 0.0);
    MirrorOptions o;
    o.mNormal = {{1.0, 0.0, 0.0}};
    const MirrorResult r = meshioplusplus::mirror(hc, o);
    for (const double a : mir_areas_z(r.mMesh))
        EXPECT_GT(a, 0.0);
    const auto& r_in = hc.Region(hc.FindRegion("boundary", RegionKind::Side));
    const auto& r_out = r.mMesh.Region(r.mMesh.FindRegion("boundary", RegionKind::Side));
    EXPECT_EQ(mir_side_nodes(hc, r_in), mir_side_nodes(r.mMesh, r_out));

    const Mesh poly = meshioplusplus::reference_cells({"polyhedron"});
    const MirrorResult rp = meshioplusplus::mirror(poly, o);
    for (const double v : mir_volumes(rp.mMesh))
        EXPECT_GT(v, 0.0);
    EXPECT_NEAR(mir_sum(mir_volumes(rp.mMesh)), mir_sum(mir_volumes(poly)), 1e-12);
}

// A half box mirrored and welded: one conforming box, the plane's region empty.
TEST(Mirror, KeepOriginalAndWeldMakeAWholeModel) {
    const Mesh half = meshioplusplus::box({{2, 2, 2}}, {{1.0, 1.0, 1.0}}, {{0.0, 0.0, 0.0}}, "tetra");
    MirrorOptions o;
    o.mKeepOriginal = true;
    const MirrorResult r = meshioplusplus::mirror(half, o);
    const std::int64_t n = static_cast<std::int64_t>(half.NumPoints());
    EXPECT_EQ(r.mNumPointsOnPlane, 9);
    EXPECT_EQ(static_cast<std::int64_t>(r.mMesh.NumPoints()), 2 * n - 9);
    const std::vector<double> vols = mir_volumes(r.mMesh);
    EXPECT_EQ(vols.size(), 2 * mir_volumes(half).size());
    for (const double v : vols)
        EXPECT_GT(v, 0.0);
    EXPECT_NEAR(mir_sum(vols), 2.0, 1e-12);
    const auto conf = meshioplusplus::check_conformity(r.mMesh);
    EXPECT_TRUE(conf.mConforming);
    EXPECT_EQ(conf.mNumUnweldedFacets, 0);
    // The [-1, 1] x [0, 1]^2 box's skin, two triangles per hexahedron face:
    // 2 x 4 on the x sides, 2 x 16 on the y sides and 2 x 16 on the z sides.
    EXPECT_EQ(conf.mNumBoundaryFacets, 8 + 8 + 16 + 16 + 16 + 16);
    const std::size_t xmin = r.mMesh.FindRegion("xmin", RegionKind::Side);
    ASSERT_NE(xmin, Mesh::npos);
    EXPECT_EQ(r.mMesh.Region(xmin).NumEntries(), 0u);
    EXPECT_EQ(r.mNumDroppedSides, 8);
    EXPECT_FALSE(r.mWarnings.empty());
    // xmax holds both copies: the half's x = 1 face and its image at x = -1.
    const auto& r_xmax = r.mMesh.Region(r.mMesh.FindRegion("xmax", RegionKind::Side));
    EXPECT_EQ(r_xmax.NumEntries(), 16u);
    // Without welding, the plane's points are duplicated and nothing is dropped.
    o.mWeld = false;
    const MirrorResult rn = meshioplusplus::mirror(half, o);
    EXPECT_EQ(static_cast<std::int64_t>(rn.mMesh.NumPoints()), 2 * n);
    EXPECT_EQ(rn.mNumDroppedSides, 0);
}

// The hemisphere probe: a half sphere's surface mirrored onto its open rim
// closes into a sphere whose every edge is shared by two consistently wound
// triangles.
TEST(Mirror, HemisphereClosesIntoASphere) {
    // A UV sphere with an even ring count has its equator on z = 0.
    const Mesh sphere = meshioplusplus::uv_sphere(1.0, 16, 8);
    // The cells of the upper half (all corners at z >= 0).
    std::vector<std::vector<std::int64_t>> rows;
    const auto cb = sphere.Cells(0);
    std::vector<std::vector<double>> pts;
    for (std::size_t p = 0; p < sphere.NumPoints(); ++p) {
        const Vec3 x = mir_point(sphere, static_cast<std::int64_t>(p));
        pts.push_back({x[0], x[1], std::abs(x[2]) < 1e-12 ? 0.0 : x[2]});
    }
    for (std::size_t c = 0; c < cb.NumCells(); ++c) {
        std::vector<std::int64_t> row;
        bool upper = true;
        for (std::size_t k = 0; k < 3; ++k) {
            row.push_back(meshioplusplus::detail::read_int(cb.Conn(), c * 3 + k));
            upper = upper && pts[static_cast<std::size_t>(row.back())][2] >= 0.0;
        }
        if (upper)
            rows.push_back(row);
    }
    const Mesh hemi = mt::make_mesh(pts, "triangle", rows);
    MirrorOptions o;
    o.mNormal = {{0.0, 0.0, 1.0}};
    o.mKeepOriginal = true;
    const MirrorResult r = meshioplusplus::mirror(hemi, o);
    std::set<std::pair<std::int64_t, std::int64_t>> directed;
    const auto ob = r.mMesh.Cells(0);
    double volume = 0.0;
    for (std::size_t c = 0; c < ob.NumCells(); ++c) {
        std::array<std::int64_t, 3> t{};
        for (std::size_t k = 0; k < 3; ++k)
            t[k] = meshioplusplus::detail::read_int(ob.Conn(), c * 3 + k);
        for (std::size_t k = 0; k < 3; ++k)
            EXPECT_TRUE(directed.insert({t[k], t[(k + 1) % 3]}).second);
        const Vec3 a = mir_point(r.mMesh, t[0]), b = mir_point(r.mMesh, t[1]),
                   z = mir_point(r.mMesh, t[2]);
        volume += (a[0] * (b[1] * z[2] - b[2] * z[1]) - a[1] * (b[0] * z[2] - b[2] * z[0]) +
                   a[2] * (b[0] * z[1] - b[1] * z[0])) /
                  6.0;
    }
    for (const auto& [u, v] : directed)
        EXPECT_TRUE(directed.count({v, u})) << u << " " << v;
    EXPECT_GT(volume, 0.0);
}

TEST(Mirror, RegionsAndDataFollowTheCopies) {
    Mesh m = meshioplusplus::box({{1, 1, 2}}, {{1.0, 1.0, 2.0}}, {{1.0, 0.0, 0.0}}, "hexahedron");
    NDArray vec(meshioplusplus::DType::Float64, {m.NumPoints(), 3});
    for (std::size_t p = 0; p < m.NumPoints(); ++p) {
        vec.As<double>()[3 * p] = 1.0;
        vec.As<double>()[3 * p + 1] = 2.0;
        vec.As<double>()[3 * p + 2] = 3.0;
    }
    m.AddPointData("v", vec);
    NDArray tag(meshioplusplus::DType::Int32, {2});
    tag.As<std::int32_t>()[0] = 7;
    tag.As<std::int32_t>()[1] = 9;
    m.AddCellData("tag", {tag});
    NDArray cells(meshioplusplus::DType::Int64, {1});
    cells.As<std::int64_t>()[0] = 1;
    m.AddRegion(meshioplusplus::Region("top_cell", RegionKind::Cell, cells));
    MirrorOptions o;
    o.mKeepOriginal = true;
    o.mReflectVectorData = true;
    const MirrorResult r = meshioplusplus::mirror(m, o);
    // No point is on x = 0, so nothing is welded.
    EXPECT_EQ(r.mMesh.NumPoints(), 2 * m.NumPoints());
    const NDArray& v = r.mMesh.PointData("v");
    EXPECT_EQ(v.As<double>()[0], 1.0);
    EXPECT_EQ(v.As<double>()[3 * m.NumPoints()], -1.0);
    EXPECT_EQ(v.As<double>()[3 * m.NumPoints() + 1], 2.0);
    const NDArray& t = r.mMesh.CellData("tag", 0);
    ASSERT_EQ(t.Size(), 4u);
    EXPECT_EQ(meshioplusplus::detail::read_int(t, 0), 7);
    EXPECT_EQ(meshioplusplus::detail::read_int(t, 3), 9);
    const auto& r_cell = r.mMesh.Region(r.mMesh.FindRegion("top_cell", RegionKind::Cell));
    EXPECT_EQ(mir_entries(r_cell), (std::vector<std::int64_t>{1, 3}));
    // Mirroring the image back gives the input's connectivity bit for bit.
    MirrorOptions once;
    once.mNormal = {{0.0, 1.0, 0.0}};
    const Mesh twice =
        meshioplusplus::mirror(meshioplusplus::mirror(m, once).mMesh, once).mMesh;
    for (std::size_t i = 0; i < m.Cells(0).Conn().Size(); ++i)
        EXPECT_EQ(meshioplusplus::detail::read_int(twice.Cells(0).Conn(), i),
                  meshioplusplus::detail::read_int(m.Cells(0).Conn(), i));
}

TEST(Mirror, Refusals) {
    MirrorOptions o;
    o.mNormal = {{0.0, 0.0, 0.0}};
    const Mesh tet = mt::tet_mesh();
    EXPECT_THROW(meshioplusplus::mirror(tet, o), std::invalid_argument);
    o.mNormal = {{0.0, 0.0, 1.0}};
    o.mTolerance = -1.0;
    EXPECT_THROW(meshioplusplus::mirror(tet, o), std::invalid_argument);
    o.mTolerance = 1e-9;
    // A planar mesh cannot be reflected out of its plane.
    EXPECT_THROW(meshioplusplus::mirror(mt::tri_mesh_2d(), o),
                 std::invalid_argument);
    o.mNormal = {{0.0, 1.0, 0.0}};
    const MirrorResult r2 = meshioplusplus::mirror(mt::tri_mesh_2d(), o);
    EXPECT_EQ(r2.mMesh.PointDim(), 2u);
    for (const double a : mir_areas_z(r2.mMesh))
        EXPECT_GT(a, 0.0);
    std::vector<std::vector<double>> pts;
    for (int i = 0; i < 10; ++i)
        pts.push_back({0.1 * i, 0.05 * i * i, 0.0});
    const Mesh t10 = mt::make_mesh(
        pts, "triangle10", {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}});
    EXPECT_THROW(meshioplusplus::mirror(t10, o), meshioplusplus::Unsupported);
}
