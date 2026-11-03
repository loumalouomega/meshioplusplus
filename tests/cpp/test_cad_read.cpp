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
 * @file test_cad_read.cpp
 * @brief STEP and IGES reading through the optional OpenCASCADE backend: a
 *        faceted STEP solid read back as triangles and NURBS faces, written as
 *        Kratos CAD JSON and sampled again; and, without the backend, the
 *        named `Unsupported`.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/formats/iges.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::BrepInfo;
using meshioplusplus::Mesh;

// The closed unit cube as six outward quads.
Mesh unit_cube() {
    return mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                          {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
                         "quad",
                         {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                          {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}});
}

}  // namespace

TEST(CadRead, CapabilityMatchesTheBuild) {
#ifdef MESHIOPLUSPLUS_HAS_OCCT
    EXPECT_TRUE(meshioplusplus::cad_read_has_occt());
#else
    EXPECT_FALSE(meshioplusplus::cad_read_has_occt());
#endif
    EXPECT_EQ(meshioplusplus::resolve_format("part.igs", ""), "iges");
    EXPECT_EQ(meshioplusplus::resolve_format("part.stp", ""), "step");
    EXPECT_TRUE(meshioplusplus::registry_readers().count("step"));
    EXPECT_TRUE(meshioplusplus::registry_readers().count("iges"));
}

#ifndef MESHIOPLUSPLUS_HAS_OCCT

TEST(CadRead, AbsentBackendIsNamed) {
    const std::string path = mt::temp_path(".step");
    meshioplusplus::write_step(path, unit_cube());
    for (bool iges : {false, true}) {
        try {
            if (iges)
                meshioplusplus::read_iges(path);
            else
                meshioplusplus::read_step(path);
            FAIL() << "expected Unsupported";
        } catch (const meshioplusplus::Unsupported& rErr) {
            EXPECT_NE(std::string(rErr.what()).find("MESHIOPLUSPLUS_WITH_OCCT"),
                      std::string::npos);
        }
    }
}

#else

TEST(CadRead, StepCubeReadsAsSixFaces) {
    const std::string path = mt::temp_path(".step");
    meshioplusplus::write_step(path, unit_cube());
    BrepInfo info;
    const Mesh mesh = meshioplusplus::read_step(path, info);
    ASSERT_EQ(info.mFaces.size(), 6u);
    EXPECT_EQ(info.mEdges.size(), 12u);
    EXPECT_EQ(info.mPatches.size(), 6u);
    for (const auto& f : info.mFaces) {
        EXPECT_TRUE(f.mTrimmed);
        ASSERT_EQ(f.mLoops.size(), 1u);
        EXPECT_EQ(f.mLoops[0].mType, "outer");
        EXPECT_EQ(f.mLoops[0].mTrims.size(), 4u);
        EXPECT_EQ(f.mControlPointIds.Size(), info.mPatches[static_cast<std::size_t>(f.mPatch)].mWeights.Size());
    }
    for (const auto& e : info.mEdges) {
        EXPECT_TRUE(e.mHasCurve);
        EXPECT_EQ(e.mTopology.size(), 2u);  // each cube edge bounds two faces
    }
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "triangle");
    EXPECT_GE(mesh.Cells(0).NumCells(), 12u);
    ASSERT_TRUE(mesh.HasCellData("cad:id"));
    // The file is in millimetres (no unit stated), and the mesh says so.
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::mesh_length_unit(mesh), 1e-3);
    // Every triangulated point lies on the cube.
    const auto* p = mesh.Points().As<double>();
    for (std::size_t i = 0; i < mesh.NumPoints(); ++i)
        for (int d = 0; d < 3; ++d) {
            EXPECT_GE(p[i * 3 + static_cast<std::size_t>(d)], -1e-9);
            EXPECT_LE(p[i * 3 + static_cast<std::size_t>(d)], 1.0 + 1e-9);
        }

    // As Kratos CAD JSON, read back and sampled: the same six faces.
    const std::string cad = mt::temp_path(".cad.json");
    meshioplusplus::write_cad_json(cad, mesh, info);
    BrepInfo back;
    const Mesh sampled = meshioplusplus::read_cad_json(cad, back);
    EXPECT_EQ(back.mFaces.size(), 6u);
    EXPECT_GT(sampled.NumPoints(), 0u);
    const auto* q = sampled.Points().As<double>();
    for (std::size_t i = 0; i < sampled.NumPoints(); ++i)
        for (int d = 0; d < 3; ++d) {
            EXPECT_GE(q[i * 3 + static_cast<std::size_t>(d)], -1e-9);
            EXPECT_LE(q[i * 3 + static_cast<std::size_t>(d)], 1.0 + 1e-9);
        }
}

TEST(CadRead, ReadsAreInMillimetresWhateverTheFileUnit) {
    // A cube of 1 m, of 1 inch and of 1 mm: OpenCASCADE converts each to
    // millimetres and the mesh states 0.001, so the physical size survives.
    for (const double unit : {1.0, 0.0254, 1e-3}) {
        Mesh cube = unit_cube();
        meshioplusplus::detail::set_mesh_length_unit(cube, unit);
        const std::string path = mt::temp_path(".step");
        meshioplusplus::write_step(path, cube);
        const Mesh mesh = meshioplusplus::read_step(path);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::mesh_length_unit(mesh), 1e-3) << unit;
        double hi = 0.0;
        const auto* p = mesh.Points().As<double>();
        for (std::size_t i = 0; i < mesh.NumPoints() * 3; ++i)
            hi = std::max(hi, p[i]);
        EXPECT_NEAR(hi, unit * 1000.0, 1e-6 * unit * 1000.0) << unit;
    }
}

namespace {

// Two planar 100 x 100 mm quads side by side whose shared edge is written
// twice, the second copy `Gap` away from the first (roadmap 7.7.1): each face
// has its own four vertices and edges, in one open shell.
std::string cr_two_face_step(double Gap) {
    const double pts[8][3] = {{0, 0, 0},     {100, 0, 0},          {100, 100, 0},
                              {0, 100, 0},   {100 + Gap, 0, Gap},  {200, 0, 0},
                              {200, 100, 0}, {100 + Gap, 100, Gap}};
    std::string out =
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('two faces'),'2;1');\n"
        "FILE_NAME('two','1970-01-01T00:00:00',(''),(''),'','','');\n"
        "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n";
    int n = 0;
    const auto ent = [&](const std::string& rBody) {
        out += "#" + std::to_string(++n) + "=" + rBody + ";\n";
        return n;
    };
    const auto ref = [](int Id) { return "#" + std::to_string(Id); };
    const auto num = [](double V) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.17g", V);
        std::string s = buf;
        if (s.find_first_of(".eE") == std::string::npos)
            s += ".";
        return s;
    };
    const auto triple = [&](const double* pV) {
        return "(" + num(pV[0]) + "," + num(pV[1]) + "," + num(pV[2]) + ")";
    };
    const int app =
        ent("APPLICATION_CONTEXT('core data for automotive mechanical design "
            "processes')");
    ent("APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000," +
        ref(app) + ")");
    const int pctx = ent("PRODUCT_CONTEXT(''," + ref(app) + ",'mechanical')");
    const int pdctx =
        ent("PRODUCT_DEFINITION_CONTEXT('part definition'," + ref(app) + ",'design')");
    const int prod = ent("PRODUCT('two','two','',(" + ref(pctx) + "))");
    const int form = ent("PRODUCT_DEFINITION_FORMATION('',''," + ref(prod) + ")");
    const int pdef = ent("PRODUCT_DEFINITION('design',''," + ref(form) + "," + ref(pdctx) + ")");
    const int pshape = ent("PRODUCT_DEFINITION_SHAPE(''," + ref(pdef) + ")");
    const int len = ent("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))");
    const int ang = ent("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))");
    const int sol = ent("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())");
    const int unc = ent("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07)," + ref(len) +
                        ",'distance_accuracy_value','confusion accuracy')");
    const int gctx =
        ent("(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_"
            "CONTEXT((" +
            ref(unc) + "))GLOBAL_UNIT_ASSIGNED_CONTEXT((" + ref(len) + "," + ref(ang) + "," +
            ref(sol) + "))REPRESENTATION_CONTEXT('Context #1','3D Context'))");
    const int o_pt = ent("CARTESIAN_POINT('',(0.,0.,0.))");
    const int o_z = ent("DIRECTION('',(0.,0.,1.))");
    const int o_x = ent("DIRECTION('',(1.,0.,0.))");
    const int origin =
        ent("AXIS2_PLACEMENT_3D(''," + ref(o_pt) + "," + ref(o_z) + "," + ref(o_x) + ")");
    int point[8], vertex[8];
    for (int i = 0; i < 8; ++i) {
        point[i] = ent("CARTESIAN_POINT(''," + triple(pts[i]) + ")");
        vertex[i] = ent("VERTEX_POINT(''," + ref(point[i]) + ")");
    }
    std::string faces;
    for (int f = 0; f < 2; ++f) {
        std::string loop;
        for (int k = 0; k < 4; ++k) {
            const int a = 4 * f + k, b = 4 * f + (k + 1) % 4;
            double d[3] = {pts[b][0] - pts[a][0], pts[b][1] - pts[a][1], pts[b][2] - pts[a][2]};
            const double ln = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            for (double& rV : d)
                rV /= ln;
            const int dir = ent("DIRECTION(''," + triple(d) + ")");
            const int vec = ent("VECTOR(''," + ref(dir) + "," + num(ln) + ")");
            const int line = ent("LINE(''," + ref(point[a]) + "," + ref(vec) + ")");
            const int edge = ent("EDGE_CURVE(''," + ref(vertex[a]) + "," + ref(vertex[b]) + "," +
                                 ref(line) + ",.T.)");
            const int oe = ent("ORIENTED_EDGE('',*,*," + ref(edge) + ",.T.)");
            loop += (k ? "," : "") + ref(oe);
        }
        const int el = ent("EDGE_LOOP('',(" + loop + "))");
        const int bound = ent("FACE_OUTER_BOUND(''," + ref(el) + ",.T.)");
        const int place = ent("AXIS2_PLACEMENT_3D(''," + ref(point[4 * f]) + "," + ref(o_z) + "," +
                              ref(o_x) + ")");
        const int plane = ent("PLANE(''," + ref(place) + ")");
        const int face = ent("ADVANCED_FACE('',(" + ref(bound) + ")," + ref(plane) + ",.T.)");
        faces += (f ? "," : "") + ref(face);
    }
    const int shell = ent("OPEN_SHELL('',(" + faces + "))");
    const int model = ent("SHELL_BASED_SURFACE_MODEL('',(" + ref(shell) + "))");
    const int rep = ent("MANIFOLD_SURFACE_SHAPE_REPRESENTATION('',(" + ref(origin) + "," +
                        ref(model) + ")," + ref(gctx) + ")");
    ent("SHAPE_DEFINITION_REPRESENTATION(" + ref(pshape) + "," + ref(rep) + ")");
    out += "ENDSEC;\nEND-ISO-10303-21;\n";
    return out;
}

// Edges of the B-rep that bound exactly one face.
std::size_t cr_open_edges(const BrepInfo& rInfo) {
    std::size_t open = 0;
    for (const auto& r_edge : rInfo.mEdges)
        if (r_edge.mTopology.size() == 1)
            ++open;
    return open;
}

}  // namespace

TEST(CadRead, HealSewsADuplicatedSharedEdge) {
    const std::string path = mt::temp_path(".step");
    {
        std::ofstream out(path);
        out << cr_two_face_step(1e-4);
    }
    BrepInfo raw;
    meshioplusplus::read_step(path, raw);
    ASSERT_EQ(raw.mFaces.size(), 2u);
    EXPECT_EQ(cr_open_edges(raw), 8u);  // the shared edge is two open edges
    EXPECT_EQ(raw.mExtra.find("heal"), std::string::npos);

    meshioplusplus::CadReadOptions heal;
    heal.mHeal = true;
    BrepInfo healed;
    const Mesh mesh = meshioplusplus::read_step(path, healed, heal);
    ASSERT_EQ(healed.mFaces.size(), 2u);
    EXPECT_EQ(cr_open_edges(healed), 6u);  // only the outer boundary is open
    EXPECT_GT(mesh.NumPoints(), 0u);
    EXPECT_NE(healed.mExtra.find("\"heal\":{\"healed\":true"), std::string::npos) << healed.mExtra;
    EXPECT_NE(healed.mExtra.find("\"free_edges_before\":8"), std::string::npos);
    EXPECT_NE(healed.mExtra.find("\"free_edges_after\":6"), std::string::npos);
    EXPECT_NE(healed.mExtra.find("\"version_number\":1"), std::string::npos);

    // A positive tolerance sews without heal.
    meshioplusplus::CadReadOptions sew_only;
    sew_only.mSewTolerance = 1e-3;
    BrepInfo sewn;
    meshioplusplus::read_step(path, sewn, sew_only);
    EXPECT_EQ(cr_open_edges(sewn), 6u);
    EXPECT_NE(sewn.mExtra.find("\"healed\":false"), std::string::npos) << sewn.mExtra;

    meshioplusplus::CadReadOptions bad;
    bad.mSewTolerance = -1.0;
    EXPECT_THROW(meshioplusplus::read_step(path, healed, bad), std::invalid_argument);
}

TEST(CadRead, MalformedFilesAreReadErrors) {
    const std::string path = mt::temp_path(".step");
    {
        std::ofstream out(path);
        out << "not a step file\n";
    }
    BrepInfo info;
    EXPECT_THROW(meshioplusplus::read_step(path, info), meshioplusplus::ReadError);
    EXPECT_THROW(meshioplusplus::read_iges(path, info), meshioplusplus::ReadError);
    meshioplusplus::CadReadOptions bad;
    bad.mAngularDeflection = 0.0;
    EXPECT_THROW(meshioplusplus::read_step(path, info, bad), std::invalid_argument);
}

#endif
