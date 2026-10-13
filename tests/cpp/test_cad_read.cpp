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
#include <cmath>
#include <fstream>
#include <string>

// Project includes
#include "mesh_fixtures.hpp"
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
