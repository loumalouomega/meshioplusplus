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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// Unit tests for detail/external_mesher.cpp (roadmap 3.11.5): exercises TetGen,
// Gmsh, fTetWild, QuadriFlow and Instant Meshes command construction, mesh I/O,
// shell/hole handling, point matching and error paths via mock shell drivers.

#if !defined(_WIN32) && !defined(__EMSCRIPTEN__)

// System includes
#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/providers.hpp"

namespace meshioplusplus::providers {
std::vector<TetMesherProvider> make_external_tet_meshers();
std::vector<QuadRemesherProvider> make_external_quad_remeshers();
}  // namespace meshioplusplus::providers

namespace {

using meshioplusplus::QuadRemeshRequest;
using meshioplusplus::QuadRemeshResult;
using meshioplusplus::TetMeshRegionRequest;
using meshioplusplus::TetMeshRegionResult;

struct EmEnvGuard {
    std::string mKey;
    bool mHadOld = false;
    std::string mOldVal;

    EmEnvGuard(const char* pKey, const std::string& rVal) : mKey(pKey) {
        if (const char* p_old = std::getenv(pKey)) {
            mHadOld = true;
            mOldVal = p_old;
        }
        ::setenv(pKey, rVal.c_str(), 1);
    }
    ~EmEnvGuard() {
        if (mHadOld)
            ::setenv(mKey.c_str(), mOldVal.c_str(), 1);
        else
            ::unsetenv(mKey.c_str());
    }
    EmEnvGuard(const EmEnvGuard&) = delete;
    EmEnvGuard& operator=(const EmEnvGuard&) = delete;
};

struct EmTempDir {
    std::filesystem::path mDir;
    explicit EmTempDir(const std::string& rTag) {
        std::string tmpl = (std::filesystem::temp_directory_path() / (rTag + "_XXXXXX")).string();
        std::vector<char> buf(tmpl.begin(), tmpl.end());
        buf.push_back('\0');
        char* p_res = ::mkdtemp(buf.data());
        if (!p_res)
            throw std::runtime_error("mkdtemp failed");
        mDir = p_res;
    }
    ~EmTempDir() {
        std::error_code ec;
        std::filesystem::remove_all(mDir, ec);
    }
    std::string WriteScript(const std::string& rName, const std::string& rBody) const {
        const auto path = mDir / rName;
        std::ofstream out(path);
        out << "#!/bin/sh\nset -eu\n" << rBody << "\n";
        out.close();
        ::chmod(path.c_str(), 0755);
        return path.string();
    }
};

// Single positively-oriented tetrahedron surface (outward-wound faces).
TetMeshRegionRequest em_test_single_tet_request() {
    TetMeshRegionRequest req;
    req.mPoints = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    req.mTriangles = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
    req.mMaxSize = 0.5;
    req.mOptimizeSteps = 2;
    return req;
}

// Two concentric shells (outer tetrahedron + small inner hole tetrahedron) to
// exercise em_shells() hole detection and .smesh hole injection.
TetMeshRegionRequest em_test_hollow_request() {
    TetMeshRegionRequest req;
    req.mPoints = {
        {0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0},
        {1.0, 1.0, 1.0}, {2.0, 1.0, 1.0},  {1.0, 2.0, 1.0},  {1.0, 1.0, 2.0},
    };
    req.mTriangles = {
        {0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {4, 5, 6}, {4, 7, 5}, {5, 7, 6}, {6, 7, 4},
    };
    req.mMaxSize = 2.0;
    req.mOptimizeSteps = 0;
    return req;
}

const char* kMsh22SingleTet = R"MSH($MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
1
1 4 2 0 1 1 2 3 4
$EndElements
)MSH";

const char* kMsh22Empty = R"MSH($MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
0
$EndElements
)MSH";

TEST(ExternalMesher, DisabledByEnvVar) {
    EmTempDir tmp("em_off");
    const std::string fake = tmp.WriteScript("fake_tetgen", "exit 0");
    EmEnvGuard g_off("MESHIOPLUSPLUS_NO_EXTERNAL_MESHERS", "1");
    EmEnvGuard g_tg("MESHIOPLUSPLUS_TETGEN", fake);
    EXPECT_TRUE(meshioplusplus::providers::make_external_tet_meshers().empty());
    EXPECT_TRUE(meshioplusplus::providers::make_external_quad_remeshers().empty());
}

TEST(ExternalMesher, MockTetgenSuccessWithHolesAndErrorPaths) {
    EmTempDir tmp("em_tetgen");
    const std::string ok_script = tmp.WriteScript("mock_tetgen",
                                                  "cat > in.1.node <<'EOF'\n"
                                                  "5 3 0 0\n"
                                                  "0 0 0 0\n"
                                                  "1 1 0 0\n"
                                                  "2 0 1 0\n"
                                                  "3 0 0 1\n"
                                                  "4 0.25 0.25 0.25\n"
                                                  "EOF\n"
                                                  "cat > in.1.ele <<'EOF'\n"
                                                  "1 4 0\n"
                                                  "0 0 1 2 3\n"
                                                  "EOF\n");
    EmEnvGuard g_off("MESHIOPLUSPLUS_NO_EXTERNAL_MESHERS", "0");
    EmEnvGuard g_tg("MESHIOPLUSPLUS_TETGEN", ok_script);
    EmEnvGuard g_gm("MESHIOPLUSPLUS_GMSH", "");
    EmEnvGuard g_ft("MESHIOPLUSPLUS_FTETWILD", "");

    const auto providers = meshioplusplus::providers::make_external_tet_meshers();
    const meshioplusplus::TetMesherProvider* p_tg = nullptr;
    for (const auto& r_p : providers)
        if (r_p.mName == "tetgen")
            p_tg = &r_p;
    ASSERT_NE(p_tg, nullptr);

    const TetMeshRegionResult res = p_tg->mMeshRegion(em_test_single_tet_request());
    EXPECT_EQ(res.mPoints.size(), 5u);
    EXPECT_EQ(res.mTets.size(), 1u);
    EXPECT_TRUE(res.mBoundaryPreserved);

    // Hollow request exercises em_shells() hole injection into .smesh.
    const TetMeshRegionResult hollow_res = p_tg->mMeshRegion(em_test_hollow_request());
    EXPECT_EQ(hollow_res.mTets.size(), 1u);
    EXPECT_FALSE(hollow_res.mBoundaryPreserved);

    // Failing exit code surfaces log tail.
    const std::string fail_script =
        tmp.WriteScript("fail_tetgen", "echo 'synthetic tetgen failure' >&2\nexit 3");
    EmEnvGuard g_tg_fail("MESHIOPLUSPLUS_TETGEN", fail_script);
    const auto fail_providers = meshioplusplus::providers::make_external_tet_meshers();
    for (const auto& r_p : fail_providers)
        if (r_p.mName == "tetgen")
            EXPECT_THROW(r_p.mMeshRegion(em_test_single_tet_request()), std::runtime_error);
}

TEST(ExternalMesher, MockGmshAndFtetwildSuccessAndEmptyTetErrors) {
    EmTempDir tmp("em_gmsh_ftw");
    const std::string gmsh_ok = tmp.WriteScript(
        "mock_gmsh",
        std::string("out=''\n"
                    "while [ $# -gt 0 ]; do\n"
                    "  if [ \"$1\" = '-o' ]; then out=\"$2\"; shift 2; else shift; fi\n"
                    "done\n"
                    "cat > \"$out\" <<'EOF'\n") +
            kMsh22SingleTet + "EOF\n");
    const std::string ftw_ok = tmp.WriteScript(
        "mock_ftw", std::string("out=''\n"
                                "while [ $# -gt 0 ]; do\n"
                                "  if [ \"$1\" = '-o' ]; then out=\"$2\"; shift 2; else shift; fi\n"
                                "done\n"
                                "cat > \"$out\" <<'EOF'\n") +
                        kMsh22SingleTet + "EOF\n");

    EmEnvGuard g_off("MESHIOPLUSPLUS_NO_EXTERNAL_MESHERS", "0");
    EmEnvGuard g_gm("MESHIOPLUSPLUS_GMSH", gmsh_ok);
    EmEnvGuard g_ft("MESHIOPLUSPLUS_FTETWILD", ftw_ok);

    const auto providers = meshioplusplus::providers::make_external_tet_meshers();
    const meshioplusplus::TetMesherProvider* p_gm = nullptr;
    const meshioplusplus::TetMesherProvider* p_ft = nullptr;
    for (const auto& r_p : providers) {
        if (r_p.mName == "gmsh")
            p_gm = &r_p;
        if (r_p.mName == "ftetwild")
            p_ft = &r_p;
    }
    ASSERT_NE(p_gm, nullptr);
    ASSERT_NE(p_ft, nullptr);

    const TetMeshRegionResult gm_res = p_gm->mMeshRegion(em_test_single_tet_request());
    EXPECT_EQ(gm_res.mTets.size(), 1u);
    EXPECT_TRUE(gm_res.mBoundaryPreserved);

    auto ftw_req = em_test_single_tet_request();
    ftw_req.mOptimizeSteps = 0;
    const TetMeshRegionResult ft_res = p_ft->mMeshRegion(ftw_req);
    EXPECT_EQ(ft_res.mTets.size(), 1u);
    EXPECT_FALSE(ft_res.mBoundaryPreserved);

    // Empty tet output triggers named error for both gmsh and ftetwild.
    const std::string empty_out = tmp.WriteScript(
        "mock_empty",
        std::string("out=''\n"
                    "while [ $# -gt 0 ]; do\n"
                    "  if [ \"$1\" = '-o' ]; then out=\"$2\"; shift 2; else shift; fi\n"
                    "done\n"
                    "cat > \"$out\" <<'EOF'\n") +
            kMsh22Empty + "EOF\n");
    EmEnvGuard g_gm2("MESHIOPLUSPLUS_GMSH", empty_out);
    EmEnvGuard g_ft2("MESHIOPLUSPLUS_FTETWILD", empty_out);
    for (const auto& r_p : meshioplusplus::providers::make_external_tet_meshers())
        if (r_p.mName == "gmsh" || r_p.mName == "ftetwild")
            EXPECT_THROW(r_p.mMeshRegion(em_test_single_tet_request()), std::runtime_error);
}

TEST(ExternalMesher, MockQuadRemeshersAndObjErrorBranches) {
    EmTempDir tmp("em_quad");
    const std::string obj_ok =
        tmp.WriteScript("mock_quad",
                        "out=''\n"
                        "while [ $# -gt 0 ]; do\n"
                        "  if [ \"$1\" = '-o' ]; then out=\"$2\"; shift 2; else shift; fi\n"
                        "done\n"
                        "cat > \"$out\" <<'EOF'\n"
                        "# comment\n"
                        "v 0 0 0\n"
                        "v 1 0 0\n"
                        "v 1 1 0\n"
                        "v 0 1 0\n"
                        "f 1/1/1 2/2 3 -1\n"
                        "f 1 2 3\n"
                        "EOF\n");

    EmEnvGuard g_off("MESHIOPLUSPLUS_NO_EXTERNAL_MESHERS", "0");
    EmEnvGuard g_qf("MESHIOPLUSPLUS_QUADRIFLOW", obj_ok);
    EmEnvGuard g_im("MESHIOPLUSPLUS_INSTANT_MESHES", obj_ok);

    QuadRemeshRequest req;
    req.mPoints = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    req.mTriangles = {{0, 1, 2}, {0, 2, 3}};
    req.mTargetFaces = 2;
    req.mSharp = true;
    req.mBoundary = true;
    req.mAdaptive = true;
    req.mDominant = true;

    const auto providers = meshioplusplus::providers::make_external_quad_remeshers();
    ASSERT_GE(providers.size(), 2u);
    for (const auto& r_p : providers) {
        if (r_p.mName != "quadriflow" && r_p.mName != "instant_meshes")
            continue;
        const QuadRemeshResult res = r_p.mRemesh(req);
        EXPECT_EQ(res.mPoints.size(), 4u);
        EXPECT_EQ(res.mQuads.size(), 1u);
        EXPECT_EQ(res.mTriangles.size(), 1u);
    }

    // Exercise OBJ error branches: missing output, bad vertex, polygon >4, bad face corner, empty.
    const std::vector<std::string> bad_bodies = {
        "exit 0",  // writes no out.obj
        "out=''; while [ $# -gt 0 ]; do if [ \"$1\" = '-o' ]; then out=\"$2\"; shift 2; else "
        "shift; fi; done; printf 'v 1 2\\n' > \"$out\"",
        "out=''; while [ $# -gt 0 ]; do if [ \"$1\" = '-o' ]; then out=\"$2\"; shift 2; else "
        "shift; fi; done; printf 'v 0 0 0\\nv 1 0 0\\nv 1 1 0\\nv 0 1 0\\nv 0.5 0.5 1\\nf 1 2 3 4 "
        "5\\n' > \"$out\"",
        "out=''; while [ $# -gt 0 ]; do if [ \"$1\" = '-o' ]; then out=\"$2\"; shift 2; else "
        "shift; fi; done; printf 'v 0 0 0\\nf 0 1 2\\n' > \"$out\"",
        "out=''; while [ $# -gt 0 ]; do if [ \"$1\" = '-o' ]; then out=\"$2\"; shift 2; else "
        "shift; fi; done; printf 'v 0 0 0\\nv 1 0 0\\nv 0 1 0\\nf 1 2 9\\n' > \"$out\"",
        "out=''; while [ $# -gt 0 ]; do if [ \"$1\" = '-o' ]; then out=\"$2\"; shift 2; else "
        "shift; fi; done; printf 'v 0 0 0\\n' > \"$out\"",
    };
    for (std::size_t i = 0; i < bad_bodies.size(); ++i) {
        const std::string bad = tmp.WriteScript("bad_qf_" + std::to_string(i), bad_bodies[i]);
        EmEnvGuard g_bad("MESHIOPLUSPLUS_QUADRIFLOW", bad);
        for (const auto& r_p : meshioplusplus::providers::make_external_quad_remeshers())
            if (r_p.mName == "quadriflow")
                EXPECT_THROW(r_p.mRemesh(req), std::runtime_error);
    }
}

}  // namespace

#endif  // !defined(_WIN32) && !defined(__EMSCRIPTEN__)
