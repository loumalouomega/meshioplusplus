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
// fit_cad golden digests: the STEP bytes, the returned mesh and the report of
// every fixture in fit_cad_fixtures.hpp, pinned so a change that is meant to
// leave the default output alone (every gated skin-to-CAD stage) provably
// does. The native backend's digests hold on every pinned platform; the
// OpenCASCADE backend's depend on the kernel's version too, so its table is
// keyed by `fit_cad_occt_version()` and an unknown version skips by name.
//
// MESHIOPLUSPLUS_RECORD_GOLDENS=1 prints every row as source to paste here.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/fit_cad.hpp"
#include "fit_cad_fixtures.hpp"
#include "golden_fixtures.hpp"

using namespace meshioplusplus;

// The fixtures call libm (sin, cos) and the fits much more, so the digests are
// pinned where golden_fixtures.hpp pins and the C library is glibc too.
#if defined(GOLDEN_PINNED) && defined(__linux__) && defined(__GLIBC__)
#define FCG_PINNED 1
#endif

namespace {

struct FcgRow {
    const char* mOcct;  // "native" for the native backend
    const char* mFixture;
    std::uint64_t mStep;
    std::uint64_t mMesh;
    std::uint64_t mReport;
};

// The fixtures themselves, as golden::ft_mesh digests (gcc, glibc). A fixture
// built from libm calls can come out in other bits under another compiler
// (clang rounds the rounded box's corners differently): its rows are then
// skipped with a note, since they pin fit_cad and not the generator.
const std::vector<std::pair<const char*, std::uint64_t>> kFcgInputs = {
    // clang-format off
    {"box", 0x6275968d717718fcull},
    {"capped_cylinder", 0xf0977bfa91cf5b29ull},
    {"turned_part", 0x0f960d340fc274ddull},
    {"filleted_plate", 0xe67eb44a33284280ull},
    {"saddle_box", 0x36476b6cbb93e463ull},
    {"wavy_capped", 0xd9fdf0cace82383bull},
    {"threaded_rod", 0xc1b588ca9b8e1519ull},
    {"rounded_box", 0x2c285c2d8152cb7full},
    {"noisy_turned_part", 0xf4578cd2374eddb2ull},
    // clang-format on
};

// Native backend (every build), recorded with MESHIOPLUSPLUS_RECORD_GOLDENS=1.
const std::vector<FcgRow> kFcgNative = {
    // clang-format off
    {"native", "box", 0x43696b9756cf8851ull, 0x041bc38bff060516ull, 0x4d89933e194a3081ull},
    {"native", "capped_cylinder", 0xb1314c08c3b7a82aull, 0xc85679404821d2c9ull, 0xa1416e12e3b60778ull},
    {"native", "turned_part", 0xf031f071b85fdf14ull, 0x5b6613001c187538ull, 0x6394f6a1719495f7ull},
    {"native", "filleted_plate", 0xbee00b6d57b28b8full, 0xa1fbfd1fb6894ee9ull, 0x84bca80af59d7a4aull},
    {"native", "saddle_box", 0x7b0385e0b64e69a9ull, 0x28b7caf6375d5072ull, 0x6bf83615b0b506acull},
    {"native", "wavy_capped", 0xb6b35691854fc2bcull, 0x5f9aa9b2e7379ba8ull, 0x1b2f18e6d0fdbc57ull},
    {"native", "threaded_rod", 0x23b963888e76fff8ull, 0x7180517275d1fd5bull, 0xc2bfba07afac00b2ull},
    {"native", "rounded_box", 0xd1eacc08a6a62994ull, 0xb9bb8a6bf508af15ull, 0xef43d32d4d262246ull},
    {"native", "noisy_turned_part", 0x6067531b52fe88eeull, 0xe4ae427a6fc0e40bull, 0xdf51b734edfe0278ull},
    // clang-format on
};

// OpenCASCADE backend, keyed by fit_cad_occt_version().
const std::vector<FcgRow> kFcgOcct = {
    // clang-format off
    {"7.6.3", "box", 0xed8bdb22b4ce64f1ull, 0xe073b156cc1154ecull, 0x73c6e572e1989b8full},
    {"7.6.3", "capped_cylinder", 0xe33194378779ddb7ull, 0xb9e9056abf4d1125ull, 0xacf88b6ee280f77cull},
    {"7.6.3", "turned_part", 0x79562744dbecd568ull, 0xbfb8c1cfaffd675eull, 0xe44910bdd1dc1d25ull},
    {"7.6.3", "filleted_plate", 0xa2295aae6f4930e5ull, 0x61ccc19d6b87f045ull, 0x31113b57e28ad516ull},
    {"7.6.3", "saddle_box", 0x6eafb487a5e02434ull, 0xb70e5093a40e3605ull, 0x9e4088243820e858ull},
    {"7.6.3", "wavy_capped", 0x5929da73cd7f7e69ull, 0x95dd65e32746a786ull, 0x9868991cd0ab0304ull},
    {"7.6.3", "threaded_rod", 0x22435a2b7dacfebcull, 0xe1bdaeea9124bbc7ull, 0xb9cd0dfb82c69aeeull},
    {"7.6.3", "rounded_box", 0x8d27c21191658323ull, 0x99d00ac337bfa4dcull, 0x9ea826ae8fb9ae2full},
    {"7.6.3", "noisy_turned_part", 0x18df6a7b3bd131b5ull, 0x6bac3ed931f55ea7ull, 0x61cc8c40a6897004ull},
    {"7.9.3", "box", 0xc83eed5f7e37ebc6ull, 0xe073b156cc1154ecull, 0x73c6e572e1989b8full},
    {"7.9.3", "capped_cylinder", 0x8abacb628e6e6e76ull, 0xb9e9056abf4d1125ull, 0xacf88b6ee280f77cull},
    {"7.9.3", "turned_part", 0xb693bb93eeceb1a1ull, 0x9ce1738499348461ull, 0xe44910bdd1dc1d25ull},
    {"7.9.3", "filleted_plate", 0xe3dd58b0603d52dfull, 0xf07fbce0db843999ull, 0x31113b57e28ad516ull},
    {"7.9.3", "saddle_box", 0x599c75f73e419a43ull, 0x1bc93b2174b358beull, 0xe85f3073487172aeull},
    {"7.9.3", "wavy_capped", 0x3f08dd859f6cc7aeull, 0x95dd65e32746a786ull, 0x9868991cd0ab0304ull},
    {"7.9.3", "threaded_rod", 0xa7aff9a64cf50849ull, 0x32b791b5e27de7daull, 0xb9cd0dfb82c69aeeull},
    {"7.9.3", "rounded_box", 0x9efd222ffa06b308ull, 0x9aff1debc503a57cull, 0x9ea826ae8fb9ae2full},
    {"7.9.3", "noisy_turned_part", 0x61d58d03515e9c54ull, 0x2aac06dc4866d59eull, 0x61cc8c40a6897004ull},
    // clang-format on
};

bool fcg_recording() {
    const char* v = std::getenv("MESHIOPLUSPLUS_RECORD_GOLDENS");
    return v && *v && std::strcmp(v, "0") != 0;
}

std::string fcg_read(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

template <class T>
void fcg_bits(golden::MeshDigest& rD, T V) {
    rD.Bytes(&V, sizeof V);
}

// Every scalar of the report, as raw bits, in declaration order.
std::uint64_t fcg_report(const FitCadResult& rR) {
    golden::MeshDigest d;
    fcg_bits(d, rR.mValid);
    for (double v : {rR.mTolerance, rR.mVolume, rR.mSkinVolume, rR.mVolumeRatio, rR.mMaxDeviation})
        fcg_bits(d, v);
    for (std::int64_t v : {rR.mNumFaces,         rR.mNumPlane,          rR.mNumSphere,
                           rR.mNumCylinder,      rR.mNumCone,           rR.mNumTorus,
                           rR.mNumFacetedFaces,  rR.mNumPatches,        rR.mNumDemoted,
                           rR.mNumExactEdges,    rR.mNumApproxEdges,    rR.mNumSolids,
                           rR.mNumBspline,       rR.mNumSplitPatches,   rR.mNumFillet,
                           rR.mNumHelical,       rR.mNumVariableFillet, rR.mNumBlend,
                           rR.mNumRoundedThread, rR.mNumVertexBlend,    rR.mNumChamfer,
                           rR.mNumConicBlend,    rR.mNumTaperedThread,  rR.mNumVariablePitchThread})
        fcg_bits(d, v);
    fcg_bits(d, static_cast<std::int32_t>(rR.mBackend));
    return d.Value();
}

struct FcgDigest {
    std::uint64_t mStep = 0;
    std::uint64_t mMesh = 0;
    std::uint64_t mReport = 0;
};

// The fixture through fit_cad, its STEP written under a fixed basename (the
// header names the file), in a directory of the running test's own, since
// ctest runs these tests in parallel.
FcgDigest fcg_run(const Mesh& rSkin, FitCadBackend Backend) {
    FitCadOptions o;
    o.mBackend = Backend;
    const ::testing::TestInfo* info = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::filesystem::path dir =
        std::filesystem::path(::testing::TempDir()) /
        (std::string("fcg_") + (info ? info->name() : "run"));
    std::filesystem::create_directories(dir);
    const std::string path = (dir / "fit_cad_golden.step").string();
    const FitCadResult r = fit_cad(rSkin, path, o);
    const std::string bytes = fcg_read(path);
    std::remove(path.c_str());
    FcgDigest g;
    golden::MeshDigest s;
    s.Bytes(bytes.data(), bytes.size());
    g.mStep = s.Value();
    g.mMesh = golden::ft_mesh(r.mMesh);
    g.mReport = fcg_report(r);
    return g;
}

void fcg_check(const char* pKey, FitCadBackend Backend, const std::vector<FcgRow>& rRows) {
    const bool record = fcg_recording();
    std::size_t matched = 0;
    for (const fcfix::Fixture& f : fcfix::golden_fixtures()) {
        SCOPED_TRACE(f.mName);
        const std::uint64_t input = golden::ft_mesh(f.mMesh);
        bool same_input = true;
        for (const auto& in : kFcgInputs)
            if (f.mName == in.first)
                same_input = input == in.second;
        if (!same_input && !record) {
            std::printf(
                "fit_cad goldens: fixture %s is built in other bits by this compiler "
                "(0x%016llx); its rows are not compared\n",
                f.mName.c_str(), static_cast<unsigned long long>(input));
            continue;
        }
        const FcgDigest g = fcg_run(f.mMesh, Backend);
        if (record)
            std::printf("    {\"%s\", \"%s\", 0x%016llxull, 0x%016llxull, 0x%016llxull},\n", pKey,
                        f.mName.c_str(), static_cast<unsigned long long>(g.mStep),
                        static_cast<unsigned long long>(g.mMesh),
                        static_cast<unsigned long long>(g.mReport));
        for (const FcgRow& row : rRows) {
            if (std::strcmp(row.mOcct, pKey) != 0 || f.mName != row.mFixture)
                continue;
            ++matched;
            if (record)
                continue;
            EXPECT_EQ(g.mStep, row.mStep) << "STEP bytes moved";
            EXPECT_EQ(g.mMesh, row.mMesh) << "returned mesh moved";
            EXPECT_EQ(g.mReport, row.mReport) << "report moved";
        }
    }
    if (!record && matched == 0)
        GTEST_SKIP() << "no golden rows for " << pKey
                     << "; record them with MESHIOPLUSPLUS_RECORD_GOLDENS=1";
}

}  // namespace

// The same input, run twice in one process, gives the same digests on every
// platform and backend this build carries.
TEST(FitCadGoldens, RepeatRunsAgree) {
    const Mesh skin = fcfix::fitted(analytic::capped_cylinder());
    for (FitCadBackend b : {FitCadBackend::Native, FitCadBackend::Auto}) {
        const FcgDigest x = fcg_run(skin, b), y = fcg_run(skin, b);
        EXPECT_EQ(x.mStep, y.mStep);
        EXPECT_EQ(x.mMesh, y.mMesh);
        EXPECT_EQ(x.mReport, y.mReport);
    }
}

TEST(FitCadGoldens, NativeBackendMatches) {
#ifdef FCG_PINNED
    fcg_check("native", FitCadBackend::Native, kFcgNative);
#else
    GTEST_SKIP() << "fit_cad digests are pinned on x86-64 Linux (glibc) only";
#endif
}

TEST(FitCadGoldens, OcctMatchesForThisVersion) {
    if (!fit_cad_has_occt()) {
        EXPECT_STREQ(fit_cad_occt_version(), "");
        GTEST_SKIP() << "built without OpenCASCADE";
    }
    EXPECT_NE(std::string(fit_cad_occt_version()), "");
#ifdef FCG_PINNED
    fcg_check(fit_cad_occt_version(), FitCadBackend::Occt, kFcgOcct);
#else
    GTEST_SKIP() << "fit_cad digests are pinned on x86-64 Linux (glibc) only";
#endif
}
