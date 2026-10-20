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
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/fit_cad.hpp"
#include "fit_cad_fixtures.hpp"
#include "golden_fixtures.hpp"

using namespace meshioplusplus;

namespace {

struct FcgRow {
    const char* mOcct;  // "native" for the native backend
    const char* mFixture;
    std::uint64_t mStep;
    std::uint64_t mMesh;
    std::uint64_t mReport;
};

// Native backend (every build), recorded with MESHIOPLUSPLUS_RECORD_GOLDENS=1.
const std::vector<FcgRow> kFcgNative = {
    // clang-format off
    // clang-format on
};

// OpenCASCADE backend, keyed by fit_cad_occt_version().
const std::vector<FcgRow> kFcgOcct = {
    // clang-format off
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
// header names the file).
FcgDigest fcg_run(const Mesh& rSkin, FitCadBackend Backend) {
    FitCadOptions o;
    o.mBackend = Backend;
    const std::string path = ::testing::TempDir() + "fit_cad_golden.step";
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
#ifdef GOLDEN_PINNED
    fcg_check("native", FitCadBackend::Native, kFcgNative);
#else
    GTEST_SKIP() << "golden digests are pinned on x86-64 GCC/Clang only";
#endif
}

TEST(FitCadGoldens, OcctMatchesForThisVersion) {
    if (!fit_cad_has_occt()) {
        EXPECT_STREQ(fit_cad_occt_version(), "");
        GTEST_SKIP() << "built without OpenCASCADE";
    }
    EXPECT_NE(std::string(fit_cad_occt_version()), "");
#ifdef GOLDEN_PINNED
    fcg_check(fit_cad_occt_version(), FitCadBackend::Occt, kFcgOcct);
#else
    GTEST_SKIP() << "golden digests are pinned on x86-64 GCC/Clang only";
#endif
}
