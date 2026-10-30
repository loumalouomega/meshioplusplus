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
 * @file test_providers.cpp
 * @brief The runtime provider registry: a fake CAD kernel is reached by
 *        `read_step`, `read_iges`, `fit_cad` and the capability queries, and
 *        removing it restores the named `Unsupported`.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <string>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/iges.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/providers.hpp"

namespace {

using namespace meshioplusplus;

/// Restores whatever kernel was registered when the test started.
class KernelGuard {
public:
    KernelGuard() : mSaved(cad_kernel()) {}
    ~KernelGuard() {
        unregister_cad_kernel();
        if (mSaved)
            register_cad_kernel(*mSaved);
    }

private:
    std::shared_ptr<const CadKernelProvider> mSaved;
};

CadKernelProvider fake_kernel(int& rCalls) {
    CadKernelProvider provider;
    provider.mName = "fake";
    provider.mVersion = "9.9.9";
    provider.mReadStep = [&rCalls](const std::string& rPath, BrepInfo&, const CadReadOptions&) {
        ++rCalls;
        (void)rPath;
        return mt::tri_mesh();
    };
    provider.mReadIges = [](const std::string&, BrepInfo&, const CadReadOptions&) -> Mesh {
        throw ReadError("fake iges");
    };
    provider.mFitCad = [&rCalls](const Mesh& rMesh, const std::string&, const FitCadOptions&) {
        ++rCalls;
        FitCadResult result;
        result.mMesh = rMesh;
        return result;
    };
    return provider;
}

}  // namespace

TEST(Providers, RegisteredKernelIsUsed) {
    KernelGuard guard;
    int calls = 0;
    register_cad_kernel(fake_kernel(calls));
    EXPECT_TRUE(cad_read_has_occt());
    EXPECT_TRUE(fit_cad_has_occt());
    EXPECT_STREQ(fit_cad_occt_version(), "9.9.9");
    EXPECT_EQ(cad_kernel()->mName, "fake");

    const Mesh mesh = read_step("any.step");
    EXPECT_EQ(calls, 1);
    EXPECT_EQ(mesh.NumPoints(), mt::tri_mesh().NumPoints());
    EXPECT_THROW(read_iges("any.iges"), ReadError);

    FitCadOptions options;
    options.mBackend = FitCadBackend::Occt;
    const FitCadResult result = fit_cad(mesh, "", options);
    EXPECT_EQ(calls, 2);
    EXPECT_EQ(result.mMesh.NumPoints(), mesh.NumPoints());
}

TEST(Providers, NoKernelNamesTheAddOn) {
    KernelGuard guard;
    unregister_cad_kernel();
    EXPECT_FALSE(cad_read_has_occt());
    EXPECT_FALSE(fit_cad_has_occt());
    EXPECT_STREQ(fit_cad_occt_version(), "");
    for (bool iges : {false, true}) {
        try {
            if (iges)
                read_iges("any.iges");
            else
                read_step("any.step");
            FAIL() << "expected Unsupported";
        } catch (const Unsupported& rErr) {
            EXPECT_NE(std::string(rErr.what()).find("meshioplusplus-occt"), std::string::npos);
        }
    }
    FitCadOptions options;
    options.mBackend = FitCadBackend::Occt;
    EXPECT_THROW(fit_cad(Mesh(), "", options), Unsupported);
}

TEST(Providers, ReplacingKeepsOldVersionStringValid) {
    KernelGuard guard;
    int calls = 0;
    register_cad_kernel(fake_kernel(calls));
    const char* old_version = fit_cad_occt_version();
    CadKernelProvider other = fake_kernel(calls);
    other.mVersion = "1.0";
    register_cad_kernel(std::move(other));
    EXPECT_STREQ(old_version, "9.9.9");
    EXPECT_STREQ(fit_cad_occt_version(), "1.0");
}
