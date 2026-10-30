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
 * @file test_plugins.cpp
 * @brief Provider plugins (providers.hpp, "plugins"): a shared library loaded at
 *        run time registers a provider through its two entry points, resolving
 *        the core's functions from this executable; one built for another ABI
 *        version is refused, and a plugin never replaces a provider that is
 *        already registered.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <string>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/providers.hpp"

#if defined(MESHIOPLUSPLUS_TEST_PLUGIN)

namespace {

using namespace meshioplusplus;

/// Restores whatever CAD kernel was registered when the test started.
class PluginKernelGuard {
public:
    PluginKernelGuard() : mSaved(cad_kernel()) {}
    ~PluginKernelGuard() {
        unregister_cad_kernel();
        if (mSaved)
            register_cad_kernel(*mSaved);
    }

private:
    std::shared_ptr<const CadKernelProvider> mSaved;
};

}  // namespace

TEST(Plugins, ALoadedPluginFillsAnEmptySlot) {
    PluginKernelGuard guard;
    unregister_cad_kernel();
    std::string error;
    ASSERT_TRUE(load_plugin(MESHIOPLUSPLUS_TEST_PLUGIN, &error)) << error;
    ASSERT_NE(cad_kernel(), nullptr);
    EXPECT_EQ(cad_kernel()->mName, "plugin-fake");
    EXPECT_TRUE(cad_read_has_occt());
    EXPECT_TRUE(fit_cad_has_occt());  // what mio_has_occt() reports
    EXPECT_STREQ(fit_cad_occt_version(), "7.7.7");
    // The plugin's code runs: its reader throws its own message.
    try {
        read_step("any.step");
        FAIL() << "expected the plugin's ReadError";
    } catch (const ReadError& rErr) {
        EXPECT_NE(std::string(rErr.what()).find("from the plugin"), std::string::npos);
    }
}

TEST(Plugins, AnotherAbiVersionIsRefused) {
    PluginKernelGuard guard;
    unregister_cad_kernel();
    std::string error;
    EXPECT_FALSE(load_plugin(MESHIOPLUSPLUS_TEST_PLUGIN_BAD_ABI, &error));
    EXPECT_NE(error.find("ABI version 1"), std::string::npos) << error;
    EXPECT_EQ(cad_kernel(), nullptr);
}

TEST(Plugins, NeverReplacesARegisteredProvider) {
    PluginKernelGuard guard;
    CadKernelProvider own;
    own.mName = "already-here";
    own.mVersion = "1";
    register_cad_kernel(std::move(own));
    EXPECT_TRUE(load_plugin(MESHIOPLUSPLUS_TEST_PLUGIN));
    EXPECT_EQ(cad_kernel()->mName, "already-here");
}

TEST(Plugins, MissingFilesAndLibrariesWithoutEntryPointsAreReported) {
    std::string error;
    EXPECT_FALSE(load_plugin("/no/such/meshioplusplus_plugin_x.so", &error));
    EXPECT_FALSE(error.empty());
    // The test executable exports no plugin entry points.
    error.clear();
    EXPECT_FALSE(load_plugin(MESHIOPLUSPLUS_TEST_PLUGIN_BAD_ABI "x", &error));
    EXPECT_FALSE(error.empty());
}

TEST(Plugins, ADirectoryAndAPathListAreSearched) {
    PluginKernelGuard guard;
    unregister_cad_kernel();
    EXPECT_EQ(load_plugins(MESHIOPLUSPLUS_TEST_PLUGIN_DIR), 1);
    EXPECT_EQ(cad_kernel()->mName, "plugin-fake");

    unregister_cad_kernel();
    const std::string list = std::string("/no/such/dir:") + MESHIOPLUSPLUS_TEST_PLUGIN_DIR + ":";
    EXPECT_EQ(load_plugins(list), 1);
    EXPECT_EQ(load_plugins(""), 0);
}

#else

TEST(Plugins, NotBuiltOnThisPlatform) {
    GTEST_SKIP() << "plugins need a platform where a library resolves its host's symbols";
}

#endif
