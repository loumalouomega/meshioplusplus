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

// A plugin for the tests: registers a fake CAD kernel through the two entry points
// providers.hpp describes, without linking the core (its calls resolve from the
// test executable, which is linked -rdynamic). `FAKE_PLUGIN_ABI` overrides the ABI
// version it reports, to build the one the loader must refuse.

#include "meshioplusplus/abi_version.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/providers.hpp"

#ifdef _WIN32
#define FAKE_EXPORT __declspec(dllexport)
#else
#define FAKE_EXPORT __attribute__((visibility("default")))
#endif

#ifndef FAKE_PLUGIN_ABI
#define FAKE_PLUGIN_ABI MESHIOPLUSPLUS_ABI_VERSION
#endif

extern "C" {

FAKE_EXPORT int meshioplusplus_plugin_abi_version(void) {
    return FAKE_PLUGIN_ABI;
}

FAKE_EXPORT void meshioplusplus_register_plugin(void) {
    meshioplusplus::CadKernelProvider provider;
    provider.mName = "plugin-fake";
    provider.mVersion = "7.7.7";
    provider.mReadStep = [](const std::string&, meshioplusplus::BrepInfo&,
                            const meshioplusplus::CadReadOptions&) -> meshioplusplus::Mesh {
        throw meshioplusplus::ReadError("from the plugin");
    };
    meshioplusplus::register_cad_kernel(std::move(provider));
}

}  // extern "C"
