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

// The runtime provider registry. See providers.hpp.

// System includes
#include <mutex>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/providers.hpp"

namespace meshioplusplus {

namespace detail {

#ifdef MESHIOPLUSPLUS_HAS_OCCT
/// The OpenCASCADE provider this build compiled in (cad_read.cpp).
CadKernelProvider builtin_occt_provider();
#endif

}  // namespace detail

namespace {

struct ProvidersState {
    std::mutex mMutex;
    std::shared_ptr<const CadKernelProvider> mCad;
    // Every provider ever registered: a version string handed out as a
    // `const char*` must outlive a later replacement. A handful at most.
    std::vector<std::shared_ptr<const CadKernelProvider>> mRetained;
    bool mBuiltinDone = false;
};

ProvidersState& providers_state() {
    // Constructed on first use and never destroyed, so a query during static
    // destruction stays safe.
    static ProvidersState* const s_state = new ProvidersState;
    return *s_state;
}

}  // namespace

void register_cad_kernel(CadKernelProvider Provider) {
    ProvidersState& r_state = providers_state();
    auto p_new = std::make_shared<const CadKernelProvider>(std::move(Provider));
    const std::lock_guard<std::mutex> lock(r_state.mMutex);
    r_state.mRetained.push_back(p_new);
    r_state.mCad = std::move(p_new);
}

void unregister_cad_kernel() noexcept {
    ProvidersState& r_state = providers_state();
    const std::lock_guard<std::mutex> lock(r_state.mMutex);
    r_state.mBuiltinDone = true;
    r_state.mCad.reset();
}

void register_builtin_providers() noexcept {
    ProvidersState& r_state = providers_state();
    const std::lock_guard<std::mutex> lock(r_state.mMutex);
    if (r_state.mBuiltinDone)
        return;
    r_state.mBuiltinDone = true;
#ifdef MESHIOPLUSPLUS_HAS_OCCT
    if (!r_state.mCad) {
        try {
            r_state.mCad =
                std::make_shared<const CadKernelProvider>(detail::builtin_occt_provider());
            r_state.mRetained.push_back(r_state.mCad);
        } catch (...) {
        }
    }
#endif
}

std::shared_ptr<const CadKernelProvider> cad_kernel() noexcept {
    register_builtin_providers();
    ProvidersState& r_state = providers_state();
    const std::lock_guard<std::mutex> lock(r_state.mMutex);
    return r_state.mCad;
}

}  // namespace meshioplusplus
