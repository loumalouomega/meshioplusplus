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
#include <type_traits>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/providers.hpp"

namespace meshioplusplus {

namespace providers {

#ifdef MESHIOPLUSPLUS_HAS_OCCT
/// The OpenCASCADE provider this build compiled in (providers/occt/).
CadKernelProvider make_occt_provider();
#endif
#ifdef MESHIOPLUSPLUS_HAS_NETGEN
/// The Netgen provider this build compiled in (providers/netgen/).
TetMesherProvider make_netgen_provider();
#endif
#ifdef MESHIOPLUSPLUS_HAS_MMG
/// The MMG provider this build compiled in (providers/mmg/).
RemesherProvider make_mmg_provider();
#endif

}  // namespace providers

namespace {

/// One kind of provider: the registered one, whether its builtin has had its
/// turn (also set by an explicit unregister, which must stick), and every one
/// ever registered, because a version string handed out as a `const char*`
/// must outlive a later replacement. A handful at most.
template <class T>
struct ProviderSlot {
    std::shared_ptr<const T> mCurrent;
    std::vector<std::shared_ptr<const T>> mRetained;
    bool mBuiltinDone = false;

    void Set(std::shared_ptr<const T> pNew) {
        mRetained.push_back(pNew);
        mCurrent = std::move(pNew);
    }
};

struct ProvidersState {
    std::mutex mMutex;
    ProviderSlot<CadKernelProvider> mCad;
    ProviderSlot<TetMesherProvider> mTet;
    ProviderSlot<RemesherProvider> mRemesher;
};

ProvidersState& providers_state() {
    // Constructed on first use and never destroyed, so a query during static
    // destruction stays safe.
    static ProvidersState* const s_state = new ProvidersState;
    return *s_state;
}

template <class T>
void slot_register(ProviderSlot<T> ProvidersState::*pSlot, T Provider) {
    ProvidersState& r_state = providers_state();
    auto p_new = std::make_shared<const T>(std::move(Provider));
    const std::lock_guard<std::mutex> lock(r_state.mMutex);
    (r_state.*pSlot).Set(std::move(p_new));
}

template <class T>
void slot_unregister(ProviderSlot<T> ProvidersState::*pSlot) noexcept {
    ProvidersState& r_state = providers_state();
    const std::lock_guard<std::mutex> lock(r_state.mMutex);
    (r_state.*pSlot).mBuiltinDone = true;
    (r_state.*pSlot).mCurrent.reset();
}

template <class T>
std::shared_ptr<const T> slot_get(ProviderSlot<T> ProvidersState::*pSlot) noexcept {
    register_builtin_providers();
    ProvidersState& r_state = providers_state();
    const std::lock_guard<std::mutex> lock(r_state.mMutex);
    return (r_state.*pSlot).mCurrent;
}

}  // namespace

void register_cad_kernel(CadKernelProvider Provider) {
    slot_register(&ProvidersState::mCad, std::move(Provider));
}
void unregister_cad_kernel() noexcept {
    slot_unregister(&ProvidersState::mCad);
}
std::shared_ptr<const CadKernelProvider> cad_kernel() noexcept {
    return slot_get(&ProvidersState::mCad);
}

void register_tet_mesher(TetMesherProvider Provider) {
    slot_register(&ProvidersState::mTet, std::move(Provider));
}
void unregister_tet_mesher() noexcept {
    slot_unregister(&ProvidersState::mTet);
}
std::shared_ptr<const TetMesherProvider> tet_mesher() noexcept {
    return slot_get(&ProvidersState::mTet);
}

void register_remesher(RemesherProvider Provider) {
    slot_register(&ProvidersState::mRemesher, std::move(Provider));
}
void unregister_remesher() noexcept {
    slot_unregister(&ProvidersState::mRemesher);
}
std::shared_ptr<const RemesherProvider> remesher() noexcept {
    return slot_get(&ProvidersState::mRemesher);
}

void register_builtin_providers() noexcept {
    ProvidersState& r_state = providers_state();
    const std::lock_guard<std::mutex> lock(r_state.mMutex);
    // A slot's builtin goes in once, and never over a provider already there.
    const auto install = [](auto& rSlot, auto&& Make) {
        if (rSlot.mBuiltinDone)
            return;
        rSlot.mBuiltinDone = true;
        if (rSlot.mCurrent)
            return;
        try {
            using Provider = std::decay_t<decltype(Make())>;
            rSlot.Set(std::make_shared<const Provider>(Make()));
        } catch (...) {
        }
    };
    (void)install;
#ifdef MESHIOPLUSPLUS_HAS_OCCT
    install(r_state.mCad, [] { return providers::make_occt_provider(); });
#else
    r_state.mCad.mBuiltinDone = true;
#endif
#ifdef MESHIOPLUSPLUS_HAS_NETGEN
    install(r_state.mTet, [] { return providers::make_netgen_provider(); });
#else
    r_state.mTet.mBuiltinDone = true;
#endif
#ifdef MESHIOPLUSPLUS_HAS_MMG
    install(r_state.mRemesher, [] { return providers::make_mmg_provider(); });
#else
    r_state.mRemesher.mBuiltinDone = true;
#endif
}

}  // namespace meshioplusplus
