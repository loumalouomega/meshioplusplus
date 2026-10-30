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
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <mutex>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/abi_version.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/providers.hpp"

// Project includes (private, not installed)
#include "detail/dynlib.hpp"

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

// Set while this thread runs a plugin's registration function: a plugin fills an
// empty slot but never replaces a provider that is already there.
thread_local bool t_loading_plugin = false;

template <class T>
void slot_register(ProviderSlot<T> ProvidersState::*pSlot, T Provider) {
    ProvidersState& r_state = providers_state();
    auto p_new = std::make_shared<const T>(std::move(Provider));
    const std::lock_guard<std::mutex> lock(r_state.mMutex);
    if (t_loading_plugin && (r_state.*pSlot).mCurrent)
        return;
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

namespace {

/// The scan of MESHIOPLUSPLUS_PLUGIN_PATH runs once. Recursive, so that a plugin
/// asking for a provider while it registers does not wait on itself; the flag
/// makes that nested call return at once.
std::recursive_mutex& plugin_scan_mutex() {
    static std::recursive_mutex* const s_mutex = new std::recursive_mutex;
    return *s_mutex;
}
bool g_plugins_scan_started = false;        // guarded by the scan mutex
std::atomic<bool> g_plugins_scan_done{false};  // the lock-free fast path

using PluginAbiFn = int (*)();
using PluginRegisterFn = void (*)();

bool is_plugin_file(const std::filesystem::path& rPath) {
    const std::string name = rPath.filename().string();
    return name.rfind("meshioplusplus_plugin_", 0) == 0 ||
           name.rfind("libmeshioplusplus_plugin_", 0) == 0;
}

}  // namespace

bool load_plugin(const std::string& rPath, std::string* pError) noexcept {
    const auto fail = [&](const std::string& rWhy) {
        if (pError)
            *pError = rWhy;
        log::warn("plugin '{}' not loaded: {}", rPath, rWhy);
        return false;
    };
    try {
        if (!detail::dynlib_available())
            return fail("loading shared libraries is not available in this build");
        std::string error;
        detail::dynlib_export_host_symbols();
        void* p_lib = detail::dynlib_open(rPath, error);
        if (!p_lib)
            return fail(error);
        const auto p_abi = reinterpret_cast<PluginAbiFn>(
            detail::dynlib_symbol(p_lib, "meshioplusplus_plugin_abi_version"));
        const auto p_register = reinterpret_cast<PluginRegisterFn>(
            detail::dynlib_symbol(p_lib, "meshioplusplus_register_plugin"));
        if (!p_abi || !p_register) {
            detail::dynlib_close(p_lib);
            return fail("it exports no meshioplusplus_plugin_abi_version / "
                        "meshioplusplus_register_plugin");
        }
        const int abi = p_abi();
        if (abi != MESHIOPLUSPLUS_ABI_VERSION) {
            detail::dynlib_close(p_lib);
            return fail("it was built for ABI version " + std::to_string(abi) +
                        ", this library is " + std::to_string(MESHIOPLUSPLUS_ABI_VERSION));
        }
        // The library stays loaded: its providers hold code from it.
        t_loading_plugin = true;
        try {
            p_register();
        } catch (...) {
            t_loading_plugin = false;
            return fail("its registration threw");
        }
        t_loading_plugin = false;
        return true;
    } catch (...) {
        return false;
    }
}

int load_plugins(const std::string& rPaths) noexcept {
    int loaded = 0;
    try {
#if defined(_WIN32)
        constexpr char kSeparator = ';';
#else
        constexpr char kSeparator = ':';
#endif
        std::size_t start = 0;
        while (start <= rPaths.size()) {
            std::size_t end = rPaths.find(kSeparator, start);
            if (end == std::string::npos)
                end = rPaths.size();
            const std::string entry = rPaths.substr(start, end - start);
            start = end + 1;
            if (entry.empty())
                continue;
            std::error_code ec;
            const std::filesystem::path path(entry);
            if (std::filesystem::is_directory(path, ec)) {
                std::vector<std::string> files;
                for (const auto& r_item : std::filesystem::directory_iterator(path, ec))
                    if (r_item.is_regular_file(ec) && is_plugin_file(r_item.path()) &&
                        r_item.path().extension() == detail::dynlib_extension())
                        files.push_back(r_item.path().string());
                std::sort(files.begin(), files.end());  // a stable order
                for (const std::string& r_file : files)
                    if (load_plugin(r_file))
                        ++loaded;
            } else if (load_plugin(entry)) {
                ++loaded;
            }
        }
    } catch (...) {
    }
    return loaded;
}

static void register_compiled_in_providers() noexcept {
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

void register_builtin_providers() noexcept {
    register_compiled_in_providers();
    // Then the plugins, once and outside the registry lock (a plugin registers
    // through it). Concurrent callers wait here until the scan is done, so none
    // sees a slot a plugin is about to fill; a nested call on the scanning
    // thread passes straight through.
    if (g_plugins_scan_done.load(std::memory_order_acquire))
        return;
    const std::lock_guard<std::recursive_mutex> scan_lock(plugin_scan_mutex());
    if (g_plugins_scan_started)
        return;  // done meanwhile, or this thread's own plugin is asking
    g_plugins_scan_started = true;
    const char* p_paths = std::getenv("MESHIOPLUSPLUS_PLUGIN_PATH");
    if (p_paths && *p_paths)
        load_plugins(p_paths);
    g_plugins_scan_done.store(true, std::memory_order_release);
}

}  // namespace meshioplusplus
