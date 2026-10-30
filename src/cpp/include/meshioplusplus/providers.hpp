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
#pragma once

// The runtime provider registry. The CAD kernel (OpenCASCADE) is reached through
// a registered provider rather than a compile-time switch, so a build that does
// not carry it can still gain it at run time: the Python add-on wheel
// `meshioplusplus-occt` registers one when it is imported. A monolithic build
// compiles its provider into the core and registers it with
// `register_builtin_providers()`, which every capability query calls first.
//
// Registration is explicit and lock-protected; there are no static
// initializers, so there is no initialization-order fragility. The core and a
// provider must come from the same release (the add-on wheels pin the base
// wheel's exact version), so the types here need only be stable within one.

// System includes
#include <functional>
#include <memory>
#include <string>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"

namespace meshioplusplus {

/// A CAD kernel: STEP and IGES reading and `fit_cad`'s kernel backend.
struct CadKernelProvider {
    /// The kernel's name ("occt").
    std::string mName;
    /// The kernel's version string (`OCC_VERSION_COMPLETE`).
    std::string mVersion;
    /// `read_step` with a B-rep; throws `ReadError` on a bad file.
    std::function<Mesh(const std::string&, BrepInfo&, const CadReadOptions&)> mReadStep;
    /// `read_iges` with a B-rep.
    std::function<Mesh(const std::string&, BrepInfo&, const CadReadOptions&)> mReadIges;
    /// The kernel backend of `fit_cad`.
    std::function<FitCadResult(const Mesh&, const std::string&, const FitCadOptions&)> mFitCad;
};

/// Registers `Provider` as the CAD kernel, replacing any earlier one.
MESHIOPLUSPLUS_API void register_cad_kernel(CadKernelProvider Provider);

/// Removes the registered CAD kernel (builtin included; the builtin is not
/// re-registered afterwards). For tests and plugin teardown.
MESHIOPLUSPLUS_API void unregister_cad_kernel() noexcept;

/// The registered CAD kernel, or null when there is none. The returned
/// provider stays valid for as long as the caller holds it, whatever is
/// registered meanwhile.
MESHIOPLUSPLUS_API std::shared_ptr<const CadKernelProvider> cad_kernel() noexcept;

/// Registers the providers this build compiled in (`-DMESHIOPLUSPLUS_WITH_OCCT`).
/// Idempotent, cheap after the first call, and called by every capability
/// query; it never overrides a provider that is already registered.
MESHIOPLUSPLUS_API void register_builtin_providers() noexcept;

}  // namespace meshioplusplus
