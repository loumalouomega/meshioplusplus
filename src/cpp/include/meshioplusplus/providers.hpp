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
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/operations/mmg.hpp"

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

/// One closed region handed to a tetrahedral mesher: its points and its
/// triangles (outward, 0-based indices into `mPoints`).
struct TetMeshRegionRequest {
    std::vector<std::array<double, 3>> mPoints;
    std::vector<std::array<std::int64_t, 3>> mTriangles;
    /// The target edge length; at most 0 lets the mesher choose.
    double mMaxSize = 0.0;
    double mGrading = 0.3;
    int mOptimizeSteps = 0;
};

/// What a tetrahedral mesher returns for one region.
struct TetMeshRegionResult {
    /// The request's points first, as the mesher holds them (so the caller can
    /// check they were kept), then the points it added.
    std::vector<std::array<double, 3>> mPoints;
    /// Tetrahedra, 0-based indices into `mPoints`.
    std::vector<std::array<std::int64_t, 4>> mTets;
    /// Whether the mesher kept the request's triangles as given.
    bool mBoundaryPreserved = true;
};

/// A tetrahedral mesher that keeps the boundary: `tetrahedralize`'s Netgen backend.
struct TetMesherProvider {
    /// The mesher's name ("netgen"), reported as the backend.
    std::string mName;
    std::string mVersion;
    /// Meshes one region; throws `std::runtime_error` when the mesher fails.
    std::function<TetMeshRegionResult(const TetMeshRegionRequest&)> mMeshRegion;
};

/// A remesher: `adapt`, `discretize_levelset`, `mesh_domain` and `move` (MMG).
struct RemesherProvider {
    std::string mName;
    std::string mVersion;
    /// Whether `move` is available (the library was built with ELAS).
    bool mElas = false;
    /// The width of the library's integer indices, 32 or 64.
    int mIndexBits = 32;
    std::function<AdaptResult(const Mesh&, const AdaptOptions&)> mAdapt;
    std::function<DiscretizeLevelsetResult(const Mesh&, const DiscretizeLevelsetOptions&)>
        mDiscretizeLevelset;
    std::function<MeshDomainResult(const Mesh&, const MeshDomainOptions&)> mMeshDomain;
    std::function<AdaptResult(const Mesh&, const MoveOptions&)> mMove;
};

/// Registers `Provider` as the tetrahedral mesher, replacing any earlier one.
MESHIOPLUSPLUS_API void register_tet_mesher(TetMesherProvider Provider);
/// Removes the registered tetrahedral mesher (the builtin is not re-registered).
MESHIOPLUSPLUS_API void unregister_tet_mesher() noexcept;
/// The registered tetrahedral mesher, or null.
MESHIOPLUSPLUS_API std::shared_ptr<const TetMesherProvider> tet_mesher() noexcept;

/// Registers `Provider` as the remesher, replacing any earlier one.
MESHIOPLUSPLUS_API void register_remesher(RemesherProvider Provider);
/// Removes the registered remesher (the builtin is not re-registered).
MESHIOPLUSPLUS_API void unregister_remesher() noexcept;
/// The registered remesher, or null.
MESHIOPLUSPLUS_API std::shared_ptr<const RemesherProvider> remesher() noexcept;

/// Registers `Provider` as the CAD kernel, replacing any earlier one.
MESHIOPLUSPLUS_API void register_cad_kernel(CadKernelProvider Provider);

/// Removes the registered CAD kernel (builtin included; the builtin is not
/// re-registered afterwards). For tests and plugin teardown.
MESHIOPLUSPLUS_API void unregister_cad_kernel() noexcept;

/// The registered CAD kernel, or null when there is none. The returned
/// provider stays valid for as long as the caller holds it, whatever is
/// registered meanwhile.
MESHIOPLUSPLUS_API std::shared_ptr<const CadKernelProvider> cad_kernel() noexcept;

/// Registers the providers this build compiled in (`-DMESHIOPLUSPLUS_WITH_OCCT`,
/// `_WITH_NETGEN`, `_WITH_MMG`).
/// Idempotent, cheap after the first call, and called by every capability
/// query; it never overrides a provider that is already registered.
MESHIOPLUSPLUS_API void register_builtin_providers() noexcept;

}  // namespace meshioplusplus
