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
 *        `read_step`, `read_iges`, `fit_cad` and the capability queries, a
 *        fake tetrahedral mesher by `tetrahedralize`, a fake remesher by
 *        `adapt` and the MMG capability queries, and removing each restores
 *        the named `Unsupported`.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <algorithm>
#include <atomic>
#include <stdexcept>
#include <string>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/iges.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/operations/mmg.hpp"
#include "meshioplusplus/operations/tetrahedralize.hpp"
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
        result.mMesh = detail::clone_mesh(rMesh);
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

// --- the tetrahedral mesher ------------------------------------------------------------

namespace {

/// Restores whatever tetrahedral mesher was registered when the test started.
class MesherGuard {
public:
    MesherGuard() : mSaved(tet_mesher()) {}
    ~MesherGuard() {
        unregister_tet_mesher();
        if (mSaved)
            register_tet_mesher(*mSaved);
    }

private:
    std::shared_ptr<const TetMesherProvider> mSaved;
};

/// A tetrahedron's four faces, outward.
Mesh tetra_surface() {
    return mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "triangle",
                         {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}});
}

/// Fills a region with the one tetrahedron of its four points, adding one point.
TetMesherProvider fake_mesher(int& rCalls, bool KeepBoundary = true) {
    TetMesherProvider provider;
    provider.mName = "fake-mesher";
    provider.mVersion = "1.2.3";
    provider.mMeshRegion = [&rCalls, KeepBoundary](const TetMeshRegionRequest& rRequest) {
        ++rCalls;
        TetMeshRegionResult result;
        result.mPoints = rRequest.mPoints;
        result.mPoints.push_back({0.25, 0.25, 0.25});
        result.mTets = {{0, 1, 2, 4}, {0, 1, 4, 3}, {0, 4, 2, 3}, {4, 1, 2, 3}};
        result.mBoundaryPreserved = KeepBoundary;
        return result;
    };
    return provider;
}

}  // namespace

TEST(Providers, RegisteredMesherIsUsedByTetrahedralize) {
    MesherGuard guard;
    int calls = 0;
    register_tet_mesher(fake_mesher(calls));
    EXPECT_TRUE(tetrahedralize_has_netgen());
    EXPECT_EQ(tet_mesher()->mVersion, "1.2.3");

    TetrahedralizeOptions options;
    options.mBackend = TetrahedralizeBackend::Netgen;
    const TetrahedralizeResult result = tetrahedralize(tetra_surface(), options);
    EXPECT_EQ(calls, 1);
    EXPECT_EQ(result.mBackend, "fake-mesher");
    EXPECT_EQ(result.mNumTetrahedra, 4);
    EXPECT_TRUE(result.mBoundaryPreserved);
    EXPECT_EQ(result.mMesh.NumPoints(), 5u);
}

TEST(Providers, MesherThatMovesTheBoundaryIsRefused) {
    MesherGuard guard;
    int calls = 0;
    register_tet_mesher(fake_mesher(calls, false));
    TetrahedralizeOptions options;
    options.mBackend = TetrahedralizeBackend::Netgen;
    EXPECT_THROW(tetrahedralize(tetra_surface(), options), std::runtime_error);
}

TEST(Providers, NamedMesherIsChosenByNameAndNotTheDefault) {
    MesherGuard guard;
    int calls = 0;
    TetMesherProvider named = fake_mesher(calls);
    named.mName = "fake-named";
    add_tet_mesher(std::move(named));
    unregister_tet_mesher();
    EXPECT_FALSE(tetrahedralize_has_netgen());  // the default slot stays empty
    ASSERT_NE(tet_mesher("fake-named"), nullptr);
    const std::vector<std::string> names = tetrahedralize_backends();
    EXPECT_NE(std::find(names.begin(), names.end(), "fake-named"), names.end());
    EXPECT_NE(std::find(names.begin(), names.end(), "octree"), names.end());

    TetrahedralizeOptions options;
    options.mBackend = TetrahedralizeBackend::Named;
    options.mProvider = "fake-named";
    const TetrahedralizeResult result = tetrahedralize(tetra_surface(), options);
    EXPECT_EQ(calls, 1);
    EXPECT_EQ(result.mBackend, "fake-named");
    EXPECT_TRUE(result.mBoundaryPreserved);

    // Auto never picks a named mesher: with no default it is the octree.
    options.mBackend = TetrahedralizeBackend::Auto;
    EXPECT_EQ(tetrahedralize(tetra_surface(), options).mBackend, "octree");
    EXPECT_EQ(calls, 1);

    // `octree` by name is the native backend.
    options.mBackend = TetrahedralizeBackend::Named;
    options.mProvider = "octree";
    EXPECT_EQ(tetrahedralize(tetra_surface(), options).mBackend, "octree");

    remove_tet_mesher("fake-named");
    EXPECT_EQ(tet_mesher("fake-named"), nullptr);
}

TEST(Providers, UnknownMesherNameListsTheAvailable) {
    MesherGuard guard;
    TetrahedralizeOptions options;
    options.mBackend = TetrahedralizeBackend::Named;
    options.mProvider = "no-such-mesher";
    try {
        tetrahedralize(tetra_surface(), options);
        FAIL() << "expected Unsupported";
    } catch (const Unsupported& rErr) {
        const std::string what = rErr.what();
        EXPECT_NE(what.find("'no-such-mesher'"), std::string::npos) << what;
        EXPECT_NE(what.find("octree"), std::string::npos) << what;
    }
    options.mProvider.clear();
    EXPECT_THROW(tetrahedralize(tetra_surface(), options), std::invalid_argument);
}

TEST(Providers, ABoundaryClaimIsCheckedNotTrusted) {
    MesherGuard guard;
    TetMesherProvider liar;
    liar.mName = "fake-liar";
    liar.mMeshRegion = [](const TetMeshRegionRequest& rRequest) {
        TetMeshRegionResult result;
        result.mPoints = rRequest.mPoints;
        result.mPoints[0][0] += 1e-3;  // moved, yet claims the surface was kept
        result.mTets = {{0, 1, 2, 3}};
        result.mBoundaryPreserved = true;
        return result;
    };
    add_tet_mesher(std::move(liar));
    TetrahedralizeOptions options;
    options.mBackend = TetrahedralizeBackend::Named;
    options.mProvider = "fake-liar";
    try {
        tetrahedralize(tetra_surface(), options);
        FAIL() << "expected a refusal";
    } catch (const std::runtime_error& rErr) {
        EXPECT_NE(std::string(rErr.what()).find("fake-liar mesher moved"), std::string::npos)
            << rErr.what();
    }
    remove_tet_mesher("fake-liar");
}

TEST(Providers, MesherFailureNamesTheRegion) {
    MesherGuard guard;
    TetMesherProvider provider;
    provider.mName = "failing";
    provider.mMeshRegion = [](const TetMeshRegionRequest&) -> TetMeshRegionResult {
        throw std::runtime_error("no luck");
    };
    register_tet_mesher(std::move(provider));
    TetrahedralizeOptions options;
    options.mBackend = TetrahedralizeBackend::Netgen;
    try {
        tetrahedralize(tetra_surface(), options);
        FAIL() << "expected the mesher's error";
    } catch (const std::runtime_error& rErr) {
        EXPECT_NE(std::string(rErr.what()).find("region 0: no luck"), std::string::npos);
    }
}

TEST(Providers, NoMesherNamesTheAddOn) {
    MesherGuard guard;
    unregister_tet_mesher();
    EXPECT_FALSE(tetrahedralize_has_netgen());
    TetrahedralizeOptions options;
    options.mBackend = TetrahedralizeBackend::Netgen;
    try {
        tetrahedralize(tetra_surface(), options);
        FAIL() << "expected Unsupported";
    } catch (const Unsupported& rErr) {
        const std::string what = rErr.what();
        EXPECT_NE(what.find("meshioplusplus-netgen"), std::string::npos);
        EXPECT_NE(what.find("MESHIOPLUSPLUS_WITH_NETGEN"), std::string::npos);
    }
}

// --- the remesher ----------------------------------------------------------------------

namespace {

class RemesherGuard {
public:
    RemesherGuard() : mSaved(remesher()) {}
    ~RemesherGuard() {
        unregister_remesher();
        if (mSaved)
            register_remesher(*mSaved);
    }

private:
    std::shared_ptr<const RemesherProvider> mSaved;
};

}  // namespace

TEST(Providers, RegisteredRemesherIsUsed) {
    RemesherGuard guard;
    int calls = 0;
    RemesherProvider provider;
    provider.mName = "fake-remesher";
    provider.mVersion = "5.8.0";
    provider.mElas = true;
    provider.mIndexBits = 64;
    provider.mAdapt = [&calls](const Mesh& rMesh, const AdaptOptions&) {
        ++calls;
        AdaptResult result;
        result.mMesh = detail::clone_mesh(rMesh);
        return result;
    };
    register_remesher(std::move(provider));
    EXPECT_TRUE(mmg_available());
    EXPECT_TRUE(mmg_has_elas());
    EXPECT_EQ(mmg_index_bits(), 64);

    const AdaptResult result = adapt(tetra_surface(), AdaptOptions{});
    EXPECT_EQ(calls, 1);
    EXPECT_EQ(result.mMesh.NumPoints(), 4u);
    // An operation the provider does not supply still names the add-on.
    EXPECT_THROW(move(tetra_surface(), MoveOptions{}), Unsupported);
}

TEST(Providers, NoRemesherNamesTheAddOn) {
    RemesherGuard guard;
    unregister_remesher();
    EXPECT_FALSE(mmg_available());
    EXPECT_FALSE(mmg_has_elas());
    EXPECT_EQ(mmg_index_bits(), 0);
    try {
        adapt(tetra_surface(), AdaptOptions{});
        FAIL() << "expected Unsupported";
    } catch (const Unsupported& rErr) {
        const std::string what = rErr.what();
        EXPECT_NE(what.find("meshioplusplus-mmg"), std::string::npos);
        EXPECT_NE(what.find("MESHIOPLUSPLUS_WITH_MMG"), std::string::npos);
    }
}

TEST(Providers, UnregisteringOneKindLeavesTheOthers) {
    KernelGuard cad_guard;
    MesherGuard mesher_guard;
    int calls = 0;
    register_cad_kernel(fake_kernel(calls));
    register_tet_mesher(fake_mesher(calls));
    unregister_tet_mesher();
    EXPECT_TRUE(cad_read_has_occt());
    EXPECT_FALSE(tetrahedralize_has_netgen());
}
