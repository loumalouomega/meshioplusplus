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
 * @file test_facet_tables.cpp
 * @brief Pins the output of every operation built on a facet or edge table.
 *
 * `extract_surface`, `extract_skin`, `smooth`'s boundary pass, `refine`,
 * `convert_cells` (elevate), `decimate` and `build_global_faces` each number or
 * count facets or edges keyed by their node ids. Roadmap §3 replaces their
 * single-threaded hash maps with one sort-based table; these digests, taken
 * from the hash-map implementation, prove the replacement changes no byte.
 * The repeated-run checks hold everywhere; the golden digests are pinned for
 * x86-64 GCC/Clang only, where no compiler contracts a*b+c into an FMA.
 */

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <vector>

#include <gtest/gtest.h>

#include "golden_fixtures.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/skin.hpp"

namespace {

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::bench::MeshDigest;

using namespace golden;

template <class T>
void ft_vec(MeshDigest& rD, const std::vector<T>& rV) {
    rD.U64(rV.size());
    rD.Bytes(rV.data(), rV.size() * sizeof(T));
}

std::uint64_t ft_faces(const Mesh& rM) {
    const meshioplusplus::detail::GlobalFaces f = meshioplusplus::detail::build_global_faces(rM);
    MeshDigest d;
    ft_vec(d, f.mFaceNodes);
    ft_vec(d, f.mFaceStart);
    ft_vec(d, f.mOwner);
    ft_vec(d, f.mNeighbour);
    ft_vec(d, f.mCellFaces);
    ft_vec(d, f.mCellFaceStart);
    ft_vec(d, f.mCellToGlobal);
    d.U64(static_cast<std::uint64_t>(f.mNumFlipped));
    d.U64(static_cast<std::uint64_t>(f.mNumUnorientable));
    d.U64(static_cast<std::uint64_t>(f.mNumNonManifold));
    // The lookup finds every face under any rotation of its corners.
    const meshioplusplus::detail::FaceLookup lookup(f);
    for (std::size_t i = 0; i < f.NumFaces(); ++i) {
        std::vector<std::int64_t> ring(f.Face(i), f.Face(i) + f.FaceSize(i));
        std::rotate(ring.begin(), ring.begin() + 1, ring.end());
        EXPECT_EQ(lookup.Find(ring.data(), ring.size()), static_cast<std::int64_t>(i));
    }
    const std::int64_t absent[3] = {0, 1, 1000000};
    EXPECT_EQ(lookup.Find(absent, 3), -1);
    return d.Value();
}

struct FtCase {
    const char* mName;
    std::uint64_t (*mRun)();
    std::uint64_t mGolden;
};

std::uint64_t ft_surface_tets() {
    return ft_mesh(meshioplusplus::extract_surface(ft_tet_cube(6, 0.2), true));
}
std::uint64_t ft_surface_mixed() {
    return ft_mesh(meshioplusplus::extract_surface(ft_hex_grid(5, 0.2, true), true));
}
std::uint64_t ft_surface_pentagons() {
    return ft_mesh(meshioplusplus::extract_surface(ft_pentagon_prisms(), true));
}
std::uint64_t ft_skin_hexes() {
    return ft_mesh(meshioplusplus::extract_skin(ft_hex_grid(5, 0.2, false)));
}
std::uint64_t ft_smooth_tets() {
    meshioplusplus::SmoothOptions o;
    o.mIterations = 4;
    return ft_mesh(meshioplusplus::smooth(ft_tet_cube(6, 0.3), o).mMesh);
}
std::uint64_t ft_smooth_mixed() {
    meshioplusplus::SmoothOptions o;
    o.mMethod = meshioplusplus::SmoothMethod::Laplacian;
    o.mIterations = 3;
    return ft_mesh(meshioplusplus::smooth(ft_hex_grid(5, 0.3, true), o).mMesh);
}
std::uint64_t ft_refine_tets() {
    return ft_mesh(meshioplusplus::refine(ft_tet_cube(4, 0.2)).mMesh);
}
std::uint64_t ft_refine_hexes() {
    return ft_mesh(meshioplusplus::refine(ft_hex_grid(4, 0.2, false)).mMesh);
}
std::uint64_t ft_refine_green() {
    meshioplusplus::RefineOptions o;
    o.mCells = {0, 7, 20, 51, 100, 222};
    o.mRecordHierarchy = true;
    return ft_mesh(meshioplusplus::refine(ft_tet_cube(4, 0.2), o).mMesh);
}
std::uint64_t ft_elevate(const Mesh& rM) {
    meshioplusplus::ConvertCellsOptions o;
    o.mMode = meshioplusplus::ConvertCellsMode::Elevate;
    auto r = meshioplusplus::convert_cells(rM, o);
    MeshDigest d;
    d.Of(r.mMesh);
    d.Array(r.mPointMap);
    d.Arrays(r.mCellMaps);
    return d.Value();
}
std::uint64_t ft_elevate_tets() {
    return ft_elevate(ft_tet_cube(5, 0.2));
}
std::uint64_t ft_elevate_hexes() {
    return ft_elevate(ft_hex_grid(4, 0.2, false));
}
std::uint64_t ft_decimate() {
    meshioplusplus::DecimateOptions o;
    o.mTargetRatio = 0.5;
    auto r = meshioplusplus::decimate(meshioplusplus::extract_surface(ft_tet_cube(8, 0.2)), o);
    MeshDigest d;
    d.Of(r.mMesh);
    d.Array(r.mPointMap);
    return d.Value();
}
std::uint64_t ft_faces_mixed() {
    return ft_faces(ft_hex_grid(5, 0.2, true));
}
std::uint64_t ft_faces_tets() {
    return ft_faces(ft_tet_cube(5, 0.2));
}
std::uint64_t ft_faces_pentagons() {
    return ft_faces(ft_pentagon_prisms());
}

const FtCase kFtCases[] = {
    {"surface_tets", ft_surface_tets, 0xa4037e5af12463f2ull},
    {"surface_mixed", ft_surface_mixed, 0x33b5cebfec51daa5ull},
    {"surface_pentagons", ft_surface_pentagons, 0xf6a322eea57e5d9dull},
    {"skin_hexes", ft_skin_hexes, 0x29742b403e1530daull},
    {"smooth_tets", ft_smooth_tets, 0x8b07f97e5039bee8ull},
    {"smooth_mixed", ft_smooth_mixed, 0xca367aebf2b0dffdull},
    {"refine_tets", ft_refine_tets, 0xab469d9e686a8e74ull},
    {"refine_hexes", ft_refine_hexes, 0x3511f2592c5bdc43ull},
    {"refine_green", ft_refine_green, 0xbfc2e13ea448af96ull},
    {"elevate_tets", ft_elevate_tets, 0xbf57f48a1c48cc5full},
    {"elevate_hexes", ft_elevate_hexes, 0xa4de69ef83591013ull},
    {"decimate", ft_decimate, 0x100c0b7e82c53f50ull},
    {"faces_mixed", ft_faces_mixed, 0xbe00299381c43fd8ull},
    {"faces_tets", ft_faces_tets, 0xdf9c379133c18f32ull},
    {"faces_pentagons", ft_faces_pentagons, 0xc494ba2c38f47b98ull},
};

TEST(FacetTables, ResultsAreStableAcrossRepeatedRuns) {
    for (const FtCase& c : kFtCases) {
        const std::uint64_t first = c.mRun();
        EXPECT_EQ(c.mRun(), first) << c.mName;
    }
}

#ifdef GOLDEN_PINNED
TEST(FacetTables, ResultsMatchTheHashMapImplementation) {
    for (const FtCase& c : kFtCases) {
        const std::uint64_t got = c.mRun();
        EXPECT_EQ(got, c.mGolden) << c.mName << ": 0x" << std::hex << got;
    }
}
#endif

}  // namespace
