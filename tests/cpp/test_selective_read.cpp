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

// External includes
#include <gtest/gtest.h>

// System includes
#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>
#include <vector>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/xdmf.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/registry.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

// Selective/partial reads (read_options.hpp) and the read_metadata summary.
// Uniform-mesh-API only, so this runs under every mesh backend.

namespace {

/** @brief Writes a fixture mesh to a fresh temp VTU and returns its path. */
std::string selread_write_temp(const Mesh& rMesh) {
    const std::string path = mt::temp_path(".vtu");
    write_vtu(path, rMesh, /*binary=*/true, /*zlib=*/false);
    return path;
}

/** @brief Best-effort cleanup, mirroring `mt::roundtrip`'s. */
void selread_remove(const std::string& rPath) {
    std::error_code ec;
    std::filesystem::remove(rPath, ec);
}

}  // namespace

TEST(ReadOptions, DefaultsReadEverything) {
    const ReadOptions opts;
    EXPECT_FALSE(opts.mPointsOnly);
    EXPECT_FALSE(opts.mMetadataOnly);
    EXPECT_FALSE(opts.mDataArrays.has_value());
    EXPECT_TRUE(opts.WantsAnyData());
    EXPECT_TRUE(opts.WantsArray("anything"));
}

TEST(ReadOptions, NulloptAndEmptyVectorDiffer) {
    // The distinction that motivates std::optional: "every array" vs "no array".
    ReadOptions all;
    EXPECT_TRUE(all.WantsAnyData());
    EXPECT_TRUE(all.WantsArray("u"));

    ReadOptions none;
    none.mDataArrays = std::vector<std::string>{};
    EXPECT_FALSE(none.WantsAnyData());
    EXPECT_FALSE(none.WantsArray("u"));

    ReadOptions some;
    some.mDataArrays = std::vector<std::string>{"u"};
    EXPECT_TRUE(some.WantsAnyData());
    EXPECT_TRUE(some.WantsArray("u"));
    EXPECT_FALSE(some.WantsArray("v"));
}

TEST(ReadOptions, PointsOnlyAndMetadataOnlySuppressData) {
    ReadOptions points_only;
    points_only.mPointsOnly = true;
    EXPECT_FALSE(points_only.WantsAnyData());

    ReadOptions metadata_only;
    metadata_only.mMetadataOnly = true;
    EXPECT_FALSE(metadata_only.WantsAnyData());
}

TEST(MetadataFromMesh, ReportsShapeAndNames) {
    const Mesh mesh = mt::tri_quad_mesh();
    const MeshMetadata meta = metadata_from_mesh(mesh);

    EXPECT_EQ(meta.mNumPoints, mesh.NumPoints());
    EXPECT_EQ(meta.mPointDim, mesh.PointDim());
    ASSERT_EQ(meta.mCellBlocks.size(), mesh.NumCellBlocks());

    std::size_t total = 0;
    for (std::size_t i = 0; i < mesh.NumCellBlocks(); ++i) {
        const auto block = mesh.Cells(i);
        EXPECT_EQ(meta.mCellBlocks[i].mType, block.Type());
        EXPECT_EQ(meta.mCellBlocks[i].mNumCells, block.NumCells());
        EXPECT_EQ(meta.mCellBlocks[i].mRagged, block.IsRagged());
        total += block.NumCells();
    }
    EXPECT_EQ(meta.NumCells(), total);
}

TEST(MetadataFromMesh, ReportsDataNamesAndBBox) {
    const Mesh mesh = mt::data_mesh();
    const MeshMetadata meta = metadata_from_mesh(mesh);

    EXPECT_EQ(meta.mPointDataNames, mesh.PointDataNames());
    EXPECT_EQ(meta.mCellDataNames, mesh.CellDataNames());
    EXPECT_EQ(meta.mFieldDataNames, mesh.FieldDataNames());

    // The mesh is already loaded, so the bbox is cheap and is filled in.
    ASSERT_TRUE(meta.mHasBBox);
    for (std::size_t d = 0; d < 3 && d < meta.mPointDim; ++d)
        EXPECT_LE(meta.mBBoxMin[d], meta.mBBoxMax[d]);
}

TEST(MetadataFromMesh, EmptyMeshHasNoBBox) {
    Mesh empty;
    const MeshMetadata meta = metadata_from_mesh(empty);
    EXPECT_EQ(meta.mNumPoints, 0u);
    EXPECT_TRUE(meta.mCellBlocks.empty());
    EXPECT_FALSE(meta.mHasBBox);
}

TEST(MetadataFromMesh, ReportsRegionsWithoutEntries) {
    // The mesh is already in memory, so regions cost nothing extra to report:
    // metadata_from_mesh must summarize every region already on the mesh.
    Mesh mesh = mt::tri_quad_mesh();

    NDArray cell_entries = NDArray::Uninit(DType::Int64, {2});
    cell_entries.As<std::int64_t>()[0] = 0;
    cell_entries.As<std::int64_t>()[1] = 1;
    mesh.AddRegion(Region("solid", RegionKind::Cell, 2, 7, std::move(cell_entries)));

    NDArray point_entries = NDArray::Uninit(DType::Int64, {1});
    point_entries.As<std::int64_t>()[0] = 0;
    mesh.AddRegion(Region("fixed", RegionKind::Point, std::move(point_entries)));

    const MeshMetadata meta = metadata_from_mesh(mesh);
    ASSERT_EQ(meta.mRegions.size(), 2u);

    const RegionSummary* solid = nullptr;
    const RegionSummary* fixed = nullptr;
    for (const auto& r : meta.mRegions) {
        if (r.mName == "solid")
            solid = &r;
        else if (r.mName == "fixed")
            fixed = &r;
    }
    ASSERT_NE(solid, nullptr);
    ASSERT_NE(fixed, nullptr);
    EXPECT_EQ(solid->mKind, RegionKind::Cell);
    EXPECT_EQ(solid->mDim, 2);
    EXPECT_EQ(solid->mTag, 7);
    EXPECT_EQ(solid->mNumEntries, 2u);
    EXPECT_EQ(fixed->mKind, RegionKind::Point);
    EXPECT_EQ(fixed->mDim, -1);
    EXPECT_EQ(fixed->mTag, -1);
    EXPECT_EQ(fixed->mNumEntries, 1u);
}

TEST(MetadataFromMesh, NoRegionsIsEmptyNotMissing) {
    const Mesh mesh = mt::tri_quad_mesh();
    const MeshMetadata meta = metadata_from_mesh(mesh);
    EXPECT_TRUE(meta.mRegions.empty());
}

TEST(RegistryReadMetadata, MatchesAFullReadSummary) {
    const Mesh source = mt::data_mesh();
    const std::string path = selread_write_temp(source);

    const MeshMetadata meta = registry_read_metadata(path, "vtu", ReadOptions{});
    const MeshMetadata reference = metadata_from_mesh(read_vtu(path));

    EXPECT_EQ(meta.mFormat, "vtu");
    EXPECT_EQ(meta.mNumPoints, reference.mNumPoints);
    ASSERT_EQ(meta.mCellBlocks.size(), reference.mCellBlocks.size());
    for (std::size_t i = 0; i < meta.mCellBlocks.size(); ++i) {
        EXPECT_EQ(meta.mCellBlocks[i].mType, reference.mCellBlocks[i].mType);
        EXPECT_EQ(meta.mCellBlocks[i].mNumCells, reference.mCellBlocks[i].mNumCells);
    }
    EXPECT_EQ(meta.mPointDataNames, reference.mPointDataNames);
    EXPECT_EQ(meta.mCellDataNames, reference.mCellDataNames);

    selread_remove(path);
}

TEST(RegistryReadMetadata, FlagsTheFullReadFallback) {
    const std::string path = selread_write_temp(mt::tri_mesh());

    const MeshMetadata meta = registry_read_metadata(path, "vtu", ReadOptions{});
    // The flag must agree with the table rather than being independently guessed.
    EXPECT_EQ(meta.mFellBackToFullRead, !registry_reader_supports_options("vtu") ||
                                            registry_metadata_readers().count("vtu") == 0);

    selread_remove(path);
}

TEST(RegistryRead, DefaultOptionsMatchAPlainRead) {
    const Mesh source = mt::data_mesh();
    const std::string path = selread_write_temp(source);

    const Mesh via_options = registry_read(path, "vtu", ReadOptions{});
    const Mesh via_plain = read_vtu(path);
    mt::expect_same_geometry(via_options, via_plain);
    EXPECT_EQ(via_options.PointDataNames(), via_plain.PointDataNames());
    EXPECT_EQ(via_options.CellDataNames(), via_plain.CellDataNames());

    selread_remove(path);
}

// ---------------------------------------------------------------------------
// Native VTU/VTP selective + metadata paths
// ---------------------------------------------------------------------------

TEST(SelectiveReadVtu, VtuAndVtpHaveNativePaths) {
    // If these regress to a fallback the tests below still pass but stop
    // testing anything, so assert the wiring itself.
    EXPECT_TRUE(registry_reader_supports_options("vtu"));
    EXPECT_TRUE(registry_reader_supports_options("vtp"));
    EXPECT_EQ(registry_metadata_readers().count("vtu"), 1u);
    EXPECT_EQ(registry_metadata_readers().count("vtp"), 1u);
}

TEST(SelectiveReadVtu, PointsOnlyKeepsGeometryDropsData) {
    const Mesh source = mt::data_mesh();
    ASSERT_FALSE(source.PointDataNames().empty());
    const std::string path = selread_write_temp(source);

    ReadOptions opts;
    opts.mPointsOnly = true;
    const Mesh got = read_vtu(path, opts);

    // Geometry survives bit-identically; only the data is gone.
    mt::expect_same_geometry(got, read_vtu(path));
    EXPECT_EQ(got.NumPointData(), 0u);
    EXPECT_EQ(got.NumCellData(), 0u);

    selread_remove(path);
}

TEST(SelectiveReadVtu, DataArraysSubsetReturnsExactlyThose) {
    const Mesh source = mt::data_mesh();
    const std::vector<std::string> all = source.PointDataNames();
    ASSERT_GE(all.size(), 2u) << "fixture must have >= 2 point_data arrays to subset";
    const std::string keep = all.front();
    const std::string path = selread_write_temp(source);

    ReadOptions opts;
    opts.mDataArrays = std::vector<std::string>{keep};
    const Mesh got = read_vtu(path, opts);

    ASSERT_EQ(got.PointDataNames().size(), 1u);
    EXPECT_EQ(got.PointDataNames().front(), keep);
    mt::expect_same_geometry(got, read_vtu(path));

    selread_remove(path);
}

TEST(SelectiveReadVtu, UnknownRequestedNamesAreIgnoredNotFatal) {
    const std::string path = selread_write_temp(mt::data_mesh());

    ReadOptions opts;
    opts.mDataArrays = std::vector<std::string>{"no-such-array"};
    const Mesh got = read_vtu(path, opts);  // must not throw
    EXPECT_EQ(got.NumPointData(), 0u);

    selread_remove(path);
}

TEST(SelectiveReadVtu, EmptyArrayListDropsAllDataButKeepsGeometry) {
    const std::string path = selread_write_temp(mt::data_mesh());

    ReadOptions opts;
    opts.mDataArrays = std::vector<std::string>{};  // "none", not "all"
    const Mesh got = read_vtu(path, opts);
    EXPECT_EQ(got.NumPointData(), 0u);
    EXPECT_GT(got.NumPoints(), 0u);

    selread_remove(path);
}

TEST(SelectiveReadVtu, DefaultOptionsStillReadEverything) {
    const Mesh source = mt::data_mesh();
    const std::string path = selread_write_temp(source);

    const Mesh got = read_vtu(path, ReadOptions{});
    EXPECT_EQ(got.PointDataNames(), read_vtu(path).PointDataNames());
    EXPECT_EQ(got.NumPointData(), source.NumPointData());

    selread_remove(path);
}

TEST(SelectiveReadVtu, MetadataMatchesAFullReadAndIsNotAFallback) {
    const Mesh source = mt::data_mesh();
    const std::string path = selread_write_temp(source);

    const MeshMetadata meta = registry_read_metadata(path, "vtu", ReadOptions{});
    const MeshMetadata reference = metadata_from_mesh(read_vtu(path));

    EXPECT_FALSE(meta.mFellBackToFullRead) << "vtu has a native metadata path";
    EXPECT_EQ(meta.mNumPoints, reference.mNumPoints);
    EXPECT_EQ(meta.mPointDim, reference.mPointDim);
    EXPECT_EQ(meta.mPointDataNames, reference.mPointDataNames);
    EXPECT_EQ(meta.mCellDataNames, reference.mCellDataNames);
    // The bbox is deliberately not computed on the native path.
    EXPECT_FALSE(meta.mHasBBox);
    EXPECT_TRUE(reference.mHasBBox);

    selread_remove(path);
}

// This is what keeps summarize_cells honest: it duplicates reconstruct_cells'
// run-grouping, so it is pinned against a real reconstruction rather than a
// comment asking future readers to keep them in sync.
TEST(SelectiveReadVtu, SummarizedBlocksMatchReconstructedBlocks) {
    const std::vector<Mesh (*)()> fixtures = {
        &mt::tri_mesh, &mt::quad_mesh, &mt::tet_mesh, &mt::hex_mesh, &mt::tri_quad_mesh,
    };
    for (auto make : fixtures) {
        const Mesh source = make();
        const std::string path = selread_write_temp(source);

        const MeshMetadata meta = read_vtu_metadata(path);
        const MeshMetadata reference = metadata_from_mesh(read_vtu(path));

        ASSERT_EQ(meta.mCellBlocks.size(), reference.mCellBlocks.size());
        for (std::size_t i = 0; i < meta.mCellBlocks.size(); ++i) {
            EXPECT_EQ(meta.mCellBlocks[i].mType, reference.mCellBlocks[i].mType);
            EXPECT_EQ(meta.mCellBlocks[i].mNumCells, reference.mCellBlocks[i].mNumCells);
            EXPECT_EQ(meta.mCellBlocks[i].mNodesPerCell, reference.mCellBlocks[i].mNodesPerCell);
        }
        selread_remove(path);
    }
}

TEST(SelectiveReadVtp, PointsOnlyAndMetadataWork) {
    const Mesh source = mt::tri_mesh();
    const std::string path = mt::temp_path(".vtp");
    write_vtp(path, source, /*binary=*/true, /*zlib=*/false);

    ReadOptions opts;
    opts.mPointsOnly = true;
    const Mesh got = read_vtp(path, opts);
    EXPECT_EQ(got.NumPointData(), 0u);
    EXPECT_EQ(got.NumPoints(), source.NumPoints());

    const MeshMetadata meta = read_vtp_metadata(path);
    const MeshMetadata reference = metadata_from_mesh(read_vtp(path));
    EXPECT_EQ(meta.mNumPoints, reference.mNumPoints);
    ASSERT_EQ(meta.mCellBlocks.size(), reference.mCellBlocks.size());
    for (std::size_t i = 0; i < meta.mCellBlocks.size(); ++i) {
        EXPECT_EQ(meta.mCellBlocks[i].mType, reference.mCellBlocks[i].mType);
        EXPECT_EQ(meta.mCellBlocks[i].mNumCells, reference.mCellBlocks[i].mNumCells);
    }

    selread_remove(path);
}

TEST(SelectiveReadVtu, MetadataRejectsWhatTheReaderRejects) {
    // A summary must never claim a file is fine when read_vtu would refuse it.
    const std::string path = mt::temp_path(".vtu");
    {
        std::ofstream os(path);
        os << "<VTKFile type=\"UnstructuredGrid\"><UnstructuredGrid>"
           << "<Piece NumberOfPoints=\"1\"/><Piece NumberOfPoints=\"1\"/>"
           << "</UnstructuredGrid></VTKFile>";
    }
    EXPECT_THROW(read_vtu(path), ReadError);
    EXPECT_THROW(read_vtu_metadata(path), ReadError);
    selread_remove(path);
}

TEST(RegistryRead, UnknownFormatThrows) {
    EXPECT_THROW(registry_read("x.zzz", "definitely-not-a-format", ReadOptions{}), ReadError);
    EXPECT_THROW(registry_read_metadata("x.zzz", "definitely-not-a-format", ReadOptions{}),
                 ReadError);
}

// ---------------------------------------------------------------------------
// XDMF: the cheapest metadata path -- DataItem @Dimensions gives exact counts
// without reading any payload (or opening the sibling .h5 on the HDF path).
// ---------------------------------------------------------------------------

TEST(SelectiveReadXdmf, MetadataMatchesFullReadWithoutPayload) {
    // Single cell type -> a plain <Topology>, which the summary can describe
    // from @Dimensions alone. (Multi-block meshes become Mixed; see below.)
    const Mesh source = mt::tri_mesh();
    const std::string path = mt::temp_path(".xdmf");
    write_xdmf(path, source, "XML");

    const MeshMetadata meta = read_xdmf_metadata(path);
    const MeshMetadata reference = metadata_from_mesh(read_xdmf(path));

    EXPECT_EQ(meta.mNumPoints, reference.mNumPoints);
    EXPECT_EQ(meta.mPointDim, reference.mPointDim);
    ASSERT_EQ(meta.mCellBlocks.size(), reference.mCellBlocks.size());
    for (std::size_t i = 0; i < meta.mCellBlocks.size(); ++i) {
        EXPECT_EQ(meta.mCellBlocks[i].mType, reference.mCellBlocks[i].mType);
        EXPECT_EQ(meta.mCellBlocks[i].mNumCells, reference.mCellBlocks[i].mNumCells);
        EXPECT_EQ(meta.mCellBlocks[i].mNodesPerCell, reference.mCellBlocks[i].mNodesPerCell);
    }
    EXPECT_EQ(meta.mPointDataNames, reference.mPointDataNames);
    EXPECT_EQ(meta.mCellDataNames, reference.mCellDataNames);
    EXPECT_FALSE(meta.mHasBBox);

    selread_remove(path);
}

// A native summary is allowed to decline a construct it cannot describe
// cheaply -- but declining must cost a slower answer, never a failed one.
TEST(SelectiveReadXdmf, MixedTopologyFallsBackInsteadOfThrowing) {
    const Mesh source = mt::tri_quad_mesh();  // two blocks -> Mixed topology
    const std::string path = mt::temp_path(".xdmf");
    write_xdmf(path, source, "XML");

    // The native path itself declines...
    EXPECT_THROW(read_xdmf_metadata(path), ReadError);

    // ...but the registry still answers, and says the answer was not cheap.
    const MeshMetadata meta = registry_read_metadata(path, "xdmf", ReadOptions{});
    EXPECT_TRUE(meta.mFellBackToFullRead);
    EXPECT_EQ(meta.mFormat, "xdmf");
    EXPECT_EQ(meta.mNumPoints, source.NumPoints());
    EXPECT_EQ(meta.mCellBlocks.size(), source.NumCellBlocks());
    EXPECT_TRUE(meta.mHasBBox) << "the full-read fallback can afford a bbox";

    selread_remove(path);
}

TEST(SelectiveReadXdmf, PointsOnlyAndSubsetSkipAttributes) {
    const Mesh source = mt::data_mesh();
    const std::vector<std::string> all = source.PointDataNames();
    ASSERT_GE(all.size(), 2u);
    const std::string path = mt::temp_path(".xdmf");
    write_xdmf(path, source, "XML");

    ReadOptions points_only;
    points_only.mPointsOnly = true;
    const Mesh bare = read_xdmf(path, points_only);
    EXPECT_EQ(bare.NumPointData(), 0u);
    EXPECT_EQ(bare.NumCellData(), 0u);
    mt::expect_same_geometry(bare, read_xdmf(path));

    ReadOptions subset;
    subset.mDataArrays = std::vector<std::string>{all.front()};
    const Mesh some = read_xdmf(path, subset);
    ASSERT_EQ(some.PointDataNames().size(), 1u);
    EXPECT_EQ(some.PointDataNames().front(), all.front());

    selread_remove(path);
}

// ---------------------------------------------------------------------------
// Gmsh: $NodeData/$ElementData carry their name above the values, so an
// unwanted section is skipped with skip_to_end instead of being parsed.
// ---------------------------------------------------------------------------

TEST(SelectiveReadGmsh, PointsOnlyAndSubsetSkipDataSections) {
    const Mesh source = mt::data_mesh();
    const std::vector<std::string> all = source.PointDataNames();
    ASSERT_GE(all.size(), 2u);
    const std::string path = mt::temp_path(".msh");
    write_gmsh41(path, source, /*binary=*/false);

    ReadOptions points_only;
    points_only.mPointsOnly = true;
    const Mesh bare = read_gmsh(path, points_only);
    EXPECT_EQ(bare.NumPointData(), 0u);
    mt::expect_same_geometry(bare, read_gmsh(path));

    ReadOptions subset;
    subset.mDataArrays = std::vector<std::string>{all.front()};
    const Mesh some = read_gmsh(path, subset);
    ASSERT_EQ(some.PointDataNames().size(), 1u);
    EXPECT_EQ(some.PointDataNames().front(), all.front());
    mt::expect_same_geometry(some, read_gmsh(path));

    selread_remove(path);
}

TEST(SelectiveReadGmsh, DefaultOptionsStillReadEverything) {
    const Mesh source = mt::data_mesh();
    const std::string path = mt::temp_path(".msh");
    write_gmsh41(path, source, /*binary=*/false);

    EXPECT_EQ(read_gmsh(path, ReadOptions{}).PointDataNames(), read_gmsh(path).PointDataNames());

    selread_remove(path);
}

// The gmsh 4.1 metadata scanner walks block headers and skips payloads by byte
// arithmetic (binary) or line counts (ascii). Both are easy to get subtly wrong
// and still "work", so they are pinned against a real read, per variant.
TEST(SelectiveReadGmsh, NativeMetadataMatchesAFullRead41) {
    const std::vector<Mesh (*)()> fixtures = {
        &mt::tri_mesh, &mt::tet_mesh, &mt::hex_mesh, &mt::tri_quad_mesh, &mt::data_mesh,
    };
    for (bool binary : {false, true}) {
        for (auto make : fixtures) {
            const Mesh source = make();
            const std::string path = mt::temp_path(".msh");
            write_gmsh41(path, source, binary);

            const MeshMetadata meta = read_gmsh_metadata(path);
            const MeshMetadata reference = metadata_from_mesh(read_gmsh(path));

            EXPECT_EQ(meta.mNumPoints, reference.mNumPoints) << (binary ? "binary" : "ascii");
            ASSERT_EQ(meta.mCellBlocks.size(), reference.mCellBlocks.size())
                << (binary ? "binary" : "ascii");
            for (std::size_t i = 0; i < meta.mCellBlocks.size(); ++i) {
                EXPECT_EQ(meta.mCellBlocks[i].mType, reference.mCellBlocks[i].mType);
                EXPECT_EQ(meta.mCellBlocks[i].mNumCells, reference.mCellBlocks[i].mNumCells);
                EXPECT_EQ(meta.mCellBlocks[i].mNodesPerCell,
                          reference.mCellBlocks[i].mNodesPerCell);
            }
            EXPECT_EQ(meta.mPointDataNames, reference.mPointDataNames)
                << (binary ? "binary" : "ascii");
            EXPECT_EQ(meta.mCellDataNames, reference.mCellDataNames)
                << (binary ? "binary" : "ascii");
            EXPECT_FALSE(meta.mHasBBox);

            selread_remove(path);
        }
    }
}

TEST(SelectiveReadGmsh, NativeMetadataIsUsedByTheRegistry) {
    const std::string path = mt::temp_path(".msh");
    write_gmsh41(path, mt::tet_mesh(), /*binary=*/false);
    const MeshMetadata meta = registry_read_metadata(path, "gmsh", ReadOptions{});
    EXPECT_FALSE(meta.mFellBackToFullRead) << "gmsh 4.1 now has a native path";
    selread_remove(path);
}

// Format 2.2 has a per-element type, so there is no cheap summary to give. It
// must decline and fall back -- correct but slower -- never return a wrong one.
TEST(SelectiveReadGmsh, Format22DeclinesAndFallsBack) {
    const Mesh source = mt::tri_mesh();
    const std::string path = mt::temp_path(".msh");
    write_gmsh22(path, source, /*binary=*/false);

    EXPECT_THROW(read_gmsh_metadata(path), ReadError);

    const MeshMetadata meta = registry_read_metadata(path, "gmsh", ReadOptions{});
    EXPECT_TRUE(meta.mFellBackToFullRead);
    EXPECT_EQ(meta.mNumPoints, source.NumPoints());
    EXPECT_EQ(meta.mCellBlocks.size(), source.NumCellBlocks());

    selread_remove(path);
}
