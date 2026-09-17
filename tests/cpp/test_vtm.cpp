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
// VTK XML MultiBlock (roadmap §1 tier B4, v11.6.0, part 3 of 3).

// System includes
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtm.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/region.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::read_vtm;
using meshioplusplus::read_vtm_metadata;
using meshioplusplus::read_vtu;
using meshioplusplus::ReadError;
using meshioplusplus::RegionKind;
using meshioplusplus::write_vtm;

namespace {

// Remove the index file and its sibling pieces directory (`<stem>/`).
void vtm_cleanup(const std::string& rPath) {
    std::error_code ec;
    const std::filesystem::path p(rPath);
    std::filesystem::remove(p, ec);
    std::filesystem::remove_all(p.parent_path() / p.stem(), ec);
}

}  // namespace

// The roadmap's own stated probe: write two (here three) blocks to .vtm, read
// them back. `tri_quad_mesh()` has three blocks -- triangle, quad, triangle --
// sharing points across blocks, so this also exercises point duplication
// across pieces and same-type block consolidation on the way back in.
TEST(Vtm, RoundTripsThreeBlocksAsThreeNamedCellRegions) {
    const Mesh m = mt::tri_quad_mesh();
    const std::string path = mt::temp_path(".vtm");
    write_vtm(path, m, /*binary=*/false, /*zlib=*/false);
    const Mesh back = read_vtm(path);

    // Each piece is pruned to only the points its own block references:
    // {0,1,5,6} + {1,2,4,5} + {2,3,4} = 4 + 4 + 3, shared points duplicated.
    EXPECT_EQ(back.NumPoints(), 11u);
    // The two "triangle" pieces (index 0 and 2) consolidate into one block,
    // same as `merge()` would for any two same-typed inputs.
    ASSERT_EQ(back.NumCellBlocks(), 2u);
    EXPECT_EQ(back.Cells(0).Type(), "triangle");
    EXPECT_EQ(back.Cells(0).NumCells(), 3u);
    EXPECT_EQ(back.Cells(1).Type(), "quad");
    EXPECT_EQ(back.Cells(1).NumCells(), 1u);

    ASSERT_EQ(back.NumRegions(), 3u);
    std::size_t total_entries = 0;
    for (std::size_t i = 0; i < back.NumRegions(); ++i) {
        const meshioplusplus::Region& r = back.Region(i);
        EXPECT_EQ(r.mKind, RegionKind::Cell);
        EXPECT_EQ(r.mName, "block_" + std::to_string(i));
        total_entries += r.NumEntries();
    }
    EXPECT_EQ(total_entries, 4u);  // 2 + 1 + 1 original cells

    vtm_cleanup(path);
}

TEST(Vtm, PiecesAreIndependentlyReadableVtuFiles) {
    const Mesh m = mt::tri_quad_mesh();
    const std::string path = mt::temp_path(".vtm");
    write_vtm(path, m, /*binary=*/true, /*zlib=*/false);

    const std::filesystem::path idx(path);
    const std::filesystem::path dir = idx.parent_path() / idx.stem();
    const std::string stem = idx.stem().string();
    const std::size_t expected_cells[3] = {2, 1, 1};
    for (int i = 0; i < 3; ++i) {
        const std::string piece = (dir / (stem + "_" + std::to_string(i) + ".vtu")).string();
        ASSERT_TRUE(std::filesystem::exists(piece)) << piece;
        const Mesh p = read_vtu(piece);
        ASSERT_EQ(p.NumCellBlocks(), 1u);
        EXPECT_EQ(p.Cells(0).NumCells(), expected_cells[i]);
    }

    vtm_cleanup(path);
}

TEST(Vtm, MetadataAgreesWithARealRead) {
    const Mesh m = mt::tri_quad_mesh();
    const std::string path = mt::temp_path(".vtm");
    write_vtm(path, m, true, false);

    const meshioplusplus::MeshMetadata meta = read_vtm_metadata(path);
    const Mesh back = read_vtm(path);
    EXPECT_EQ(meta.mNumPoints, back.NumPoints());
    ASSERT_EQ(meta.mCellBlocks.size(), back.NumCellBlocks());
    for (std::size_t i = 0; i < meta.mCellBlocks.size(); ++i) {
        EXPECT_EQ(meta.mCellBlocks[i].mType, back.Cells(i).Type());
        EXPECT_EQ(meta.mCellBlocks[i].mNumCells, back.Cells(i).NumCells());
    }

    vtm_cleanup(path);
}

TEST(Vtm, EmptyIndexReadsAsAnEmptyMesh) {
    const std::string path = mt::temp_path(".vtm");
    std::ofstream os(path, std::ios::binary);
    os << "<?xml version=\"1.0\"?>\n"
          "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">\n"
          "<vtkMultiBlockDataSet>\n"
          "<Block index=\"0\"/>\n"
          "</vtkMultiBlockDataSet>\n"
          "</VTKFile>\n";
    os.close();

    const Mesh back = read_vtm(path);
    EXPECT_EQ(back.NumPoints(), 0u);
    EXPECT_EQ(back.NumCellBlocks(), 0u);

    const meshioplusplus::MeshMetadata meta = read_vtm_metadata(path);
    EXPECT_EQ(meta.mNumPoints, 0u);
    EXPECT_TRUE(meta.mCellBlocks.empty());

    std::filesystem::remove(path);
}

TEST(Vtm, DeclinesAPieceThatIsNeitherVtuNorVtp) {
    const std::string path = mt::temp_path(".vtm");
    std::ofstream os(path, std::ios::binary);
    os << "<?xml version=\"1.0\"?>\n"
          "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">\n"
          "<vtkMultiBlockDataSet>\n"
          "<Block index=\"0\">\n"
          "<DataSet index=\"0\" name=\"a\" file=\"piece.obj\"/>\n"
          "</Block>\n"
          "</vtkMultiBlockDataSet>\n"
          "</VTKFile>\n";
    os.close();

    EXPECT_THROW(read_vtm(path), ReadError);
    std::filesystem::remove(path);
}
