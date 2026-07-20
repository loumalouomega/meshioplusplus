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

// Project includes
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/read_options.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;
using meshioplusplus::detail::FileSource;

namespace {

void mmap_remove(const std::string& rPath) {
    std::error_code ec;
    std::filesystem::remove(rPath, ec);
}

std::string mmap_write_text(const std::string& rSuffix, const std::string& rText) {
    const std::string path = mt::temp_path(rSuffix);
    std::ofstream os(path, std::ios::binary);
    os.write(rText.data(), static_cast<std::streamsize>(rText.size()));
    return path;
}

}  // namespace

TEST(FileSource, BufferedReadReturnsWholeFile) {
    const std::string text = "hello\nmeshio++\n";
    const std::string path = mmap_write_text(".txt", text);

    const FileSource source(path, FileSource::Mode::Buffered);
    EXPECT_FALSE(source.IsMapped());
    EXPECT_EQ(source.Size(), text.size());
    EXPECT_EQ(std::string(source.View()), text);

    mmap_remove(path);
}

TEST(FileSource, MappedAndBufferedAgreeByteForByte) {
    // Big enough to clear the Auto threshold and to be worth mapping at all.
    std::string text;
    text.reserve(1u << 21);
    for (int i = 0; i < 60000; ++i)
        text += "0.125 3.5 42\n";
    const std::string path = mmap_write_text(".txt", text);

    const FileSource buffered(path, FileSource::Mode::Buffered);
    const FileSource mapped(path, FileSource::Mode::Mmap);

    EXPECT_EQ(buffered.Size(), mapped.Size());
    EXPECT_EQ(std::string(buffered.View()), std::string(mapped.View()));

    mmap_remove(path);
}

TEST(FileSource, AutoLeavesSmallFilesBuffered) {
    const std::string path = mmap_write_text(".txt", "tiny\n");
    const FileSource source(path, FileSource::Mode::Auto);
    // Below the documented threshold the copy is cheap; mapping is not used.
    EXPECT_FALSE(source.IsMapped());
    EXPECT_EQ(std::string(source.View()), "tiny\n");
    mmap_remove(path);
}

TEST(FileSource, EmptyFileFallsBackAndIsEmpty) {
    // A zero-length file cannot be mapped; it must still read cleanly.
    const std::string path = mmap_write_text(".txt", "");
    const FileSource source(path, FileSource::Mode::Mmap);
    EXPECT_FALSE(source.IsMapped());
    EXPECT_EQ(source.Size(), 0u);
    EXPECT_NE(source.Data(), nullptr) << "Data() must never be null";
    mmap_remove(path);
}

TEST(FileSource, MissingFileThrows) {
    EXPECT_THROW(FileSource(mt::temp_path(".does-not-exist")), ReadError);
}

TEST(FileSource, PageMultipleSizedFilesAreNotMapped) {
    // A mapped buffer behaves as NUL-terminated only because the kernel
    // zero-fills the final partial page. An exact page multiple has no such
    // slack, so FileSource must decline to map it -- even when asked to.
    const std::size_t page = 4096;
    const std::string path = mmap_write_text(".bin", std::string(page * 4, 'x'));

    const FileSource source(path, FileSource::Mode::Mmap);
    EXPECT_FALSE(source.IsMapped()) << "page-multiple file must fall back to buffered";
    EXPECT_EQ(source.Size(), page * 4);

    mmap_remove(path);
}

TEST(FileSource, MoveTransfersOwnership) {
    const std::string path = mmap_write_text(".txt", "movable\n");
    FileSource a(path, FileSource::Mode::Buffered);
    FileSource b(std::move(a));
    EXPECT_EQ(std::string(b.View()), "movable\n");
    mmap_remove(path);
}

// ---------------------------------------------------------------------------
// Reader integration
// ---------------------------------------------------------------------------

TEST(MmapRead, GmshMappedEqualsGmshBuffered) {
    const Mesh source = mt::data_mesh();
    const std::string path = mt::temp_path(".msh");
    write_gmsh41(path, source, /*binary=*/false);

    ReadOptions on;
    on.mMmap = MmapMode::On;
    ReadOptions off;
    off.mMmap = MmapMode::Off;

    const Mesh mapped = read_gmsh(path, on);
    const Mesh buffered = read_gmsh(path, off);

    mt::expect_same_geometry(mapped, buffered);
    EXPECT_EQ(mapped.PointDataNames(), buffered.PointDataNames());
    EXPECT_EQ(mapped.CellDataNames(), buffered.CellDataNames());

    mmap_remove(path);
}

// The single rule that makes FileSource's function-local lifetime safe:
// nothing derived from View() may survive in the returned mesh. Deleting and
// overwriting the file must not disturb an already-read mesh.
TEST(MmapRead, MeshOutlivesTheFileItWasReadFrom) {
    const Mesh source = mt::data_mesh();
    const std::string path = mt::temp_path(".msh");
    write_gmsh41(path, source, /*binary=*/false);

    ReadOptions on;
    on.mMmap = MmapMode::On;
    const Mesh mesh = read_gmsh(path, on);

    const std::size_t num_points = mesh.NumPoints();
    const std::size_t num_blocks = mesh.NumCellBlocks();

    // Destroy the backing file, then scribble a different file over the path.
    mmap_remove(path);
    {
        std::ofstream os(path, std::ios::binary);
        os << std::string(64 * 1024, 'Z');
    }

    EXPECT_EQ(mesh.NumPoints(), num_points);
    ASSERT_EQ(mesh.NumCellBlocks(), num_blocks);
    mt::expect_same_geometry(mesh, source);

    mmap_remove(path);
}

TEST(MmapRead, DefaultOptionsStillRead) {
    const Mesh source = mt::tri_mesh();
    const std::string path = mt::temp_path(".msh");
    write_gmsh41(path, source, /*binary=*/false);

    mt::expect_same_geometry(read_gmsh(path), source);

    mmap_remove(path);
}
