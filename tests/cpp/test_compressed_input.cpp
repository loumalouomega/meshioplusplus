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
// Compressed input files (roadmap 7.1.1): the codecs, the single-member zip
// reader, and registry reads of `cube.msh.gz` and friends.

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
#include <zlib.h>
#endif

#include "meshioplusplus/detail/compressed_input.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/registry.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;
using detail::Compression;

namespace {

std::string cmp_slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void cmp_spit(const std::string& rPath, const std::string& rBytes) {
    std::ofstream out(rPath, std::ios::binary);
    out << rBytes;
}

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
// Deflate @p rIn: gzip wrapper (WindowBits 31) or raw (-15, as zip stores it).
std::string cmp_deflate(const std::string& rIn, int WindowBits) {
    z_stream z{};
    deflateInit2(&z, Z_BEST_COMPRESSION, Z_DEFLATED, WindowBits, 8, Z_DEFAULT_STRATEGY);
    std::string out(deflateBound(&z, static_cast<uLong>(rIn.size())) + 32, '\0');
    z.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(rIn.data()));
    z.avail_in = static_cast<uInt>(rIn.size());
    z.next_out = reinterpret_cast<Bytef*>(&out[0]);
    z.avail_out = static_cast<uInt>(out.size());
    deflate(&z, Z_FINISH);
    out.resize(z.total_out);
    deflateEnd(&z);
    return out;
}
#endif

void cmp_put16(std::string& rOut, std::uint32_t V) {
    rOut += static_cast<char>(V & 0xFF);
    rOut += static_cast<char>((V >> 8) & 0xFF);
}

void cmp_put32(std::string& rOut, std::uint32_t V) {
    cmp_put16(rOut, V & 0xFFFF);
    cmp_put16(rOut, V >> 16);
}

// A zip archive of @p rMembers (name, stored data, method, raw size). The CRC
// is left zero: the reader does not check it.
struct CmpMember {
    std::string mName;
    std::string mData;
    std::uint32_t mMethod;
    std::uint32_t mSize;
};

std::string cmp_zip(const std::vector<CmpMember>& rMembers) {
    std::string out, central;
    for (const CmpMember& m : rMembers) {
        const std::uint32_t offset = static_cast<std::uint32_t>(out.size());
        cmp_put32(out, 0x04034b50u);
        cmp_put16(out, 20);
        cmp_put16(out, 0);
        cmp_put16(out, m.mMethod);
        cmp_put32(out, 0);  // time, date
        cmp_put32(out, 0);  // crc
        cmp_put32(out, static_cast<std::uint32_t>(m.mData.size()));
        cmp_put32(out, m.mSize);
        cmp_put16(out, static_cast<std::uint32_t>(m.mName.size()));
        cmp_put16(out, 0);
        out += m.mName;
        out += m.mData;
        cmp_put32(central, 0x02014b50u);
        cmp_put16(central, 20);
        cmp_put16(central, 20);
        cmp_put16(central, 0);
        cmp_put16(central, m.mMethod);
        cmp_put32(central, 0);
        cmp_put32(central, 0);
        cmp_put32(central, static_cast<std::uint32_t>(m.mData.size()));
        cmp_put32(central, m.mSize);
        cmp_put16(central, static_cast<std::uint32_t>(m.mName.size()));
        cmp_put16(central, 0);
        cmp_put16(central, 0);
        cmp_put16(central, 0);
        cmp_put16(central, 0);
        cmp_put32(central, 0);
        cmp_put32(central, offset);
        central += m.mName;
    }
    const std::uint32_t cd_offset = static_cast<std::uint32_t>(out.size());
    out += central;
    cmp_put32(out, 0x06054b50u);
    cmp_put16(out, 0);
    cmp_put16(out, 0);
    cmp_put16(out, static_cast<std::uint32_t>(rMembers.size()));
    cmp_put16(out, static_cast<std::uint32_t>(rMembers.size()));
    cmp_put32(out, static_cast<std::uint32_t>(central.size()));
    cmp_put32(out, cd_offset);
    cmp_put16(out, 0);
    return out;
}

// A fresh directory holding `cube.msh`, the Gmsh file of the tet fixture.
struct CmpDir {
    std::filesystem::path mDir;
    std::string mMsh;
    CmpDir() {
        mDir = std::filesystem::path(mt::temp_path("_cmp"));
        std::filesystem::create_directories(mDir);
        const std::string path = (mDir / "cube.msh").string();
        registry_writers().at("gmsh")(path, mt::tet_mesh());
        mMsh = cmp_slurp(path);
    }
    ~CmpDir() {
        std::error_code ec;
        std::filesystem::remove_all(mDir, ec);
    }
    std::string Path(const std::string& rName) const { return (mDir / rName).string(); }
};

// How many of the decompressed copies' directories exist right now.
std::size_t cmp_temp_dirs() {
    std::size_t n = 0;
    for (const auto& entry :
         std::filesystem::directory_iterator(std::filesystem::temp_directory_path()))
        if (entry.path().filename().string().rfind("meshioplusplus_decompress_", 0) == 0)
            ++n;
    return n;
}

void cmp_expect_tets(const Mesh& rMesh) {
    mt::expect_same_geometry(mt::tet_mesh(), rMesh);
}

}  // namespace

TEST(CompressedInput, SuffixesAndNames) {
    EXPECT_EQ(detail::compression_from_suffix("a/cube.msh.gz"), Compression::Gzip);
    EXPECT_EQ(detail::compression_from_suffix("CUBE.VTU.BZ2"), Compression::Bzip2);
    EXPECT_EQ(detail::compression_from_suffix("x.inp.xz"), Compression::Xz);
    EXPECT_EQ(detail::compression_from_suffix("x.stl.zst"), Compression::Zstd);
    EXPECT_EQ(detail::compression_from_suffix("case.zip"), Compression::Zip);
    EXPECT_EQ(detail::compression_from_suffix("cube.msh"), Compression::None);
    EXPECT_EQ(detail::compression_from_suffix(".gz"), Compression::None);
    EXPECT_EQ(detail::decompressed_name("dir/cube.msh.gz", Compression::Gzip), "dir/cube.msh");
    EXPECT_TRUE(detail::compression_available(Compression::None));
    EXPECT_STREQ(detail::compression_option(Compression::Zstd), "MESHIOPLUSPLUS_WITH_ZSTD");
}

TEST(CompressedInput, ResolveFormatStripsTheSuffix) {
    EXPECT_EQ(resolve_format("cube.vtu.gz", ""), "vtu");
    EXPECT_EQ(resolve_format("part.stl.bz2", ""), "stl");
    EXPECT_EQ(resolve_format("run.inp.xz", ""), "abaqus");
    // A compound extension a format owns is not stripped.
    EXPECT_EQ(resolve_format("mesh.vol.gz", ""), "netgen");
    EXPECT_EQ(resolve_format("mesh.xda.gz", ""), "libmesh");
    EXPECT_THROW(resolve_format("data.gz", ""), ReadError);
}

TEST(CompressedInput, StoredZip) {
    CmpDir d;
    const std::string zip =
        cmp_zip({{"inner/cube.msh", d.mMsh, 0, static_cast<std::uint32_t>(d.mMsh.size())}});
    EXPECT_EQ(detail::zip_member_name(zip), "inner/cube.msh");
    EXPECT_EQ(detail::decompress_bytes(Compression::Zip, zip), d.mMsh);
    cmp_spit(d.Path("case.zip"), zip);
    EXPECT_EQ(detail::decompressed_name(d.Path("case.zip"), Compression::Zip), d.Path("cube.msh"));
    EXPECT_EQ(resolve_format(d.Path("case.zip"), ""), "gmsh");
    cmp_expect_tets(registry_read(d.Path("case.zip"), "gmsh", {}));
    // Two files, a directory entry ignored, an unknown method: refused.
    const std::string two = cmp_zip({{"a.msh", d.mMsh, 0, 1}, {"b.msh", d.mMsh, 0, 1}});
    EXPECT_THROW(detail::decompress_bytes(Compression::Zip, two), ReadError);
    const std::string dir_and_file =
        cmp_zip({{"inner/", "", 0, 0},
                 {"inner/cube.msh", d.mMsh, 0, static_cast<std::uint32_t>(d.mMsh.size())}});
    EXPECT_EQ(detail::decompress_bytes(Compression::Zip, dir_and_file), d.mMsh);
    EXPECT_THROW(detail::decompress_bytes(Compression::Zip, cmp_zip({{"a", "x", 12, 1}})),
                 ReadError);
    EXPECT_THROW(detail::decompress_bytes(Compression::Zip, "not a zip at all, no"), ReadError);
}

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
TEST(CompressedInput, GzipAndDeflatedZipThroughTheRegistry) {
    CmpDir d;
    const std::size_t before = cmp_temp_dirs();
    const std::string gz = cmp_deflate(d.mMsh, 31);
    EXPECT_EQ(detail::decompress_bytes(Compression::Gzip, gz), d.mMsh);
    // Concatenated members decompress to the concatenation, as `gzip -d` does.
    EXPECT_EQ(detail::decompress_bytes(Compression::Gzip, gz + cmp_deflate("tail", 31)),
              d.mMsh + "tail");
    EXPECT_THROW(detail::decompress_bytes(Compression::Gzip, gz.substr(0, gz.size() / 2)),
                 ReadError);

    cmp_spit(d.Path("cube.msh.gz"), gz);
    cmp_expect_tets(
        registry_read(d.Path("cube.msh.gz"), resolve_format(d.Path("cube.msh.gz"), ""), {}));
    // Every plain reader of the registry takes the compressed path too.
    cmp_expect_tets(registry_readers().at("gmsh")(d.Path("cube.msh.gz")));
    const MeshMetadata meta = registry_read_metadata(d.Path("cube.msh.gz"), "gmsh", {});
    EXPECT_EQ(meta.mNumPoints, 5u);

    const std::string zip = cmp_zip(
        {{"cube.msh", cmp_deflate(d.mMsh, -15), 8, static_cast<std::uint32_t>(d.mMsh.size())}});
    cmp_spit(d.Path("cube.msh.zip"), zip);
    cmp_expect_tets(
        registry_read(d.Path("cube.msh.zip"), resolve_format(d.Path("cube.msh.zip"), ""), {}));

    // The temporary copies are gone afterwards.
    EXPECT_EQ(cmp_temp_dirs(), before);
}

TEST(CompressedInput, MultiFileFormatsAreRefusedByName) {
    CmpDir d;
    cmp_spit(d.Path("case.pvtu.gz"), cmp_deflate("<VTKFile/>", 31));
    try {
        registry_read(d.Path("case.pvtu.gz"), "pvtu", {});
        FAIL() << "a compressed pvtu was read";
    } catch (const ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("pvtu"), std::string::npos) << e.what();
        EXPECT_NE(std::string(e.what()).find("several files"), std::string::npos) << e.what();
    }
}
#endif

TEST(CompressedInput, MissingCodecsNameTheirOption) {
    for (const Compression kind : {Compression::Bzip2, Compression::Xz, Compression::Zstd}) {
        if (detail::compression_available(kind))
            continue;
        try {
            detail::decompress_bytes(kind, "whatever");
            FAIL() << "decompressed without the codec";
        } catch (const ReadError& e) {
            EXPECT_NE(std::string(e.what()).find(detail::compression_option(kind)),
                      std::string::npos)
                << e.what();
        }
    }
}
