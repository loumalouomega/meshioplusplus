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

// The chunked parallel number parse (roadmap 5.7.2): equal to the serial
// prefix-parser loop at every thread count, refusing what that loop would not
// take whole, and the readers built on it equal to their serial selves.

// System includes
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/write_options.hpp"
#include "thread_sweep.hpp"

#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "../../src/cpp/src/detail/chunked_parse.hpp"

namespace {

namespace mio = meshioplusplus;
namespace det = meshioplusplus::detail;

// Numbers in many spellings, separated by every classic blank.
std::string cp_text(std::size_t n, std::uint64_t Seed) {
    static const char* kBlanks[] = {" ", "  ", "\t", "\n", "\r\n", " \v ", "\f"};
    std::string out = "  \n";
    std::uint64_t s = Seed;
    char buf[64];
    for (std::size_t i = 0; i < n; ++i) {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        const double v = static_cast<double>(static_cast<std::int64_t>(s % 2000001) - 1000000) /
                         static_cast<double>(1 + (s >> 40) % 997);
        switch ((s >> 20) % 6) {
            case 0:
                std::snprintf(buf, sizeof buf, "%.17g", v);
                break;
            case 1:
                std::snprintf(buf, sizeof buf, "%.3e", v);
                break;
            case 2:
                std::snprintf(buf, sizeof buf, "+%.6f", v < 0 ? -v : v);
                break;
            case 3:
                std::snprintf(buf, sizeof buf, "%lld", static_cast<long long>(v));
                break;
            case 4:
                std::snprintf(buf, sizeof buf, "%a", v);  // hexadecimal: the fallback parser
                break;
            default:
                std::snprintf(buf, sizeof buf, "%.40f", v);  // a long token
                break;
        }
        out += buf;
        out += kBlanks[(s >> 30) % 7];
    }
    return out + "\nTAIL 1 2 3";
}

// The serial loop the readers keep: DoublePrefix value after value.
std::vector<double> cp_serial(std::string_view Text, std::size_t Count) {
    det::TextCursor tc(Text);
    std::vector<double> out(Count);
    for (double& v : out)
        EXPECT_TRUE(tc.DoublePrefix(v));
    return out;
}

std::size_t cp_parallel(std::string_view Text, std::size_t Count, std::vector<double>& rOut) {
    rOut.assign(Count, 0.0);
    return det::parse_tokens(Text, 0, Count, [&](std::string_view Token, std::size_t i) {
        det::TextCursor tc(Token);
        double x = 0.0;
        if (!tc.DoublePrefix(x) || tc.Pos() != Token.size())
            return false;
        rOut[i] = x;
        return true;
    });
}

// A mesh made of `n` points and `n - 3` tetrahedra over consecutive points,
// with a point field: big enough that its ASCII files take the parallel paths.
mio::Mesh cp_big_mesh(std::size_t n) {
    mio::Mesh big;
    mio::NDArray pts(mio::DType::Float64, {n, std::size_t{3}});
    for (std::size_t i = 0; i < 3 * n; ++i)
        pts.As<double>()[i] = 0.1 * static_cast<double>(i % 977) - 1.0 / (1.0 + i % 13);
    big.AssignPoints(std::move(pts));
    mio::NDArray conn(mio::DType::Int64, {n - 3, std::size_t{4}});
    for (std::size_t c = 0; c + 3 < n; ++c)
        for (std::size_t k = 0; k < 4; ++k)
            conn.As<std::int64_t>()[4 * c + k] = static_cast<std::int64_t>(c + k);
    big.AddCellBlock("tetra", std::move(conn));
    mio::NDArray u(mio::DType::Float64, {n});
    for (std::size_t i = 0; i < n; ++i)
        u.As<double>()[i] = 1e-3 * static_cast<double>(i) - 7.25;
    big.AddPointData("u", std::move(u));
    return big;
}

std::string cp_slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), {});
}

void cp_spit(const std::string& rPath, const std::string& rText) {
    std::ofstream(rPath, std::ios::binary) << rText;
}

std::uint64_t cp_digest(const mio::Mesh& rMesh) {
    mio::bench::MeshDigest d;
    d.Of(rMesh);
    return d.Value();
}

struct CpThreads {
    explicit CpThreads(int Count) : mSaved(mio_test::max_threads()) {
        mio_test::set_threads(Count);
    }
    ~CpThreads() { mio_test::set_threads(mSaved); }
    int mSaved;
};

}  // namespace

TEST(ChunkedParse, EqualsTheSerialLoopAtEveryThreadCount) {
    const std::size_t n = 300000;
    const std::string text = cp_text(n, 42);
    const std::vector<double> want = cp_serial(text, n);
    EXPECT_EQ(det::count_tokens(text), n + 4);  // and "TAIL 1 2 3"
    for (int threads : {1, 2, 4, 7}) {
        CpThreads guard(threads);
        std::vector<double> got;
        const std::size_t end = cp_parallel(text, n, got);
        ASSERT_NE(end, std::string_view::npos) << threads;
        // Just past the last number: only blanks, then the tail.
        const std::size_t tail = text.find("TAIL", end);
        EXPECT_TRUE(end > 0 && !det::text_is_blank(text[end - 1])) << threads;
        for (std::size_t i = end; i < tail; ++i)
            ASSERT_TRUE(det::text_is_blank(text[i])) << threads;
        ASSERT_EQ(got.size(), want.size());
        for (std::size_t i = 0; i < n; ++i)
            ASSERT_EQ(std::memcmp(&got[i], &want[i], sizeof(double)), 0) << i << " @" << threads;
    }
}

TEST(ChunkedParse, RefusesWhatTheSerialLoopWouldNotTakeWhole) {
    CpThreads guard(4);
    std::string text = cp_text(200000, 7);
    // A glued pair: the serial loop reads 1.5 then -2.5; a token parse
    // cannot, so the whole parse is refused and the caller goes serial.
    const std::size_t at = text.size() / 2;
    const std::size_t blank = text.find(' ', at);
    text.insert(blank + 1, "1.5-2.5 ");
    std::vector<double> got;
    EXPECT_EQ(cp_parallel(text, 200001, got), std::string_view::npos);
    // Too few tokens: refused too.
    EXPECT_EQ(cp_parallel(cp_text(1000, 3), 5000, got), std::string_view::npos);
    // Nothing to parse is the start itself.
    EXPECT_EQ(cp_parallel(text, 0, got), 0u);
}

// The readers that parse in parallel give what their serial loop gives: the
// same mesh at one thread and at four, from files large enough to take the
// parallel path.
TEST(ChunkedParse, AsciiReadersEqualTheirSerialSelves) {
    const mio::Mesh big = cp_big_mesh(40000);
    struct Case {
        const char* mFormat;
        const char* mExt;
    };
    for (const Case c : {Case{"vtk", ".vtk"}, Case{"vtu", ".vtu"}, Case{"gmsh", ".msh"}}) {
        const std::string path = mt::temp_path(c.mExt);
        mio::WriteOptions o;
        o.mEncoding = mio::WriteEncoding::Ascii;
        mio::registry_write_ex(path, big, c.mFormat, o);
        mio::Mesh serial, parallel;
        {
            CpThreads guard(1);
            serial = mio::registry_read(path, c.mFormat, mio::ReadOptions{});
        }
        {
            CpThreads guard(4);
            parallel = mio::registry_read(path, c.mFormat, mio::ReadOptions{});
        }
        mt::expect_same_geometry(parallel, serial);
        EXPECT_EQ(cp_digest(parallel), cp_digest(serial)) << c.mFormat;
        ASSERT_EQ(parallel.NumPoints(), big.NumPoints()) << c.mFormat;
        std::filesystem::remove(path);
    }
    const std::string path = mt::temp_path(".msh");
    mio::write_gmsh22(path, big, /*binary=*/false);
    mio::Mesh serial, parallel;
    {
        CpThreads guard(1);
        serial = mio::read_gmsh(path);
    }
    {
        CpThreads guard(4);
        parallel = mio::read_gmsh(path);
    }
    mt::expect_same_geometry(parallel, serial);
    std::filesystem::remove(path);
}

namespace {

// Reads `rPath` with `rRead` at one thread and at four, and expects the same
// mesh, data and all; returns the serial read.
template <class Read>
mio::Mesh cp_same_at_one_and_four(const std::string& rPath, Read&& rRead) {
    mio::Mesh serial, parallel;
    {
        CpThreads guard(1);
        serial = rRead(rPath);
    }
    {
        CpThreads guard(4);
        parallel = rRead(rPath);
    }
    mt::expect_same_geometry(parallel, serial);
    EXPECT_EQ(cp_digest(parallel), cp_digest(serial)) << rPath;
    return serial;
}

}  // namespace

// Whole lines, about a chunk each, covering the text: every cut but the last
// lands just after a newline, and the chunks' lines are the text's lines.
TEST(ChunkedParse, ChunkLinesCutsAfterNewlines) {
    std::string text;
    for (int i = 0; i < 200000; ++i)
        text += std::to_string(i) + (i % 3 ? " x y\n" : "\r\n");
    text += "no newline at the end";
    const std::vector<std::size_t> cuts = det::chunk_lines(text, 3, text.size());
    ASSERT_GE(cuts.size(), 3u);
    EXPECT_EQ(cuts.front(), 3u);
    EXPECT_EQ(cuts.back(), text.size());
    std::vector<std::string_view> joined;
    for (std::size_t c = 0; c + 1 < cuts.size(); ++c) {
        if (c + 2 < cuts.size())
            EXPECT_EQ(text[cuts[c + 1] - 1], '\n');
        for (std::string_view l :
             det::split_lines(std::string_view(text).substr(cuts[c], cuts[c + 1] - cuts[c])))
            joined.push_back(l);
    }
    EXPECT_EQ(joined, det::split_lines(std::string_view(text).substr(3)));
    EXPECT_EQ(det::chunk_lines(text, 7, 7), std::vector<std::size_t>{7});
}

// Medit's tables, parsed in chunks, give the serial read; and its comments are
// not blanks to the chunked parse: a `#` anywhere in a table sends that table
// to the serial loop, which skips the comment and reads the same mesh.
TEST(ChunkedParse, MeditCommentsInsideTablesFallBack) {
    const std::string path = mt::temp_path(".mesh");
    mio::write_medit_ascii(path, cp_big_mesh(40000));
    const mio::Mesh plain = cp_same_at_one_and_four(path, [](const std::string& rP) {
        return mio::read_medit_ascii(rP);
    });
    std::string text = cp_slurp(path), out;
    std::size_t line = 0;
    for (std::size_t p = 0; p < text.size();) {
        std::size_t e = text.find('\n', p);
        e = e == std::string::npos ? text.size() : e + 1;
        out.append(text, p, e - p);
        if (++line % 7919 == 0)
            out += "# a comment 1 2 3\n";
        if (line == 20011 && out.size() > 1 && out[out.size() - 1] == '\n') {
            out.pop_back();
            out += " #trailing 4 5\n";
        }
        p = e;
    }
    cp_spit(path, out);
    const mio::Mesh m = cp_same_at_one_and_four(path, [](const std::string& rP) {
        return mio::read_medit_ascii(rP);
    });
    EXPECT_EQ(m.NumPoints(), 40000u);
    EXPECT_EQ(cp_digest(m), cp_digest(plain));
    std::filesystem::remove(path);
}

// Gmsh 2.2's elements: blocks split where the type or the tag count changes,
// CRLF and indented lines, an extra token after a row and a row whose extra
// token is not an integer (the chunked pass refuses the file; the serial loop
// reads it) -- the same mesh at one thread and at four either way.
TEST(ChunkedParse, Gmsh22ElementBlocksEqualTheSerialLoop) {
    const std::size_t n = 20000;
    for (const bool odd_token : {false, true}) {
        std::string t = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n" + std::to_string(n) + "\n";
        for (std::size_t i = 0; i < n; ++i)
            t += std::to_string(i + 1) + " " + std::to_string(0.5 * i) + " " +
                 std::to_string(i % 17) + " " + std::to_string(-0.25 * (i % 5)) + "\n";
        t += "$EndNodes\n$Elements\n";
        std::string rows;
        std::size_t count = 0;
        for (std::size_t c = 0; c + 3 < n; ++c) {
            const std::size_t run = (c / 2500) % 4;
            std::string row = std::to_string(++count);
            if (run == 0)  // tetra, two tags
                row += " 4 2 7 " + std::to_string(c % 9);
            else if (run == 1)  // tetra, three tags
                row += " 4 3 7 1 " + std::to_string(c % 4);
            else if (run == 2)  // triangle, two tags
                row += " 2 2 3 3";
            else  // line, no tags
                row += " 1 0";
            const std::size_t nodes = run < 2 ? 4 : run == 2 ? 3 : 2;
            for (std::size_t k = 0; k < nodes; ++k)
                row += " " + std::to_string(c + k + 1);
            if (c % 1013 == 0)
                row += " 99";  // ignored by the serial loop
            if (odd_token && c == 7777)
                row += " 1.5";  // the serial loop stops at it; the chunked one refuses
            rows += (c % 3 == 0 ? "  " : "") + row + (c % 5 == 0 ? "\r\n" : "\n");
        }
        t += std::to_string(count) + "\n" + rows + "$EndElements\n";
        const std::string path = mt::temp_path(".msh");
        cp_spit(path, t);
        const mio::Mesh m =
            cp_same_at_one_and_four(path, [](const std::string& rP) { return mio::read_gmsh(rP); });
        EXPECT_EQ(m.NumPoints(), n) << odd_token;
        std::filesystem::remove(path);
    }
}

// Gmsh 4.0's ASCII element blocks (int32 tags) on the chunked path.
TEST(ChunkedParse, Gmsh40ElementBlocksEqualTheSerialLoop) {
    const std::size_t n = 40000;
    std::string t = "$MeshFormat\n4.0 0 8\n$EndMeshFormat\n$Entities\n0 0 0 1\n"
                    "1 0 0 0 1 1 1 0 0\n$EndEntities\n$Nodes\n1 " +
                    std::to_string(n) + "\n1 3 0 " + std::to_string(n) + "\n";
    for (std::size_t i = 0; i < n; ++i)
        t += std::to_string(i + 1) + " " + std::to_string(0.5 * i) + " " +
             std::to_string(i % 17) + " 0\n";
    t += "$EndNodes\n$Elements\n1 " + std::to_string(n - 3) + "\n1 3 4 " +
         std::to_string(n - 3) + "\n";
    for (std::size_t c = 0; c + 3 < n; ++c)
        t += std::to_string(c + 1) + " " + std::to_string(c + 1) + " " + std::to_string(c + 2) +
             " " + std::to_string(c + 3) + " " + std::to_string(c + 4) + "\n";
    t += "$EndElements\n";
    const std::string path = mt::temp_path(".msh");
    cp_spit(path, t);
    const mio::Mesh m =
        cp_same_at_one_and_four(path, [](const std::string& rP) { return mio::read_gmsh(rP); });
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).NumCells(), n - 3);
    // A node tag below 1 is still the serial loop's error.
    const std::string first_row = "\n1 1 2 3 4\n";
    std::string bad = t;
    bad.replace(bad.find(first_row), first_row.size(), "\n1 0 2 3 4\n");
    cp_spit(path, bad);
    for (const int threads : {1, 4}) {
        CpThreads guard(threads);
        EXPECT_THROW(mio::read_gmsh(path), mio::ReadError) << threads;
    }
    std::filesystem::remove(path);
}

// ASCII STL read in chunks of whole lines gives the serial read's mesh.
TEST(ChunkedParse, AsciiStlEqualsItsSerialSelf) {
    const std::size_t n = 40000;
    mio::Mesh surf;
    {
        mio::NDArray pts(mio::DType::Float64, {n, std::size_t{3}});
        for (std::size_t i = 0; i < 3 * n; ++i)
            pts.As<double>()[i] = 0.1 * static_cast<double>(i % 977) - 1.0 / (1.0 + i % 13);
        surf.AssignPoints(std::move(pts));
        mio::NDArray conn(mio::DType::Int64, {n - 2, std::size_t{3}});
        for (std::size_t c = 0; c + 2 < n; ++c)
            for (std::size_t k = 0; k < 3; ++k)
                conn.As<std::int64_t>()[3 * c + k] = static_cast<std::int64_t>(c + k);
        surf.AddCellBlock("triangle", std::move(conn));
    }
    const std::string path = mt::temp_path(".stl");
    mio::WriteOptions o;
    o.mEncoding = mio::WriteEncoding::Ascii;
    mio::registry_write_ex(path, surf, "stl", o);
    ASSERT_GT(std::filesystem::file_size(path), std::size_t{1} << 20);
    const mio::Mesh m =
        cp_same_at_one_and_four(path, [](const std::string& rP) { return mio::read_stl(rP); });
    EXPECT_EQ(m.Cells(0).NumCells(), n - 2);
    // A truncated file is malformed at every thread count.
    std::string text = cp_slurp(path);
    text.resize(text.rfind("vertex"));
    cp_spit(path, text);
    for (const int threads : {1, 4}) {
        CpThreads guard(threads);
        EXPECT_THROW(mio::read_stl(path), mio::ReadError) << threads;
    }
    std::filesystem::remove(path);
}
