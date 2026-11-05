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
#include <string>
#include <string_view>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/write_options.hpp"
#include "thread_sweep.hpp"

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

struct CpThreads {
    explicit CpThreads(int Count) : mSaved(mio_test::max_threads()) { mio_test::set_threads(Count); }
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
    mio::Mesh big;
    {
        const std::size_t n = 40000;
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
    }
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
