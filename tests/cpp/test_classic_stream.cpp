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

// Locale independence of the *stream* number path: `detail/classic_stream.hpp`.
//
// `fast_number.hpp` covers `strtod`/`snprintf`, which `setlocale` drives. A C++
// stream is driven by `std::locale::global` instead, which no `pytest` can reach
// (`locale.setlocale` moves the C locale, not the C++ one) -- hence a C++ test.
//
// The hostile locale is synthesized rather than looked up by name: a bare CI
// image ships no `de_DE`, and a test that skips wherever it cannot find one is a
// test that does not run. A `numpunct` with a comma radix and 3-digit grouping
// is also *more* hostile than a real German locale, because the grouping is
// guaranteed active. It has no name, so `std::locale::global` leaves the C
// locale alone: this isolates the stream exposure with `fast_number.hpp` out of
// the picture.

// External includes
#include <gtest/gtest.h>

// System includes
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/formats/vtu.hpp"

namespace {

using meshioplusplus::detail::imbue_classic;
using meshioplusplus::detail::make_classic_ifstream;
using meshioplusplus::detail::make_classic_istringstream;
using meshioplusplus::detail::make_classic_ofstream;
using meshioplusplus::detail::make_classic_ostringstream;
using meshioplusplus::detail::make_classic_stringstream;

/// Comma radix, dot thousands separator, groups of three: `1234567.5` -> `1.234.567,5`.
struct CsCommaGroupingPunct : std::numpunct<char> {
    char do_decimal_point() const override { return ','; }
    char do_thousands_sep() const override { return '.'; }
    std::string do_grouping() const override { return "\3"; }
};

/// RAII: installs the hostile global locale and restores the previous one.
/// Mandatory -- every gtest in this binary shares one process, and a leaked
/// global locale would corrupt every later test's streams.
class CsHostileGlobalLocale {
public:
    CsHostileGlobalLocale()
        : mSaved(
              std::locale::global(std::locale(std::locale::classic(), new CsCommaGroupingPunct))) {}
    ~CsHostileGlobalLocale() { std::locale::global(mSaved); }
    CsHostileGlobalLocale(const CsHostileGlobalLocale&) = delete;
    CsHostileGlobalLocale& operator=(const CsHostileGlobalLocale&) = delete;

private:
    std::locale mSaved;
};

std::string cs_slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

/// A triangulated `n x n` grid: 1600 points at n = 40, so vertex indices pass 1000
/// (where an integer `operator<<` grows a separator under grouping) and the
/// coordinates are non-integers (where the radix matters).
mt::Mesh cs_grid_mesh(int n = 40) {
    std::vector<std::vector<double>> pts;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            pts.push_back({0.37 * i + 0.125, 0.53 * j + 0.0625, 0.0});
    std::vector<std::vector<std::int64_t>> tris;
    for (int j = 0; j + 1 < n; ++j) {
        for (int i = 0; i + 1 < n; ++i) {
            const std::int64_t a = static_cast<std::int64_t>(j) * n + i;
            tris.push_back({a, a + 1, a + n});
            tris.push_back({a + 1, a + n + 1, a + n});
        }
    }
    return mt::make_mesh(std::move(pts), "triangle", std::move(tris));
}

}  // namespace

// --- the hostile locale really is hostile (the tests below are not vacuous) --

TEST(ClassicStream, PlainStreamsAreCorruptedByTheHostileLocale) {
    CsHostileGlobalLocale hostile;

    std::ostringstream plain_out;
    plain_out << 1234567 << ' ' << 1.5;
    EXPECT_EQ(plain_out.str(), "1.234.567 1,5");

    std::istringstream plain_in("1.5");
    double x = 0.0;
    plain_in >> x;
    EXPECT_NE(x, 1.5);
}

// --- the factories ------------------------------------------------------------

TEST(ClassicStream, FactoriesIgnoreTheGlobalLocale) {
    CsHostileGlobalLocale hostile;

    auto out = make_classic_ostringstream();
    out << 1234567 << ' ' << 1.5;
    EXPECT_EQ(out.str(), "1234567 1.5");

    auto in = make_classic_istringstream("1.5 1234567");
    double x = 0.0;
    std::int64_t n = 0;
    in >> x >> n;
    EXPECT_EQ(x, 1.5);
    EXPECT_EQ(n, 1234567);

    auto both = make_classic_stringstream("2.5");
    double y = 0.0;
    both >> y;
    EXPECT_EQ(y, 2.5);
}

TEST(ClassicStream, FactoriesAreClassicInTheDefaultConfigurationToo) {
    // The skip-if-already-classic path in imbue_classic must leave a correct stream.
    auto out = make_classic_ostringstream();
    out << 1234567 << ' ' << 1.5;
    EXPECT_EQ(out.str(), "1234567 1.5");
    EXPECT_EQ(out.getloc(), std::locale::classic());
}

TEST(ClassicStream, ImbueClassicSurvivesTheFactoryReturnMove) {
    CsHostileGlobalLocale hostile;
    const std::string path = mt::temp_path(".txt");
    {
        auto f = make_classic_ofstream(path, std::ios::binary);
        ASSERT_TRUE(f.is_open());
        EXPECT_EQ(f.getloc(), std::locale::classic());
        f << 1234567;
    }
    {
        auto f = make_classic_ifstream(path, std::ios::binary);
        ASSERT_TRUE(f.is_open());
        EXPECT_EQ(f.getloc(), std::locale::classic());
        std::int64_t n = 0;
        f >> n;
        EXPECT_EQ(n, 1234567);
    }
    EXPECT_EQ(cs_slurp(path), "1234567");
    std::remove(path.c_str());
}

TEST(ClassicStream, ImbueClassicRepairsAStreamBuiltUnderTheHostileLocale) {
    std::ostringstream plain;
    {
        CsHostileGlobalLocale hostile;
        std::ostringstream inner;
        imbue_classic(inner);
        inner << 1234567;
        EXPECT_EQ(inner.str(), "1234567");
    }
    imbue_classic(plain);
    EXPECT_EQ(plain.getloc(), std::locale::classic());
}

TEST(ClassicStream, FileStreamFailureStillSetsFailbit) {
    // The factories open() after construction rather than in the constructor;
    // the observable contract (`!stream` on failure) must be unchanged.
    auto in = make_classic_ifstream("/nonexistent-dir/definitely/missing.txt");
    EXPECT_FALSE(in.is_open());
    EXPECT_TRUE(!in);
    auto out = make_classic_ofstream("/nonexistent-dir/definitely/missing.txt", std::ios::binary);
    EXPECT_FALSE(out.is_open());
    EXPECT_TRUE(!out);
}

// --- end to end: writers and readers -------------------------------------------

namespace {

/// Write @p rMesh with @p write under the default locale and again under the
/// hostile one; the two files must be byte-identical.
template <typename WriteFn>
void cs_expect_writer_locale_independent(const mt::Mesh& rMesh, const std::string& rSuffix,
                                         WriteFn write) {
    const std::string a = mt::temp_path(rSuffix);
    const std::string b = mt::temp_path(rSuffix);
    write(a, rMesh);
    {
        CsHostileGlobalLocale hostile;
        write(b, rMesh);
    }
    const std::string bytes_a = cs_slurp(a);
    const std::string bytes_b = cs_slurp(b);
    EXPECT_FALSE(bytes_a.empty());
    // Not EXPECT_EQ on the strings: on failure gtest would dump both whole files.
    std::size_t first_diff = 0;
    while (first_diff < bytes_a.size() && first_diff < bytes_b.size() &&
           bytes_a[first_diff] == bytes_b[first_diff])
        ++first_diff;
    EXPECT_TRUE(bytes_a == bytes_b) << rSuffix << ": output changed under a comma+grouping locale"
                                    << " (first difference at byte " << first_diff << ")";
    std::remove(a.c_str());
    std::remove(b.c_str());
}

}  // namespace

TEST(ClassicStream, ObjWriterIsLocaleIndependent) {
    cs_expect_writer_locale_independent(
        cs_grid_mesh(), ".obj",
        [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_obj(p, m); });
}

TEST(ClassicStream, OffWriterIsLocaleIndependent) {
    cs_expect_writer_locale_independent(
        cs_grid_mesh(), ".off",
        [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_off(p, m); });
}

TEST(ClassicStream, MeditWriterIsLocaleIndependent) {
    cs_expect_writer_locale_independent(
        cs_grid_mesh(), ".mesh",
        [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_medit_ascii(p, m); });
}

TEST(ClassicStream, VtuAsciiWriterIsLocaleIndependent) {
    cs_expect_writer_locale_independent(
        cs_grid_mesh(), ".vtu", [](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_vtu(p, m, /*binary=*/false, /*zlib=*/false);
        });
}

TEST(ClassicStream, ReadersRoundTripUnderTheHostileLocale) {
    // The read half: `off.cpp`'s `in >> coordinate` and every `istringstream >> x`
    // fail here without the factories.
    const mt::Mesh mesh = cs_grid_mesh();
    const std::string off = mt::temp_path(".off");
    const std::string obj = mt::temp_path(".obj");
    meshioplusplus::write_off(off, mesh);
    meshioplusplus::write_obj(obj, mesh);
    {
        CsHostileGlobalLocale hostile;
        mt::expect_mesh_eq(mesh, meshioplusplus::read_off(off));
        mt::expect_mesh_eq(mesh, meshioplusplus::read_obj(obj));
    }
    std::remove(off.c_str());
    std::remove(obj.c_str());
}
