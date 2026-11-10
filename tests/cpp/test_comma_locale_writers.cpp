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

// Every registered writer writes the same bytes under a comma-decimal
// LC_NUMERIC as under "C" (roadmap 5.7.1): number formatting goes through
// `snprintf_c`/`CNumber`, never the process locale. The ASCII flavour of the
// formats that have one is written too.

// System includes
#include <clocale>
#include <filesystem>
#include <functional>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "comma_locale.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/write_options.hpp"

namespace {

namespace mio = meshioplusplus;
namespace fs = std::filesystem;

// Tetrahedra with fractional coordinates and point and cell fields, so every
// writer has numbers with a fractional part to format.
mio::Mesh clw_mesh() {
    mio::Mesh m = mt::make_mesh({{0.125, 0, 0},
                                 {1.5, 0.25, 0},
                                 {1.25, 1.75, 0.0625},
                                 {0, 1.5, 0},
                                 {0.5, 0.5, 0.75},
                                 {0.625, 0.375, -1.25}},
                                "tetra", {{0, 1, 2, 4}, {0, 2, 3, 4}, {0, 2, 1, 5}});
    mio::NDArray u(mio::DType::Float64, {6});
    for (std::size_t i = 0; i < 6; ++i)
        u.As<double>()[i] = 0.1 * static_cast<double>(i) - 0.35;
    m.AddPointData("u", std::move(u));
    mio::NDArray c(mio::DType::Float64, {3});
    for (std::size_t i = 0; i < 3; ++i)
        c.As<double>()[i] = 2.5 + 0.3 * static_cast<double>(i);
    m.AddCellData("c", {std::move(c)});
    return m;
}

// Every file under `rDir`, by relative path.
std::map<std::string, std::string> clw_files(const fs::path& rDir) {
    std::map<std::string, std::string> out;
    for (const auto& e : fs::recursive_directory_iterator(rDir)) {
        if (!e.is_regular_file())
            continue;
        std::ifstream in(e.path(), std::ios::binary);
        std::ostringstream ss;
        ss << in.rdbuf();
        out[fs::relative(e.path(), rDir).string()] = ss.str();
    }
    return out;
}

using ClwWrite = std::function<void(const std::string&)>;

// Writes through `rWrite` into a fresh directory; empty when the writer
// declines the mesh.
std::map<std::string, std::string> clw_write(const ClwWrite& rWrite, const std::string& rName) {
    const fs::path dir = fs::temp_directory_path() / ("mio_clw_" + rName);
    fs::remove_all(dir);
    fs::create_directories(dir);
    std::map<std::string, std::string> out;
    // A fresh per-thread provenance record, as a new write starts with: the
    // PLY writer notes its int64 cast after rendering its header, so without
    // this the note of one write surfaces in the next one's header.
    mio::detail::provenance_begin_write();
    try {
        rWrite((dir / ("mesh_" + rName)).string());
        out = clw_files(dir);
    } catch (const std::exception&) {
        out.clear();
    }
    fs::remove_all(dir);
    return out;
}

// A file extension the registry maps to `rFormat` (the first in sorted
// order), or ".dat".
std::string clw_extension(const std::string& rFormat) {
    for (const auto& [ext, format] : mio::registry_extension_defaults())
        if (format == rFormat)
            return ext;
    return ".dat";
}

}  // namespace

TEST(CommaLocaleWriters, EveryWriterIgnoresTheProcessLocale) {
    mio_test::LcNumericGuard guard;
    if (mio_test::find_comma_locale().empty())
        GTEST_SKIP() << "no comma-decimal locale available on this machine";
    std::setlocale(LC_NUMERIC, "C");
    const mio::Mesh mesh = clw_mesh();
    std::vector<std::pair<std::string, ClwWrite>> writers;
    for (const auto& [name, fn] : mio::registry_writers()) {
        const std::string ext = clw_extension(name);
        writers.emplace_back(name,
                             [&mesh, fn, ext](const std::string& rPath) { fn(rPath + ext, mesh); });
        writers.emplace_back(name + "_ascii", [&mesh, name, ext](const std::string& rPath) {
            mio::WriteOptions o;
            o.mEncoding = mio::WriteEncoding::Ascii;
            mio::registry_write_ex(rPath + ext, mesh, name, o);
        });
    }
    writers.emplace_back("gmsh22_ascii", [&](const std::string& rPath) {
        mio::write_gmsh22(rPath + ".msh", mesh, /*binary=*/false);
    });
    std::size_t compared = 0;
    for (const auto& [name, write] : writers) {
        std::setlocale(LC_NUMERIC, "C");
        const auto in_c = clw_write(write, name);
        const std::string comma = mio_test::find_comma_locale();
        const auto in_comma = clw_write(write, name);
        std::setlocale(LC_NUMERIC, "C");
        if (in_c.empty() && in_comma.empty())
            continue;
        ++compared;
        SCOPED_TRACE("writer " + name + " under LC_NUMERIC=" + comma);
        if (in_c.size() != in_comma.size()) {
            // clw_write swallows a throw, so a writer that fails only under the
            // comma locale shows up as a size mismatch: report why it threw.
            std::setlocale(LC_NUMERIC, comma.c_str());
            std::string why = "no exception";
            try {
                write((fs::temp_directory_path() / ("mio_clw_err_" + name)).string());
            } catch (const std::exception& rErr) {
                why = rErr.what();
            }
            std::setlocale(LC_NUMERIC, "C");
            FAIL() << name << ": " << in_c.size() << " file(s) under C, " << in_comma.size()
                   << " under the comma locale; the comma write says: " << why;
        }
        for (const auto& [file, bytes] : in_c) {
            const auto it = in_comma.find(file);
            ASSERT_NE(it, in_comma.end()) << name << ": " << file << " missing under comma";
            if (it->second == bytes)
                continue;
            std::size_t at = 0;
            while (at < bytes.size() && at < it->second.size() && bytes[at] == it->second[at])
                ++at;
            ADD_FAILURE() << name << ": " << file << " differs at byte " << at << " (C: \""
                          << bytes.substr(at, 24) << "\", comma: \"" << it->second.substr(at, 24)
                          << "\")";
        }
    }
    EXPECT_GT(compared, 20u);
}
