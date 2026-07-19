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
 * @file main.cpp
 * @brief Standalone native command-line interface for meshio++.
 *
 * A Python-free CLI over `meshioplusplus_core_obj`: it mirrors the verb set of
 * the Python console script (`src/meshioplusplus/_cli/`) but links only the C++
 * core, so it builds and runs without Python or pybind11. Format dispatch reuses
 * the shared registry (`registry.hpp`); the per-format binary/ASCII/compression
 * variants call the format writer free functions directly (the registry bakes a
 * single default per format and cannot express those flags); operation verbs use
 * the backend-agnostic operations layer (`operations/*.hpp`).
 *
 * Point/cell *sets* and the `convert -s/-d` sets<->data conversions never cross
 * into the C++ core (they live only in the Python `Mesh`), so those are reported
 * as unsupported here -- matching the documented C-API/WASM limitations.
 */

// System includes
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/operations/diff.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/crop.hpp"
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/operations/stats.hpp"
// Per-format writers for the ASCII/binary/compress variants (the registry bakes
// one default each, so these are called directly).
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/ply.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/formats/xdmf.hpp"
#ifdef MESHIOPLUSPLUS_HAS_HDF5
#include "meshioplusplus/formats/cgns.hpp"
#include "meshioplusplus/formats/h5m.hpp"
#endif

#ifndef MESHIOPLUSPLUS_CLI_VERSION
#define MESHIOPLUSPLUS_CLI_VERSION "unknown"
#endif

namespace {

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::ReadError;
using meshioplusplus::WriteError;

// --------------------------------------------------------------------------
// Small argument parser
// --------------------------------------------------------------------------

/// One option/flag spec: a canonical long name, its aliases, and whether it
/// consumes a following value.
struct cli_opt_spec {
    std::string canonical;
    std::vector<std::string> aliases;
    bool takes_value;
};

/// Parsed result of a verb's argument list.
struct cli_parsed {
    std::vector<std::string> positionals;
    std::unordered_map<std::string, std::string> values;  ///< canonical -> value
    std::unordered_set<std::string> flags;                ///< canonical present
};

/// Parse `args` against `specs`. Throws std::runtime_error on an unknown option
/// or a value-option missing its argument. `--` stops option parsing.
cli_parsed cli_parse(const std::vector<std::string>& rArgs,
                     const std::vector<cli_opt_spec>& rSpecs) {
    cli_parsed out;
    bool positional_only = false;
    for (std::size_t i = 0; i < rArgs.size(); ++i) {
        const std::string& tok = rArgs[i];
        if (positional_only || tok.empty() || tok[0] != '-' || tok == "-") {
            out.positionals.push_back(tok);
            continue;
        }
        if (tok == "--") {
            positional_only = true;
            continue;
        }
        // Split "--name=value".
        std::string name = tok;
        std::string inline_value;
        bool has_inline = false;
        auto eq = tok.find('=');
        if (tok.rfind("--", 0) == 0 && eq != std::string::npos) {
            name = tok.substr(0, eq);
            inline_value = tok.substr(eq + 1);
            has_inline = true;
        }
        const cli_opt_spec* match = nullptr;
        for (const auto& s : rSpecs) {
            if (name == "--" + s.canonical) {
                match = &s;
                break;
            }
            for (const auto& a : s.aliases) {
                if (name == a) {
                    match = &s;
                    break;
                }
            }
            if (match)
                break;
        }
        if (!match)
            throw std::runtime_error("unknown option '" + name + "'");
        if (match->takes_value) {
            std::string value;
            if (has_inline) {
                value = inline_value;
            } else if (i + 1 < rArgs.size()) {
                value = rArgs[++i];
            } else {
                throw std::runtime_error("option '" + name + "' requires a value");
            }
            out.values[match->canonical] = value;
        } else {
            out.flags.insert(match->canonical);
        }
    }
    return out;
}

std::string opt_value(const cli_parsed& rP, const std::string& rName,
                      const std::string& rDefault = "") {
    auto it = rP.values.find(rName);
    return it == rP.values.end() ? rDefault : it->second;
}
bool has_flag(const cli_parsed& rP, const std::string& rName) {
    return rP.flags.count(rName) != 0;
}

// --------------------------------------------------------------------------
// Format dispatch (mirrors bindings_c / bindings_js)
// --------------------------------------------------------------------------

std::string compiled_out_hint(const std::string& rFormat) {
    const char* dep = meshioplusplus::registry_compiled_out(rFormat);
    if (dep)
        return " (format '" + rFormat + "' is not available in this build; requires " + dep + ")";
    return "";
}

Mesh read_mesh_cli(const std::string& rPath, const std::string& rFormat) {
    std::string fmt;
    try {
        fmt = meshioplusplus::resolve_format(rPath, rFormat);
    } catch (const ReadError&) {
        fmt = meshioplusplus::sniff_format(rPath);
        if (fmt.empty())
            throw;
    }
    const auto& readers = meshioplusplus::registry_readers();
    auto it = readers.find(fmt);
    if (it == readers.end())
        throw ReadError("no reader for format '" + fmt + "'" + compiled_out_hint(fmt));
    return it->second(rPath);
}

void write_mesh_cli(const std::string& rPath, const Mesh& rMesh, const std::string& rFormat) {
    std::string fmt = meshioplusplus::resolve_format(rPath, rFormat);
    const auto& writers = meshioplusplus::registry_writers();
    auto it = writers.find(fmt);
    if (it == writers.end())
        throw WriteError("no writer for format '" + fmt + "'" + compiled_out_hint(fmt));
    it->second(rPath, rMesh);
}

/// Write `rMesh` in the ASCII (`binary=false`) or binary (`binary=true`) variant
/// of `rFormat`, bypassing the registry's baked-in default. Returns false when
/// the format has no such variant (the caller reports the error).
bool write_binary_variant(const std::string& rPath, const Mesh& rMesh, const std::string& rFormat,
                          bool binary, const std::string& rFloatFmt) {
    const std::string& ff = rFloatFmt.empty() ? std::string(".16e") : rFloatFmt;
    if (rFormat == "ansys") {
        meshioplusplus::write_ansys(rPath, rMesh, binary);
    } else if (rFormat == "flac3d") {
        meshioplusplus::write_flac3d(rPath, rMesh, ff, binary);
    } else if (rFormat == "gmsh") {
        meshioplusplus::write_gmsh41(rPath, rMesh, binary);
    } else if (rFormat == "ply") {
        meshioplusplus::write_ply(rPath, rMesh, binary, /*skin=*/true);
    } else if (rFormat == "stl") {
        meshioplusplus::write_stl(rPath, rMesh, binary, /*skin=*/true);
    } else if (rFormat == "vtk") {
        meshioplusplus::write_vtk(rPath, rMesh, binary, /*v51=*/true);
    } else if (rFormat == "vtu") {
        meshioplusplus::write_vtu(rPath, rMesh, binary, /*zlib=*/binary);
    } else if (rFormat == "xdmf") {
        meshioplusplus::write_xdmf(rPath, rMesh, binary ? "HDF" : "XML");
    } else {
        return false;
    }
    return true;
}

// --------------------------------------------------------------------------
// Helpers
// --------------------------------------------------------------------------

std::string file_size_mb(const std::string& rPath) {
    std::error_code ec;
    auto n = std::filesystem::file_size(rPath, ec);
    double mb = ec ? 0.0 : static_cast<double>(n) / (1024.0 * 1024.0);
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.2f MB", mb);
    return buf;
}

/// Invoke `f(v)` for every element of an integer connectivity array as int64.
template <class F>
void each_conn_value(const NDArray& rArr, F&& f) {
    std::size_t n = rArr.Size();
    switch (rArr.Dtype()) {
        case DType::Int64:
            for (std::size_t i = 0; i < n; ++i)
                f(static_cast<std::int64_t>(rArr.As<std::int64_t>()[i]));
            break;
        case DType::Int32:
            for (std::size_t i = 0; i < n; ++i)
                f(static_cast<std::int64_t>(rArr.As<std::int32_t>()[i]));
            break;
        case DType::UInt64:
            for (std::size_t i = 0; i < n; ++i)
                f(static_cast<std::int64_t>(rArr.As<std::uint64_t>()[i]));
            break;
        case DType::UInt32:
            for (std::size_t i = 0; i < n; ++i)
                f(static_cast<std::int64_t>(rArr.As<std::uint32_t>()[i]));
            break;
        default:
            break;  // non-integer connectivity: skip the diagnostic
    }
}

/// Iterate every node index referenced by a cell block (rectangular or ragged).
template <class CellView, class F>
void each_block_index(const CellView& rCv, F&& f) {
    if (!rCv.IsRagged()) {
        each_conn_value(rCv.Conn(), f);
        return;
    }
    for (std::size_t c = 0; c < rCv.NumCells(); ++c) {
        if (rCv.IsPolyhedron()) {
            for (std::size_t face = 0; face < rCv.NumFaces(c); ++face) {
                auto [ptr, len] = rCv.Face(c, face);
                for (std::size_t k = 0; k < len; ++k)
                    f(ptr[k]);
            }
        } else {
            const std::int64_t* ptr = rCv.Row(c);
            std::size_t len = rCv.RowSize(c);
            for (std::size_t k = 0; k < len; ++k)
                f(ptr[k]);
        }
    }
}

std::string join(const std::vector<std::string>& rNames) {
    std::string out;
    for (std::size_t i = 0; i < rNames.size(); ++i) {
        if (i)
            out += ", ";
        out += rNames[i];
    }
    return out;
}

// --------------------------------------------------------------------------
// Verbs
// --------------------------------------------------------------------------

const char* kVersionText = "meshio++ " MESHIOPLUSPLUS_CLI_VERSION
                           " [C++ native]\n"
                           "Copyright (c) 2015-2021 Nico Schlomer et al. (as meshio)\n"
                           "Copyright (c) 2025 Vicente Mataix Ferrandiz\n"
                           "Copyright (c) 2026 the meshio++ contributors";

void print_usage(std::ostream& os) {
    os << "usage: meshioplusplus <command> [options]\n\n"
          "Mesh input/output tools (native C++ CLI, no Python required).\n\n"
          "commands:\n"
          "  convert (c)             Convert between mesh formats\n"
          "  info (i)                Print mesh info\n"
          "  ascii (a)               Rewrite a file in its ASCII variant (in place)\n"
          "  binary (b)              Rewrite a file in its binary variant (in place)\n"
          "  compress                Compress a mesh file (in place)\n"
          "  decompress              Decompress a mesh file (in place)\n"
          "  quality (q)             Print mesh quality metrics\n"
          "  extract-surface (surface)  Extract the boundary surface/edges\n"
          "  reorder                 Renumber nodes/elements (RCM / Morton / Hilbert)\n"
          "  diff                    Compare two meshes (nonzero exit if different)\n"
          "  merge                   Merge two or more meshes into one\n"
          "  transform               Affine transform (translate/scale/rotate/matrix/units)\n"
          "  clean                   Weld / prune / de-dup a mesh\n"
          "  crop                    Subset by bounding box or half-space\n"
          "  split                   Partition into multiple files (type/region/component)\n"
          "  stats                   Print geometric statistics (bbox/area/volume)\n\n"
          "  -v, --version           Display version information\n"
          "  -h, --help              Show this message\n\n"
          "notes: point/cell sets and 'convert -s/-d' are unavailable in the native\n"
          "       CLI (they live only in the Python Mesh); use the Python CLI for those.\n";
}

int cmd_convert(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"float-format", {"-f"}, true},
                                  {"ascii", {"-a"}, false},
                                  {"sets-to-int-data", {"-s"}, false},
                                  {"int-data-to-sets", {"-d"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("convert requires exactly INFILE and OUTFILE");
    if (has_flag(p, "sets-to-int-data") || has_flag(p, "int-data-to-sets"))
        throw std::runtime_error(
            "the -s/--sets-to-int-data and -d/--int-data-to-sets options are not "
            "supported by the native CLI (sets live only in the Python Mesh); use the "
            "Python CLI for these");

    const std::string& infile = p.positionals[0];
    const std::string& outfile = p.positionals[1];
    std::string in_fmt = opt_value(p, "input-format");
    std::string out_fmt = opt_value(p, "output-format");
    std::string float_fmt = opt_value(p, "float-format");
    bool ascii = has_flag(p, "ascii");

    Mesh mesh = read_mesh_cli(infile, in_fmt);

    if (ascii) {
        std::string fmt = meshioplusplus::resolve_format(outfile, out_fmt);
        if (!write_binary_variant(outfile, mesh, fmt, /*binary=*/false, float_fmt))
            throw std::runtime_error("format '" + fmt + "' has no ASCII variant");
    } else {
        if (!float_fmt.empty())
            std::cerr << "note: --float-format only affects ASCII output (--ascii)\n";
        write_mesh_cli(outfile, mesh, out_fmt);
    }
    return 0;
}

int cmd_info(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {{"input-format", {"-i"}, true}});
    if (p.positionals.size() != 1)
        throw std::runtime_error("info requires exactly INFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    // Mirror Mesh.__repr__ (point/cell sets are unavailable in the C++ core).
    std::cout << "<meshio++ mesh object>\n";
    std::cout << "  Number of points: " << mesh.NumPoints() << "\n";
    if (mesh.NumCellBlocks() > 0) {
        std::cout << "  Number of cells:\n";
        for (const auto cb : mesh.CellRange()) {
            std::string type = cb.Type();
            if ((type == "polygon" || type == "polyhedron") && !cb.IsRagged())
                type += "(" + std::to_string(cb.NodesPerCell()) + ")";
            std::cout << "    " << type << ": " << cb.NumCells() << "\n";
        }
    } else {
        std::cout << "  No cells.\n";
    }
    if (mesh.NumPointData() > 0)
        std::cout << "  Point data: " << join(mesh.PointDataNames()) << "\n";
    if (mesh.NumCellData() > 0)
        std::cout << "  Cell data: " << join(mesh.CellDataNames()) << "\n";
    if (mesh.NumFieldData() > 0)
        std::cout << "  Field data: " << join(mesh.FieldDataNames()) << "\n";

    // Consistency checks (Python _info.py).
    const std::int64_t num_points = static_cast<std::int64_t>(mesh.NumPoints());
    bool is_consistent = true;
    for (const auto cb : mesh.CellRange()) {
        bool bad = false;
        each_block_index(cb, [&](std::int64_t v) {
            if (v > num_points)
                bad = true;
        });
        if (bad) {
            std::cerr << "Warning: Inconsistent mesh. Cells refer to nonexistent points.\n";
            is_consistent = false;
            break;
        }
    }
    if (is_consistent && num_points > 0) {
        std::vector<char> used(static_cast<std::size_t>(num_points), 0);
        for (const auto cb : mesh.CellRange())
            each_block_index(cb, [&](std::int64_t v) {
                if (v >= 0 && v < num_points)
                    used[static_cast<std::size_t>(v)] = 1;
            });
        bool any_unused = false;
        for (char u : used)
            if (!u) {
                any_unused = true;
                break;
            }
        if (any_unused)
            std::cerr << "Warning: Some points are not part of any cell.\n";
    }
    return 0;
}

int cmd_ascii_binary(const std::vector<std::string>& rArgs, bool binary) {
    auto p = cli_parse(rArgs, {{"input-format", {"-i"}, true}});
    if (p.positionals.size() != 1)
        throw std::runtime_error(std::string(binary ? "binary" : "ascii") +
                                 " requires exactly INFILE");
    const std::string& infile = p.positionals[0];
    std::cout << "File size before: " << file_size_mb(infile) << "\n";
    Mesh mesh = read_mesh_cli(infile, opt_value(p, "input-format"));
    std::string fmt = meshioplusplus::resolve_format(infile, opt_value(p, "input-format"));
    if (!write_binary_variant(infile, mesh, fmt, binary, ""))
        throw std::runtime_error("don't know how to write '" + fmt + "' as " +
                                 (binary ? "binary" : "ASCII"));
    std::cout << "File size after: " << file_size_mb(infile) << "\n";
    return 0;
}

int cmd_compress(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {{"input-format", {"-i"}, true}, {"max", {"-max"}, false}});
    if (p.positionals.size() != 1)
        throw std::runtime_error("compress requires exactly INFILE");
    const std::string& infile = p.positionals[0];
    bool max = has_flag(p, "max");
    int gzip = max ? 9 : 4;
    std::cout << "File size before: " << file_size_mb(infile) << "\n";
    Mesh mesh = read_mesh_cli(infile, opt_value(p, "input-format"));
    std::string fmt = meshioplusplus::resolve_format(infile, opt_value(p, "input-format"));

    if (fmt == "ansys" || fmt == "gmsh" || fmt == "ply" || fmt == "stl" || fmt == "vtk") {
        write_binary_variant(infile, mesh, fmt, /*binary=*/true, "");
    } else if (fmt == "vtu") {
        if (max)
            std::cerr << "note: the native CLI has no lzma variant; using zlib\n";
        meshioplusplus::write_vtu(infile, mesh, /*binary=*/true, /*zlib=*/true);
    } else if (fmt == "xdmf") {
        meshioplusplus::write_xdmf(infile, mesh, "HDF", gzip);
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    } else if (fmt == "cgns") {
        meshioplusplus::write_cgns(infile, mesh, gzip);
    } else if (fmt == "h5m") {
        meshioplusplus::write_h5m(infile, mesh, /*add_global_ids=*/true, gzip);
#endif
    } else {
        throw std::runtime_error("don't know how to compress '" + fmt + "'" +
                                 compiled_out_hint(fmt));
    }
    std::cout << "File size after: " << file_size_mb(infile) << "\n";
    return 0;
}

int cmd_decompress(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {{"input-format", {"-i"}, true}});
    if (p.positionals.size() != 1)
        throw std::runtime_error("decompress requires exactly INFILE");
    const std::string& infile = p.positionals[0];
    std::cout << "File size before: " << file_size_mb(infile) << "\n";
    Mesh mesh = read_mesh_cli(infile, opt_value(p, "input-format"));
    std::string fmt = meshioplusplus::resolve_format(infile, opt_value(p, "input-format"));

    if (fmt == "vtu") {
        meshioplusplus::write_vtu(infile, mesh, /*binary=*/true, /*zlib=*/false);
    } else if (fmt == "xdmf") {
        meshioplusplus::write_xdmf(infile, mesh, "HDF", /*gzip_level=*/-1);
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    } else if (fmt == "cgns") {
        meshioplusplus::write_cgns(infile, mesh, /*gzip_level=*/-1);
    } else if (fmt == "h5m") {
        meshioplusplus::write_h5m(infile, mesh, /*add_global_ids=*/true, /*gzip_level=*/-1);
#endif
    } else {
        throw std::runtime_error("don't know how to decompress '" + fmt + "'" +
                                 compiled_out_hint(fmt));
    }
    std::cout << "File size after: " << file_size_mb(infile) << "\n";
    return 0;
}

std::string quality_value(double v) {
    if (v != v)  // NaN
        return "  n/a ";
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%8.4f", v);
    return buf;
}
std::string right(const std::string& rS, int w) {
    if (static_cast<int>(rS.size()) >= w)
        return rS;
    return std::string(w - rS.size(), ' ') + rS;
}
std::string left(const std::string& rS, int w) {
    if (static_cast<int>(rS.size()) >= w)
        return rS;
    return rS + std::string(w - rS.size(), ' ');
}

int cmd_quality(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {{"input-format", {"-i"}, true}, {"output", {"-o"}, true}});
    if (p.positionals.size() != 1)
        throw std::runtime_error("quality requires exactly INFILE");
    const std::string& infile = p.positionals[0];
    Mesh mesh = read_mesh_cli(infile, opt_value(p, "input-format"));
    auto report = meshioplusplus::compute_quality(mesh);

    std::cout << "Mesh quality report for " << infile << "\n";
    std::cout << "  cells: " << report.mNumCells << "   inverted: " << report.mNumInverted
              << "   degenerate: " << report.mNumDegenerate << "\n";
    std::cout << "  " << left("metric", 24) << right("min", 10) << right("mean", 10)
              << right("max", 10) << right("count", 10) << "\n";
    for (const auto& [name, s] : report.mMetrics) {
        if (s.mCount == 0)
            continue;
        std::cout << "  " << left(name, 24) << right(quality_value(s.mMin), 10)
                  << right(quality_value(s.mMean), 10) << right(quality_value(s.mMax), 10)
                  << right(std::to_string(s.mCount), 10) << "\n";
    }

    std::string output = opt_value(p, "output");
    if (!output.empty())
        write_mesh_cli(output, meshioplusplus::attach_quality(mesh), "");
    return 0;
}

int cmd_extract_surface(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"parent-ids", {"-p"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("extract-surface requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    Mesh surface =
        meshioplusplus::extract_surface(mesh, /*recordParentIds=*/has_flag(p, "parent-ids"));
    write_mesh_cli(p.positionals[1], surface, opt_value(p, "output-format"));
    return 0;
}

std::vector<double> parse_doubles(const std::string& rText) {
    std::vector<double> out;
    std::size_t start = 0;
    while (start <= rText.size()) {
        std::size_t comma = rText.find(',', start);
        std::string tok =
            rText.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        if (!tok.empty())
            out.push_back(std::stod(tok));
        if (comma == std::string::npos)
            break;
        start = comma + 1;
    }
    return out;
}

int cmd_transform(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"translate", {}, true},
                                  {"scale", {}, true},
                                  {"rotate", {}, true},
                                  {"matrix", {}, true},
                                  {"scale-units", {}, true},
                                  {"rotate-data", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("transform requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::AffineTransform xf;
    int given = 0;
    if (p.values.count("translate")) {
        auto v = parse_doubles(opt_value(p, "translate"));
        if (v.size() != 3)
            throw std::runtime_error("transform: --translate expects 'x,y,z'");
        xf = meshioplusplus::transform_translation(v[0], v[1], v[2]);
        ++given;
    }
    if (p.values.count("scale")) {
        auto v = parse_doubles(opt_value(p, "scale"));
        if (v.size() == 1)
            xf = meshioplusplus::transform_scale(v[0], v[0], v[0]);
        else if (v.size() == 3)
            xf = meshioplusplus::transform_scale(v[0], v[1], v[2]);
        else
            throw std::runtime_error("transform: --scale expects 's' or 'sx,sy,sz'");
        ++given;
    }
    if (p.values.count("rotate")) {
        auto v = parse_doubles(opt_value(p, "rotate"));
        const std::string& raw = opt_value(p, "rotate");
        double ax = 0, ay = 0, az = 0, deg = 0;
        if (v.size() == 4) {
            ax = v[0];
            ay = v[1];
            az = v[2];
            deg = v[3];
        } else if (v.size() == 2 && !raw.empty() &&
                   (raw[0] == 'x' || raw[0] == 'y' || raw[0] == 'z')) {
            ax = raw[0] == 'x' ? 1 : 0;
            ay = raw[0] == 'y' ? 1 : 0;
            az = raw[0] == 'z' ? 1 : 0;
            deg = v[0];  // parse_doubles skips the non-numeric axis token
        } else {
            throw std::runtime_error(
                "transform: --rotate expects 'axis,deg' (x|y|z) or "
                "'nx,ny,nz,deg'");
        }
        xf = meshioplusplus::transform_rotation(ax, ay, az, deg * 3.14159265358979323846 / 180.0);
        ++given;
    }
    if (p.values.count("matrix")) {
        auto v = parse_doubles(opt_value(p, "matrix"));
        if (v.size() != 16)
            throw std::runtime_error("transform: --matrix expects 16 values (row-major 4x4)");
        xf = meshioplusplus::transform_from_matrix(v.data());
        ++given;
    }
    if (p.values.count("scale-units")) {
        xf = meshioplusplus::transform_units(std::stod(opt_value(p, "scale-units")));
        ++given;
    }
    if (given != 1)
        throw std::runtime_error(
            "transform: give exactly one of "
            "--translate/--scale/--rotate/--matrix/--scale-units");

    Mesh out = meshioplusplus::transform(mesh, xf, has_flag(p, "rotate-data"));
    write_mesh_cli(p.positionals[1], out, opt_value(p, "output-format"));
    return 0;
}

int cmd_clean(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"atol", {}, true},
                                  {"weld", {}, false},
                                  {"remove-orphans", {}, false},
                                  {"drop-degenerate", {}, false},
                                  {"drop-duplicates", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("clean requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::CleanOptions opts;
    const bool any_flag = has_flag(p, "weld") || has_flag(p, "remove-orphans") ||
                          has_flag(p, "drop-degenerate") || has_flag(p, "drop-duplicates");
    if (any_flag) {
        opts.weld = has_flag(p, "weld");
        opts.remove_orphans = has_flag(p, "remove-orphans");
        opts.drop_degenerate = has_flag(p, "drop-degenerate");
        opts.drop_duplicate_cells = has_flag(p, "drop-duplicates");
    } else {
        opts.weld = false;
        opts.remove_orphans = true;
        opts.drop_degenerate = true;
        opts.drop_duplicate_cells = true;
    }
    if (p.values.count("atol"))
        opts.atol = std::stod(opt_value(p, "atol"));

    auto r = meshioplusplus::clean(mesh, opts);
    if (!has_flag(p, "quiet")) {
        std::cout << "cleaned mesh\n";
        std::cout << "  points welded:             " << r.mPointsWelded << "\n";
        std::cout << "  points removed (orphan):   " << r.mPointsRemovedOrphan << "\n";
        std::cout << "  cells dropped (degenerate): " << r.mCellsDroppedDegenerate << "\n";
        std::cout << "  cells dropped (duplicate):  " << r.mCellsDroppedDuplicate << "\n";
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_crop(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"bbox", {}, true},
                                  {"plane", {}, true},
                                  {"mode", {}, true},
                                  {"record-ids", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("crop requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    const bool has_bbox = p.values.count("bbox") != 0;
    const bool has_plane = p.values.count("plane") != 0;
    if (has_bbox == has_plane)
        throw std::runtime_error("crop: give exactly one of --bbox or --plane");
    std::string mode_s = opt_value(p, "mode", "all");
    if (mode_s != "all" && mode_s != "any")
        throw std::runtime_error("crop: --mode must be 'all' or 'any'");
    meshioplusplus::CropMode mode =
        mode_s == "any" ? meshioplusplus::CropMode::Any : meshioplusplus::CropMode::All;
    const bool record_ids = has_flag(p, "record-ids");

    meshioplusplus::CropResult r;
    if (has_bbox) {
        auto v = parse_doubles(opt_value(p, "bbox"));
        if (v.size() != 6)
            throw std::runtime_error("crop: --bbox expects 'xmin,ymin,zmin,xmax,ymax,zmax'");
        double lo[3] = {v[0], v[1], v[2]}, hi[3] = {v[3], v[4], v[5]};
        r = meshioplusplus::crop_bbox(mesh, lo, hi, mode, record_ids);
    } else {
        auto v = parse_doubles(opt_value(p, "plane"));
        if (v.size() != 6)
            throw std::runtime_error("crop: --plane expects 'px,py,pz,nx,ny,nz'");
        double point[3] = {v[0], v[1], v[2]}, normal[3] = {v[3], v[4], v[5]};
        r = meshioplusplus::crop_halfspace(mesh, point, normal, mode, record_ids);
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

std::string replace_key(const std::string& rPattern, const std::string& rKey) {
    const std::string token = "{key}";
    std::string out = rPattern;
    std::size_t pos = out.find(token);
    while (pos != std::string::npos) {
        out.replace(pos, token.size(), rKey);
        pos = out.find(token, pos + rKey.size());
    }
    return out;
}

int cmd_split(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"by", {}, true},
                                  {"tag", {}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("split requires exactly INFILE and OUTPATTERN (with {key})");
    if (p.positionals[1].find("{key}") == std::string::npos)
        throw std::runtime_error("split: output pattern must contain '{key}' (e.g. out_{key}.vtu)");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    std::string by = opt_value(p, "by", "type");
    auto crit = meshioplusplus::split_by_from_name(by);
    auto result = meshioplusplus::split(mesh, crit, opt_value(p, "tag"));

    std::string out_fmt = opt_value(p, "output-format");
    if (!has_flag(p, "quiet"))
        std::cout << "split into " << result.mPieces.size() << " piece(s) by " << by << "\n";
    for (auto& piece : result.mPieces) {
        std::string path = replace_key(p.positionals[1], piece.mKey);
        std::int64_t ncells = 0;
        for (const auto cb : piece.mMesh.CellRange())
            ncells += static_cast<std::int64_t>(cb.NumCells());
        if (!has_flag(p, "quiet"))
            std::cout << "  " << piece.mKey << ": " << piece.mMesh.NumPoints() << " points, "
                      << ncells << " cells -> " << path << "\n";
        write_mesh_cli(path, piece.mMesh, out_fmt);
    }
    return 0;
}

std::string stats_g6(double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.6g", v);
    return buf;
}

int cmd_stats(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"json", {}, false},
                              });
    if (p.positionals.size() != 1)
        throw std::runtime_error("stats requires exactly INFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    meshioplusplus::StatsReport s = meshioplusplus::compute_stats(mesh);

    auto vec3 = [&](const double* v) {
        return "[" + stats_g6(v[0]) + ", " + stats_g6(v[1]) + ", " + stats_g6(v[2]) + "]";
    };

    if (has_flag(p, "json")) {
        std::cout << "{\n";
        std::cout << "  \"num_points\": " << s.mNumPoints << ",\n";
        std::cout << "  \"num_cells\": " << s.mNumCells << ",\n";
        std::cout << "  \"bbox_min\": " << vec3(s.mBBoxMin) << ",\n";
        std::cout << "  \"bbox_max\": " << vec3(s.mBBoxMax) << ",\n";
        std::cout << "  \"extent\": " << vec3(s.mExtent) << ",\n";
        std::cout << "  \"centroid\": " << vec3(s.mCentroid) << ",\n";
        std::cout << "  \"cell_type_counts\": {";
        for (std::size_t i = 0; i < s.mCellTypeCounts.size(); ++i)
            std::cout << (i ? ", " : "") << "\"" << s.mCellTypeCounts[i].first
                      << "\": " << s.mCellTypeCounts[i].second;
        std::cout << "},\n";
        std::cout << "  \"total_area\": " << stats_g6(s.mTotalArea) << ",\n";
        std::cout << "  \"signed_volume\": " << stats_g6(s.mSignedVolume) << ",\n";
        std::cout << "  \"unsigned_volume\": " << stats_g6(s.mUnsignedVolume) << ",\n";
        std::cout << "  \"num_inverted\": " << s.mNumInverted << "\n";
        std::cout << "}\n";
        return 0;
    }

    std::cout << "<meshio++ geometric stats>\n";
    std::cout << "  points:          " << s.mNumPoints << "\n";
    std::cout << "  cells:           " << s.mNumCells << "\n";
    std::cout << "  bbox min:        " << vec3(s.mBBoxMin) << "\n";
    std::cout << "  bbox max:        " << vec3(s.mBBoxMax) << "\n";
    std::cout << "  extent:          " << vec3(s.mExtent) << "\n";
    std::cout << "  centroid:        " << vec3(s.mCentroid) << "\n";
    std::cout << "  cell types:\n";
    for (const auto& kv : s.mCellTypeCounts)
        std::cout << "    " << kv.first << ": " << kv.second << "\n";
    std::cout << "  total area:      " << stats_g6(s.mTotalArea) << "\n";
    std::cout << "  signed volume:   " << stats_g6(s.mSignedVolume) << "\n";
    std::cout << "  unsigned volume: " << stats_g6(s.mUnsignedVolume) << "\n";
    std::cout << "  inverted cells:  " << s.mNumInverted << "\n";
    return 0;
}

int cmd_reorder(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"method", {"-m"}, true},
                                  {"report", {"-r"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("reorder requires exactly INFILE and OUTFILE");
    std::string method_name = opt_value(p, "method", "rcm");
    auto method = meshioplusplus::reorder_method_from_name(method_name);
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    std::int64_t before = 0;
    bool report = has_flag(p, "report");
    if (report)
        before = meshioplusplus::compute_bandwidth(mesh);
    auto result = meshioplusplus::reorder(mesh, method);
    if (report) {
        std::int64_t after = meshioplusplus::compute_bandwidth(result.mMesh);
        std::cout << "bandwidth (" << method_name << "): " << before << " -> " << after << "\n";
    }
    write_mesh_cli(p.positionals[1], result.mMesh, opt_value(p, "output-format"));
    return 0;
}

std::string sci3(double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return buf;
}

void print_data_diff(const std::string& rSection, const meshioplusplus::DataDiff& rDd) {
    if (!rDd.mOnlyInA.empty())
        std::cout << rSection << ": only in A: " << join(rDd.mOnlyInA) << "\n";
    if (!rDd.mOnlyInB.empty())
        std::cout << rSection << ": only in B: " << join(rDd.mOnlyInB) << "\n";
    for (const auto& ad : rDd.mShared) {
        if (ad.mShapeMismatch)
            std::cout << rSection << "['" << ad.mName << "']: shape mismatch\n";
        else if (!ad.mExact)
            std::cout << rSection << "['" << ad.mName << "']: max_abs=" << sci3(ad.mMaxAbsError)
                      << " max_rel=" << sci3(ad.mMaxRelError) << " exceeding=" << ad.mNumExceeding
                      << "\n";
    }
}

int cmd_diff(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format-a", {}, true},
                                  {"input-format-b", {}, true},
                                  {"atol", {}, true},
                                  {"rtol", {}, true},
                                  {"unordered", {}, false},
                                  {"exact", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("diff requires exactly two mesh files A and B");
    Mesh a = read_mesh_cli(p.positionals[0], opt_value(p, "input-format-a"));
    Mesh b = read_mesh_cli(p.positionals[1], opt_value(p, "input-format-b"));

    meshioplusplus::DiffOptions opts;
    if (p.values.count("atol"))
        opts.atol = std::stod(opt_value(p, "atol"));
    if (p.values.count("rtol"))
        opts.rtol = std::stod(opt_value(p, "rtol"));
    opts.unordered = has_flag(p, "unordered");
    auto report = meshioplusplus::diff(a, b, opts);

    if (!has_flag(p, "quiet")) {
        std::cout << "verdict: " << meshioplusplus::diff_verdict_name(report.mVerdict) << "\n";
        for (const auto& msg : report.mMessages)
            std::cout << "  note: " << msg << "\n";
        if (report.mPointCountMismatch) {
            std::cout << "points: count/dim differ (" << report.mNumPointsA << " vs "
                      << report.mNumPointsB << ")\n";
        } else if (report.mPoints.mShapeMismatch) {
            std::cout << "points: shape mismatch\n";
        } else if (!report.mPoints.mExact) {
            std::cout << "points: max_abs=" << sci3(report.mPoints.mMaxAbsError)
                      << " max_rel=" << sci3(report.mPoints.mMaxRelError)
                      << " exceeding=" << report.mPoints.mNumExceeding
                      << " worst_index=" << report.mPoints.mWorstIndex << "\n";
        }
        if (report.mBlockCountMismatch)
            std::cout << "cells: block count differs (" << report.mNumBlocksA << " vs "
                      << report.mNumBlocksB << ")\n";
        for (const auto& bd : report.mBlocks) {
            if (bd.mTypeMismatch)
                std::cout << "cells[" << bd.mBlock << "]: type " << bd.mTypeA << " vs " << bd.mTypeB
                          << "\n";
            else if (bd.mCountMismatch)
                std::cout << "cells[" << bd.mBlock << "] (" << bd.mTypeA << "): count "
                          << bd.mCountA << " vs " << bd.mCountB << "\n";
            else if (bd.mConnMismatchCount > 0) {
                std::cout << "cells[" << bd.mBlock << "] (" << bd.mTypeA
                          << "): " << bd.mConnMismatchCount << " connectivity mismatch(es)";
                if (!bd.mFirstMismatches.empty()) {
                    std::cout << " [first: ";
                    for (std::size_t i = 0; i < bd.mFirstMismatches.size(); ++i)
                        std::cout << (i ? ", " : "") << bd.mFirstMismatches[i];
                    std::cout << "]";
                }
                std::cout << "\n";
            }
        }
        print_data_diff("point_data", report.mPointData);
        print_data_diff("cell_data", report.mCellData);
        print_data_diff("field_data", report.mFieldData);
    }

    bool exact = has_flag(p, "exact");
    bool equal = report.mVerdict == meshioplusplus::DiffVerdict::Identical ||
                 (report.mVerdict == meshioplusplus::DiffVerdict::EqualWithinTolerance && !exact);
    return equal ? 0 : 1;
}

std::int64_t total_cells(const Mesh& rMesh) {
    std::int64_t n = 0;
    for (const auto cb : rMesh.CellRange())
        n += static_cast<std::int64_t>(cb.NumCells());
    return n;
}

int cmd_merge(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"atol", {}, true},
                                  {"data-policy", {}, true},
                                  {"weld", {}, false},
                                  {"drop-duplicate-cells", {}, false},
                                  {"no-source-tag", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() < 2)
        throw std::runtime_error(
            "merge: need at least one input mesh and an output file "
            "(e.g. `meshioplusplus merge a.vtu b.vtu out.vtu`)");
    std::string in_fmt = opt_value(p, "input-format");
    std::string out_fmt = opt_value(p, "output-format");
    std::string out_path = p.positionals.back();

    std::vector<Mesh> meshes;
    meshes.reserve(p.positionals.size() - 1);
    for (std::size_t i = 0; i + 1 < p.positionals.size(); ++i)
        meshes.push_back(read_mesh_cli(p.positionals[i], in_fmt));

    std::int64_t total_in_points = 0, total_in_cells = 0;
    for (const auto& m : meshes) {
        total_in_points += static_cast<std::int64_t>(m.NumPoints());
        total_in_cells += total_cells(m);
    }

    meshioplusplus::MergeOptions opts;
    opts.weld = has_flag(p, "weld");
    if (p.values.count("atol"))
        opts.atol = std::stod(opt_value(p, "atol"));
    opts.source_tag = !has_flag(p, "no-source-tag");
    opts.drop_duplicate_cells = has_flag(p, "drop-duplicate-cells");
    std::string policy = opt_value(p, "data-policy", "intersection");
    if (policy == "fill")
        opts.data_policy = meshioplusplus::MergeDataPolicy::Fill;
    else if (policy == "intersection")
        opts.data_policy = meshioplusplus::MergeDataPolicy::Intersection;
    else
        throw std::runtime_error("merge: --data-policy must be 'intersection' or 'fill'");

    std::vector<const Mesh*> ptrs;
    ptrs.reserve(meshes.size());
    for (const auto& m : meshes)
        ptrs.push_back(&m);
    auto result = meshioplusplus::merge(ptrs, opts);

    std::int64_t out_points = static_cast<std::int64_t>(result.mMesh.NumPoints());
    std::int64_t out_cells = total_cells(result.mMesh);
    if (!has_flag(p, "quiet")) {
        std::cout << "merged " << meshes.size() << " meshes\n";
        std::cout << "  points in:  " << total_in_points << "\n";
        std::cout << "  cells in:   " << total_in_cells << "\n";
        if (opts.weld)
            std::cout << "  points welded: " << (total_in_points - out_points) << "\n";
        std::cout << "  points out: " << out_points << "\n";
        std::cout << "  cells out:  " << out_cells << "\n";
    }
    write_mesh_cli(out_path, result.mMesh, out_fmt);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> args(argv + (argc > 0 ? 1 : 0), argv + argc);

    if (args.empty()) {
        print_usage(std::cerr);
        return 2;
    }
    std::string cmd = args[0];
    std::vector<std::string> rest(args.begin() + 1, args.end());

    if (cmd == "-v" || cmd == "--version") {
        std::cout << kVersionText << "\n";
        return 0;
    }
    if (cmd == "-h" || cmd == "--help") {
        print_usage(std::cout);
        return 0;
    }

    try {
        if (cmd == "convert" || cmd == "c")
            return cmd_convert(rest);
        if (cmd == "info" || cmd == "i")
            return cmd_info(rest);
        if (cmd == "ascii" || cmd == "a")
            return cmd_ascii_binary(rest, /*binary=*/false);
        if (cmd == "binary" || cmd == "b")
            return cmd_ascii_binary(rest, /*binary=*/true);
        if (cmd == "compress")
            return cmd_compress(rest);
        if (cmd == "decompress")
            return cmd_decompress(rest);
        if (cmd == "quality" || cmd == "q")
            return cmd_quality(rest);
        if (cmd == "extract-surface" || cmd == "surface")
            return cmd_extract_surface(rest);
        if (cmd == "reorder")
            return cmd_reorder(rest);
        if (cmd == "diff")
            return cmd_diff(rest);
        if (cmd == "merge")
            return cmd_merge(rest);
        if (cmd == "transform")
            return cmd_transform(rest);
        if (cmd == "clean")
            return cmd_clean(rest);
        if (cmd == "crop")
            return cmd_crop(rest);
        if (cmd == "split")
            return cmd_split(rest);
        if (cmd == "stats")
            return cmd_stats(rest);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }

    std::cerr << "error: unknown command '" << cmd << "'\n\n";
    print_usage(std::cerr);
    return 2;
}
