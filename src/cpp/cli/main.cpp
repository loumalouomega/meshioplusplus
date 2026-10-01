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
 * the Python console script (`src/python/meshioplusplus/_cli/`) but links only the C++
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
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iterator>
#include <cstring>
#include <utility>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/formats/pvtp.hpp"
#include "meshioplusplus/formats/pvtu.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/write_options.hpp"

#include "polyscope_view.hpp"
#include "view_payload.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/pipeline.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/data_calc.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/data_condition.hpp"
#include "meshioplusplus/operations/tensor_invariants.hpp"
#include "meshioplusplus/operations/data_info.hpp"
#include "meshioplusplus/operations/data_integrate.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/operations/data_manage.hpp"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/decimate_volume.hpp"
#include "meshioplusplus/operations/remesh.hpp"
#include "meshioplusplus/operations/remesh_volume.hpp"
#include "meshioplusplus/operations/optimize_volume.hpp"
#include "meshioplusplus/operations/conservative_interpolate.hpp"
#include "meshioplusplus/operations/diff.hpp"
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/operations/voxelize.hpp"
#include "meshioplusplus/operations/error.hpp"
#include "meshioplusplus/operations/gradient.hpp"
#include "meshioplusplus/operations/hessian.hpp"
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/curvature.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/operations/feature_edges.hpp"
#include "meshioplusplus/operations/interfaces.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/operations/periodic.hpp"
#include "meshioplusplus/operations/quality_gate.hpp"
#include "meshioplusplus/operations/region_ops.hpp"
#include "meshioplusplus/operations/repair.hpp"
#include "meshioplusplus/operations/shrinkwrap.hpp"
#include "meshioplusplus/operations/sobolev_deform.hpp"
#include "meshioplusplus/operations/crop.hpp"
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/operations/subdivide.hpp"
#include "meshioplusplus/operations/undo_green.hpp"
// Per-format writers for the ASCII/binary/compress variants (the registry bakes
// one default each, so these are called directly).
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/pcd.hpp"
#include "meshioplusplus/formats/ply.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/gltf.hpp"
#include "meshioplusplus/formats/svg.hpp"
#include "meshioplusplus/formats/tikz.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/formats/xdmf.hpp"
#ifdef MESHIOPLUSPLUS_HAS_HDF5
#include "meshioplusplus/formats/cgns.hpp"
#include "meshioplusplus/formats/h5m.hpp"
#endif

#ifndef MESHIOPLUSPLUS_CLI_VERSION
#define MESHIOPLUSPLUS_CLI_VERSION "unknown"
#endif

// CLI-private
#include "json_out.hpp"
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
    std::unordered_map<std::string, std::string> values;  ///< canonical -> value (last wins)
    std::unordered_set<std::string> flags;                ///< canonical present
    /// Every occurrence of each value-option, in order. Always populated
    /// alongside `values`, which keeps its last-wins semantics so no existing
    /// verb changes behaviour; the `data` verbs read this instead because they
    /// accept repeated `--point`/`--cell`/`--field`.
    std::unordered_map<std::string, std::vector<std::string>> multi;
    /// Every value-option occurrence as (canonical name, value), in
    /// command-line order across options -- what `regions`' edit flags need.
    std::vector<std::pair<std::string, std::string>> ordered;
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
            out.multi[match->canonical].push_back(value);
            out.ordered.emplace_back(match->canonical, value);
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
/// Every occurrence of a repeatable value-option, in order (empty if absent).
const std::vector<std::string>& opt_values(const cli_parsed& rP, const std::string& rName) {
    static const std::vector<std::string> empty;
    auto it = rP.multi.find(rName);
    return it == rP.multi.end() ? empty : it->second;
}
/// Whether a value-option was supplied at all (distinct from "supplied empty").
bool has_opt(const cli_parsed& rP, const std::string& rName) {
    return rP.values.count(rName) != 0;
}

// --------------------------------------------------------------------------
// Format dispatch (mirrors bindings/c / bindings/wasm)
// --------------------------------------------------------------------------

std::string compiled_out_hint(const std::string& rFormat) {
    const char* dep = meshioplusplus::registry_compiled_out(rFormat);
    if (dep)
        return " (format '" + rFormat + "' is not available in this build; requires " + dep + ")";
    return "";
}

Mesh read_mesh_cli(const std::string& rPath, const std::string& rFormat,
                   const meshioplusplus::ReadOptions& rOpts = {}) {
    std::string fmt;
    try {
        fmt = meshioplusplus::resolve_format(rPath, rFormat);
    } catch (const ReadError&) {
        fmt = meshioplusplus::sniff_format(rPath);
        if (fmt.empty())
            throw;
    }
    if (!meshioplusplus::registry_readers().count(fmt) &&
        !meshioplusplus::registry_reader_supports_options(fmt))
        throw ReadError("no reader for format '" + fmt + "'" + compiled_out_hint(fmt));
    // registry_read honours the options where the format supports them and
    // falls back to a full read where it does not.
    return meshioplusplus::registry_read(rPath, fmt, rOpts);
}

void write_mesh_cli(const std::string& rPath, const Mesh& rMesh, const std::string& rFormat) {
    std::string fmt = meshioplusplus::resolve_write_format(rPath, rFormat);
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
    meshioplusplus::WriteOptions opts;
    opts.mEncoding =
        binary ? meshioplusplus::WriteEncoding::Binary : meshioplusplus::WriteEncoding::Ascii;
    // Only pass the float format to a writer that takes one: `--float-format`
    // on a format without one has always been ignored rather than rejected,
    // and registry_write_supports() would now turn that into an error.
    std::string why;
    meshioplusplus::WriteOptions probe;
    probe.mFloatFormat = ".16e";
    if (!rFloatFmt.empty() && meshioplusplus::registry_write_supports(rFormat, probe, why))
        opts.mFloatFormat = rFloatFmt;

    if (!meshioplusplus::registry_write_supports(rFormat, opts, why))
        return false;
    meshioplusplus::registry_write_ex(rPath, rMesh, rFormat, opts);
    return true;
}

/// Rewrite a PCD file in place with the `data` encoding, at the mesh's own coordinate precision:
/// the in-place verbs (`ascii`, `binary`, `compress`, `decompress`) must not narrow a float64
/// cloud, unlike `convert`, whose default is PCL's float32.
void write_pcd_in_place(const std::string& rPath, const Mesh& rMesh, meshioplusplus::PcdData data) {
    meshioplusplus::write_pcd(rPath, rMesh, data,
                              rMesh.Points().Dtype() == meshioplusplus::DType::Float64);
}

/// The formats `--color-by` applies to. Everything else errors rather than
/// silently ignoring the flags.
bool cli_is_colorable_format(const std::string& rFormat) {
    return rFormat == "svg" || rFormat == "tikz" || rFormat == "gltf";
}

/// Write `rMesh` as a data-coloured SVG/TikZ figure, bypassing the registry
/// (whose `(path, mesh)` writer lambdas cannot carry parameters) exactly as
/// `write_binary_variant` above does for the ASCII/binary variants. Every
/// pre-camera and camera argument keeps the writer's own default.
/// Returns false when the format has no coloured variant.
bool write_colored_variant(const std::string& rPath, const Mesh& rMesh, const std::string& rFormat,
                           const std::string& rColorBy, const std::optional<int>& rComponent,
                           const std::string& rCmap, const std::optional<double>& rVMin,
                           const std::optional<double>& rVMax, const std::string& rNanColor,
                           bool colorbar) {
    if (rFormat == "svg") {
        meshioplusplus::write_svg(rPath, rMesh, ".3f", std::nullopt, 100.0, "#c8c5bd", "#000080",
                                  45.0, 35.264389682754654, 0.0, rColorBy, rComponent, rCmap, rVMin,
                                  rVMax, rNanColor.empty() ? "#808080" : rNanColor, colorbar);
    } else if (rFormat == "tikz") {
        meshioplusplus::write_tikz(rPath, rMesh, ".6f", true, std::nullopt, "gray!30", "black",
                                   std::nullopt, 45.0, 35.264389682754654, 0.0, rColorBy,
                                   rComponent, rCmap, rVMin, rVMax,
                                   rNanColor.empty() ? "gray" : rNanColor, colorbar);
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
    meshioplusplus::detail::snprintf_c(buf, sizeof(buf), "%.2f MB", mb);
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
          "                            --points-only / --arrays a,b narrow what is read\n"
          "                            --time-step=N picks a step of a multi-step file\n"
          "                            --piece=N keeps one piece of a partitioned file\n"
          "                            --lenient skips constructs the reader cannot represent\n"
          "                            --drop-ghosts removes a .pvtu/.pvd halo (vtkGhostType)\n"
          "                            --appended writes a .vtu's arrays as raw binary in\n"
          "                            one <AppendedData> section (no base64)\n"
          "                            'out_*.vtu' (quoted) or repeated --input fans a\n"
          "                            sequence IN; an out_{step}.vtu output fans one OUT\n"
          "                            (--times, --time-from, --sequence/--no-sequence)\n"
          "                            --color-by NAME colours svg/tikz output by a data\n"
          "                            array (--component --cmap --vmin --vmax\n"
          "                            --nan-color --colorbar)\n"
          "  info (i)                Print mesh info (--fast summarizes from the header)\n"
          "  ascii (a)               Rewrite a file in its ASCII variant (in place)\n"
          "  binary (b)              Rewrite a file in its binary variant (in place)\n"
          "  compress                Compress a mesh file (in place)\n"
          "                            --codec zlib|lz4|zstd for vti/vtu/vtp\n"
          "  decompress              Decompress a mesh file (in place)\n"
          "  quality (q)             Print mesh quality metrics\n"
          "  check                   Quality gate: exit 1 when cells break thresholds\n"
          "                            (--require 'scaled_jacobian >= 0.2; ...', --gate FILE,\n"
          "                            --max-inverted/--max-degenerate; exit 2: not checked)\n"
          "  curvature               Per-vertex mean/Gaussian curvature of a surface\n"
          "  normals                 Point/cell normals of a surface, optionally split at creases\n"
           "  feature-edges           Sharp/open/non-manifold edges of a surface as a line mesh\n"
           "  region-adjacency        Shared facets between Cell regions or cell blocks\n"
           "  find-interface         Match boundary facets of two Cell regions\n"
           "  contact-pairs          Project slave Point regions onto master facets\n"
           "  split-interface        Split point fans along a Side region\n"
          "  hausdorff               Hausdorff distance between two surfaces\n"
          "                            (--max D: nonzero exit when it is exceeded)\n"
          "  periodic                Match the nodes of two regions a transform maps together\n"
          "  repair                  Fix a surface's orientation, holes and bowties\n"
          "  shrinkwrap              Project a mesh's points onto a target surface\n"
          "  sobolev-deform          Filter a displacement field and apply it\n"
          "  extract-surface (surface)  Extract the boundary surface/edges\n"
          "  reorder                 Renumber nodes/elements (RCM / Morton / Hilbert)\n"
          "  diff                    Compare two meshes (nonzero exit if different)\n"
          "  merge                   Merge two or more meshes into one\n"
          "  transform               Affine transform (translate/scale/rotate/matrix/units)\n"
          "  clean                   Weld / prune / de-dup a mesh\n"
          "  crop                    Subset by bounding box, half-space or a data\n"
          "                            predicate: --bbox / --plane / --where 'N < V'\n"
          "  slice                   Planar cross-section (volume->surface, surface->lines)\n"
          "  voxelize                Regular hexahedron grid around a mesh\n"
          "                            exactly one of --resolution/--cell-size;\n"
          "                            --fill all|surface|inside\n"
          "  isosurface              Level set of a scalar point_data field (contours)\n"
          "                            --array NAME --values v1,v2 [--component I]\n"
          "  sdf                     Signed distance field: a grid over a surface,\n"
          "                            filled -- --structure voxel|octree\n"
          "  split                   Partition into multiple files "
          "(type/region/regions/component)\n"
          "                            --by regions is one piece per named Cell region\n"
          "                            (not a partition -- overlapping regions overlap)\n"
          "  regions                 List a mesh's named regions (name/kind/dim/tag/entries),\n"
          "                            or edit them into OUTFILE: --union/--intersection/\n"
          "                            --difference OUT=A,B --rename OLD=NEW\n"
          "                            --retag NAME=TAG[:DIM] --delete NAME\n"
          "  convert-cells           Convert elements (linearize/simplexify/elevate)\n"
          "  subdivide               Polyhedrally refine: one polyhedral child per 3D\n"
          "                            cell face, connected to a new interior point\n"
          "  agglomerate             Polyhedrally coarsen: merge groups of cells into\n"
          "                            single larger polyhedral cells\n"
          "  resample                Resample a sequence onto new times (--times\n"
          "                            START:STOP:STEP | --times-from PATTERN)\n"
          "  refine                  Subdivide cells into same-type children (all, or a\n"
          "                          selected subset with a conforming closure)\n"
          "  undo-green              Restore a transitional (green) cell to its coarse\n"
          "                            parent, read verbatim from COARSE FINE OUTFILE\n"
          "  decimate                Reduce a surface mesh's face count (QEM edge collapse)\n"
          "                            exactly one of --ratio/--target-faces/--max-error\n"
          "  decimate-volume         Reduce a tet mesh's cell count (QEM tet-edge collapse)\n"
          "                            exactly one of --ratio/--target-cells/--max-error\n"
          "  remesh                  Replace a surface's triangulation with a new,\n"
          "                            well-shaped one (ACVD clustering); --num-clusters\n"
          "                            is required, --metric isotropic|quadric|anisotropic\n"
          "  remesh-volume           Retetrahedralize a volume (or closed surface) by\n"
          "                            isosurface stuffing; --cell-size or --resolution\n"
          "                            is required, --warp-fraction trades boundary\n"
          "                            quality for a small chance of non-manifold edges\n"
          "  optimize-volume         ODT-remesh a tetrahedral mesh: relocate vertices and\n"
          "                            flip connectivity (2-3/3-2) to raise element quality\n"
          "  partition               Decompose into N balanced parts (SFC / KaHIP)\n"
          "                            OUT pattern needs {part}, or is a .pvtu/.pvtp index\n"
          "                            over every part; --labels-only writes one file\n"
          "                            with the partition:part cell_data instead\n"
          "  smooth                  Relax node positions (Laplacian / Taubin / ODT)\n"
          "  interpolate             Sample data arrays from a source mesh onto a target\n"
          "                            (nearest / barycentric; --arrays a,b names them)\n"
          "  conservative-interpolate Mass-preserving (overlap-measure weighted) transfer\n"
          "                            from a source mesh onto a target mesh\n"
          "  stats                   Print geometric statistics (bbox/area/volume)\n"
          "  view                    Open a mesh in an interactive viewer\n"
          "  screenshot              Render a mesh to a PNG without a window\n"
          "  data <verb>             Inspect / rename / average / compute on data arrays\n"
          "  pipeline                Run a settings.json operation chain (read -> ops ->\n"
          "                          write; see doc/pipeline.md). --input/--output\n"
          "                          override the paths in the file; --json for a\n"
          "                          machine-readable report. A Pattern/Paths Input\n"
          "                          or a {step} Output runs the chain per step over a\n"
          "                          whole transient dataset (see doc/sequences.md)\n\n"
          "  -v, --version           Display version information\n"
          "  -h, --help              Show this message\n\n"
          "notes: point/cell sets and 'convert -s/-d' are unavailable in the native\n"
          "       CLI (they live only in the Python Mesh); use the Python CLI for those.\n"
          "       view/screenshot need a build with -DMESHIOPLUSPLUS_WITH_POLYSCOPE=ON;\n"
          "       they are listed in every build but otherwise report that.\n";
}

/// "a,b" -> {"a", "b"}, skipping empty entries.
std::vector<std::string> data_split_names(const std::string& rValue) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start <= rValue.size()) {
        const std::size_t comma = rValue.find(',', start);
        const std::size_t end = comma == std::string::npos ? rValue.size() : comma;
        std::string part = rValue.substr(start, end - start);
        // trim
        const std::size_t b = part.find_first_not_of(" \t");
        const std::size_t e = part.find_last_not_of(" \t");
        if (b != std::string::npos)
            out.push_back(part.substr(b, e - b + 1));
        if (comma == std::string::npos)
            break;
        start = comma + 1;
    }
    return out;
}

// ---- JSON shapes shared with the Python CLI (tests/python/test_cli_json.py) --

void cli_regions_json(meshioplusplus::cli::JsonOut& rJson, const std::string& rName,
                      meshioplusplus::RegionKind Kind, int Dim, std::int64_t Tag,
                      std::size_t NumEntries) {
    rJson.BeginObject();
    rJson.Field("name", rName);
    rJson.Field("kind", meshioplusplus::region_kind_name(Kind));
    rJson.Field("dim", Dim);
    rJson.Field("tag", Tag);
    rJson.Field("num_entries", NumEntries);
    rJson.EndObject();
}

void cli_names_json(meshioplusplus::cli::JsonOut& rJson, const char* pKey,
                    const std::vector<std::string>& rNames) {
    rJson.Key(pKey);
    rJson.BeginArray();
    for (const std::string& n : rNames)
        rJson.String(n);
    rJson.EndArray();
}

// `info --fast --json`: the keys of the Python `read_metadata` dict.
void cli_metadata_json(meshioplusplus::cli::JsonOut& rJson,
                       const meshioplusplus::MeshMetadata& rMeta) {
    std::size_t ncells = 0;
    for (const auto& b : rMeta.mCellBlocks)
        ncells += b.mNumCells;
    rJson.BeginObject();
    rJson.Field("num_points", rMeta.mNumPoints);
    rJson.Field("point_dim", rMeta.mPointDim);
    rJson.Field("num_cells", ncells);
    rJson.Key("cell_blocks");
    rJson.BeginArray();
    for (const auto& b : rMeta.mCellBlocks) {
        rJson.BeginObject();
        rJson.Field("type", b.mType);
        rJson.Field("num_cells", b.mNumCells);
        rJson.Field("nodes_per_cell", b.mNodesPerCell);
        rJson.Field("ragged", b.mRagged);
        rJson.EndObject();
    }
    rJson.EndArray();
    cli_names_json(rJson, "point_data_names", rMeta.mPointDataNames);
    cli_names_json(rJson, "cell_data_names", rMeta.mCellDataNames);
    cli_names_json(rJson, "field_data_names", rMeta.mFieldDataNames);
    rJson.Field("fell_back_to_full_read", rMeta.mFellBackToFullRead);
    rJson.Field("format", rMeta.mFormat);
    rJson.Key("time_values");
    rJson.BeginArray();
    for (double t : rMeta.mTimeValues)
        rJson.Number(t);
    rJson.EndArray();
    cli_names_json(rJson, "provenance", rMeta.mProvenance);
    rJson.Field("provenance_recognised", rMeta.mProvenanceRecognised);
    if (rMeta.mHasBBox) {
        for (const auto& [key, v] :
             {std::pair<const char*, const double*>{"bbox_min", rMeta.mBBoxMin},
              {"bbox_max", rMeta.mBBoxMax}}) {
            rJson.Key(key);
            rJson.BeginArray();
            for (int k = 0; k < 3; ++k)
                rJson.Number(v[k]);
            rJson.EndArray();
        }
    }
    rJson.Key("regions");
    rJson.BeginArray();
    for (const auto& r : rMeta.mRegions)
        cli_regions_json(rJson, r.mName, r.mKind, r.mDim, r.mTag, r.mNumEntries);
    rJson.EndArray();
    rJson.EndObject();
}

// `info --json` after a full read: the Python CLI's `_mesh_summary`.
void cli_mesh_summary_json(meshioplusplus::cli::JsonOut& rJson, const Mesh& rMesh) {
    const std::int64_t n = static_cast<std::int64_t>(rMesh.NumPoints());
    std::size_t ncells = 0;
    bool consistent = true;
    std::vector<char> used(static_cast<std::size_t>(n), 0);
    for (const auto cb : rMesh.CellRange()) {
        ncells += cb.NumCells();
        each_block_index(cb, [&](std::int64_t v) {
            if (v > n)
                consistent = false;
            else if (v >= 0 && v < n)
                used[static_cast<std::size_t>(v)] = 1;
        });
    }
    rJson.BeginObject();
    rJson.Field("num_points", n);
    rJson.Field("point_dim", n ? rMesh.PointDim() : std::size_t{0});
    rJson.Field("num_cells", ncells);
    rJson.Key("cell_blocks");
    rJson.BeginArray();
    for (const auto cb : rMesh.CellRange()) {
        rJson.BeginObject();
        rJson.Field("type", std::string(cb.Type()));
        rJson.Field("num_cells", cb.NumCells());
        rJson.Field("nodes_per_cell", cb.IsRagged() ? std::size_t{0} : cb.NodesPerCell());
        rJson.Field("ragged", cb.IsRagged());
        rJson.EndObject();
    }
    rJson.EndArray();
    cli_names_json(rJson, "point_data_names", rMesh.PointDataNames());
    cli_names_json(rJson, "cell_data_names", rMesh.CellDataNames());
    cli_names_json(rJson, "field_data_names", rMesh.FieldDataNames());
    rJson.Key("regions");
    rJson.BeginArray();
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        cli_regions_json(rJson, r.mName, r.mKind, r.mDim, r.mTag, r.NumEntries());
    }
    rJson.EndArray();
    rJson.Field("consistent", consistent);
    rJson.Key("num_unused_points");
    if (consistent)
        rJson.Int(static_cast<std::int64_t>(std::count(used.begin(), used.end(), char{0})));
    else
        rJson.Null();
    rJson.EndObject();
}

/// The transient half of `convert`: fan-in, fan-out or N->N. Kept separate so
/// the single-file path above is physically untouched.
int convert_sequence(const cli_parsed& rParsed, const std::string& rInfile,
                     const std::string& rOutfile, const std::string& rInFmt,
                     const std::string& rOutFmt, const meshioplusplus::ReadOptions& rOpts,
                     bool Ascii, bool Appended, const std::string& rFloatFmt,
                     const std::vector<std::string>& rExtraInputs) {
    meshioplusplus::SequenceInput in;
    // A pre-expanded argv (`--input a --input b`) and a quoted pattern
    // ('out_*.vtu') must reach exactly the same code.
    if (!rExtraInputs.empty()) {
        in.mPaths.push_back(rInfile);
        for (const std::string& extra : rExtraInputs)
            in.mPaths.push_back(extra);
    } else if (rInfile.find('*') != std::string::npos || rInfile.find('?') != std::string::npos) {
        in.mPattern = rInfile;
    } else {
        in.mPaths.push_back(rInfile);
    }
    in.mFormat = rInFmt;
    in.mOptions = rOpts;
    if (has_opt(rParsed, "times"))
        for (const std::string& t : data_split_names(opt_value(rParsed, "times")))
            in.mTimes.push_back(meshioplusplus::detail::stod_c(t));
    in.mTimeFrom = meshioplusplus::sequence_time_from_name(opt_value(rParsed, "time-from"));

    meshioplusplus::SequenceOutput out;
    out.mPath = rOutfile;
    out.mFormat = rOutFmt;
    if (Ascii)
        out.mOptions.mEncoding = meshioplusplus::WriteEncoding::Ascii;
    if (Appended)
        out.mOptions.mEncoding = meshioplusplus::WriteEncoding::RawAppended;
    out.mOptions.mFloatFormat = rFloatFmt;

    const std::vector<meshioplusplus::SequenceEntry> entries = meshioplusplus::sequence_expand(in);
    const meshioplusplus::SequenceMode mode =
        meshioplusplus::sequence_resolve_mode(entries, out, meshioplusplus::SequenceMode::Auto);
    if (mode == meshioplusplus::SequenceMode::FanIn) {
        meshioplusplus::sequence_to_timeseries(in, out);
        if (has_flag(rParsed, "json")) {
            meshioplusplus::cli::JsonOut json(std::cout);
            json.BeginObject();
            json.Field("mode", "fan-in");
            json.Field("num_steps", entries.size());
            json.Field("num_files", std::int64_t{1});
            json.EndObject();
        } else {
            std::cout << "fan-in: " << entries.size() << " step(s) -> 1 file\n";
        }
        return 0;
    }
    // Sequence / fan-out: one output file per entry, through the shared driver.
    meshioplusplus::SequencePipeline pipeline;
    pipeline.mInput = in;
    pipeline.mOutput = out;
    meshioplusplus::run_sequence_pipeline(pipeline);
    const char* what = mode == meshioplusplus::SequenceMode::FanOut ? "fan-out" : "sequence";
    if (has_flag(rParsed, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Field("mode", what);
        json.Field("num_steps", entries.size());
        json.Field("num_files", entries.size());
        json.EndObject();
    } else {
        std::cout << what << ": " << entries.size() << " step(s) -> " << entries.size()
                  << " file(s)\n";
    }
    return 0;
}

int cmd_convert(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"float-format", {"-f"}, true},
                                  {"ascii", {"-a"}, false},
                                  {"appended", {}, false},
                                  {"sets-to-int-data", {"-s"}, false},
                                  {"int-data-to-sets", {"-d"}, false},
                                  {"points-only", {}, false},
                                  {"arrays", {}, true},
                                  {"time-step", {}, true},
                                  {"piece", {}, true},
                                  {"lenient", {}, false},
                                  {"drop-ghosts", {}, false},
                                  {"color-by", {}, true},
                                  {"component", {}, true},
                                  {"cmap", {}, true},
                                  {"vmin", {}, true},
                                  {"vmax", {}, true},
                                  {"nan-color", {}, true},
                                  {"colorbar", {}, false},
                                  {"split-angle", {}, true},
                                  {"up-axis", {}, true},
                                  {"input", {}, true},
                                  {"times", {}, true},
                                  {"time-from", {}, true},
                                  {"sequence", {}, false},
                                  {"no-sequence", {}, false},
                                  {"json", {}, false},
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
    const bool appended = has_flag(p, "appended");
    if (ascii && appended)
        throw std::runtime_error("--ascii and --appended are mutually exclusive");

    // Selective read: --points-only drops every data array, --arrays keeps only
    // the named ones. Connectivity is kept either way, so the output is still a
    // usable mesh -- only the data is narrowed.
    meshioplusplus::ReadOptions opts;
    opts.mPointsOnly = has_flag(p, "points-only");
    if (has_opt(p, "arrays")) {
        if (opts.mPointsOnly)
            throw std::runtime_error("--points-only and --arrays are mutually exclusive");
        opts.mDataArrays = data_split_names(opt_value(p, "arrays"));
    }
    // --time-step is not a narrowing option: no post-filter can recover a step
    // that was never read, so it goes to the reader or nowhere.
    if (has_opt(p, "time-step"))
        opts.mTimeStep = std::stoi(opt_value(p, "time-step"));
    // --piece keeps one partition / composite block of a partitioned file (VTKHDF)
    // instead of the merged mesh. Like --time-step it cannot be emulated after the
    // fact, so it goes to the reader or fails there.
    if (has_opt(p, "piece")) {
        opts.mPiece = std::stoll(opt_value(p, "piece"));
        opts.mPieceSet = true;
    }
    // --lenient downgrades "this reader cannot represent construct X" to a
    // warning plus a skip. Not "ignore all errors": a malformed file still
    // fails. There is no Python fallback here, so this is what makes a
    // production .mdpa readable at all from the native CLI.
    opts.mLenient = has_flag(p, "lenient");
    // --drop-ghosts removes the ghost cells (halo) of a partitioned .pvtu/.pvtp/.pvd
    // and the points only they used; every other reader ignores it, like --lenient.
    if (has_flag(p, "drop-ghosts"))
        opts.mGhosts = meshioplusplus::GhostPolicy::Drop;

    // Transient sequences: a quoted glob, repeated --input, or a {step}/{index}
    // output. `--input` is read through the parser's `multi` map (added for the
    // `data` verbs) so a shell that already expanded a glob is served without
    // changing convert's two-positional shape.
    const std::vector<std::string>& extra_inputs = opt_values(p, "input");
    if (has_flag(p, "sequence") && has_flag(p, "no-sequence"))
        throw std::runtime_error("--sequence and --no-sequence are mutually exclusive");
    bool sequence = has_flag(p, "sequence");
    if (!has_flag(p, "no-sequence") && !sequence) {
        sequence = !extra_inputs.empty() || infile.find('*') != std::string::npos ||
                   infile.find('?') != std::string::npos ||
                   meshioplusplus::sequence_pattern_has_token(outfile);
        // A multi-step input aimed at a single-step output must not quietly
        // become step 0: route it here so the driver refuses by name. An
        // explicit --time-step IS a deliberate single-step selection.
        if (!sequence && !has_opt(p, "time-step"))
            sequence = meshioplusplus::sequence_num_steps(infile, in_fmt) > 1;
    }
    if (sequence)
        return convert_sequence(p, infile, outfile, in_fmt, out_fmt, opts, ascii, appended,
                                float_fmt, extra_inputs);

    // Data-driven colouring (svg/tikz/gltf) and the glTF-only options. Validated
    // before the read so a bad flag combination fails immediately rather than
    // after loading a big mesh.
    const bool color = has_opt(p, "color-by");
    const bool gltf_flags = has_opt(p, "split-angle") || has_opt(p, "up-axis");
    std::string target_fmt;
    if (color || gltf_flags)
        target_fmt = meshioplusplus::resolve_write_format(outfile, out_fmt);
    if (gltf_flags && target_fmt != "gltf")
        throw std::runtime_error("--split-angle/--up-axis only apply to glTF output, not '" +
                                 target_fmt + "'");
    if (!color) {
        for (const char* flag : {"component", "cmap", "vmin", "vmax", "nan-color"})
            if (has_opt(p, flag))
                throw std::runtime_error(std::string("--") + flag + " requires --color-by");
        if (has_flag(p, "colorbar"))
            throw std::runtime_error("--colorbar requires --color-by");
    } else {
        if (!cli_is_colorable_format(target_fmt))
            throw std::runtime_error(
                "--color-by is only supported for svg/tikz/gltf output, not '" + target_fmt + "'");
        if (ascii || appended)
            throw std::runtime_error(std::string(ascii ? "--ascii" : "--appended") +
                                     " has no meaning for " + target_fmt + " output");
        if (target_fmt == "gltf" && has_flag(p, "colorbar"))
            throw std::runtime_error("--colorbar has no meaning for gltf output (svg/tikz only)");
    }

    Mesh mesh = read_mesh_cli(infile, in_fmt, opts);

    std::optional<int> component;
    std::optional<double> vmin;
    std::optional<double> vmax;
    if (color) {
        if (has_opt(p, "component"))
            component = std::stoi(opt_value(p, "component"));
        if (has_opt(p, "vmin"))
            vmin = meshioplusplus::detail::stod_c(opt_value(p, "vmin"));
        if (has_opt(p, "vmax"))
            vmax = meshioplusplus::detail::stod_c(opt_value(p, "vmax"));
    }
    const std::string cmap = has_opt(p, "cmap") ? opt_value(p, "cmap") : "viridis";

    if (target_fmt == "gltf" && (color || gltf_flags)) {
        meshioplusplus::GltfWriteOptions gltf;
        if (has_opt(p, "split-angle"))
            gltf.mSplitAngle = meshioplusplus::detail::stod_c(opt_value(p, "split-angle"));
        if (has_opt(p, "up-axis"))
            gltf.mUpAxis = meshioplusplus::gltf_up_axis_from_name(opt_value(p, "up-axis"));
        if (color) {
            gltf.mColorBy = opt_value(p, "color-by");
            gltf.mComponent = component;
            gltf.mCmap = cmap;
            gltf.mVMin = vmin;
            gltf.mVMax = vmax;
            if (has_opt(p, "nan-color"))
                gltf.mNanColor = opt_value(p, "nan-color");
        }
        meshioplusplus::write_gltf(outfile, mesh, gltf);
    } else if (color) {
        write_colored_variant(outfile, mesh, target_fmt, opt_value(p, "color-by"), component, cmap,
                              vmin, vmax, opt_value(p, "nan-color"), has_flag(p, "colorbar"));
    } else if (ascii) {
        std::string fmt = meshioplusplus::resolve_write_format(outfile, out_fmt);
        if (!write_binary_variant(outfile, mesh, fmt, /*binary=*/false, float_fmt))
            throw std::runtime_error("format '" + fmt + "' has no ASCII variant");
    } else if (appended) {
        std::string fmt = meshioplusplus::resolve_format(outfile, out_fmt);
        meshioplusplus::WriteOptions wopts;
        wopts.mEncoding = meshioplusplus::WriteEncoding::RawAppended;
        std::string why;
        if (!meshioplusplus::registry_write_supports(fmt, wopts, why))
            throw std::runtime_error("--appended: " + why);
        if (!float_fmt.empty())
            std::cerr << "note: --float-format only affects ASCII output (--ascii)\n";
        meshioplusplus::registry_write_ex(outfile, mesh, fmt, wopts);
    } else {
        if (!float_fmt.empty())
            std::cerr << "note: --float-format only affects ASCII output (--ascii)\n";
        write_mesh_cli(outfile, mesh, out_fmt);
    }
    if (has_flag(p, "json")) {
        // The Python CLI's convert summary.
        std::size_t ncells = 0;
        for (const auto cb : mesh.CellRange())
            ncells += cb.NumCells();
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Field("input", infile);
        json.Field("output", outfile);
        json.Field("num_points", mesh.NumPoints());
        json.Field("num_cells", ncells);
        json.Key("cell_blocks");
        json.BeginArray();
        for (const auto cb : mesh.CellRange()) {
            json.BeginObject();
            json.Field("type", std::string(cb.Type()));
            json.Field("num_cells", cb.NumCells());
            json.EndObject();
        }
        json.EndArray();
        cli_names_json(json, "point_data_names", mesh.PointDataNames());
        cli_names_json(json, "cell_data_names", mesh.CellDataNames());
        cli_names_json(json, "field_data_names", mesh.FieldDataNames());
        json.Field("num_regions", mesh.NumRegions());
        json.EndObject();
    }
    return 0;
}

/// Print the header-only summary behind `info --fast`.
void print_metadata_summary(const meshioplusplus::MeshMetadata& rMeta) {
    std::cout << "<meshio++ mesh summary>\n";
    std::cout << "  Format: " << (rMeta.mFormat.empty() ? "unknown" : rMeta.mFormat) << "\n";
    std::cout << "  Number of points: " << rMeta.mNumPoints << "\n";
    if (rMeta.mCellBlocks.empty()) {
        std::cout << "  No cells.\n";
    } else {
        std::cout << "  Number of cells:\n";
        for (const auto& block : rMeta.mCellBlocks)
            std::cout << "    " << block.mType << ": " << block.mNumCells << "\n";
    }
    const std::pair<const char*, const std::vector<std::string>*> sections[] = {
        {"Point data", &rMeta.mPointDataNames},
        {"Cell data", &rMeta.mCellDataNames},
        {"Field data", &rMeta.mFieldDataNames},
    };
    for (const auto& section : sections) {
        if (section.second->empty())
            continue;
        std::cout << "  " << section.first << ": ";
        for (std::size_t i = 0; i < section.second->size(); ++i)
            std::cout << (i ? ", " : "") << (*section.second)[i];
        std::cout << "\n";
    }
    // Only worth printing when there is a choice to make: a single-step file
    // gives a caller nothing to pass to --time-step.
    if (rMeta.mTimeValues.size() > 1) {
        std::cout << "  Time steps: " << rMeta.mTimeValues.size() << " [";
        const std::size_t shown = std::min<std::size_t>(rMeta.mTimeValues.size(), 8);
        for (std::size_t i = 0; i < shown; ++i)
            std::cout << (i ? ", " : "") << rMeta.mTimeValues[i];
        if (rMeta.mTimeValues.size() > shown)
            std::cout << ", ...";
        std::cout << "]\n";
    }
    if (!rMeta.mRegions.empty()) {
        std::cout << "  Regions (" << rMeta.mRegions.size() << "):\n";
        for (const auto& r : rMeta.mRegions) {
            std::cout << "    " << r.mName << " (" << meshioplusplus::region_kind_name(r.mKind)
                      << ", " << r.mNumEntries << " entries";
            if (r.mTag >= 0)
                std::cout << ", tag=" << r.mTag;
            std::cout << ")\n";
        }
    }
    if (!rMeta.mProvenance.empty()) {
        // Say whose block it is: a file can carry a comment that merely looks
        // like a header, and reporting that as meshio++'s would be a lie.
        std::cout << "  Provenance"
                  << (rMeta.mProvenanceRecognised ? "" : " (not written by meshio++)") << ":\n";
        for (const std::string& line : rMeta.mProvenance)
            std::cout << "    " << line << "\n";
    }
    if (rMeta.mHasBBox) {
        std::cout << "  Bounding box: [" << rMeta.mBBoxMin[0] << ", " << rMeta.mBBoxMin[1] << ", "
                  << rMeta.mBBoxMin[2] << "] - [" << rMeta.mBBoxMax[0] << ", " << rMeta.mBBoxMax[1]
                  << ", " << rMeta.mBBoxMax[2] << "]\n";
    }
    // Say plainly when "fast" was not, rather than implying a saving that did
    // not happen.
    if (rMeta.mFellBackToFullRead)
        std::cout << "  (no header-only path for this format; the file was read in full)\n";
}

int cmd_info(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs,
                       {{"input-format", {"-i"}, true}, {"fast", {}, false}, {"json", {}, false}});
    if (has_flag(p, "fast")) {
        if (p.positionals.size() != 1)
            throw std::runtime_error("info requires exactly one INFILE");
        std::string fmt;
        try {
            fmt = meshioplusplus::resolve_format(p.positionals[0], opt_value(p, "input-format"));
        } catch (const ReadError&) {
            fmt = meshioplusplus::sniff_format(p.positionals[0]);
            if (fmt.empty())
                throw;
        }
        const auto meta = meshioplusplus::registry_read_metadata(p.positionals[0], fmt, {});
        if (has_flag(p, "json")) {
            meshioplusplus::cli::JsonOut json(std::cout);
            cli_metadata_json(json, meta);
        } else {
            print_metadata_summary(meta);
        }
        return 0;
    }
    if (p.positionals.size() != 1)
        throw std::runtime_error("info requires exactly INFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        cli_mesh_summary_json(json, mesh);
        return 0;
    }

    // Mirror Mesh.__repr__, named regions included (doc/regions.md).
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
    // Named groups, one line per kind, matching Mesh.__repr__'s order.
    {
        std::vector<std::string> point_sets, cell_sets, side_sets;
        for (std::size_t i = 0; i < mesh.NumRegions(); ++i) {
            const meshioplusplus::Region& r = mesh.Region(i);
            switch (r.mKind) {
                case meshioplusplus::RegionKind::Point:
                    point_sets.push_back(r.mName);
                    break;
                case meshioplusplus::RegionKind::Cell:
                    cell_sets.push_back(r.mName);
                    break;
                case meshioplusplus::RegionKind::Side:
                    side_sets.push_back(r.mName);
                    break;
            }
        }
        if (!point_sets.empty())
            std::cout << "  Point sets: " << join(point_sets) << "\n";
        if (!cell_sets.empty())
            std::cout << "  Cell sets: " << join(cell_sets) << "\n";
        if (!side_sets.empty())
            std::cout << "  Side sets: " << join(side_sets) << "\n";
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

// One region edit from the CLI spelling (the Python CLI's `parse_edit`):
// union/intersection/difference OUT=A,B[,...]; rename OLD=NEW; retag
// NAME=TAG[:DIM]; delete NAME. A name may be prefixed point:/cell:/side:.
meshioplusplus::RegionEdit cli_region_edit(const std::string& rOp, const std::string& rText,
                                           bool KeepInputs) {
    auto selector = [](const std::string& rSpec) {
        meshioplusplus::RegionSelector sel;
        const auto colon = rSpec.find(':');
        const std::string kind = colon == std::string::npos ? "" : rSpec.substr(0, colon);
        if (kind == "point" || kind == "cell" || kind == "side") {
            sel.mKind = static_cast<std::int32_t>(meshioplusplus::region_kind_from_name(kind));
            sel.mName = rSpec.substr(colon + 1);
        } else {
            sel.mName = rSpec;
        }
        return sel;
    };
    meshioplusplus::RegionEdit e;
    e.mOp = meshioplusplus::region_op_from_name(rOp);
    if (e.mOp == meshioplusplus::RegionOp::Delete) {
        e.mInputs.push_back(selector(rText));
        return e;
    }
    const auto eq = rText.find('=');
    if (eq == std::string::npos || eq == 0 || eq + 1 == rText.size())
        throw std::runtime_error("regions: --" + rOp + " expects NAME=VALUE, got '" + rText + "'");
    const std::string left = rText.substr(0, eq);
    const std::string right = rText.substr(eq + 1);
    switch (e.mOp) {
        case meshioplusplus::RegionOp::Union:
        case meshioplusplus::RegionOp::Intersection:
        case meshioplusplus::RegionOp::Difference: {
            std::size_t start = 0;
            while (start <= right.size()) {
                const auto comma = right.find(',', start);
                const std::string name = right.substr(start, comma - start);
                if (!name.empty())
                    e.mInputs.push_back(selector(name));
                if (comma == std::string::npos)
                    break;
                start = comma + 1;
            }
            e.mOutputName = left;
            e.mKeepInputs = KeepInputs;
            break;
        }
        case meshioplusplus::RegionOp::Rename:
            e.mInputs.push_back(selector(left));
            e.mOutputName = right;
            break;
        default: {
            e.mInputs.push_back(selector(left));
            const auto colon = right.find(':');
            const std::string tag = right.substr(0, colon);
            try {
                if (!tag.empty())
                    e.mOutputTag = std::stoll(tag);
                if (colon != std::string::npos)
                    e.mOutputDim = std::stoll(right.substr(colon + 1));
            } catch (const std::exception&) {
                throw std::runtime_error("regions: --retag expects NAME=TAG[:DIM], got '" + rText +
                                         "'");
            }
        }
    }
    return e;
}

void cli_print_regions_json(meshioplusplus::cli::JsonOut& rJson, const Mesh& rMesh) {
    rJson.BeginArray();
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        rJson.BeginObject();
        rJson.Field("name", r.mName);
        rJson.Field("kind", meshioplusplus::region_kind_name(r.mKind));
        rJson.Field("dim", r.mDim);
        rJson.Field("tag", r.mTag);
        rJson.Field("num_entries", r.NumEntries());
        rJson.EndObject();
    }
    rJson.EndArray();
}

int cmd_regions(const std::vector<std::string>& rArgs) {
    static const char* const kEdits[] = {"union",  "intersection", "difference",
                                         "rename", "retag",        "delete"};
    std::vector<cli_opt_spec> specs = {{"input-format", {"-i"}, true},
                                       {"output-format", {"-o"}, true},
                                       {"json", {}, false},
                                       {"drop-inputs", {}, false}};
    for (const char* e : kEdits)
        specs.push_back({e, {}, true});
    auto p = cli_parse(rArgs, specs);
    if (!p.ordered.empty() && std::any_of(p.ordered.begin(), p.ordered.end(), [](const auto& rKv) {
            return rKv.first != "input-format" && rKv.first != "output-format";
        })) {
        if (p.positionals.size() != 2)
            throw std::runtime_error("regions: edits need exactly INFILE and OUTFILE");
        std::vector<meshioplusplus::RegionEdit> edits;
        for (const auto& [name, value] : p.ordered)
            if (name != "input-format" && name != "output-format")
                edits.push_back(cli_region_edit(name, value, !has_flag(p, "drop-inputs")));
        Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
        Mesh out = meshioplusplus::edit_regions(mesh, edits);
        write_mesh_cli(p.positionals[1], out, opt_value(p, "output-format"));
        if (has_flag(p, "json")) {
            meshioplusplus::cli::JsonOut json(std::cout);
            cli_print_regions_json(json, out);
            return 0;
        }
        std::cout << "<meshio++ mesh regions> (" << out.NumRegions() << ") after " << edits.size()
                  << " edit(s)\n";
        for (std::size_t i = 0; i < out.NumRegions(); ++i) {
            const meshioplusplus::Region& r = out.Region(i);
            std::cout << "  " << r.mName << " (" << meshioplusplus::region_kind_name(r.mKind)
                      << ", " << r.NumEntries() << " entries";
            if (r.mDim >= 0)
                std::cout << ", dim=" << r.mDim;
            if (r.mTag >= 0)
                std::cout << ", tag=" << r.mTag;
            std::cout << ")\n";
        }
        return 0;
    }
    if (p.positionals.size() != 1)
        throw std::runtime_error("regions requires exactly INFILE (and OUTFILE with edits)");

    // Goes through registry_read_metadata rather than a full read: a native
    // metadata path costs nothing extra to report regions from an
    // already-read mesh (Exodus, and every fallback path), and this stays
    // cheap on any format -- nothing here needs the connectivity.
    std::string fmt;
    try {
        fmt = meshioplusplus::resolve_format(p.positionals[0], opt_value(p, "input-format"));
    } catch (const ReadError&) {
        fmt = meshioplusplus::sniff_format(p.positionals[0]);
        if (fmt.empty())
            throw;
    }
    const meshioplusplus::MeshMetadata meta =
        meshioplusplus::registry_read_metadata(p.positionals[0], fmt, {});

    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginArray();
        for (const auto& r : meta.mRegions) {
            json.BeginObject();
            json.Field("name", r.mName);
            json.Field("kind", meshioplusplus::region_kind_name(r.mKind));
            json.Field("dim", r.mDim);
            json.Field("tag", r.mTag);
            json.Field("num_entries", r.mNumEntries);
            json.EndObject();
        }
        json.EndArray();
        return 0;
    }

    if (meta.mRegions.empty()) {
        std::cout << "<meshio++ mesh regions>\n  No regions.\n";
        if (meta.mFellBackToFullRead)
            std::cout << "  (the full mesh was read; this format may simply carry none)\n";
        return 0;
    }
    std::cout << "<meshio++ mesh regions> (" << meta.mRegions.size() << ")\n";
    for (const auto& r : meta.mRegions) {
        std::cout << "  " << r.mName << " (" << meshioplusplus::region_kind_name(r.mKind) << ", "
                  << r.mNumEntries << " entries";
        if (r.mDim >= 0)
            std::cout << ", dim=" << r.mDim;
        if (r.mTag >= 0)
            std::cout << ", tag=" << r.mTag;
        std::cout << ")\n";
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
    if (fmt == "pcd")
        write_pcd_in_place(
            infile, mesh,
            binary ? meshioplusplus::PcdData::Binary : meshioplusplus::PcdData::Ascii);
    else if (!write_binary_variant(infile, mesh, fmt, binary, ""))
        throw std::runtime_error("don't know how to write '" + fmt + "' as " +
                                 (binary ? "binary" : "ASCII"));
    std::cout << "File size after: " << file_size_mb(infile) << "\n";
    return 0;
}

/// `--codec` value -> VtkCodec, rejecting anything that is not a block codec.
meshioplusplus::detail::VtkCodec codec_from_name(const std::string& rName) {
    using meshioplusplus::detail::VtkCodec;
    if (rName == "zlib")
        return VtkCodec::Zlib;
    if (rName == "lz4")
        return VtkCodec::LZ4;
    if (rName == "zstd")
        return VtkCodec::ZSTD;
    throw std::runtime_error("unknown --codec '" + rName + "' (expected zlib, lz4 or zstd)");
}

int cmd_compress(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(
        rArgs, {{"input-format", {"-i"}, true}, {"max", {"-max"}, false}, {"codec", {}, true}});
    if (p.positionals.size() != 1)
        throw std::runtime_error("compress requires exactly INFILE");
    const std::string& infile = p.positionals[0];
    bool max = has_flag(p, "max");
    int gzip = max ? 9 : 4;
    std::cout << "File size before: " << file_size_mb(infile) << "\n";
    Mesh mesh = read_mesh_cli(infile, opt_value(p, "input-format"));
    std::string fmt = meshioplusplus::resolve_format(infile, opt_value(p, "input-format"));

    // --codec only means something where a block codec is actually chosen.
    // Accepting and ignoring it elsewhere would be the worst outcome: the user
    // would believe they got zstd and silently get gzip (or plain binary).
    const bool has_codec = has_opt(p, "codec");
    if (has_codec && fmt != "vti" && fmt != "vtu" && fmt != "vtp")
        throw std::runtime_error("--codec is not applicable to '" + fmt +
                                 "'; it selects the VTK XML block codec and only "
                                 "vti/vtu/vtp have one");

    if (fmt == "ansys" || fmt == "gmsh" || fmt == "ply" || fmt == "stl" || fmt == "vtk") {
        write_binary_variant(infile, mesh, fmt, /*binary=*/true, "");
    } else if (fmt == "pcd") {
        // LZF over the struct-of-arrays layout: the one compression PCD defines.
        write_pcd_in_place(infile, mesh, meshioplusplus::PcdData::BinaryCompressed);
    } else if (fmt == "vtu" || fmt == "vtp") {
        meshioplusplus::detail::VtkCodec codec = meshioplusplus::detail::VtkCodec::Zlib;
        if (has_codec)
            codec = codec_from_name(opt_value(p, "codec"));
        else if (max)
            std::cerr << "note: the native CLI has no lzma variant; using zlib\n";
        if (fmt == "vtu")
            meshioplusplus::write_vtu_codec(infile, mesh, /*binary=*/true, codec);
        else
            meshioplusplus::write_vtp_codec(infile, mesh, /*binary=*/true, codec);
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

    if (fmt == "pcd") {
        write_pcd_in_place(infile, mesh, meshioplusplus::PcdData::Binary);
    } else if (fmt == "vtu") {
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
    meshioplusplus::detail::snprintf_c(buf, sizeof(buf), "%8.4f", v);
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
    auto p = cli_parse(
        rArgs, {{"input-format", {"-i"}, true}, {"output", {"-o"}, true}, {"json", {}, false}});
    if (p.positionals.size() != 1)
        throw std::runtime_error("quality requires exactly INFILE");
    const std::string& infile = p.positionals[0];
    Mesh mesh = read_mesh_cli(infile, opt_value(p, "input-format"));
    auto report = meshioplusplus::compute_quality(mesh);
    if (has_flag(p, "json")) {
        const std::string output = opt_value(p, "output");
        if (!output.empty())
            write_mesh_cli(output, meshioplusplus::attach_quality(mesh), "");
        // compute_quality's Python dict, without the per-cell arrays.
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Field("num_cells", report.mNumCells);
        json.Field("num_inverted", report.mNumInverted);
        json.Field("num_degenerate", report.mNumDegenerate);
        json.Key("metrics");
        json.BeginObject();
        for (const auto& [name, m] : report.mMetrics) {
            json.Key(name);
            json.BeginObject();
            json.Field("min", m.mMin);
            json.Field("max", m.mMax);
            json.Field("mean", m.mMean);
            json.Field("count", m.mCount);
            json.Key("histogram");
            json.BeginArray();
            for (auto h : m.mHistogram)
                json.Int(h);
            json.EndArray();
            json.EndObject();
        }
        json.EndObject();
        json.EndObject();
        return 0;
    }

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

int cmd_curvature(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"no-mean", {}, false},
                                  {"no-gaussian", {}, false},
                                  {"dual-area", {}, true},
                                  {"include-boundary", {}, false},
                                  {"record-area", {}, false},
                                  {"record-principal", {}, false},
                                  {"region", {}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("curvature requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::CurvatureOptions options;
    options.mMean = !has_flag(p, "no-mean");
    options.mGaussian = !has_flag(p, "no-gaussian");
    const std::string dual = opt_value(p, "dual-area");
    options.mDualArea =
        meshioplusplus::curvature_dual_area_from_name(dual.empty() ? "mixed-voronoi" : dual);
    options.mIncludeBoundary = has_flag(p, "include-boundary");
    options.mRecordArea = has_flag(p, "record-area");
    options.mRecordPrincipal = has_flag(p, "record-principal");
    options.mRegion = opt_value(p, "region");

    meshioplusplus::CurvatureResult r = meshioplusplus::compute_curvature(mesh, options);

    if (!has_flag(p, "quiet")) {
        std::cout << "curvature (" << meshioplusplus::curvature_dual_area_name(options.mDualArea)
                  << ")\n";
        // The Gauss-Bonnet oracle, printed beside its expected value: on a
        // closed surface the defects sum to 2*pi*chi whatever the
        // tessellation, so a reader can check the result without knowing the
        // estimator.
        std::cout << "  total angle defect:       " << std::fixed << std::setprecision(6)
                  << r.mTotalAngleDefect << "   (4*pi = " << (4.0 * 3.14159265358979323846)
                  << " for a closed genus-0 surface)\n";
        std::cout << "  boundary vertices (NaN):  " << r.mNumBoundary << "\n";
        std::cout << "  isolated vertices (NaN):  " << r.mNumIsolated << "\n";
        std::cout << "  degenerate triangles:     " << r.mNumDegenerate << "\n";
        std::cout << "  surface: ";
        if (r.mQuality.mWatertight)
            std::cout << "watertight\n";
        else
            std::cout << r.mQuality.mBoundaryEdges << " boundary / " << r.mQuality.mNonManifoldEdges
                      << " non-manifold / " << r.mQuality.mInconsistentPairs
                      << " inconsistent edge(s)\n";
        std::cout.unsetf(std::ios::floatfield);
    }

    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_normals(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"cell", {}, false},
                                  {"no-point", {}, false},
                                  {"weight", {}, true},
                                  {"split-angle", {}, true},
                                  {"record-parent-ids", {}, false},
                                  {"region", {}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("normals requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::NormalsOptions options;
    options.mPointNormals = !has_flag(p, "no-point");
    options.mCellNormals = has_flag(p, "cell");
    const std::string weight = opt_value(p, "weight");
    options.mWeight = meshioplusplus::sdf_weight_from_name(weight.empty() ? "angle" : weight);
    if (has_opt(p, "split-angle")) {
        options.mSplit = true;
        options.mSplitAngle = meshioplusplus::detail::stod_c(opt_value(p, "split-angle"));
    }
    options.mRecordParentIds = has_flag(p, "record-parent-ids");
    options.mRegion = opt_value(p, "region");

    meshioplusplus::NormalsResult r = meshioplusplus::compute_normals(mesh, options);

    if (!has_flag(p, "quiet")) {
        std::cout << "normals (" << (weight.empty() ? "angle" : weight) << ")\n";
        std::cout << "  points added by the split: " << r.mNumAddedPoints << "\n";
        std::cout << "  points split:              " << r.mNumSplitPoints << "\n";
        std::cout << "  isolated points (NaN):     " << r.mNumIsolated << "\n";
        std::cout << "  undefined points (NaN):    " << r.mNumUndefined << "\n";
        std::cout << "  degenerate triangles:      " << r.mNumDegenerate << "\n";
        std::cout << "  surface: ";
        if (r.mQuality.mWatertight)
            std::cout << "watertight\n";
        else
            std::cout << r.mQuality.mBoundaryEdges << " boundary / " << r.mQuality.mNonManifoldEdges
                      << " non-manifold / " << r.mQuality.mInconsistentPairs
                      << " inconsistent edge(s)\n";
    }

    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

void cli_quality_gate_json(meshioplusplus::cli::JsonOut& rJson,
                           const meshioplusplus::QualityGateResult& rR) {
    rJson.BeginObject();
    rJson.Field("passed", rR.mPassed);
    rJson.Field("num_cells", rR.mReport.mNumCells);
    rJson.Field("num_inverted", rR.mReport.mNumInverted);
    rJson.Field("num_degenerate", rR.mReport.mNumDegenerate);
    rJson.Key("checks");
    rJson.BeginArray();
    for (const meshioplusplus::QualityCheck& c : rR.mChecks) {
        rJson.BeginObject();
        rJson.Field("name", c.mName);
        rJson.Field("metric", c.mMetric);
        rJson.Field("min", c.mMin);
        rJson.Field("max", c.mMax);
        rJson.Field("max_fraction", c.mMaxFraction);
        rJson.Field("evaluated", c.mEvaluated);
        rJson.Field("violations", c.mViolations);
        rJson.Field("fraction", c.mFraction);
        rJson.Field("worst", c.mWorst);
        rJson.Field("worst_cell", c.mWorstCell);
        rJson.Field("passed", c.mPassed);
        rJson.EndObject();
    }
    rJson.EndArray();
    rJson.EndObject();
}

int cmd_check(const std::vector<std::string>& rArgs) {
    // 0 pass, 1 fail, 2 could not check (unreadable file, bad specification):
    // a CI job must tell "the mesh is bad" from "the check did not run".
    meshioplusplus::QualityGateResult r;
    bool json = false;
    try {
        auto p = cli_parse(rArgs, {
                                      {"input-format", {"-i"}, true},
                                      {"require", {"-r"}, true},
                                      {"gate", {}, true},
                                      {"max-inverted", {}, true},
                                      {"max-degenerate", {}, true},
                                      {"json", {}, false},
                                  });
        if (p.positionals.size() != 1)
            throw std::runtime_error("check requires exactly INFILE");
        json = has_flag(p, "json");
        meshioplusplus::QualityGateOptions options;
        std::vector<std::string> specs = opt_values(p, "require");
        if (has_opt(p, "gate")) {
            auto in = meshioplusplus::detail::make_classic_ifstream(opt_value(p, "gate"));
            if (!in)
                throw std::runtime_error("check: cannot open '" + opt_value(p, "gate") + "'");
            specs.emplace_back(std::istreambuf_iterator<char>(in),
                               std::istreambuf_iterator<char>());
        }
        for (const std::string& spec : specs) {
            const auto t = meshioplusplus::parse_quality_thresholds(spec);
            options.mThresholds.insert(options.mThresholds.end(), t.begin(), t.end());
        }
        if (has_opt(p, "max-inverted"))
            options.mMaxInverted = std::stoll(opt_value(p, "max-inverted"));
        if (has_opt(p, "max-degenerate"))
            options.mMaxDegenerate = std::stoll(opt_value(p, "max-degenerate"));
        Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
        r = meshioplusplus::check_quality(mesh, options);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 2;
    }
    if (json) {
        meshioplusplus::cli::JsonOut out(std::cout);
        cli_quality_gate_json(out, r);
    } else {
        std::cout << meshioplusplus::quality_gate_summary(r);
    }
    return r.mPassed ? 0 : 1;
}

int cmd_feature_edges(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"angle", {}, true},
                                  {"no-feature", {}, false},
                                  {"no-boundary", {}, false},
                                  {"no-non-manifold", {}, false},
                                  {"no-inconsistent", {}, false},
                                  {"region", {}, true},
                                  {"json", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("feature-edges requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    meshioplusplus::FeatureEdgeOptions options;
    if (has_opt(p, "angle"))
        options.mFeatureAngleDeg = meshioplusplus::detail::stod_c(opt_value(p, "angle"));
    options.mFeature = !has_flag(p, "no-feature");
    options.mBoundary = !has_flag(p, "no-boundary");
    options.mNonManifold = !has_flag(p, "no-non-manifold");
    options.mInconsistent = !has_flag(p, "no-inconsistent");
    options.mRegion = opt_value(p, "region");
    meshioplusplus::FeatureEdgeResult r = meshioplusplus::feature_edges(mesh, options);
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    const std::size_t edges = r.mMesh.NumCellBlocks() ? r.mMesh.Cells(0).NumCells() : 0;
    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Field("edges", edges);
        json.Field("num_feature", r.mNumFeature);
        json.Field("num_boundary", r.mNumBoundary);
        json.Field("num_non_manifold", r.mNumNonManifold);
        json.Field("num_inconsistent", r.mNumInconsistent);
        json.EndObject();
    } else if (!has_flag(p, "quiet")) {
        char deg[32];
        meshioplusplus::detail::snprintf_c(deg, sizeof(deg), "%g", options.mFeatureAngleDeg);
        std::cout << "feature edges (" << deg << " degrees): " << edges << " written\n";
        std::cout << "  sharp:        " << r.mNumFeature << "\n";
        std::cout << "  boundary:     " << r.mNumBoundary << "\n";
        std::cout << "  non-manifold: " << r.mNumNonManifold << "\n";
        std::cout << "  inconsistent: " << r.mNumInconsistent << "\n";
    }
    return 0;
}

int cmd_region_adjacency(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"regions", {}, true},
                                  {"json", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("region-adjacency requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    std::vector<meshioplusplus::RegionSelector> selectors;
    for (const std::string& value : opt_values(p, "regions")) {
        std::size_t start = 0;
        while (start <= value.size()) {
            const std::size_t comma = value.find(',', start);
            const std::string name = value.substr(start, comma - start);
            if (!name.empty()) {
                meshioplusplus::RegionSelector selector;
                selector.mName = name;
                selector.mKind = static_cast<std::int32_t>(meshioplusplus::RegionKind::Cell);
                selectors.push_back(std::move(selector));
            }
            if (comma == std::string::npos)
                break;
            start = comma + 1;
        }
    }
    const Mesh out = meshioplusplus::region_adjacency(mesh, selectors);
    write_mesh_cli(p.positionals[1], out, opt_value(p, "output-format"));
    std::int64_t facets = 0;
    double area = 0.0;
    for (std::size_t b = 0; b < out.NumCellBlocks(); ++b) {
        const NDArray& measures = out.CellData("interface:measure", b);
        facets += static_cast<std::int64_t>(out.Cells(b).NumCells());
        for (std::size_t i = 0; i < measures.Size(); ++i)
            area += meshioplusplus::detail::read_double(measures, i);
    }
    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Field("facets", facets);
        json.Field("measure", area);
        json.EndObject();
    } else if (!has_flag(p, "quiet")) {
        char measure[64];
        meshioplusplus::detail::snprintf_c(measure, sizeof(measure), "%g", area);
        std::cout << "region adjacency: " << facets << " shared facets, measure " << measure
                  << " written\n";
    }
    return 0;
}

void cli_json_side_entries(meshioplusplus::cli::JsonOut& rJson, const NDArray& rEntries) {
    rJson.BeginArray();
    const std::size_t count = rEntries.Shape().empty() ? 0 : rEntries.Shape()[0];
    for (std::size_t i = 0; i < count; ++i) {
        rJson.BeginArray();
        rJson.Int(meshioplusplus::detail::read_int(rEntries, 2 * i));
        rJson.Int(meshioplusplus::detail::read_int(rEntries, 2 * i + 1));
        rJson.EndArray();
    }
    rJson.EndArray();
}

int cmd_find_interface(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"format-b", {}, true},
                                  {"output-format", {"-o"}, true},
                                  {"region-a", {}, true},
                                  {"region-b", {}, true},
                                  {"mesh-b", {}, true},
                                  {"mode", {}, true},
                                  {"master", {}, true},
                                  {"gap-tolerance", {}, true},
                                  {"angle-tolerance", {}, true},
                                  {"overlap-tolerance", {}, true},
                                  {"json", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("find-interface requires exactly INFILE and OUTFILE");
    if (!has_opt(p, "region-a") || !has_opt(p, "region-b"))
        throw std::runtime_error("find-interface requires --region-a and --region-b");
    Mesh mesh_a = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    const std::string mesh_b_path = opt_value(p, "mesh-b");
    Mesh mesh_b;
    if (!mesh_b_path.empty())
        mesh_b = read_mesh_cli(mesh_b_path, opt_value(p, "format-b"));
    meshioplusplus::RegionSelector region_a, region_b;
    region_a.mName = opt_value(p, "region-a");
    region_a.mKind = static_cast<std::int32_t>(meshioplusplus::RegionKind::Cell);
    region_b.mName = opt_value(p, "region-b");
    region_b.mKind = static_cast<std::int32_t>(meshioplusplus::RegionKind::Cell);
    meshioplusplus::FindInterfaceOptions options;
    const std::string mode = opt_value(p, "mode", "conforming");
    if (mode == "conforming")
        options.mMode = meshioplusplus::InterfaceMode::Conforming;
    else if (mode == "proximity")
        options.mMode = meshioplusplus::InterfaceMode::Proximity;
    else
        throw std::runtime_error("find-interface: --mode must be conforming or proximity");
    const std::string master = opt_value(p, "master", "a");
    if (master == "a")
        options.mMaster = meshioplusplus::InterfaceMaster::A;
    else if (master == "b")
        options.mMaster = meshioplusplus::InterfaceMaster::B;
    else
        throw std::runtime_error("find-interface: --master must be a or b");
    if (has_opt(p, "gap-tolerance"))
        options.mGapTolerance = meshioplusplus::detail::stod_c(opt_value(p, "gap-tolerance"));
    if (has_opt(p, "angle-tolerance"))
        options.mAngleTolerance = meshioplusplus::detail::stod_c(opt_value(p, "angle-tolerance"));
    if (has_opt(p, "overlap-tolerance"))
        options.mOverlapTolerance = meshioplusplus::detail::stod_c(opt_value(p, "overlap-tolerance"));
    auto result = mesh_b_path.empty()
                      ? meshioplusplus::find_interface(mesh_a, region_a, region_b, options)
                      : meshioplusplus::find_interface(mesh_a, region_a, mesh_b, region_b, options);
    write_mesh_cli(p.positionals[1], result.mMesh, opt_value(p, "output-format"));
    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Field("num_pairs", result.mReport.mNumPairs);
        json.Field("area", result.mReport.mArea);
        json.Field("max_gap", result.mReport.mMaxGap);
        json.Field("unmatched_a", result.mReport.mUnmatchedA);
        json.Field("unmatched_b", result.mReport.mUnmatchedB);
        json.Key("side_a");
        cli_json_side_entries(json, result.mSideA.mEntries);
        json.Key("side_b");
        cli_json_side_entries(json, result.mSideB.mEntries);
        json.EndObject();
    } else {
        char area[40], gap[40];
        meshioplusplus::detail::snprintf_c(area, sizeof(area), "%g", result.mReport.mArea);
        meshioplusplus::detail::snprintf_c(gap, sizeof(gap), "%g", result.mReport.mMaxGap);
        std::cout << "find interface: " << result.mReport.mNumPairs << " facet pair(s), area "
                  << area << ", max gap " << gap << " written\n";
    }
    return 0;
}

void cli_json_xyz(meshioplusplus::cli::JsonOut& rJson, const NDArray& rValues,
                  std::size_t Row) {
    rJson.BeginArray();
    for (std::size_t d = 0; d < 3; ++d)
        rJson.Number(meshioplusplus::detail::read_double(rValues, 3 * Row + d));
    rJson.EndArray();
}

std::string cli_csv_number(double Value) {
    char buffer[48];
    meshioplusplus::detail::snprintf_c(buffer, sizeof(buffer), "%g", Value);
    return buffer;
}

void cli_contact_pairs_json(meshioplusplus::cli::JsonOut& rJson,
                            const meshioplusplus::ContactPairsResult& rResult) {
    const std::size_t count = rResult.mSlavePoint.Shape()[0];
    rJson.BeginObject();
    rJson.Field("num_pairs", count);
    rJson.Key("pairs");
    rJson.BeginArray();
    for (std::size_t i = 0; i < count; ++i) {
        rJson.BeginObject();
        rJson.Field("slave_point", meshioplusplus::detail::read_int(rResult.mSlavePoint, i));
        rJson.Field("master_cell", meshioplusplus::detail::read_int(rResult.mMasterCell, i));
        rJson.Field("master_facet", meshioplusplus::detail::read_int(rResult.mMasterFacet, i));
        rJson.Field("master_subfacet", meshioplusplus::detail::read_int(rResult.mMasterSubfacet, i));
        rJson.Key("local_coordinates");
        cli_json_xyz(rJson, rResult.mLocalCoordinates, i);
        rJson.Key("closest_point");
        cli_json_xyz(rJson, rResult.mClosestPoint, i);
        rJson.Field("gap", meshioplusplus::detail::read_double(rResult.mGap, i));
        rJson.Key("normal");
        cli_json_xyz(rJson, rResult.mNormal, i);
        rJson.EndObject();
    }
    rJson.EndArray();
    rJson.Key("unmatched");
    rJson.BeginArray();
    for (std::size_t i = 0; i < rResult.mUnmatched.Size(); ++i)
        rJson.Int(meshioplusplus::detail::read_int(rResult.mUnmatched, i));
    rJson.EndArray();
    rJson.EndObject();
}

void cli_contact_pairs_csv(std::ostream& rStream,
                           const meshioplusplus::ContactPairsResult& rResult, bool Full) {
    if (Full)
        rStream << "slave_point,master_cell,master_facet,master_subfacet,u,v,w,x,y,z,gap,nx,ny,nz\n";
    else
        rStream << "slave_point,master_cell,master_facet,master_subfacet,gap\n";
    const std::size_t count = rResult.mSlavePoint.Shape()[0];
    for (std::size_t i = 0; i < count; ++i) {
        rStream << std::to_string(meshioplusplus::detail::read_int(rResult.mSlavePoint, i)) << ','
                << std::to_string(meshioplusplus::detail::read_int(rResult.mMasterCell, i)) << ','
                << std::to_string(meshioplusplus::detail::read_int(rResult.mMasterFacet, i)) << ','
                << std::to_string(meshioplusplus::detail::read_int(rResult.mMasterSubfacet, i));
        if (Full) {
            for (const NDArray* array : {&rResult.mLocalCoordinates, &rResult.mClosestPoint})
                for (std::size_t d = 0; d < 3; ++d)
                    rStream << ',' << cli_csv_number(
                        meshioplusplus::detail::read_double(*array, 3 * i + d));
            rStream << ',' << cli_csv_number(
                meshioplusplus::detail::read_double(rResult.mGap, i));
            for (std::size_t d = 0; d < 3; ++d)
                rStream << ',' << cli_csv_number(
                    meshioplusplus::detail::read_double(rResult.mNormal, 3 * i + d));
        } else {
            rStream << ',' << cli_csv_number(
                meshioplusplus::detail::read_double(rResult.mGap, i));
        }
        rStream << '\n';
    }
}

int cmd_contact_pairs(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"slave-region", {}, true},
                                  {"master-region", {}, true},
                                  {"slave-format", {}, true},
                                  {"master-format", {}, true},
                                  {"output", {"-o"}, true},
                                  {"tolerance", {}, true},
                                  {"require-complete", {}, false},
                                  {"json", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("contact-pairs requires exactly SLAVE_FILE and MASTER_FILE");
    if (!has_opt(p, "slave-region") || !has_opt(p, "master-region"))
        throw std::runtime_error("contact-pairs requires --slave-region and --master-region");
    Mesh slave = read_mesh_cli(p.positionals[0], opt_value(p, "slave-format"));
    Mesh master = read_mesh_cli(p.positionals[1], opt_value(p, "master-format"));
    meshioplusplus::RegionSelector slave_region, master_region;
    slave_region.mName = opt_value(p, "slave-region");
    slave_region.mKind = static_cast<std::int32_t>(meshioplusplus::RegionKind::Point);
    master_region.mName = opt_value(p, "master-region");
    master_region.mKind = static_cast<std::int32_t>(meshioplusplus::RegionKind::Cell);
    meshioplusplus::ContactPairsOptions options;
    if (has_opt(p, "tolerance"))
        options.mTolerance = meshioplusplus::detail::stod_c(opt_value(p, "tolerance"));
    options.mRequireComplete = has_flag(p, "require-complete");
    const auto result = meshioplusplus::contact_pairs(slave, slave_region, master,
                                                     master_region, options);
    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        cli_contact_pairs_json(json, result);
        return 0;
    }
    const std::string output = opt_value(p, "output");
    if (output.empty()) {
        cli_contact_pairs_csv(std::cout, result, false);
    } else {
        auto file = meshioplusplus::detail::make_classic_ofstream(output);
        if (!file)
            throw WriteError("meshio++: contact-pairs: cannot open output '" + output + "'");
        cli_contact_pairs_csv(file, result, true);
    }
    return 0;
}

int cmd_split_interface(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"side-region", {}, true},
                                  {"side-entries", {}, true},
                                  {"add-cohesive", {}, false},
                                  {"json", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("split-interface requires exactly INFILE and OUTFILE");
    if (has_opt(p, "side-region") == has_opt(p, "side-entries"))
        throw std::runtime_error("split-interface requires exactly one of --side-region and --side-entries");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    meshioplusplus::RegionSelector side;
    side.mKind = static_cast<std::int32_t>(meshioplusplus::RegionKind::Side);
    if (has_opt(p, "side-entries")) {
        const std::string& text = opt_value(p, "side-entries");
        std::vector<std::int64_t> values;
        std::size_t pos = 0;
        auto skip_space = [&]() {
            while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
        };
        auto expect = [&](char c) {
            skip_space();
            if (pos >= text.size() || text[pos] != c)
                throw std::runtime_error("split-interface: --side-entries must be a JSON list of [cell, facet] pairs");
            ++pos;
        };
        auto integer = [&]() {
            skip_space();
            const std::size_t start = pos;
            if (pos < text.size() && (text[pos] == '-' || text[pos] == '+')) ++pos;
            const std::size_t digits = pos;
            while (pos < text.size() && text[pos] >= '0' && text[pos] <= '9') ++pos;
            if (digits == pos)
                throw std::runtime_error("split-interface: --side-entries contains a non-integer entry");
            return std::stoll(text.substr(start, pos - start));
        };
        expect('[');
        skip_space();
        if (pos < text.size() && text[pos] != ']') {
            while (true) {
                expect('[');
                values.push_back(integer());
                expect(',');
                values.push_back(integer());
                expect(']');
                skip_space();
                if (pos < text.size() && text[pos] == ',') {
                    ++pos;
                    continue;
                }
                break;
            }
        }
        expect(']');
        skip_space();
        if (pos != text.size())
            throw std::runtime_error("split-interface: trailing text in --side-entries");
        Mesh with_side = meshioplusplus::detail::clone_mesh(mesh);
        NDArray entries = NDArray::Uninit(DType::Int64, {values.size() / 2, 2});
        std::copy(values.begin(), values.end(), entries.As<std::int64_t>());
        with_side.AddRegion(meshioplusplus::Region("cli:interface",
            meshioplusplus::RegionKind::Side, std::move(entries)));
        mesh = std::move(with_side);
        side.mName = "cli:interface";
    } else {
        side.mName = opt_value(p, "side-region");
    }
    meshioplusplus::SplitInterfaceOptions options;
    options.mAddCohesive = has_flag(p, "add-cohesive");
    const auto result = meshioplusplus::split_interface(mesh, side, options);
    write_mesh_cli(p.positionals[1], result.mMesh, opt_value(p, "output-format"));
    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Field("num_duplicated_points", result.mNumDuplicatedPoints);
        json.Field("num_cohesive_cells", result.mNumCohesiveCells);
        json.EndObject();
    } else {
        std::cout << "split interface: duplicated " << result.mNumDuplicatedPoints
                  << " point(s), added " << result.mNumCohesiveCells << " cohesive cell(s)\n";
    }
    return 0;
}

int cmd_hausdorff(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format-a", {}, true},
                                  {"input-format-b", {}, true},
                                  {"region-a", {}, true},
                                  {"region-b", {}, true},
                                  {"face-samples", {}, true},
                                  {"max", {}, true},
                                  {"json", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("hausdorff requires exactly INFILE_A and INFILE_B");
    Mesh a = read_mesh_cli(p.positionals[0], opt_value(p, "input-format-a"));
    Mesh b = read_mesh_cli(p.positionals[1], opt_value(p, "input-format-b"));
    meshioplusplus::HausdorffOptions options;
    options.mRegionA = opt_value(p, "region-a");
    options.mRegionB = opt_value(p, "region-b");
    if (has_opt(p, "face-samples"))
        options.mFaceSamples = std::stoll(opt_value(p, "face-samples"));
    const meshioplusplus::HausdorffResult r = meshioplusplus::hausdorff_distance(a, b, options);
    const bool has_max = has_opt(p, "max");
    const double limit = has_max ? meshioplusplus::detail::stod_c(opt_value(p, "max")) : 0.0;
    const bool failed = has_max && r.mDistance > limit;
    // "%.9g" / "%g", the Python CLI's own formats.
    auto g = [](double v) {
        char buf[32];
        meshioplusplus::detail::snprintf_c(buf, sizeof(buf), "%.9g", v);
        return std::string(buf);
    };
    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Field("distance", r.mDistance);
        json.Field("a_to_b", r.mAtoB);
        json.Field("b_to_a", r.mBtoA);
        json.Field("mean_a_to_b", r.mMeanAtoB);
        json.Field("rms_a_to_b", r.mRmsAtoB);
        json.Field("mean_b_to_a", r.mMeanBtoA);
        json.Field("rms_b_to_a", r.mRmsBtoA);
        json.Field("num_samples_a", r.mNumSamplesA);
        json.Field("num_samples_b", r.mNumSamplesB);
        json.Key("worst_point_a");
        json.BeginArray();
        for (double v : r.mWorstPointA)
            json.Number(v);
        json.EndArray();
        json.Key("worst_point_b");
        json.BeginArray();
        for (double v : r.mWorstPointB)
            json.Number(v);
        json.EndArray();
        if (has_max) {
            json.Field("max", limit);
            json.Field("passed", !failed);
        }
        json.EndObject();
    } else {
        std::cout << "Hausdorff distance: " << g(r.mDistance) << "\n";
        std::cout << "  A -> B: max " << g(r.mAtoB) << ", mean " << g(r.mMeanAtoB) << ", rms "
                  << g(r.mRmsAtoB) << " (" << r.mNumSamplesA << " samples)\n";
        std::cout << "  B -> A: max " << g(r.mBtoA) << ", mean " << g(r.mMeanBtoA) << ", rms "
                  << g(r.mRmsBtoA) << " (" << r.mNumSamplesB << " samples)\n";
        if (has_max) {
            char lim[32];
            meshioplusplus::detail::snprintf_c(lim, sizeof(lim), "%g", limit);
            std::cout << "  " << (failed ? "FAIL" : "pass") << ": limit " << lim << "\n";
        }
    }
    return failed ? 1 : 0;
}

// Comma-separated numbers, exactly `Count` of them (any count when 0).
std::vector<double> cli_numbers(const std::string& rText, std::size_t Count, const char* pFlag) {
    std::vector<double> out;
    std::size_t start = 0;
    while (true) {
        const auto comma = rText.find(',', start);
        out.push_back(meshioplusplus::detail::stod_c(rText.substr(start, comma - start)));
        if (comma == std::string::npos)
            break;
        start = comma + 1;
    }
    if (Count != 0 && out.size() != Count)
        throw std::runtime_error(std::string("periodic: ") + pFlag + " expects " +
                                 std::to_string(Count) + " numbers, got " +
                                 std::to_string(out.size()));
    return out;
}

int cmd_periodic(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"slave", {}, true},
                                  {"master", {}, true},
                                  {"translate", {}, true},
                                  {"rotate", {}, true},
                                  {"origin", {}, true},
                                  {"matrix", {}, true},
                                  {"atol", {}, true},
                                  {"allow-incomplete", {}, false},
                                  {"output", {"-o"}, true},
                                  {"json", {}, false},
                              });
    if (p.positionals.size() != 1)
        throw std::runtime_error("periodic requires exactly INFILE");
    if (!has_opt(p, "slave") || !has_opt(p, "master"))
        throw std::runtime_error("periodic requires --slave and --master");
    if (!has_opt(p, "matrix") && !has_opt(p, "translate") && !has_opt(p, "rotate"))
        throw std::runtime_error("periodic: give a transform: --translate, --rotate or --matrix");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::PeriodicOptions options;
    if (has_opt(p, "matrix")) {
        const std::vector<double> m = cli_numbers(opt_value(p, "matrix"), 16, "--matrix");
        options.mTransform = meshioplusplus::transform_from_matrix(m.data());
    } else {
        meshioplusplus::AffineTransform xf;
        if (has_opt(p, "rotate")) {
            const std::string text = opt_value(p, "rotate");
            double ax = 0.0, ay = 0.0, az = 0.0, deg = 0.0;
            const auto comma = text.find(',');
            const std::string head = text.substr(0, comma);
            if (comma != std::string::npos && (head == "x" || head == "y" || head == "z")) {
                (head == "x" ? ax : head == "y" ? ay : az) = 1.0;
                deg = meshioplusplus::detail::stod_c(text.substr(comma + 1));
            } else {
                const std::vector<double> v = cli_numbers(text, 4, "--rotate");
                ax = v[0];
                ay = v[1];
                az = v[2];
                deg = v[3];
            }
            std::vector<double> o = {0.0, 0.0, 0.0};
            if (has_opt(p, "origin"))
                o = cli_numbers(opt_value(p, "origin"), 3, "--origin");
            const double rad = deg * 3.14159265358979323846 / 180.0;
            xf = meshioplusplus::transform_compose(
                meshioplusplus::transform_translation(o[0], o[1], o[2]),
                meshioplusplus::transform_compose(
                    meshioplusplus::transform_rotation(ax, ay, az, rad),
                    meshioplusplus::transform_translation(-o[0], -o[1], -o[2])));
        }
        if (has_opt(p, "translate")) {
            const std::vector<double> t = cli_numbers(opt_value(p, "translate"), 3, "--translate");
            xf = meshioplusplus::transform_compose(
                meshioplusplus::transform_translation(t[0], t[1], t[2]), xf);
        }
        options.mTransform = xf;
    }
    if (has_opt(p, "atol"))
        options.mAtol = meshioplusplus::detail::stod_c(opt_value(p, "atol"));
    options.mRequireComplete = !has_flag(p, "allow-incomplete");
    meshioplusplus::RegionSelector slave;
    slave.mName = opt_value(p, "slave");
    meshioplusplus::RegionSelector master;
    master.mName = opt_value(p, "master");
    const meshioplusplus::PeriodicPairs r =
        meshioplusplus::match_periodic_nodes(mesh, slave, master, options);
    const std::size_t n = r.mSlave.Size();
    const std::int64_t* ps = n ? r.mSlave.As<std::int64_t>() : nullptr;
    const std::int64_t* pm = n ? r.mMaster.As<std::int64_t>() : nullptr;

    if (has_opt(p, "output")) {
        auto out = meshioplusplus::detail::make_classic_ofstream(opt_value(p, "output"));
        if (!out)
            throw std::runtime_error("periodic: cannot open '" + opt_value(p, "output") + "'");
        out << "# slave,master\n";
        for (std::size_t i = 0; i < n; ++i)
            out << ps[i] << "," << pm[i] << "\n";
    }
    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Field("num_pairs", n);
        json.Key("pairs");
        json.BeginArray();
        for (std::size_t i = 0; i < n; ++i) {
            json.BeginArray();
            json.Int(ps[i]);
            json.Int(pm[i]);
            json.EndArray();
        }
        json.EndArray();
        json.Key("unmatched");
        json.BeginArray();
        for (std::size_t i = 0; i < r.mUnmatched.Size(); ++i)
            json.Int(r.mUnmatched.As<std::int64_t>()[i]);
        json.EndArray();
        json.Field("num_fixed", r.mNumFixed);
        json.Field("max_residual", r.mMaxResidual);
        json.EndObject();
        return 0;
    }
    std::cout << "periodic pairs " << slave.mName << " -> " << master.mName << ": " << n << "\n";
    std::cout << "  unmatched slave nodes: " << r.mUnmatched.Size() << "\n";
    std::cout << "  fixed points:          " << r.mNumFixed << "\n";
    char res[32];
    meshioplusplus::detail::snprintf_c(res, sizeof(res), "%.3g", r.mMaxResidual);
    std::cout << "  max residual:          " << res << "\n";
    if (!has_opt(p, "output")) {
        for (std::size_t i = 0; i < n && i < 20; ++i)
            std::cout << "  " << ps[i] << " -> " << pm[i] << "\n";
        if (n > 20)
            std::cout << "  ... (" << n - 20 << " more; use --output or --json)\n";
    }
    return 0;
}

int cmd_repair(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"no-fix-orientation", {}, false},
                                  {"no-orient-outward", {}, false},
                                  {"no-fill-holes", {}, false},
                                  {"no-split-non-manifold", {}, false},
                                  {"max-hole-edges", {}, true},
                                  {"weld-tolerance", {}, true},
                                  {"record-provenance", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("repair requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::RepairOptions options;
    options.mFixOrientation = !has_flag(p, "no-fix-orientation");
    options.mOrientOutward = !has_flag(p, "no-orient-outward");
    options.mFillHoles = !has_flag(p, "no-fill-holes");
    options.mSplitNonManifold = !has_flag(p, "no-split-non-manifold");
    if (has_opt(p, "max-hole-edges"))
        options.mMaxHoleEdges = std::stoll(opt_value(p, "max-hole-edges"));
    if (has_opt(p, "weld-tolerance"))
        options.mWeldTolerance = meshioplusplus::detail::stod_c(opt_value(p, "weld-tolerance"));
    options.mRecordProvenance = has_flag(p, "record-provenance");

    meshioplusplus::RepairResult r = meshioplusplus::repair(mesh, options);

    if (!has_flag(p, "quiet")) {
        auto quality_line = [](const meshioplusplus::SurfaceQuality& q) {
            if (q.mWatertight)
                return std::string("watertight");
            return std::to_string(q.mBoundaryEdges) + " boundary / " +
                   std::to_string(q.mNonManifoldEdges) + " non-manifold / " +
                   std::to_string(q.mInconsistentPairs) + " inconsistent edge(s), " +
                   std::to_string(q.mDegenerateTriangles) + " degenerate triangle(s)";
        };
        std::cout << "repair\n";
        std::cout << "  before:              " << quality_line(r.mQualityBefore) << "\n";
        std::cout << "  after:               " << quality_line(r.mQualityAfter) << "\n";
        std::cout << "  triangles rewound:   " << r.mNumFlipped << "\n";
        std::cout << "  components:          " << r.mNumComponents << " (largest "
                  << r.mLargestComponent << ", " << r.mNumOrientedOutward << " oriented outward, "
                  << r.mNumUnorientable << " unorientable)\n";
        std::cout << "  vertices split:      " << r.mNumVerticesSplit << "\n";
        std::cout << "  holes:               " << r.mNumHolesFilled << " filled, "
                  << r.mNumHolesSkipped << " skipped (of " << r.mNumHolesDetected << " detected)\n";
        std::cout << "  added:               " << r.mNumFacesAdded << " face(s), "
                  << r.mNumPointsAdded << " point(s)\n";
        if (r.mPointsWelded > 0)
            std::cout << "  points welded:       " << r.mPointsWelded << "\n";
    }

    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_shrinkwrap(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"target-format", {}, true},
                                  {"offset", {}, true},
                                  {"max-distance", {}, true},
                                  {"weights", {}, true},
                                  {"target-region", {}, true},
                                  {"normal-weight", {}, true},
                                  {"record-distance", {}, false},
                                  {"record-closest-cell", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 3)
        throw std::runtime_error("shrinkwrap requires exactly INFILE TARGET and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    Mesh target = read_mesh_cli(p.positionals[1], opt_value(p, "target-format"));

    meshioplusplus::ShrinkwrapOptions options;
    if (has_opt(p, "offset"))
        options.mOffset = meshioplusplus::detail::stod_c(opt_value(p, "offset"));
    if (has_opt(p, "max-distance"))
        options.mMaxDistance = meshioplusplus::detail::stod_c(opt_value(p, "max-distance"));
    options.mWeights = opt_value(p, "weights");
    options.mTargetRegion = opt_value(p, "target-region");
    const std::string weight = opt_value(p, "normal-weight");
    options.mNormalWeight = meshioplusplus::sdf_weight_from_name(weight.empty() ? "angle" : weight);
    options.mRecordDistance = has_flag(p, "record-distance");
    options.mRecordClosestCell = has_flag(p, "record-closest-cell");

    meshioplusplus::ShrinkwrapResult r = meshioplusplus::shrinkwrap(mesh, target, options);

    if (!has_flag(p, "quiet")) {
        std::cout << "shrinkwrap\n";
        std::cout << "  points projected:    " << r.mNumProjected << "\n";
        std::cout << "  points missed:       " << r.mNumMissed << "\n";
        std::cout << "  points skipped:      " << r.mNumSkipped << "\n";
        std::cout << "  max displacement:    " << r.mMaxDisplacement << "\n";
        std::cout << "  target: ";
        if (r.mQuality.mWatertight)
            std::cout << "watertight\n";
        else
            std::cout << r.mQuality.mBoundaryEdges << " boundary / " << r.mQuality.mNonManifoldEdges
                      << " non-manifold / " << r.mQuality.mInconsistentPairs
                      << " inconsistent edge(s)\n";
    }

    write_mesh_cli(p.positionals[2], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_sobolev_deform(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"array", {}, true},
                                  {"length-scale", {}, true},
                                  {"fixed-points-array", {}, true},
                                  {"fix-boundary", {}, false},
                                  {"record-filtered", {}, false},
                                  {"max-iterations", {}, true},
                                  {"tolerance", {}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("sobolev-deform requires exactly INFILE and OUTFILE");
    if (!has_opt(p, "array"))
        throw std::runtime_error("sobolev-deform requires --array NAME");
    if (!has_opt(p, "length-scale"))
        throw std::runtime_error("sobolev-deform requires --length-scale L");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::SobolevOptions options;
    options.mArrayName = opt_value(p, "array");
    options.mLengthScale = meshioplusplus::detail::stod_c(opt_value(p, "length-scale"));
    options.mFixedPointsArray = opt_value(p, "fixed-points-array");
    options.mFixBoundary = has_flag(p, "fix-boundary");
    options.mRecordFiltered = has_flag(p, "record-filtered");
    if (has_opt(p, "max-iterations"))
        options.mMaxIterations = std::stoi(opt_value(p, "max-iterations"));
    if (has_opt(p, "tolerance"))
        options.mTolerance = meshioplusplus::detail::stod_c(opt_value(p, "tolerance"));

    meshioplusplus::SobolevResult r = meshioplusplus::sobolev_deform(mesh, options);

    if (!has_flag(p, "quiet")) {
        std::cout << "sobolev_deform\n";
        std::cout << "  iterations:          " << r.mNumIterations
                  << (r.mConverged ? " (converged)" : " (NOT converged)") << "\n";
        std::cout << "  relative residual:   " << r.mResidual << "\n";
        std::cout << "  points fixed:        " << r.mNumFixed << "\n";
        std::cout << "  points isolated:     " << r.mNumIsolated << "\n";
        std::cout << "  max displacement:    " << r.mMaxDisplacement << "\n";
    }

    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
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
            out.push_back(meshioplusplus::detail::stod_c(tok));
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
        xf = meshioplusplus::transform_units(
            meshioplusplus::detail::stod_c(opt_value(p, "scale-units")));
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
        opts.atol = meshioplusplus::detail::stod_c(opt_value(p, "atol"));

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

// `NAME OP VALUE`, the spelling `crop --where` takes.
//
// Scanned longest-operator-first so "<=" is not read as "<" plus a stray "=",
// and from the start of the string, because array names routinely contain a
// colon (`sdf:distance`, `quality:scaled_jacobian`) but never an operator. The
// Python CLI's `_parse_where` implements the identical rule.
struct cli_where {
    std::string mName;
    meshioplusplus::RefineCompare mOp = meshioplusplus::RefineCompare::Less;
    double mValue = 0.0;
};

cli_where cli_parse_where(const std::string& rText, const char* pVerb) {
    static const char* ops[] = {"<=", ">=", "==", "!=", "<", ">"};
    for (const char* op : ops) {
        const std::size_t idx = rText.find(op);
        if (idx == std::string::npos || idx == 0)
            continue;
        std::string name = rText.substr(0, idx);
        std::string value = rText.substr(idx + std::strlen(op));
        const auto trim = [](std::string& s) {
            while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front())))
                s.erase(s.begin());
            while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())))
                s.pop_back();
        };
        trim(name);
        trim(value);
        if (name.empty() || value.empty())
            continue;
        cli_where out;
        out.mName = name;
        out.mOp = meshioplusplus::refine_compare_from_name(op);
        out.mValue = meshioplusplus::detail::stod_c(value);
        return out;
    }
    throw std::runtime_error(std::string(pVerb) +
                             ": --where expects 'NAME OP VALUE', with OP one of "
                             "<=, >=, ==, !=, <, >");
}

int cmd_crop(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"bbox", {}, true},
                                  {"plane", {}, true},
                                  {"where", {}, true},
                                  {"mode", {}, true},
                                  {"record-ids", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("crop requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    const bool has_bbox = p.values.count("bbox") != 0;
    const bool has_plane = p.values.count("plane") != 0;
    const bool has_where = p.values.count("where") != 0;
    if (static_cast<int>(has_bbox) + static_cast<int>(has_plane) + static_cast<int>(has_where) != 1)
        throw std::runtime_error("crop: give exactly one of --bbox, --plane or --where");
    if (has_where && p.values.count("mode"))
        throw std::runtime_error(
            "crop: --mode applies to --bbox and --plane, which test points; a --where "
            "predicate is already one value per cell and has nothing to reduce");
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
    } else if (has_plane) {
        auto v = parse_doubles(opt_value(p, "plane"));
        if (v.size() != 6)
            throw std::runtime_error("crop: --plane expects 'px,py,pz,nx,ny,nz'");
        double point[3] = {v[0], v[1], v[2]}, normal[3] = {v[3], v[4], v[5]};
        r = meshioplusplus::crop_halfspace(mesh, point, normal, mode, record_ids);
    } else {
        const auto pred = cli_parse_where(opt_value(p, "where"), "crop");
        r = meshioplusplus::crop_predicate(mesh, pred.mName, pred.mOp, pred.mValue, record_ids);
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_slice(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"origin", {}, true},
                                  {"normal", {}, true},
                                  {"record-parent-ids", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("slice requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    // Negatives need the --origin=/--normal= form (the parser rule shared with
    // --bbox / --mu).
    auto o = parse_doubles(opt_value(p, "origin", "0,0,0"));
    auto n = parse_doubles(opt_value(p, "normal", "0,0,1"));
    if (o.size() != 3)
        throw std::runtime_error("slice: --origin expects 'x,y,z'");
    if (n.size() != 3)
        throw std::runtime_error("slice: --normal expects 'x,y,z'");

    meshioplusplus::SliceOptions options;
    options.mOrigin = {o[0], o[1], o[2]};
    options.mNormal = {n[0], n[1], n[2]};
    options.mRecordParentIds = has_flag(p, "record-parent-ids");
    Mesh out = meshioplusplus::slice(mesh, options);

    write_mesh_cli(p.positionals[1], out, opt_value(p, "output-format"));
    return 0;
}

int cmd_voxelize(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"resolution", {}, true},
                                  {"cell-size", {}, true},
                                  {"bounds", {}, true},
                                  {"padding", {}, true},
                                  {"padding-relative", {}, true},
                                  {"fill", {}, true},
                                  {"sign", {}, true},
                                  {"attach-occupancy", {}, false},
                                  {"max-cells", {}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("voxelize requires exactly INFILE and OUTFILE");

    const bool has_res = p.values.count("resolution") != 0;
    const bool has_cell = p.values.count("cell-size") != 0;
    if (has_res == has_cell)
        throw std::runtime_error("voxelize: give exactly one of --resolution and --cell-size");

    meshioplusplus::VoxelOptions options;
    if (has_res) {
        // parse_doubles rather than the int64 parser, which is defined further
        // down this file; the values are small counts either way.
        auto v = parse_doubles(opt_value(p, "resolution"));
        if (v.size() != 3)
            throw std::runtime_error("voxelize: --resolution expects 'nx,ny,nz'");
        options.mResolution = std::array<std::int64_t, 3>{{static_cast<std::int64_t>(v[0]),
                                                           static_cast<std::int64_t>(v[1]),
                                                           static_cast<std::int64_t>(v[2])}};
    } else {
        options.mCellSize = meshioplusplus::detail::stod_c(opt_value(p, "cell-size"));
    }
    // Negatives need the --bounds= form (the parser rule shared with --bbox).
    if (p.values.count("bounds")) {
        auto v = parse_doubles(opt_value(p, "bounds"));
        if (v.size() != 6)
            throw std::runtime_error("voxelize: --bounds expects 'xlo,ylo,zlo,xhi,yhi,zhi'");
        options.mBounds = std::array<double, 6>{{v[0], v[1], v[2], v[3], v[4], v[5]}};
    }
    if (p.values.count("padding"))
        options.mPadding = meshioplusplus::detail::stod_c(opt_value(p, "padding"));
    if (p.values.count("padding-relative"))
        options.mPaddingRelative = meshioplusplus::detail::stod_c(opt_value(p, "padding-relative"));
    options.mFill = meshioplusplus::voxel_fill_from_name(opt_value(p, "fill", "all"));
    options.mDistance.mSign =
        meshioplusplus::sdf_sign_from_name(opt_value(p, "sign", "pseudonormal"));
    options.mAttachOccupancy = has_flag(p, "attach-occupancy");
    if (p.values.count("max-cells"))
        options.mMaxCells = std::stoll(opt_value(p, "max-cells"));
    // Only the inside fill depends on the surface being closed, so only it warns.
    options.mDistance.mWatertightCheck = options.mFill == meshioplusplus::VoxelFill::Inside
                                             ? meshioplusplus::SdfWatertightCheck::Warn
                                             : meshioplusplus::SdfWatertightCheck::Off;

    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    meshioplusplus::VoxelResult r = meshioplusplus::voxelize(mesh, options);

    if (!has_flag(p, "quiet")) {
        std::cout << "voxelized (" << opt_value(p, "fill", "all") << ")\n";
        std::cout << "  grid:           " << r.mDims[0] << " x " << r.mDims[1] << " x "
                  << r.mDims[2] << "\n";
        std::cout << "  cell size:      " << r.mSpacing[0] << ", " << r.mSpacing[1] << ", "
                  << r.mSpacing[2] << "\n";
        std::cout << "  cells kept:     " << r.mNumOccupied << "\n";
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_sdf(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"structure", {}, true},
                                  {"resolution", {}, true},
                                  {"cell-size", {}, true},
                                  {"bounds", {}, true},
                                  {"padding", {}, true},
                                  {"padding-relative", {}, true},
                                  {"root-resolution", {}, true},
                                  {"max-depth", {}, true},
                                  {"band-cells", {}, true},
                                  {"sign", {}, true},
                                  {"location", {}, true},
                                  {"band", {}, true},
                                  {"watertight-check", {}, true},
                                  {"max-cells", {}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("sdf requires exactly INFILE and OUTFILE");

    meshioplusplus::SdfOptions options;
    const std::string structure = opt_value(p, "structure", "voxel");
    options.mStructure = meshioplusplus::sdf_structure_from_name(structure);
    // resolution/cell_size size a voxel grid; the octree's finest cell is
    // root/2^depth and is therefore already determined. The core refuses the
    // combination by name, so this verb just forwards whatever was given.
    if (p.values.count("resolution")) {
        auto v = parse_doubles(opt_value(p, "resolution"));
        if (v.size() != 3)
            throw std::runtime_error("sdf: --resolution expects 'nx,ny,nz'");
        options.mResolution = std::array<std::int64_t, 3>{{static_cast<std::int64_t>(v[0]),
                                                           static_cast<std::int64_t>(v[1]),
                                                           static_cast<std::int64_t>(v[2])}};
    }
    if (p.values.count("cell-size"))
        options.mCellSize = meshioplusplus::detail::stod_c(opt_value(p, "cell-size"));
    // Negatives need the --bounds= form (the parser rule shared with --bbox).
    if (p.values.count("bounds")) {
        auto v = parse_doubles(opt_value(p, "bounds"));
        if (v.size() != 6)
            throw std::runtime_error("sdf: --bounds expects 'xlo,ylo,zlo,xhi,yhi,zhi'");
        options.mBounds = std::array<double, 6>{{v[0], v[1], v[2], v[3], v[4], v[5]}};
    }
    if (p.values.count("padding"))
        options.mPadding = meshioplusplus::detail::stod_c(opt_value(p, "padding"));
    if (p.values.count("padding-relative"))
        options.mPaddingRelative = meshioplusplus::detail::stod_c(opt_value(p, "padding-relative"));
    if (p.values.count("root-resolution"))
        options.mRootResolution = std::stoll(opt_value(p, "root-resolution"));
    if (p.values.count("max-depth"))
        options.mMaxDepth = std::stoll(opt_value(p, "max-depth"));
    if (p.values.count("band-cells"))
        options.mBandCells = meshioplusplus::detail::stod_c(opt_value(p, "band-cells"));
    if (p.values.count("max-cells"))
        options.mMaxCells = std::stoll(opt_value(p, "max-cells"));
    options.mDistance.mSign =
        meshioplusplus::sdf_sign_from_name(opt_value(p, "sign", "pseudonormal"));
    options.mDistance.mLocation =
        meshioplusplus::sdf_location_from_name(opt_value(p, "location", "corner"));
    if (p.values.count("band"))
        options.mDistance.mBand = meshioplusplus::detail::stod_c(opt_value(p, "band"));
    options.mDistance.mWatertightCheck =
        meshioplusplus::sdf_watertight_check_from_name(opt_value(p, "watertight-check", "warn"));

    Mesh surface = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    meshioplusplus::SdfResult r = meshioplusplus::compute_sdf(surface, options);

    if (!has_flag(p, "quiet")) {
        std::size_t cells = 0;
        for (const auto cb : r.mMesh.CellRange())
            cells += cb.NumCells();
        std::cout << "signed distance field (" << structure << ")\n";
        std::cout << "  root grid:      " << r.mDims[0] << " x " << r.mDims[1] << " x "
                  << r.mDims[2] << "\n";
        std::cout << "  finest cell:    " << r.mSpacing[0] << ", " << r.mSpacing[1] << ", "
                  << r.mSpacing[2] << "\n";
        if (r.mMaxDepth != 0)
            std::cout << "  octree depth:   " << r.mMaxDepth << "\n";
        std::cout << "  cells:          " << cells << "\n";
        if (options.mDistance.mBand > 0.0)
            std::cout << "  banded:         " << r.mNumBanded << "\n";
        if (!r.mQuality.mWatertight)
            std::cout << "  surface:        NOT watertight (" << r.mQuality.mBoundaryEdges
                      << " boundary, " << r.mQuality.mNonManifoldEdges << " non-manifold, "
                      << r.mQuality.mInconsistentPairs << " inconsistent, "
                      << r.mQuality.mDegenerateTriangles << " degenerate)\n";
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_isosurface(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"array", {}, true},
                                  {"values", {}, true},
                                  {"component", {}, true},
                                  {"record-parent-ids", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("isosurface requires exactly INFILE and OUTFILE");
    const std::string array = opt_value(p, "array");
    if (array.empty())
        throw std::runtime_error("isosurface: --array NAME is required");
    const std::string values = opt_value(p, "values");
    if (values.empty())
        throw std::runtime_error("isosurface: --values 'v1,v2,...' is required");

    // Negatives need the --values= form (the parser rule shared with --bbox /
    // --origin / --mu).
    auto vals = parse_doubles(values);
    if (vals.empty())
        throw std::runtime_error("isosurface: --values expects at least one isovalue");

    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::IsosurfaceOptions options;
    options.mArrayName = array;
    options.mIsovalues = std::move(vals);
    const std::string comp = opt_value(p, "component");
    if (!comp.empty()) {
        const int c = std::stoi(comp);
        if (c >= 0)
            options.mComponent = c;
    }
    options.mRecordParentIds = has_flag(p, "record-parent-ids");
    Mesh out = meshioplusplus::isosurface(mesh, options);

    write_mesh_cli(p.positionals[1], out, opt_value(p, "output-format"));
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

int cmd_convert_cells(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"mode", {}, true},
                                  {"record-parent-ids", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("convert-cells requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::ConvertCellsOptions options;
    options.mMode = meshioplusplus::convert_cells_mode_from_name(opt_value(p, "mode", "linearize"));
    options.mRecordParentIds = has_flag(p, "record-parent-ids");

    auto result = meshioplusplus::convert_cells(mesh, options);
    write_mesh_cli(p.positionals[1], result.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_subdivide(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"record-parent-ids", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("subdivide requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::SubdivideOptions options;
    options.mRecordParentIds = has_flag(p, "record-parent-ids");

    auto result = meshioplusplus::subdivide(mesh, options);
    write_mesh_cli(p.positionals[1], result.mMesh, opt_value(p, "output-format"));
    return 0;
}

std::int64_t total_cells(const Mesh& rMesh) {
    std::int64_t n = 0;
    for (const auto cb : rMesh.CellRange())
        n += static_cast<std::int64_t>(cb.NumCells());
    return n;
}

int cmd_agglomerate(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"target-group-size", {}, true},
                                  {"merge-coplanar-faces", {}, false},
                                  {"coplanar-angle", {}, true},
                                  {"min-sphericity", {}, true},
                                  {"json", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("agglomerate requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::AgglomerateOptions options;
    options.mTargetGroupSize =
        static_cast<std::size_t>(std::stoull(opt_value(p, "target-group-size", "8")));
    options.mMergeCoplanarFaces = has_flag(p, "merge-coplanar-faces");
    if (has_opt(p, "coplanar-angle"))
        options.mCoplanarAngleDeg = meshioplusplus::detail::stod_c(opt_value(p, "coplanar-angle"));
    if (has_opt(p, "min-sphericity"))
        options.mMinSphericity = meshioplusplus::detail::stod_c(opt_value(p, "min-sphericity"));

    auto result = meshioplusplus::agglomerate(mesh, options);
    write_mesh_cli(p.positionals[1], result.mMesh, opt_value(p, "output-format"));
    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Field("cells_in", total_cells(mesh));
        json.Field("cells_out", total_cells(result.mMesh));
        json.Field("num_faces_merged", result.mNumFacesMerged);
        json.Field("num_rejected", result.mNumRejected);
        json.EndObject();
    }
    return 0;
}

int cmd_resample(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input", {}, true},
                                  {"input-format", {"-i"}, true},
                                  {"time-from", {}, true},
                                  {"times", {}, true},
                                  {"times-from", {}, true},
                                  {"method", {}, true},
                                  {"clamp", {}, false},
                                  {"blend-points", {}, false},
                                  {"json", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("resample requires exactly INFILE and OUTFILE");
    if (has_opt(p, "times") == has_opt(p, "times-from"))
        throw std::runtime_error("resample needs exactly one of --times and --times-from");
    meshioplusplus::SequencePipeline pipeline;
    meshioplusplus::SequenceInput& in = pipeline.mInput;
    const std::string& infile = p.positionals[0];
    const std::vector<std::string>& extra_inputs = opt_values(p, "input");
    if (!extra_inputs.empty()) {
        in.mPaths.push_back(infile);
        for (const std::string& extra : extra_inputs)
            in.mPaths.push_back(extra);
    } else if (infile.find('*') != std::string::npos || infile.find('?') != std::string::npos) {
        in.mPattern = infile;
    } else {
        in.mPaths.push_back(infile);
    }
    in.mFormat = opt_value(p, "input-format");
    in.mTimeFrom = meshioplusplus::sequence_time_from_name(opt_value(p, "time-from"));
    pipeline.mOutput.mPath = p.positionals[1];

    // The targets: START:STOP:STEP or a comma list, or another sequence's times.
    meshioplusplus::SequenceResample rs;
    if (has_opt(p, "times")) {
        const std::string spec = opt_value(p, "times");
        if (std::count(spec.begin(), spec.end(), ':') == 2) {
            const std::size_t a = spec.find(':');
            const std::size_t b = spec.find(':', a + 1);
            rs.mTimes = meshioplusplus::resample_times_range(
                meshioplusplus::detail::stod_c(spec.substr(0, a)),
                meshioplusplus::detail::stod_c(spec.substr(a + 1, b - a - 1)),
                meshioplusplus::detail::stod_c(spec.substr(b + 1)));
        } else {
            for (const std::string& t : data_split_names(spec))
                rs.mTimes.push_back(meshioplusplus::detail::stod_c(t));
        }
    } else {
        meshioplusplus::SequenceInput other;
        other.mPattern = opt_value(p, "times-from");
        other.mTimeFrom = in.mTimeFrom;  // as the pipeline's TimesFrom resolves it
        rs.mTimes = meshioplusplus::sequence_times(other);
    }
    const std::string method = opt_value(p, "method", "linear");
    rs.mMethod = meshioplusplus::resample_method_from_name(method);
    if (has_flag(p, "clamp"))
        rs.mExtrapolate = meshioplusplus::ResampleExtrapolate::Clamp;
    rs.mBlendPoints = has_flag(p, "blend-points");
    const std::vector<double> targets = rs.mTimes;
    pipeline.mResample = rs;
    meshioplusplus::run_sequence_pipeline(pipeline);

    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Field("method", method);
        json.Field("num_steps", targets.size());
        json.Key("times");
        json.BeginArray();
        for (double t : targets)
            json.Number(t);
        json.EndArray();
        json.EndObject();
    } else {
        std::cout << "resample (" << method << "): " << targets.size() << " step(s) -> "
                  << pipeline.mOutput.mPath << "\n";
    }
    return 0;
}

// The comma-separated Int64 list `--cells` takes; parse_doubles' twin.
std::vector<std::int64_t> refine_parse_int64s(const std::string& rText) {
    std::vector<std::int64_t> out;
    std::size_t start = 0;
    while (start <= rText.size()) {
        const std::size_t comma = rText.find(',', start);
        const std::string tok =
            rText.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        if (!tok.empty())
            out.push_back(static_cast<std::int64_t>(std::stoll(tok)));
        if (comma == std::string::npos)
            break;
        start = comma + 1;
    }
    return out;
}

std::string refine_trim(const std::string& rText) {
    const std::size_t first = rText.find_first_not_of(" \t");
    if (first == std::string::npos)
        return "";
    return rText.substr(first, rText.find_last_not_of(" \t") - first + 1);
}

int cmd_refine(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"levels", {}, true},
                                  {"record-parent-ids", {}, false},
                                  {"cells", {}, true},
                                  {"region", {}, true},
                                  {"where", {}, true},
                                  {"closure", {}, true},
                                  {"record-levels", {}, false},
                                  {"record-hierarchy", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("refine requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::RefineOptions options;
    options.mLevels = std::stoi(opt_value(p, "levels", "1"));
    options.mRecordParentIds = has_flag(p, "record-parent-ids");
    options.mRecordLevels = has_flag(p, "record-levels");
    options.mRecordHierarchy = has_flag(p, "record-hierarchy");
    options.mClosure = meshioplusplus::refine_closure_from_name(opt_value(p, "closure"));
    if (has_opt(p, "cells"))
        options.mCells = refine_parse_int64s(opt_value(p, "cells"));
    options.mRegion = opt_value(p, "region");
    if (has_opt(p, "where")) {
        // The operator is the maximal run of <>=! characters, located from the
        // first of them: array names routinely contain ':' and '.' but never a
        // comparison character. Deliberately not a second data_calc grammar.
        const std::string where = opt_value(p, "where");
        const std::size_t start = where.find_first_of("<>=!");
        std::size_t end = start;
        while (end != std::string::npos && end < where.size() &&
               std::string("<>=!").find(where[end]) != std::string::npos)
            ++end;
        if (start == std::string::npos || start == 0 || end >= where.size())
            throw std::runtime_error("refine: cannot parse --where '" + where +
                                     "' (expected 'NAME OP VALUE', e.g. "
                                     "'quality:scaled_jacobian < 0.3')");
        options.mPredicateArray = refine_trim(where.substr(0, start));
        options.mPredicateOp =
            meshioplusplus::refine_compare_from_name(where.substr(start, end - start));
        options.mPredicateValue = meshioplusplus::detail::stod_c(refine_trim(where.substr(end)));
    }

    auto result = meshioplusplus::refine(mesh, options);
    write_mesh_cli(p.positionals[1], result.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_undo_green(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 3)
        throw std::runtime_error("undo-green requires exactly COARSE FINE and OUTFILE");
    Mesh coarse = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    Mesh fine = read_mesh_cli(p.positionals[1], opt_value(p, "input-format"));

    auto result = meshioplusplus::undo_green(coarse, fine);
    if (!has_flag(p, "quiet"))
        std::cout << "undid " << result.mNumGroupsUndone << " green group(s), removing "
                  << result.mNumCellsRemoved << " cell(s)\n";
    write_mesh_cli(p.positionals[2], result.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_decimate(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"ratio", {}, true},
                                  {"target-faces", {}, true},
                                  {"max-error", {}, true},
                                  {"placement", {}, true},
                                  {"feature-angle", {}, true},
                                  {"no-preserve-boundary", {}, false},
                                  {"no-preserve-features", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("decimate requires exactly INFILE and OUTFILE");
    const int num_set = (has_opt(p, "ratio") ? 1 : 0) + (has_opt(p, "target-faces") ? 1 : 0) +
                        (has_opt(p, "max-error") ? 1 : 0);
    if (num_set != 1)
        throw std::runtime_error(
            "decimate: give exactly one of --ratio, --target-faces or --max-error");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::DecimateOptions options;
    // Negative values are the "unset" sentinels the operation validates.
    options.mTargetRatio = meshioplusplus::detail::stod_c(opt_value(p, "ratio", "-1"));
    options.mTargetFaces = std::stoll(opt_value(p, "target-faces", "-1"));
    options.mMaxError = meshioplusplus::detail::stod_c(opt_value(p, "max-error", "-1"));
    options.mPlacement =
        meshioplusplus::decimate_placement_from_name(opt_value(p, "placement", "optimal"));
    options.mFeatureAngleDeg = meshioplusplus::detail::stod_c(opt_value(p, "feature-angle", "30"));
    options.mPreserveBoundary = !has_flag(p, "no-preserve-boundary");
    options.mPreserveFeatures = !has_flag(p, "no-preserve-features");

    auto r = meshioplusplus::decimate(mesh, options);
    if (!has_flag(p, "quiet")) {
        std::size_t faces_out = 0;
        for (const auto cb : r.mMesh.CellRange())
            faces_out += cb.NumCells();
        std::cout << "decimated to " << faces_out << " faces\n";
        std::cout << "  faces removed:            " << r.mFacesRemoved << "\n";
        std::cout << "  points removed:           " << r.mPointsRemoved << "\n";
        std::cout << "  collapses rejected:       " << r.mCollapsesRejected << "\n";
        std::cout << "  max error applied:        " << r.mMaxErrorApplied << "\n";
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_decimate_volume(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"ratio", {}, true},
                                  {"target-cells", {}, true},
                                  {"max-error", {}, true},
                                  {"placement", {}, true},
                                  {"feature-angle", {}, true},
                                  {"preserve-boundary", {}, false},
                                  {"no-preserve-features", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("decimate-volume requires exactly INFILE and OUTFILE");
    const int num_set = (has_opt(p, "ratio") ? 1 : 0) + (has_opt(p, "target-cells") ? 1 : 0) +
                        (has_opt(p, "max-error") ? 1 : 0);
    if (num_set != 1)
        throw std::runtime_error(
            "decimate-volume: give exactly one of --ratio, --target-cells or --max-error");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::DecimateVolumeOptions options;
    options.mTargetRatio = meshioplusplus::detail::stod_c(opt_value(p, "ratio", "-1"));
    options.mTargetCells = std::stoll(opt_value(p, "target-cells", "-1"));
    options.mMaxError = meshioplusplus::detail::stod_c(opt_value(p, "max-error", "-1"));
    options.mPlacement =
        meshioplusplus::decimate_placement_from_name(opt_value(p, "placement", "optimal"));
    options.mFeatureAngleDeg = meshioplusplus::detail::stod_c(opt_value(p, "feature-angle", "30"));
    // Unlike decimate, boundary vertices participate by default here.
    options.mPreserveBoundary = has_flag(p, "preserve-boundary");
    options.mPreserveFeatures = !has_flag(p, "no-preserve-features");

    auto r = meshioplusplus::decimate_volume(mesh, options);
    if (!has_flag(p, "quiet")) {
        std::size_t tets_out = 0;
        for (const auto cb : r.mMesh.CellRange())
            tets_out += cb.NumCells();
        std::cout << "decimated to " << tets_out << " tets\n";
        std::cout << "  tets removed:             " << r.mTetsRemoved << "\n";
        std::cout << "  points removed:           " << r.mPointsRemoved << "\n";
        std::cout << "  collapses rejected:       " << r.mCollapsesRejected << "\n";
        std::cout << "  max error applied:        " << r.mMaxErrorApplied << "\n";
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_remesh(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"num-clusters", {}, true},
                                  {"subdivide", {}, true},
                                  {"subsample-ratio", {}, true},
                                  {"max-subdivide", {}, true},
                                  {"iterations", {}, true},
                                  {"repair-passes", {}, true},
                                  {"metric", {}, true},
                                  {"gradation", {}, true},
                                  {"no-preserve-boundary", {}, false},
                                  {"max-anisotropy", {}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("remesh requires exactly INFILE and OUTFILE");
    if (!has_opt(p, "num-clusters"))
        throw std::runtime_error("remesh: --num-clusters is required");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::RemeshOptions options;
    options.mNumClusters = std::stoll(opt_value(p, "num-clusters"));
    options.mSubdivide = std::stoi(opt_value(p, "subdivide", "-1"));
    options.mSubsampleRatio = meshioplusplus::detail::stod_c(opt_value(p, "subsample-ratio", "10"));
    options.mMaxSubdivide = std::stoi(opt_value(p, "max-subdivide", "4"));
    options.mMaxIterations = std::stoi(opt_value(p, "iterations", "100"));
    options.mMaxRepairPasses = std::stoi(opt_value(p, "repair-passes", "10"));
    options.mMetric = meshioplusplus::remesh_metric_from_name(opt_value(p, "metric", "isotropic"));
    options.mGradation = meshioplusplus::detail::stod_c(opt_value(p, "gradation", "0"));
    options.mPreserveBoundary = !has_flag(p, "no-preserve-boundary");
    options.mMaxAnisotropy = meshioplusplus::detail::stod_c(opt_value(
        p, "max-anisotropy", std::to_string(meshioplusplus::kRemeshDefaultMaxAnisotropy)));

    auto r = meshioplusplus::remesh(mesh, options);
    if (!has_flag(p, "quiet")) {
        std::cout << "remeshed to " << r.mNumClusters << " clusters\n";
        std::cout << "  iterations:               " << r.mNumIterations << "\n";
        std::cout << "  subdivide passes applied: " << r.mSubdivideApplied << "\n";
        std::cout << "  isolated clusters:        " << r.mNumIsolatedClusters << "\n";
        std::cout << "  non-manifold vertices:    " << r.mNumNonManifoldVertices << "\n";
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_remesh_volume(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"resolution", {}, true},
                                  {"cell-size", {}, true},
                                  {"bounds", {}, true},
                                  {"padding", {}, true},
                                  {"padding-relative", {}, true},
                                  {"max-cells", {}, true},
                                  {"max-tets", {}, true},
                                  {"warp-fraction", {}, true},
                                  {"sign", {}, true},
                                  {"watertight-check", {}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("remesh-volume requires exactly INFILE and OUTFILE");

    const bool has_res = p.values.count("resolution") != 0;
    const bool has_cell = p.values.count("cell-size") != 0;
    if (has_res == has_cell)
        throw std::runtime_error("remesh-volume: give exactly one of --resolution and --cell-size");

    meshioplusplus::RemeshVolumeOptions options;
    if (has_res) {
        auto v = parse_doubles(opt_value(p, "resolution"));
        if (v.size() != 3)
            throw std::runtime_error("remesh-volume: --resolution expects 'nx,ny,nz'");
        options.mResolution = std::array<std::int64_t, 3>{{static_cast<std::int64_t>(v[0]),
                                                           static_cast<std::int64_t>(v[1]),
                                                           static_cast<std::int64_t>(v[2])}};
    } else {
        options.mCellSize = meshioplusplus::detail::stod_c(opt_value(p, "cell-size"));
    }
    // Negatives need the --bounds= form, voxelize's own parser rule.
    if (p.values.count("bounds")) {
        auto v = parse_doubles(opt_value(p, "bounds"));
        if (v.size() != 6)
            throw std::runtime_error("remesh-volume: --bounds expects 'xlo,ylo,zlo,xhi,yhi,zhi'");
        options.mBounds = std::array<double, 6>{{v[0], v[1], v[2], v[3], v[4], v[5]}};
    }
    if (p.values.count("padding"))
        options.mPadding = meshioplusplus::detail::stod_c(opt_value(p, "padding"));
    if (p.values.count("padding-relative"))
        options.mPaddingRelative = meshioplusplus::detail::stod_c(opt_value(p, "padding-relative"));
    if (p.values.count("max-cells"))
        options.mMaxCells = std::stoll(opt_value(p, "max-cells"));
    if (p.values.count("max-tets"))
        options.mMaxTets = std::stoll(opt_value(p, "max-tets"));
    // A negative value is rejected by the operation itself, naming the
    // option; --warp-fraction=-0.1 form is needed for a caller to trigger it.
    options.mWarpFraction = meshioplusplus::detail::stod_c(opt_value(p, "warp-fraction", "0.35"));
    options.mDistance.mSign =
        meshioplusplus::sdf_sign_from_name(opt_value(p, "sign", "pseudonormal"));
    options.mDistance.mWatertightCheck =
        meshioplusplus::sdf_watertight_check_from_name(opt_value(p, "watertight-check", "warn"));

    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    meshioplusplus::RemeshVolumeResult r = meshioplusplus::remesh_volume(mesh, options);

    if (!has_flag(p, "quiet")) {
        std::cout << "retetrahedralized\n";
        std::cout << "  tets emitted:             " << r.mNumTets << "\n";
        std::cout << "  vertices warped:          " << r.mNumVerticesWarped << "\n";
        std::cout << "  tets rejected (degen.):   " << r.mNumTetsRejected << "\n";
        std::cout << "  non-manifold edges:       " << r.mNumNonManifoldEdges << "\n";
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_optimize_volume(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"max-iterations", {}, true},
                                  {"no-relocate", {}, false},
                                  {"no-flip", {}, false},
                                  {"no-preserve-boundary", {}, false},
                                  {"min-improvement", {}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("optimize-volume requires exactly INFILE and OUTFILE");

    meshioplusplus::OptimizeVolumeOptions options;
    if (p.values.count("max-iterations"))
        options.mMaxIterations = static_cast<int>(std::stoll(opt_value(p, "max-iterations")));
    options.mRelocate = !has_flag(p, "no-relocate");
    options.mFlip = !has_flag(p, "no-flip");
    options.mPreserveBoundary = !has_flag(p, "no-preserve-boundary");
    // A negative value is rejected/ignored by the operation; --min-improvement=
    // form is needed for a caller to pass one.
    options.mMinImprovement =
        meshioplusplus::detail::stod_c(opt_value(p, "min-improvement", "0.000001"));

    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    meshioplusplus::OptimizeVolumeResult r = meshioplusplus::optimize_volume(mesh, options);

    if (!has_flag(p, "quiet")) {
        std::cout << "ODT-remeshed\n";
        std::cout << "  flips (2-3 / 3-2):        " << r.mNumFlips << " (" << r.mNum23Flips << " / "
                  << r.mNum32Flips << ")\n";
        std::cout << "  vertices moved:           " << r.mNumVerticesMoved << "\n";
        std::cout << "  tets:                     " << r.mNumTets << "\n";
        std::cout << "  min quality before/after: " << r.mMinQualityBefore << " / "
                  << r.mMinQualityAfter << "\n";
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_smooth(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"method", {}, true},
                                  {"iterations", {}, true},
                                  {"lambda", {}, true},
                                  {"mu", {}, true},
                                  {"feature-angle", {}, true},
                                  {"no-fix-boundary", {}, false},
                                  {"no-preserve-features", {}, false},
                                  {"no-guard-inversion", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("smooth requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::SmoothOptions options;
    options.mMethod = meshioplusplus::smooth_method_from_name(opt_value(p, "method", "taubin"));
    options.mIterations = std::stoi(opt_value(p, "iterations", "10"));
    // Negative lambda is the sentinel for "this method's own default"
    // (0.5 Laplacian / 0.33 Taubin), so it is what we pass when unset.
    options.mLambda = meshioplusplus::detail::stod_c(opt_value(p, "lambda", "-1"));
    options.mMu = meshioplusplus::detail::stod_c(opt_value(p, "mu", "-0.34"));
    options.mFeatureAngleDeg = meshioplusplus::detail::stod_c(opt_value(p, "feature-angle", "30"));
    options.mFixBoundary = !has_flag(p, "no-fix-boundary");
    options.mPreserveFeatures = !has_flag(p, "no-preserve-features");
    options.mGuardInversion = !has_flag(p, "no-guard-inversion");

    auto r = meshioplusplus::smooth(mesh, options);
    if (!has_flag(p, "quiet")) {
        std::cout << "smoothed mesh\n";
        std::cout << "  nodes moved:              " << r.mNumNodesMoved << "\n";
        std::cout << "  max displacement:         " << r.mMaxDisplacement << "\n";
        std::cout << "  skipped (inversion):      " << r.mNumSkippedInversion << "\n";
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_interpolate(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"method", {}, true},
                                  {"arrays", {}, true},
                                  {"extrapolate", {}, false},
                                  {"default-value", {}, true},
                                  {"on-conflict", {}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 3)
        throw std::runtime_error("interpolate requires exactly SOURCE TARGET and OUTFILE");
    Mesh source = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    Mesh target = read_mesh_cli(p.positionals[1], opt_value(p, "input-format"));

    meshioplusplus::InterpolateOptions options;
    options.mMethod =
        meshioplusplus::interpolate_method_from_name(opt_value(p, "method", "nearest"));
    if (has_opt(p, "arrays"))
        options.mArrays = data_split_names(opt_value(p, "arrays"));
    options.mExtrapolate = has_flag(p, "extrapolate");
    // A negative value needs the --default-value=-1 form (the parser rule
    // shared with --mu and --bbox).
    options.mDefaultValue = meshioplusplus::detail::stod_c(opt_value(p, "default-value", "0"));
    options.mOnConflict =
        meshioplusplus::interpolate_conflict_from_name(opt_value(p, "on-conflict", "error"));

    Mesh out = meshioplusplus::interpolate(source, target, options);
    if (!has_flag(p, "quiet")) {
        const std::size_t n =
            options.mArrays.empty() ? source.PointDataNames().size() : options.mArrays.size();
        std::cout << "interpolated " << n << " array(s) onto " << out.NumPoints()
                  << " target points\n";
    }
    write_mesh_cli(p.positionals[2], out, opt_value(p, "output-format"));
    return 0;
}

int cmd_conservative_interpolate(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"arrays", {}, true},
                                  {"default-value", {}, true},
                                  {"on-conflict", {}, true},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 3)
        throw std::runtime_error(
            "conservative-interpolate requires exactly SOURCE TARGET and OUTFILE");
    Mesh source = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    Mesh target = read_mesh_cli(p.positionals[1], opt_value(p, "input-format"));

    meshioplusplus::ConservativeInterpolateOptions options;
    if (has_opt(p, "arrays"))
        options.mArrays = data_split_names(opt_value(p, "arrays"));
    // A negative value needs the --default-value=-1 form (the parser rule
    // shared with --mu and --bbox).
    options.mDefaultValue = meshioplusplus::detail::stod_c(opt_value(p, "default-value", "0"));
    options.mOnConflict = meshioplusplus::conservative_interpolate_conflict_from_name(
        opt_value(p, "on-conflict", "error"));

    Mesh out = meshioplusplus::conservative_interpolate(source, target, options);
    if (!has_flag(p, "quiet")) {
        const std::size_t n = options.mArrays.empty()
                                  ? source.PointDataNames().size() + source.CellDataNames().size()
                                  : options.mArrays.size();
        std::cout << "conservatively interpolated " << n << " array(s) onto " << out.NumPoints()
                  << " target points\n";
    }
    write_mesh_cli(p.positionals[2], out, opt_value(p, "output-format"));
    return 0;
}

// The {part} analogue of replace_key (kept separate so split stays untouched).
std::string partition_replace_part(const std::string& rPattern, int part) {
    const std::string token = "{part}";
    const std::string value = std::to_string(part);
    std::string out = rPattern;
    std::size_t pos = out.find(token);
    while (pos != std::string::npos) {
        out.replace(pos, token.size(), value);
        pos = out.find(token, pos + value.size());
    }
    return out;
}

/// `"pvtu"` / `"pvtp"` when the output (an explicit `--output-format`, else the
/// path's extension) is a parallel index, otherwise empty.
std::string partition_index_kind(const std::string& rPath, const std::string& rFormat) {
    if (!rFormat.empty())
        return rFormat == "pvtu" || rFormat == "pvtp" ? rFormat : std::string();
    std::string ext = std::filesystem::path(rPath).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (ext == ".pvtu")
        return "pvtu";
    if (ext == ".pvtp")
        return "pvtp";
    return {};
}

int cmd_partition(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input-format", {"-i"}, true},
                                  {"output-format", {"-o"}, true},
                                  {"nparts", {"-n"}, true},
                                  {"method", {}, true},
                                  {"imbalance", {}, true},
                                  {"mode", {}, true},
                                  {"seed", {}, true},
                                  {"weights", {}, true},
                                  {"ghost-layers", {}, true},
                                  {"record-ids", {}, false},
                                  {"labels-only", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error(
            "partition requires exactly INFILE and OUTPATTERN (with {part}, a .pvtu/.pvtp "
            "index, or a plain path with --labels-only)");
    if (opt_value(p, "nparts").empty())
        throw std::runtime_error("partition: --nparts N is required");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::PartitionOptions options;
    options.mNParts = std::stoi(opt_value(p, "nparts"));
    options.mMethod = meshioplusplus::partition_method_from_name(opt_value(p, "method", "auto"));
    options.mImbalance = meshioplusplus::detail::stod_c(opt_value(p, "imbalance", "0.03"));
    options.mMode = meshioplusplus::partition_mode_from_name(opt_value(p, "mode", "eco"));
    options.mSeed = std::stoi(opt_value(p, "seed", "0"));
    options.mRecordIds = has_flag(p, "record-ids");
    options.mGhostLayers = std::stoi(opt_value(p, "ghost-layers", "0"));
    options.mWeightsKey = opt_value(p, "weights");

    std::string out_fmt = opt_value(p, "output-format");
    if (has_flag(p, "labels-only")) {
        std::vector<meshioplusplus::NDArray> labels =
            meshioplusplus::partition_labels(mesh, options);
        mesh.AddCellData("partition:part", std::move(labels));
        write_mesh_cli(p.positionals[1], mesh, out_fmt);
        return 0;
    }

    const bool has_token = p.positionals[1].find("{part}") != std::string::npos;
    const std::string index_kind = partition_index_kind(p.positionals[1], out_fmt);
    if (!has_token && index_kind.empty())
        throw std::runtime_error(
            "partition: output pattern must contain '{part}' (e.g. out_{part}.vtu), or be a "
            ".pvtu/.pvtp index (one file per part plus an index), or pass --labels-only");
    auto result = meshioplusplus::partition(mesh, options);
    if (!has_flag(p, "quiet"))
        std::cout << "partitioned into " << result.mPieces.size() << " piece(s)\n";
    if (!has_token) {
        // One index over every piece, so the halo layers survive as vtkGhostType.
        std::vector<const Mesh*> ptrs;
        for (const auto& piece : result.mPieces)
            ptrs.push_back(&piece.mMesh);
        const auto codec =
            meshioplusplus::detail::vtk_codec_available(meshioplusplus::detail::VtkCodec::Zlib)
                ? meshioplusplus::detail::VtkCodec::Zlib
                : meshioplusplus::detail::VtkCodec::None;
        if (index_kind == "pvtu")
            meshioplusplus::write_pvtu_pieces_codec(p.positionals[1], ptrs, /*binary=*/true, codec);
        else
            meshioplusplus::write_pvtp_pieces_codec(p.positionals[1], ptrs, /*binary=*/true, codec);
        if (!has_flag(p, "quiet"))
            std::cout << "  " << ptrs.size() << " pieces -> " << p.positionals[1] << "\n";
        return 0;
    }
    for (auto& piece : result.mPieces) {
        std::string path = partition_replace_part(p.positionals[1], piece.mPartId);
        std::int64_t ncells = 0;
        for (const auto cb : piece.mMesh.CellRange())
            ncells += static_cast<std::int64_t>(cb.NumCells());
        if (!has_flag(p, "quiet"))
            std::cout << "  " << piece.mPartId << ": " << piece.mMesh.NumPoints() << " points, "
                      << ncells << " cells -> " << path << "\n";
        write_mesh_cli(path, piece.mMesh, out_fmt);
    }
    return 0;
}

std::string stats_g6(double v) {
    char buf[32];
    meshioplusplus::detail::snprintf_c(buf, sizeof(buf), "%.6g", v);
    return buf;
}

/// Shared options for `view` and `screenshot`, which differ only in output.
std::vector<cli_opt_spec> view_common_specs() {
    return {
        {"input-format", {"-i"}, true},
        {"kind", {}, true},
        {"color-by", {}, true},
        {"name", {}, true},
    };
}

int cmd_view(const std::vector<std::string>& rArgs) {
    auto specs = view_common_specs();
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 1)
        throw std::runtime_error("view requires exactly INFILE");

    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    meshioplusplus::cli::view_mesh(
        mesh, meshioplusplus::cli::view_kind_from_name(opt_value(p, "kind", "auto")),
        opt_value(p, "color-by"), opt_value(p, "name", "mesh"));
    return 0;
}

int cmd_screenshot(const std::vector<std::string>& rArgs) {
    auto specs = view_common_specs();
    specs.push_back({"size", {}, true});
    specs.push_back({"transparent", {}, false});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 2)
        throw std::runtime_error("screenshot requires exactly INFILE and OUTFILE");

    // WxH rather than two positionals, so the parser needs no multi-value
    // option and OUTFILE stays unambiguous.
    int width = 1280;
    int height = 960;
    const std::string size = opt_value(p, "size");
    if (!size.empty()) {
        const auto x = size.find('x');
        if (x == std::string::npos)
            throw std::runtime_error("--size expects WIDTHxHEIGHT, e.g. --size 1600x1200");
        width = std::stoi(size.substr(0, x));
        height = std::stoi(size.substr(x + 1));
    }

    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    meshioplusplus::cli::screenshot_mesh(
        mesh, meshioplusplus::cli::view_kind_from_name(opt_value(p, "kind", "auto")),
        opt_value(p, "color-by"), opt_value(p, "name", "mesh"), p.positionals[1], width, height,
        has_flag(p, "transparent"));
    return 0;
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
        meshioplusplus::cli::JsonOut json(std::cout);
        auto vec = [&](const char* pKey, const double* pV) {
            json.Key(pKey);
            json.BeginArray();
            for (int k = 0; k < 3; ++k)
                json.Number(pV[k]);
            json.EndArray();
        };
        json.BeginObject();
        json.Field("num_points", static_cast<std::int64_t>(s.mNumPoints));
        json.Field("num_cells", static_cast<std::int64_t>(s.mNumCells));
        vec("bbox_min", s.mBBoxMin);
        vec("bbox_max", s.mBBoxMax);
        vec("extent", s.mExtent);
        vec("centroid", s.mCentroid);
        json.Key("cell_type_counts");
        json.BeginObject();
        for (const auto& kv : s.mCellTypeCounts)
            json.Field(kv.first, static_cast<std::int64_t>(kv.second));
        json.EndObject();
        json.Field("total_area", s.mTotalArea);
        json.Field("signed_volume", s.mSignedVolume);
        json.Field("unsigned_volume", s.mUnsignedVolume);
        json.Field("num_inverted", static_cast<std::int64_t>(s.mNumInverted));
        json.EndObject();
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

// --------------------------------------------------------------------------
// `data` — the nested verb group over the data operations
//
// These act on the mesh's point/cell/field data arrays; the geometry is never
// modified. Mirrors the Python CLI's nine verbs one-for-one, including the
// colon-splitting rules, so the two surfaces stay interchangeable.
// --------------------------------------------------------------------------

/// The three location flags every data verb accepts, in report order.
const char* const DATA_LOCATION_FLAGS[3] = {"point", "cell", "field"};

meshioplusplus::DataLocation data_location_of_flag(int i) {
    switch (i) {
        case 0:
            return meshioplusplus::DataLocation::Point;
        case 1:
            return meshioplusplus::DataLocation::Cell;
        default:
            return meshioplusplus::DataLocation::Field;
    }
}

/// "OLD:NEW" split on the LAST colon, because data names routinely contain
/// colons (`gmsh:physical`). Identical rule to the Python CLI and doc/cli.md.
std::pair<std::string, std::string> data_split_rename(const std::string& rValue) {
    const std::size_t pos = rValue.rfind(':');
    if (pos == std::string::npos || pos == 0 || pos + 1 >= rValue.size())
        throw std::runtime_error("data rename: expected OLD:NEW, got '" + rValue + "'");
    return {rValue.substr(0, pos), rValue.substr(pos + 1)};
}

/// "NAME = EXPR" split on the FIRST '='.
std::pair<std::string, std::string> data_split_assignment(const std::string& rValue) {
    const std::size_t pos = rValue.find('=');
    if (pos == std::string::npos)
        throw std::runtime_error("data calc: expected 'NAME = EXPRESSION', got '" + rValue + "'");
    auto trim = [](std::string s) {
        const std::size_t b = s.find_first_not_of(" \t");
        const std::size_t e = s.find_last_not_of(" \t");
        return b == std::string::npos ? std::string() : s.substr(b, e - b + 1);
    };
    const std::string name = trim(rValue.substr(0, pos));
    const std::string expr = trim(rValue.substr(pos + 1));
    if (name.empty() || expr.empty())
        throw std::runtime_error("data calc: expected 'NAME = EXPRESSION', got '" + rValue + "'");
    return {name, expr};
}

/// The IO options every data verb shares.
std::vector<cli_opt_spec> data_io_specs(bool with_output) {
    std::vector<cli_opt_spec> specs = {{"input-format", {"-i"}, true}};
    if (with_output)
        specs.push_back({"output-format", {"-o"}, true});
    return specs;
}

std::string data_g6(double v) {
    char buf[32];
    meshioplusplus::detail::snprintf_c(buf, sizeof(buf), "%.6g", v);
    return buf;
}

int cmd_data_info(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(false);
    specs.push_back({"json", {}, false});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 1)
        throw std::runtime_error("data info requires exactly INFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    meshioplusplus::DataInfoReport r = meshioplusplus::data_info(mesh);

    if (has_flag(p, "json")) {
        // The keys of the Python `data_info` report, in its order.
        meshioplusplus::cli::JsonOut json(std::cout);
        auto numbers = [&](const char* pKey, const std::vector<double>& rV) {
            json.Key(pKey);
            json.BeginArray();
            for (double v : rV)
                json.Number(v);
            json.EndArray();
        };
        json.BeginArray();
        for (const auto& a : r.mArrays) {
            json.BeginObject();
            json.Field("location", meshioplusplus::data_location_name(a.mLocation));
            json.Field("name", a.mName);
            json.Field("dtype", meshioplusplus::dtype_numpy_str(a.mDtype));
            json.Key("shape");
            json.BeginArray();
            for (auto d : a.mShape)
                json.Int(static_cast<std::int64_t>(d));
            json.EndArray();
            json.Field("num_blocks", static_cast<std::int64_t>(a.mNumBlocks));
            json.Field("num_entries", static_cast<std::int64_t>(a.mNumEntries));
            json.Field("num_components", static_cast<std::int64_t>(a.mNumComponents));
            json.Field("num_values", static_cast<std::int64_t>(a.mNumValues));
            json.Field("min", a.mMin);
            json.Field("max", a.mMax);
            json.Field("mean", a.mMean);
            numbers("min_per_component", a.mMinPerComponent);
            numbers("max_per_component", a.mMaxPerComponent);
            numbers("mean_per_component", a.mMeanPerComponent);
            json.Field("num_nan", static_cast<std::int64_t>(a.mNumNan));
            json.Field("num_inf", static_cast<std::int64_t>(a.mNumInf));
            json.Field("num_finite", static_cast<std::int64_t>(a.mNumFinite));
            json.Field("inconsistent_blocks", a.mInconsistentBlocks);
            json.EndObject();
        }
        json.EndArray();
        return 0;
    }

    std::cout << "<meshio++ data summary>\n";
    if (r.mArrays.empty()) {
        std::cout << "  (the mesh carries no data arrays)\n";
        return 0;
    }
    std::cout << "  location    name                 dtype  comp  entries          min"
                 "          max         mean   nan   inf\n";
    std::cout << "  " << std::string(104, '-') << "\n";
    for (const auto& a : r.mArrays) {
        char line[512];
        std::snprintf(
            line, sizeof(line), "  %-11s %-20s %-6s %4lld %8lld %12s %12s %12s %5lld %5lld",
            meshioplusplus::data_location_name(a.mLocation), a.mName.c_str(),
            meshioplusplus::dtype_numpy_str(a.mDtype), static_cast<long long>(a.mNumComponents),
            static_cast<long long>(a.mNumEntries), data_g6(a.mMin).c_str(), data_g6(a.mMax).c_str(),
            data_g6(a.mMean).c_str(), static_cast<long long>(a.mNumNan),
            static_cast<long long>(a.mNumInf));
        std::cout << line << "\n";
        if (a.mInconsistentBlocks)
            std::cout << "              (warning: cell blocks disagree in components)\n";
    }
    return 0;
}

int cmd_data_rename(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(true);
    for (const char* loc : DATA_LOCATION_FLAGS)
        specs.push_back({loc, {}, true});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 2)
        throw std::runtime_error("data rename requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    int total = 0;
    for (int i = 0; i < 3; ++i) {
        for (const std::string& spec : opt_values(p, DATA_LOCATION_FLAGS[i])) {
            auto [old_name, new_name] = data_split_rename(spec);
            mesh = meshioplusplus::data_rename(mesh, data_location_of_flag(i), old_name, new_name);
            std::cout << "renamed " << meshioplusplus::data_location_name(data_location_of_flag(i))
                      << " '" << old_name << "' -> '" << new_name << "'\n";
            ++total;
        }
    }
    if (total == 0)
        std::cout << "data rename: nothing to do (pass --point/--cell/--field OLD:NEW)\n";
    write_mesh_cli(p.positionals[1], mesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_data_drop(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(true);
    for (const char* loc : DATA_LOCATION_FLAGS)
        specs.push_back({loc, {}, true});
    specs.push_back({"ignore-missing", {}, false});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 2)
        throw std::runtime_error("data drop requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    const bool ignore = has_flag(p, "ignore-missing");

    for (int i = 0; i < 3; ++i) {
        if (!has_opt(p, DATA_LOCATION_FLAGS[i]))
            continue;
        const std::vector<std::string> names =
            data_split_names(opt_value(p, DATA_LOCATION_FLAGS[i]));
        if (names.empty())
            continue;
        mesh = meshioplusplus::data_drop(mesh, data_location_of_flag(i), names, ignore);
        std::cout << "dropped " << meshioplusplus::data_location_name(data_location_of_flag(i))
                  << ": ";
        for (std::size_t k = 0; k < names.size(); ++k)
            std::cout << (k ? ", " : "") << names[k];
        std::cout << "\n";
    }
    write_mesh_cli(p.positionals[1], mesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_data_keep(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(true);
    for (const char* loc : DATA_LOCATION_FLAGS)
        specs.push_back({loc, {}, true});
    specs.push_back({"ignore-missing", {}, false});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 2)
        throw std::runtime_error("data keep requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));
    const bool ignore = has_flag(p, "ignore-missing");

    for (int i = 0; i < 3; ++i) {
        // A location the user did not name is left alone; naming it with an
        // empty list means "keep nothing there".
        if (!has_opt(p, DATA_LOCATION_FLAGS[i]))
            continue;
        const std::vector<std::string> names =
            data_split_names(opt_value(p, DATA_LOCATION_FLAGS[i]));
        mesh = meshioplusplus::data_keep(mesh, data_location_of_flag(i), names, ignore);
        std::cout << "kept " << meshioplusplus::data_location_name(data_location_of_flag(i))
                  << ": ";
        if (names.empty()) {
            std::cout << "(nothing)";
        } else {
            for (std::size_t k = 0; k < names.size(); ++k)
                std::cout << (k ? ", " : "") << names[k];
        }
        std::cout << "\n";
    }
    write_mesh_cli(p.positionals[1], mesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_data_to_cell(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(true);
    specs.push_back({"keys", {}, true});
    specs.push_back({"target-suffix", {}, true});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 2)
        throw std::runtime_error("data to-cell requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::DataAverageOptions opts;
    opts.names = data_split_names(opt_value(p, "keys"));
    opts.suffix = opt_value(p, "target-suffix");
    const std::size_t n = opts.names.empty() ? mesh.NumPointData() : opts.names.size();
    Mesh out = meshioplusplus::point_data_to_cell_data(mesh, opts);
    std::cout << "averaged " << n << " point_data array(s) onto the cells\n";
    write_mesh_cli(p.positionals[1], out, opt_value(p, "output-format"));
    return 0;
}

int cmd_data_to_point(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(true);
    specs.push_back({"keys", {}, true});
    specs.push_back({"weighted", {}, false});
    specs.push_back({"target-suffix", {}, true});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 2)
        throw std::runtime_error("data to-point requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::DataAverageOptions opts;
    opts.names = data_split_names(opt_value(p, "keys"));
    opts.weight = has_flag(p, "weighted") ? meshioplusplus::CellPointWeight::Measure
                                          : meshioplusplus::CellPointWeight::Uniform;
    opts.suffix = opt_value(p, "target-suffix");
    const std::size_t n = opts.names.empty() ? mesh.NumCellData() : opts.names.size();
    Mesh out = meshioplusplus::cell_data_to_point_data(mesh, opts);
    std::cout << "averaged " << n << " cell_data array(s) onto the points ("
              << (has_flag(p, "weighted") ? "measure-weighted" : "unweighted") << ")\n";
    write_mesh_cli(p.positionals[1], out, opt_value(p, "output-format"));
    return 0;
}

int cmd_data_calc(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(true);
    for (const char* loc : DATA_LOCATION_FLAGS)
        specs.push_back({loc, {}, true});
    specs.push_back({"overwrite", {}, false});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 2)
        throw std::runtime_error("data calc requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    int total = 0;
    for (int i = 0; i < 3; ++i) {
        for (const std::string& spec : opt_values(p, DATA_LOCATION_FLAGS[i])) {
            auto [name, expr] = data_split_assignment(spec);
            meshioplusplus::DataCalcOptions opts;
            opts.location = data_location_of_flag(i);
            opts.output = name;
            opts.overwrite = has_flag(p, "overwrite");
            mesh = meshioplusplus::data_calc(mesh, expr, opts);
            std::cout << "computed " << meshioplusplus::data_location_name(opts.location) << " '"
                      << name << "' = " << expr << "\n";
            ++total;
        }
    }
    if (total == 0)
        std::cout << "data calc: nothing to do (pass --point/--cell/--field 'NAME = EXPR')\n";
    write_mesh_cli(p.positionals[1], mesh, opt_value(p, "output-format"));
    return 0;
}

/// Shared body of `data clamp` and `data normalize`.
int cmd_data_condition_impl(const std::vector<std::string>& rArgs, bool clamp) {
    auto specs = data_io_specs(true);
    for (const char* loc : DATA_LOCATION_FLAGS)
        specs.push_back({loc, {}, true});
    specs.push_back({"magnitude", {}, false});
    specs.push_back({"nan", {}, true});
    specs.push_back({"nan-value", {}, true});
    specs.push_back({"suffix", {}, true});
    if (clamp) {
        specs.push_back({"min", {}, true});
        specs.push_back({"max", {}, true});
    } else {
        specs.push_back({"to", {}, true});
        specs.push_back({"zero-mean", {}, false});
    }
    auto p = cli_parse(rArgs, specs);
    const char* verb = clamp ? "data clamp" : "data normalize";
    if (p.positionals.size() != 2)
        throw std::runtime_error(std::string(verb) + " requires exactly INFILE and OUTFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    const bool zero_mean = !clamp && has_flag(p, "zero-mean");
    double lo = 0.0;
    double hi = 1.0;
    if (clamp) {
        if (!has_opt(p, "min") || !has_opt(p, "max"))
            throw std::runtime_error("data clamp requires --min and --max");
        lo = meshioplusplus::detail::stod_c(opt_value(p, "min"));
        hi = meshioplusplus::detail::stod_c(opt_value(p, "max"));
    } else if (!zero_mean) {
        const std::vector<std::string> parts = data_split_names(opt_value(p, "to", "0,1"));
        if (parts.size() != 2)
            throw std::runtime_error("data normalize: --to expects LO,HI");
        lo = meshioplusplus::detail::stod_c(parts[0]);
        hi = meshioplusplus::detail::stod_c(parts[1]);
    }

    int touched = 0;
    for (int i = 0; i < 3; ++i) {
        if (!has_opt(p, DATA_LOCATION_FLAGS[i]))
            continue;
        meshioplusplus::DataConditionOptions opts;
        opts.location = data_location_of_flag(i);
        opts.names = data_split_names(opt_value(p, DATA_LOCATION_FLAGS[i]));
        opts.mode = clamp ? meshioplusplus::ConditionMode::Clamp
                          : (zero_mean ? meshioplusplus::ConditionMode::Standardize
                                       : meshioplusplus::ConditionMode::Normalize);
        opts.scope = has_flag(p, "magnitude") ? meshioplusplus::ConditionScope::Magnitude
                                              : meshioplusplus::ConditionScope::Component;
        opts.lo = lo;
        opts.hi = hi;
        opts.nan_policy = meshioplusplus::nan_policy_from_name(opt_value(p, "nan", "ignore"));
        opts.nan_replacement = meshioplusplus::detail::stod_c(opt_value(p, "nan-value", "0"));
        opts.suffix = opt_value(p, "suffix");
        mesh = meshioplusplus::data_condition(mesh, opts);
        std::cout << (clamp ? "clamped " : "normalized ")
                  << meshioplusplus::data_location_name(opts.location) << " ";
        if (opts.names.empty()) {
            std::cout << "(all)";
        } else {
            for (std::size_t k = 0; k < opts.names.size(); ++k)
                std::cout << (k ? ", " : "") << opts.names[k];
        }
        if (clamp)
            std::cout << " to [" << data_g6(lo) << ", " << data_g6(hi) << "]\n";
        else if (zero_mean)
            std::cout << " to zero mean / unit std\n";
        else
            std::cout << " to [" << data_g6(lo) << ", " << data_g6(hi) << "]\n";
        ++touched;
    }
    if (touched == 0)
        std::cout << verb << ": nothing to do (pass --point/--cell/--field NAME)\n";
    write_mesh_cli(p.positionals[1], mesh, opt_value(p, "output-format"));
    return 0;
}

// `data invariants` derives von Mises / principal / hydrostatic / deviatoric
// fields from a symmetric (6-component) or general 3x3 (9-component) tensor
// array. See doc/tensor_invariants.md.
int cmd_data_invariants(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(true);
    specs.push_back({"point", {}, true});
    specs.push_back({"cell", {}, true});
    specs.push_back({"outputs", {}, true});
    specs.push_back({"prefix", {}, true});
    specs.push_back({"suffix", {}, true});
    specs.push_back({"no-overwrite", {}, false});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 2)
        throw std::runtime_error("data invariants requires exactly INFILE and OUTFILE");
    if (has_opt(p, "point") == has_opt(p, "cell"))
        throw std::runtime_error("data invariants requires exactly one of --point NAMES or --cell NAMES");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::TensorInvariantsOptions opts;
    const bool cell = has_opt(p, "cell");
    opts.location = cell ? meshioplusplus::DataLocation::Cell : meshioplusplus::DataLocation::Point;
    opts.names = data_split_names(opt_value(p, cell ? "cell" : "point"));
    opts.outputs = has_opt(p, "outputs") ? meshioplusplus::tensor_invariant_from_name(opt_value(p, "outputs"))
                                         : meshioplusplus::TensorInvariant::All;
    opts.prefix = opt_value(p, "prefix");
    opts.suffix = opt_value(p, "suffix");
    opts.overwrite = !has_flag(p, "no-overwrite");

    mesh = meshioplusplus::tensor_invariants(mesh, opts);
    std::cout << "computed invariants for " << meshioplusplus::data_location_name(opts.location)
              << " ";
    if (opts.names.empty()) {
        std::cout << "(all tensor-shaped)";
    } else {
        for (std::size_t k = 0; k < opts.names.size(); ++k)
            std::cout << (k ? ", " : "") << opts.names[k];
    }
    std::cout << "\n";
    write_mesh_cli(p.positionals[1], mesh, opt_value(p, "output-format"));
    return 0;
}

int cmd_data_clamp(const std::vector<std::string>& rArgs) {
    return cmd_data_condition_impl(rArgs, /*clamp=*/true);
}

int cmd_data_normalize(const std::vector<std::string>& rArgs) {
    return cmd_data_condition_impl(rArgs, /*clamp=*/false);
}

// `data gradient` is a slight departure from the rest of this group: it is a
// mesh operation (it reads geometry and topology), not one of the data_* family.
// It lives here because it consumes and produces data arrays, which is where a
// user looks for it. See doc/gradient.md.
int cmd_data_gradient(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(true);
    specs.push_back({"array", {}, true});
    specs.push_back({"op", {}, true});
    specs.push_back({"method", {}, true});
    specs.push_back({"location", {}, true});
    specs.push_back({"output", {}, true});
    specs.push_back({"component", {}, true});
    specs.push_back({"overwrite", {}, false});
    specs.push_back({"quiet", {"-q"}, false});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 2)
        throw std::runtime_error("data gradient requires exactly INFILE and OUTFILE");
    if (!has_opt(p, "array"))
        throw std::runtime_error("data gradient requires --array NAME");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::GradientOptions opts;
    opts.mArrayName = opt_value(p, "array");
    // The defaults here must stay string-identical to the Python CLI's.
    opts.mOperator = meshioplusplus::gradient_operator_from_name(opt_value(p, "op", "gradient"));
    opts.mMethod = meshioplusplus::gradient_method_from_name(opt_value(p, "method", "green-gauss"));
    opts.mLocation = meshioplusplus::data_location_from_name(opt_value(p, "location", "cell"));
    opts.mOutputName = opt_value(p, "output");
    if (has_opt(p, "component"))
        opts.mComponent = std::stoi(opt_value(p, "component"));
    opts.mOverwrite = has_flag(p, "overwrite");

    meshioplusplus::GradientResult r = meshioplusplus::gradient(mesh, opts);
    if (!has_flag(p, "quiet")) {
        const std::string name =
            opts.mOutputName.empty()
                ? opts.mArrayName + std::string(opt_value(p, "op", "gradient") == "divergence"
                                                    ? meshioplusplus::kDivergenceSuffix
                                                    : (opt_value(p, "op", "gradient") == "curl"
                                                           ? meshioplusplus::kCurlSuffix
                                                           : meshioplusplus::kGradientSuffix))
                : opts.mOutputName;
        std::cout << "wrote " << meshioplusplus::data_location_name(opts.mLocation) << " '" << name
                  << "' (" << opt_value(p, "method", "green-gauss") << ")\n";
        std::cout << "  cells skipped (NaN):      " << r.mNumSkipped << "\n";
        std::cout << "  cells fell back to GG:    " << r.mNumFallback << "\n";
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

// `data hessian` is `data gradient`'s companion one order further: a
// composition of TWO gradient() calls, not a new kernel. See
// doc/hessian.md.
int cmd_data_hessian(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(true);
    specs.push_back({"array", {}, true});
    specs.push_back({"method", {}, true});
    specs.push_back({"location", {}, true});
    specs.push_back({"output", {}, true});
    specs.push_back({"overwrite", {}, false});
    specs.push_back({"quiet", {"-q"}, false});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 2)
        throw std::runtime_error("data hessian requires exactly INFILE and OUTFILE");
    if (!has_opt(p, "array"))
        throw std::runtime_error("data hessian requires --array NAME");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::HessianOptions opts;
    opts.mArrayName = opt_value(p, "array");
    // The defaults here must stay string-identical to the Python CLI's.
    opts.mMethod = meshioplusplus::gradient_method_from_name(opt_value(p, "method", "green-gauss"));
    opts.mLocation = meshioplusplus::data_location_from_name(opt_value(p, "location", "cell"));
    opts.mOutputName = opt_value(p, "output");
    opts.mOverwrite = has_flag(p, "overwrite");

    meshioplusplus::HessianResult r = meshioplusplus::hessian(mesh, opts);
    if (!has_flag(p, "quiet")) {
        const std::string name = opts.mOutputName.empty()
                                     ? opts.mArrayName + std::string(meshioplusplus::kHessianSuffix)
                                     : opts.mOutputName;
        std::cout << "wrote " << meshioplusplus::data_location_name(opts.mLocation) << " '" << name
                  << "' (" << opt_value(p, "method", "green-gauss") << ")\n";
        std::cout << "  cells skipped (NaN):      " << r.mNumSkipped << "\n";
        std::cout << "  cells fell back to GG:    " << r.mNumFallback << "\n";
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

// `data integrate` is `data gradient`'s companion the other direction: a
// **mesh** operation (it reads cell measures), living here because it reduces
// data arrays, not because it belongs to the no-geometry data_* bundle. See
// doc/field_integration.md.
int cmd_data_integrate(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(false);
    specs.push_back({"array", {}, true});
    specs.push_back({"json", {}, false});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 1)
        throw std::runtime_error("data integrate requires exactly INFILE");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::DataIntegrateOptions opts;
    opts.mArrayNames = opt_values(p, "array");
    meshioplusplus::DataIntegrateReport r = meshioplusplus::data_integrate(mesh, opts);

    auto print_region = [](const meshioplusplus::FieldIntegralRegion& reg) {
        const std::string label = reg.mName.empty() ? "(whole mesh)" : reg.mName;
        char line[128];
        std::snprintf(line, sizeof(line), "    %-20s cells=%-6lld skipped=%-4lld", label.c_str(),
                      static_cast<long long>(reg.mNumCells),
                      static_cast<long long>(reg.mNumSkipped));
        std::cout << line << " total=[";
        for (std::size_t k = 0; k < reg.mTotalPerComponent.size(); ++k)
            std::cout << (k ? ", " : "") << data_g6(reg.mTotalPerComponent[k]);
        std::cout << "] mean=[";
        for (std::size_t k = 0; k < reg.mMeanPerComponent.size(); ++k)
            std::cout << (k ? ", " : "") << data_g6(reg.mMeanPerComponent[k]);
        std::cout << "]\n";
    };

    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        auto region = [&](const meshioplusplus::FieldIntegralRegion& reg) {
            auto numbers = [&](const char* pKey, const std::vector<double>& rV) {
                json.Key(pKey);
                json.BeginArray();
                for (double v : rV)
                    json.Number(v);
                json.EndArray();
            };
            json.BeginObject();
            json.Field("name", reg.mName);
            json.Field("num_cells", static_cast<std::int64_t>(reg.mNumCells));
            json.Field("num_skipped", static_cast<std::int64_t>(reg.mNumSkipped));
            numbers("domain_measure_per_component", reg.mDomainMeasurePerComponent);
            numbers("total_per_component", reg.mTotalPerComponent);
            numbers("mean_per_component", reg.mMeanPerComponent);
            json.Key("num_nan_per_component");
            json.BeginArray();
            for (auto v : reg.mNumNanPerComponent)
                json.Int(static_cast<std::int64_t>(v));
            json.EndArray();
            json.EndObject();
        };
        json.BeginArray();
        for (const auto& a : r.mArrays) {
            json.BeginObject();
            json.Field("name", a.mName);
            json.Field("num_components", static_cast<std::int64_t>(a.mNumComponents));
            json.Key("domain");
            region(a.mDomain);
            json.Key("regions");
            json.BeginArray();
            for (const auto& reg : a.mRegions)
                region(reg);
            json.EndArray();
            json.EndObject();
        }
        json.EndArray();
        return 0;
    }

    std::cout << "<meshio++ field integration>\n";
    if (r.mArrays.empty()) {
        std::cout << "  (the mesh carries no cell_data arrays)\n";
        return 0;
    }
    for (const auto& a : r.mArrays) {
        std::cout << "  " << a.mName << " (" << a.mNumComponents << " component(s))\n";
        print_region(a.mDomain);
        for (const auto& reg : a.mRegions)
            print_region(reg);
    }
    return 0;
}

// `data estimate-error` is a slight departure from the rest of this group too,
// for the same reason `data gradient` is: a mesh operation (it reads geometry
// and topology), living here because it consumes and produces data arrays. It
// composes `gradient` with the point<->cell averaging round trip to produce
// the ZZ recovery-based error indicator `refine`'s own `--where` consumes. See
// doc/error.md.
int cmd_data_estimate_error(const std::vector<std::string>& rArgs) {
    auto specs = data_io_specs(true);
    specs.push_back({"array", {}, true});
    specs.push_back({"method", {}, true});
    specs.push_back({"marking", {}, true});
    specs.push_back({"marking-value", {}, true});
    specs.push_back({"output", {}, true});
    specs.push_back({"marked", {}, true});
    specs.push_back({"overwrite", {}, false});
    specs.push_back({"quiet", {"-q"}, false});
    auto p = cli_parse(rArgs, specs);
    if (p.positionals.size() != 2)
        throw std::runtime_error("data estimate-error requires exactly INFILE and OUTFILE");
    if (!has_opt(p, "array"))
        throw std::runtime_error("data estimate-error requires --array NAME");
    Mesh mesh = read_mesh_cli(p.positionals[0], opt_value(p, "input-format"));

    meshioplusplus::ErrorOptions opts;
    opts.mArrayName = opt_value(p, "array");
    // The defaults here must stay string-identical to the Python CLI's.
    opts.mMethod = meshioplusplus::error_method_from_name(opt_value(p, "method", "zz"));
    opts.mMarking = meshioplusplus::error_marking_from_name(opt_value(p, "marking", "none"));
    if (has_opt(p, "marking-value"))
        opts.mMarkingValue = meshioplusplus::detail::stod_c(opt_value(p, "marking-value"));
    opts.mOutputName = opt_value(p, "output");
    opts.mMarkedName = opt_value(p, "marked");
    opts.mOverwrite = has_flag(p, "overwrite");

    meshioplusplus::ErrorResult r = meshioplusplus::estimate_error(mesh, opts);
    if (!has_flag(p, "quiet")) {
        const std::string name = opts.mOutputName.empty() ? "error:zz" : opts.mOutputName;
        std::cout << "wrote cell_data '" << name << "' (" << opt_value(p, "method", "zz") << ")\n";
        std::cout << "  global error:             " << r.mGlobalError << "\n";
        std::cout << "  cells skipped (NaN):      " << r.mNumSkipped << "\n";
        const std::string marking = opt_value(p, "marking", "none");
        if (marking != "none") {
            const std::string marked_name =
                opts.mMarkedName.empty() ? "error:marked" : opts.mMarkedName;
            std::cout << "  wrote cell_data '" << marked_name << "' (" << marking << ")\n";
            std::cout << "  cells marked:             " << r.mNumMarked << "\n";
        }
    }
    write_mesh_cli(p.positionals[1], r.mMesh, opt_value(p, "output-format"));
    return 0;
}

/// One report counter, printed as an integer when it is one (most are counts;
/// only e.g. Smooth's MaxDisplacement is genuinely fractional).
void pipeline_print_counter(std::ostream& rOut, double value) {
    if (value == std::floor(value) && std::abs(value) < 1e15)
        rOut << static_cast<long long>(value);
    else
        rOut << value;
}

int cmd_pipeline(const std::vector<std::string>& rArgs) {
    auto p = cli_parse(rArgs, {
                                  {"input", {}, true},
                                  {"output", {}, true},
                                  {"json", {}, false},
                                  {"quiet", {"-q"}, false},
                              });
    if (p.positionals.size() != 1)
        throw std::runtime_error("pipeline requires exactly SETTINGS.json");
    // Needs a build with the JSON parser (-DMESHIOPLUSPLUS_WITH_JSON=ON, the
    // default when the submodule is checked out); otherwise this throws naming
    // the flag -- the view/screenshot contract, the verb always exists.
    // One parse for both document shapes: a sequence document (a glob/list
    // Input, a {step} Output, Mode/Parallel/Workers) runs through the sequence
    // driver, and a plain single-file one takes the physically unchanged
    // v9.11.0 path -- so nobody has to know which kind of document they hold.
    meshioplusplus::SequencePipeline pipeline =
        meshioplusplus::parse_sequence_file(p.positionals[0]);
    if (has_opt(p, "input")) {
        pipeline.mInput.mPaths = {opt_value(p, "input")};
        pipeline.mInput.mPattern.clear();
    }
    if (has_opt(p, "output"))
        pipeline.mOutput.mPath = opt_value(p, "output");
    const meshioplusplus::PipelineReport report = meshioplusplus::run_sequence_pipeline(pipeline);

    if (has_flag(p, "json")) {
        meshioplusplus::cli::JsonOut json(std::cout);
        json.BeginObject();
        json.Key("steps");
        json.BeginArray();
        for (const auto& step : report.mSteps) {
            json.BeginObject();
            json.Field("op", step.mOp);
            for (const auto& counter : step.mCounters) {
                // A counter is a count unless it has a fractional part.
                const double v = counter.second;
                json.Key(counter.first);
                if (std::isfinite(v) && v == std::floor(v) && std::fabs(v) < 9.0e15)
                    json.Int(static_cast<std::int64_t>(v));
                else
                    json.Number(v);
            }
            json.EndObject();
        }
        json.EndArray();
        json.Key("warnings");
        json.BeginArray();
        for (const auto& warning : report.mWarnings)
            json.String(warning);
        json.EndArray();
        json.EndObject();
    } else if (!has_flag(p, "quiet")) {
        for (std::size_t i = 0; i < report.mSteps.size(); ++i) {
            const auto& step = report.mSteps[i];
            std::cout << "step " << (i + 1) << ": " << step.mOp;
            for (std::size_t c = 0; c < step.mCounters.size(); ++c) {
                std::cout << (c ? ", " : " (") << step.mCounters[c].first << "=";
                pipeline_print_counter(std::cout, step.mCounters[c].second);
                if (c + 1 == step.mCounters.size())
                    std::cout << ")";
            }
            std::cout << "\n";
        }
        for (const auto& warning : report.mWarnings)
            std::cout << "warning: " << warning << "\n";
        std::cout << "wrote " << pipeline.mOutput.mPath << "\n";
    }
    return 0;
}

void print_data_usage(std::ostream& rOut) {
    rOut << "usage: meshioplusplus data <subcommand> [options]\n\n"
            "Operations on a mesh's data arrays (the geometry is never modified).\n\n"
            "subcommands:\n"
            "  info        Summarize every data array (dtype/shape/min/max/mean/NaN)\n"
            "  rename      Rename data arrays (--point OLD:NEW, split on the last ':')\n"
            "  drop        Drop data arrays by name (--point A,B)\n"
            "  keep        Keep only the named data arrays (--point T,p)\n"
            "  to-cell     Average point_data onto the cells\n"
            "  to-point    Average cell_data onto the points (--weighted)\n"
            "  calc        Derive an array from an expression (--point 'n = norm(v)')\n"
            "  clamp       Clamp values into [--min, --max]\n"
            "  normalize   Rescale to --to LO,HI (or --zero-mean)\n"
            "  gradient    Differentiate a point_data field (--array NAME --op "
            "gradient|divergence|curl)\n"
            "  hessian     Second derivative of a scalar point_data field "
            "(--array NAME), gradient's companion one order further\n"
            "  integrate   Cell-measure-weighted total/mean of a cell_data field, "
            "whole mesh and per named Cell region (--array NAME, repeatable)\n"
            "  estimate-error  ZZ recovery-based error indicator (--array NAME "
            "--marking none|absolute|fraction|dorfler)\n\n";
}

int cmd_data(const std::vector<std::string>& rArgs) {
    if (rArgs.empty()) {
        print_data_usage(std::cerr);
        return 2;
    }
    const std::string sub = rArgs[0];
    const std::vector<std::string> rest(rArgs.begin() + 1, rArgs.end());
    if (sub == "-h" || sub == "--help") {
        print_data_usage(std::cout);
        return 0;
    }
    if (sub == "info")
        return cmd_data_info(rest);
    if (sub == "rename")
        return cmd_data_rename(rest);
    if (sub == "drop")
        return cmd_data_drop(rest);
    if (sub == "keep")
        return cmd_data_keep(rest);
    if (sub == "to-cell")
        return cmd_data_to_cell(rest);
    if (sub == "to-point")
        return cmd_data_to_point(rest);
    if (sub == "calc")
        return cmd_data_calc(rest);
    if (sub == "clamp")
        return cmd_data_clamp(rest);
    if (sub == "normalize")
        return cmd_data_normalize(rest);
    if (sub == "invariants")
        return cmd_data_invariants(rest);
    if (sub == "gradient")
        return cmd_data_gradient(rest);
    if (sub == "hessian")
        return cmd_data_hessian(rest);
    if (sub == "integrate")
        return cmd_data_integrate(rest);
    if (sub == "estimate-error")
        return cmd_data_estimate_error(rest);
    std::cerr << "error: unknown data subcommand '" << sub << "'\n\n";
    print_data_usage(std::cerr);
    return 2;
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
    meshioplusplus::detail::snprintf_c(buf, sizeof(buf), "%.3e", v);
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

// `diff --json`: the Python `diff` report (the pybind dict), with `regions`;
// the object is left open so the caller can add `equal`.
void cli_diff_json(meshioplusplus::cli::JsonOut& rJson, const meshioplusplus::DiffReport& rRep) {
    auto array_diff = [&](const meshioplusplus::ArrayDiff& rAd) {
        rJson.BeginObject();
        rJson.Field("name", rAd.mName);
        rJson.Field("shape_mismatch", rAd.mShapeMismatch);
        rJson.Field("size_a", static_cast<std::int64_t>(rAd.mSizeA));
        rJson.Field("size_b", static_cast<std::int64_t>(rAd.mSizeB));
        rJson.Field("max_abs_error", rAd.mMaxAbsError);
        rJson.Field("max_rel_error", rAd.mMaxRelError);
        rJson.Field("worst_index", static_cast<std::int64_t>(rAd.mWorstIndex));
        rJson.Field("num_exceeding", static_cast<std::int64_t>(rAd.mNumExceeding));
        rJson.Field("exact", rAd.mExact);
        rJson.EndObject();
    };
    auto data_diff = [&](const char* pKey, const meshioplusplus::DataDiff& rDd) {
        rJson.Key(pKey);
        rJson.BeginObject();
        cli_names_json(rJson, "only_in_a", rDd.mOnlyInA);
        cli_names_json(rJson, "only_in_b", rDd.mOnlyInB);
        rJson.Key("shared");
        rJson.BeginArray();
        for (const auto& ad : rDd.mShared)
            array_diff(ad);
        rJson.EndArray();
        rJson.EndObject();
    };
    rJson.BeginObject();
    rJson.Field("verdict", meshioplusplus::diff_verdict_name(rRep.mVerdict));
    rJson.Field("unordered", rRep.mUnordered);
    rJson.Field("correspondence_failed", rRep.mCorrespondenceFailed);
    rJson.Field("point_count_mismatch", rRep.mPointCountMismatch);
    rJson.Field("num_points_a", static_cast<std::int64_t>(rRep.mNumPointsA));
    rJson.Field("num_points_b", static_cast<std::int64_t>(rRep.mNumPointsB));
    rJson.Key("points");
    array_diff(rRep.mPoints);
    rJson.Field("block_count_mismatch", rRep.mBlockCountMismatch);
    rJson.Field("num_blocks_a", static_cast<std::int64_t>(rRep.mNumBlocksA));
    rJson.Field("num_blocks_b", static_cast<std::int64_t>(rRep.mNumBlocksB));
    rJson.Key("blocks");
    rJson.BeginArray();
    for (const auto& bd : rRep.mBlocks) {
        rJson.BeginObject();
        rJson.Field("block", static_cast<std::int64_t>(bd.mBlock));
        rJson.Field("type_a", bd.mTypeA);
        rJson.Field("type_b", bd.mTypeB);
        rJson.Field("count_a", static_cast<std::int64_t>(bd.mCountA));
        rJson.Field("count_b", static_cast<std::int64_t>(bd.mCountB));
        rJson.Field("type_mismatch", bd.mTypeMismatch);
        rJson.Field("count_mismatch", bd.mCountMismatch);
        rJson.Field("conn_mismatch_count", static_cast<std::int64_t>(bd.mConnMismatchCount));
        rJson.Key("first_mismatches");
        rJson.BeginArray();
        for (auto f : bd.mFirstMismatches)
            rJson.Int(static_cast<std::int64_t>(f));
        rJson.EndArray();
        rJson.EndObject();
    }
    rJson.EndArray();
    data_diff("point_data", rRep.mPointData);
    data_diff("cell_data", rRep.mCellData);
    data_diff("field_data", rRep.mFieldData);
    rJson.Key("regions");
    rJson.BeginObject();
    cli_names_json(rJson, "only_in_a", rRep.mRegions.mOnlyInA);
    cli_names_json(rJson, "only_in_b", rRep.mRegions.mOnlyInB);
    cli_names_json(rJson, "changed", rRep.mRegions.mChanged);
    rJson.EndObject();
    cli_names_json(rJson, "messages", rRep.mMessages);
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
                                  {"json", {}, false},
                              });
    if (p.positionals.size() != 2)
        throw std::runtime_error("diff requires exactly two mesh files A and B");
    Mesh a = read_mesh_cli(p.positionals[0], opt_value(p, "input-format-a"));
    Mesh b = read_mesh_cli(p.positionals[1], opt_value(p, "input-format-b"));

    meshioplusplus::DiffOptions opts;
    if (p.values.count("atol"))
        opts.atol = meshioplusplus::detail::stod_c(opt_value(p, "atol"));
    if (p.values.count("rtol"))
        opts.rtol = meshioplusplus::detail::stod_c(opt_value(p, "rtol"));
    opts.unordered = has_flag(p, "unordered");
    auto report = meshioplusplus::diff(a, b, opts);
    if (has_flag(p, "json")) {
        const bool equal = report.mVerdict == meshioplusplus::DiffVerdict::Identical ||
                           (report.mVerdict == meshioplusplus::DiffVerdict::EqualWithinTolerance &&
                            !has_flag(p, "exact"));
        meshioplusplus::cli::JsonOut json(std::cout);
        cli_diff_json(json, report);
        // cli_diff_json leaves the object open for the CLI's own key.
        json.Field("equal", equal);
        json.EndObject();
        return equal ? 0 : 1;
    }

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

        // Named regions (doc/regions.md). This used to be the Python shim's
        // job alone -- sets never reached the C++ core, so the native CLI's
        // diff simply omitted them.
        const meshioplusplus::RegionDiff& rd = report.mRegions;
        if (rd.Differs()) {
            std::cout << "regions:\n";
            if (!rd.mOnlyInA.empty())
                std::cout << "  only in A: " << join(rd.mOnlyInA) << "\n";
            if (!rd.mOnlyInB.empty())
                std::cout << "  only in B: " << join(rd.mOnlyInB) << "\n";
            if (!rd.mChanged.empty())
                std::cout << "  differing membership: " << join(rd.mChanged) << "\n";
        }
    }

    bool exact = has_flag(p, "exact");
    bool equal = report.mVerdict == meshioplusplus::DiffVerdict::Identical ||
                 (report.mVerdict == meshioplusplus::DiffVerdict::EqualWithinTolerance && !exact);
    return equal ? 0 : 1;
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
        opts.atol = meshioplusplus::detail::stod_c(opt_value(p, "atol"));
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
        if (cmd == "repair")
            return cmd_repair(rest);
        if (cmd == "shrinkwrap")
            return cmd_shrinkwrap(rest);
        if (cmd == "sobolev-deform")
            return cmd_sobolev_deform(rest);
        if (cmd == "curvature")
            return cmd_curvature(rest);
        if (cmd == "normals")
            return cmd_normals(rest);
        if (cmd == "check")
            return cmd_check(rest);
        if (cmd == "feature-edges")
            return cmd_feature_edges(rest);
        if (cmd == "region-adjacency")
            return cmd_region_adjacency(rest);
        if (cmd == "find-interface")
            return cmd_find_interface(rest);
        if (cmd == "contact-pairs")
            return cmd_contact_pairs(rest);
        if (cmd == "split-interface")
            return cmd_split_interface(rest);
        if (cmd == "hausdorff")
            return cmd_hausdorff(rest);
        if (cmd == "periodic")
            return cmd_periodic(rest);
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
        if (cmd == "slice")
            return cmd_slice(rest);
        if (cmd == "isosurface")
            return cmd_isosurface(rest);
        if (cmd == "voxelize")
            return cmd_voxelize(rest);
        if (cmd == "sdf")
            return cmd_sdf(rest);
        if (cmd == "split")
            return cmd_split(rest);
        if (cmd == "regions")
            return cmd_regions(rest);
        if (cmd == "convert-cells")
            return cmd_convert_cells(rest);
        if (cmd == "subdivide")
            return cmd_subdivide(rest);
        if (cmd == "agglomerate")
            return cmd_agglomerate(rest);
        if (cmd == "resample")
            return cmd_resample(rest);
        if (cmd == "refine")
            return cmd_refine(rest);
        if (cmd == "undo-green")
            return cmd_undo_green(rest);
        if (cmd == "decimate")
            return cmd_decimate(rest);
        if (cmd == "decimate-volume")
            return cmd_decimate_volume(rest);
        if (cmd == "remesh")
            return cmd_remesh(rest);
        if (cmd == "remesh-volume")
            return cmd_remesh_volume(rest);
        if (cmd == "optimize-volume")
            return cmd_optimize_volume(rest);
        if (cmd == "smooth")
            return cmd_smooth(rest);
        if (cmd == "interpolate")
            return cmd_interpolate(rest);
        if (cmd == "conservative-interpolate")
            return cmd_conservative_interpolate(rest);
        if (cmd == "partition")
            return cmd_partition(rest);
        if (cmd == "data")
            return cmd_data(rest);
        if (cmd == "pipeline")
            return cmd_pipeline(rest);
        if (cmd == "stats")
            return cmd_stats(rest);
        if (cmd == "view")
            return cmd_view(rest);
        if (cmd == "screenshot")
            return cmd_screenshot(rest);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }

    std::cerr << "error: unknown command '" << cmd << "'\n\n";
    print_usage(std::cerr);
    return 2;
}
