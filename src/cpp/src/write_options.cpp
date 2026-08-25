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
 * @file write_options.cpp
 * @brief `registry_write_ex()` -- the single owner of "write this format with
 * these parameters", shared by the CLI and the flat bindings.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/write_options.hpp"

#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/mdpa.hpp"
#include "meshioplusplus/formats/ply.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/vti.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/formats/xdmf.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/detail/provenance.hpp"

namespace meshioplusplus {
namespace {

bool wopt_is_text_only(const std::string& rFormat);

/// Formats with both an ASCII and a binary variant reachable from here.
bool wopt_has_encoding_variants(const std::string& rFormat) {
    return rFormat == "ansys" || rFormat == "flac3d" || rFormat == "gmsh" || rFormat == "ply" ||
           rFormat == "stl" || rFormat == "vtk" || rFormat == "vti" || rFormat == "vtu" ||
           rFormat == "vtp" || rFormat == "xdmf" || wopt_is_text_only(rFormat);
}

/// Text-only formats that still accept an explicit ASCII request.
///
/// "Write this as ASCII" is a perfectly sensible thing to ask of a text format
/// -- it is just the normal write -- so it succeeds, while a BINARY request
/// fails by name. Kept to the set the Python CLI's `ascii` verb also accepts
/// (`_cli/_ascii.py`), so the two CLIs agree on which files the verb handles.
bool wopt_is_text_only(const std::string& rFormat) {
    return rFormat == "mdpa";
}

/// Formats with a VTK-XML block codec.
bool wopt_has_codec(const std::string& rFormat) {
    return rFormat == "vti" || rFormat == "vtu" || rFormat == "vtp";
}

/// Formats whose ASCII writer takes a float format string.
bool wopt_has_float_format(const std::string& rFormat) {
    return rFormat == "flac3d";
}

}  // namespace

bool registry_write_supports(const std::string& rFormat, const WriteOptions& rOptions,
                             std::string& rWhy) {
    if (rOptions.mEncoding != WriteEncoding::Default && !wopt_has_encoding_variants(rFormat)) {
        rWhy = "format '" + rFormat + "' has no ASCII/binary variant to select";
        return false;
    }
    if (rOptions.mCodecSet && !wopt_has_codec(rFormat)) {
        rWhy = "format '" + rFormat + "' has no block compression codec (only vti/vtu/vtp do)";
        return false;
    }
    if (!rOptions.mFloatFormat.empty() && !wopt_has_float_format(rFormat)) {
        rWhy = "format '" + rFormat + "' takes no float-format string";
        return false;
    }
    return true;
}

void registry_write_ex(const std::string& rPath, const Mesh& rMesh, const std::string& rFormat,
                       const WriteOptions& rOptions) {
    const std::string fmt = resolve_format(rPath, rFormat);
    // Bound scope-less notes to this write -- see provenance_begin_write().
    detail::provenance_begin_write();

    // All-defaults goes straight through the registry, so the common path is
    // byte-for-byte what it always was (and picks up formats this file has no
    // parameterized entry for).
    const bool defaulted = rOptions.mEncoding == WriteEncoding::Default && !rOptions.mCodecSet &&
                           rOptions.mFloatFormat.empty();
    if (defaulted) {
        const auto& writers = registry_writers();
        const auto it = writers.find(fmt);
        if (it == writers.end()) {
            if (const char* p_missing = registry_compiled_out(fmt))
                throw WriteError("meshio++: cannot write '" + fmt + "': " + p_missing);
            throw WriteError("meshio++: no writer for format '" + fmt + "'");
        }
        it->second(rPath, rMesh);
        return;
    }

    std::string why;
    if (!registry_write_supports(fmt, rOptions, why))
        throw WriteError("meshio++: " + why);

    const bool binary = rOptions.mEncoding == WriteEncoding::Binary;
    const std::string& r_ff =
        rOptions.mFloatFormat.empty() ? std::string(".16e") : rOptions.mFloatFormat;

    if (fmt == "ansys") {
        write_ansys(rPath, rMesh, binary);
    } else if (fmt == "flac3d") {
        write_flac3d(rPath, rMesh, r_ff, binary);
    } else if (fmt == "gmsh") {
        write_gmsh41(rPath, rMesh, binary);
    } else if (fmt == "ply") {
        write_ply(rPath, rMesh, binary, /*skin=*/true);
    } else if (fmt == "stl") {
        write_stl(rPath, rMesh, binary, /*skin=*/true);
    } else if (fmt == "vtk") {
        write_vtk(rPath, rMesh, binary, /*v51=*/true);
    } else if (fmt == "vti") {
        const detail::VtkCodec codec =
            rOptions.mCodecSet ? rOptions.mCodec
                               : (binary ? detail::VtkCodec::Zlib : detail::VtkCodec::None);
        write_vti_codec(rPath, rMesh, binary, codec);
    } else if (fmt == "vtu") {
        // Compression only exists in the binary encoding; an explicit codec on
        // an ASCII request would otherwise be silently dropped.
        const detail::VtkCodec codec =
            rOptions.mCodecSet ? rOptions.mCodec
                               : (binary ? detail::VtkCodec::Zlib : detail::VtkCodec::None);
        write_vtu_codec(rPath, rMesh, binary, codec);
    } else if (fmt == "vtp") {
        const detail::VtkCodec codec =
            rOptions.mCodecSet ? rOptions.mCodec
                               : (binary ? detail::VtkCodec::Zlib : detail::VtkCodec::None);
        write_vtp_codec(rPath, rMesh, binary, codec);
    } else if (fmt == "xdmf") {
        write_xdmf(rPath, rMesh, binary ? "HDF" : "XML");
    } else if (wopt_is_text_only(fmt)) {
        if (binary)
            throw WriteError("meshio++: format '" + fmt +
                             "' is text-only and has no binary variant");
        registry_write_ex(rPath, rMesh, fmt, WriteOptions{});  // the plain (ASCII) write
    } else {
        // registry_write_supports() vetted the combination, so this is a bug
        // in this file rather than user error.
        throw WriteError("meshio++: internal: no parameterized writer for '" + fmt + "'");
    }
}

}  // namespace meshioplusplus
