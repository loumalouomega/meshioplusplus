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
#pragma once

/**
 * @file write_options.hpp
 * @brief `WriteOptions`, the write-side counterpart of `read_options.hpp`.
 *
 * The dispatch registry stores writers as `void(const std::string&, const
 * Mesh&)` lambdas with one set of parameters baked in per format
 * (`registry.cpp`), so nothing that goes through `registry_writers()` can pick
 * ASCII vs binary, a compression codec, or a float format. Callers that need
 * to were expected to call the per-format free function directly -- which the
 * CLI did, with a private variant table, and which the flat bindings (C API,
 * Fortran, WASM) could not do at all.
 *
 * `registry_write_ex()` is that table, hoisted into the core so there is one
 * owner instead of one per consumer. `registry_write(path, mesh, format)` with
 * default options is exactly the registry's own writer, so existing behaviour
 * is unchanged.
 *
 * Not every option applies to every format. The rule is deliberate and matches
 * the CLI's: an option a format cannot honour is an **error**, never silently
 * ignored -- `--codec zstd` on a Gmsh file is a mistake worth reporting, not a
 * no-op. `registry_write_supports()` answers the question up front.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief ASCII/binary selection for formats that have both variants. */
enum class WriteEncoding {
    Default,  ///< the format's registry default (unchanged behaviour)
    Ascii,
    Binary,
};

/**
 * @brief Parameters for `registry_write_ex()`.
 *
 * Default-constructed reproduces `registry_write()` exactly, which is what lets
 * this be added without touching a single existing caller.
 */
struct WriteOptions {
    /** @brief ASCII vs binary. Errors for a format with only one variant. */
    WriteEncoding mEncoding = WriteEncoding::Default;

    /**
     * @brief Block compression codec for the VTK-XML formats (vtu/vtp).
     *
     * `Zlib` is the default those formats already use when binary; `None`
     * disables compression. Errors for any other format -- no other format in
     * meshio++ has a block codec, so naming one is a mistake.
     */
    detail::VtkCodec mCodec = detail::VtkCodec::Zlib;
    bool mCodecSet = false;  ///< whether mCodec was chosen explicitly

    /**
     * @brief `printf`-style float format for ASCII writers that take one
     * (e.g. `".16e"`, the default). Empty keeps the writer's own default.
     */
    std::string mFloatFormat;
};

/**
 * @brief Whether @p rFormat can honour @p rOptions.
 * @param rFormat a resolved format name (see `resolve_format`).
 * @param rOptions the options to check.
 * @param rWhy set to a human-readable reason when the answer is false.
 */
MESHIOPLUSPLUS_API bool registry_write_supports(const std::string& rFormat,
                                                const WriteOptions& rOptions, std::string& rWhy);

/**
 * @brief Write @p rMesh honouring @p rOptions.
 *
 * Resolves the format like `registry_read` (extension default when @p rFormat
 * is empty), then either calls the parameterized per-format writer or, when the
 * options are all defaults, the registry's own writer.
 *
 * @throws WriteError if the format cannot honour an option, or on any writer
 *         failure.
 */
MESHIOPLUSPLUS_API void registry_write_ex(const std::string& rPath, const Mesh& rMesh,
                                          const std::string& rFormat, const WriteOptions& rOptions);

}  // namespace meshioplusplus
