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
 * @file registry.hpp
 * @brief C++-level format dispatch registry shared by the flat bindings
 *        (WASM/JS in `bindings/wasm/`, the C API in `bindings/c/`).
 *
 * The pybind11 binding (`bindings/python/_core.cpp`) exposes one function per
 * format and leaves extension dispatch entirely to Python (`_helpers.py`); the
 * flat bindings instead need a C++-side `format name -> read/write function`
 * table plus an `extension -> default format` map. Those tables originally
 * lived in `bindings/wasm/js_bindings.cpp`; they are hoisted here (compiled into
 * `meshioplusplus_core_obj` via `src/cpp/src/registry.cpp`) so the JS and C
 * bindings share one copy that cannot drift.
 *
 * Parameterized writers get a fixed default here (documented per entry in
 * registry.cpp, matching each format's own Python reference default);
 * per-call overrides are a possible future API addition, deliberately out of
 * scope for v1 of both flat bindings.
 *
 * HDF5-backed formats (cgns, h5m, hmf, med, plus XDMF's HDF data path) and
 * netCDF-backed ones (exodus) are registered only when the corresponding
 * `MESHIOPLUSPLUS_HAS_*` macro is defined -- never under Emscripten, so the
 * WASM format set is unchanged by this refactor. Their extensions are mapped
 * unconditionally so that a build without the dependency reports "format
 * compiled out" (see registry_compiled_out()) instead of the misleading
 * "cannot infer format".
 */

// System includes
#include <functional>
#include <map>
#include <string>
#include <unordered_map>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

using ReadFn = std::function<Mesh(const std::string&)>;
using WriteFn = std::function<void(const std::string&, const Mesh&)>;

/** @brief A reader that can act on `ReadOptions` (selective/partial reads). */
using ReadExFn = std::function<Mesh(const std::string&, const ReadOptions&)>;

/** @brief A reader that can summarize a file without loading its heavy arrays. */
using MetadataFn = std::function<MeshMetadata(const std::string&, const ReadOptions&)>;

/** @brief `format name -> reader` for every format readable in this build. */
MESHIOPLUSPLUS_API const std::map<std::string, ReadFn>& registry_readers();

/** @brief `format name -> writer` for every format writable in this build
 *         (read-only formats like `openfoam` have no entry). */
MESHIOPLUSPLUS_API const std::map<std::string, WriteFn>& registry_writers();

/**
 * @brief `extension (with leading dot) -> default format name`.
 *
 * Ambiguous extensions get this repo's own import-order default (`.msh` ->
 * gmsh, `.inp` -> abaqus); pass an explicit format to select ansys/freefem
 * (.msh) or ansysinp (.inp) instead. Extensions of optional-dependency
 * formats are present even when the format itself is compiled out.
 */
MESHIOPLUSPLUS_API const std::map<std::string, std::string>& registry_extension_defaults();

/**
 * @brief Resolve the effective format: `rFormat` if non-empty, else the
 *        extension default for `rPath`.
 * @throws ReadError if `rFormat` is empty and the extension is unknown.
 */
MESHIOPLUSPLUS_API std::string resolve_format(const std::string& rPath, const std::string& rFormat);

/**
 * @brief The optional dependency a known-but-absent format was compiled out
 *        with, or `nullptr`.
 * @return `"HDF5"` / `"netCDF"` when `rFormat` names a format this build
 *         excluded for lack of that dependency; `nullptr` for formats that
 *         are present or simply unknown. Lets bindings say "format 'med' is
 *         not available in this build (requires HDF5)" instead of "unknown
 *         format".
 */
MESHIOPLUSPLUS_API const char* registry_compiled_out(const std::string& rFormat);

/**
 * @name Selective reads and metadata
 *
 * These tables are deliberately **sparse**: only formats with a native
 * options-aware reader appear in them. `ReadFn`/`registry_readers()`/
 * `registry_writers()` are untouched, so `mio_read`, the WASM `readMesh`, and
 * every existing caller keep working unchanged -- a format's absence here costs
 * a full read, never an error.
 *
 * Prefer the `registry_read`/`registry_read_metadata` free functions over the
 * tables; they own the fallback so each binding surface does not reimplement it.
 * @{
 */

/**
 * @brief `format name -> options-aware reader`, for the formats that have one.
 *
 * Unordered: these are pure keyed lookups, never iterated for output, so they
 * follow the repo default rather than the `std::map` of the older tables above
 * (whose type is kept only because it is baked into the flat bindings).
 */
MESHIOPLUSPLUS_API const std::unordered_map<std::string, ReadExFn>& registry_readers_ex();

/** @brief `format name -> metadata reader`, for the formats that have one. */
MESHIOPLUSPLUS_API const std::unordered_map<std::string, MetadataFn>& registry_metadata_readers();

/**
 * @brief Whether @p rFormat can act on `ReadOptions` rather than being read
 *        whole and filtered afterwards.
 */
MESHIOPLUSPLUS_API bool registry_reader_supports_options(const std::string& rFormat);

/**
 * @brief Read @p rPath honouring @p rOptions, falling back to a full read.
 *
 * Formats without a native options-aware reader are read whole; the result is
 * correct either way, just not faster. Prefer this over indexing the tables.
 * @throws ReadError if the format is unknown or compiled out.
 */
MESHIOPLUSPLUS_API Mesh registry_read(const std::string& rPath, const std::string& rFormat,
                   const ReadOptions& rOptions);

/**
 * @brief Summarize @p rPath without loading its heavy arrays where possible.
 *
 * Formats without a native metadata path are read in full and summarized via
 * `metadata_from_mesh`, with `MeshMetadata::mFellBackToFullRead` set so the
 * caller can tell "fast" from merely "worked".
 * @throws ReadError if the format is unknown or compiled out.
 */
MESHIOPLUSPLUS_API MeshMetadata registry_read_metadata(const std::string& rPath, const std::string& rFormat,
                                    const ReadOptions& rOptions);
/** @} */

}  // namespace meshioplusplus
