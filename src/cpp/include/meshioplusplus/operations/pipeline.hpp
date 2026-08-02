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
 * @file operations/pipeline.hpp
 * @brief The declarative operation pipeline: read a mesh, apply a sequence of
 * mesh/data operations, write the result -- described by a typed `Pipeline`
 * value or, through the JSON front-end, by a `settings.json` document.
 *
 * The engine is split in two layers, and the split is load-bearing:
 *
 *  - The **typed layer** (`PipelineStep`, `apply_pipeline_step`,
 *    `run_pipeline_steps`, `run_pipeline`, `validate_pipeline_step`) is
 *    JSON-independent and always compiles. It is the single owner of the
 *    step vocabulary -- the WASM `convertSurfaceOps` binding dispatches
 *    through it, so the browser viewer's pipeline and a `settings.json`
 *    cannot drift apart.
 *  - The **JSON front-end** (`parse_pipeline_json`/`parse_pipeline_file` and
 *    the `run_pipeline_json`/`run_pipeline_file` conveniences) parses the
 *    settings document with the vendored nlohmann/json (a git submodule at
 *    `src/cpp/third_party/json`, like Eigen). When the submodule is absent
 *    or `-DMESHIOPLUSPLUS_WITH_JSON=OFF`, these entry points still exist and
 *    throw naming the option -- never a link error, never a silent no-op
 *    (the `partition_kahip_parts` contract). `pipeline_has_json()` reports
 *    which build this is.
 *
 * **Vocabulary**: op names and parameter keys are PascalCase (`"Op":
 * "ConvertCells"`, `"RemoveOrphans": true`); enum *values* keep the exact
 * lowercase spellings the `*_from_name` helpers accept (`"simplexify"`,
 * `"redgreen"`, `"cell"`). Refine's comparison key is `Compare`, never `Op` --
 * `Op` is the step discriminant. Parsing is **strict**: an unknown op, an
 * unknown key on a step, or an unknown top-level key is an error naming it,
 * never silently ignored (the `registry_write_supports` rule).
 *
 * **Scope (v1)**: a single-mesh chain. Multi-input (`Merge`, `Interpolate`)
 * and multi-output (`Split`) operations are rejected with an error naming the
 * CLI verb to use instead; `Partition` as a step attaches the `partition:part`
 * labels rather than splitting into pieces (exactly what the browser viewer's
 * partition chip always did).
 *
 * No installed header names an nlohmann type -- the parser is confined to
 * `pipeline.cpp` (the pugixml rule).
 */

// System includes
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <variant>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/write_options.hpp"

namespace meshioplusplus {

/**
 * @brief One parameter value on a pipeline step.
 *
 * The closed set of shapes a step parameter can take: scalars, and homogeneous
 * arrays of them. Deliberately no nested object -- a step is flat, which is
 * what keeps the vocabulary expressible from every binding (a JS object, a
 * settings.json entry, a Python dict of plain values).
 */
using PipelineValue = std::variant<bool, std::int64_t, double, std::string, std::vector<double>,
                                   std::vector<std::int64_t>, std::vector<std::string>>;

/**
 * @brief One operation in the chain: the PascalCase op name plus its
 * parameters, keyed by PascalCase parameter name.
 *
 * `std::map` rather than `unordered_map` deliberately: validation errors and
 * report echoes enumerate the keys, and that order is observable.
 */
struct PipelineStep {
    std::string mOp;
    std::map<std::string, PipelineValue> mParams;
};

/** @brief Where the pipeline reads from, and how (`Input` in settings.json). */
struct PipelineInput {
    std::string mPath;
    /// Empty resolves from the extension, with the `sniff_format` fallback.
    std::string mFormat;
    ReadOptions mOptions;
};

/** @brief Where the pipeline writes to, and how (`Output` in settings.json). */
struct PipelineOutput {
    std::string mPath;
    /// Empty resolves from the extension.
    std::string mFormat;
    /// Honoured through `registry_write_ex`: an option the format cannot
    /// honour is an error, never silently ignored.
    WriteOptions mOptions;
};

/** @brief A whole settings document: input, the operation chain, output. */
struct Pipeline {
    int mVersion = 1;
    PipelineInput mInput;
    std::vector<PipelineStep> mSteps;
    PipelineOutput mOutput;
};

/**
 * @brief What one executed step reports back: its op name and the counters the
 * underlying operation returned (PascalCase keys, e.g. `PointsWelded`).
 *
 * Counters are doubles because every consumer of this shape (JS, JSON, Python)
 * speaks doubles; the underlying `std::int64_t` counters are far below 2^53.
 */
struct PipelineStepReport {
    std::string mOp;
    std::vector<std::pair<std::string, double>> mCounters;
};

/** @brief The report of a whole run: one entry per step, plus warnings. */
struct PipelineReport {
    std::vector<PipelineStepReport> mSteps;
    std::vector<std::string> mWarnings;
};

/**
 * @brief The step vocabulary itself: op name -> its parameter keys, in table
 * order (aliases like `Section` included).
 *
 * Exposed so the pure-Python runner's transcribed table can be pinned against
 * this one by a test (`_core.pipeline_op_table`, the `refine_mask_table` /
 * `colormap_table` precedent) instead of drifting silently.
 */
MESHIOPLUSPLUS_API std::vector<std::pair<std::string, std::vector<std::string>>>
pipeline_op_table();

/**
 * @brief Checks @p rStep against the step vocabulary.
 *
 * Strict: an unknown op throws listing the known ops; a key the op does not
 * take throws naming the key and listing the op's keys. The excluded
 * multi-mesh ops (`Merge`, `Interpolate`, `Split`, `Diff`) throw an error
 * naming the CLI verb to use instead.
 *
 * @throws std::invalid_argument as above.
 */
MESHIOPLUSPLUS_API void validate_pipeline_step(const PipelineStep& rStep);

/**
 * @brief Applies one step to @p rMesh, appending its entry to @p rReport.
 *
 * The single owner of the step dispatch -- the WASM `convertSurfaceOps`
 * binding and the JSON pipeline both go through here. Parameter defaults are
 * the ones the browser viewer's pipeline always used (e.g. `Clean` welds only
 * on request, `Decimate` with no criterion keeps half the faces).
 *
 * Takes the mesh by value (move it in): most ops return a new mesh anyway,
 * and the KRATOS backend's `Mesh` is not copy-constructible.
 *
 * @throws std::invalid_argument on an unknown op or a mis-typed parameter,
 *         plus whatever the underlying operation throws.
 */
MESHIOPLUSPLUS_API Mesh apply_pipeline_step(Mesh mesh, const PipelineStep& rStep,
                                            PipelineReport& rReport);

/**
 * @brief Applies @p rSteps in order (validating each first), moving the mesh
 * through the chain.
 */
MESHIOPLUSPLUS_API Mesh run_pipeline_steps(Mesh mesh, const std::vector<PipelineStep>& rSteps,
                                           PipelineReport& rReport);

/**
 * @brief Runs a whole pipeline: `registry_read` the input, apply the steps,
 * `registry_write_ex` the output.
 *
 * Every step is validated before the input is read, so a typo fails fast
 * rather than after an expensive read.
 */
MESHIOPLUSPLUS_API PipelineReport run_pipeline(const Pipeline& rPipeline);

/**
 * @brief Parses a settings.json document (see the file comment for the
 * schema).
 * @throws std::runtime_error naming `-DMESHIOPLUSPLUS_WITH_JSON=ON` when the
 *         parser is compiled out; std::invalid_argument on any schema error.
 */
MESHIOPLUSPLUS_API Pipeline parse_pipeline_json(const std::string& rText);

/** @brief `parse_pipeline_json` over the contents of @p rPath. */
MESHIOPLUSPLUS_API Pipeline parse_pipeline_file(const std::string& rPath);

/** @brief `parse_pipeline_json` + `run_pipeline`. */
MESHIOPLUSPLUS_API PipelineReport run_pipeline_json(const std::string& rText);

/** @brief `parse_pipeline_file` + `run_pipeline`. */
MESHIOPLUSPLUS_API PipelineReport run_pipeline_file(const std::string& rPath);

/** @brief Whether this build carries the JSON front-end. */
MESHIOPLUSPLUS_API bool pipeline_has_json() noexcept;

// --------------------------------------------------------------------------
// The Output/Input option spellings, owned here so every front-end (the JSON
// parser, the WASM settings-object converter, the CLIs) accepts exactly the
// same words. All always compiled (JSON-independent).
// --------------------------------------------------------------------------

/** @brief `""`/`"default"`/`"ascii"`/`"binary"`; anything else throws. */
MESHIOPLUSPLUS_API WriteEncoding pipeline_encoding_from_name(const std::string& rName);

/** @brief `"none"`/`"zlib"`/`"lz4"`/`"zstd"`; anything else (incl. empty) throws. */
MESHIOPLUSPLUS_API detail::VtkCodec pipeline_codec_from_name(const std::string& rName);

/** @brief `""`/`"auto"`/`"on"`/`"off"`; anything else throws. */
MESHIOPLUSPLUS_API MmapMode pipeline_mmap_from_name(const std::string& rName);

}  // namespace meshioplusplus
