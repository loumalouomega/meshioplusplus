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
 * @file operations/sequence.hpp
 * @brief Multi-file / transient datasets: treat a *set* of files (or the steps
 * inside one multi-step file) as one ordered logical sequence.
 *
 * This is how transient solver output actually arrives -- `out_0000.vtu …
 * out_0500.vtu` -- and how most of the 42 formats have to express time, since
 * only a minority carry several steps natively.
 *
 * **It is a driver, not a new mesh operation.** Everything here reads and
 * writes through the existing registry (`registry_read`, `registry_write_ex`)
 * and runs operation chains through the existing typed pipeline layer
 * (`run_pipeline_steps`). The typed layer stays the single owner of the step
 * dispatch -- there is deliberately no second `if (op == ...)` chain in this
 * file, so a sequence document and the browser viewer's `convertSurfaceOps`
 * still cannot drift apart.
 *
 * Three shapes, distinguished by how many files the input expands to and
 * whether the output path carries a `{step}`/`{index}` token:
 *
 *  - **FanOut** — one multi-step file -> `out_{step}.vtu`. The answer for every
 *    format that cannot express time, which is most of them.
 *  - **FanIn** — N single-step files -> one multi-step file (XDMF today).
 *  - **Sequence** — N files -> N files, the operation chain applied per step.
 *
 * @section sequence_ordering Ordering
 *
 * A sequence has a defined order and it is **not lexicographic**: `out_10.vtu`
 * must follow `out_9.vtu`. `sequence_natural_less` implements natural-numeric
 * ordering; see its documentation for the exact rule, which is a documented
 * contract and not an implementation detail.
 *
 * @section sequence_streaming The streaming invariant
 *
 * **At most one `Mesh` is alive at any point** inside `sequence_to_timeseries`,
 * `timeseries_to_sequence` and `run_sequence_pipeline`, and `sequence_read_step`
 * hands back exactly one. This is a contract, not an optimization -- the whole
 * feature exists so that a 500-step dataset is traversable on a laptop. No
 * implementation may buffer the sequence, including "just for sorting" or "just
 * to compute the time range": `sequence_expand` returns the *plan* (paths, step
 * indices and times), never meshes. `tests/cpp/test_sequence.cpp` pins this
 * through the `BufferAllocator` hook, asserting the peak is O(1) in the step
 * count rather than merely small.
 *
 * @section sequence_json The JSON front-end
 *
 * `parse_sequence_json`/`_file` and `run_sequence_json`/`_file` extend the
 * `settings.json` schema with `Mode` and `Input.Pattern`/`Paths`/`Times`/
 * `TimeFrom`. They live behind the same guard as the pipeline's own parser:
 * with `-DMESHIOPLUSPLUS_WITH_JSON=OFF` they exist and throw naming the option,
 * never a link error and never a silent no-op. No installed header names an
 * nlohmann type -- the parser is confined to `pipeline.cpp` (the pugixml rule).
 */

// System includes
#include <cstddef>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/pipeline.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/write_options.hpp"

namespace meshioplusplus {

/**
 * @brief Where an entry's time value came from.
 *
 * Reported per entry rather than inferred by the caller: "this time is the
 * integer index because the file said nothing" and "this time is 0.25 because
 * the file said so" are different facts, and a user debugging a plot needs to
 * know which they have.
 */
enum class SequenceTimeSource {
    Explicit,  ///< From the caller's own list (`SequenceInput::mTimes`).
    File,      ///< From the file: a step's time value, or `field_data["meshio:time"]`.
    Filename,  ///< Parsed from the last digit run of the file's stem.
    Index,     ///< Nothing said; the integer position in the sequence.
};

/// Which time source to use. `Auto` walks the documented precedence.
enum class SequenceTimeFrom {
    Auto,      ///< Explicit -> File -> Filename -> Index (the default).
    File,      ///< Only the file; falls back to the index if it says nothing.
    Filename,  ///< Only the filename; falls back to the index.
    Index,     ///< Always the integer position.
};

/// Which of the three sequence shapes a run takes.
enum class SequenceMode {
    Auto,      ///< Infer from the expanded input and the output pattern (the default).
    Sequence,  ///< N inputs -> N outputs.
    FanIn,     ///< N inputs -> one multi-step output.
    FanOut,    ///< One multi-step input -> N outputs.
};

/**
 * @brief One logical step of a sequence: a file, plus which step *inside* it.
 *
 * `mStep` is always 0 for a single-step file, so a directory of `.vtu`s and one
 * multi-step `.xdmf` produce the same shape and the drivers need no special
 * case for either.
 */
struct SequenceEntry {
    std::string mPath;
    std::size_t mStep = 0;  ///< Step index within `mPath`; 0 for a single-step file.
    double mTime = 0.0;
    SequenceTimeSource mTimeSource = SequenceTimeSource::Index;
};

/**
 * @brief Where a sequence reads from, and how (`Input` in settings.json).
 *
 * Exactly one of `mPaths` and `mPattern` is set; both empty, or both set, is an
 * error naming the conflict.
 */
struct SequenceInput {
    /// An explicit, ordered list. Not re-sorted unless `mSortExplicit`, because
    /// a caller-supplied list is a stated order.
    std::vector<std::string> mPaths;
    /// A glob. See `sequence_glob_match` for the (deliberately narrow) pattern
    /// language; the directory component is taken literally. Always sorted
    /// natural-numerically, since a directory listing has no meaningful order.
    std::string mPattern;
    /// Empty resolves per file from the extension, with the `sniff_format`
    /// fallback -- the same rule `run_pipeline` uses for its single input.
    std::string mFormat;
    ReadOptions mOptions;
    /// Explicit per-entry times. When non-empty its size must equal the entry
    /// count, or expansion throws naming both counts.
    std::vector<double> mTimes;
    SequenceTimeFrom mTimeFrom = SequenceTimeFrom::Auto;
    /// Sort `mPaths` natural-numerically too. Ignored for `mPattern`, which is
    /// always sorted.
    bool mSortExplicit = false;
};

/** @brief Where a sequence writes to, and how (`Output` in settings.json). */
struct SequenceOutput {
    /// Either a literal path (one multi-step file) or a pattern carrying
    /// `{step}` / `{index}` (one file per step). See `sequence_expand_pattern`.
    std::string mPath;
    /// Empty resolves from the extension.
    std::string mFormat;
    /// Honoured through `registry_write_ex`: an option the format cannot
    /// honour is an error, never silently ignored.
    WriteOptions mOptions;
};

/** @brief A whole sequence settings document. */
struct SequencePipeline {
    int mVersion = 1;
    SequenceInput mInput;
    /// Applied to **every** step, through `run_pipeline_steps`.
    std::vector<PipelineStep> mSteps;
    SequenceOutput mOutput;
    /// `Auto` infers; anything else *asserts* the inference and errors naming
    /// both on a mismatch (see `run_sequence_pipeline`).
    SequenceMode mMode = SequenceMode::Auto;
    /// Run the steps over several files at once. **A Python-driver feature**:
    /// the C++ engine carries the request so a settings document round-trips,
    /// validates it, and then runs serially with a warning -- every operation
    /// already parallelizes internally, so a step-level parallel region nested
    /// over them would oversubscribe on every backend.
    bool mParallel = false;
    /// Worker count for `mParallel`; 0 means "as many as there are cores".
    /// Ignored, with `mParallel`, by the C++ engine.
    int mWorkers = 0;
};

/**
 * @brief The `field_data` key carrying a single mesh's own time value.
 *
 * A length-1 Float64 array, generalizing the `exodus:time` convention the
 * Exodus reader/writer already round-trips ("the one value this mesh is a
 * snapshot at"). Colon-namespaced like `partition:part` and `iso:value`,
 * because a bare `"time"` would collide with a solver's own field far too often.
 *
 * `timeseries_to_sequence` attaches it to every file it writes. **It only
 * survives where the target format carries `field_data` at all**, which today
 * means Exodus, Gmsh and MDPA -- VTU, VTK and most others carry none in either
 * direction. So a fan-out to `out_{step}.vtu` followed by a fan-in recovers the
 * step *index* from the filename, not the original time value, and
 * `SequenceEntry::mTimeSource` reports `Filename` so a caller can see that.
 * Supply `SequenceInput::mTimes` on the way back in to close the round trip on
 * time for such a format; a fan-in to XDMF needs no help, since the series
 * records its own step times.
 */
inline constexpr const char* kSequenceTimeKey = "meshio:time";

// --------------------------------------------------------------------------
// Pure units. No I/O, no registry, no mesh -- these are the pieces the whole
// feature's correctness rests on, and they are testable in isolation.
// --------------------------------------------------------------------------

/**
 * @brief Natural-numeric ordering, so `out_9.vtu` sorts before `out_10.vtu`.
 *
 * The rule, which is a documented contract:
 *
 *  1. Both strings are split into maximal runs of digits and maximal runs of
 *     non-digits, and compared run by run.
 *  2. Two non-digit runs compare byte by byte as `unsigned char` (plain `char`
 *     signedness is implementation-defined, and would order the same UTF-8
 *     paths differently on ARM than on x86).
 *  3. Two digit runs compare *numerically*: leading zeros are stripped, a
 *     shorter stripped run is less, and equal lengths compare
 *     lexicographically. This is done on the digits themselves and never
 *     through `stoull`, so a 40-digit hash-named file cannot overflow.
 *  4. A digit run at the same position as a non-digit run sorts first.
 *  5. If every run compares equal, the *unstripped* strings are compared
 *     byte-wise, so `out_1` < `out_01` deterministically.
 *  6. The whole path is compared, not the basename, so `a/out_9` < `a/out_10`
 *     < `b/out_1`.
 *
 * Rule 5 is what makes this a **strict weak ordering**: without it `out_1` and
 * `out_01` are mutually "not less" yet not equivalent, and `std::sort` would be
 * undefined behaviour on a directory mixing padded and unpadded names.
 * `tests/cpp/test_sequence.cpp` brute-forces the four axioms over a table that
 * includes exactly that case.
 */
MESHIOPLUSPLUS_API bool sequence_natural_less(const std::string& rA, const std::string& rB);

/**
 * @brief Glob-match @p rName against @p rPattern.
 *
 * The meshio++ pattern language is exactly `*` (any run, possibly empty) and
 * `?` (exactly one character) -- **no** `**`, no `[set]`, no brace expansion,
 * and no special treatment of a leading dot. It is deliberately narrower than
 * POSIX `glob(3)` and than Python's `glob`/`fnmatch`, so that the C++ matcher
 * and its Python twin cannot accept different things; `[abc]` matches those
 * three characters literally rather than as a set.
 *
 * Matching is iterative with backtracking (no recursion, no allocation), so a
 * pathological `a*a*a*a*b` cannot blow the stack.
 */
MESHIOPLUSPLUS_API bool sequence_glob_match(const std::string& rPattern, const std::string& rName);

/**
 * @brief Expand `{step}` / `{index}` in an output pattern.
 *
 * `{index}` is the plain decimal index; `{step}` is that index zero-padded to
 * `max(4, digits(Count - 1))` -- so a 12-step run writes `out_0000 …
 * out_0011`, and a 20000-step run widens to `out_00000 …` rather than sorting
 * wrongly in a directory listing.
 *
 * Substring replace-all, **not** `std::format`/`str.format` semantics: an
 * unrelated `{` in the path is a literal here. This deliberately matches the
 * native CLI's existing `replace_key`/`partition_replace_part` helpers, and
 * deliberately differs from the Python CLI's `str.format`-based `{key}`/`{part}`
 * expansion, which raises on a stray brace. The asymmetry is pre-existing;
 * these tokens inherit it rather than introducing a third convention.
 *
 * @param rPattern the output pattern; must contain `{step}` or `{index}`.
 * @param Index the 0-based position in the sequence.
 * @param Count the total number of entries, used only for the padding width.
 */
MESHIOPLUSPLUS_API std::string sequence_expand_pattern(const std::string& rPattern,
                                                       std::size_t Index, std::size_t Count);

/** @brief Whether @p rPath contains a `{step}` or `{index}` token. */
MESHIOPLUSPLUS_API bool sequence_pattern_has_token(const std::string& rPath);

// --------------------------------------------------------------------------
// Registry-derived capability queries.
// --------------------------------------------------------------------------

/**
 * @brief How many time steps @p rPath carries.
 *
 * Derived from the registry rather than a hardcoded per-format table:
 * `registry_read_metadata(...).mTimeValues.size()`, clamped to at least 1 so
 * that a single-step file is always a valid degenerate sequence. A format whose
 * metadata reader does not fill `mTimeValues` therefore reports 1, which is the
 * truthful answer for every format that cannot express time.
 *
 * Today that means XDMF and Exodus report real counts. MED honours
 * `ReadOptions::mTimeStep` but has no metadata reader, so a multi-step `.med`
 * reports 1; that is a recorded gap in MED's metadata support and not a special
 * case here -- the moment `read_med_metadata` fills `mTimeValues`, MED fan-out
 * starts working with no change to this file.
 *
 * Never throws for an unreadable file: an unreadable path reports 1 and the
 * failure surfaces from the actual read, with its own diagnostics.
 */
MESHIOPLUSPLUS_API std::size_t sequence_num_steps(const std::string& rPath,
                                                  const std::string& rFormat);

/**
 * @brief Whether @p rFormat's writer can emit a **multi-step** file.
 *
 * There is no file to probe on the write side, so unlike `sequence_num_steps`
 * this is a predicate in the exact style of `registry_write_supports`: a small
 * owned set, with @p rWhy filled with a message naming the format and the
 * remedy on failure.
 *
 * The anti-drift mechanism is a test rather than the table: a gtest iterates
 * `registry_writers()` and asserts this predicate agrees with whether a real
 * two-step `sequence_to_timeseries` to that format actually succeeds, so a
 * format that grows a series writer without updating the predicate turns CI red
 * naming itself.
 */
MESHIOPLUSPLUS_API bool sequence_write_supports_time(const std::string& rFormat, std::string& rWhy);

// --------------------------------------------------------------------------
// The driver.
// --------------------------------------------------------------------------

/**
 * @brief Expand an input to its ordered entry list, with each entry's time and
 * the source that time came from.
 *
 * Reads no heavy data: at most one `registry_read_metadata` per file, and only
 * when a time actually has to be derived from a multi-step file. A single-step
 * file whose time would come from `field_data` is left provisionally on the
 * filename or index source; the drivers upgrade it to `File` when the mesh they
 * were going to read anyway turns out to carry `meshio:time`, so no file is
 * ever read twice.
 *
 * Time-value precedence under `SequenceTimeFrom::Auto`:
 *
 *  1. **Explicit** — `SequenceInput::mTimes`.
 *  2. **File** — a multi-step file's own step time, or a single-step file's
 *     `field_data["meshio:time"]` (or a length-1 `field_data["exodus:time"]`).
 *  3. **Filename** — the last maximal digit run in the stem, so
 *     `run17/out_0042.vtu` gives 42.
 *  4. **Index** — the integer position, recorded so the caller can see it was
 *     a fallback.
 *
 * @throws ReadError if a pattern matches nothing (never a silently empty
 *         sequence), or if a listed path does not exist.
 * @throws std::invalid_argument if both/neither of `mPaths`/`mPattern` are set,
 *         or `mTimes` has the wrong length (naming both counts).
 */
MESHIOPLUSPLUS_API std::vector<SequenceEntry> sequence_expand(const SequenceInput& rInput);

/**
 * @brief Read exactly one entry.
 *
 * This is the whole lazy read surface: a caller loops over the entries and
 * never holds two meshes. See @ref sequence_streaming.
 *
 * @param rEntries the plan from `sequence_expand`.
 * @param Index which entry to read.
 * @param rFormat empty resolves per file, as `SequenceInput::mFormat` does.
 * @param rOptions read options; the entry's own step index overrides
 *        `ReadOptions::mTimeStep`, which is what makes fan-out work.
 * @throws std::out_of_range if @p Index is past the end.
 */
MESHIOPLUSPLUS_API Mesh sequence_read_step(const std::vector<SequenceEntry>& rEntries,
                                           std::size_t Index, const std::string& rFormat,
                                           const ReadOptions& rOptions);

/**
 * @brief Fan-in: N entries -> one multi-step file, streaming.
 *
 * Reads step i, writes step i, releases it. At most one `Mesh` is alive.
 *
 * @throws WriteError naming the format when the target cannot hold a series
 *         (see `sequence_write_supports_time`) -- never a silent truncation to
 *         the first step.
 */
MESHIOPLUSPLUS_API void sequence_to_timeseries(const SequenceInput& rInput,
                                               const SequenceOutput& rOutput);

/**
 * @brief Fan-out: one multi-step file -> `out_{step}.vtu`, streaming.
 *
 * Each written mesh carries `field_data["meshio:time"]`, so the round trip back
 * through `sequence_to_timeseries` recovers the times as well as the geometry.
 *
 * @throws std::invalid_argument if `rOutput.mPath` carries no `{step}`/`{index}`.
 */
MESHIOPLUSPLUS_API void timeseries_to_sequence(const std::string& rInPath,
                                               const std::string& rInFormat,
                                               const ReadOptions& rOptions,
                                               const SequenceOutput& rOutput);

/**
 * @brief Whether a document that used no sequence-only key must still be run by
 * the sequence driver.
 *
 * True when the output is a `{step}`/`{index}` pattern, or when the input turns
 * out to carry **several steps** — because writing the first one and calling it
 * a conversion is exactly the silent truncation this layer exists to prevent.
 * An explicit `ReadOptions::mTimeStep` **is** a deliberate single-step
 * selection and opts out.
 *
 * The step probe is gated on the formats that can carry time at all, so it
 * costs nothing for the 39 that cannot (`read_metadata` on a format with no
 * native metadata reader is a full read). Shared by every front-end — the JSON
 * one, the WASM settings converter, both CLIs and the Python twin — so they
 * cannot disagree about which documents are transient.
 */
MESHIOPLUSPLUS_API bool sequence_input_needs_driver(const SequenceInput& rInput,
                                                    const SequenceOutput& rOutput);

/**
 * @brief Resolve @p rPipeline's mode: infer it, then check a stated one agrees.
 *
 * The inference, from the expanded entry count and whether the output carries a
 * token:
 *
 * | entries | output has `{step}`/`{index}` | mode |
 * |---|---|---|
 * | 1 | no | `Sequence` (degenerate; identical to a plain single-file run) |
 * | >1 (one file, several steps) | yes | `FanOut` |
 * | >1 (one file, several steps) | no | **error** — refuses to truncate |
 * | >1 files | no | `FanIn` |
 * | >1 files | yes | `Sequence` |
 * | 1 | yes | `Sequence` with one entry |
 *
 * A non-`Auto` `mMode` never *changes* the run; it asserts the inference and
 * throws naming both on a mismatch. That matters because the inference depends
 * on how many files a glob happened to match: a pattern that matched exactly
 * one file would otherwise quietly take the single-file path.
 *
 * @throws std::invalid_argument on a mismatch, or when a multi-step input is
 *         aimed at a single-step output.
 */
MESHIOPLUSPLUS_API SequenceMode sequence_resolve_mode(const std::vector<SequenceEntry>& rEntries,
                                                      const SequenceOutput& rOutput,
                                                      SequenceMode Stated);

/**
 * @brief Run a whole sequence: expand, then per step read -> `run_pipeline_steps`
 * -> write, streaming throughout.
 *
 * The operation chain is validated once, before anything is read. The report
 * carries one entry per (step, op), each op's name suffixed with nothing --
 * the step is identifiable by position, and a `Step` counter is attached to
 * every entry so a 500-step run stays machine-readable.
 */
MESHIOPLUSPLUS_API PipelineReport run_sequence_pipeline(const SequencePipeline& rPipeline);

/**
 * @brief Parse a name into a `SequenceMode`: `"auto"`, `"sequence"`, `"fan-in"`
 * or `"fan-out"`.
 *
 * Lowercase and hyphenated, following the repo-wide rule that settings *keys*
 * are PascalCase but enum *values* keep their `*_from_name` spelling -- as in
 * gradient's `"green-gauss"` and convert_cells' `"simplexify"`.
 */
MESHIOPLUSPLUS_API SequenceMode sequence_mode_from_name(const std::string& rName);

/** @brief The canonical name of a `SequenceMode`, for errors and reports. */
MESHIOPLUSPLUS_API const char* sequence_mode_name(SequenceMode Mode);

/** @brief Parse a name into a `SequenceTimeFrom`: `"auto"`, `"file"`, `"filename"`, `"index"`. */
MESHIOPLUSPLUS_API SequenceTimeFrom sequence_time_from_name(const std::string& rName);

/** @brief The name of a `SequenceTimeSource`, for reports (`"explicit"`, ...). */
MESHIOPLUSPLUS_API const char* sequence_time_source_name(SequenceTimeSource Source);

// --------------------------------------------------------------------------
// The JSON front-end. Same guard as the pipeline's own parser: these exist in
// every build and throw naming -DMESHIOPLUSPLUS_WITH_JSON=ON when compiled out.
// --------------------------------------------------------------------------

/**
 * @brief Parse a sequence settings document.
 *
 * The v9.11.0 schema plus `Mode`, `Parallel`, `Workers` at the top level and
 * `Pattern`, `Paths`, `Times`, `TimeFrom` under `Input`. A document using none
 * of them parses to a one-entry sequence and behaves exactly as
 * `parse_pipeline_json` + `run_pipeline` always did.
 *
 * `Parallel`/`Workers` are accepted and validated here but are a **Python
 * driver** feature: the C++ engine runs steps serially, because every operation
 * already parallelizes internally and a step-level parallel region nested over
 * them would oversubscribe. `run_sequence_pipeline` records a warning when a
 * document asks for them.
 *
 * @throws std::runtime_error naming `-DMESHIOPLUSPLUS_WITH_JSON=ON` when the
 *         parser is compiled out; std::invalid_argument on any schema error.
 */
MESHIOPLUSPLUS_API SequencePipeline parse_sequence_json(const std::string& rText);

/** @brief `parse_sequence_json` over the contents of @p rPath. */
MESHIOPLUSPLUS_API SequencePipeline parse_sequence_file(const std::string& rPath);

/** @brief `parse_sequence_json` + `run_sequence_pipeline`. */
MESHIOPLUSPLUS_API PipelineReport run_sequence_json(const std::string& rText);

/** @brief `parse_sequence_file` + `run_sequence_pipeline`. */
MESHIOPLUSPLUS_API PipelineReport run_sequence_file(const std::string& rPath);

}  // namespace meshioplusplus
