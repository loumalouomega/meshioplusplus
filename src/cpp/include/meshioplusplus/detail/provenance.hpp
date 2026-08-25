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
 * @file detail/provenance.hpp
 * @brief The provenance record every writer renders into its own header slot.
 *
 * ### Two layers
 *
 * `kProvenanceTag` (v10.15.0) is the unconditional one-line credit: it is
 * always emitted, by both engines, character-identically. Nothing about it
 * changed here.
 *
 * The `Record`/`Scope` surface (v10.16.0) is the **opt-in** richer block --
 * source, target, the operation chain and the conversion assumptions accepted
 * on the way. `Mode::Off` is the default, so a caller who asks for nothing
 * gets byte-for-byte what v10.15.0 wrote.
 *
 * ### Why a scoped context rather than `WriteOptions`
 *
 * `registry_write_ex` is the documented "single owner of write-this-format-
 * with-these-parameters", so `WriteOptions` looks like the obvious home. It is
 * not: only four write paths reach it (the native CLI, `mio_write_ex`, the
 * pipeline and the sequence driver). The Python path
 * (`_helpers.write` -> the format shim -> `_core.<fmt>_write` -> the per-format
 * free function), the WASM path (the `registry_writers()` lambdas) and plain
 * `mio_write` all bypass it entirely -- so a `WriteOptions`-borne record would
 * be invisible from the primary user surface. Growing `WriteOptions` is also a
 * Tier A ABI change, and growing `mio_write_opts` only preserves `sizeof`
 * where `sizeof(void*) == 8`.
 *
 * A context solves all of that at once, because writers read it exactly where
 * they already read `kProvenanceTag` -- inside their own bodies -- so no
 * signature anywhere has to change.
 *
 * **Thread-local, not process-global.** This is the `set_buffer_allocator`
 * shape (`ndarray.hpp`) with one deliberate difference: writes are not
 * serialized against each other, so a process-global record would let one
 * write's provenance leak into a concurrent write on another thread. The
 * scope is RAII and restores the previous record on destruction, so nesting
 * works and an exception cannot strand it.
 *
 * ### Slot tiers
 *
 * Formats differ in what their header slot can physically hold, and the
 * difference is not a preference -- a binary STL header is 80 bytes and the
 * VTK legacy title is one line. `SlotTier` classifies it, and `render()`
 * degrades accordingly. `Mode::Required` turns that degradation into an error
 * instead, which is how `write_options.hpp`'s standing rule ("an option a
 * format cannot honour is an error, never silently ignored") is honoured
 * without making a provenance request fail on a third of the format table by
 * default.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/version.hpp"

namespace meshioplusplus {
namespace detail {

/// The canonical one-line provenance tag every writer emits, wrapped in each
/// format's own comment syntax at each writer's existing header position.
inline constexpr const char* kProvenanceTag = "Written by meshio++ v" MESHIOPLUSPLUS_VERSION_STRING;

/// How much provenance a writer should render.
enum class ProvenanceMode : std::uint8_t {
    /// Only `kProvenanceTag` -- exactly what v10.15.0 wrote.
    Off = 0,
    /// The full block where the slot allows it, degrading silently otherwise.
    /// The process default -- see `default_provenance_mode()`.
    BestEffort = 1,
    /// The full block, or a `WriteError` naming the format whose slot cannot
    /// hold one.
    Required = 2,
};

/// What a format's header slot can physically hold.
enum class SlotTier : std::uint8_t {
    /// No free-text slot at all (`wkt`, `ugrid`, the HDF5 containers, ...).
    None = 0,
    /// One line, hard-capped in bytes (binary STL's 80, EnSight's `str80`).
    Bounded = 1,
    /// Exactly one line (the VTK legacy title, Tecplot's `TITLE=`, ...).
    SingleLine = 2,
    /// Arbitrarily many lines (every `#`/`!`/`*`/`$`/`%`/`//` comment syntax,
    /// PLY `comment`, XML `<!--`, Abaqus `*HEADING`, the OpenFOAM banner).
    Block = 3,
};

/// One conversion assumption accepted on the way to the output file.
///
/// `mCategory` is a short stable slug (`"regions-dropped"`, `"node-order"`,
/// `"dtype"`, ...) so a consumer can filter without parsing prose;
/// `mDetail` is the human-readable specifics.
struct ProvenanceNote {
    std::string mCategory;
    std::string mDetail;
};

/**
 * @brief Everything the record can carry.
 *
 * Deliberately **not** carrying which engine rendered it. The roadmap asks
 * for that diagnostic, but it contradicts the harder guarantee that the C++
 * core and its numpy twins emit character-identical bytes -- an engine marker
 * is precisely the `(C++ core)`-vs-`v{version}` drift v10.15.0 removed. It is
 * reported through `current_provenance()` instead, where it costs nothing.
 */
struct ProvenanceRecord {
    /// Where the mesh was read from, when a single call spans read and write.
    std::string mSourcePath;
    /// The source format, with its detected sub-version where one exists
    /// (`"gmsh 4.1"`, `"med 4"`).
    std::string mSourceFormat;

    /// The resolved target format.
    std::string mTargetFormat;
    /// The encoding actually used, not the one requested (`"ascii"`/`"binary"`).
    std::string mEncoding;
    /// The block codec actually used, where the format has one.
    std::string mCodec;
    /// The float format actually used, where the writer takes one.
    std::string mFloatFormat;

    /// The operation chain, in order, each rendered as `op(key=value, ...)`.
    std::vector<std::string> mOperations;

    /// The conversion assumptions accepted, in the order they were made.
    std::vector<ProvenanceNote> mNotes;

    /// ISO-8601 UTC to the second, or empty when suppressed.
    std::string mTimestamp;
};

/**
 * @brief Installs @p rRecord as the active record for the current thread.
 *
 * RAII: the previous record (if any) is restored on destruction, so scopes
 * nest and an exception in the middle of a write cannot strand one. Move-only
 * -- a copied scope would restore twice.
 */
class MESHIOPLUSPLUS_API ProvenanceScope {
public:
    ProvenanceScope(ProvenanceMode mode, ProvenanceRecord record);
    explicit ProvenanceScope(ProvenanceMode mode) : ProvenanceScope(mode, ProvenanceRecord{}) {}
    ~ProvenanceScope();

    ProvenanceScope(const ProvenanceScope&) = delete;
    ProvenanceScope& operator=(const ProvenanceScope&) = delete;
    ProvenanceScope(ProvenanceScope&&) = delete;
    ProvenanceScope& operator=(ProvenanceScope&&) = delete;

    /// The record as it stands now, including every `provenance_note()` made
    /// since the scope opened. This is where a caller reads back what was
    /// recorded -- including facts (like which engine ran) deliberately kept
    /// out of the file.
    const ProvenanceRecord& Get() const;

private:
    bool mHadPrevious;
    ProvenanceRecord mPrevious;
    ProvenanceMode mPreviousMode;
};

/// The active mode for this thread (`Off` when no scope is open).
MESHIOPLUSPLUS_API ProvenanceMode current_provenance_mode();

/// The active record for this thread; empty when no scope is open.
MESHIOPLUSPLUS_API const ProvenanceRecord& current_provenance();

/**
 * @brief Records one conversion assumption against the active scope.
 *
 * A **no-op outside a scope**, which is what lets it be called
 * unconditionally from a writer or an operation without either of them
 * knowing whether anyone is listening. Duplicate `(category, detail)` pairs
 * are collapsed: a per-cell warning must not produce a per-cell record.
 */
MESHIOPLUSPLUS_API void provenance_note(std::string_view category, std::string_view detail);

/// Sets a field on the active record. No-ops outside a scope.
///@{
MESHIOPLUSPLUS_API void provenance_set_source(std::string_view path, std::string_view format);
MESHIOPLUSPLUS_API void provenance_set_target(std::string_view format, std::string_view encoding,
                                              std::string_view codec,
                                              std::string_view float_format);
MESHIOPLUSPLUS_API void provenance_add_operation(std::string_view rendered);
///@}

/**
 * @brief The lines a writer should emit, without any comment prefix.
 *
 * Always at least `kProvenanceTag`, as line 0. What follows depends on both
 * @p tier and the active `ProvenanceMode`:
 *
 *  - `Off` (the default, no scope open): always just the tag. This is the
 *    whole reason every writer touched here keeps producing v10.15.0's exact
 *    bytes until a caller opts in.
 *  - `BestEffort`: the full block (source, target, operations, notes,
 *    timestamp -- one per line) when `tier == Block`; otherwise just the tag,
 *    silently -- a `SingleLine`/`Bounded` slot structurally cannot hold more
 *    than one line, and a `None` slot cannot hold even that (so a writer with
 *    no slot at all simply never calls this function, the v10.15.0 rule).
 *  - `Required`: identical to `BestEffort` for `Block`/`SingleLine`/`Bounded`
 *    -- degrading to the tag on a slot that structurally cannot hold more is
 *    not a failure, it is the honest maximum that slot can ever offer, and
 *    the caller already knows that from `doc/formats.md`'s tier column.
 *    `SlotTier::None` is the one genuine failure: there is nothing to write
 *    at all, so this throws.
 *
 * @param tier what the calling format's slot can hold.
 * @throws WriteError under `Mode::Required` when @p tier is `SlotTier::None`.
 */
MESHIOPLUSPLUS_API std::vector<std::string> provenance_lines(SlotTier tier);

/**
 * @brief `provenance_lines(tier)`, each line prefixed with @p prefix and
 * newline-terminated, ready to stream to a writer's `std::ostream`.
 *
 * Covers the common case -- every `#`/`!`/`*`/`$`/`//`-style line comment,
 * PLY's `comment ` keyword, and Abaqus's unprefixed `*HEADING` body (pass an
 * empty prefix). Formats whose slot wraps differently (XML's `<!--...-->`,
 * Ansys's split `(0 "...")`/`(1 "...")` records, Nastran's sentinel needing
 * to stay first, Exodus's single netCDF attribute) call `provenance_lines`
 * directly instead.
 */
MESHIOPLUSPLUS_API std::string provenance_render_lines(SlotTier tier, std::string_view prefix);

/**
 * @brief `provenance_lines(tier)` wrapped as one XML comment.
 *
 * A single line renders as `<!--line-->` -- byte-identical to what every
 * `<!--...-->` writer emitted before this existed. Multiple lines render as
 * one comment spanning them (`<!--\nline\nline\n-->`) rather than one
 * `<!--...-->` per line, since a document with a repeated warning-worthy
 * substring is otherwise indistinguishable from several unrelated comments.
 */
MESHIOPLUSPLUS_API std::string provenance_render_xml_comment(SlotTier tier);

/**
 * @brief The timestamp `ProvenanceScope` stamps a record with by default.
 *
 * Honours `SOURCE_DATE_EPOCH` (the reproducible-builds convention) so a
 * byte-comparable build can pin it; returns empty when
 * `MESHIOPLUSPLUS_PROVENANCE_TIMESTAMP=off`. Host and user are never
 * recorded -- these files get shared.
 */
MESHIOPLUSPLUS_API std::string provenance_timestamp();

/**
 * @brief The mode writers use when no `ProvenanceScope` is open.
 *
 * **`BestEffort` by default** -- provenance is on, not opt-in. With nothing
 * scoped a writer still renders any conversion assumptions recorded while it
 * ran, which is the half of the record with no substitute elsewhere. The
 * source/target/operation-chain fields need a caller or driver to set them and
 * so stay absent, and no timestamp is added (a scope sets that), which is what
 * keeps default output deterministic.
 *
 * Overridden by the `MESHIOPLUSPLUS_PROVENANCE` environment variable
 * (`off`/`none`/`0`, or `required`), read once on first use.
 */
MESHIOPLUSPLUS_API ProvenanceMode default_provenance_mode();

/// Sets what `default_provenance_mode()` returns, overriding the environment.
/// Process-wide: this is configuration, not per-write state.
MESHIOPLUSPLUS_API void set_default_provenance_mode(ProvenanceMode mode);

/**
 * @brief Marks the start of a write, bounding scope-less notes to it.
 *
 * With provenance on by default a writer-side `provenance_note` fires during
 * an ordinary `write()` with nothing scoped, and those notes need a lifetime
 * or they attach to whatever file is written *next* -- including an unrelated
 * one. (Reproduced before this existed: `extract_surface(A)` dropping a
 * region put `Note [regions-dropped]: extract_surface: ...` into unrelated
 * mesh B's header.) This gives that lifetime: the write entry points call it,
 * and it resets the ambient record so only notes raised by *this* write can be
 * rendered.
 *
 * **A no-op when a `ProvenanceScope` is open** -- the caller's scope owns the
 * lifetime then, and the whole point of opening one is to span more than a
 * single write.
 *
 * The trade this makes deliberately: a note raised *before* the write begins
 * (the operation-side sites, e.g. `warn_regions_dropped`) is discarded rather
 * than misattributed. Capturing those needs an explicit scope spanning the
 * operations and the write, which is exactly what `ProvenanceScope` is for.
 */
MESHIOPLUSPLUS_API void provenance_begin_write();

/// What `read_provenance_lines` found in a file.
struct ProvenanceReadResult {
    /// The block's lines as found, comment punctuation stripped, in file
    /// order. Empty when the file carries none.
    std::vector<std::string> mLines;
    /// Whether the first line is our own tag format (`"Written by meshio++
    /// v<something>"`). False for a file with no block at all, and for one
    /// whose block does not start the way this library writes one.
    bool mRecognised = false;
};

/**
 * @brief Recovers the provenance block a writer left in @p rPath.
 *
 * **One scanner, not 44 parsers.** `render_block`'s output lines are
 * format-independent by construction -- only the comment punctuation wrapping
 * them differs -- so recovering them is a matter of reading the file's head,
 * stripping whatever leading marker each format uses (`#`, `!`, `*`, `$`,
 * `%`, `//`, PLY's `comment `, an XML `<!--`) plus any trailing box
 * structure, and keeping the runs that match the shapes this file writes.
 *
 * **Deliberately returns raw lines rather than a re-parsed
 * `ProvenanceRecord`.** A block can have been hand-edited, truncated by a
 * `Bounded` slot, or written by a later release carrying fields this build
 * has never heard of; handing back what is actually there, plus a flag for
 * "this at least starts like ours", cannot silently drop or mis-attribute
 * any of those. Callers wanting structure can match on the documented line
 * prefixes themselves.
 *
 * Two formats are **not** served by this scanner and are handled by their own
 * readers instead: exodus keeps the block in a netCDF `title` attribute
 * rather than the head bytes (`read_exodus_metadata` reads it, having already
 * opened the file for `mTimeValues`), and the HDF5 containers carry no block
 * at all (`SlotTier::None`).
 *
 * Never throws for a missing or unreadable file -- an absent block and an
 * unopenable path are both simply "nothing found", since this is a
 * best-effort enrichment of a summary, not a read whose failure should
 * fail the summary.
 *
 * @param rPath the file to scan.
 * @param max_bytes how much of the head to read; the block is always at or
 *        near the top of every slot that can hold one.
 */
MESHIOPLUSPLUS_API ProvenanceReadResult read_provenance_lines(const std::string& rPath,
                                                              std::size_t max_bytes = 8192);

/**
 * @brief `read_provenance_lines`' scanner, over bytes already in hand.
 *
 * The half exodus needs: it has the block as one newline-joined netCDF
 * attribute string rather than a file to open, and re-opening the file to
 * scan it would be both wasteful and wrong (the attribute is not in the head
 * bytes). Same stripping and matching rules, so the two paths cannot
 * disagree about what counts as a block.
 */
MESHIOPLUSPLUS_API ProvenanceReadResult scan_provenance_text(std::string_view text);

}  // namespace detail
}  // namespace meshioplusplus
