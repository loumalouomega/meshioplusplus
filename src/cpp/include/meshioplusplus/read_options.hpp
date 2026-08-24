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
 * @file read_options.hpp
 * @brief Per-call reader options (selective/partial reads, memory mapping) and
 *        the lightweight mesh summary returned by `read_metadata`.
 *
 * Readers are all-or-nothing by default: every `<DataArray>`, every section,
 * fully decoded. `ReadOptions` lets a caller ask for less -- geometry only, a
 * named subset of data arrays, or just the header summary -- without
 * materializing the rest. **Every field defaults to "read everything"**, so a
 * default-constructed `ReadOptions` reproduces the historical behaviour exactly;
 * that is what lets the option be threaded through a reader without touching any
 * existing caller or test.
 *
 * Only some readers can act on the options (currently VTU, VTP, XDMF and Gmsh).
 * The rest fall back to a full read, which is correct but not faster -- and
 * always *says so*, via `MeshMetadata::mFellBackToFullRead`. A partial read that
 * silently wasn't partial would be worse than no feature at all.
 *
 * What `mMetadataOnly` actually buys is format-dependent and worth stating
 * plainly: for Gmsh it is close to free (the `$Nodes`/`$Elements` headers carry
 * the counts, so the bodies are skipped outright), while for VTU/VTP the whole
 * file is still read and parsed as XML -- what is skipped is base64 decoding,
 * decompression, allocation and byte-swapping. That is a solid multiple, not an
 * asymptotic change. Genuinely O(1) VTU metadata would need the *appended* data
 * format with `offset=` attributes, which this reader does not accept.
 */

// System includes
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

/**
 * @brief Whether a reader should memory-map its input.
 *
 * `Auto` maps only regular files at or above `mmap_auto_threshold_bytes`, where
 * avoiding the whole-file copy is worth the mapping setup; smaller files are
 * read normally. Mapping is advisory throughout -- an unmappable input (a pipe,
 * a zero-length file, a platform without support, WASM) silently falls back to
 * buffered reading rather than failing, so `On` is a preference, not a demand.
 */
enum class MmapMode { Auto, On, Off };

/** @brief File size at or above which `MmapMode::Auto` maps instead of copying. */
inline constexpr std::size_t mmap_auto_threshold_bytes = 16u * 1024u * 1024u;

/**
 * @brief Per-call options narrowing what a reader materializes.
 *
 * Defaults reproduce a full read exactly. The three narrowing options are
 * ordered by how much they skip: `mDataArrays` (a named subset) is the least
 * aggressive, `mPointsOnly` drops all data, `mMetadataOnly` drops the heavy
 * arrays entirely.
 */
struct ReadOptions {
    /** @brief Read geometry (points + connectivity) but no data arrays at all. */
    bool mPointsOnly = false;

    /**
     * @brief Read only the header/summary, skipping the heavy arrays.
     *
     * Used by `read_metadata`; a reader honouring this may return a `Mesh` whose
     * arrays are absent or empty, so it is not meaningful on the `read` path.
     */
    bool mMetadataOnly = false;

    /**
     * @brief Restrict point/cell/field data to these names.
     *
     * `std::nullopt` means every array (the default); an **empty vector** means
     * none. The distinction is load-bearing -- a bare `std::vector` could not
     * express both, and the ambiguity would be unfixable once it reached the C
     * ABI. Names not present in the file are ignored, not an error: a caller
     * asking for `{"u", "v"}` across a directory of meshes should not have to
     * know which files happen to carry which.
     */
    std::optional<std::vector<std::string>> mDataArrays;

    /** @brief Memory-mapping preference; honoured where the reader supports it. */
    MmapMode mMmap = MmapMode::Auto;

    /**
     * @brief Which time step to materialize from a multi-step file.
     *
     * Formats carrying a time series (Exodus `time_whole`, XDMF temporal
     * collections, CGNS, MED) hold one value per array *per step*, but a `Mesh`
     * holds exactly one. Before this option every such reader silently took the
     * first step, which is correct for the overwhelmingly common single-step
     * file and quietly wrong for the rest.
     *
     * `0` (the default) is that historical behaviour -- the first step -- so a
     * default-constructed `ReadOptions` still reproduces it exactly. Negative
     * values count from the end, `-1` being the last step, which is what a
     * caller wanting "the final state of the simulation" actually means and
     * cannot express without knowing the step count up front. Out of range is an
     * error naming the available count, never a silent clamp: quietly handing
     * back step 0 when step 7 was requested is the failure mode this option
     * exists to remove.
     *
     * A reader for a format with no time concept ignores this field.
     */
    int mTimeStep = 0;

    /**
     * @brief Downgrade "this reader cannot represent construct X" to a warning.
     *
     * A reader that meets a construct it has no way to express normally throws
     * `ReadError` naming it, so a caller with a Python fallback can take it and a
     * caller without one at least learns what was in the file. That is the right
     * default, but it makes formats whose *optional* sections are common
     * -- MDPA's `Table`/`Geometries`/`Constraints` blocks, say -- unreadable in
     * practice for a caller that only wants the mesh.
     *
     * `true` turns those specific rejections into a `log::warn` plus a skip.
     * `false` (the default) reproduces the historical behaviour exactly.
     *
     * This is deliberately **not** "ignore all errors". It only applies where a
     * reader can skip a construct and still return a *correct* mesh: a malformed
     * file, a truncated block, a bad node reference or an unknown element type
     * still throw, because continuing past those would hand back a mesh that is
     * quietly wrong rather than merely incomplete.
     *
     * Currently honoured by `mdpa` only; every other reader ignores it.
     */
    bool mLenient = false;

    /**
     * @brief Resolve `mTimeStep` against an actual step count.
     *
     * @param NumSteps Number of steps the file carries.
     * @return The 0-based step index.
     * @throws ReadError when the request is out of range.
     */
    MESHIOPLUSPLUS_API std::size_t ResolveTimeStep(std::size_t NumSteps) const;

    /** @brief Whether @p rName survives the `mDataArrays` filter. */
    bool WantsArray(const std::string& rName) const {
        if (!mDataArrays.has_value())
            return true;
        for (const std::string& candidate : *mDataArrays)
            if (candidate == rName)
                return true;
        return false;
    }

    /** @brief Whether any data array at all should be read. */
    bool WantsAnyData() const {
        if (mPointsOnly || mMetadataOnly)
            return false;
        return !mDataArrays.has_value() || !mDataArrays->empty();
    }
};

/** @brief One cell block's shape, without its connectivity. */
struct CellBlockInfo {
    std::string mType;              ///< meshio++ cell type name, e.g. `"tetra10"`.
    std::size_t mNumCells = 0;      ///< Number of cells in the block.
    std::size_t mNodesPerCell = 0;  ///< Nodes per cell; 0 when ragged.
    bool mRagged = false;           ///< Whether the block is a polygon/polyhedron block.
};

/**
 * @brief One named region's shape, without its entries.
 *
 * The region complement to `CellBlockInfo`: enough to enumerate and identify a
 * mesh's regions (build a SubModelPart tree, say) without the cost of a full
 * read just to learn how many there are and what they're called.
 */
struct RegionSummary {
    std::string mName;                     ///< The region's name.
    RegionKind mKind = RegionKind::Point;  ///< Point / Cell / Side.
    int mDim = -1;                         ///< Topological dimension, or -1 if unspecified.
    std::int64_t mTag = -1;                ///< Format-native integer id, or -1 if none.
    std::size_t mNumEntries = 0;  ///< Number of grouped entities (not the entries themselves).
};

/**
 * @brief A mesh summary: what a file contains, without its contents.
 *
 * The topological complement to `compute_stats` (geometry) and `data_info`
 * (data values) -- this one is about *shape and names*, and is cheap enough to
 * run on a file you have no intention of loading.
 */
struct MeshMetadata {
    std::size_t mNumPoints = 0;  ///< Number of points.
    std::size_t mPointDim = 0;   ///< Coordinate dimension (2 or 3).

    std::vector<CellBlockInfo> mCellBlocks;  ///< Per-block shape, in file order.

    /** @name Available array names, sorted (matching the uniform API's guarantee). */
    ///@{
    std::vector<std::string> mPointDataNames;
    std::vector<std::string> mCellDataNames;
    std::vector<std::string> mFieldDataNames;
    ///@}

    /**
     * @brief Whether the bounding box below was computed.
     *
     * False on a native metadata path: the bounding box requires decoding the
     * point coordinates, usually the single largest array in the file, which
     * would defeat the purpose. It is only filled in when a full read happened
     * anyway (see `mFellBackToFullRead`).
     */
    bool mHasBBox = false;
    double mBBoxMin[3] = {0.0, 0.0, 0.0};
    double mBBoxMax[3] = {0.0, 0.0, 0.0};

    /**
     * @brief Whether the whole file had to be read to produce this summary.
     *
     * True for every format without a native metadata path. The summary is still
     * correct -- it just wasn't cheap. Exposed everywhere (Python, C, Fortran,
     * JS, both CLIs) so a caller can tell "fast" from "worked".
     */
    bool mFellBackToFullRead = false;

    /** @brief The resolved format name the summary was read as. */
    std::string mFormat;

    /**
     * @brief The file's time-series values, or empty when it carries none.
     *
     * The companion to `ReadOptions::mTimeStep`: `mTimeValues.size()` is the
     * number of steps a caller may ask for, and the values themselves are the
     * simulation times. Empty means either "this format has no time concept" or
     * "this file is single-step with no recorded time" -- the two are not worth
     * distinguishing, since neither leaves a caller anything to choose between.
     */
    std::vector<double> mTimeValues;

    /**
     * @brief The file's named regions, without their entries.
     *
     * Populated whenever the mesh producing this summary was already fully in
     * memory (i.e. whenever `metadata_from_mesh` ran) -- regions are read
     * alongside geometry by every region-capable reader, so there is nothing
     * left to save by skipping them once a read has happened anyway. A native
     * metadata path that declines a full read (VTU/VTP/XDMF/Gmsh 4.1) reports
     * no regions rather than guessing; none of those formats map regions yet
     * regardless (see `doc/regions.md`), so this is never a wrong answer, only
     * an incomplete one on formats already known to fall back for other
     * reasons.
     */
    std::vector<RegionSummary> mRegions;

    /**
     * @brief The provenance block the file carries, or empty when it has none.
     *
     * The lines as found, comment punctuation stripped, in file order --
     * deliberately not re-parsed into a `ProvenanceRecord`, since a block can
     * have been hand-edited, truncated by a fixed-width slot, or written by a
     * later release carrying fields this build has never heard of. See
     * `detail/provenance.hpp`.
     *
     * Filled from the file's bytes by `registry_read_metadata`, on both the
     * native and the fall-back path -- unlike everything else here it cannot
     * come from `metadata_from_mesh`, since a `Mesh` does not carry it.
     */
    std::vector<std::string> mProvenance;

    /**
     * @brief Whether `mProvenance`'s first line is meshio++'s own tag format.
     *
     * False both for a file with no block at all and for one whose block does
     * not start the way this library writes one -- the honest distinction
     * between "meshio++ wrote this" and "something left a comment here".
     */
    bool mProvenanceRecognised = false;

    /** @brief Total cells across every block. */
    std::size_t NumCells() const {
        std::size_t total = 0;
        for (const CellBlockInfo& block : mCellBlocks)
            total += block.mNumCells;
        return total;
    }
};

/**
 * @brief Summarize an already-loaded mesh.
 *
 * The single fallback implementation behind `read_metadata` for every format
 * lacking a native metadata path -- one function rather than one per format.
 * Takes the mesh by const reference and returns a plain aggregate: it must never
 * copy the mesh, because the KRATOS backend's `Mesh` is not copy-constructible.
 *
 * Does **not** set `mFellBackToFullRead`; that is the caller's business, since
 * this is also useful on a mesh that was going to be read regardless.
 */
MESHIOPLUSPLUS_API MeshMetadata metadata_from_mesh(const Mesh& rMesh);

}  // namespace meshioplusplus
