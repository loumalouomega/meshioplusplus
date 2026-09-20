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
 * @file formats/pvd.hpp
 * @brief ParaView collection (`.pvd`): a time-indexed list of VTK XML files
 * (v15.0.0).
 *
 * A `.pvd` is `<VTKFile type="Collection"><Collection><DataSet timestep="t"
 * part="p" group="g" file="..."/>...`. Entries name serial or parallel XML files
 * (`.vtu`, `.vtp`, `.vtm`, `.pvtu`, `.pvtp`), never legacy `.vtk`, so
 * `.pvd` -> `.pvtu` -> `.vtu` is the ordinary layout of a partitioned transient
 * run and nests without special cases.
 *
 * ### Two axes
 *
 * `timestep` selects the *step* (`ReadOptions::mTimeStep`) and, within a step,
 * `part` selects the *piece* (`ReadOptions::mPiece`). The steps are the distinct
 * `timestep` values in ascending order; a missing `timestep` or `part` is 0, as
 * ParaView reads it. `read_pvd` returns step 0 with every part merged (one
 * `RegionKind::Cell` region per entry, named from `name=`, else `group/part_<p>`,
 * else `part_<p>`); a step with a single entry is returned as that file reads, so
 * a step that is one `.pvtu` keeps its own `piece_<i>` regions. The chosen step's
 * time is attached as `field_data["meshio:time"]`.
 *
 * `read_pvd_metadata` reports every step's time from the index alone, without
 * opening a piece; the rest of the summary is step 0's pieces'.
 *
 * ### Write
 *
 * `write_pvd` writes a one-step collection: one `.vtu` in a sibling directory
 * named after the index's stem (`<stem>/<stem>_0000.vtu`, relative paths) whose
 * `timestep` is the mesh's `meshio:time` field data when present, else 0.
 * `PvdSeriesWriter` streams many steps: one mesh alive at a time, and the index
 * is rewritten after every step, so a run that is killed leaves a collection
 * ParaView opens covering every finished step.
 */

// System includes
#include <memory>
#include <string>

// Project includes
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/formats/pvtu.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Streams `(time, mesh)` steps into a `.pvd` plus one `.vtu` per step.
 *
 * Move-only. A moved-from writer's `Write` throws `WriteError`; `Finalize` is
 * idempotent.
 */
class MESHIOPLUSPLUS_API PvdSeriesWriter {
public:
    /**
     * @param rPath the `.pvd` index; pieces land in `<stem>/` beside it.
     * @param binary base64-encode each piece's arrays instead of writing text.
     * @param codec the block compressor; `None` writes uncompressed base64.
     * @throws WriteError when @p codec is not in this build, or the piece
     *         directory cannot be created.
     */
    explicit PvdSeriesWriter(const std::string& rPath, bool binary = true,
                             detail::VtkCodec codec = detail::VtkCodec::Zlib);
    ~PvdSeriesWriter();

    PvdSeriesWriter(const PvdSeriesWriter&) = delete;
    PvdSeriesWriter& operator=(const PvdSeriesWriter&) = delete;
    PvdSeriesWriter(PvdSeriesWriter&&) noexcept;
    PvdSeriesWriter& operator=(PvdSeriesWriter&&) noexcept;

    /**
     * @brief Write one step: its piece file, then the index listing every step so far.
     * @param Time the step's time; must be finite.
     * @throws WriteError on a non-finite time, an unwritable path, or a moved-from writer.
     */
    void Write(double Time, const Mesh& rMesh);

    /// The number of steps written so far.
    std::size_t NumSteps() const noexcept;

    /// Write the final index. Idempotent. @throws WriteError when no step was written.
    void Finalize();

private:
    struct Impl;
    std::unique_ptr<Impl> mImpl;
};

/**
 * @brief Write a one-step `.pvd` collection.
 * @param rPath the `.pvd` path.
 * @param rMesh the mesh; its `field_data["meshio:time"]` (one value) is the step's time.
 * @param binary base64-encode the piece's arrays instead of writing text.
 * @param zlib compress the piece's binary blocks. Ignored when @p binary is false.
 */
MESHIOPLUSPLUS_API void write_pvd(const std::string& rPath, const Mesh& rMesh, bool binary = true,
                                  bool zlib = true);

/// `write_pvd` with an explicit block codec.
MESHIOPLUSPLUS_API void write_pvd_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                                        detail::VtkCodec codec);

/**
 * @brief Read one step of a `.pvd`, its parts merged (or one part).
 * @param rOpts `mTimeStep` selects the step (`ResolveTimeStep`: negative counts
 *        from the end, out of range names the step count); `mPieceSet` selects one
 *        `part` of that step; the narrowing options and `mGhosts` reach every piece.
 * @throws ReadError on an unparsable index, a missing piece file, or an entry
 *         that is not one of `.vtu`/`.vtp`/`.vtm`/`.pvtu`/`.pvtp`.
 */
MESHIOPLUSPLUS_API Mesh read_pvd(const std::string& rPath, const ReadOptions& rOpts = {});

/// Every step's time (`mTimeValues`) from the index alone, plus step 0's pieces' summary.
MESHIOPLUSPLUS_API MeshMetadata read_pvd_metadata(const std::string& rPath,
                                                  const ReadOptions& rOpts = {});

}  // namespace meshioplusplus
