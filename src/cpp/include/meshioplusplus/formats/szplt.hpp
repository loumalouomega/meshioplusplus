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
 * @file szplt.hpp
 * @brief Tecplot SZL (`.szplt`) reader through a user-installed TecIO, the
 *        `szplt` format.
 *
 * SZL is Tecplot's undocumented subzone-loadable format; TecIO, Tecplot's own
 * library, is its only reader. meshio++ never redistributes it: the reader is
 * compiled only with `MESHIOPLUSPLUS_WITH_TECIO=ON` and `TECIO_ROOT` pointing at
 * a TecIO the user already has (a Tecplot 360 install, or a build of the TecIO
 * source Tecplot distributes).
 *
 * TecIO's reader fills the same zone model the ASCII and `.plt` readers decode
 * into (see tecplot.hpp), so a `.szplt` reads exactly as the same data saved as
 * `.plt`: one cell block and Cell region per zone of the selected step, point
 * or cell data by each variable's location, sharing and passive variables
 * resolved, `ReadOptions::mTimeStep` over the distinct solution times.
 * FEPOLYGON/FEPOLYHEDRON zones are refused by name.
 * See doc/formats/szplt.md.
 */

#ifdef MESHIOPLUSPLUS_HAS_TECIO

// System includes
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Read one step of a `.szplt` file through TecIO.
 * @throws ReadError if TecIO cannot open the file, a zone is face-based, `X`
 *         or `Y` is missing, or `mTimeStep` is out of range
 */
MESHIOPLUSPLUS_API Mesh read_szplt(const std::string& rPath, const ReadOptions& rOptions = {});

/** @brief The distinct solution times, ascending; empty for a static file. */
MESHIOPLUSPLUS_API std::vector<double> szplt_time_values(const std::string& rPath);

/** @brief Metadata of step 0 (a points-only read) and `mTimeValues`. */
MESHIOPLUSPLUS_API MeshMetadata read_szplt_metadata(const std::string& rPath,
                                                    const ReadOptions& rOptions = {});

/**
 * @brief Write a `.szplt` file through TecIO: the zones of the ASCII and `.plt`
 *        writers (one per cell block, a cell-centred variable passive where a
 *        block has none), doubles throughout, wedges and pyramids as degenerate
 *        bricks. TecIO shares a variable between zones only with their
 *        connectivity, so a mesh of one block is one zone with every point and a
 *        mesh of several gives each zone the points its cells use (a point
 *        several blocks use is repeated, one no cell uses is dropped). A static
 *        data set. The Python twin writes the same zones through a shared TecIO.
 * @throws WriteError when TecIO fails or the mesh has no cell type SZL can hold
 *         (polygons and polyhedra are skipped with a warning)
 */
MESHIOPLUSPLUS_API void write_szplt(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief A strand of zones per step in one `.szplt` file; the later steps share
 *        the first step's coordinates and connectivity, so every step needs the
 *        first's points, cells and variables. The Python twin is
 *        `meshioplusplus.szplt.SeriesWriter`.
 */
class MESHIOPLUSPLUS_API SzpltSeriesWriter {
public:
    explicit SzpltSeriesWriter(const std::string& rPath);
    ~SzpltSeriesWriter();

    SzpltSeriesWriter(const SzpltSeriesWriter&) = delete;
    SzpltSeriesWriter& operator=(const SzpltSeriesWriter&) = delete;
    SzpltSeriesWriter(SzpltSeriesWriter&&) noexcept;
    SzpltSeriesWriter& operator=(SzpltSeriesWriter&&) noexcept;

    /// Add one step. @throws WriteError on a mismatch with the first step or a TecIO failure.
    void Write(double Time, const Mesh& rMesh);

    /// The number of steps added so far.
    std::size_t NumSteps() const noexcept;

    /// Close the file. Idempotent.
    void Finalize();

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
};

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_TECIO
