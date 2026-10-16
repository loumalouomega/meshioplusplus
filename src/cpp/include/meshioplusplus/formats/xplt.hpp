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
 * @file xplt.hpp
 * @brief FEBio plot file (`.xplt`) C++ reader: the results FEBio writes.
 *
 * An `.xplt` file is a tree of tagged chunks (`u32 id, u32 size, payload`)
 * after a 4-byte magic: a root holding the header and the dictionary of plot
 * variables, the mesh, then one state per saved time. With compression on,
 * every chunk after the first mesh is its own zlib stream. Plot versions
 * 0x0030 and later (FEBio 3 and 4) are read, either byte order.
 *
 *  - Each domain is a cell block and a `Cell` region (named after the domain,
 *    else the element set holding exactly its elements, else its part), tagged
 *    with its part id. hex27 goes through the `"febio"` node-order table.
 *  - Node sets, element sets and surfaces become `Point`, `Cell` and `Side`
 *    regions, as in `.feb`. Each data surface is also a block of facet cells
 *    per facet type, after the domains, with a `Cell` region named after it
 *    and `cell_data["xplt:surface"]` (its 1-based id, 0 on the domains).
 *  - A remeshed run writes a new mesh before the states that use it; each
 *    state is read on its own mesh. Metadata describes the first.
 *  - `ReadOptions::mTimeStep` picks the state. Its time is
 *    `field_data["meshio:time"]`, its index `"xplt:step"`, its status
 *    `"xplt:status"`. Nodal variables are point data; per-element variables
 *    (`FMT_ITEM`) cell data, NaN on domains without them; per-region ones
 *    (`FMT_REGION`) cell data repeated over the domain; per-element-node ones
 *    (`FMT_NODE`, `FMT_MULT`) point data, averaged over the elements sharing a
 *    node. Global variables are field data. Surface variables live on the facet
 *    cells the same way (per facet, per surface, or averaged per facet node);
 *    one on no surface of the mesh gives no array. Edge and material-point
 *    variables are not read (a warning names them).
 *  - Symmetric tensors keep FEBio's 6 components (xx, yy, zz, xy, yz, xz).
 *
 * See doc/formats/xplt.md.
 */

// System includes
#include <cstddef>
#include <memory>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Read one state of an FEBio plot file.
 * @param rPath filesystem path to read
 * @param rOptions `mTimeStep` picks the state (negative counts from the end);
 *        `mPointsOnly`/`mDataArrays` narrow the data read; `mLenient`
 *        downgrades tet5/tet15 domains
 * @return the mesh with the chosen state's data
 * @throws ReadError for a bad magic, an unsupported plot version, a
 *         truncated mesh, or an out-of-range time step
 */
MESHIOPLUSPLUS_API Mesh read_xplt(const std::string& rPath, const ReadOptions& rOptions = {});

/**
 * @brief Summarise an FEBio plot file, its state times included.
 */
MESHIOPLUSPLUS_API MeshMetadata read_xplt_metadata(const std::string& rPath,
                                                   const ReadOptions& rOptions = {});

/**
 * @brief Write an FEBio plot file holding one state (plot version 0x35,
 *        little-endian, uncompressed): the mesh and its data at
 *        `field_data["meshio:time"]` (0 when absent), status `"xplt:status"`.
 *        Cell blocks become domains (a region that is exactly one block names
 *        it and tags its part), point regions node sets and cell regions
 *        element sets; point data are nodal variables, cell data per-element
 *        variables (`FMT_ITEM`) and numeric field data global ones. Vertices,
 *        polygons, side regions and surfaces are not written. The Python twin
 *        writes the same bytes.
 * @throws WriteError when the file cannot be written or a cell names a missing point.
 */
MESHIOPLUSPLUS_API void write_xplt(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief A plot file with a state per step: the mesh and the variables of the
 *        first step, which every later step must share. The file is written by
 *        `Finalize`. The Python twin is `meshioplusplus.xplt.SeriesWriter`.
 */
class MESHIOPLUSPLUS_API XpltSeriesWriter {
public:
    /// @throws WriteError when @p rPath cannot be opened.
    explicit XpltSeriesWriter(const std::string& rPath);
    ~XpltSeriesWriter();

    XpltSeriesWriter(const XpltSeriesWriter&) = delete;
    XpltSeriesWriter& operator=(const XpltSeriesWriter&) = delete;
    XpltSeriesWriter(XpltSeriesWriter&&) noexcept;
    XpltSeriesWriter& operator=(XpltSeriesWriter&&) noexcept;

    /// Add one step. @throws WriteError when its points, cells or variables differ
    /// from the first step's, or on a moved-from writer.
    void Write(double Time, const Mesh& rMesh);

    /// The number of steps added so far.
    std::size_t NumSteps() const noexcept;

    /// Write the file. Idempotent. @throws WriteError when no step was added.
    void Finalize();

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
};

}  // namespace meshioplusplus
