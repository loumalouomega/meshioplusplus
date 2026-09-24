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
 * @file tecplot.hpp
 * @brief Tecplot reader (ASCII `.dat`/`.tec` and binary `.plt`) and ASCII
 *        writer.
 *
 * Reading: both encodings decode into one zone model. The binary form is the
 * Data Format Guide's appendix A (`#!TDV112`, either byte order); a file is
 * binary when it starts with `#!TDV`, whatever its extension. Every zone of the
 * selected step becomes its own cell block and a Cell region named after the
 * zone title (`zone_<k>` otherwise), with `cell_data["tecplot:zone"]` holding
 * the zone index. FE zones map LINESEG/TRIANGLE/QUADRILATERAL/TETRAHEDRON/BRICK
 * to line/triangle/quad/tetra/hexahedron; ORDERED zones become lines, quads or
 * hexahedra over their dimensions longer than one (VTK corner order, `i`
 * fastest). FEPOLYGON/FEPOLYHEDRON zones are refused with a `ReadError`.
 * VARSHARELIST, PASSIVEVARLIST (NaN) and CONNECTIVITYSHAREZONE are resolved;
 * a zone reuses an earlier zone's points only when all its nodal variables
 * are shared from it. `X`/`Y` (and `Z`) are the coordinates; every other
 * variable is point data where a zone stores it at the nodes and cell data
 * where it stores it cell-centred. `ReadOptions::mTimeStep` selects one
 * distinct SOLUTIONTIME of a transient file (the zones sharing the first
 * zone's strand); a file without SOLUTIONTIME is one step holding every zone.
 *
 * Writing: ASCII only, one FEBLOCK zone per cell block, later zones sharing
 * the first zone's coordinates and nodal fields through VARSHARELIST;
 * pyramid/wedge/hexahedron degrade to an 8-node FEBRICK (pyramid:
 * `[0,1,2,3,4,4,4,4]`; wedge: `[0,1,4,3,2,2,5,5]`).
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as Tecplot ASCII, one FE zone per supported cell block.
 *
 * Emits `VARIABLES`/`ZONE` headers, FEBLOCK-packed coordinate/point_data/
 * cell_data columns (wrapped at 20 values per line) and 1-based
 * connectivity; later zones share the first zone's coordinates and nodal
 * fields through VARSHARELIST.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError if no cell block has a Tecplot FE type
 */
MESHIOPLUSPLUS_API void write_tecplot(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a Tecplot ASCII or binary (`.plt`) file's first step.
 *
 * See the file comment for how zones become cell blocks, regions and data.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh
 * @throws ReadError if `X` or `Y` is missing, a zone is polygonal/polyhedral,
 *         a `.plt` is not version 112 or is truncated, or the file otherwise
 *         does not parse
 * @note point_data/cell_data keys are the raw Tecplot variable names (no
 *       prefix); `X`/`Y`/`Z` are reserved for coordinates.
 */
MESHIOPLUSPLUS_API Mesh read_tecplot(const std::string& rPath);

/**
 * @brief `read_tecplot` with read options — selects one zone of a transient
 *        file's timeline instead of always the first.
 *
 * @param rPath filesystem path to read
 * @param rOptions read options; only `mTimeStep` is consulted
 * @return the read Mesh, as the plain overload
 * @throws ReadError as the plain overload, plus a `mTimeStep` out of range
 *         naming the step count
 */
MESHIOPLUSPLUS_API Mesh read_tecplot(const std::string& rPath, const ReadOptions& rOptions);

/**
 * @brief Summarize a Tecplot file's shape and available time steps without
 *        decoding any zone's data body.
 *
 * A native metadata path (`MeshMetadata::mFellBackToFullRead` is `false`):
 * every zone header is scanned (the ASCII token-budget walk, or the binary
 * header and data-section offsets) without decoding any values. `mCellBlocks`/`mNumPoints` describe
 * the timeline's first zone (transient zones typically share topology); `mTimeValues` is the
 * resolved timeline's `SOLUTIONTIME`s in the same sorted order `mTimeStep` indexes into -- empty
 * when no zone carries one.
 *
 * @param rPath filesystem path to read
 * @param rOptions unused (metadata carries no timestep of its own to select)
 * @return the file's shape and time values
 * @throws ReadError as `read_tecplot`.
 */
MESHIOPLUSPLUS_API MeshMetadata read_tecplot_metadata(const std::string& rPath,
                                                      const ReadOptions& rOptions);

}  // namespace meshioplusplus
