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
 * @brief Tecplot ASCII finite-element (.dat/.tec) C++ reader/writer,
 *        single-zone only.
 *
 * A Tecplot FE file is a `VARIABLES = "X" "Y" "Z" ...` list followed by one
 * or more `ZONE T="..." N=<nodes> E=<elements> F=FEPOINT|FEBLOCK
 * ET=TRIANGLE|FEQUADRILATERAL|FETETRAHEDRON|FEBRICK [SOLUTIONTIME=<t>]
 * [STRANDID=<id>] [VARLOCATION=(...)]` blocks. meshio++ writes a **single**
 * FE zone; on read, `ReadOptions::mTimeStep` (since v11.3.0) selects one zone
 * of a transient file's timeline -- every zone sharing the first zone's
 * `STRANDID` (or, absent one, every zone with a `SOLUTIONTIME` at all),
 * sorted by `SOLUTIONTIME`. With no `SOLUTIONTIME` anywhere, only the first
 * zone is read, as before (several static, non-transient zones are a
 * documented roadmap remainder, not concatenated or errored on).
 * `VARLOCATION=([a-b]=CELLCENTERED)` (1-based,
 * inclusive ranges) marks cell-centered variables, otherwise cell-centered-
 * ness is inferred from `NV=`. `FEBLOCK` packing reads one variable's full
 * array before the next; `FEPOINT` reads one full-variable-tuple row per
 * node. `X`/`x` (and optional `Y`/`Z`) become point coordinates; everything
 * else becomes point_data or cell_data, keyed by the raw variable name (no
 * `tecplot:` prefix).
 *
 * Zone type -> meshio++ type: LINESEG/FELINESEG->line,
 * TRIANGLE/FETRIANGLE->triangle, QUADRILATERAL/FEQUADRILATERAL->quad,
 * TETRAHEDRON/FETETRAHEDRON->tetra, BRICK/FEBRICK->hexahedron. On write,
 * pyramid/wedge/hexahedron all degrade to an 8-node FEBRICK zone, padding
 * with duplicated corner nodes (pyramid: `[0,1,2,3,4,4,4,4]`; wedge:
 * `[0,1,4,3,2,2,5,5]`).
 *
 * The multi-cell-type write path (Python degrades everything into a single
 * FEQUADRILATERAL/FEBRICK zone via "order_2" padding tables) exists **only
 * in the Python writer**: the C++ writer throws WriteError as soon as more
 * than one distinct cell type is present, forcing the Python fallback.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as a single Tecplot FE zone.
 *
 * Emits `VARIABLES`/`ZONE` headers for the one supported cell type present
 * (line/triangle/quad/tetra/hexahedron, with pyramid/wedge padded into
 * FEBRICK), then FEBLOCK-packed coordinate/point_data/cell_data columns
 * (data wrapped at 20 values per line) and 1-based connectivity.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError if the mesh contains **more than one** distinct cell
 *         type (the Python fallback handles that case by degrading
 *         everything into one FEQUADRILATERAL/FEBRICK zone)
 */
MESHIOPLUSPLUS_API void write_tecplot(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a Tecplot ASCII file's first FE zone.
 *
 * Parses the `VARIABLES` list and every `ZONE` header (tolerating multi-line
 * continuation and a quoted `T="..."` title), then the selected zone's
 * FEBLOCK or FEPOINT data body and 1-based connectivity.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh
 * @throws ReadError if `X`/`x` is missing, a zone header uses an unsupported
 *         `F=`/`ZONETYPE=` combination, or the header/data otherwise doesn't
 *         parse (e.g. an adversarial zone title that is literally the string
 *         `"VARLOCATION"`) — the shim then falls back to the more tolerant
 *         Python reader.
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
 * every `ZONE` header is scanned for its `N=`/`E=`/`SOLUTIONTIME=`/
 * `STRANDID=`, using the same token-budget walk `read_tecplot` uses to skip
 * from one zone's header to the next, but without decoding the tokens
 * themselves. `mCellBlocks`/`mNumPoints` describe the timeline's first zone
 * (transient zones typically share topology); `mTimeValues` is the resolved
 * timeline's `SOLUTIONTIME`s in the same sorted order `mTimeStep` indexes
 * into -- empty when no zone carries one.
 *
 * @param rPath filesystem path to read
 * @param rOptions unused (metadata carries no timestep of its own to select)
 * @return the file's shape and time values
 * @throws ReadError as `read_tecplot`.
 */
MESHIOPLUSPLUS_API MeshMetadata read_tecplot_metadata(const std::string& rPath,
                                                      const ReadOptions& rOptions);

}  // namespace meshioplusplus
