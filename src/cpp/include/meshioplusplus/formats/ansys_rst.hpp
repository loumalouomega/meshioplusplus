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
 * @file ansys_rst.hpp
 * @brief Ansys MAPDL binary results (`.rst`, `.rth`) C++ reader.
 *
 * A results file is a sequence of Fortran-style records -- `i32 length` (in
 * 4-byte words), `i32 flags`, the data, a trailer word -- addressed by word
 * offsets from the start of the file. Integer, int16, float32, float64,
 * bit-sparse and windowed-sparse records are read; zlib-compressed ones
 * (`/FCOMP`) are refused. The layout follows Ansys's `fdresu.inc` as the open
 * reader pymapdl-reader (MIT) documents it.
 *
 *  - The geometry records give the nodes, element types, elements and
 *    components, which become cells, `ansys:*` cell data and point/cell regions
 *    exactly as a `.cdb` deck's do (`detail::ansys_build_mesh`).
 *  - `ReadOptions::mTimeStep` picks the result set. Its time (a frequency in a
 *    modal analysis) is `field_data["meshio:time"]`; `"ansys:load_step"`,
 *    `"ansys:substep"` and `"ansys:cumulative"` place it.
 *  - The set's nodal DOF solution is point data: `U`, `ROT`, `A` and `V` as
 *    3-component vectors, rotated from each node's coordinate system to the
 *    global one; other DOFs (`TEMP`, `PRES`, `VOLT` ...) as scalars. A node the
 *    set has no solution for, and MAPDL's undefined value, are NaN.
 *  - A partial file of a distributed solve is refused (read the combined file);
 *    of a cyclic-symmetry model only the base sector is read.
 * See doc/formats/ansys_rst.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Read one result set of an Ansys MAPDL results file.
 * @param rPath filesystem path to read
 * @param rOptions `mTimeStep` picks the set (negative counts from the end);
 *        `mPointsOnly`/`mDataArrays` narrow the data read; `mLenient` skips
 *        elements whose type has no meshio++ cell
 * @return the mesh with the chosen set's nodal solution
 * @throws ReadError for a file that is not a MAPDL results file, a compressed
 *         or truncated record, a distributed partial file, or an out-of-range
 *         time step
 */
MESHIOPLUSPLUS_API Mesh read_ansys_rst(const std::string& rPath, const ReadOptions& rOptions = {});

/**
 * @brief Summarise an Ansys MAPDL results file, its set times included.
 */
MESHIOPLUSPLUS_API MeshMetadata read_ansys_rst_metadata(const std::string& rPath,
                                                        const ReadOptions& rOptions = {});

}  // namespace meshioplusplus
