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
 *  - Reaction forces are point data `RF` and `RMOM` (rotated to the global
 *    axes) and `RF_<DOF label>`; NaN where a node has none.
 *  - Element nodal stresses and strains (`S`, `EPEL`, `EPPL`, `EPCR`, `EPTH`:
 *    `xx yy zz xy yz xz`, rotated from the element system by its Euler angles)
 *    are cell data per element node, `(cells, nodes * 6)` flattened
 *    point-major with the widest block's node count in every block (NaN at nodes
 *    that carry none, midside nodes, and past a narrower cell's nodes;
 *    `field_data["ansys:layout:<name>"]` is `[nodes, 6]`), and point data
 *    averaged over the elements at each corner node; a layered shell's top
 *    surface is `<name>@top`. Element nodal forces are cell data `ENF`,
 *    `(cells, nodes * DOFs)` in the set's DOF order. Line and point elements
 *    carry none.
 *  - The main file of a distributed solve (`<job>0.rst`) reads its partial
 *    files (`<job>1.rst` ...) with it, merged by node number; another partial
 *    file is refused. Of a cyclic-symmetry model only the base sector is read;
 *    `read_ansys_rst_cyclic` expands a static one to the full rotor.
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
 *         or truncated record, a distributed partial file other than the main
 *         one (or a main one whose partial files are missing), or an
 *         out-of-range time step
 */
MESHIOPLUSPLUS_API Mesh read_ansys_rst(const std::string& rPath, const ReadOptions& rOptions = {});

/**
 * @brief Summarise an Ansys MAPDL results file, its set times included.
 */
MESHIOPLUSPLUS_API MeshMetadata read_ansys_rst_metadata(const std::string& rPath,
                                                        const ReadOptions& rOptions = {});

/**
 * @brief Read one result set of a static cyclic-symmetry model as the full rotor.
 *
 * The base sector's cells (element numbers up to the model's `csEls`) and their
 * points are repeated round the cyclic axis (global Z, or the Z axis of the
 * local coordinate system the model names), each copy's vectors and tensors
 * rotated with it; `"ansys:sector"` cell data numbers the copies and
 * `"ansys:sectors"` field data counts them. Coincident nodes on the sector
 * boundaries are not merged. Otherwise as read_ansys_rst().
 * @throws ReadError for a model that is not cyclic, or not a static analysis
 */
MESHIOPLUSPLUS_API Mesh read_ansys_rst_cyclic(const std::string& rPath,
                                              const ReadOptions& rOptions = {});

/**
 * @brief Summarise a static cyclic-symmetry model's full rotor, its set times included.
 */
MESHIOPLUSPLUS_API MeshMetadata read_ansys_rst_cyclic_metadata(const std::string& rPath,
                                                               const ReadOptions& rOptions = {});

}  // namespace meshioplusplus
