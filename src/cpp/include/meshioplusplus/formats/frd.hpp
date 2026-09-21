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
 * @file frd.hpp
 * @brief CalculiX result file (`.frd`) C++ reader.
 *
 * The ASCII file `ccx` writes and `cgx` reads: fixed-column records keyed by their
 * first columns (`1C`/`1U` header, `2C` nodes, `3C` elements, one `100C` block per
 * result per increment, `9999`). Values are `E12.5` with no separator, so every field
 * is sliced by column. The short (`I5` ids, flag 0) and long (`I10` ids, flag 1)
 * layouts are read; the binary layout (flag 2) is refused.
 *
 *  - The mesh is the one `ccx` wrote, not the `.inp` mesh: shells and beams are
 *    expanded into solids. The twelve cgx element types map to `hexahedron`, `wedge`,
 *    `tetra`, `hexahedron20`, `wedge15`, `tetra10`, `triangle`, `triangle6`, `quad`,
 *    `quad8`, `line` and `line3`; the he20/pe15 mid-edge groups and the be3 mid node
 *    are permuted to the Abaqus/meshio++ order. Element group and material become
 *    the integer cell data `frd:group` and `frd:material`.
 *  - Each `100C` increment is one step for the sequence engine: `ReadOptions::mTimeStep`
 *    selects it (negative counts from the end), and its value, step number and
 *    analysis type (0 static, 1 time, 2 frequency, 3 buckling, 4 user) are the field
 *    data `meshio:time`, `frd:step` and `frd:analysis`.
 *  - Every result block becomes a point data array named as the file names it, NaN at
 *    the nodes the block has no value for; components the file marks as calculated
 *    (the `ALL` of `DISP`) are not data. Symmetric tensors keep the file's order
 *    `xx yy zz xy yz zx`.
 *  - `FrdReadOptions::mDerived` adds `<NAME>_mises` and `<NAME>_principal` (ascending
 *    min, mid, max) beside each `STRESS`/`TOSTRAIN`/`MESTRAIN` tensor.
 *
 * See doc/formats/frd.md for the record layouts and the limits.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/** @brief Format-specific options of `read_frd`. */
struct FrdReadOptions {
    /** @brief Add `<NAME>_mises` and `<NAME>_principal` beside each 6-component tensor. */
    bool mDerived = false;
};

/**
 * @brief Read a CalculiX `.frd` result file.
 *
 * @param rPath filesystem path to read
 * @param rOpts `mTimeStep` selects the increment (`ResolveTimeStep`); `mPointsOnly` and
 *        `mDataArrays` narrow the result blocks that are parsed
 * @return the mesh, with the selected increment's results as point data
 * @throws ReadError if the file can't be read, is binary, a field is malformed, an
 *         element or a result refers to an undefined node, or the step is out of range
 */
MESHIOPLUSPLUS_API Mesh read_frd(const std::string& rPath, const ReadOptions& rOpts = {});

/** @brief `read_frd` with the format-specific options. */
MESHIOPLUSPLUS_API Mesh read_frd(const std::string& rPath, const ReadOptions& rOpts,
                                 const FrdReadOptions& rFrdOpts);

/**
 * @brief Summarize a `.frd` file: the mesh, step 0's result names, and every
 *        increment's value as `mTimeValues`.
 */
MESHIOPLUSPLUS_API MeshMetadata read_frd_metadata(const std::string& rPath,
                                                  const ReadOptions& rOpts = {});

}  // namespace meshioplusplus
