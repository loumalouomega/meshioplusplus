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
 * @file formats/plot3d.hpp
 * @brief PLOT3D multi-block structured grids, with `.q` solutions and `.f` function
 * files (roadmap 2.4.2).
 *
 * A PLOT3D grid does not state its variant, so the reader infers it: ASCII,
 * Fortran sequential unformatted (record markers, read through
 * `detail/fortran_records.hpp`) or raw binary (C stream); single or multi-block
 * (a leading block count); 2-D (`ni nj`, `x y`) or 3-D (`ni nj nk`, `x y z`);
 * single or double precision; with or without an `iblank` array; either byte order.
 * Every hypothesis is checked against the file's exact size (the record sizes of a
 * Fortran file); one match is read, none is a `ReadError`, and several are a
 * `ReadError` naming them and the options that pick one.
 *
 * Each grid block becomes one cell block (hexahedra, quads or lines, by how many
 * of its dimensions exceed one), its points in PLOT3D's order (`i` fastest) after
 * the previous block's, and a cell region `block_<b>`. The dimensions go to
 * `field_data["plot3d:dims"]` (`(nblocks, 3)`, `nk = 1` in 2-D) and `iblank` to the
 * Int32 point data `plot3d:iblank` (overset files use negative values, so it is not
 * a `uint8`). A `.q` solution (`solution=`) adds the conserved variables as point
 * data (`density`, `momentum` as a vector, `energy`) and its four reference values
 * per block to `field_data["plot3d:q_reference"]` (mach, alpha, Reynolds number,
 * time); a `.f` function file (`function=`) adds `plot3d:f` (`(npoints, nvar)`).
 *
 * The writer emits the Fortran unformatted multi-block double-precision form
 * (little-endian, 4-byte markers), 3-D, with `iblank` when the mesh carries
 * `plot3d:iblank`; it needs `plot3d:dims` and refuses any other mesh by name.
 * The Python twin is `meshioplusplus.plot3d`.
 * @see doc/formats/plot3d.md
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief The choices a PLOT3D file does not state, `0`/`-1` meaning "infer". */
struct Plot3dReadOptions {
    std::string mSolution;  ///< a `.q` file to read as point data ("" = none)
    std::string mFunction;  ///< a `.f` function file to read as point data ("" = none)
    int mPrecision = 0;     ///< 4 or 8 bytes per real; 0 infers
    int mMultiblock = -1;   ///< 1 multi-block, 0 single; -1 infers
    int mIblank = -1;       ///< 1 with iblank, 0 without; -1 infers
    int mDimensions = 0;    ///< 2 or 3; 0 infers
};

/** @brief Read a PLOT3D grid (and optionally its solution and function files). */
MESHIOPLUSPLUS_API Mesh read_plot3d(const std::string& rPath, const Plot3dReadOptions& rOptions);

/** @brief Read a PLOT3D grid, inferring its variant. */
MESHIOPLUSPLUS_API Mesh read_plot3d(const std::string& rPath);

/** @brief Write a mesh carrying `plot3d:dims` as an unformatted multi-block grid. */
MESHIOPLUSPLUS_API void write_plot3d(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus
