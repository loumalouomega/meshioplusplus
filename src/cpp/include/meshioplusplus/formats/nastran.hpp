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
 * @file nastran.hpp
 * @brief MSC/NX Nastran and Altair OptiStruct bulk-data (.bdf/.fem/.nas)
 *        reader and writer.
 *
 * Bulk data is card text (`GRID`, `CTRIA3`, `CTETRA`, `CHEXA`, …) between a
 * `BEGIN BULK` line and `ENDDATA`, in small-field (10 x 8 columns),
 * large-field (8 + 4x16 + 8 columns, keyword suffixed `*`) or free
 * (comma-separated) layout, with `+`/`*` or implicit blank continuations.
 * The reader accepts any such deck; the Python reference
 * (`nastran/_nastran.py`) is its twin and both give the same Mesh.
 *
 * HyperMesh/OptiStruct decks keep their components in comment cards:
 * `$HMMOVE <id>` followed by `$` lines of element ids (`a THRU b` ranges)
 * and `$HMNAME COMP <id>"name"`, which may come after the elements. Each
 * component becomes a cell Region tagged with its id. OptiStruct
 * `SET,<id>,GRID|ELEM,LIST,…` cards become point or cell Regions, named
 * by `$HMSET`. Cards meshio++ does not read are skipped; the ones that are
 * not properties, materials, loads, constraints, coordinate systems or
 * solution parameters are named in one warning (optimization, contact,
 * unsupported elements such as `CONM2` or `RBE2`).
 *
 * `CTETRA`/`CPYRA`/`CPENTA`/`CHEXA` become their 10/13/15/20-node
 * quadratic forms when a card lists that many nodes. Node orders:
 * `CTRIAX6`/`CTRAX6` list corner, mid, corner… (to meshio `[0,2,4,1,3,5]`);
 * `hexahedron20` and `wedge15` put the vertical mid-edges before the top
 * ring (involutions `[0..11,16..19,12..15]`, `[0..8,12,13,14,9,10,11]`).
 * Other card fields past the nodes (THETA, ZOFFS, orientation vectors) are
 * ignored. See doc/formats/nastran.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to a Nastran bulk-data file (large-field `GRID*`,
 *        small-field element cards).
 *
 * `GRID*` coordinates use the shortest 16-character string that round-trips.
 * `nastran:ref` point/cell data fills the CP/PID field (0 is written blank).
 * Disjoint cell regions are written as HyperMesh `$HMMOVE`/`$HMNAME COMP`
 * comment blocks; point, side and overlapping regions are dropped with a
 * warning. The first line is the comment `meshioplusplus-cpp-nastran`,
 * which releases before 16.1 required to read the file back.
 *
 * @param rPath filesystem path to the .bdf/.fem/.nas file to create/overwrite
 * @param rMesh the mesh to write
 * @throws WriteError on an unsupported cell type
 */
MESHIOPLUSPLUS_API void write_nastran(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a Nastran/OptiStruct bulk-data file into a Mesh.
 *
 * Produces `nastran:ref` point data (GRID CP field) and cell data (element
 * PID field) when any card fills that field, blanks reading as 0; HyperMesh
 * components and OptiStruct SETs as Regions.
 *
 * @param rPath filesystem path to the .bdf/.fem/.nas file to read
 * @return the read Mesh
 * @throws ReadError when there is no `BEGIN BULK`, on a malformed field, on
 *         an element with a wrong node count, or on an undefined grid
 */
MESHIOPLUSPLUS_API Mesh read_nastran(const std::string& rPath);

}  // namespace meshioplusplus
