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
 * @file unv.hpp
 * @brief I-DEAS Universal (.unv) C++ reader/writer — datasets 2411 (nodes)
 *        and 2412 (elements) only.
 *
 * A UNV file is a sequence of datasets, each delimited by a line containing
 * only `-1`, followed by a numeric dataset-id line and the dataset body.
 * Dataset **2411**: two-line node records (`label CS1 CS2 color` then a
 * coordinate line, Fortran `D`/`d` exponents normalized before parsing);
 * node labels are arbitrary integers, so a `label -> 0-based index` map is
 * built while reading. Dataset **2412**: a 6-integer record (`label fedesc
 * pid ... ... num_nodes`) selects the meshio++ type from the FE-descriptor
 * id (11/21->line, 22/24->line3, 41/81/91->triangle, 42/82/92->triangle6,
 * 44/84/94/122->quad, 45/85/95->quad8, 111->tetra, 118->tetra10,
 * 112->wedge, 115->hexahedron, 116->hexahedron20), followed by an extra
 * discarded 3-integer orientation record for beam descriptors (11/21/22/24)
 * — beam orientation is a genuinely lossy round-trip, always rewritten as
 * `0 0 0` on write — then the node-label records themselves.
 *
 * Parabolic (second-order) types use the Salome/Code-Aster mid-node
 * "sandwich" ordering (corner, mid-node, corner, mid-node, ...), converted
 * to meshio++'s "all corners then all edge nodes" convention via a fixed
 * permutation table per type (line3 `[0,2,1]`, triangle6
 * `[0,3,1,4,2,5]`, quad8 `[0,4,1,5,2,6,3,7]`, tetra10
 * `[0,4,1,5,2,6,7,8,9,3]`, hexahedron20 20-entry table) — applied directly
 * on read and inverted on write.
 *
 * Dataset 2467/2477 (permanent groups -> point_sets/cell_sets) and field
 * datasets (2414/55/56/57) are **not implemented in the C++ core**: the
 * reader throws ReadError as soon as it meets one, and the writer shim
 * routes any mesh carrying point_sets/cell_sets to the Python writer
 * (groups have no C++ writer path).
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as a UNV file (datasets 2411 + 2412 only).
 *
 * Emits node records (dataset 2411, labels = 1-based row index) and element
 * records (dataset 2412), choosing one canonical FE descriptor per meshio++
 * type (line->21, line3->24, triangle->91, triangle6->92, quad->94,
 * quad8->95, tetra->111, tetra10->118, wedge->112, hexahedron->115,
 * hexahedron20->116), applying the inverse sandwich permutation for
 * parabolic types, and always writing a placeholder `0 0 0` beam
 * orientation record for line/line3 elements.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError if the mesh carries `point_sets`/`cell_sets` (no
 *         dataset-2467 writer in C++ — the shim falls back to Python) or
 *         contains an unsupported cell type
 * @note reads `cell_data["unv:pid"]` for the per-element property id
 *       (defaults to `1` if absent).
 */
void write_unv(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a UNV file's node (2411) and element (2412) datasets.
 *
 * Splits the file into datasets on `-1` delimiter lines, builds a node
 * label->index map from dataset 2411, then decodes dataset 2412 element
 * records into typed cell blocks using the FE-descriptor table and the
 * sandwich-order permutation for parabolic types.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh
 * @throws ReadError if the file contains a 2467/2477 (permanent group) or a
 *         field (2414/55/56/57) dataset — the shim then falls back to the
 *         Python reader, which does support those (as point_sets/cell_sets).
 * @note cell_data key produced: `"unv:pid"` (element property id, dataset-
 *       2412 record-1 field 2).
 */
Mesh read_unv(const std::string& rPath);

}  // namespace meshioplusplus
