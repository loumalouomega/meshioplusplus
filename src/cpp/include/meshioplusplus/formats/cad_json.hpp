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
 * @file formats/cad_json.hpp
 * @brief Kratos CAD JSON (`*.cad.json`): the NURBS boundary representation
 * Kratos Multiphysics' IGA reads (`CadJsonInput`, `CadIoModeler`) and writes
 * (`CadJsonOutput`) -- trimmed NURBS faces, their boundary loops of trimming
 * curves, edges with 3-D curves and the coupling topology between faces, in
 * the JSON schema of Teschemacher et al. (AMSES 2018). The reader also takes
 * the `Surfaces`/`Lines` dialect of Kratos's `ImportNurbsSbmModeler`.
 *
 * ### The mesh
 *
 * Every knot span of every face becomes one Lagrange quadrilateral of the
 * face's higher degree (`quad`, `quad9`, else `VTK_LAGRANGE_QUADRILATERAL`),
 * its nodes the exact surface at the cell's lattice; `mSubdivisions` splits
 * each span further (by default only on a trimmed face, so its trim shows). A trimmed face keeps the cells whose parametric centre
 * lies inside its boundary loops (even-odd over all of its loops, sampled at
 * `mTrimSamples` points per knot span), so the trimmed edge of the mesh is
 * staircase-accurate to a cell while the trimming curves stay exact in the
 * side channel. Each edge with a 3-D curve becomes Lagrange line cells. Cell
 * data `cad:kind`, `cad:index`, `cad:id`, `cad:span` and `cad:trimmed` link
 * every cell back (see `BrepInfo`). Points are shared within a face or edge,
 * not between them.
 *
 * ### The side channel
 *
 * The overloads taking a `BrepInfo` return the faces, loops, edges and
 * vertices exactly, and `write_cad_json` with a non-empty `BrepInfo` writes
 * them back, in the layout of Kratos's own `PrettyPrintJsonString`; a file
 * meshio++ writes reads back to the same `BrepInfo` and rewrites to the same
 * bytes. Without one (the registry, a plain write) the writer turns every
 * `quad` cell into an untrimmed bilinear face and every `line` cell into an
 * edge, and drops other cells with a warning.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// How a B-rep is sampled into a mesh.
struct CadJsonReadOptions {
    /// Cells per knot span and direction, 1 to 64; 0 (automatic) keeps one per
    /// span except on a clipped face, refined to about 16 cells per direction.
    std::int32_t mSubdivisions = 0;
    /// Points per knot span a trimming curve is sampled at for the clip, 2 to 4096.
    std::int32_t mTrimSamples = 16;
    /// Whether edges with a 3-D curve become line cells.
    bool mEdges = true;
    /// Whether trimmed faces are clipped by their loops (false keeps every span).
    bool mTrim = true;
};

/// Read a CAD JSON file as its sampled mesh.
MESHIOPLUSPLUS_API Mesh read_cad_json(const std::string& rPath);

/// Read a CAD JSON file as its sampled mesh and its B-rep.
/// @throws ReadError on a malformed file or an inconsistent NURBS.
MESHIOPLUSPLUS_API Mesh read_cad_json(const std::string& rPath, BrepInfo& rInfo,
                                      const CadJsonReadOptions& rOptions = {});

/// Read only the B-rep of a CAD JSON file (no sampling).
MESHIOPLUSPLUS_API BrepInfo read_cad_json_info(const std::string& rPath);

/**
 * @brief Sample a B-rep into a mesh, as `read_cad_json` does.
 * @throws std::invalid_argument on an option out of range or an inconsistent
 * `BrepInfo` (an index out of range, a control net that does not match its
 * knots).
 */
MESHIOPLUSPLUS_API Mesh sample_brep(const BrepInfo& rInfo, const CadJsonReadOptions& rOptions = {});

/// Write a mesh as CAD JSON: `quad` cells as bilinear faces, `line` cells as edges.
MESHIOPLUSPLUS_API void write_cad_json(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Write CAD JSON from a B-rep; an empty @p rInfo writes @p rMesh as the
 * one-argument overload does.
 * @throws WriteError on an inconsistent `BrepInfo` or an unwritable path.
 */
MESHIOPLUSPLUS_API void write_cad_json(const std::string& rPath, const Mesh& rMesh,
                                       const BrepInfo& rInfo);

}  // namespace meshioplusplus
