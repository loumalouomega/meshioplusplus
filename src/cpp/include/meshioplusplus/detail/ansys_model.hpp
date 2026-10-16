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
 * @file ansys_model.hpp
 * @brief The mesh of an Ansys MAPDL model: element categories, degenerate shapes
 *        and components, shared by the `.cdb` and `.rst` readers.
 *
 * An element becomes a cell by its routine's category (pymapdl-reader's table;
 * MESH200 by KEYOPT(1)): a brick with repeated nodes is a wedge, pyramid or
 * tetrahedron, a shell with K == L a triangle, and a missing midside node (node 0,
 * or a row cut short) is placed at its edge midpoint. The Python twin is `_build`
 * in `ansysInp/_ansysInp.py`.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/formats/ansysinp.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// How an element routine's nodes become a cell: a point, a line (quadratic when
/// its third node is set), a shell or plane, a degenerate brick, a native
/// tetrahedron, or a line whose extra nodes only orient it. `Skip` has no cell.
enum class AnsysCategory { Skip, Point, Line, Shell, Brick, Tet, LinearLine };

/// The category of element routine `Routine` (186 for SOLID186 ...).
MESHIOPLUSPLUS_API AnsysCategory ansys_category(int Routine);

/// MESH200's category, which its KEYOPT(1) chooses.
MESHIOPLUSPLUS_API AnsysCategory ansys_mesh200_category(int Keyopt1);

/// One element: its ET slot, attributes, number and node numbers.
struct AnsysElement {
    int mSlot = 0;
    std::int64_t mMat = 0, mReal = 0, mSecnum = 0, mId = 0;
    std::vector<std::int64_t> mNodes;
};

/// A NODE (`mNodes`) or ELEM component: its name and node or element numbers.
struct AnsysComponent {
    std::string mName;
    bool mNodes = false;
    std::vector<std::int64_t> mIds;
};

/// A model as a reader collects it.
struct AnsysModel {
    std::map<int, int> mRoutine;                // ET slot -> element routine
    std::map<int, std::map<int, int>> mKeyopt;  // ET slot -> KEYOPT k -> value
    std::vector<std::int64_t> mNodeIds;
    std::vector<double> mCoords;  // three per node
    std::vector<AnsysElement> mElements;
    std::vector<AnsysComponent> mComponents;
};

/**
 * @brief The mesh of `rModel`: cells grouped by type in first-seen order,
 *        `ansys:element`/`ansys:type`/`ansys:mat`/`ansys:real`/`ansys:secnum`
 *        cell data, components as point/cell regions (and in `rInfo`).
 *
 * `Label` prefixes messages (`"Ansys .cdb"`). An element whose type has no cell
 * is a ReadError, or skipped with a warning when `Lenient`. `rNodeIndex` receives
 * the node-number -> point-index map.
 */
MESHIOPLUSPLUS_API Mesh
ansys_build_mesh(const AnsysModel& rModel, bool Lenient, std::string_view Label, AnsysInfo& rInfo,
                 std::unordered_map<std::int64_t, std::int64_t>& rNodeIndex);

/// Where one element of an AnsysModel landed: its block and row in the built
/// mesh, and which element node became each of the cell's nodes (`mSlots[j]`
/// for cell node `j`). `mBlock` is -1 for an element that has no cell.
struct AnsysCellLocation {
    std::int64_t mBlock = -1;
    std::int64_t mRow = -1;
    std::vector<int> mSlots;
};

/**
 * @brief As above, also filling `rCells` (parallel to `rModel.mElements`) with
 *        each element's cell, so element results can be placed on it.
 */
MESHIOPLUSPLUS_API Mesh ansys_build_mesh(const AnsysModel& rModel, bool Lenient,
                                         std::string_view Label, AnsysInfo& rInfo,
                                         std::unordered_map<std::int64_t, std::int64_t>& rNodeIndex,
                                         std::vector<AnsysCellLocation>& rCells);

/**
 * @brief The element nodes of cell type `Type` under an element of `Category`:
 *        entry `c` is the cell node element node `c` takes (a degenerate brick
 *        repeats cell nodes); empty when the category has no such layout.
 */
MESHIOPLUSPLUS_API std::vector<int> ansys_cell_layout(std::string_view Type,
                                                      AnsysCategory Category);

/// The element routine a cell type is written as when nothing better is known
/// (185 for a hexahedron ...), or -1 when Ansys has no element for it.
MESHIOPLUSPLUS_API int ansys_default_routine(std::string_view Type);

/// The corner nodes of a cell type (8 for `"hexahedron20"`); 0 when unknown.
MESHIOPLUSPLUS_API std::size_t ansys_num_corners(std::string_view Type);

/// The Ansys element of every cell of a mesh: `mCellEtype[block][row]` is
/// `(routine, ET slot)` and `mSlotRoutine` maps each ET slot to its routine.
struct AnsysElementTypes {
    std::map<int, int> mSlotRoutine;
    std::vector<std::vector<std::pair<int, int>>> mCellEtype;
};

/**
 * @brief Choose the Ansys element of every cell: its `ansys:element` routine and
 *        `ansys:type` slot where the element's node layout fits the cell, else
 *        the default routine of its type in a fresh slot.
 * @param What Names the writer in messages (`"Ansys .cdb writer"`).
 * @throws WriteError for a cell type (or a ragged block) with no element type.
 */
MESHIOPLUSPLUS_API AnsysElementTypes ansys_element_types(const Mesh& rMesh, std::string_view What);

}  // namespace detail
}  // namespace meshioplusplus
