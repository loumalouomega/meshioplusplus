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
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
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

}  // namespace detail
}  // namespace meshioplusplus
