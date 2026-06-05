#pragma once
//
// VTK cell-type metadata shared by the VTU and VTK formats, ported from
// src/meshio/_vtk_common.py (meshio_to_vtk_type and the cell-ordering quirks).

#include <string>
#include <unordered_map>
#include <vector>

namespace meshio {

// meshio cell type -> VTK cell type id (inverse of vtk_to_meshio_type).
inline const std::unordered_map<std::string, int>& meshio_to_vtk_type() {
    static const std::unordered_map<std::string, int> m = {
        {"empty", 0},
        {"vertex", 1},
        {"line", 3},
        {"triangle", 5},
        {"polygon", 7},
        {"pixel", 8},
        {"quad", 9},
        {"tetra", 10},
        {"hexahedron", 12},
        {"wedge", 13},
        {"pyramid", 14},
        {"penta_prism", 15},
        {"hexa_prism", 16},
        {"line3", 21},
        {"triangle6", 22},
        {"quad8", 23},
        {"tetra10", 24},
        {"hexahedron20", 25},
        {"wedge15", 26},
        {"pyramid13", 27},
        {"quad9", 28},
        {"hexahedron27", 29},
        {"quad6", 30},
        {"wedge12", 31},
        {"wedge18", 32},
        {"hexahedron24", 33},
        {"triangle7", 34},
        {"line4", 35},
        {"polyhedron", 42},
        {"VTK_LAGRANGE_CURVE", 68},
        {"VTK_LAGRANGE_TRIANGLE", 69},
        {"VTK_LAGRANGE_QUADRILATERAL", 70},
        {"VTK_LAGRANGE_TETRAHEDRON", 71},
        {"VTK_LAGRANGE_HEXAHEDRON", 72},
        {"VTK_LAGRANGE_WEDGE", 73},
        {"VTK_LAGRANGE_PYRAMID", 74},
        {"VTK_BEZIER_CURVE", 75},
        {"VTK_BEZIER_TRIANGLE", 76},
        {"VTK_BEZIER_QUADRILATERAL", 77},
        {"VTK_BEZIER_TETRAHEDRON", 78},
        {"VTK_BEZIER_HEXAHEDRON", 79},
        {"VTK_BEZIER_WEDGE", 80},
        {"VTK_BEZIER_PYRAMID", 81},
    };
    return m;
}

// Node reordering when going from meshio to VTK. Only the linear wedge differs
// (meshio/gmsh prism ordering vs. vtkWedge). Returns empty for identity.
inline std::vector<int> meshio_to_vtk_order(const std::string& meshio_type) {
    if (meshio_type == "wedge") return {0, 2, 1, 3, 5, 4};
    return {};
}

}  // namespace meshio
