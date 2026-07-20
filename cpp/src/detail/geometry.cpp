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

// Project includes
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {
namespace detail {

Vec3 read_point(const NDArray& rPoints, std::size_t pointDim, std::int64_t nodeId) {
    Vec3 p = {0.0, 0.0, 0.0};
    const std::size_t base = static_cast<std::size_t>(nodeId) * pointDim;
    const std::size_t k = pointDim < 3 ? pointDim : 3;
    for (std::size_t c = 0; c < k; ++c)
        p[c] = read_double(rPoints, base + c);
    return p;
}

void read_corner_coords(const NDArray& rPoints, std::size_t pointDim, const NDArray& rConn,
                        std::size_t rowOffset, std::size_t n, std::vector<Vec3>& rOut) {
    rOut.clear();
    rOut.reserve(n);
    for (std::size_t k = 0; k < n; ++k)
        rOut.push_back(read_point(rPoints, pointDim, read_int(rConn, rowOffset + k)));
}

int cell_corner_count(CellType type) {
    switch (type) {
        case CellType::Vertex:
            return 1;
        case CellType::Line:
        case CellType::Line3:
        case CellType::Line4:
        case CellType::Line5:
        case CellType::Line6:
        case CellType::Line7:
        case CellType::Line8:
        case CellType::Line9:
        case CellType::Line10:
        case CellType::Line11:
            return 2;
        case CellType::Triangle:
        case CellType::Triangle6:
        case CellType::Triangle10:
        case CellType::Triangle15:
        case CellType::Triangle21:
        case CellType::Triangle28:
        case CellType::Triangle36:
        case CellType::Triangle45:
        case CellType::Triangle55:
        case CellType::Triangle66:
            return 3;
        case CellType::Quad:
        case CellType::Quad8:
        case CellType::Quad9:
        case CellType::Quad16:
        case CellType::Quad25:
        case CellType::Quad36:
        case CellType::Quad49:
        case CellType::Quad64:
        case CellType::Quad81:
        case CellType::Quad100:
        case CellType::Quad121:
            return 4;
        case CellType::Tetra:
        case CellType::Tetra10:
        case CellType::Tetra20:
        case CellType::Tetra35:
        case CellType::Tetra56:
        case CellType::Tetra84:
        case CellType::Tetra120:
        case CellType::Tetra165:
        case CellType::Tetra220:
        case CellType::Tetra286:
            return 4;
        case CellType::Hexahedron:
        case CellType::Hexahedron20:
        case CellType::Hexahedron24:
        case CellType::Hexahedron27:
        case CellType::Hexahedron64:
        case CellType::Hexahedron125:
        case CellType::Hexahedron216:
        case CellType::Hexahedron343:
        case CellType::Hexahedron512:
        case CellType::Hexahedron729:
        case CellType::Hexahedron1000:
        case CellType::Hexahedron1331:
            return 8;
        case CellType::Wedge:
        case CellType::Wedge15:
        case CellType::Wedge18:
        case CellType::Wedge40:
        case CellType::Wedge75:
        case CellType::Wedge126:
        case CellType::Wedge196:
        case CellType::Wedge288:
        case CellType::Wedge405:
        case CellType::Wedge550:
            return 6;
        case CellType::Pyramid:
        case CellType::Pyramid13:
        case CellType::Pyramid14:
            return 5;
        default:  // Polygon, Polyhedron, VtkLagrange*, Custom
            return 0;
    }
}

}  // namespace detail
}  // namespace meshioplusplus
