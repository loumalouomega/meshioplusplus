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

// System includes
#include <stdexcept>
#include <string>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "tet_scope.hpp"

namespace meshioplusplus {
namespace detail {

void require_tet_only_blocks(const Mesh& rMesh, const char* pPrefix) {
    const std::string prefix(pPrefix);
    bool has_tet = false;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            throw std::invalid_argument(
                prefix +
                "operates on tet-only meshes, but the mesh contains a polyhedron cell block -- "
                "run convert_cells(mode='simplexify') first");
        if (cb.IsRagged())
            throw std::invalid_argument(
                prefix + "operates on tet-only meshes, but the mesh contains ragged cell block '" +
                std::string(cb.Type()) + "'");
        const CellType ct = cell_type_from_name(cb.Type());
        if (ct == CellType::Tetra) {
            has_tet = true;
            continue;
        }
        const int dim = cell_type_dimension(ct);
        if (dim == 3)
            throw std::invalid_argument(
                prefix + "operates on tet-only meshes, but the mesh contains 3D cell block '" +
                std::string(cb.Type()) +
                "' that is not linear tetra -- run convert_cells(mode='simplexify') first");
        throw std::invalid_argument(
            prefix + "operates on tet-only meshes, but the mesh contains non-3D cell block '" +
            std::string(cb.Type()) +
            "' alongside its tets (its nodes would dangle; drop it first, e.g. via split)");
    }
    if (!has_tet)
        throw std::invalid_argument(prefix + "requires at least one tetra cell block");
}

}  // namespace detail
}  // namespace meshioplusplus
