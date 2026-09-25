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
 * @file formats/abaqus_face.hpp
 * @brief Abaqus element face identifiers (`S1`..`S6`, `SPOS`/`SNEG`) and
 *        meshio++'s local facets, shared by the `.inp` and `.fil` readers.
 */

// System includes
#include <cctype>
#include <cstdint>
#include <string>

namespace meshioplusplus::detail {

/**
 * @brief Abaqus face identifier (`S1`..`S6`) -> the meshio++ local facet index
 * of `detail/cell_faces.hpp` / `cell_edges.hpp`.
 *
 * The two numberings genuinely differ, so this table is not the identity. It is
 * derived by matching node sets: Abaqus C3D8 `S1` is the 1-2-3-4 face, i.e.
 * local nodes {0,1,2,3}, which is meshio++'s face 4 (`{0,3,2,1}`) — same face,
 * different slot and winding. Getting this wrong yields a plausible-looking
 * side set pointing at the wrong faces, so each row is spelled out.
 *
 * Shell elements use `SPOS`/`SNEG` for their two sides rather than a facet;
 * there is no facet to name, so they map to 0 and 1 respectively and are
 * documented as such in doc/regions.md.
 *
 * @param rCellType The meshio++ cell type of the element.
 * @param rFace The Abaqus face identifier, upper-cased (e.g. `"S3"`).
 * @return The local facet index, or -1 when the pair has no mapping.
 */
inline int abaqus_face_index(const std::string& rCellType, const std::string& rFace) {
    if (rFace == "SPOS")
        return 0;
    if (rFace == "SNEG")
        return 1;
    if (rFace.size() < 2 || rFace[0] != 'S')
        return -1;
    int n = 0;
    for (std::size_t k = 1; k < rFace.size(); ++k) {
        if (!std::isdigit(static_cast<unsigned char>(rFace[k])))
            return -1;
        n = n * 10 + (rFace[k] - '0');
    }
    if (n < 1)
        return -1;
    const int s = n - 1;  // 0-based Abaqus face number

    // tetra: Abaqus S1=1-2-3, S2=1-2-4, S3=2-3-4, S4=1-3-4
    //        meshio++ 0={0,1,3} 1={1,2,3} 2={2,0,3} 3={0,2,1}
    static const int tetra[4] = {3, 0, 1, 2};
    // hexahedron: Abaqus S1=1-2-3-4, S2=5-8-7-6, S3=1-5-6-2,
    //                    S4=2-6-7-3,  S5=3-7-8-4, S6=4-8-5-1
    //             meshio++ 0={0,4,7,3} 1={1,2,6,5} 2={0,1,5,4}
    //                      3={3,7,6,2} 4={0,3,2,1} 5={4,5,6,7}
    static const int hexa[6] = {4, 5, 2, 1, 3, 0};
    // wedge: Abaqus S1=1-2-3, S2=4-5-6, S3=1-2-5-4, S4=2-3-6-5, S5=3-1-4-6
    //        meshio++ 0={0,2,1} 1={3,4,5} 2={0,1,4,3} 3={1,2,5,4} 4={2,0,3,5}
    static const int wedge[5] = {0, 1, 2, 3, 4};
    // 2-D elements: Abaqus numbers the edges 1-2, 2-3, ... in the same order
    // detail/cell_edges.hpp does.
    static const int identity[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};

    const int* table = nullptr;
    int count = 0;
    if (rCellType == "tetra" || rCellType == "tetra10") {
        table = tetra;
        count = 4;
    } else if (rCellType == "hexahedron" || rCellType == "hexahedron20") {
        table = hexa;
        count = 6;
    } else if (rCellType == "wedge" || rCellType == "wedge15") {
        table = wedge;
        count = 5;
    } else if (rCellType == "triangle" || rCellType == "triangle6") {
        table = identity;
        count = 3;
    } else if (rCellType == "quad" || rCellType == "quad8" || rCellType == "quad9") {
        table = identity;
        count = 4;
    }
    if (table == nullptr || s >= count)
        return -1;
    return table[s];
}

/// The inverse of `abaqus_face_index`, for the writer.
inline std::string abaqus_face_name(const std::string& rCellType, std::int64_t Facet) {
    for (int n = 1; n <= 6; ++n) {
        const std::string face = "S" + std::to_string(n);
        if (abaqus_face_index(rCellType, face) == static_cast<int>(Facet))
            return face;
    }
    return {};
}

}  // namespace meshioplusplus::detail
