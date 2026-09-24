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
// System includes
#include <cctype>
#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/abaqus_types.hpp"
#include "meshioplusplus/cell_type.hpp"

namespace meshioplusplus {
namespace detail {

// (abaqus type, meshio type) in source order; the meshio->abaqus inverse keeps
// the last entry per meshio type (matching the Python dict comprehension).
const std::vector<std::pair<std::string, std::string>>& abaqus_type_table() {
    static const std::vector<std::pair<std::string, std::string>> t = {

        {"T2D2", "line"},
        {"T2D2H", "line"},
        {"T2D3", "line3"},
        {"T2D3H", "line3"},
        {"T3D2", "line"},
        {"T3D2H", "line"},
        {"T3D3", "line3"},
        {"T3D3H", "line3"},
        {"B21", "line"},
        {"B21H", "line"},
        {"B22", "line3"},
        {"B22H", "line3"},
        {"B31", "line"},
        {"B31H", "line"},
        {"B32", "line3"},
        {"B32H", "line3"},
        {"B33", "line3"},
        {"B33H", "line3"},
        {"CPS4", "quad"},
        {"CPS4R", "quad"},
        {"S4", "quad"},
        {"S4R", "quad"},
        {"S4RS", "quad"},
        {"S4RSW", "quad"},
        {"S4R5", "quad"},
        {"S8R", "quad8"},
        {"S8R5", "quad8"},
        {"S9R5", "quad9"},
        {"CPS3", "triangle"},
        {"STRI3", "triangle"},
        {"S3", "triangle"},
        {"S3R", "triangle"},
        {"S3RS", "triangle"},
        {"R3D3", "triangle"},
        {"STRI65", "triangle6"},
        {"C3D8", "hexahedron"},
        {"C3D8H", "hexahedron"},
        {"C3D8I", "hexahedron"},
        {"C3D8IH", "hexahedron"},
        {"C3D8R", "hexahedron"},
        {"C3D8RH", "hexahedron"},
        {"C3D20", "hexahedron20"},
        {"C3D20H", "hexahedron20"},
        {"C3D20R", "hexahedron20"},
        {"C3D20RH", "hexahedron20"},
        // C3D4H before C3D4: the writer's inverse keeps the last name per type.
        {"C3D4H", "tetra"},
        {"C3D4", "tetra"},
        {"C3D10", "tetra10"},
        {"C3D10H", "tetra10"},
        {"C3D10I", "tetra10"},
        {"C3D10M", "tetra10"},
        {"C3D10MH", "tetra10"},
        {"C3D6", "wedge"},
        {"C3D15", "wedge15"},
        {"CAX4P", "quad"},
        {"CPE6", "triangle6"},
    };
    return t;
}

namespace {

bool abq_starts_with(std::string_view Name, std::string_view Prefix) {
    return Name.substr(0, Prefix.size()) == Prefix;
}

// "S" followed by a digit: S3, S4R, S8R5, S9R5 (not SAX1, SC8R, SPRINGA).
bool abq_is_shell(std::string_view Name) {
    return (Name.size() > 1 && Name[0] == 'S' &&
            std::isdigit(static_cast<unsigned char>(Name[1]))) ||
           abq_starts_with(Name, "STRI");
}

}  // namespace

std::string abaqus_cell_type(std::string_view Name, std::size_t NodeCount) {
    static const std::unordered_map<std::string, std::string> table = [] {
        std::unordered_map<std::string, std::string> m;
        for (const auto& [abq, type] : abaqus_type_table())
            m.emplace(abq, type);
        return m;
    }();
    const auto it = table.find(std::string(Name));
    if (it != table.end() &&
        cell_type_num_nodes(cell_type_from_name(it->second)) == static_cast<int>(NodeCount))
        return it->second;

    static const char* const solid[] = {"C3D", "DC3D", "AC3D", "COH3D", "SC", "CCL"};
    static const char* const planar[] = {"CPE",   "CPS",   "CAX", "CGAX", "DC2D",  "DCAX",
                                         "COH2D", "COHAX", "M3D", "R3D",  "SFM3D", "CPEG"};
    static const char* const line[] = {"T2D",  "T3D",  "B2",   "B3",  "PIPE", "R2D",
                                       "RB2D", "RB3D", "DC1D", "SAX", "FRAME"};
    auto any_prefix = [&](const auto& rList) {
        for (const char* p : rList)
            if (abq_starts_with(Name, p))
                return true;
        return false;
    };
    if (any_prefix(solid)) {
        switch (NodeCount) {
            case 4:
                return "tetra";
            case 5:
                return "pyramid";
            case 6:
                return "wedge";
            case 8:
                return "hexahedron";
            case 10:
                return "tetra10";
            case 13:
                return "pyramid13";
            case 15:
                return "wedge15";
            case 20:
                return "hexahedron20";
            default:
                return {};
        }
    }
    if (any_prefix(planar) || abq_is_shell(Name)) {
        switch (NodeCount) {
            case 3:
                return "triangle";
            case 4:
                return "quad";
            case 6:
                return "triangle6";
            case 8:
                return "quad8";
            case 9:
                return "quad9";
            default:
                return {};
        }
    }
    if (any_prefix(line)) {
        if (NodeCount == 2)
            return "line";
        if (NodeCount == 3)
            return "line3";
    }
    return {};
}

}  // namespace detail
}  // namespace meshioplusplus
