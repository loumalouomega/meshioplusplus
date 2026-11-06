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
 * @file detail/nastran_write_common.hpp
 * @brief The writer-side pieces the Nastran HDF5 and OP2 writers share: the
 * `nastran:card` number to card mapping, which cell block holds which card, the
 * cell type to card table, a series' signature, the GRID ids and the PCOMP plies.
 *
 * A **core-private** header: no installed header names it, and it adds nothing
 * to the API or the ABI. It needs no HDF5.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// The card a line or vertex takes from its `nastran:card` number, or nullptr when the
/// writers hold no cell type for it.
inline const char* nastran_line_card(std::int64_t Code) {
    switch (Code) {
        case 30:
            return "CROD";
        case 16:
            return "CONROD";
        case 37:
            return "CTUBE";
        case 6:
            return "CELAS1";
        case 2:
            return "CDAMP1";
        case 15:
            return "CONM2";
        case 10:
            return "CMASS1";
        case 11:
            return "CMASS2";
        case 14:
            return "CONM1";
        case 24:
            return "CBAR";
        case 54:
            return "CBEAM";
        case 26:
            return "CBUSH";
        case 3:
            return "CDAMP2";
        case 7:
            return "CELAS2";
        case 39:
            return "CVISC";
        case 52:
            return "PLOTEL";
        case 19:
            return "CGAP";
        default:
            return nullptr;
    }
}

/// Whether a cell block of @p rType may hold an element of card @p Card.
inline bool nastran_block_holds(const std::string& rType, std::string_view Card) {
    if (rType == "vertex")
        return Card == "CONM2" || Card == "CELAS1" || Card == "CDAMP1" || Card == "CELAS2" ||
               Card == "CDAMP2" || Card == "CMASS1" || Card == "CMASS2" || Card == "CONM1";
    if (rType == "line")
        return Card == "CROD" || Card == "CONROD" || Card == "CTUBE" || Card == "CELAS1" ||
               Card == "CDAMP1" || Card == "CBAR" || Card == "CBEAM" || Card == "CBUSH" ||
               Card == "CELAS2" || Card == "CDAMP2" || Card == "CVISC" || Card == "PLOTEL" ||
               Card == "CGAP" || Card == "CMASS1" || Card == "CMASS2";
    return false;
}

/// The card a `nastran:card` number names (its GEOM2 type), for the warnings.
std::string nastran_card_name(std::int64_t Code);

/// The card and default property id a cell block of a meshio type writes as.
struct NastranTypeSpec {
    const char* mType;
    const char* mCard;  // null: a line takes CROD, CTUBE or CONROD by its property
    std::size_t mNodes;
    std::int64_t mPid;  // the property id of a block with no 'nastran:pid'
};

/// The entry for @p rType, or nullptr when no card holds it.
const NastranTypeSpec* nastran_type_spec(const std::string& rType);

/// What every step of a series shares: the points' count and the cell blocks'
/// types and shapes.
using NastranSignature = std::vector<std::tuple<std::string, std::size_t, std::size_t>>;

NastranSignature nastran_signature(const Mesh& rMesh);

/// The GRID id of each point: `point_data["nastran:nid"]`, else 1 to n. Throws a
/// `WriteError` prefixed with @p pWho and ": ".
std::vector<std::int64_t> nastran_node_ids(const Mesh& rMesh, std::size_t N, const char* pWho);

/// The plies of each PCOMP id: `nastran:PCOMP_<id>:PLY`, `(n, 4)` rows of MID, T, THETA, SOUT.
using NastranPlies = std::map<std::int64_t, std::vector<std::array<double, 4>>>;

/// Throws a `WriteError` prefixed with @p pWho and ": ".
NastranPlies nastran_plies(const Mesh& rMesh, const char* pWho);

/// The names joined with ", ", in the set's order.
std::string nastran_join(const std::set<std::string>& rNames);

}  // namespace detail
}  // namespace meshioplusplus
