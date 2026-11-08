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
#include <initializer_list>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

#include "typed_view.hpp"

namespace meshioplusplus {
namespace detail {

/// The card a line or vertex takes from its `nastran:card` number, or nullptr when the
/// writers hold no cell type for it.
/// Throws a `WriteError` "<pWho>: <message>".
[[noreturn]] void nastran_write_fail(const char* pWho, const std::string& rMessage);

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

/// The id in "<PTYPE>_<digits>" when the name starts with @p rPtype and "_",
/// else -1 (a region naming a property, e.g. `PTUBE_7`).
std::int64_t nastran_ptype_pid(const std::string& rName, const std::string& rPtype);

/// Where @p Name sits in a writer's card table, or `N` when it is not there.
template <class TCard, std::size_t N>
std::size_t nastran_card_index(const TCard (&rCards)[N], std::string_view Name) {
    for (std::size_t i = 0; i < N; ++i)
        if (Name == rCards[i].mName)
            return i;
    return N;
}

/// The mesh block (and cells) behind one element row, as both writers build it.
struct NastranWriteBlock {
    std::size_t mIndex = 0;  // the mesh's block
    std::size_t mCount = 0;
    std::size_t mNodes = 0;
    std::vector<std::int64_t> mConn;  // 0-based, mCount x mNodes
    std::vector<std::int64_t> mEids;
    std::vector<std::int64_t> mPids;
    std::vector<std::size_t> mCards;   // card index per cell
    std::vector<std::int64_t> mComps;  // two per cell: a spring's or damper's components
    std::vector<double> mBar;          // the writer's bar width per cell: a CBAR's fields, a
                                       // CELAS2's K, GE and S, a CDAMP2's B, a CONM1's terms
};

/// One element row: its id and the block and cell it comes from.
struct NastranWriteRow {
    std::int64_t mEid;
    std::size_t mBlock;  // into the writer's blocks
    std::size_t mRow;
};

/// `{(card, id): {member: value}}`.
using NastranValues = std::map<std::pair<std::string, std::int64_t>, std::map<std::string, double>>;

/// `{(card, id): {member: value}}` from the `nastran:<CARD>_<id>:<MEMBER>` field
/// data, for the cards and members of @p rCards (each with `mCard` and
/// `mMembers`, a member with `mName` and `mInt`); keys that name no card or
/// member there are left alone. A key holding more than one value throws a
/// `WriteError` prefixed with @p pWho and ": ".
template <class TCards>
NastranValues nastran_property_values(const Mesh& rMesh, const TCards& rCards, const char* pWho) {
    NastranValues out;
    for (const std::string& name : rMesh.FieldDataNames()) {
        if (name.rfind("nastran:", 0) != 0)
            continue;
        const std::size_t colon = name.rfind(':');
        const std::size_t under = name.rfind('_', colon);
        if (colon <= 8 || under == std::string::npos || under < 8 || under + 1 >= colon)
            continue;
        const std::string card = name.substr(8, under - 8);
        const std::string member = name.substr(colon + 1);
        const std::string digits = name.substr(under + 1, colon - under - 1);
        if (digits.find_first_not_of("0123456789") != std::string::npos || digits.size() > 15)
            continue;
        for (const auto& c : rCards) {
            if (card != c.mCard)
                continue;
            for (const auto& m : c.mMembers)
                if (member == m.mName) {
                    const NDArray& a = rMesh.FieldData(name);
                    if (a.Size() != 1)
                        nastran_write_fail(pWho, "'" + name + "' must hold one value");
                    out[{card, std::stoll(digits)}][member] =
                        m.mInt ? static_cast<double>(Int64View(a).Data()[0])
                               : DoubleView(a).Data()[0];
                }
        }
    }
    return out;
}

/// `nastran:<CARD>_<id>:<MEMBER>` field data of one card (PBEAM, PBUSH, ...), by id,
/// each member's values as doubles.
using NastranArrays = std::map<std::int64_t, std::map<std::string, std::vector<double>>>;

/// An int64 or double array's values as doubles.
std::vector<double> nastran_as_doubles(const NDArray& rA);

NastranArrays nastran_named_arrays(const Mesh& rMesh, const std::string& rCard);

/// A spring's or damper's components into @p rBlk's `mComps`, two per cell:
/// `nastran:c1` and `nastran:c2` of mesh block @p B, else 1 at each end; a
/// vertex's second end is ground, so 0.
void nastran_read_components(const Mesh& rMesh, std::size_t B, const char* pWho,
                             NastranWriteBlock& rBlk);

/// Cell data @p pName of mesh block @p B, @p Width values per cell, into columns
/// @p Lo.. of the rows of @p rBlk's `mBar` (@p BarWidth per cell) whose card is
/// one of @p Cards: a row takes only its own card's fields.
void nastran_read_field(const Mesh& rMesh, std::size_t B, const char* pName,
                        std::initializer_list<std::size_t> Cards, std::size_t Lo, std::size_t Width,
                        std::size_t BarWidth, const char* pWho, NastranWriteBlock& rBlk);

/// Drops the cells of @p rBlk whose @p rKeep flag is clear (@p Kept remain).
void nastran_compact_block(const std::vector<bool>& rKeep, std::size_t Kept, std::size_t BarWidth,
                           NastranWriteBlock& rBlk);

/// Each card's rows (@p NumCards of them), by element id; equal ids keep block order.
std::vector<std::vector<NastranWriteRow>> nastran_rows_by_card(
    const std::vector<NastranWriteBlock>& rBlocks, std::size_t NumCards);

/// The ids of the regions named "<rPtype>_<id>" (the PTUBE properties).
std::set<std::int64_t> nastran_region_pids(const Mesh& rMesh, const std::string& rPtype);

/// Fills @p rBlk from mesh block @p B, which @p rSpec writes: its count and nodes
/// per cell, the 0-based connectivity checked against @p NumPoints,
/// `nastran:eid` (or ids from @p rSequential on, which advances past the block
/// either way) and `nastran:pid` (or the spec's default). Throws a `WriteError`
/// prefixed with @p pWho and ": ".
void nastran_fill_block(const Mesh& rMesh, std::size_t B, const NastranTypeSpec& rSpec,
                        std::size_t NumPoints, std::int64_t& rSequential, const char* pWho,
                        NastranWriteBlock& rBlk);

}  // namespace detail
}  // namespace meshioplusplus
