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
 * @file nastran_model.hpp
 * @brief The Nastran model half shared by the result readers (`nastran_h5`,
 *        `nastran_op2`): element cards to cell blocks, element ids to cells,
 *        property ids to regions.
 *
 * Both readers get their model as rows per card (`EID`, `PID`, the `G`
 * connectivity) and a GRID id -> point index map, and build the same mesh from
 * them: the linear or quadratic cell type of each element (a quadratic one only
 * when every mid-side node is given), `cell_data["nastran:eid"]` and
 * `["nastran:pid"]`, and one cell region per property id named
 * `<PTYPE>_<pid>` (`PID_<pid>` when the property card is unknown). An element
 * on a scalar point is dropped; one on an undefined GRID is an error. GRIDs
 * given in local coordinate systems (CP) are moved to the basic system, and
 * results in a GRID's output system (CD) can be rotated to basic.
 *
 * Python twin: `src/python/meshioplusplus/nastran/_model.py`.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace detail {

/** @brief One element card with a cell type: its linear and (optional) quadratic shape. */
struct NastranCardSpec {
    const char* mCard;
    const char* mLinear;
    std::size_t mLinearNodes;
    const char* mQuadratic;  ///< nullptr: no quadratic variant
    std::size_t mQuadraticNodes;
    const int* mPermutation;  ///< quadratic connectivity: conn[k] = G[perm[k]]; nullptr = identity
};

/** @brief The spec of an element card, or nullptr when it has no cell type. */
const NastranCardSpec* nastran_card_spec(std::string_view Card);

/** @brief The elements of one card: rows of `mWidth` grid ids each. */
struct NastranCardRows {
    std::string mCard;
    std::vector<std::int64_t> mEid;
    std::vector<std::int64_t> mPid;  ///< -1 where the card has none (CONM2)
    std::vector<std::int64_t> mNodes;
    std::size_t mWidth = 0;
};

/** @brief How elements map to the cells a Nastran model's blocks hold. */
struct NastranCells {
    /// EID -> global cell. A CONM2 has none: MSC accepts one sharing its id
    /// with a structural element, so it stays out.
    std::unordered_map<std::int64_t, std::size_t> mCellIndex;
    std::vector<std::size_t> mOffsets;  ///< first global cell of each block
    std::vector<std::size_t> mSizes;    ///< cells per block
    std::size_t mNumCells = 0;
};

/**
 * @brief Add the cell blocks, `nastran:eid`/`nastran:pid` and property regions
 *        of a model to @p rMesh, whose points are the GRIDs.
 * @param rCards the cards with a spec, in the order their blocks should take
 * @param rGridIndex GRID id -> point index
 * @param rScalarPoints SPOINT/EPOINT ids (an element on one is dropped)
 * @param rPtype property id -> property card name, for the region names
 * @param rWho the reader's name, prefixed to its warnings and errors
 * @throws ReadError when an element references an undefined GRID or no node
 */
NastranCells nastran_add_cells(Mesh& rMesh, const std::vector<NastranCardRows>& rCards,
                               const std::unordered_map<std::int64_t, std::size_t>& rGridIndex,
                               const std::unordered_set<std::int64_t>& rScalarPoints,
                               const std::map<std::int64_t, std::string>& rPtype,
                               const std::string& rWho);

/** @brief The arrays nastran_add_cells makes for the blocks it built: the
 *         `nastran:card`, `nastran:c1` and `nastran:c2` of every block and, when some
 *         block holds a CBAR (`mHasBar`), the CBAR's `nastran:x`, `nastran:g0`,
 *         `nastran:pin` and `nastran:w`; when some block holds a CELAS2 (`mHasSpring`),
 *         its `nastran:k`, `nastran:ge` and `nastran:s`; when some block holds a CDAMP2
 *         (`mHasDamper`), its `nastran:b`. */
struct NastranCellIds {
    std::vector<NDArray> mCard;  ///< the card's type number (its GEOM2 record's second key)
    std::vector<NDArray> mC1;    ///< a spring or damper's component at its first end
    std::vector<NDArray> mC2;    ///< its component at its second end (0 at a vertex)
    bool mHasBar = false;        ///< some block holds a CBAR: the four arrays below are filled
    std::vector<NDArray> mX;     ///< a CBAR's orientation vector X1..X3, `(n, 3)`
    std::vector<NDArray> mG0;    ///< a CBAR's G0 (0 with an orientation vector), `(n,)`
    std::vector<NDArray> mPin;   ///< a CBAR's pins PA and PB, `(n, 2)`
    std::vector<NDArray> mW;     ///< a CBAR's offsets W1A..W3A and W1B..W3B, `(n, 6)`
    bool mHasSpring = false;     ///< some block holds a CELAS2: the three arrays below are filled
    std::vector<NDArray> mK;     ///< a CELAS2's stiffness K, `(n,)`
    std::vector<NDArray> mGe;    ///< a CELAS2's damping coefficient GE, `(n,)`
    std::vector<NDArray> mS;     ///< a CELAS2's stress coefficient S, `(n,)`
    bool mHasDamper = false;     ///< some block holds a CDAMP2: the array below is filled
    std::vector<NDArray> mB;     ///< a CDAMP2's damping B, `(n,)`
};

/**
 * @brief nastran_add_cells, which also reports each block's card, spring components
 *        and CBAR fields.
 * @param rComps per card of @p rCards, two `int64` per element (its components
 *               at the first and at the second end), or empty for none
 * @param rBar per card of @p rCards, twelve doubles per element: for a CBAR, its
 *             orientation vector X1..X3, its G0 (0 with a vector), its pins PA and PB,
 *             and its offsets W1A..W3A and W1B..W3B; for a CELAS2, its K, GE and S in
 *             the first three; for a CDAMP2, its B in the first; empty for none
 * @param pIds receives the blocks' arrays, which the caller adds to @p rMesh; nullptr
 *        leaves them out
 */
NastranCells nastran_add_cells(Mesh& rMesh, const std::vector<NastranCardRows>& rCards,
                               const std::unordered_map<std::int64_t, std::size_t>& rGridIndex,
                               const std::unordered_set<std::int64_t>& rScalarPoints,
                               const std::map<std::int64_t, std::string>& rPtype,
                               const std::string& rWho,
                               const std::vector<std::vector<std::int64_t>>& rComps,
                               const std::vector<std::vector<double>>& rBar, NastranCellIds* pIds);

/** @brief One coordinate system card: CORD1R/C/S (three GRIDs) or CORD2R/C/S. */
struct NastranCoordCard {
    std::int64_t mCid = 0;
    int mType = 1;          ///< 1 rectangular, 2 cylindrical, 3 spherical
    bool mByGrids = false;  ///< CORD1*: `mGrids`; CORD2*: `mRid` and `mAbc`
    std::int64_t mRid = 0;
    double mAbc[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};  ///< A (origin), B (on +z), C (in the xz plane)
    std::int64_t mGrids[3] = {0, 0, 0};            ///< origin, +z and xz-plane GRIDs
};

/**
 * @brief Resolved coordinate systems: each an origin and axes in the basic
 *        system, and its type.
 *
 * Local coordinates are rectangular (x, y, z), cylindrical (r, theta, z) or
 * spherical (r, theta, phi), angles in degrees, theta of a spherical system
 * measured from its z axis.
 */
class NastranCoordSystems {
public:
    struct System {
        int mType = 1;
        double mOrigin[3] = {0, 0, 0};
        double mAxes[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};  ///< rows: x, y, z axis in basic
    };

    NastranCoordSystems();
    bool Has(std::int64_t Cid) const { return mSystems.count(Cid) != 0; }
    /** @brief The resolved system @p Cid, which must exist. */
    const System& Get(std::int64_t Cid) const { return mSystems.at(Cid); }
    /** @brief Local coordinates of system @p Cid (which must exist) to basic. */
    void ToBasic(std::int64_t Cid, const double* pLocal, double* pBasic) const;
    /**
     * @brief A vector given in the components of system @p Cid at the point
     *        @p pBasicPoint (the local basis of a cylindrical or spherical
     *        system depends on it) to basic components, in place.
     */
    void VectorToBasic(std::int64_t Cid, const double* pBasicPoint, double* pVector) const;
    void Add(std::int64_t Cid, const System& rSystem) { mSystems[Cid] = rSystem; }

private:
    std::map<std::int64_t, System> mSystems;
};

/**
 * @brief Move GRID points to the basic system and keep their frames.
 *
 * Resolves @p rCards (in any order: a CORD2 through its reference system, a
 * CORD1 through its GRIDs, which may themselves be in local systems), then
 * turns every point with CP != 0 into basic coordinates. `nastran:cp` and
 * `nastran:cd` are added as point data when any value is non-zero. A GRID
 * whose system cannot be resolved (undefined, circular or degenerate) keeps
 * its coordinates as written, with a warning.
 *
 * @param rMesh a mesh whose points are the GRIDs, as written, in @p rIds order
 * @return the resolved systems, for the CD components of nodal results
 */
NastranCoordSystems nastran_apply_frames(Mesh& rMesh, const std::vector<NastranCoordCard>& rCards,
                                         const std::vector<std::int64_t>& rIds,
                                         const std::vector<std::int64_t>& rCp,
                                         const std::vector<std::int64_t>& rCd,
                                         const std::string& rWho);

/**
 * @brief Rotate a point's result vectors from its CD system to basic, in place.
 *
 * @param pValues `Count` consecutive triplets (e.g. translations then rotations)
 * @return false (values untouched) when @p Cd is 0 or not a resolved system
 */
bool nastran_rotate_to_basic(const NastranCoordSystems& rSystems, std::int64_t Cd,
                             const double* pBasicPoint, double* pValues, std::size_t Count);

/**
 * @brief Read a bulk-data deck as `read_nastran` does, also returning the GRID
 *        ids in point order and the element ids in global cell order (the OP2
 *        reader's sibling-deck route).
 */
Mesh nastran_read_deck(const std::string& rPath, std::vector<std::int64_t>& rGridIds,
                       std::vector<std::int64_t>& rCellIds);

}  // namespace detail
}  // namespace meshioplusplus
