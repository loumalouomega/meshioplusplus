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
 * @file face_mesh.hpp
 * @brief A mesh's volume cells re-expressed as a globally deduplicated face
 * list with owner/neighbour pairing.
 *
 * Two formats are *defined* in these terms rather than in terms of cells:
 *
 *  - **OpenFOAM `polyMesh`** is literally this structure on disk -- `faces`,
 *    `owner`, `neighbour` -- and never stores a cell's node list at all.
 *  - **CGNS `NGON_n` + `NFACE_n`**: `NGON_n` is the face list, `NFACE_n` gives
 *    each cell as a list of *signed* face ids, the sign meaning "traverse this
 *    face reversed".
 *
 * Building it twice would mean two implementations of the one genuinely
 * error-prone step (recognising that a hexahedron and a polyhedron meeting on a
 * face are meeting on *one* face), so `mCellFaces` is stored in CGNS's own
 * signed-1-based encoding: the CGNS writer copies a row with an element-id
 * offset and no per-entry branch, and OpenFOAM reads only the sign.
 *
 * ### Cells are numbered in a compact space, not the block-major global one
 *
 * `detail/cell_index.hpp`'s `block_bases` numbering counts *every* cell block,
 * including 2D ones. But every mesh `read_openfoam` produces carries its
 * boundary faces as `triangle`/`quad`/`polygon<N>` blocks, so in that numbering
 * "cell 8" is routinely a boundary quad rather than a cell. `mOwner`,
 * `mNeighbour` and `mCellFaceStart` therefore index a compact space containing
 * only the blocks that bound a volume, and `mCellToGlobal` bridges back to it
 * for `cell_data` lookups. Confusing the two is the bug that member exists to
 * prevent.
 *
 * ### Winding is repaired, never assumed
 *
 * `cell_faces.hpp`'s rows are outward-wound *on the reference element*. An
 * inverted (negative-Jacobian) hexahedron therefore yields six inward-pointing
 * normals, and OpenFOAM would report every one of its faces as misoriented. So
 * every cell -- tabulated and polyhedral alike -- goes through
 * `detail::orient_rings`, which rewinds and then flips on the enclosed signed
 * volume. `mNumFlipped` reports how often that mattered.
 *
 * ### Determinism
 *
 * Face ids are assigned by a **serial** first-seen sweep in ascending
 * (compact cell, local face) order over a buffer filled in parallel --
 * `surface.cpp`'s phase split. No hash-map iteration order reaches the output,
 * so the result is byte-identical across thread counts and across the three
 * mesh backends.
 *
 * See `doc/polyhedra.md`.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief A mesh's volume cells as a deduplicated face list plus the pairing
 * that tells you which cells each face separates.
 */
struct GlobalFaces {
    // ---- the faces themselves (OpenFOAM `faces`, CGNS `NGON_n`) ----

    /// CSR over the face node rings: global node ids, **corners only**, wound
    /// so that the normal points *out of* `mOwner[f]`.
    std::vector<std::int64_t> mFaceNodes;
    std::vector<std::int64_t> mFaceStart;  ///< `NumFaces() + 1` entries

    // ---- the pairing (OpenFOAM `owner` / `neighbour`) ----

    std::vector<std::int64_t> mOwner;      ///< compact cell id; never -1
    std::vector<std::int64_t> mNeighbour;  ///< compact cell id, or -1 = boundary

    // ---- the cells (CGNS `NFACE_n`) ----

    /// CSR over compact cell ids. Each entry is a **signed, 1-based** face id:
    /// `+(f+1)` means "traverse `Face(f)` as stored", `-(f+1)` means "reversed".
    /// This is CGNS `NFACE_n`'s own encoding, deliberately.
    ///
    /// Entry `k` of cell `c`'s row corresponds to face `k` of that cell as
    /// `detail::cell_rings` enumerated it, so "which face of which cell, and
    /// was it reversed" is answered by the position in the row plus the sign --
    /// no parallel side table to keep in step.
    std::vector<std::int64_t> mCellFaces;
    std::vector<std::int64_t> mCellFaceStart;  ///< `NumCells() + 1` entries

    /// Compact cell id -> the mesh's global (block-major) cell index, the
    /// numbering `detail/cell_index.hpp` and `cell_data` use. See the file
    /// header: these two numberings genuinely differ on any mesh with 2D blocks.
    std::vector<std::int64_t> mCellToGlobal;

    /// Block indices that contributed no cells: 1D/2D blocks, 1-level ragged
    /// (polygon) blocks, and 3D types with no `cell_faces` row. The **caller**
    /// decides which are acceptable -- a 2D boundary block is expected, a
    /// skipped 3D block is a silently dropped solid.
    std::vector<std::size_t> mNonCellBlocks;

    std::int64_t mNumFlipped = 0;       ///< cells whose winding needed repair
    std::int64_t mNumUnorientable = 0;  ///< cells that are not closed+orientable
    std::int64_t mNumNonManifold = 0;   ///< faces used by three or more cells

    /** @brief Number of distinct faces. */
    std::size_t NumFaces() const { return mFaceStart.empty() ? 0 : mFaceStart.size() - 1; }

    /** @brief Number of volume cells (the compact space). */
    std::size_t NumCells() const {
        return mCellFaceStart.empty() ? 0 : mCellFaceStart.size() - 1;
    }

    /** @brief Corner count of face @p Face. */
    std::size_t FaceSize(std::size_t Face) const {
        return static_cast<std::size_t>(mFaceStart[Face + 1] - mFaceStart[Face]);
    }

    /** @brief Global node ids of face @p Face, `FaceSize(Face)` of them. */
    const std::int64_t* Face(std::size_t Face) const {
        return mFaceNodes.data() + mFaceStart[Face];
    }

    /** @brief Number of faces bounding compact cell @p Cell. */
    std::size_t NumCellFaces(std::size_t Cell) const {
        return static_cast<std::size_t>(mCellFaceStart[Cell + 1] - mCellFaceStart[Cell]);
    }

    /** @brief Signed 1-based face ids of compact cell @p Cell. */
    const std::int64_t* CellFaces(std::size_t Cell) const {
        return mCellFaces.data() + mCellFaceStart[Cell];
    }
};

/**
 * @brief Build the global face list of @p rMesh's volume cells.
 *
 * Deterministic: face ids are assigned in ascending (compact cell, local face)
 * order, independent of thread count, mesh backend and hash order.
 *
 * Quadratic cells contribute their **corner** faces only (neither target format
 * has quadratic faces), so a `hexahedron20` mesh's mid-edge nodes end up
 * unreferenced by any face. That is reported by the caller rather than fixed
 * here: pruning them would renumber points for no benefit.
 */
MESHIOPLUSPLUS_API GlobalFaces build_global_faces(const Mesh& rMesh);

/**
 * @brief `build_global_faces` restricted to the cell blocks in @p rBlocks.
 *
 * Every other block is reported in `mNonCellBlocks` exactly as a 2D one would
 * be, so a caller cannot silently lose a solid: the distinction between "this
 * block has no volume" and "you did not ask for this block" is the caller's to
 * make, and CGNS's is the case that needs it. Writing a mesh that mixes
 * hexahedra with polyhedra emits ordinary `HEXA_8` sections for the former and
 * an `NGON_n`/`NFACE_n` pair for the latter, and the face list must then contain
 * the polyhedra's faces **only** -- a hexahedron's faces would be `NGON_n`
 * elements no `NFACE_n` cell ever references.
 *
 * Blocks keep their relative order; ids out of range are ignored.
 */
MESHIOPLUSPLUS_API GlobalFaces build_global_faces(const Mesh& rMesh,
                                                  const std::vector<std::size_t>& rBlocks);

/**
 * @brief Find a global face by its (unordered) corner ids.
 *
 * How a caller maps a 2D cell block -- an OpenFOAM boundary patch, say -- back
 * onto the face list. Built once, then queried per 2D cell.
 */
class MESHIOPLUSPLUS_API FaceLookup {
public:
    /** @brief Index @p rFaces. The reference need not outlive the lookup. */
    explicit FaceLookup(const GlobalFaces& rFaces);

    /**
     * @brief Look up the face whose corner set is exactly @p pIds.
     * @return the global face id, or -1 when no face has those corners.
     */
    std::int64_t Find(const std::int64_t* pIds, std::size_t N) const;

private:
    std::unordered_map<FacetKey, std::int64_t, FacetKeyHash> mMap;
};

}  // namespace detail
}  // namespace meshioplusplus
