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
 * @file polyhedron.hpp
 * @brief The geometric kernel for cells bounded by arbitrary polygonal faces:
 * volume, centroid, surface area, face normals and areas, and the winding
 * repair every one of those depends on.
 *
 * ### Why this is a separate header
 *
 * `detail/geometry.hpp` is inline `Vec3` arithmetic included nearly everywhere.
 * This kernel needs `mesh.hpp` and `cell_faces.hpp`, so folding it in there
 * would drag heavy includes into every consumer. Keeping the two apart also
 * keeps `geometry.hpp` byte-identical across this release.
 *
 * ### The impedance mismatch this resolves
 *
 * `cell_faces.hpp` describes a tabulated type's faces as **local** node indices
 * into a rectangular connectivity row. A polyhedron block has no row: it has
 * `NumFaces(cell)` faces of **global** node ids. `CellRings` is the single
 * shape both reduce to, so one measure kernel serves a `hexahedron` and a
 * `polyhedron12` alike, and every existing loop that wrote
 * `coords[face.mNodes[i]]` ports across unchanged.
 *
 * ### Winding: repaired, not required
 *
 * Faces *should* wind so their right-hand normal points out of the cell. Real
 * meshes do not always oblige -- this repository's own long-standing
 * `polyhedron5` test fixture has two faces traversing a shared edge in the
 * same direction -- so `orient_rings` fixes the winding per cell (BFS over the
 * faces' shared-edge dual, then a global flip if the enclosed volume came out
 * negative) rather than trusting or rejecting it. A face set that is not a
 * closed orientable surface is reported as `Unorientable`, never guessed at.
 *
 * ### Non-planar faces: the corner-average fan
 *
 * A face with four or more non-coplanar corners bounds no unique volume. This
 * kernel resolves it as the fan of triangles from each of the face's edges to
 * that face's **corner average**, each paired with the cell's own corner
 * average. Five independent reasons, because any one alone would read as
 * taste:
 *
 *  1. **The apex is a function of the FACE ALONE, so a shared face is
 *     triangulated identically from both sides whatever each cell's local node
 *     order happens to be.** Two cells meeting on a warped quad therefore agree
 *     on it exactly, and their boundary integrals cancel there. A fan about the
 *     ring's first node is a function of the cell's *storage*, so the two sides
 *     can pick different diagonals and the sum over a closed region drifts.
 *     (A structured generator numbers shared faces consistently and so happens
 *     to escape this -- which is precisely why the test for it uses two cells
 *     that store the shared ring differently, and why a mesh that looks fine
 *     can hide the defect until it is read from a file that does not.)
 *     `PolyhedronKernel.CellsSharingAWarpedFaceAgreeOnIt` pins it.
 *  2. **Invariance to the ring's start node**, which is the same property seen
 *     from one cell. OpenFOAM's `faces`, MED's `INN` and VTU's `faces` stream
 *     all choose that start node arbitrarily, so two readers of one mesh could
 *     otherwise legitimately disagree about a cell's volume. A quantity that
 *     moves under a rotation the file format does not constrain is a trap, not
 *     a rounding difference.
 *  3. **It is already the forced choice elsewhere.** `operations/gradient.cpp`'s
 *     Green-Gauss integration fans about the face centre because the corner
 *     average is the only apex whose value is known *exactly* for a linear
 *     field. Sharing the apex means gradient and stats measure the same
 *     surface instead of disagreeing ten lines apart.
 *  4. **It is OpenFOAM's own cell volume** (`primitiveMesh`'s face-centre to
 *     cell-centre pyramid decomposition), so a round-tripped case reports the
 *     volumes its solver would.
 *  5. **The convention is already in the repo**: `test_skin.cpp`'s
 *     `SkinFaceTables.OutwardWindingAndNodeCounts` asserts geometrically that a
 *     quad9's mid-face node sits at the face's corner average.
 *
 * See `doc/polyhedra.md`.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief One cell's boundary as node-indexed face rings -- the uniform shape a
 * polyhedron (`CellView::Face`, global ids) and a tabulated type
 * (`cell_faces.hpp`, local indices into a row) both reduce to.
 *
 * Corner nodes only: mid-edge and mid-face nodes are dropped, matching every
 * existing measure kernel in the repo. Reuse one instance across cells via
 * `Clear()`, which keeps the vectors' capacity.
 */
struct CellRings {
    std::vector<std::int64_t> mNodes;       ///< global node ids, first-seen order
    std::vector<std::uint32_t> mFaceNodes;  ///< local indices into `mNodes`
    std::vector<std::uint32_t> mFaceStart;  ///< `NumFaces() + 1` entries

    /** @brief Number of faces bounding the cell. */
    std::size_t NumFaces() const { return mFaceStart.empty() ? 0 : mFaceStart.size() - 1; }

    /** @brief Corner count of face @p Face. */
    std::size_t FaceSize(std::size_t Face) const { return mFaceStart[Face + 1] - mFaceStart[Face]; }

    /** @brief Local node indices of face @p Face, `FaceSize(Face)` of them. */
    const std::uint32_t* Face(std::size_t Face) const {
        return mFaceNodes.data() + mFaceStart[Face];
    }

    /** @brief Mutable form of `Face`, for the winding repair. */
    std::uint32_t* FaceMutable(std::size_t Face) { return mFaceNodes.data() + mFaceStart[Face]; }

    /** @brief Reset to empty, keeping the allocated capacity. */
    void Clear() {
        mNodes.clear();
        mFaceNodes.clear();
        mFaceStart.clear();
    }
};

/**
 * @brief Fill @p rOut with cell @p Cell's face rings and @p rCoords with the
 * matching coordinates, one entry per `rOut.mNodes` id.
 *
 * Works for a polyhedron block (faces come from `CellView::Face`) and for any
 * rectangular 3D type with a `cell_faces.hpp` row (faces come from the table,
 * corners only).
 *
 * @return false when the block has no face topology at all -- a 1D/2D block, a
 *   1-level ragged (polygon) block, or a 3D type with no `cell_faces` row (the
 *   3D Lagrange family). The caller then reports NaN and counts the cell rather
 *   than guessing, which is `compute_quality`'s standing convention.
 */
MESHIOPLUSPLUS_API bool cell_rings(const Mesh::CellView& rBlock, std::size_t Cell,
                                   const NDArray& rPoints, std::size_t PointDim, CellRings& rOut,
                                   std::vector<Vec3>& rCoords);

/** @brief Outcome of `orient_rings`. */
enum class RingOrientation : std::uint8_t {
    Consistent,   ///< already consistently wound outward; nothing changed
    Repaired,     ///< rewound (and/or globally flipped) to point outward
    Unorientable  ///< not a closed orientable surface; measures are undefined
};

/**
 * @brief Make every ring wind consistently, then flip them all if the enclosed
 * signed volume came out negative -- so on return the normals point *out*.
 *
 * BFS over the faces' shared-edge dual: two faces sharing an undirected edge
 * agree iff they traverse it in *opposite* directions. Returns `Unorientable`
 * (leaving @p rRings untouched) when some undirected edge is not used exactly
 * twice, which is what an open or non-manifold face set looks like.
 */
MESHIOPLUSPLUS_API RingOrientation orient_rings(CellRings& rRings, const Vec3* pCoords);

/** @brief What `poly_measure` computes about one cell. */
struct PolyMeasure {
    double mVolume = 0.0;  ///< signed; positive once the rings wind outward
    double mSurfaceArea = 0.0;
    Vec3 mCentroid{0.0, 0.0, 0.0};  ///< volume centroid
    bool mClosed = false;           ///< every undirected edge used exactly twice
    bool mOrientable = false;       ///< a consistent winding exists
};

/**
 * @brief Volume, volume centroid and surface area of the cell bounded by
 * @p rRings, by the corner-average fan documented at the top of this file.
 *
 * Call `orient_rings` first if the winding is not known to be outward; this
 * function measures whatever it is given, so inward-wound input yields a
 * negative volume rather than an error. Coordinates are recentred on the cell's
 * corner average before any arithmetic, which is what keeps the volume's
 * `sum(x . A) / 3` telescoping accurate for a mesh far from the origin.
 */
MESHIOPLUSPLUS_API PolyMeasure poly_measure(const CellRings& rRings, const Vec3* pCoords);

/**
 * @brief Newell area vector of one ring: `|result|` is the area, its direction
 * the normal. Translation-invariant and independent of which corner the ring
 * starts at, which is exactly why Newell rather than a fan is used here.
 */
MESHIOPLUSPLUS_API Vec3 polygon_area_vector(const Vec3* pCoords, const std::uint32_t* pRing,
                                            std::size_t N);

/** @brief Area of one ring (the magnitude of `polygon_area_vector`). */
MESHIOPLUSPLUS_API double polygon_area(const Vec3* pCoords, const std::uint32_t* pRing,
                                       std::size_t N);

/** @brief Corner average of one ring -- the fan apex, and the face "centre". */
MESHIOPLUSPLUS_API Vec3 polygon_centroid(const Vec3* pCoords, const std::uint32_t* pRing,
                                         std::size_t N);

/**
 * @brief `polygon_area_vector` for coordinates already stored in ring order,
 * i.e. with an implicit identity ring `0, 1, ... N-1`.
 *
 * This is the shape a rectangular 2D cell arrives in (`read_corner_coords`
 * hands back the corners in connectivity order), so it saves materializing an
 * identity index array once per cell in a hot loop.
 */
MESHIOPLUSPLUS_API Vec3 polygon_area_vector(const Vec3* pCoords, std::size_t N);

/** @brief `polygon_area` for coordinates already in ring order. */
MESHIOPLUSPLUS_API double polygon_area(const Vec3* pCoords, std::size_t N);

/**
 * @brief Signed volume of a tabulated 3D cell whose CORNER coordinates are
 * given in connectivity order, by the same corner-average fan as
 * `poly_measure`.
 *
 * For callers that already hold the coordinates and cannot go through
 * `cell_rings` -- `clean`, for instance, measures a cell built from *welded*
 * node representatives rather than from the mesh's own connectivity. Positive
 * when the cell is correctly oriented, since `cell_faces.hpp`'s rows are
 * outward-wound.
 *
 * @return NaN when @p Type has no `cell_faces` row.
 */
MESHIOPLUSPLUS_API double cell_volume_from_corners(const Vec3* pCoords, CellType Type);

/**
 * @brief The sorted corner ids of one facet, of any arity.
 *
 * `surface.cpp`'s existing key is a fixed `std::array<std::int64_t, 4>`, which
 * cannot hold a general polygonal face. This keeps that four-entry inline fast
 * path (so a triangle or quad still costs no allocation) and spills to the heap
 * beyond it.
 */
struct FacetKey {
    /** @brief An empty key, so a `FacetKey` can sit in a pre-sized buffer that
     *  a parallel pass then fills slot by slot (surface.cpp's phase split). */
    FacetKey() = default;

    /** @brief Build from @p N unsorted node ids. */
    FacetKey(const std::int64_t* pIds, std::size_t N);

    std::size_t Size() const { return mN; }
    const std::int64_t* Data() const { return mN <= kInline ? mInline.data() : mHeap.data(); }

    bool operator==(const FacetKey& rOther) const;

    static constexpr std::size_t kInline = 4;

private:
    std::size_t mN = 0;
    std::array<std::int64_t, kInline> mInline{};
    std::vector<std::int64_t> mHeap;
};

/** @brief Hash for `FacetKey`, for `unordered_map`/`unordered_set` keying. */
struct FacetKeyHash {
    std::size_t operator()(const FacetKey& rKey) const;
};

}  // namespace detail
}  // namespace meshioplusplus
