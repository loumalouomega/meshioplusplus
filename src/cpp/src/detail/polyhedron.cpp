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
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <utility>

// Project includes
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/cell_type.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// A directed ring edge, keyed undirected. The BFS below asks "which other face
// uses this same undirected edge, and in which direction".
using PolyEdge = std::pair<std::int64_t, std::int64_t>;

struct PolyEdgeHash {
    std::size_t operator()(const PolyEdge& rE) const {
        // The surface.cpp mixing constant; the pair is already small ints.
        const std::size_t a = static_cast<std::size_t>(rE.first);
        const std::size_t b = static_cast<std::size_t>(rE.second);
        return a * 0x9E3779B97F4A7C15ULL + b;
    }
};

PolyEdge poly_undirected(std::uint32_t a, std::uint32_t b) {
    return a < b ? PolyEdge{a, b} : PolyEdge{b, a};
}

// Append `id` to rOut.mNodes if new, and return its local index. Cells have a
// handful of nodes, so the linear scan beats a hash map comfortably.
std::uint32_t poly_intern(CellRings& rOut, std::vector<Vec3>& rCoords, std::int64_t Id,
                          const NDArray& rPoints, std::size_t PointDim) {
    for (std::size_t i = 0; i < rOut.mNodes.size(); ++i)
        if (rOut.mNodes[i] == Id)
            return static_cast<std::uint32_t>(i);
    rOut.mNodes.push_back(Id);
    rCoords.push_back(read_point(rPoints, PointDim, Id));
    return static_cast<std::uint32_t>(rOut.mNodes.size() - 1);
}

}  // namespace

bool cell_rings(const Mesh::CellView& rBlock, std::size_t Cell, const NDArray& rPoints,
                std::size_t PointDim, CellRings& rOut, std::vector<Vec3>& rCoords) {
    rOut.Clear();
    rCoords.clear();

    if (rBlock.IsPolyhedron()) {
        const std::size_t nf = rBlock.NumFaces(Cell);
        if (nf == 0)
            return false;
        rOut.mFaceStart.push_back(0);
        for (std::size_t f = 0; f < nf; ++f) {
            const auto face = rBlock.Face(Cell, f);
            for (std::size_t k = 0; k < face.second; ++k)
                rOut.mFaceNodes.push_back(
                    poly_intern(rOut, rCoords, face.first[k], rPoints, PointDim));
            rOut.mFaceStart.push_back(static_cast<std::uint32_t>(rOut.mFaceNodes.size()));
        }
        return true;
    }
    if (rBlock.IsRagged())
        return false;  // 1-level (polygon): a surface cell, no enclosed volume

    const CellType ct = cell_type_from_name(rBlock.Type());
    const std::vector<CellFaceDef>& faces = cell_faces(ct);
    if (faces.empty())
        return false;
    const int corners = cell_corner_count(ct);
    if (corners <= 0)
        return false;

    const NDArray& conn = rBlock.Conn();
    const std::size_t npc = rBlock.NodesPerCell();
    // Intern the corners in CONNECTIVITY order first, so a tabulated cell's
    // local ring indices are exactly `cell_faces`' own `mNodes` values. That is
    // a contract, not an accident: `gradient` carries per-corner field values in
    // connectivity order and indexes them with the ring, and
    // `cell_volume_from_corners` hands in corner coordinates the same way.
    // Interning lazily per face would instead number nodes in face-traversal
    // order and silently permute both.
    rOut.mNodes.reserve(static_cast<std::size_t>(corners));
    rCoords.reserve(static_cast<std::size_t>(corners));
    for (int i = 0; i < corners; ++i) {
        const std::int64_t id = read_int(conn, Cell * npc + static_cast<std::size_t>(i));
        rOut.mNodes.push_back(id);
        rCoords.push_back(read_point(rPoints, PointDim, id));
    }
    rOut.mFaceStart.push_back(0);
    for (const CellFaceDef& f : faces) {
        // Corners only, matching every other measure kernel in the repo.
        for (int i = 0; i < f.mNumCorners; ++i)
            rOut.mFaceNodes.push_back(f.mNodes[i]);
        rOut.mFaceStart.push_back(static_cast<std::uint32_t>(rOut.mFaceNodes.size()));
    }
    return true;
}

Vec3 polygon_area_vector(const Vec3* pCoords, const std::uint32_t* pRing, std::size_t N) {
    // Newell: sum of cross products around the ring. Translation-invariance and
    // start-node independence both follow from the sum telescoping, which is
    // what makes it the right primitive for a face whose node order the file
    // format does not constrain.
    Vec3 n{0.0, 0.0, 0.0};
    if (N < 3)
        return n;
    for (std::size_t i = 0; i < N; ++i) {
        const Vec3& a = pCoords[pRing[i]];
        const Vec3& b = pCoords[pRing[(i + 1) % N]];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    return vec3_scale(n, 0.5);
}

double polygon_area(const Vec3* pCoords, const std::uint32_t* pRing, std::size_t N) {
    return vec3_norm(polygon_area_vector(pCoords, pRing, N));
}

Vec3 polygon_area_vector(const Vec3* pCoords, std::size_t N) {
    // Same arithmetic as the ring form, with the indirection removed rather
    // than an identity ring materialized per call.
    Vec3 n{0.0, 0.0, 0.0};
    if (N < 3)
        return n;
    for (std::size_t i = 0; i < N; ++i) {
        const Vec3& a = pCoords[i];
        const Vec3& b = pCoords[(i + 1) % N];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    return vec3_scale(n, 0.5);
}

double polygon_area(const Vec3* pCoords, std::size_t N) {
    return vec3_norm(polygon_area_vector(pCoords, N));
}

double cell_volume_from_corners(const Vec3* pCoords, CellType Type) {
    const std::vector<CellFaceDef>& faces = cell_faces(Type);
    if (faces.empty())
        return std::nan("");
    const int corners = cell_corner_count(Type);
    if (corners <= 0)
        return std::nan("");
    // The coordinates are already the cell's corners in order, so the rings are
    // just cell_faces' local indices with an identity node mapping.
    CellRings rings;
    rings.mNodes.resize(static_cast<std::size_t>(corners));
    for (int i = 0; i < corners; ++i)
        rings.mNodes[static_cast<std::size_t>(i)] = i;
    rings.mFaceStart.push_back(0);
    for (const CellFaceDef& f : faces) {
        for (int i = 0; i < f.mNumCorners; ++i)
            rings.mFaceNodes.push_back(f.mNodes[i]);
        rings.mFaceStart.push_back(static_cast<std::uint32_t>(rings.mFaceNodes.size()));
    }
    return poly_measure(rings, pCoords).mVolume;
}

Vec3 polygon_centroid(const Vec3* pCoords, const std::uint32_t* pRing, std::size_t N) {
    Vec3 c{0.0, 0.0, 0.0};
    if (N == 0)
        return c;
    for (std::size_t i = 0; i < N; ++i)
        c = vec3_add(c, pCoords[pRing[i]]);
    return vec3_scale(c, 1.0 / static_cast<double>(N));
}

RingOrientation orient_rings(CellRings& rRings, const Vec3* pCoords) {
    const std::size_t nf = rRings.NumFaces();
    if (nf == 0)
        return RingOrientation::Unorientable;

    // Map each undirected edge to the (face, directed-as-stored) uses of it. A
    // closed surface uses every undirected edge exactly twice.
    struct EdgeUse {
        std::size_t mFace;
        bool mForward;  // traversed low->high as stored
    };
    std::unordered_map<PolyEdge, std::vector<EdgeUse>, PolyEdgeHash> uses;
    uses.reserve(nf * 4);
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t n = rRings.FaceSize(f);
        if (n < 3)
            return RingOrientation::Unorientable;
        const std::uint32_t* ring = rRings.Face(f);
        for (std::size_t i = 0; i < n; ++i) {
            const std::uint32_t a = ring[i], b = ring[(i + 1) % n];
            if (a == b)
                return RingOrientation::Unorientable;
            uses[poly_undirected(a, b)].push_back({f, a < b});
        }
    }
    for (const auto& kv : uses)
        if (kv.second.size() != 2)
            return RingOrientation::Unorientable;

    // BFS the face dual. Two faces sharing an undirected edge are consistently
    // wound iff they traverse it in OPPOSITE directions, so an equal
    // `mForward` means the neighbour must be reversed.
    std::vector<int> flip(nf, -1);  // -1 unvisited, 0 keep, 1 reverse
    std::vector<std::size_t> stack;
    flip[0] = 0;
    stack.push_back(0);
    std::size_t visited = 0;
    while (!stack.empty()) {
        const std::size_t f = stack.back();
        stack.pop_back();
        ++visited;
        const std::size_t n = rRings.FaceSize(f);
        const std::uint32_t* ring = rRings.Face(f);
        for (std::size_t i = 0; i < n; ++i) {
            const std::uint32_t a = ring[i], b = ring[(i + 1) % n];
            const auto& u = uses[poly_undirected(a, b)];
            const EdgeUse& mine = u[0].mFace == f && u[0].mForward == (a < b) ? u[0] : u[1];
            const EdgeUse& other = (&mine == &u[0]) ? u[1] : u[0];
            if (other.mFace == f)
                return RingOrientation::Unorientable;  // a face using its own edge twice
            // Effective direction of this edge in each face after its own flip.
            const bool mine_fwd = mine.mForward != (flip[f] == 1);
            const int want = (other.mForward == mine_fwd) ? 1 : 0;
            if (flip[other.mFace] == -1) {
                flip[other.mFace] = want;
                stack.push_back(other.mFace);
            } else if (flip[other.mFace] != want) {
                return RingOrientation::Unorientable;  // non-orientable (Mobius-like)
            }
        }
    }
    if (visited != nf)
        return RingOrientation::Unorientable;  // disconnected face set

    bool changed = false;
    for (std::size_t f = 0; f < nf; ++f) {
        if (flip[f] == 1) {
            std::uint32_t* ring = rRings.FaceMutable(f);
            std::reverse(ring, ring + rRings.FaceSize(f));
            changed = true;
        }
    }

    // Consistent now, but possibly consistently INWARD. One signed-volume probe
    // decides, and a global reverse fixes it.
    if (poly_measure(rRings, pCoords).mVolume < 0.0) {
        for (std::size_t f = 0; f < nf; ++f) {
            std::uint32_t* ring = rRings.FaceMutable(f);
            std::reverse(ring, ring + rRings.FaceSize(f));
        }
        changed = true;
    }
    return changed ? RingOrientation::Repaired : RingOrientation::Consistent;
}

PolyMeasure poly_measure(const CellRings& rRings, const Vec3* pCoords) {
    PolyMeasure out;
    const std::size_t nf = rRings.NumFaces();
    if (nf == 0 || rRings.mNodes.empty())
        return out;

    // Recentre on the cell's corner average before any arithmetic. This is
    // load-bearing, not an optimization: V = sum(x . A) / 3 only telescopes
    // because sum(A) == 0 over a closed surface, so at |x| ~ 1e8 the raw form
    // loses eight digits to cancellation and then divides by the result.
    Vec3 origin{0.0, 0.0, 0.0};
    const std::size_t nn = rRings.mNodes.size();
    for (std::size_t i = 0; i < nn; ++i)
        origin = vec3_add(origin, pCoords[i]);
    origin = vec3_scale(origin, 1.0 / static_cast<double>(nn));

    std::vector<Vec3> local(nn);
    for (std::size_t i = 0; i < nn; ++i)
        local[i] = vec3_sub(pCoords[i], origin);

    double vol6 = 0.0;           // 6 * volume
    Vec3 moment{0.0, 0.0, 0.0};  // sum over tets of (vol6_i * 4*centroid_i)
    double area = 0.0;
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t n = rRings.FaceSize(f);
        const std::uint32_t* ring = rRings.Face(f);
        if (n < 3)
            continue;
        // The corner-average fan: apex at the face's own corner average, one
        // triangle per edge of the ring. See the header for why this apex and
        // not the ring's first node.
        const Vec3 fc = polygon_centroid(local.data(), ring, n);
        for (std::size_t i = 0; i < n; ++i) {
            const Vec3& a = local[ring[i]];
            const Vec3& b = local[ring[(i + 1) % n]];
            const Vec3 cr = vec3_cross(vec3_sub(a, fc), vec3_sub(b, fc));
            area += 0.5 * vec3_norm(cr);
            // Signed volume of the tet (origin, fc, a, b), x6.
            const double t = triple_product(fc, a, b);
            vol6 += t;
            // Its centroid is (0 + fc + a + b)/4; accumulate 4*centroid*vol6.
            const Vec3 s = vec3_add(fc, vec3_add(a, b));
            moment = vec3_add(moment, vec3_scale(s, t));
        }
    }

    out.mVolume = vol6 / 6.0;
    out.mSurfaceArea = area;
    out.mClosed = true;
    out.mOrientable = true;
    if (std::abs(vol6) > 0.0) {
        // centroid = sum(v_i * c_i) / sum(v_i), with c_i = s_i/4 and
        // v_i = t_i/6 -- the constant factors cancel except the 1/4.
        out.mCentroid = vec3_add(origin, vec3_scale(moment, 0.25 / vol6));
    } else {
        out.mCentroid = origin;
    }
    return out;
}

FacetKey::FacetKey(const std::int64_t* pIds, std::size_t N) : mN(N) {
    if (N <= kInline) {
        for (std::size_t i = 0; i < N; ++i)
            mInline[i] = pIds[i];
        std::sort(mInline.begin(), mInline.begin() + static_cast<std::ptrdiff_t>(N));
    } else {
        mHeap.assign(pIds, pIds + N);
        std::sort(mHeap.begin(), mHeap.end());
    }
}

bool FacetKey::operator==(const FacetKey& rOther) const {
    if (mN != rOther.mN)
        return false;
    return std::equal(Data(), Data() + mN, rOther.Data());
}

std::size_t FacetKeyHash::operator()(const FacetKey& rKey) const {
    std::size_t h = rKey.Size() * 0x9E3779B97F4A7C15ULL;
    const std::int64_t* p = rKey.Data();
    for (std::size_t i = 0; i < rKey.Size(); ++i)
        h = h * 0x100000001B3ULL + static_cast<std::size_t>(p[i]);
    return h;
}

}  // namespace detail
}  // namespace meshioplusplus
