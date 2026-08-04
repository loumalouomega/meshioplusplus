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
// Field differential operators (gradient / divergence / curl) of a point_data
// field, by Green-Gauss over the cell's own faces or by a least-squares fit over
// its node-sharing neighbours. Built entirely through the uniform mesh API so it
// compiles under every mesh backend. See operations/gradient.hpp for the
// contract, the exactness argument and the naming/shape rules.
//
// FP-order discipline: every accumulation inside a cell runs in a fixed order
// (faces in cell_faces() order, fan triangles in ring order, neighbours in
// ascending global cell index) because addition is not associative, and the
// numpy twin in src/python/meshioplusplus/_gradient.py replicates that order
// token for token.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/gradient.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_adjacency.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/data_manage.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

constexpr const char* kGradPrefix = "meshio++: gradient: ";

// --- per-cell geometry -------------------------------------------------------

/// One cell, reduced to what both methods need: its corner coordinates and the
/// matching field values, both recentred on the corner average.
///
/// Recentring is load-bearing, not an optimization. `V = (1/3) sum x_j . A_j`
/// only telescopes to the cell's volume because `sum A_j == 0` over the closed
/// fan surface; in floating point, on a mesh at x ~ 1e6 with cells of size 1,
/// every term is ~1e6 * area and they cancel to ~area, losing about six decimal
/// digits of a quantity we then divide by. The same applies to the numerator
/// when the field carries a large constant offset (Kelvin, Pascals). Subtracting
/// the corner means changes nothing mathematically and removes the cancellation.
struct GradCell {
    std::vector<Vec3> mCorners;   ///< Corner coords, relative to `mOrigin`.
    std::vector<double> mValues;  ///< Corner field values, relative to `mValue0`, `ncomp` each.
    Vec3 mOrigin{0.0, 0.0, 0.0};  ///< Arithmetic mean of the raw corner coords.
    std::vector<double> mValue0;  ///< Arithmetic mean of the raw corner values, `ncomp` long.
    bool mSupported = false;      ///< Whether the cell can be differentiated at all.
    CellType mType = CellType::Custom;
    int mDim = -1;
    /// The cell's boundary as face rings, indexing `mCorners`/`mValues`. Empty
    /// for a 2D/1D cell. This is what makes Green-Gauss polyhedral: a
    /// `polyhedron12` and a `hexahedron` differ only in what fills this.
    detail::CellRings mRings;
};

/// Reads the leading `NumCorners` connectivity entries of one cell.
void grad_corner_nodes(const Mesh::CellView& rBlock, std::size_t Cell, std::size_t NumCorners,
                       std::vector<std::int64_t>& rOut) {
    rOut.clear();
    const NDArray& conn = rBlock.Conn();
    const std::size_t npc = rBlock.NodesPerCell();
    for (std::size_t k = 0; k < NumCorners; ++k)
        rOut.push_back(detail::read_int(conn, Cell * npc + k));
}

/// Gathers a cell's corner coordinates and field values, recentred.
///
/// `rOut.mSupported` is false when the block is ragged/polyhedron, has no corner
/// count, or references an out-of-range node — the caller then emits NaN and
/// counts the cell rather than guessing.
void grad_load_cell(const Mesh::CellView& rBlock, std::size_t Cell, const NDArray& rPoints,
                    std::size_t PointDim, const NDArray& rField, std::size_t NumComp,
                    std::size_t NumPoints, std::size_t NumCorners, GradCell& rOut) {
    rOut.mSupported = false;
    rOut.mCorners.clear();
    rOut.mValues.clear();
    rOut.mValue0.assign(NumComp, 0.0);
    rOut.mRings.Clear();

    static thread_local std::vector<std::int64_t> nodes;
    static thread_local std::vector<Vec3> ring_coords;
    // The rings are filled for every 3D cell, tabulated or polyhedral, and are
    // what Green-Gauss integrates over. `cell_rings` interns a tabulated cell's
    // nodes in CONNECTIVITY order, so the ring indices line up with
    // mCorners/mValues below either way -- that is the contract this relies on.
    const bool has_rings =
        detail::cell_rings(rBlock, Cell, rPoints, PointDim, rOut.mRings, ring_coords);
    if (rBlock.IsPolyhedron()) {
        if (!has_rings)
            return;
        nodes = rOut.mRings.mNodes;
    } else {
        // 2D/1D cells legitimately have no rings; they take the corner-ring
        // path further down and never reach grad_green_gauss_3d.
        if (rBlock.IsRagged() || NumCorners == 0)
            return;
        grad_corner_nodes(rBlock, Cell, NumCorners, nodes);
    }
    for (std::int64_t nid : nodes)
        if (nid < 0 || static_cast<std::size_t>(nid) >= NumPoints)
            return;
    const std::size_t ncorners = nodes.size();

    rOut.mCorners.reserve(ncorners);
    rOut.mValues.reserve(ncorners * NumComp);
    for (std::int64_t nid : nodes) {
        rOut.mCorners.push_back(detail::read_point(rPoints, PointDim, nid));
        for (std::size_t k = 0; k < NumComp; ++k)
            rOut.mValues.push_back(
                detail::read_double(rField, static_cast<std::size_t>(nid) * NumComp + k));
    }

    // Ascending corner order, left to right — the numpy twin sums identically.
    Vec3 origin{0.0, 0.0, 0.0};
    for (const Vec3& r_p : rOut.mCorners)
        origin = detail::vec3_add(origin, r_p);
    const double inv = 1.0 / static_cast<double>(ncorners);
    origin = detail::vec3_scale(origin, inv);
    for (std::size_t k = 0; k < NumComp; ++k) {
        double s = 0.0;
        for (std::size_t i = 0; i < ncorners; ++i)
            s += rOut.mValues[i * NumComp + k];
        rOut.mValue0[k] = s * inv;
    }
    for (Vec3& r_p : rOut.mCorners)
        r_p = detail::vec3_sub(r_p, origin);
    for (std::size_t i = 0; i < ncorners; ++i)
        for (std::size_t k = 0; k < NumComp; ++k)
            rOut.mValues[i * NumComp + k] -= rOut.mValue0[k];
    rOut.mOrigin = origin;
    rOut.mSupported = true;
}

// --- Green-Gauss -------------------------------------------------------------

/// Green-Gauss over a 3D cell's faces. Writes `NumComp * 3` values to `pOut`
/// (`[component][derivative]`), or returns false when the cell's fan volume is
/// degenerate relative to its own size.
bool grad_green_gauss_3d(const GradCell& rCell, std::size_t NumComp, double* pOut) {
    // Face rings rather than `cell_faces(mType)`: that is the ONLY difference
    // between differentiating a hexahedron and a `polyhedron12`, which is what
    // makes Green-Gauss naturally polyhedral. The arithmetic below is unchanged.
    const std::size_t nfaces = rCell.mRings.NumFaces();
    if (nfaces == 0)
        return false;

    std::vector<double> num(NumComp * 3, 0.0);
    double volume = 0.0;
    double area_scale = 0.0;

    for (std::size_t fi = 0; fi < nfaces; ++fi) {
        const std::uint32_t* p_ring = rCell.mRings.Face(fi);
        const std::size_t m = rCell.mRings.FaceSize(fi);
        // Face centre and face-centre value: the fan apex. For linear f the
        // corner average is the ONLY apex whose value is known exactly without
        // shape functions -- see the header. Do not switch to the area centroid.
        Vec3 c{0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < m; ++i)
            c = detail::vec3_add(c, rCell.mCorners[p_ring[i]]);
        const double inv_m = 1.0 / static_cast<double>(m);
        c = detail::vec3_scale(c, inv_m);

        static thread_local std::vector<double> fc;
        fc.assign(NumComp, 0.0);
        for (std::size_t k = 0; k < NumComp; ++k) {
            double s = 0.0;
            for (std::size_t i = 0; i < m; ++i)
                s += rCell.mValues[static_cast<std::size_t>(p_ring[i]) * NumComp + k];
            fc[k] = s * inv_m;
        }

        for (std::size_t i = 0; i < m; ++i) {
            const std::size_t a = p_ring[i];
            const std::size_t b = p_ring[(i + 1) % m];
            const Vec3& r_pa = rCell.mCorners[a];
            const Vec3& r_pb = rCell.mCorners[b];
            const Vec3 av = detail::vec3_scale(
                detail::vec3_cross(detail::vec3_sub(r_pa, c), detail::vec3_sub(r_pb, c)), 0.5);
            const Vec3 xj =
                detail::vec3_scale(detail::vec3_add(detail::vec3_add(c, r_pa), r_pb), 1.0 / 3.0);
            volume += detail::vec3_dot(xj, av) * (1.0 / 3.0);
            area_scale += detail::vec3_norm(av);
            for (std::size_t k = 0; k < NumComp; ++k) {
                const double fj =
                    (fc[k] + rCell.mValues[a * NumComp + k] + rCell.mValues[b * NumComp + k]) *
                    (1.0 / 3.0);
                num[k * 3 + 0] += fj * av[0];
                num[k * 3 + 1] += fj * av[1];
                num[k * 3 + 2] += fj * av[2];
            }
        }
    }

    // Relative degeneracy test: a volume must be large next to the cell's own
    // (area)^(3/2) scale, never against an absolute epsilon.
    const double scale = area_scale * std::sqrt(area_scale);
    if (!(std::abs(volume) > 1e-12 * scale))
        return false;
    const double inv_v = 1.0 / volume;
    for (std::size_t i = 0; i < NumComp * 3; ++i)
        pOut[i] = num[i] * inv_v;
    return true;
}

/// The recentred Newell area vector of a 2D cell's corner ring.
///
/// Deliberately the recentred form rather than `dataops_polygon_area`'s
/// `sum p_i x p_{i+1}`: algebraically identical (the cross terms telescope) but
/// far better conditioned. Do not "unify" the two.
Vec3 grad_ring_area_vector(const std::vector<Vec3>& rRing) {
    const std::size_t n = rRing.size();
    Vec3 av{0.0, 0.0, 0.0};
    for (std::size_t i = 1; i + 1 < n; ++i)
        av = detail::vec3_add(av, detail::vec3_cross(detail::vec3_sub(rRing[i], rRing[0]),
                                                     detail::vec3_sub(rRing[i + 1], rRing[0])));
    return detail::vec3_scale(av, 0.5);
}

/// Green-Gauss over a 2D cell's corner ring, using the in-plane outward normal
/// `t x N`. Invariant under reversal and cyclic rotation of the ring.
bool grad_green_gauss_2d(const GradCell& rCell, std::size_t NumComp, double* pOut) {
    const std::size_t n = rCell.mCorners.size();
    if (n < 3)
        return false;
    const Vec3 av = grad_ring_area_vector(rCell.mCorners);
    const double area = detail::vec3_norm(av);
    double edge_scale = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        edge_scale +=
            detail::vec3_norm(detail::vec3_sub(rCell.mCorners[(i + 1) % n], rCell.mCorners[i]));
    if (!(area > 1e-12 * (edge_scale * edge_scale)))
        return false;
    const Vec3 nrm = detail::vec3_scale(av, 1.0 / area);

    std::vector<double> num(NumComp * 3, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t a = i;
        const std::size_t b = (i + 1) % n;
        const Vec3 t = detail::vec3_sub(rCell.mCorners[b], rCell.mCorners[a]);
        const Vec3 mds = detail::vec3_cross(t, nrm);
        for (std::size_t k = 0; k < NumComp; ++k) {
            const double fbar =
                (rCell.mValues[a * NumComp + k] + rCell.mValues[b * NumComp + k]) * 0.5;
            num[k * 3 + 0] += fbar * mds[0];
            num[k * 3 + 1] += fbar * mds[1];
            num[k * 3 + 2] += fbar * mds[2];
        }
    }
    const double inv_a = 1.0 / area;
    for (std::size_t i = 0; i < NumComp * 3; ++i)
        pOut[i] = num[i] * inv_a;
    return true;
}

/// Green-Gauss on a 1D cell: the two endpoints are its "faces", with outward
/// normals -/+ t_hat and "volume" |p1 - p0|.
bool grad_green_gauss_1d(const GradCell& rCell, std::size_t NumComp, double* pOut) {
    if (rCell.mCorners.size() < 2)
        return false;
    const Vec3 t = detail::vec3_sub(rCell.mCorners[1], rCell.mCorners[0]);
    const double len = detail::vec3_norm(t);
    if (!(len > 0.0))
        return false;
    const Vec3 dir = detail::vec3_scale(t, 1.0 / len);
    for (std::size_t k = 0; k < NumComp; ++k) {
        const double slope =
            (rCell.mValues[1 * NumComp + k] - rCell.mValues[0 * NumComp + k]) / len;
        pOut[k * 3 + 0] = slope * dir[0];
        pOut[k * 3 + 1] = slope * dir[1];
        pOut[k * 3 + 2] = slope * dir[2];
    }
    return true;
}

bool grad_green_gauss(const GradCell& rCell, int Dim, std::size_t NumComp, double* pOut) {
    if (Dim == 3)
        return grad_green_gauss_3d(rCell, NumComp, pOut);
    if (Dim == 2)
        return grad_green_gauss_2d(rCell, NumComp, pOut);
    if (Dim == 1)
        return grad_green_gauss_1d(rCell, NumComp, pOut);
    return false;
}

// --- least squares -----------------------------------------------------------

/// Solves the symmetric 3x3 system `M g = b` in closed cofactor form.
///
/// The conditioning bound is transcribed from `decimate.cpp`'s optimal-placement
/// solve, including the `scale * scale * scale` multiplication form: `std::pow`
/// (and `np.power` in the twin) rounds differently and would break parity.
bool grad_solve3(const double* pM, const double* pB, double* pOut) {
    const double c00 = pM[4] * pM[8] - pM[5] * pM[5];
    const double c01 = pM[2] * pM[5] - pM[1] * pM[8];
    const double c02 = pM[1] * pM[5] - pM[2] * pM[4];
    const double det = pM[0] * c00 + pM[1] * c01 + pM[2] * c02;
    double scale = std::abs(pM[0]);
    scale = std::max(scale, std::abs(pM[1]));
    scale = std::max(scale, std::abs(pM[2]));
    scale = std::max(scale, std::abs(pM[4]));
    scale = std::max(scale, std::abs(pM[5]));
    scale = std::max(scale, std::abs(pM[8]));
    if (!(std::abs(det) > 1e-12 * (scale * scale * scale)))
        return false;
    const double c11 = pM[0] * pM[8] - pM[2] * pM[2];
    const double c12 = pM[1] * pM[2] - pM[0] * pM[5];
    const double c22 = pM[0] * pM[4] - pM[1] * pM[1];
    const double inv = 1.0 / det;
    pOut[0] = (c00 * pB[0] + c01 * pB[1] + c02 * pB[2]) * inv;
    pOut[1] = (c01 * pB[0] + c11 * pB[1] + c12 * pB[2]) * inv;
    pOut[2] = (c02 * pB[0] + c12 * pB[1] + c22 * pB[2]) * inv;
    return true;
}

/// One cell's least-squares stencil: the neighbour offsets and value deltas,
/// already recentred on this cell.
struct GradStencil {
    std::vector<Vec3> mOffsets;
    std::vector<double> mDeltas;  ///< `NumComp` per offset.
};

/// Least-squares fit over a prepared stencil. `pNormal` is null in 3D and the
/// cell's unit plane normal in 2D.
bool grad_least_squares(const GradStencil& rStencil, std::size_t NumComp, int Dim,
                        const Vec3* pNormal, double* pOut) {
    double m[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> b(NumComp * 3, 0.0);
    const std::size_t n = rStencil.mOffsets.size();
    for (std::size_t j = 0; j < n; ++j) {
        Vec3 d = rStencil.mOffsets[j];
        if (pNormal != nullptr)
            d = detail::vec3_sub(d, detail::vec3_scale(*pNormal, detail::vec3_dot(d, *pNormal)));
        const double d2 = detail::vec3_norm_sq(d);
        // A zero offset would give an infinite weight and poison the whole
        // matrix with NaN; two node-sharing cells really can have coincident
        // corner means (a collapsed cell, a line lying inside a triangle).
        if (!(d2 > 0.0))
            continue;
        const double w = 1.0 / d2;
        for (std::size_t r = 0; r < 3; ++r)
            for (std::size_t c = 0; c < 3; ++c)
                m[r * 3 + c] += w * d[r] * d[c];
        for (std::size_t k = 0; k < NumComp; ++k) {
            const double wf = w * rStencil.mDeltas[j * NumComp + k];
            b[k * 3 + 0] += wf * d[0];
            b[k * 3 + 1] += wf * d[1];
            b[k * 3 + 2] += wf * d[2];
        }
    }

    if (Dim == 3 || Dim == 2) {
        if (pNormal != nullptr) {
            // Rank-1 regularization instead of a tangent-plane basis: adding
            // s * N N^T makes the in-plane rank-2 matrix invertible without
            // touching b (every offset was projected, so b . N == 0), so the
            // solution automatically satisfies g . N == 0 and the SAME 3x3
            // solver serves both dimensions. A Gram-Schmidt basis would need an
            // arbitrary tie-break transcribed bit-for-bit into the numpy twin.
            double s = (m[0] + m[4] + m[8]) * 0.5;
            if (!(s > 0.0))
                s = 1.0;
            const Vec3& r_nrm = *pNormal;
            for (std::size_t r = 0; r < 3; ++r)
                for (std::size_t c = 0; c < 3; ++c)
                    m[r * 3 + c] += s * r_nrm[r] * r_nrm[c];
        }
        for (std::size_t k = 0; k < NumComp; ++k)
            if (!grad_solve3(m, b.data() + k * 3, pOut + k * 3))
                return false;
        return true;
    }

    // 1D: a single scalar equation per component. The generic bound would
    // degenerate to |M00| <= 1e-12 * |M00|, which is only true at exactly zero,
    // so say that outright rather than shipping an expression that reads as a bug.
    if (m[0] == 0.0 && m[4] == 0.0 && m[8] == 0.0)
        return false;
    const double trace = m[0] + m[4] + m[8];
    if (trace == 0.0)
        return false;
    for (std::size_t k = 0; k < NumComp; ++k) {
        // The offsets are collinear, so M = w |d|^2 d_hat d_hat^T has rank 1 and
        // the fit reduces to b projected onto that direction over the trace.
        pOut[k * 3 + 0] = b[k * 3 + 0] / trace;
        pOut[k * 3 + 1] = b[k * 3 + 1] / trace;
        pOut[k * 3 + 2] = b[k * 3 + 2] / trace;
    }
    return true;
}

// --- operator application ----------------------------------------------------

/// Number of output components for one operator.
std::size_t grad_out_components(GradientOperator Op, std::size_t NumComp) {
    if (Op == GradientOperator::Divergence)
        return 1;
    if (Op == GradientOperator::Curl)
        return 3;
    return NumComp * 3;
}

/// Reduces a full `[component][derivative]` block to the requested operator.
void grad_apply_operator(GradientOperator Op, const double* pGrad, std::size_t NumComp,
                         double* pOut) {
    if (Op == GradientOperator::Gradient) {
        for (std::size_t i = 0; i < NumComp * 3; ++i)
            pOut[i] = pGrad[i];
        return;
    }
    // Divergence and curl read a 3x3 block; a 2-component field is (u, v, 0),
    // the same padding convention read_point applies to 2D points.
    double g[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (std::size_t k = 0; k < NumComp && k < 3; ++k)
        for (std::size_t j = 0; j < 3; ++j)
            g[k * 3 + j] = pGrad[k * 3 + j];
    if (Op == GradientOperator::Divergence) {
        pOut[0] = g[0 * 3 + 0] + g[1 * 3 + 1] + g[2 * 3 + 2];
        return;
    }
    pOut[0] = g[2 * 3 + 1] - g[1 * 3 + 2];
    pOut[1] = g[0 * 3 + 2] - g[2 * 3 + 0];
    pOut[2] = g[1 * 3 + 0] - g[0 * 3 + 1];
}

/// The default output name for an operator.
std::string grad_default_name(const std::string& rArray, GradientOperator Op) {
    if (Op == GradientOperator::Divergence)
        return rArray + kDivergenceSuffix;
    if (Op == GradientOperator::Curl)
        return rArray + kCurlSuffix;
    return rArray + kGradientSuffix;
}

/// `{rows}` for a scalar, `{rows, ncomp}` otherwise — data_average's rule, so a
/// scalar result stays 1-D rather than becoming an (n, 1) column.
std::vector<std::size_t> grad_shape(std::size_t Rows, std::size_t NumComp) {
    if (NumComp <= 1)
        return {Rows};
    return {Rows, NumComp};
}

/// Max topological dimension over the mesh's blocks; ragged polygon blocks count
/// as 2 and polyhedron blocks as 3, matching partition's dual-dimension rule.
int grad_mesh_dimension(const Mesh& rMesh) {
    int dim = -1;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.NumCells() == 0)
            continue;
        int d;
        if (cb.IsPolyhedron())
            d = 3;
        else if (cb.IsRagged())
            d = 2;
        else
            d = cell_type_dimension(cell_type_from_name(cb.Type()));
        dim = std::max(dim, d);
    }
    return dim;
}

}  // namespace

GradientOperator gradient_operator_from_name(const std::string& rName) {
    if (rName == "gradient" || rName == "grad" || rName.empty())
        return GradientOperator::Gradient;
    if (rName == "divergence" || rName == "div")
        return GradientOperator::Divergence;
    if (rName == "curl" || rName == "rot")
        return GradientOperator::Curl;
    throw std::invalid_argument(std::string(kGradPrefix) + "unknown operator '" + rName +
                                "' (expected 'gradient', 'divergence', or 'curl')");
}

GradientMethod gradient_method_from_name(const std::string& rName) {
    if (rName == "green-gauss" || rName == "green_gauss" || rName == "gg" || rName.empty())
        return GradientMethod::GreenGauss;
    if (rName == "least-squares" || rName == "least_squares" || rName == "lsq")
        return GradientMethod::LeastSquares;
    throw std::invalid_argument(std::string(kGradPrefix) + "unknown method '" + rName +
                                "' (expected 'green-gauss' or 'least-squares')");
}

GradientResult gradient(const Mesh& rMesh, const GradientOptions& rOptions) {
    // --- validation, all before any work ------------------------------------
    if (rOptions.mArrayName.empty())
        throw std::invalid_argument(std::string(kGradPrefix) + "an array name is required");
    if (rOptions.mLocation == DataLocation::Field)
        throw std::invalid_argument(std::string(kGradPrefix) +
                                    "field_data has no location on the mesh; the result must be "
                                    "'point' or 'cell'");
    if (!rMesh.HasPointData(rOptions.mArrayName)) {
        // A cell_data field is piecewise constant and has no derivative. Name
        // the fix, exactly as isosurface does for the same reason.
        if (rMesh.HasCellData(rOptions.mArrayName))
            throw std::invalid_argument(
                std::string(kGradPrefix) + "'" + rOptions.mArrayName +
                "' is cell_data, which is piecewise constant and has no derivative; convert it "
                "first with cell_data_to_point_data (CLI: meshioplusplus data to-point)");
        throw std::invalid_argument(
            std::string(kGradPrefix) +
            data_unknown_key_message(rMesh, DataLocation::Point, rOptions.mArrayName));
    }

    const NDArray& field = rMesh.PointData(rOptions.mArrayName);
    const std::size_t field_comp = data_num_components(field);
    if (field_comp == 0 || detail::rows(field) == 0)
        throw std::invalid_argument(std::string(kGradPrefix) + "'" + rOptions.mArrayName +
                                    "' carries no values");
    if (detail::rows(field) != rMesh.NumPoints())
        throw std::invalid_argument(std::string(kGradPrefix) + "'" + rOptions.mArrayName +
                                    "' has " + std::to_string(detail::rows(field)) +
                                    " rows on a mesh of " + std::to_string(rMesh.NumPoints()) +
                                    " points");
    if (rOptions.mComponent.has_value()) {
        if (rOptions.mOperator != GradientOperator::Gradient)
            throw std::invalid_argument(std::string(kGradPrefix) +
                                        "a component selection applies to the gradient only; "
                                        "divergence and curl consume the whole vector");
        const int c = *rOptions.mComponent;
        if (c < 0 || static_cast<std::size_t>(c) >= field_comp)
            throw std::invalid_argument(std::string(kGradPrefix) + "component " +
                                        std::to_string(c) + " is out of range for '" +
                                        rOptions.mArrayName + "', which has " +
                                        std::to_string(field_comp) + " component(s)");
    }
    if (rOptions.mOperator != GradientOperator::Gradient && field_comp != 2 && field_comp != 3)
        throw std::invalid_argument(
            std::string(kGradPrefix) + "'" + rOptions.mArrayName + "' has " +
            std::to_string(field_comp) +
            " components; divergence and curl need a vector field of 2 or 3");

    const std::string out_name = rOptions.mOutputName.empty()
                                     ? grad_default_name(rOptions.mArrayName, rOptions.mOperator)
                                     : rOptions.mOutputName;
    if (!rOptions.mOverwrite) {
        const bool taken = rOptions.mLocation == DataLocation::Point ? rMesh.HasPointData(out_name)
                                                                     : rMesh.HasCellData(out_name);
        if (taken)
            throw std::invalid_argument(
                std::string(kGradPrefix) + "'" + out_name + "' already exists in " +
                data_location_name(rOptions.mLocation) + " (pass overwrite=true to replace it)");
    }

    // The field components actually differentiated: one when a component was
    // selected, all of them otherwise.
    const std::size_t work_comp = rOptions.mComponent.has_value() ? 1 : field_comp;
    const std::size_t comp_base =
        rOptions.mComponent.has_value() ? static_cast<std::size_t>(*rOptions.mComponent) : 0;
    const std::size_t out_comp = grad_out_components(rOptions.mOperator, work_comp);

    // Geometry, connectivity, regions, property sets and every existing array
    // ride through untouched. No warn_regions_dropped here -- unlike every
    // restructuring operation, this one does not renumber anything.
    Mesh out = detail::clone_mesh(rMesh);
    const std::size_t nblocks = rMesh.NumCellBlocks();
    if (nblocks == 0)
        return GradientResult{std::move(out), 0, 0};

    const int dim = grad_mesh_dimension(rMesh);
    if (dim <= 0)
        throw std::invalid_argument(
            std::string(kGradPrefix) +
            "the mesh has no cells of topological dimension 1 or more to differentiate over");

    // A compacted view of the field holding only the components being worked on,
    // so every downstream loop reads a dense `work_comp`-wide row.
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    const std::size_t npoints = rMesh.NumPoints();

    NDArray work(DType::Float64, grad_shape(npoints, work_comp));
    {
        double* pw = work.As<double>();
        parallel_for_bw(npoints, [&](std::size_t i) {
            for (std::size_t k = 0; k < work_comp; ++k)
                pw[i * work_comp + k] = detail::read_double(field, i * field_comp + comp_base + k);
        });
    }

    // Per-block corner counts and eligibility, resolved once.
    std::vector<std::size_t> corners(nblocks, 0);
    std::vector<std::uint8_t> eligible(nblocks, 0);
    std::vector<CellType> types(nblocks, CellType::Custom);
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.IsPolyhedron()) {
            // Polyhedra are 3D by definition and carry their own faces, so the
            // corner-count and face-table checks below have nothing to say
            // about them. `corners` stays 0: grad_load_cell reads the count
            // from the cell's own rings, which vary per cell.
            if (dim != 3)
                continue;
            types[b] = CellType::Polyhedron;
            eligible[b] = 1;
            continue;
        }
        if (cb.IsRagged())
            continue;
        const CellType t = cell_type_from_name(cb.Type());
        types[b] = t;
        if (cell_type_dimension(t) != dim)
            continue;
        const int nc = detail::cell_corner_count(t);
        if (nc <= 0)
            continue;
        // 3D additionally needs a face table; the 3D Lagrange family has none
        // and is skipped rather than guessed at. 2D walks the corner ring
        // directly, which is exactly what cell_edges would have given.
        if (dim == 3 && detail::cell_faces(t).empty())
            continue;
        if (dim == 2 && nc < 3)
            continue;
        if (dim == 1 && nc < 2)
            continue;
        corners[b] = static_cast<std::size_t>(nc);
        eligible[b] = 1;
    }

    // The least-squares stencil needs cell-centred state for the WHOLE mesh
    // before any cell can be fitted, so it is prepared up front.
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t total = static_cast<std::size_t>(detail::total_cells(bases));
    detail::CellIncidence cell_nodes, node_cells;
    std::vector<double> cell_x, cell_f;
    std::vector<std::uint8_t> cell_ok;
    if (rOptions.mMethod == GradientMethod::LeastSquares) {
        cell_nodes = detail::build_cell_node_incidence(rMesh, npoints);
        node_cells = detail::invert_cell_node_incidence(cell_nodes, npoints, total);
        cell_x.assign(total * 3, 0.0);
        cell_f.assign(total * work_comp, 0.0);
        cell_ok.assign(total, 0);
        for (std::size_t b = 0; b < nblocks; ++b) {
            if (!eligible[b])
                continue;
            const auto cb = rMesh.Cells(b);
            const std::size_t base = static_cast<std::size_t>(bases[b]);
            const std::size_t ncells = cb.NumCells();
            parallel_for(ncells, [&](std::size_t c) {
                GradCell cell;
                grad_load_cell(cb, c, points, pdim, work, work_comp, npoints, corners[b], cell);
                if (!cell.mSupported)
                    return;
                const std::size_t g = base + c;
                cell_x[g * 3 + 0] = cell.mOrigin[0];
                cell_x[g * 3 + 1] = cell.mOrigin[1];
                cell_x[g * 3 + 2] = cell.mOrigin[2];
                for (std::size_t k = 0; k < work_comp; ++k)
                    cell_f[g * work_comp + k] = cell.mValue0[k];
                cell_ok[g] = 1;
            });
        }
    }

    // --- per-cell evaluation -------------------------------------------------
    std::vector<NDArray> blocks;
    blocks.reserve(nblocks);
    std::int64_t skipped = 0;
    std::int64_t fallback = 0;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t ncells = cb.NumCells();
        NDArray dst(DType::Float64, grad_shape(ncells, out_comp));
        double* pdst = dst.As<double>();
        // Per-cell flags reduced serially afterwards: an atomic counter would
        // work but this keeps the counts independent of scheduling and mirrors
        // data_average's "apply the policy serially" idiom.
        std::vector<std::uint8_t> flag_skipped(ncells, 0);
        std::vector<std::uint8_t> flag_fallback(ncells, 0);

        if (!eligible[b]) {
            for (std::size_t i = 0; i < ncells * out_comp; ++i)
                pdst[i] = std::nan("");
            skipped += static_cast<std::int64_t>(ncells);
            blocks.push_back(std::move(dst));
            continue;
        }

        const std::size_t base = static_cast<std::size_t>(bases[b]);
        const CellType type = types[b];
        const std::size_t ncorners = corners[b];
        parallel_for(ncells, [&](std::size_t c) {
            std::vector<double> grad(work_comp * 3, 0.0);
            GradCell cell;
            grad_load_cell(cb, c, points, pdim, work, work_comp, npoints, ncorners, cell);
            cell.mType = type;
            cell.mDim = dim;

            bool ok = false;
            if (cell.mSupported) {
                if (rOptions.mMethod == GradientMethod::LeastSquares) {
                    static thread_local std::vector<std::int64_t> nbrs;
                    detail::cell_node_neighbors(cell_nodes, node_cells, base + c, nbrs);
                    GradStencil stencil;
                    stencil.mOffsets.reserve(nbrs.size());
                    stencil.mDeltas.reserve(nbrs.size() * work_comp);
                    // Ascending global cell index: the neighbour order is part
                    // of cell_node_neighbors' contract precisely so this sum is
                    // reproducible.
                    for (std::int64_t nb : nbrs) {
                        const std::size_t g = static_cast<std::size_t>(nb);
                        if (!cell_ok[g])
                            continue;
                        stencil.mOffsets.push_back(Vec3{cell_x[g * 3 + 0] - cell.mOrigin[0],
                                                        cell_x[g * 3 + 1] - cell.mOrigin[1],
                                                        cell_x[g * 3 + 2] - cell.mOrigin[2]});
                        for (std::size_t k = 0; k < work_comp; ++k)
                            stencil.mDeltas.push_back(cell_f[g * work_comp + k] - cell.mValue0[k]);
                    }
                    Vec3 nrm{0.0, 0.0, 0.0};
                    const Vec3* p_normal = nullptr;
                    if (dim == 2) {
                        const Vec3 av = grad_ring_area_vector(cell.mCorners);
                        const double a = detail::vec3_norm(av);
                        if (a > 0.0) {
                            nrm = detail::vec3_scale(av, 1.0 / a);
                            p_normal = &nrm;
                        }
                    }
                    ok = grad_least_squares(stencil, work_comp, dim, p_normal, grad.data());
                    if (!ok) {
                        // Degenerate neighbourhood: an isolated cell, a
                        // collinear strip. Green-Gauss still has a usable
                        // answer, so take it and say so rather than emit NaN.
                        ok = grad_green_gauss(cell, dim, work_comp, grad.data());
                        if (ok)
                            flag_fallback[c] = 1;
                    }
                } else {
                    ok = grad_green_gauss(cell, dim, work_comp, grad.data());
                }
            }

            if (!ok) {
                for (std::size_t i = 0; i < out_comp; ++i)
                    pdst[c * out_comp + i] = std::nan("");
                flag_skipped[c] = 1;
                return;
            }
            grad_apply_operator(rOptions.mOperator, grad.data(), work_comp, pdst + c * out_comp);
        });

        for (std::size_t c = 0; c < ncells; ++c) {
            skipped += flag_skipped[c];
            fallback += flag_fallback[c];
        }
        blocks.push_back(std::move(dst));
    }

    out.AddCellData(out_name, std::move(blocks));

    if (rOptions.mLocation == DataLocation::Point) {
        // Compose the existing averaging rather than reimplementing a scatter:
        // its accumulation is already deliberately serial, so the point values
        // stay reproducible. All contributing cell values are equal for a linear
        // field, so the result is exact to within rounding -- but summing n
        // copies of g and dividing by n is not bit-exact in IEEE, which is why
        // the point-located tests use a tolerance.
        DataAverageOptions avg;
        avg.names = {out_name};
        avg.weight = CellPointWeight::Uniform;
        avg.overwrite = true;
        Mesh with_points = cell_data_to_point_data(out, avg);
        out = data_drop(with_points, DataLocation::Cell, {out_name});
    }

    return GradientResult{std::move(out), skipped, fallback};
}

}  // namespace meshioplusplus
