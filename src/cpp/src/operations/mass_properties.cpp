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
// Mass, centre of mass and inertia: see operations/mass_properties.hpp for the
// contract. Every cell contributes an exact integral of a constant density over
// tetrahedra (or triangles, for a shell), taken about one global reference
// point so the sums stay accurate far from the origin.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/mass_properties.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/sym3_eigen.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "../detail/cell_region_mask.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

// One cell's contribution about the reference point `o`:
// [0] mass, [1] signed measure, [2..4] first moment (mass * centroid),
// [5..10] second moment about o in xx yy zz xy yz zx order.
using MpContribution = std::array<double, 11>;

constexpr std::size_t kMpNoBlock = static_cast<std::size_t>(-1);

void mp_add_second_moment(MpContribution& rC, double Weight, const Vec3& rSum,
                          const std::array<Vec3, 4>& rP, std::size_t N) {
    // Weight * (sum_a p_a p_a^T + s s^T), the exact second moment of a simplex
    // for a constant density once Weight carries its mass over (N + 1) * (N + 2).
    double xx = rSum[0] * rSum[0];
    double yy = rSum[1] * rSum[1];
    double zz = rSum[2] * rSum[2];
    double xy = rSum[0] * rSum[1];
    double yz = rSum[1] * rSum[2];
    double zx = rSum[2] * rSum[0];
    for (std::size_t a = 0; a < N; ++a) {
        xx += rP[a][0] * rP[a][0];
        yy += rP[a][1] * rP[a][1];
        zz += rP[a][2] * rP[a][2];
        xy += rP[a][0] * rP[a][1];
        yz += rP[a][1] * rP[a][2];
        zx += rP[a][2] * rP[a][0];
    }
    rC[5] += Weight * xx;
    rC[6] += Weight * yy;
    rC[7] += Weight * zz;
    rC[8] += Weight * xy;
    rC[9] += Weight * yz;
    rC[10] += Weight * zx;
}

// A signed tetrahedron (rA, rB, rC, rD) of density Rho, added to rOut.
void mp_add_tet(MpContribution& rOut, double Rho, const Vec3& rA, const Vec3& rB, const Vec3& rC,
                const Vec3& rD) {
    const double vol = detail::triple_product(detail::vec3_sub(rB, rA), detail::vec3_sub(rC, rA),
                                              detail::vec3_sub(rD, rA)) /
                       6.0;
    const double mass = Rho * vol;
    const std::array<Vec3, 4> p{rA, rB, rC, rD};
    const Vec3 s = detail::vec3_add(detail::vec3_add(rA, rB), detail::vec3_add(rC, rD));
    rOut[0] += mass;
    rOut[1] += vol;
    for (int k = 0; k < 3; ++k)
        rOut[2 + k] += mass * s[static_cast<std::size_t>(k)] / 4.0;
    mp_add_second_moment(rOut, mass / 20.0, s, p, 4);
}

// A triangle of area `Area` and areal density Sigma (a shell).
void mp_add_triangle(MpContribution& rOut, double Sigma, const Vec3& rA, const Vec3& rB,
                     const Vec3& rC) {
    const double area = 0.5 * detail::vec3_norm(detail::vec3_cross(detail::vec3_sub(rB, rA),
                                                                   detail::vec3_sub(rC, rA)));
    const double mass = Sigma * area;
    const std::array<Vec3, 4> p{rA, rB, rC, Vec3{0.0, 0.0, 0.0}};
    const Vec3 s = detail::vec3_add(detail::vec3_add(rA, rB), rC);
    rOut[0] += mass;
    rOut[1] += area;
    for (int k = 0; k < 3; ++k)
        rOut[2 + k] += mass * s[static_cast<std::size_t>(k)] / 3.0;
    mp_add_second_moment(rOut, mass / 12.0, s, p, 3);
}

// What a block can offer: the cell dimension that integrates (3 solid, 2
// surface), or 0 when it cannot (lines, points, ragged polygons, higher-order
// 3D types without a face table are found per cell).
int mp_block_dimension(const Mesh::CellView& rBlock) {
    if (rBlock.IsPolyhedron())
        return 3;
    if (rBlock.IsRagged())
        return 0;
    const CellType ct = cell_type_from_name(rBlock.Type());
    if (detail::cell_corner_count(ct) <= 0)
        return 0;
    const int dim = cell_type_dimension(ct);
    return (dim == 2 || dim == 3) ? dim : 0;
}

struct MpBody {
    const Mesh* mpMesh = nullptr;
    std::vector<std::int64_t> mBases;
    std::vector<int> mBlockDim;
    Vec3 mOrigin{0.0, 0.0, 0.0};
    int mDim = 0;  // 3 solid cells, 2 surface cells
    bool mShell = false;
    double mThickness = 0.0;
    double mDensity = 1.0;
    std::vector<const NDArray*> mDensityArrays;  // one per block, or empty
    // Per global cell: contribution, and whether it exists.
    std::vector<MpContribution> mContrib;
    std::vector<std::uint8_t> mUsable;
};

double mp_cell_density(const MpBody& rBody, std::size_t Block, std::size_t Row) {
    if (rBody.mDensityArrays.empty())
        return rBody.mDensity;
    return detail::read_double(*rBody.mDensityArrays[Block], Row);
}

// A solid cell: the corner-average fan, oriented outward.
bool mp_solid_cell(const MpBody& rBody, std::size_t Block, std::size_t Row, MpContribution& rOut) {
    const Mesh& mesh = *rBody.mpMesh;
    const auto cb = mesh.Cells(Block);
    detail::CellRings rings;
    std::vector<Vec3> coords;
    if (!detail::cell_rings(cb, Row, mesh.Points(), mesh.PointDim(), rings, coords))
        return false;
    if (detail::orient_rings(rings, coords.data()) == detail::RingOrientation::Unorientable)
        return false;
    const std::size_t nn = rings.mNodes.size();
    if (nn == 0)
        return false;
    Vec3 apex{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < nn; ++i) {
        coords[i] = detail::vec3_sub(coords[i], rBody.mOrigin);
        apex = detail::vec3_add(apex, coords[i]);
    }
    apex = detail::vec3_scale(apex, 1.0 / static_cast<double>(nn));
    const double rho = mp_cell_density(rBody, Block, Row);
    rOut.fill(0.0);
    for (std::size_t f = 0; f < rings.NumFaces(); ++f) {
        const std::size_t n = rings.FaceSize(f);
        const std::uint32_t* ring = rings.Face(f);
        if (n < 3)
            continue;
        const Vec3 fc = detail::polygon_centroid(coords.data(), ring, n);
        for (std::size_t i = 0; i < n; ++i)
            mp_add_tet(rOut, rho, apex, fc, coords[ring[i]], coords[ring[(i + 1) % n]]);
    }
    return true;
}

// A surface cell: a shell triangle fan, or (closed body) signed tets to the origin.
bool mp_surface_cell(const MpBody& rBody, std::size_t Block, std::size_t Row,
                     MpContribution& rOut) {
    const Mesh& mesh = *rBody.mpMesh;
    const auto cb = mesh.Cells(Block);
    const CellType ct = cell_type_from_name(cb.Type());
    const std::size_t corners = static_cast<std::size_t>(detail::cell_corner_count(ct));
    std::vector<Vec3> c;
    detail::read_corner_coords(mesh.Points(), mesh.PointDim(), cb.Conn(), Row * cb.NodesPerCell(),
                               corners, c);
    for (Vec3& p : c)
        p = detail::vec3_sub(p, rBody.mOrigin);
    rOut.fill(0.0);
    if (corners < 3)
        return false;
    Vec3 centre{0.0, 0.0, 0.0};  // the corner average, the fan apex of poly_measure
    for (const Vec3& p : c)
        centre = detail::vec3_add(centre, p);
    centre = detail::vec3_scale(centre, 1.0 / static_cast<double>(corners));
    if (rBody.mShell) {
        const double sigma = mp_cell_density(rBody, Block, Row) * rBody.mThickness;
        for (std::size_t i = 0; i < corners; ++i)
            mp_add_triangle(rOut, sigma, centre, c[i], c[(i + 1) % corners]);
    } else {
        const double rho = mp_cell_density(rBody, Block, Row);
        const Vec3 zero{0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < corners; ++i)
            mp_add_tet(rOut, rho, zero, centre, c[i], c[(i + 1) % corners]);
    }
    return true;
}

// Whether the selected surface cells form a closed, consistently wound shell:
// every directed edge once, with its reverse once.
bool mp_closed_and_consistent(const MpBody& rBody, const std::vector<std::uint8_t>& rSel) {
    const Mesh& mesh = *rBody.mpMesh;
    std::vector<std::pair<std::int64_t, std::int64_t>> fwd;
    for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b) {
        if (rBody.mBlockDim[b] != 2)
            continue;
        const auto cb = mesh.Cells(b);
        const std::size_t corners =
            static_cast<std::size_t>(detail::cell_corner_count(cell_type_from_name(cb.Type())));
        const std::size_t stride = cb.NodesPerCell();
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            const std::size_t g = static_cast<std::size_t>(rBody.mBases[b]) + r;
            if (!rSel[g] || !rBody.mUsable[g])
                continue;
            for (std::size_t i = 0; i < corners; ++i) {
                const std::int64_t a = detail::read_int(cb.Conn(), r * stride + i);
                const std::int64_t d = detail::read_int(cb.Conn(), r * stride + (i + 1) % corners);
                fwd.emplace_back(a, d);
            }
        }
    }
    if (fwd.empty())
        return false;
    std::vector<std::pair<std::int64_t, std::int64_t>> rev;
    rev.reserve(fwd.size());
    for (const auto& e : fwd)
        rev.emplace_back(e.second, e.first);
    std::sort(fwd.begin(), fwd.end());
    std::sort(rev.begin(), rev.end());
    if (std::adjacent_find(fwd.begin(), fwd.end()) != fwd.end())
        return false;
    return fwd == rev;
}

void mp_inertia_to_tensor(const double* pS, double* pI) {
    // Second moment S -> inertia: Ixx = Syy + Szz, ..., Ixy = -Sxy, ...
    pI[0] = pS[1] + pS[2];
    pI[1] = pS[0] + pS[2];
    pI[2] = pS[0] + pS[1];
    pI[3] = -pS[3];
    pI[4] = -pS[4];
    pI[5] = -pS[5];
}

// Fold the selected cells' contributions (serial, in cell order).
MassProperties mp_reduce(const MpBody& rBody, const std::vector<std::uint8_t>& rSel,
                         MassAbout About, std::string Name) {
    MpContribution sum{};
    std::int64_t used = 0;
    for (std::size_t g = 0; g < rSel.size(); ++g) {
        if (!rSel[g] || !rBody.mUsable[g])
            continue;
        for (std::size_t k = 0; k < sum.size(); ++k)
            sum[k] += rBody.mContrib[g][k];
        ++used;
    }
    if (rBody.mDim == 2 && !rBody.mShell && sum[1] < 0.0)
        for (double& v : sum)
            v = -v;  // an inward-wound closed shell: turn it outward

    MassProperties out;
    out.mName = std::move(Name);
    out.mNumCells = used;
    out.mMeasure = sum[1];
    out.mMass = sum[0];
    const double nan = std::nan("");
    const Vec3 m1{sum[2], sum[3], sum[4]};
    const Vec3 o = rBody.mOrigin;
    // Centre of mass.
    if (sum[0] != 0.0)
        out.mCentre = detail::vec3_add(o, detail::vec3_scale(m1, 1.0 / sum[0]));
    else
        out.mCentre = {nan, nan, nan};
    // Second moment about the chosen point, from the one about o.
    double s[6] = {sum[5], sum[6], sum[7], sum[8], sum[9], sum[10]};
    if (About == MassAbout::Centroid) {
        if (sum[0] != 0.0) {
            const double inv = 1.0 / sum[0];
            s[0] -= m1[0] * m1[0] * inv;
            s[1] -= m1[1] * m1[1] * inv;
            s[2] -= m1[2] * m1[2] * inv;
            s[3] -= m1[0] * m1[1] * inv;
            s[4] -= m1[1] * m1[2] * inv;
            s[5] -= m1[2] * m1[0] * inv;
        }
    } else {
        // About the coordinate origin: S0 = S_o + o m1^T + m1 o^T + m o o^T.
        s[0] += 2.0 * o[0] * m1[0] + sum[0] * o[0] * o[0];
        s[1] += 2.0 * o[1] * m1[1] + sum[0] * o[1] * o[1];
        s[2] += 2.0 * o[2] * m1[2] + sum[0] * o[2] * o[2];
        s[3] += o[0] * m1[1] + o[1] * m1[0] + sum[0] * o[0] * o[1];
        s[4] += o[1] * m1[2] + o[2] * m1[1] + sum[0] * o[1] * o[2];
        s[5] += o[2] * m1[0] + o[0] * m1[2] + sum[0] * o[2] * o[0];
    }
    mp_inertia_to_tensor(s, out.mInertia.data());
    detail::sym3_eigen(out.mInertia.data(), out.mPrincipalMoments.data(),
                       out.mPrincipalAxes.data());
    for (std::size_t k = 0; k < 3; ++k) {
        const double ratio = sum[0] != 0.0 ? out.mPrincipalMoments[k] / sum[0] : nan;
        out.mGyrationRadii[k] = ratio >= 0.0 ? std::sqrt(ratio) : nan;
    }
    return out;
}

}  // namespace

MassAbout mass_about_from_name(const std::string& rName) {
    if (rName == "centroid" || rName == "centre" || rName == "center")
        return MassAbout::Centroid;
    if (rName == "origin")
        return MassAbout::Origin;
    throw std::invalid_argument("meshio++: mass_properties: unknown reference point '" + rName +
                                "' (expected 'centroid' or 'origin')");
}

MassPropertiesReport mass_properties(const Mesh& rMesh, const MassPropertiesOptions& rOptions) {
    static const char* const kPrefix = "meshio++: mass_properties: ";
    if (!std::isfinite(rOptions.mDensity))
        throw std::invalid_argument(std::string(kPrefix) + "density must be finite");
    if (!(rOptions.mThickness >= 0.0) || !std::isfinite(rOptions.mThickness))
        throw std::invalid_argument(std::string(kPrefix) +
                                    "thickness must be a finite number, zero or above");

    MassPropertiesReport report;
    MpBody body;
    body.mpMesh = &rMesh;
    body.mBases = detail::block_bases(rMesh);
    const std::size_t nblocks = rMesh.NumCellBlocks();
    const std::size_t ncells = static_cast<std::size_t>(detail::total_cells(body.mBases));
    body.mDensity = rOptions.mDensity;
    body.mThickness = rOptions.mThickness;

    // Which cells take part, and at what dimension.
    const std::vector<std::uint8_t> selection =
        detail::cell_region_mask(rMesh, rOptions.mRegion, ncells, kPrefix);
    body.mBlockDim.assign(nblocks, 0);
    int dim = 0;
    for (std::size_t b = 0; b < nblocks; ++b) {
        body.mBlockDim[b] = mp_block_dimension(rMesh.Cells(b));
        bool any = false;
        for (std::size_t r = 0; r < rMesh.Cells(b).NumCells() && !any; ++r)
            any = selection[static_cast<std::size_t>(body.mBases[b]) + r] != 0;
        if (any)
            dim = std::max(dim, body.mBlockDim[b]);
    }
    if (dim == 0)
        throw std::invalid_argument(std::string(kPrefix) +
                                    "no triangle, polygon or volume cell to integrate");
    body.mDim = dim;
    body.mShell = dim == 2 && rOptions.mThickness > 0.0;
    if (dim == 3 && rOptions.mThickness > 0.0)
        report.mWarnings.push_back("thickness ignored: the mesh has volume cells");
    report.mKind = (dim == 2 && body.mShell) ? "shell" : "solid";

    if (!rOptions.mDensityArray.empty()) {
        if (dim == 2 && !body.mShell)
            throw std::invalid_argument(std::string(kPrefix) +
                                        "a density array needs volume cells or a shell: a closed "
                                        "surface has no cell of its own to carry a density");
        if (!rMesh.HasCellData(rOptions.mDensityArray))
            throw std::invalid_argument(std::string(kPrefix) + "no cell_data array named '" +
                                        rOptions.mDensityArray + "'");
        if (rMesh.CellDataNumBlocks(rOptions.mDensityArray) != nblocks)
            throw std::invalid_argument(std::string(kPrefix) + "cell_data '" +
                                        rOptions.mDensityArray +
                                        "' does not have one array per cell block");
        for (std::size_t b = 0; b < nblocks; ++b) {
            const NDArray& a = rMesh.CellData(rOptions.mDensityArray, b);
            if (a.Size() != rMesh.Cells(b).NumCells())
                throw std::invalid_argument(std::string(kPrefix) + "density array '" +
                                            rOptions.mDensityArray +
                                            "' must hold one scalar per cell");
            body.mDensityArrays.push_back(&a);
        }
    }

    // Reference point: the mean of the points, so offsets stay small.
    {
        const NDArray& pts = rMesh.Points();
        const std::size_t n = rMesh.NumPoints();
        const std::size_t pd = rMesh.PointDim();
        Vec3 sum{0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < n; ++i)
            sum = detail::vec3_add(sum, detail::read_point(pts, pd, static_cast<std::int64_t>(i)));
        if (n > 0)
            body.mOrigin = detail::vec3_scale(sum, 1.0 / static_cast<double>(n));
        if (!std::isfinite(body.mOrigin[0] + body.mOrigin[1] + body.mOrigin[2]))
            throw std::invalid_argument(std::string(kPrefix) + "a point coordinate is not finite");
    }

    // Per-cell contributions, in parallel; every cell is independent.
    body.mContrib.assign(ncells, MpContribution{});
    body.mUsable.assign(ncells, 0);
    std::int64_t skipped = 0;
    std::int64_t lower = 0;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t base = static_cast<std::size_t>(body.mBases[b]);
        const std::size_t n = cb.NumCells();
        if (body.mBlockDim[b] != dim) {
            for (std::size_t r = 0; r < n; ++r)
                lower += selection[base + r] ? 1 : 0;
            continue;
        }
        parallel_for(n, [&](std::size_t r) {
            if (!selection[base + r])
                return;
            const bool ok = dim == 3 ? mp_solid_cell(body, b, r, body.mContrib[base + r])
                                     : mp_surface_cell(body, b, r, body.mContrib[base + r]);
            body.mUsable[base + r] = ok ? 1 : 0;
        });
        for (std::size_t r = 0; r < n; ++r)
            if (selection[base + r] && !body.mUsable[base + r])
                ++skipped;
    }
    if (lower > 0)
        report.mWarnings.push_back(std::to_string(lower) + " cell(s) of lower dimension ignored");
    if (skipped > 0)
        report.mWarnings.push_back(std::to_string(skipped) +
                                   " cell(s) skipped: no face table or an unorientable polyhedron");
    if (std::none_of(body.mUsable.begin(), body.mUsable.end(),
                     [](std::uint8_t u) { return u != 0; }))
        throw std::invalid_argument(std::string(kPrefix) + "no cell could be integrated");

    const auto check_closed = [&](const std::vector<std::uint8_t>& rSel, const std::string& rWhat) {
        if (dim == 2 && !body.mShell && !mp_closed_and_consistent(body, rSel))
            throw std::invalid_argument(
                std::string(kPrefix) + rWhat +
                " is not a closed, consistently wound surface; close it with `repair` or give a "
                "thickness to treat it as a shell");
    };
    check_closed(selection,
                 rOptions.mRegion.empty() ? "the surface" : "region '" + rOptions.mRegion + "'");
    report.mTotal = mp_reduce(body, selection, rOptions.mAbout, "");

    if (rOptions.mRegion.empty()) {
        const std::vector<std::int64_t>& bases = body.mBases;
        for (const std::string& name : rMesh.RegionNames()) {
            const std::size_t idx = rMesh.FindRegion(name, RegionKind::Cell);
            if (idx == Mesh::npos)
                continue;
            std::vector<std::uint8_t> sel(ncells, 0);
            const NDArray& entries = rMesh.Region(idx).mEntries;
            for (std::size_t e = 0; e < entries.Size(); ++e) {
                const std::int64_t g = detail::read_int(entries, e);
                if (g >= 0 && g < bases.back())
                    sel[static_cast<std::size_t>(g)] = 1;
            }
            // A region whose cells all have another dimension has nothing here.
            bool any = false;
            for (std::size_t g = 0; g < ncells && !any; ++g)
                any = sel[g] && body.mUsable[g];
            if (!any)
                continue;
            check_closed(sel, "region '" + name + "'");
            report.mRegions.push_back(mp_reduce(body, sel, rOptions.mAbout, name));
        }
    }
    return report;
}

}  // namespace meshioplusplus
