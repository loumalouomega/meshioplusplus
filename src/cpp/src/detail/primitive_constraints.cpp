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
// Fit related primitives together. See detail/primitive_constraints.hpp.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/segment.hpp"

// Project includes (private, not installed)
#include "lsq_fit.hpp"
#include "primitive_constraints.hpp"
#include "primitive_geom.hpp"
#include "surface_rings.hpp"

namespace meshioplusplus {
namespace detail {
namespace {

constexpr const char* kPcPrefix = "fit_cad: ";
// The most samples a joint fit reads; a larger group is read with a stride.
constexpr std::size_t kPcMaxJointSamples = 20000;
// How far, in multiples of a patch's tolerance, a shared axis or normal may
// move the patch for it to be tried.
constexpr double kPcReach = 2.0;

double pc_dot(const double* pA, const double* pB) {
    return pA[0] * pB[0] + pA[1] * pB[1] + pA[2] * pB[2];
}

void pc_cross(const double* pA, const double* pB, double* pOut) {
    pOut[0] = pA[1] * pB[2] - pA[2] * pB[1];
    pOut[1] = pA[2] * pB[0] - pA[0] * pB[2];
    pOut[2] = pA[0] * pB[1] - pA[1] * pB[0];
}

bool pc_normalize(double* pV) {
    const double n = std::sqrt(pc_dot(pV, pV));
    if (!(n > 0.0) || !std::isfinite(n))
        return false;
    for (int k = 0; k < 3; ++k)
        pV[k] /= n;
    return true;
}

// The sign convention of `fit:patch_params`: the largest-magnitude component
// (the first on a tie) positive.
void pc_canonical(double* pD) {
    int k = 0;
    for (int j = 1; j < 3; ++j)
        if (std::fabs(pD[j]) > std::fabs(pD[k]))
            k = j;
    if (pD[k] < 0.0)
        for (int j = 0; j < 3; ++j)
            pD[j] = -pD[j];
}

// The sine of the angle between two unit directions, sign ignored.
double pc_sin(const double* pA, const double* pB) {
    double c[3];
    pc_cross(pA, pB, c);
    return std::sqrt(pc_dot(c, c));
}

struct PcStats {
    double mRms = 0.0;
    double mMax = 0.0;
    double mFraction = 0.0;
};

struct PcPatch {
    std::int64_t mType = 0;
    double mParams[kPrimitiveParams] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double mOriginal[kPrimitiveParams] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> mPts;  // the patch's points, three per point
    double mCentroid[3] = {0.0, 0.0, 0.0};
    double mExtent = 0.0;  // the diagonal of the points' bounding box
    double mTol = 0.0;
    PcStats mStats;  // of mParams

    std::size_t NumPts() const { return mPts.size() / 3; }
    bool Axial() const { return mType == 3 || mType == 4 || mType == 5; }
};

PcStats pc_stats(std::int64_t Type, const double* pParams, const std::vector<double>& rPts,
                 double Tol) {
    PcStats s;
    const std::size_t n = rPts.size() / 3;
    if (n == 0)
        return s;
    double sum = 0.0;
    std::size_t within = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double d = primitive_point_distance(Type, pParams, &rPts[i * 3]);
        sum += d * d;
        s.mMax = std::max(s.mMax, d);
        if (d <= Tol)
            ++within;
    }
    s.mRms = std::sqrt(sum / static_cast<double>(n));
    s.mFraction = static_cast<double>(within) / static_cast<double>(n);
    return s;
}

// Whether a constrained fit still passes: the test the patch was accepted
// under (its tolerance, and no smaller share of points within it than the
// default or the patch had), and a largest distance no writer demotes at
// unless the patch's own fit already exceeded it.
bool pc_accept(const PcPatch& rPatch, const double* pParams, double FaceTol, double MinFraction,
               PcStats& rOut) {
    for (int k = 0; k < kPrimitiveParams; ++k)
        if (!std::isfinite(pParams[k]))
            return false;
    rOut = pc_stats(rPatch.mType, pParams, rPatch.mPts, rPatch.mTol);
    if (!std::isfinite(rOut.mRms) || !std::isfinite(rOut.mMax))
        return false;
    return rOut.mRms <= rPatch.mTol &&
           rOut.mFraction >= std::min(MinFraction, rPatch.mStats.mFraction) &&
           (rOut.mMax <= FaceTol || rOut.mMax <= rPatch.mStats.mMax);
}

// An axial primitive's axis point (a cylinder's axis point, a cone's apex, a
// torus's or sphere's centre) moved along its axis to the height of the
// patch's centroid.
void pc_axis_point(const PcPatch& rPatch, double* pOut) {
    const double* x = rPatch.mParams;
    const double* d = rPatch.mParams + 3;
    const double v[3] = {rPatch.mCentroid[0] - x[0], rPatch.mCentroid[1] - x[1],
                         rPatch.mCentroid[2] - x[2]};
    const double h = pc_dot(v, d);
    for (int k = 0; k < 3; ++k)
        pOut[k] = x[k] + h * d[k];
}

// The distance from pP to the line through pX along the unit pD.
double pc_line_distance(const double* pX, const double* pD, const double* pP) {
    const double v[3] = {pP[0] - pX[0], pP[1] - pX[1], pP[2] - pX[2]};
    const double h = pc_dot(v, pD);
    const double r[3] = {v[0] - h * pD[0], v[1] - h * pD[1], v[2] - h * pD[2]};
    return std::sqrt(pc_dot(r, r));
}

// Whether a patch's primitive can take the axis (pX, pD): an axial one whose
// axis is parallel to it and lies on it, a sphere whose centre is on it, each
// within what moves the patch by kPcReach of its tolerance.
bool pc_on_axis(const PcPatch& rPatch, const double* pX, const double* pD) {
    const double reach = kPcReach * rPatch.mTol;
    if (rPatch.mType == 2)
        return pc_line_distance(pX, pD, rPatch.mParams) <= reach;
    if (!rPatch.Axial())
        return false;
    if (pc_sin(pD, rPatch.mParams + 3) * rPatch.mExtent > reach)
        return false;
    double p[3];
    pc_axis_point(rPatch, p);
    return pc_line_distance(pX, pD, p) <= reach;
}

// One member of a coaxial group, parametrised along the shared axis.
struct PcMember {
    std::size_t mRow = 0;
    double mSign = 1.0;  // a cone's axis along (+1) or against (-1) the shared one
    std::size_t mFirst = 0;  // its first parameter in the joint vector
};

std::size_t pc_member_params(std::int64_t Type) {
    switch (Type) {
        case 2:
            return 2;  // centre height, radius
        case 3:
            return 1;  // radius
        case 4:
            return 2;  // apex height, half-angle
        case 5:
            return 3;  // centre height, major and minor radius
        default:
            return 0;
    }
}

// The primitive of a member given the shared axis (pP, pD) and its own
// parameters pTheta.
void pc_member_build(const PcPatch& rPatch, const PcMember& rMember, const double* pP,
                     const double* pD, const double* pTheta, double* pOut) {
    std::fill(pOut, pOut + kPrimitiveParams, 0.0);
    double axis[3] = {pD[0], pD[1], pD[2]};
    switch (rPatch.mType) {
        case 2:
            for (int k = 0; k < 3; ++k)
                pOut[k] = pP[k] + pTheta[0] * pD[k];
            pOut[6] = pTheta[1];
            return;
        case 3: {
            const double v[3] = {rPatch.mCentroid[0] - pP[0], rPatch.mCentroid[1] - pP[1],
                                 rPatch.mCentroid[2] - pP[2]};
            const double h = pc_dot(v, pD);
            pc_canonical(axis);
            for (int k = 0; k < 3; ++k) {
                pOut[k] = pP[k] + h * pD[k];
                pOut[3 + k] = axis[k];
            }
            pOut[6] = pTheta[0];
            return;
        }
        case 4:
            for (int k = 0; k < 3; ++k) {
                pOut[k] = pP[k] + pTheta[0] * pD[k];
                pOut[3 + k] = rMember.mSign * pD[k];
            }
            pOut[8] = pTheta[1];
            return;
        case 5:
            pc_canonical(axis);
            for (int k = 0; k < 3; ++k) {
                pOut[k] = pP[k] + pTheta[0] * pD[k];
                pOut[3 + k] = axis[k];
            }
            pOut[6] = pTheta[1];
            pOut[7] = pTheta[2];
            return;
        default:
            return;
    }
}

// A member's own parameters read off its primitive, its centre or apex
// projected on the axis (pP, pD).
void pc_member_start(const PcPatch& rPatch, const double* pP, const double* pD, double* pTheta) {
    const double* x = rPatch.mParams;
    const double v[3] = {x[0] - pP[0], x[1] - pP[1], x[2] - pP[2]};
    const double h = pc_dot(v, pD);
    switch (rPatch.mType) {
        case 2:
            pTheta[0] = h;
            pTheta[1] = x[6];
            return;
        case 3:
            pTheta[0] = x[6];
            return;
        case 4:
            pTheta[0] = h;
            pTheta[1] = x[8];
            return;
        case 5:
            pTheta[0] = h;
            pTheta[1] = x[6];
            pTheta[2] = x[7];
            return;
        default:
            return;
    }
}

// Fit the rows of a coaxial group jointly: one axis, each member's own
// parameters along it. Writes each member's constrained primitive to rOut.
void pc_fit_group(const std::vector<PcPatch>& rPatches, const std::vector<std::size_t>& rRows,
                  std::vector<std::vector<double>>& rOut) {
    // The start: the point-weighted mean of the axial members' axes, signed
    // alike, through the weighted mean of their axis points.
    double d0[3] = {0.0, 0.0, 0.0}, p0[3] = {0.0, 0.0, 0.0}, w = 0.0;
    const double* seed = nullptr;
    double extent = 0.0;
    for (std::size_t r : rRows) {
        const PcPatch& p = rPatches[r];
        extent = std::max(extent, p.mExtent);
        if (!p.Axial())
            continue;
        const double* a = p.mParams + 3;
        if (!seed)
            seed = a;
        const double s = pc_dot(a, seed) < 0.0 ? -1.0 : 1.0;
        const double wi = static_cast<double>(p.NumPts());
        double x[3];
        pc_axis_point(p, x);
        for (int k = 0; k < 3; ++k) {
            d0[k] += wi * s * a[k];
            p0[k] += wi * x[k];
        }
        w += wi;
    }
    rOut.assign(rRows.size(), std::vector<double>());
    if (!seed || !(w > 0.0) || !pc_normalize(d0))
        return;
    for (int k = 0; k < 3; ++k)
        p0[k] /= w;
    double e1[3], e2[3];
    lsq_frame(d0, e1, e2);

    std::vector<PcMember> members(rRows.size());
    std::vector<double> theta = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> typical = {1.0, 1.0, extent > 0.0 ? extent : 1.0,
                                   extent > 0.0 ? extent : 1.0};
    for (std::size_t i = 0; i < rRows.size(); ++i) {
        const PcPatch& p = rPatches[rRows[i]];
        PcMember& m = members[i];
        m.mRow = rRows[i];
        m.mSign = (p.mType == 4 && pc_dot(p.mParams + 3, d0) < 0.0) ? -1.0 : 1.0;
        m.mFirst = theta.size();
        const std::size_t np = pc_member_params(p.mType);
        double start[3] = {0.0, 0.0, 0.0};
        pc_member_start(p, p0, d0, start);
        for (std::size_t k = 0; k < np; ++k) {
            theta.push_back(start[k]);
            // A cone's half-angle is in degrees; everything else a length.
            const bool angle = p.mType == 4 && k == 1;
            typical.push_back(angle ? 10.0 : (extent > 0.0 ? extent : 1.0));
        }
    }

    // The samples: every point of every member, read with one stride when the
    // group holds more than kPcMaxJointSamples.
    std::size_t total = 0;
    for (std::size_t r : rRows)
        total += rPatches[r].NumPts();
    const std::size_t stride =
        total > kPcMaxJointSamples ? (total + kPcMaxJointSamples - 1) / kPcMaxJointSamples : 1;
    std::vector<std::pair<std::size_t, std::size_t>> samples;  // (member, point)
    for (std::size_t i = 0; i < rRows.size(); ++i)
        for (std::size_t j = 0; j < rPatches[rRows[i]].NumPts(); j += stride)
            samples.emplace_back(i, j);

    const auto axis_of = [&](const std::vector<double>& rTheta, double* pP, double* pD) {
        for (int k = 0; k < 3; ++k) {
            pD[k] = d0[k] + rTheta[0] * e1[k] + rTheta[1] * e2[k];
            pP[k] = p0[k] + rTheta[2] * e1[k] + rTheta[3] * e2[k];
        }
        pc_normalize(pD);
    };
    const LsqResidual residual = [&](const std::vector<double>& rTheta, double* pOut) {
        double p[3], d[3];
        axis_of(rTheta, p, d);
        std::vector<std::vector<double>> prim(rRows.size(), std::vector<double>(kPrimitiveParams));
        for (std::size_t i = 0; i < rRows.size(); ++i)
            pc_member_build(rPatches[rRows[i]], members[i], p, d, &rTheta[members[i].mFirst],
                            prim[i].data());
        for (std::size_t s = 0; s < samples.size(); ++s) {
            const std::size_t i = samples[s].first;
            const PcPatch& patch = rPatches[rRows[i]];
            pOut[s] = primitive_point_distance(patch.mType, prim[i].data(),
                                               &patch.mPts[samples[s].second * 3]);
        }
    };
    lsq_levenberg_marquardt(theta, typical, samples.size(), residual, 50);

    double p[3], d[3];
    axis_of(theta, p, d);
    for (std::size_t i = 0; i < rRows.size(); ++i) {
        rOut[i].assign(kPrimitiveParams, 0.0);
        pc_member_build(rPatches[rRows[i]], members[i], p, d, &theta[members[i].mFirst],
                        rOut[i].data());
    }
}

}  // namespace

PrimitiveConstraintReport constrain_primitives(Mesh& rMesh, double FaceTolerance) {
    PrimitiveConstraintReport report;
    for (const char* name : {"fit:patch_id", "fit:patch_type", "fit:patch_params",
                             "fit:patch_rms", "fit:patch_max", "fit:patch_fraction",
                             "fit:patch_tolerance"})
        if (!rMesh.HasFieldData(name))
            return report;
    if (!rMesh.HasCellData(kSegmentPatchName))
        return report;

    const NDArray& ids = rMesh.FieldData("fit:patch_id");
    const NDArray& types = rMesh.FieldData("fit:patch_type");
    const NDArray& params = rMesh.FieldData("fit:patch_params");
    const NDArray& tols = rMesh.FieldData("fit:patch_tolerance");
    const NDArray& rms = rMesh.FieldData("fit:patch_rms");
    const NDArray& maxs = rMesh.FieldData("fit:patch_max");
    const NDArray& fracs = rMesh.FieldData("fit:patch_fraction");
    const std::size_t m = ids.Size();
    if (types.Size() != m || params.Size() != m * kPrimitiveParams || tols.Size() != m ||
        rms.Size() != m || maxs.Size() != m || fracs.Size() != m)
        return report;

    std::map<std::int64_t, std::size_t> row_of;
    for (std::size_t r = 0; r < m; ++r)
        row_of.emplace(read_int(ids, r), r);

    // Each row's patch: its points (each once), centroid and extent.
    const std::vector<std::int64_t> bases = block_bases(rMesh);
    const std::size_t ncells = static_cast<std::size_t>(total_cells(bases));
    std::vector<std::int64_t> cell_row(ncells, -1);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const NDArray& a = rMesh.CellData(kSegmentPatchName, b);
        const std::size_t nc = rMesh.Cells(b).NumCells();
        const std::size_t base = static_cast<std::size_t>(bases[b]);
        for (std::size_t c = 0; c < nc && c < a.Size(); ++c) {
            auto it = row_of.find(read_int(a, c));
            if (it != row_of.end())
                cell_row[base + c] = static_cast<std::int64_t>(it->second);
        }
    }
    const std::vector<std::uint8_t> all(ncells, 1);
    const SurfaceRings rings = surface_rings(rMesh, all, kPcPrefix);
    const std::vector<double> xyz = surface_xyz(rMesh);
    std::vector<std::vector<std::int64_t>> verts(m);
    for (std::size_t f = 0; f < rings.NumFaces(); ++f) {
        const std::int64_t r = cell_row[static_cast<std::size_t>(rings.mCell[f])];
        if (r < 0)
            continue;
        for (std::int64_t k = rings.mStart[f]; k < rings.mStart[f + 1]; ++k)
            verts[static_cast<std::size_t>(r)].push_back(rings.mNodes[static_cast<std::size_t>(k)]);
    }

    std::vector<PcPatch> patches(m);
    for (std::size_t r = 0; r < m; ++r) {
        PcPatch& p = patches[r];
        p.mType = read_int(types, r);
        bool finite = true;
        for (int k = 0; k < kPrimitiveParams; ++k) {
            p.mParams[k] = p.mOriginal[k] =
                read_double(params, r * kPrimitiveParams + static_cast<std::size_t>(k));
            finite = finite && std::isfinite(p.mParams[k]);
        }
        p.mTol = read_double(tols, r);
        p.mStats.mRms = read_double(rms, r);
        p.mStats.mMax = read_double(maxs, r);
        p.mStats.mFraction = read_double(fracs, r);
        if (p.mType < 1 || p.mType > 5 || !finite || !(p.mTol > 0.0)) {
            p.mType = 0;  // nothing to constrain
            continue;
        }
        std::vector<std::int64_t>& v = verts[r];
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        p.mPts.resize(v.size() * 3);
        double lo[3] = {0.0, 0.0, 0.0}, hi[3] = {0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < v.size(); ++i)
            for (std::size_t k = 0; k < 3; ++k) {
                const double x = xyz[static_cast<std::size_t>(v[i]) * 3 + k];
                p.mPts[i * 3 + k] = x;
                p.mCentroid[k] += x;
                lo[k] = i == 0 ? x : std::min(lo[k], x);
                hi[k] = i == 0 ? x : std::max(hi[k], x);
            }
        if (v.empty()) {
            p.mType = 0;
            continue;
        }
        for (std::size_t k = 0; k < 3; ++k)
            p.mCentroid[k] /= static_cast<double>(v.size());
        const double ext[3] = {hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
        p.mExtent = std::sqrt(pc_dot(ext, ext));
    }
    const double min_fraction = FitPrimitivesOptions().mMinInlierFraction;
    const auto keep = [&](PcPatch& rPatch, const double* pParams) {
        PcStats s;
        if (!pc_accept(rPatch, pParams, FaceTolerance, min_fraction, s))
            return false;
        std::copy(pParams, pParams + kPrimitiveParams, rPatch.mParams);
        rPatch.mStats = s;
        return true;
    };

    // 1. Coaxial groups, seeded in row order by the axial primitives.
    std::vector<std::uint8_t> grouped(m, 0);
    std::vector<std::vector<double>> references;  // unit directions, axes first
    for (std::size_t a = 0; a < m; ++a) {
        if (grouped[a] || !patches[a].Axial())
            continue;
        grouped[a] = 1;
        double x[3];
        pc_axis_point(patches[a], x);
        const double* d = patches[a].mParams + 3;
        std::vector<std::size_t> rows = {a};
        // Every other row not yet grouped, in order: a sphere before the
        // seed too (an earlier axial row seeded its own group already).
        for (std::size_t b = 0; b < m; ++b)
            if (b != a && !grouped[b] && pc_on_axis(patches[b], x, d))
                rows.push_back(b);
        // Fit jointly; a member that no longer passes leaves and the rest are
        // fitted again, until every member passes or one is left.
        while (rows.size() >= 2) {
            std::vector<std::vector<double>> fitted;
            pc_fit_group(patches, rows, fitted);
            std::vector<std::size_t> passing;
            for (std::size_t i = 0; i < rows.size(); ++i) {
                PcStats s;
                if (!fitted[i].empty() &&
                    pc_accept(patches[rows[i]], fitted[i].data(), FaceTolerance, min_fraction, s))
                    passing.push_back(i);
            }
            if (passing.size() == rows.size()) {
                for (std::size_t i = 0; i < rows.size(); ++i) {
                    keep(patches[rows[i]], fitted[i].data());
                    grouped[rows[i]] = 1;
                }
                ++report.mNumGroups;
                break;
            }
            // Only the seed's own fit can leave the group unformed.
            if (passing.empty() || passing.front() != 0)
                break;
            std::vector<std::size_t> next;
            for (std::size_t i : passing)
                next.push_back(rows[i]);
            rows = std::move(next);
        }
        double axis[3] = {patches[a].mParams[3], patches[a].mParams[4], patches[a].mParams[5]};
        pc_canonical(axis);
        references.push_back({axis[0], axis[1], axis[2]});
    }

    // 2. Planes against the references: the axes, then the planes seen.
    std::vector<std::uint8_t> family(references.size(), 0);
    for (std::size_t r = 0; r < m; ++r) {
        PcPatch& p = patches[r];
        if (p.mType != 1)
            continue;
        const double budget = p.mExtent > 0.0 ? kPcReach * p.mTol / p.mExtent : 0.0;
        const double* n = p.mParams + 3;
        double target[3];
        bool have = false;
        std::size_t parallel = references.size();
        for (std::size_t i = 0; i < references.size() && !have; ++i) {
            const double* ref = references[i].data();
            if (pc_sin(n, ref) <= budget) {
                const double s = pc_dot(n, ref) < 0.0 ? -1.0 : 1.0;
                for (int k = 0; k < 3; ++k)
                    target[k] = s * ref[k];
                have = true;
                parallel = i;
            }
        }
        if (!have) {
            std::vector<std::size_t> perp;
            for (std::size_t i = 0; i < references.size(); ++i)
                if (std::fabs(pc_dot(n, references[i].data())) <= budget)
                    perp.push_back(i);
            if (!perp.empty()) {
                const double* r1 = references[perp[0]].data();
                const double* r2 = nullptr;
                for (std::size_t j = 1; j < perp.size() && !r2; ++j)
                    if (pc_sin(r1, references[perp[j]].data()) > 1e-6)
                        r2 = references[perp[j]].data();
                if (r2) {
                    pc_cross(r1, r2, target);
                } else {
                    const double h = pc_dot(n, r1);
                    for (int k = 0; k < 3; ++k)
                        target[k] = n[k] - h * r1[k];
                }
                have = pc_normalize(target);
                if (have && pc_dot(target, n) < 0.0)
                    for (int k = 0; k < 3; ++k)
                        target[k] = -target[k];
            }
        }
        bool snapped = false;
        if (have && (target[0] != n[0] || target[1] != n[1] || target[2] != n[2])) {
            // The offset refitted: the mean height of the points along it.
            double off = 0.0;
            for (std::size_t i = 0; i < p.NumPts(); ++i)
                off += pc_dot(target, &p.mPts[i * 3]);
            off /= static_cast<double>(p.NumPts());
            const double h = off - pc_dot(target, p.mCentroid);
            double cand[kPrimitiveParams] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            for (int k = 0; k < 3; ++k) {
                cand[k] = p.mCentroid[k] + h * target[k];
                cand[3 + k] = target[k];
            }
            snapped = keep(p, cand);
        }
        if (snapped && parallel < references.size())
            family[parallel] = 1;
        else if (parallel == references.size()) {
            references.push_back({p.mParams[3], p.mParams[4], p.mParams[5]});
            family.push_back(0);
        }
    }
    for (std::uint8_t f : family)
        report.mNumGroups += f;

    // 3. Equal radii: cylinders and spheres by their radius, tori by their
    // minor radius, each against the first of its kind in row order.
    const auto radius_slot = [](std::int64_t Type) { return Type == 5 ? 7 : 6; };
    std::vector<std::uint8_t> sized(m, 0);
    for (std::size_t a = 0; a < m; ++a) {
        const std::int64_t t = patches[a].mType;
        if (sized[a] || (t != 2 && t != 3 && t != 5))
            continue;
        sized[a] = 1;
        const int slot = radius_slot(t);
        std::vector<std::size_t> rows = {a};
        for (std::size_t b = a + 1; b < m; ++b)
            if (!sized[b] && patches[b].mType == t &&
                std::fabs(patches[b].mParams[slot] - patches[a].mParams[slot]) <=
                    std::min(patches[a].mTol, patches[b].mTol))
                rows.push_back(b);
        while (rows.size() >= 2) {
            double sum = 0.0, w = 0.0;
            for (std::size_t r : rows) {
                const double wi = static_cast<double>(patches[r].NumPts());
                sum += wi * patches[r].mParams[slot];
                w += wi;
            }
            const double radius = sum / w;
            std::vector<std::size_t> passing;
            for (std::size_t r : rows) {
                double cand[kPrimitiveParams];
                std::copy(patches[r].mParams, patches[r].mParams + kPrimitiveParams, cand);
                cand[slot] = radius;
                PcStats s;
                if (pc_accept(patches[r], cand, FaceTolerance, min_fraction, s))
                    passing.push_back(r);
            }
            if (passing.size() == rows.size()) {
                bool changed = false;
                for (std::size_t r : rows) {
                    double cand[kPrimitiveParams];
                    std::copy(patches[r].mParams, patches[r].mParams + kPrimitiveParams, cand);
                    changed = changed || cand[slot] != radius;
                    cand[slot] = radius;
                    keep(patches[r], cand);
                    sized[r] = 1;
                }
                if (changed)
                    ++report.mNumGroups;
                break;
            }
            rows = std::move(passing);
        }
    }

    // Write back the rows that changed.
    std::vector<std::uint8_t> changed(m, 0);
    for (std::size_t r = 0; r < m; ++r)
        for (int k = 0; k < kPrimitiveParams && patches[r].mType != 0; ++k)
            if (patches[r].mParams[k] != patches[r].mOriginal[k])
                changed[r] = 1;
    for (std::uint8_t c : changed)
        report.mNumConstrained += c;
    if (report.mNumConstrained == 0)
        return report;

    NDArray new_params = NDArray::Uninit(
        DType::Float64, {m, static_cast<std::size_t>(kPrimitiveParams)});
    NDArray new_rms = NDArray::Uninit(DType::Float64, {m});
    NDArray new_max = NDArray::Uninit(DType::Float64, {m});
    NDArray new_frac = NDArray::Uninit(DType::Float64, {m});
    for (std::size_t r = 0; r < m; ++r) {
        for (int k = 0; k < kPrimitiveParams; ++k) {
            const std::size_t i = r * kPrimitiveParams + static_cast<std::size_t>(k);
            new_params.As<double>()[i] =
                changed[r] ? patches[r].mParams[k] : read_double(params, i);
        }
        new_rms.As<double>()[r] = changed[r] ? patches[r].mStats.mRms : read_double(rms, r);
        new_max.As<double>()[r] = changed[r] ? patches[r].mStats.mMax : read_double(maxs, r);
        new_frac.As<double>()[r] =
            changed[r] ? patches[r].mStats.mFraction : read_double(fracs, r);
    }
    if (rMesh.HasCellData(kFitRmsName) && rMesh.HasCellData(kFitMaxName)) {
        std::vector<NDArray> rms_blocks, max_blocks;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const NDArray& old_rms = rMesh.CellData(kFitRmsName, b);
            const NDArray& old_max = rMesh.CellData(kFitMaxName, b);
            const std::size_t nc = old_rms.Size();
            const std::size_t base = static_cast<std::size_t>(bases[b]);
            NDArray br = NDArray::Uninit(DType::Float64, {nc});
            NDArray bm = NDArray::Uninit(DType::Float64, {nc});
            for (std::size_t c = 0; c < nc; ++c) {
                const std::int64_t row = base + c < ncells ? cell_row[base + c] : -1;
                const bool ch = row >= 0 && changed[static_cast<std::size_t>(row)];
                br.As<double>()[c] =
                    ch ? patches[static_cast<std::size_t>(row)].mStats.mRms : read_double(old_rms, c);
                bm.As<double>()[c] = ch ? patches[static_cast<std::size_t>(row)].mStats.mMax
                                        : (c < old_max.Size() ? read_double(old_max, c) : 0.0);
            }
            rms_blocks.push_back(std::move(br));
            max_blocks.push_back(std::move(bm));
        }
        rMesh.AddCellData(kFitRmsName, std::move(rms_blocks));
        rMesh.AddCellData(kFitMaxName, std::move(max_blocks));
    }
    rMesh.AddFieldData("fit:patch_params", std::move(new_params));
    rMesh.AddFieldData("fit:patch_rms", std::move(new_rms));
    rMesh.AddFieldData("fit:patch_max", std::move(new_max));
    rMesh.AddFieldData("fit:patch_fraction", std::move(new_frac));
    return report;
}

}  // namespace detail
}  // namespace meshioplusplus
