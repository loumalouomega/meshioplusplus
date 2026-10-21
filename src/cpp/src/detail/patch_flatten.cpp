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
#include <map>
#include <numeric>
#include <utility>

// Project includes
#include "patch_flatten.hpp"
#include "lsq_fit.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kPfPi = 3.14159265358979323846;

double pf_dot(const PatchVec& rA, const PatchVec& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

PatchVec pf_cross(const PatchVec& rA, const PatchVec& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}

double pf_norm(const PatchVec& rA) {
    return std::sqrt(pf_dot(rA, rA));
}

PatchVec pf_sub(const PatchVec& rA, const PatchVec& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}

PatchVec pf_add(const PatchVec& rA, const PatchVec& rB) {
    return {rA[0] + rB[0], rA[1] + rB[1], rA[2] + rB[2]};
}

PatchVec pf_scale(const PatchVec& rA, double S) {
    return {rA[0] * S, rA[1] * S, rA[2] * S};
}

void pf_normalize(PatchVec& rA) {
    const double d = pf_norm(rA);
    rA = {rA[0] / d, rA[1] / d, rA[2] / d};
}

// Two unit vectors completing @p rA to a right-handed frame.
void pf_frame(const PatchVec& rA, PatchVec& rE1, PatchVec& rE2) {
    const PatchVec seed =
        std::fabs(rA[0]) < 0.9 ? PatchVec{1.0, 0.0, 0.0} : PatchVec{0.0, 1.0, 0.0};
    rE1 = pf_sub(seed, pf_scale(rA, pf_dot(seed, rA)));
    pf_normalize(rE1);
    rE2 = pf_cross(rA, rE1);
}

// Local points and triangles of a subset of a patch's triangles.
void pf_local(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
              const std::vector<std::size_t>& rSubset, std::vector<PatchVec>& rLocalPts,
              std::vector<PatchTri>& rLocalTris) {
    std::map<std::size_t, std::size_t> local;
    for (std::size_t t : rSubset)
        for (std::size_t v : rTris[t])
            local.emplace(v, 0);
    rLocalPts.clear();
    for (auto& kv : local) {
        kv.second = rLocalPts.size();
        rLocalPts.push_back(rPts[kv.first]);
    }
    rLocalTris.clear();
    for (std::size_t t : rSubset)
        rLocalTris.push_back({local.at(rTris[t][0]), local.at(rTris[t][1]), local.at(rTris[t][2])});
}

// The edge-connected components of a subset (each ascending), in the order of
// their lowest triangle.
std::vector<std::vector<std::size_t>> pf_components(const std::vector<PatchTri>& rTris,
                                                    const std::vector<std::size_t>& rSubset) {
    std::vector<std::size_t> parent(rSubset.size());
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    auto find = [&](std::size_t X) {
        while (parent[X] != X) {
            parent[X] = parent[parent[X]];
            X = parent[X];
        }
        return X;
    };
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> first;
    for (std::size_t i = 0; i < rSubset.size(); ++i) {
        const PatchTri& t = rTris[rSubset[i]];
        for (int k = 0; k < 3; ++k) {
            const std::size_t a = t[static_cast<std::size_t>(k)];
            const std::size_t b = t[static_cast<std::size_t>((k + 1) % 3)];
            auto ins = first.emplace(std::make_pair(std::min(a, b), std::max(a, b)), i);
            if (!ins.second) {
                const std::size_t ra = find(i), rb = find(ins.first->second);
                if (ra != rb)
                    parent[std::max(ra, rb)] = std::min(ra, rb);
            }
        }
    }
    std::map<std::size_t, std::vector<std::size_t>> by_root;
    for (std::size_t i = 0; i < rSubset.size(); ++i)
        by_root[find(i)].push_back(rSubset[i]);
    std::vector<std::vector<std::size_t>> out;
    for (auto& kv : by_root)
        out.push_back(std::move(kv.second));
    std::sort(out.begin(), out.end(),
              [](const auto& rA, const auto& rB) { return rA.front() < rB.front(); });
    return out;
}

void pf_split(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
              const std::vector<PatchVec>& rNormals, const std::vector<double>& rAreas,
              const std::vector<std::size_t>& rSubset, int Depth, const PatchSplitOptions& rOptions,
              std::vector<std::vector<std::size_t>>& rPieces) {
    if (rSubset.size() < rOptions.mMinTriangles)
        return;
    std::vector<PatchVec> lp;
    std::vector<PatchTri> lt;
    pf_local(rPts, rTris, rSubset, lp, lt);
    if (patch_flattenable(lp, lt)) {
        rPieces.push_back(rSubset);
        return;
    }
    if (Depth >= rOptions.mMaxDepth)
        return;
    // The direction the normals spread along most.
    double m[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, w[3], v[9];
    double total = 0.0;
    for (std::size_t t : rSubset) {
        const PatchVec& n = rNormals[t];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m[i * 3 + j] +=
                    rAreas[t] * n[static_cast<std::size_t>(i)] * n[static_cast<std::size_t>(j)];
        total += rAreas[t];
    }
    lsq_sym3_eigen(m, w, v);
    const PatchVec e{v[6], v[7], v[8]};
    std::vector<std::pair<double, std::size_t>> keyed;
    keyed.reserve(rSubset.size());
    for (std::size_t t : rSubset)
        keyed.emplace_back(pf_dot(rNormals[t], e), t);
    std::sort(keyed.begin(), keyed.end());
    // The area-weighted median; both halves keep at least one triangle.
    std::size_t k = 0;
    double acc = 0.0;
    while (k + 1 < keyed.size() && acc + rAreas[keyed[k].second] <= 0.5 * total) {
        acc += rAreas[keyed[k].second];
        ++k;
    }
    k = std::max<std::size_t>(k, 1);
    std::vector<std::size_t> lo, hi;
    for (std::size_t i = 0; i < keyed.size(); ++i)
        (i < k ? lo : hi).push_back(keyed[i].second);
    for (std::vector<std::size_t>* half : {&lo, &hi}) {
        std::sort(half->begin(), half->end());
        for (const auto& comp : pf_components(rTris, *half))
            pf_split(rPts, rTris, rNormals, rAreas, comp, Depth + 1, rOptions, rPieces);
    }
}

}  // namespace

PatchVec PatchFlat::At(double S, double T, double H) const {
    if (!mCylinder)
        return pf_add(pf_add(pf_add(mOrigin, pf_scale(mE1, S)), pf_scale(mE2, T)),
                      pf_scale(mAxis, H));
    const double th = mTheta0 + S / mRadius;
    return pf_add(
        pf_add(mOrigin, pf_scale(mAxis, T)),
        pf_scale(pf_add(pf_scale(mE1, std::cos(th)), pf_scale(mE2, std::sin(th))), mRadius + H));
}

bool patch_flatten(int Kind, const std::vector<PatchVec>& rPts,
                   const std::vector<PatchVec>& rNormals, const std::vector<double>& rAreas,
                   PatchFlat& rOut) {
    const std::size_t n = rPts.size();
    rOut = PatchFlat{};
    if (n == 0)
        return false;
    PatchVec c{0.0, 0.0, 0.0};
    for (const PatchVec& p : rPts)
        c = pf_add(c, p);
    c = {c[0] / static_cast<double>(n), c[1] / static_cast<double>(n),
         c[2] / static_cast<double>(n)};
    double m[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, w[3], v[9];
    if (Kind == 0) {
        for (const PatchVec& p : rPts) {
            const PatchVec d = pf_sub(p, c);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    m[i * 3 + j] += d[static_cast<std::size_t>(i)] * d[static_cast<std::size_t>(j)];
        }
        lsq_sym3_eigen(m, w, v);
        rOut.mOrigin = c;
        rOut.mAxis = {v[0], v[1], v[2]};  // the smallest spread: the normal
        rOut.mE1 = {v[6], v[7], v[8]};
        rOut.mE2 = pf_cross(rOut.mAxis, rOut.mE1);
        for (const PatchVec& p : rPts) {
            const PatchVec d = pf_sub(p, c);
            rOut.mS.push_back(pf_dot(d, rOut.mE1));
            rOut.mT.push_back(pf_dot(d, rOut.mE2));
            rOut.mH.push_back(pf_dot(d, rOut.mAxis));
        }
        return true;
    }
    for (std::size_t f = 0; f < rNormals.size(); ++f)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m[i * 3 + j] += rAreas[f] * rNormals[f][static_cast<std::size_t>(i)] *
                                rNormals[f][static_cast<std::size_t>(j)];
    lsq_sym3_eigen(m, w, v);
    PatchVec axis{v[0], v[1], v[2]};
    if (pf_norm(axis) < 0.5)
        return false;
    pf_normalize(axis);
    PatchVec e1, e2;
    pf_frame(axis, e1, e2);
    // The circle the points' projections fit (Kasa): x^2 + y^2 + D x + E y + F = 0.
    double a[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, b[3] = {0, 0, 0};
    for (const PatchVec& p : rPts) {
        const PatchVec d = pf_sub(p, c);
        const double x = pf_dot(d, e1), y = pf_dot(d, e2), r[3] = {x, y, 1.0};
        const double rhs = -(x * x + y * y);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j)
                a[i * 3 + j] += r[i] * r[j];
            b[i] += r[i] * rhs;
        }
    }
    // Solve the 3x3 normal equations by Cramer's rule.
    auto det3 = [](const double* q) {
        return q[0] * (q[4] * q[8] - q[5] * q[7]) - q[1] * (q[3] * q[8] - q[5] * q[6]) +
               q[2] * (q[3] * q[7] - q[4] * q[6]);
    };
    const double d0 = det3(a);
    if (!(std::fabs(d0) > 0.0))
        return false;
    double sol[3];
    for (int k = 0; k < 3; ++k) {
        double q[9];
        std::copy(a, a + 9, q);
        for (int i = 0; i < 3; ++i)
            q[i * 3 + k] = b[i];
        sol[k] = det3(q) / d0;
    }
    const double cx = -sol[0] / 2.0, cy = -sol[1] / 2.0;
    const double r2 = cx * cx + cy * cy - sol[2];
    if (!(r2 > 0.0) || !std::isfinite(r2))
        return false;
    const double radius = std::sqrt(r2);
    const PatchVec origin = pf_add(pf_add(c, pf_scale(e1, cx)), pf_scale(e2, cy));
    std::vector<double> theta(n);
    for (std::size_t i = 0; i < n; ++i) {
        const PatchVec d = pf_sub(rPts[i], origin);
        theta[i] = std::atan2(pf_dot(d, e2), pf_dot(d, e1));
    }
    // Cut the circle in the widest gap between the points' angles.
    std::vector<double> sorted = theta;
    std::sort(sorted.begin(), sorted.end());
    double gap = sorted.front() + 2.0 * kPfPi - sorted.back(), cut = sorted.back();
    for (std::size_t i = 1; i < sorted.size(); ++i)
        if (sorted[i] - sorted[i - 1] > gap) {
            gap = sorted[i] - sorted[i - 1];
            cut = sorted[i - 1];
        }
    if (!(gap > 1e-6))
        return false;
    rOut.mCylinder = true;
    rOut.mOrigin = origin;
    rOut.mAxis = axis;
    rOut.mE1 = e1;
    rOut.mE2 = e2;
    rOut.mRadius = radius;
    rOut.mTheta0 = cut + gap / 2.0;  // the parameter starts just past the gap's middle
    for (std::size_t i = 0; i < n; ++i) {
        double th = theta[i] - rOut.mTheta0;
        th -= 2.0 * kPfPi * std::floor(th / (2.0 * kPfPi));
        const PatchVec d = pf_sub(rPts[i], origin);
        rOut.mS.push_back(radius * th);
        rOut.mT.push_back(pf_dot(d, axis));
        const PatchVec radial = pf_sub(d, pf_scale(axis, pf_dot(d, axis)));
        rOut.mH.push_back(pf_norm(radial) - radius);
    }
    return true;
}

bool patch_flat_folds(const PatchFlat& rFlat, const std::vector<PatchTri>& rTris) {
    const std::size_t nt = rTris.size();
    double max_area = 0.0, pos = 0.0, neg = 0.0;
    std::vector<double> signed_area(nt);
    for (std::size_t t = 0; t < nt; ++t) {
        const PatchTri& tr = rTris[t];
        const double ax = rFlat.mS[tr[1]] - rFlat.mS[tr[0]], ay = rFlat.mT[tr[1]] - rFlat.mT[tr[0]];
        const double bx = rFlat.mS[tr[2]] - rFlat.mS[tr[0]], by = rFlat.mT[tr[2]] - rFlat.mT[tr[0]];
        signed_area[t] = 0.5 * (ax * by - ay * bx);
        max_area = std::max(max_area, std::fabs(signed_area[t]));
        (signed_area[t] > 0.0 ? pos : neg) += std::fabs(signed_area[t]);
    }
    const double sign = pos >= neg ? 1.0 : -1.0;
    bool folds = !(max_area > 0.0);
    for (double sa : signed_area)
        folds = folds || sa * sign < -1e-9 * max_area;
    return folds;
}

void patch_normals(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                   std::vector<PatchVec>& rNormals, std::vector<double>& rAreas) {
    const std::size_t nt = rTris.size();
    rNormals.assign(nt, PatchVec{0.0, 0.0, 0.0});
    rAreas.assign(nt, 0.0);
    for (std::size_t t = 0; t < nt; ++t) {
        const PatchVec e1 = pf_sub(rPts[rTris[t][1]], rPts[rTris[t][0]]);
        const PatchVec e2 = pf_sub(rPts[rTris[t][2]], rPts[rTris[t][0]]);
        PatchVec nrm = pf_cross(e1, e2);
        rAreas[t] = 0.5 * pf_norm(nrm);
        if (pf_norm(nrm) > 0.0)
            pf_normalize(nrm);
        rNormals[t] = nrm;
    }
}

bool patch_flattenable(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris) {
    if (rPts.size() < 4 || rTris.size() < 2)
        return false;
    std::vector<PatchVec> normals;
    std::vector<double> areas;
    patch_normals(rPts, rTris, normals, areas);
    for (int kind = 0; kind < 2; ++kind) {
        PatchFlat flat;
        if (patch_flatten(kind, rPts, normals, areas, flat) && !patch_flat_folds(flat, rTris))
            return true;
    }
    return false;
}

std::vector<std::int64_t> patch_bisect(const std::vector<PatchVec>& rPts,
                                       const std::vector<PatchTri>& rTris,
                                       std::size_t& rNumPieces) {
    rNumPieces = 0;
    std::vector<std::int64_t> out(rTris.size(), -1);
    if (rTris.size() < 2)
        return out;
    std::vector<PatchVec> normals;
    std::vector<double> areas;
    patch_normals(rPts, rTris, normals, areas);
    double m[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, w[3], v[9];
    double total = 0.0;
    for (std::size_t t = 0; t < rTris.size(); ++t) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m[i * 3 + j] += areas[t] * normals[t][static_cast<std::size_t>(i)] *
                                normals[t][static_cast<std::size_t>(j)];
        total += areas[t];
    }
    lsq_sym3_eigen(m, w, v);
    const PatchVec e{v[6], v[7], v[8]};
    std::vector<std::pair<double, std::size_t>> keyed;
    for (std::size_t t = 0; t < rTris.size(); ++t)
        keyed.emplace_back(pf_dot(normals[t], e), t);
    std::sort(keyed.begin(), keyed.end());
    std::size_t k = 0;
    double acc = 0.0;
    while (k + 1 < keyed.size() && acc + areas[keyed[k].second] <= 0.5 * total) {
        acc += areas[keyed[k].second];
        ++k;
    }
    k = std::max<std::size_t>(k, 1);
    std::vector<std::size_t> lo, hi;
    for (std::size_t i = 0; i < keyed.size(); ++i)
        (i < k ? lo : hi).push_back(keyed[i].second);
    std::vector<std::vector<std::size_t>> pieces;
    for (std::vector<std::size_t>* half : {&lo, &hi}) {
        std::sort(half->begin(), half->end());
        for (auto& comp : pf_components(rTris, *half))
            pieces.push_back(std::move(comp));
    }
    std::sort(pieces.begin(), pieces.end(),
              [](const auto& rA, const auto& rB) { return rA.front() < rB.front(); });
    for (std::size_t p = 0; p < pieces.size(); ++p)
        for (std::size_t t : pieces[p])
            out[t] = static_cast<std::int64_t>(p);
    rNumPieces = pieces.size();
    return out;
}

std::vector<std::int64_t> patch_split_until_flat(const std::vector<PatchVec>& rPts,
                                                 const std::vector<PatchTri>& rTris,
                                                 const PatchSplitOptions& rOptions,
                                                 std::size_t& rNumPieces) {
    std::vector<PatchVec> normals;
    std::vector<double> areas;
    patch_normals(rPts, rTris, normals, areas);
    std::vector<std::size_t> all(rTris.size());
    std::iota(all.begin(), all.end(), std::size_t{0});
    std::vector<std::vector<std::size_t>> pieces;
    for (const auto& comp : pf_components(rTris, all))
        pf_split(rPts, rTris, normals, areas, comp, 0, rOptions, pieces);
    std::sort(pieces.begin(), pieces.end(),
              [](const auto& rA, const auto& rB) { return rA.front() < rB.front(); });
    std::vector<std::int64_t> out(rTris.size(), -1);
    for (std::size_t p = 0; p < pieces.size(); ++p)
        for (std::size_t t : pieces[p])
            out[t] = static_cast<std::int64_t>(p);
    rNumPieces = pieces.size();
    return out;
}

std::vector<std::int64_t> patch_bisect(const std::vector<PatchVec>& rPts,
                                       const std::vector<PatchTri>& rTris,
                                       std::size_t MinTriangles, std::size_t& rNumPieces) {
    rNumPieces = 0;
    const std::size_t nt = rTris.size();
    if (nt < 2 || nt < 2 * MinTriangles)
        return {};
    std::vector<PatchVec> normals, centroids(nt);
    std::vector<double> areas;
    patch_normals(rPts, rTris, normals, areas);
    PatchVec mean{0.0, 0.0, 0.0};
    double total = 0.0;
    for (std::size_t t = 0; t < nt; ++t) {
        for (int j = 0; j < 3; ++j) {
            double c = 0.0;
            for (std::size_t v : rTris[t])
                c += rPts[v][static_cast<std::size_t>(j)];
            centroids[t][static_cast<std::size_t>(j)] = c / 3.0;
            mean[static_cast<std::size_t>(j)] += areas[t] * c / 3.0;
        }
        total += areas[t];
    }
    if (!(total > 0.0))
        return {};
    for (double& x : mean)
        x /= total;
    // The direction the centroids spread along most.
    double m[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, w[3], v[9];
    for (std::size_t t = 0; t < nt; ++t) {
        const PatchVec d = pf_sub(centroids[t], mean);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m[i * 3 + j] += areas[t] * d[static_cast<std::size_t>(i)] * d[static_cast<std::size_t>(j)];
    }
    lsq_sym3_eigen(m, w, v);
    const PatchVec e{v[6], v[7], v[8]};
    std::vector<std::pair<double, std::size_t>> keyed;
    keyed.reserve(nt);
    for (std::size_t t = 0; t < nt; ++t)
        keyed.emplace_back(pf_dot(pf_sub(centroids[t], mean), e), t);
    std::sort(keyed.begin(), keyed.end());
    std::size_t k = 0;
    double acc = 0.0;
    while (k + 1 < keyed.size() && acc + areas[keyed[k].second] <= 0.5 * total) {
        acc += areas[keyed[k].second];
        ++k;
    }
    k = std::max<std::size_t>(k, 1);
    std::vector<std::size_t> lo, hi;
    for (std::size_t i = 0; i < keyed.size(); ++i)
        (i < k ? lo : hi).push_back(keyed[i].second);
    std::vector<std::vector<std::size_t>> pieces;
    for (std::vector<std::size_t>* half : {&lo, &hi}) {
        std::sort(half->begin(), half->end());
        for (auto& comp : pf_components(rTris, *half))
            if (comp.size() >= MinTriangles)
                pieces.push_back(std::move(comp));
    }
    std::sort(pieces.begin(), pieces.end(),
              [](const auto& rA, const auto& rB) { return rA.front() < rB.front(); });
    std::vector<std::int64_t> out(nt, -1);
    for (std::size_t p = 0; p < pieces.size(); ++p)
        for (std::size_t t : pieces[p])
            out[t] = static_cast<std::int64_t>(p);
    rNumPieces = pieces.size();
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
