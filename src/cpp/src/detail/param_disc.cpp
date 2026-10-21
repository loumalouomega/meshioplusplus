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
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

// Project includes (private, not installed)
#include "param_disc.hpp"
#include "sparse_solve.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kPdPi = 3.14159265358979323846;

PatchVec pd_sub(const PatchVec& rA, const PatchVec& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}
double pd_dot(const PatchVec& rA, const PatchVec& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}
PatchVec pd_cross(const PatchVec& rA, const PatchVec& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}
double pd_norm(const PatchVec& rA) {
    return std::sqrt(pd_dot(rA, rA));
}

// The vertices the triangles use, ascending.
std::vector<std::size_t> pd_used(const std::vector<PatchTri>& rTris) {
    std::vector<std::size_t> v;
    for (const PatchTri& t : rTris)
        v.insert(v.end(), t.begin(), t.end());
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

}  // namespace

std::vector<std::vector<std::size_t>> param_boundary_loops(const std::vector<PatchTri>& rTris,
                                                           bool& rOk) {
    rOk = true;
    // Directed edges; a boundary edge is one whose reverse no triangle uses.
    std::map<std::pair<std::size_t, std::size_t>, int> directed;
    for (const PatchTri& t : rTris)
        for (int k = 0; k < 3; ++k)
            ++directed[{t[static_cast<std::size_t>(k)], t[static_cast<std::size_t>((k + 1) % 3)]}];
    std::map<std::size_t, std::size_t> next;
    for (const auto& kv : directed) {
        if (kv.second > 1)
            rOk = false;  // an edge used twice the same way: inconsistent or non-manifold
        if (directed.count({kv.first.second, kv.first.first}))
            continue;
        if (next.count(kv.first.first))
            rOk = false;  // a pinched boundary vertex
        next[kv.first.first] = kv.first.second;
    }
    std::vector<std::vector<std::size_t>> loops;
    std::set<std::size_t> done;
    for (const auto& kv : next) {
        if (done.count(kv.first))
            continue;
        std::vector<std::size_t> loop;
        std::size_t v = kv.first;
        while (!done.count(v)) {
            done.insert(v);
            loop.push_back(v);
            const auto it = next.find(v);
            if (it == next.end()) {
                rOk = false;
                break;
            }
            v = it->second;
        }
        if (v != kv.first)
            rOk = false;
        loops.push_back(std::move(loop));
    }
    return loops;
}

bool param_is_disc(const std::vector<PatchTri>& rTris) {
    if (rTris.empty())
        return false;
    bool ok = true;
    const std::vector<std::vector<std::size_t>> loops = param_boundary_loops(rTris, ok);
    if (!ok || loops.size() != 1)
        return false;
    std::set<std::pair<std::size_t, std::size_t>> edges;
    for (const PatchTri& t : rTris)
        for (int k = 0; k < 3; ++k) {
            std::size_t a = t[static_cast<std::size_t>(k)],
                        b = t[static_cast<std::size_t>((k + 1) % 3)];
            if (a > b)
                std::swap(a, b);
            edges.insert({a, b});
        }
    const long long v = static_cast<long long>(pd_used(rTris).size());
    const long long e = static_cast<long long>(edges.size());
    const long long f = static_cast<long long>(rTris.size());
    return v - e + f == 1;
}

std::size_t param_flips(const std::vector<PatchTri>& rTris, const std::vector<double>& rU,
                        const std::vector<double>& rV) {
    std::size_t pos = 0, neg = 0, zero = 0;
    for (const PatchTri& t : rTris) {
        const double a = (rU[t[1]] - rU[t[0]]) * (rV[t[2]] - rV[t[0]]) -
                         (rU[t[2]] - rU[t[0]]) * (rV[t[1]] - rV[t[0]]);
        if (a > 0.0)
            ++pos;
        else if (a < 0.0)
            ++neg;
        else
            ++zero;
    }
    return std::min(pos, neg) + zero;
}

DiscParam param_lscm(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris) {
    DiscParam out;
    bool ok = true;
    const std::vector<std::vector<std::size_t>> loops = param_boundary_loops(rTris, ok);
    if (!ok || loops.empty())
        return out;
    const std::size_t n = rPts.size();
    // Pin the two boundary vertices farthest apart (ties to the lowest ids).
    const std::vector<std::size_t>& loop = loops.front();
    std::size_t pa = loop.front(), pb = loop.front();
    double far = -1.0;
    for (std::size_t i = 0; i < loop.size(); ++i)
        for (std::size_t j = i + 1; j < loop.size(); ++j) {
            const double d = pd_norm(pd_sub(rPts[loop[i]], rPts[loop[j]]));
            if (d > far) {
                far = d;
                pa = loop[i];
                pb = loop[j];
            }
        }
    if (!(far > 0.0))
        return out;
    // Unknowns: u then v of every vertex but the pinned ones.
    std::vector<std::ptrdiff_t> col(n, -1);
    std::size_t nfree = 0;
    for (std::size_t i : pd_used(rTris))
        if (i != pa && i != pb)
            col[i] = static_cast<std::ptrdiff_t>(nfree++);
    const double pin_u[2] = {0.0, far}, pin_v[2] = {0.0, 0.0};
    const std::size_t m = 2 * nfree;
    std::vector<SparseTriplet> ata;
    std::vector<double> atb(m, 0.0);
    for (const PatchTri& t : rTris) {
        const PatchVec e1 = pd_sub(rPts[t[1]], rPts[t[0]]);
        const PatchVec e2 = pd_sub(rPts[t[2]], rPts[t[0]]);
        const PatchVec nrm = pd_cross(e1, e2);
        const double twice = pd_norm(nrm);
        const double l1 = pd_norm(e1);
        if (!(twice > 0.0) || !(l1 > 0.0))
            continue;
        // The triangle in its own plane, counter-clockwise.
        const PatchVec ex = {e1[0] / l1, e1[1] / l1, e1[2] / l1};
        const PatchVec nz = {nrm[0] / twice, nrm[1] / twice, nrm[2] / twice};
        const PatchVec ey = pd_cross(nz, ex);
        const double x[3] = {0.0, l1, pd_dot(e2, ex)};
        const double y[3] = {0.0, 0.0, pd_dot(e2, ey)};
        // The gradient of each hat function: the opposite edge turned 90 degrees.
        double gx[3], gy[3];
        for (int j = 0; j < 3; ++j) {
            const int k = (j + 1) % 3, l = (j + 2) % 3;
            gx[j] = -(y[l] - y[k]) / twice;
            gy[j] = (x[l] - x[k]) / twice;
        }
        // Residual R90 grad u - grad v, two rows, weight sqrt(area).
        const double w = std::sqrt(0.5 * twice);
        // Row 0: -(grad u)_y - (grad v)_x; row 1: (grad u)_x - (grad v)_y.
        for (int row = 0; row < 2; ++row) {
            double cu[3], cv[3];
            for (int j = 0; j < 3; ++j) {
                cu[j] = w * (row == 0 ? -gy[j] : gx[j]);
                cv[j] = w * (row == 0 ? -gx[j] : -gy[j]);
            }
            // The pinned vertices' contribution moves to the right-hand side.
            double rhs = 0.0;
            for (int j = 0; j < 3; ++j) {
                const std::size_t vtx = t[static_cast<std::size_t>(j)];
                if (vtx == pa || vtx == pb) {
                    const int p = vtx == pa ? 0 : 1;
                    rhs -= cu[j] * pin_u[p] + cv[j] * pin_v[p];
                }
            }
            for (int a = 0; a < 3; ++a) {
                const std::ptrdiff_t ca = col[t[static_cast<std::size_t>(a)]];
                if (ca < 0)
                    continue;
                const std::size_t iu = static_cast<std::size_t>(ca);
                const std::size_t iv = nfree + iu;
                atb[iu] += cu[a] * rhs;
                atb[iv] += cv[a] * rhs;
                for (int b = 0; b < 3; ++b) {
                    const std::ptrdiff_t cb = col[t[static_cast<std::size_t>(b)]];
                    if (cb < 0)
                        continue;
                    const std::size_t ju = static_cast<std::size_t>(cb);
                    const std::size_t jv = nfree + ju;
                    ata.push_back({iu, ju, cu[a] * cu[b]});
                    ata.push_back({iu, jv, cu[a] * cv[b]});
                    ata.push_back({iv, ju, cv[a] * cu[b]});
                    ata.push_back({iv, jv, cv[a] * cv[b]});
                }
            }
        }
    }
    const SparseCsr a = csr_from_triplets(m, m, std::move(ata));
    std::vector<double> x(m, 0.0);
    const SparseSolveReport rep = csr_pcg(a, atb.data(), x.data(), 1e-12, 20 * m + 100);
    if (!(rep.mRelativeResidual <= 1e-8))
        return out;
    out.mU.assign(n, 0.0);
    out.mV.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        if (col[i] >= 0) {
            out.mU[i] = x[static_cast<std::size_t>(col[i])];
            out.mV[i] = x[nfree + static_cast<std::size_t>(col[i])];
        }
    out.mU[pa] = pin_u[0];
    out.mV[pa] = pin_v[0];
    out.mU[pb] = pin_u[1];
    out.mV[pb] = pin_v[1];
    // Keep the 3-D orientation: the majority of triangles counter-clockwise.
    std::size_t pos = 0, neg = 0;
    for (const PatchTri& t : rTris) {
        const double s = (out.mU[t[1]] - out.mU[t[0]]) * (out.mV[t[2]] - out.mV[t[0]]) -
                         (out.mU[t[2]] - out.mU[t[0]]) * (out.mV[t[1]] - out.mV[t[0]]);
        (s >= 0.0 ? pos : neg) += 1;
    }
    if (neg > pos)
        for (double& v : out.mV)
            v = -v;
    out.mNumFlips = param_flips(rTris, out.mU, out.mV);
    out.mOk = true;
    return out;
}

DiscParam param_mean_value(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris) {
    DiscParam out;
    bool ok = true;
    const std::vector<std::vector<std::size_t>> loops = param_boundary_loops(rTris, ok);
    if (!ok || loops.size() != 1)
        return out;
    const std::size_t n = rPts.size();
    const std::vector<std::size_t>& loop = loops.front();
    // The boundary on a circle by arc length.
    std::vector<double> cum(loop.size() + 1, 0.0);
    for (std::size_t i = 0; i < loop.size(); ++i)
        cum[i + 1] = cum[i] + pd_norm(pd_sub(rPts[loop[(i + 1) % loop.size()]], rPts[loop[i]]));
    const double len = cum.back();
    if (!(len > 0.0))
        return out;
    const double radius = len / (2.0 * kPdPi);
    out.mU.assign(n, 0.0);
    out.mV.assign(n, 0.0);
    std::vector<std::uint8_t> on_boundary(n, 0);
    for (std::size_t i = 0; i < loop.size(); ++i) {
        const double a = 2.0 * kPdPi * cum[i] / len;
        out.mU[loop[i]] = radius * std::cos(a);
        out.mV[loop[i]] = radius * std::sin(a);
        on_boundary[loop[i]] = 1;
    }
    std::vector<std::ptrdiff_t> col(n, -1);
    std::size_t nfree = 0;
    for (std::size_t i : pd_used(rTris))
        if (!on_boundary[i])
            col[i] = static_cast<std::ptrdiff_t>(nfree++);
    if (nfree > 0) {
        // Floater's weights: w_ij = (tan(a/2) + tan(b/2)) / |p_i - p_j|, the
        // angles at p_i in the two triangles beside the edge.
        std::map<std::pair<std::size_t, std::size_t>, double> w;
        for (const PatchTri& t : rTris)
            for (int c = 0; c < 3; ++c) {
                const std::size_t i = t[static_cast<std::size_t>(c)];
                const std::size_t j = t[static_cast<std::size_t>((c + 1) % 3)];
                const std::size_t k = t[static_cast<std::size_t>((c + 2) % 3)];
                const PatchVec a = pd_sub(rPts[j], rPts[i]);
                const PatchVec b = pd_sub(rPts[k], rPts[i]);
                const double la = pd_norm(a), lb = pd_norm(b);
                if (!(la > 0.0) || !(lb > 0.0))
                    continue;
                const double cosv = std::min(1.0, std::max(-1.0, pd_dot(a, b) / (la * lb)));
                const double half = std::tan(0.5 * std::acos(cosv));
                w[{i, j}] += half / la;
                w[{i, k}] += half / lb;
            }
        std::vector<SparseTriplet> trip;
        std::vector<double> bu(nfree, 0.0), bv(nfree, 0.0), diag(nfree, 0.0);
        for (const auto& kv : w) {
            const std::size_t i = kv.first.first, j = kv.first.second;
            if (col[i] < 0)
                continue;
            const std::size_t r = static_cast<std::size_t>(col[i]);
            diag[r] += kv.second;
            if (col[j] >= 0) {
                trip.push_back({r, static_cast<std::size_t>(col[j]), -kv.second});
            } else {
                bu[r] += kv.second * out.mU[j];
                bv[r] += kv.second * out.mV[j];
            }
        }
        for (std::size_t r = 0; r < nfree; ++r)
            trip.push_back({r, r, diag[r]});
        const SparseCsr a = csr_from_triplets(nfree, nfree, std::move(trip));
        std::vector<double> xu(nfree, 0.0), xv(nfree, 0.0);
        const SparseSolveReport ru = csr_bicgstab(a, bu.data(), xu.data(), 1e-12, 20 * nfree + 100);
        const SparseSolveReport rv = csr_bicgstab(a, bv.data(), xv.data(), 1e-12, 20 * nfree + 100);
        if (!(ru.mRelativeResidual <= 1e-8) || !(rv.mRelativeResidual <= 1e-8))
            return out;
        for (std::size_t i = 0; i < n; ++i)
            if (col[i] >= 0) {
                out.mU[i] = xu[static_cast<std::size_t>(col[i])];
                out.mV[i] = xv[static_cast<std::size_t>(col[i])];
            }
    }
    // The loop runs the way the triangles do, so they come out counter-
    // clockwise when the loop is: flip v otherwise.
    std::size_t pos = 0, neg = 0;
    for (const PatchTri& t : rTris) {
        const double s = (out.mU[t[1]] - out.mU[t[0]]) * (out.mV[t[2]] - out.mV[t[0]]) -
                         (out.mU[t[2]] - out.mU[t[0]]) * (out.mV[t[1]] - out.mV[t[0]]);
        (s >= 0.0 ? pos : neg) += 1;
    }
    if (neg > pos)
        for (double& v : out.mV)
            v = -v;
    out.mNumFlips = param_flips(rTris, out.mU, out.mV);
    out.mOk = true;
    return out;
}

namespace {

// Turn a parametrisation (a conformal map is one up to a rotation) so its
// bounding box -- the spline's domain -- is the smallest: the boundary's box
// over rotations in half-degree steps, the first smallest kept.
void pd_align(const std::vector<PatchTri>& rTris, DiscParam& rP) {
    bool ok = true;
    const std::vector<std::vector<std::size_t>> loops = param_boundary_loops(rTris, ok);
    if (!ok || loops.empty())
        return;
    const std::vector<std::size_t>& loop = loops.front();
    double best_area = std::numeric_limits<double>::infinity();
    double best = 0.0;
    for (int k = 0; k < 180; ++k) {
        const double a = static_cast<double>(k) * kPdPi / 360.0;
        const double c = std::cos(a), s = std::sin(a);
        double lo[2] = {std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
        double hi[2] = {-lo[0], -lo[1]};
        for (std::size_t i : loop) {
            const double x = c * rP.mU[i] + s * rP.mV[i];
            const double y = -s * rP.mU[i] + c * rP.mV[i];
            lo[0] = std::min(lo[0], x);
            hi[0] = std::max(hi[0], x);
            lo[1] = std::min(lo[1], y);
            hi[1] = std::max(hi[1], y);
        }
        const double area = (hi[0] - lo[0]) * (hi[1] - lo[1]);
        if (area < best_area * (1.0 - 1e-12)) {
            best_area = area;
            best = a;
        }
    }
    const double c = std::cos(best), s = std::sin(best);
    for (std::size_t i = 0; i < rP.mU.size(); ++i) {
        const double u = rP.mU[i], v = rP.mV[i];
        rP.mU[i] = c * u + s * v;
        rP.mV[i] = -s * u + c * v;
    }
}

}  // namespace

DiscParam param_disc(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris) {
    if (!param_is_disc(rTris))
        return DiscParam{};
    DiscParam p = param_lscm(rPts, rTris);
    if (p.mOk && p.mNumFlips == 0) {
        pd_align(rTris, p);
        return p;
    }
    p = param_mean_value(rPts, rTris);
    if (p.mOk && p.mNumFlips == 0)
        return p;
    p.mOk = false;
    return p;
}

}  // namespace detail
}  // namespace meshioplusplus
