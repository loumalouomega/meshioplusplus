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
//                   SPDX-License-Identifier: MIT
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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "parallel_workers.hpp"
#include "sparse_solve.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// Vector loops this short run serially: a parallel region costs more than
// they do, and far more on a loaded machine.
constexpr std::size_t kSpsGrain = std::size_t{1} << 15;

// y = a x + y, element-wise in parallel (independent entries).
void sps_axpy(double A, const double* pX, double* pY, std::size_t N) {
    parallel_for(N, [&](std::size_t i) { pY[i] += A * pX[i]; }, kSpsGrain);
}

// The relative residual ||b - A x|| / ||b|| (0 for a zero b), into pR.
double sps_residual(const SparseCsr& rA, const double* pB, const double* pX, double* pR,
                    double BNorm) {
    const std::size_t n = rA.mRows;
    csr_multiply(rA, pX, pR);
    parallel_for(n, [&](std::size_t i) { pR[i] = pB[i] - pR[i]; }, kSpsGrain);
    const double rn = std::sqrt(sparse_dot(pR, pR, n));
    return BNorm > 0.0 ? rn / BNorm : rn;
}

// The inverse diagonal, 1 where the diagonal is zero or not finite.
std::vector<double> sps_jacobi(const SparseCsr& rA) {
    std::vector<double> d = csr_diagonal(rA);
    for (double& v : d)
        v = (v != 0.0 && std::isfinite(v)) ? 1.0 / v : 1.0;
    return d;
}

void sps_check_square(const SparseCsr& rA, const char* pWho) {
    if (rA.mRows != rA.mCols)
        throw std::invalid_argument(std::string("meshio++: ") + pWho +
                                    ": the matrix is not square");
}

// Core preconditioned conjugate gradient body taking a preconditioner functor
// `rApplyPrecond(const double* pR, double* pZ)` and caller-owned scratch.
template <class TApply>
SparseSolveReport sps_pcg_core(const SparseCsr& rA, const double* pB, double* pX,
                               const TApply& rApplyPrecond, double* pR, double* pZ, double* pP,
                               double* pAp, double RelativeTolerance, std::size_t MaxIterations) {
    const std::size_t n = rA.mRows;
    SparseSolveReport rep;
    if (n == 0) {
        rep.mConverged = true;
        return rep;
    }
    if (MaxIterations == 0)
        MaxIterations = 2 * n + 10;
    const double bnorm = std::sqrt(sparse_dot(pB, pB, n));
    rep.mRelativeResidual = sps_residual(rA, pB, pX, pR, bnorm);
    if (bnorm == 0.0) {
        // b = 0: the solution of an SPD system is x = 0.
        std::fill(pX, pX + n, 0.0);
        rep.mConverged = true;
        rep.mRelativeResidual = 0.0;
        return rep;
    }
    if (rep.mRelativeResidual <= RelativeTolerance) {
        rep.mConverged = true;
        return rep;
    }
    rApplyPrecond(pR, pZ);
    std::copy(pZ, pZ + n, pP);
    double rz = sparse_dot(pR, pZ, n);
    for (std::size_t it = 1; it <= MaxIterations; ++it) {
        csr_multiply(rA, pP, pAp);
        const double pap = sparse_dot(pP, pAp, n);
        if (!(pap > 0.0) || !std::isfinite(pap))
            break;
        const double alpha = rz / pap;
        sps_axpy(alpha, pP, pX, n);
        sps_axpy(-alpha, pAp, pR, n);
        rep.mIterations = it;
        rep.mRelativeResidual = std::sqrt(sparse_dot(pR, pR, n)) / bnorm;
        if (rep.mRelativeResidual <= RelativeTolerance) {
            rep.mConverged = true;
            break;
        }
        rApplyPrecond(pR, pZ);
        const double rz_new = sparse_dot(pR, pZ, n);
        const double beta = rz_new / rz;
        rz = rz_new;
        parallel_for(n, [&](std::size_t i) { pP[i] = pZ[i] + beta * pP[i]; }, kSpsGrain);
    }
    // The recurrence's residual drifts from the true one; report the true one.
    rep.mRelativeResidual = sps_residual(rA, pB, pX, pR, bnorm);
    rep.mConverged = rep.mRelativeResidual <= 10.0 * RelativeTolerance;
    return rep;
}

// Multi-RHS PCG streaming `rA` once per iteration across all active columns.
// Every floating-point operation for each column is bit-identical to running
// `sps_pcg_core` independently for that column.
template <class TApply>
std::vector<SparseSolveReport> sps_pcg_block_core(const SparseCsr& rA, std::size_t NumRhs,
                                                  const double* const* ppB, double* const* ppX,
                                                  const TApply& rApplyPrecond,
                                                  double RelativeTolerance,
                                                  std::size_t MaxIterations) {
    const std::size_t n = rA.mRows;
    std::vector<SparseSolveReport> reps(NumRhs);
    if (NumRhs == 0)
        return reps;
    if (n == 0) {
        for (SparseSolveReport& rep : reps)
            rep.mConverged = true;
        return reps;
    }
    if (MaxIterations == 0)
        MaxIterations = 2 * n + 10;

    // On threaded builds with small-to-medium n, running independent RHS
    // columns in parallel gives full multi-core scaling with zero cross-column
    // synchronization and identical arithmetic to the serial solver.
    if (NumRhs > 1 && n <= kSpsGrain && parallel_workers() > 1) {
        parallel_for(
            NumRhs,
            [&](std::size_t j) {
                std::vector<double> r(n), z(n), p(n), ap(n);
                reps[j] = sps_pcg_core(rA, ppB[j], ppX[j], rApplyPrecond, r.data(), z.data(),
                                       p.data(), ap.data(), RelativeTolerance, MaxIterations);
            },
            std::size_t{1});
        return reps;
    }

    // Stream `rA` once per iteration across all active columns.
    std::vector<double> r_buf(NumRhs * n);
    std::vector<double> z_buf(NumRhs * n);
    std::vector<double> p_buf(NumRhs * n);
    std::vector<double> ap_buf(NumRhs * n);
    std::vector<double> bnorm(NumRhs, 0.0);
    std::vector<double> rz(NumRhs, 0.0);
    std::vector<std::size_t> active;
    active.reserve(NumRhs);

    for (std::size_t j = 0; j < NumRhs; ++j) {
        double* pr = r_buf.data() + j * n;
        double* pz = z_buf.data() + j * n;
        double* pp = p_buf.data() + j * n;
        bnorm[j] = std::sqrt(sparse_dot(ppB[j], ppB[j], n));
        reps[j].mRelativeResidual = sps_residual(rA, ppB[j], ppX[j], pr, bnorm[j]);
        if (bnorm[j] == 0.0) {
            std::fill(ppX[j], ppX[j] + n, 0.0);
            reps[j].mConverged = true;
            reps[j].mRelativeResidual = 0.0;
            continue;
        }
        if (reps[j].mRelativeResidual <= RelativeTolerance) {
            reps[j].mConverged = true;
            continue;
        }
        rApplyPrecond(pr, pz);
        std::copy(pz, pz + n, pp);
        rz[j] = sparse_dot(pr, pz, n);
        active.push_back(j);
    }

    std::vector<std::size_t> next_active;
    next_active.reserve(NumRhs);
    std::vector<double> row_sums(NumRhs);
    std::vector<const double*> p_ptrs(NumRhs);
    std::vector<double*> ap_ptrs(NumRhs);

    for (std::size_t it = 1; it <= MaxIterations && !active.empty(); ++it) {
        const std::size_t ka = active.size();
        for (std::size_t a = 0; a < ka; ++a) {
            p_ptrs[a] = p_buf.data() + active[a] * n;
            ap_ptrs[a] = ap_buf.data() + active[a] * n;
        }
        if (n <= kSpsGrain) {
            for (std::size_t row = 0; row < n; ++row) {
                for (std::size_t a = 0; a < ka; ++a)
                    row_sums[a] = 0.0;
                const std::size_t beg = rA.mStart[row], end = rA.mStart[row + 1];
                for (std::size_t e = beg; e < end; ++e) {
                    const std::size_t col = rA.mCol[e];
                    const double val = rA.mValue[e];
                    for (std::size_t a = 0; a < ka; ++a)
                        row_sums[a] += val * p_ptrs[a][col];
                }
                for (std::size_t a = 0; a < ka; ++a)
                    ap_ptrs[a][row] = row_sums[a];
            }
        } else {
            parallel_for(
                n,
                [&](std::size_t row) {
                    static thread_local std::vector<double> s;
                    s.assign(ka, 0.0);
                    const std::size_t beg = rA.mStart[row], end = rA.mStart[row + 1];
                    for (std::size_t e = beg; e < end; ++e) {
                        const std::size_t col = rA.mCol[e];
                        const double val = rA.mValue[e];
                        for (std::size_t a = 0; a < ka; ++a)
                            s[a] += val * p_ptrs[a][col];
                    }
                    for (std::size_t a = 0; a < ka; ++a)
                        ap_ptrs[a][row] = s[a];
                },
                kSpsGrain);
        }

        next_active.clear();
        for (std::size_t a = 0; a < ka; ++a) {
            const std::size_t j = active[a];
            double* pr = r_buf.data() + j * n;
            double* pz = z_buf.data() + j * n;
            double* pp = p_buf.data() + j * n;
            const double* pap = ap_ptrs[a];
            const double pap_dot = sparse_dot(pp, pap, n);
            if (!(pap_dot > 0.0) || !std::isfinite(pap_dot))
                continue;
            const double alpha = rz[j] / pap_dot;
            sps_axpy(alpha, pp, ppX[j], n);
            sps_axpy(-alpha, pap, pr, n);
            reps[j].mIterations = it;
            reps[j].mRelativeResidual = std::sqrt(sparse_dot(pr, pr, n)) / bnorm[j];
            if (reps[j].mRelativeResidual <= RelativeTolerance) {
                reps[j].mConverged = true;
                continue;
            }
            rApplyPrecond(pr, pz);
            const double rz_new = sparse_dot(pr, pz, n);
            const double beta = rz_new / rz[j];
            rz[j] = rz_new;
            parallel_for(n, [&](std::size_t i) { pp[i] = pz[i] + beta * pp[i]; }, kSpsGrain);
            next_active.push_back(j);
        }
        active.swap(next_active);
    }

    for (std::size_t j = 0; j < NumRhs; ++j) {
        if (bnorm[j] == 0.0 || reps[j].mIterations == 0)
            continue;
        double* pr = r_buf.data() + j * n;
        reps[j].mRelativeResidual = sps_residual(rA, ppB[j], ppX[j], pr, bnorm[j]);
        reps[j].mConverged = reps[j].mRelativeResidual <= 10.0 * RelativeTolerance;
    }
    return reps;
}

// Deterministic BFS-level-set Nested Dissection ordering on the symmetrised
// adjacency graph of `rA`, for O(n log n) fill on 2-D meshes (and reduced fill
// on 3-D meshes) in `CsrLdlt`.
void sps_nd_order(const SparseCsr& rA, std::vector<std::size_t>& rPerm,
                  std::vector<std::size_t>& rInvPerm) {
    const std::size_t n = rA.mRows;
    rPerm.clear();
    rPerm.reserve(n);
    rInvPerm.assign(n, 0);
    if (n == 0)
        return;

    // Build symmetrised loopless adjacency list (CSR).
    std::vector<std::size_t> adj_start(n + 1, 0);
    for (std::size_t r = 0; r < n; ++r) {
        for (std::size_t k = rA.mStart[r]; k < rA.mStart[r + 1]; ++k) {
            const std::size_t c = rA.mCol[k];
            if (c != r) {
                ++adj_start[r + 1];
                ++adj_start[c + 1];
            }
        }
    }
    for (std::size_t i = 0; i < n; ++i)
        adj_start[i + 1] += adj_start[i];
    std::vector<std::size_t> adj_col(adj_start[n]);
    {
        std::vector<std::size_t> cur = adj_start;
        for (std::size_t r = 0; r < n; ++r) {
            for (std::size_t k = rA.mStart[r]; k < rA.mStart[r + 1]; ++k) {
                const std::size_t c = rA.mCol[k];
                if (c != r) {
                    adj_col[cur[r]++] = c;
                    adj_col[cur[c]++] = r;
                }
            }
        }
    }
    // Deduplicate each row in place so adjacency is clean and compact.
    std::vector<std::size_t> sym_start(n + 1, 0);
    std::size_t write = 0;
    for (std::size_t r = 0; r < n; ++r) {
        auto beg = adj_col.begin() + static_cast<std::ptrdiff_t>(adj_start[r]);
        auto end = adj_col.begin() + static_cast<std::ptrdiff_t>(adj_start[r + 1]);
        std::sort(beg, end);
        end = std::unique(beg, end);
        for (auto it = beg; it != end; ++it)
            adj_col[write++] = *it;
        sym_start[r + 1] = write;
    }
    adj_col.resize(write);

    std::vector<std::uint32_t> sub_tag(n, 0);
    std::vector<std::uint32_t> vis_tag(n, 0);
    std::vector<std::int32_t> lvl(n, -1);
    std::vector<std::size_t> deg(n, 0);
    std::vector<std::size_t> bfs;
    bfs.reserve(n);
    std::uint32_t cur_sub = 0;
    std::uint32_t cur_vis = 0;

    constexpr std::size_t kMinSplit = 32;
    constexpr int kMaxDepth = 48;

    const auto next_sub = [&]() -> std::uint32_t {
        if (++cur_sub == 0) {
            std::fill(sub_tag.begin(), sub_tag.end(), 0);
            cur_sub = 1;
        }
        return cur_sub;
    };
    const auto next_vis = [&]() -> std::uint32_t {
        if (++cur_vis == 0) {
            std::fill(vis_tag.begin(), vis_tag.end(), 0);
            cur_vis = 1;
        }
        return cur_vis;
    };

    // Order a small leaf component by local RCM so its bandwidth stays tight.
    const auto order_leaf = [&](const std::vector<std::size_t>& rComp) {
        if (rComp.size() <= 4) {
            rPerm.insert(rPerm.end(), rComp.begin(), rComp.end());
            return;
        }
        const std::uint32_t stag = next_sub();
        for (std::size_t u : rComp)
            sub_tag[u] = stag;
        std::size_t seed = rComp.front();
        std::size_t min_d = n + 1;
        for (std::size_t u : rComp) {
            std::size_t d = 0;
            for (std::size_t k = sym_start[u]; k < sym_start[u + 1]; ++k)
                if (sub_tag[adj_col[k]] == stag)
                    ++d;
            deg[u] = d;
            if (d < min_d || (d == min_d && u < seed)) {
                min_d = d;
                seed = u;
            }
        }
        const std::uint32_t vtag = next_vis();
        bfs.clear();
        std::vector<std::size_t> nbrs;
        for (std::size_t start_u : rComp) {
            if (vis_tag[start_u] == vtag)
                continue;
            const std::size_t s = vis_tag[seed] != vtag ? seed : start_u;
            const std::size_t comp_beg = bfs.size();
            vis_tag[s] = vtag;
            bfs.push_back(s);
            for (std::size_t h = comp_beg; h < bfs.size(); ++h) {
                const std::size_t u = bfs[h];
                nbrs.clear();
                for (std::size_t k = sym_start[u]; k < sym_start[u + 1]; ++k) {
                    const std::size_t v = adj_col[k];
                    if (sub_tag[v] == stag && vis_tag[v] != vtag) {
                        vis_tag[v] = vtag;
                        nbrs.push_back(v);
                    }
                }
                std::sort(nbrs.begin(), nbrs.end(), [&](std::size_t a, std::size_t b) {
                    return deg[a] != deg[b] ? deg[a] < deg[b] : a < b;
                });
                for (std::size_t v : nbrs)
                    bfs.push_back(v);
            }
            for (std::size_t idx = bfs.size(); idx-- > comp_beg;)
                rPerm.push_back(bfs[idx]);
        }
    };

    const auto recurse = [&](auto&& rSelf, std::vector<std::size_t> rNodes, int Depth) -> void {
        if (rNodes.size() <= kMinSplit || Depth >= kMaxDepth) {
            order_leaf(rNodes);
            return;
        }
        const std::uint32_t stag = next_sub();
        for (std::size_t u : rNodes)
            sub_tag[u] = stag;

        // Split into connected components within `rNodes`.
        const std::uint32_t cvis = next_vis();
        std::vector<std::vector<std::size_t>> comps;
        for (std::size_t root : rNodes) {
            if (vis_tag[root] == cvis)
                continue;
            std::vector<std::size_t> comp;
            vis_tag[root] = cvis;
            comp.push_back(root);
            for (std::size_t h = 0; h < comp.size(); ++h) {
                const std::size_t u = comp[h];
                for (std::size_t k = sym_start[u]; k < sym_start[u + 1]; ++k) {
                    const std::size_t v = adj_col[k];
                    if (sub_tag[v] == stag && vis_tag[v] != cvis) {
                        vis_tag[v] = cvis;
                        comp.push_back(v);
                    }
                }
            }
            comps.push_back(std::move(comp));
        }

        if (comps.size() > 1) {
            rNodes.clear();
            rNodes.shrink_to_fit();
            for (auto& comp : comps)
                rSelf(rSelf, std::move(comp), Depth);
            return;
        }

        // Single connected component: bisect via pseudo-peripheral BFS level set.
        std::size_t seed = rNodes.front();
        std::size_t min_d = n + 1;
        for (std::size_t u : rNodes) {
            std::size_t d = 0;
            for (std::size_t k = sym_start[u]; k < sym_start[u + 1]; ++k)
                if (sub_tag[adj_col[k]] == stag)
                    ++d;
            deg[u] = d;
            if (d < min_d || (d == min_d && u < seed)) {
                min_d = d;
                seed = u;
            }
        }
        for (int pass = 0; pass < 2; ++pass) {
            const std::uint32_t vtag = next_vis();
            bfs.clear();
            vis_tag[seed] = vtag;
            lvl[seed] = 0;
            bfs.push_back(seed);
            std::size_t far = seed;
            std::int32_t max_l = 0;
            for (std::size_t h = 0; h < bfs.size(); ++h) {
                const std::size_t u = bfs[h];
                const std::int32_t lu = lvl[u];
                if (lu > max_l ||
                    (lu == max_l && (deg[u] < deg[far] || (deg[u] == deg[far] && u < far)))) {
                    max_l = lu;
                    far = u;
                }
                for (std::size_t k = sym_start[u]; k < sym_start[u + 1]; ++k) {
                    const std::size_t v = adj_col[k];
                    if (sub_tag[v] == stag && vis_tag[v] != vtag) {
                        vis_tag[v] = vtag;
                        lvl[v] = lu + 1;
                        bfs.push_back(v);
                    }
                }
            }
            seed = far;
        }

        // Final BFS from `seed` to get level structure.
        const std::uint32_t vtag = next_vis();
        bfs.clear();
        vis_tag[seed] = vtag;
        lvl[seed] = 0;
        bfs.push_back(seed);
        std::int32_t max_lvl = 0;
        for (std::size_t h = 0; h < bfs.size(); ++h) {
            const std::size_t u = bfs[h];
            const std::int32_t lu = lvl[u];
            if (lu > max_lvl)
                max_lvl = lu;
            for (std::size_t k = sym_start[u]; k < sym_start[u + 1]; ++k) {
                const std::size_t v = adj_col[k];
                if (sub_tag[v] == stag && vis_tag[v] != vtag) {
                    vis_tag[v] = vtag;
                    lvl[v] = lu + 1;
                    bfs.push_back(v);
                }
            }
        }
        if (max_lvl < 2) {
            order_leaf(rNodes);
            return;
        }
        std::int32_t mid = lvl[bfs[bfs.size() / 2]];
        if (mid == 0)
            mid = 1;
        if (mid == max_lvl)
            mid = max_lvl - 1;

        std::vector<std::size_t> part_a, part_b, sep;
        part_a.reserve(bfs.size() / 2);
        part_b.reserve(bfs.size() / 2);
        for (std::size_t u : rNodes) {
            const std::int32_t lu = lvl[u];
            if (lu < mid) {
                part_a.push_back(u);
            } else if (lu > mid) {
                part_b.push_back(u);
            } else {
                bool touches_b = false;
                for (std::size_t k = sym_start[u]; k < sym_start[u + 1]; ++k) {
                    const std::size_t v = adj_col[k];
                    if (sub_tag[v] == stag && lvl[v] > mid) {
                        touches_b = true;
                        break;
                    }
                }
                if (touches_b)
                    sep.push_back(u);
                else
                    part_a.push_back(u);
            }
        }
        if (part_a.empty() || part_b.empty()) {
            order_leaf(rNodes);
            return;
        }
        rNodes.clear();
        rNodes.shrink_to_fit();
        rSelf(rSelf, std::move(part_a), Depth + 1);
        rSelf(rSelf, std::move(part_b), Depth + 1);
        rPerm.insert(rPerm.end(), sep.begin(), sep.end());
    };

    std::vector<std::size_t> all(n);
    for (std::size_t i = 0; i < n; ++i)
        all[i] = i;
    recurse(recurse, std::move(all), 0);

    for (std::size_t i = 0; i < n; ++i)
        rInvPerm[rPerm[i]] = i;
}

}  // namespace

SparseCsr csr_from_triplets(std::size_t Rows, std::size_t Cols,
                            std::vector<SparseTriplet> rTriplets) {
    for (const SparseTriplet& t : rTriplets)
        if (t.mRow >= Rows || t.mCol >= Cols)
            throw std::invalid_argument(
                "meshio++: csr_from_triplets: an entry lies outside the matrix");
    std::stable_sort(rTriplets.begin(), rTriplets.end(),
                     [](const SparseTriplet& rL, const SparseTriplet& rR) {
                         return rL.mRow != rR.mRow ? rL.mRow < rR.mRow : rL.mCol < rR.mCol;
                     });
    SparseCsr a;
    a.mRows = Rows;
    a.mCols = Cols;
    a.mStart.assign(Rows + 1, 0);
    a.mCol.reserve(rTriplets.size());
    a.mValue.reserve(rTriplets.size());
    for (std::size_t k = 0; k < rTriplets.size(); ++k) {
        const SparseTriplet& t = rTriplets[k];
        if (k && t.mRow == rTriplets[k - 1].mRow && t.mCol == rTriplets[k - 1].mCol) {
            a.mValue.back() += t.mValue;
            continue;
        }
        a.mCol.push_back(t.mCol);
        a.mValue.push_back(t.mValue);
        ++a.mStart[t.mRow + 1];
    }
    for (std::size_t r = 0; r < Rows; ++r)
        a.mStart[r + 1] += a.mStart[r];
    return a;
}

void csr_multiply(const SparseCsr& rA, const double* pX, double* pY) {
    parallel_for(
        rA.mRows,
        [&](std::size_t r) {
            double s = 0.0;
            for (std::size_t k = rA.mStart[r]; k < rA.mStart[r + 1]; ++k)
                s += rA.mValue[k] * pX[rA.mCol[k]];
            pY[r] = s;
        },
        kSpsGrain);
}

std::vector<double> csr_diagonal(const SparseCsr& rA) {
    std::vector<double> d(std::min(rA.mRows, rA.mCols), 0.0);
    for (std::size_t r = 0; r < d.size(); ++r)
        for (std::size_t k = rA.mStart[r]; k < rA.mStart[r + 1]; ++k)
            if (rA.mCol[k] == r)
                d[r] = rA.mValue[k];
    return d;
}

double sparse_dot(const double* pA, const double* pB, std::size_t N) {
    const auto chunk = [&](std::size_t Begin, std::size_t End) {
        double s = 0.0;
        for (std::size_t i = Begin; i < End; ++i)
            s += pA[i] * pB[i];
        return s;
    };
    if (N <= kSpsGrain) {
        // The same chunks, summed the same way, without a parallel region.
        double acc = 0.0;
        for (std::size_t b = 0; b < N; b += kSparseChunk)
            acc = acc + chunk(b, std::min(N, b + kSparseChunk));
        return acc;
    }
    return parallel_reduce(N, kSparseChunk, 0.0, chunk, [](double L, double R) { return L + R; });
}

void CsrPcgWorkspace::Prepare(const SparseCsr& rA) {
    sps_check_square(rA, "CsrPcgWorkspace::Prepare");
    const std::size_t n = rA.mRows;
    mRows = n;
    mMinv = sps_jacobi(rA);
    mR.resize(n);
    mZ.resize(n);
    mP.resize(n);
    mAp.resize(n);
}

SparseSolveReport csr_pcg(const SparseCsr& rA, const double* pB, double* pX,
                          CsrPcgWorkspace& rWorkspace, double RelativeTolerance,
                          std::size_t MaxIterations) {
    sps_check_square(rA, "csr_pcg");
    const std::size_t n = rA.mRows;
    if (rWorkspace.mRows != n || rWorkspace.mMinv.size() != n)
        rWorkspace.Prepare(rA);
    const double* minv = rWorkspace.mMinv.data();
    return sps_pcg_core(
        rA, pB, pX,
        [minv, n](const double* pR, double* pZ) {
            parallel_for(n, [&](std::size_t i) { pZ[i] = minv[i] * pR[i]; }, kSpsGrain);
        },
        rWorkspace.mR.data(), rWorkspace.mZ.data(), rWorkspace.mP.data(), rWorkspace.mAp.data(),
        RelativeTolerance, MaxIterations);
}

SparseSolveReport csr_pcg(const SparseCsr& rA, const double* pB, double* pX,
                          double RelativeTolerance, std::size_t MaxIterations) {
    sps_check_square(rA, "csr_pcg");
    if (rA.mRows == 0) {
        SparseSolveReport rep;
        rep.mConverged = true;
        return rep;
    }
    CsrPcgWorkspace ws;
    ws.Prepare(rA);
    return csr_pcg(rA, pB, pX, ws, RelativeTolerance, MaxIterations);
}

std::vector<SparseSolveReport> csr_pcg_block(const SparseCsr& rA, std::size_t NumRhs,
                                             const double* const* ppB, double* const* ppX,
                                             CsrPcgWorkspace& rWorkspace, double RelativeTolerance,
                                             std::size_t MaxIterations) {
    sps_check_square(rA, "csr_pcg_block");
    if (NumRhs == 0)
        return {};
    if (NumRhs == 1)
        return {csr_pcg(rA, ppB[0], ppX[0], rWorkspace, RelativeTolerance, MaxIterations)};
    const std::size_t n = rA.mRows;
    if (rWorkspace.mRows != n || rWorkspace.mMinv.size() != n)
        rWorkspace.Prepare(rA);
    const double* minv = rWorkspace.mMinv.data();
    return sps_pcg_block_core(
        rA, NumRhs, ppB, ppX,
        [minv, n](const double* pR, double* pZ) {
            parallel_for(n, [&](std::size_t i) { pZ[i] = minv[i] * pR[i]; }, kSpsGrain);
        },
        RelativeTolerance, MaxIterations);
}

bool CsrIc0::Build(const SparseCsr& rA) {
    sps_check_square(rA, "CsrIc0::Build");
    const std::size_t n = rA.mRows;
    mRows = n;
    mJacobiInv = sps_jacobi(rA);
    mLStart.assign(n + 1, 0);
    mLCol.clear();
    mLVal.clear();
    mDiagInv.assign(n, 1.0);
    mValid = false;
    if (n == 0) {
        mValid = true;
        return true;
    }

    // Extract lower triangle (j <= i) with the diagonal always placed last in
    // each row.
    std::vector<double> a_lower;
    mLCol.reserve(rA.NumNonzeros() / 2 + n);
    a_lower.reserve(rA.NumNonzeros() / 2 + n);
    for (std::size_t i = 0; i < n; ++i) {
        double diag = 0.0;
        for (std::size_t k = rA.mStart[i]; k < rA.mStart[i + 1]; ++k) {
            const std::size_t col = rA.mCol[k];
            if (col < i) {
                mLCol.push_back(col);
                a_lower.push_back(rA.mValue[k]);
            } else if (col == i) {
                diag = rA.mValue[k];
            }
        }
        mLCol.push_back(i);
        a_lower.push_back(diag);
        mLStart[i + 1] = mLCol.size();
    }

    mLVal.resize(a_lower.size());
    std::vector<double> work(n, 0.0);
    const double shifts[] = {0.0, 1e-3, 1e-2, 1e-1, 1.0};
    for (const double alpha : shifts) {
        std::copy(a_lower.begin(), a_lower.end(), mLVal.begin());
        bool ok = true;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t beg = mLStart[i];
            const std::size_t end = mLStart[i + 1] - 1;  // diagonal slot
            for (std::size_t p = beg; p < end; ++p)
                work[mLCol[p]] = mLVal[p];
            double diag = mLVal[end] * (1.0 + alpha);
            for (std::size_t p = beg; p < end; ++p) {
                const std::size_t j = mLCol[p];
                double s = work[j];
                const std::size_t j_beg = mLStart[j];
                const std::size_t j_end = mLStart[j + 1] - 1;
                for (std::size_t q = j_beg; q < j_end; ++q)
                    s -= mLVal[q] * work[mLCol[q]];
                const double lij = s * mDiagInv[j];
                mLVal[p] = lij;
                work[j] = lij;
                diag -= lij * lij;
            }
            for (std::size_t p = beg; p < end; ++p)
                work[mLCol[p]] = 0.0;
            if (!(diag > 0.0) || !std::isfinite(diag)) {
                ok = false;
                break;
            }
            const double lii = std::sqrt(diag);
            mLVal[end] = lii;
            mDiagInv[i] = 1.0 / lii;
        }
        if (ok) {
            mValid = true;
            return true;
        }
    }
    mValid = false;
    return false;
}

void CsrIc0::Apply(const double* pR, double* pZ) const {
    const std::size_t n = mRows;
    if (!mValid) {
        const double* minv = mJacobiInv.data();
        parallel_for(n, [&](std::size_t i) { pZ[i] = minv[i] * pR[i]; }, kSpsGrain);
        return;
    }
    // Forward solve: L y = r (into pZ).
    for (std::size_t i = 0; i < n; ++i) {
        double s = pR[i];
        const std::size_t beg = mLStart[i];
        const std::size_t end = mLStart[i + 1] - 1;
        for (std::size_t k = beg; k < end; ++k)
            s -= mLVal[k] * pZ[mLCol[k]];
        pZ[i] = s * mDiagInv[i];
    }
    // Backward solve: L^T z = y (in-place on pZ, using row-oriented L as
    // column-oriented L^T).
    for (std::size_t i = n; i-- > 0;) {
        const double zi = pZ[i] * mDiagInv[i];
        pZ[i] = zi;
        const std::size_t beg = mLStart[i];
        const std::size_t end = mLStart[i + 1] - 1;
        for (std::size_t k = beg; k < end; ++k)
            pZ[mLCol[k]] -= mLVal[k] * zi;
    }
}

SparseSolveReport csr_ic0_pcg(const SparseCsr& rA, const CsrIc0& rPrecond, const double* pB,
                              double* pX, CsrPcgWorkspace& rWorkspace, double RelativeTolerance,
                              std::size_t MaxIterations) {
    sps_check_square(rA, "csr_ic0_pcg");
    const std::size_t n = rA.mRows;
    if (rWorkspace.mR.size() != n) {
        rWorkspace.mRows = n;
        rWorkspace.mR.resize(n);
        rWorkspace.mZ.resize(n);
        rWorkspace.mP.resize(n);
        rWorkspace.mAp.resize(n);
    }
    return sps_pcg_core(
        rA, pB, pX, [&](const double* pR, double* pZ) { rPrecond.Apply(pR, pZ); },
        rWorkspace.mR.data(), rWorkspace.mZ.data(), rWorkspace.mP.data(), rWorkspace.mAp.data(),
        RelativeTolerance, MaxIterations);
}

SparseSolveReport csr_ic0_pcg(const SparseCsr& rA, const double* pB, double* pX,
                              double RelativeTolerance, std::size_t MaxIterations) {
    sps_check_square(rA, "csr_ic0_pcg");
    if (rA.mRows == 0) {
        SparseSolveReport rep;
        rep.mConverged = true;
        return rep;
    }
    CsrIc0 precond;
    precond.Build(rA);
    CsrPcgWorkspace ws;
    return csr_ic0_pcg(rA, precond, pB, pX, ws, RelativeTolerance, MaxIterations);
}

std::vector<SparseSolveReport> csr_ic0_pcg_block(const SparseCsr& rA, const CsrIc0& rPrecond,
                                                 std::size_t NumRhs, const double* const* ppB,
                                                 double* const* ppX, CsrPcgWorkspace& rWorkspace,
                                                 double RelativeTolerance,
                                                 std::size_t MaxIterations) {
    sps_check_square(rA, "csr_ic0_pcg_block");
    if (NumRhs == 0)
        return {};
    if (NumRhs == 1)
        return {csr_ic0_pcg(rA, rPrecond, ppB[0], ppX[0], rWorkspace, RelativeTolerance,
                            MaxIterations)};
    return sps_pcg_block_core(
        rA, NumRhs, ppB, ppX, [&](const double* pR, double* pZ) { rPrecond.Apply(pR, pZ); },
        RelativeTolerance, MaxIterations);
}

bool CsrLdlt::Factor(const SparseCsr& rA, std::size_t MaxLNonzeros) {
    sps_check_square(rA, "CsrLdlt::Factor");
    const std::size_t n = rA.mRows;
    mRows = n;
    mValid = false;
    mLp.clear();
    mLi.clear();
    mLx.clear();
    mDInv.clear();
    if (n == 0) {
        mPerm.clear();
        mInvPerm.clear();
        mLp.assign(1, 0);
        mValid = true;
        return true;
    }

    sps_nd_order(rA, mPerm, mInvPerm);

    // Build upper triangle of P A P^T in column-compressed format (row <= col).
    std::vector<std::size_t> u_start(n + 1, 0);
    for (std::size_t r = 0; r < n; ++r) {
        const std::size_t pr = mInvPerm[r];
        for (std::size_t k = rA.mStart[r]; k < rA.mStart[r + 1]; ++k) {
            const std::size_t pc = mInvPerm[rA.mCol[k]];
            if (pr <= pc)
                ++u_start[pc + 1];
        }
    }
    for (std::size_t c = 0; c < n; ++c)
        u_start[c + 1] += u_start[c];
    std::vector<std::size_t> u_row(u_start[n]);
    std::vector<double> u_val(u_start[n]);
    {
        std::vector<std::size_t> cursor = u_start;
        for (std::size_t r = 0; r < n; ++r) {
            const std::size_t pr = mInvPerm[r];
            for (std::size_t k = rA.mStart[r]; k < rA.mStart[r + 1]; ++k) {
                const std::size_t pc = mInvPerm[rA.mCol[k]];
                if (pr <= pc) {
                    const std::size_t dst = cursor[pc]++;
                    u_row[dst] = pr;
                    u_val[dst] = rA.mValue[k];
                }
            }
        }
    }

    // Symbolic analysis: elimination tree and column counts of L (Davis's LDL).
    std::vector<std::size_t> parent(n, n);
    std::vector<std::size_t> flag(n, n);
    std::vector<std::size_t> lnz(n, 0);
    for (std::size_t k = 0; k < n; ++k) {
        parent[k] = n;
        flag[k] = k;
        lnz[k] = 0;
        for (std::size_t p = u_start[k]; p < u_start[k + 1]; ++p) {
            std::size_t i = u_row[p];
            if (i < k) {
                for (; flag[i] != k; i = parent[i]) {
                    if (parent[i] == n)
                        parent[i] = k;
                    ++lnz[i];
                    flag[i] = k;
                }
            }
        }
    }

    mLp.assign(n + 1, 0);
    for (std::size_t k = 0; k < n; ++k)
        mLp[k + 1] = mLp[k] + lnz[k];
    if (mLp[n] > MaxLNonzeros) {
        mLp.clear();
        return false;
    }

    // Numeric up-looking LDL^T factorisation.
    mLi.resize(mLp[n]);
    mLx.resize(mLp[n]);
    mDInv.resize(n);
    std::vector<double> y(n, 0.0);
    std::vector<std::size_t> pattern(n, 0);
    std::fill(flag.begin(), flag.end(), n);
    std::fill(lnz.begin(), lnz.end(), 0);

    for (std::size_t k = 0; k < n; ++k) {
        y[k] = 0.0;
        std::size_t top = n;
        flag[k] = k;
        lnz[k] = 0;
        for (std::size_t p = u_start[k]; p < u_start[k + 1]; ++p) {
            std::size_t i = u_row[p];
            if (i <= k) {
                y[i] += u_val[p];
                std::size_t len = 0;
                for (; flag[i] != k; i = parent[i]) {
                    pattern[len++] = i;
                    flag[i] = k;
                }
                while (len > 0)
                    pattern[--top] = pattern[--len];
            }
        }
        double d = y[k];
        y[k] = 0.0;
        for (; top < n; ++top) {
            const std::size_t i = pattern[top];
            const double yi = y[i];
            y[i] = 0.0;
            const std::size_t p2 = mLp[i] + lnz[i];
            for (std::size_t p = mLp[i]; p < p2; ++p)
                y[mLi[p]] -= mLx[p] * yi;
            const double l_ki = yi * mDInv[i];
            d -= l_ki * yi;
            mLi[p2] = k;
            mLx[p2] = l_ki;
            ++lnz[i];
        }
        if (!(d > 0.0) || !std::isfinite(d)) {
            mLp.clear();
            mLi.clear();
            mLx.clear();
            mDInv.clear();
            return false;
        }
        mDInv[k] = 1.0 / d;
    }
    mValid = true;
    return true;
}

void CsrLdlt::Solve(const double* pB, double* pX, std::vector<double>& rWork) const {
    const std::size_t n = mRows;
    if (n == 0)
        return;
    rWork.resize(n);
    double* w = rWork.data();
    const std::size_t* perm = mPerm.data();
    const std::size_t* lp = mLp.data();
    const std::size_t* li = mLi.data();
    const double* lx = mLx.data();
    const double* dinv = mDInv.data();
    for (std::size_t i = 0; i < n; ++i)
        w[i] = pB[perm[i]];
    // Forward solve: L w = P b (unit lower, column-oriented).
    for (std::size_t j = 0; j < n; ++j) {
        const double wj = w[j];
        const std::size_t beg = lp[j], end = lp[j + 1];
        for (std::size_t p = beg; p < end; ++p)
            w[li[p]] -= lx[p] * wj;
    }
    // Diagonal solve: D w = w.
    for (std::size_t j = 0; j < n; ++j)
        w[j] *= dinv[j];
    // Backward solve: L^T w = w.
    for (std::size_t j = n; j-- > 0;) {
        double wj = w[j];
        const std::size_t beg = lp[j], end = lp[j + 1];
        for (std::size_t p = beg; p < end; ++p)
            wj -= lx[p] * w[li[p]];
        w[j] = wj;
    }
    for (std::size_t i = 0; i < n; ++i)
        pX[perm[i]] = w[i];
}

void CsrLdlt::Solve(const double* pB, double* pX) const {
    std::vector<double> work;
    Solve(pB, pX, work);
}

void CsrLdlt::SolveBlock(std::size_t NumRhs, const double* const* ppB, double* const* ppX) const {
    const std::size_t n = mRows;
    if (NumRhs == 0 || n == 0)
        return;
    if (NumRhs == 1) {
        Solve(ppB[0], ppX[0]);
        return;
    }
    const std::size_t* perm = mPerm.data();
    const std::size_t* lp = mLp.data();
    const std::size_t* li = mLi.data();
    const double* lx = mLx.data();
    const double* dinv = mDInv.data();

    // Process up to 4 right-hand sides per chunk so each L nonzero is loaded
    // once per 4 columns while keeping every column bit-identical to `Solve`.
    const std::size_t num_chunks = (NumRhs + 3) / 4;
    parallel_for(
        num_chunks,
        [&](std::size_t chunk) {
            const std::size_t c0 = chunk * 4;
            const std::size_t cnt = std::min<std::size_t>(4, NumRhs - c0);
            static thread_local std::vector<double> work;
            if (cnt == 4) {
                work.resize(4 * n);
                double* w = work.data();
                const double* b0 = ppB[c0 + 0];
                const double* b1 = ppB[c0 + 1];
                const double* b2 = ppB[c0 + 2];
                const double* b3 = ppB[c0 + 3];
                for (std::size_t i = 0; i < n; ++i) {
                    const std::size_t pi = perm[i];
                    w[4 * i + 0] = b0[pi];
                    w[4 * i + 1] = b1[pi];
                    w[4 * i + 2] = b2[pi];
                    w[4 * i + 3] = b3[pi];
                }
                for (std::size_t j = 0; j < n; ++j) {
                    const double w0 = w[4 * j + 0];
                    const double w1 = w[4 * j + 1];
                    const double w2 = w[4 * j + 2];
                    const double w3 = w[4 * j + 3];
                    const std::size_t beg = lp[j], end = lp[j + 1];
                    for (std::size_t p = beg; p < end; ++p) {
                        const std::size_t r4 = 4 * li[p];
                        const double lxp = lx[p];
                        w[r4 + 0] -= lxp * w0;
                        w[r4 + 1] -= lxp * w1;
                        w[r4 + 2] -= lxp * w2;
                        w[r4 + 3] -= lxp * w3;
                    }
                }
                for (std::size_t j = 0; j < n; ++j) {
                    const double di = dinv[j];
                    w[4 * j + 0] *= di;
                    w[4 * j + 1] *= di;
                    w[4 * j + 2] *= di;
                    w[4 * j + 3] *= di;
                }
                for (std::size_t j = n; j-- > 0;) {
                    double w0 = w[4 * j + 0];
                    double w1 = w[4 * j + 1];
                    double w2 = w[4 * j + 2];
                    double w3 = w[4 * j + 3];
                    const std::size_t beg = lp[j], end = lp[j + 1];
                    for (std::size_t p = beg; p < end; ++p) {
                        const std::size_t r4 = 4 * li[p];
                        const double lxp = lx[p];
                        w0 -= lxp * w[r4 + 0];
                        w1 -= lxp * w[r4 + 1];
                        w2 -= lxp * w[r4 + 2];
                        w3 -= lxp * w[r4 + 3];
                    }
                    w[4 * j + 0] = w0;
                    w[4 * j + 1] = w1;
                    w[4 * j + 2] = w2;
                    w[4 * j + 3] = w3;
                }
                double* x0 = ppX[c0 + 0];
                double* x1 = ppX[c0 + 1];
                double* x2 = ppX[c0 + 2];
                double* x3 = ppX[c0 + 3];
                for (std::size_t i = 0; i < n; ++i) {
                    const std::size_t pi = perm[i];
                    x0[pi] = w[4 * i + 0];
                    x1[pi] = w[4 * i + 1];
                    x2[pi] = w[4 * i + 2];
                    x3[pi] = w[4 * i + 3];
                }
            } else {
                for (std::size_t k = 0; k < cnt; ++k)
                    Solve(ppB[c0 + k], ppX[c0 + k], work);
            }
        },
        std::size_t{1});
}

SparseSolveReport csr_bicgstab(const SparseCsr& rA, const double* pB, double* pX,
                               double RelativeTolerance, std::size_t MaxIterations) {
    sps_check_square(rA, "csr_bicgstab");
    const std::size_t n = rA.mRows;
    SparseSolveReport rep;
    if (n == 0) {
        rep.mConverged = true;
        return rep;
    }
    if (MaxIterations == 0)
        MaxIterations = 2 * n + 10;
    const std::vector<double> minv = sps_jacobi(rA);
    const double bnorm = std::sqrt(sparse_dot(pB, pB, n));
    std::vector<double> r(n), r0(n), p(n, 0.0), v(n, 0.0), s(n), t(n), ph(n), sh(n);
    rep.mRelativeResidual = sps_residual(rA, pB, pX, r.data(), bnorm);
    if (rep.mRelativeResidual <= RelativeTolerance) {
        rep.mConverged = true;
        return rep;
    }
    r0 = r;
    double rho = 1.0, alpha = 1.0, omega = 1.0;
    std::vector<double> best(pX, pX + n);
    double best_res = rep.mRelativeResidual;
    for (std::size_t it = 1; it <= MaxIterations; ++it) {
        const double rho_new = sparse_dot(r0.data(), r.data(), n);
        if (rho_new == 0.0 || !std::isfinite(rho_new) || omega == 0.0)
            break;
        const double beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;
        parallel_for(
            n, [&](std::size_t i) { p[i] = r[i] + beta * (p[i] - omega * v[i]); }, kSpsGrain);
        parallel_for(n, [&](std::size_t i) { ph[i] = minv[i] * p[i]; }, kSpsGrain);
        csr_multiply(rA, ph.data(), v.data());
        const double r0v = sparse_dot(r0.data(), v.data(), n);
        if (r0v == 0.0 || !std::isfinite(r0v))
            break;
        alpha = rho / r0v;
        parallel_for(n, [&](std::size_t i) { s[i] = r[i] - alpha * v[i]; }, kSpsGrain);
        rep.mIterations = it;
        const double snorm = std::sqrt(sparse_dot(s.data(), s.data(), n));
        if (snorm / bnorm <= RelativeTolerance) {
            sps_axpy(alpha, ph.data(), pX, n);
            break;
        }
        parallel_for(n, [&](std::size_t i) { sh[i] = minv[i] * s[i]; }, kSpsGrain);
        csr_multiply(rA, sh.data(), t.data());
        const double tt = sparse_dot(t.data(), t.data(), n);
        if (!(tt > 0.0) || !std::isfinite(tt))
            break;
        omega = sparse_dot(t.data(), s.data(), n) / tt;
        parallel_for(n, [&](std::size_t i) { pX[i] += alpha * ph[i] + omega * sh[i]; }, kSpsGrain);
        parallel_for(n, [&](std::size_t i) { r[i] = s[i] - omega * t[i]; }, kSpsGrain);
        const double res = std::sqrt(sparse_dot(r.data(), r.data(), n)) / bnorm;
        if (res < best_res) {
            best_res = res;
            std::copy(pX, pX + n, best.begin());
        }
        if (res <= RelativeTolerance)
            break;
    }
    std::vector<double> tmp(n);
    double res = sps_residual(rA, pB, pX, tmp.data(), bnorm);
    if (best_res < res) {
        std::copy(best.begin(), best.end(), pX);
        res = best_res;
    }
    rep.mRelativeResidual = res;
    rep.mConverged = res <= 10.0 * RelativeTolerance;
    return rep;
}

}  // namespace detail
}  // namespace meshioplusplus
