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
// Mesh comparison ("diff"): compare two meshes elementwise (ordered) or up to a
// spatial node correspondence (unordered) and produce a structured report with
// an overall verdict. Built entirely through the uniform mesh API so it
// compiles under every mesh backend. See operations/diff.hpp for the contract.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/operations/diff.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// --- elementwise array reduction --------------------------------------------

// Accumulated comparison of one array (or a chunk of it): worst absolute /
// relative error, its flat index, count over tolerance, and whether everything
// was bitwise-equal. Combine is associative + deterministic (ties -> lower
// index) so the result is identical across thread counts and backends.
struct DiffAcc {
    double mMaxAbs = 0.0;
    double mMaxRel = 0.0;
    std::int64_t mWorstIndex = -1;
    std::int64_t mNumExceeding = 0;
    bool mExact = true;
};

void diff_acc_consider(DiffAcc& rAcc, double a, double b, std::int64_t flat, double atol,
                       double rtol) {
    const double abs_err = std::fabs(a - b);
    if (a != b)
        rAcc.mExact = false;
    if (abs_err > rAcc.mMaxAbs || (abs_err == rAcc.mMaxAbs && rAcc.mWorstIndex == -1)) {
        rAcc.mMaxAbs = abs_err;
        rAcc.mWorstIndex = flat;
    }
    const double denom = std::fabs(a);
    const double rel_err = denom > 0.0 ? abs_err / denom : 0.0;
    if (rel_err > rAcc.mMaxRel)
        rAcc.mMaxRel = rel_err;
    if (abs_err > atol + rtol * denom)
        ++rAcc.mNumExceeding;
}

void diff_acc_combine(DiffAcc& rInto, const DiffAcc& rFrom) {
    if (rFrom.mWorstIndex != -1 &&
        (rFrom.mMaxAbs > rInto.mMaxAbs ||
         (rFrom.mMaxAbs == rInto.mMaxAbs &&
          (rInto.mWorstIndex == -1 || rFrom.mWorstIndex < rInto.mWorstIndex)))) {
        rInto.mMaxAbs = rFrom.mMaxAbs;
        rInto.mWorstIndex = rFrom.mWorstIndex;
    }
    rInto.mMaxRel = std::max(rInto.mMaxRel, rFrom.mMaxRel);
    rInto.mNumExceeding += rFrom.mNumExceeding;
    rInto.mExact = rInto.mExact && rFrom.mExact;
}

// Compare two arrays row-wise. When `pRowMapA` is non-null (unordered mode) row
// r of B is compared against row `pRowMapA[r]` of A; otherwise row r vs row r.
// Fills `rOut` (shape flags + reduction). The reduction runs in parallel over
// row chunks, then combines serially (deterministic).
ArrayDiff diff_compare_array(const NDArray& rA, const NDArray& rB, const std::int64_t* pRowMapA,
                             double atol, double rtol) {
    ArrayDiff out;
    out.mSizeA = static_cast<std::int64_t>(rA.Size());
    out.mSizeB = static_cast<std::int64_t>(rB.Size());
    const std::size_t cols_a = detail::cols(rA);
    const std::size_t cols_b = detail::cols(rB);
    const std::size_t rows_a = detail::rows(rA);
    const std::size_t rows_b = detail::rows(rB);
    // Column count must match; row count must match unless a correspondence maps
    // every B row into A.
    if (cols_a != cols_b || (pRowMapA == nullptr && rows_a != rows_b)) {
        out.mShapeMismatch = true;
        out.mExact = false;
        return out;
    }
    const std::size_t nrows = pRowMapA ? rows_b : rows_a;
    const std::size_t ncols = cols_a;
    if (nrows == 0 || ncols == 0)
        return out;

    // Chunk over rows; grain 1 so even a handful of chunks dispatch in parallel.
    constexpr std::size_t kRowsPerChunk = 4096;
    const std::size_t nchunks = (nrows + kRowsPerChunk - 1) / kRowsPerChunk;
    std::vector<DiffAcc> partials(nchunks);
    parallel_for(
        nchunks,
        [&](std::size_t chunk) {
            DiffAcc acc;
            const std::size_t r0 = chunk * kRowsPerChunk;
            const std::size_t r1 = std::min(r0 + kRowsPerChunk, nrows);
            for (std::size_t r = r0; r < r1; ++r) {
                const std::size_t ar = pRowMapA ? static_cast<std::size_t>(pRowMapA[r]) : r;
                for (std::size_t c = 0; c < ncols; ++c) {
                    const double a = detail::read_double(rA, ar * cols_a + c);
                    const double b = detail::read_double(rB, r * cols_b + c);
                    diff_acc_consider(acc, a, b, static_cast<std::int64_t>(r * ncols + c), atol,
                                      rtol);
                }
            }
            partials[chunk] = acc;
        },
        1);

    DiffAcc total;
    for (const DiffAcc& p : partials)
        diff_acc_combine(total, p);
    out.mMaxAbsError = total.mMaxAbs;
    out.mMaxRelError = total.mMaxRel;
    out.mWorstIndex = total.mWorstIndex;
    out.mNumExceeding = total.mNumExceeding;
    out.mExact = total.mExact;
    return out;
}

// --- cell block comparison --------------------------------------------------

// Collect the node ids of one cell (rectangular / polygon / polyhedron) into
// `rOut`, optionally remapping each id through `pNodeMap` (unordered mode maps B
// node ids into A's index space).
void diff_cell_nodes(const Mesh::CellView& rCb, std::size_t c, const std::int64_t* pNodeMap,
                     std::size_t nNodeMap, std::vector<std::int64_t>& rOut) {
    rOut.clear();
    auto push = [&](std::int64_t id) {
        if (pNodeMap && id >= 0 && static_cast<std::size_t>(id) < nNodeMap)
            rOut.push_back(pNodeMap[static_cast<std::size_t>(id)]);
        else
            rOut.push_back(id);
    };
    if (rCb.IsPolyhedron()) {
        for (std::size_t f = 0; f < rCb.NumFaces(c); ++f) {
            std::pair<const std::int64_t*, std::size_t> face = rCb.Face(c, f);
            for (std::size_t k = 0; k < face.second; ++k)
                push(face.first[k]);
        }
    } else if (rCb.IsRagged()) {
        const std::int64_t* row = rCb.Row(c);
        for (std::size_t k = 0; k < rCb.RowSize(c); ++k)
            push(row[k]);
    } else {
        const NDArray& conn = rCb.Conn();
        const std::size_t npc = rCb.NodesPerCell();
        for (std::size_t k = 0; k < npc; ++k)
            push(detail::read_int(conn, c * npc + k));
    }
}

// Compare cell block `block` of A and B (matched by position). When `pNodeMap`
// is non-null (unordered), B's node ids are remapped into A's space first.
BlockDiff diff_compare_block(const Mesh::CellView& rA, const Mesh::CellView& rB, std::size_t block,
                             const std::int64_t* pNodeMap, std::size_t nNodeMap,
                             std::int64_t max_reported, std::vector<std::string>& rMessages) {
    BlockDiff out;
    out.mBlock = static_cast<std::int64_t>(block);
    out.mTypeA = rA.Type();
    out.mTypeB = rB.Type();
    out.mCountA = static_cast<std::int64_t>(rA.NumCells());
    out.mCountB = static_cast<std::int64_t>(rB.NumCells());
    out.mTypeMismatch = out.mTypeA != out.mTypeB;
    out.mCountMismatch = out.mCountA != out.mCountB;
    if (out.mTypeMismatch || out.mCountMismatch)
        return out;

    const std::size_t nc = rA.NumCells();
    std::vector<char> mism(nc, 0);
    parallel_for(nc, [&](std::size_t c) {
        std::vector<std::int64_t> na, nb;
        diff_cell_nodes(rA, c, nullptr, 0, na);
        diff_cell_nodes(rB, c, pNodeMap, nNodeMap, nb);
        mism[c] = (na != nb) ? 1 : 0;
    });

    for (std::size_t c = 0; c < nc; ++c) {
        if (mism[c]) {
            ++out.mConnMismatchCount;
            if (static_cast<std::int64_t>(out.mFirstMismatches.size()) < max_reported)
                out.mFirstMismatches.push_back(static_cast<std::int64_t>(c));
        }
    }
    if (out.mConnMismatchCount > 0 && (rA.IsRagged() || rB.IsRagged()))
        rMessages.push_back("block " + std::to_string(block) + " (" + out.mTypeA +
                            "): ragged connectivity differs");
    return out;
}

// --- data section comparison ------------------------------------------------

// Sorted-key set difference into `rDiff.mOnlyInA` / `rDiff.mOnlyInB`, returning
// the shared keys (both name lists are already sorted).
std::vector<std::string> diff_key_diff(const std::vector<std::string>& rNamesA,
                                       const std::vector<std::string>& rNamesB, DataDiff& rDiff) {
    std::vector<std::string> shared;
    std::set_difference(rNamesA.begin(), rNamesA.end(), rNamesB.begin(), rNamesB.end(),
                        std::back_inserter(rDiff.mOnlyInA));
    std::set_difference(rNamesB.begin(), rNamesB.end(), rNamesA.begin(), rNamesA.end(),
                        std::back_inserter(rDiff.mOnlyInB));
    std::set_intersection(rNamesA.begin(), rNamesA.end(), rNamesB.begin(), rNamesB.end(),
                          std::back_inserter(shared));
    return shared;
}

// Fold an ArrayDiff `add` into the running per-key summary `rInto` (used to
// aggregate cell_data's per-block arrays into one summary per key).
void diff_fold_array(ArrayDiff& rInto, const ArrayDiff& rAdd, std::int64_t index_base) {
    rInto.mSizeA += rAdd.mSizeA;
    rInto.mSizeB += rAdd.mSizeB;
    rInto.mShapeMismatch = rInto.mShapeMismatch || rAdd.mShapeMismatch;
    rInto.mNumExceeding += rAdd.mNumExceeding;
    rInto.mExact = rInto.mExact && rAdd.mExact;
    rInto.mMaxRelError = std::max(rInto.mMaxRelError, rAdd.mMaxRelError);
    if (rAdd.mWorstIndex != -1 && rAdd.mMaxAbsError > rInto.mMaxAbsError) {
        rInto.mMaxAbsError = rAdd.mMaxAbsError;
        rInto.mWorstIndex = index_base + rAdd.mWorstIndex;
    }
}

// --- unordered node correspondence ------------------------------------------

// A packed integer key for a grid cell (three 21-bit signed-ish buckets).
std::int64_t diff_bucket_key(std::int64_t ix, std::int64_t iy, std::int64_t iz) {
    const std::uint64_t ux = static_cast<std::uint64_t>(ix + (1 << 20)) & 0x1fffff;
    const std::uint64_t uy = static_cast<std::uint64_t>(iy + (1 << 20)) & 0x1fffff;
    const std::uint64_t uz = static_cast<std::uint64_t>(iz + (1 << 20)) & 0x1fffff;
    return static_cast<std::int64_t>(ux | (uy << 21) | (uz << 42));
}

// Match every point of B to a unique point of A by spatial proximity within
// tolerance, using a bucket grid over A's bounding box (no O(N^2) scan). On
// success returns a full injective map `pointA = out[pointB]` and sets
// `rMaxResidual` to the worst matched coordinate error; on failure returns an
// empty vector and sets `rFailReason`.
std::vector<std::int64_t> diff_match_points(const Mesh& rA, const Mesh& rB, double atol,
                                            double rtol, double& rMaxResidual,
                                            std::string& rFailReason) {
    rMaxResidual = 0.0;
    const std::size_t na = rA.NumPoints();
    const std::size_t nb = rB.NumPoints();
    if (na != nb) {
        rFailReason = "point counts differ (" + std::to_string(na) + " vs " + std::to_string(nb) +
                      "); no one-to-one correspondence";
        return {};
    }
    const std::size_t pdim = rA.PointDim();
    if (rB.PointDim() != pdim) {
        rFailReason = "point dimensions differ; no correspondence";
        return {};
    }
    const std::size_t gdim = std::min<std::size_t>(pdim, 3);
    const NDArray& pa = rA.Points();
    const NDArray& pb = rB.Points();

    // Bounding box + per-axis grid cell size (>= the match tolerance so a match
    // within tolerance always lands in an adjacent cell).
    double lo[3] = {0, 0, 0}, h[3] = {1, 1, 1};
    for (std::size_t d = 0; d < gdim; ++d) {
        double mn = std::numeric_limits<double>::infinity();
        double mx = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < na; ++i) {
            const double c = detail::read_double(pa, i * pdim + d);
            mn = std::min(mn, c);
            mx = std::max(mx, c);
        }
        if (!std::isfinite(mn)) {
            mn = 0.0;
            mx = 0.0;
        }
        lo[d] = mn;
        const double maxabs = std::max(std::fabs(mn), std::fabs(mx));
        double cell = atol + rtol * maxabs;
        if (!(cell > 0.0))
            cell = (mx > mn) ? (mx - mn) : 1.0;
        h[d] = cell;
    }
    auto bucket_of = [&](const NDArray& p, std::size_t i, std::int64_t out[3]) {
        for (std::size_t d = 0; d < 3; ++d) {
            if (d < gdim) {
                const double c = detail::read_double(p, i * pdim + d);
                out[d] = static_cast<std::int64_t>(std::floor((c - lo[d]) / h[d]));
            } else {
                out[d] = 0;
            }
        }
    };

    std::unordered_map<std::int64_t, std::vector<std::int64_t>> grid;
    grid.reserve(na * 2);
    for (std::size_t i = 0; i < na; ++i) {
        std::int64_t b[3];
        bucket_of(pa, i, b);
        grid[diff_bucket_key(b[0], b[1], b[2])].push_back(static_cast<std::int64_t>(i));
    }

    auto within_tol = [&](std::size_t ai, std::size_t bi, double& rWorst) -> bool {
        double worst = 0.0;
        for (std::size_t d = 0; d < pdim; ++d) {
            const double a = detail::read_double(pa, ai * pdim + d);
            const double bb = detail::read_double(pb, bi * pdim + d);
            const double e = std::fabs(a - bb);
            if (e > atol + rtol * std::fabs(a))
                return false;
            worst = std::max(worst, e);
        }
        rWorst = worst;
        return true;
    };

    std::vector<std::int64_t> a_of_b(nb, -1);
    std::vector<char> used(na, 0);
    for (std::size_t j = 0; j < nb; ++j) {
        std::int64_t bb[3];
        bucket_of(pb, j, bb);
        std::int64_t best = -1;
        double best_err = std::numeric_limits<double>::infinity();
        for (std::int64_t dx = -1; dx <= 1; ++dx) {
            for (std::int64_t dy = -1; dy <= 1; ++dy) {
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    auto it = grid.find(diff_bucket_key(bb[0] + dx, bb[1] + dy, bb[2] + dz));
                    if (it == grid.end())
                        continue;
                    for (std::int64_t ai : it->second) {
                        if (used[static_cast<std::size_t>(ai)])
                            continue;
                        double err;
                        if (within_tol(static_cast<std::size_t>(ai), j, err) && err < best_err) {
                            best_err = err;
                            best = ai;
                        }
                    }
                }
            }
        }
        if (best < 0) {
            rFailReason = "point " + std::to_string(j) +
                          " of the second mesh has no unmatched neighbour within tolerance";
            return {};
        }
        used[static_cast<std::size_t>(best)] = 1;
        a_of_b[j] = best;
        rMaxResidual = std::max(rMaxResidual, best_err);
    }
    return a_of_b;
}

// --- verdict ----------------------------------------------------------------

void diff_bump(DiffVerdict& rV, DiffVerdict candidate) {
    if (static_cast<int>(candidate) > static_cast<int>(rV))
        rV = candidate;
}

DiffVerdict diff_array_verdict(const ArrayDiff& rA) {
    if (rA.mShapeMismatch || rA.mNumExceeding > 0)
        return DiffVerdict::Different;
    if (!rA.mExact)
        return DiffVerdict::EqualWithinTolerance;
    return DiffVerdict::Identical;
}

void diff_compare_data_section(const std::vector<std::string>& rNamesA,
                               const std::vector<std::string>& rNamesB, DataDiff& rOut,
                               const Mesh& rA, const Mesh& rB, const std::int64_t* pRowMapA,
                               bool cell_data, double atol, double rtol) {
    std::vector<std::string> shared = diff_key_diff(rNamesA, rNamesB, rOut);
    for (const std::string& name : shared) {
        ArrayDiff summary;
        summary.mName = name;
        if (cell_data) {
            const std::size_t nba = rA.CellDataNumBlocks(name);
            const std::size_t nbb = rB.CellDataNumBlocks(name);
            if (nba != nbb) {
                summary.mShapeMismatch = true;
                summary.mExact = false;
            } else {
                std::int64_t base = 0;
                for (std::size_t b = 0; b < nba; ++b) {
                    const NDArray& aa = rA.CellData(name, b);
                    const NDArray& bb = rB.CellData(name, b);
                    // Cell blocks are not reordered, so cell_data compares
                    // index-by-index even in unordered mode.
                    ArrayDiff part = diff_compare_array(aa, bb, nullptr, atol, rtol);
                    diff_fold_array(summary, part, base);
                    base += static_cast<std::int64_t>(aa.Size());
                }
            }
        } else {
            summary =
                diff_compare_array(rA.PointData(name), rB.PointData(name), pRowMapA, atol, rtol);
            summary.mName = name;
        }
        rOut.mShared.push_back(std::move(summary));
    }
}

}  // namespace

// --- public API -------------------------------------------------------------

const std::string& diff_verdict_name(DiffVerdict verdict) {
    static const std::string kIdentical = "identical";
    static const std::string kEqual = "equal within tolerance";
    static const std::string kDifferent = "different";
    switch (verdict) {
        case DiffVerdict::Identical:
            return kIdentical;
        case DiffVerdict::EqualWithinTolerance:
            return kEqual;
        case DiffVerdict::Different:
            return kDifferent;
    }
    return kDifferent;
}

// Compare two meshes' named regions. Identity is (name, kind): the same name
// may legitimately label a node set and an element set. Membership comparison
// is exact -- canonical entries make two equal groups bitwise equal, so no
// tolerance is involved.
void diff_compare_regions(const Mesh& rA, const Mesh& rB, RegionDiff& rOut) {
    auto label = [](const Region& r) { return r.mName + " (" + region_kind_name(r.mKind) + ")"; };
    auto find = [](const Mesh& rMesh, const Region& r) -> const Region* {
        const std::size_t i = rMesh.FindRegion(r.mName, r.mKind);
        return i == Mesh::npos ? nullptr : &rMesh.Region(i);
    };

    for (std::size_t i = 0; i < rA.NumRegions(); ++i) {
        const Region& a = rA.Region(i);
        const Region* b = find(rB, a);
        if (b == nullptr)
            rOut.mOnlyInA.push_back(label(a));
        else if (!regions_equal(a, *b))
            rOut.mChanged.push_back(label(a));
    }
    for (std::size_t i = 0; i < rB.NumRegions(); ++i) {
        const Region& b = rB.Region(i);
        if (find(rA, b) == nullptr)
            rOut.mOnlyInB.push_back(label(b));
    }
    std::sort(rOut.mOnlyInA.begin(), rOut.mOnlyInA.end());
    std::sort(rOut.mOnlyInB.begin(), rOut.mOnlyInB.end());
    std::sort(rOut.mChanged.begin(), rOut.mChanged.end());
}

DiffReport diff(const Mesh& rA, const Mesh& rB, const DiffOptions& rOpts) {
    DiffReport rep;
    rep.mUnordered = rOpts.unordered;
    rep.mNumPointsA = static_cast<std::int64_t>(rA.NumPoints());
    rep.mNumPointsB = static_cast<std::int64_t>(rB.NumPoints());
    rep.mNumBlocksA = static_cast<std::int64_t>(rA.NumCellBlocks());
    rep.mNumBlocksB = static_cast<std::int64_t>(rB.NumCellBlocks());

    const double atol = rOpts.atol;
    const double rtol = rOpts.rtol;

    // Unordered: build the node correspondence up front. On failure, everything
    // downstream is Different and we say why.
    std::vector<std::int64_t> a_of_b;
    const std::int64_t* row_map = nullptr;
    const std::int64_t* node_map = nullptr;
    if (rOpts.unordered) {
        double residual = 0.0;
        std::string reason;
        a_of_b = diff_match_points(rA, rB, atol, rtol, residual, reason);
        if (a_of_b.empty()) {
            rep.mCorrespondenceFailed = true;
            rep.mVerdict = DiffVerdict::Different;
            rep.mMessages.push_back("unordered matching failed: " + reason);
            rep.mPointCountMismatch = rep.mNumPointsA != rep.mNumPointsB;
            return rep;
        }
        row_map = a_of_b.data();
        node_map = a_of_b.data();
        // Points are matched within tolerance by construction; reflect the
        // worst residual so a nonzero drift downgrades the verdict.
        rep.mPoints.mMaxAbsError = residual;
        rep.mPoints.mExact = residual == 0.0;
        rep.mMessages.push_back("unordered: matched " + std::to_string(rep.mNumPointsB) +
                                " points by proximity");
    }

    // --- points ---
    if (rA.NumPoints() != rB.NumPoints() || rA.PointDim() != rB.PointDim()) {
        rep.mPointCountMismatch = true;
        rep.mPoints.mShapeMismatch = true;
        rep.mPoints.mExact = false;
        rep.mPoints.mSizeA = static_cast<std::int64_t>(rA.Points().Size());
        rep.mPoints.mSizeB = static_cast<std::int64_t>(rB.Points().Size());
    } else if (!rOpts.unordered) {
        rep.mPoints = diff_compare_array(rA.Points(), rB.Points(), nullptr, atol, rtol);
    } else {
        // Unordered: compare matched coordinate rows to record the full summary.
        rep.mPoints = diff_compare_array(rA.Points(), rB.Points(), row_map, atol, rtol);
    }

    // --- cells ---
    rep.mBlockCountMismatch = rep.mNumBlocksA != rep.mNumBlocksB;
    const std::size_t nblocks = std::min(rA.NumCellBlocks(), rB.NumCellBlocks());
    const std::size_t nnode = a_of_b.size();
    for (std::size_t b = 0; b < nblocks; ++b) {
        rep.mBlocks.push_back(diff_compare_block(rA.Cells(b), rB.Cells(b), b, node_map, nnode,
                                                 rOpts.max_reported, rep.mMessages));
    }

    // --- data ---
    diff_compare_data_section(rA.PointDataNames(), rB.PointDataNames(), rep.mPointData, rA, rB,
                              row_map, /*cell_data=*/false, atol, rtol);
    diff_compare_data_section(rA.CellDataNames(), rB.CellDataNames(), rep.mCellData, rA, rB,
                              nullptr, /*cell_data=*/true, atol, rtol);
    diff_compare_data_section(rA.FieldDataNames(), rB.FieldDataNames(), rep.mFieldData, rA, rB,
                              nullptr, /*cell_data=*/false, atol, rtol);

    diff_compare_regions(rA, rB, rep.mRegions);

    // --- verdict ---
    DiffVerdict v = DiffVerdict::Identical;
    if (rep.mPointCountMismatch || rep.mBlockCountMismatch)
        v = DiffVerdict::Different;
    diff_bump(v, diff_array_verdict(rep.mPoints));
    for (const BlockDiff& bd : rep.mBlocks) {
        if (bd.mTypeMismatch || bd.mCountMismatch || bd.mConnMismatchCount > 0)
            v = DiffVerdict::Different;
    }
    for (const DataDiff* dd : {&rep.mPointData, &rep.mCellData, &rep.mFieldData}) {
        if (!dd->mOnlyInA.empty() || !dd->mOnlyInB.empty())
            v = DiffVerdict::Different;
        for (const ArrayDiff& ad : dd->mShared)
            diff_bump(v, diff_array_verdict(ad));
    }
    rep.mVerdictIgnoringRegions = v;
    // A named group differing is a structural difference, not a numerical one,
    // so it goes straight to Different rather than through diff_bump.
    if (rep.mRegions.Differs())
        v = DiffVerdict::Different;
    rep.mVerdict = v;
    return rep;
}

bool meshes_equal(const Mesh& rA, const Mesh& rB, double atol, double rtol) {
    DiffOptions opts;
    opts.atol = atol;
    opts.rtol = rtol;
    return diff(rA, rB, opts).mVerdictIgnoringRegions != DiffVerdict::Different;
}

}  // namespace meshioplusplus
