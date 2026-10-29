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
// Reflect a mesh across a plane with every cell right side out (roadmap
// 3.2.5). See operations/mirror.hpp for the contract.

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
#include "../detail/cell_flip.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/mirror.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

// How one block's cells are flipped: a fixed table, a polygon ring or a
// polyhedron's face rings.
struct MirBlockFlip {
    const std::vector<int>* pPerm = nullptr;
    const std::vector<int>* pFacets = nullptr;
    std::vector<int> mPolygonPerm;  // a rectangular polygon block's own table
    bool mRaggedPolygon = false;
    bool mPolyhedron = false;
};

template <class CellViewT>
MirBlockFlip mir_block_flip(const CellViewT& rCells) {
    MirBlockFlip f;
    const std::string type(rCells.Type());
    const CellType ct = cell_type_from_name(type);
    if (rCells.IsPolyhedron()) {
        f.mPolyhedron = true;
        return f;
    }
    if (ct == CellType::Polygon) {
        if (rCells.IsRagged()) {
            f.mRaggedPolygon = true;
            return f;
        }
        // A polygon with rectangular storage: every row has the same size.
        const int npc = static_cast<int>(rCells.NodesPerCell());
        f.mPolygonPerm.push_back(0);
        for (int i = npc - 1; i >= 1; --i)
            f.mPolygonPerm.push_back(i);
        return f;
    }
    f.pPerm = rCells.IsRagged() ? nullptr : detail::cell_flip_permutation(ct);
    f.pFacets = detail::cell_flip_facet_map(ct);
    if (f.pPerm == nullptr ||
        static_cast<std::size_t>(cell_type_num_nodes(ct)) != rCells.NodesPerCell())
        throw Unsupported("meshio++: mirror: no orientation flip for cell type '" + type +
                          "' (linear and quadratic types, polygons and polyhedra only)");
    return f;
}

// The local facet a flipped cell gives facet F of a cell with N facets (a
// polygon's N is its node count), or -1 when F is not a facet of the cell.
std::int64_t mir_flip_facet(const MirBlockFlip& rFlip, std::int64_t F, std::int64_t N) {
    if (F < 0)
        return -1;
    if (rFlip.mPolyhedron)
        return F < N ? F : -1;
    if (rFlip.mRaggedPolygon || !rFlip.mPolygonPerm.empty())
        return F < N ? N - 1 - F : -1;
    if (rFlip.pFacets == nullptr || F >= static_cast<std::int64_t>(rFlip.pFacets->size()))
        return -1;
    return (*rFlip.pFacets)[static_cast<std::size_t>(F)];
}

// Reflect rows [Start, Start + Rows) of a float array of 3 (v' = H v) or 9
// (A' = H A H, row-major) columns in place.
void mir_reflect_rows(NDArray& rA, std::size_t Start, std::size_t Rows, std::size_t Cols,
                      const std::array<double, 9>& rH) {
    detail::dispatch_dtype(rA.Dtype(), [&]<class T>() {
        T* p = rA.As<T>();
        parallel_for_bw(Rows, [&](std::size_t r) {
            const std::size_t b = (Start + r) * Cols;
            if (Cols == 3) {
                const double v[3] = {static_cast<double>(p[b]), static_cast<double>(p[b + 1]),
                                     static_cast<double>(p[b + 2])};
                for (std::size_t i = 0; i < 3; ++i)
                    p[b + i] =
                        static_cast<T>(rH[i * 3] * v[0] + rH[i * 3 + 1] * v[1] + rH[i * 3 + 2] * v[2]);
                return;
            }
            double a[9], t[9];
            for (std::size_t k = 0; k < 9; ++k)
                a[k] = static_cast<double>(p[b + k]);
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j)
                    t[i * 3 + j] = rH[i * 3] * a[j] + rH[i * 3 + 1] * a[3 + j] + rH[i * 3 + 2] * a[6 + j];
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j)
                    p[b + i * 3 + j] = static_cast<T>(t[i * 3] * rH[j * 3] +
                                                      t[i * 3 + 1] * rH[j * 3 + 1] +
                                                      t[i * 3 + 2] * rH[j * 3 + 2]);
        });
    });
}

// Gather rows `rIdx` of a data array and reflect those from `Start` on, when
// asked and the array is a float vector or tensor.
NDArray mir_carry(const NDArray& rSrc, std::size_t SrcRows, const std::vector<std::int64_t>& rIdx,
                  std::size_t Start, bool Reflect, const std::array<double, 9>& rH) {
    NDArray out = detail::subset_gather_rows(rSrc, rIdx);
    const std::size_t cols = SrcRows > 0 ? rSrc.Size() / SrcRows : 0;
    if (Reflect && detail::is_float_dtype(rSrc.Dtype()) && (cols == 3 || cols == 9) &&
        rIdx.size() > Start)
        mir_reflect_rows(out, Start, rIdx.size() - Start, cols, rH);
    return out;
}

std::vector<std::int64_t> mir_iota(std::size_t N) {
    std::vector<std::int64_t> v(N);
    for (std::size_t i = 0; i < N; ++i)
        v[i] = static_cast<std::int64_t>(i);
    return v;
}

NDArray mir_entries(const std::vector<std::int64_t>& rFlat, bool Side) {
    const std::size_t n = Side ? rFlat.size() / 2 : rFlat.size();
    NDArray a = Side ? NDArray(DType::Int64, {n, 2}) : NDArray(DType::Int64, {n});
    std::copy(rFlat.begin(), rFlat.end(), a.As<std::int64_t>());
    return a;
}

}  // namespace

MirrorResult mirror(const Mesh& rMesh, const MirrorOptions& rOptions) {
    const MirrorOptions& o = rOptions;
    if (!std::isfinite(o.mTolerance) || o.mTolerance < 0.0)
        throw std::invalid_argument("meshio++: mirror: tolerance must be finite and non-negative");
    double len2 = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        if (!std::isfinite(o.mOrigin[d]) || !std::isfinite(o.mNormal[d]))
            throw std::invalid_argument("meshio++: mirror: origin and normal must be finite");
        len2 += o.mNormal[d] * o.mNormal[d];
    }
    if (!(len2 > 0.0))
        throw std::invalid_argument("meshio++: mirror: the normal must not be zero");
    const double len = std::sqrt(len2);
    const std::array<double, 3> nrm{{o.mNormal[0] / len, o.mNormal[1] / len, o.mNormal[2] / len}};

    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    for (std::size_t d = std::min<std::size_t>(dim, 3); d < 3; ++d)
        if (nrm[d] != 0.0)
            throw std::invalid_argument(
                "meshio++: mirror: the points have " + std::to_string(dim) +
                " coordinates, so the normal's component " + std::to_string(d) +
                " must be zero (the image would leave their space)");
    const std::size_t sdim = std::min<std::size_t>(dim, 3);
    const std::array<double, 9> H{{1.0 - 2.0 * nrm[0] * nrm[0], -2.0 * nrm[0] * nrm[1],
                                   -2.0 * nrm[0] * nrm[2], -2.0 * nrm[1] * nrm[0],
                                   1.0 - 2.0 * nrm[1] * nrm[1], -2.0 * nrm[1] * nrm[2],
                                   -2.0 * nrm[2] * nrm[0], -2.0 * nrm[2] * nrm[1],
                                   1.0 - 2.0 * nrm[2] * nrm[2]}};

    // Every block's flip, checked before anything is built.
    std::vector<MirBlockFlip> flips;
    for (const auto cb : rMesh.CellRange())
        flips.push_back(mir_block_flip(cb));

    // Signed distances to the plane and the points on it.
    const NDArray& pts = rMesh.Points();
    std::vector<double> sdist(n, 0.0);
    std::array<double, 3> lo{{0.0, 0.0, 0.0}}, hi{{0.0, 0.0, 0.0}};
    for (std::size_t p = 0; p < n; ++p) {
        double s = 0.0;
        for (std::size_t d = 0; d < sdim; ++d) {
            const double x = detail::read_double(pts, p * dim + d);
            s += (x - o.mOrigin[d]) * nrm[d];
            lo[d] = p == 0 ? x : std::min(lo[d], x);
            hi[d] = p == 0 ? x : std::max(hi[d], x);
        }
        sdist[p] = s;
    }
    double diag2 = 0.0;
    for (std::size_t d = 0; d < sdim; ++d)
        diag2 += (hi[d] - lo[d]) * (hi[d] - lo[d]);
    const double tol = o.mTolerance * std::sqrt(diag2);
    std::vector<char> on_plane(n, 0);
    MirrorResult res;
    for (std::size_t p = 0; p < n; ++p)
        if (std::abs(sdist[p]) <= tol) {
            on_plane[p] = 1;
            ++res.mNumPointsOnPlane;
        }

    const bool keep = o.mKeepOriginal;
    const bool weld = keep && o.mWeld;

    // The image of each point, and the rows the output points come from.
    std::vector<std::int64_t> image(n);
    std::vector<std::int64_t> point_rows = mir_iota(n);
    for (std::size_t p = 0; p < n; ++p) {
        if (!keep || (weld && on_plane[p])) {
            image[p] = static_cast<std::int64_t>(p);
            continue;
        }
        image[p] = static_cast<std::int64_t>(point_rows.size());
        point_rows.push_back(static_cast<std::int64_t>(p));
    }
    const std::size_t first_image = keep ? n : 0;  // the first reflected output point
    {
        NDArray out_pts = detail::subset_gather_rows(pts, point_rows);
        const std::size_t nimg = point_rows.size() - first_image;
        detail::dispatch_dtype(out_pts.Dtype(), [&]<class T>() {
            T* q = out_pts.As<T>();
            parallel_for_bw(nimg, [&](std::size_t i) {
                const std::size_t row = first_image + i;
                const double s = sdist[static_cast<std::size_t>(point_rows[row])];
                for (std::size_t d = 0; d < sdim; ++d)
                    q[row * dim + d] =
                        static_cast<T>(static_cast<double>(q[row * dim + d]) - 2.0 * s * nrm[d]);
            });
        });
        res.mMesh.AssignPoints(std::move(out_pts));
    }

    // The cells: per block, its own cells when kept, then the reflected ones;
    // a reflected cell lying on the plane is dropped when welding.
    const std::size_t nb = rMesh.NumCellBlocks();
    std::vector<std::vector<std::int64_t>> mirrored(nb);  // per block, the kept sources
    std::vector<std::vector<std::int64_t>> rank(nb);      // per cell, its copy's rank or -1
    for (std::size_t b = 0; b < nb; ++b) {
        const auto cb = rMesh.Cells(b);
        const MirBlockFlip& r_flip = flips[b];
        const std::size_t nc = cb.NumCells();
        const auto all_on = [&](const std::int64_t* pNodes, std::size_t K) {
            for (std::size_t k = 0; k < K; ++k)
                if (!on_plane[static_cast<std::size_t>(pNodes[k])])
                    return false;
            return K > 0;
        };
        std::vector<std::int64_t>& r_rank = rank[b];
        r_rank.assign(nc, -1);
        std::vector<std::int64_t>& r_kept = mirrored[b];
        const std::string type(cb.Type());
        if (r_flip.mPolyhedron) {
            std::vector<std::int64_t> flat, rows{0}, faces{0};
            const auto put = [&](std::size_t c, bool Flip) {
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    const auto face = cb.Face(c, f);
                    for (std::size_t k = 0; k < face.second; ++k) {
                        const std::size_t at = Flip && k > 0 ? face.second - k : k;
                        const std::int64_t v = face.first[at];
                        flat.push_back(Flip ? image[static_cast<std::size_t>(v)] : v);
                    }
                    rows.push_back(static_cast<std::int64_t>(flat.size()));
                }
                faces.push_back(static_cast<std::int64_t>(rows.size() - 1));
            };
            for (std::size_t c = 0; c < nc; ++c) {
                bool on = weld;
                for (std::size_t f = 0; f < cb.NumFaces(c) && on; ++f) {
                    const auto face = cb.Face(c, f);
                    on = all_on(face.first, face.second);
                }
                if (!on) {
                    r_rank[c] = static_cast<std::int64_t>(r_kept.size());
                    r_kept.push_back(static_cast<std::int64_t>(c));
                }
            }
            if (keep)
                for (std::size_t c = 0; c < nc; ++c)
                    put(c, false);
            for (const std::int64_t c : r_kept)
                put(static_cast<std::size_t>(c), true);
            res.mMesh.AddPolyhedronBlock(type, std::move(flat), std::move(rows), std::move(faces));
        } else if (r_flip.mRaggedPolygon) {
            std::vector<std::int64_t> flat, rows{0};
            for (std::size_t c = 0; c < nc; ++c)
                if (!(weld && all_on(cb.Row(c), cb.RowSize(c)))) {
                    r_rank[c] = static_cast<std::int64_t>(r_kept.size());
                    r_kept.push_back(static_cast<std::int64_t>(c));
                }
            if (keep)
                for (std::size_t c = 0; c < nc; ++c) {
                    flat.insert(flat.end(), cb.Row(c), cb.Row(c) + cb.RowSize(c));
                    rows.push_back(static_cast<std::int64_t>(flat.size()));
                }
            for (const std::int64_t cc : r_kept) {
                const auto c = static_cast<std::size_t>(cc);
                const std::int64_t* row = cb.Row(c);
                const std::size_t k = cb.RowSize(c);
                for (std::size_t i = 0; i < k; ++i)
                    flat.push_back(image[static_cast<std::size_t>(row[i == 0 ? 0 : k - i])]);
                rows.push_back(static_cast<std::int64_t>(flat.size()));
            }
            res.mMesh.AddPolygonBlock(type, std::move(flat), std::move(rows));
        } else {
            const std::vector<int>& perm =
                r_flip.mPolygonPerm.empty() ? *r_flip.pPerm : r_flip.mPolygonPerm;
            const std::size_t npc = cb.NodesPerCell();
            const NDArray& conn = cb.Conn();
            std::vector<std::int64_t> row(npc);
            for (std::size_t c = 0; c < nc; ++c) {
                bool on = weld;
                for (std::size_t k = 0; k < npc && on; ++k)
                    on = on_plane[static_cast<std::size_t>(detail::read_int(conn, c * npc + k))] != 0;
                if (!on || npc == 0) {
                    r_rank[c] = static_cast<std::int64_t>(r_kept.size());
                    r_kept.push_back(static_cast<std::int64_t>(c));
                }
            }
            const std::size_t base = keep ? nc : 0;
            NDArray out(DType::Int64, {base + r_kept.size(), npc});
            std::int64_t* q = out.As<std::int64_t>();
            parallel_for_bw(base, [&](std::size_t c) {
                for (std::size_t k = 0; k < npc; ++k)
                    q[c * npc + k] = detail::read_int(conn, c * npc + k);
            });
            parallel_for_bw(r_kept.size(), [&](std::size_t r) {
                const auto c = static_cast<std::size_t>(r_kept[r]);
                for (std::size_t k = 0; k < npc; ++k) {
                    const std::int64_t v =
                        detail::read_int(conn, c * npc + static_cast<std::size_t>(perm[k]));
                    q[(base + r) * npc + k] = image[static_cast<std::size_t>(v)];
                }
            });
            res.mMesh.AddCellBlock(type, std::move(out));
        }
        res.mNumFlippedCells += static_cast<std::int64_t>(r_kept.size());
        res.mNumDroppedCells += static_cast<std::int64_t>(nc - r_kept.size());
    }

    // Data, row by row.
    const bool reflect = o.mReflectVectorData;
    for (const std::string& name : rMesh.PointDataNames())
        res.mMesh.AddPointData(
            name, mir_carry(rMesh.PointData(name), n, point_rows, first_image, reflect, H));
    for (const std::string& name : rMesh.CellDataNames()) {
        std::vector<NDArray> blocks;
        for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name) && b < nb; ++b) {
            const std::size_t nc = rMesh.Cells(b).NumCells();
            std::vector<std::int64_t> rows = keep ? mir_iota(nc) : std::vector<std::int64_t>{};
            rows.insert(rows.end(), mirrored[b].begin(), mirrored[b].end());
            blocks.push_back(
                mir_carry(rMesh.CellData(name, b), nc, rows, keep ? nc : 0, reflect, H));
        }
        res.mMesh.AddCellData(name, std::move(blocks));
    }
    for (const std::string& name : rMesh.FieldDataNames()) {
        NDArray a = rMesh.FieldData(name);
        a.MakeOwned();
        res.mMesh.AddFieldData(name, std::move(a));
    }
    for (std::size_t i = 0; i < rMesh.NumPropertySets(); ++i)
        res.mMesh.AddPropertySet(rMesh.GetPropertySet(i));

    // Regions. A cell (b, c) is out_base[b] + c when kept, and its copy
    // out_base[b] + (kept ? n_b : 0) + rank.
    const std::vector<std::int64_t> in_bases = detail::block_bases(rMesh);
    std::vector<std::int64_t> out_bases(nb + 1, 0);
    for (std::size_t b = 0; b < nb; ++b)
        out_bases[b + 1] = out_bases[b] +
                           (keep ? static_cast<std::int64_t>(rMesh.Cells(b).NumCells()) : 0) +
                           static_cast<std::int64_t>(mirrored[b].size());
    const auto original_of = [&](std::size_t B, std::int64_t C) { return out_bases[B] + C; };
    const auto copy_of = [&](std::size_t B, std::int64_t C) -> std::int64_t {
        const std::int64_t r = rank[B][static_cast<std::size_t>(C)];
        if (r < 0)
            return -1;
        return out_bases[B] +
               (keep ? static_cast<std::int64_t>(rMesh.Cells(B).NumCells()) : 0) + r;
    };
    std::vector<std::int64_t> facet;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r_in = rMesh.Region(i);
        const std::int64_t* e = r_in.Entries();
        const std::size_t ne = r_in.NumEntries();
        std::vector<std::int64_t> flat;
        std::int64_t dropped = 0, invalid = 0;
        if (r_in.mKind == RegionKind::Point) {
            for (std::size_t k = 0; k < ne; ++k) {
                if (e[k] < 0 || static_cast<std::size_t>(e[k]) >= n) {
                    ++invalid;
                    continue;
                }
                flat.push_back(e[k]);
                if (keep)
                    flat.push_back(image[static_cast<std::size_t>(e[k])]);
            }
        } else if (r_in.mKind == RegionKind::Cell) {
            for (std::size_t k = 0; k < ne; ++k) {
                const auto [b, c] = detail::global_to_block_row(in_bases, e[k]);
                if (b >= nb) {
                    ++invalid;
                    continue;
                }
                if (keep)
                    flat.push_back(original_of(b, c));
                const std::int64_t m = copy_of(b, c);
                if (m >= 0)
                    flat.push_back(m);
            }
        } else {
            for (std::size_t k = 0; k < ne; ++k) {
                const std::int64_t g = e[2 * k], f = e[2 * k + 1];
                const auto [b, c] = detail::global_to_block_row(in_bases, g);
                if (b >= nb) {
                    ++invalid;
                    continue;
                }
                const auto cb = rMesh.Cells(b);
                const auto uc = static_cast<std::size_t>(c);
                const std::int64_t nfac =
                    flips[b].mPolyhedron
                        ? static_cast<std::int64_t>(cb.NumFaces(uc))
                        : (flips[b].mRaggedPolygon
                               ? static_cast<std::int64_t>(cb.RowSize(uc))
                               : static_cast<std::int64_t>(cb.NodesPerCell()));
                const std::int64_t nf = mir_flip_facet(flips[b], f, nfac);
                if (nf < 0) {
                    ++invalid;
                    continue;
                }
                if (weld) {
                    CellType ft = CellType::Custom;
                    bool on = detail::facet_nodes(rMesh, g, f, ft, facet) && !facet.empty();
                    for (const std::int64_t v : facet)
                        on = on && on_plane[static_cast<std::size_t>(v)];
                    if (on) {
                        ++dropped;
                        continue;
                    }
                }
                if (keep) {
                    flat.push_back(original_of(b, c));
                    flat.push_back(f);
                }
                const std::int64_t m = copy_of(b, c);
                if (m >= 0) {
                    flat.push_back(m);
                    flat.push_back(nf);
                }
            }
        }
        res.mNumDroppedSides += dropped;
        if (dropped > 0)
            res.mWarnings.push_back("region '" + r_in.mName + "': " + std::to_string(dropped) +
                                    " facet(s) on the mirror plane are interior and were dropped");
        if (invalid > 0)
            res.mWarnings.push_back("region '" + r_in.mName + "': " + std::to_string(invalid) +
                                    " entr(ies) naming no entity of the mesh were dropped");
        Region r_out(r_in.mName, r_in.mKind, r_in.mDim, r_in.mTag,
                     mir_entries(flat, r_in.mKind == RegionKind::Side));
        res.mMesh.AddRegion(std::move(r_out));
    }
    return res;
}

}  // namespace meshioplusplus
