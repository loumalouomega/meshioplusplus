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
#include <cstring>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "crease_edges.hpp"
#include "decimate_shared.hpp"
#include "ndarray_utils.hpp"

namespace meshioplusplus {
namespace detail {

DecimStop decim_resolve_stop(double TargetRatio, std::int64_t TargetCount, double MaxError,
                             const std::string& rPrefix, const char* pCountName) {
    const int num_set =
        (TargetRatio >= 0.0 ? 1 : 0) + (TargetCount >= 0 ? 1 : 0) + (MaxError >= 0.0 ? 1 : 0);
    if (num_set != 1)
        throw std::invalid_argument(rPrefix + "exactly one of target_ratio, " + pCountName +
                                    ", max_error must be set; got " + std::to_string(num_set));
    DecimStop stop;
    if (TargetRatio >= 0.0) {
        if (!(TargetRatio > 0.0 && TargetRatio <= 1.0))
            throw std::invalid_argument(rPrefix + "target_ratio must lie in (0, 1]; got " +
                                        std::to_string(TargetRatio));
        stop.mUseTarget = true;  // the count is filled in once it is known
    } else if (TargetCount >= 0) {
        stop.mUseTarget = true;
        stop.mTarget = TargetCount;
    } else {
        stop.mUseMaxError = true;
        stop.mMaxError = MaxError;
    }
    return stop;
}

void decim_check_frozen(std::size_t FrozenSize, std::size_t NumPoints, const std::string& rPrefix) {
    if (FrozenSize != 0 && FrozenSize != NumPoints)
        throw std::invalid_argument(rPrefix + "frozen mask has " + std::to_string(FrozenSize) +
                                    " entries but the mesh has " + std::to_string(NumPoints) +
                                    " points");
}

std::int64_t decim_ratio_target(double Ratio, std::size_t Count) {
    return std::max<std::int64_t>(
        1, static_cast<std::int64_t>(Ratio * static_cast<double>(Count) + 0.5));
}

std::vector<DecimFloatData> decim_float_working_copies(const Mesh& rMesh, std::size_t n) {
    std::vector<DecimFloatData> float_data;
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (!detail::is_float_dtype(a.Dtype()) || detail::rows(a) != n || n == 0)
            continue;
        DecimFloatData fd;
        fd.mName = name;
        fd.mNumComponents = a.Size() / n;
        fd.mValues.resize(a.Size());
        parallel_for_bw(a.Size(),
                        [&](std::size_t i) { fd.mValues[i] = detail::read_double(a, i); });
        float_data.push_back(std::move(fd));
    }
    return float_data;
}

std::vector<std::uint8_t> decim_pin_mask(const std::vector<std::uint8_t>& rFrozen,
                                         const std::vector<std::uint8_t>& rBoundary,
                                         bool PreserveBoundary, bool PreserveFeatures,
                                         double FeatureAngleDeg,
                                         const std::vector<std::int64_t>& rCorners,
                                         const std::vector<double>& rNormals, std::size_t n) {
    std::vector<std::uint8_t> pinned(n, 0);
    if (!rFrozen.empty())
        for (std::size_t i = 0; i < n; ++i)
            pinned[i] = rFrozen[i] ? 1 : 0;
    if (PreserveBoundary)
        for (std::size_t i = 0; i < n; ++i)
            if (rBoundary[i])
                pinned[i] = 1;
    if (PreserveFeatures) {
        // The shared per-edge crease test (v16.23.0): only the two faces that
        // share an edge are compared, and non-manifold edges pin too.
        detail::pin_crease_endpoints(
            detail::crease_edges_triangles(rCorners, rNormals, FeatureAngleDeg), pinned);
    }
    return pinned;
}

bool decim_any_normal_flip(const std::vector<double>& rXyz,
                           const std::vector<std::int64_t>& rCorners,
                           const std::vector<std::int64_t>& rFa,
                           const std::vector<std::int64_t>& rFb,
                           const std::vector<std::int64_t>& rShared, std::int64_t a, std::int64_t b,
                           const double* pPlaced, std::vector<std::int64_t>& rMerged) {
    rMerged.clear();
    std::merge(rFa.begin(), rFa.end(), rFb.begin(), rFb.end(), std::back_inserter(rMerged));
    std::size_t si = 0;
    for (const std::int64_t f : rMerged) {
        while (si < rShared.size() && rShared[si] < f)
            ++si;
        if (si < rShared.size() && rShared[si] == f)
            continue;  // dies in the collapse: no constraint
        const std::int64_t* c = rCorners.data() + static_cast<std::size_t>(f) * 3;
        double n0[3];
        double n1[3];
        decim_face_normal(rXyz, c, a, b, nullptr, n0);
        if (n0[0] == 0.0 && n0[1] == 0.0 && n0[2] == 0.0)
            continue;
        decim_face_normal(rXyz, c, a, b, pPlaced, n1);
        if (n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2] <= 0.0)
            return true;
    }
    return false;
}

void decim_commit_blend(std::vector<double>& rXyz, std::vector<double>& rQuads,
                        std::vector<DecimFloatData>& rFloatData, const double* pPlaced,
                        std::int64_t s, std::int64_t r) {
    // Blend parameter: project the placed point onto the edge, clamped to
    // [0, 1]. Exact for Midpoint/Endpoint, an approximation for Optimal
    // (the optimal point need not lie on the edge). A pinned survivor has
    // placed == its own position, so t == 0 and the data provably keeps the
    // survivor's values.
    {
        const double* xs = rXyz.data() + static_cast<std::size_t>(s) * 3;
        const double* xr = rXyz.data() + static_cast<std::size_t>(r) * 3;
        const double dx = xr[0] - xs[0];
        const double dy = xr[1] - xs[1];
        const double dz = xr[2] - xs[2];
        const double den = dx * dx + dy * dy + dz * dz;
        double t;
        if (den == 0.0) {
            t = 0.5;
        } else {
            t = ((pPlaced[0] - xs[0]) * dx + (pPlaced[1] - xs[1]) * dy +
                 (pPlaced[2] - xs[2]) * dz) /
                den;
            if (t < 0.0)
                t = 0.0;
            if (t > 1.0)
                t = 1.0;
        }
        for (DecimFloatData& fd : rFloatData) {
            double* vals = fd.mValues.data();
            const std::size_t so = static_cast<std::size_t>(s) * fd.mNumComponents;
            const std::size_t ro = static_cast<std::size_t>(r) * fd.mNumComponents;
            for (std::size_t k = 0; k < fd.mNumComponents; ++k)
                vals[so + k] = vals[so + k] + t * (vals[ro + k] - vals[so + k]);
        }
    }

    rXyz[static_cast<std::size_t>(s) * 3] = pPlaced[0];
    rXyz[static_cast<std::size_t>(s) * 3 + 1] = pPlaced[1];
    rXyz[static_cast<std::size_t>(s) * 3 + 2] = pPlaced[2];
    {
        double* qs = rQuads.data() + static_cast<std::size_t>(s) * 10;
        const double* qr = rQuads.data() + static_cast<std::size_t>(r) * 10;
        for (int i = 0; i < 10; ++i)
            qs[i] += qr[i];
    }
}

NDArray decim_emit_points(const NDArray& rPoints, const std::vector<double>& rXyz,
                          const std::vector<char>& rUsed, std::size_t n, std::size_t num_used,
                          std::size_t dim) {
    NDArray new_points = NDArray::Uninit(rPoints.Dtype(), {num_used, dim});
    detail::dispatch_dtype(rPoints.Dtype(), [&]<class T>() {
        T* dst = new_points.As<T>();
        std::size_t w = 0;
        for (std::size_t i = 0; i < n; ++i) {
            if (!rUsed[i])
                continue;
            for (std::size_t d = 0; d < dim && d < 3; ++d)
                dst[w * dim + d] = static_cast<T>(rXyz[i * 3 + d]);
            ++w;
        }
    });
    return new_points;
}

void decim_emit_data(const Mesh& rSrc, const std::vector<std::uint8_t>& rAlive,
                     const std::vector<std::int64_t>& rBlockBase, const char* pSkipName,
                     const std::vector<DecimFloatData>& rFloatData, const std::vector<char>& rUsed,
                     std::size_t n, std::size_t num_used, Mesh& rOut) {
    const std::size_t nblocks = rSrc.NumCellBlocks();

    // cell_data: the surviving cells' own rows.
    for (const std::string& name : rSrc.CellDataNames()) {
        if (pSkipName != nullptr && name == pSkipName)
            continue;
        const std::size_t ndata = rSrc.CellDataNumBlocks(name);
        std::vector<NDArray> blocks;
        blocks.reserve(ndata);
        std::size_t bi = 0;
        for (const auto cb : rSrc.CellRange()) {
            if (bi >= ndata)
                break;
            const NDArray& a = rSrc.CellData(name, bi);
            const std::size_t nc = cb.NumCells();
            if (ndata != nblocks || detail::rows(a) != nc || nc == 0) {
                blocks.push_back(detail::data_owned_copy(a));
                ++bi;
                continue;
            }
            const std::size_t base = static_cast<std::size_t>(rBlockBase[bi]);
            std::size_t alive_in_block = 0;
            for (std::size_t c = 0; c < nc; ++c)
                if (rAlive[base + c])
                    ++alive_in_block;
            std::vector<std::size_t> shape = a.Shape();
            shape[0] = alive_in_block;
            const std::size_t row_bytes = a.Nbytes() / nc;
            NDArray sub = NDArray::Uninit(a.Dtype(), std::move(shape));
            const std::byte* src = a.Data();
            std::byte* dst = sub.Data();
            std::size_t w = 0;
            for (std::size_t c = 0; c < nc; ++c)
                if (rAlive[base + c])
                    std::memcpy(dst + (w++) * row_bytes, src + c * row_bytes, row_bytes);
            blocks.push_back(std::move(sub));
            ++bi;
        }
        rOut.AddCellData(name, std::move(blocks));
    }

    // point_data: float-kind arrays from the blended working copies (cast back
    // to the source dtype), everything else as the survivor's own byte rows.
    for (const std::string& name : rSrc.PointDataNames()) {
        const NDArray& a = rSrc.PointData(name);
        if (detail::rows(a) != n) {
            rOut.AddPointData(name, detail::data_owned_copy(a));
            continue;
        }
        std::vector<std::size_t> shape = a.Shape();
        shape[0] = num_used;
        if (detail::is_float_dtype(a.Dtype())) {
            const DecimFloatData* p_fd = nullptr;
            for (const DecimFloatData& fd : rFloatData)
                if (fd.mName == name) {
                    p_fd = &fd;
                    break;
                }
            const std::size_t ncomp = p_fd != nullptr ? p_fd->mNumComponents : 0;
            NDArray b = NDArray::Uninit(a.Dtype(), std::move(shape));
            detail::dispatch_dtype(a.Dtype(), [&]<class T>() {
                T* dst = b.As<T>();
                std::size_t w = 0;
                for (std::size_t i = 0; i < n; ++i) {
                    if (!rUsed[i])
                        continue;
                    for (std::size_t k = 0; k < ncomp; ++k)
                        dst[w * ncomp + k] = static_cast<T>(p_fd->mValues[i * ncomp + k]);
                    ++w;
                }
            });
            rOut.AddPointData(name, std::move(b));
        } else {
            const std::size_t row_bytes = n == 0 ? 0 : a.Nbytes() / n;
            NDArray b = NDArray::Uninit(a.Dtype(), std::move(shape));
            const std::byte* src = a.Data();
            std::byte* dst = b.Data();
            std::size_t w = 0;
            for (std::size_t i = 0; i < n; ++i)
                if (rUsed[i])
                    std::memcpy(dst + (w++) * row_bytes, src + i * row_bytes, row_bytes);
            rOut.AddPointData(name, std::move(b));
        }
    }
    for (const std::string& name : rSrc.FieldDataNames())
        rOut.AddFieldData(name, detail::data_owned_copy(rSrc.FieldData(name)));
}

NDArray decim_point_map(const std::vector<std::int64_t>& rSuccessor,
                        const std::vector<std::int64_t>& rRemap, std::size_t n) {
    NDArray pm = NDArray::Uninit(DType::Int64, {n});
    std::int64_t* dst = pm.As<std::int64_t>();
    for (std::size_t i = 0; i < n; ++i) {
        std::int64_t v = static_cast<std::int64_t>(i);
        while (rSuccessor[static_cast<std::size_t>(v)] >= 0)
            v = rSuccessor[static_cast<std::size_t>(v)];
        dst[i] = rRemap[static_cast<std::size_t>(v)];
    }
    return pm;
}

template <int K>
void decim_emit_blocks(const Mesh& rIn, const Mesh& rSrc, const char* pParentName,
                       const std::vector<std::uint8_t>& rAlive,
                       const std::vector<std::int64_t>& rConn,
                       const std::vector<std::int64_t>& rBases,
                       const std::vector<std::int64_t>& rRemap, CellType Type, Mesh& rOut,
                       std::vector<NDArray>& rCellMaps) {
    constexpr std::size_t kCorners = K;
    const std::size_t nblocks = rSrc.NumCellBlocks();
    const bool has_parent = pParentName != nullptr && rSrc.HasCellData(pParentName) &&
                            rSrc.CellDataNumBlocks(pParentName) == nblocks;
    std::size_t bi = 0;
    std::vector<std::size_t> in_cells;
    for (const auto cb : rIn.CellRange())
        in_cells.push_back(cb.NumCells());
    for (const auto cb : rSrc.CellRange()) {
        const std::size_t nc = cb.NumCells();
        const std::size_t base = static_cast<std::size_t>(rBases[bi]);
        std::size_t alive_in_block = 0;
        for (std::size_t c = 0; c < nc; ++c)
            if (rAlive[base + c])
                ++alive_in_block;
        NDArray block = NDArray::Uninit(DType::Int64, {alive_in_block, kCorners});
        std::int64_t* dst = block.As<std::int64_t>();
        const NDArray* p_parent = has_parent ? &rSrc.CellData(pParentName, bi) : nullptr;
        std::vector<std::int64_t> cell_map(bi < in_cells.size() ? in_cells[bi] : nc, -1);
        std::size_t w = 0;
        for (std::size_t c = 0; c < nc; ++c) {
            if (!rAlive[base + c])
                continue;
            for (std::size_t k = 0; k < kCorners; ++k)
                dst[w * kCorners + k] =
                    rRemap[static_cast<std::size_t>(rConn[(base + c) * kCorners + k])];
            const std::int64_t parent =
                p_parent != nullptr ? read_int(*p_parent, c) : static_cast<std::int64_t>(c);
            if (parent >= 0 && static_cast<std::size_t>(parent) < cell_map.size() &&
                cell_map[static_cast<std::size_t>(parent)] < 0)
                cell_map[static_cast<std::size_t>(parent)] = static_cast<std::int64_t>(w);
            ++w;
        }
        rOut.AddCellBlock(cell_type_name(Type), std::move(block));
        rCellMaps.push_back(NDArrayUtils::Int64Vector(cell_map));
        ++bi;
    }
}

template void decim_emit_blocks<3>(const Mesh&, const Mesh&, const char*,
                                   const std::vector<std::uint8_t>&,
                                   const std::vector<std::int64_t>&,
                                   const std::vector<std::int64_t>&,
                                   const std::vector<std::int64_t>&, CellType, Mesh&,
                                   std::vector<NDArray>&);
template void decim_emit_blocks<4>(const Mesh&, const Mesh&, const char*,
                                   const std::vector<std::uint8_t>&,
                                   const std::vector<std::int64_t>&,
                                   const std::vector<std::int64_t>&,
                                   const std::vector<std::int64_t>&, CellType, Mesh&,
                                   std::vector<NDArray>&);

void decim_remap_regions(const Mesh& rIn, Mesh& rOut, const NDArray& rPointMap,
                         const std::vector<NDArray>& rCellMaps, CellMapKind Kind,
                         const char* pOpName) {
    RegionRemap rmap;
    rmap.pPointMap = &rPointMap;
    rmap.mCellMapKind = Kind;
    rmap.pCellMaps = &rCellMaps;
    rmap.mOpName = pOpName;
    remap_regions(rIn, rOut, rmap);
}

}  // namespace detail
}  // namespace meshioplusplus
