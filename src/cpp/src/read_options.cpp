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
#include <cstddef>
#include <limits>
#include <string>

// Project includes
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

std::size_t ReadOptions::ResolveTimeStep(std::size_t NumSteps) const {
    // A file with no recorded steps still has "the data", conceptually step 0 --
    // asking for it must not become an error just because the format omitted a
    // time array. Anything else on such a file is genuinely unanswerable.
    if (NumSteps == 0) {
        if (mTimeStep == 0 || mTimeStep == -1)
            return 0;
        throw ReadError("meshio++: time step " + std::to_string(mTimeStep) +
                        " requested, but this file carries no time steps");
    }

    const long long n = static_cast<long long>(NumSteps);
    // Negative counts from the end (-1 = last), which is the only way to say
    // "the final state" without knowing the count up front.
    const long long resolved = mTimeStep < 0 ? n + mTimeStep : mTimeStep;
    if (resolved < 0 || resolved >= n)
        throw ReadError("meshio++: time step " + std::to_string(mTimeStep) +
                        " is out of range: this file has " + std::to_string(NumSteps) +
                        (NumSteps == 1 ? " step" : " steps"));
    return static_cast<std::size_t>(resolved);
}

MeshMetadata metadata_from_mesh(const Mesh& rMesh) {
    MeshMetadata meta;
    meta.mNumPoints = rMesh.NumPoints();
    meta.mPointDim = rMesh.PointDim();

    meta.mCellBlocks.reserve(rMesh.NumCellBlocks());
    for (const auto block : rMesh.CellRange()) {
        CellBlockInfo info;
        info.mType = block.Type();
        info.mNumCells = block.NumCells();
        info.mRagged = block.IsRagged();
        info.mNodesPerCell = info.mRagged ? 0 : block.NodesPerCell();
        meta.mCellBlocks.push_back(std::move(info));
    }

    // The uniform API already guarantees sorted names, so no re-sort here --
    // re-sorting would be harmless but would imply the guarantee is in doubt.
    meta.mPointDataNames = rMesh.PointDataNames();
    meta.mCellDataNames = rMesh.CellDataNames();
    meta.mFieldDataNames = rMesh.FieldDataNames();

    // The mesh is already in memory, so this is nearly free -- every
    // region-capable reader builds regions alongside geometry, not lazily.
    meta.mRegions.reserve(rMesh.NumRegions());
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        RegionSummary summary;
        summary.mName = r.mName;
        summary.mKind = r.mKind;
        summary.mDim = r.mDim;
        summary.mTag = r.mTag;
        summary.mNumEntries = r.NumEntries();
        meta.mRegions.push_back(std::move(summary));
    }

    // The mesh is already in memory, so the bounding box is nearly free here --
    // unlike on a native metadata path, where it would force decoding the point
    // coordinates and defeat the whole point.
    const std::size_t dim = std::min<std::size_t>(meta.mPointDim, 3);
    if (meta.mNumPoints > 0 && dim > 0) {
        const NDArray& points = rMesh.Points();
        double lo[3] = {std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
        double hi[3] = {-std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity()};
        const std::size_t stride = detail::cols(points);
        for (std::size_t i = 0; i < meta.mNumPoints; ++i) {
            for (std::size_t d = 0; d < dim; ++d) {
                const double v = detail::read_double(points, i * stride + d);
                lo[d] = std::min(lo[d], v);
                hi[d] = std::max(hi[d], v);
            }
        }
        for (std::size_t d = 0; d < dim; ++d) {
            meta.mBBoxMin[d] = lo[d];
            meta.mBBoxMax[d] = hi[d];
        }
        meta.mHasBBox = true;
    }

    return meta;
}

}  // namespace meshioplusplus
