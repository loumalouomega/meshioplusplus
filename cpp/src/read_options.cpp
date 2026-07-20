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

// Project includes
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {

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
