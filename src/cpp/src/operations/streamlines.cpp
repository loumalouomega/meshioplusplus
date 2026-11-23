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
//
// Streamlines of a point vector field as a mesh of line segments (see
// operations/streamlines.hpp). The tracing is detail/streamlines.cpp's; this
// file reads the field, calls it and builds the mesh.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/streamlines.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "../detail/streamlines.hpp"

namespace meshioplusplus {

namespace {

[[noreturn]] void sl_fail(const std::string& rMessage) {
    throw std::invalid_argument("meshio++: streamlines: " + rMessage);
}

}  // namespace

StreamlinesResult streamlines(const Mesh& rMesh, const StreamlinesOptions& rOptions) {
    if (rOptions.mArray.empty())
        sl_fail("an array name is required");
    if (!rMesh.HasPointData(rOptions.mArray))
        sl_fail(data_unknown_key_message(rMesh, DataLocation::Point, rOptions.mArray));
    if (rOptions.mSeeds < 0 || rOptions.mMaxSteps < 0)
        sl_fail("the seed and step counts must not be negative");
    if (!(rOptions.mMaxLength > 0.0) || rOptions.mStep < 0.0 || rOptions.mTerminalSpeed < 0.0 ||
        !(rOptions.mTolerance > 0.0))
        sl_fail(
            "the length and the tolerance must be positive, the step and the terminal speed "
            "not negative");
    if (rOptions.mSeedPoints.size() % 3 != 0)
        sl_fail("the seed points are 3 coordinates each");

    const NDArray& field = rMesh.PointData(rOptions.mArray);
    const std::size_t comps = data_num_components(field);
    if (comps != 2 && comps != 3)
        sl_fail("'" + rOptions.mArray + "' has " + std::to_string(comps) +
                " component(s); a vector field of 2 or 3 is needed");
    const std::size_t npoints = rMesh.NumPoints();
    if (detail::rows(field) != npoints)
        sl_fail("'" + rOptions.mArray + "' has " + std::to_string(detail::rows(field)) +
                " rows on a mesh of " + std::to_string(npoints) + " points");

    // Points and the field as 3 doubles each (a 2D point or field gets z = 0).
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    std::vector<double> xyz(npoints * 3, 0.0);
    std::vector<double> vec(npoints * 3, 0.0);
    for (std::size_t i = 0; i < npoints; ++i) {
        for (std::size_t c = 0; c < pdim && c < 3; ++c)
            xyz[i * 3 + c] = detail::read_double(points, i * pdim + c);
        for (std::size_t c = 0; c < comps; ++c)
            vec[i * 3 + c] = detail::read_double(field, i * comps + c);
    }

    detail::StreamlineTraceOptions to;
    to.mSeeds = static_cast<std::size_t>(rOptions.mSeeds);
    to.mSeedPoints = rOptions.mSeedPoints;
    to.mBackward = rOptions.mDirection != StreamlineDirection::Forward;
    to.mForward = rOptions.mDirection != StreamlineDirection::Backward;
    to.mAdaptive = rOptions.mIntegrator == StreamlineIntegrator::RK45;
    to.mStep = rOptions.mStep;
    to.mMaxLength = rOptions.mMaxLength;
    to.mMaxSteps = static_cast<std::size_t>(rOptions.mMaxSteps);
    to.mTerminalSpeed = rOptions.mTerminalSpeed;
    to.mTolerance = rOptions.mTolerance;
    const detail::StreamlineTraces traces =
        detail::trace_streamlines_timed(rMesh, xyz.data(), vec.data(), to);

    StreamlinesResult result;
    result.mDim = traces.mDim;
    result.mNumSeedsOutside = traces.mNumSeedsOutside;
    result.mNumLines = static_cast<std::int64_t>(traces.mLines.size());

    std::size_t total_points = 0;
    std::size_t total_segments = 0;
    for (const detail::StreamlineTrace& line : traces.mLines) {
        total_points += line.mTime.size();
        total_segments += line.mTime.size() - 1;
    }
    NDArray out_points = NDArray::Uninit(DType::Float64, {total_points, std::size_t{3}});
    NDArray out_time = NDArray::Uninit(DType::Float64, {total_points});
    NDArray out_speed = NDArray::Uninit(DType::Float64, {total_points});
    NDArray out_seed = NDArray::Uninit(DType::Int64, {total_points});
    NDArray out_field = NDArray::Uninit(DType::Float64, {total_points, comps});
    NDArray conn = NDArray::Uninit(DType::Int64, {total_segments, std::size_t{2}});
    NDArray seg_seed = NDArray::Uninit(DType::Int64, {total_segments});
    double* pp = out_points.As<double>();
    double* pt = out_time.As<double>();
    double* ps = out_speed.As<double>();
    std::int64_t* pi = out_seed.As<std::int64_t>();
    double* pf = out_field.As<double>();
    std::int64_t* pc = conn.As<std::int64_t>();
    std::int64_t* pg = seg_seed.As<std::int64_t>();
    std::size_t base = 0;
    std::size_t seg = 0;
    for (const detail::StreamlineTrace& line : traces.mLines) {
        const std::size_t n = line.mTime.size();
        for (std::size_t k = 0; k < n; ++k) {
            for (std::size_t c = 0; c < 3; ++c)
                pp[(base + k) * 3 + c] = line.mXyz[k * 3 + c];
            pt[base + k] = line.mTime[k];
            const double* v = line.mVelocity.data() + k * 3;
            ps[base + k] = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            pi[base + k] = line.mSeedId;
            for (std::size_t c = 0; c < comps; ++c)
                pf[(base + k) * comps + c] = v[c];
            if (k + 1 < n) {
                pc[seg * 2] = static_cast<std::int64_t>(base + k);
                pc[seg * 2 + 1] = static_cast<std::int64_t>(base + k + 1);
                pg[seg] = line.mSeedId;
                ++seg;
            }
        }
        base += n;
        result.mSeedIds.push_back(line.mSeedId);
        result.mEndBackward.push_back(static_cast<std::int32_t>(line.mBackwardEnd));
        result.mEndForward.push_back(static_cast<std::int32_t>(line.mForwardEnd));
    }
    Mesh& out = result.mMesh;
    out.AssignPoints(std::move(out_points));
    if (total_segments > 0) {
        out.AddCellBlock("line", std::move(conn));
        std::vector<NDArray> seg_blocks;
        seg_blocks.push_back(std::move(seg_seed));
        out.AddCellData(kStreamlineSeedIdName, std::move(seg_blocks));
    }
    out.AddPointData(kStreamlineTimeName, std::move(out_time));
    out.AddPointData(kStreamlineSpeedName, std::move(out_speed));
    out.AddPointData(kStreamlineSeedIdName, std::move(out_seed));
    out.AddPointData(rOptions.mArray, std::move(out_field));
    return result;
}

}  // namespace meshioplusplus
