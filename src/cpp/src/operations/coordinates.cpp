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
// Cylindrical and spherical components about an axis. See
// operations/coordinates.hpp for the contract.

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
#include "meshioplusplus/operations/coordinates.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/node_adjacency.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {

namespace {

using CoordVec = std::array<double, 3>;

CoordVec coord_unit(const CoordVec& rA) {
    const double n = std::sqrt(detail::vec3_dot(rA, rA));
    return {rA[0] / n, rA[1] / n, rA[2] / n};
}

// The frame every point shares: the unit axis, and the two directions of angle
// 0 and pi/2 about it.
struct CoordFrame {
    CoordVec mAxis, mE1, mE2, mOrigin;
};

// The coordinates of position @p rP, and the local basis (rows) there; true
// when the radial direction was undefined and e1 stood in.
bool coord_local(const CoordFrame& rF, CoordinateSystem System, const CoordVec& rP,
                 CoordVec& rCoords, std::array<CoordVec, 3>& rBasis) {
    const CoordVec d{rP[0] - rF.mOrigin[0], rP[1] - rF.mOrigin[1], rP[2] - rF.mOrigin[2]};
    const double z = detail::vec3_dot(d, rF.mAxis);
    const CoordVec rho{d[0] - z * rF.mAxis[0], d[1] - z * rF.mAxis[1], d[2] - z * rF.mAxis[2]};
    const double r_cyl = std::sqrt(detail::vec3_dot(rho, rho));
    const double scale = std::sqrt(detail::vec3_dot(d, d));
    const bool on_axis = !(r_cyl > 1e-14 * std::max(scale, 1.0));
    const CoordVec e_rho = on_axis ? rF.mE1 : coord_unit(rho);
    const double angle =
        on_axis ? 0.0 : std::atan2(detail::vec3_dot(rho, rF.mE2), detail::vec3_dot(rho, rF.mE1));
    const CoordVec e_phi = detail::vec3_cross(rF.mAxis, e_rho);
    if (System == CoordinateSystem::Cylindrical) {
        rCoords = {r_cyl, angle, z};
        rBasis = {e_rho, e_phi, rF.mAxis};
        return on_axis;
    }
    const double r = scale;
    const bool at_centre = !(r > 0.0);
    const double polar = at_centre ? 0.0 : std::atan2(r_cyl, z);
    const double c = std::cos(polar), s = std::sin(polar);
    const CoordVec e_r{s * e_rho[0] + c * rF.mAxis[0], s * e_rho[1] + c * rF.mAxis[1],
                       s * e_rho[2] + c * rF.mAxis[2]};
    const CoordVec e_theta{c * e_rho[0] - s * rF.mAxis[0], c * e_rho[1] - s * rF.mAxis[1],
                           c * e_rho[2] - s * rF.mAxis[2]};
    rCoords = {r, polar, angle};
    rBasis = {e_r, e_theta, e_phi};
    return on_axis || at_centre;
}

// v' = R v (forward) or R^T v (inverse), R's rows the basis.
void coord_vector(const std::array<CoordVec, 3>& rB, bool Inverse, const double* pIn,
                  double* pOut) {
    for (std::size_t i = 0; i < 3; ++i) {
        double acc = 0.0;
        for (std::size_t k = 0; k < 3; ++k)
            acc += (Inverse ? rB[k][i] : rB[i][k]) * pIn[k];
        pOut[i] = acc;
    }
}

// T' = R T R^T (forward) or R^T T R (inverse), row-major 3x3.
void coord_tensor(const std::array<CoordVec, 3>& rB, bool Inverse, const double* pIn,
                  double* pOut) {
    const auto m = [&](std::size_t i, std::size_t k) { return Inverse ? rB[k][i] : rB[i][k]; };
    double tmp[9];
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            double acc = 0.0;
            for (std::size_t k = 0; k < 3; ++k)
                acc += m(i, k) * pIn[k * 3 + j];
            tmp[i * 3 + j] = acc;
        }
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            double acc = 0.0;
            for (std::size_t k = 0; k < 3; ++k)
                acc += tmp[i * 3 + k] * m(j, k);
            pOut[i * 3 + j] = acc;
        }
}

// Rewrite one array of @p Rows rows at @p rPositions; returns the new array.
NDArray coord_convert_array(const NDArray& rIn, const std::string& rName,
                            const std::vector<CoordVec>& rPositions, const CoordFrame& rF,
                            const CoordinateOptions& rOptions) {
    const std::vector<std::size_t>& shape = rIn.Shape();
    const std::size_t rows = rPositions.size();
    std::size_t width = 0;
    if (shape.size() == 2 && shape[0] == rows && (shape[1] == 3 || shape[1] == 9))
        width = shape[1];
    else if (shape.size() == 3 && shape[0] == rows && shape[1] == 3 && shape[2] == 3)
        width = 9;
    if (width == 0)
        throw std::invalid_argument("convert_coordinates: '" + rName +
                                    "' is neither a vector (n, 3) nor a tensor (n, 9) or "
                                    "(n, 3, 3) array");
    NDArray out(DType::Float64, shape);
    double* po = out.As<double>();
    parallel_for(rows, [&](std::size_t i) {
        CoordVec coords;
        std::array<CoordVec, 3> basis;
        coord_local(rF, rOptions.mSystem, rPositions[i], coords, basis);
        double in[9];
        for (std::size_t k = 0; k < width; ++k)
            in[k] = detail::read_double(rIn, i * width + k);
        if (width == 3)
            coord_vector(basis, rOptions.mInverse, in, po + i * 3);
        else
            coord_tensor(basis, rOptions.mInverse, in, po + i * 9);
    });
    return out;
}

}  // namespace

CoordinateSystem coordinate_system_from_name(const std::string& rName) {
    if (rName == "cylindrical")
        return CoordinateSystem::Cylindrical;
    if (rName == "spherical")
        return CoordinateSystem::Spherical;
    throw std::invalid_argument("convert_coordinates: unknown system '" + rName +
                                "' (expected 'cylindrical' or 'spherical')");
}

CoordinateResult convert_coordinates(const Mesh& rMesh, const CoordinateOptions& rOptions) {
    CoordFrame frame;
    const double axis_len = std::sqrt(detail::vec3_dot(rOptions.mAxis, rOptions.mAxis));
    if (!(axis_len > 0.0) || !std::isfinite(axis_len))
        throw std::invalid_argument("convert_coordinates: the axis must be a nonzero vector");
    frame.mAxis = coord_unit(rOptions.mAxis);
    frame.mOrigin = rOptions.mOrigin;
    // e1: the reference projected off the axis, or x / y when it lies along it.
    for (const CoordVec& candidate :
         {rOptions.mReference, CoordVec{1.0, 0.0, 0.0}, CoordVec{0.0, 1.0, 0.0}}) {
        const double along = detail::vec3_dot(candidate, frame.mAxis);
        const CoordVec off{candidate[0] - along * frame.mAxis[0],
                           candidate[1] - along * frame.mAxis[1],
                           candidate[2] - along * frame.mAxis[2]};
        const double len = std::sqrt(detail::vec3_dot(off, off));
        const double clen = std::sqrt(detail::vec3_dot(candidate, candidate));
        if (len > 1e-9 * clen) {
            frame.mE1 = coord_unit(off);
            break;
        }
    }
    frame.mE2 = detail::vec3_cross(frame.mAxis, frame.mE1);

    // Positions: the points, and the cell centroids (distinct nodes).
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    std::vector<CoordVec> points(n, CoordVec{0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < dim && k < 3; ++k)
            points[i][k] = detail::read_double(rMesh.Points(), i * dim + k);

    CoordinateResult result;
    result.mMesh = detail::clone_mesh(rMesh);
    if (rOptions.mRecordCoordinates) {
        NDArray coords(DType::Float64, {n, std::size_t{3}});
        double* pc = coords.As<double>();
        std::int64_t singular = 0;
        for (std::size_t i = 0; i < n; ++i) {
            CoordVec c;
            std::array<CoordVec, 3> basis;
            if (coord_local(frame, rOptions.mSystem, points[i], c, basis))
                ++singular;
            pc[i * 3 + 0] = c[0];
            pc[i * 3 + 1] = c[1];
            pc[i * 3 + 2] = c[2];
        }
        result.mNumSingular = singular;
        result.mMesh.AddPointData(rOptions.mSystem == CoordinateSystem::Cylindrical
                                      ? kCylindricalCoordsName
                                      : kSphericalCoordsName,
                                  std::move(coords));
    }

    std::vector<CoordVec> centroids;
    bool have_centroids = false;
    const auto cell_positions = [&]() -> const std::vector<CoordVec>& {
        if (!have_centroids) {
            std::vector<std::int64_t> ids;
            for (const auto cb : rMesh.CellRange())
                for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                    detail::cell_node_ids(cb, c, n, ids);
                    std::sort(ids.begin(), ids.end());
                    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
                    CoordVec m{0.0, 0.0, 0.0};
                    for (const std::int64_t id : ids)
                        for (std::size_t k = 0; k < 3; ++k)
                            m[k] += points[static_cast<std::size_t>(id)][k];
                    if (!ids.empty())
                        for (std::size_t k = 0; k < 3; ++k)
                            m[k] /= static_cast<double>(ids.size());
                    centroids.push_back(m);
                }
            have_centroids = true;
        }
        return centroids;
    };

    for (const std::string& name : rOptions.mArrays) {
        bool found = false;
        if (rMesh.HasPointData(name)) {
            found = true;
            result.mMesh.AddPointData(
                name, coord_convert_array(rMesh.PointData(name), name, points, frame, rOptions));
            ++result.mNumArrays;
        }
        if (rMesh.HasCellData(name)) {
            found = true;
            const std::vector<CoordVec>& positions = cell_positions();
            std::vector<NDArray> blocks;
            std::size_t base = 0;
            std::size_t block = 0;
            for (const auto cb : rMesh.CellRange()) {
                const std::size_t count = cb.NumCells();
                const std::vector<CoordVec> slice(positions.begin() + base,
                                                  positions.begin() + base + count);
                blocks.push_back(
                    coord_convert_array(rMesh.CellData(name, block), name, slice, frame, rOptions));
                base += count;
                ++block;
            }
            result.mMesh.AddCellData(name, std::move(blocks));
            ++result.mNumArrays;
        }
        if (!found)
            throw std::invalid_argument("convert_coordinates: no point or cell data array named '" +
                                        name + "'");
    }
    return result;
}

}  // namespace meshioplusplus
