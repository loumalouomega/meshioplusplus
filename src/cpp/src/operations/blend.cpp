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

// blend_steps. See operations/blend.hpp for the contract.

// System includes
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/blend.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kBlPrefix = "meshio++: blend_steps: ";

bool bl_is_float(DType Dt) {
    return Dt == DType::Float32 || Dt == DType::Float64;
}

// (1 - w) a + w b over two arrays of one shape; integers from the nearer one.
NDArray bl_blend(const NDArray& rA, const NDArray& rB, double W, const std::string& rWhat) {
    if (rA.Shape() != rB.Shape())
        throw std::invalid_argument(std::string(kBlPrefix) + rWhat +
                                    " has different shapes in the two steps");
    if (!bl_is_float(rA.Dtype()) || !bl_is_float(rB.Dtype()))
        return detail::data_owned_copy(W < 0.5 ? rA : rB);
    NDArray out = NDArray::Uninit(rA.Dtype(), rA.Shape());
    const double wa = 1.0 - W;
    parallel_for_bw(rA.Size(), [&](std::size_t i) {
        detail::write_double(out, i,
                             wa * detail::read_double(rA, i) + W * detail::read_double(rB, i));
    });
    return out;
}

void bl_check_topology(const Mesh& rA, const Mesh& rB) {
    if (rA.NumPoints() != rB.NumPoints() || rA.PointDim() != rB.PointDim())
        throw std::invalid_argument(std::string(kBlPrefix) + "the steps have " +
                                    std::to_string(rA.NumPoints()) + " and " +
                                    std::to_string(rB.NumPoints()) +
                                    " points; blending needs one topology across the steps");
    if (rA.NumCellBlocks() != rB.NumCellBlocks())
        throw std::invalid_argument(std::string(kBlPrefix) +
                                    "the steps have different numbers of cell blocks");
    for (std::size_t b = 0; b < rA.NumCellBlocks(); ++b) {
        const auto ca = rA.Cells(b);
        const auto cb = rB.Cells(b);
        if (std::string(ca.Type()) != std::string(cb.Type()) || ca.NumCells() != cb.NumCells())
            throw std::invalid_argument(std::string(kBlPrefix) + "cell block " + std::to_string(b) +
                                        " differs between the steps (" + std::string(ca.Type()) +
                                        " x" + std::to_string(ca.NumCells()) + " vs " +
                                        std::string(cb.Type()) + " x" +
                                        std::to_string(cb.NumCells()) + ")");
    }
    auto same_names = [&](const std::vector<std::string>& rX, const std::vector<std::string>& rY,
                          const char* pWhat) {
        if (rX != rY)
            throw std::invalid_argument(std::string(kBlPrefix) + "the steps carry different " +
                                        pWhat + " arrays");
    };
    same_names(rA.PointDataNames(), rB.PointDataNames(), "point_data");
    same_names(rA.CellDataNames(), rB.CellDataNames(), "cell_data");
    same_names(rA.FieldDataNames(), rB.FieldDataNames(), "field_data");
}

}  // namespace

Mesh blend_steps(const Mesh& rA, const Mesh& rB, double W, const BlendOptions& rOptions) {
    bl_check_topology(rA, rB);
    Mesh out = detail::clone_geometry(rA);
    if (rOptions.mBlendPoints)
        out.AssignPoints(bl_blend(rA.Points(), rB.Points(), W, "the point array"));
    for (const std::string& name : rA.PointDataNames())
        out.AddPointData(
            name, bl_blend(rA.PointData(name), rB.PointData(name), W, "point_data '" + name + "'"));
    for (const std::string& name : rA.CellDataNames()) {
        if (rA.CellDataNumBlocks(name) != rB.CellDataNumBlocks(name))
            throw std::invalid_argument(std::string(kBlPrefix) + "cell_data '" + name +
                                        "' has different block counts in the two steps");
        std::vector<NDArray> blocks;
        for (std::size_t b = 0; b < rA.CellDataNumBlocks(name); ++b)
            blocks.push_back(bl_blend(rA.CellData(name, b), rB.CellData(name, b), W,
                                      "cell_data '" + name + "'"));
        out.AddCellData(name, std::move(blocks));
    }
    for (const std::string& name : rA.FieldDataNames())
        out.AddFieldData(
            name, bl_blend(rA.FieldData(name), rB.FieldData(name), W, "field_data '" + name + "'"));
    return out;
}

}  // namespace meshioplusplus
