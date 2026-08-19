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
// The Hessian (second derivative) of a scalar point_data field. Deliberately
// a COMPOSITION of two gradient() calls, not a new numerical kernel; see
// operations/hessian.hpp for the full contract and the exactness argument.
//
// The two internal gradient() results are private working state under fixed
// internal names -- never returned to the caller -- so a name collision with
// a real array in the input mesh cannot lose data: the actual result mesh is
// built separately from a fresh clone of the ORIGINAL input, exactly
// error.cpp's own isolation discipline.

// System includes
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/hessian.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/data_manage.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kHessPrefix = "meshio++: hessian: ";
// Private working names: never reach the returned mesh (see the file
// comment above), so any clash with a user array only affects internal
// intermediate copies, not the caller's data.
constexpr const char* kHessRawGradName = "__hessian_raw_gradient__";
constexpr const char* kHessRawName = "__hessian_raw__";

}  // namespace

HessianResult hessian(const Mesh& rMesh, const HessianOptions& rOptions) {
    // --- validation, all before any work -------------------------------------
    if (rOptions.mArrayName.empty())
        throw std::invalid_argument(std::string(kHessPrefix) + "an array name is required");
    if (rOptions.mLocation == DataLocation::Field)
        throw std::invalid_argument(std::string(kHessPrefix) +
                                    "field_data has no location on the mesh; the result must be "
                                    "'point' or 'cell'");
    if (!rMesh.HasPointData(rOptions.mArrayName)) {
        // Same restriction gradient states, for the same reason.
        if (rMesh.HasCellData(rOptions.mArrayName))
            throw std::invalid_argument(
                std::string(kHessPrefix) + "'" + rOptions.mArrayName +
                "' is cell_data, which is piecewise constant and has no derivative; convert it "
                "first with cell_data_to_point_data (CLI: meshioplusplus data to-point)");
        throw std::invalid_argument(
            std::string(kHessPrefix) +
            data_unknown_key_message(rMesh, DataLocation::Point, rOptions.mArrayName));
    }
    const NDArray& field = rMesh.PointData(rOptions.mArrayName);
    const std::size_t field_comp = data_num_components(field);
    if (field_comp != 1)
        throw std::invalid_argument(
            std::string(kHessPrefix) + "'" + rOptions.mArrayName + "' has " +
            std::to_string(field_comp) +
            " components; hessian currently supports scalar fields only -- call it once per "
            "component of a vector field");

    const std::string out_name =
        rOptions.mOutputName.empty() ? rOptions.mArrayName + kHessianSuffix : rOptions.mOutputName;
    if (!rOptions.mOverwrite) {
        const bool taken = rOptions.mLocation == DataLocation::Point ? rMesh.HasPointData(out_name)
                                                                     : rMesh.HasCellData(out_name);
        if (taken)
            throw std::invalid_argument(
                std::string(kHessPrefix) + "'" + out_name + "' already exists in " +
                data_location_name(rOptions.mLocation) + " (pass overwrite=true to replace it)");
    }

    // --- two composed gradient() passes -----------------------------------
    // The first pass MUST be Point-located: the second pass' own validation
    // requires a point_data input, which is why this can't stay at gradient's
    // Cell default the way estimate_error's own single pass does.
    GradientOptions grad1;
    grad1.mArrayName = rOptions.mArrayName;
    grad1.mLocation = DataLocation::Point;
    grad1.mMethod = rOptions.mMethod;
    grad1.mOutputName = kHessRawGradName;
    grad1.mOverwrite = true;
    const GradientResult r1 = gradient(rMesh, grad1);

    // The second pass differentiates the first pass' (n,3) gradient with the
    // default Gradient operator, which is generic over the input's own
    // component count -- feeding a 3-component field back in yields (n,9),
    // the flattened row-major 3x3 Hessian. Always computed at Cell first,
    // mirroring gradient()'s own unconditional cell-first step.
    GradientOptions grad2;
    grad2.mArrayName = kHessRawGradName;
    grad2.mLocation = DataLocation::Cell;
    grad2.mMethod = rOptions.mMethod;
    grad2.mOutputName = kHessRawName;
    grad2.mOverwrite = true;
    const GradientResult r2 = gradient(r1.mMesh, grad2);

    // --- build the actual result from a fresh clone of the ORIGINAL input ---
    Mesh out = detail::clone_mesh(rMesh);
    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<NDArray> hess_blocks;
    hess_blocks.reserve(nblocks);
    for (std::size_t b = 0; b < nblocks; ++b)
        hess_blocks.push_back(r2.mMesh.CellData(kHessRawName, b));  // deep copy: NDArray owns its buffer
    out.AddCellData(out_name, std::move(hess_blocks));

    if (rOptions.mLocation == DataLocation::Point) {
        // Compose the existing averaging rather than reimplementing a scatter
        // -- gradient()'s own Point-location handling, applied to the final
        // Hessian array.
        DataAverageOptions avg;
        avg.names = {out_name};
        avg.weight = CellPointWeight::Uniform;
        avg.overwrite = true;
        Mesh with_points = cell_data_to_point_data(out, avg);
        out = data_drop(with_points, DataLocation::Cell, {out_name});
    }

    return HessianResult{std::move(out), r2.mNumSkipped, r1.mNumFallback + r2.mNumFallback};
}

}  // namespace meshioplusplus
