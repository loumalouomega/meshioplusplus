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
#pragma once

/**
 * @file operations/tensor_invariants.hpp
 * @brief `tensor_invariants`: von Mises, principal values, and the hydrostatic
 * / deviatoric split of a symmetric-tensor data array.
 *
 * Every format that stores a stress, strain or other rank-2 tensor field
 * needs these; before this operation the only route was the CalculiX `.frd`
 * reader's format-specific `derived=` option (`formats/frd.hpp`), which
 * neither the generic `read`, the CLIs, MCP nor the flat bindings carried.
 * `.frd`'s `derived=` now calls this operation internally and keeps its
 * existing output names; see `doc/formats/frd.md`.
 *
 * ### Input layout
 *
 * An array is treated as one tensor per row, selected by its trailing
 * component count:
 *
 *  - **6 components**: a symmetric tensor `xx yy zz xy yz zx`, the order used
 *    throughout this codebase (`doc/mesh_data_model.md`).
 *  - **9 components**: a general 3x3 tensor, row-major (`xx xy xz yx yy yz zx
 *    zy zz`), the layout `operations/gradient.hpp` and `operations/hessian.hpp`
 *    produce. `Mises` and `Principal` are only defined for a symmetric tensor,
 *    so both are computed from the symmetric part `0.5*(T + T^T)`.
 *
 * Any other component count is a `std::invalid_argument`.
 *
 * ### Outputs
 *
 * Requested through the `mOutputs` bitmask, each written as `rOpts.prefix +
 * name + suffix`:
 *
 *  - `Mises`: `sqrt(0.5*((xx-yy)^2+(yy-zz)^2+(zz-xx)^2+6*(xy^2+yz^2+zx^2)))`,
 *    using the symmetric part for a 9-component input. Shape `(n,)`.
 *  - `Principal`: eigenvalues of the symmetric part, ascending. Shape `(n, 3)`.
 *  - `Hydrostatic`: `(xx+yy+zz)/3`, the mean of the tensor's own diagonal
 *    (unsymmetrized). Shape `(n,)`.
 *  - `Deviatoric`: the input with `Hydrostatic` subtracted from its three
 *    diagonal entries, off-diagonal entries unchanged; same shape and
 *    component count as the input.
 *
 * A row with any non-finite input component produces NaN in every output for
 * that row rather than throwing (mirrors `.frd`'s `derived=` and the general
 * "non-finite in, non-finite out" rule of `operations/data_common.hpp`, except
 * that here there is no reduction to exclude the value from).
 *
 * Points, connectivity and block order are left bit-identical — this is a
 * data-only operation and reuses `detail/clone_mesh`.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus {

/// Which invariants to compute; combine with `|`.
enum class TensorInvariant : unsigned {
    Mises = 1u << 0,
    Principal = 1u << 1,
    Hydrostatic = 1u << 2,
    Deviatoric = 1u << 3,
    All = Mises | Principal | Hydrostatic | Deviatoric,
};

/// `a | b` for `TensorInvariant` flags.
inline TensorInvariant operator|(TensorInvariant a, TensorInvariant b) {
    return static_cast<TensorInvariant>(static_cast<unsigned>(a) | static_cast<unsigned>(b));
}

/// `a & b` for `TensorInvariant` flags, as a plain bitmask (not a `TensorInvariant`
/// itself, so it can be tested with `!= 0` without another cast).
inline unsigned operator&(TensorInvariant a, TensorInvariant b) {
    return static_cast<unsigned>(a) & static_cast<unsigned>(b);
}

/// Options for `tensor_invariants`.
struct TensorInvariantsOptions {
    /// Which data map the array lives in. `Field` is rejected (no per-row
    /// tensor to reduce).
    DataLocation location = DataLocation::Point;
    /// Names to process; empty means every 6- or 9-component array at
    /// `location` (arrays with another component count are skipped, not an
    /// error, when the list is empty).
    std::vector<std::string> names;
    /// Which invariants to compute.
    TensorInvariant outputs = TensorInvariant::All;
    /// Prefix prepended to every output array's name.
    std::string prefix;
    /// Suffix appended after the invariant's own name segment
    /// (`prefix + name + "_mises" + suffix`, etc.).
    std::string suffix;
    /// Throw instead of silently overwriting an existing array of the target name.
    bool overwrite = true;
};

/**
 * @brief Computes von Mises, principal, hydrostatic and/or deviatoric fields
 * for the selected tensor arrays.
 * @param rMesh the source mesh (unmodified).
 * @param rOpts which arrays, which outputs, and how to name them.
 * @return a new mesh with the source data plus the requested invariant arrays;
 *         geometry is bit-identical to @p rMesh.
 * @throws std::invalid_argument on an unknown name, `location ==
 *         DataLocation::Field`, an array whose trailing component count is
 *         neither 6 nor 9 when named explicitly, a `cell_data` array whose
 *         block count disagrees with the mesh, or a name collision when
 *         `overwrite` is `false`.
 */
MESHIOPLUSPLUS_API Mesh tensor_invariants(const Mesh& rMesh, const TensorInvariantsOptions& rOpts);

/**
 * @brief Parses a comma-separated list of invariant names into a bitmask.
 * @param rName one or more of `mises`, `principal`, `hydrostatic`,
 *        `deviatoric`, `all`, separated by commas.
 * @return the matching `TensorInvariant` bitmask.
 * @throws std::invalid_argument on an unknown name or an empty list.
 */
MESHIOPLUSPLUS_API TensorInvariant tensor_invariant_from_name(const std::string& rName);

}  // namespace meshioplusplus
