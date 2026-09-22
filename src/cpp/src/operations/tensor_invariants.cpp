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
// von Mises / principal / hydrostatic / deviatoric invariants of a symmetric
// (6-component) or general (9-component row-major) tensor data array. See
// operations/tensor_invariants.hpp for the contract; the eigensolver lives in
// detail/sym3_eigen.hpp, shared with the CalculiX .frd reader's derived= option.

// System includes
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/tensor_invariants.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/sym3_eigen.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

/// Reads row @p r (0-based) of a 6- or 9-component array into the
/// `xx yy zz xy yz zx` order `detail::sym3_*` expects, symmetrizing the
/// off-diagonal terms for a 9-component input. Returns false (sym left
/// untouched) when any source component is non-finite.
bool tinv_read_symmetric(const NDArray& rArray, std::size_t Ncomp, std::size_t r, double* pSym) {
    double t[9];
    for (std::size_t k = 0; k < Ncomp; ++k)
        t[k] = detail::read_double(rArray, r * Ncomp + k);
    for (std::size_t k = 0; k < Ncomp; ++k)
        if (!std::isfinite(t[k]))
            return false;
    if (Ncomp == 6) {
        for (int k = 0; k < 6; ++k)
            pSym[k] = t[k];
        return true;
    }
    // 9-component row-major: xx xy xz yx yy yz zx zy zz.
    pSym[0] = t[0];
    pSym[1] = t[4];
    pSym[2] = t[8];
    pSym[3] = 0.5 * (t[1] + t[3]);
    pSym[4] = 0.5 * (t[5] + t[7]);
    pSym[5] = 0.5 * (t[2] + t[6]);
    return true;
}

/// The (unsymmetrized) hydrostatic mean of the tensor's own diagonal.
double tinv_hydrostatic(const NDArray& rArray, std::size_t Ncomp, std::size_t r) {
    const double xx = detail::read_double(rArray, r * Ncomp + 0);
    const double yy = detail::read_double(rArray, r * Ncomp + (Ncomp == 6 ? 1 : 4));
    const double zz = detail::read_double(rArray, r * Ncomp + (Ncomp == 6 ? 2 : 8));
    if (!std::isfinite(xx) || !std::isfinite(yy) || !std::isfinite(zz))
        return std::numeric_limits<double>::quiet_NaN();
    return (xx + yy + zz) / 3.0;
}

/// Number of rows of a Ncomp-wide array.
std::size_t tinv_nrows(const NDArray& rArray, std::size_t Ncomp) {
    return Ncomp > 0 ? rArray.Size() / Ncomp : 0;
}

NDArray tinv_mises(const NDArray& rSrc, std::size_t Ncomp) {
    const std::size_t nrows = tinv_nrows(rSrc, Ncomp);
    NDArray dst(DType::Float64, {nrows});
    double* out = dst.As<double>();
    parallel_for(nrows, [&](std::size_t r) {
        double sym[6];
        out[r] = tinv_read_symmetric(rSrc, Ncomp, r, sym)
                     ? detail::sym3_mises(sym)
                     : std::numeric_limits<double>::quiet_NaN();
    });
    return dst;
}

NDArray tinv_principal(const NDArray& rSrc, std::size_t Ncomp) {
    const std::size_t nrows = tinv_nrows(rSrc, Ncomp);
    NDArray dst(DType::Float64, {nrows, std::size_t{3}});
    double* out = dst.As<double>();
    parallel_for(nrows, [&](std::size_t r) {
        double sym[6];
        if (tinv_read_symmetric(rSrc, Ncomp, r, sym))
            detail::sym3_principal(sym, out + r * 3);
        else
            for (int k = 0; k < 3; ++k)
                out[r * 3 + k] = std::numeric_limits<double>::quiet_NaN();
    });
    return dst;
}

NDArray tinv_hydrostatic_array(const NDArray& rSrc, std::size_t Ncomp) {
    const std::size_t nrows = tinv_nrows(rSrc, Ncomp);
    NDArray dst(DType::Float64, {nrows});
    double* out = dst.As<double>();
    parallel_for(nrows, [&](std::size_t r) { out[r] = tinv_hydrostatic(rSrc, Ncomp, r); });
    return dst;
}

NDArray tinv_deviatoric(const NDArray& rSrc, std::size_t Ncomp) {
    const std::size_t nrows = tinv_nrows(rSrc, Ncomp);
    NDArray dst = NDArray::Uninit(DType::Float64, rSrc.Shape());
    double* out = dst.As<double>();
    const std::size_t yy_idx = Ncomp == 6 ? 1 : 4;
    const std::size_t zz_idx = Ncomp == 6 ? 2 : 8;
    parallel_for(nrows, [&](std::size_t r) {
        const double h = tinv_hydrostatic(rSrc, Ncomp, r);
        for (std::size_t k = 0; k < Ncomp; ++k) {
            const double v = detail::read_double(rSrc, r * Ncomp + k);
            const bool diag = (k == 0 || k == yy_idx || k == zz_idx);
            out[r * Ncomp + k] = (diag && std::isfinite(v) && std::isfinite(h)) ? v - h : v;
        }
    });
    return dst;
}

/// Whether an array's trailing component count is one `tensor_invariants` can
/// process.
bool tinv_is_tensor_shaped(const NDArray& rArray) {
    const std::size_t nc = data_num_components(rArray);
    return nc == 6 || nc == 9;
}

}  // namespace

TensorInvariant tensor_invariant_from_name(const std::string& rName) {
    unsigned mask = 0;
    std::size_t start = 0;
    bool any = false;
    while (start <= rName.size()) {
        const std::size_t comma = rName.find(',', start);
        const std::string tok = rName.substr(start, comma == std::string::npos ? std::string::npos
                                                                                : comma - start);
        if (tok == "mises")
            mask |= static_cast<unsigned>(TensorInvariant::Mises);
        else if (tok == "principal")
            mask |= static_cast<unsigned>(TensorInvariant::Principal);
        else if (tok == "hydrostatic")
            mask |= static_cast<unsigned>(TensorInvariant::Hydrostatic);
        else if (tok == "deviatoric")
            mask |= static_cast<unsigned>(TensorInvariant::Deviatoric);
        else if (tok == "all")
            mask |= static_cast<unsigned>(TensorInvariant::All);
        else
            throw std::invalid_argument(
                "meshio++: unknown tensor invariant '" + tok +
                "' (expected 'mises', 'principal', 'hydrostatic', 'deviatoric' or 'all')");
        any = true;
        if (comma == std::string::npos)
            break;
        start = comma + 1;
    }
    if (!any)
        throw std::invalid_argument("meshio++: tensor_invariant_from_name: empty list");
    return static_cast<TensorInvariant>(mask);
}

Mesh tensor_invariants(const Mesh& rMesh, const TensorInvariantsOptions& rOpts) {
    if (rOpts.location == DataLocation::Field)
        throw std::invalid_argument(
            "meshio++: tensor_invariants: field_data has no per-row tensor to reduce");

    std::vector<std::string> names = rOpts.names;
    if (names.empty()) {
        for (const std::string& n : data_names(rMesh, rOpts.location)) {
            const NDArray& src = rOpts.location == DataLocation::Point ? rMesh.PointData(n)
                                                                        : rMesh.CellData(n, 0);
            if (rOpts.location == DataLocation::Cell && rMesh.CellDataNumBlocks(n) == 0)
                continue;
            if (tinv_is_tensor_shaped(src))
                names.push_back(n);
        }
    } else {
        for (const std::string& n : names) {
            if (!data_has(rMesh, rOpts.location, n))
                throw std::invalid_argument(data_unknown_key_message(rMesh, rOpts.location, n));
            const NDArray& src = rOpts.location == DataLocation::Point ? rMesh.PointData(n)
                                                                        : rMesh.CellData(n, 0);
            if (!tinv_is_tensor_shaped(src))
                throw std::invalid_argument(
                    "meshio++: tensor_invariants: '" + n + "' has " +
                    std::to_string(data_num_components(src)) +
                    " component(s); expected 6 (symmetric) or 9 (general 3x3)");
        }
    }

    // Only ever called for the Point branch below (Field is rejected up front,
    // Cell has its own per-block helper), so it only needs AddPointData.
    auto add_point = [&](Mesh& rOut, const std::string& rTarget, NDArray Value) {
        if (!rOpts.overwrite && data_has(rOut, DataLocation::Point, rTarget))
            throw std::invalid_argument("meshio++: tensor_invariants: '" + rTarget +
                                        "' already exists (overwrite=false)");
        rOut.AddPointData(rTarget, std::move(Value));
    };

    Mesh out = detail::clone_mesh(rMesh);

    for (const std::string& name : names) {
        const std::string base = rOpts.prefix + name;

        if (rOpts.location == DataLocation::Cell) {
            const std::size_t nblocks = rMesh.NumCellBlocks();
            if (rMesh.CellDataNumBlocks(name) != nblocks)
                throw std::invalid_argument(
                    "meshio++: tensor_invariants: cell_data '" + name + "' has " +
                    std::to_string(rMesh.CellDataNumBlocks(name)) + " block(s) but the mesh has " +
                    std::to_string(nblocks) + " cell block(s)");
            if (nblocks == 0)
                continue;
            const std::size_t ncomp = data_num_components(rMesh.CellData(name, 0));

            auto add_cell = [&](const std::string& rSuffix,
                                NDArray (*pFn)(const NDArray&, std::size_t)) {
                const std::string target = base + rSuffix + rOpts.suffix;
                if (!rOpts.overwrite && data_has(out, DataLocation::Cell, target))
                    throw std::invalid_argument("meshio++: tensor_invariants: '" + target +
                                                "' already exists (overwrite=false)");
                std::vector<NDArray> blocks;
                blocks.reserve(nblocks);
                for (std::size_t b = 0; b < nblocks; ++b)
                    blocks.push_back(pFn(rMesh.CellData(name, b), ncomp));
                out.AddCellData(target, std::move(blocks));
            };
            if ((rOpts.outputs & TensorInvariant::Mises) != 0)
                add_cell("_mises", tinv_mises);
            if ((rOpts.outputs & TensorInvariant::Principal) != 0)
                add_cell("_principal", tinv_principal);
            if ((rOpts.outputs & TensorInvariant::Hydrostatic) != 0)
                add_cell("_hydrostatic", tinv_hydrostatic_array);
            if ((rOpts.outputs & TensorInvariant::Deviatoric) != 0)
                add_cell("_deviatoric", tinv_deviatoric);
            continue;
        }

        const NDArray& src = rMesh.PointData(name);
        const std::size_t ncomp = data_num_components(src);
        if ((rOpts.outputs & TensorInvariant::Mises) != 0)
            add_point(out, base + "_mises" + rOpts.suffix, tinv_mises(src, ncomp));
        if ((rOpts.outputs & TensorInvariant::Principal) != 0)
            add_point(out, base + "_principal" + rOpts.suffix, tinv_principal(src, ncomp));
        if ((rOpts.outputs & TensorInvariant::Hydrostatic) != 0)
            add_point(out, base + "_hydrostatic" + rOpts.suffix, tinv_hydrostatic_array(src, ncomp));
        if ((rOpts.outputs & TensorInvariant::Deviatoric) != 0)
            add_point(out, base + "_deviatoric" + rOpts.suffix, tinv_deviatoric(src, ncomp));
    }
    return out;
}

}  // namespace meshioplusplus
