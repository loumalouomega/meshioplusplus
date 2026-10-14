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
 * @file operations/iga.hpp
 * @brief The isogeometric layer: evaluate, refine and sample the NURBS
 * models a reader keeps beside its sampled mesh -- `NurbsInfo` (MFEM NURBS
 * meshes) and `BrepInfo` (Kratos CAD JSON, STEP, IGES, `fit_cad`) -- without
 * going through the sampling.
 *
 * - **Evaluation** (`iga_evaluate`): the points, first derivatives and unit
 *   normals of a patch (or a trimmed face) at parametric coordinates, rational
 *   in 1, 2 or 3 parametric dimensions; a `NurbsInfo` patch also evaluates
 *   its NURBS fields, a `BrepInfo` face reports whether each point lies
 *   inside its trimming loops.
 * - **Refinement** (`iga_refine`): exact knot insertion (h) and degree
 *   elevation (p) of every patch, k-refinement being the two in that order
 *   (elevate, then insert, which keeps the highest continuity). The geometry
 *   and the parametrisation do not change, so a trimmed face's trimming
 *   curves stay exact; a `NurbsInfo`'s fields refine with the geometry. The
 *   patch-level kernels (`iga_insert_knots`, `iga_elevate_degree`) take one
 *   direction at a time.
 * - **Sampling** (`iga_sample`): the model as Lagrange cells again.
 *
 * Every direction beyond a patch's parametric dimension is ignored. A control
 * net is rational when any of its weights is not exactly 1; a face whose
 * `mRational` is false is a polynomial B-spline whatever its weights are, and
 * refines to a patch of unit weights. The Python twin (`meshioplusplus._iga`)
 * repeats the arithmetic, so both engines give the same bits.
 * See doc/iga.md.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/nurbs.hpp"

namespace meshioplusplus {

/// Cell data (Int64) a sampled `NurbsInfo` carries: each cell's patch.
inline constexpr const char* kIgaPatchName = "iga:patch";
/// Cell data (Int64, 3 columns): each cell's knot span per direction, counted
/// over the spans of nonzero length (0 for an unused direction).
inline constexpr const char* kIgaSpanName = "iga:span";

/// How `iga_refine` refines every patch of a model.
struct IgaRefineOptions {
    /// Equal parts every nonzero knot span is split into, per direction (1: none).
    std::array<std::int32_t, 3> mSubdivisions{1, 1, 1};
    /// Degree elevation per direction (0: none), applied before the knots are inserted.
    std::array<std::int32_t, 3> mElevate{0, 0, 0};
};

/// How `iga_sample` turns a model into cells.
struct IgaSampleOptions {
    /// Cells per knot span and direction; 0 is automatic (one per span; on a
    /// trimmed face about 16 per direction so its trim shows; for a `NurbsInfo`
    /// with a patch topology, the MFEM reader's own sampling).
    std::int32_t mSubdivisions = 0;
    /// `BrepInfo` only: points per knot span a trimming curve is sampled at, 2 to 4096.
    std::int32_t mTrimSamples = 16;
    /// `BrepInfo` only: keep only the cells inside the trimming loops.
    bool mTrim = true;
    /// `BrepInfo` only: sample the edges that have a 3-D curve as line cells.
    bool mEdges = true;
};

/// What `iga_evaluate` returns; the arrays a request does not ask for are empty.
struct MESHIOPLUSPLUS_API IgaEvaluation {
    /// Float64 `(N, space dimension)`.
    NDArray mPoints{DType::Float64, {0, 0}};
    /// Float64 `(N, parametric dimension, space dimension)`: `d x / d u_k`.
    NDArray mDerivatives{DType::Float64, {0}};
    /// Float64 `(N, 3)`: the unit normal `d1 x d2 / |d1 x d2|` of a surface in
    /// 3-D (zeros where the tangents are parallel).
    NDArray mNormals{DType::Float64, {0}};
    /// One Float64 `(N, components)` per `NurbsInfo::mFieldNames` entry.
    std::vector<NDArray> mFields;
    /// UInt8 `(N)`: 1 where a face's trimming loops contain the point (a face
    /// that is not trimmed contains every point).
    NDArray mInside{DType::UInt8, {0}};
};

/**
 * @brief Insert knots into direction @p Direction of one patch, exactly.
 * @param rKnots The knots, ascending, strictly inside the parameter range
 * (repeating a knot raises its multiplicity); fields refine with the net.
 * @param Rational Evaluate the weights (false: a polynomial B-spline whose
 * weights are written back as 1).
 * @throws std::invalid_argument for a bad direction or a knot out of range.
 */
MESHIOPLUSPLUS_API NurbsPatch iga_insert_knots(const NurbsPatch& rPatch, std::int32_t Direction,
                                               const std::vector<double>& rKnots,
                                               bool Rational = true);

/**
 * @brief Raise the degree of direction @p Direction of one patch by @p Times,
 * exactly, keeping every knot's continuity.
 * @throws std::invalid_argument for a bad direction, `Times < 0` or a degree above 32.
 */
MESHIOPLUSPLUS_API NurbsPatch iga_elevate_degree(const NurbsPatch& rPatch, std::int32_t Direction,
                                                 std::int32_t Times, bool Rational = true);

/**
 * @brief Refine every patch of @p rInfo: degree elevation, then knot insertion.
 *
 * The patch topology and field names are kept; the cell-to-patch map and the
 * boundary cell map are cleared (the sampling they described is gone; sample
 * again with `iga_sample`).
 */
MESHIOPLUSPLUS_API NurbsInfo iga_refine(const NurbsInfo& rInfo, const IgaRefineOptions& rOptions);

/**
 * @brief Refine every face surface of @p rInfo.
 *
 * Trimming curves, edges and vertices are untouched (the parametrisation is);
 * a face's control-point ids are cleared.
 * @throws std::invalid_argument when a patch is shared by faces that disagree
 * on `mRational`.
 */
MESHIOPLUSPLUS_API BrepInfo iga_refine(const BrepInfo& rInfo, const IgaRefineOptions& rOptions);

/**
 * @brief Evaluate patch @p Patch of @p rInfo, and its fields, at the parameters
 * @p rParams (Float64 `(N, parametric dimension)`).
 * @param Derivatives 0 or 1: with 1, also the first derivatives and, for a
 * surface in 3-D, the unit normals.
 * @throws std::invalid_argument for a bad patch, parameters or a parameter
 * outside the patch's knot range (by more than rounding).
 */
MESHIOPLUSPLUS_API IgaEvaluation iga_evaluate(const NurbsInfo& rInfo, std::size_t Patch,
                                              const NDArray& rParams, std::int32_t Derivatives = 0);

/**
 * @brief Evaluate face @p Face of @p rInfo, with `mInside` from its trimming
 * loops (each trim sampled at @p TrimSamples points per knot span).
 */
MESHIOPLUSPLUS_API IgaEvaluation iga_evaluate(const BrepInfo& rInfo, std::size_t Face,
                                              const NDArray& rParams, std::int32_t Derivatives = 0,
                                              std::int32_t TrimSamples = 16);

/**
 * @brief Sample every patch of @p rInfo into Lagrange cells of the highest
 * degree (`line`/`line3`/`VTK_LAGRANGE_CURVE`, the quad family or the
 * hexahedron family by parametric dimension), their nodes the exact
 * (rational) geometry at the cell's lattice, the NURBS fields evaluated as
 * point data, `iga:patch` and `iga:span` as cell data.
 *
 * Points are shared within a patch, not between patches. A model with its
 * MFEM patch topology (`mTopology`) and `mSubdivisions == 0` samples as
 * `sample_mfem_nurbs` does instead -- exactly the mesh `read_mfem` gives, with
 * its points shared across patches and its boundary cells -- so a refined
 * model writes back through `write_mfem`.
 * @param pMapped When not null, receives @p rInfo with `mCellPatch` and
 * `mCellSpan` filled from the cells and `mBoundaryCellPatch` cleared.
 */
MESHIOPLUSPLUS_API Mesh iga_sample(const NurbsInfo& rInfo, const IgaSampleOptions& rOptions = {},
                                   NurbsInfo* pMapped = nullptr);

/// Sample a B-rep as `sample_brep` does (`mSubdivisions`, `mTrimSamples`, `mTrim`, `mEdges`).
MESHIOPLUSPLUS_API Mesh iga_sample(const BrepInfo& rInfo, const IgaSampleOptions& rOptions = {});

}  // namespace meshioplusplus
