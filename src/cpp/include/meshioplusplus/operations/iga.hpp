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
    /// Float64 `(N, parametric dimension, parametric dimension, space
    /// dimension)`: the full symmetric second derivatives `d2 x / d u_k d u_l`,
    /// filled only by `Derivatives == 2` (roadmap 9.8.1).
    NDArray mSecondDerivatives{DType::Float64, {0}};
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
 * @param Derivatives 0, 1 or 2: with 1, also the first derivatives and, for a
 * surface in 3-D, the unit normals; with 2, also the second derivatives
 * (`mSecondDerivatives`). The first derivatives of order 2 come from the
 * second-order basis recurrence and may differ in the last bits from those of
 * order 1.
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

/// What `iga_remove_knots` did, per distinct knot value, ascending.
struct IgaRemoveKnotsReport {
    /// The knot values tried.
    std::vector<double> mKnots;
    /// Copies of each that were asked for.
    std::vector<std::int32_t> mRequested;
    /// Copies of each that were removed.
    std::vector<std::int32_t> mRemoved;
    /// The sum of `mRemoved`.
    std::int32_t mNumRemoved = 0;
};

/**
 * @brief Knot removal (roadmap 9.8.3): the NURBS Book's A5.8 on one patch.
 *
 * Removes interior knots from direction @p Direction of patch @p Patch while
 * the geometry (and the NURBS fields) deviates by at most @p Tolerance: the
 * book's bound `TOL = d * wmin / (1 + Pmax)` on the homogeneous net, `wmin` the
 * smallest weight and `Pmax` the largest control-point norm, so the curve or
 * surface moves by at most `Tolerance`. A tolerance of 0 removes only what is
 * exact up to rounding (`1e-12` relative to the net's size), which is what
 * undoes a knot insertion. A knot of a surface or volume goes only when it
 * goes from every line of the net.
 *
 * @param rKnots The knots to remove, each listing one copy (a value listed
 *        twice asks for two); empty asks for every copy of every interior knot.
 * @param pReport Receives, per distinct value, the copies asked for and removed.
 * @return the model with that patch replaced; the cell maps cleared as
 *         `iga_refine` does.
 * @throws std::invalid_argument for a bad patch or direction, a negative or
 *         non-finite tolerance, a knot that is not an interior knot of the
 *         direction (or listed more often than its multiplicity), or a model
 *         whose patch topology joins several patches (they share knot vectors).
 */
MESHIOPLUSPLUS_API NurbsInfo iga_remove_knots(const NurbsInfo& rInfo, std::size_t Patch,
                                              std::int32_t Direction,
                                              const std::vector<double>& rKnots, double Tolerance,
                                              IgaRemoveKnotsReport* pReport = nullptr);

/**
 * @brief Bezier extraction of one patch (roadmap 9.8.2): the patch as VTK Bezier cells.
 *
 * Every interior knot is inserted up to multiplicity `p` in each direction
 * (exactly, on the homogeneous net, so a rational patch stays exact), and the
 * refined net is cut at its nonzero knot spans: one `VTK_BEZIER_CURVE`,
 * `VTK_BEZIER_QUADRILATERAL` or `VTK_BEZIER_HEXAHEDRON` cell per span product,
 * with `(p+1)^dim` points in VTK's Lagrange node order, in Cartesian
 * coordinates (the homogeneous values divided by the weight). Points are not
 * shared between cells. Cell data `HigherOrderDegrees` (Float64, `(cells, 3)`,
 * unused directions 0), `iga:patch` and `iga:span` (as `iga_sample` writes
 * them); point data `RationalWeights` when the patch is rational (any weight
 * other than 1) and one array per NURBS field.
 * @throws std::invalid_argument for a bad patch index or patch, or unequal
 *         degrees in the patch's directions (elevate with `iga_elevate_degree`).
 */
MESHIOPLUSPLUS_API Mesh iga_bezier_extract(const NurbsInfo& rInfo, std::size_t Patch);

/// Sample a B-rep as `sample_brep` does (`mSubdivisions`, `mTrimSamples`, `mTrim`, `mEdges`).
MESHIOPLUSPLUS_API Mesh iga_sample(const BrepInfo& rInfo, const IgaSampleOptions& rOptions = {});

}  // namespace meshioplusplus
