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
 * @file detail/mesh_cad_size.hpp
 * @brief The size field of `mesh_cad_surface` (roadmap 6.3): one function of
 * space for the whole model, so faces and the edges they share agree on it by
 * construction.
 *
 * Its sources are point samples `(x_i, h_i)` taken on a first, uniform mesh:
 * the curvature size `2 sqrt(2 d / k - d^2)` for a chordal error `d` on a
 * surface of largest principal curvature `k`, and the gap size `t / n` for a
 * local feature size `t` (a ray along the normal to the nearest face that is
 * neither the point's own nor a neighbour) and `n` elements across it. The
 * field is their gradation-limited envelope
 *
 *     h(x) = min(H, min_i (h_i + s max(0, |x - x_i| - r))),
 *
 * whose gradient is at most `s = grading - 1` everywhere: the continuous
 * form of gradation (Persson 2006), exact for point sources. Each sample
 * holds its size within `r`, the reach of the first mesh's points (they are
 * about `H` apart), since the size it measured is the size of the patch of
 * surface around it, not of the point alone. A minimum of
 * doubles does not depend on the order it is taken in, so the field is the
 * same over thread counts.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// How much longer than the local size an edge may come out: refinement
/// splits a triangle whose circumradius exceeds this times an equilateral
/// one's of that size. A curvature or gap size is a bound, so the field's
/// target is the bound over this.
inline constexpr double kMeshCadSizeSlack = 1.5;

/// A point of the first mesh on one face: where it is, the face, the
/// surface's unit normal there (outward, with the face's orientation), and
/// the largest absolute principal curvature.
struct MeshCadSizeSample {
    Vec3 mX{0.0, 0.0, 0.0};
    std::int64_t mFace = -1;
    Vec3 mNormal{0.0, 0.0, 0.0};
    double mCurvature = 0.0;
};

class MeshCadBackgroundSize;

/// `h(x)` from its sources; see the file comment.
class MeshCadSizeField {
public:
    /// @param H the target size, the field's ceiling.
    /// @param Floor the field's floor.
    /// @param Slope the largest gradient, `grading - 1`.
    /// @param Reach how far a sample's own size holds.
    /// @param rX, rSize the samples whose size is below `H` (others are dropped).
    MeshCadSizeField(double H, double Floor, double Slope, double Reach,
                     const std::vector<Vec3>& rX, const std::vector<double>& rSize);
    double At(const Vec3& rX) const;
    /// Also take a background's size at every point, as given (its own
    /// gradation is the user's), within the field's floor and ceiling.
    void SetBackground(const MeshCadBackgroundSize* pBackground) { mpBackground = pBackground; }
    bool HasBackground() const { return mpBackground != nullptr; }
    /// Scale every size the field gives by `S` (roadmap 6.8.3: the
    /// all-quad route meshes at twice the size and splits each quad in
    /// four, so it asks the field for `2 h(x)`).
    void SetScale(double S) { mScale = S; }
    double Ceiling() const { return mScale * mH; }
    double Floor() const { return mScale * mFloor; }
    /// The smallest sample size (the ceiling with none).
    double Smallest() const { return mScale * mSmallest; }
    std::size_t NumSamples() const { return mX.size(); }

private:
    double mH = 0.0;
    double mFloor = 0.0;
    double mSlope = 0.0;
    double mReach = 0.0;
    double mSmallest = 0.0;
    double mScale = 1.0;
    double mCell = 1.0;
    Vec3 mLo{0.0, 0.0, 0.0};
    std::vector<Vec3> mX;
    std::vector<double> mSize;
    std::unordered_map<std::int64_t, std::vector<std::int32_t>> mCells;
    std::array<std::int64_t, 3> mDims{1, 1, 1};
    const MeshCadBackgroundSize* mpBackground = nullptr;

    std::int64_t Key(std::int64_t I, std::int64_t J, std::int64_t K) const;
};

/// The size each sample asks for: the curvature size for chordal error
/// @p Deviation (0: off) and the gap size for @p CellsAcrossGaps elements
/// across (0: off), measured by rays on @p rFirst's triangles (`cad:face`
/// names each one's face); @p rAdjacent lists each face's neighbours. Never
/// above @p H nor below @p Floor.
std::vector<double> mesh_cad_sample_sizes(const std::vector<MeshCadSizeSample>& rSamples,
                                          const Mesh& rFirst,
                                          const std::vector<std::vector<std::int64_t>>& rAdjacent,
                                          double H, double Floor, double Deviation,
                                          double CellsAcrossGaps);

/// A background mesh's size (roadmap 6.3.1): its one-component point_data
/// interpolated linearly at the nearest point of its surface (a volume mesh
/// is skinned first). Built once; `At` is safe to call from many threads.
class MeshCadBackgroundSize {
public:
    MeshCadBackgroundSize(const Mesh& rBackground, const std::string& rArray);
    ~MeshCadBackgroundSize();
    /// Not positive or not finite there, or no surface: infinity.
    double At(const Vec3& rX) const;

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
};

/**
 * @brief An anisotropic metric from a field's Hessian (roadmap 6.10.2):
 * per point of a background mesh (a volume mesh is skinned first), `M =
 * R |Λ| R^T / (8 ε)` from the symmetrised Hessian `R Λ R^T` in its 9-component
 * point_data and the target interpolation error `ε` -- the metric under
 * which linear interpolation of a quadratic errs by at most `ε` on a unit
 * edge (A. Loseille and F. Alauzet, "Continuous mesh framework part I",
 * SIAM J. Numer. Anal. 49(1), 2011; the bound `h <= sqrt(8 ε / |λ|)` along
 * an eigenvector). Interpolated linearly at the nearest point of the
 * background's surface (a convex combination keeps it positive
 * semi-definite), projected on that surface's tangent plane (`P M P`), then
 * clamped so that every direction's size lies between a floor and the local
 * isotropic size: the metric only refines. Built once; `At` is safe to call
 * from many threads.
 */
class MeshCadMetricField {
public:
    MeshCadMetricField(const Mesh& rBackground, const std::string& rArray, double Error,
                       double Floor);
    ~MeshCadMetricField();
    /// The metric at `rX` (row-major 3x3) with its sizes clamped to
    /// [`Floor`, `Size`], and its square root in `pRoot` when given; off the
    /// background (no surface) the isotropic `I / Size^2`.
    void At(const Vec3& rX, double Size, std::array<double, 9>& rM,
            std::array<double, 9>* pRoot = nullptr) const;
    /// Every size the metric asks for, and its floor, times `Scale`: the
    /// metric over `Scale^2` (all-quad meshes at twice the size, 6.8.3).
    void SetScale(double Scale);

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
};

}  // namespace detail
}  // namespace meshioplusplus
