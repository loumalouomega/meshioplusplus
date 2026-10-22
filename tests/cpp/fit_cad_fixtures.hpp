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
 * @file fit_cad_fixtures.hpp
 * @brief The skins `fit_cad`'s tests share (`test_fit_cad.cpp`,
 *        `test_fit_cad_goldens.cpp`): generated in code from
 *        `analytic_surfaces.hpp`, segmented and fitted the way the tests feed
 *        `fit_cad`, so no fixture file (LFS or otherwise) is involved.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "analytic_surfaces.hpp"

namespace fcfix {

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

/// `segment` (angle method, the given curvature test and minimum patch) then
/// `fit_primitives` at `Tolerance`.
inline Mesh fitted(const Mesh& rSkin, double CurvatureTolerance = 0.0, double Tolerance = 1e-3,
                   std::int64_t MinPatchCells = 1) {
    meshioplusplus::SegmentOptions so;
    so.mCurvatureTolerance = CurvatureTolerance;
    so.mMinPatchCells = MinPatchCells;
    meshioplusplus::FitPrimitivesOptions fo;
    fo.mTolerance = Tolerance;
    return meshioplusplus::fit_primitives(meshioplusplus::segment(rSkin, so).mMesh, fo).mMesh;
}

/// A closed box whose top is the saddle z = 1 + 0.2 (x - 1/2)(y - 1/2) * 4: no
/// primitive explains the top.
inline Mesh saddle_box(int n) {
    Mesh box = analytic::box_skin(n);
    const NDArray& p = box.Points();
    NDArray q = NDArray::Uninit(DType::Float64, {box.NumPoints(), 3});
    for (std::size_t i = 0; i < box.NumPoints(); ++i) {
        const double x = p.As<double>()[i * 3], y = p.As<double>()[i * 3 + 1];
        double z = p.As<double>()[i * 3 + 2];
        z *= 1.0 + 0.8 * (x - 0.5) * (y - 0.5);
        q.As<double>()[i * 3] = x;
        q.As<double>()[i * 3 + 1] = y;
        q.As<double>()[i * 3 + 2] = z;
    }
    box.AssignPoints(std::move(q));
    return box;
}

/// A closed box whose top waves as z = 1 + A sin(2 pi x) sin(2 pi y): the
/// sides stay planes (the factor is 1 on them), the top is no primitive.
inline Mesh wavy_top_box(int n, double A = 0.05) {
    Mesh box = analytic::box_skin(n);
    const NDArray& p = box.Points();
    NDArray q = NDArray::Uninit(DType::Float64, {box.NumPoints(), 3});
    constexpr double kTwoPi = 6.28318530717958647692;
    for (std::size_t i = 0; i < box.NumPoints(); ++i) {
        const double x = p.As<double>()[i * 3], y = p.As<double>()[i * 3 + 1];
        q.As<double>()[i * 3] = x;
        q.As<double>()[i * 3 + 1] = y;
        q.As<double>()[i * 3 + 2] =
            p.As<double>()[i * 3 + 2] * (1.0 + A * std::sin(kTwoPi * x) * std::sin(kTwoPi * y));
    }
    box.AssignPoints(std::move(q));
    return box;
}

/// Every point of `rSkin` moved by up to `Amplitude` per coordinate, by an
/// integer LCG (the same bits on every platform).
inline Mesh jittered(Mesh rSkin, double Amplitude, std::uint64_t Seed = 7) {
    const NDArray& p = rSkin.Points();
    NDArray q = NDArray::Uninit(DType::Float64, {rSkin.NumPoints(), 3});
    std::uint64_t state = Seed;
    for (std::size_t i = 0; i < rSkin.NumPoints() * 3; ++i) {
        state = state * 6364136223846793005ull + 1442695040888963407ull;
        const double u = static_cast<double>(state >> 11) / 9007199254740992.0;
        q.As<double>()[i] = p.As<double>()[i] + Amplitude * (2.0 * u - 1.0);
    }
    rSkin.AssignPoints(std::move(q));
    return rSkin;
}

/// One named input of the golden set.
struct Fixture {
    std::string mName;
    Mesh mMesh;
};

/**
 * @brief The golden set: one skin per path through `fit_cad` (analytic faces,
 * exact and approximated edges, a raw skin segmented inside, freeform and
 * split freeform faces, thread flanks, conic and vertex blends, a noisy skin).
 * Names are stable: they key the golden tables.
 */
inline std::vector<Fixture> golden_fixtures() {
    std::vector<Fixture> f;
    f.push_back({"box", fitted(analytic::box_skin(3))});
    f.push_back({"capped_cylinder", fitted(analytic::capped_cylinder())});
    f.push_back({"turned_part", fitted(analytic::turned_part(), 1.0, 2e-3)});
    f.push_back({"filleted_plate", analytic::filleted_plate()});
    f.push_back({"saddle_box", fitted(saddle_box(24))});
    f.push_back({"wavy_capped", fitted(analytic::wavy_capped())});
    f.push_back({"threaded_rod", fitted(analytic::threaded_rod())});
    f.push_back({"rounded_box", meshioplusplus::fit_primitives(analytic::rounded_box()).mMesh});
    f.push_back(
        {"noisy_turned_part", fitted(jittered(analytic::turned_part(), 1e-3), 1.0, 2e-3, 20)});
    return f;
}

/// One input of the freeform-mode set (v17.12.0): a fitted skin and the
/// absolute tolerance its patches were fitted at and `fit_cad` runs at.
struct FreeformFixture {
    std::string mName;
    Mesh mMesh;
    double mTolerance = 0.0;
};

/// The bounding-box diagonal of a mesh's points.
inline double diagonal(const Mesh& rMesh) {
    const NDArray& p = rMesh.Points();
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i)
        for (int j = 0; j < 3; ++j) {
            lo[j] = std::min(lo[j], p.As<double>()[i * 3 + j]);
            hi[j] = std::max(hi[j], p.As<double>()[i * 3 + j]);
        }
    return std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                     (hi[2] - lo[2]) * (hi[2] - lo[2]));
}

/// `segment` by angles, then `fit_primitives` at the absolute `Tolerance`.
inline Mesh fitted_abs(const Mesh& rSkin, double Tolerance, std::int64_t MinPatchCells = 1) {
    meshioplusplus::SegmentOptions so;
    so.mMinPatchCells = MinPatchCells;
    meshioplusplus::FitPrimitivesOptions fo;
    fo.mTolerance = Tolerance;
    fo.mRelative = false;
    return meshioplusplus::fit_primitives(meshioplusplus::segment(rSkin, so).mMesh, fo).mMesh;
}

/**
 * @brief A segmentation that missed a crease: `rFitted` (a box from
 * `fitted_abs`) with its planar side at the smallest x merged into its type-0
 * (freeform) patch, the side's fit-table row dropped.
 */
inline Mesh merged_side(Mesh rFitted) {
    using meshioplusplus::detail::read_double;
    using meshioplusplus::detail::read_int;
    const double* p = rFitted.Points().As<double>();
    double lo = 1e300;
    for (std::size_t i = 0; i < rFitted.NumPoints(); ++i)
        lo = std::min(lo, p[3 * i]);
    const NDArray& ids = rFitted.FieldData("fit:patch_id");
    const NDArray& types = rFitted.FieldData("fit:patch_type");
    std::int64_t freeform = -1, side = -1;
    for (std::size_t r = 0; r < ids.Size(); ++r)
        if (read_int(types, r) == 0)
            freeform = read_int(ids, r);
    std::vector<NDArray> blocks;
    for (std::size_t b = 0; b < rFitted.NumCellBlocks(); ++b) {
        const NDArray& label = rFitted.CellData("segment:patch", b);
        const NDArray& conn = rFitted.Cells(b).Conn();
        const std::size_t k = rFitted.Cells(b).NodesPerCell(), nc = rFitted.Cells(b).NumCells();
        NDArray out = NDArray::Uninit(DType::Int64, {nc});
        for (std::size_t c = 0; c < nc; ++c) {
            bool on = true;
            for (std::size_t j = 0; j < k; ++j)
                on = on && p[3 * read_int(conn, c * k + j)] == lo;
            out.As<std::int64_t>()[c] = on ? freeform : read_int(label, c);
            if (on)
                side = read_int(label, c);
        }
        blocks.push_back(std::move(out));
    }
    rFitted.AddCellData("segment:patch", std::move(blocks));
    const NDArray& params = rFitted.FieldData("fit:patch_params");
    const NDArray& tols = rFitted.FieldData("fit:patch_tolerance");
    const std::size_t rows = ids.Size(), np = params.Size() / rows;
    NDArray nid = NDArray::Uninit(DType::Int64, {rows - 1});
    NDArray ntype = NDArray::Uninit(DType::Int64, {rows - 1});
    NDArray nparam = NDArray::Uninit(DType::Float64, {rows - 1, np});
    NDArray ntol = NDArray::Uninit(DType::Float64, {rows - 1});
    for (std::size_t r = 0, w = 0; r < rows; ++r) {
        if (read_int(ids, r) == side)
            continue;
        nid.As<std::int64_t>()[w] = read_int(ids, r);
        ntype.As<std::int64_t>()[w] = read_int(types, r);
        for (std::size_t k = 0; k < np; ++k)
            nparam.As<double>()[w * np + k] = read_double(params, r * np + k);
        ntol.As<double>()[w] = read_double(tols, r);
        ++w;
    }
    rFitted.AddFieldData("fit:patch_id", std::move(nid));
    rFitted.AddFieldData("fit:patch_type", std::move(ntype));
    rFitted.AddFieldData("fit:patch_params", std::move(nparam));
    rFitted.AddFieldData("fit:patch_tolerance", std::move(ntol));
    return rFitted;
}

/**
 * @brief The freeform-mode set: skins whose faces no primitive explains (a
 * monkey saddle, a Gaussian bump, an S band and a hood over planar boxes, a
 * genus-1 organic torus, a sphere with a handle, a genus-2 double torus,
 * Catmull-Clark limit surfaces with extraordinary vertices of valence 3, 5
 * and 6, and a loaf whose top and sides meet G1). With `Variants`, each also
 * with Gaussian noise of 0.1 and 0.5 times the tolerance (`_noise01`,
 * `_noise05`) and the boxes on a graded, anisotropic triangulation
 * (`_warped`). The tolerance is `2e-3` of each skin's diagonal; noisy skins
 * segment with a minimum patch of 20 cells, as a scan would. Names are
 * stable: they key the facet-area tables in doc/fit_cad.md.
 */
inline std::vector<FreeformFixture> freeform_fixtures(bool Variants = true) {
    struct Base {
        const char* mName;
        Mesh (*mMake)(bool Warp);
        bool mBox;
    };
    const Base bases[] = {
        {"monkey_saddle", [](bool W) { return analytic::monkey_saddle_box(24, W); }, true},
        {"gaussian_bump", [](bool W) { return analytic::gaussian_bump_box(24, W); }, true},
        {"s_band", [](bool W) { return analytic::s_band_box(24, W); }, true},
        {"hood", [](bool W) { return analytic::hood_box(24, W); }, true},
        {"organic_torus", [](bool) { return analytic::organic_torus(); }, false},
        {"sphere_handle", [](bool) { return analytic::sphere_with_handle(); }, false},
        {"double_torus", [](bool) { return analytic::double_torus(); }, false},
        {"cc_valence3", [](bool) { return analytic::cc_cage(3); }, false},
        {"cc_valence5", [](bool) { return analytic::cc_cage(5); }, false},
        {"cc_valence6", [](bool) { return analytic::cc_cage(6); }, false},
        {"g1_loaf", [](bool) { return analytic::g1_loaf(); }, false},
    };
    std::vector<FreeformFixture> out;
    for (const Base& b : bases) {
        const Mesh skin = b.mMake(false);
        const double tol = 2e-3 * diagonal(skin);
        out.push_back({b.mName, fitted_abs(skin, tol), tol});
        if (!Variants)
            continue;
        out.push_back({std::string(b.mName) + "_noise01",
                       fitted_abs(analytic::gaussian_noise(skin, 0.1 * tol), tol, 20), tol});
        out.push_back({std::string(b.mName) + "_noise05",
                       fitted_abs(analytic::gaussian_noise(skin, 0.5 * tol), tol, 20), tol});
        if (b.mBox)
            out.push_back({std::string(b.mName) + "_warped", fitted_abs(b.mMake(true), tol), tol});
    }
    return out;
}

}  // namespace fcfix
