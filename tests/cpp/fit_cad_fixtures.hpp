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

#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

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

}  // namespace fcfix
