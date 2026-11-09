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
 * @file detail/hausdorff_samples.hpp
 * @brief Every sample distance `hausdorff_distance` reduces, kept rather than
 * reduced: the distance of each sample of A to B, then of each sample of B to
 * A, sampled and measured exactly as `hausdorff_distance` does (so the largest
 * is its Hausdorff distance). `decimate`'s `max_deviation` takes a quantile of
 * them (roadmap 7.8.5).
 *
 * A **core-private** header: no installed header names it.
 */

// System includes
#include <cstdint>
#include <memory>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus::detail {

/// The unsigned sample distances A to B then B to A (`FaceSamples` as
/// `HausdorffOptions::mFaceSamples`).
std::vector<double> hausdorff_sample_distances(const Mesh& rA, const Mesh& rB,
                                               std::int64_t FaceSamples);

/**
 * @brief The same against one fixed reference, its samples and its distance
 * query built once: `Distances(rOther)` is
 * `hausdorff_sample_distances(rOther, rReference, FaceSamples)` when
 * `ReferenceFaceSamples == FaceSamples`; a dense reference (a mesh being
 * decimated) is sampled well enough at its vertices (`ReferenceFaceSamples`
 * 0) while the coarse other side keeps its face samples.
 */
class HausdorffSampler {
public:
    HausdorffSampler(const Mesh& rReference, std::int64_t ReferenceFaceSamples,
                     std::int64_t FaceSamples);
    ~HausdorffSampler();
    HausdorffSampler(const HausdorffSampler&) = delete;
    HausdorffSampler& operator=(const HausdorffSampler&) = delete;

    std::vector<double> Distances(const Mesh& rOther) const;

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
};

/// The @p Coverage quantile of @p Values (the smallest value at least that
/// share of them do not exceed); 0 for none.
double coverage_quantile(std::vector<double> Values, double Coverage);

}  // namespace meshioplusplus::detail
