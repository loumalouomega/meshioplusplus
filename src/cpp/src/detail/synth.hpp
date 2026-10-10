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
 * @file detail/synth.hpp
 * @brief A small, deterministic synthesizer for the viewer's optional
 * soundtrack: a slow minor-key loop as 16-bit mono PCM, and the WAV container
 * around it.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * Everything is integer arithmetic on phase accumulators: waveforms come from
 * the accumulators' top bits, the sine is a fixed-point polynomial, the filter
 * is a one-pole integer low-pass and the noise is an xorshift generator. No
 * floating-point value and no libm call touches a sample, so the same
 * `SynthOptions` give the same bytes on every compiler, platform and
 * optimisation level. The loop is original: nothing is sampled, copied from or
 * modelled on an existing piece, and no audio file ships with the library.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace meshioplusplus {
namespace detail {

struct SynthOptions {
    /// Picks the chord progression and its variation.
    std::uint32_t mSeed = 1;
    /// Beats per minute, 60 to 140 (the 1980s sit around 90 to 110).
    double mTempo = 100.0;
    /// The tonic as a semitone above C (0 to 11); the default is A.
    int mKey = 9;
    /// Length of the loop in bars of four beats.
    int mBars = 8;
    /// Samples per second.
    int mSampleRate = 22050;
    /// Loudness of the loudest sample as a fraction of full scale, in (0, 0.9].
    double mGain = 0.3;
};

/// Samples in one beat for these options.
std::size_t synth_samples_per_beat(const SynthOptions& rOptions);

/// Samples in the whole loop: `mBars * 4 * synth_samples_per_beat`.
std::size_t synth_length(const SynthOptions& rOptions);

/// The loop as 16-bit mono PCM. Throws std::invalid_argument on a tempo, key,
/// length, rate or gain out of range.
std::vector<std::int16_t> synth_synthwave(const SynthOptions& rOptions);

/// A RIFF/WAVE file (PCM, mono, 16-bit, little-endian) around the samples.
std::string synth_wav(const std::vector<std::int16_t>& rSamples, int SampleRate);

}  // namespace detail
}  // namespace meshioplusplus
