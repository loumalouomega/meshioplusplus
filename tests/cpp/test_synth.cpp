// SPDX-License-Identifier: MIT
// The soundtrack synthesizer (detail/synth.*): determinism, length, headroom and the loop seam.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/synth.hpp"

using namespace meshioplusplus::detail;

namespace {

std::uint32_t le32(const std::string& rS, std::size_t At) {
    std::uint32_t v = 0;
    for (int k = 3; k >= 0; --k)
        v = (v << 8) | static_cast<unsigned char>(rS[At + static_cast<std::size_t>(k)]);
    return v;
}

}  // namespace

TEST(Synth, TheSameOptionsGiveTheSameBytes) {
    SynthOptions o;
    const auto a = synth_synthwave(o);
    const auto b = synth_synthwave(o);
    EXPECT_EQ(a, b);
    EXPECT_EQ(synth_wav(a, o.mSampleRate), synth_wav(b, o.mSampleRate));
    SynthOptions other = o;
    other.mSeed = 2;
    EXPECT_NE(synth_synthwave(other), a);
    other = o;
    other.mKey = 0;
    EXPECT_NE(synth_synthwave(other), a);
}

TEST(Synth, ThePinnedLoopDoesNotDrift) {
    // A checksum of the default loop, so a change to the synthesis is a visible,
    // deliberate edit of this number (no floating point touches a sample, so it is
    // the same on every platform).
    SynthOptions o;
    const auto samples = synth_synthwave(o);
    std::uint64_t h = 1469598103934665603ull;
    for (std::int16_t s : samples) {
        h ^= static_cast<std::uint16_t>(s);
        h *= 1099511628211ull;
    }
    EXPECT_EQ(h, 5088833146301154402ull) << "update the pinned checksum to " << h;
}

TEST(Synth, TheLengthIsBarsTimesBeatsTimesSamplesPerBeat) {
    SynthOptions o;
    o.mBars = 4;
    o.mTempo = 90.0;
    o.mSampleRate = 22050;
    const std::size_t spb = synth_samples_per_beat(o);
    EXPECT_EQ(spb, 14700u);  // 22050 * 60 / 90
    EXPECT_EQ(synth_length(o), 4u * 4u * spb);
    EXPECT_EQ(synth_synthwave(o).size(), synth_length(o));
}

TEST(Synth, ThePeakStaysBelowFullScaleAndFollowsTheGain) {
    SynthOptions o;
    for (double gain : {0.1, 0.3, 0.9}) {
        o.mGain = gain;
        const auto s = synth_synthwave(o);
        int peak = 0;
        for (std::int16_t v : s)
            peak = std::max(peak, std::abs(static_cast<int>(v)));
        EXPECT_LT(peak, 32767);
        EXPECT_NEAR(static_cast<double>(peak) / 32767.0, gain, 0.01);
    }
}

TEST(Synth, TheLoopJoinsWithoutAStep) {
    // The last sample and the first are no further apart than the largest step
    // inside the loop, so repeating it makes no click of its own.
    for (std::uint32_t seed : {1u, 2u, 3u, 7u}) {
        SynthOptions o;
        o.mSeed = seed;
        const auto s = synth_synthwave(o);
        int largest = 0;
        for (std::size_t i = 1; i < s.size(); ++i)
            largest = std::max(largest, std::abs(s[i] - s[i - 1]));
        EXPECT_LE(std::abs(s.front() - s.back()), largest) << "seed " << seed;
    }
}

TEST(Synth, TheSoundHasEnergyAndIsNotSilent) {
    SynthOptions o;
    const auto s = synth_synthwave(o);
    double sum = 0.0;
    for (std::int16_t v : s)
        sum += static_cast<double>(v) * v;
    EXPECT_GT(std::sqrt(sum / static_cast<double>(s.size())), 500.0);
}

TEST(Synth, ARiffWaveContainerWrapsTheSamples) {
    SynthOptions o;
    o.mBars = 1;
    const auto s = synth_synthwave(o);
    const std::string wav = synth_wav(s, o.mSampleRate);
    ASSERT_EQ(wav.size(), 44u + 2u * s.size());
    EXPECT_EQ(wav.substr(0, 4), "RIFF");
    EXPECT_EQ(wav.substr(8, 8), "WAVEfmt ");
    EXPECT_EQ(le32(wav, 4), 36u + 2u * s.size());
    EXPECT_EQ(le32(wav, 24), static_cast<std::uint32_t>(o.mSampleRate));
    EXPECT_EQ(wav.substr(36, 4), "data");
    EXPECT_EQ(le32(wav, 40), 2u * s.size());
    EXPECT_EQ(static_cast<unsigned char>(wav[44]), static_cast<std::uint16_t>(s[0]) & 0xFFu);
}

TEST(Synth, OutOfRangeArgumentsAreRefusedByName) {
    SynthOptions o;
    o.mTempo = 30.0;
    EXPECT_THROW(synth_synthwave(o), std::invalid_argument);
    o = {};
    o.mKey = 12;
    EXPECT_THROW(synth_synthwave(o), std::invalid_argument);
    o = {};
    o.mBars = 0;
    EXPECT_THROW(synth_synthwave(o), std::invalid_argument);
    o = {};
    o.mSampleRate = 100;
    EXPECT_THROW(synth_synthwave(o), std::invalid_argument);
    o = {};
    o.mGain = 1.0;
    EXPECT_THROW(synth_synthwave(o), std::invalid_argument);
    try {
        o = {};
        o.mTempo = 200.0;
        synth_synthwave(o);
        FAIL();
    } catch (const std::invalid_argument& rErr) {
        EXPECT_NE(std::string(rErr.what()).find("tempo"), std::string::npos);
    }
}
