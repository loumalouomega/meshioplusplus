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

// The soundtrack synthesizer. See detail/synth.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes (private, not installed)
#include "synth.hpp"

namespace meshioplusplus {
namespace detail {
namespace {

// --- Fixed-point oscillators: a 32-bit phase, waveforms from its top bits ---

int synth_tri(std::uint32_t Phase) {
    const int t = static_cast<int>(Phase >> 16);  // 0 .. 65535
    return t < 32768 ? t * 2 - 32767 : 98303 - t * 2;
}

int synth_saw(std::uint32_t Phase) {
    return static_cast<int>(Phase >> 16) - 32768;
}

int synth_pulse(std::uint32_t Phase, std::uint32_t Duty) {
    return Phase < Duty ? 24000 : -24000;
}

// sin(pi x / 2) ~ x (3 - x^2) / 2 on the triangle wave: a smooth, odd, integer sine.
int synth_sine(std::uint32_t Phase) {
    const std::int64_t x = synth_tri(Phase);
    const std::int64_t x2 = (x * x) >> 15;
    return static_cast<int>((x * (98304 - x2)) >> 16);
}

// The twelve pitch classes of octave 0 in Q8 hertz (C0 = 16.35 Hz, A0 = 27.5 Hz).
constexpr std::array<std::uint32_t, 12> kPitchQ8 = {4186, 4435, 4699, 4978, 5274, 5588,
                                                    5920, 6272, 6645, 7040, 7459, 7902};

// The phase step per sample of a MIDI-like note number (C0 = 0).
std::uint32_t synth_step(int Note, int Rate) {
    const std::uint64_t q8 =
        static_cast<std::uint64_t>(kPitchQ8[static_cast<std::size_t>(Note % 12)]) << (Note / 12);
    return static_cast<std::uint32_t>((q8 << 24) / static_cast<std::uint64_t>(Rate));
}

std::uint32_t synth_noise(std::uint32_t& rState) {
    rState ^= rState << 13;
    rState ^= rState >> 17;
    rState ^= rState << 5;
    return rState;
}

// Linear decay from 32767 to 0 over `Length` samples.
int synth_decay(std::int64_t Age, std::int64_t Length) {
    return Age < 0 || Age >= Length ? 0 : static_cast<int>(32767 - 32767 * Age / Length);
}

struct SynthChord {
    int mRoot;  // semitones above the tonic
    bool mMajor;
};

constexpr std::array<std::array<SynthChord, 4>, 4> kProgressions = {{
    {{{0, false}, {8, true}, {3, true}, {10, true}}},   // i  VI  III VII
    {{{0, false}, {5, false}, {8, true}, {7, false}}},  // i  iv  VI  v
    {{{0, false}, {3, true}, {10, true}, {8, true}}},   // i  III VII VI
    {{{0, false}, {10, true}, {8, true}, {10, true}}},  // i  VII VI  VII
}};

// Sixteenth-note hit masks (bit k = step k of the bar) for the bass.
constexpr std::array<std::uint32_t, 4> kBassHits = {0x9292, 0x8888 | 0x0020, 0xA4A4, 0x9249};
// Arpeggio orders over the chord's three tones.
constexpr std::array<std::array<int, 4>, 3> kArp = {{{0, 1, 2, 1}, {0, 2, 1, 2}, {2, 1, 0, 1}}};

}  // namespace

std::size_t synth_samples_per_beat(const SynthOptions& rOptions) {
    return static_cast<std::size_t>(
        static_cast<double>(rOptions.mSampleRate) * 60.0 / rOptions.mTempo + 0.5);
}

std::size_t synth_length(const SynthOptions& rOptions) {
    return static_cast<std::size_t>(rOptions.mBars) * 4 * synth_samples_per_beat(rOptions);
}

std::vector<std::int16_t> synth_synthwave(const SynthOptions& rOptions) {
    const SynthOptions& o = rOptions;
    if (!(o.mTempo >= 60.0 && o.mTempo <= 140.0))
        throw std::invalid_argument("meshio++: synth: the tempo must lie in [60, 140] BPM");
    if (o.mKey < 0 || o.mKey > 11)
        throw std::invalid_argument("meshio++: synth: the key must lie in [0, 11]");
    if (o.mBars < 1 || o.mBars > 64)
        throw std::invalid_argument("meshio++: synth: the length must lie in [1, 64] bars");
    if (o.mSampleRate < 8000 || o.mSampleRate > 48000)
        throw std::invalid_argument("meshio++: synth: the sample rate must lie in [8000, 48000]");
    if (!(o.mGain > 0.0 && o.mGain <= 0.9))
        throw std::invalid_argument("meshio++: synth: the gain must lie in (0, 0.9]");

    const std::size_t spb = synth_samples_per_beat(o);
    const std::size_t bar = 4 * spb;
    const std::size_t total = o.mBars * bar;
    const std::size_t tail = 6 * spb;
    const int rate = o.mSampleRate;
    const auto& progression = kProgressions[o.mSeed % 4];
    const auto& arp = kArp[(o.mSeed / 4) % 3];
    const std::uint32_t bass_hits = kBassHits[(o.mSeed / 16) % 4];
    const std::size_t delay_len = 3 * spb / 4;

    std::vector<std::int32_t> mix(total + tail, 0);
    std::vector<std::int32_t> ring(delay_len, 0);
    std::size_t ring_at = 0;
    std::uint32_t noise = 0x9E3779B9u ^ (o.mSeed * 2654435761u);
    std::array<std::uint32_t, 6> pad_phase = {0,          0x2AAAAAAA, 0x55555555,
                                              0x80000000, 0xAAAAAAAA, 0xD5555555};
    std::uint32_t arp_phase = 0;
    std::uint32_t bass_phase = 0;
    std::uint32_t kick_phase = 0;
    int lp_state = 0;
    int hat_prev = 0;

    for (std::size_t i = 0; i < total + tail; ++i) {
        std::int32_t dry_bed = 0;  // pad and arpeggio: what the delay hears
        std::int32_t drums = 0;
        if (i < total) {
            const std::size_t bar_index = i / bar;
            const std::size_t in_bar = i % bar;
            const std::size_t sixteenth = (in_bar * 16) / bar;
            const std::size_t sixteenth_start = (sixteenth * bar) / 16;
            const std::size_t step_len = ((sixteenth + 1) * bar) / 16 - sixteenth_start;
            const std::int64_t in_step = static_cast<std::int64_t>(in_bar - sixteenth_start);
            const SynthChord chord = progression[bar_index % 4];
            const int root = o.mKey + chord.mRoot;
            const int tones[3] = {root, root + (chord.mMajor ? 4 : 3), root + 7};

            // Pad: two detuned saws per chord tone, swelling in and out of each bar.
            std::int64_t pad = 0;
            for (int v = 0; v < 3; ++v) {
                const std::uint32_t step = synth_step(tones[v] + 36, rate);
                pad_phase[static_cast<std::size_t>(2 * v)] += step;
                pad_phase[static_cast<std::size_t>(2 * v + 1)] += step + (step >> 8);
                pad += synth_saw(pad_phase[static_cast<std::size_t>(2 * v)]) +
                       synth_saw(pad_phase[static_cast<std::size_t>(2 * v + 1)]);
            }
            const std::int64_t attack =
                std::min<std::int64_t>(32767, static_cast<std::int64_t>(in_bar) * 32767 /
                                                  static_cast<std::int64_t>(bar / 4));
            const std::int64_t release =
                std::min<std::int64_t>(32767, static_cast<std::int64_t>(bar - in_bar) * 32767 /
                                                  static_cast<std::int64_t>(bar / 8));
            const std::int64_t swell = std::min(attack, release);
            pad = (pad / 6) * swell >> 15;
            // A slow low-pass sweep across the whole loop.
            const int lfo = synth_tri(
                static_cast<std::uint32_t>((static_cast<std::uint64_t>(i) << 32) / total));
            const int coeff = 2500 + ((lfo < 0 ? -lfo : lfo) * 9000 >> 15);
            lp_state += static_cast<int>((static_cast<std::int64_t>(pad - lp_state) * coeff) >> 15);
            dry_bed += lp_state * 5 / 10;

            // Arpeggio: a quarter-width pulse stepping through the chord on sixteenths.
            const int arp_note = tones[arp[sixteenth % 4]] + 60;
            arp_phase += synth_step(arp_note, rate);
            const int arp_env = synth_decay(in_step, static_cast<std::int64_t>(step_len) * 7 / 10);
            dry_bed += (synth_pulse(arp_phase, 0x40000000u) * arp_env >> 15) * 3 / 10;

            // Bass: a sine on the root, two octaves below the pad, on the hit mask.
            bass_phase += synth_step(root + 24, rate);
            const bool hit = ((bass_hits >> sixteenth) & 1u) != 0;
            const int bass_env =
                hit ? synth_decay(in_step, static_cast<std::int64_t>(step_len) * 2) : 0;
            drums += (synth_sine(bass_phase) * bass_env >> 15) * 6 / 10;

            // Kick: a sine that falls from 140 Hz to 45 Hz, on every beat.
            const std::int64_t in_beat = static_cast<std::int64_t>(i % spb);
            const std::int64_t kick_len = static_cast<std::int64_t>(rate) / 4;
            const std::int64_t sweep = static_cast<std::int64_t>(rate) / 8;
            const std::int64_t f_q8 =
                (in_beat < sweep ? 140 * 256 - (140 - 45) * 256 * in_beat / sweep : 45 * 256);
            kick_phase += static_cast<std::uint32_t>(((static_cast<std::uint64_t>(f_q8)) << 24) /
                                                     static_cast<std::uint64_t>(rate));
            drums += (synth_sine(kick_phase) * synth_decay(in_beat, kick_len) >> 15) * 9 / 10;

            // Hat: high-passed noise on the off-beat eighths, quieter on the other sixteenths.
            const int n = static_cast<int>(synth_noise(noise) >> 17) - 16384;
            const int hp = n - hat_prev;
            hat_prev = n;
            const bool off_beat = (sixteenth % 4) == 2;
            const int hat_env =
                synth_decay(in_step, static_cast<std::int64_t>(rate) / 25) / (off_beat ? 1 : 4);
            drums += (hp * hat_env >> 15) * 2 / 10;
        }
        // Feedback delay of a dotted eighth, over the bed.
        const std::int32_t echo = ring[ring_at];
        ring[ring_at] = dry_bed + (echo * 45 / 100);
        ring_at = (ring_at + 1) % delay_len;
        mix[i] = dry_bed + (echo * 40 / 100) + drums;
    }

    // Fold the delay's tail onto the start so the loop joins without a click.
    for (std::size_t i = 0; i < tail; ++i)
        mix[i] += mix[total + i];
    mix.resize(total);

    std::int32_t peak = 1;
    for (std::int32_t v : mix)
        peak = std::max<std::int32_t>(peak, v < 0 ? -v : v);
    const std::int64_t target = static_cast<std::int64_t>(o.mGain * 32767.0);
    std::vector<std::int16_t> out(total);
    for (std::size_t i = 0; i < total; ++i)
        out[i] = static_cast<std::int16_t>(static_cast<std::int64_t>(mix[i]) * target / peak);
    return out;
}

std::string synth_wav(const std::vector<std::int16_t>& rSamples, int SampleRate) {
    auto le32 = [](std::string& rOut, std::uint32_t v) {
        for (int k = 0; k < 4; ++k)
            rOut.push_back(static_cast<char>((v >> (8 * k)) & 0xFF));
    };
    auto le16 = [](std::string& rOut, std::uint16_t v) {
        rOut.push_back(static_cast<char>(v & 0xFF));
        rOut.push_back(static_cast<char>(v >> 8));
    };
    const std::uint32_t bytes = static_cast<std::uint32_t>(rSamples.size() * 2);
    std::string out;
    out.reserve(44 + bytes);
    out += "RIFF";
    le32(out, 36 + bytes);
    out += "WAVEfmt ";
    le32(out, 16);
    le16(out, 1);  // PCM
    le16(out, 1);  // mono
    le32(out, static_cast<std::uint32_t>(SampleRate));
    le32(out, static_cast<std::uint32_t>(SampleRate) * 2);
    le16(out, 2);
    le16(out, 16);
    out += "data";
    le32(out, bytes);
    for (std::int16_t s : rSamples)
        le16(out, static_cast<std::uint16_t>(s));
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
