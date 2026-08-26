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

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/gid.hpp"
#include "meshioplusplus/formats/xdmf_time_series.hpp"
#include "meshioplusplus/operations/sniff.hpp"

namespace meshioplusplus {

namespace {

// --------------------------------------------------------------------------
// Every helper here is `seq_`-prefixed: the amalgamation concatenates every
// src/cpp/src/*.cpp into one translation unit, so an anonymous-namespace name
// must be unique across the whole tree.
// --------------------------------------------------------------------------

bool seq_is_digit(char c) {
    return c >= '0' && c <= '9';
}

/// The digits of the run starting at `pos`, with leading zeros stripped, and
/// `pos` advanced past the run. Returns a view into `rS`.
std::string_view seq_digit_run(const std::string& rS, std::size_t& rPos) {
    const std::size_t begin = rPos;
    while (rPos < rS.size() && seq_is_digit(rS[rPos]))
        ++rPos;
    std::size_t first = begin;
    // Strip leading zeros, but keep one digit for a run that is all zeros.
    while (first + 1 < rPos && rS[first] == '0')
        ++first;
    return std::string_view(rS).substr(first, rPos - first);
}

}  // namespace

bool sequence_natural_less(const std::string& rA, const std::string& rB) {
    std::size_t ia = 0;
    std::size_t ib = 0;
    while (ia < rA.size() && ib < rB.size()) {
        const bool da = seq_is_digit(rA[ia]);
        const bool db = seq_is_digit(rB[ib]);
        if (da != db)
            return da;  // rule 4: a digit run sorts before a non-digit run
        if (da) {
            // Rule 3: numeric comparison done on the digits themselves, never
            // through stoull -- a 40-digit hash-named file must not overflow.
            const std::string_view va = seq_digit_run(rA, ia);
            const std::string_view vb = seq_digit_run(rB, ib);
            if (va.size() != vb.size())
                return va.size() < vb.size();
            if (va != vb)
                return va < vb;
        } else {
            // Rule 2: as unsigned char. Plain `char` signedness is
            // implementation-defined, and would order the same UTF-8 paths
            // differently on ARM than on x86.
            const unsigned char ca = static_cast<unsigned char>(rA[ia]);
            const unsigned char cb = static_cast<unsigned char>(rB[ib]);
            if (ca != cb)
                return ca < cb;
            ++ia;
            ++ib;
        }
    }
    if (ia < rA.size() || ib < rB.size())
        return ib < rB.size();  // the exhausted string is the shorter one
    // Rule 5: every run compared equal, so fall back to a plain byte-wise
    // comparison of the ORIGINAL strings. This is what makes the comparator a
    // strict weak ordering: without it `out_1` and `out_01` would be mutually
    // "not less" yet not equivalent, and std::sort would be UB on a directory
    // mixing padded and unpadded names.
    return rA < rB;
}

bool sequence_glob_match(const std::string& rPattern, const std::string& rName) {
    // Iterative backtracking: `star` remembers the last '*' seen and `mark`
    // where in the name it was matched, so a pathological `a*a*a*a*b` cannot
    // blow the stack the way a recursive matcher would.
    std::size_t p = 0;
    std::size_t n = 0;
    std::size_t star = std::string::npos;
    std::size_t mark = 0;
    while (n < rName.size()) {
        if (p < rPattern.size() && (rPattern[p] == '?' || rPattern[p] == rName[n])) {
            ++p;
            ++n;
        } else if (p < rPattern.size() && rPattern[p] == '*') {
            star = p++;
            mark = n;
        } else if (star != std::string::npos) {
            p = star + 1;
            n = ++mark;
        } else {
            return false;
        }
    }
    while (p < rPattern.size() && rPattern[p] == '*')
        ++p;
    return p == rPattern.size();
}

bool sequence_pattern_has_token(const std::string& rPath) {
    return rPath.find("{step}") != std::string::npos || rPath.find("{index}") != std::string::npos;
}

namespace {

/// Replace every occurrence of `rToken` in `rIn` with `rValue`, advancing past
/// each substitution so a value containing the token cannot loop. The native
/// CLI's `replace_key`/`partition_replace_part` shape, deliberately: these
/// tokens inherit the existing substring-replace semantics rather than
/// introducing a third convention. (The Python twin uses `str.format`, which
/// raises on a stray brace; that asymmetry is pre-existing.)
std::string seq_replace_all(const std::string& rIn, const std::string& rToken,
                            const std::string& rValue) {
    std::string out = rIn;
    std::size_t pos = out.find(rToken);
    while (pos != std::string::npos) {
        out.replace(pos, rToken.size(), rValue);
        pos = out.find(rToken, pos + rValue.size());
    }
    return out;
}

}  // namespace

std::string sequence_expand_pattern(const std::string& rPattern, std::size_t Index,
                                    std::size_t Count) {
    const std::string plain = std::to_string(Index);
    // Pad to at least 4, and wider once the count needs it, so a directory
    // listing of a 20000-step run still sorts correctly for a naive tool.
    std::size_t width = 4;
    if (Count > 1) {
        const std::string last = std::to_string(Count - 1);
        width = std::max(width, last.size());
    }
    std::string padded = plain;
    if (padded.size() < width)
        padded.insert(padded.begin(), width - padded.size(), '0');
    return seq_replace_all(seq_replace_all(rPattern, "{step}", padded), "{index}", plain);
}

namespace {

/// resolve_format with the sniff_format fallback -- the read-path rule
/// everywhere (the CLIs, mio_read, the wasm read_mesh, run_pipeline).
std::string seq_resolve_read_format(const std::string& rPath, const std::string& rFormat) {
    try {
        return resolve_format(rPath, rFormat);
    } catch (const ReadError&) {
        const std::string sniffed = sniff_format(rPath);
        if (sniffed.empty())
            throw;
        return sniffed;
    }
}

}  // namespace

namespace {

/// Whether `rFormat`'s reader can return more than one step at all.
///
/// Consulted BEFORE `registry_read_metadata`, which for a format with no
/// native metadata reader costs a full read -- so this is what keeps the step
/// probe free for the 38 formats that cannot carry time. Same shape and same
/// anti-drift discipline as `sequence_write_supports_time`: a small owned set,
/// cross-checked by a gtest against which readers actually honour
/// `ReadOptions::mTimeStep`.
///
/// **MED is deliberately absent.** It honours `ReadOptions::mTimeStep`, but has
/// no entry in `registry_metadata_readers()`, so there is no count to read:
/// probing it would cost a full read and still report one step. That is a
/// recorded gap in MED's metadata support, and it closes here for free the
/// moment `read_med_metadata` fills `mTimeValues`.
bool seq_format_may_have_steps(const std::string& rFormat) {
    // gid joined in v10.19.0: its reader has always honoured mTimeStep, but
    // read_gid_metadata never opened the results sibling where steps live, so
    // it reported one step and this predicate had nothing to gate on.
    return rFormat == "xdmf" || rFormat == "exodus" || rFormat == "gid";
}

}  // namespace

std::size_t sequence_num_steps(const std::string& rPath, const std::string& rFormat) {
    // Registry-derived, never a per-format table: a format whose metadata
    // reader does not fill mTimeValues reports 1, which is the truthful answer
    // for every format that cannot express time. An unreadable file also
    // reports 1 -- the failure belongs to the real read, with its own
    // diagnostics, not to a capability query.
    try {
        const std::string fmt = seq_resolve_read_format(rPath, rFormat);
        if (!seq_format_may_have_steps(fmt))
            return 1u;
        ReadOptions opts;
        opts.mMetadataOnly = true;
        const MeshMetadata meta = registry_read_metadata(rPath, fmt, opts);
        return meta.mTimeValues.empty() ? 1u : meta.mTimeValues.size();
    } catch (const std::exception&) {
        return 1u;
    }
}

bool sequence_write_supports_time(const std::string& rFormat, std::string& rWhy) {
    // The one multi-step writer in the repo. Unlike sequence_num_steps there is
    // no file to probe, so this is a predicate in registry_write_supports'
    // style. It is kept honest by a gtest that cross-checks it against what
    // sequence_to_timeseries actually accepts, over every registry_writers()
    // entry -- a format that grows a series writer without updating this turns
    // CI red naming itself.
    if (rFormat == "xdmf" || rFormat == "gid") {
        rWhy.clear();
        return true;
    }
    rWhy = "meshio++: sequence: format '" + rFormat +
           "' cannot hold a multi-step series (only 'xdmf' and 'gid' can); write one file per "
           "step with an Output path containing '{step}' instead";
    return false;
}

namespace {

/// The last maximal digit run of a path's stem, as a double. Returns false when
/// the stem has no digits at all.
bool seq_time_from_filename(const std::string& rPath, double& rOut) {
    std::string stem;
    try {
        stem = std::filesystem::path(rPath).stem().string();
    } catch (const std::exception&) {
        stem = rPath;
    }
    // Walk backwards to the last run: `run17/out_0042.vtu` must give 42, not 17
    // (and the directory is already excluded by taking the stem).
    std::size_t end = stem.size();
    while (end > 0 && !seq_is_digit(stem[end - 1]))
        --end;
    if (end == 0)
        return false;
    std::size_t begin = end;
    while (begin > 0 && seq_is_digit(stem[begin - 1]))
        --begin;
    try {
        rOut = static_cast<double>(std::stoll(stem.substr(begin, end - begin)));
    } catch (const std::exception&) {
        // A run too long for long long: the value is not a plausible time
        // anyway, so decline rather than saturate.
        return false;
    }
    return true;
}

/// A mesh's own single time value: `meshio:time`, or a length-1 `exodus:time`
/// (the convention this generalizes -- "the one value this mesh is a snapshot
/// at"). A longer `exodus:time` is the writer's whole-history vector and is
/// deliberately NOT read as a per-file scalar.
bool seq_time_from_mesh(const Mesh& rMesh, double& rOut) {
    for (const char* key : {kSequenceTimeKey, "exodus:time"}) {
        if (!rMesh.HasFieldData(key))
            continue;
        const NDArray& a = rMesh.FieldData(key);
        if (a.Size() != 1)
            continue;
        rOut = detail::read_double(a, 0);
        return true;
    }
    return false;
}

/// Split a glob into (directory, basename pattern). The directory component is
/// taken literally: a `*` there is an error naming the restriction, rather than
/// a recursive walk nobody asked for.
void seq_split_pattern(const std::string& rPattern, std::string& rDir, std::string& rBase) {
    const std::size_t slash = rPattern.find_last_of("/\\");
    if (slash == std::string::npos) {
        rDir = ".";
        rBase = rPattern;
    } else {
        rDir = rPattern.substr(0, slash);
        rBase = rPattern.substr(slash + 1);
        if (rDir.empty())
            rDir = "/";
    }
    if (rDir.find('*') != std::string::npos || rDir.find('?') != std::string::npos)
        throw std::invalid_argument(
            "meshio++: sequence: the directory part of a pattern is taken literally, so '" + rDir +
            "' cannot contain '*' or '?'; glob one directory at a time");
}

std::vector<std::string> seq_glob(const std::string& rPattern) {
    std::string dir;
    std::string base;
    seq_split_pattern(rPattern, dir, base);

    std::error_code ec;
    std::vector<std::string> out;
    std::filesystem::directory_iterator it(dir, ec);
    if (ec)
        throw ReadError("meshio++: sequence: cannot list directory '" + dir + "' for pattern '" +
                        rPattern + "': " + ec.message());
    for (const std::filesystem::directory_entry& entry : it) {
        if (!entry.is_regular_file(ec))
            continue;
        const std::string name = entry.path().filename().string();
        if (sequence_glob_match(base, name))
            out.push_back(entry.path().string());
    }
    std::sort(out.begin(), out.end(), sequence_natural_less);
    if (out.empty())
        throw ReadError("meshio++: sequence: pattern '" + rPattern + "' matched no files");
    return out;
}

}  // namespace

std::vector<SequenceEntry> sequence_expand(const SequenceInput& rInput) {
    const bool has_paths = !rInput.mPaths.empty();
    const bool has_pattern = !rInput.mPattern.empty();
    if (has_paths == has_pattern)
        throw std::invalid_argument(
            has_paths ? "meshio++: sequence: set either an explicit path list or a pattern, "
                        "not both"
                      : "meshio++: sequence: no input; set an explicit path list or a pattern");

    std::vector<std::string> files = has_pattern ? seq_glob(rInput.mPattern) : rInput.mPaths;
    if (has_paths && rInput.mSortExplicit)
        std::sort(files.begin(), files.end(), sequence_natural_less);

    // One file may carry several steps, so the entry list is not the file list.
    // This is the only place a metadata read happens, and only when a time
    // could come from a multi-step file at all.
    std::vector<SequenceEntry> entries;
    for (const std::string& path : files) {
        if (!std::filesystem::exists(path))
            throw ReadError("meshio++: sequence: input file not found: '" + path + "'");
        std::vector<double> times;
        std::size_t steps = 1;
        if (rInput.mTimeFrom == SequenceTimeFrom::Auto ||
            rInput.mTimeFrom == SequenceTimeFrom::File) {
            try {
                const std::string fmt = seq_resolve_read_format(path, rInput.mFormat);
                // Gated: for a format that cannot carry time this metadata read
                // would be a full read, and a plain single-file run must not
                // pay for a probe whose answer is known to be 1.
                if (seq_format_may_have_steps(fmt)) {
                    ReadOptions meta_opts;
                    meta_opts.mMetadataOnly = true;
                    times = registry_read_metadata(path, fmt, meta_opts).mTimeValues;
                }
            } catch (const std::exception&) {
                times.clear();  // no metadata reader, or an unreadable file
            }
            steps = times.empty() ? 1u : times.size();
        }
        for (std::size_t k = 0; k < steps; ++k) {
            SequenceEntry e;
            e.mPath = path;
            e.mStep = k;
            if (k < times.size()) {
                e.mTime = times[k];
                e.mTimeSource = SequenceTimeSource::File;
            }
            entries.push_back(std::move(e));
        }
    }

    // Precedence, applied over the whole expanded list.
    if (!rInput.mTimes.empty()) {
        if (rInput.mTimes.size() != entries.size())
            throw std::invalid_argument(
                "meshio++: sequence: " + std::to_string(rInput.mTimes.size()) +
                " explicit time value(s) for " + std::to_string(entries.size()) +
                " sequence entr(ies); the counts must match");
        for (std::size_t i = 0; i < entries.size(); ++i) {
            entries[i].mTime = rInput.mTimes[i];
            entries[i].mTimeSource = SequenceTimeSource::Explicit;
        }
        return entries;
    }

    for (std::size_t i = 0; i < entries.size(); ++i) {
        SequenceEntry& e = entries[i];
        if (rInput.mTimeFrom == SequenceTimeFrom::Index) {
            e.mTime = static_cast<double>(i);
            e.mTimeSource = SequenceTimeSource::Index;
            continue;
        }
        if (e.mTimeSource == SequenceTimeSource::File)
            continue;  // a multi-step file already told us
        double t = 0.0;
        if (rInput.mTimeFrom != SequenceTimeFrom::File && seq_time_from_filename(e.mPath, t)) {
            e.mTime = t;
            e.mTimeSource = SequenceTimeSource::Filename;
            continue;
        }
        // Provisional. A single-step file's `meshio:time` needs a real read, so
        // the drivers upgrade this to `File` from the mesh they were going to
        // read anyway -- no file is ever read twice just to find its time.
        e.mTime = static_cast<double>(i);
        e.mTimeSource = SequenceTimeSource::Index;
    }
    return entries;
}

Mesh sequence_read_step(const std::vector<SequenceEntry>& rEntries, std::size_t Index,
                        const std::string& rFormat, const ReadOptions& rOptions) {
    if (Index >= rEntries.size())
        throw std::out_of_range("meshio++: sequence: step index " + std::to_string(Index) +
                                " is past the end of a " + std::to_string(rEntries.size()) +
                                "-step sequence");
    const SequenceEntry& e = rEntries[Index];
    const std::string fmt = seq_resolve_read_format(e.mPath, rFormat);
    ReadOptions opts = rOptions;
    // The entry's own step wins: that is what makes fan-out work at all.
    opts.mTimeStep = static_cast<int>(e.mStep);
    return registry_read(e.mPath, fmt, opts);
}

namespace {

/// Attach `field_data["meshio:time"]`, so fan-out -> fan-in closes on time as
/// well as on geometry.
void seq_attach_time(Mesh& rMesh, double Time) {
    NDArray t(DType::Float64, {1});
    *reinterpret_cast<double*>(t.Data()) = Time;
    rMesh.AddFieldData(kSequenceTimeKey, std::move(t));
}

/// The write format of a sequence output, resolved once. For a pattern the
/// extension of the *expanded* first name is what matters, since `{step}`
/// sits before it.
std::string seq_resolve_write_format(const SequenceOutput& rOutput, std::size_t Count) {
    if (!rOutput.mFormat.empty())
        return rOutput.mFormat;
    const std::string probe = sequence_pattern_has_token(rOutput.mPath)
                                  ? sequence_expand_pattern(rOutput.mPath, 0, Count)
                                  : rOutput.mPath;
    return resolve_format(probe, "");
}

/// The transient writer's data-format choice. An explicit request always
/// wins -- `Binary` asks for `"HDF"` and lets `XdmfTimeSeriesWriter`'s own
/// constructor throw naming the missing build flag when this core has no
/// HDF5, which is the existing "reject what cannot be honoured" rule.  Left
/// at `Default`, this follows the build exactly like the registry's own
/// xdmf entry does (`registry.cpp`: "HDF when HDF5 is available, XML
/// otherwise") -- an HDF-format series is unreadable by a core with no
/// HDF5 support, including the very core that would have just written it.
std::string seq_resolve_data_format(const WriteOptions& rOptions) {
    if (rOptions.mEncoding == WriteEncoding::Ascii)
        return "XML";
    if (rOptions.mEncoding == WriteEncoding::Binary)
        return "HDF";
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    return "HDF";
#else
    return "XML";
#endif
}

/// The transient writer bypasses `registry_write_ex` (it drives
/// `XdmfTimeSeriesWriter` directly, not a `(path, mesh)` registry entry), so
/// it must apply the write_options.hpp rule -- "an option the writer cannot
/// honour is an error, never silently ignored" -- itself. Only `Encoding`
/// (XML vs HDF) has anywhere to go; `Codec`/`FloatFormat` do not.
void seq_check_series_write_options(const WriteOptions& rOptions) {
    if (rOptions.mCodecSet)
        throw WriteError("meshio++: sequence: the transient XDMF writer does not support Codec");
    if (!rOptions.mFloatFormat.empty())
        throw WriteError(
            "meshio++: sequence: the transient XDMF writer does not support FloatFormat");
}

}  // namespace

void sequence_to_timeseries(const SequenceInput& rInput, const SequenceOutput& rOutput) {
    if (rOutput.mPath.empty())
        throw std::invalid_argument("meshio++: sequence: an output path is required");
    const std::vector<SequenceEntry> entries = sequence_expand(rInput);

    const std::string ofmt = seq_resolve_write_format(rOutput, entries.size());
    std::string why;
    if (!sequence_write_supports_time(ofmt, why))
        throw WriteError(why);
    seq_check_series_write_options(rOutput.mOptions);

    // Streaming: one mesh enters scope per iteration and leaves it. There is
    // deliberately no std::vector<Mesh> anywhere in this file.
    if (ofmt == "gid") {
        // gid's series writer is a free function taking a provider rather than
        // a stateful class -- the transient surface meshio++ exposes IS this
        // layer, which already carries every binding, so it needs none of its
        // own. The lambda keeps the same one-mesh-alive property the loop
        // below has.
        write_gid_series(rOutput.mPath, [&](std::size_t i, double& rTime, Mesh& rMesh) {
            if (i >= entries.size())
                return false;
            rMesh = sequence_read_step(entries, i, rInput.mFormat, rInput.mOptions);
            rTime = entries[i].mTime;
            if (entries[i].mTimeSource == SequenceTimeSource::Index &&
                rInput.mTimeFrom != SequenceTimeFrom::Index) {
                double t = 0.0;
                if (seq_time_from_mesh(rMesh, t))
                    rTime = t;
            }
            return true;
        });
        return;
    }

    XdmfTimeSeriesWriter writer(rOutput.mPath, seq_resolve_data_format(rOutput.mOptions));
    bool grid_written = false;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        Mesh mesh = sequence_read_step(entries, i, rInput.mFormat, rInput.mOptions);
        double time = entries[i].mTime;
        if (entries[i].mTimeSource == SequenceTimeSource::Index &&
            rInput.mTimeFrom != SequenceTimeFrom::Index) {
            // The provisional-index upgrade: the mesh is already in hand, so
            // reading its own `meshio:time` costs nothing extra.
            double t = 0.0;
            if (seq_time_from_mesh(mesh, t))
                time = t;
        }
        if (!grid_written) {
            writer.WritePointsCells(mesh);
            grid_written = true;
        }
        writer.WriteData(time, mesh);
    }
    writer.Finalize();
}

void timeseries_to_sequence(const std::string& rInPath, const std::string& rInFormat,
                            const ReadOptions& rOptions, const SequenceOutput& rOutput) {
    if (!sequence_pattern_has_token(rOutput.mPath))
        throw std::invalid_argument(
            "meshio++: sequence: fanning a multi-step file out needs an output pattern "
            "containing '{step}' or '{index}' (e.g. 'out_{step}.vtu'), not '" +
            rOutput.mPath + "'");

    SequenceInput in;
    in.mPaths = {rInPath};
    in.mFormat = rInFormat;
    in.mOptions = rOptions;
    const std::vector<SequenceEntry> entries = sequence_expand(in);

    for (std::size_t i = 0; i < entries.size(); ++i) {
        Mesh mesh = sequence_read_step(entries, i, rInFormat, rOptions);
        seq_attach_time(mesh, entries[i].mTime);
        registry_write_ex(sequence_expand_pattern(rOutput.mPath, i, entries.size()), mesh,
                          rOutput.mFormat, rOutput.mOptions);
    }
}

bool sequence_input_needs_driver(const SequenceInput& rInput, const SequenceOutput& rOutput) {
    if (sequence_pattern_has_token(rOutput.mPath))
        return true;
    // An explicit step request is a deliberate single-step selection.
    if (rInput.mOptions.mTimeStep != 0)
        return false;
    if (!rInput.mPattern.empty() || rInput.mPaths.size() != 1)
        return true;  // a pattern or a list is a sequence by construction
    return sequence_num_steps(rInput.mPaths[0], rInput.mFormat) > 1;
}

SequenceMode sequence_resolve_mode(const std::vector<SequenceEntry>& rEntries,
                                   const SequenceOutput& rOutput, SequenceMode Stated) {
    const bool token = sequence_pattern_has_token(rOutput.mPath);
    // "One file with several steps" and "several files" are different inputs
    // and produce different errors, so count the distinct paths.
    std::size_t files = 0;
    for (std::size_t i = 0; i < rEntries.size(); ++i)
        if (i == 0 || rEntries[i].mPath != rEntries[i - 1].mPath)
            ++files;

    SequenceMode inferred;
    if (token) {
        inferred =
            (files == 1 && rEntries.size() > 1) ? SequenceMode::FanOut : SequenceMode::Sequence;
    } else if (rEntries.size() > 1) {
        if (files == 1)
            throw std::invalid_argument(
                "meshio++: sequence: the input has " + std::to_string(rEntries.size()) +
                " time steps but the output path '" + rOutput.mPath +
                "' names a single file; add '{step}' to write one file per step, or select "
                "a single step (--time-step on the CLI, Input.Options.TimeStep in a "
                "settings document)");
        inferred = SequenceMode::FanIn;
    } else {
        inferred = SequenceMode::Sequence;
    }

    if (Stated != SequenceMode::Auto && Stated != inferred)
        throw std::invalid_argument(
            std::string("meshio++: sequence: Mode says '") + sequence_mode_name(Stated) +
            "' but the input (" + std::to_string(files) + " file(s), " +
            std::to_string(rEntries.size()) + " step(s)) and output '" + rOutput.mPath +
            "' describe '" + sequence_mode_name(inferred) + "'");
    return inferred;
}

const char* sequence_time_source_name(SequenceTimeSource Source) {
    switch (Source) {
        case SequenceTimeSource::Explicit:
            return "explicit";
        case SequenceTimeSource::File:
            return "file";
        case SequenceTimeSource::Filename:
            return "filename";
        case SequenceTimeSource::Index:
            return "index";
    }
    return "index";
}

const char* sequence_mode_name(SequenceMode Mode) {
    switch (Mode) {
        case SequenceMode::Sequence:
            return "sequence";
        case SequenceMode::FanIn:
            return "fan-in";
        case SequenceMode::FanOut:
            return "fan-out";
        case SequenceMode::Auto:
            break;
    }
    return "auto";
}

SequenceMode sequence_mode_from_name(const std::string& rName) {
    // Lowercase, hyphenated -- the repo-wide rule that settings *keys* are
    // PascalCase but enum *values* keep their `*_from_name` spelling, as in
    // gradient's "green-gauss" and convert_cells' "simplexify".
    if (rName.empty() || rName == "auto")
        return SequenceMode::Auto;
    if (rName == "sequence")
        return SequenceMode::Sequence;
    if (rName == "fan-in")
        return SequenceMode::FanIn;
    if (rName == "fan-out")
        return SequenceMode::FanOut;
    throw std::invalid_argument(
        "meshio++: sequence: Mode must be 'sequence', 'fan-in' or 'fan-out', not '" + rName + "'");
}

SequenceTimeFrom sequence_time_from_name(const std::string& rName) {
    if (rName.empty() || rName == "auto")
        return SequenceTimeFrom::Auto;
    if (rName == "file")
        return SequenceTimeFrom::File;
    if (rName == "filename")
        return SequenceTimeFrom::Filename;
    if (rName == "index")
        return SequenceTimeFrom::Index;
    throw std::invalid_argument(
        "meshio++: sequence: Input.TimeFrom must be 'auto', 'file', 'filename' or 'index', "
        "not '" +
        rName + "'");
}

PipelineReport run_sequence_pipeline(const SequencePipeline& rPipeline) {
    if (rPipeline.mVersion != 1)
        throw std::invalid_argument("meshio++: sequence: unsupported Version " +
                                    std::to_string(rPipeline.mVersion) + " (this build knows 1)");
    if (rPipeline.mOutput.mPath.empty())
        throw std::invalid_argument("meshio++: sequence: Output.Path is required");

    // Validate the whole chain before anything is read -- a typo in step 7 must
    // not cost reading 500 files first.
    for (const PipelineStep& step : rPipeline.mSteps)
        validate_pipeline_step(step);

    const std::vector<SequenceEntry> entries = sequence_expand(rPipeline.mInput);
    const SequenceMode mode = sequence_resolve_mode(entries, rPipeline.mOutput, rPipeline.mMode);

    PipelineReport report;
    std::size_t index_fallbacks = 0;
    for (const SequenceEntry& e : entries)
        if (e.mTimeSource == SequenceTimeSource::Index)
            ++index_fallbacks;
    if (index_fallbacks > 0 && rPipeline.mInput.mTimeFrom != SequenceTimeFrom::Index)
        report.mWarnings.push_back(
            "sequence: no time value found for " + std::to_string(index_fallbacks) + " of " +
            std::to_string(entries.size()) + " entr(ies); using the integer index for those");
    if (rPipeline.mParallel)
        report.mWarnings.push_back(
            "sequence: Parallel is a Python-driver feature; this engine runs the steps "
            "serially (every operation already parallelizes internally)");

    if (mode == SequenceMode::FanIn) {
        const std::string ofmt = seq_resolve_write_format(rPipeline.mOutput, entries.size());
        std::string why;
        if (!sequence_write_supports_time(ofmt, why))
            throw WriteError(why);
        seq_check_series_write_options(rPipeline.mOutput.mOptions);
        XdmfTimeSeriesWriter writer(rPipeline.mOutput.mPath,
                                    seq_resolve_data_format(rPipeline.mOutput.mOptions));
        bool grid_written = false;
        for (std::size_t i = 0; i < entries.size(); ++i) {
            Mesh mesh =
                sequence_read_step(entries, i, rPipeline.mInput.mFormat, rPipeline.mInput.mOptions);
            double time = entries[i].mTime;
            if (entries[i].mTimeSource == SequenceTimeSource::Index &&
                rPipeline.mInput.mTimeFrom != SequenceTimeFrom::Index) {
                double t = 0.0;
                if (seq_time_from_mesh(mesh, t))
                    time = t;
            }
            // The single owner of step dispatch, unchanged: sequences are a
            // driver AROUND run_pipeline_steps, never a second dispatch path.
            mesh = run_pipeline_steps(std::move(mesh), rPipeline.mSteps, report);
            if (!grid_written) {
                writer.WritePointsCells(mesh);
                grid_written = true;
            }
            writer.WriteData(time, mesh);
        }
        writer.Finalize();
        return report;
    }

    // Sequence / FanOut: both write one file per entry, and differ only in
    // where the entries came from -- which sequence_expand already settled.
    for (std::size_t i = 0; i < entries.size(); ++i) {
        Mesh mesh =
            sequence_read_step(entries, i, rPipeline.mInput.mFormat, rPipeline.mInput.mOptions);
        double time = entries[i].mTime;
        if (entries[i].mTimeSource == SequenceTimeSource::Index &&
            rPipeline.mInput.mTimeFrom != SequenceTimeFrom::Index) {
            double t = 0.0;
            if (seq_time_from_mesh(mesh, t))
                time = t;
        }
        mesh = run_pipeline_steps(std::move(mesh), rPipeline.mSteps, report);
        // `meshio:time` labels a file that is one step OF A SERIES, so it is
        // attached only when the output is a pattern. A document with a single
        // plain output path must produce byte-identical output to the
        // single-file pipeline it degenerates to.
        const bool token = sequence_pattern_has_token(rPipeline.mOutput.mPath);
        if (token)
            seq_attach_time(mesh, time);
        const std::string path =
            token ? sequence_expand_pattern(rPipeline.mOutput.mPath, i, entries.size())
                  : rPipeline.mOutput.mPath;
        registry_write_ex(path, mesh, rPipeline.mOutput.mFormat, rPipeline.mOutput.mOptions);
    }
    return report;
}

}  // namespace meshioplusplus
