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
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <mutex>
#include <string_view>

// Project includes
#include "meshioplusplus/detail/provenance.hpp"

#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

/// Thread-local active state. A raw pointer (never owning) so `ProvenanceScope`
/// -- which does own the record -- can point it at its own storage while it is
/// alive and null it on destruction; that is what makes "no scope open" free
/// to check from every writer without touching a mutex.
thread_local ProvenanceRecord* g_active_record = nullptr;
thread_local ProvenanceMode g_active_mode = ProvenanceMode::Off;

const ProvenanceRecord& empty_record() {
    static const ProvenanceRecord r;
    return r;
}

}  // namespace

ProvenanceMode current_provenance_mode() {
    return g_active_mode;
}

const ProvenanceRecord& current_provenance() {
    return g_active_record ? *g_active_record : empty_record();
}

ProvenanceScope::ProvenanceScope(ProvenanceMode mode, ProvenanceRecord record)
    : mHadPrevious(g_active_record != nullptr),
      mPrevious(mHadPrevious ? *g_active_record : ProvenanceRecord{}),
      mPreviousMode(g_active_mode) {
    // Stamp the timestamp now, once, rather than lazily at render time: a
    // caller reading Get() mid-write should see the same instant the file
    // will (if it records one at all), and a lazy stamp taken once per
    // render call would let a single write's block disagree with itself
    // across formats sharing one scope (never happens today, but "the
    // timestamp is fixed for the life of the scope" is the simpler contract
    // to keep, and it costs one clock read).
    if (record.mTimestamp.empty())
        record.mTimestamp = provenance_timestamp();
    // The record this scope owns lives for the scope's lifetime, heap-backed
    // so the thread-local pointer stays valid across the constructor's
    // return; every note()/set_* call afterwards mutates it in place rather
    // than copying.
    g_active_record = new ProvenanceRecord(std::move(record));
    g_active_mode = mode;
}

ProvenanceScope::~ProvenanceScope() {
    delete g_active_record;
    if (mHadPrevious) {
        g_active_record = new ProvenanceRecord(std::move(mPrevious));
    } else {
        g_active_record = nullptr;
    }
    g_active_mode = mPreviousMode;
}

const ProvenanceRecord& ProvenanceScope::Get() const {
    return current_provenance();
}

void provenance_note(std::string_view category, std::string_view detail) {
    if (!g_active_record)
        return;
    for (const auto& n : g_active_record->mNotes)
        if (n.mCategory == category && n.mDetail == detail)
            return;  // collapse duplicates -- a per-cell warning must not
                     // produce a per-cell record.
    g_active_record->mNotes.push_back(ProvenanceNote{std::string(category), std::string(detail)});
}

void provenance_set_source(std::string_view path, std::string_view format) {
    if (!g_active_record)
        return;
    g_active_record->mSourcePath = std::string(path);
    g_active_record->mSourceFormat = std::string(format);
}

void provenance_set_target(std::string_view format, std::string_view encoding,
                           std::string_view codec, std::string_view float_format) {
    if (!g_active_record)
        return;
    g_active_record->mTargetFormat = std::string(format);
    g_active_record->mEncoding = std::string(encoding);
    g_active_record->mCodec = std::string(codec);
    g_active_record->mFloatFormat = std::string(float_format);
}

void provenance_add_operation(std::string_view rendered) {
    if (!g_active_record)
        return;
    g_active_record->mOperations.emplace_back(rendered);
}

namespace {

std::string render_block(const ProvenanceRecord& rRecord) {
    std::vector<std::string> lines;
    if (!rRecord.mSourcePath.empty() || !rRecord.mSourceFormat.empty()) {
        std::string line = "Converted from ";
        line += rRecord.mSourcePath.empty() ? "(in-memory)" : rRecord.mSourcePath;
        if (!rRecord.mSourceFormat.empty())
            line += " (" + rRecord.mSourceFormat + ")";
        lines.push_back(std::move(line));
    }
    if (!rRecord.mTargetFormat.empty()) {
        std::string line = "Written as " + rRecord.mTargetFormat;
        std::vector<std::string> extras;
        if (!rRecord.mEncoding.empty())
            extras.push_back(rRecord.mEncoding);
        if (!rRecord.mCodec.empty())
            extras.push_back("codec=" + rRecord.mCodec);
        if (!rRecord.mFloatFormat.empty())
            extras.push_back("float_format=" + rRecord.mFloatFormat);
        if (!extras.empty()) {
            line += " (";
            for (std::size_t i = 0; i < extras.size(); ++i) {
                if (i)
                    line += ", ";
                line += extras[i];
            }
            line += ")";
        }
        lines.push_back(std::move(line));
    }
    for (const auto& op : rRecord.mOperations)
        lines.push_back("Operation: " + op);
    for (const auto& note : rRecord.mNotes)
        lines.push_back("Note [" + note.mCategory + "]: " + note.mDetail);
    if (!rRecord.mTimestamp.empty())
        lines.push_back("Timestamp: " + rRecord.mTimestamp);

    std::string out;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        if (i)
            out += '\n';
        out += lines[i];
    }
    return out;
}

}  // namespace

std::vector<std::string> provenance_lines(SlotTier tier) {
    std::vector<std::string> out{kProvenanceTag};

    const ProvenanceMode mode = g_active_record ? g_active_mode : ProvenanceMode::Off;
    if (mode == ProvenanceMode::Off)
        return out;

    if (tier == SlotTier::None) {
        if (mode == ProvenanceMode::Required)
            throw WriteError(
                "meshio++: provenance: this format's header slot cannot hold any "
                "provenance -- Mode::Required cannot be honoured here");
        return out;  // BestEffort on a slot-less format: nothing to add.
    }

    if (tier != SlotTier::Block)
        return out;  // SingleLine/Bounded: the tag is the honest maximum, and
                     // Required does not treat that as a failure (see the
                     // header doc comment).

    const ProvenanceRecord& record = current_provenance();
    std::string block = render_block(record);
    if (block.empty())
        return out;  // nothing was ever recorded beyond the tag itself
    // Split render_block's newline-joined text back into individual lines so
    // each one can be comment-prefixed independently by the caller.
    std::size_t start = 0;
    while (start <= block.size()) {
        std::size_t nl = block.find('\n', start);
        if (nl == std::string::npos) {
            out.push_back(block.substr(start));
            break;
        }
        out.push_back(block.substr(start, nl - start));
        start = nl + 1;
    }
    return out;
}

std::string provenance_render_lines(SlotTier tier, std::string_view prefix) {
    std::string out;
    for (const auto& line : provenance_lines(tier)) {
        out += prefix;
        out += line;
        out += '\n';
    }
    return out;
}

std::string provenance_render_xml_comment(SlotTier tier) {
    auto lines = provenance_lines(tier);
    if (lines.size() == 1)
        return "<!--" + lines[0] + "-->";
    std::string out = "<!--\n";
    for (const auto& line : lines) {
        out += line;
        out += '\n';
    }
    out += "-->";
    return out;
}

std::string provenance_timestamp() {
    if (const char* off = std::getenv("MESHIOPLUSPLUS_PROVENANCE_TIMESTAMP");
        off != nullptr && std::string_view(off) == "off") {
        // SOURCE_DATE_EPOCH, when also set, is a more specific ask than the
        // blanket off switch and wins -- a reproducible-build pipeline that
        // exports both wants the pinned epoch, not silence.
        if (std::getenv("SOURCE_DATE_EPOCH") == nullptr)
            return {};
    }

    std::time_t t;
    if (const char* epoch = std::getenv("SOURCE_DATE_EPOCH")) {
        char* end = nullptr;
        long long v = std::strtoll(epoch, &end, 10);
        t = (end != epoch) ? static_cast<std::time_t>(v) : std::time(nullptr);
    } else {
        t = std::time(nullptr);
    }

    std::tm tm_utc{};
#if defined(_WIN32)
    gmtime_s(&tm_utc, &t);
#else
    gmtime_r(&t, &tm_utc);
#endif
    char buf[32];
    std::size_t n = std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm_utc);
    return std::string(buf, n);
}

namespace {

/// The line shapes `render_block` emits, plus the tag itself. A line counts as
/// provenance iff it starts with one of these once its comment punctuation is
/// stripped -- which is what makes one scanner serve every format.
constexpr const char* kProvenancePrefixes[] = {
    "Written by meshio++ v", "Converted from ", "Written as ",
    "Operation: ",           "Note [",          "Timestamp: ",
};

bool prov_starts_with(std::string_view s, std::string_view prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

/// Trims ASCII whitespace from both ends.
std::string_view prov_trim(std::string_view s) {
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t' || s.front() == '\r'))
        s.remove_prefix(1);
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r'))
        s.remove_suffix(1);
    return s;
}

bool prov_is_provenance_line(std::string_view s) {
    for (const char* p : kProvenancePrefixes)
        if (prov_starts_with(s, p))
            return true;
    return false;
}

/// Strips one leading comment marker, if present.
///
/// Deliberately tolerant rather than per-format: the marker sets across the 44
/// formats are small and none of them collides with the content prefixes
/// above, so stripping any of them can only expose a real provenance line or
/// leave a non-matching one to be rejected by the prefix test anyway.
std::string_view prov_strip_marker(std::string_view line) {
    line = prov_trim(line);
    if (prov_starts_with(line, "<!--"))
        line.remove_prefix(4);
    else if (prov_starts_with(line, "//"))
        line.remove_prefix(2);
    else if (prov_starts_with(line, "comment "))
        line.remove_prefix(8);
    else if (!line.empty() && (line.front() == '#' || line.front() == '!' || line.front() == '*' ||
                               line.front() == '$' || line.front() == '%'))
        line.remove_prefix(1);
    return prov_trim(line);
}

/// Yields every substring of @p rRaw that could be a provenance line, best
/// candidate first, and calls @p rFn with each until one is accepted.
///
/// Three shapes have to be tried rather than one, and each earns its place:
/// the plain comment-stripped line (almost every format); the body of a
/// keyword slot that *wraps* the text in a quote or paren (Tecplot's
/// `TITLE = "..."`, Ansys's `(1 "...")`) -- where the closing punctuation must
/// be removed only because this opener put it there, never unconditionally,
/// since `Converted from x (fmt)` and `Operation: Clean(Weld=true)` both end
/// in a legitimate `)`; and each `|`-delimited cell of a box-drawn banner
/// (OpenFOAM), whose credit sits in an interior cell rather than at the start
/// of the line.
template <class F>
void prov_for_each_candidate(std::string_view rRaw, F&& rFn) {
    std::string_view stripped = prov_strip_marker(rRaw);
    if (rFn(stripped))
        return;

    for (std::string_view lead : {std::string_view("TITLE = \""), std::string_view("TITLE=\""),
                                  std::string_view("(1 \"")}) {
        if (!prov_starts_with(stripped, lead))
            continue;
        std::string_view body = stripped.substr(lead.size());
        // Drop exactly the closer this opener implies, nothing more.
        while (!body.empty() && (body.back() == ')' || body.back() == '"'))
            body.remove_suffix(1);
        if (rFn(prov_trim(body)))
            return;
    }

    if (rRaw.find('|') == std::string_view::npos)
        return;
    std::size_t pos = 0;
    while (pos <= rRaw.size()) {
        std::size_t bar = rRaw.find('|', pos);
        std::string_view cell =
            rRaw.substr(pos, bar == std::string_view::npos ? std::string_view::npos : bar - pos);
        if (rFn(prov_trim(cell)))
            return;
        if (bar == std::string_view::npos)
            break;
        pos = bar + 1;
    }
}

}  // namespace

ProvenanceReadResult scan_provenance_text(std::string_view text) {
    ProvenanceReadResult out;
    std::size_t pos = 0;
    while (pos <= text.size()) {
        std::size_t nl = text.find('\n', pos);
        std::string_view raw =
            text.substr(pos, nl == std::string_view::npos ? std::string_view::npos : nl - pos);

        prov_for_each_candidate(raw, [&](std::string_view candidate) {
            if (!prov_is_provenance_line(candidate))
                return false;
            // An XML block's closing "-->" can share the last content line.
            while (candidate.size() >= 3 && candidate.compare(candidate.size() - 3, 3, "-->") == 0)
                candidate = prov_trim(candidate.substr(0, candidate.size() - 3));
            out.mLines.emplace_back(candidate);
            return true;
        });

        if (nl == std::string_view::npos)
            break;
        pos = nl + 1;
    }
    out.mRecognised =
        !out.mLines.empty() && prov_starts_with(out.mLines[0], "Written by meshio++ v");
    return out;
}

ProvenanceReadResult read_provenance_lines(const std::string& rPath, std::size_t max_bytes) {
    // Best-effort: an unopenable path is "nothing found", never a throw. This
    // enriches a summary; it must not be able to fail one.
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        return {};
    std::string head(max_bytes, '\0');
    in.read(head.data(), static_cast<std::streamsize>(max_bytes));
    head.resize(static_cast<std::size_t>(in.gcount()));

    // Binary slots (STL's 80-byte header, EnSight's str80 records) NUL-pad
    // their text. Treating NUL as a line break is what lets the same
    // line-oriented scan find those without a separate code path.
    for (char& c : head)
        if (c == '\0')
            c = '\n';

    return scan_provenance_text(head);
}

}  // namespace detail
}  // namespace meshioplusplus
