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
// A ranked format detection report: the content sniffer, the extension table and a
// metadata read per extension candidate. See operations/detect.hpp.

// System includes
#include <algorithm>
#include <cctype>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/detect.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/registry.hpp"

namespace meshioplusplus {

namespace {

/// Extensions several formats register, in the Python registry's order (the twin's
/// `extension_to_filetypes`); every other extension names its registry default only.
struct DetectShared {
    const char* mSuffix;
    std::vector<const char*> mFormats;
};

const std::vector<DetectShared>& detect_shared_extensions() {
    static const std::vector<DetectShared> table = {
        {".inp", {"abaqus", "ansysInp"}}, {".msh", {"ansys", "gmsh", "freefem"}},
        {".dat", {"marc", "tecplot"}},    {".mesh", {"medit", "mfem"}},
        {".ele", {"tetgen", "triangle"}}, {".node", {"tetgen", "triangle"}},
    };
    return table;
}

std::string detect_lower(std::string Text) {
    std::transform(Text.begin(), Text.end(), Text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return Text;
}

/// The extension candidates of @p rPath and the suffix that named them (longest first,
/// as the registry resolves it), or none.
std::pair<std::string, std::vector<std::string>> detect_extension_candidates(
    const std::string& rPath) {
    std::vector<std::string> out;
    std::string named;
    try {
        named = resolve_write_format(rPath, "");
    } catch (const std::exception&) {
        return {};
    }
    const std::string base = detect_lower(std::filesystem::path(rPath).filename().string());
    std::string suffix;
    const auto& defaults = registry_extension_defaults();
    for (std::size_t pos = base.find('.'); pos != std::string::npos;
         pos = base.find('.', pos + 1)) {
        if (defaults.count(base.substr(pos))) {
            suffix = base.substr(pos);
            break;
        }
    }
    for (const DetectShared& shared : detect_shared_extensions())
        if (suffix == shared.mSuffix)
            for (const char* fmt : shared.mFormats)
                if (registry_readers().count(fmt))
                    out.emplace_back(fmt);
    if (std::find(out.begin(), out.end(), named) == out.end())
        out.insert(out.begin(), named);
    return {suffix.empty() ? base : suffix, out};
}

int detect_rank(const std::string& rConfidence) {
    return rConfidence == "certain" ? 0 : rConfidence == "likely" ? 1 : 2;
}

}  // namespace

std::vector<DetectCandidate> detect_format(const std::string& rPath) {
    std::vector<DetectCandidate> out;
    std::string sniffed;
    try {
        sniffed = sniff_format(rPath);
    } catch (const std::exception&) {
        sniffed.clear();
    }
    if (!sniffed.empty())
        out.push_back({sniffed, "certain", "content signature"});
    const auto [suffix, candidates] = detect_extension_candidates(rPath);
    for (const std::string& fmt : candidates) {
        if (fmt == sniffed)
            continue;
        bool reads = false;
        if (registry_readers().count(fmt)) {
            try {
                (void)registry_read_metadata(rPath, fmt, ReadOptions{});
                reads = true;
            } catch (const std::exception&) {
                reads = false;
            }
        }
        if (reads)
            out.push_back({fmt, "likely", "extension " + suffix + " and the file reads as " + fmt});
        else
            out.push_back({fmt, "guess", "extension " + suffix + " only"});
    }
    std::stable_sort(out.begin(), out.end(),
                     [](const DetectCandidate& rA, const DetectCandidate& rB) {
                         return detect_rank(rA.mConfidence) < detect_rank(rB.mConfidence);
                     });
    return out;
}

std::string detect_format_text(const std::string& rPath) {
    std::string text;
    for (const DetectCandidate& c : detect_format(rPath))
        text += c.mFormat + "\t" + c.mConfidence + "\t" + c.mReason + "\n";
    return text;
}

}  // namespace meshioplusplus
