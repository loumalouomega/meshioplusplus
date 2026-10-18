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
#include <array>
#include <charconv>
#include <bit>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <deque>
#include <memory>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/openfoam.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/ragged_csr.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "face_cells_common.hpp"
#include "../detail/row_writer.hpp"
#include "../detail/text_cursor.hpp"
#include "../detail/typed_view.hpp"

namespace fs = std::filesystem;

namespace meshioplusplus {

namespace {

using Face = std::vector<std::int64_t>;

struct FoamFormat {
    bool mBinary = false;
    int mLabelBytes = 8;
    int mScalarBytes = 8;
};

/**
 * @brief Whole-file access, mapped where that pays (detail/file_source.hpp).
 *
 * This replaces an `ostringstream` + `.str()` slurp, which paid for **two**
 * extra full-file copies on top of the read -- by far the worst of the
 * whole-file readers, and the reason this one benefits most from mapping.
 * Returns the source itself so the caller controls its lifetime; everything
 * below takes a view into it.
 */
detail::FileSource read_whole(const std::string& rPath) {
    try {
        return detail::FileSource(rPath);
    } catch (const ReadError&) {
        throw ReadError("Could not open OpenFOAM file: " + rPath);
    }
}

std::string openfoam_strip(const std::string& rS) {
    std::size_t a = rS.find_first_not_of(" \t\r\n");
    if (a == std::string::npos)
        return "";
    std::size_t b = rS.find_last_not_of(" \t\r\n");
    return rS.substr(a, b - a + 1);
}

// `openfoam_strip` over a view: the same " \t\r\n" set, no copy.
std::string_view foam_trim(std::string_view S) {
    const std::size_t a = S.find_first_not_of(" \t\r\n");
    if (a == std::string_view::npos)
        return {};
    const std::size_t b = S.find_last_not_of(" \t\r\n");
    return S.substr(a, b - a + 1);
}

// Parse the FoamFile header for format/arch (label/scalar byte widths) from the
// file's text, which the caller keeps mapped.
FoamFormat detect_format(std::string_view rRaw) {
    FoamFormat fmt;
    detail::TextCursor lines(rRaw);
    while (!lines.AtEnd()) {
        const std::string_view s = foam_trim(lines.Line());
        // format <word>;
        if (s.starts_with("format")) {
            std::string_view rest = foam_trim(s.substr(6));
            if (!rest.empty() && rest.back() == ';')
                rest.remove_suffix(1);
            rest = foam_trim(rest);
            if (rest == "binary")
                fmt.mBinary = true;
            else if (rest == "ascii")
                fmt.mBinary = false;
        }
        if (s.starts_with("arch")) {
            // OpenFOAM's own arch strings are "LSB;label=32;scalar=64" (the
            // "BSB" spelling is what a big-endian host's files would carry).
            // Binary bytes this reader decodes are always little-endian, so a
            // file naming anything else is refused by name rather than
            // silently misread.
            if (s.find("BSB") != std::string_view::npos)
                throw ReadError(
                    "OpenFOAM: big-endian ('BSB') binary files are not supported, only "
                    "little-endian ('LSB')");
            const std::string arch(s);
            std::size_t lp = arch.find("label=");
            if (lp != std::string::npos) {
                int bits = std::atoi(arch.c_str() + lp + 6);
                if (bits)
                    fmt.mLabelBytes = bits / 8;
            }
            std::size_t sp = arch.find("scalar=");
            if (sp != std::string::npos) {
                int bits = std::atoi(arch.c_str() + sp + 7);
                if (bits)
                    fmt.mScalarBytes = bits / 8;
            }
            // Only the widths the binary decoders implement: any other value would size a list
            // by a width it is not decoded with (a zero width divides by zero in the extent check).
            if ((fmt.mLabelBytes != 4 && fmt.mLabelBytes != 8) ||
                (fmt.mScalarBytes != 4 && fmt.mScalarBytes != 8))
                throw ReadError(
                    "OpenFOAM: unsupported 'arch' widths (label and scalar must be 32 "
                    "or 64 bits): " +
                    arch);
        }
        if (s == "}")
            break;
    }
    return fmt;
}

// Strip C-style /* */ and // comments and drop the FoamFile { ... } block.
std::string strip_comments_and_header(std::string_view rText) {
    std::string out;
    out.reserve(rText.size());
    // remove /* */ and //; text between slashes is copied in runs
    for (std::size_t i = 0; i < rText.size();) {
        const std::size_t slash = rText.find('/', i);
        if (slash == std::string_view::npos) {
            out.append(rText.substr(i));
            break;
        }
        out.append(rText.substr(i, slash - i));
        i = slash;
        if (i + 1 < rText.size() && rText[i + 1] == '*') {
            std::size_t e = rText.find("*/", i + 2);
            i = (e == std::string::npos) ? rText.size() : e + 2;
        } else if (i + 1 < rText.size() && rText[i + 1] == '/') {
            std::size_t e = rText.find('\n', i + 2);
            i = (e == std::string::npos) ? rText.size() : e;
        } else {
            out.push_back(rText[i++]);
        }
    }
    // drop FoamFile { ... }
    detail::TextCursor lines(out);
    std::string result;
    result.reserve(out.size() + 1);
    bool in_header = false;
    int depth = 0;
    while (!lines.AtEnd()) {
        const std::string_view line = lines.Line();
        if (line.find("FoamFile") != std::string_view::npos)
            in_header = true;
        if (in_header) {
            for (char c : line) {
                if (c == '{')
                    ++depth;
                else if (c == '}')
                    --depth;
            }
            if (depth <= 0)
                in_header = false;
            continue;
        }
        result.append(line);
        result.push_back('\n');
    }
    return result;
}

// ---- ASCII parsers ----

// The capacity to reserve from a list's count line: the count, capped by what
// the body could hold at `MinBytes` bytes per entry, so a wrong count only
// sizes the reservation and never fails a read that succeeded before.
std::size_t openfoam_count_hint(std::string_view rCount, std::size_t BodySize,
                                std::size_t MinBytes) {
    std::uint64_t n = 0;
    const auto r = std::from_chars(rCount.data(), rCount.data() + rCount.size(), n);
    if (r.ec != std::errc{})
        return 0;
    return static_cast<std::size_t>(std::min<std::uint64_t>(n, BodySize / MinBytes));
}

std::vector<std::array<double, 3>> parse_points_ascii(std::string_view rBody) {
    std::vector<std::array<double, 3>> pts;
    detail::TextCursor lines(rBody);
    std::string t;  // one point's text with its parentheses blanked, reused
    bool in_block = false;
    bool have_n = false;
    while (!lines.AtEnd()) {
        const std::string_view s = foam_trim(lines.Line());
        if (s.empty())
            continue;
        if (!have_n && s.find_first_not_of("0123456789") == std::string_view::npos) {
            have_n = true;
            pts.reserve(openfoam_count_hint(s, rBody.size(), 7));
            continue;
        }
        if (s == "(" && have_n) {
            in_block = true;
            continue;
        }
        if (s == ")" && in_block)
            break;
        if (in_block) {
            // extract up to 3 numbers from within parentheses
            t.assign(s);
            for (char& c : t)
                if (c == '(' || c == ')')
                    c = ' ';
            detail::TextStream ns{std::string_view(t)};
            double a, b, c;
            if (ns >> a >> b >> c)
                pts.push_back({a, b, c});
        }
    }
    return pts;
}

std::vector<Face> parse_faces_ascii(std::string_view rBody) {
    std::vector<Face> faces;
    detail::TextCursor lines(rBody);
    bool in_block = false, have_n = false;
    while (!lines.AtEnd()) {
        const std::string_view s = foam_trim(lines.Line());
        if (s.empty())
            continue;
        if (!have_n && s.find_first_not_of("0123456789") == std::string_view::npos) {
            have_n = true;
            faces.reserve(openfoam_count_hint(s, rBody.size(), 8));
            continue;
        }
        if (s == "(" && have_n) {
            in_block = true;
            continue;
        }
        if (s == ")" && in_block)
            break;
        if (in_block) {
            // form: <count>(<ids...>)
            std::size_t lp = s.find('(');
            std::size_t rp = s.find(')', lp);
            if (lp == std::string_view::npos || rp == std::string_view::npos)
                continue;
            detail::TextStream ns(s.substr(lp + 1, rp - lp - 1));
            Face f;
            std::int64_t v;
            while (ns >> v)
                f.push_back(v);
            faces.push_back(std::move(f));
        }
    }
    return faces;
}

std::vector<std::int64_t> parse_int_list_ascii(std::string_view rBody) {
    std::vector<std::int64_t> out;
    detail::TextCursor lines(rBody);
    bool in_block = false, have_n = false;
    while (!lines.AtEnd()) {
        const std::string_view s = foam_trim(lines.Line());
        if (s.empty())
            continue;
        if (!have_n && s.find_first_not_of("0123456789") == std::string_view::npos) {
            have_n = true;
            out.reserve(openfoam_count_hint(s, rBody.size(), 2));
            continue;
        }
        if (s == "(") {
            in_block = true;
            continue;
        }
        if (s == ")")
            break;
        if (in_block) {
            detail::TextStream ns(s);
            std::int64_t v;
            while (ns >> v)
                out.push_back(v);
        }
    }
    return out;
}

// Boundary patch descriptor. `mNFaces`/`mStartFace` deliberately mirror
// OpenFOAM's own on-disk `boundary` field names (`nFaces`/`startFace`).
struct Patch {
    std::string mName;
    std::string mType;  ///< the `type` entry; empty when the file omitted it
    std::int64_t mNFaces = 0;
    std::int64_t mStartFace = 0;
};

/**
 * @brief Read the value of key @p pKey from a `boundary` sub-dictionary body.
 *
 * `nFaces`/`startFace` are read with `atoll`, but `type` is a word, so it needs
 * real tokenising. The word-boundary guard matters: a bare `find("type")` also
 * matches `physicalType` and `patchType`, both of which are legal entries in the
 * same dictionary and neither of which is the patch's type.
 *
 * @return the value token, or "" when the key is absent.
 */
std::string openfoam_dict_word(const std::string& rBlock, const char* pKey) {
    const std::size_t klen = std::strlen(pKey);
    std::size_t p = 0;
    while ((p = rBlock.find(pKey, p)) != std::string::npos) {
        const bool left_ok = p == 0 || std::isspace(static_cast<unsigned char>(rBlock[p - 1])) ||
                             rBlock[p - 1] == ';';
        const std::size_t after = p + klen;
        const bool right_ok =
            after < rBlock.size() && std::isspace(static_cast<unsigned char>(rBlock[after]));
        if (!left_ok || !right_ok) {
            p = after;
            continue;
        }
        std::size_t a = after;
        while (a < rBlock.size() && std::isspace(static_cast<unsigned char>(rBlock[a])))
            ++a;
        std::size_t b = a;
        while (b < rBlock.size() && !std::isspace(static_cast<unsigned char>(rBlock[b])) &&
               rBlock[b] != ';')
            ++b;
        return rBlock.substr(a, b - a);
    }
    return "";
}

/**
 * @brief Scan top-level `name { ... }` blocks, matching braces by DEPTH.
 *
 * Shared by `parse_boundary` and the zone-file parsers below: both formats
 * are a flat `N ( name { ... } name { ... } ... )` list of named
 * sub-dictionaries. Depth matching (not "the first `}`") matters here too --
 * a patch's `transform`/`sample` sub-block would otherwise truncate it.
 *
 * @return `(name, block body)` pairs, in file order.
 */
std::vector<std::pair<std::string, std::string>> foam_named_blocks(const std::string& rBody) {
    std::vector<std::pair<std::string, std::string>> blocks;
    std::size_t i = 0, n = rBody.size();
    auto skip_ws = [&](std::size_t& p) {
        while (p < n && std::isspace(static_cast<unsigned char>(rBody[p])))
            ++p;
    };
    while (i < n) {
        skip_ws(i);
        std::size_t start = i;
        while (i < n && !std::isspace(static_cast<unsigned char>(rBody[i])) && rBody[i] != '{' &&
               rBody[i] != '(' && rBody[i] != ')')
            ++i;
        std::string name = rBody.substr(start, i - start);
        skip_ws(i);
        if (i < n && rBody[i] == '{') {
            std::size_t close = std::string::npos;
            int depth = 0;
            for (std::size_t p = i; p < n; ++p) {
                if (rBody[p] == '{') {
                    ++depth;
                } else if (rBody[p] == '}') {
                    if (--depth == 0) {
                        close = p;
                        break;
                    }
                }
            }
            if (close == std::string::npos)
                break;
            if (!name.empty())
                blocks.emplace_back(name, rBody.substr(i + 1, close - i - 1));
            i = close + 1;
        } else if (i < n && (rBody[i] == '(' || rBody[i] == ')')) {
            ++i;  // skip list delimiters
        } else if (name.empty()) {
            ++i;
        }
    }
    return blocks;
}

std::vector<Patch> parse_boundary(const std::string& rBody) {
    std::vector<Patch> patches;
    for (const auto& [name, block] : foam_named_blocks(rBody)) {
        Patch pt;
        pt.mName = name;
        pt.mType = openfoam_dict_word(block, "type");
        bool has_n = false, has_s = false;
        std::size_t np = block.find("nFaces");
        if (np != std::string::npos) {
            pt.mNFaces = std::atoll(block.c_str() + np + 6);
            has_n = true;
        }
        std::size_t sp = block.find("startFace");
        if (sp != std::string::npos) {
            pt.mStartFace = std::atoll(block.c_str() + sp + 9);
            has_s = true;
        }
        if (has_n && has_s)
            patches.push_back(pt);
    }
    return patches;
}

/// One named `cellZone`/`faceZone`/`pointZone` entry: a name plus its member
/// ids (cell/face/point ids, in the file's own numbering).
struct Zone {
    std::string mName;
    std::vector<std::int64_t> mIds;
};

/**
 * @brief Read the `List<label>` value of key @p pKey out of a zone block.
 *
 * Zone blocks look like `type cellZone; cellLabels List<label> 3(0 5 9);` --
 * `flipMap` (a `faceZone`-only `List<bool>`) is deliberately never read: a
 * flip only matters for a zone consumer that walks faces directionally
 * (cyclic AMI construction, e.g.), and `Region`'s `Side` entries carry no
 * orientation bit to hold it in. See doc/formats/openfoam.md.
 */
std::vector<std::int64_t> foam_zone_label_list(const std::string& rBlock, const char* pKey) {
    std::size_t p = rBlock.find(pKey);
    if (p == std::string::npos)
        return {};
    std::size_t lp = rBlock.find('(', p);
    if (lp == std::string::npos)
        return {};
    std::size_t rp = lp + 1;
    int depth = 1;
    while (rp < rBlock.size() && depth > 0) {
        if (rBlock[rp] == '(')
            ++depth;
        else if (rBlock[rp] == ')')
            --depth;
        ++rp;
    }
    const std::string inside = rBlock.substr(lp + 1, rp - lp - 2);
    detail::TextStream ss(inside);
    std::vector<std::int64_t> out;
    std::int64_t v;
    while (ss >> v)
        out.push_back(v);
    return out;
}

/// Parse an **ASCII** `cellZones`/`faceZones`/`pointZones` file body (the
/// comment-and-header-stripped text). `parse_zone_file_binary`, further
/// down, is the binary counterpart -- it cannot reuse this one, since
/// `strip_comments_and_header`'s comment-removal pass is unsafe to run over
/// a binary body.
std::vector<Zone> parse_zone_file(const std::string& rBody, const char* pLabelKey) {
    std::vector<Zone> zones;
    for (const auto& [name, block] : foam_named_blocks(rBody))
        zones.push_back({name, foam_zone_label_list(block, pLabelKey)});
    return zones;
}

// ---- binary parsers ----

// Return (N, offset just after the outer '(').
std::pair<std::int64_t, std::size_t> data_start(std::string_view rRaw) {
    std::size_t end = rRaw.find('}');
    if (end == std::string::npos)
        throw ReadError("OpenFOAM: no FoamFile header");
    std::size_t lp = rRaw.find('(', end);
    if (lp == std::string::npos)
        throw ReadError("OpenFOAM: no data list '('");
    // last integer between end and lp
    std::int64_t n = 0;
    bool found = false;
    std::size_t i = end;
    while (i < lp) {
        if (std::isdigit(static_cast<unsigned char>(rRaw[i]))) {
            std::int64_t v = 0;
            while (i < lp && std::isdigit(static_cast<unsigned char>(rRaw[i])))
                v = v * 10 + (rRaw[i++] - '0');
            n = v;
            found = true;
        } else {
            ++i;
        }
    }
    if (!found)
        throw ReadError("OpenFOAM: no element count before '('");
    return {n, lp + 1};
}

template <typename T>
T read_le(const char* pP) {
    T v;
    std::memcpy(&v, pP, sizeof(T));
    return v;
}

/// A binary list of @p Count entries of @p EntryBytes each, starting at
/// @p Start, must lie inside the file: the count comes from the file's own
/// text, so a truncated or edited file would otherwise be read past its end.
void foam_check_binary_extent(std::string_view rRaw, std::size_t Start, std::int64_t Count,
                              std::size_t EntryBytes, const char* pWhat) {
    if (Count < 0 || Start > rRaw.size() ||
        static_cast<std::uint64_t>(Count) > (rRaw.size() - Start) / EntryBytes)
        throw ReadError(detail::format_compat(
            "OpenFOAM: binary {} list of {} entries runs past the end of the file", pWhat, Count));
}

std::vector<std::array<double, 3>> read_binary_points(std::string_view rRaw, int scalar_bytes) {
    auto [n, start] = data_start(rRaw);
    foam_check_binary_extent(rRaw, start, n, 3 * static_cast<std::size_t>(scalar_bytes), "points");
    std::vector<std::array<double, 3>> pts(static_cast<std::size_t>(n));
    const char* base = rRaw.data() + start;
    for (std::int64_t i = 0; i < n; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::size_t off =
                (static_cast<std::size_t>(i) * 3 + j) * static_cast<std::size_t>(scalar_bytes);
            pts[i][j] = scalar_bytes == 4 ? static_cast<double>(read_le<float>(base + off))
                                          : read_le<double>(base + off);
        }
    }
    return pts;
}

std::vector<std::int64_t> read_binary_labels(std::string_view rRaw, int label_bytes) {
    auto [n, start] = data_start(rRaw);
    foam_check_binary_extent(rRaw, start, n, static_cast<std::size_t>(label_bytes), "label");
    std::vector<std::int64_t> out(static_cast<std::size_t>(n));
    const char* base = rRaw.data() + start;
    for (std::int64_t i = 0; i < n; ++i) {
        std::size_t off = static_cast<std::size_t>(i) * static_cast<std::size_t>(label_bytes);
        out[i] = label_bytes == 4 ? static_cast<std::int64_t>(read_le<std::int32_t>(base + off))
                                  : read_le<std::int64_t>(base + off);
    }
    return out;
}

std::vector<Face> read_binary_faces(std::string_view rRaw, int label_bytes) {
    auto [nfaces, pos] = data_start(rRaw);
    foam_check_binary_extent(rRaw, pos, nfaces, 1, "faces");
    std::vector<Face> faces(static_cast<std::size_t>(nfaces));
    std::size_t p = pos;
    for (std::int64_t i = 0; i < nfaces; ++i) {
        std::size_t lp = rRaw.find('(', p);
        if (lp == std::string::npos)
            throw ReadError("OpenFOAM: missing '(' in faces");
        std::int64_t count = std::atoll(std::string(rRaw.substr(p, lp - p)).c_str());
        std::size_t blob = lp + 1;
        foam_check_binary_extent(rRaw, blob, count, static_cast<std::size_t>(label_bytes),
                                 "face");
        Face f(static_cast<std::size_t>(count));
        for (std::int64_t j = 0; j < count; ++j) {
            std::size_t off =
                blob + static_cast<std::size_t>(j) * static_cast<std::size_t>(label_bytes);
            f[j] = label_bytes == 4
                       ? static_cast<std::int64_t>(read_le<std::int32_t>(rRaw.data() + off))
                       : read_le<std::int64_t>(rRaw.data() + off);
        }
        faces[i] = std::move(f);
        p = blob + static_cast<std::size_t>(count) * static_cast<std::size_t>(label_bytes) + 1;
    }
    return faces;
}

/**
 * @brief Binary counterpart of `parse_zone_file`.
 *
 * A zone file's *structure* (zone count, names, `{`/`}`, the `type ...;`
 * line, the `List<label>` keyword and its own decimal count) stays plain
 * text even under `format binary;` -- only each zone's id payload is raw
 * bytes. That payload can legitimately contain a byte equal to `{`, `}` or
 * `/` (a small id's low byte routinely does, e.g. `123` as little-endian
 * `int32` starts with `0x7B` == `{`), so this walks the text directly on
 * the **raw, unstripped** file bytes -- never through `strip_comments_and_header`,
 * whose comment-removal pass scans the whole body and would desync on
 * exactly those bytes -- and explicitly skips `count * LabelBytes` bytes as
 * one opaque unit wherever it recognizes @p pLabelKey, instead of ever
 * scanning byte-by-byte across a blob for a delimiter.
 */
std::vector<Zone> parse_zone_file_binary(std::string_view rRaw, const char* pLabelKey,
                                         int LabelBytes) {
    std::vector<Zone> zones;
    auto [nzones, pos] = data_start(rRaw);
    std::size_t p = pos;
    for (std::int64_t z = 0; z < nzones; ++z) {
        while (p < rRaw.size() && std::isspace(static_cast<unsigned char>(rRaw[p])))
            ++p;
        const std::size_t name_start = p;
        while (p < rRaw.size() && !std::isspace(static_cast<unsigned char>(rRaw[p])))
            ++p;
        const std::string name(rRaw.substr(name_start, p - name_start));

        const std::size_t brace = rRaw.find('{', p);
        if (brace == std::string_view::npos)
            throw ReadError("OpenFOAM: zone '" + name + "' has no opening '{'");
        const std::size_t key_pos = rRaw.find(pLabelKey, brace + 1);
        if (key_pos == std::string_view::npos)
            throw ReadError("OpenFOAM: zone '" + name + "' has no '" + pLabelKey + "'");
        const std::size_t lparen = rRaw.find('(', key_pos);
        if (lparen == std::string_view::npos)
            throw ReadError("OpenFOAM: zone '" + name + "' has no '(' after '" + pLabelKey + "'");

        std::int64_t count = 0;
        bool found = false;
        for (std::size_t i = key_pos; i < lparen; ++i)
            if (std::isdigit(static_cast<unsigned char>(rRaw[i]))) {
                std::int64_t v = 0;
                while (i < lparen && std::isdigit(static_cast<unsigned char>(rRaw[i])))
                    v = v * 10 + (rRaw[i++] - '0');
                count = v;
                found = true;
            }
        if (!found)
            throw ReadError("OpenFOAM: zone '" + name + "' has no count before '('");

        foam_check_binary_extent(rRaw, lparen + 1, count, static_cast<std::size_t>(LabelBytes),
                                 "zone");
        std::vector<std::int64_t> ids(static_cast<std::size_t>(count));
        const char* base = rRaw.data() + lparen + 1;
        for (std::int64_t i = 0; i < count; ++i) {
            const std::size_t off =
                static_cast<std::size_t>(i) * static_cast<std::size_t>(LabelBytes);
            ids[static_cast<std::size_t>(i)] =
                LabelBytes == 4 ? static_cast<std::int64_t>(read_le<std::int32_t>(base + off))
                                : read_le<std::int64_t>(base + off);
        }
        zones.push_back({name, std::move(ids)});

        // Resume the text scan only *after* the raw blob -- everything from
        // here to this zone's own closing '}' (its `);` terminator, the
        // closing brace) is plain text again.
        const std::size_t blob_end =
            lparen + 1 + static_cast<std::size_t>(count) * static_cast<std::size_t>(LabelBytes);
        const std::size_t zone_close = rRaw.find('}', blob_end);
        if (zone_close == std::string_view::npos)
            throw ReadError("OpenFOAM: zone '" + name + "' has no closing '}'");
        p = zone_close + 1;
    }
    return zones;
}

// ---- dispatch readers ----

std::vector<std::array<double, 3>> read_points(const fs::path& rPath) {
    const detail::FileSource source = read_whole(rPath.string());
    const std::string_view raw = source.View();
    const FoamFormat fmt = detect_format(raw);
    if (fmt.mBinary)
        return read_binary_points(raw, fmt.mScalarBytes);
    return parse_points_ascii(strip_comments_and_header(raw));
}

std::vector<Face> read_faces(const fs::path& rPath) {
    const detail::FileSource source = read_whole(rPath.string());
    const std::string_view raw = source.View();
    const FoamFormat fmt = detect_format(raw);
    if (fmt.mBinary)
        return read_binary_faces(raw, fmt.mLabelBytes);
    return parse_faces_ascii(strip_comments_and_header(raw));
}

std::vector<std::int64_t> read_int_list(const fs::path& rPath) {
    const detail::FileSource source = read_whole(rPath.string());
    const std::string_view raw = source.View();
    const FoamFormat fmt = detect_format(raw);
    if (fmt.mBinary)
        return read_binary_labels(raw, fmt.mLabelBytes);
    return parse_int_list_ascii(strip_comments_and_header(raw));
}

// ---- geometry ----

// The face-to-cell kernel (triple, match_top, build_*, reconstruct_cell) is
// shared with the Fluent reader: formats/face_cells_common.hpp.
using face_cells::P3;
using face_cells::reconstruct_cell;
using face_cells::unique_node_count;

// ---- decomposed (processorN) cases (v11.4.0, roadmap §1 tier B2) ----

/// A single `polyMesh`'s raw, un-reconstructed file contents -- what
/// `read_openfoam` used to read inline before it could also come from
/// `reconstruct_decomposed`.
struct RawPolyMesh {
    P3 mPoints;
    std::vector<Face> mFaces;
    std::vector<std::int64_t> mOwner;
    std::vector<std::int64_t> mNeighbour;  ///< only the internal faces
    std::vector<Patch> mBoundary;
    /// Decomposed cases only: per patch, per face, the processor (index into the processor list)
    /// and the face's offset within that processor's copy of the patch -- where a field's patch
    /// values live.
    std::vector<std::vector<std::pair<std::size_t, std::int64_t>>> mBoundarySource;
    /// Decomposed cases only: the decomposed case's global face id (the `faceProcAddressing`
    /// number, minus one) -> this reconstruction's face index; -1 for a face not reconstructed
    /// (a `processor` patch face).
    std::vector<std::int64_t> mFaceOfGlobal;
};

/// Refuse a polyMesh whose lists disagree, before anything indexes one list
/// with another's values: a face naming a point that is not there, an owner or
/// neighbour list longer than the faces, a negative owner, a cell id the lists
/// cannot account for, a patch reaching past the last face. A short points
/// list (one line skipped for having fewer than three numbers, or a list the
/// parser read nothing from) used to leave the cell reconstruction reading past
/// the points.
void foam_check_polymesh(const RawPolyMesh& rRaw) {
    const std::size_t n_points = rRaw.mPoints.size();
    const std::size_t n_faces = rRaw.mFaces.size();
    for (std::size_t f = 0; f < n_faces; ++f)
        for (const std::int64_t id : rRaw.mFaces[f])
            if (id < 0 || static_cast<std::uint64_t>(id) >= n_points)
                throw ReadError(detail::format_compat(
                    "OpenFOAM: face {} names point {}, but the mesh has {} points", f, id,
                    n_points));
    if (rRaw.mOwner.size() > n_faces)
        throw ReadError(detail::format_compat(
            "OpenFOAM: owner lists {} faces, but the mesh has {}", rRaw.mOwner.size(), n_faces));
    if (rRaw.mNeighbour.size() > n_faces)
        throw ReadError(detail::format_compat("OpenFOAM: neighbour lists {} faces, but the mesh has {}",
                                              rRaw.mNeighbour.size(), n_faces));
    const std::uint64_t n_ids = rRaw.mOwner.size() + rRaw.mNeighbour.size();
    for (std::size_t f = 0; f < rRaw.mOwner.size(); ++f) {
        const std::int64_t id = rRaw.mOwner[f];
        if (id < 0 || static_cast<std::uint64_t>(id) >= n_ids)
            throw ReadError(detail::format_compat("OpenFOAM: face {} has owner {}", f, id));
    }
    for (std::size_t f = 0; f < rRaw.mNeighbour.size(); ++f) {
        const std::int64_t id = rRaw.mNeighbour[f];
        if (id < -1 || (id >= 0 && static_cast<std::uint64_t>(id) >= n_ids))
            throw ReadError(detail::format_compat("OpenFOAM: face {} has neighbour {}", f, id));
    }
    for (const Patch& patch : rRaw.mBoundary)
        if (patch.mStartFace < 0 || patch.mNFaces < 0 ||
            static_cast<std::uint64_t>(patch.mStartFace) + static_cast<std::uint64_t>(patch.mNFaces) >
                n_faces)
            throw ReadError(detail::format_compat(
                "OpenFOAM: patch '{}' covers faces {} to {}, but the mesh has {} faces",
                patch.mName, patch.mStartFace, patch.mStartFace + patch.mNFaces, n_faces));
}

RawPolyMesh read_raw_polymesh(const fs::path& rPoly) {
    RawPolyMesh raw;
    raw.mPoints = read_points(rPoly / "points");
    raw.mFaces = read_faces(rPoly / "faces");
    raw.mOwner = read_int_list(rPoly / "owner");
    if (fs::exists(rPoly / "neighbour"))
        raw.mNeighbour = read_int_list(rPoly / "neighbour");
    if (fs::exists(rPoly / "boundary"))
        raw.mBoundary = parse_boundary(
            strip_comments_and_header(read_whole((rPoly / "boundary").string()).View()));
    foam_check_polymesh(raw);
    return raw;
}

/// One processor's `polyMesh` plus the addressing lists that map its local
/// ids back onto the undecomposed case's global ones.
struct ProcMesh {
    RawPolyMesh mRaw;
    std::vector<std::int64_t> mPointAddr;     ///< local point -> global point id
    std::vector<std::int64_t> mCellAddr;      ///< local cell -> global cell id
    std::vector<std::int64_t> mFaceAddr;      ///< local face -> signed (global face id + 1)
    std::vector<std::int64_t> mBoundaryAddr;  ///< local patch -> global patch id, -1 if none
};

/// The processor-local face (and its orientation) claiming a given global
/// face id -- one entry for a face interior to one processor or an original
/// external boundary face, two for a face split by decomposition.
struct DecompFaceClaim {
    std::int64_t mProc = -1;
    std::int64_t mLocalFace = -1;
    bool mFlipped = false;
};

/// The local patch containing local face @p LocalFace, or `npos`.
std::size_t foam_patch_of_local_face(const std::vector<Patch>& rBoundary, std::int64_t LocalFace) {
    for (std::size_t p = 0; p < rBoundary.size(); ++p)
        if (LocalFace >= rBoundary[p].mStartFace &&
            LocalFace < rBoundary[p].mStartFace + rBoundary[p].mNFaces)
            return p;
    return static_cast<std::size_t>(-1);
}

/// The sorted processor indices of a decomposed case's `<root>/processorN/`
/// directories that carry a `constant/polyMesh` -- not assumed contiguous.
std::vector<std::size_t> foam_processor_ids(const fs::path& rCaseRoot) {
    std::vector<std::size_t> ids;
    std::error_code ec;
    for (const auto& entry : fs::directory_iterator(rCaseRoot, ec)) {
        const std::string name = entry.path().filename().string();
        if (name.rfind("processor", 0) != 0)
            continue;
        const std::string digits = name.substr(std::strlen("processor"));
        if (digits.empty() || digits.find_first_not_of("0123456789") != std::string::npos)
            continue;
        if (fs::exists(entry.path() / "constant" / "polyMesh"))
            ids.push_back(static_cast<std::size_t>(std::atoll(digits.c_str())));
    }
    std::sort(ids.begin(), ids.end());
    return ids;
}

/**
 * @brief Reassemble a decomposed case's `processorN` directories into one
 * global `RawPolyMesh`, mirroring what `reconstructParMesh` does on disk.
 *
 * Points/cells/faces are placed at the global ids their processor's own
 * `*ProcAddressing` files (plain `labelList`s, read like `owner`) name.
 * A global face claimed by exactly one processor is either an original
 * external boundary face or one interior to that processor alone; claimed by
 * two, it is an internal face `decomposePar` split at a processor boundary --
 * the positive-signed `faceProcAddressing` entry names the true owner side,
 * the negative-signed one the neighbour side (the sign meaning "this local
 * copy is stored reversed relative to the global orientation"). A boundary
 * face's global patch comes from its owning processor's `boundaryProcAddressing`
 * (missing/negative marks a `processor*` inter-rank patch, dropped -- it has no
 * counterpart in the original case).
 */
RawPolyMesh reconstruct_decomposed(const fs::path& rCaseRoot,
                                   const std::vector<std::size_t>& rProcIds) {
    std::vector<ProcMesh> procs(rProcIds.size());
    std::int64_t max_point = -1, max_cell = -1, max_face = -1;
    for (std::size_t k = 0; k < rProcIds.size(); ++k) {
        const fs::path poly =
            rCaseRoot / ("processor" + std::to_string(rProcIds[k])) / "constant" / "polyMesh";
        ProcMesh& pm = procs[k];
        pm.mRaw = read_raw_polymesh(poly);
        pm.mPointAddr = read_int_list(poly / "pointProcAddressing");
        pm.mCellAddr = read_int_list(poly / "cellProcAddressing");
        pm.mFaceAddr = read_int_list(poly / "faceProcAddressing");
        if (fs::exists(poly / "boundaryProcAddressing"))
            pm.mBoundaryAddr = read_int_list(poly / "boundaryProcAddressing");
        for (std::int64_t v : pm.mPointAddr)
            max_point = std::max(max_point, v);
        for (std::int64_t v : pm.mCellAddr)
            max_cell = std::max(max_cell, v);
        for (std::int64_t v : pm.mFaceAddr)
            max_face = std::max(max_face, std::abs(v) - 1);
    }
    const std::size_t n_points = static_cast<std::size_t>(max_point + 1);
    const std::size_t n_faces = static_cast<std::size_t>(max_face + 1);

    P3 points(n_points);
    for (const ProcMesh& pm : procs)
        for (std::size_t i = 0; i < pm.mPointAddr.size(); ++i)
            points[static_cast<std::size_t>(pm.mPointAddr[i])] = pm.mRaw.mPoints[i];

    // Every processor-local face that claims a given global id.
    std::vector<std::array<DecompFaceClaim, 2>> claims(n_faces);
    std::vector<std::uint8_t> n_claims(n_faces, 0);
    for (std::size_t k = 0; k < procs.size(); ++k) {
        const auto& fa = procs[k].mFaceAddr;
        for (std::size_t i = 0; i < fa.size(); ++i) {
            const std::size_t g = static_cast<std::size_t>(std::abs(fa[i]) - 1);
            const std::uint8_t slot = n_claims[g]++;
            if (slot < 2)
                claims[g][slot] = {static_cast<std::int64_t>(k), static_cast<std::int64_t>(i),
                                   fa[i] < 0};
        }
    }

    std::vector<Face> internal_faces;
    std::vector<std::int64_t> internal_owner, internal_neighbour, internal_gids;
    std::map<std::int64_t, std::vector<std::int64_t>> patch_gids;
    // Global patch id -> its (name, type) plus the member faces' node rings
    // and owner cells, in the order they are found.
    std::map<std::int64_t, Patch> patch_table;
    std::map<std::int64_t, std::vector<Face>> patch_faces;
    std::map<std::int64_t, std::vector<std::int64_t>> patch_owners;
    std::map<std::int64_t, std::vector<std::pair<std::size_t, std::int64_t>>> patch_sources;
    std::size_t n_dropped_processor_faces = 0, n_dropped_unclaimed = 0;

    for (std::size_t g = 0; g < n_faces; ++g) {
        if (n_claims[g] == 0) {
            ++n_dropped_unclaimed;  // an id `faceProcAddressing` never actually used
            continue;
        }
        const DecompFaceClaim* owner_claim = nullptr;
        const DecompFaceClaim* neigh_claim = nullptr;
        for (std::uint8_t k = 0; k < std::min<std::uint8_t>(n_claims[g], 2); ++k) {
            const DecompFaceClaim& c = claims[g][k];
            (c.mFlipped ? neigh_claim : owner_claim) = &c;
        }
        if (!owner_claim)
            owner_claim = &claims[g][0];  // defensive: both flipped should not happen

        const ProcMesh& op = procs[static_cast<std::size_t>(owner_claim->mProc)];
        const Face& lf = op.mRaw.mFaces[static_cast<std::size_t>(owner_claim->mLocalFace)];
        Face gf(lf.size());
        for (std::size_t k = 0; k < lf.size(); ++k)
            gf[k] = op.mPointAddr[static_cast<std::size_t>(lf[k])];
        const std::int64_t owner_cell = op.mCellAddr[static_cast<std::size_t>(
            op.mRaw.mOwner[static_cast<std::size_t>(owner_claim->mLocalFace)])];

        if (neigh_claim) {
            const ProcMesh& np = procs[static_cast<std::size_t>(neigh_claim->mProc)];
            const std::int64_t neigh_cell = np.mCellAddr[static_cast<std::size_t>(
                np.mRaw.mOwner[static_cast<std::size_t>(neigh_claim->mLocalFace)])];
            internal_faces.push_back(std::move(gf));
            internal_owner.push_back(owner_cell);
            internal_neighbour.push_back(neigh_cell);
            internal_gids.push_back(static_cast<std::int64_t>(g));
            continue;
        }

        // A face in the owning processor's internal range with no second claim lies inside
        // that processor: both cells are its own.
        if (static_cast<std::size_t>(owner_claim->mLocalFace) < op.mRaw.mNeighbour.size()) {
            const std::int64_t local_neigh =
                op.mRaw.mNeighbour[static_cast<std::size_t>(owner_claim->mLocalFace)];
            if (local_neigh >= 0 && static_cast<std::size_t>(local_neigh) < op.mCellAddr.size()) {
                internal_faces.push_back(std::move(gf));
                internal_owner.push_back(owner_cell);
                internal_neighbour.push_back(op.mCellAddr[static_cast<std::size_t>(local_neigh)]);
                internal_gids.push_back(static_cast<std::int64_t>(g));
                continue;
            }
        }

        // A genuine boundary face: resolve its global patch via the owning
        // processor's own local patch + `boundaryProcAddressing`.
        const std::size_t local_patch =
            foam_patch_of_local_face(op.mRaw.mBoundary, owner_claim->mLocalFace);
        std::int64_t global_patch = -1;
        if (local_patch != static_cast<std::size_t>(-1) && local_patch < op.mBoundaryAddr.size())
            global_patch = op.mBoundaryAddr[local_patch];
        const bool is_processor_patch =
            global_patch < 0 || (local_patch != static_cast<std::size_t>(-1) &&
                                 op.mRaw.mBoundary[local_patch].mType.rfind("processor", 0) == 0);
        if (is_processor_patch) {
            ++n_dropped_processor_faces;
            continue;
        }
        if (!patch_table.count(global_patch)) {
            Patch p = op.mRaw.mBoundary[local_patch];  // name/type only; counts recomputed below
            p.mNFaces = 0;
            p.mStartFace = 0;
            patch_table[global_patch] = p;
        }
        patch_faces[global_patch].push_back(std::move(gf));
        patch_owners[global_patch].push_back(owner_cell);
        patch_gids[global_patch].push_back(static_cast<std::int64_t>(g));
        patch_sources[global_patch].emplace_back(
            static_cast<std::size_t>(owner_claim->mProc),
            owner_claim->mLocalFace - op.mRaw.mBoundary[local_patch].mStartFace);
    }
    if (n_dropped_processor_faces > 0)
        log::info("OpenFOAM: reconstructed case drops {} inter-processor patch face(s)",
                  n_dropped_processor_faces);
    if (n_dropped_unclaimed > 0)
        log::warn("OpenFOAM: {} face id(s) in *ProcAddressing were never claimed",
                  n_dropped_unclaimed);

    RawPolyMesh out;
    out.mPoints = std::move(points);
    out.mFaces = std::move(internal_faces);
    out.mOwner = internal_owner;
    out.mNeighbour = std::move(internal_neighbour);
    out.mFaceOfGlobal.assign(n_faces, -1);
    for (std::size_t i = 0; i < internal_gids.size(); ++i)
        out.mFaceOfGlobal[static_cast<std::size_t>(internal_gids[i])] =
            static_cast<std::int64_t>(i);
    for (auto& kv : patch_table) {
        kv.second.mStartFace = static_cast<std::int64_t>(out.mFaces.size());
        for (std::size_t j = 0; j < patch_gids[kv.first].size(); ++j)
            out.mFaceOfGlobal[static_cast<std::size_t>(patch_gids[kv.first][j])] =
                kv.second.mStartFace + static_cast<std::int64_t>(j);
        for (Face& f : patch_faces[kv.first])
            out.mFaces.push_back(std::move(f));
        for (std::int64_t c : patch_owners[kv.first])
            out.mOwner.push_back(c);
        kv.second.mNFaces = static_cast<std::int64_t>(out.mFaces.size()) - kv.second.mStartFace;
        out.mBoundary.push_back(kv.second);
        out.mBoundarySource.push_back(std::move(patch_sources[kv.first]));
    }
    return out;
}

/// One `cellZones`/`faceZones`/`pointZones` file; empty when the file is absent.
std::vector<Zone> foam_read_zone_file(const fs::path& rPath, const char* pKey) {
    if (!fs::exists(rPath))
        return {};
    const detail::FileSource source = read_whole(rPath.string());
    const FoamFormat fmt = detect_format(source.View());
    // strip_comments_and_header's comment-removal pass scans the whole body and is unsafe over
    // raw id bytes -- see parse_zone_file_binary's own doc comment.
    if (fmt.mBinary)
        return parse_zone_file_binary(source.View(), pKey, fmt.mLabelBytes);
    return parse_zone_file(strip_comments_and_header(source.View()), pKey);
}

/// A decomposed case's zones of one kind, back in the global numbering: each processor's
/// `<file>` ids go through that processor's `<kind>ProcAddressing` (a face's, through the signed
/// `faceProcAddressing` and then @p pFaceOfGlobal), and a zone's members from every processor are
/// merged, sorted and made unique -- a face on a processor boundary is listed by both. Zones keep
/// the order they first appear in. An id no addressing entry covers is dropped, as is a face that
/// the reconstruction did not keep (a `processor` patch face).
std::vector<Zone> foam_read_decomposed_zones(const fs::path& rCaseRoot,
                                             const std::vector<std::size_t>& rProcIds,
                                             const char* pFile, const char* pKey,
                                             const char* pAddressing,
                                             const std::vector<std::int64_t>* pFaceOfGlobal) {
    std::vector<std::string> order;
    std::map<std::string, std::set<std::int64_t>> merged;
    for (const std::size_t proc : rProcIds) {
        const fs::path poly =
            rCaseRoot / ("processor" + std::to_string(proc)) / "constant" / "polyMesh";
        const std::vector<Zone> zones = foam_read_zone_file(poly / pFile, pKey);
        if (zones.empty())
            continue;
        const std::vector<std::int64_t> addr = read_int_list(poly / pAddressing);
        for (const Zone& z : zones) {
            if (!merged.count(z.mName))
                order.push_back(z.mName);
            std::set<std::int64_t>& out = merged[z.mName];
            for (const std::int64_t local : z.mIds) {
                if (local < 0 || static_cast<std::size_t>(local) >= addr.size())
                    continue;
                std::int64_t id = addr[static_cast<std::size_t>(local)];
                if (pFaceOfGlobal) {
                    const std::int64_t g = std::abs(id) - 1;
                    if (g < 0 || static_cast<std::size_t>(g) >= pFaceOfGlobal->size())
                        continue;
                    id = (*pFaceOfGlobal)[static_cast<std::size_t>(g)];
                    if (id < 0)
                        continue;
                }
                out.insert(id);
            }
        }
    }
    std::vector<Zone> result;
    for (const std::string& name : order)
        result.push_back(
            {name, std::vector<std::int64_t>(merged[name].begin(), merged[name].end())});
    return result;
}

// ---- time-directory fields (v11.4.0, roadmap §1 tier B2) ----

/// One field's shape: number of components, and whether it belongs on
/// points (`pointScalarField`/`pointVectorField`) rather than cells.
struct FoamFieldClass {
    int mComponents = 0;
    bool mIsPoint = false;
    bool mSupported = true;
};

FoamFieldClass foam_field_class(const std::string& rClass) {
    if (rClass == "volScalarField")
        return {1, false, true};
    if (rClass == "volVectorField")
        return {3, false, true};
    if (rClass == "volSymmTensorField")
        return {6, false, true};
    if (rClass == "volTensorField")
        return {9, false, true};
    if (rClass == "pointScalarField")
        return {1, true, true};
    if (rClass == "pointVectorField")
        return {3, true, true};
    return {0, false, false};  // e.g. surfaceScalarField: no cell/point home
}

/// A field value list: `mFlat` holds `mComponents` entries (uniform) or
/// `mCount * mComponents` (nonuniform), read by `foam_read_internal_field`.
struct FoamField {
    bool mUniform = false;
    std::int64_t mCount = 0;
    std::vector<double> mFlat;
};

/// Parse `uniform <value>` (`rText` starting AT the `uniform` keyword).
std::vector<double> foam_scan_uniform_value(std::string_view rText, int components) {
    std::size_t p = std::strlen("uniform");
    while (p < rText.size() && std::isspace(static_cast<unsigned char>(rText[p])))
        ++p;
    std::vector<double> out;
    if (components == 1) {
        out.push_back(detail::parse_double(std::string(rText.substr(p))));
        return out;
    }
    const std::size_t lp = rText.find('(', p);
    const std::size_t rp = rText.find(')', lp);
    if (lp == std::string::npos || rp == std::string::npos)
        return out;
    detail::TextStream ss(rText.substr(lp + 1, rp - lp - 1));
    double v;
    while (ss >> v)
        out.push_back(v);
    return out;
}

/// Parse `nonuniform List<T>\n<N>\n(\n<entries>\n)` (`rText` starting AT the
/// `nonuniform` keyword) -- ASCII only; the binary variant is read directly
/// via `data_start` in `foam_read_internal_field`, which needs the whole raw
/// buffer rather than a text view.
FoamField foam_scan_nonuniform_list(std::string_view rText, int components) {
    FoamField out;
    detail::TextCursor lines(rText);
    bool have_n = false;
    std::int64_t n = 0;
    while (!lines.AtEnd()) {
        const std::string_view s = foam_trim(lines.Line());
        if (s.empty())
            continue;
        if (!have_n) {
            if (s.find_first_not_of("0123456789") == std::string_view::npos) {
                n = detail::strtoll_token(s);
                have_n = true;
            }
            continue;
        }
        if (s == "(")
            break;
    }
    out.mCount = n;
    out.mFlat.reserve(std::min<std::size_t>(
        static_cast<std::size_t>(n) * static_cast<std::size_t>(components), rText.size()));
    std::string t;  // one entry's text with its parentheses blanked, reused
    for (std::int64_t i = 0; i < n && !lines.AtEnd();) {
        const std::string_view s = foam_trim(lines.Line());
        if (s.empty())
            continue;
        if (components == 1) {
            out.mFlat.push_back(detail::parse_double_prefix(s));
        } else {
            t.assign(s);
            for (char& c : t)
                if (c == '(' || c == ')')
                    c = ' ';
            detail::TextStream ls{std::string_view(t)};
            double v;
            while (ls >> v)
                out.mFlat.push_back(v);
        }
        ++i;
    }
    return out;
}

/**
 * @brief Read a field file's `internalField`.
 *
 * `uniform` is always plain text, in both ASCII and binary field files (a
 * single small value is never worth binary-encoding), so it is scanned the
 * same way regardless of `FoamFormat`. `nonuniform` follows `points`/`faces`/
 * `owner`'s own dispatch: ASCII is line-scanned, binary reuses `data_start` --
 * safe here because nothing between the FoamFile header and `internalField`'s
 * own data list can introduce a stray `(` (`dimensions` uses `[...]`).
 */
FoamField foam_read_internal_field(const fs::path& rPath, int components) {
    const detail::FileSource source = read_whole(rPath.string());
    const std::string_view raw = source.View();
    const FoamFormat fmt = detect_format(raw);
    const std::size_t kp = raw.find("internalField");
    if (kp == std::string::npos)
        throw ReadError("OpenFOAM: field file has no internalField: " + rPath.string());
    std::size_t p = kp + std::strlen("internalField");
    while (p < raw.size() && std::isspace(static_cast<unsigned char>(raw[p])))
        ++p;
    if (raw.compare(p, 7, "uniform") == 0) {
        FoamField out;
        out.mUniform = true;
        out.mFlat = foam_scan_uniform_value(raw.substr(p), components);
        return out;
    }
    if (raw.compare(p, 10, "nonuniform") != 0)
        throw ReadError("OpenFOAM: internalField is neither uniform nor nonuniform: " +
                        rPath.string());
    if (!fmt.mBinary)
        return foam_scan_nonuniform_list(raw.substr(p), components);

    auto [n, start] = data_start(raw);
    foam_check_binary_extent(raw, start, n,
                             static_cast<std::size_t>(components) *
                                 static_cast<std::size_t>(fmt.mScalarBytes),
                             "field");
    FoamField out;
    out.mCount = n;
    out.mFlat.resize(static_cast<std::size_t>(n) * static_cast<std::size_t>(components));
    const char* base = raw.data() + start;
    for (std::int64_t i = 0; i < n; ++i)
        for (int c = 0; c < components; ++c) {
            const std::size_t off =
                (static_cast<std::size_t>(i) * static_cast<std::size_t>(components) +
                 static_cast<std::size_t>(c)) *
                static_cast<std::size_t>(fmt.mScalarBytes);
            out.mFlat[static_cast<std::size_t>(i) * static_cast<std::size_t>(components) +
                      static_cast<std::size_t>(c)] =
                fmt.mScalarBytes == 4 ? static_cast<double>(read_le<float>(base + off))
                                      : read_le<double>(base + off);
        }
    return out;
}

/// One patch's `value` entry from a field's `boundaryField`: empty (`mFlat`) when the patch has no
/// value of its own (`zeroGradient`, `empty`, `$internalField`, ...).
struct FoamPatchField {
    std::string mName;
    FoamField mValue;
};

/**
 * @brief Read the per-patch `value` entries of a field file's `boundaryField`.
 *
 * A small dictionary scanner rather than the line scan `internalField` gets: patch values are
 * short lists OpenFOAM writes inline (`4(1 2 3 4)`), as `N{value}`, or as packed binary, and any
 * other entry of a patch (`inletValue`, `uniformValue`, nested dictionaries) has to be skipped
 * without mistaking its bytes or braces for structure. It never throws: a file it cannot follow
 * ends the scan and the patches read up to there are returned.
 */
class FoamBoundaryScanner {
public:
    FoamBoundaryScanner(std::string_view Raw, const FoamFormat& rFormat, int Components)
        : mRaw(Raw), mFormat(rFormat), mComponents(Components) {}

    std::vector<FoamPatchField> Scan() {
        std::vector<FoamPatchField> out;
        while (!mFailed) {
            SkipSpace();
            if (mPos >= mRaw.size())
                break;
            const std::string key = ReadToken();
            if (key.empty()) {
                mFailed = true;
                break;
            }
            SkipSpace();
            if (mPos < mRaw.size() && mRaw[mPos] == '{') {
                if (key == "boundaryField")
                    ScanPatches(out);
                else
                    ScanBody(1, nullptr);
            } else {
                SkipValue();
            }
        }
        return out;
    }

private:
    static constexpr int kMaxDepth = 32;

    static bool IsSpace(char c) { return std::isspace(static_cast<unsigned char>(c)) != 0; }
    static bool IsStop(char c) {
        return IsSpace(c) || c == ';' || c == '{' || c == '}' || c == '(' || c == ')';
    }

    void SkipSpace() {
        while (mPos < mRaw.size()) {
            if (IsSpace(mRaw[mPos])) {
                ++mPos;
            } else if (mRaw.compare(mPos, 2, "//") == 0) {
                while (mPos < mRaw.size() && mRaw[mPos] != '\n')
                    ++mPos;
            } else if (mRaw.compare(mPos, 2, "/*") == 0) {
                const std::size_t end = mRaw.find("*/", mPos + 2);
                mPos = end == std::string_view::npos ? mRaw.size() : end + 2;
            } else {
                break;
            }
        }
    }

    /// A word or a quoted string (quotes dropped); empty at punctuation or the end.
    std::string ReadToken() {
        SkipSpace();
        std::string tok;
        if (mPos < mRaw.size() && mRaw[mPos] == '"') {
            ++mPos;
            while (mPos < mRaw.size() && mRaw[mPos] != '"') {
                if (mRaw[mPos] == '\\' && mPos + 1 < mRaw.size())
                    ++mPos;
                tok += mRaw[mPos++];
            }
            if (mPos < mRaw.size())
                ++mPos;
            else
                mFailed = true;
            return tok;
        }
        while (mPos < mRaw.size() && !IsStop(mRaw[mPos]))
            tok += mRaw[mPos++];
        return tok;
    }

    /// Past a `( ... )` group (the current character is the `(`), parentheses balanced.
    void SkipParens() {
        int depth = 0;
        while (mPos < mRaw.size()) {
            const char c = mRaw[mPos];
            if (c == '"') {
                ReadToken();
                continue;
            }
            ++mPos;
            if (c == '(')
                ++depth;
            else if (c == ')' && --depth == 0)
                return;
        }
        mFailed = true;
    }

    /// The entries of a dictionary, from its `{` to its `}`; @p pValue receives a `value` entry.
    void ScanBody(int Depth, FoamField* pValue) {
        if (Depth > kMaxDepth) {
            mFailed = true;
            return;
        }
        ++mPos;  // the '{'
        while (!mFailed) {
            SkipSpace();
            if (mPos >= mRaw.size()) {
                mFailed = true;
                return;
            }
            if (mRaw[mPos] == '}') {
                ++mPos;
                return;
            }
            const std::string key = ReadToken();
            if (key.empty()) {
                mFailed = true;
                return;
            }
            SkipSpace();
            if (mPos < mRaw.size() && mRaw[mPos] == '{')
                ScanBody(Depth + 1, nullptr);
            else if (pValue != nullptr && key == "value")
                ReadValue(*pValue);
            else
                SkipValue();
        }
    }

    void ScanPatches(std::vector<FoamPatchField>& rOut) {
        ++mPos;  // the '{'
        while (!mFailed) {
            SkipSpace();
            if (mPos >= mRaw.size()) {
                mFailed = true;
                return;
            }
            if (mRaw[mPos] == '}') {
                ++mPos;
                return;
            }
            FoamPatchField pf;
            pf.mName = ReadToken();
            if (pf.mName.empty()) {
                mFailed = true;
                return;
            }
            SkipSpace();
            if (mPos < mRaw.size() && mRaw[mPos] == '{') {
                ScanBody(2, &pf.mValue);
                if (!mFailed)
                    rOut.push_back(std::move(pf));
            } else {
                SkipValue();
            }
        }
    }

    /// The rest of an entry, through its `;` (or up to a closing `}`).
    void SkipValue() {
        while (!mFailed) {
            SkipSpace();
            if (mPos >= mRaw.size()) {
                mFailed = true;
                return;
            }
            const char c = mRaw[mPos];
            if (c == ';') {
                ++mPos;
                return;
            }
            if (c == '}')
                return;
            if (c == '{') {
                ScanBody(kMaxDepth, nullptr);  // an entry never holds a bare dictionary
                return;
            }
            if (c == '(') {
                SkipParens();
                continue;
            }
            const std::string tok = ReadToken();
            if (tok.empty()) {
                mFailed = true;
                return;
            }
            if (tok == "nonuniform")
                ReadNonuniform(nullptr);
        }
    }

    /// Entries per list element for a `List<T>` type name, -1 for a label, 0 when unknown.
    static int TypeComponents(const std::string& rType) {
        if (rType == "scalar" || rType == "sphericalTensor")
            return 1;
        if (rType == "vector")
            return 3;
        if (rType == "symmTensor")
            return 6;
        if (rType == "tensor")
            return 9;
        if (rType == "label")
            return -1;
        return 0;
    }

    /// Numbers of an ASCII list body, the opening `(` or `{` already consumed, through its
    /// closing @p Close.
    void ReadAsciiNumbers(char Close, std::vector<double>* pOut) {
        int depth = 0;
        while (mPos < mRaw.size()) {
            SkipSpace();
            if (mPos >= mRaw.size())
                break;
            const char c = mRaw[mPos];
            if (c == Close && depth == 0) {
                ++mPos;
                return;
            }
            if (c == '(') {
                ++depth;
                ++mPos;
            } else if (c == ')') {
                if (--depth < 0)
                    break;
                ++mPos;
            } else if (c == ';' || c == '{' || c == '}') {
                break;
            } else {
                const std::size_t start = mPos;
                while (mPos < mRaw.size() && !IsStop(mRaw[mPos]))
                    ++mPos;
                if (pOut != nullptr)
                    pOut->push_back(detail::parse_double_prefix(mRaw.substr(start, mPos - start)));
            }
        }
        mFailed = true;
    }

    /// After the `nonuniform` keyword: past the list; into @p pOut when it holds the field's own
    /// type.
    void ReadNonuniform(FoamField* pOut) {
        std::string tok = ReadToken();
        int comps = mComponents;  // `nonuniform 0()` names no type
        if (tok.rfind("List<", 0) == 0) {
            const std::size_t gt = tok.find('>');
            comps = TypeComponents(tok.substr(5, gt == std::string::npos ? gt : gt - 5));
            tok = ReadToken();
        }
        if (tok.empty() || tok.find_first_not_of("0123456789") != std::string::npos) {
            mFailed = true;
            return;
        }
        const std::int64_t count = detail::strtoll_token(tok);
        SkipSpace();
        const bool want = pOut != nullptr && comps == mComponents;
        if (mPos < mRaw.size() && mRaw[mPos] == '{') {  // N{value}: one entry repeated
            ++mPos;
            std::vector<double> one;
            ReadAsciiNumbers('}', &one);
            if (!mFailed && want && one.size() == static_cast<std::size_t>(mComponents)) {
                pOut->mUniform = true;
                pOut->mCount = count;
                pOut->mFlat = std::move(one);
            }
            return;
        }
        if (mPos >= mRaw.size() || mRaw[mPos] != '(') {
            mFailed = true;
            return;
        }
        ++mPos;
        if (!mFormat.mBinary || count == 0) {
            std::vector<double> flat;
            ReadAsciiNumbers(')', want ? &flat : nullptr);
            if (!mFailed && want) {
                pOut->mUniform = false;
                pOut->mCount = count;
                pOut->mFlat = std::move(flat);
            }
            return;
        }
        if (comps == 0) {  // a binary list of a type this scanner cannot size
            mFailed = true;
            return;
        }
        const std::size_t entry = comps < 0 ? static_cast<std::size_t>(mFormat.mLabelBytes)
                                            : static_cast<std::size_t>(comps) *
                                                  static_cast<std::size_t>(mFormat.mScalarBytes);
        if (static_cast<std::uint64_t>(count) > (mRaw.size() - mPos) / entry) {
            mFailed = true;
            return;
        }
        if (want) {
            pOut->mUniform = false;
            pOut->mCount = count;
            pOut->mFlat.resize(static_cast<std::size_t>(count) * static_cast<std::size_t>(comps));
            const char* base = mRaw.data() + mPos;
            for (std::size_t i = 0; i < pOut->mFlat.size(); ++i)
                pOut->mFlat[i] = mFormat.mScalarBytes == 4
                                     ? static_cast<double>(read_le<float>(base + i * 4))
                                     : read_le<double>(base + i * 8);
        }
        mPos += static_cast<std::size_t>(count) * entry;
        if (mPos < mRaw.size() && mRaw[mPos] == ')')
            ++mPos;
        else
            mFailed = true;
    }

    /// A patch's `value` entry, through its `;`.
    void ReadValue(FoamField& rOut) {
        const std::string tok = ReadToken();
        if (tok == "uniform") {
            SkipSpace();
            std::vector<double> flat;
            if (mComponents == 1) {
                flat.push_back(detail::parse_double_prefix(ReadToken()));
            } else if (mPos < mRaw.size() && mRaw[mPos] == '(') {
                ++mPos;
                ReadAsciiNumbers(')', &flat);
            }
            if (!mFailed && flat.size() == static_cast<std::size_t>(mComponents)) {
                rOut.mUniform = true;
                rOut.mFlat = std::move(flat);
            }
        } else if (tok == "nonuniform") {
            ReadNonuniform(&rOut);
        }
        SkipValue();
    }

    std::string_view mRaw;
    FoamFormat mFormat;
    int mComponents;
    std::size_t mPos = 0;
    bool mFailed = false;
};

/// The per-patch values of a field file's `boundaryField` (see `FoamBoundaryScanner`).
std::vector<FoamPatchField> foam_read_boundary_field(const fs::path& rPath, int Components) {
    const detail::FileSource source = read_whole(rPath.string());
    const std::string_view raw = source.View();
    return FoamBoundaryScanner(raw, detect_format(raw), Components).Scan();
}

/// Parse a time directory's name as an OpenFOAM time value, requiring the
/// WHOLE name to be consumed (so `"0.1_backup"` is correctly not a time dir).
bool foam_parse_time_dir_name(const std::string& rName, double& rValue) {
    if (rName.empty())
        return false;
    const char* end = nullptr;
    const double v = detail::parse_double(rName.c_str(), end);
    if (end != rName.c_str() + rName.size())
        return false;
    rValue = v;
    return true;
}

/// One `<case_root>/<numeric>/` time directory: its parsed value and its
/// own (exact, on-disk) name -- kept together because a value re-formatted
/// back to text ("0.1" vs "0.100000") is not reliably the same string.
struct FoamTimeDir {
    double mValue = 0.0;
    std::string mName;
};

/// Time directories holding at least one regular file (a field), sorted
/// ascending by value. `0` is included only when it holds fields -- the
/// historical (field-free) `read_openfoam` never depended on a `0/`
/// directory existing at all.
std::vector<FoamTimeDir> foam_time_dirs(const fs::path& rCaseRoot, const std::string& rRegion = {}) {
    std::vector<FoamTimeDir> dirs;
    std::error_code ec;
    for (const auto& entry : fs::directory_iterator(rCaseRoot, ec)) {
        if (!entry.is_directory())
            continue;
        const std::string name = entry.path().filename().string();
        double t = 0.0;
        if (!foam_parse_time_dir_name(name, t))
            continue;
        std::error_code ec2;
        bool has_field = false;
        // a region's fields live in `<time>/<region>/`
        const fs::path fields_dir = rRegion.empty() ? entry.path() : entry.path() / rRegion;
        for (const auto& f : fs::directory_iterator(fields_dir, ec2)) {
            if (f.is_regular_file()) {
                has_field = true;
                break;
            }
        }
        if (has_field)
            dirs.push_back({t, name});
    }
    std::sort(dirs.begin(), dirs.end(),
              [](const FoamTimeDir& a, const FoamTimeDir& b) { return a.mValue < b.mValue; });
    return dirs;
}

/// Field file names directly inside a time directory (regular files only,
/// non-recursive -- `uniform/`, `polyMesh/` and other sub-directories a
/// moving-mesh case may carry there are not field files).
std::vector<std::string> foam_field_files(const fs::path& rTimeDir) {
    std::vector<std::string> names;
    std::error_code ec;
    for (const auto& entry : fs::directory_iterator(rTimeDir, ec))
        if (entry.is_regular_file())
            names.push_back(entry.path().filename().string());
    std::sort(names.begin(), names.end());
    return names;
}

/// One field of a decomposed case: each processor's `internalField` scattered to the global cell
/// (or point) numbering its `*ProcAddressing` names; NaN where no processor has a value.
FoamField foam_read_decomposed_field(const fs::path& rCaseRoot,
                                     const std::vector<std::size_t>& rProcs,
                                     const std::string& rTime, const std::string& rName,
                                     int Components, bool IsPoint, std::int64_t NumGlobal) {
    FoamField out;
    out.mCount = NumGlobal;
    const auto ncomp = static_cast<std::size_t>(Components);
    out.mFlat.assign(static_cast<std::size_t>(NumGlobal) * ncomp,
                     std::numeric_limits<double>::quiet_NaN());
    for (const std::size_t proc : rProcs) {
        const fs::path pdir = rCaseRoot / ("processor" + std::to_string(proc));
        const fs::path file = pdir / rTime / rName;
        if (!fs::exists(file))
            continue;
        const FoamField local = foam_read_internal_field(file, Components);
        const std::vector<std::int64_t> addr = read_int_list(
            pdir / "constant" / "polyMesh" / (IsPoint ? "pointProcAddressing" : "cellProcAddressing"));
        if (!local.mUniform && static_cast<std::size_t>(local.mCount) != addr.size())
            continue;
        for (std::size_t i = 0; i < addr.size(); ++i) {
            if (addr[i] < 0 || addr[i] >= NumGlobal)
                continue;
            for (std::size_t c = 0; c < ncomp; ++c)
                out.mFlat[static_cast<std::size_t>(addr[i]) * ncomp + c] =
                    local.mUniform ? local.mFlat[c] : local.mFlat[i * ncomp + c];
        }
    }
    return out;
}

/// The `class` entry of a field file's `FoamFile` header (`volScalarField`,
/// …), read the same cheap line-scan way `detect_format` reads `format`/
/// `arch` -- no full parse needed just to classify the field.
std::string foam_field_file_class(const fs::path& rPath) {
    auto f = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!f)
        return {};
    std::string line;
    while (std::getline(f, line)) {
        const std::string s = openfoam_strip(line);
        if (s.rfind("class", 0) == 0) {
            std::string rest = openfoam_strip(s.substr(std::strlen("class")));
            if (!rest.empty() && rest.back() == ';')
                rest.pop_back();
            return openfoam_strip(rest);
        }
        if (s == "}")
            break;
    }
    return {};
}

}  // namespace

Mesh read_openfoam(const std::string& rPathIn, OpenFoamInfo& rInfo) {
    return read_openfoam(rPathIn, ReadOptions{}, rInfo);
}

Mesh read_openfoam(const std::string& rPathIn, const ReadOptions& rOptions, OpenFoamInfo& rInfo) {
    // resolve polyMesh directory
    fs::path path(rPathIn);
    fs::path poly;
    if (path.extension() == ".foam") {
        fs::path c = path.parent_path() / "constant" / "polyMesh";
        if (fs::exists(c))
            poly = c;
    }
    if (poly.empty() && path.filename() == "polyMesh" && fs::is_directory(path))
        poly = path;
    if (poly.empty()) {
        for (const fs::path& c : {path / "constant" / "polyMesh", path / "polyMesh"}) {
            if (fs::exists(c)) {
                poly = c;
                break;
            }
        }
    }
    // Multi-region case (v11.4.0, roadmap §1 tier B2): no single
    // `constant/polyMesh`, but `<case>/constant/<region>/polyMesh` per
    // region. `case_root` is the case directory regardless of which of the
    // three `rPathIn` forms was given.
    const fs::path case_root = path.extension() == ".foam" ? path.parent_path() : path;
    if (poly.empty() && !rInfo.mRegion.empty()) {
        const fs::path c = case_root / "constant" / rInfo.mRegion / "polyMesh";
        if (fs::exists(c))
            poly = c;
        else
            throw ReadError(detail::format_compat("OpenFOAM: region '{}' has no {}", rInfo.mRegion,
                                                  c.string()));
    }
    if (poly.empty()) {
        // Regions are the `constant/<name>/polyMesh` directories; `regionProperties` (the solver's
        // region families) is not needed to find them.
        std::vector<std::string> regions;
        std::error_code ec;
        for (const auto& entry : fs::directory_iterator(case_root / "constant", ec)) {
            if (entry.is_directory() && fs::exists(entry.path() / "polyMesh"))
                regions.push_back(entry.path().filename().string());
        }
        if (!regions.empty() || fs::exists(case_root / "constant" / "regionProperties")) {
            if (regions.empty()) {
                // regionProperties without region directories: name the regions it lists
                try {
                    for (const OpenFoamRegionGroup& group :
                         read_openfoam_region_properties(case_root.string()))
                        regions.insert(regions.end(), group.mRegions.begin(), group.mRegions.end());
                } catch (const ReadError&) {
                }
            }
            std::sort(regions.begin(), regions.end());
            std::string joined;
            for (std::size_t i = 0; i < regions.size(); ++i)
                joined += (i ? ", " : "") + regions[i];
            throw ReadError(detail::format_compat(
                "'{}' is a multi-region case; set OpenFoamInfo::mRegion to one of: {}", rPathIn,
                joined));
        }
    }

    // Decomposed case (v11.4.0, roadmap §1 tier B2): no single
    // `constant/polyMesh` and no `<region>` selected, but `processorN`
    // directories each carrying their own `constant/polyMesh`.
    bool decomposed = false;
    std::vector<std::size_t> proc_ids;
    if (poly.empty()) {
        proc_ids = foam_processor_ids(case_root);
        decomposed = !proc_ids.empty();
        if (decomposed)
            poly = case_root;  // for the "Reading polyMesh from" log line below only
    }
    if (poly.empty())
        throw ReadError(detail::format_compat(
            "Could not locate polyMesh from '{}'. Expected <case>/constant/polyMesh/.", rPathIn));
    log::info("Reading polyMesh from {}", poly.string());

    RawPolyMesh raw;
    if (decomposed) {
        log::info("OpenFOAM: reconstructing {} from {} processor director{}", case_root.string(),
                  proc_ids.size(), proc_ids.size() == 1 ? "y" : "ies");
        raw = reconstruct_decomposed(case_root, proc_ids);
    } else {
        raw = read_raw_polymesh(poly);
    }
    P3 points = std::move(raw.mPoints);
    std::vector<Face> faces = std::move(raw.mFaces);
    std::vector<std::int64_t> owner = std::move(raw.mOwner);
    std::vector<std::int64_t> neighbour = std::move(raw.mNeighbour);
    std::vector<Patch> boundary = std::move(raw.mBoundary);
    const std::vector<std::int64_t> raw_face_of_global = std::move(raw.mFaceOfGlobal);

    std::int64_t owner_max = -1, neigh_max = -1;
    for (std::int64_t v : owner)
        owner_max = std::max(owner_max, v);
    for (std::int64_t v : neighbour)
        neigh_max = std::max(neigh_max, v);
    std::int64_t n_cells = owner.empty() ? 0 : std::max(owner_max, neigh_max) + 1;
    log::info("{} points, {} faces, {} cells, {} patches", points.size(), faces.size(), n_cells,
              boundary.size());

    // cell -> face ids
    std::vector<std::vector<std::int64_t>> cell_faces(static_cast<std::size_t>(n_cells));
    for (std::size_t fid = 0; fid < owner.size(); ++fid)
        cell_faces[static_cast<std::size_t>(owner[fid])].push_back(static_cast<std::int64_t>(fid));
    for (std::size_t fid = 0; fid < neighbour.size(); ++fid)
        if (neighbour[fid] >= 0)
            cell_faces[static_cast<std::size_t>(neighbour[fid])].push_back(
                static_cast<std::int64_t>(fid));

    // reconstruct volume cells
    std::vector<std::string> vol_order;
    std::map<std::string, std::vector<Face>> vol_buckets;
    // polyhedra grouped by unique node count -> "polyhedron<N>"
    std::vector<std::string> poly_order;
    std::map<std::string, std::vector<std::vector<Face>>> poly_buckets;

    // Per-cell geometric reconstruction is the expensive part and every cell
    // only reads faces/owner/points -> compute all cells in parallel into a
    // pre-sized result array, then do the (ordered) bucket grouping
    // sequentially.
    struct CellResult {
        std::string mType;         // "" = degenerate (skipped)
        Face mConn;                // named types
        std::vector<Face> mFaces;  // oriented faces, polyhedra only
    };
    std::vector<CellResult> results(static_cast<std::size_t>(n_cells));
    parallel_for(static_cast<std::size_t>(n_cells), [&](std::size_t cs) {
        const std::int64_t cid = static_cast<std::int64_t>(cs);
        std::vector<Face> oriented;
        for (std::int64_t fid : cell_faces[cs]) {
            Face f = faces[static_cast<std::size_t>(fid)];
            if (owner[static_cast<std::size_t>(fid)] != cid)
                std::reverse(f.begin(), f.end());
            oriented.push_back(std::move(f));
        }
        auto [mtype, conn] = reconstruct_cell(oriented, points);
        if (mtype == "polyhedron") {
            results[cs] = {"polyhedron", {}, std::move(oriented)};
        } else if (conn.empty()) {
            results[cs] = {};  // degenerate topology
        } else {
            results[cs] = {std::move(mtype), std::move(conn), {}};
        }
    });

    // Original OpenFOAM cell id -> (is_poly, bucket key, row within that
    // bucket), captured as cells are bucketed so zone regions can later
    // recover each cell's position in the final Mesh (see `orig_cell_to_global`
    // below). `row == npos` marks a skipped (degenerate) cell.
    constexpr std::size_t npos = static_cast<std::size_t>(-1);
    std::vector<std::tuple<bool, std::string, std::size_t>> placement(
        static_cast<std::size_t>(n_cells),
        std::tuple<bool, std::string, std::size_t>{false, "", npos});

    std::size_t n_skipped = 0;
    std::size_t n_polyhedra = 0;
    for (std::size_t cid = 0; cid < results.size(); ++cid) {
        auto& res = results[cid];
        if (res.mType == "polyhedron") {
            std::size_t nn = unique_node_count(res.mFaces);
            std::string key = "polyhedron" + std::to_string(nn);
            if (!poly_buckets.count(key))
                poly_order.push_back(key);
            placement[cid] = {true, key, poly_buckets[key].size()};
            poly_buckets[key].push_back(std::move(res.mFaces));
            ++n_polyhedra;
        } else if (res.mType.empty()) {
            ++n_skipped;
        } else {
            if (!vol_buckets.count(res.mType))
                vol_order.push_back(res.mType);
            placement[cid] = {false, res.mType, vol_buckets[res.mType].size()};
            vol_buckets[res.mType].push_back(std::move(res.mConn));
        }
    }
    if (n_skipped > 0)
        log::warn("{} cell(s) skipped (degenerate topology).", n_skipped);
    if (n_polyhedra > 0)
        log::info("{} general polyhedron cell(s) found.", n_polyhedra);

    // Block index of each bucket in the order blocks are about to be added
    // (volume types first, in `vol_order`, then polyhedron buckets in
    // `poly_order`) -- the boundary 2D blocks added further down come after
    // both, so this table stays valid for `detail::block_bases` once those
    // volume/polyhedron blocks are on the mesh.
    std::unordered_map<std::string, std::size_t> vol_block_index, poly_block_index;
    for (std::size_t k = 0; k < vol_order.size(); ++k)
        vol_block_index[vol_order[k]] = k;
    for (std::size_t k = 0; k < poly_order.size(); ++k)
        poly_block_index[poly_order[k]] = vol_order.size() + k;

    Mesh mesh;
    std::size_t npts = points.size();
    {
        NDArray pts(DType::Float64, {npts, 3});
        double* pdst = pts.As<double>();
        parallel_for(npts, [&](std::size_t i) {
            for (std::size_t j = 0; j < 3; ++j)
                pdst[i * 3 + j] = points[i][j];
        });
        mesh.AssignPoints(std::move(pts));
    }

    std::vector<NDArray> cell_tags;  // one per block, in final block order

    // rectangular volume blocks
    for (const std::string& t : vol_order) {
        const auto& rows = vol_buckets[t];
        std::size_t nc = rows.size();
        std::size_t k = nc ? rows[0].size() : 0;
        NDArray data(DType::Int64, {nc, k});
        std::int64_t* dp = data.As<std::int64_t>();
        parallel_for(nc, [&](std::size_t r) {
            for (std::size_t c = 0; c < k; ++c)
                dp[r * k + c] = rows[r][c];
        });
        mesh.AddCellBlock(t, std::move(data));
        cell_tags.emplace_back(DType::Int64, std::vector<std::size_t>{nc});  // zeros
    }
    // ragged polyhedron blocks
    for (const std::string& key : poly_order) {
        // Flattened once into the CSR the mesh stores (no per-face copy).
        auto& bucket = poly_buckets[key];
        const std::size_t nc = bucket.size();
        std::vector<std::int64_t> flat, row_offsets, face_offsets;
        detail::csr_from_cells(bucket, flat, row_offsets, face_offsets);
        decltype(poly_buckets)::mapped_type().swap(bucket);
        mesh.AddPolyhedronBlock(key, std::move(flat), std::move(row_offsets),
                                std::move(face_offsets));
        cell_tags.emplace_back(DType::Int64, std::vector<std::size_t>{nc});  // zeros
    }

    // Original OpenFOAM cell id -> global (block-major) cell index, used by
    // both the zones-as-regions block below and time-directory field
    // attachment further down. `detail::block_bases(mesh)` only needs to be
    // right for the blocks added so far (volume + polyhedron); boundary (2D)
    // blocks are always appended after, so their prefix sums never change
    // what is computed here.
    std::vector<std::int64_t> orig_cell_to_global(static_cast<std::size_t>(n_cells), -1);
    {
        const std::vector<std::int64_t> bases = detail::block_bases(mesh);
        for (std::size_t cid = 0; cid < placement.size(); ++cid) {
            const auto& [is_poly, key, row] = placement[cid];
            if (row == npos)
                continue;
            const std::size_t block = is_poly ? poly_block_index.at(key) : vol_block_index.at(key);
            orig_cell_to_global[cid] =
                detail::block_row_to_global(bases, block, static_cast<std::int64_t>(row));
        }
    }

    // ---- zones as named regions (cellZones/faceZones/pointZones) ----
    {
        // Local facet index of face `fid` within its owner cell `cid`,
        // matching the order `AddPolyhedronBlock`/`AddCellBlock` end up
        // storing: for a polyhedron, that is simply `fid`'s position in
        // `cell_faces[cid]` (the same order the per-cell reconstruction loop
        // above walked to build `oriented`, and that `AddPolyhedronBlock`
        // preserves verbatim); a named type's reconstruction instead rewinds
        // into `cell_faces.hpp`'s canonical per-type order, so the position
        // is recovered by matching each canonical face's corner *set*
        // against `fid`'s -- orientation- and start-point-independent.
        auto local_facet = [&](std::size_t cid, std::int64_t fid) -> int {
            const auto& [is_poly, key, row] = placement[cid];
            if (row == npos)
                return -1;
            if (is_poly) {
                const auto& cf = cell_faces[cid];
                for (std::size_t k = 0; k < cf.size(); ++k)
                    if (cf[k] == fid)
                        return static_cast<int>(k);
                return -1;
            }
            const Face& conn = vol_buckets.at(key)[row];
            const auto& facedefs = detail::cell_faces(cell_type_from_name(key));
            const Face& fnodes = faces[static_cast<std::size_t>(fid)];
            const std::unordered_set<std::int64_t> target(fnodes.begin(), fnodes.end());
            for (std::size_t k = 0; k < facedefs.size(); ++k) {
                std::unordered_set<std::int64_t> cand;
                for (int c = 0; c < facedefs[k].mNumCorners; ++c)
                    cand.insert(conn[facedefs[k].mNodes[c]]);
                if (cand == target)
                    return static_cast<int>(k);
            }
            return -1;
        };

        // A decomposed case keeps its zones per processor, in processor-local ids.
        auto read_zone_file = [&](const char* pFile, const char* pKey) {
            if (!decomposed)
                return foam_read_zone_file(poly / pFile, pKey);
            const std::string file = pFile;
            return foam_read_decomposed_zones(case_root, proc_ids, pFile, pKey,
                                              file == "cellZones"    ? "cellProcAddressing"
                                              : file == "pointZones" ? "pointProcAddressing"
                                                                     : "faceProcAddressing",
                                              file == "faceZones" ? &raw_face_of_global : nullptr);
        };

        for (const Zone& z : read_zone_file("cellZones", "cellLabels")) {
            std::size_t n_dropped = 0;
            std::vector<std::int64_t> entries;
            entries.reserve(z.mIds.size());
            for (std::int64_t cid : z.mIds) {
                if (cid < 0 || cid >= n_cells ||
                    orig_cell_to_global[static_cast<std::size_t>(cid)] < 0) {
                    ++n_dropped;
                    continue;
                }
                entries.push_back(orig_cell_to_global[static_cast<std::size_t>(cid)]);
            }
            if (n_dropped > 0)
                log::warn("OpenFOAM: cellZone '{}' drops {} entr{} (degenerate cell)", z.mName,
                          n_dropped, n_dropped == 1 ? "y" : "ies");
            NDArray arr = NDArray::Uninit(DType::Int64, {entries.size()});
            std::copy(entries.begin(), entries.end(), arr.As<std::int64_t>());
            mesh.AddRegion(Region(z.mName, RegionKind::Cell, std::move(arr)));
        }

        for (const Zone& z : read_zone_file("pointZones", "pointLabels")) {
            std::vector<std::int64_t> entries;
            entries.reserve(z.mIds.size());
            for (std::int64_t pid : z.mIds)
                if (pid >= 0 && static_cast<std::size_t>(pid) < npts)
                    entries.push_back(pid);
            NDArray arr = NDArray::Uninit(DType::Int64, {entries.size()});
            std::copy(entries.begin(), entries.end(), arr.As<std::int64_t>());
            mesh.AddRegion(Region(z.mName, RegionKind::Point, std::move(arr)));
        }

        for (const Zone& z : read_zone_file("faceZones", "faceLabels")) {
            std::size_t n_dropped = 0;
            std::vector<std::int64_t> pairs;
            pairs.reserve(z.mIds.size() * 2);
            for (std::int64_t fid : z.mIds) {
                if (fid < 0 || static_cast<std::size_t>(fid) >= faces.size()) {
                    ++n_dropped;
                    continue;
                }
                const std::int64_t cid = owner[static_cast<std::size_t>(fid)];
                const std::int64_t global = orig_cell_to_global[static_cast<std::size_t>(cid)];
                const int facet = local_facet(static_cast<std::size_t>(cid), fid);
                if (global < 0 || facet < 0) {
                    ++n_dropped;
                    continue;
                }
                pairs.push_back(global);
                pairs.push_back(facet);
            }
            if (n_dropped > 0)
                log::warn("OpenFOAM: faceZone '{}' drops {} entr{} (degenerate owner cell)",
                          z.mName, n_dropped, n_dropped == 1 ? "y" : "ies");
            NDArray arr = NDArray::Uninit(DType::Int64, {pairs.size() / 2, 2});
            std::copy(pairs.begin(), pairs.end(), arr.As<std::int64_t>());
            mesh.AddRegion(Region(z.mName, RegionKind::Side, std::move(arr)));
        }
    }

    // boundary cells grouped by size, with patch family tags
    std::map<int, std::vector<Face>> bysize;  // 3 -> triangles, 4 -> quads
    std::map<int, std::vector<std::int64_t>> tagsize;
    std::vector<Face> poly_faces;
    std::vector<std::int64_t> poly_tags;
    // Where each boundary-file face went: its vertex count (which names its block) and its row
    // in that block, so a field's per-patch values can be attached to the right cell.
    struct BoundaryFaceDest {
        std::size_t mSize = 0;
        // Not `npos`: MSVC rejects a local class whose initializer names a local of the
        // enclosing function (C2672 on `std::construct_at`).
        std::size_t mRow = static_cast<std::size_t>(-1);
    };
    std::vector<std::vector<BoundaryFaceDest>> patch_dest(boundary.size());
    std::map<std::size_t, std::size_t> dest_rows;  // vertex count -> rows handed out so far
    for (std::size_t pidx = 0; pidx < boundary.size(); ++pidx) {
        std::int64_t fam = -(static_cast<std::int64_t>(pidx) + 1);
        rInfo.mCellTags[fam] = {boundary[pidx].mName};
        if (!boundary[pidx].mType.empty())
            rInfo.mPatchTypes[fam] = boundary[pidx].mType;
        if (boundary[pidx].mNFaces > 0)
            patch_dest[pidx].resize(
                std::min(static_cast<std::size_t>(boundary[pidx].mNFaces), faces.size()));
        for (std::int64_t fid = boundary[pidx].mStartFace;
             fid < boundary[pidx].mStartFace + boundary[pidx].mNFaces; ++fid) {
            if (fid < 0 || static_cast<std::size_t>(fid) >= faces.size())
                continue;
            const Face& f = faces[static_cast<std::size_t>(fid)];
            const auto offset = static_cast<std::size_t>(fid - boundary[pidx].mStartFace);
            if (offset < patch_dest[pidx].size())
                patch_dest[pidx][offset] = {f.size(), dest_rows[f.size()]++};
            if (f.size() == 3) {
                bysize[3].push_back(f);
                tagsize[3].push_back(fam);
            } else if (f.size() == 4) {
                bysize[4].push_back(f);
                tagsize[4].push_back(fam);
            } else {
                poly_faces.push_back(f);
                poly_tags.push_back(fam);
            }
        }
    }
    auto add_boundary_block = [&](const std::string& type, const std::vector<Face>& rows,
                                  const std::vector<std::int64_t>& tags) {
        std::size_t nc = rows.size();
        std::size_t k = nc ? rows[0].size() : 0;
        NDArray data(DType::Int64, {nc, k});
        NDArray tag(DType::Int64, {nc});
        std::int64_t* dp = data.As<std::int64_t>();
        std::int64_t* tp = tag.As<std::int64_t>();
        parallel_for(nc, [&](std::size_t r) {
            for (std::size_t c = 0; c < k; ++c)
                dp[r * k + c] = rows[r][c];
            tp[r] = tags[r];
        });
        mesh.AddCellBlock(type, std::move(data));
        cell_tags.push_back(std::move(tag));
    };
    if (!bysize[3].empty())
        add_boundary_block("triangle", bysize[3], tagsize[3]);
    if (!bysize[4].empty())
        add_boundary_block("quad", bysize[4], tagsize[4]);
    if (!poly_faces.empty()) {
        // group boundary polygons by vertex count -> polygon<N>
        std::map<std::size_t, std::vector<Face>> by_n;
        std::map<std::size_t, std::vector<std::int64_t>> tag_n;
        for (std::size_t i = 0; i < poly_faces.size(); ++i) {
            by_n[poly_faces[i].size()].push_back(poly_faces[i]);
            tag_n[poly_faces[i].size()].push_back(poly_tags[i]);
        }
        for (auto& kv : by_n)
            add_boundary_block("polygon" + std::to_string(kv.first), kv.second, tag_n[kv.first]);
    }

    if (!cell_tags.empty())
        mesh.AddCellData("cell_tags", std::move(cell_tags));

    // Boundary-face block of each vertex count, and each patch by name (first wins).
    std::map<std::size_t, std::size_t> size_block;
    for (std::size_t b = vol_order.size() + poly_order.size(); b < mesh.NumCellBlocks(); ++b) {
        const std::string type(mesh.Cells(b).Type());
        size_block[type == "triangle" ? 3
                   : type == "quad"   ? 4
                                      : static_cast<std::size_t>(std::atoll(type.c_str() + 7))] = b;
    }
    std::unordered_map<std::string, std::size_t> patch_by_name;
    for (std::size_t pidx = 0; pidx < boundary.size(); ++pidx)
        patch_by_name.emplace(boundary[pidx].mName, pidx);

    // ---- time-directory fields (v11.4.0, roadmap §1 tier B2) ----
    if (rOptions.WantsAnyData()) {
        // A decomposed case lists its times under the first processor; a region's fields live
        // in `<time>/<region>/`.
        const fs::path time_root =
            decomposed ? case_root / ("processor" + std::to_string(proc_ids.front())) : case_root;
        const std::string time_region = decomposed ? std::string() : rInfo.mRegion;
        const std::vector<FoamTimeDir> time_dirs = foam_time_dirs(time_root, time_region);
        if (!time_dirs.empty()) {
            const std::size_t step = rOptions.ResolveTimeStep(time_dirs.size());
            const fs::path time_dir = time_region.empty()
                                          ? time_root / time_dirs[step].mName
                                          : time_root / time_dirs[step].mName / time_region;
            for (const std::string& field_name : foam_field_files(time_dir)) {
                if (!rOptions.WantsArray(field_name))
                    continue;
                const fs::path field_path = time_dir / field_name;
                const std::string cls = foam_field_file_class(field_path);
                const FoamFieldClass fc = foam_field_class(cls);
                if (!fc.mSupported) {
                    log::warn(
                        "OpenFOAM: field '{}' has class '{}', which has no point/cell home; "
                        "skipped",
                        field_name, cls.empty() ? "?" : cls);
                    continue;
                }
                const FoamField values =
                    decomposed
                        ? foam_read_decomposed_field(case_root, proc_ids, time_dirs[step].mName,
                                                     field_name, fc.mComponents, fc.mIsPoint,
                                                     fc.mIsPoint ? static_cast<std::int64_t>(npts)
                                                                 : n_cells)
                        : foam_read_internal_field(field_path, fc.mComponents);
                const std::size_t ncomp = static_cast<std::size_t>(fc.mComponents);
                // A value list shorter than its declared count (a line that
                // gave too few numbers, a truncated file) is refused here, not
                // read past: the checks below compare the declared count only.
                if (values.mUniform ? values.mFlat.size() < ncomp
                                    : values.mFlat.size() !=
                                          static_cast<std::size_t>(values.mCount) * ncomp) {
                    log::warn("OpenFOAM: field '{}' holds {} number(s), which does not match its "
                              "declaration; skipped",
                              field_name, values.mFlat.size());
                    continue;
                }

                if (fc.mIsPoint) {
                    if (!values.mUniform && static_cast<std::size_t>(values.mCount) != npts) {
                        log::warn(
                            "OpenFOAM: field '{}' has {} value(s), expected {} point(s); skipped",
                            field_name, values.mCount, npts);
                        continue;
                    }
                    NDArray arr = NDArray::Uninit(DType::Float64,
                                                  {npts, static_cast<std::size_t>(fc.mComponents)});
                    double* dst = arr.As<double>();
                    if (values.mUniform) {
                        for (std::size_t i = 0; i < npts; ++i)
                            for (int c = 0; c < fc.mComponents; ++c)
                                dst[i * static_cast<std::size_t>(fc.mComponents) +
                                    static_cast<std::size_t>(c)] =
                                    values.mFlat[static_cast<std::size_t>(c)];
                    } else {
                        std::copy(values.mFlat.begin(), values.mFlat.end(), dst);
                    }
                    mesh.AddPointData(field_name, std::move(arr));
                    continue;
                }

                if (!values.mUniform && values.mCount != n_cells) {
                    log::warn("OpenFOAM: field '{}' has {} value(s), expected {} cell(s); skipped",
                              field_name, values.mCount, n_cells);
                    continue;
                }
                std::vector<NDArray> blocks(mesh.NumCellBlocks());
                for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b) {
                    NDArray arr = NDArray::Uninit(
                        DType::Float64,
                        {mesh.Cells(b).NumCells(), static_cast<std::size_t>(fc.mComponents)});
                    double* dst = arr.As<double>();
                    std::fill(dst, dst + arr.Size(), std::numeric_limits<double>::quiet_NaN());
                    blocks[b] = std::move(arr);
                }
                for (std::size_t cid = 0; cid < placement.size(); ++cid) {
                    const auto& [is_poly, key, row] = placement[cid];
                    if (row == npos)
                        continue;
                    const std::size_t block =
                        is_poly ? poly_block_index.at(key) : vol_block_index.at(key);
                    double* dst = blocks[block].As<double>();
                    for (int c = 0; c < fc.mComponents; ++c)
                        dst[row * static_cast<std::size_t>(fc.mComponents) +
                            static_cast<std::size_t>(c)] =
                            values.mUniform
                                ? values.mFlat[static_cast<std::size_t>(c)]
                                : values.mFlat[cid * static_cast<std::size_t>(fc.mComponents) +
                                               static_cast<std::size_t>(c)];
                }

                // Boundary-face blocks: the `value` of each patch's `boundaryField` entry (a patch
                // without one, or a processor patch, keeps the NaN). A decomposed case reads each
                // processor's copy of the patch face by face.
                auto set_face = [&](std::size_t pidx, std::size_t face, const FoamField& rValue,
                                    std::size_t index) {
                    const BoundaryFaceDest& d = patch_dest[pidx][face];
                    if (d.mRow == npos)
                        return;
                    double* dst = blocks[size_block.at(d.mSize)].As<double>();
                    for (std::size_t c = 0; c < ncomp; ++c)
                        dst[d.mRow * ncomp + c] =
                            rValue.mUniform ? rValue.mFlat[c] : rValue.mFlat[index * ncomp + c];
                };
                auto usable = [&](const FoamField& rValue) {
                    return rValue.mUniform ? rValue.mFlat.size() == ncomp
                                           : rValue.mFlat.size() ==
                                                 static_cast<std::size_t>(rValue.mCount) * ncomp;
                };
                if (!patch_dest.empty() && !size_block.empty()) {
                    if (!decomposed) {
                        for (const FoamPatchField& pf :
                             foam_read_boundary_field(field_path, fc.mComponents)) {
                            const auto it = patch_by_name.find(pf.mName);
                            if (it == patch_by_name.end() || pf.mValue.mFlat.empty() ||
                                !usable(pf.mValue))
                                continue;
                            const std::size_t pidx = it->second;
                            if (!pf.mValue.mUniform && pf.mValue.mCount != boundary[pidx].mNFaces) {
                                log::warn(
                                    "OpenFOAM: field '{}' patch '{}' has {} value(s), expected {} "
                                    "face(s); skipped",
                                    field_name, pf.mName, pf.mValue.mCount, boundary[pidx].mNFaces);
                                continue;
                            }
                            for (std::size_t j = 0; j < patch_dest[pidx].size(); ++j)
                                set_face(pidx, j, pf.mValue, j);
                        }
                    } else {
                        std::vector<std::vector<FoamPatchField>> per_proc(proc_ids.size());
                        for (std::size_t k = 0; k < proc_ids.size(); ++k) {
                            const fs::path file = case_root /
                                                  ("processor" + std::to_string(proc_ids[k])) /
                                                  time_dirs[step].mName / field_name;
                            if (fs::exists(file))
                                per_proc[k] = foam_read_boundary_field(file, fc.mComponents);
                        }
                        for (std::size_t pidx = 0;
                             pidx < boundary.size() && pidx < raw.mBoundarySource.size(); ++pidx) {
                            std::vector<const FoamPatchField*> local(proc_ids.size(), nullptr);
                            for (std::size_t k = 0; k < proc_ids.size(); ++k)
                                for (const FoamPatchField& pf : per_proc[k])
                                    if (pf.mName == boundary[pidx].mName && usable(pf.mValue) &&
                                        !pf.mValue.mFlat.empty()) {
                                        local[k] = &pf;
                                        break;
                                    }
                            const auto& sources = raw.mBoundarySource[pidx];
                            for (std::size_t j = 0;
                                 j < sources.size() && j < patch_dest[pidx].size(); ++j) {
                                const FoamPatchField* pf = local[sources[j].first];
                                const std::int64_t at = sources[j].second;
                                if (pf == nullptr || at < 0 ||
                                    (!pf->mValue.mUniform && at >= pf->mValue.mCount))
                                    continue;
                                set_face(pidx, j, pf->mValue, static_cast<std::size_t>(at));
                            }
                        }
                    }
                }
                mesh.AddCellData(field_name, std::move(blocks));
            }
        }
    }

    return mesh;
}

MeshMetadata read_openfoam_metadata(const std::string& rPathIn, const ReadOptions& rOptions) {
    const fs::path path(rPathIn);
    const fs::path case_root = path.extension() == ".foam" ? path.parent_path() : path;

    std::vector<double> times;
    const std::vector<std::size_t> procs = foam_processor_ids(case_root);
    const bool decomposed_case = !fs::exists(case_root / "constant" / "polyMesh") && !procs.empty();
    for (const FoamTimeDir& t : foam_time_dirs(
             decomposed_case ? case_root / ("processor" + std::to_string(procs.front())) : case_root))
        times.push_back(t.mValue);

    OpenFoamInfo info;
    MeshMetadata meta = metadata_from_mesh(read_openfoam(rPathIn, rOptions, info));
    meta.mTimeValues = std::move(times);
    meta.mFellBackToFullRead = true;
    return meta;
}

// ==========================================================================
//                                  WRITER
// ==========================================================================

namespace {

/// One boundary patch as it will be written.
struct FoamPatchOut {
    std::string mName;
    std::string mType = "patch";
    std::int64_t mNFaces = 0;
    std::int64_t mStartFace = 0;
    std::string mExtra;  ///< further dictionary lines (a processor patch's neighbour)
};

/// The written face order plus the patch table describing its tail.
struct FoamFaceOrder {
    std::vector<std::int64_t> mNewToOld;  ///< written face id -> GlobalFaces id
    std::int64_t mNumInternal = 0;
    std::vector<FoamPatchOut> mPatches;  ///< ascending mStartFace
};

/**
 * @brief Patch types that survive a round trip unchanged.
 *
 * Everything else needs companion dictionary entries `OpenFoamInfo` does not
 * carry -- `cyclic`/`cyclicAMI` need `neighbourPatch`, `processor` needs
 * `myProcNo`/`neighbProcNo`, `mapped*` needs `sample*` -- and OpenFOAM refuses
 * to *load* a case whose patch declares such a type without them. Downgrading
 * to `patch` yields a case that loads and solves with visibly wrong boundary
 * conditions, which is strictly better than one that does not open.
 */
bool foam_type_is_self_contained(const std::string& rType) {
    return rType == "patch" || rType == "wall" || rType == "symmetry" || rType == "symmetryPlane" ||
           rType == "empty" || rType == "wedge";
}

/// Resolve the polyMesh directory. One function for both directions, so the
/// reader's resolution and the writer's cannot drift apart.
fs::path foam_polymesh_dir(const fs::path& rPath, bool ForWrite) {
    if (rPath.extension() == ".foam") {
        const fs::path c = rPath.parent_path() / "constant" / "polyMesh";
        if (ForWrite || fs::exists(c))
            return c;
    }
    if (rPath.filename() == "polyMesh" && (ForWrite || fs::is_directory(rPath)))
        return rPath;
    if (!ForWrite) {
        for (const fs::path& c : {rPath / "constant" / "polyMesh", rPath / "polyMesh"}) {
            if (fs::exists(c))
                return c;
        }
        return {};
    }
    return rPath / "constant" / "polyMesh";
}

/// Standard FoamFile header. `detect_format` reads only `format` and `arch`,
/// but the rest is what makes the file legible to OpenFOAM itself.
void foam_write_header(std::ostream& rOs, const std::string& rClass, const std::string& rObject,
                       const OpenFoamWriteOptions& rOpts,
                       const char* pLocation = "constant/polyMesh") {
    // The credit cell is fixed-width (48 chars before the closing box edge) so
    // the banner stays aligned regardless of how long the release string is.
    std::string credit = detail::provenance_lines(detail::SlotTier::Bounded)[0];
    if (credit.size() < 48)
        credit.append(48 - credit.size(), ' ');
    rOs << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
           "| =========                 |                                                 |\n"
           "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
           "|  \\\\    /   O peration     |                                                 |\n"
           "|   \\\\  /    A nd           | "
        << credit
        << "|\n"
           "|    \\\\/     M anipulation  |                                                 |\n"
           "\\*---------------------------------------------------------------------------*/\n"
           "FoamFile\n"
           "{\n"
           "    version     2.0;\n"
           "    format      "
        << (rOpts.mBinary ? "binary" : "ascii") << ";\n";
    if (rOpts.mBinary)
        rOs << "    arch        \"LSB;label=" << rOpts.mLabelBits << ";scalar=" << rOpts.mScalarBits
            << "\";\n";
    rOs << "    class       " << rClass
        << ";\n"
           "    location    \""
        << pLocation
        << "\";\n"
           "    object      "
        << rObject
        << ";\n"
           "}\n"
           "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
}

// ---- binary body writers ----
// Raw little-endian bytes, matching read_binary_points/read_binary_labels/
// read_binary_faces's own expectations exactly (this file's read side, not
// invented independently) -- see foam_write_header's arch line for the
// widths a reader needs to parse these back.

void foam_write_binary_label(std::ostream& rOs, std::int64_t v, int LabelBits) {
    if (LabelBits == 32) {
        const std::int32_t v32 = static_cast<std::int32_t>(v);
        rOs.write(reinterpret_cast<const char*>(&v32), sizeof(v32));
    } else {
        rOs.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }
}

void foam_write_binary_scalar(std::ostream& rOs, double v, int ScalarBits) {
    if (ScalarBits == 32) {
        const float f = static_cast<float>(v);
        rOs.write(reinterpret_cast<const char*>(&f), sizeof(f));
    } else {
        rOs.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }
}

/**
 * @brief Assign every boundary face to a patch.
 *
 * @return per-`GlobalFaces`-face patch index, `-1` for internal or unassigned.
 */
struct FoamPatchAssignment {
    std::vector<std::int64_t> mFacePatch;
    /// per face, the global (block-major) index of the boundary-face cell that claimed it, or -1
    /// -- where a patch face's cell-data values come from
    std::vector<std::int64_t> mFaceSource;
    std::vector<FoamPatchOut> mPatches;
    std::int64_t mNumOrphan = 0;          ///< 2D cell matching no face at all
    std::int64_t mNumInternalTagged = 0;  ///< 2D cell matching an INTERNAL face
};

FoamPatchAssignment foam_assign_patches(const Mesh& rMesh, const detail::GlobalFaces& rFaces,
                                        const OpenFoamInfo& rInfo) {
    FoamPatchAssignment out;
    out.mFacePatch.assign(rFaces.NumFaces(), -1);
    out.mFaceSource.assign(rFaces.NumFaces(), -1);
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);

    // Family ids in the reader's own order (ascending -fam == ascending patch
    // index), so an OpenFOAM round trip preserves the boundary file's order.
    std::vector<std::int64_t> fams;
    fams.reserve(rInfo.mCellTags.size());
    for (const auto& kv : rInfo.mCellTags)
        fams.push_back(kv.first);
    std::sort(fams.begin(), fams.end(), [](std::int64_t a, std::int64_t b) { return a > b; });

    std::unordered_map<std::int64_t, std::size_t> fam_to_patch;
    for (std::int64_t fam : fams) {
        FoamPatchOut p;
        const auto& names = rInfo.mCellTags.at(fam);
        p.mName = names.empty() ? ("patch" + std::to_string(-fam)) : names.front();
        const auto it = rInfo.mPatchTypes.find(fam);
        if (it != rInfo.mPatchTypes.end() && !it->second.empty()) {
            if (foam_type_is_self_contained(it->second)) {
                p.mType = it->second;
            } else {
                log::warn(
                    "OpenFOAM: patch '{}' has type '{}', which needs dictionary entries meshio++ "
                    "does not carry; writing 'patch' instead so the case still loads",
                    p.mName, it->second);
                detail::provenance_note(
                    "data-dropped", "patch '" + p.mName + "' downgraded from type '" + it->second +
                                        "' to 'patch' -- its companion dictionary entries "
                                        "are not carried");
            }
        }
        fam_to_patch[fam] = out.mPatches.size();
        out.mPatches.push_back(std::move(p));
    }

    const detail::FaceLookup lookup(rFaces);
    const bool have_tags = rMesh.HasCellData("cell_tags");

    for (std::size_t block : rFaces.mNonCellBlocks) {
        const auto cb = rMesh.Cells(block);
        if (cell_type_dimension(cell_type_from_name(cb.Type())) != 2)
            continue;
        const NDArray* tags =
            have_tags && rMesh.CellDataNumBlocks("cell_tags") == rMesh.NumCellBlocks()
                ? &rMesh.CellData("cell_tags", block)
                : nullptr;

        std::vector<std::int64_t> ids;
        std::optional<detail::Int64View> conn_values;
        if (!cb.IsRagged())
            conn_values.emplace(cb.Conn());
        std::optional<detail::Int64View> tag_values;
        if (tags)
            tag_values.emplace(*tags);
        for (std::size_t i = 0; i < cb.NumCells(); ++i) {
            ids.clear();
            if (cb.IsRagged()) {
                const std::size_t n = cb.RowSize(i);
                const std::int64_t* row = cb.Row(i);
                ids.assign(row, row + n);
            } else {
                const std::size_t npc = cb.NodesPerCell();
                for (std::size_t k = 0; k < npc; ++k)
                    ids.push_back((*conn_values)[i * npc + k]);
            }
            const std::int64_t fid = lookup.Find(ids.data(), ids.size());
            if (fid < 0) {
                ++out.mNumOrphan;
                continue;
            }
            if (rFaces.mNeighbour[static_cast<std::size_t>(fid)] >= 0) {
                // A face between two cells cannot also be a patch member --
                // OpenFOAM would refuse the case. Legitimate input (an interior
                // baffle), so drop it rather than throw.
                ++out.mNumInternalTagged;
                continue;
            }
            if (out.mFacePatch[static_cast<std::size_t>(fid)] >= 0)
                continue;  // first claim wins
            // A claim by a cell the patch table cannot place (no tag) still carries its values
            // to `defaultFaces`; a tagged claim below replaces it.
            std::int64_t& rSource = out.mFaceSource[static_cast<std::size_t>(fid)];
            if (rSource < 0)
                rSource = detail::block_row_to_global(bases, block, static_cast<std::int64_t>(i));
            std::int64_t fam = 0;
            if (tag_values && i < tags->Shape()[0])
                fam = (*tag_values)[i];
            const auto it = fam_to_patch.find(fam);
            if (fam < 0 && it != fam_to_patch.end()) {
                out.mFacePatch[static_cast<std::size_t>(fid)] =
                    static_cast<std::int64_t>(it->second);
                rSource = detail::block_row_to_global(bases, block, static_cast<std::int64_t>(i));
            }
        }
    }

    // Everything still unassigned joins `defaultFaces` -- blockMesh's own name
    // for exactly this, so the result is a valid single-patch case rather than
    // an error or an invented decomposition.
    std::int64_t n_unassigned = 0;
    for (std::size_t f = 0; f < rFaces.NumFaces(); ++f)
        if (rFaces.mNeighbour[f] < 0 && out.mFacePatch[f] < 0)
            ++n_unassigned;
    if (n_unassigned > 0) {
        const std::int64_t idx = static_cast<std::int64_t>(out.mPatches.size());
        FoamPatchOut p;
        p.mName = "defaultFaces";
        p.mType = "patch";
        out.mPatches.push_back(std::move(p));
        for (std::size_t f = 0; f < rFaces.NumFaces(); ++f)
            if (rFaces.mNeighbour[f] < 0 && out.mFacePatch[f] < 0)
                out.mFacePatch[f] = idx;
    }

    // An empty patch is legal but checkMesh flags it, and its emptiness always
    // means something went wrong upstream.
    std::vector<std::int64_t> counts(out.mPatches.size(), 0);
    for (std::size_t f = 0; f < rFaces.NumFaces(); ++f)
        if (out.mFacePatch[f] >= 0)
            ++counts[static_cast<std::size_t>(out.mFacePatch[f])];
    std::vector<FoamPatchOut> kept;
    std::vector<std::int64_t> remap(out.mPatches.size(), -1);
    for (std::size_t p = 0; p < out.mPatches.size(); ++p) {
        if (counts[p] == 0) {
            log::warn("OpenFOAM: patch '{}' has no faces and is not written",
                      out.mPatches[p].mName);
            detail::provenance_note("regions-dropped", "patch '" + out.mPatches[p].mName +
                                                           "' dropped -- it has no faces");
            continue;
        }
        remap[p] = static_cast<std::int64_t>(kept.size());
        kept.push_back(out.mPatches[p]);
    }
    for (std::size_t f = 0; f < rFaces.NumFaces(); ++f)
        if (out.mFacePatch[f] >= 0)
            out.mFacePatch[f] = remap[static_cast<std::size_t>(out.mFacePatch[f])];
    out.mPatches = std::move(kept);
    return out;
}

/// Order the faces the way OpenFOAM requires.
FoamFaceOrder foam_order_faces(const detail::GlobalFaces& rFaces,
                               const FoamPatchAssignment& rAssign) {
    FoamFaceOrder order;
    order.mPatches = rAssign.mPatches;

    std::vector<std::int64_t> internal, boundary;
    for (std::size_t f = 0; f < rFaces.NumFaces(); ++f)
        (rFaces.mNeighbour[f] >= 0 ? internal : boundary).push_back(static_cast<std::int64_t>(f));

    // Upper-triangular order. The face id is a THIRD key, not decoration: two
    // cells can share two distinct faces, so (owner, neighbour) alone is not a
    // strict weak ordering and std::sort would be undefined behaviour on it.
    std::sort(internal.begin(), internal.end(), [&](std::int64_t a, std::int64_t b) {
        const std::size_t ia = static_cast<std::size_t>(a), ib = static_cast<std::size_t>(b);
        if (rFaces.mOwner[ia] != rFaces.mOwner[ib])
            return rFaces.mOwner[ia] < rFaces.mOwner[ib];
        if (rFaces.mNeighbour[ia] != rFaces.mNeighbour[ib])
            return rFaces.mNeighbour[ia] < rFaces.mNeighbour[ib];
        return a < b;
    });
    // Boundary faces contiguous per patch, patches in table order.
    std::stable_sort(boundary.begin(), boundary.end(), [&](std::int64_t a, std::int64_t b) {
        return rAssign.mFacePatch[static_cast<std::size_t>(a)] <
               rAssign.mFacePatch[static_cast<std::size_t>(b)];
    });

    order.mNumInternal = static_cast<std::int64_t>(internal.size());
    order.mNewToOld = std::move(internal);
    order.mNewToOld.insert(order.mNewToOld.end(), boundary.begin(), boundary.end());

    std::int64_t start = order.mNumInternal;
    for (std::size_t p = 0; p < order.mPatches.size(); ++p) {
        std::int64_t n = 0;
        for (std::int64_t f : boundary)
            if (rAssign.mFacePatch[static_cast<std::size_t>(f)] == static_cast<std::int64_t>(p))
                ++n;
        order.mPatches[p].mStartFace = start;
        order.mPatches[p].mNFaces = n;
        start += n;
    }
    return order;
}

/**
 * @brief Check every clause of the polyMesh ordering contract.
 *
 * Runs in release builds too, deliberately: release is exactly where someone
 * writes a ten-million-cell case, the cost is a handful of flops per face
 * against ASCII formatting that costs far more, and a failure means we were
 * about to hand a solver a corrupt mesh.
 *
 * @return "" when valid, else the first violated clause, named.
 */
std::string foam_validate_order(const detail::GlobalFaces& rFaces, const FoamFaceOrder& rOrder,
                                const FoamPatchAssignment& rAssign, const NDArray& rPoints,
                                std::size_t PointDim, std::size_t NumPoints) {
    const std::size_t nf = rOrder.mNewToOld.size();
    if (nf != rFaces.NumFaces())
        return "C0: the written face list does not cover every face";

    // Cell centroids, for the two normal-direction clauses.
    std::vector<detail::Vec3> centroid(rFaces.NumCells(), detail::Vec3{0, 0, 0});
    for (std::size_t c = 0; c < rFaces.NumCells(); ++c) {
        detail::Vec3 acc{0, 0, 0};
        std::size_t n = 0;
        for (std::size_t k = 0; k < rFaces.NumCellFaces(c); ++k) {
            const std::size_t f = static_cast<std::size_t>(std::abs(rFaces.CellFaces(c)[k]) - 1);
            for (std::size_t j = 0; j < rFaces.FaceSize(f); ++j) {
                const detail::Vec3 p = detail::read_point(rPoints, PointDim, rFaces.Face(f)[j]);
                acc[0] += p[0];
                acc[1] += p[1];
                acc[2] += p[2];
                ++n;
            }
        }
        if (n)
            for (int a = 0; a < 3; ++a)
                acc[a] /= static_cast<double>(n);
        centroid[c] = acc;
    }

    std::vector<detail::Vec3> ring;
    for (std::size_t i = 0; i < nf; ++i) {
        const std::size_t f = static_cast<std::size_t>(rOrder.mNewToOld[i]);
        const bool is_internal = i < static_cast<std::size_t>(rOrder.mNumInternal);

        // C2: internal faces first.
        if (is_internal != (rFaces.mNeighbour[f] >= 0))
            return detail::format_compat("C2: face {} is {} but sits in the {} range", i,
                                         rFaces.mNeighbour[f] >= 0 ? "internal" : "boundary",
                                         is_internal ? "internal" : "boundary");

        // C7: node ids in range, ring big enough to bound an area.
        if (rFaces.FaceSize(f) < 3)
            return detail::format_compat("C7: face {} has fewer than three nodes", i);
        for (std::size_t k = 0; k < rFaces.FaceSize(f); ++k) {
            const std::int64_t id = rFaces.Face(f)[k];
            if (id < 0 || static_cast<std::size_t>(id) >= NumPoints)
                return detail::format_compat("C7: face {} references node {}, out of range", i, id);
        }

        ring.clear();
        for (std::size_t k = 0; k < rFaces.FaceSize(f); ++k)
            ring.push_back(detail::read_point(rPoints, PointDim, rFaces.Face(f)[k]));
        const detail::Vec3 nrm = detail::polygon_area_vector(ring.data(), ring.size());
        detail::Vec3 fc{0, 0, 0};
        for (const detail::Vec3& p : ring)
            for (int a = 0; a < 3; ++a)
                fc[a] += p[a] / static_cast<double>(ring.size());

        if (is_internal) {
            // C1: owner < neighbour.
            if (!(rFaces.mOwner[f] < rFaces.mNeighbour[f]))
                return detail::format_compat("C1: face {} has owner {} >= neighbour {}", i,
                                             rFaces.mOwner[f], rFaces.mNeighbour[f]);
            // C3: strictly increasing (owner, neighbour).
            if (i > 0) {
                const std::size_t g = static_cast<std::size_t>(rOrder.mNewToOld[i - 1]);
                const bool ok = rFaces.mOwner[g] < rFaces.mOwner[f] ||
                                (rFaces.mOwner[g] == rFaces.mOwner[f] &&
                                 rFaces.mNeighbour[g] <= rFaces.mNeighbour[f]);
                if (!ok)
                    return detail::format_compat(
                        "C3: face {} has (owner,neighbour)=({},{}) after ({},{})", i,
                        rFaces.mOwner[f], rFaces.mNeighbour[f], rFaces.mOwner[g],
                        rFaces.mNeighbour[g]);
            }
            // C4: normal points owner -> neighbour.
            const detail::Vec3& co = centroid[static_cast<std::size_t>(rFaces.mOwner[f])];
            const detail::Vec3& cn = centroid[static_cast<std::size_t>(rFaces.mNeighbour[f])];
            const double d =
                nrm[0] * (cn[0] - co[0]) + nrm[1] * (cn[1] - co[1]) + nrm[2] * (cn[2] - co[2]);
            if (!(d > 0.0))
                return detail::format_compat(
                    "C4: internal face {} does not point from owner to neighbour", i);
        } else {
            // C5: boundary normal points out of the domain.
            const detail::Vec3& co = centroid[static_cast<std::size_t>(rFaces.mOwner[f])];
            const double d =
                nrm[0] * (fc[0] - co[0]) + nrm[1] * (fc[1] - co[1]) + nrm[2] * (fc[2] - co[2]);
            if (!(d > 0.0))
                return detail::format_compat("C5: boundary face {} is wound inward", i);
        }
    }

    // C6: patches partition the boundary range exactly -- AND every face in a
    // patch's range really belongs to that patch. Checking only that the
    // start/count table tiles the range is not enough: the table is built by
    // counting, so it describes a contiguity the face order may simply not
    // have, and the resulting case is silently wrong.
    std::int64_t expect = rOrder.mNumInternal;
    for (std::size_t p = 0; p < rOrder.mPatches.size(); ++p) {
        const FoamPatchOut& patch = rOrder.mPatches[p];
        if (patch.mStartFace != expect)
            return detail::format_compat("C6: patch '{}' starts at {}, expected {}", patch.mName,
                                         patch.mStartFace, expect);
        for (std::int64_t i = patch.mStartFace; i < patch.mStartFace + patch.mNFaces; ++i) {
            const std::size_t f =
                static_cast<std::size_t>(rOrder.mNewToOld[static_cast<std::size_t>(i)]);
            if (rAssign.mFacePatch[f] != static_cast<std::int64_t>(p))
                return detail::format_compat(
                    "C6: face {} sits in patch '{}'s range but belongs to patch {}", i, patch.mName,
                    rAssign.mFacePatch[f]);
        }
        expect += patch.mNFaces;
    }
    if (expect != static_cast<std::int64_t>(nf))
        return detail::format_compat("C6: patches cover {} faces, expected {}",
                                     expect - rOrder.mNumInternal,
                                     static_cast<std::int64_t>(nf) - rOrder.mNumInternal);
    return "";
}

auto foam_open(const fs::path& rPath) {
    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("OpenFOAM: could not open for writing: " + rPath.string());
    return f;
}

/// One zone as it will be written: a name plus its member ids, already
/// converted to the OpenFOAM numbering (compact cell id / point id / written
/// face id -- see the three `foam_collect_*_zones` callers).
using FoamZoneOut = std::pair<std::string, std::vector<std::int64_t>>;

void foam_write_zone_file(const fs::path& rPath, const char* pClass, const char* pObject,
                          const char* pZoneType, const char* pLabelKey,
                          const std::vector<FoamZoneOut>& rZones,
                          const OpenFoamWriteOptions& rOpts) {
    auto f = foam_open(rPath);
    foam_write_header(f, pClass, pObject, rOpts);
    f << rZones.size() << "\n(\n";
    for (const auto& [name, ids] : rZones) {
        f << name << "\n{\n";
        f << "    type " << pZoneType << ";\n";
        f << "    " << pLabelKey << " List<label>\n    " << ids.size() << "\n    (";
        if (rOpts.mBinary) {
            for (std::int64_t id : ids)
                foam_write_binary_label(f, id, rOpts.mLabelBits);
            f << ");\n";
        } else {
            f << "\n";
            for (std::int64_t id : ids)
                f << "    " << id << "\n";
            f << "    );\n";
        }
        f << "}\n";
    }
    f << ")\n";
}

/// `Region`s of kind `Cell` -> `cellZones` entries: a global cell index maps
/// 1:1 onto a written OpenFOAM cell id via `rG2C` (`GlobalFaces::mCellToGlobal`
/// inverted) -- cells are never reordered on write, unlike faces.
std::vector<FoamZoneOut> foam_collect_cell_zones(
    const Mesh& rMesh, const std::unordered_map<std::int64_t, std::int64_t>& rG2C) {
    std::vector<FoamZoneOut> zones;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Cell)
            continue;
        std::vector<std::int64_t> ids;
        const std::int64_t* e = r.Entries();
        std::size_t n_dropped = 0;
        for (std::size_t k = 0; k < r.NumEntries(); ++k) {
            const auto it = rG2C.find(e[k]);
            if (it == rG2C.end()) {
                ++n_dropped;
                continue;
            }
            ids.push_back(it->second);
        }
        if (n_dropped > 0)
            log::warn("OpenFOAM: cellZone '{}' drops {} entr{} outside the volume mesh", r.mName,
                      n_dropped, n_dropped == 1 ? "y" : "ies");
        zones.emplace_back(r.mName, std::move(ids));
    }
    return zones;
}

/// `Region`s of kind `Point` -> `pointZones` entries: a point index needs no
/// conversion, since points are never reordered on write either.
std::vector<FoamZoneOut> foam_collect_point_zones(const Mesh& rMesh) {
    std::vector<FoamZoneOut> zones;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Point)
            continue;
        const std::int64_t* e = r.Entries();
        zones.emplace_back(r.mName, std::vector<std::int64_t>(e, e + r.NumEntries()));
    }
    return zones;
}

/// `Region`s of kind `Side` -> `faceZones` entries. A `(global cell, local
/// facet)` pair becomes a written face id via `rG2C` (global -> compact cell),
/// `GlobalFaces::CellFaces` (compact cell + local facet -> signed GlobalFaces
/// face id) and `rOldToNew` (`FoamFaceOrder::mNewToOld` inverted). `flipMap` is
/// never written -- see `foam_zone_label_list`'s doc comment on the read side.
std::vector<FoamZoneOut> foam_collect_face_zones(
    const Mesh& rMesh, const detail::GlobalFaces& rFaces,
    const std::unordered_map<std::int64_t, std::int64_t>& rG2C,
    const std::vector<std::int64_t>& rOldToNew) {
    std::vector<FoamZoneOut> zones;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Side)
            continue;
        std::vector<std::int64_t> ids;
        const std::int64_t* e = r.Entries();
        std::size_t n_dropped = 0;
        for (std::size_t k = 0; k < r.NumEntries(); ++k) {
            const std::int64_t global_cell = e[2 * k];
            const std::int64_t facet = e[2 * k + 1];
            const auto it = rG2C.find(global_cell);
            if (it == rG2C.end()) {
                ++n_dropped;
                continue;
            }
            const std::size_t compact = static_cast<std::size_t>(it->second);
            if (facet < 0 || static_cast<std::size_t>(facet) >= rFaces.NumCellFaces(compact)) {
                ++n_dropped;
                continue;
            }
            const std::int64_t signed_face = rFaces.CellFaces(compact)[facet];
            const std::size_t old_face = static_cast<std::size_t>(std::abs(signed_face) - 1);
            ids.push_back(rOldToNew[old_face]);
        }
        if (n_dropped > 0)
            log::warn("OpenFOAM: faceZone '{}' drops {} entr{} outside the volume mesh", r.mName,
                      n_dropped, n_dropped == 1 ? "y" : "ies");
        zones.emplace_back(r.mName, std::move(ids));
    }
    return zones;
}

/// One field to write into the time directory: a name, its OpenFOAM class, and
/// `mComponents` values per entry in OpenFOAM's own entry order.
struct FoamFieldOut {
    std::string mName;
    const char* mClass;
    const char* mListType;  // scalar, vector, symmTensor or tensor
    std::size_t mComponents;
    bool mIsPoint;
    std::vector<double> mValues;  // entries x components
    /// A boundary patch's face values, taken from the boundary-face cell blocks.
    struct PatchValues {
        std::size_t mPatch;           ///< index into the written patch table
        std::vector<double> mValues;  ///< the patch's faces in written order x components
    };
    std::vector<PatchValues> mPatches;  ///< only patches whose every face carries a finite value
};

const char* foam_field_class_name(std::size_t Components, bool IsPoint, const char** pListType) {
    switch (Components) {
        case 1:
            *pListType = "scalar";
            return IsPoint ? "pointScalarField" : "volScalarField";
        case 3:
            *pListType = "vector";
            return IsPoint ? "pointVectorField" : "volVectorField";
        case 6:
            *pListType = "symmTensor";
            return IsPoint ? nullptr : "volSymmTensorField";
        case 9:
            *pListType = "tensor";
            return IsPoint ? nullptr : "volTensorField";
        default:
            return nullptr;
    }
}

/// Patch types whose field patches are fixed by OpenFOAM (`empty`, `wedge`, `symmetry`,
/// `symmetryPlane`): a `calculated` entry with values there is refused, so they get none.
bool foam_is_constraint_patch(const std::string& rType) {
    return rType == "empty" || rType == "wedge" || rType == "symmetry" || rType == "symmetryPlane";
}

/// `point_data` and `cell_data` as OpenFOAM fields, the inverse of the reader: volume cell data
/// in OpenFOAM's compact cell order, point data in point order. The values of the
/// boundary-face cell blocks become per-patch `boundaryField` values, for every patch whose
/// faces all carry a finite value (so the NaN the reader fills a field without patch values
/// with writes nothing) and that is not a constraint patch. Arrays with no OpenFOAM field
/// class (component counts other than 1/3/6/9, vectors-only on points) are skipped with a
/// warning.
std::vector<FoamFieldOut> foam_collect_fields(const Mesh& rMesh, const detail::GlobalFaces& rFaces,
                                              const FoamFaceOrder& rOrder,
                                              const FoamPatchAssignment& rAssign) {
    std::vector<FoamFieldOut> out;
    auto usable_name = [](const std::string& rName) {
        return !rName.empty() && rName != "cell_tags" &&
               rName.find_first_of("/\\") == std::string::npos && rName[0] != '.';
    };
    auto components_of = [](const NDArray& rArr) {
        return rArr.Shape().size() >= 2 ? detail::cols(rArr) : std::size_t{1};
    };
    for (const std::string& name : rMesh.PointDataNames()) {
        if (!usable_name(name))
            continue;
        const NDArray& arr = rMesh.PointData(name);
        const std::size_t nc = components_of(arr);
        FoamFieldOut f{name, nullptr, nullptr, nc, true, {}};
        f.mClass = foam_field_class_name(nc, true, &f.mListType);
        if (f.mClass == nullptr) {
            log::warn("OpenFOAM: point data '{}' has {} component(s); OpenFOAM point fields hold "
                      "1 or 3. Not written",
                      name, nc);
            continue;
        }
        const detail::DoubleView values(arr);
        f.mValues.assign(values.Data(), values.Data() + arr.Size());
        out.push_back(std::move(f));
    }
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    for (const std::string& name : rMesh.CellDataNames()) {
        if (!usable_name(name))
            continue;
        if (rMesh.CellDataNumBlocks(name) != rMesh.NumCellBlocks()) {
            log::warn("OpenFOAM: cell data '{}' is not present on every cell block. Not written",
                      name);
            continue;
        }
        const std::size_t nc = components_of(rMesh.CellData(name, 0));
        FoamFieldOut f{name, nullptr, nullptr, nc, false, {}};
        f.mClass = foam_field_class_name(nc, false, &f.mListType);
        if (f.mClass == nullptr) {
            log::warn("OpenFOAM: cell data '{}' has {} component(s); OpenFOAM cell fields hold "
                      "1, 3, 6 or 9. Not written",
                      name, nc);
            continue;
        }
        std::deque<detail::DoubleView> views;  // views are neither copyable nor movable
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            views.emplace_back(rMesh.CellData(name, b));
        f.mValues.resize(rFaces.NumCells() * nc);
        for (std::size_t cell = 0; cell < rFaces.NumCells(); ++cell) {
            const auto [block, row] = detail::global_to_block_row(bases, rFaces.mCellToGlobal[cell]);
            for (std::size_t c = 0; c < nc; ++c)
                f.mValues[cell * nc + c] = views[block][row * nc + c];
        }
        for (std::size_t pi = 0; pi < rOrder.mPatches.size(); ++pi) {
            const FoamPatchOut& patch = rOrder.mPatches[pi];
            if (patch.mNFaces <= 0 || foam_is_constraint_patch(patch.mType))
                continue;
            FoamFieldOut::PatchValues pv{pi, {}};
            pv.mValues.reserve(static_cast<std::size_t>(patch.mNFaces) * nc);
            bool complete = true;
            for (std::int64_t k = 0; k < patch.mNFaces && complete; ++k) {
                const auto gf = static_cast<std::size_t>(
                    rOrder.mNewToOld[static_cast<std::size_t>(patch.mStartFace + k)]);
                const std::int64_t source =
                    gf < rAssign.mFaceSource.size() ? rAssign.mFaceSource[gf] : -1;
                if (source < 0) {
                    complete = false;
                    break;
                }
                const auto [block, row] = detail::global_to_block_row(bases, source);
                for (std::size_t c = 0; c < nc; ++c) {
                    const double v = views[block][row * nc + c];
                    if (!std::isfinite(v))
                        complete = false;
                    pv.mValues.push_back(v);
                }
            }
            if (complete)
                f.mPatches.push_back(std::move(pv));
        }
        out.push_back(std::move(f));
    }
    return out;
}

/// Write fields into `<rDir>` (a time directory); `rLocation` is the header's `location`.
/// A `nonuniform List<T>` value, from the type keyword to the closing `)` plus @p pClose: the
/// count on a line of its own, then the entries (binary: packed scalars, ASCII: one per line).
void foam_write_value_list(std::ostream& rOs, const FoamFieldOut& rField,
                           const std::vector<double>& rValues, const OpenFoamWriteOptions& rOpts,
                           const char* pClose) {
    const std::size_t entries = rValues.size() / rField.mComponents;
    rOs << "nonuniform List<" << rField.mListType << ">\n" << entries << "\n(";
    if (rOpts.mBinary) {
        for (const double v : rValues)
            foam_write_binary_scalar(rOs, v, rOpts.mScalarBits);
        rOs << ")" << pClose;
        return;
    }
    rOs << "\n";
    const detail::CNumber num;
    std::string buf;
    for (std::size_t i = 0; i < entries; ++i) {
        if (rField.mComponents > 1)
            buf += '(';
        for (std::size_t c = 0; c < rField.mComponents; ++c) {
            if (c)
                buf += ' ';
            num.Append(buf, "%.16g", rValues[i * rField.mComponents + c]);
        }
        buf += rField.mComponents > 1 ? ")\n" : "\n";
    }
    rOs << buf << ")" << pClose;
}

void foam_emit_fields(const fs::path& rDir, const std::string& rLocation,
                      const std::vector<FoamFieldOut>& rFields,
                      const std::vector<FoamPatchOut>& rPatches,
                      const OpenFoamWriteOptions& rOpts) {
    if (rFields.empty())
        return;
    const std::vector<FoamFieldOut>& out = rFields;
    const fs::path dir = rDir;
    const std::string& location = rLocation;
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec && !fs::is_directory(dir))
        throw WriteError("OpenFOAM: could not create directory " + dir.string() + ": " +
                         ec.message());
    for (const FoamFieldOut& field : out) {
        auto f = foam_open(dir / field.mName);
        foam_write_header(f, field.mClass, field.mName, rOpts, location.c_str());
        f << "dimensions      [0 0 0 0 0 0 0];\n\n";
        f << "internalField   ";
        foam_write_value_list(f, field, field.mValues, rOpts, "\n;\n");
        // Patches whose boundary-face blocks carry values get them as explicit entries; a
        // catch-all patch value keeps the case loadable for every other patch, and a solver
        // replaces these with its real conditions.
        f << "\nboundaryField\n{\n";
        for (const FoamFieldOut::PatchValues& pv : field.mPatches) {
            f << "    " << rPatches[pv.mPatch].mName
              << "\n    {\n        type            calculated;\n        value           ";
            foam_write_value_list(f, field, pv.mValues, rOpts, "\n        ;\n");
            f << "    }\n";
        }
        f << "    \".*\"\n    {\n        type            calculated;\n"
             "        value           uniform ";
        if (field.mComponents == 1) {
            f << "0";
        } else {
            f << "(0";
            for (std::size_t c = 1; c < field.mComponents; ++c)
                f << " 0";
            f << ")";
        }
        f << ";\n    }\n}\n";
    }
}

void foam_write_fields(const fs::path& rCaseRoot, const Mesh& rMesh,
                       const detail::GlobalFaces& rFaces, const FoamFaceOrder& rOrder,
                       const FoamPatchAssignment& rAssign, const OpenFoamWriteOptions& rOpts,
                       const std::string& rTime = "0", const std::string& rRegion = {}) {
    const std::vector<FoamFieldOut> out = foam_collect_fields(rMesh, rFaces, rOrder, rAssign);
    const fs::path dir = rRegion.empty() ? rCaseRoot / rTime : rCaseRoot / rTime / rRegion;
    foam_emit_fields(dir, rRegion.empty() ? rTime : rTime + "/" + rRegion, out, rOrder.mPatches,
                     rOpts);
}

}  // namespace

namespace {

/// What a write touches: the whole case, or only the fields of one time directory (the later
/// steps of a series, over a mesh the first step already wrote).
enum class FoamWriteMode { Full, FieldsOnly, Decomposed };

/// Where a write goes: an empty region is the single-region layout.
struct FoamTarget {
    std::string mRegion;
    std::string mTime = "0";
    const std::vector<std::int64_t>* mParts = nullptr;  ///< Decomposed: a processor per compact cell
};

fs::path foam_region_polymesh_dir(const fs::path& rPath, const std::string& rRegion) {
    if (rRegion.empty())
        return foam_polymesh_dir(rPath, /*ForWrite=*/true);
    if (rRegion.find_first_of("/\\") != std::string::npos || rRegion[0] == '.')
        throw WriteError("OpenFOAM: invalid region name '" + rRegion + "'");
    if (rPath.filename() == "polyMesh")
        throw WriteError(
            "OpenFOAM: a region needs a case path (a .foam file or a case directory), not a "
            "polyMesh directory");
    const fs::path root = rPath.extension() == ".foam" ? rPath.parent_path() : rPath;
    return root / "constant" / rRegion / "polyMesh";
}

/// The case directory above `constant/[<region>/]polyMesh`.
fs::path foam_case_root(const fs::path& rPoly, const std::string& rRegion) {
    fs::path root = rPoly.parent_path().parent_path();
    return rRegion.empty() ? root : root.parent_path();
}


/// One face of a processor's mesh: the written (serial) face it is, and whether the processor's
/// copy is the reversed one (its cell is the serial face's neighbour).
struct FoamLocalFace {
    std::int64_t mSerial = 0;
    bool mFlipped = false;
};

/// The zones of a mesh in the serial case's numbering (compact cell ids, point ids, written face
/// ids), which a decomposed write splits per processor.
struct FoamZonesOut {
    std::vector<FoamZoneOut> mCell, mPoint, mFace;
};

void foam_emit_label_list(const fs::path& rFile, const char* pObject,
                          const std::vector<std::int64_t>& rValues,
                          const OpenFoamWriteOptions& rOpts) {
    auto f = foam_open(rFile);
    foam_write_header(f, "labelList", pObject, rOpts);
    if (rOpts.mBinary) {
        f << rValues.size() << "\n(";
        for (const std::int64_t v : rValues)
            foam_write_binary_label(f, v, rOpts.mLabelBits);
        f << ")\n";
    } else {
        f << rValues.size() << "\n(\n";
        for (const std::int64_t v : rValues)
            f << v << "\n";
        f << ")\n";
    }
}

/// A decomposed case, `decomposePar`'s layout: for each processor `processorN/constant/polyMesh`
/// holding its cells' faces (internal faces, then the original patches, then a `processor` patch
/// per neighbour) plus the four addressing lists the reader reassembles the global case from, and
/// its share of the fields in `processorN/<time>/`.
void foam_write_decomposed(const std::string& rPath, const Mesh& rMesh,
                           const detail::GlobalFaces& rFaces, const FoamFaceOrder& rOrder,
                           const FoamPatchAssignment& rAssign, const OpenFoamWriteOptions& rOpts,
                           const std::vector<std::int64_t>& rParts, const std::string& rTime,
                           const FoamZonesOut& rZones) {
    const fs::path given(rPath);
    const fs::path root = given.extension() == ".foam" ? given.parent_path() : given;
    if (fs::exists(root / "constant" / "polyMesh"))
        log::warn(
            "OpenFOAM: {} already holds a serial constant/polyMesh, which a reader prefers to the "
            "processor directories; write a decomposed case into a fresh directory",
            root.string());
    std::error_code ec;
    fs::create_directories(root, ec);
    if (ec && !fs::is_directory(root))
        throw WriteError("OpenFOAM: could not create directory " + root.string() + ": " +
                         ec.message());
    if (given.extension() == ".foam")
        auto marker = detail::make_classic_ofstream(rPath, std::ios::binary);

    const std::size_t ncells = rFaces.NumCells();
    std::set<std::int64_t> proc_set(rParts.begin(), rParts.end());
    if (!proc_set.empty() && *proc_set.begin() < 0)
        throw WriteError("OpenFOAM: partition:part holds a negative processor number");
    const std::size_t nfaces = rOrder.mNewToOld.size();
    std::vector<std::int64_t> face_patch(nfaces, -1);
    for (std::size_t pi = 0; pi < rOrder.mPatches.size(); ++pi)
        for (std::int64_t k = 0; k < rOrder.mPatches[pi].mNFaces; ++k)
            face_patch[static_cast<std::size_t>(rOrder.mPatches[pi].mStartFace + k)] =
                static_cast<std::int64_t>(pi);
    const std::vector<FoamFieldOut> fields = foam_collect_fields(rMesh, rFaces, rOrder, rAssign);
    const NDArray& pts = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();

    for (const std::int64_t proc : proc_set) {
        const fs::path pdir = root / ("processor" + std::to_string(proc));
        const fs::path poly = pdir / "constant" / "polyMesh";
        fs::create_directories(poly, ec);
        if (ec && !fs::is_directory(poly))
            throw WriteError("OpenFOAM: could not create directory " + poly.string() + ": " +
                             ec.message());
        std::vector<std::int64_t> cell_addr, cell_local(ncells, -1);
        for (std::size_t c = 0; c < ncells; ++c)
            if (rParts[c] == proc) {
                cell_local[c] = static_cast<std::int64_t>(cell_addr.size());
                cell_addr.push_back(static_cast<std::int64_t>(c));
            }

        std::vector<FoamLocalFace> internal;
        std::vector<std::vector<FoamLocalFace>> by_patch(rOrder.mPatches.size());
        std::map<std::int64_t, std::vector<FoamLocalFace>> by_proc;
        for (std::size_t sf = 0; sf < nfaces; ++sf) {
            const auto gf = static_cast<std::size_t>(rOrder.mNewToOld[sf]);
            const std::int64_t o = rFaces.mOwner[gf];
            const std::int64_t n = rFaces.mNeighbour[gf];
            const auto serial = static_cast<std::int64_t>(sf);
            if (static_cast<std::int64_t>(sf) < rOrder.mNumInternal) {
                const std::int64_t po = rParts[static_cast<std::size_t>(o)];
                const std::int64_t pn = rParts[static_cast<std::size_t>(n)];
                if (po == proc && pn == proc)
                    internal.push_back({serial, false});
                else if (po == proc)
                    by_proc[pn].push_back({serial, false});
                else if (pn == proc)
                    by_proc[po].push_back({serial, true});
            } else if (rParts[static_cast<std::size_t>(o)] == proc && face_patch[sf] >= 0) {
                by_patch[static_cast<std::size_t>(face_patch[sf])].push_back({serial, false});
            }
        }

        // the processor's faces in written order
        std::vector<FoamLocalFace> all = internal;
        std::vector<FoamPatchOut> patches;
        std::vector<std::int64_t> boundary_addr;
        for (std::size_t pi = 0; pi < rOrder.mPatches.size(); ++pi) {
            FoamPatchOut p = rOrder.mPatches[pi];
            p.mStartFace = static_cast<std::int64_t>(all.size());
            p.mNFaces = static_cast<std::int64_t>(by_patch[pi].size());
            p.mExtra.clear();
            all.insert(all.end(), by_patch[pi].begin(), by_patch[pi].end());
            patches.push_back(std::move(p));
            boundary_addr.push_back(static_cast<std::int64_t>(pi));
        }
        for (const auto& [other, list] : by_proc) {
            FoamPatchOut p;
            p.mName = "procBoundary" + std::to_string(proc) + "to" + std::to_string(other);
            p.mType = "processor";
            p.mStartFace = static_cast<std::int64_t>(all.size());
            p.mNFaces = static_cast<std::int64_t>(list.size());
            p.mExtra =
                "        inGroups        List<word> 1(processor);\n"
                "        matchTolerance  0.0001;\n"
                "        transform       unknown;\n"
                "        myProcNo        " +
                std::to_string(proc) + ";\n        neighbProcNo    " + std::to_string(other) + ";\n";
            all.insert(all.end(), list.begin(), list.end());
            patches.push_back(std::move(p));
            boundary_addr.push_back(-1);
        }

        // local points: the nodes of those faces, in global order
        std::vector<std::int64_t> point_local(rMesh.NumPoints(), -1);
        for (const FoamLocalFace& lf : all) {
            const auto gf = static_cast<std::size_t>(rOrder.mNewToOld[static_cast<std::size_t>(lf.mSerial)]);
            for (std::size_t k = 0; k < rFaces.FaceSize(gf); ++k)
                point_local[static_cast<std::size_t>(rFaces.Face(gf)[k])] = 0;
        }
        std::vector<std::int64_t> point_addr;
        for (std::size_t i = 0; i < point_local.size(); ++i)
            if (point_local[i] == 0) {
                point_local[i] = static_cast<std::int64_t>(point_addr.size());
                point_addr.push_back(static_cast<std::int64_t>(i));
            }

        {
            auto f = foam_open(poly / "points");
            foam_write_header(f, "vectorField", "points", rOpts);
            if (rOpts.mBinary) {
                f << point_addr.size() << "\n(";
                for (const std::int64_t gp : point_addr)
                    for (const double c : detail::read_point(pts, dim, gp))
                        foam_write_binary_scalar(f, c, rOpts.mScalarBits);
                f << ")\n";
            } else {
                f << point_addr.size() << "\n(\n";
                const detail::CNumber num;
                std::string buf;
                for (const std::int64_t gp : point_addr) {
                    const detail::Vec3 p = detail::read_point(pts, dim, gp);
                    buf += '(';
                    num.Append(buf, "%.16g", p[0]);
                    buf += ' ';
                    num.Append(buf, "%.16g", p[1]);
                    buf += ' ';
                    num.Append(buf, "%.16g", p[2]);
                    buf += ")\n";
                }
                f << buf << ")\n";
            }
        }
        std::vector<std::int64_t> owner, neighbour, face_addr;
        {
            auto f = foam_open(poly / "faces");
            foam_write_header(f, "faceList", "faces", rOpts);
            f << all.size() << (rOpts.mBinary ? "\n(" : "\n(\n");
            for (std::size_t li = 0; li < all.size(); ++li) {
                const FoamLocalFace& lf = all[li];
                const auto gf =
                    static_cast<std::size_t>(rOrder.mNewToOld[static_cast<std::size_t>(lf.mSerial)]);
                const std::size_t size = rFaces.FaceSize(gf);
                std::vector<std::int64_t> ring(size);
                for (std::size_t k = 0; k < size; ++k) {
                    // OpenFOAM's reverseFace: the first point stays, the rest run backwards
                    const std::size_t src = lf.mFlipped ? (k == 0 ? 0 : size - k) : k;
                    ring[k] = point_local[static_cast<std::size_t>(rFaces.Face(gf)[src])];
                }
                if (rOpts.mBinary) {
                    f << size << "(";
                    for (const std::int64_t v : ring)
                        foam_write_binary_label(f, v, rOpts.mLabelBits);
                    f << ")";
                } else {
                    f << size << "(";
                    for (std::size_t k = 0; k < size; ++k)
                        f << (k ? " " : "") << ring[k];
                    f << ")\n";
                }
                const std::int64_t o = rFaces.mOwner[gf];
                const std::int64_t n = rFaces.mNeighbour[gf];
                owner.push_back(cell_local[static_cast<std::size_t>(lf.mFlipped ? n : o)]);
                if (static_cast<std::int64_t>(li) < static_cast<std::int64_t>(internal.size()))
                    neighbour.push_back(cell_local[static_cast<std::size_t>(n)]);
                face_addr.push_back(lf.mFlipped ? -(lf.mSerial + 1) : lf.mSerial + 1);
            }
            f << (rOpts.mBinary ? ")\n" : ")\n");
        }
        foam_emit_label_list(poly / "owner", "owner", owner, rOpts);
        foam_emit_label_list(poly / "neighbour", "neighbour", neighbour, rOpts);
        {
            auto f = foam_open(poly / "boundary");
            foam_write_header(f, "polyBoundaryMesh", "boundary", rOpts);
            f << patches.size() << "\n(\n";
            for (const FoamPatchOut& p : patches) {
                f << "    " << p.mName << "\n    {\n";
                f << "        type            " << p.mType << ";\n";
                f << p.mExtra;
                f << "        nFaces          " << p.mNFaces << ";\n";
                f << "        startFace       " << p.mStartFace << ";\n";
                f << "    }\n";
            }
            f << ")\n";
        }
        foam_emit_label_list(poly / "pointProcAddressing", "pointProcAddressing", point_addr, rOpts);
        foam_emit_label_list(poly / "cellProcAddressing", "cellProcAddressing", cell_addr, rOpts);
        foam_emit_label_list(poly / "faceProcAddressing", "faceProcAddressing", face_addr, rOpts);
        foam_emit_label_list(poly / "boundaryProcAddressing", "boundaryProcAddressing",
                             boundary_addr, rOpts);

        // This processor's zones, as decomposePar writes them: every zone on every processor
        // (an empty list where it has no member), with processor-local ids. A face on a
        // processor boundary is in both processors' zones.
        {
            std::vector<std::int64_t> face_local(nfaces, -1);
            for (std::size_t li = 0; li < all.size(); ++li)
                face_local[static_cast<std::size_t>(all[li].mSerial)] =
                    static_cast<std::int64_t>(li);
            auto localize = [&](const std::vector<FoamZoneOut>& rGlobal,
                                const std::vector<std::int64_t>& rMap) {
                std::vector<FoamZoneOut> local;
                for (const auto& [name, ids] : rGlobal) {
                    FoamZoneOut z{name, {}};
                    for (const std::int64_t g : ids) {
                        const std::int64_t l = rMap[static_cast<std::size_t>(g)];
                        if (l >= 0)
                            z.second.push_back(l);
                    }
                    local.push_back(std::move(z));
                }
                return local;
            };
            auto emit = [&](const char* pFile, const std::vector<FoamZoneOut>& rLocal,
                            const char* pClass, const char* pType, const char* pKey) {
                if (rLocal.empty()) {
                    std::error_code rc;
                    fs::remove(poly / pFile, rc);
                    return;
                }
                foam_write_zone_file(poly / pFile, pClass, pFile, pType, pKey, rLocal, rOpts);
            };
            emit("cellZones", localize(rZones.mCell, cell_local), "cellZoneList", "cellZone",
                 "cellLabels");
            emit("pointZones", localize(rZones.mPoint, point_local), "pointZoneList", "pointZone",
                 "pointLabels");
            emit("faceZones", localize(rZones.mFace, face_local), "faceZoneList", "faceZone",
                 "faceLabels");
        }

        // this processor's share of the fields
        std::vector<FoamFieldOut> local;
        for (const FoamFieldOut& field : fields) {
            FoamFieldOut part{field.mName, field.mClass, field.mListType, field.mComponents,
                              field.mIsPoint, {}};
            const std::vector<std::int64_t>& addr = field.mIsPoint ? point_addr : cell_addr;
            for (const std::int64_t g : addr)
                for (std::size_t c = 0; c < field.mComponents; ++c)
                    part.mValues.push_back(
                        field.mValues[static_cast<std::size_t>(g) * field.mComponents + c]);
            // A processor keeps the values of the faces of each original patch it holds; the
            // `processor` patches that follow them (index >= the original count) get none.
            for (const FoamFieldOut::PatchValues& pv : field.mPatches) {
                if (by_patch[pv.mPatch].empty())
                    continue;
                FoamFieldOut::PatchValues lpv{pv.mPatch, {}};
                for (const FoamLocalFace& lf : by_patch[pv.mPatch]) {
                    const auto at = static_cast<std::size_t>(lf.mSerial -
                                                             rOrder.mPatches[pv.mPatch].mStartFace);
                    for (std::size_t c = 0; c < field.mComponents; ++c)
                        lpv.mValues.push_back(pv.mValues[at * field.mComponents + c]);
                }
                part.mPatches.push_back(std::move(lpv));
            }
            local.push_back(std::move(part));
        }
        foam_emit_fields(pdir / rTime, rTime, local, patches, rOpts);
    }
    log::info("Wrote {} processor directories under {}", proc_set.size(), root.string());
}

void foam_write_impl(const std::string& rPath, const Mesh& rMesh, const OpenFoamInfo& rInfo,
                     const OpenFoamWriteOptions& rOpts, const FoamTarget& rTarget,
                     FoamWriteMode Mode) {
    if (rOpts.mBinary && std::endian::native == std::endian::big)
        throw WriteError(
            "OpenFOAM: binary write requested on a big-endian host; OpenFOAM binary files are "
            "little-endian only, and writing big-endian bytes under a 'LSB' arch header would "
            "silently corrupt every value read back");
    if (rOpts.mBinary && rOpts.mLabelBits != 32 && rOpts.mLabelBits != 64)
        throw WriteError("OpenFOAM: label_bits must be 32 or 64");
    if (rOpts.mBinary && rOpts.mScalarBits != 32 && rOpts.mScalarBits != 64)
        throw WriteError("OpenFOAM: scalar_bits must be 32 or 64");

    const detail::GlobalFaces faces = detail::build_global_faces(rMesh);

    if (Mode == FoamWriteMode::FieldsOnly) {
        // The mesh is already on disk; its patch table is rebuilt (quietly: the first step
        // reported any dropped cell) only to know which patch each boundary face went to.
        FoamPatchAssignment assign;
        FoamFaceOrder order;
        if (!faces.mNonCellBlocks.empty() && faces.NumCells() > 0) {
            assign = foam_assign_patches(rMesh, faces, rInfo);
            order = foam_order_faces(faces, assign);
        }
        const fs::path only = foam_region_polymesh_dir(fs::path(rPath), rTarget.mRegion);
        foam_write_fields(foam_case_root(only, rTarget.mRegion), rMesh, faces, order, assign, rOpts,
                          rTarget.mTime, rTarget.mRegion);
        return;
    }

    if (faces.NumCells() == 0)
        throw WriteError(
            "OpenFOAM: the mesh has no volume cells; a polyMesh needs at least one "
            "tetra/pyramid/wedge/hexahedron/polyhedron cell");
    // A 3D block we could not turn into faces would be a silently dropped solid.
    for (std::size_t block : faces.mNonCellBlocks) {
        const auto cb = rMesh.Cells(block);
        const std::string type(cb.Type());
        if (cell_type_dimension(cell_type_from_name(type)) == 3)
            throw WriteError(detail::format_compat(
                "OpenFOAM: cell type '{}' is 3D but has no face topology in meshio++, so writing "
                "it would silently drop those cells",
                type));
    }
    if (faces.mNumNonManifold > 0)
        throw WriteError(detail::format_compat(
            "OpenFOAM: {} face(s) are shared by three or more cells; a polyMesh face has at most "
            "an owner and one neighbour",
            faces.mNumNonManifold));
    if (faces.mNumUnorientable > 0)
        log::warn(
            "OpenFOAM: {} cell(s) are not closed orientable surfaces; their faces are "
            "written with the winding they arrived with",
            faces.mNumUnorientable);
    if (faces.mNumFlipped > 0)
        log::info("OpenFOAM: rewound {} inverted cell(s) so their faces point outward",
                  faces.mNumFlipped);

    const FoamPatchAssignment assign = foam_assign_patches(rMesh, faces, rInfo);
    if (assign.mNumOrphan > 0) {
        log::warn("OpenFOAM: {} boundary cell(s) match no cell face and are not written",
                  assign.mNumOrphan);
        detail::provenance_note("cells-dropped",
                                std::to_string(assign.mNumOrphan) +
                                    " boundary cell(s) dropped -- they match no cell face");
    }
    if (assign.mNumInternalTagged > 0) {
        log::warn(
            "OpenFOAM: {} boundary cell(s) coincide with an INTERNAL face; OpenFOAM cannot "
            "put such a face on a patch, so they are not written",
            assign.mNumInternalTagged);
        detail::provenance_note("cells-dropped",
                                std::to_string(assign.mNumInternalTagged) +
                                    " boundary cell(s) dropped -- they coincide with an internal "
                                    "face, which OpenFOAM cannot put on a patch");
    }

    const FoamFaceOrder order = foam_order_faces(faces, assign);
    const std::string bad = foam_validate_order(faces, order, assign, rMesh.Points(),
                                                rMesh.PointDim(), rMesh.NumPoints());
    if (!bad.empty())
        throw WriteError("OpenFOAM: internal error, the written face order violates " + bad);

    if (Mode == FoamWriteMode::Decomposed) {
        std::unordered_map<std::int64_t, std::int64_t> global_to_compact;
        for (std::size_t c = 0; c < faces.mCellToGlobal.size(); ++c)
            global_to_compact[faces.mCellToGlobal[c]] = static_cast<std::int64_t>(c);
        std::vector<std::int64_t> old_to_new(faces.NumFaces());
        for (std::size_t i = 0; i < order.mNewToOld.size(); ++i)
            old_to_new[static_cast<std::size_t>(order.mNewToOld[i])] = static_cast<std::int64_t>(i);
        FoamZonesOut zones;
        zones.mCell = foam_collect_cell_zones(rMesh, global_to_compact);
        zones.mPoint = foam_collect_point_zones(rMesh);
        zones.mFace = foam_collect_face_zones(rMesh, faces, global_to_compact, old_to_new);
        foam_write_decomposed(rPath, rMesh, faces, order, assign, rOpts, *rTarget.mParts,
                              rTarget.mTime, zones);
        return;
    }

    const fs::path poly = foam_region_polymesh_dir(fs::path(rPath), rTarget.mRegion);
    std::error_code ec;
    fs::create_directories(poly, ec);
    if (ec && !fs::is_directory(poly))
        throw WriteError("OpenFOAM: could not create directory " + poly.string() + ": " +
                         ec.message());
    if (fs::path(rPath).extension() == ".foam") {
        // The marker file is what makes the case openable by ParaView and by
        // this reader's own `.foam` branch.
        auto marker = detail::make_classic_ofstream(rPath, std::ios::binary);
    }

    // Companion files this writer does not produce but OpenFOAM would read.
    // Leaving a stale one behind corrupts the case, so remove exactly these --
    // never the whole directory, which may hold a user's own files.
    // `cellZones`/`faceZones`/`pointZones` are handled separately below: this
    // writer produces them once the mesh carries the matching `Region` kind,
    // and only deletes a stale one when it no longer does (so an old zone
    // file is not left behind once its region is removed from the mesh).
    for (const char* name : {"meshModifiers", "boundaryProcAddressing", "cellProcAddressing",
                             "faceProcAddressing", "pointProcAddressing", "cellLevel", "pointLevel",
                             "level0Edge", "refinementHistory", "surfaceIndex"}) {
        std::error_code rc;
        if (fs::remove(poly / name, rc))
            log::info("OpenFOAM: removed stale {}", name);
    }

    const NDArray& pts = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t np = rMesh.NumPoints();

    {
        auto f = foam_open(poly / "points");
        foam_write_header(f, "vectorField", "points", rOpts);
        // The count MUST be on a line of its own: every ASCII parser here takes
        // "the first line that is entirely digits" as the count, so `8(` would
        // be read as data and the list would come back EMPTY, not as an error.
        if (rOpts.mBinary) {
            f << np << "\n(";
            for (std::size_t i = 0; i < np; ++i) {
                const detail::Vec3 p = detail::read_point(pts, dim, static_cast<std::int64_t>(i));
                for (double c : p)
                    foam_write_binary_scalar(f, c, rOpts.mScalarBits);
            }
            f << ")\n";
        } else {
            f << np << "\n(\n";
            // `f << std::setprecision(16) << x` on this classic-locale stream,
            // with no float field set, prints exactly printf's "%.16g"; rows are
            // formatted that way in parallel chunks (row_writer.hpp).
            const detail::CNumber num;
            detail::write_row_chunks(
                f, np, [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                    for (std::size_t i = First; i < Last; ++i) {
                        const detail::Vec3 p =
                            detail::read_point(pts, dim, static_cast<std::int64_t>(i));
                        rBuf += '(';
                        num.Append(rBuf, "%.16g", p[0]);
                        rBuf += ' ';
                        num.Append(rBuf, "%.16g", p[1]);
                        rBuf += ' ';
                        num.Append(rBuf, "%.16g", p[2]);
                        rBuf += ")\n";
                    }
                });
            f << ")\n";
        }
    }
    {
        auto f = foam_open(poly / "faces");
        foam_write_header(f, "faceList", "faces", rOpts);
        // Non-contiguous, so each face is its own length-prefixed labelList,
        // in binary exactly as in ASCII -- the same shape read_binary_faces
        // (this file's own reader) expects, not CompactListList.
        if (rOpts.mBinary) {
            f << order.mNewToOld.size() << "\n(";
            for (std::int64_t old : order.mNewToOld) {
                const std::size_t fi = static_cast<std::size_t>(old);
                f << faces.FaceSize(fi) << "(";
                for (std::size_t k = 0; k < faces.FaceSize(fi); ++k)
                    foam_write_binary_label(f, faces.Face(fi)[k], rOpts.mLabelBits);
                f << ")";
            }
            f << ")\n";
        } else {
            f << order.mNewToOld.size() << "\n(\n";
            for (std::int64_t old : order.mNewToOld) {
                const std::size_t fi = static_cast<std::size_t>(old);
                f << faces.FaceSize(fi) << "(";
                for (std::size_t k = 0; k < faces.FaceSize(fi); ++k)
                    f << (k ? " " : "") << faces.Face(fi)[k];
                f << ")\n";
            }
            f << ")\n";
        }
    }
    {
        auto f = foam_open(poly / "owner");
        foam_write_header(f, "labelList", "owner", rOpts);
        if (rOpts.mBinary) {
            f << order.mNewToOld.size() << "\n(";
            for (std::int64_t old : order.mNewToOld)
                foam_write_binary_label(f, faces.mOwner[static_cast<std::size_t>(old)],
                                        rOpts.mLabelBits);
            f << ")\n";
        } else {
            f << order.mNewToOld.size() << "\n(\n";
            for (std::int64_t old : order.mNewToOld)
                f << faces.mOwner[static_cast<std::size_t>(old)] << "\n";
            f << ")\n";
        }
    }
    {
        // Always written, even with zero entries: a stale `neighbour` left from
        // a previous, larger case is one of the nastiest ways to corrupt one.
        // OpenFOAM's `neighbour` holds ONLY internal faces -- our own reader
        // also accepts a -1-padded full-length list, which is exactly why a
        // round trip through it is a weak oracle for this writer.
        auto f = foam_open(poly / "neighbour");
        foam_write_header(f, "labelList", "neighbour", rOpts);
        if (rOpts.mBinary) {
            f << order.mNumInternal << "\n(";
            for (std::int64_t i = 0; i < order.mNumInternal; ++i)
                foam_write_binary_label(f,
                                        faces.mNeighbour[static_cast<std::size_t>(
                                            order.mNewToOld[static_cast<std::size_t>(i)])],
                                        rOpts.mLabelBits);
            f << ")\n";
        } else {
            f << order.mNumInternal << "\n(\n";
            for (std::int64_t i = 0; i < order.mNumInternal; ++i)
                f << faces.mNeighbour[static_cast<std::size_t>(
                         order.mNewToOld[static_cast<std::size_t>(i)])]
                  << "\n";
            f << ")\n";
        }
    }
    {
        auto f = foam_open(poly / "boundary");
        foam_write_header(f, "polyBoundaryMesh", "boundary", rOpts);
        f << order.mPatches.size() << "\n(\n";
        for (const FoamPatchOut& p : order.mPatches) {
            f << "    " << p.mName << "\n    {\n";
            f << "        type            " << p.mType << ";\n";
            f << "        nFaces          " << p.mNFaces << ";\n";
            f << "        startFace       " << p.mStartFace << ";\n";
            f << "    }\n";
        }
        f << ")\n";
    }

    // ---- zones from named regions (cellZones/faceZones/pointZones) ----
    {
        std::unordered_map<std::int64_t, std::int64_t> global_to_compact;
        for (std::size_t c = 0; c < faces.mCellToGlobal.size(); ++c)
            global_to_compact[faces.mCellToGlobal[c]] = static_cast<std::int64_t>(c);
        std::vector<std::int64_t> old_to_new(faces.NumFaces());
        for (std::size_t i = 0; i < order.mNewToOld.size(); ++i)
            old_to_new[static_cast<std::size_t>(order.mNewToOld[i])] = static_cast<std::int64_t>(i);

        auto write_or_remove = [&](const char* pFile, const std::vector<FoamZoneOut>& rZones,
                                   const char* pClass, const char* pZoneType,
                                   const char* pLabelKey) {
            if (rZones.empty()) {
                std::error_code rc;
                if (fs::remove(poly / pFile, rc))
                    log::info("OpenFOAM: removed stale {}", pFile);
                return;
            }
            foam_write_zone_file(poly / pFile, pClass, pFile, pZoneType, pLabelKey, rZones, rOpts);
        };

        write_or_remove("cellZones", foam_collect_cell_zones(rMesh, global_to_compact),
                        "cellZoneList", "cellZone", "cellLabels");
        write_or_remove("pointZones", foam_collect_point_zones(rMesh), "pointZoneList", "pointZone",
                        "pointLabels");
        write_or_remove("faceZones",
                        foam_collect_face_zones(rMesh, faces, global_to_compact, old_to_new),
                        "faceZoneList", "faceZone", "faceLabels");
    }

    foam_write_fields(foam_case_root(poly, rTarget.mRegion), rMesh, faces, order, assign, rOpts,
                      rTarget.mTime, rTarget.mRegion);

    log::info("Wrote polyMesh to {} ({} cells, {} faces, {} internal, {} patches)", poly.string(),
              faces.NumCells(), order.mNewToOld.size(), order.mNumInternal, order.mPatches.size());
}

/// A time directory's name for a time value: the shortest %.15g form, the way OpenFOAM
/// names them (`0.5`, `1e-05`).
std::string foam_time_name(double Time) {
    const detail::CNumber num;
    std::string name;
    num.Append(name, "%.15g", Time);
    return name;
}

}  // namespace

void write_openfoam(const std::string& rPath, const Mesh& rMesh, const OpenFoamInfo& rInfo,
                    const OpenFoamWriteOptions& rOpts) {
    foam_write_impl(rPath, rMesh, rInfo, rOpts, FoamTarget{}, FoamWriteMode::Full);
}

void write_openfoam(const std::string& rPath, const Mesh& rMesh, const OpenFoamInfo& rInfo) {
    write_openfoam(rPath, rMesh, rInfo, OpenFoamWriteOptions{});
}

void write_openfoam_region(const std::string& rPath, const Mesh& rMesh, const OpenFoamInfo& rInfo,
                           const OpenFoamWriteOptions& rOpts, const std::string& rRegion) {
    FoamTarget target;
    target.mRegion = rRegion;
    foam_write_impl(rPath, rMesh, rInfo, rOpts, target, FoamWriteMode::Full);
}

void write_openfoam_decomposed(const std::string& rPath, const Mesh& rMesh,
                               const OpenFoamInfo& rInfo, const OpenFoamWriteOptions& rOptions) {
    if (!rMesh.HasCellData("partition:part") ||
        rMesh.CellDataNumBlocks("partition:part") != rMesh.NumCellBlocks())
        throw WriteError(
            "OpenFOAM: a decomposed case needs cell_data[\"partition:part\"] on every cell "
            "block (meshioplusplus.partition_labels)");
    const detail::GlobalFaces faces = detail::build_global_faces(rMesh);
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    std::vector<std::int64_t> parts(faces.NumCells());
    for (std::size_t c = 0; c < parts.size(); ++c) {
        const auto [block, row] = detail::global_to_block_row(bases, faces.mCellToGlobal[c]);
        parts[c] = detail::read_int(rMesh.CellData("partition:part", block), row);
    }
    FoamTarget target;
    target.mParts = &parts;
    foam_write_impl(rPath, rMesh, rInfo, rOptions, target, FoamWriteMode::Decomposed);
}

namespace {

/// `constant/regionProperties` of the case @p rPath names (a `.foam` file or a case directory).
fs::path foam_region_properties_path(const std::string& rPath) {
    const fs::path given(rPath);
    const fs::path root = given.extension() == ".foam" ? given.parent_path() : given;
    return root / "constant" / "regionProperties";
}

/// A name that can sit in a `regionProperties` list unquoted: one OpenFOAM word.
bool foam_region_properties_word_ok(const std::string& rWord) {
    if (rWord.empty())
        return false;
    for (const char c : rWord)
        if (static_cast<unsigned char>(c) <= ' ' || static_cast<unsigned char>(c) >= 127 ||
            std::strchr("(){};\"'#$/\\", c) != nullptr)
            return false;
    return true;
}

/// The words of @p rText, with `(`, `)` and `;` as words of their own.
std::vector<std::string> foam_region_properties_tokens(const std::string& rText) {
    std::vector<std::string> tokens;
    std::string word;
    const auto flush = [&] {
        if (!word.empty())
            tokens.push_back(std::move(word));
        word.clear();
    };
    for (const char c : rText) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            flush();
        } else if (c == '(' || c == ')' || c == ';') {
            flush();
            tokens.emplace_back(1, c);
        } else {
            word.push_back(c);
        }
    }
    flush();
    return tokens;
}

}  // namespace

void write_openfoam_region_properties(const std::string& rPath,
                                      const std::vector<OpenFoamRegionGroup>& rGroups) {
    if (rGroups.empty())
        throw WriteError("OpenFOAM: regionProperties needs at least one family of regions");
    if (fs::path(rPath).filename() == "polyMesh")
        throw WriteError(
            "OpenFOAM: regionProperties belongs to a case (a .foam file or a case directory), not "
            "a polyMesh directory");
    std::set<std::string> families, regions;
    for (const OpenFoamRegionGroup& group : rGroups) {
        if (!foam_region_properties_word_ok(group.mFamily))
            throw WriteError("OpenFOAM: invalid regionProperties family name '" + group.mFamily +
                             "'");
        if (!families.insert(group.mFamily).second)
            throw WriteError("OpenFOAM: regionProperties lists the family '" + group.mFamily +
                             "' twice");
        if (group.mRegions.empty())
            throw WriteError("OpenFOAM: regionProperties family '" + group.mFamily +
                             "' has no region");
        for (const std::string& region : group.mRegions) {
            if (!foam_region_properties_word_ok(region) || region[0] == '.')
                throw WriteError("OpenFOAM: invalid region name '" + region + "'");
            if (!regions.insert(region).second)
                throw WriteError("OpenFOAM: regionProperties lists the region '" + region +
                                 "' twice");
        }
    }
    const fs::path file = foam_region_properties_path(rPath);
    std::error_code ec;
    fs::create_directories(file.parent_path(), ec);
    if (ec && !fs::is_directory(file.parent_path()))
        throw WriteError("OpenFOAM: could not create directory " + file.parent_path().string() +
                         ": " + ec.message());
    auto f = foam_open(file);
    foam_write_header(f, "dictionary", "regionProperties", OpenFoamWriteOptions{}, "constant");
    f << "regions\n(\n";
    for (const OpenFoamRegionGroup& group : rGroups) {
        f << "    " << group.mFamily << " (";
        for (std::size_t i = 0; i < group.mRegions.size(); ++i)
            f << (i ? " " : "") << group.mRegions[i];
        f << ")\n";
    }
    f << ");\n";
}

std::vector<OpenFoamRegionGroup> read_openfoam_region_properties(const std::string& rPath) {
    const fs::path file = foam_region_properties_path(rPath);
    if (!fs::exists(file))
        throw ReadError("OpenFOAM: no " + file.string());
    const std::vector<std::string> tokens =
        foam_region_properties_tokens(strip_comments_and_header(read_whole(file.string()).View()));
    const auto bad = [&](const std::string& rWhy) {
        return ReadError("OpenFOAM: " + file.string() + ": " + rWhy);
    };
    std::size_t i = 0;
    while (i < tokens.size() && tokens[i] != "regions")
        ++i;
    if (i + 1 >= tokens.size() || tokens[i + 1] != "(")
        throw bad("no `regions ( family (names) ... )` list");
    std::vector<OpenFoamRegionGroup> groups;
    for (i += 2;;) {
        if (i >= tokens.size())
            throw bad("the `regions` list is not closed");
        if (tokens[i] == ")")
            break;
        OpenFoamRegionGroup group;
        group.mFamily = tokens[i++];
        if (group.mFamily == "(" || group.mFamily == ";" || i >= tokens.size() || tokens[i] != "(")
            throw bad("expected `family (regions)`, found '" + group.mFamily + "'");
        for (++i;; ++i) {
            if (i >= tokens.size() || tokens[i] == "(" || tokens[i] == ";")
                throw bad("the region list of '" + group.mFamily + "' is not closed");
            if (tokens[i] == ")") {
                ++i;
                break;
            }
            group.mRegions.push_back(tokens[i]);
        }
        groups.push_back(std::move(group));
    }
    return groups;
}

struct OpenFoamSeriesWriter::Impl {
    std::string mPath;
    OpenFoamInfo mInfo;
    OpenFoamWriteOptions mOptions;
    std::string mRegion;
    std::size_t mSteps = 0;
    std::size_t mPoints = 0;
    std::vector<std::pair<std::string, std::size_t>> mBlocks;  // type, cells
    std::vector<std::string> mTimes;
};

OpenFoamSeriesWriter::OpenFoamSeriesWriter(const std::string& rPath, const OpenFoamInfo& rInfo,
                                           const OpenFoamWriteOptions& rOptions,
                                           const std::string& rRegion)
    : mpImpl(std::make_unique<Impl>()) {
    mpImpl->mPath = rPath;
    mpImpl->mInfo = rInfo;
    mpImpl->mOptions = rOptions;
    mpImpl->mRegion = rRegion;
}

OpenFoamSeriesWriter::~OpenFoamSeriesWriter() = default;
OpenFoamSeriesWriter::OpenFoamSeriesWriter(OpenFoamSeriesWriter&&) noexcept = default;
OpenFoamSeriesWriter& OpenFoamSeriesWriter::operator=(OpenFoamSeriesWriter&&) noexcept = default;

void OpenFoamSeriesWriter::Write(double Time, const Mesh& rMesh) {
    if (!mpImpl)
        throw WriteError("OpenFOAM: series: the writer was moved from");
    Impl& s = *mpImpl;
    std::vector<std::pair<std::string, std::size_t>> blocks;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        blocks.emplace_back(std::string(rMesh.Cells(b).Type()), rMesh.Cells(b).NumCells());
    const std::string name = foam_time_name(Time);
    if (std::find(s.mTimes.begin(), s.mTimes.end(), name) != s.mTimes.end())
        throw WriteError("OpenFOAM: series: two steps share the time directory '" + name + "'");
    FoamTarget target;
    target.mRegion = s.mRegion;
    target.mTime = name;
    if (s.mSteps == 0) {
        s.mPoints = rMesh.NumPoints();
        s.mBlocks = blocks;
        foam_write_impl(s.mPath, rMesh, s.mInfo, s.mOptions, target, FoamWriteMode::Full);
    } else {
        if (rMesh.NumPoints() != s.mPoints || blocks != s.mBlocks)
            throw WriteError("OpenFOAM: series: step " + std::to_string(s.mSteps) +
                             "'s points or cells differ from the first step's; a case holds one "
                             "polyMesh -- write one case per step with '{step}'");
        foam_write_impl(s.mPath, rMesh, s.mInfo, s.mOptions, target, FoamWriteMode::FieldsOnly);
    }
    s.mTimes.push_back(name);
    ++s.mSteps;
}

std::size_t OpenFoamSeriesWriter::NumSteps() const noexcept {
    return mpImpl ? mpImpl->mSteps : 0;
}

void OpenFoamSeriesWriter::Finalize() {
    if (mpImpl && mpImpl->mSteps == 0)
        throw WriteError("OpenFOAM: series: no step was written");
}

}  // namespace meshioplusplus
