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
 * @file detail/drawing.hpp
 * @brief What the `svg` and `tikz` drawings share: the marker comment their
 * writers emit and their readers require, and the builder that turns the
 * drawing's paths back into a planar mesh.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed header
 * names it, and it adds nothing to the API or the ABI. Its Python twin is
 * `meshioplusplus/_drawing.py`.
 *
 * ### The marker
 *
 * Both drawings carry a comment right after the provenance tag
 * (`Written by meshio++ v...`): one line `Drawing: planar` or
 * `Drawing: projected azimuth=A elevation=E roll=R`, and for `svg`
 * `scale=S offset=OX,OY` after it. Numbers are `%.17g`, so the reader restores
 * them exactly. The reader refuses a file without the line unless it is told
 * `ReadOptions::mLenient` (`strict=False`), because a drawing is not a mesh: only
 * the structure meshio++'s own writers emit is read back, and the marker is how
 * a reader knows it is looking at one.
 *
 * ### The model, as the reader restores it
 *
 *     model_x = svg_x / S + OX      model_y = OY - svg_y / S      (svg)
 *
 * with `S` the scale `image_width` applied (1 when none was) and `OY` the
 * `max_y + min_y` the writer flipped y about. Without a marker `S = 1`, `OX = 0`
 * and `OY = 2 * viewBox_y + viewBox_height`, which undoes the flip in drawing
 * units.
 */

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/mesh.hpp"
#include "text_cursor.hpp"

namespace meshioplusplus {
namespace detail {

/// The prefix of the marker line, shared by both drawings.
inline constexpr const char* kDrawingMarkerPrefix = "Drawing: ";

/// File text for an error message: anything but printable ASCII becomes '?', so a
/// message built from hostile bytes is always valid UTF-8 where it is shown.
inline std::string drawing_safe(std::string_view Text) {
    std::string out(Text);
    for (char& c : out)
        if (static_cast<unsigned char>(c) < 0x20 || static_cast<unsigned char>(c) > 0x7e)
            c = '?';
    return out;
}

/// `%.17g` -- what the marker prints numbers with, so the reader restores them.
inline std::string drawing_num17(double Value) {
    char buf[40];
    snprintf_c(buf, sizeof(buf), "%.17g", Value);
    return buf;
}

/// `%g` -- what the read warning prints numbers with.
inline std::string drawing_numg(double Value) {
    char buf[40];
    snprintf_c(buf, sizeof(buf), "%g", Value);
    return buf;
}

/// The marker's `Drawing: ...` line. @p pScaleOffset is null for tikz.
inline std::string drawing_marker_line(bool Projected, double Azimuth, double Elevation,
                                       double Roll, const double* pScaleOffset) {
    std::string line = kDrawingMarkerPrefix;
    if (Projected) {
        line += "projected azimuth=" + drawing_num17(Azimuth) +
                " elevation=" + drawing_num17(Elevation) + " roll=" + drawing_num17(Roll);
    } else {
        line += "planar";
    }
    if (pScaleOffset != nullptr) {
        line += " scale=" + drawing_num17(pScaleOffset[0]) +
                " offset=" + drawing_num17(pScaleOffset[1]) + "," + drawing_num17(pScaleOffset[2]);
    }
    return line;
}

/// The provenance lines of a drawing's comment, the marker line right after the
/// tag (so it is there under every provenance mode).
inline std::vector<std::string> drawing_comment_lines(const std::string& rMarkerLine) {
    std::vector<std::string> lines = provenance_lines(SlotTier::Block);
    lines.insert(lines.begin() + 1, rMarkerLine);
    return lines;
}

/// What a drawing's marker said.
struct DrawingMarker {
    bool mFound = false;
    bool mProjected = false;
    double mAzimuth = 0.0;
    double mElevation = 0.0;
    double mRoll = 0.0;
    double mScale = 1.0;
    double mOffsetX = 0.0;
    double mOffsetY = 0.0;
};

/**
 * @brief Scans an SVG/TikZ number at @p Pos in @p rText: `[+-]? (digits [. digits?]
 * | . digits) ([eE] [+-]? digits)?`.
 *
 * @return the token's length, or 0 when no number starts there.
 */
inline std::size_t drawing_scan_number(std::string_view Text, std::size_t Pos) {
    auto digit = [&](std::size_t i) { return i < Text.size() && Text[i] >= '0' && Text[i] <= '9'; };
    std::size_t i = Pos;
    if (i < Text.size() && (Text[i] == '+' || Text[i] == '-'))
        ++i;
    if (digit(i)) {
        while (digit(i))
            ++i;
        if (i < Text.size() && Text[i] == '.') {
            ++i;
            while (digit(i))
                ++i;
        }
    } else if (i < Text.size() && Text[i] == '.' && digit(i + 1)) {
        ++i;
        while (digit(i))
            ++i;
    } else {
        return 0;
    }
    if (i < Text.size() && (Text[i] == 'e' || Text[i] == 'E')) {
        std::size_t j = i + 1;
        if (j < Text.size() && (Text[j] == '+' || Text[j] == '-'))
            ++j;
        if (digit(j)) {
            while (digit(j))
                ++j;
            i = j;
        }
    }
    return i - Pos;
}

namespace drawing_detail {

inline std::string_view trim(std::string_view Text) {
    while (!Text.empty() && text_is_blank(Text.front()))
        Text.remove_prefix(1);
    while (!Text.empty() && text_is_blank(Text.back()))
        Text.remove_suffix(1);
    return Text;
}

inline double marker_number(std::string_view Token, const char* pFormat, const char* pKey) {
    double v = 0.0;
    if (drawing_scan_number(Token, 0) != Token.size() || !parse_double_token(Token, v) ||
        !std::isfinite(v))
        throw ReadError(std::string(pFormat) + ": malformed drawing marker, " + pKey +
                        " is not a number");
    return v;
}

}  // namespace drawing_detail

/**
 * @brief The marker among @p rCommentLines (the comment text of a drawing's
 * head, one entry per line).
 *
 * Found when both the provenance tag and a `Drawing: ` line are present. A
 * `Drawing: ` line that is malformed is a `ReadError` (the file claims to be
 * ours and is not); keys this version does not know are skipped.
 *
 * @param NeedScale whether the format requires `scale=`/`offset=` (svg).
 */
template <class TLines>
DrawingMarker parse_drawing_marker(const TLines& rCommentLines, const char* pFormat,
                                   bool NeedScale) {
    bool tag = false;
    std::string_view drawing;
    bool have_drawing = false;
    for (const auto& raw : rCommentLines) {
        const std::string_view line = drawing_detail::trim(std::string_view(raw));
        if (line.rfind("Written by meshio++ v", 0) == 0)
            tag = true;
        else if (!have_drawing && line.rfind(kDrawingMarkerPrefix, 0) == 0) {
            drawing = line;
            have_drawing = true;
        }
    }
    DrawingMarker marker;
    if (!tag || !have_drawing)
        return marker;
    drawing.remove_prefix(std::string_view(kDrawingMarkerPrefix).size());
    const std::vector<std::string_view> tokens = split_blanks(drawing);
    if (tokens.empty())
        throw ReadError(std::string(pFormat) + ": malformed drawing marker, empty");
    if (tokens[0] == "planar")
        marker.mProjected = false;
    else if (tokens[0] == "projected")
        marker.mProjected = true;
    else
        throw ReadError(std::string(pFormat) + ": malformed drawing marker, unknown kind '" +
                        drawing_safe(tokens[0]) + "'");
    bool has_az = false, has_el = false, has_roll = false, has_scale = false, has_offset = false;
    for (std::size_t i = 1; i < tokens.size(); ++i) {
        const std::size_t eq = tokens[i].find('=');
        if (eq == std::string_view::npos)
            continue;
        const std::string_view key = tokens[i].substr(0, eq);
        const std::string_view value = tokens[i].substr(eq + 1);
        if (key == "azimuth") {
            marker.mAzimuth = drawing_detail::marker_number(value, pFormat, "azimuth");
            has_az = true;
        } else if (key == "elevation") {
            marker.mElevation = drawing_detail::marker_number(value, pFormat, "elevation");
            has_el = true;
        } else if (key == "roll") {
            marker.mRoll = drawing_detail::marker_number(value, pFormat, "roll");
            has_roll = true;
        } else if (key == "scale") {
            marker.mScale = drawing_detail::marker_number(value, pFormat, "scale");
            has_scale = true;
        } else if (key == "offset") {
            const std::size_t comma = value.find(',');
            if (comma == std::string_view::npos)
                throw ReadError(std::string(pFormat) +
                                ": malformed drawing marker, offset needs two numbers");
            marker.mOffsetX =
                drawing_detail::marker_number(value.substr(0, comma), pFormat, "offset");
            marker.mOffsetY =
                drawing_detail::marker_number(value.substr(comma + 1), pFormat, "offset");
            has_offset = true;
        }
    }
    if (marker.mProjected && !(has_az && has_el && has_roll))
        throw ReadError(std::string(pFormat) +
                        ": malformed drawing marker, a projection needs azimuth, elevation, roll");
    if (NeedScale && !(has_scale && has_offset))
        throw ReadError(std::string(pFormat) +
                        ": malformed drawing marker, scale and offset are missing");
    if (has_scale && !(marker.mScale > 0.0))
        throw ReadError(std::string(pFormat) + ": malformed drawing marker, scale must be > 0");
    marker.mFound = true;
    return marker;
}

/// The one warning a drawing read raises: what the file was, and what is lost.
inline void drawing_read_warning(const char* pFormat, const DrawingMarker& rMarker) {
    if (!rMarker.mFound) {
        log::warn(
            "{}: no meshio++ drawing marker; read as a plain drawing in its own units "
            "(depth, scale and offset are unknown), z = 0 assumed",
            pFormat);
    } else if (rMarker.mProjected) {
        log::warn(
            "{}: projected from 3-D (azimuth {}, elevation {}, roll {}), depth lost; z = 0 "
            "assumed",
            pFormat, drawing_numg(rMarker.mAzimuth), drawing_numg(rMarker.mElevation),
            drawing_numg(rMarker.mRoll));
    } else {
        log::warn("{}: planar drawing, z = 0 assumed", pFormat);
    }
}

/// The refusal for a drawing without a marker.
inline ReadError drawing_unmarked_error(const char* pFormat) {
    return ReadError(
        std::string(pFormat) +
        ": no meshio++ drawing marker. Only drawings written by meshio++ are read back (a "
        "drawing is not a mesh); pass strict=False (ReadOptions::mLenient) to read this file "
        "as a plain drawing anyway");
}

/// An `#rgb` or `#rrggbb` colour, or nullopt.
inline std::optional<std::array<int, 3>> drawing_hex_color(std::string_view Text) {
    Text = drawing_detail::trim(Text);
    if (Text.size() != 4 && Text.size() != 7)
        return std::nullopt;
    if (Text[0] != '#')
        return std::nullopt;
    auto hex = [](char c) -> int {
        if (c >= '0' && c <= '9')
            return c - '0';
        if (c >= 'a' && c <= 'f')
            return c - 'a' + 10;
        if (c >= 'A' && c <= 'F')
            return c - 'A' + 10;
        return -1;
    };
    std::array<int, 3> rgb{};
    for (std::size_t k = 0; k < 3; ++k) {
        int hi, lo;
        if (Text.size() == 4) {
            hi = lo = hex(Text[1 + k]);
        } else {
            hi = hex(Text[1 + 2 * k]);
            lo = hex(Text[2 + 2 * k]);
        }
        if (hi < 0 || lo < 0)
            return std::nullopt;
        rgb[k] = hi * 16 + lo;
    }
    return rgb;
}

/**
 * @brief Collects a drawing's paths and builds the planar mesh.
 *
 * Points are welded on the exact text of their coordinates, in order of first
 * appearance; consecutive cells of one type and size form one block. A fill is
 * kept as an `int32` `(n, 3)` cell array named by the caller, `-1 -1 -1` for a
 * cell without an RGB fill, and only when some cell has one.
 */
class DrawingBuilder {
public:
    DrawingBuilder(const char* pFormat, const char* pFillName)
        : mpFormat(pFormat), mpFillName(pFillName) {}

    /// Whether the fill array is kept (a `ReadOptions` narrowing it away says no).
    void SetWantFill(bool Want) { mWantFill = Want; }

    /// The welded index of the vertex whose coordinate texts are @p rX/@p rY.
    /// @p rParse returns the `{x, y}` values and is called for a new vertex only.
    template <class TParse>
    std::int64_t Weld(std::string_view X, std::string_view Y, TParse&& rParse) {
        std::string key;
        key.reserve(X.size() + Y.size() + 1);
        key.append(X);
        key.push_back(' ');
        key.append(Y);
        const auto it = mIndex.find(key);
        if (it != mIndex.end())
            return it->second;
        const std::pair<double, double> xy = rParse();
        const auto id = static_cast<std::int64_t>(mX.size());
        mX.push_back(xy.first);
        mY.push_back(xy.second);
        mIndex.emplace(std::move(key), id);
        return id;
    }

    /// One cell. @p Closed is a closed path (a face), otherwise a line.
    void AddCell(const std::vector<std::int64_t>& rIds, bool Closed,
                 const std::optional<std::array<int, 3>>& rFill) {
        const std::size_t n = rIds.size();
        const char* type = nullptr;
        if (!Closed) {
            if (n != 2)
                throw ReadError(std::string(mpFormat) + ": an open path of " + std::to_string(n) +
                                " vertices is not a line (only 2-vertex open paths are read)");
            type = "line";
        } else if (n == 3) {
            type = "triangle";
        } else if (n == 4) {
            type = "quad";
        } else if (n > 4) {
            type = "polygon";
        } else {
            throw ReadError(std::string(mpFormat) + ": a closed path of " + std::to_string(n) +
                            " vertices is not a face");
        }
        if (mBlocks.empty() || mBlocks.back().mType != type || mBlocks.back().mWidth != n) {
            mBlocks.push_back(Block{type, n, {}, {}});
        }
        Block& block = mBlocks.back();
        block.mConn.insert(block.mConn.end(), rIds.begin(), rIds.end());
        if (rFill && mWantFill) {
            mAnyFill = true;
            block.mFill.insert(block.mFill.end(), rFill->begin(), rFill->end());
        } else {
            block.mFill.insert(block.mFill.end(), {-1, -1, -1});
        }
    }

    /**
     * @brief The mesh, points mapped back to model units.
     *
     * With @p Flip: `x / Scale + OffsetX`, `OffsetY - y / Scale`. Without it and
     * with the identity scale/offset, the coordinates are left bit-exact.
     */
    Mesh Finish(bool Flip, double Scale, double OffsetX, double OffsetY) {
        const std::size_t n = mX.size();
        NDArray pts(DType::Float64, {n, static_cast<std::size_t>(3)});
        double* pp = pts.As<double>();
        const bool identity = !Flip && Scale == 1.0 && OffsetX == 0.0 && OffsetY == 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double x = mX[i], y = mY[i];
            if (!identity) {
                x = x / Scale + OffsetX;
                y = Flip ? OffsetY - y / Scale : y / Scale + OffsetY;
            }
            if (!std::isfinite(x) || !std::isfinite(y))
                throw ReadError(std::string(mpFormat) + ": a coordinate is not finite");
            pp[i * 3 + 0] = x;
            pp[i * 3 + 1] = y;
            pp[i * 3 + 2] = 0.0;
        }
        Mesh mesh;
        mesh.AssignPoints(std::move(pts));
        std::vector<NDArray> fills;
        for (const Block& block : mBlocks) {
            const std::size_t cells = block.mConn.size() / block.mWidth;
            NDArray conn(DType::Int64, {cells, block.mWidth});
            std::copy(block.mConn.begin(), block.mConn.end(), conn.As<std::int64_t>());
            mesh.AddCellBlock(block.mType, std::move(conn));
            if (mAnyFill) {
                NDArray fill(DType::Int32, {cells, static_cast<std::size_t>(3)});
                std::copy(block.mFill.begin(), block.mFill.end(), fill.As<std::int32_t>());
                fills.push_back(std::move(fill));
            }
        }
        if (mAnyFill)
            mesh.AddCellData(mpFillName, std::move(fills));
        return mesh;
    }

private:
    struct Block {
        const char* mType;
        std::size_t mWidth;
        std::vector<std::int64_t> mConn;
        std::vector<std::int32_t> mFill;
    };

    const char* mpFormat;
    const char* mpFillName;
    std::unordered_map<std::string, std::int64_t> mIndex;
    std::vector<double> mX, mY;
    std::vector<Block> mBlocks;
    bool mAnyFill = false;
    bool mWantFill = true;
};

/// A scanned number's value, or `ReadError`.
inline double drawing_number_value(std::string_view Token, const char* pFormat) {
    double v = 0.0;
    if (!parse_double_token(Token, v) || !std::isfinite(v))
        throw ReadError(std::string(pFormat) + ": bad number '" + drawing_safe(Token) + "'");
    return v;
}

}  // namespace detail
}  // namespace meshioplusplus
